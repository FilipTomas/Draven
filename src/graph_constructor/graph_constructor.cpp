
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <cassert>
#include "graph.hpp"
#include "graph_constructor.h"
#include "parser.h"
#include "overlap.h"
#include "overlap_parser.h"
#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "biosoup/timer.hpp"
#include "extended_overlap.h"
#include "overlap_helpers.h"
#include "annotation_helpers.h"
#include "ram_overlap.h"
#include "overlap_selection.h"
#include "catboost_model.hpp"
#include <ATen/Parallel.h>
#include <torch/script.h>
#include <atomic>
#include <mutex>
#include <iostream>
#include <numeric>
#include "cnn_preprocessing.hpp"
// extern "C" {
//   #include "/home/ftomas/github/draven/build/_deps/wfa2-src/wavefront/wavefront_align.h"
// }
#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;
namespace raven {

edlib_align wfa_wrapper(const std::string& query, const std::string& target) {
    edlib_align result;

    thread_local wavefront_aligner_t* wf_aligner = []() {
        wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;

        attr.distance_metric = gap_affine;
        attr.affine_penalties.match = 0;
        attr.affine_penalties.mismatch = 4;
        attr.affine_penalties.gap_opening = 6;
        attr.affine_penalties.gap_extension = 2;

        attr.alignment_scope = compute_alignment;

        // Usually fastest if memory remains acceptable.
        attr.memory_mode = wavefront_memory_med;

        attr.heuristic.strategy = wf_heuristic_wfadaptive;
        attr.heuristic.min_wavefront_length = 10;
        attr.heuristic.max_distance_threshold = 100;
        attr.heuristic.steps_between_cutoffs = 5;

        return wavefront_aligner_new(&attr);
    }();

    wavefront_align(wf_aligner,
                    query.c_str(), query.size(),
                    target.c_str(), target.size());

    cigar_t* cigar = wf_aligner->cigar;
    std::string cigar_result;
    int count = 0;
    char last_op = 0;
    result.matches = 0;
    result.block_length = 0;

    for (int i = cigar->begin_offset; i < cigar->end_offset; i++) {
        char op = cigar->operations[i];
        char std_op;
        switch (op) {
            case 'M': case 'X': std_op = 'M'; break;
            case 'I': std_op = 'I'; break;
            case 'D': std_op = 'D'; break;
            default:  std_op = 'M'; break;
        }
        if (std_op == last_op) {
            count++;
        } else {
            if (last_op != 0) cigar_result += std::to_string(count) + last_op;
            last_op = std_op;
            count = 1;
        }
        if (std_op == 'M') result.matches++;
        result.block_length++;
    }
    if (last_op != 0) cigar_result += std::to_string(count) + last_op;

    result.cigar = cigar_result;
    result.edit_distance = cigar->score;

    // no wavefront_aligner_delete — reused across calls on same thread
    return result;
}

static inline ram_overlap
start_group(const biosoup::Overlap& ov) {
    extended_overlap eo = {ov, {}, 0, 0};
    ram_overlap ro;
    ro.lhs_id = ov.lhs_id;          // adapt to your field names
    ro.rhs_id = ov.rhs_id;

    ro.lhs_begin = ov.lhs_begin;    // first fragment → starts here
    ro.lhs_end   = ov.lhs_end;      // will be updated as we append
    ro.rhs_begin = ov.rhs_begin;
    ro.rhs_end   = ov.rhs_end;
    ro.n_fragments = 1;

    ro.overlap_fragments.push_back(eo);
    return ro;
}

std::vector<biosoup::Overlap> remove_duplicates(const std::vector<biosoup::Overlap>& ovlps){
    std::vector<biosoup::Overlap> result;
    result.reserve(ovlps.size());

    for (const auto& ov : ovlps) {
        bool is_duplicate = false;
        for (const auto& res_ov : result) {
            if (ov.lhs_id == res_ov.lhs_id &&
                ov.rhs_id == res_ov.rhs_id) {
                is_duplicate = true;
                break;
            }
        }
        if (!is_duplicate) {
            result.push_back(ov);
        }
    }
    return result;
}

// Merge Map's overlap fragments into ONE overlap per pair. A pair's chain is
// often broken into several fragments by indel drift; keeping only the longest
// and extrapolating its diagonal to the read ends mis-places the endpoints once
// there are errors. Instead we span the OUTERMOST anchors of the best collinear
// cluster, so both endpoints are real matched minimizers with the errors in
// between spanned (not extrapolated). Fragments join a cluster only if they are
// same-strand and on a consistent diagonal (within DIAG_TOL) -- this follows
// indel drift but rejects repeat-induced off-diagonal fragments that would
// otherwise fuse into a bogus giant overlap.
std::vector<biosoup::Overlap>
merge_fragments(std::vector<biosoup::Overlap> ovlps) {
    if (ovlps.empty()) {
        return {};
    }

    // Diagonal along which a collinear chain runs: fwd overlaps keep rhs - lhs
    // ~constant, reverse overlaps keep rhs + lhs ~constant.
    auto diagonal = [](const biosoup::Overlap& o) -> std::int64_t {
        return o.strand
            ? static_cast<std::int64_t>(o.rhs_begin) - static_cast<std::int64_t>(o.lhs_begin)
            : static_cast<std::int64_t>(o.rhs_begin) + static_cast<std::int64_t>(o.lhs_begin);
    };
    constexpr std::int64_t DIAG_TOL = 2000;  // allow indel drift, reject repeats

    std::sort(ovlps.begin(), ovlps.end(),
        [](const biosoup::Overlap& a, const biosoup::Overlap& b) {
            if (a.lhs_id != b.lhs_id) return a.lhs_id < b.lhs_id;
            if (a.rhs_id != b.rhs_id) return a.rhs_id < b.rhs_id;
            if (a.strand != b.strand) return a.strand < b.strand;
            return a.lhs_begin < b.lhs_begin;
        });

    std::vector<biosoup::Overlap> merged;
    std::size_t g = 0;
    while (g < ovlps.size()) {
        std::size_t h = g;
        while (h < ovlps.size() &&
               ovlps[h].lhs_id == ovlps[g].lhs_id &&
               ovlps[h].rhs_id == ovlps[g].rhs_id) {
            ++h;
        }

        // Within the pair [g, h) (sorted by strand then lhs_begin), greedily
        // grow collinear same-strand clusters and keep the longest-query-span one.
        biosoup::Overlap best{};
        std::int64_t best_span = -1;

        biosoup::Overlap cur = ovlps[g];
        std::int64_t cur_diag = diagonal(ovlps[g]);
        std::uint64_t cur_score = ovlps[g].score;

        auto flush = [&]() {
            cur.score = static_cast<std::uint32_t>(
                std::min<std::uint64_t>(cur_score, 0xFFFFFFFFULL));
            const std::int64_t span = static_cast<std::int64_t>(cur.lhs_end) -
                                      static_cast<std::int64_t>(cur.lhs_begin);
            if (span > best_span) { best = cur; best_span = span; }
        };

        for (std::size_t i = g + 1; i < h; ++i) {
            std::int64_t dd = diagonal(ovlps[i]) - cur_diag;
            if (dd < 0) dd = -dd;
            const bool collinear =
                (ovlps[i].strand == cur.strand) && (dd <= DIAG_TOL);
            if (collinear) {
                cur.lhs_begin = std::min(cur.lhs_begin, ovlps[i].lhs_begin);
                cur.lhs_end   = std::max(cur.lhs_end,   ovlps[i].lhs_end);
                cur.rhs_begin = std::min(cur.rhs_begin, ovlps[i].rhs_begin);
                cur.rhs_end   = std::max(cur.rhs_end,   ovlps[i].rhs_end);
                cur_diag = diagonal(ovlps[i]);  // follow the drift
                cur_score += ovlps[i].score;
            } else {
                flush();
                cur = ovlps[i];
                cur_diag = diagonal(ovlps[i]);
                cur_score = ovlps[i].score;
            }
        }
        flush();
        merged.push_back(best);
        g = h;
    }
    return merged;
}

// Anchors for aligning `kept` (the one overlap kept for a read pair)
// piecewise. The pair's chain is often split into several fragments, so start
// from the longest fragment inside `kept` and walk outwards along lhs, taking
// fragments that continue its diagonal; the anchors then also cover the parts
// `kept` was extrapolated over. The union is reduced to a strictly collinear
// chain with anchors at least `spacing` apart.
//
// Every anchor pins the alignment path, so one exact k-mer match from another
// repeat copy (e.g. a HOR period away) forces a detour that full NW would not
// take. Diagonal shifts are therefore only allowed at the rate indel drift
// produces (~1/8 of the distance), and anchors off their neighbours' median
// diagonal are dropped.
std::vector<std::pair<std::uint32_t, std::uint32_t>> collect_anchors(
    const biosoup::Overlap &kept,
    const std::vector<biosoup::Overlap> &fragments,
    const std::vector<ram::MinimizerEngine::Anchors> &fragment_anchors,
    std::uint32_t spacing) {
    using Anchor = std::pair<std::uint32_t, std::uint32_t>;
    const bool fwd = kept.strand;
    // fwd chains keep rhs - lhs ~constant, reverse chains keep rhs + lhs
    auto diagonal = [&](const Anchor& a) -> std::int64_t {
        return fwd
            ? static_cast<std::int64_t>(a.second) - static_cast<std::int64_t>(a.first)
            : static_cast<std::int64_t>(a.second) + static_cast<std::int64_t>(a.first);
    };
    // largest diagonal shift plausible from indel drift over `distance` bp
    auto max_drift = [](std::int64_t distance) -> std::int64_t {
        return 32 + distance / 8;
    };
    // fragments are joined across ram's chain gaps, allow a little more
    auto continues = [&](const Anchor& from, const Anchor& to) -> bool {
        const std::int64_t distance = static_cast<std::int64_t>(to.first) - from.first;
        return std::abs(diagonal(to) - diagonal(from)) <= 32 + max_drift(distance);
    };
    constexpr std::size_t MEDIAN_RADIUS = 8;  // neighbours on each side
    constexpr std::int64_t MEDIAN_TOL = 20;   // bp off the local median diagonal

    std::vector<std::size_t> pair;  // fragments of this pair and strand
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        if (fragments[f].rhs_id == kept.rhs_id &&
            fragments[f].strand == kept.strand &&
            !fragment_anchors[f].empty()) {
            pair.push_back(f);
        }
    }
    std::sort(pair.begin(), pair.end(), [&](std::size_t a, std::size_t b) {
        return fragment_anchors[a].front().first < fragment_anchors[b].front().first;
    });

    std::size_t seed = pair.size();
    std::uint32_t seed_len = 0;
    for (std::size_t p = 0; p < pair.size(); ++p) {
        const auto& o = fragments[pair[p]];
        if (o.lhs_begin < kept.lhs_begin || o.lhs_end > kept.lhs_end ||
            o.rhs_begin < kept.rhs_begin || o.rhs_end > kept.rhs_end) {
            continue;
        }
        if (seed == pair.size() || o.lhs_end - o.lhs_begin > seed_len) {
            seed = p;
            seed_len = o.lhs_end - o.lhs_begin;
        }
    }
    if (seed == pair.size()) {
        return {};
    }

    std::vector<std::size_t> chain{pair[seed]};
    for (std::size_t p = seed + 1, last = pair[seed]; p < pair.size(); ++p) {
        const auto& a = fragment_anchors[pair[p]];
        const auto& b = fragment_anchors[last];
        if (a.front().first <= b.back().first || !continues(b.back(), a.front())) {
            continue;
        }
        chain.push_back(pair[p]);
        last = pair[p];
    }
    for (std::size_t p = seed, first = pair[seed]; p-- > 0;) {
        const auto& a = fragment_anchors[pair[p]];
        const auto& b = fragment_anchors[first];
        if (a.back().first >= b.front().first || !continues(a.back(), b.front())) {
            continue;
        }
        chain.push_back(pair[p]);
        first = pair[p];
    }
    std::sort(chain.begin(), chain.end(), [&](std::size_t a, std::size_t b) {
        return fragment_anchors[a].front().first < fragment_anchors[b].front().first;
    });

    std::vector<Anchor> dense;
    for (std::size_t f : chain) {
        dense.insert(dense.end(), fragment_anchors[f].begin(), fragment_anchors[f].end());
    }

    std::vector<std::int64_t> window;
    window.reserve(2 * MEDIAN_RADIUS + 1);
    std::vector<Anchor> dst;
    for (std::size_t n = 0; n < dense.size(); ++n) {
        const Anchor& a = dense[n];
        if (!dst.empty() && a.first < dst.back().first + spacing) {
            continue;
        }
        window.clear();
        for (std::size_t m = n > MEDIAN_RADIUS ? n - MEDIAN_RADIUS : 0;
             m < std::min(dense.size(), n + MEDIAN_RADIUS + 1); ++m) {
            window.push_back(diagonal(dense[m]));
        }
        std::nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
        if (std::abs(diagonal(a) - window[window.size() / 2]) > MEDIAN_TOL) {
            continue;
        }
        if (!dst.empty()) {
            const Anchor& b = dst.back();
            const std::int64_t dl = static_cast<std::int64_t>(a.first) - b.first;
            const std::int64_t dr = fwd
                ? static_cast<std::int64_t>(a.second) - b.second
                : static_cast<std::int64_t>(b.second) - a.second;
            // not collinear, or more indel imbalance than drift explains
            if (dr <= 0 || std::abs(dl - dr) > max_drift(dl)) {
                continue;
            }
        }
        dst.push_back(a);
    }
    return dst;
}

Graph_Constructor::Graph_Constructor(Graph &graph, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
  : graph_(graph), thread_pool_(thread_pool ?
                                thread_pool :
                                std::make_shared<thread_pool::ThreadPool>(1)) {
}

void Graph_Constructor::LoadAndCount(
    const std::string &sequence_path,
    Program_Parameters &param,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences) {
  raven::Parser parser;
  auto sparser = parser.CreateParser(sequence_path);
  if (sparser == nullptr) {
    return;
  }

  // Bloom sizing. Total distinct k-mers isn't known before parsing, so for the
  // auto case estimate an upper bound from the file size (bytes >= bases), which
  // safely over-sizes the bloom (extra bloom memory, lower false-positive rate).
  std::uint64_t bloom_hint;
  if (param.bloom_gb < 0.0) {
    std::uint64_t fsize = 0;
    try { fsize = std::filesystem::file_size(sequence_path); } catch (...) {}
    bloom_hint = fsize;  // ~1 bit per input byte after next-pow2 rounding
    std::cerr << "[raven::] internal counter bloom: auto (by file size)" << std::endl;
  } else if (param.bloom_gb == 0.0) {
    bloom_hint = 0;  // disable bloom -> exact count + prune
    std::cerr << "[raven::] internal counter bloom: disabled (exact count)" << std::endl;
  } else {
    bloom_hint = static_cast<std::uint64_t>(
        param.bloom_gb * static_cast<double>(1ULL << 33));  // GB -> bits
    std::cerr << "[raven::] internal counter bloom: " << param.bloom_gb << " GB" << std::endl;
  }

  internal_counter_ = std::make_unique<ram::KmerCounter>(
      param.kmer_len, thread_pool_, /*min_count=*/2, /*bloom_capacity_hint=*/bloom_hint);

  biosoup::Timer timer{};
  timer.Start();

  // Parse one chunk at a time (single-threaded decode); dispatch each chunk's
  // counting to the pool WITHOUT blocking, so the count overlaps parsing the
  // next chunk. Reads are appended to `sequences`; the NucleicAcid objects are
  // heap-allocated so the raw pointers handed to CountAsync stay valid across
  // the moves/reallocations.
  const std::uint64_t chunk_bytes = 1ULL << 28;  // 256 MB of input per chunk
  std::uint64_t total = 0;
  while (true) {
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> chunk;
    try {
      chunk = sparser->Parse(chunk_bytes);
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
      break;
    }
    if (chunk.empty()) {
      break;
    }
    for (const auto &read : chunk) {
      internal_counter_->CountAsync(read.get());
    }
    for (auto &read : chunk) {
      sequences.emplace_back(std::move(read));
    }
    total += chunk.size();
  }

  internal_counter_->Finalize();
  std::cerr << "[raven::] loaded + counted " << total << " sequences ("
            << internal_counter_->num_distinct() << " distinct k-mers >= 2) "
            << std::fixed << timer.Stop() << "s" << std::endl;
}

void Graph_Constructor::Construct(
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
  Program_Parameters &param) {

  disagreement_ = param.disagreement;

  if (sequences.empty()) {
    return;
  }

  if (graph_.state_manager_.state() != GraphState::Construct_Graph) {
    return;
  }

  std::vector<std::vector<extended_overlap>> extended_overlaps(sequences.size());
  biosoup::Timer timer{};

  // checkpoint test
  if (graph_.use_checkpoints()) {
    graph_.Store(param.cereal_filename);
  }


  if(param.gt_overlaps.empty()){
    ConstructOverlaps(sequences, extended_overlaps, timer, param);
  } else {
    ConstructOverlapsFromGT(sequences, extended_overlaps, timer, param);
  };

  graph_.state_manager_.advance_state();
  if (graph_.use_checkpoints()) {
    timer.Start();
    graph_.Store(param.cereal_filename);
    std::cerr << "[raven::Graph::Construct] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }
  ConstructOnlyBackboneGraph(sequences, extended_overlaps, timer, param);
  ConstructBackboneGraph(sequences, extended_overlaps, timer, param);
  //exit(0);
  // ConstructOverlapGraph(sequences, extended_overlaps, timer, param);
  // ConstructAssemblyGraph(sequences, extended_overlaps, timer, param);
  graph_.PrintOverlaps(extended_overlaps, sequences, false, "after_backbone.paf");
  exit(0);
  graph_.state_manager_.advance_state();
  if (graph_.use_checkpoints()) {
    timer.Start();
    graph_.Store(param.cereal_filename);
    std::cerr << "[raven::Graph::Construct] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  std::cerr << "[raven::Graph::Construct] "
            << std::fixed << timer.elapsed_time() << "s"
            << std::endl;
}

// Called once all output files are closed. Ends the process while the overlap
// piles are still in scope, so neither they nor the reads, k-mer counter and
// minimizer index are destroyed: that teardown is single-threaded and took
// ~190 s of a 21 min chr18 run.
[[noreturn]] static void ExitWithoutTeardown() {
  std::cout.flush();
  std::cerr.flush();
  std::fflush(nullptr);
  std::_Exit(0);
}

// Log what the overlap piles hold: entries, and bytes of the structs (by
// capacity), their anchors and their CIGARs.
static void LogOverlapMemory(
    const std::vector<std::vector<extended_overlap>> &piles, const char *when) {
  std::uint64_t entries = 0, structs = 0, anchors = 0, cigars = 0;
  for (const auto &pile : piles) {
    entries += pile.size();
    structs += pile.capacity() * sizeof(extended_overlap);
    for (const auto &it : pile) {
      anchors += it.anchors.capacity() * sizeof(it.anchors[0]) +
                 it.gap_positions.capacity() * sizeof(it.gap_positions[0]);
      if (it.edlib_alignment.cigar.capacity() > 15) {  // past the in-place buffer
        cigars += it.edlib_alignment.cigar.capacity();
      }
    }
  }
  const double gb = 1024.0 * 1024.0 * 1024.0;
  std::cerr << std::fixed << std::setprecision(2)
            << "[raven::] overlaps " << when << ": " << entries << " entries, "
            << structs / gb << " GB structs + " << anchors / gb << " GB anchors + "
            << cigars / gb << " GB CIGARs" << std::endl;
}

void Graph_Constructor::ConstructCorrectionOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, Program_Parameters &param){
    std::vector<std::vector<extended_overlap>> extended_overlaps(sequences.size());
    biosoup::Timer timer{};
    // Overlaps (and their features) are computed here unless we only load an
    // externally selected PAF (LoadOverlapsAndAlign).
    const bool compute_overlaps = param.dump_overlaps_only ||
                                  param.dump_all_overlaps ||
                                  !param.rf_model_path.empty();

    if (compute_overlaps) {
      ConstructOverlapsForCorrection(sequences, extended_overlaps, timer, param);

      if (param.dump_all_overlaps) {
        // Raw overlaps were already dumped inside ConstructOverlapsForCorrection.
        ExitWithoutTeardown();
      }

      if (param.dump_overlaps_only) {
        // Debug: dump the correction overlaps (backbone, non strong-contained
        // reads) with full features BEFORE alignment, then stop -- skips the
        // costly alignment step.
        std::cerr << "Printing pre-alignment backbone correction overlaps with features!"
                  << std::endl;
        graph_.PrintBackboneOverlapsWithFeatures(
            extended_overlaps, sequences, "correction_overlaps_features.paf");
        ExitWithoutTeardown();
      }
    }

    // correction.paf is written pile by pile as piles finish aligning, and
    // each pile is freed with its CIGARs right away, instead of holding every
    // CIGAR until the end. Piles land in completion order, but each pile's
    // lines (all alignments of one target read) stay together and in order,
    // which is all HERRO reads per target.
    std::ofstream paf("correction.paf");
    if (!paf) {
      std::cerr << "[raven::] error: cannot open correction.paf" << std::endl;
      ExitWithoutTeardown();
    }
    std::mutex paf_mutex;
    std::uint64_t written = 0;
    std::uint64_t skipped = 0;
    auto write_pile = [&](std::uint32_t i) -> void {
      CorrectionPafText text =
          graph_.FormatCorrectionOverlaps(extended_overlaps[i], sequences);
      std::vector<extended_overlap>().swap(extended_overlaps[i]);
      std::lock_guard<std::mutex> lock(paf_mutex);
      paf << text.paf;
      std::cerr << text.log;
      written += text.written;
      skipped += text.skipped;
    };

    if (compute_overlaps) {
      LogOverlapMemory(extended_overlaps, "before selection");
      SelectCorrectionOverlaps(sequences, extended_overlaps, param);
      LogOverlapMemory(extended_overlaps, "after selection");
      AlignCorrectionOverlaps(sequences, extended_overlaps, param.kmer_len,
                          param.max_edit_rate, write_pile);
      LogOverlapMemory(extended_overlaps, "after alignment");
    } else {
      LoadOverlapsAndAlign(param, sequences, extended_overlaps, write_pile);
    }
    paf.close();
    //AlignBackboneReads(sequences, extended_overlaps, timer, param);
    //ConstructOnlyBackboneGraph(sequences, extended_overlaps, timer, param);
    std::cerr << "[raven::Graph::PrintCorrectionOverlaps] printed "
              << written << " overlaps, skipped " << skipped
              << " (written while aligning)" << std::endl;
    std::cerr << "[raven::Graph::Construct] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
    ExitWithoutTeardown();
}

// In-process replacement for mm2_draven_match.ipynb: score every backbone
// correction overlap with the exported random forest, keep the top-K per
// target window, and reset the kept overlaps to the state
// LoadOverlapsAndAlign would load them in (plain coordinates, not aligned).
void Graph_Constructor::SelectCorrectionOverlaps(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
    std::vector<std::vector<extended_overlap>> &extended_overlaps,
    Program_Parameters &param) {
  biosoup::Timer timer{};
  timer.Start();

  overlap_selection::RandomForest forest;
  forest.Load(param.rf_model_path);
  std::cerr << "[raven::SelectCorrectionOverlaps] loaded " << forest.num_trees()
            << " trees from " << param.rf_model_path << std::endl;

  // Same reads PrintBackboneOverlapsWithFeatures dumps for the notebook.
  auto is_corrected = [&](std::uint32_t i) -> bool {
    return i < graph_.piles_.size() && graph_.piles_[i]->is_backbone() &&
           !graph_.piles_[i]->is_strong_contained();
  };

  std::vector<std::vector<double>> probabilities(extended_overlaps.size());
  std::vector<std::vector<bool>> selected(extended_overlaps.size());

  auto score_pile = [&](std::uint32_t i) -> void {
    const auto &pile = extended_overlaps[i];
    if (!is_corrected(i) || pile.empty()) return;

    std::vector<std::uint32_t> query_lengths(pile.size());
    for (std::size_t k = 0; k < pile.size(); ++k) {
      query_lengths[k] = sequences[pile[k].overlap.rhs_id]->inflated_len;
    }

    const std::uint32_t target_length = sequences[i]->inflated_len;
    const auto features = overlap_selection::ComputeFeatures(
        pile, target_length, query_lengths);
    const std::size_t n_features = overlap_selection::FeatureNames().size();

    auto &probs = probabilities[i];
    probs.resize(pile.size());
    for (std::size_t k = 0; k < pile.size(); ++k) {
      probs[k] = forest.PredictProba(&features[k * n_features]);
    }

    selected[i] = overlap_selection::SelectWindowTopK(
        pile, probs, target_length, param.rf_window_size, param.rf_top_k,
        param.rf_require_full_window);
  };

  std::vector<std::future<void>> futures;
  futures.reserve(extended_overlaps.size());
  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    futures.emplace_back(thread_pool_->Submit(score_pile, i));
  }
  for (const auto &it : futures) {
    it.wait();
  }

  if (param.dump_rf_scores) {
    // PAF key columns (query = rhs, target = lhs) + probability, joinable
    // with the notebook's selected_top30_per_window.tsv.
    std::ofstream os("rf_selection.tsv");
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "qn\tql\tqs\tqe\tstrand\ttn\ttl\tts\tte\tpred_proba\tselected\n";
    for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
      for (std::size_t k = 0; k < probabilities[i].size(); ++k) {
        const auto &o = extended_overlaps[i][k].overlap;
        os << sequences[o.rhs_id]->name << '\t'
           << sequences[o.rhs_id]->inflated_len << '\t'
           << o.rhs_begin << '\t' << o.rhs_end << '\t'
           << (o.strand ? '+' : '-') << '\t'
           << sequences[o.lhs_id]->name << '\t'
           << sequences[o.lhs_id]->inflated_len << '\t'
           << o.lhs_begin << '\t' << o.lhs_end << '\t'
           << probabilities[i][k] << '\t'
           << (selected[i][k] ? 1 : 0) << '\n';
      }
    }
    std::cerr << "[raven::SelectCorrectionOverlaps] wrote rf_selection.tsv"
              << std::endl;
  }

  std::uint64_t candidates = 0;
  std::uint64_t kept = 0;
  std::uint64_t skipped = 0;
  // kept overlaps whose support read is strongly contained: in the target,
  // or overlapping it otherwise (dovetail, or containing it)
  std::uint64_t kept_contained_support = 0;
  std::uint64_t kept_dovetail_contained_support = 0;

  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    std::vector<extended_overlap> pile;
    pile.swap(extended_overlaps[i]);
    if (selected[i].empty()) continue;

    candidates += pile.size();
    for (std::size_t k = 0; k < pile.size(); ++k) {
      if (!selected[i][k]) continue;

      const auto &o = pile[k].overlap;
      // same validation as LoadOverlapsAndAlign
      if (o.lhs_end <= o.lhs_begin || o.rhs_end <= o.rhs_begin ||
          o.lhs_end > sequences[o.lhs_id]->inflated_len ||
          o.rhs_end > sequences[o.rhs_id]->inflated_len) {
        ++skipped;
        continue;
      }

      biosoup::Overlap overlap{o.lhs_id, o.lhs_begin, o.lhs_end,
                               o.rhs_id, o.rhs_begin, o.rhs_end,
                               255, o.strand};
      edlib_align tmp = {};
      extended_overlap total_ovlp{overlap, tmp, 0, 0};
      total_ovlp.correction_overlap = false;
      total_ovlp.anchors = std::move(pile[k].anchors);
      total_ovlp.lhs_begin_original = o.lhs_begin;
      total_ovlp.lhs_end_original = o.lhs_end;
      total_ovlp.rhs_begin_original = o.rhs_begin;
      total_ovlp.rhs_end_original = o.rhs_end;

      extended_overlaps[i].emplace_back(std::move(total_ovlp));
      ++kept;
      if (graph_.piles_[o.rhs_id]->is_strong_contained()) {
        const std::uint32_t type = overlap_type_extended(
            o, graph_, sequences[o.lhs_id]->inflated_len,
            sequences[o.rhs_id]->inflated_len);
        ++(type == 2 ? kept_contained_support : kept_dovetail_contained_support);
      }
    }
  }

  std::cerr << "[raven::SelectCorrectionOverlaps] kept " << kept << " / "
            << candidates << " candidate overlaps ("
            << (candidates ? 100.0 * kept / candidates : 0.0) << "%), "
            << "window " << param.rf_window_size << ", top " << param.rf_top_k
            << (param.rf_require_full_window ? ", full-window" : ", any-overlap")
            << ", skipped " << skipped << " invalid, "
            << std::fixed << timer.Stop() << "s" << std::endl;
  std::cerr << "[raven::SelectCorrectionOverlaps] kept overlaps with a strongly "
               "contained support read: " << kept_contained_support
            << " contained in the target, " << kept_dovetail_contained_support
            << " otherwise" << std::endl;
}

// void Graph_Constructor::LoadOverlapsAndAlign(Program_Parameters &param, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps){
//   std::ifstream file("/mnt/share1_Jabba/ftomas/centromere_fix/overlap_selection/selected_top30_per_window.paf");
//   if (!file.is_open()) {
//     throw std::runtime_error("Error opening file: " + param.input_paf_path);
//   }
//   std::cerr << "Got to here" << std::endl;
//   std::map<std::string, std::uint32_t> sequence_name_to_seq_id;
//   for (std::uint32_t i = 0; i < sequences.size(); ++i) {
//     sequence_name_to_seq_id[sequences[i]->name] = sequences[i]->id;
//   }

//   std::cerr << "[raven::Graph::LoadHerroSNPs] loading overlaps from: " << param.input_paf_path << std::endl;
//   std::string line;
//   while (std::getline(file, line)) {
//     std::istringstream iss(line);
//     std::string item;
//     std::vector<std::string> items;
//     std::uint32_t lhs_seq_id;
//     std::uint32_t rhs_seq_id;

//     while (std::getline(iss, item, '\t')) {
//       items.push_back(item);
//     };

//     // lhs_seq_id = get_read_id(items[0], sequences);
//     lhs_seq_id = sequence_name_to_seq_id[items[0]];
//     // rhs_seq_id = get_read_id(items[5], sequences);
//     rhs_seq_id = sequence_name_to_seq_id[items[5]];

//     if (lhs_seq_id == (std::uint32_t )-1 || rhs_seq_id == (std::uint32_t )-1) {
//       continue;
//     } else {
//       biosoup::Overlap overlap{ lhs_seq_id, (std::uint32_t)std::stoi(items[2]), (std::uint32_t)std::stoi(items[3]),
//                                 rhs_seq_id, (std::uint32_t)std::stoi(items[7]), (std::uint32_t)std::stoi(items[8]),
//                                 255, items[4] == "+" ? true : false};
//       edlib_align tmp = {};
//       extended_overlap total_ovlp{ overlap, tmp, 0, 0 };
//       total_ovlp.correction_overlap = true;
//       extended_overlaps[lhs_seq_id].emplace_back(total_ovlp);
//     }
//   }
//   std::cerr << "[raven::LoadAndAlign] loaded overlaps from: " << param.input_paf_path << std::endl;

//   auto resolve_alignments = [&](std::uint32_t i) -> void {
//     std::string lhs_full = sequences[i]->InflateData();
//     for (auto& overlap : extended_overlaps[i]) {
//       std::string lhs = lhs_full.substr(
//           overlap.overlap.lhs_begin,
//           overlap.overlap.lhs_end - overlap.overlap.lhs_begin);

//       biosoup::NucleicAcid rhs_{ "",
//                                 sequences[overlap.overlap.rhs_id]->InflateData(
//                                     overlap.overlap.rhs_begin,
//                                     overlap.overlap.rhs_end - overlap.overlap.rhs_begin) };

//       if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

//       auto rhs = rhs_.InflateData();

//       overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
//     }
//   };

//   std::vector<std::future<void>> alignment_futures;
//   alignment_futures.reserve(sequences.size());
//   for (std::uint32_t i = 0; i < sequences.size(); ++i) {
//     alignment_futures.emplace_back(thread_pool_->Submit(resolve_alignments, i));
//   }
//   for (const auto& it : alignment_futures) {
//     it.wait();
//   }
// };

void Graph_Constructor::LoadOverlapsAndAlign(
    Program_Parameters &param,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
    std::vector<std::vector<extended_overlap>> &extended_overlaps,
    const std::function<void(std::uint32_t)> &pile_aligned) {

  // const std::string overlaps_path =
  //     param.input_paf_path;

  std::ifstream file("/storage4/ftomas/jabba_trans/overlaps_features_fix_chr18/selected_top30_per_window.paf");

  if (!file.is_open()) {
    throw std::runtime_error(
        "Error opening file: ");
  }


  // ========================================================================
  // READ NAME -> SEQUENCE ID
  // ========================================================================

  std::map<std::string, std::uint32_t>
      sequence_name_to_seq_id;

  for (std::uint32_t i = 0;
       i < sequences.size();
       ++i) {

    sequence_name_to_seq_id[
        sequences[i]->name
    ] = sequences[i]->id;
  }


  std::cerr
      << "[raven::LoadAndAlign] loading overlaps from: "
      << std::endl;


  // ========================================================================
  // LOAD PAF
  //
  // Input PAF orientation:
  //
  //   query  = supporting read
  //   target = backbone/read being corrected
  //
  // Internal orientation:
  //
  //   lhs = backbone
  //   rhs = support
  // ========================================================================

  std::string line;

  std::uint64_t loaded = 0;
  std::uint64_t skipped_load = 0;


  while (std::getline(file, line)) {

    if (line.empty()) {
      continue;
    }


    std::istringstream iss(line);

    std::string item;
    std::vector<std::string> items;


    while (std::getline(
        iss,
        item,
        '\t')) {

      items.emplace_back(
          std::move(item));
    }


    // Mandatory PAF fields.
    if (items.size() < 12) {

      ++skipped_load;
      continue;
    }


    // ----------------------------------------------------------------------
    // PAF target -> lhs / backbone
    // ----------------------------------------------------------------------

    const auto lhs_it =
        sequence_name_to_seq_id.find(
            items[5]);


    // ----------------------------------------------------------------------
    // PAF query -> rhs / supporting read
    // ----------------------------------------------------------------------

    const auto rhs_it =
        sequence_name_to_seq_id.find(
            items[0]);


    if (lhs_it ==
            sequence_name_to_seq_id.end() ||
        rhs_it ==
            sequence_name_to_seq_id.end()) {

      ++skipped_load;
      continue;
    }


    const std::uint32_t lhs_seq_id =
        lhs_it->second;

    const std::uint32_t rhs_seq_id =
        rhs_it->second;


    // ----------------------------------------------------------------------
    // Coordinates
    //
    // PAF target coordinates -> lhs/backbone
    // PAF query coordinates  -> rhs/support
    //
    // These remain forward-coordinate intervals even for '-' overlaps.
    // ----------------------------------------------------------------------

    const std::uint32_t lhs_begin =
        static_cast<std::uint32_t>(
            std::stoul(items[7]));

    const std::uint32_t lhs_end =
        static_cast<std::uint32_t>(
            std::stoul(items[8]));


    const std::uint32_t rhs_begin =
        static_cast<std::uint32_t>(
            std::stoul(items[2]));

    const std::uint32_t rhs_end =
        static_cast<std::uint32_t>(
            std::stoul(items[3]));


    const bool strand =
        items[4] == "+";


    // ----------------------------------------------------------------------
    // Validate coordinates
    // ----------------------------------------------------------------------

    if (lhs_end <= lhs_begin ||
        rhs_end <= rhs_begin ||
        lhs_end >
            sequences[lhs_seq_id]->inflated_len ||
        rhs_end >
            sequences[rhs_seq_id]->inflated_len) {

      ++skipped_load;
      continue;
    }


    biosoup::Overlap overlap{
        lhs_seq_id,
        lhs_begin,
        lhs_end,

        rhs_seq_id,
        rhs_begin,
        rhs_end,

        255,
        strand
    };


    edlib_align tmp = {};


    extended_overlap total_ovlp{
        overlap,
        tmp,
        0,
        0
    };


    // Alignment has not been calculated yet.
    total_ovlp.correction_overlap =
        false;


    // Keep original input coordinates around if these fields are used
    // elsewhere.
    total_ovlp.lhs_begin_original =
        lhs_begin;

    total_ovlp.lhs_end_original =
        lhs_end;

    total_ovlp.rhs_begin_original =
        rhs_begin;

    total_ovlp.rhs_end_original =
        rhs_end;


    extended_overlaps[
        lhs_seq_id
    ].emplace_back(
        std::move(total_ovlp));


    ++loaded;
  }


  file.close();


  std::cerr
      << "[raven::LoadAndAlign] loaded "
      << loaded
      << " overlaps, skipped "
      << skipped_load
      << std::endl;

  AlignCorrectionOverlaps(sequences, extended_overlaps, param.kmer_len,
                          param.max_edit_rate, pile_aligned);
}

void Graph_Constructor::AlignCorrectionOverlaps(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
    std::vector<std::vector<extended_overlap>> &extended_overlaps,
    std::uint32_t kmer_len,
    double max_edit_rate,
    const std::function<void(std::uint32_t)> &pile_aligned) {

  // ========================================================================
  // TRIM TERMINAL EDLIB GAPS
  //
  // Edlib orientation:
  //
  //   query  = lhs = backbone
  //   target = rhs = support
  //
  // Therefore in the INTERNAL Edlib CIGAR:
  //
  //   I consumes lhs
  //   D consumes rhs
  //
  //
  // RHS orientation:
  //
  // '+' :
  //
  //   rhs_begin ----------------> rhs_end
  //
  // '-' :
  //
  //   rhs_begin <---------------- rhs_end
  //
  // because the extracted RHS interval is reverse-complemented.
  // ========================================================================

  auto trim_terminal_gaps =
      [](extended_overlap &ovlp) -> bool {

    auto &alignment =
        ovlp.edlib_alignment;

    auto &overlap =
        ovlp.overlap;


    if (alignment.cigar.empty()) {
      return false;
    }


    // ======================================================================
    // PARSE CIGAR
    // ======================================================================

    std::vector<
        std::pair<std::uint32_t, char>
    > ops;


    std::uint32_t num = 0;


    for (char c :
         alignment.cigar) {

      if (std::isdigit(
              static_cast<unsigned char>(c))) {

        num =
            num * 10 +
            static_cast<std::uint32_t>(
                c - '0');

        continue;
      }


      if (num == 0) {
        return false;
      }


      // Edlib standard CIGAR should contain only M/I/D.
      if (c != 'M' &&
          c != 'I' &&
          c != 'D') {

        return false;
      }


      ops.emplace_back(
          num,
          c);

      num = 0;
    }


    if (ops.empty()) {
      return false;
    }


    std::size_t lo = 0;
    std::size_t hi = ops.size();


    std::uint64_t removed_gap_bases = 0;


    // ======================================================================
    // REMOVE LEADING I/D
    // ======================================================================

    while (lo < hi &&
           (ops[lo].second == 'I' ||
            ops[lo].second == 'D')) {

      const std::uint32_t len =
          ops[lo].first;

      const char op =
          ops[lo].second;


      if (op == 'I') {

        // --------------------------------------------------------------
        // Internal I consumes lhs/backbone.
        //
        // Leading lhs sequence corresponds to lhs_begin.
        // --------------------------------------------------------------

        if (len >
            overlap.lhs_end -
            overlap.lhs_begin) {

          return false;
        }


        overlap.lhs_begin +=
            len;


      } else {

        // --------------------------------------------------------------
        // Internal D consumes rhs/support.
        //
        // '+':
        //      leading RHS -> rhs_begin += len
        //
        // '-':
        //      RHS was reverse-complemented.
        //      Leading aligned sequence corresponds to the high
        //      forward-coordinate end.
        //
        //      leading RHS -> rhs_end -= len
        // --------------------------------------------------------------

        if (len >
            overlap.rhs_end -
            overlap.rhs_begin) {

          return false;
        }


        if (overlap.strand) {

          overlap.rhs_begin +=
              len;

        } else {

          overlap.rhs_end -=
              len;
        }
      }


      removed_gap_bases +=
          len;

      ++lo;
    }


    // ======================================================================
    // REMOVE TRAILING I/D
    // ======================================================================

    while (hi > lo &&
           (ops[hi - 1].second == 'I' ||
            ops[hi - 1].second == 'D')) {

      const std::uint32_t len =
          ops[hi - 1].first;

      const char op =
          ops[hi - 1].second;


      if (op == 'I') {

        // --------------------------------------------------------------
        // Internal I consumes lhs/backbone.
        //
        // Trailing lhs sequence corresponds to lhs_end.
        // --------------------------------------------------------------

        if (len >
            overlap.lhs_end -
            overlap.lhs_begin) {

          return false;
        }


        overlap.lhs_end -=
            len;


      } else {

        // --------------------------------------------------------------
        // Internal D consumes rhs/support.
        //
        // '+':
        //      trailing RHS -> rhs_end -= len
        //
        // '-':
        //      trailing reverse-complemented RHS corresponds to the low
        //      forward-coordinate end.
        //
        //      trailing RHS -> rhs_begin += len
        // --------------------------------------------------------------

        if (len >
            overlap.rhs_end -
            overlap.rhs_begin) {

          return false;
        }


        if (overlap.strand) {

          overlap.rhs_end -=
              len;

        } else {

          overlap.rhs_begin +=
              len;
        }
      }


      removed_gap_bases +=
          len;

      --hi;
    }


    // ======================================================================
    // VALIDATE REMAINING ALIGNMENT
    // ======================================================================

    if (lo >= hi) {
      return false;
    }


    if (overlap.lhs_end <=
            overlap.lhs_begin ||
        overlap.rhs_end <=
            overlap.rhs_begin) {

      return false;
    }


    // After removing terminal gaps we expect the alignment to start
    // and end with aligned sequence.
    if (ops[lo].second != 'M' ||
        ops[hi - 1].second != 'M') {

      return false;
    }


    // ======================================================================
    // RECONSTRUCT TRIMMED CIGAR
    // ======================================================================

    std::string trimmed_cigar;

    trimmed_cigar.reserve(
        alignment.cigar.size());


    std::uint64_t new_block_length =
        0;


    for (std::size_t j = lo;
         j < hi;
         ++j) {

      trimmed_cigar +=
          std::to_string(
              ops[j].first);

      trimmed_cigar +=
          ops[j].second;


      // M, I and D all contribute alignment columns to the Edlib
      // alignment/block length.
      new_block_length +=
          ops[j].first;
    }


    alignment.cigar =
        std::move(
            trimmed_cigar);


    // ======================================================================
    // UPDATE ALIGNMENT STATISTICS
    // ======================================================================

    alignment.block_length =
        static_cast<decltype(
            alignment.block_length)>(
                new_block_length);


    // Every terminal I/D base contributes exactly one edit.
    if (removed_gap_bases >
        static_cast<std::uint64_t>(
            alignment.edit_distance)) {

      return false;
    }


    alignment.edit_distance -=
        static_cast<decltype(
            alignment.edit_distance)>(
                removed_gap_bases);


    // alignment.matches does NOT change because only gaps were removed.


    // ======================================================================
    // SAVE EXACT COORDINATES CORRESPONDING TO THE TRIMMED CIGAR
    // ======================================================================

    ovlp.correction_lhs_begin =
        overlap.lhs_begin;

    ovlp.correction_lhs_end =
        overlap.lhs_end;

    ovlp.correction_rhs_begin =
        overlap.rhs_begin;

    ovlp.correction_rhs_end =
        overlap.rhs_end;


    return true;
  };


  // ========================================================================
  // ALIGN EACH PILE
  // ========================================================================

  biosoup::Timer timer{};
  timer.Start();

  std::atomic<std::uint64_t> anchored_overlaps{0};
  std::atomic<std::uint64_t> anchors_used{0};
  std::atomic<std::uint64_t> anchors_rejected{0};
  std::atomic<std::uint64_t> over_edit_cap{0};
  std::atomic<std::uint64_t> aligned{0};
  std::atomic<std::uint64_t> rejected{0};

  auto resolve_alignments =
      [&](std::uint32_t i) -> void {


    const std::string lhs_full =
        sequences[i]->InflateData();


    for (auto &overlap :
         extended_overlaps[i]) {


      const std::uint32_t rhs_id =
          overlap.overlap.rhs_id;


      // ====================================================================
      // EXTRACT LHS / BACKBONE
      // ====================================================================

      const std::string lhs =
          lhs_full.substr(
              overlap.overlap.lhs_begin,
              overlap.overlap.lhs_end -
                  overlap.overlap.lhs_begin);


      // ====================================================================
      // EXTRACT RHS / SUPPORT
      // ====================================================================

      biosoup::NucleicAcid rhs_{
          "",
          sequences[rhs_id]->InflateData(
              overlap.overlap.rhs_begin,
              overlap.overlap.rhs_end -
                  overlap.overlap.rhs_begin)
      };


      // PAF coordinates remain forward-coordinate intervals.
      // Only the sequence itself is reverse-complemented.
      if (!overlap.overlap.strand) {

        rhs_.ReverseAndComplement();
      }


      const std::string rhs =
          rhs_.InflateData();


      if (lhs.empty() ||
          rhs.empty()) {

        overlap.correction_overlap =
            false;

        continue;
      }


      // ====================================================================
      // GLOBAL EDLIB ALIGNMENT
      //
      // With chain anchors: NW between consecutive anchors, mapped into the
      // extracted (on '-' reverse-complemented) lhs / rhs coordinates.
      // ====================================================================

      // Give up on overlaps of unrelated sequence (~0.5 edits per base, real
      // overlaps ~0.08) instead of aligning them in full.
      const std::int64_t max_edit_distance =
          max_edit_rate > 0 ?
              static_cast<std::int64_t>(
                  max_edit_rate * std::max(lhs.size(), rhs.size())) :
              -1;

      if (overlap.anchors.empty()) {

        overlap.edlib_alignment =
            edlib_wrapper(
                lhs,
                rhs,
                max_edit_distance);

      } else {

        std::vector<std::pair<std::uint32_t, std::uint32_t>> local;
        local.reserve(overlap.anchors.size());

        for (const auto &a : overlap.anchors) {
          if (a.first < overlap.overlap.lhs_begin ||
              a.first + kmer_len > overlap.overlap.lhs_end ||
              a.second < overlap.overlap.rhs_begin ||
              a.second + kmer_len > overlap.overlap.rhs_end) {
            continue;
          }
          local.emplace_back(
              a.first - overlap.overlap.lhs_begin,
              overlap.overlap.strand ?
                  a.second - overlap.overlap.rhs_begin :
                  overlap.overlap.rhs_end - (a.second + kmer_len));
        }

        std::uint32_t used = 0;
        std::uint32_t rejected = 0;
        overlap.edlib_alignment =
            anchored_edlib_wrapper(
                lhs,
                rhs,
                local,
                kmer_len,
                max_edit_distance,
                &used,
                &rejected);

        anchored_overlaps += 1;
        anchors_used += used;
        anchors_rejected += rejected + (overlap.anchors.size() - local.size());

        std::vector<std::pair<std::uint32_t, std::uint32_t>>().swap(
            overlap.anchors);
      }

      if (max_edit_distance >= 0 &&
          overlap.edlib_alignment.cigar.empty()) {
        over_edit_cap += 1;
      }


      // ====================================================================
      // REMOVE TERMINAL I/D
      //
      // This modifies both:
      //
      //   overlap.edlib_alignment.cigar
      //
      // and
      //
      //   overlap.overlap.{lhs,rhs}_{begin,end}
      //
      // so they describe exactly the same trimmed alignment.
      // ====================================================================

      if (!trim_terminal_gaps(
              overlap)) {

        overlap.correction_overlap =
            false;

        continue;
      }


      overlap.correction_overlap =
          true;
    }

    std::uint64_t pile_aligned_count = 0;
    for (const auto &overlap : extended_overlaps[i]) {
      pile_aligned_count += overlap.correction_overlap;
    }
    aligned += pile_aligned_count;
    rejected += extended_overlaps[i].size() - pile_aligned_count;

    if (pile_aligned) {
      pile_aligned(i);  // e.g. write the pile out and free it
    }
  };


  // ========================================================================
  // RUN ALIGNMENTS IN PARALLEL
  // ========================================================================

  std::vector<std::future<void>>
      alignment_futures;


  alignment_futures.reserve(
      sequences.size());


  // Largest piles first, so a few expensive ones (repeats aligned without
  // anchors) do not run alone at the end.
  std::vector<std::uint64_t> pile_bases(sequences.size(), 0);
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    for (const auto &overlap : extended_overlaps[i]) {
      pile_bases[i] += overlap.overlap.lhs_end - overlap.overlap.lhs_begin;
    }
  }
  std::vector<std::uint32_t> pile_order(sequences.size());
  std::iota(pile_order.begin(), pile_order.end(), 0);
  std::sort(pile_order.begin(), pile_order.end(),
            [&](std::uint32_t a, std::uint32_t b) {
              return pile_bases[a] > pile_bases[b];
            });


  for (std::uint32_t i : pile_order) {

    alignment_futures.emplace_back(
        thread_pool_->Submit(
            resolve_alignments,
            i));
  }


  for (const auto &it :
       alignment_futures) {

    it.wait();
  }


  alignment_futures.clear();


  // ========================================================================
  // DIAGNOSTICS
  // ========================================================================

  std::cerr
      << "[raven::LoadAndAlign] successfully aligned/trimmed "
      << aligned
      << " overlaps, rejected "
      << rejected
      << " (" << over_edit_cap << " above --max-edit-rate " << max_edit_rate << ")"
      << "; " << anchored_overlaps << " aligned between "
      << anchors_used << " anchors ("
      << anchors_rejected << " anchors unusable) "
      << std::fixed << timer.Stop() << "s"
      << std::endl;
}

void Graph_Constructor::ConstructOverlapsForCorrection(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                       std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                       biosoup::Timer &timer,
                                       Program_Parameters &param) {

  graph_.annotations_.resize(sequences.size());
                                      
  // The correction path reads only k-mer types and flags from piles, so they
  // go without the per-base coverage.
  for (const auto &it : sequences) {
    graph_.piles_.emplace_back(new Pile(it->id, it->inflated_len, /*coverage=*/false));
  }

  // Containment votes are applied as mapping stores overlaps, so piles of
  // reads that become strongly contained are freed during mapping.
  if (!param.dump_all_overlaps) {
    containment_votes_.assign(sequences.size(), 0);
    containment_threshold_ = 3;
  }

  MapSequencesFast(sequences, extended_overlaps, timer, param);

  // After mapping, nothing on the correction path reads k-mer counts.
  internal_counter_.reset();

  // Piles grew one overlap at a time; give back their spare capacity.
  {
    std::vector<std::future<void>> futures;
    futures.reserve(extended_overlaps.size());
    for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit([&](std::uint32_t j) -> void {
        extended_overlaps[j].shrink_to_fit();
        for (auto &it : extended_overlaps[j]) {
          it.anchors.shrink_to_fit();
        }
      }, i));
    }
    for (const auto &it : futures) {
      it.wait();
    }
  }
  LogOverlapMemory(extended_overlaps, "after mapping");

  if (param.dump_all_overlaps) {
    // Debug: raw ram overlaps (all fragments, no dedup/merge/features) are now
    // in extended_overlaps[i] with read i as query. Dump and stop before any
    // contained-read resolution or classification touches them.
    std::cerr << "Dumping ALL raw overlaps (no dedup/merge/features) to "
                 "all_overlaps.paf" << std::endl;
    graph_.PrintAllRawOverlaps(extended_overlaps, sequences, "all_overlaps.paf");
    return;
  }

  // graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_initial_overlaps_filename);
  // std::cerr << "[raven::Graph::Construct] initial overlaps printed"
  //           << std::endl;
  // exit(0);

  ResolveContainedReads(sequences, extended_overlaps, timer, containment_threshold_);

  std::cerr << "Resolved Contained!" << std::endl;
  //graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_contained_filename);
}

void Graph_Constructor::ConstructOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                       std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                       biosoup::Timer &timer,
                                       Program_Parameters &param) {

  graph_.annotations_.resize(sequences.size());
                                      
  for (const auto &it : sequences) {
    graph_.piles_.emplace_back(new Pile(it->id, it->inflated_len));
  }

  bool load_cigar = false;
  if (!param.load_paf.empty()) {
    LoadOverlapsFromPaf(sequences, extended_overlaps, load_cigar, param);
  } else {
    MapSequencesFast(sequences, extended_overlaps, timer, param);
    //MapSequences(sequences, extended_overlaps, timer, param);
  }
  graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_initial_overlaps_filename);
  std::cerr << "[raven::Graph::Construct] initial overlaps printed"
            << std::endl;

  std::cerr << "[raven::Graph::Construct] initial overlaps resolved"
          << std::endl;

//   if (!load_cigar) {
//     std::vector<std::future<void>> void_futures;
//     for (int i = 0; i < (int)sequences.size(); i++) {
//       void_futures.emplace_back(thread_pool_->Submit(
//         [&](std::uint32_t i) -> void {
//           find_pairwise_alignment(i, extended_overlaps[i], sequences, graph_);
//         },
//         i));
//     };
//     for (const auto &it : void_futures) {
//       it.wait();
//     }
//     void_futures.clear();
//   }

//   if (param.herro_snps_path == "") {
//     LoadAnnotations(sequences, extended_overlaps, param);
//   } else {
//     LoadHerroSNPs(param.herro_snps_path, sequences);
//   }

//  // graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_snp_filename);
//   ResolveSnps(sequences, extended_overlaps, timer, param);

//   //ResolveOverlapType(sequences, extended_overlaps, timer, param);
//   graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_snp_filename);

  ResolveContainedReads(sequences, extended_overlaps, timer, 3);

  std::cerr << "Resolved Contained!" << std::endl;
  graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_contained_filename);
  //exit(0);

 // ResolveChimericSequences(sequences, extended_overlaps, timer);
  //graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_chimeric_filename);
}

void Graph_Constructor::ConstructOverlapsFromGT(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                        std::vector<std::vector<extended_overlap>> &extended_overlaps,
                        biosoup::Timer &timer,
                        Program_Parameters &param){
  
  std::cerr << "[raven::Graph::Construct] Constructing overlaps from GT"
            << std::endl;
  graph_.annotations_.resize(sequences.size());

  for (const auto &it : sequences) {
    graph_.piles_.emplace_back(new Pile(it->id, it->inflated_len));
  }
  std::ofstream piles_tmp("piles_2.txt");
  for (const auto &it : graph_.piles_) {
    piles_tmp << it->id() << "\t" 
              << sequences[it->id()]->name << "\t"
              << it->get_data().size() << "\t"
              << it->length() << std::endl;
  }
  LoadGTOverlaps(param.gt_overlaps, sequences, extended_overlaps, false);
  graph_.PrintOverlaps(extended_overlaps, sequences, true, "gt_overlaps.paf");
  TrimAndAnnotatePiles(sequences, extended_overlaps, timer, param);
  PrintPiles(sequences);
  graph_.PrintOverlaps(extended_overlaps, sequences, true, "beforeContainedGT.paf");
  ResolveContainedReadsGT(sequences, extended_overlaps, timer);
  graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_contained_filename);
  
};

void Graph_Constructor::LoadOverlapsFromPaf(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                         std::vector<std::vector<extended_overlap>> &extended_overlaps,
                         bool load_cigar,
                         Program_Parameters &param) {

  LoadOverlaps(param.load_paf, sequences, extended_overlaps, load_cigar);

  std::vector<std::future<void>> extended_layers_futures;
  std::uint16_t counter = 0;
  for (const auto &it : graph_.piles_) {
    counter += 1;
    //std::cerr << counter << std::endl;
    if (extended_overlaps[it->id()].empty()) {
      continue;
    }
    extended_layers_futures.emplace_back(thread_pool_->Submit(
      [&]() -> void {
        it->AddExtendedLayers(
          extended_overlaps[it->id()].begin(),
          extended_overlaps[it->id()].end());

      }));
  }

  for (const auto &it : extended_layers_futures) {
    it.wait();
  }

  extended_layers_futures.clear();
}

std::vector<std::pair<std::uint64_t, std::uint16_t>> window_min(const std::vector<std::pair<std::uint64_t, std::uint16_t>>& dst, std::uint16_t w_size = 64, bool keep_tail = true) {
  std::vector<std::pair<std::uint64_t, std::uint16_t>> mins;
  mins.reserve(dst.size() / w_size + 1);

  for (std::size_t i = 0; i < dst.size(); i += w_size) {
      std::size_t block_end = std::min<std::size_t>(i + w_size, dst.size());
      if (!keep_tail && block_end - i < w_size) break;   // drop short tail

      auto min_pair = *std::min_element(
          dst.begin() + i,
          dst.begin() + block_end,
          [](const auto& a, const auto& b) {
              return a.second < b.second;
          }
      );
      mins.push_back(min_pair);
  }
  return mins;   // one value per w_size window
}

void Graph_Constructor::MapSequencesFast(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                     std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                     biosoup::Timer &timer,
                                     Program_Parameters &param){
  std::size_t bytes = 0;

  ram::MinimizerEngine minimizer_engine{
    thread_pool_,
    param.kmer_len,
    param.window_len,
    param.bandwidth,
    param.chain_n,
    param.match_n,
    param.gap_size,
    param.fraction,
    param.coverage,
    0.5f,
    param.fastk_path,
    param.minimizers,
    100000,
    50,
    500U,
    4U,
    80U,
    10000U,
    50};
  const std::string cnnont_path =
    //"/home/ftomas/github/draven/third_party/pytorch_model/CNNONT_model_multi_v1.int8.ts.pt";
    "/storage4/ftomas/jabba_trans/CNNONT_model_multi_uncapped.ts.pt";
  // /home/ftomas/github/CNNONT/mlruns/4/f1546b13a1e645cfa8e5d9340c82894d/artifacts/checkpoints/cnnont_model_best.pth
  at::set_num_threads(1);
  at::set_num_interop_threads(1);
  std::vector<std::vector<float>> all_features {};

  // Set below: sketching also builds each read's minimizers (fused path).
  // They are kept for the index and mapping, or, for a two-pass index, only
  // counted into minimizer_counts and freed.
  bool store_minimizers = false;
  std::unique_ptr<ram::KmerCounter> minimizer_counts;

  auto sketch_sequence = [&](std::uint32_t i) -> void {
    std::vector<std::uint64_t> ids;
    std::vector<float> sketch;
    std::vector<std::uint32_t> avg_base_qualities;
    std::vector<std::uint32_t> min_base_qualities;

    if (store_minimizers) {
      // Fused: also builds and stores this read's minimizer sketch.
      minimizer_engine.SketchAndMinimizeFastK(
          sequences[i], 1U, ids, sketch, avg_base_qualities, min_base_qualities);
    } else {
      minimizer_engine.FastKSketchReadInto(
          sequences[i], 1U, ids, sketch, avg_base_qualities, min_base_qualities);
    }
    if (minimizer_counts) {
      minimizer_engine.CountMinimizers(sequences[i], minimizer_counts.get());
    }

    if (!sketch.empty()) {
      int over_30k = std::count_if(sketch.begin(), sketch.end(), [](float v) { return v > 30000.f; });
      float hor_pct = 100.f * over_30k / (float)sketch.size();
      graph_.piles_[i]->set_hor_percent(hor_pct);
      if (over_30k >= (int)(sketch.size() * 0.05)) {
        graph_.piles_[i]->set_is_hor();
      }
    }

    for(int l = 0; l < sketch.size(); l++){
      sketch[l] = std::min<float>(sketch[l], 200.f);
    };


    // Store per-position k-mer ids only if requested; otherwise they are
    // recomputed on demand during the Jaccard step (saves ~8 B/position).
    if (param.store_kmer_ids) {
      graph_.piles_[i]->set_k_kmer_ids(ids);
    }


    WindowedCnnInput cnn_input = ComputeWindowFeaturesInference(
        sketch,
        avg_base_qualities,
        min_base_qualities,
        10  // window size
    );

    const int64_t L = cnn_input.num_windows;
    if (!cnn_input.defined || L == 0) {
      graph_.piles_[i]->clear_k_mer_types();
      return;
    }

    const float g0 = 50.0f;
    const float g1 = 25.0f;

    torch::InferenceMode guard;

    // ---------------------------------------------------------
    // Thread-local model
    // ---------------------------------------------------------
    thread_local torch::jit::script::Module cnnont_local =
        torch::jit::load(cnnont_path, torch::kCPU);

    thread_local bool cnnont_inited = false;
    if (!cnnont_inited) {
      cnnont_local.eval();
      cnnont_inited = true;
    }

    // Thread-local buffers
    thread_local torch::Tensor g_buf;
    thread_local std::vector<torch::jit::IValue> inputs;

    // Input tensor: [1, 12, L]
    auto x = torch::from_blob(
        cnn_input.X.data(),
        {1, 12, L},
        torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU)
    );

    if (!g_buf.defined()) {
      g_buf = torch::empty(
          {1, 2},
          torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU));
    }
    g_buf[0][0] = g0;
    g_buf[0][1] = g1;

    inputs.clear();
    inputs.push_back(x);
    inputs.push_back(g_buf);

    torch::Tensor logits = cnnont_local.forward(inputs).toTensor();

    // ---------------------------------------------------------
    // Convert logits -> class IDs
    // Expected: [1, 3, L]
    // ---------------------------------------------------------
    torch::Tensor pred = logits.argmax(1).squeeze(0).to(torch::kCPU);  // [L]
    auto pred_acc = pred.accessor<int64_t, 1>();

    auto class_to_type = [](int64_t c) -> kMerType {
      switch (c) {
        case 0: return kMerType::Diploid;  // D
        case 1: return kMerType::Haploid;  // H
        case 2: return kMerType::Error;    // E
        default: return kMerType::None;
      }
    };

    // One class per window of 10 positions; positions past the last window
    // take its class. The pile keeps them per window (1 B per 10 positions).
    std::vector<kMerType> window_types(L);
    for (int64_t w = 0; w < L; ++w) {
      window_types[w] = class_to_type(pred_acc[w]);
    }
    graph_.piles_[i]->set_k_mer_type_windows(
        std::move(window_types), 10, ids.size());

    // Prevent potential optimization removal
    volatile float sink = logits.reshape({-1})[0].item<float>();
    (void)sink;
  };

  // Jaccard over two SORTED, DEDUPED vectors (same semantics as the previous
  // std::set version, but cache-friendly and allocation-light).
  auto calc_jaccard =[&](const std::vector<std::uint64_t>& a,
                         const std::vector<std::uint64_t>& b) -> float{
    std::size_t i = 0, j = 0;
    std::size_t intersection = 0;
    std::size_t union_count = 0;

    while(i < a.size() && j < b.size()){
      if(a[i] == b[j]){
        ++intersection;
        ++union_count;
        ++i;
        ++j;
      } else if(a[i] < b[j]){
        ++union_count;
        ++i;
      } else {
        ++union_count;
        ++j;
      };
    };

    union_count += (a.size() - i) + (b.size() - j);

    if (union_count == 0) return 0.0f;

    return static_cast<float>(intersection) / static_cast<float>(union_count);
  };

  // Gather the k-mer ids of the given regions of a pile into one sorted, deduped
  // vector (replaces per-region std::set construction + merging).
  auto collect_sorted_kmers =[&](std::uint32_t pile_id,
      const std::vector<std::pair<std::uint32_t, std::uint32_t>>& positions)
      -> std::vector<std::uint64_t>{
    std::vector<std::uint64_t> kmers;
    for(const auto& p : positions){
      if (param.store_kmer_ids) {
        graph_.piles_[pile_id]->append_kmers_in_region(p.first, p.second, kmers);
      } else {
        // recompute the region's k-mer ids from the read (same ids as stored)
        minimizer_engine.KmerIdsInRegion(*sequences[pile_id], p.first, p.second, kmers);
      }
    }
    std::sort(kmers.begin(), kmers.end());
    kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
    return kmers;
  };

  auto find_similarity_for_region =
  [&](std::uint32_t lhs_id,
      std::uint32_t rhs_id,
      kMerType k_mer_type,
      std::uint32_t lhs_start,
      std::uint32_t rhs_start,
      std::uint32_t lhs_end,
      std::uint32_t rhs_end,
      std::uint32_t k_mer_size) {

      const std::uint32_t lhs_from = std::min(lhs_start, lhs_end);
      const std::uint32_t lhs_to   = std::max(lhs_start, lhs_end);

      const std::uint32_t rhs_from = std::min(rhs_start, rhs_end);
      const std::uint32_t rhs_to   = std::max(rhs_start, rhs_end);

      std::vector<std::pair<std::uint32_t, std::uint32_t>> lhs_positions =
          graph_.piles_[lhs_id]->find_region_positions(
              k_mer_type,
              lhs_from > k_mer_size ? lhs_from - k_mer_size : 0,
              lhs_to   > k_mer_size ? lhs_to   - k_mer_size : 0);

      std::vector<std::pair<std::uint32_t, std::uint32_t>> rhs_positions =
          graph_.piles_[rhs_id]->find_region_positions(
              k_mer_type,
              rhs_from > k_mer_size ? rhs_from - k_mer_size : 0,
              rhs_to   > k_mer_size ? rhs_to   - k_mer_size : 0);

      if (!lhs_positions.empty() && !rhs_positions.empty()) {
          std::vector<std::uint64_t> lhs_kmers = collect_sorted_kmers(lhs_id, lhs_positions);
          std::vector<std::uint64_t> rhs_kmers = collect_sorted_kmers(rhs_id, rhs_positions);
          return calc_jaccard(lhs_kmers, rhs_kmers);
      } else {
          return 0.0f;
      }
  };

    auto safe_region = [&](uint32_t a, uint32_t b, uint32_t k) {
      uint32_t from = std::min(a, b);
      uint32_t to   = std::max(a, b);

      if (to <= k) return std::pair<uint32_t,uint32_t>(0,0);

      from = from > k ? from - k : 0;
      to   = to   > k ? to   - k : 0;

      return std::make_pair(from, to);
  };

  // Chain anchors are only used to align the RF-selected correction overlaps.
  const bool keep_anchors = param.correction_overlaps &&
                            !param.rf_model_path.empty() &&
                            param.anchor_spacing > 0 &&
                            !param.dump_all_overlaps &&
                            !param.dump_overlaps_only;

  // Correction path: containment votes are applied as overlaps are stored
  // (StoreCorrectionCandidate), so a read at the threshold is final.
  const bool tracking_containment = !containment_votes_.empty();
  auto is_contained = [&](std::uint32_t id) -> bool {
    return __atomic_load_n(&containment_votes_[id], __ATOMIC_RELAXED) >=
           containment_threshold_;
  };

  auto map_sequences = [&](std::uint32_t i) -> std::vector<extended_overlap> { // map sequences
  //auto map_sequences = [&](std::uint32_t i) -> std::vector<std::vector<float>> { // map sequences
    std::vector<biosoup::Overlap> ovlps;
    std::vector<ram::MinimizerEngine::Anchors> ovlp_anchors;
    auto *anchors_out = keep_anchors ? &ovlp_anchors : nullptr;
    if(graph_.piles_[i]->is_hor()){
     ovlps = minimizer_engine.MapRepetitive(sequences[i], true, true,
                                                               false, nullptr, anchors_out);
    } else {
      ovlps = minimizer_engine.Map(sequences[i], true, true,
                                                                false, nullptr, anchors_out);
    }
    // ovlps is reordered below; keep the fragments ovlp_anchors is indexed by
    std::vector<biosoup::Overlap> fragments;
    if (keep_anchors) {
      fragments = ovlps;
    }

    if (param.dump_all_overlaps) {
      // Debug: keep EVERY raw fragment ram found (no dedup, no merge, no
      // feature extraction, no alignment) so raw ram coordinates can be
      // compared against minimap2. Query is read i (lhs); coordinates are the
      // raw anchor spans exactly as ram returned them.
      std::vector<extended_overlap> raw;
      raw.reserve(ovlps.size());
      for (const auto &ovlp : ovlps) {
        extended_overlap eo{ ovlp, {}, 0, 0 };
        eo.lhs_begin_original = ovlp.lhs_begin;
        eo.lhs_end_original   = ovlp.lhs_end;
        eo.rhs_begin_original = ovlp.rhs_begin;
        eo.rhs_end_original   = ovlp.rhs_end;
        eo.found_matches      = ovlp.score;
        eo.found_length       = std::max(ovlp.lhs_end - ovlp.lhs_begin,
                                         ovlp.rhs_end - ovlp.rhs_begin);
        raw.push_back(eo);
      }
      return raw;
    }

    if (!ovlps.empty()) {
      std::vector<extended_overlap> ovlps_final;
      std::vector<extended_overlap> ovlps_tmp;
      std::vector<biosoup::Overlap> no_dups_overlaps;
      std::vector<std::vector<float>> features {};
        
      if (param.merge_fragments) {
        // Merge collinear fragments per pair -> outermost-anchor span.
        no_dups_overlaps = merge_fragments(ovlps);
      } else {
        // Keep only the longest fragment per pair, then extrapolate it.
        std::sort(ovlps.begin(), ovlps.end(),
            [&](const biosoup::Overlap &lhs,
                const biosoup::Overlap &rhs) -> bool {
              return overlap_length(lhs) > overlap_length(rhs);
            });
        no_dups_overlaps = remove_duplicates(ovlps);
      }

      std::vector<biosoup::Overlap> tmp;

      for (auto &ovlp : no_dups_overlaps) {

        // Both reads already strongly contained: neither pile is corrected and
        // the overlap's votes and backbone marks go to contained reads only.
        if (tracking_containment && is_contained(i) && is_contained(ovlp.rhs_id)) {
          continue;
        }

        // before extrapolation: anchors are looked up by the found span
        std::vector<std::pair<std::uint32_t, std::uint32_t>> anchors;
        if (keep_anchors) {
          anchors = collect_anchors(ovlp, fragments, ovlp_anchors,
                                    param.anchor_spacing);
        }

        std::uint32_t left_overhang = 0;
        std::uint32_t right_overhang = 0;
        const std::uint32_t lhs_len = sequences[i]->inflated_len;
        const std::uint32_t rhs_len = sequences[ovlp.rhs_id]->inflated_len;

        std::uint32_t lhs_begin_original = ovlp.lhs_begin;
        std::uint32_t lhs_end_original   = ovlp.lhs_end;
        std::uint32_t rhs_begin_original = ovlp.rhs_begin;
        std::uint32_t rhs_end_original   = ovlp.rhs_end;

        if (ovlp.strand) {
            // Forward/forward overlap
            left_overhang = std::min(ovlp.lhs_begin, ovlp.rhs_begin);
            right_overhang = std::min(lhs_len - ovlp.lhs_end,
                                      rhs_len - ovlp.rhs_end);

            ovlp.lhs_begin -= left_overhang;
            ovlp.rhs_begin -= left_overhang;
            ovlp.lhs_end   += right_overhang;
            ovlp.rhs_end   += right_overhang;

        } else {
            // Forward/reverse overlap
            left_overhang = std::min(ovlp.lhs_begin,
                                    rhs_len - ovlp.rhs_end);

            right_overhang = std::min(lhs_len - ovlp.lhs_end,
                                      ovlp.rhs_begin);

            ovlp.lhs_begin -= left_overhang;
            ovlp.rhs_end   += left_overhang;
            ovlp.lhs_end   += right_overhang;
            ovlp.rhs_begin -= right_overhang;
        }
          

        extended_overlap total_ovlp{ ovlp, {}, 0, 0 };
        total_ovlp.anchors = std::move(anchors);
        total_ovlp.lhs_begin_original = lhs_begin_original;
        total_ovlp.lhs_end_original = lhs_end_original;
        total_ovlp.rhs_begin_original = rhs_begin_original;
        total_ovlp.rhs_end_original = rhs_end_original;
        total_ovlp.found_length = std::max(lhs_end_original - lhs_begin_original,
                                          rhs_end_original - rhs_begin_original);
        total_ovlp.extended_length = overlap_length(ovlp);
        total_ovlp.found_matches = ovlp.score;

        std::vector<std::pair<std::uint32_t, std::uint32_t>> lhs_hap_positions = graph_.piles_[ovlp.lhs_id]->find_region_positions(kMerType::Haploid, ovlp.lhs_begin > param.kmer_len ? ovlp.lhs_begin - param.kmer_len : 0, ovlp.lhs_end - param.kmer_len); // hardcoded to 21 ## TODO
        std::vector<std::pair<std::uint32_t, std::uint32_t>> rhs_hap_positions = graph_.piles_[ovlp.rhs_id]->find_region_positions(kMerType::Haploid, ovlp.rhs_begin > param.kmer_len ? ovlp.rhs_begin - param.kmer_len : 0, ovlp.rhs_end - param.kmer_len);
        
        if(!lhs_hap_positions.empty() && !rhs_hap_positions.empty()){
          total_ovlp.hap_regions = true;
          std::vector<std::uint64_t> lhs_hap_kmers = collect_sorted_kmers(ovlp.lhs_id, lhs_hap_positions);
          std::vector<std::uint64_t> rhs_hap_kmers = collect_sorted_kmers(ovlp.rhs_id, rhs_hap_positions);
          total_ovlp.jaccard_index = calc_jaccard(lhs_hap_kmers, rhs_hap_kmers);
        };
        total_ovlp.lhs_hap = graph_.piles_[total_ovlp.overlap.lhs_id]->return_haploid(
          ovlp.lhs_begin > param.kmer_len ? ovlp.lhs_begin - param.kmer_len : 0,
          ovlp.lhs_end - param.kmer_len);

        total_ovlp.rhs_hap = graph_.piles_[total_ovlp.overlap.rhs_id]->return_haploid(
          ovlp.rhs_begin > param.kmer_len ? ovlp.rhs_begin - param.kmer_len : 0,
          ovlp.rhs_end - param.kmer_len); 

        total_ovlp.score_to_length = static_cast<float>(total_ovlp.overlap.score) /
                                      static_cast<float>(overlap_length(total_ovlp.overlap));
        total_ovlp.found_to_extended_length = static_cast<float>(total_ovlp.found_length) /
                                      static_cast<float>(overlap_length(total_ovlp.overlap));
        total_ovlp.hap_ratio = total_ovlp.lhs_hap / (total_ovlp.rhs_hap + 0.0001f);

        total_ovlp.q_hor = graph_.piles_[total_ovlp.overlap.lhs_id]->is_hor();
        total_ovlp.t_hor = graph_.piles_[total_ovlp.overlap.rhs_id]->is_hor();

        

 
        bool has_hap_signal = (total_ovlp.lhs_hap > 0.0005) && (total_ovlp.rhs_hap > 0.0005);
        bool both_hor = total_ovlp.q_hor && total_ovlp.t_hor;
        bool low_hap = !has_hap_signal;

        if (both_hor || low_hap) {
            std::vector<std::pair<std::uint32_t, std::uint32_t>> lhs_dip_positions = graph_.piles_[ovlp.lhs_id]->find_region_positions(kMerType::Diploid,
                ovlp.lhs_begin > param.kmer_len ? ovlp.lhs_begin - param.kmer_len : 0,
                ovlp.lhs_end - param.kmer_len);
            std::vector<std::pair<std::uint32_t, std::uint32_t>> rhs_dip_positions = graph_.piles_[ovlp.rhs_id]->find_region_positions(kMerType::Diploid,
                ovlp.rhs_begin > param.kmer_len ? ovlp.rhs_begin - param.kmer_len : 0,
                ovlp.rhs_end - param.kmer_len);

            float dip_jaccard = 0.0f;
            if (!lhs_dip_positions.empty() && !rhs_dip_positions.empty()) {
                std::vector<std::uint64_t> lhs_dip_kmers = collect_sorted_kmers(ovlp.lhs_id, lhs_dip_positions);
                std::vector<std::uint64_t> rhs_dip_kmers = collect_sorted_kmers(ovlp.rhs_id, rhs_dip_positions);
                dip_jaccard = calc_jaccard(lhs_dip_kmers, rhs_dip_kmers);
            }
            if(both_hor){
              total_ovlp.diploid_jaccard_full_overlap = dip_jaccard;
              total_ovlp.classification_label = (dip_jaccard > 0.33) ? 2 : 0;
            }else {
              total_ovlp.diploid_jaccard_full_overlap = dip_jaccard;
              total_ovlp.classification_label = (dip_jaccard > 0.66) ? 2 : 0;
          }
        } else {
            total_ovlp.classification_label = ((total_ovlp.jaccard_index > 0.05) && has_hap_signal) ? 1 : 0;
            if(total_ovlp.classification_label != 1 && total_ovlp.jaccard_index > 0.025){
              total_ovlp.classification_label = 2;
            }
        }
        ovlps_final.emplace_back(total_ovlp);

      }
      return ovlps_final;
    }

    std::vector<extended_overlap> total_ovlps{};
    return total_ovlps;
  };

  
  // In the FastK path we fuse per-read feature extraction with minimizer
  // construction (one pass, one FastK lookup per base) and reuse the stored
  // minimizers for both index building and mapping. The minimizers-mode path
  // is left byte-identical.
  const bool use_fused_fastk = !param.minimizers && !param.no_fuse;

  // Count source: in-process internal counter (--internal-count) or the
  // precomputed FastK table. The counter must outlive every engine call below.
  std::unique_ptr<ram::KmerCounter> internal_counter;  // local fallback
  if (param.internal_count) {
    if (internal_counter_ != nullptr) {
      // Preloaded by LoadAndCount (counting was overlapped with parsing).
      minimizer_engine.set_internal_counter(internal_counter_.get());
    } else {
      // Not preloaded: count here (not overlapped with the read parse).
      // Bloom sizing (bloom_capacity_hint is a minimum bit count):
      //   --bloom-gb < 0  -> auto: size by total k-mer instances
      //   --bloom-gb == 0 -> disable bloom: exact count + prune
      //   --bloom-gb  > 0 -> that many GB of bloom bits (1 GB = 2^33 bits)
      std::uint64_t bloom_hint;
      if (param.bloom_gb < 0.0) {
        bloom_hint = 0;
        for (const auto& s : sequences) {
          bloom_hint += s->inflated_len;
        }
        std::cerr << "[raven::] internal counter bloom: auto (by k-mer instances)"
                  << std::endl;
      } else if (param.bloom_gb == 0.0) {
        bloom_hint = 0;  // disable bloom -> exact
        std::cerr << "[raven::] internal counter bloom: disabled (exact count)"
                  << std::endl;
      } else {
        bloom_hint = static_cast<std::uint64_t>(
            param.bloom_gb * static_cast<double>(1ULL << 33));  // GB -> bits
        std::cerr << "[raven::] internal counter bloom: " << param.bloom_gb
                  << " GB" << std::endl;
      }
      internal_counter = std::make_unique<ram::KmerCounter>(
          param.kmer_len, thread_pool_, /*min_count=*/2,
          /*bloom_capacity_hint=*/bloom_hint);
      timer.Start();
      internal_counter->Count(sequences.begin(), sequences.end());
      std::cerr << "[raven::] internal k-mer count: "
                << internal_counter->num_distinct() << " distinct (>= 2) "
                << std::fixed << timer.Stop() << "s" << std::endl;
      minimizer_engine.set_internal_counter(internal_counter.get());
    }
  } else if(!param.minimizers){
    minimizer_engine.LoadFastK();
  }

  if (keep_anchors) {
    // A unique minimizer occurs about once per covering read; repeat copies
    // multiply that, and their matches may pin the alignment to another copy.
    const ram::KmerCounter *counter =
        internal_counter_ ? internal_counter_.get() : internal_counter.get();
    biosoup::Timer peak_timer{};
    peak_timer.Start();
    const std::uint32_t peak = counter ? counter->EstimatePeak() : 0;
    const std::uint32_t max_occurrence = 2 * (peak ? peak : param.coverage);
    minimizer_engine.set_anchor_max_occurrence(max_occurrence);
    std::cerr << "[raven::] alignment anchors: minimizers with > "
              << max_occurrence << " occurrences left out (coverage peak "
              << peak << ") " << std::fixed << peak_timer.Stop() << "s" << std::endl;
  }

  // One minimizer index over all reads, and every read mapped against it once
  // (Map keeps pairs with query id <= target id, so each pair is found once).
  // Up to param.two_pass_index_bases the fused path keeps every read's
  // minimizers from sketching for the index and mapping. Above it keeping
  // them would cost ~16 B per minimizer of the whole read set, so sketching
  // only counts them, the index is filled in a second minimizing pass and Map
  // recomputes each read's.
  std::uint64_t total_bases = 0;
  for (const auto &it : sequences) {
    total_bases += it->inflated_len;
  }
  const bool two_pass = total_bases > param.two_pass_index_bases;
  store_minimizers = use_fused_fastk;
  if (store_minimizers) {
    minimizer_engine.InitReadStore(sequences.size());
  }
  if (two_pass) {
    // bloom sized at ~4x the minimizer instances (~14 bits per distinct one)
    const std::uint64_t instances = 2 * total_bases / (param.window_len + 1);
    minimizer_counts = std::make_unique<ram::KmerCounter>(
        32, thread_pool_, /*min_count=*/2, /*bloom_capacity_hint=*/4 * instances);
  }

  // ---- sketch every read first: mapping computes features from the k-mer
  //      types of both reads of an overlap ----
  timer.Start();
  std::cerr << "Starting sketching" << std::endl;
  {
    std::vector<std::future<void>> sketch_futures;
    sketch_futures.reserve(sequences.size());
    for (std::uint32_t k = 0; k < sequences.size(); ++k) {
      sketch_futures.emplace_back(thread_pool_->Submit(sketch_sequence, k));
    }
    for (const auto &it : sketch_futures) {
      it.wait();
    }
  }
  std::cerr << "[raven::Graph::Construct] sketched sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  graph_.PrintPiles();

  // Correcting, qualities are read only by sketching (CNN features); free
  // them (1 B per base). Assembly keeps them.
  if (param.correction_overlaps) {
    for (auto &it : sequences) {
      std::vector<std::uint8_t>().swap(it->quality);
    }
  }

  // ---- build the minimizer index ----
  timer.Start();
  if (two_pass) {
    minimizer_counts->Finalize();
    std::cerr << "[raven::Graph::Construct] counted "
              << minimizer_counts->num_distinct()
              << " minimizers occurring more than once" << std::endl;
    minimizer_engine.BuildIndexFromCounts(
      *minimizer_counts, sequences.begin(), sequences.end());
    minimizer_counts.reset();
  } else if (store_minimizers) {
    minimizer_engine.BuildIndexFromStored(sequences.begin(), sequences.end());
  } else {
    minimizer_engine.Minimize(sequences.begin(), sequences.end(), true);
  }

  minimizer_engine.Filter(param.freq);

  std::cerr << "[raven::Graph::Construct] minimized "
            << sequences.size() << " sequences" << (two_pass ? " (two-pass index) " : " ")
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  // ---- map every read against it ----
  {
    std::vector<std::future<std::vector<extended_overlap>>> thread_futures;

    timer.Start();
    for (std::uint32_t k = 0; k < sequences.size(); ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(map_sequences, k));

        bytes += sequences[k]->inflated_len;
        if (k != sequences.size() - 1 && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        for (auto &it : thread_futures) {
          for (auto &jt : it.get()) {
            if (tracking_containment) {
              extended_overlap reversed = feature_overlap_reverse(jt);
              StoreCorrectionCandidate(std::move(jt), extended_overlaps, sequences);
              StoreCorrectionCandidate(std::move(reversed), extended_overlaps, sequences);
              continue;
            }
            extended_overlaps[jt.overlap.lhs_id].emplace_back(jt);
            // In raw-dump mode keep each fragment once, oriented with read i
            // (lhs) as query -- no reversed copy on the target pile.
            if (!param.dump_all_overlaps) {
              extended_overlaps[jt.overlap.rhs_id].emplace_back(feature_overlap_reverse(jt));
            }
          }
        }
        thread_futures.clear();
      std::cerr << "[raven::Graph::Construct] mapped sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
  }

  // ---- release the stored minimizers ----
  if (store_minimizers) {
    minimizer_engine.ClearReadStore();
  }
}

/*void Graph_Constructor::MapSequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                     std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                     biosoup::Timer &timer,
                                     Program_Parameters &param) {
  std::size_t bytes = 0;

  ram::MinimizerEngine minimizer_engine{
    thread_pool_,
    param.kmer_len,
    param.window_len,
    param.bandwidth,
    param.chain_n,
    param.match_n,
    param.gap_size
  };

  auto map_sequences = [&](std::uint32_t i) -> std::vector<extended_overlap> { // map sequences
    std::vector<biosoup::Overlap> ovlps = minimizer_engine.Map(sequences[i], true, true,
                                                               false);

    if (!ovlps.empty()) {
      std::vector<extended_overlap> ovlps_final;
      std::sort(ovlps.begin(), ovlps.end(),
                [&](const biosoup::Overlap &lhs,
                    const biosoup::Overlap &rhs) -> bool {
                  return overlap_length(lhs) > overlap_length(rhs);
                });
      
      if(ovlps.size() > minimizer_engine.hom_peak()*1.5){
        ovlps.resize(minimizer_engine.hom_peak()*1.5);
      }

      std::vector<biosoup::Overlap> tmp;

      for (auto &ovlp : ovlps) {
        if (overlap_length(ovlp) > sequences[i]->inflated_len*0.05) {

          ovlp.lhs_begin = ovlp.lhs_begin - (param.window_len + param.kmer_len - 1) ? ovlp.lhs_begin
            - (param.window_len + param.kmer_len - 1) : 0;
          ovlp.lhs_end =
            ovlp.lhs_end + (param.window_len + param.kmer_len - 1) < sequences[ovlp.lhs_id]->inflated_len ?
            ovlp.lhs_end + (param.window_len + param.kmer_len - 1) : sequences[ovlp.lhs_id]->inflated_len;

          ovlp.rhs_begin = ovlp.rhs_begin - (param.window_len + param.kmer_len - 1) ? ovlp.rhs_begin
            - (param.window_len + param.kmer_len - 1) : 0;
          ovlp.rhs_end =
            ovlp.rhs_end + (param.window_len + param.kmer_len - 1) < sequences[ovlp.rhs_id]->inflated_len ?
            ovlp.rhs_end + (param.window_len + param.kmer_len - 1) : sequences[ovlp.rhs_id]->inflated_len;

          auto lhs = sequences[i]->InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

          biosoup::NucleicAcid rhs_{ "",
                                     sequences[ovlp.rhs_id]->InflateData(ovlp.rhs_begin,
                                                                         ovlp.rhs_end - ovlp.rhs_begin) };

          if (!ovlp.strand) rhs_.ReverseAndComplement();

          auto rhs = rhs_.InflateData();

          edlib_align tmp = edlib_wrapper(lhs, rhs);
          if (static_cast<float>(tmp.matches) / tmp.block_length > 0.9) {
            // edlib_align tmp;
            biosoup::Overlap ovlp_tmp{ ovlp.lhs_id, ovlp.lhs_begin, ovlp.lhs_end,
                                       ovlp.rhs_id, ovlp.rhs_begin, ovlp.rhs_end,
                                       ovlp.score, ovlp.strand };

            extended_overlap total_ovlp{ ovlp_tmp, tmp, 0, 0 };
            ovlps_final.emplace_back(total_ovlp);
          }
        }

      }
      return ovlps_final;
    }
    std::vector<extended_overlap> total_ovlps;
    return total_ovlps;
  };

  minimizer_engine.Count(sequences.begin(),
                         sequences.begin() + static_cast<std::ptrdiff_t>(sequences.size()*param.fraction),
                         param.fraction,
                         false);  // count k-mers in preconstructed minimizer index
  for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
    bytes += sequences[i]->inflated_len;
    if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
      continue;
    }
    bytes = 0;

    timer.Start();

    minimizer_engine.Minimize(
      sequences.begin() + j,
      sequences.begin() + i + 1,
      true);

    minimizer_engine.Filter(param.freq);

    std::cerr << "aaaaaaaaaa" << std::endl;
    std::cerr << "[raven::Graph::Construct] minimized "
              << j << " - " << i + 1 << " / " << sequences.size() << " "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    std::vector<std::uint32_t> num_overlaps(extended_overlaps.size());
    for (std::uint32_t k = 0; k < extended_overlaps.size(); ++k) {
      num_overlaps[k] = extended_overlaps[k].size();
    }

    std::vector<std::future<std::vector<extended_overlap>>> thread_futures;

    for (std::uint32_t k = 0; k < i + 1; ++k) {
      thread_futures.emplace_back(thread_pool_->Submit(map_sequences, k));

      bytes += sequences[k]->inflated_len;
      if (k != i && bytes < (1U << 30)) {
        continue;
      }
      bytes = 0;

      for (auto &it : thread_futures) {
        for (const auto &jt : it.get()) {
          extended_overlaps[jt.overlap.lhs_id].emplace_back(jt);
          //overlaps.emplace_back(jt.overlap);
          extended_overlaps[jt.overlap.rhs_id].emplace_back(cigar_extended_overlap_reverse(jt));
          //overlaps.emplace_back(overlap_reverse(jt.overlap));
        }
      }
      thread_futures.clear();
    }

    std::vector<std::future<void>> void_futures;
    for (const auto &it : graph_.piles_) {
      if (extended_overlaps[it->id()].empty()
        || extended_overlaps[it->id()].size() == num_overlaps[it->id()]
        ) {
        continue;
      }

      void_futures.emplace_back(thread_pool_->Submit(
        [&](std::uint32_t i) -> void {

          graph_.piles_[i]->AddExtendedLayers(
            extended_overlaps[i].begin(),
            extended_overlaps[i].end());

        },
        it->id()));
    }
    for (const auto &it : void_futures) {
      it.wait();
    }

    std::cerr << "[raven::Graph::Construct] mapped sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    j = i + 1;
  }
}*/

void Graph_Constructor::TrimAndAnnotatePiles(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                             std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                             biosoup::Timer &timer,
                                             Program_Parameters &param) {

  auto thread_func = [&](std::uint32_t i) -> void {
    graph_.piles_[i]->FindValidRegion(param.valid_region_coverage_threshold, param.valid_region_length_threshold);
    if (graph_.piles_[i]->is_invalid()) { // the sequence needs to be at least 512 bases long
      std::vector<extended_overlap>().swap(extended_overlaps[i]);

    } else {
      graph_.piles_[i]->FindMedian();
      //graph_.piles_[i]->FindChimericRegions();
    }
  };

  timer.Start();
  std::vector<std::future<void>> thread_futures;
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    thread_futures.emplace_back(thread_pool_->Submit(thread_func, i));
  }

  for (const auto &it : thread_futures) {
    it.wait();
  }

  thread_futures.clear();

  std::cerr << "[raven::Graph::Construct] annotated piles "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  std::ofstream outdata3;
  outdata3.open("valid_regions.fasta");
  for (int i = 0; i < (int)graph_.piles_.size(); i++) {
    outdata3 << ">" << graph_.piles_[i]->id() << "," << graph_.piles_[i]->begin() << "," << graph_.piles_[i]->end()
             << std::endl;
    outdata3 << std::endl;
    outdata3 << sequences[graph_.piles_[i]->id()]->InflateData(graph_.piles_[i]->begin(),
                                                               graph_.piles_[i]->end() - graph_.piles_[i]->begin())
             << std::endl;
    outdata3 << std::endl;

  }

  // std::ofstream chimeric_out;
  // chimeric_out.open("chimeric_regions.txt");
  // for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
  //   if (graph_.piles_[i]->is_maybe_chimeric()) {
  //     chimeric_out << sequences[i]->name << std::endl;
  //   }
  // }
  // chimeric_out.close();
}

void Graph_Constructor::ResolveSnps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                   std::vector<std::vector<extended_overlap>> &extended_overlaps,
                   biosoup::Timer &timer,
                   Program_Parameters &param) {
  timer.Start();
  std::vector<std::future<void>> futures;

  auto snp_match_thread_func = [&](std::uint32_t i) -> void {
    std::uint32_t k = 0;
    for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
      if (!overlap_update(extended_overlaps[i][j].overlap, graph_)) {
        continue;
      }

      const auto &it = extended_overlaps[i][j];

      auto lhs_anno = annotation_extract(
        it.overlap.lhs_id,
        it.overlap.lhs_begin,
        it.overlap.lhs_end,
        sequences[it.overlap.lhs_id]->inflated_len,
        true, graph_);

      auto rhs_anno = annotation_extract(
        it.overlap.rhs_id,
        it.overlap.rhs_begin,
        it.overlap.rhs_end,
        sequences[it.overlap.rhs_id]->inflated_len,
        it.overlap.strand, graph_);

      if (!lhs_anno.empty() || !rhs_anno.empty()) {
        //std::vector<std::pair<char, int>> cigar = parse_cigar_string(it.alignment);
        std::string edlib_alignment = cigar_to_edlib_alignment(it.edlib_alignment.cigar);
        std::uint32_t lhs_pos = it.overlap.lhs_begin;
        std::uint32_t rhs_pos = it.overlap.strand ?
                                it.overlap.rhs_begin :
                                sequences[it.overlap.rhs_id]->inflated_len - it.overlap.rhs_end;

        std::uint32_t mismatches = 0;
        std::uint32_t snps = 0;

        for (int t = 0; t < static_cast<std::int32_t>(edlib_alignment.length()); t++) {
          if (lhs_anno.find(lhs_pos) != lhs_anno.end() ||
            rhs_anno.find(rhs_pos) != rhs_anno.end()) {
            ++snps;
            if (edlib_alignment[t] == 3) {
              ++mismatches;
            }
          }
          switch (edlib_alignment[t]) {
            case 0:
            case 3: {
              ++lhs_pos;
              ++rhs_pos;
              break;
            }
            case 1: {
              ++lhs_pos;
              break;
            }
            case 2: {
              ++rhs_pos;
              break;
            }
            default:break;
          }
        }
        //outdata4 << sequences[it.overlap.lhs_id]->name << " " << sequences[it.overlap.rhs_id]->name << " " << mismatches << " " << snps << std::endl;
        extended_overlaps[i][j].total_overlap_snps = snps;
        extended_overlaps[i][j].total_overlap_snp_mismatches = mismatches;
        std::float_t identity = static_cast<float>(extended_overlaps[i][j].edlib_alignment.matches) / static_cast<float>(extended_overlaps[i][j].edlib_alignment.block_length);
        std::uint16_t matches = snps - mismatches;
        std::float_t missmatch_rate = static_cast<float>(mismatches) / static_cast<float>(snps);
        std::float_t heterozygosity_rate = static_cast<float>(snps) / static_cast<float>(extended_overlaps[i][j].edlib_alignment.block_length);
        extended_overlaps[i][j].identity = identity;
        extended_overlaps[i][j].heterozygosity_rate = heterozygosity_rate;
     //   if (mismatches / static_cast<double>(snps) > disagreement_) {
       //   continue;
        //}
      }

      //extended_overlaps[i][k++] = extended_overlaps[i][j];
    }
   // extended_overlaps[i].resize(k);
  };

  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    futures.emplace_back(thread_pool_->Submit(snp_match_thread_func, i));
  }


  for (const auto &it : futures) {
    it.wait();
  }

  futures.clear();

}


void Graph_Constructor::ResolveOverlapType(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                           std::vector<std::vector<extended_overlap>> &extended_overlaps,
                           biosoup::Timer &timer,
                           Program_Parameters &param) {

  timer.Start();
  std::vector<std::future<void>> futures;

  

  auto overlap_analysis_func = [&](std::uint32_t i) -> void {
    for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
      if (!overlap_update(extended_overlaps[i][j].overlap, graph_)) {
        continue;
      }
      std::float_t identity = static_cast<float>(extended_overlaps[i][j].edlib_alignment.matches) / static_cast<float>(extended_overlaps[i][j].edlib_alignment.block_length);
      std::uint16_t snps = extended_overlaps[i][j].total_overlap_snps;
      std::uint16_t snp_mismatches = extended_overlaps[i][j].total_overlap_snp_mismatches;
      std::uint16_t matches = snps - snp_mismatches;
      std::float_t missmatch_rate = static_cast<float>(snp_mismatches) / static_cast<float>(snps);
      std::float_t heterozygosity_rate = static_cast<float>(snps) / static_cast<float>(extended_overlaps[i][j].edlib_alignment.block_length);
      extended_overlaps[i][j].identity = identity;
      extended_overlaps[i][j].heterozygosity_rate = heterozygosity_rate;
    }
  };

  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    futures.emplace_back(thread_pool_->Submit(overlap_analysis_func, i));
  }

  for (const auto &it : futures) {
    it.wait();
  }

  futures.clear();

}

/* void Graph_Constructor::ResolveContainedReads(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                           std::vector<std::vector<extended_overlap>> &extended_overlaps,
                           biosoup::Timer &timer) {

  timer.Start();


  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    std::uint32_t k = 0;
    for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
      //if (!overlap_update(extended_overlaps[i][j].overlap, graph_)) {
      if(extended_overlaps[i][j].classification_label != 1){
        continue;
      }
      std::uint32_t type = overlap_type_extended(extended_overlaps[i][j].overlap, graph_, 
                                                 sequences[extended_overlaps[i][j].overlap.lhs_id]->inflated_len, 
                                                 sequences[extended_overlaps[i][j].overlap.rhs_id]->inflated_len);
      extended_overlaps[i][j].graph_overlap_type = type;
   //   if (type == 1 && !graph_.piles_[i]->is_maybe_chimeric()) {
      if(type == 1){
        graph_.piles_[i]->set_is_contained();
   //   } else if (type == 2 && !graph_.piles_[extended_overlaps[i][j].overlap.rhs_id]->is_maybe_chimeric()) {
      } else if (type == 2){
        graph_.piles_[extended_overlaps[i][j].overlap.rhs_id]->set_is_contained();
      } 
      else {
        //if(safe_overlap(extended_overlaps[i][j])){
          extended_overlaps[i][k++] = extended_overlaps[i][j]; // might be fine to remove temporarily
        //};
       }
    }
    extended_overlaps[i].resize(k); // this with the above line might be fine to remove
  }


  std::cerr << "[raven::Graph::Construct] removed contained sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;
} */

// --------------------------------------------------
// PASS 1 of ResolveContainedReads for one overlap copy:
// Use only strong overlaps (classification_label == 1)
// to detect containment and mark backbone reads.
// Label-2 (HOR) overlaps are never marked as contained.
// Returns the read that got a containment vote, or -1.
// --------------------------------------------------
std::uint32_t Graph_Constructor::ApplyContainmentPass1(
    extended_overlap& ext_ovlp,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::uint32_t* containment_votes) {

  ext_ovlp.cat = overlapCategory::None;

  if (ext_ovlp.classification_label != 1 &&
      ext_ovlp.classification_label != 2) {
    return -1;
  }

  std::uint32_t lhs_id = ext_ovlp.overlap.lhs_id;
  std::uint32_t rhs_id = ext_ovlp.overlap.rhs_id;

  std::uint32_t type = overlap_type_extended(
      ext_ovlp.overlap,
      graph_,
      sequences[lhs_id]->inflated_len,
      sequences[rhs_id]->inflated_len);

  ext_ovlp.graph_overlap_type = type;

  if (ext_ovlp.classification_label == 2) {
    // HOR overlap — never mark as contained, always backbone
    ext_ovlp.cat = overlapCategory::Backbone;
    graph_.piles_[lhs_id]->set_backbone();
    graph_.piles_[rhs_id]->set_backbone();
    return -1;
  }
  // label == 1 — vote-based containment logic (votes are read by mapping
  // workers while the main thread adds them)
  if (type == 1) {
    __atomic_fetch_add(&containment_votes[lhs_id], 1U, __ATOMIC_RELAXED);
    ext_ovlp.cat = overlapCategory::StrongContained;
    return lhs_id;
  }
  if (type == 2) {
    __atomic_fetch_add(&containment_votes[rhs_id], 1U, __ATOMIC_RELAXED);
    ext_ovlp.cat = overlapCategory::StrongContained;
    return rhs_id;
  }
  ext_ovlp.cat = overlapCategory::Backbone;
  graph_.piles_[lhs_id]->set_backbone();
  graph_.piles_[rhs_id]->set_backbone();
  return -1;
}

void Graph_Constructor::StoreCorrectionCandidate(
    extended_overlap&& ext_ovlp,
    std::vector<std::vector<extended_overlap>>& extended_overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  const std::uint32_t voted =
      ApplyContainmentPass1(ext_ovlp, sequences, containment_votes_.data());
  // Votes only grow, so a read at the threshold stays strongly contained
  // and its pile is never corrected: free it once, store nothing more in it.
  if (voted != static_cast<std::uint32_t>(-1) &&
      containment_votes_[voted] == containment_threshold_) {
    std::vector<extended_overlap>().swap(extended_overlaps[voted]);
  }
  const std::uint32_t pile = ext_ovlp.overlap.lhs_id;
  if (containment_votes_[pile] >= containment_threshold_) {
    return;
  }
  extended_overlaps[pile].emplace_back(std::move(ext_ovlp));
}

void Graph_Constructor::ResolveContainedReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::vector<std::vector<extended_overlap>>& extended_overlaps,
    biosoup::Timer& timer,
    std::uint32_t containment_threshold) {

  timer.Start();

  // PASS 1, unless already applied to every overlap as mapping stored it
  // (then piles of contained reads are already empty).
  std::vector<std::uint32_t> containment_votes;
  if (!containment_votes_.empty()) {
    containment_votes.swap(containment_votes_);
  } else {
    containment_votes.assign(sequences.size(), 0);
    for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
      for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
        ApplyContainmentPass1(extended_overlaps[i][j], sequences,
                              containment_votes.data());
      }
    }
  }

  // Apply containment threshold
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    if (containment_votes[i] >= containment_threshold) {
      graph_.piles_[i]->set_strong_contained();
    } else if (containment_votes[i] > 0) {
      // not enough votes to be contained — treat as backbone
      graph_.piles_[i]->set_backbone();
    }
  }

  // --------------------------------------------------
  // PASS 2:
  // Classify all overlaps using node marks from PASS 1
  // --------------------------------------------------
  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {

      auto& ext_ovlp = extended_overlaps[i][j];

      std::uint32_t lhs_id = ext_ovlp.overlap.lhs_id;
      std::uint32_t rhs_id = ext_ovlp.overlap.rhs_id;

      bool lhs_strong_contained = graph_.piles_[lhs_id]->is_strong_contained();
      bool rhs_strong_contained = graph_.piles_[rhs_id]->is_strong_contained();

      bool lhs_backbone = graph_.piles_[lhs_id]->is_backbone();
      bool rhs_backbone = graph_.piles_[rhs_id]->is_backbone();

      // Any overlap touching a strongly contained read
      if (lhs_strong_contained || rhs_strong_contained) {
        ext_ovlp.cat = overlapCategory::StrongContained;
        continue;
      }

      // Strong or HOR overlap
      if (ext_ovlp.classification_label == 1 ||
          ext_ovlp.classification_label == 2) {
        ext_ovlp.cat = overlapCategory::Backbone;
        continue;
      }

      // From here on, overlap is weak
      if (lhs_backbone && rhs_backbone) {
        ext_ovlp.cat = overlapCategory::WeakBackbone;
      } else if (lhs_backbone != rhs_backbone) {
        ext_ovlp.cat = overlapCategory::WeakConnecting;
      } else if (!lhs_backbone && !rhs_backbone) {
        ext_ovlp.cat = overlapCategory::Weak;
      } else {
        ext_ovlp.cat = overlapCategory::None;
      }
    }
  }
}

void Graph_Constructor::ResolveContainedReadsGT(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                            std::vector<std::vector<extended_overlap>> &extended_overlaps,
                            biosoup::Timer &timer){
  timer.Start();

  auto safe_overlap = [&](const extended_overlap &overlap) -> bool {
    return (overlap.ground_truth);
  };

  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    std::uint32_t k = 0;
    for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
      if (!overlap_update(extended_overlaps[i][j].overlap, graph_)) {
        continue;
      }
      std::uint32_t type = overlap_type(extended_overlaps[i][j].overlap, graph_);
      extended_overlaps[i][j].graph_overlap_type = type;
      //if (type == 1 && safe_overlap(extended_overlaps[i][j]) && extended_overlaps[i][j].total_overlap_snps > 0) {
      //if (extended_overlaps[i][j].graph_overlap_type == 1 && safe_overlap(extended_overlaps[i][j]) && extended_overlaps[i][j].total_overlap_snps > 0) {
     // if (extended_overlaps[i][j].graph_overlap_type == 1 && safe_overlap(extended_overlaps[i][j])) {
      if (extended_overlaps[i][j].graph_overlap_type == 1){
        graph_.piles_[i]->set_is_contained();
      } else //if (type == 2 && safe_overlap(extended_overlaps[i][j]) && extended_overlaps[i][j].total_overlap_snps > 0) {
      //  if (extended_overlaps[i][j].graph_overlap_type == 2 && safe_overlap(extended_overlaps[i][j]) && extended_overlaps[i][j].total_overlap_snps > 0) {
    //  if (extended_overlaps[i][j].graph_overlap_type == 2 && safe_overlap(extended_overlaps[i][j])) {
      if (extended_overlaps[i][j].graph_overlap_type == 2) {
        graph_.piles_[extended_overlaps[i][j].overlap.rhs_id]->set_is_contained();
      } 
     else {
       if(safe_overlap(extended_overlaps[i][j])){
         extended_overlaps[i][k++] = extended_overlaps[i][j]; // might be fine to remove temporarily
       };
      }
    }
   extended_overlaps[i].resize(k); // this with the above line might be fine to remove
  }


  std::cerr << "[raven::Graph::Construct] removed contained sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;
}

void Graph_Constructor::ResolveChimericSequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                              std::vector<std::vector<extended_overlap>> &extended_overlaps,
                              biosoup::Timer &timer) {

  timer.Start();

  while (true) {
    auto components = connected_components(sequences, extended_overlaps, graph_); // HERE
    for (const auto &it : components) {
      std::vector<std::uint16_t> medians;
      for (const auto &jt : it) {
        medians.emplace_back(graph_.piles_[jt]->median());
      }
      std::nth_element(
        medians.begin(),
        medians.begin() + medians.size() / 2,
        medians.end());
      std::uint16_t median = medians[medians.size() / 2];

      std::vector<std::future<void>> thread_futures;
      for (const auto &jt : it) {
        thread_futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t i) -> void {
            graph_.piles_[i]->ClearChimericRegions(median);
            if (graph_.piles_[i]->is_invalid()) {
              std::vector<extended_overlap>().swap(extended_overlaps[i]);
            }
          },
          jt));
      }
      for (const auto &it : thread_futures) {
        it.wait();
      }
      thread_futures.clear();
    }

    bool is_changed = false;
    for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
      std::uint32_t k = 0;
      for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
        if (overlap_update(extended_overlaps[i][j].overlap, graph_)) {
          extended_overlaps[i][k++] = extended_overlaps[i][j];
        } else {
          is_changed = true;
        }
      }
      extended_overlaps[i].resize(k);
    }

    if (!is_changed) {
      for (const auto &it : extended_overlaps) {
        for (const auto &jt : it) {
          std::uint32_t type = overlap_type(jt.overlap, graph_);
          if (type == 1) {
            graph_.piles_[jt.overlap.lhs_id]->set_is_contained();
            graph_.piles_[jt.overlap.lhs_id]->set_is_invalid();
          } else if (type == 2) {
            graph_.piles_[jt.overlap.rhs_id]->set_is_contained();
            graph_.piles_[jt.overlap.rhs_id]->set_is_invalid();
          }
        }
      }
      //extended_overlaps.clear();
      break;
    }
  }

  std::cerr << "[raven::Graph::Construct] removed chimeric sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;
}

// void Graph_Constructor::AlignBackboneReads(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
//                                                std::vector<std::vector<extended_overlap>> &overlaps,
//                                                biosoup::Timer &timer,
//                                                Program_Parameters &param){

//   auto resolve_alignments = [&](std::uint32_t i) -> void {
//     if(graph_.piles_[i]->is_backbone() && !graph_.piles_[i]->is_strong_contained()){
//       std::string lhs_full = sequences[i]->InflateData();
//       for(auto& overlap : overlaps[i]){
//         if (param.infix_align) {
//           // Infix (HW) alignment: align read i's raw anchor region as an infix
//           // of a PADDED target window, and recover the exact target coordinates
//           // from edlib -- robust to the anchor/extension coordinate imprecision.
//           if (overlap.lhs_end_original <= overlap.lhs_begin_original) {
//             continue;
//           }
//           std::string query = lhs_full.substr(
//               overlap.lhs_begin_original,
//               overlap.lhs_end_original - overlap.lhs_begin_original);

//           const std::uint32_t margin =
//               std::max<std::uint32_t>(500u, static_cast<std::uint32_t>(query.size()) / 20u);
//           const std::uint32_t tlen = sequences[overlap.overlap.rhs_id]->inflated_len;
//           const std::uint32_t w_start =
//               overlap.rhs_begin_original > margin ? overlap.rhs_begin_original - margin : 0;
//           const std::uint32_t w_end =
//               std::min<std::uint32_t>(overlap.rhs_end_original + margin, tlen);
//           if (w_end <= w_start) {
//             continue;
//           }
//           const std::uint32_t wl = w_end - w_start;

//           biosoup::NucleicAcid tw{ "",
//               sequences[overlap.overlap.rhs_id]->InflateData(w_start, wl) };
//           if (!overlap.overlap.strand) tw.ReverseAndComplement();
//           const std::string target_window = tw.InflateData();

//           overlap.edlib_alignment = edlib_wrapper_infix(query, target_window);

//           // Recover exact target coords from the aligned infix [ts, te] (window
//           // coordinates). For reverse overlaps the window was reverse-complemented,
//           // so map positions back to the forward target.
//           const std::int32_t ts = overlap.edlib_alignment.target_start;
//           const std::int32_t te = overlap.edlib_alignment.target_end;
//           if (overlap.overlap.strand) {
//             overlap.rhs_begin_original = w_start + static_cast<std::uint32_t>(ts);
//             overlap.rhs_end_original   = w_start + static_cast<std::uint32_t>(te) + 1;
//           } else {
//             overlap.rhs_begin_original = w_start + (wl - 1 - static_cast<std::uint32_t>(te));
//             overlap.rhs_end_original   = w_start + (wl - 1 - static_cast<std::uint32_t>(ts)) + 1;
//           }
//           overlap.correction_overlap = true;
//           continue;
//         }
//         if(overlap.classification_label != 1){
//           // std::string lhs_full = sequences[i]->InflateData();
//           std::string lhs = lhs_full.substr(
//                             overlap.lhs_begin_original,
//                             overlap.lhs_end_original - overlap.lhs_begin_original);
//           // biosoup::NucleicAcid lhs_{ "",
//           //                           sequences[overlap.overlap.lhs_id]->InflateData(overlap.overlap.lhs_begin,
//           //                                                                     overlap.overlap.lhs_end - overlap.overlap.lhs_begin)};
//           biosoup::NucleicAcid rhs_{ "",
//                                     sequences[overlap.overlap.rhs_id]->InflateData(overlap.rhs_begin_original,
//                                                                               overlap.rhs_end_original - overlap.rhs_begin_original) };

//           if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

//           auto rhs = rhs_.InflateData();
//       //    auto lhs = lhs_.InflateData();
//           overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
//           overlap.correction_overlap = true;
//         } else if(overlap.classification_label == 1){
//           std::string lhs = lhs_full.substr(
//                             overlap.overlap.lhs_begin,
//                             overlap.overlap.lhs_end - overlap.overlap.lhs_begin);
//           // biosoup::NucleicAcid lhs_{ "",
//           //                           sequences[overlap.overlap.lhs_id]->InflateData(overlap.overlap.lhs_begin,
//           //                                                                     overlap.overlap.lhs_end - overlap.overlap.lhs_begin)};
//           biosoup::NucleicAcid rhs_{ "",
//                                     sequences[overlap.overlap.rhs_id]->InflateData(overlap.overlap.rhs_begin,
//                                                                               overlap.overlap.rhs_end - overlap.overlap.rhs_begin) };

//           if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

//           auto rhs = rhs_.InflateData();
//       //    auto lhs = lhs_.InflateData();
//           overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
//           overlap.correction_overlap = true;        
//         } //else {
//         //   std::string lhs = lhs_full.substr(
//         //                     overlap.overlap.lhs_begin,
//         //                     overlap.overlap.lhs_end - overlap.overlap.lhs_begin);
//         //   biosoup::NucleicAcid rhs_{ "",
//         //                             sequences[overlap.overlap.rhs_id]->InflateData(overlap.overlap.rhs_begin,
//         //                                                                       overlap.overlap.rhs_end - overlap.overlap.rhs_begin) };

//         //   if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

//         //   auto rhs = rhs_.InflateData();
//         //   // overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
//         //   overlap.correction_overlap = true;
//           // for(auto& rhs_overlap : overlaps[overlap.overlap.rhs_id]){
//           //   if(rhs_overlap.overlap.rhs_id == i){
//           //     rhs_overlap.edlib_alignment.cigar = cigar_alignment_reverse(overlap.edlib_alignment.cigar);
//           //     rhs_overlap.edlib_alignment.matches = overlap.edlib_alignment.matches;
//           //     rhs_overlap.edlib_alignment.block_length = overlap.edlib_alignment.block_length;
//           //     rhs_overlap.edlib_alignment.edit_distance = overlap.edlib_alignment.edit_distance;
//           //     rhs_overlap.correction_overlap = true;
//           //     break;
//           //   }
//           // }
//       //  }
//       }
//     }
//   };

//   std::vector<std::future<void>> alignment_futures;
//   for(int i = 0; i < sequences.size(); i++){
//     alignment_futures.emplace_back(thread_pool_->Submit(resolve_alignments, i));
//     //resolve_alignments(i);
//   }

//   for (const auto &it : alignment_futures) {
//     it.wait();
//   }
  
//   alignment_futures.clear();

// }

void Graph_Constructor::AlignBackboneReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
    std::vector<std::vector<extended_overlap>> &overlaps,
    biosoup::Timer &timer,
    Program_Parameters &param) {

  // Extra sequence given to NW around the geometrically estimated
  // overlap boundaries.
  constexpr std::uint32_t kNwPadding = 50;


  auto resolve_alignments =
      [&](std::uint32_t i) -> void {

    // Only backbone reads are corrected/aligned.
    if (!graph_.piles_[i]->is_backbone() ||
        graph_.piles_[i]->is_strong_contained()) {

      return;
    }


    const std::string lhs_full =
        sequences[i]->InflateData();

    const std::uint32_t lhs_len =
        sequences[i]->inflated_len;


    for (auto &overlap : overlaps[i]) {

      const std::uint32_t rhs_id =
          overlap.overlap.rhs_id;

      const std::uint32_t rhs_len =
          sequences[rhs_id]->inflated_len;


      // ================================================================
      // START FROM THE GEOMETRICALLY EXTENDED OVERLAP FOR ALL CLASSES
      //
      // classification_label no longer changes which coordinates are
      // used for correction/alignment.
      // ================================================================

      std::uint32_t lhs_begin =
          overlap.overlap.lhs_begin;

      std::uint32_t lhs_end =
          overlap.overlap.lhs_end;

      std::uint32_t rhs_begin =
          overlap.overlap.rhs_begin;

      std::uint32_t rhs_end =
          overlap.overlap.rhs_end;


      // ================================================================
      // VALIDATE ORIGINAL GEOMETRIC COORDINATES
      // ================================================================

      if (lhs_end <= lhs_begin ||
          rhs_end <= rhs_begin ||
          lhs_end > lhs_len ||
          rhs_end > rhs_len) {

        continue;
      }


      // ================================================================
      // ADD NW PADDING
      //
      // Give the global alignment some freedom around the geometrically
      // predicted overlap endpoints.
      //
      // Coordinates are clamped to the physical read boundaries.
      //
      // In a normal dovetail overlap, geometric extension will already
      // place one sequence at its physical end on each side. Therefore
      // padding mostly gives extra sequence to the read which still has
      // an overhang.
      // ================================================================

      lhs_begin =
          lhs_begin > kNwPadding
              ? lhs_begin - kNwPadding
              : 0;

      rhs_begin =
          rhs_begin > kNwPadding
              ? rhs_begin - kNwPadding
              : 0;


      lhs_end +=
          std::min<std::uint32_t>(
              kNwPadding,
              lhs_len - lhs_end);

      rhs_end +=
          std::min<std::uint32_t>(
              kNwPadding,
              rhs_len - rhs_end);


      // ================================================================
      // EXTRACT THE EXACT SEQUENCES THAT WILL BE ALIGNED
      // ================================================================

      const std::string lhs =
          lhs_full.substr(
              lhs_begin,
              lhs_end - lhs_begin);


      biosoup::NucleicAcid rhs_{
          "",
          sequences[rhs_id]->InflateData(
              rhs_begin,
              rhs_end - rhs_begin)
      };


      // For reverse overlaps, align against the reverse complement of
      // the selected forward-coordinate RHS interval.
      if (!overlap.overlap.strand) {
        rhs_.ReverseAndComplement();
      }


      const std::string rhs =
          rhs_.InflateData();


      if (lhs.empty() ||
          rhs.empty()) {

        continue;
      }


      // ================================================================
      // CORRECTION COORDINATES
      //
      // IMPORTANT:
      //
      // These coordinates describe EXACTLY the substrings passed to
      // Edlib below.
      //
      // PrintCorrectedOverlaps should therefore use correction_* when
      // emitting the PAF coordinates alongside the resulting cg:Z CIGAR.
      //
      // Do NOT use overlap.overlap.* for the CIGAR-containing PAF,
      // because those are the unpadded geometrical coordinates.
      // ================================================================

      overlap.correction_lhs_begin =
          lhs_begin;

      overlap.correction_lhs_end =
          lhs_end;

      overlap.correction_rhs_begin =
          rhs_begin;

      overlap.correction_rhs_end =
          rhs_end;


      // ================================================================
      // FINAL GLOBAL ALIGNMENT
      // ================================================================

      overlap.edlib_alignment =
          edlib_wrapper(
              lhs,
              rhs);


      overlap.correction_overlap =
          true;
    }
  };


  // ====================================================================
  // RUN BACKBONE ALIGNMENTS IN PARALLEL
  // ====================================================================

  std::vector<std::future<void>>
      alignment_futures;


  for (std::uint32_t i = 0;
       i < sequences.size();
       ++i) {

    alignment_futures.emplace_back(
        thread_pool_->Submit(
            resolve_alignments,
            i));
  }


  for (const auto &it :
       alignment_futures) {

    it.wait();
  }


  alignment_futures.clear();
}

void Graph_Constructor::ConstructAssemblyGraphInPhases(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                               std::vector<std::vector<extended_overlap>> &overlaps,
                                               biosoup::Timer &timer,
                                               Program_Parameters &param){
  std::ofstream outdata_invalid;
  outdata_invalid.open("invalid_reads.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_invalid()) {
      outdata_invalid << sequences[i]->name << std::endl;
    }
  }

  std::ofstream outdata_contained;
  outdata_contained.open("contained_reads.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_contained()) {
      outdata_contained << sequences[i]->name << std::endl;
    }
  }
  std::uint16_t n_phases = 5;


  //exit(0);
  Node::num_objects = 0;
  Edge::num_objects = 0;

  std::vector<std::int32_t> sequence_to_node(graph_.piles_.size(), -1);

  
  for (const auto &it : graph_.piles_) {  // create nodes
    if (it->is_invalid() || it->is_contained() || overlaps[it->id()].empty()) {
      continue;
    }

    bool any_edge = false;

    for (int j = 0; j < (int)overlaps[it->id()].size(); j++) {
      if (overlap_type(overlaps[it->id()][j].overlap, graph_) > 2) {
        if (!graph_.piles_[overlaps[it->id()][j].overlap.rhs_id]->is_invalid()) {
          any_edge = true;
          break;
        }
      }
    }

    if (!any_edge) {
      continue;
    }

    std::unordered_set<std::uint32_t> annotations;
    for (const auto &jt : graph_.annotations_[it->id()]) {
      if (it->begin() <= jt && jt < it->end()) {
        annotations.emplace(jt - it->begin());
      }
    }
    graph_.annotations_[it->id()].swap(annotations);

    auto sequence = biosoup::NucleicAcid{
      sequences[it->id()]->name,
      sequences[it->id()]->InflateData(it->begin(), it->end() - it->begin()) };  // NOLINT
    sequence.id = it->id();

    sequence_to_node[it->id()] = Node::num_objects;

    auto node = std::make_shared<Node>(sequence);
    sequence.ReverseAndComplement();
    graph_.nodes_.emplace_back(node);
    graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
    node->pair = graph_.nodes_.back().get();
    node->pair->pair = node.get();

    if (it->id() < param.split) {
      node->color = 1;
      node->pair->color = 1;
    }
  }

  std::cerr << "[raven::Graph::Construct] stored " << graph_.nodes_.size() << " nodes "  // NOLINT
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  timer.Start();

  int counter = 0;

  graph_.PrintOverlaps(overlaps, sequences, true, param.paf_before_parsing_edges_filename);

  auto return_viable_overlaps = [&](const std::vector<extended_overlap> &overlaps, raven::Graph graph_) -> std::vector<extended_overlap>{
    
  };

  for(std::uint8_t round; round < n_phases; round++){
    std::cerr << "Phase: " << round << std::endl;
    for (const auto &it: graph_.piles_){
      if (it->is_invalid() || it->is_contained() || overlaps[it->id()].empty()) {
        continue;
      }
    bool any_edge = false;

    for(int j = 0; j < (int)overlaps[it->id()].size(); j++){
      if(overlap_type(overlaps[it->id()][j].overlap, graph_) > 2, round){
        if(!graph_.piles_[overlaps[it->id()][j].overlap.rhs_id]->is_invalid()){
          continue;
        }
          any_edge = true;
          break;
      }
    }
    if (!any_edge) {
      continue;
    }
    }
  };



  for (int i = 0; i < (int)overlaps.size(); i++) {

    for (auto &it : overlaps[i]) {  // create edges

      if (!overlap_finalize(it.overlap, graph_)) {
        continue;
      }

      counter++;
      auto tail_seq_id = sequence_to_node[it.overlap.lhs_id];
      auto head_seq_id = sequence_to_node[it.overlap.rhs_id];

      if (tail_seq_id == -1 || head_seq_id == -1) {
        continue;
      }
      auto tail = graph_.nodes_[sequence_to_node[it.overlap.lhs_id]].get();
      auto head = graph_.nodes_[sequence_to_node[it.overlap.rhs_id] + 1 - it.overlap.strand].get();

      auto length = it.overlap.lhs_begin - it.overlap.rhs_begin;
      auto length_pair =
        (graph_.piles_[it.overlap.rhs_id]->length() - it.overlap.rhs_end) -
          (graph_.piles_[it.overlap.lhs_id]->length() - it.overlap.lhs_end);

      if (it.overlap.score == 4) {
        std::swap(head, tail);
        length *= -1;
        length_pair *= -1;
      }

      auto edge = std::make_shared<Edge>(tail, head, length);
      graph_.edges_.emplace_back(edge);
      graph_.edges_.emplace_back(std::make_shared<Edge>(head->pair, tail->pair, length_pair));  // NOLINT
      edge->pair = graph_.edges_.back().get();
      edge->pair->pair = edge.get();

    }
  }
  

  std::cerr << "[raven::Graph::Construct] stored " << graph_.edges_.size() << " edges "  // NOLINT
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  graph_.PrintGfa(param.gfa_after_construction_filename, false);                                                
};

void Graph_Constructor::ConstructOverlapGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                               std::vector<std::vector<extended_overlap>> &overlaps,
                                               biosoup::Timer &timer,
                                               Program_Parameters &param) {

  std::vector<std::int32_t> sequence_to_node(graph_.piles_.size(), -1);
  for (const auto &it : graph_.piles_){
    if(it->is_invalid()){
      continue;
    }
    std::unordered_set<std::uint32_t> annotations;
    for (const auto &jt : graph_.annotations_[it->id()]) {
      if (it->begin() <= jt && jt < it->end()) {
        annotations.emplace(jt - it->begin());
      }
    }
    graph_.annotations_[it->id()].swap(annotations);

    auto sequence = biosoup::NucleicAcid{
      sequences[it->id()]->name,
      sequences[it->id()]->InflateData(it->begin(), it->end() - it->begin()) };  // NOLINT
    sequence.id = it->id();

    sequence_to_node[it->id()] = Node::num_objects;

    auto node = std::make_shared<Node>(sequence);
    sequence.ReverseAndComplement();
    graph_.nodes_.emplace_back(node);
    graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
    node->pair = graph_.nodes_.back().get();
    node->pair->pair = node.get();

    if (it->id() < param.split) {
      node->color = 1;
      node->pair->color = 1;
    }
  }
  std::cerr << "[raven::Graph::ConstructOverlapGraph] stored " << graph_.nodes_.size() << " nodes "  // NOLINT
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  timer.Start();
  int counter = 0;

  for (int i = 0; i < (int)overlaps.size(); i++) {

    for (auto &it : overlaps[i]) {  // create edges

      // if (!overlap_finalize(it.overlap, graph_)) {
      //   continue;
      // }

      it.overlap.score = overlap_type(it.overlap, graph_);

      counter++;
      auto tail_seq_id = sequence_to_node[it.overlap.lhs_id];
      auto head_seq_id = sequence_to_node[it.overlap.rhs_id];

      if (tail_seq_id == -1 || head_seq_id == -1) {
        continue;
      }
      auto tail = graph_.nodes_[sequence_to_node[it.overlap.lhs_id]].get();
      auto head = graph_.nodes_[sequence_to_node[it.overlap.rhs_id] + 1 - it.overlap.strand].get();

      auto length = it.overlap.lhs_begin - it.overlap.rhs_begin;
      auto length_pair =
        (graph_.piles_[it.overlap.rhs_id]->length() - it.overlap.rhs_end) -
          (graph_.piles_[it.overlap.lhs_id]->length() - it.overlap.lhs_end);

      if (it.overlap.score == 4) {
        std::swap(head, tail);
        length *= -1;
        length_pair *= -1;
      }

      auto edge = std::make_shared<Edge>(tail, head, length);
      graph_.edges_.emplace_back(edge);
      graph_.edges_.emplace_back(std::make_shared<Edge>(head->pair, tail->pair, length_pair));  // NOLINT
      edge->pair = graph_.edges_.back().get();
      edge->pair->pair = edge.get();

    }
  }

  std::cerr << "[raven::Graph::ConstructOverlapGraph] stored " << graph_.edges_.size() << " edges "  // NOLINT
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  graph_.PrintGfa(param.gfa_after_overlap_graph_construction_filename, false);  

  };

void Graph_Constructor::ConstructOnlyBackboneGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::vector<std::vector<extended_overlap>>& overlaps,
    biosoup::Timer& timer,
    Program_Parameters& param) {

  std::ofstream outdata_invalid("invalid_reads_backbone_only.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_invalid()) {
      outdata_invalid << sequences[i]->name << std::endl;
    }
  }

  std::ofstream outdata_contained("contained_reads_backbone_only.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_strong_contained()) {
      outdata_contained << sequences[i]->name << std::endl;
    }
  }

  graph_.nodes_.clear();
  graph_.edges_.clear();

  Node::num_objects = 0;
  Edge::num_objects = 0;

  const std::uint32_t n_reads = graph_.piles_.size();

  std::vector<bool> read_is_backbone(n_reads, false);
  for (std::uint32_t i = 0; i < n_reads; ++i) {
    read_is_backbone[i] =
        !graph_.piles_[i]->is_invalid() &&
        !graph_.piles_[i]->is_strong_contained() &&
        graph_.piles_[i]->is_backbone();
  }

  timer.Start();

  // --------------------------------------------------
  // PASS 1:
  // Create nodes only for backbone reads
  // --------------------------------------------------
  std::vector<std::int32_t> sequence_to_node(graph_.piles_.size(), -1);

  for (const auto& it : graph_.piles_) {
    std::uint32_t id = it->id();

    if (!read_is_backbone[id]) {
      continue;
    }

    std::unordered_set<std::uint32_t> annotations;
    for (const auto& jt : graph_.annotations_[id]) {
      if (it->begin() <= jt && jt < it->end()) {
        annotations.emplace(jt - it->begin());
      }
    }
    graph_.annotations_[id].swap(annotations);

    auto sequence = biosoup::NucleicAcid{
        sequences[id]->name,
        sequences[id]->InflateData(it->begin(), it->end() - it->begin())};
    sequence.id = id;

    sequence_to_node[id] = Node::num_objects;

    auto node = std::make_shared<Node>(sequence);
    sequence.ReverseAndComplement();
    graph_.nodes_.emplace_back(node);
    graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
    node->pair = graph_.nodes_.back().get();
    node->pair->pair = node.get();

    if (id < param.split) {
      node->color = 1;
      node->pair->color = 1;
    }
  }

  std::cerr << "[raven::Graph::ConstructOnlyBackboneGraph] stored "
            << graph_.nodes_.size() << " nodes "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  timer.Start();

  using edge_check = std::pair<std::uint32_t, std::uint32_t>;

  struct EdgeHash {
    std::size_t operator()(const edge_check& e) const {
      return std::hash<std::uint32_t>()(e.first) ^
             (std::hash<std::uint32_t>()(e.second) << 1);
    }
  };

  std::unordered_set<edge_check, EdgeHash> used_edges;

  // Optional dump before edge construction
 // graph_.PrintOverlaps(overlaps, sequences, true, param.paf_before_parsing_edges_filename);

  // --------------------------------------------------
  // PASS 2:
  // Add only strong backbone overlaps
  // --------------------------------------------------
  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    for (auto& ext_ovlp : overlaps[i]) {
      if (ext_ovlp.cat != overlapCategory::Backbone) {
        continue;
      }

      auto& ovlp = ext_ovlp.overlap;
      std::uint32_t lhs_id = ovlp.lhs_id;
      std::uint32_t rhs_id = ovlp.rhs_id;

      if (!read_is_backbone[lhs_id] || !read_is_backbone[rhs_id]) {
        continue;
      }

      std::uint32_t type = ext_ovlp.graph_overlap_type;
      if (type <= 2) {
          type = overlap_type(ovlp, graph_);
          ext_ovlp.graph_overlap_type = type;
      }
      ext_ovlp.overlap.score = type;

      // for label-2 HOR overlaps, don't drop on contained geometry
      if (type <= 2 && ext_ovlp.classification_label != 2) {
          continue;
      }

      auto tail_seq_id = sequence_to_node[lhs_id];
      auto head_seq_id = sequence_to_node[rhs_id];

      if (tail_seq_id == -1 || head_seq_id == -1) {
        continue;
      }

      std::uint32_t a = std::min(tail_seq_id, head_seq_id);
      std::uint32_t b = std::max(tail_seq_id, head_seq_id);
      edge_check edge_ids = {a, b};

      if (used_edges.count(edge_ids)) {
        continue;
      }
      used_edges.insert(edge_ids);

      auto tail = graph_.nodes_[sequence_to_node[lhs_id]].get();
      auto head = graph_.nodes_[sequence_to_node[rhs_id] + 1 - ovlp.strand].get();

      auto length = ovlp.lhs_begin - ovlp.rhs_begin;
      auto length_pair =
          (graph_.piles_[rhs_id]->length() - ovlp.rhs_end) -
          (graph_.piles_[lhs_id]->length() - ovlp.lhs_end);

      if (type == 4) {
        std::swap(head, tail);
        length *= -1;
        length_pair *= -1;
      }

      auto edge = std::make_shared<Edge>(tail, head, length);
      graph_.edges_.emplace_back(edge);
      graph_.edges_.emplace_back(
          std::make_shared<Edge>(head->pair, tail->pair, length_pair));
      edge->pair = graph_.edges_.back().get();
      edge->pair->pair = edge.get();
    }
  }

  std::cerr << "[raven::Graph::ConstructOnlyBackboneGraph] stored "
            << graph_.edges_.size() << " edges "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  graph_.PrintGfa("backbone_only.gfa", false);
}

void Graph_Constructor::ConstructBackboneGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::vector<std::vector<extended_overlap>>& overlaps,
    biosoup::Timer& timer,
    Program_Parameters& param) {

  const std::uint32_t weak_k = 1;  // try 1 first, then 2 if needed

  std::ofstream outdata_invalid("invalid_reads_backbone.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_invalid()) {
      outdata_invalid << sequences[i]->name << std::endl;
    }
  }

  std::ofstream outdata_contained("contained_reads_backbone.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_strong_contained()) {
      outdata_contained << sequences[i]->name << std::endl;
    }
  }

  graph_.nodes_.clear();
  graph_.edges_.clear();

  Node::num_objects = 0;
  Edge::num_objects = 0;

  const std::uint32_t n_reads = graph_.piles_.size();

  std::vector<bool> read_is_backbone(n_reads, false);
  for (std::uint32_t i = 0; i < n_reads; ++i) {
    read_is_backbone[i] =
        !graph_.piles_[i]->is_invalid() &&
        !graph_.piles_[i]->is_strong_contained() &&
        graph_.piles_[i]->is_backbone();
  }

  // --------------------------------------------------
  // PASS 0:
  // Reset candidate flags
  // --------------------------------------------------
  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    for (auto& ext_ovlp : overlaps[i]) {
      ext_ovlp.candidate_overlap = false;
    }
  }

  timer.Start();

  // --------------------------------------------------
  // PASS 1:
  // Create nodes only for backbone reads
  // --------------------------------------------------
  std::vector<std::int32_t> sequence_to_node(graph_.piles_.size(), -1);

  for (const auto& it : graph_.piles_) {
    std::uint32_t id = it->id();

    if (!read_is_backbone[id]) {
      continue;
    }

    std::unordered_set<std::uint32_t> annotations;
    for (const auto& jt : graph_.annotations_[id]) {
      if (it->begin() <= jt && jt < it->end()) {
        annotations.emplace(jt - it->begin());
      }
    }
    graph_.annotations_[id].swap(annotations);

    auto sequence = biosoup::NucleicAcid{
        sequences[id]->name,
        sequences[id]->InflateData(it->begin(), it->end() - it->begin())};
    sequence.id = id;

    sequence_to_node[id] = Node::num_objects;

    auto node = std::make_shared<Node>(sequence);
    sequence.ReverseAndComplement();
    graph_.nodes_.emplace_back(node);
    graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
    node->pair = graph_.nodes_.back().get();
    node->pair->pair = node.get();

    if (id < param.split) {
      node->color = 1;
      node->pair->color = 1;
    }
  }

  std::cerr << "[raven::Graph::ConstructBackboneGraph] stored "
            << graph_.nodes_.size() << " nodes "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  timer.Start();

  using edge_check = std::pair<std::uint32_t, std::uint32_t>;

  struct EdgeHash {
    std::size_t operator()(const edge_check& e) const {
      return std::hash<std::uint32_t>()(e.first) ^
             (std::hash<std::uint32_t>()(e.second) << 1);
    }
  };

  std::unordered_set<edge_check, EdgeHash> used_edges;

  // --------------------------------------------------
  // PASS 2:
  // Determine which backbone reads already have strong edges
  // in each direction.
  //
  // type == 3:
  //   lhs uses dir3, rhs uses dir4
  //
  // type == 4:
  //   lhs uses dir4, rhs uses dir3
  // --------------------------------------------------
  std::vector<bool> has_strong_dir3(n_reads, false);
  std::vector<bool> has_strong_dir4(n_reads, false);

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    for (auto& ext_ovlp : overlaps[i]) {
      if (ext_ovlp.cat != overlapCategory::Backbone) {
        continue;
      }

      auto& ovlp = ext_ovlp.overlap;
      std::uint32_t lhs_id = ovlp.lhs_id;
      std::uint32_t rhs_id = ovlp.rhs_id;

      if (!read_is_backbone[lhs_id] || !read_is_backbone[rhs_id]) {
        continue;
      }

      std::uint32_t type = ext_ovlp.graph_overlap_type;
      if (type <= 2) {
        type = overlap_type(ovlp, graph_);
        ext_ovlp.graph_overlap_type = type;
      }
      ext_ovlp.overlap.score = type;

      if (type == 3) {
        has_strong_dir3[lhs_id] = true;
        has_strong_dir4[rhs_id] = true;
      } else if (type == 4) {
        has_strong_dir4[lhs_id] = true;
        has_strong_dir3[rhs_id] = true;
      }
    }
  }

  // Optional dump before edge construction
  graph_.PrintOverlaps(overlaps, sequences, true, param.paf_before_parsing_edges_filename);

  auto add_edge_if_possible = [&](extended_overlap& ext_ovlp) {
    auto& ovlp = ext_ovlp.overlap;

    std::uint32_t lhs_id = ovlp.lhs_id;
    std::uint32_t rhs_id = ovlp.rhs_id;

    if (!read_is_backbone[lhs_id] || !read_is_backbone[rhs_id]) {
      return;
    }

    auto tail_seq_id = sequence_to_node[lhs_id];
    auto head_seq_id = sequence_to_node[rhs_id];

    if (tail_seq_id == -1 || head_seq_id == -1) {
      return;
    }

    std::uint32_t a = std::min(tail_seq_id, head_seq_id);
    std::uint32_t b = std::max(tail_seq_id, head_seq_id);
    edge_check edge_ids = {a, b};

    if (used_edges.count(edge_ids)) {
      return;
    }
    used_edges.insert(edge_ids);

    std::uint32_t type = ext_ovlp.graph_overlap_type;
    if (type <= 2) {
      type = overlap_type(ovlp, graph_);
      ext_ovlp.graph_overlap_type = type;
    }
    ext_ovlp.overlap.score = type;

    auto tail = graph_.nodes_[sequence_to_node[lhs_id]].get();
    auto head = graph_.nodes_[sequence_to_node[rhs_id] + 1 - ovlp.strand].get();

    auto length = ovlp.lhs_begin - ovlp.rhs_begin;
    auto length_pair =
        (graph_.piles_[rhs_id]->length() - ovlp.rhs_end) -
        (graph_.piles_[lhs_id]->length() - ovlp.lhs_end);

    if (type == 4) {
      std::swap(head, tail);
      length *= -1;
      length_pair *= -1;
    }

    auto edge = std::make_shared<Edge>(tail, head, length);
    graph_.edges_.emplace_back(edge);
    graph_.edges_.emplace_back(
        std::make_shared<Edge>(head->pair, tail->pair, length_pair));
    edge->pair = graph_.edges_.back().get();
    edge->pair->pair = edge.get();
  };

  // --------------------------------------------------
  // PASS 3:
  // Add all strong backbone overlaps first
  // --------------------------------------------------
  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    for (auto& ext_ovlp : overlaps[i]) {
      if (ext_ovlp.cat != overlapCategory::Backbone) {
        continue;
      }

      std::uint32_t type = ext_ovlp.graph_overlap_type;
      if (type <= 2) {
        type = overlap_type(ext_ovlp.overlap, graph_);
        ext_ovlp.graph_overlap_type = type;
      }
      ext_ovlp.overlap.score = type;

      if (type <= 2) {
        continue;
      }

      add_edge_if_possible(ext_ovlp);
    }
  }

  // --------------------------------------------------
  // PASS 4:
  // Mark and collect candidate weak-backbone overlaps
  // per missing direction.
  // --------------------------------------------------
  std::vector<std::vector<extended_overlap*>> weak_dir3_candidates(n_reads);
  std::vector<std::vector<extended_overlap*>> weak_dir4_candidates(n_reads);

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    for (auto& ext_ovlp : overlaps[i]) {
      if (ext_ovlp.cat != overlapCategory::WeakBackbone) {
        continue;
      }

      auto& ovlp = ext_ovlp.overlap;
      std::uint32_t lhs_id = ovlp.lhs_id;
      std::uint32_t rhs_id = ovlp.rhs_id;

      if (!read_is_backbone[lhs_id] || !read_is_backbone[rhs_id]) {
        continue;
      }

      std::uint32_t type = ext_ovlp.graph_overlap_type;
      if (type <= 2) {
        type = overlap_type(ovlp, graph_);
        ext_ovlp.graph_overlap_type = type;
      }
      ext_ovlp.overlap.score = type;

      if (type <= 2) {
        continue;
      }

      bool is_candidate = false;

      if (type == 3) {
        // lhs uses dir3, rhs uses dir4
        if (!has_strong_dir3[lhs_id]) {
          weak_dir3_candidates[lhs_id].push_back(&ext_ovlp);
          is_candidate = true;
        }
        if (!has_strong_dir4[rhs_id]) {
          weak_dir4_candidates[rhs_id].push_back(&ext_ovlp);
          is_candidate = true;
        }
      } else if (type == 4) {
        // lhs uses dir4, rhs uses dir3
        if (!has_strong_dir4[lhs_id]) {
          weak_dir4_candidates[lhs_id].push_back(&ext_ovlp);
          is_candidate = true;
        }
        if (!has_strong_dir3[rhs_id]) {
          weak_dir3_candidates[rhs_id].push_back(&ext_ovlp);
          is_candidate = true;
        }
      }

      ext_ovlp.candidate_overlap = is_candidate;
    }
  }

  auto overlap_len = [](const extended_overlap* e) -> std::uint32_t {
    const auto& ov = e->overlap;
    return std::max(ov.lhs_end - ov.lhs_begin, ov.rhs_end - ov.rhs_begin);
  };

  auto sort_and_truncate = [&](std::vector<extended_overlap*>& v) {
    std::sort(v.begin(), v.end(),
              [&](const extended_overlap* a, const extended_overlap* b) {
                return overlap_len(a) > overlap_len(b);
              });

    // remove duplicate pointers
    v.erase(std::unique(v.begin(), v.end()), v.end());

    if (v.size() > weak_k) {
      v.resize(weak_k);
    }
  };

  for (std::uint32_t i = 0; i < n_reads; ++i) {
    sort_and_truncate(weak_dir3_candidates[i]);
    sort_and_truncate(weak_dir4_candidates[i]);
  }

  // --------------------------------------------------
  // PASS 5:
  // Add selected weak-backbone overlaps
  // --------------------------------------------------
  for (std::uint32_t i = 0; i < n_reads; ++i) {
    for (auto* ext_ovlp : weak_dir3_candidates[i]) {
      add_edge_if_possible(*ext_ovlp);
    }
    for (auto* ext_ovlp : weak_dir4_candidates[i]) {
      add_edge_if_possible(*ext_ovlp);
    }
  }

  std::cerr << "[raven::Graph::ConstructBackboneGraph] stored "
            << graph_.edges_.size() << " edges "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  graph_.PrintGfa("backbone.gfa", false);
}

void Graph_Constructor::ConstructAssemblyGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                               std::vector<std::vector<extended_overlap>> &overlaps,
                                               biosoup::Timer &timer,
                                               Program_Parameters &param) {

  std::ofstream outdata_invalid;
  outdata_invalid.open("invalid_reads.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_invalid()) {
      outdata_invalid << sequences[i]->name << std::endl;
    }
  }

  std::ofstream outdata_contained;
  outdata_contained.open("contained_reads.txt");
  for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
    if (graph_.piles_[i]->is_contained()) {
      outdata_contained << sequences[i]->name << std::endl;
    }
  }
  graph_.nodes_.clear();
  graph_.edges_.clear();

  //exit(0);
  Node::num_objects = 0;
  Edge::num_objects = 0;

  std::vector<std::int32_t> sequence_to_node(graph_.piles_.size(), -1);
  for (const auto &it : graph_.piles_) {  // create nodes
    if (it->is_invalid() || it->is_contained()) {
      continue;
    }

    bool any_edge = false;

    // for (int j = 0; j < (int)overlaps[it->id()].size(); j++) {
    //   if (overlap_type(overlaps[it->id()][j].overlap, graph_) > 2) {
    //     if (!graph_.piles_[overlaps[it->id()][j].overlap.rhs_id]->is_invalid()) {
    //       any_edge = true;
    //       break;
    //     }
    //   }
    // }

    // if (!any_edge) {
    //   continue;
    // }

    std::unordered_set<std::uint32_t> annotations;
    for (const auto &jt : graph_.annotations_[it->id()]) {
      if (it->begin() <= jt && jt < it->end()) {
        annotations.emplace(jt - it->begin());
      }
    }
    graph_.annotations_[it->id()].swap(annotations);

    auto sequence = biosoup::NucleicAcid{
      sequences[it->id()]->name,
      sequences[it->id()]->InflateData(it->begin(), it->end() - it->begin()) };  // NOLINT
    sequence.id = it->id();

    sequence_to_node[it->id()] = Node::num_objects;

    auto node = std::make_shared<Node>(sequence);
    sequence.ReverseAndComplement();
    graph_.nodes_.emplace_back(node);
    graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
    node->pair = graph_.nodes_.back().get();
    node->pair->pair = node.get();

    if (it->id() < param.split) {
      node->color = 1;
      node->pair->color = 1;
    }
  }

  std::cerr << "[raven::Graph::ConstructAssemblyGraph] stored " << graph_.nodes_.size() << " nodes "  // NOLINT
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  timer.Start();

  int counter = 0;
  using edge_check = std::pair<std::uint32_t, std::uint32_t>;

  struct EdgeHash {
      std::size_t operator()(const edge_check& e) const {
          return std::hash<std::uint32_t>()(e.first) ^ (std::hash<std::uint32_t>()(e.second) << 1);
      }
  };
  std::unordered_set<edge_check, EdgeHash> used_edges;

  graph_.PrintOverlaps(overlaps, sequences, true, param.paf_before_parsing_edges_filename);

  for (int i = 0; i < (int)overlaps.size(); i++) {

    for (auto &it : overlaps[i]) {  // create edges

      // if (!overlap_finalize(it.overlap, graph_)) {
      //   std::cerr << "Error in overlap finalization" << std::endl;
      //   continue;
      // }
      it.overlap.score = overlap_type(it.overlap, graph_);

      // if (MLOverlapResolve(it) != 1){
      //   continue;
      // }
      // if(it.ground_truth == false){
      //   continue;
      // }
      counter++;
      auto tail_seq_id = sequence_to_node[it.overlap.lhs_id];
      auto head_seq_id = sequence_to_node[it.overlap.rhs_id];

      if (tail_seq_id == -1 || head_seq_id == -1) {
        continue;
      }
      
      std::uint32_t a = std::min(tail_seq_id, head_seq_id);
      std::uint32_t b = std::max(tail_seq_id, head_seq_id);
      edge_check edge_ids = {a, b};

      if (used_edges.count(edge_ids)) {
          continue;
      } else {
          used_edges.insert(edge_ids);
      }


      auto tail = graph_.nodes_[sequence_to_node[it.overlap.lhs_id]].get();
      auto head = graph_.nodes_[sequence_to_node[it.overlap.rhs_id] + 1 - it.overlap.strand].get();

      auto length = it.overlap.lhs_begin - it.overlap.rhs_begin;
      auto length_pair =
        (graph_.piles_[it.overlap.rhs_id]->length() - it.overlap.rhs_end) -
          (graph_.piles_[it.overlap.lhs_id]->length() - it.overlap.lhs_end);

      if (it.overlap.score == 4) {
        std::swap(head, tail);
        length *= -1;
        length_pair *= -1;
      }

      auto edge = std::make_shared<Edge>(tail, head, length);
      graph_.edges_.emplace_back(edge);
      graph_.edges_.emplace_back(std::make_shared<Edge>(head->pair, tail->pair, length_pair));  // NOLINT
      edge->pair = graph_.edges_.back().get();
      edge->pair->pair = edge.get();

    }
  }

  std::cerr << "[raven::Graph::ConstructAssemblyGraph] stored " << graph_.edges_.size() << " edges "  // NOLINT
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  graph_.PrintGfa(param.gfa_after_construction_filename, false);
}

void Graph_Constructor::PrintPiles(const std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences) {

  std::ofstream outdata;
  outdata.open("piles.csv");
  std::cout << "Writing pile data to piles.csv" << std::endl;
  for (int i = 0; i < (int)graph_.piles_.size(); i++) {
    outdata << sequences[i].get()->name << ",";
    if (graph_.piles_[i]->get_data().empty()) {
      continue;
    }
    std::vector<uint16_t> coverages = graph_.piles_[i]->get_data();
    for (auto &element : coverages) {
      outdata << element << ";";
    }
    outdata << std::endl;

  }
}

void Graph_Constructor::LoadGTOverlaps(const std::string &overlaps_path,
                                        std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                        std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                        bool load_cigar){
  std::ifstream file(overlaps_path);
  if (!file.is_open()) {
    throw std::runtime_error("Error opening file: " + overlaps_path);
  }

  std::map<std::string, std::uint32_t> sequence_name_to_seq_id;
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    sequence_name_to_seq_id[sequences[i]->name] = sequences[i]->id;
  }

  std::cerr << "[raven::Graph::LoadHerroSNPs] loading overlaps from: " << overlaps_path << std::endl;
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string item;
    std::vector<std::string> items;
    // std::uint32_t lhs_seq_id;
    // std::uint32_t rhs_seq_id;

    while (std::getline(iss, item, '\t')) {
      items.push_back(item);
    };

    // lhs_seq_id = get_read_id(items[0], sequences);
    auto lhs_seq_id = sequence_name_to_seq_id.find(items[0]);
    // rhs_seq_id = get_read_id(items[5], sequences);
    auto rhs_seq_id = sequence_name_to_seq_id.find(items[5]);

    if (lhs_seq_id == sequence_name_to_seq_id.end() || rhs_seq_id == sequence_name_to_seq_id.end()) {
      continue;
    } else {
      biosoup::Overlap overlap{ lhs_seq_id->second, (std::uint32_t)std::stoi(items[2]), (std::uint32_t)std::stoi(items[3]) > sequences[lhs_seq_id->second]->inflated_len ? sequences[lhs_seq_id->second]->inflated_len : (std::uint32_t)std::stoi(items[3]),
                                rhs_seq_id->second, (std::uint32_t)std::stoi(items[7]), (std::uint32_t)std::stoi(items[8]) > sequences[rhs_seq_id->second]->inflated_len ? sequences[rhs_seq_id->second]->inflated_len : (std::uint32_t)std::stoi(items[8]),
                                255, items[4] == "+" ? true : false};
      edlib_align tmp;
      if (load_cigar) {
        std::stringstream ss(items[16]);
        std::string segment;
        std::vector<std::string> seglist;
        while (std::getline(ss, segment, ':')) {
          seglist.push_back(segment);
        }
        tmp = { 0, 0, seglist[2], 0 };
      } else {
        tmp = {};
      }
      extended_overlap total_ovlp{};
      total_ovlp.overlap = overlap;
      total_ovlp.edlib_alignment = tmp;
      total_ovlp.total_overlap_snps = (std::uint32_t)std::stoi(items[11]);
      total_ovlp.total_overlap_snp_mismatches = (std::uint32_t)std::stoi(items[12]);
      total_ovlp.identity = (float)std::stoi(items[13]);
      total_ovlp.heterozygosity_rate = (float)std::stoi(items[14]);
      total_ovlp.graph_overlap_type = (std::uint8_t)std::stoi(items[15]);
      //total_ovlp.ol_class = items[24] == "same_strand" ? 1 : 0;
      total_ovlp.ol_class = 1;
      total_ovlp.ground_truth = (std::uint8_t)std::stoi(items[16]) == 1 ? true : false;

      // total_ovlp.total_overlap_snps = (std::uint32_t)std::stoi(items[17]);
      // total_ovlp.total_overlap_snp_mismatches = (std::uint32_t)std::stoi(items[18]);
      // total_ovlp.identity = (float)std::stoi(items[19]);
      // total_ovlp.heterozygosity_rate = (float)std::stoi(items[20]);
      // total_ovlp.graph_overlap_type = (std::uint8_t)std::stoi(items[21]);
      // total_ovlp.ol_type = OverlapType::perfect_heterozygous_high_match;
      // total_ovlp.ol_class = items[24] == "same_strand" ? 1 : 0;
      // total_ovlp.ground_truth = (std::uint8_t)std::stoi(items[33]) == 1 ? true : false;
    // extended_overlap total_ovlp{overlap, tmp, 
    //                             (std::uint32_t)std::stoi(items[19]), (std::uint32_t)std::stoi(items[20]),
    //                             (std::uint32_t)std::stoi(items[21]), (std::uint32_t)std::stoi(items[20]),
    //                             (std::uint32_t)std::stoi(items[18]), OverlapType::perfect_heterozygous_high_match,
    //                             items[24] == "same_strand" ? 1 : 0};
      extended_overlaps[lhs_seq_id->second].emplace_back(total_ovlp);
    }
  }
  std::cerr << "[raven::Graph::LoadHerroSNPs] loaded overlaps from: " << overlaps_path << std::endl;
  graph_.PrintOverlaps(extended_overlaps, sequences, false, "gt.paf");
  std::ofstream piles_tmp("piles_2.txt");
  for (const auto &it : graph_.piles_) {
    piles_tmp << it->id() << "\t" 
              << sequences[it->id()]->name << "\t"
              << it->get_data().size() << "\t"
              << it->length() << std::endl;
  }
  std::vector<std::future<void>> extended_layers_futures;
  std::uint16_t counter = 0;
  for (const auto &it : graph_.piles_) {
      counter += 1;
      auto id = it->id();

      if (id >= extended_overlaps.size() || extended_overlaps[id].empty()) {
          continue;
      }

      auto overlaps_copy = extended_overlaps[id];
      auto pile_data = it->get_data();  // copy of the pointer (deep copy)

      extended_layers_futures.emplace_back(thread_pool_->Submit(
          [&it, overlaps_copy]() -> void {
              it->AddExtendedLayers(
                  overlaps_copy.begin(),
                  overlaps_copy.end());
          }));
      // //auto data = it->get_data();
  }


  for (const auto &it : extended_layers_futures) {
    it.wait();
  }

  extended_layers_futures.clear();
};

void Graph_Constructor::LoadOverlaps(const std::string &overlaps_path,
                                     std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                     std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                     bool load_cigar) {
  std::ifstream file(overlaps_path);
  if (!file.is_open()) {
    throw std::runtime_error("Error opening file: " + overlaps_path);
  }

  std::map<std::string, std::uint32_t> sequence_name_to_seq_id;
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    sequence_name_to_seq_id[sequences[i]->name] = sequences[i]->id;
  }

  std::cerr << "[raven::Graph::LoadHerroSNPs] loading overlaps from: " << overlaps_path << std::endl;
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string item;
    std::vector<std::string> items;
    std::uint32_t lhs_seq_id;
    std::uint32_t rhs_seq_id;

    while (std::getline(iss, item, '\t')) {
      items.push_back(item);
    };

    // lhs_seq_id = get_read_id(items[0], sequences);
    lhs_seq_id = sequence_name_to_seq_id[items[0]];
    // rhs_seq_id = get_read_id(items[5], sequences);
    rhs_seq_id = sequence_name_to_seq_id[items[5]];

    if (lhs_seq_id == (std::uint32_t )-1 || rhs_seq_id == (std::uint32_t )-1) {
      continue;
    } else {
      biosoup::Overlap overlap{ lhs_seq_id, (std::uint32_t)std::stoi(items[2]), (std::uint32_t)std::stoi(items[3]),
                                rhs_seq_id, (std::uint32_t)std::stoi(items[7]), (std::uint32_t)std::stoi(items[8]),
                                255, items[4] == "+" ? true : false};
      edlib_align tmp;
      if (load_cigar) {
        std::stringstream ss(items[16]);
        std::string segment;
        std::vector<std::string> seglist;
        while (std::getline(ss, segment, ':')) {
          seglist.push_back(segment);
        }
        tmp = { 0, 0, seglist[2], 0 };
      } else {
        tmp = {};
      }
      extended_overlap total_ovlp{ overlap, tmp, 0, 0 };
      extended_overlaps[lhs_seq_id].emplace_back(total_ovlp);
    }
  }
  std::cerr << "[raven::Graph::LoadHerroSNPs] loaded overlaps from: " << overlaps_path << std::endl;
}

void Graph_Constructor::LoadHerroSNPs(const std::string &herro_snps_path,
                                      std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences) {
  std::cerr << "Loading snps" << std::endl;

  std::ifstream file(herro_snps_path);
  if (!file.is_open()) {
    throw std::runtime_error("Error opening file: " + herro_snps_path);
  }

  std::cerr << "[raven::Graph::LoadHerroSNPs] loading snps from: " << herro_snps_path << std::endl;
  std::string line;
  std::map<std::string, std::uint32_t> sequence_name_to_seq_id;

  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    sequence_name_to_seq_id[sequences[i]->name] = sequences[i]->id;
  }

  while (std::getline(file, line)) {
    std::string single_line = line;
    std::istringstream iss(single_line);

    std::string item;
    std::uint32_t seq_id;
    bool found = false;

    std::uint32_t last_id;
    std::vector<std::string> elements(3);
    std::vector<std::string> elements2(2);

    for (int i = 0; i < 3 && std::getline(iss, item, '\t'); ++i) {
      elements[i] = item;
    }
    if (elements[0] == elements2[0] && elements[1] == elements2[1]) {
      seq_id = last_id;
      found = true;
    } else {
      auto it = sequence_name_to_seq_id.find(elements[0]);
      if (it != sequence_name_to_seq_id.end()) {
        seq_id = it->second;
        found = true;
      } else {
        auto it = sequence_name_to_seq_id.find(elements[0] + ":" + elements[1]);
        if (it != sequence_name_to_seq_id.end()) {
          seq_id = it->second;
          found = true;
        }
      }
    }
    if (found) {
      last_id = seq_id;
      std::uint32_t pos = std::stoi(elements[2]);
      graph_.annotations_[seq_id].emplace(pos);
      elements2[0] = elements[0];
      elements2[1] = elements[1];
      found = false;
    }
  }

  std::ofstream outdata;
  outdata.open("snp_annotations_check.anno");
  for (std::uint32_t i = 0; i < graph_.annotations_.size(); ++i) {
    if (graph_.annotations_[i].empty()) {
      continue;
    }
    outdata << sequences[i]->name << " ";
    for (const auto &jt : graph_.annotations_[i]) {
      outdata << " " << jt;
    }
    outdata << std::endl;
  }
}

void Graph_Constructor::LoadAnnotations(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                        std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                        Program_Parameters &param) {
  std::vector<std::future<void>> void_futures;
  for (int i = 0; i < (int)sequences.size(); i++) {
    void_futures.emplace_back(thread_pool_->Submit(
      [&](std::uint32_t i) -> void {
        call_snps(i, extended_overlaps[i], sequences, graph_);
      },
      i));
  }

  for (const auto &it : void_futures) {
    it.wait();
  }

  void_futures.clear();

  std::cerr << "[raven::Graph::Construct] snps called"
            << std::endl;

  if (!print_snp_data || param.ploidy < 2)
    return;

  std::ofstream outdata;
  outdata.open("snp_annotations.anno");

  for (std::uint32_t i = 0; i < graph_.annotations_.size(); ++i) {
    if (graph_.annotations_[i].empty()) {
      continue;
    }
    outdata << sequences[i]->name << " ";
    for (const auto &jt : graph_.annotations_[i]) {
      outdata << " " << jt;
    }
    outdata << std::endl;
  }
}

void Graph_Constructor::LoadFromGfa(const std::string &gfa_path) {
  try {
    std::string gfa_path_without_leading_whitespace;
    if (!gfa_path.empty()) {
      gfa_path_without_leading_whitespace = gfa_path.substr(1);
    }
    std::ifstream file(gfa_path_without_leading_whitespace);

    if (!file.is_open()) {
      throw std::runtime_error("Error opening file: " + gfa_path_without_leading_whitespace);
    }

    std::string line;
    Node::num_objects = 0;
    Edge::num_objects = 0;
    std::map<std::string, std::shared_ptr<Node>> sequence_to_node;

    while (std::getline(file, line)) {
      // Process each line here
      std::string single_line = line;
      std::istringstream iss(single_line);

      std::string item;
      std::string first_item;
      std::uint8_t counter = 0;
      std::uint32_t sequence_counter = 0;
      std::string seq_name;
      std::string nuc_sequence;

      while (std::getline(iss, item, '\t')) {
        if (counter == 0) {
          first_item = item;
          if (first_item == "S") {
            sequence_counter++;
            counter++;
            continue;
          } else {
            break;
          }
        }

        if (first_item == "S") {
          if (counter == 1) {
            seq_name = item;
          } else if (counter == 2) {
            nuc_sequence = item;

          }
        }
        counter++;
      }
      if (first_item == "S") {
        auto sequence = biosoup::NucleicAcid{
          seq_name,
          nuc_sequence
        };
        sequence.id = sequence_counter;

        auto node = std::make_shared<Node>(sequence);
        sequence.ReverseAndComplement();
        graph_.nodes_.emplace_back(node);
        graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
        node->pair = graph_.nodes_.back().get();
        node->pair->pair = node.get();
        sequence_to_node.emplace(seq_name, node);
      }
    }
    if (file.eof()) {
      std::cerr << "[raven::Graph::LoadFromGfa] loaded sequences from: " << gfa_path_without_leading_whitespace
                << std::endl;
    }

    std::ifstream file2(gfa_path_without_leading_whitespace);

    while (std::getline(file2, line)) {
      // Process each line here
      std::string single_line = line;
      std::istringstream iss(single_line);

      std::string item;
      std::string first_item;
      std::uint8_t counter = 0;
      std::uint32_t sequence_counter = 0;
      std::string seq_name;
      std::string nuc_sequence;
      std::string tail_node_name;
      std::string head_node_name;

      bool tail_node_strand;
      bool head_node_strand;

      std::string edge_length;
      std::string item2;
      std::string ol_length;

      while (std::getline(iss, item, '\t')) {
        if (counter == 0) {
          first_item = item;
          if (first_item == "L") {
            sequence_counter++;
            counter++;
            continue;
          } else {
            break;
          }
        }
        if (first_item == "L") {
          if (counter == 1) {
            tail_node_name = item;
          } else if (counter == 2) {
            tail_node_strand = item == "+" ? true : false;
          } else if (counter == 3) {
            head_node_name = item;
          } else if (counter == 4) {
            head_node_strand = item == "+" ? true : false;
          } else if (counter == 5) {
            std::stringstream ss(item);
            while (std::getline(ss, item2, 'M')) {
              ol_length = item2;
            }
          } else if (counter == 6) {
            std::uint8_t mini_counter = 0;
            std::stringstream ss(item);
            while (std::getline(ss, item2, ':')) {
              if (mini_counter == 2) edge_length = item2;
              mini_counter++;
            }
          }
        }
        counter++;
      }
      if (first_item == "L") {

        auto tail_node =
          tail_node_strand ? sequence_to_node[tail_node_name].get() : sequence_to_node[tail_node_name]->pair;
        auto head_node =
          head_node_strand ? sequence_to_node[head_node_name].get() : sequence_to_node[head_node_name]->pair;

        auto length = std::stoi(edge_length);
        auto length_pair = head_node->sequence.inflated_len - std::stoi(ol_length);

        auto edge = std::make_shared<Edge>(tail_node, head_node, length);
        graph_.edges_.emplace_back(edge);
        graph_.edges_.emplace_back(std::make_shared<Edge>(head_node->pair, tail_node->pair, length_pair));  // NOLINT
        edge->pair = graph_.edges_.back().get();
        edge->pair->pair = edge.get();

      }
    }
    //std::cout << line << std::endl;
    file2.close();
    if (file.eof()) {
      // File has been read successfully
      file.close();
      std::cerr << "[raven::Graph::LoadFromGfa] successfully loaded graph from: " << gfa_path_without_leading_whitespace
                << std::endl;
    } else {
      throw std::runtime_error("Error reading file: " + gfa_path_without_leading_whitespace);
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  graph_.state_manager_.set_state(GraphState::Assemble_Transitive_Edges);
}

void Graph_Constructor::LoadFromPaf(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                    const std::string &paf_path) {
  try {
    gzFile file = gzopen(paf_path.c_str(), "r");
    OverlapParser parser{ file };
    std::vector<std::unique_ptr<Overlap>> overlaps = parser.ParsePaf((std::uint64_t)-1);

    Node::num_objects = 0;
    Edge::num_objects = 0;
    std::map<std::string, std::shared_ptr<Node>> sequence_to_node;
    std::uint32_t sequence_counter = 0;

    for (std::unique_ptr<biosoup::NucleicAcid> &seq_ptr : sequences) {
      biosoup::NucleicAcid seq = *seq_ptr.get();
      seq.id = ++sequence_counter;
      std::shared_ptr<Node> node = std::make_shared<Node>(seq);
      graph_.nodes_.emplace_back(node);
      seq.ReverseAndComplement();
      graph_.nodes_.emplace_back(std::make_shared<Node>(seq));
      node->pair = graph_.nodes_.back().get();
      node->pair->pair = node.get();
      sequence_to_node.emplace(seq.name, node);
    }

    bool tail_node_strand;
    bool head_node_strand;

    for (std::unique_ptr<Overlap> &overlap_ptr : overlaps) {
      Overlap overlap = *overlap_ptr.get();
      tail_node_strand = true;
      head_node_strand = overlap.strand;

      Node
        *tail_node = tail_node_strand ? sequence_to_node[overlap.q_name].get() : sequence_to_node[overlap.q_name]->pair;
      Node
        *head_node = head_node_strand ? sequence_to_node[overlap.t_name].get() : sequence_to_node[overlap.t_name]->pair;

      uint32_t length = overlap.q_len - overlap.overlap_len;
      uint32_t length_pair = overlap.t_len - overlap.overlap_len;

      std::shared_ptr<Edge> edge = std::make_shared<Edge>(tail_node, head_node, length);
      graph_.edges_.emplace_back(edge);
      graph_.edges_.emplace_back(std::make_shared<Edge>(head_node->pair, tail_node->pair, length_pair));  // NOLINT
      edge->pair = graph_.edges_.back().get();
      edge->pair->pair = edge.get();
    }

  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  graph_.state_manager_.set_state(GraphState::Assemble_Transitive_Edges);
}
// NOLINT

} // raven
