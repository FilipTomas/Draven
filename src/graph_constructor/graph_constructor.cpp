
#include <cmath>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <cassert>
#include "graph.hpp"
#include "graph_constructor.h"
#include "overlap.h"
#include "overlap_parser.h"
#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "biosoup/timer.hpp"
#include "extended_overlap.h"
#include "overlap_helpers.h"
#include "annotation_helpers.h"
#include "ram_overlap.h"
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

std::vector<ram_overlap>
merge_fragments(const std::vector<biosoup::Overlap>& ovlps)
{
    std::vector<ram_overlap> result;
    result.reserve(ovlps.size());         // upper bound

    if (ovlps.empty()){
        return result;
    }

    // begin first group
    ram_overlap curr = start_group(ovlps.front());
    curr.n_fragments = 1;

    for (std::size_t i = 1; i < ovlps.size(); ++i) {
        const auto& ov = ovlps[i];

        bool same_pair =
            (ov.lhs_id == curr.lhs_id && ov.rhs_id == curr.rhs_id);

        if (!same_pair) {
            // commit the finished group
            result.push_back(std::move(curr));
            curr = start_group(ov);       // start new group
        } else {
            // extend current aggregate
            curr.overlap_fragments.push_back(extended_overlap {ov, {}, 0, 0});
            if (curr.lhs_begin > ov.lhs_begin){
                curr.lhs_begin = ov.lhs_begin;
            } else if (curr.lhs_end < ov.lhs_end) {
                curr.lhs_end = ov.lhs_end;
            };
            // curr.lhs_end = ov.lhs_end;    // ovlps are ordered by lhs_begin
            // curr.rhs_end = ov.rhs_end;
            curr.n_fragments++;
        }
    }
    result.push_back(std::move(curr));    // flush last group
    return result;
}

Graph_Constructor::Graph_Constructor(Graph &graph, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
  : graph_(graph), thread_pool_(thread_pool ?
                                thread_pool :
                                std::make_shared<thread_pool::ThreadPool>(1)) {
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

void Graph_Constructor::ConstructCorrectionOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, Program_Parameters &param){
    std::vector<std::vector<extended_overlap>> extended_overlaps(sequences.size());
    biosoup::Timer timer{};
    std::cerr << "Only here!" << std::endl;
    ConstructOverlapsForCorrection(sequences, extended_overlaps, timer, param);
    //LoadOverlapsAndAlign(param, sequences, extended_overlaps);
    AlignBackboneReads(sequences, extended_overlaps, timer, param);
    //ConstructOnlyBackboneGraph(sequences, extended_overlaps, timer, param);
    std::cerr << "[raven::Graph::Construct] reached checkpoint "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
    std::cerr << "Printing correction overlaps!" << std::endl;
    graph_.PrintCorrectionOverlaps(extended_overlaps, sequences, "correction.paf");
}

void Graph_Constructor::LoadOverlapsAndAlign(Program_Parameters &param, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps){
  std::ifstream file("/mnt/share1_Jabba/ftomas/centromere_fix/minimap2/minimap2_dedup.paf");
  if (!file.is_open()) {
    throw std::runtime_error("Error opening file: " + param.input_paf_path);
  }
  std::cerr << "Got to here" << std::endl;
  std::map<std::string, std::uint32_t> sequence_name_to_seq_id;
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    sequence_name_to_seq_id[sequences[i]->name] = sequences[i]->id;
  }

  std::cerr << "[raven::Graph::LoadHerroSNPs] loading overlaps from: " << param.input_paf_path << std::endl;
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
      edlib_align tmp = {};
      extended_overlap total_ovlp{ overlap, tmp, 0, 0 };
      total_ovlp.correction_overlap = true;
      extended_overlaps[lhs_seq_id].emplace_back(total_ovlp);
    }
  }
  std::cerr << "[raven::LoadAndAlign] loaded overlaps from: " << param.input_paf_path << std::endl;

  auto resolve_alignments = [&](std::uint32_t i) -> void {
    std::string lhs_full = sequences[i]->InflateData();
    for (auto& overlap : extended_overlaps[i]) {
      std::string lhs = lhs_full.substr(
          overlap.overlap.lhs_begin,
          overlap.overlap.lhs_end - overlap.overlap.lhs_begin);

      biosoup::NucleicAcid rhs_{ "",
                                sequences[overlap.overlap.rhs_id]->InflateData(
                                    overlap.overlap.rhs_begin,
                                    overlap.overlap.rhs_end - overlap.overlap.rhs_begin) };

      if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

      auto rhs = rhs_.InflateData();

      overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
    }
  };

  std::vector<std::future<void>> alignment_futures;
  alignment_futures.reserve(sequences.size());
  for (std::uint32_t i = 0; i < sequences.size(); ++i) {
    alignment_futures.emplace_back(thread_pool_->Submit(resolve_alignments, i));
  }
  for (const auto& it : alignment_futures) {
    it.wait();
  }
};

void Graph_Constructor::ConstructOverlapsForCorrection(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                       std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                       biosoup::Timer &timer,
                                       Program_Parameters &param) {

  graph_.annotations_.resize(sequences.size());
                                      
  for (const auto &it : sequences) {
    graph_.piles_.emplace_back(new Pile(it->id, it->inflated_len));
  }

  MapSequencesFast(sequences, extended_overlaps, timer, param);

  graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_initial_overlaps_filename);
  std::cerr << "[raven::Graph::Construct] initial overlaps printed"
            << std::endl;
  exit(0);

  ResolveContainedReads(sequences, extended_overlaps, timer, 3);

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
    "/mnt/share1_Jabba/ftomas/CNNONT_models/CNNONT_model_multi_uncapped.ts.pt";
  // /home/ftomas/github/CNNONT/mlruns/4/f1546b13a1e645cfa8e5d9340c82894d/artifacts/checkpoints/cnnont_model_best.pth
  at::set_num_threads(1);
  at::set_num_interop_threads(1);
  std::vector<std::vector<float>> all_features {};

  auto sketch_sequence = [&](std::uint32_t i) -> void {
    std::vector<std::uint64_t> ids;
    std::vector<float> sketch;
    std::vector<std::uint32_t> avg_base_qualities;
    std::vector<std::uint32_t> min_base_qualities;

    if (!param.minimizers && !param.no_fuse) {
      // Fused: also builds and stores this read's minimizer sketch.
      minimizer_engine.SketchAndMinimizeFastK(
          sequences[i], 1U, ids, sketch, avg_base_qualities, min_base_qualities);
    } else {
      minimizer_engine.FastKSketchReadInto(
          sequences[i], 1U, ids, sketch, avg_base_qualities, min_base_qualities);
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


    graph_.piles_[i]->set_k_kmer_ids(ids);


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

    std::vector<kMerType> expanded_kmer_types;
    expanded_kmer_types.reserve(ids.size());

    auto class_to_type = [](int64_t c) -> kMerType {
      switch (c) {
        case 0: return kMerType::Diploid;  // D
        case 1: return kMerType::Haploid;  // H
        case 2: return kMerType::Error;    // E
        default: return kMerType::None;
      }
    };

    // Expand each window prediction to 10 positions
    for (int64_t w = 0; w < L; ++w) {
      kMerType t = class_to_type(pred_acc[w]);

      for (std::uint32_t r = 0; r < 10; ++r) {
        if (expanded_kmer_types.size() >= ids.size()) break;
        expanded_kmer_types.push_back(t);
      }
    }

    // ---------------------------------------------------------
    // Tail handling: repeat LAST classification
    // ---------------------------------------------------------
    kMerType tail = expanded_kmer_types.empty()
                        ? kMerType::None
                        : expanded_kmer_types.back();

    while (expanded_kmer_types.size() < ids.size()) {
      expanded_kmer_types.push_back(tail);
    }

    graph_.piles_[i]->set_k_mer_types(expanded_kmer_types);

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
      graph_.piles_[pile_id]->append_kmers_in_region(p.first, p.second, kmers);
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

  auto map_sequences = [&](std::uint32_t i) -> std::vector<extended_overlap> { // map sequences
  //auto map_sequences = [&](std::uint32_t i) -> std::vector<std::vector<float>> { // map sequences
    std::vector<biosoup::Overlap> ovlps;
    if(graph_.piles_[i]->is_hor()){
     ovlps = minimizer_engine.MapRepetitive(sequences[i], true, true,
                                                               false);
    } else {
      ovlps = minimizer_engine.Map(sequences[i], true, true,
                                                                false);
    }
    
    if (!ovlps.empty()) {
      std::vector<extended_overlap> ovlps_final;
      std::vector<extended_overlap> ovlps_tmp;
      std::vector<biosoup::Overlap> no_dups_overlaps;
      std::vector<std::vector<float>> features {};
        
      std::sort(ovlps.begin(), ovlps.end(),
          [&](const biosoup::Overlap &lhs,
              const biosoup::Overlap &rhs) -> bool {
            return overlap_length(lhs) > overlap_length(rhs);
          });
      no_dups_overlaps = remove_duplicates(ovlps);

      std::vector<biosoup::Overlap> tmp;

      for (auto &ovlp : no_dups_overlaps) {

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
  std::unique_ptr<ram::KmerCounter> internal_counter;
  if (param.internal_count) {
    // Bloom sizing (bloom_capacity_hint is interpreted as a minimum bit count):
    //   --bloom-gb < 0  -> auto: size by total k-mer instances (safe upper
    //                      bound on distinct k-mers)
    //   --bloom-gb == 0 -> disable bloom: exact count + prune (no false
    //                      positives, holds singletons transiently)
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
  } else if(!param.minimizers){
    minimizer_engine.LoadFastK();
  }

  if (use_fused_fastk) {
    minimizer_engine.InitReadStore(sequences.size());
  }
  for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
    bytes += sequences[i]->inflated_len;
    if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
      continue;
    }
    bytes = 0;

    const std::uint32_t block_start = j;

    // ---- sketch first: feature extraction (+ fused minimizer construction,
    //      which populates the per-read minimizer store) ----
    timer.Start();
    std::cerr << "Starting sketching" << std::endl;
    std::vector<std::future<void>> sketch_futures;
    for (std::uint32_t k = block_start; k < i + 1; ++k) {
      sketch_futures.emplace_back(thread_pool_->Submit(sketch_sequence, k));
    }
    for (const auto &it : sketch_futures) {
      it.wait();
    }
    sketch_futures.clear();
    std::cerr << "[raven::Graph::Construct] sketched sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    graph_.PrintPiles();

    // ---- build the minimizer index (from stored sketches on the FastK path) ----
    timer.Start();
    if (use_fused_fastk) {
      minimizer_engine.BuildIndexFromStored(
        sequences.begin() + block_start,
        sequences.begin() + i + 1);
    } else {
      minimizer_engine.Minimize(
        sequences.begin() + block_start,
        sequences.begin() + i + 1,
        true);
    }

    minimizer_engine.Filter(param.freq);

    std::cerr << "[raven::Graph::Construct] minimized "
              << block_start << " - " << i + 1 << " / " << sequences.size() << " "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    std::vector<std::uint32_t> num_overlaps(extended_overlaps.size());
    for (std::uint32_t k = 0; k < extended_overlaps.size(); ++k) {
      num_overlaps[k] = extended_overlaps[k].size();
    }

    // ---- map ----
    std::vector<std::future<std::vector<extended_overlap>>> thread_futures;

    timer.Start();
    for (std::uint32_t k = block_start; k < i + 1; ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(map_sequences, k));

        bytes += sequences[k]->inflated_len;
        if (k != i && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        for (auto &it : thread_futures) {
          for (const auto &jt : it.get()) {
            extended_overlaps[jt.overlap.lhs_id].emplace_back(jt);
            extended_overlaps[jt.overlap.rhs_id].emplace_back(feature_overlap_reverse(jt));
          }
        }
        thread_futures.clear();
      std::cerr << "[raven::Graph::Construct] mapped sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }

    // ---- release this block's stored minimizers ----
    if (use_fused_fastk) {
      minimizer_engine.FreeReadStoreRange(
        sequences.begin() + block_start,
        sequences.begin() + i + 1);
    }

    j = i + 1;
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

void Graph_Constructor::ResolveContainedReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::vector<std::vector<extended_overlap>>& extended_overlaps,
    biosoup::Timer& timer,
    std::uint32_t containment_threshold) {

  timer.Start();

  // --------------------------------------------------
  // PASS 1:
  // Use only strong overlaps (classification_label == 1)
  // to detect containment and mark backbone reads.
  // Label-2 (HOR) overlaps are never marked as contained.
  // A read is only marked as contained if it receives
  // at least containment_threshold votes.
  // --------------------------------------------------
  std::vector<std::uint32_t> containment_votes(sequences.size(), 0);

  for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
    for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {

      auto& ext_ovlp = extended_overlaps[i][j];
      ext_ovlp.cat = overlapCategory::None;

      if (ext_ovlp.classification_label != 1 &&
          ext_ovlp.classification_label != 2) {
        continue;
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
      } else {
        // label == 1 — vote-based containment logic
        if (type == 1) {
          containment_votes[lhs_id]++;
          ext_ovlp.cat = overlapCategory::StrongContained;
        } else if (type == 2) {
          containment_votes[rhs_id]++;
          ext_ovlp.cat = overlapCategory::StrongContained;
        } else {
          ext_ovlp.cat = overlapCategory::Backbone;
          graph_.piles_[lhs_id]->set_backbone();
          graph_.piles_[rhs_id]->set_backbone();
        }
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

void Graph_Constructor::AlignBackboneReads(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                               std::vector<std::vector<extended_overlap>> &overlaps,
                                               biosoup::Timer &timer,
                                               Program_Parameters &param){

  auto resolve_alignments = [&](std::uint32_t i) -> void {
    if(graph_.piles_[i]->is_backbone() && !graph_.piles_[i]->is_strong_contained()){
      std::string lhs_full = sequences[i]->InflateData();
      for(auto& overlap : overlaps[i]){
        if(overlap.classification_label != 1){
          // std::string lhs_full = sequences[i]->InflateData();
          std::string lhs = lhs_full.substr(
                            overlap.lhs_begin_original,
                            overlap.lhs_end_original - overlap.lhs_begin_original);
          // biosoup::NucleicAcid lhs_{ "",
          //                           sequences[overlap.overlap.lhs_id]->InflateData(overlap.overlap.lhs_begin,
          //                                                                     overlap.overlap.lhs_end - overlap.overlap.lhs_begin)};
          biosoup::NucleicAcid rhs_{ "",
                                    sequences[overlap.overlap.rhs_id]->InflateData(overlap.rhs_begin_original,
                                                                              overlap.rhs_end_original - overlap.rhs_begin_original) };

          if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

          auto rhs = rhs_.InflateData();
      //    auto lhs = lhs_.InflateData();
          overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
          overlap.correction_overlap = true;
        } else if(overlap.classification_label == 1){
          std::string lhs = lhs_full.substr(
                            overlap.overlap.lhs_begin,
                            overlap.overlap.lhs_end - overlap.overlap.lhs_begin);
          // biosoup::NucleicAcid lhs_{ "",
          //                           sequences[overlap.overlap.lhs_id]->InflateData(overlap.overlap.lhs_begin,
          //                                                                     overlap.overlap.lhs_end - overlap.overlap.lhs_begin)};
          biosoup::NucleicAcid rhs_{ "",
                                    sequences[overlap.overlap.rhs_id]->InflateData(overlap.overlap.rhs_begin,
                                                                              overlap.overlap.rhs_end - overlap.overlap.rhs_begin) };

          if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

          auto rhs = rhs_.InflateData();
      //    auto lhs = lhs_.InflateData();
          overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
          overlap.correction_overlap = true;        
        } //else {
        //   std::string lhs = lhs_full.substr(
        //                     overlap.overlap.lhs_begin,
        //                     overlap.overlap.lhs_end - overlap.overlap.lhs_begin);
        //   biosoup::NucleicAcid rhs_{ "",
        //                             sequences[overlap.overlap.rhs_id]->InflateData(overlap.overlap.rhs_begin,
        //                                                                       overlap.overlap.rhs_end - overlap.overlap.rhs_begin) };

        //   if (!overlap.overlap.strand) rhs_.ReverseAndComplement();

        //   auto rhs = rhs_.InflateData();
        //   // overlap.edlib_alignment = edlib_wrapper(lhs, rhs);
        //   overlap.correction_overlap = true;
          // for(auto& rhs_overlap : overlaps[overlap.overlap.rhs_id]){
          //   if(rhs_overlap.overlap.rhs_id == i){
          //     rhs_overlap.edlib_alignment.cigar = cigar_alignment_reverse(overlap.edlib_alignment.cigar);
          //     rhs_overlap.edlib_alignment.matches = overlap.edlib_alignment.matches;
          //     rhs_overlap.edlib_alignment.block_length = overlap.edlib_alignment.block_length;
          //     rhs_overlap.edlib_alignment.edit_distance = overlap.edlib_alignment.edit_distance;
          //     rhs_overlap.correction_overlap = true;
          //     break;
          //   }
          // }
      //  }
      }
    }
  };

  std::vector<std::future<void>> alignment_futures;
  for(int i = 0; i < sequences.size(); i++){
    alignment_futures.emplace_back(thread_pool_->Submit(resolve_alignments, i));
    //resolve_alignments(i);
  }

  for (const auto &it : alignment_futures) {
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
