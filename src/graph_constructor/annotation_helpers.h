
#ifndef RAVEN_SRC_GRAPH_CONSTRUCTOR_ANNOTATION_HELPERS_H_
#define RAVEN_SRC_GRAPH_CONSTRUCTOR_ANNOTATION_HELPERS_H_

#include <fstream>
#include "../graph.hpp"
#include "overlap_helpers.h"
#include "edlib.h"

struct base_pile {
  std::uint32_t a;
  std::uint32_t c;
  std::uint32_t g;
  std::uint32_t t;
  std::uint32_t i;
  std::uint32_t d;
};

// NW alignment with path that gives up (editDistance -1) once the edit
// distance would exceed max_k; max_k < 0 means no limit. Like edlib's own
// k = -1 it starts at k = 64 and doubles, so an alignment under the cap costs
// and comes out the same (edlib computes the path from the final edit
// distance, not from k). Unrelated sequence costs ~0.5 edits per base and,
// uncapped, close to quadratic time.
inline EdlibAlignResult capped_nw_align(
  const char *query, int query_len,
  const char *target, int target_len,
  std::int64_t max_k) {
  if (max_k < 0) {
    return edlibAlign(query, query_len, target, target_len,
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  }
  for (std::int64_t k = std::min<std::int64_t>(64, max_k);; k = std::min(2 * k, max_k)) {
    EdlibAlignResult result = edlibAlign(query, query_len, target, target_len,
        edlibNewAlignConfig(static_cast<int>(k), EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
    if (result.status != EDLIB_STATUS_OK || result.editDistance >= 0 || k == max_k) {
      return result;
    }
    edlibFreeAlignResult(result);
  }
}

edlib_align edlib_wrapper(
  const std::string &lhs,
  const std::string &rhs,
  std::int64_t max_k = -1) {

  edlib_align alignment_result;

  EdlibAlignResult result = capped_nw_align(
    lhs.c_str(), lhs.size(),
    rhs.c_str(), rhs.size(),
    max_k); // align lhs and rhs

  if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0) {
    alignment_result.cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    alignment_result.edit_distance = result.editDistance;
    alignment_result.block_length = result.alignmentLength;
    alignment_result.matches = 0;
    for (int i = 0; i < result.alignmentLength; i++) {
      if (result.alignment[i] == 0) {
        ++alignment_result.matches;
      }
    }
    edlibFreeAlignResult(result);
    return alignment_result;
  } else {
    edlibFreeAlignResult(result);
    alignment_result.cigar = "";
    alignment_result.matches = 0;
    alignment_result.edit_distance = -1;
    alignment_result.block_length = 0;
    return alignment_result;
  }
}

// Global (NW) alignment pinned at exact k-mer anchors: lhs and rhs are aligned
// piece by piece between consecutive anchors (local (lhs, rhs) k-mer starts,
// increasing) and the pieces concatenated, so the cost grows with the pieces
// instead of the whole overlap. The path must also pass through both corners,
// so anchors are chained from (0, 0) and must reach (lhs.size(), rhs.size())
// with no more diagonal shift than indel drift explains; an anchor far off
// that path (e.g. from another repeat copy) would force a detour full NW does
// not take. Anchors that are out of order, out of range, not an exact k-mer
// match or off the path are skipped and counted in `rejected`. Without usable
// anchors this is edlib_wrapper. With max_edit_distance >= 0 it gives up
// (empty cigar) once the pieces' summed edit distance exceeds it: each piece
// gets the budget the earlier ones left.
inline edlib_align anchored_edlib_wrapper(
  const std::string &lhs,
  const std::string &rhs,
  const std::vector<std::pair<std::uint32_t, std::uint32_t>> &anchors,
  std::uint32_t k,
  std::int64_t max_edit_distance = -1,
  std::uint32_t *used = nullptr,
  std::uint32_t *rejected = nullptr) {

  edlib_align alignment_result;
  alignment_result.cigar = "";
  alignment_result.matches = 0;
  alignment_result.edit_distance = -1;
  alignment_result.block_length = 0;

  std::vector<unsigned char> ops;  // edlib alignment ops of all pieces
  ops.reserve(lhs.size() + lhs.size() / 8);
  std::int64_t edit_distance = 0;
  std::uint32_t lhs_prev = 0;
  std::uint32_t rhs_prev = 0;

  auto align_piece = [&](std::uint32_t lhs_to, std::uint32_t rhs_to) -> bool {
    const std::uint32_t lhs_len = lhs_to - lhs_prev;
    const std::uint32_t rhs_len = rhs_to - rhs_prev;
    if (lhs_len == 0 || rhs_len == 0) {  // pure gap
      ops.insert(ops.end(), lhs_len, EDLIB_EDOP_INSERT);
      ops.insert(ops.end(), rhs_len, EDLIB_EDOP_DELETE);
      edit_distance += lhs_len + rhs_len;
      if (max_edit_distance >= 0 && edit_distance > max_edit_distance) {
        return false;
      }
    } else {
      EdlibAlignResult result = capped_nw_align(
        lhs.data() + lhs_prev, lhs_len,
        rhs.data() + rhs_prev, rhs_len,
        max_edit_distance < 0 ? -1 : max_edit_distance - edit_distance);
      if (result.status != EDLIB_STATUS_OK || result.editDistance < 0) {
        edlibFreeAlignResult(result);
        return false;
      }
      ops.insert(ops.end(), result.alignment, result.alignment + result.alignmentLength);
      edit_distance += result.editDistance;
      edlibFreeAlignResult(result);
    }
    lhs_prev = lhs_to;
    rhs_prev = rhs_to;
    return true;
  };

  // largest diagonal shift plausible from indel drift between two path points
  auto on_path = [](std::pair<std::uint32_t, std::uint32_t> from,
                    std::pair<std::uint32_t, std::uint32_t> to) -> bool {
    const std::int64_t dl = static_cast<std::int64_t>(to.first) - from.first;
    const std::int64_t dr = static_cast<std::int64_t>(to.second) - from.second;
    return std::abs(dl - dr) <= 32 + std::max(dl, dr) / 8;
  };

  std::vector<std::pair<std::uint32_t, std::uint32_t>> path;
  path.reserve(anchors.size());
  std::pair<std::uint32_t, std::uint32_t> prev{0, 0};
  for (const auto &a : anchors) {
    if (a.first < prev.first || a.second < prev.second ||
        a.first + k > lhs.size() || a.second + k > rhs.size() ||
        !on_path(prev, a) ||
        lhs.compare(a.first, k, rhs, a.second, k) != 0) {
      continue;
    }
    path.push_back(a);
    prev = a;
  }
  const std::pair<std::uint32_t, std::uint32_t> end{lhs.size(), rhs.size()};
  while (!path.empty() && !on_path(path.back(), end)) {
    path.pop_back();
  }

  for (const auto &a : path) {
    if (!align_piece(a.first, a.second)) {
      return alignment_result;
    }
  }
  if (!align_piece(lhs.size(), rhs.size())) {
    return alignment_result;
  }
  if (used) *used = path.size();
  if (rejected) *rejected = anchors.size() - path.size();

  char *cigar = edlibAlignmentToCigar(ops.data(), ops.size(), EDLIB_CIGAR_STANDARD);
  alignment_result.cigar = cigar;
  free(cigar);
  alignment_result.edit_distance = static_cast<std::int32_t>(edit_distance);
  alignment_result.block_length = ops.size();
  alignment_result.matches = std::count(ops.begin(), ops.end(), EDLIB_EDOP_MATCH);
  return alignment_result;
}

// Infix (HW) alignment: align `query` fully as an infix of `target`; the
// target's unaligned prefix/suffix are free. Returns the cigar of query vs the
// matched infix, plus that infix's [target_start, target_end] (inclusive)
// within `target`, so the caller can recover the exact target coordinates
// regardless of imprecision in the seed region it padded.
inline edlib_align edlib_wrapper_infix(
  const std::string &query,
  const std::string &target) {

  edlib_align alignment_result;

  EdlibAlignResult result = edlibAlign(
    query.c_str(), query.size(),
    target.c_str(), target.size(),
    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0));

  if (result.status == EDLIB_STATUS_OK) {
    alignment_result.cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    alignment_result.edit_distance = result.editDistance;
    alignment_result.block_length = result.alignmentLength;
    alignment_result.matches = 0;
    for (int i = 0; i < result.alignmentLength; i++) {
      if (result.alignment[i] == 0) {
        ++alignment_result.matches;
      }
    }
    alignment_result.target_end =
        (result.numLocations > 0 && result.endLocations != nullptr)
            ? result.endLocations[0]
            : static_cast<std::int32_t>(target.size()) - 1;
    alignment_result.target_start =
        (result.numLocations > 0 && result.startLocations != nullptr)
            ? result.startLocations[0]
            : 0;
    edlibFreeAlignResult(result);
    return alignment_result;
  } else {
    edlibFreeAlignResult(result);
    alignment_result.cigar = "";
    alignment_result.matches = 0;
    alignment_result.edit_distance = -1;
    alignment_result.block_length = 0;
    return alignment_result;
  }
}

void find_pairwise_alignment(std::uint32_t i,
                             std::vector<extended_overlap> &ovlps,
                             std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                             raven::Graph &graph) {
  for (auto &it : ovlps) {
    if (!overlap_update(it.overlap, graph)) {
      continue;
    }
    std::string lhs = sequences[i]->InflateData(it.overlap.lhs_begin, it.overlap.lhs_end - it.overlap.lhs_begin);
    biosoup::NucleicAcid rhs_{ "",
                               sequences[it.overlap.rhs_id]->InflateData(it.overlap.rhs_begin,
                                                                         it.overlap.rhs_end - it.overlap.rhs_begin) };

    if (!it.overlap.strand) rhs_.ReverseAndComplement();

    auto rhs = rhs_.InflateData();

    it.edlib_alignment = edlib_wrapper(lhs, rhs);
  }
}

void call_snps(std::uint32_t i,
               std::vector<extended_overlap> ovlps_final,
               std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
               raven::Graph &graph) {
  std::uint32_t seq_inflated_len = sequences[i]->inflated_len;
  std::vector<base_pile> base_pile_tmp(seq_inflated_len);
  std::vector<std::uint32_t> cov;

  for (auto &ovlp : ovlps_final) {
    if (!(ovlp.edlib_alignment.cigar.empty())) {
      std::uint32_t lhs_pos = ovlp.overlap.lhs_begin;
      std::uint32_t rhs_pos = 0;
      biosoup::NucleicAcid rhs{ "",
                                sequences[ovlp.overlap.rhs_id]->InflateData(ovlp.overlap.rhs_begin,
                                                                            ovlp.overlap.rhs_end
                                                                              - ovlp.overlap.rhs_begin) };
      if (!ovlp.overlap.strand) rhs.ReverseAndComplement();
      std::string rhs_tmp = rhs.InflateData();

      std::string edlib_alignment = cigar_to_edlib_alignment(ovlp.edlib_alignment.cigar);
      for (auto &edlib_align : edlib_alignment) {
        switch (edlib_align) {
          case 0:
          case 3: {
            switch (rhs_tmp[rhs_pos]) {
              case 'A':++base_pile_tmp[lhs_pos].a;
                break;
              case 'C':++base_pile_tmp[lhs_pos].c;
                break;
              case 'G':++base_pile_tmp[lhs_pos].g;
                break;
              case 'T':++base_pile_tmp[lhs_pos].t;
                break;
              default:break; // if they align
            }
            ++lhs_pos;
            ++rhs_pos;
            break;
          }
          case 1: {
            ++base_pile_tmp[lhs_pos].i;
            ++lhs_pos;
            break; // insertion on the left hand side
          }
          case 2: {
            if (!(lhs_pos >= base_pile_tmp.size())) {
              ++base_pile_tmp[lhs_pos].d;
            }
            ++rhs_pos;
            break; // deletion on the left hand side
          }
          default:break;
        }
      }
    }
  }

  for (const auto &jt : base_pile_tmp) {
    //cov.emplace_back(jt.a + jt.c + jt.g + jt.t);
    cov.emplace_back(jt.a + jt.c + jt.g + jt.t + jt.d + jt.i);
  }

  std::nth_element(cov.begin(), cov.begin() + cov.size() / 2, cov.end());
  double m = cov[cov.size() / 2] * 2. / 3.;

  std::size_t j = 0;
  for (const auto &jt : base_pile_tmp) {
    std::vector<double> counts = {
      static_cast<double>(jt.a),
      static_cast<double>(jt.c),
      static_cast<double>(jt.g),
      static_cast<double>(jt.t),
      static_cast<double>(jt.d),
      static_cast<double>(jt.i)
    };

    double sum = std::accumulate(counts.begin(), counts.end(), 0);

    if (raven::use_frequencies) {
      for (auto &kt : counts) {
        kt /= sum;
      }
    }

    if (sum > m) {
      std::size_t variants = 0;
      for (const auto &it : counts) {
        if (raven::use_frequencies) {
          if (raven::freq_low_th < it && it < raven::freq_high_th) {
            ++variants;
          }
        } else {
          if (it > raven::variant_call_th) {
            ++variants;
          }
        }
      }
      if (variants > 1) graph.annotations_[i].emplace(j);
    }
    ++j;
  }
}

std::unordered_set<std::uint32_t> annotation_extract(
  std::uint32_t i,
  std::uint32_t begin,
  std::uint32_t end,
  std::uint32_t len,
  bool strand, raven::Graph &graph) {
  std::unordered_set<std::uint32_t> dst;
  if (graph.annotations_[i].empty()) {
    return dst;
  }
  for (const auto &it : graph.annotations_[i]) {
    if (begin <= it && it <= end) {
      dst.emplace(strand ? it : len - 1 - it);
    }
  }
  return dst;
}

#endif //RAVEN_SRC_GRAPH_CONSTRUCTOR_ANNOTATION_HELPERS_H_
