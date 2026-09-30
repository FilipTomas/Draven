
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

edlib_align edlib_wrapper(
  const std::string &lhs,
  const std::string &rhs) {

  edlib_align alignment_result;

  EdlibAlignResult result = edlibAlign(
    lhs.c_str(), lhs.size(),
    rhs.c_str(), rhs.size(),
    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0)); // align lhs and rhs

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
