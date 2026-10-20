#ifndef RAVEN_CORRECTION_CANDIDATE_H
#define RAVEN_CORRECTION_CANDIDATE_H

#include <cstdint>
#include <utility>
#include <vector>

#include "extended_overlap.h"

// A candidate correction overlap between mapping and RF selection: the
// fields containment, the RF features (overlap_selection) and the feature
// dump read, without the rest of extended_overlap (152 B instead of 376 B).
// lhs is the read of the pile it is stored in. Coordinates are those of the
// extended overlap; graph_overlap_type is set by containment for strong
// overlaps and is 0 for the others.
struct correction_candidate {
  std::uint32_t lhs_id = 0;
  std::uint32_t lhs_begin = 0;
  std::uint32_t lhs_end = 0;
  std::uint32_t rhs_id = 0;
  std::uint32_t rhs_begin = 0;
  std::uint32_t rhs_end = 0;
  std::uint32_t score = 0;
  std::uint32_t lhs_begin_original = 0;
  std::uint32_t lhs_end_original = 0;
  std::uint32_t rhs_begin_original = 0;
  std::uint32_t rhs_end_original = 0;
  std::uint32_t extended_length = 0;
  std::uint32_t found_length = 0;
  std::uint32_t longer_overhang_len = 0;
  std::uint32_t shorter_overhang_len = 0;
  float score_to_length = 0;
  float found_to_extended_length = 0;
  float hap_ratio = 0;
  float jaccard_index = 0;
  float jaccard_non_extended_overlap = 0;
  float diploid_jaccard_full_overlap = 0;
  float diploid_jaccard_longer = 0;
  float diploid_jaccard_shorter = 0;
  float lhs_hap = 0;
  float rhs_hap = 0;
  float lhs_err = 0;
  float rhs_err = 0;
  float lhs_rep = 0;
  float rhs_rep = 0;
  std::int8_t classification_label = 0;
  std::uint8_t graph_overlap_type = 0;
  bool strand = true;
  bool q_hor = false;
  bool t_hor = false;
  // chain anchors (lhs, rhs), see extended_overlap::anchors
  std::vector<std::pair<std::uint32_t, std::uint32_t>> anchors;
};

// The candidate fields of a mapped overlap (its anchors are moved).
// graph_overlap_type is left 0 for containment to set.
inline correction_candidate ToCorrectionCandidate(extended_overlap&& eo) {
  correction_candidate c;
  c.lhs_id = eo.overlap.lhs_id;
  c.lhs_begin = eo.overlap.lhs_begin;
  c.lhs_end = eo.overlap.lhs_end;
  c.rhs_id = eo.overlap.rhs_id;
  c.rhs_begin = eo.overlap.rhs_begin;
  c.rhs_end = eo.overlap.rhs_end;
  c.score = eo.overlap.score;
  c.strand = eo.overlap.strand;
  c.lhs_begin_original = eo.lhs_begin_original;
  c.lhs_end_original = eo.lhs_end_original;
  c.rhs_begin_original = eo.rhs_begin_original;
  c.rhs_end_original = eo.rhs_end_original;
  c.extended_length = eo.extended_length;
  c.found_length = eo.found_length;
  c.longer_overhang_len = eo.longer_overhang_len;
  c.shorter_overhang_len = eo.shorter_overhang_len;
  c.score_to_length = eo.score_to_length;
  c.found_to_extended_length = eo.found_to_extended_length;
  c.hap_ratio = eo.hap_ratio;
  c.jaccard_index = eo.jaccard_index;
  c.jaccard_non_extended_overlap = eo.jaccard_non_extended_overlap;
  c.diploid_jaccard_full_overlap = eo.diploid_jaccard_full_overlap;
  c.diploid_jaccard_longer = eo.diploid_jaccard_longer;
  c.diploid_jaccard_shorter = eo.diploid_jaccard_shorter;
  c.lhs_hap = eo.lhs_hap;
  c.rhs_hap = eo.rhs_hap;
  c.lhs_err = eo.lhs_err;
  c.rhs_err = eo.rhs_err;
  c.lhs_rep = eo.lhs_rep;
  c.rhs_rep = eo.rhs_rep;
  c.classification_label = static_cast<std::int8_t>(eo.classification_label);
  c.q_hor = eo.q_hor;
  c.t_hor = eo.t_hor;
  c.anchors = std::move(eo.anchors);
  return c;
}

// An extended_overlap with a candidate's fields (others 0), as the RF
// features and the feature dump read them. Anchors are not copied.
inline extended_overlap ToExtendedOverlap(const correction_candidate& c) {
  extended_overlap eo{};
  eo.overlap.lhs_id = c.lhs_id;
  eo.overlap.lhs_begin = c.lhs_begin;
  eo.overlap.lhs_end = c.lhs_end;
  eo.overlap.rhs_id = c.rhs_id;
  eo.overlap.rhs_begin = c.rhs_begin;
  eo.overlap.rhs_end = c.rhs_end;
  eo.overlap.score = c.score;
  eo.overlap.strand = c.strand;
  eo.lhs_begin_original = c.lhs_begin_original;
  eo.lhs_end_original = c.lhs_end_original;
  eo.rhs_begin_original = c.rhs_begin_original;
  eo.rhs_end_original = c.rhs_end_original;
  eo.extended_length = c.extended_length;
  eo.found_length = c.found_length;
  eo.longer_overhang_len = c.longer_overhang_len;
  eo.shorter_overhang_len = c.shorter_overhang_len;
  eo.score_to_length = c.score_to_length;
  eo.found_to_extended_length = c.found_to_extended_length;
  eo.hap_ratio = c.hap_ratio;
  eo.jaccard_index = c.jaccard_index;
  eo.jaccard_non_extended_overlap = c.jaccard_non_extended_overlap;
  eo.diploid_jaccard_full_overlap = c.diploid_jaccard_full_overlap;
  eo.diploid_jaccard_longer = c.diploid_jaccard_longer;
  eo.diploid_jaccard_shorter = c.diploid_jaccard_shorter;
  eo.lhs_hap = c.lhs_hap;
  eo.rhs_hap = c.rhs_hap;
  eo.lhs_err = c.lhs_err;
  eo.rhs_err = c.rhs_err;
  eo.lhs_rep = c.lhs_rep;
  eo.rhs_rep = c.rhs_rep;
  eo.classification_label = c.classification_label;
  eo.graph_overlap_type = c.graph_overlap_type;
  eo.q_hor = c.q_hor;
  eo.t_hor = c.t_hor;
  return eo;
}

#endif  // RAVEN_CORRECTION_CANDIDATE_H
