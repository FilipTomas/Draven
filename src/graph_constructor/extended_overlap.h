#ifndef RAVEN_EXTENDED_OVERLAP_H
#define RAVEN_EXTENDED_OVERLAP_H

#include <vector>
#include "biosoup/overlap.hpp"


// Backbone - strong prefix-suffx, StrongContain - strong contain, WeakBackbone - weak between 2 piles marked as backbones(those with at 
// least 1 strong prefix-suffx overlap), WeakConnecting - one backbone and one none backbone, Weak - weak between 2 non-backbone piles, None
enum class overlapCategory : std::uint8_t { Backbone = 0, StrongContained = 1, WeakBackbone = 2, WeakConnecting = 3, Weak = 4,  None = 5};

struct edlib_align {
  std::uint32_t matches;
  std::uint32_t block_length;
  std::string cigar;
  std::int32_t edit_distance;
  // For HW/infix alignment: the [start, end] positions (inclusive) of the
  // aligned infix within the target that was passed to the aligner. Unused
  // (0) for global (NW) alignment.
  std::int32_t target_start = 0;
  std::int32_t target_end = 0;
};

struct extended_overlap {
  biosoup::Overlap overlap;
  edlib_align edlib_alignment;
  float identity;
  float heterozygosity_rate;
  std::uint32_t graph_overlap_type;
  std::uint32_t total_overlap_snps;
  std::uint32_t total_overlap_snp_mismatches;
  std::uint8_t ol_class;
  bool ground_truth;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> gap_positions;

  std::uint32_t extended_length;
  std::uint32_t found_length;
  std::uint32_t found_matches;
  std::uint32_t lhs_begin_original;
  std::uint32_t lhs_end_original;
  std::uint32_t rhs_begin_original;
  std::uint32_t rhs_end_original;

  std::uint32_t correction_lhs_begin;
  std::uint32_t correction_lhs_end;
  std::uint32_t correction_rhs_begin;
  std::uint32_t correction_rhs_end;

  float lhs_hap;
  float rhs_hap;

  float lhs_err;
  float rhs_err;

  float lhs_rep;
  float rhs_rep;

  
  float score_to_length;
  float found_to_extended_length;
  float hap_ratio;

  bool q_hor;
  bool t_hor;

  float jaccard_index;
  bool hap_regions;

  float diploid_jaccard_longer;
  float diploid_jaccard_shorter;

  float error_jaccard_longer;
  float error_jaccard_shorter;
  float diploid_jaccard_full_overlap;
  float jaccard_non_extended_overlap;

  std::uint32_t longer_overhang_len;
  std::uint32_t shorter_overhang_len;

  float hap_rate_shorter_overhang_lhs;
  float hap_rate_longer_overhang_lhs;

  float hap_rate_shorter_overhang_rhs;
  float hap_rate_longer_overhang_rhs;

  float dip_rate_shorter_overhang_lhs;
  float dip_rate_longer_overhang_lhs;

  float dip_rate_shorter_overhang_rhs;
  float dip_rate_longer_overhang_rhs;

  float err_rate_shorter_overhang_lhs;
  float err_rate_longer_overhang_lhs;

  float err_rate_shorter_overhang_rhs;
  float err_rate_longer_overhang_rhs;

  bool backbone_overlap;
  overlapCategory cat;
  bool candidate_overlap;

  bool correction_overlap;

  int classification_label;
};




#endif // RAVEN_EXTENDED_OVERLAP_H


