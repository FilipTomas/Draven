#ifndef RAVEN_EXTENDED_OVERLAP_H
#define RAVEN_EXTENDED_OVERLAP_H

#include <vector>
#include "biosoup/overlap.hpp"

struct edlib_align {
  std::uint32_t matches;
  std::uint32_t block_length;
  std::string cigar;
  std::int32_t edit_distance;
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

  int classification_label;
};




#endif // RAVEN_EXTENDED_OVERLAP_H


