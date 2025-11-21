    #include <cstdint>
    #include "extended_overlap.h"
    #include <vector>
    
    
    struct edlib_align {
      std::uint32_t matches;
      std::uint32_t block_legnth;
      std::string cigar;
      std::int32_t edit_distance;
    };

    struct extended_overlap {
      biosoup::Overlap overlap;
      edlib_align edlib_alignment;
      std::uint32_t total_overlap_snps;
      std::uint32_t total_overlap_snp_mismatches;
      float identity;
      float heterozygosity_rate;
      std::uint32_t graph_overlap_type;
      std::uint8_t ol_class;
      bool ground_truth = NULL;
      std::vector<std::pair<std::uint32_t, std::uint32_t>> gap_positions; // odd entries are start positions, even entries are end positions
      std::uint32_t found_length = 0;
      std::uint32_t found_matches = 0;
      std::uint32_t lhs_begin_original;
      std::uint32_t lhs_end_original;
      std::uint32_t rhs_begin_original;
      std::uint32_t rhs_end_original;

      float lhs_hap = 0.0;
      float rhs_hap = 0.0;

      float lhs_err = 0.0;
      float rhs_err = 0.0;

      float lhs_rep = 0.0;
      float rhs_rep = 0.0;

      float score_to_length = 0.0;
      float found_to_extended_length = 0.0;
      float hap_ratio = 0.0;

      bool q_hor = false;
      bool t_hor = false;

      int classification_label = -1;
    };