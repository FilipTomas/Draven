
#ifndef RAVEN_OPTION_MANAGER_H
#define RAVEN_OPTION_MANAGER_H

#include <cstdint>
#include <string>

struct Program_Parameters {
  unsigned split = 0;

  std::uint8_t ploidy = 2;

  std::uint8_t kmer_len = 21;
  std::uint8_t window_len = 21;
  std::uint16_t bandwidth = 200;
  std::uint16_t chain_n = 4;
  std::uint16_t match_n = 200;
  std::uint16_t gap_size = 500;
  double fraction = 1; // fraction of sequences to use for kmer counting
  double freq = 1000;
  std::uint32_t coverage = 50;
  bool hpc = false;
  bool minimizers = false;
  bool no_fuse = false;  // force the un-fused 3-pass FastK path (benchmarking)
  bool internal_count = false;  // count k-mers in-process instead of loading FastK
  // internal-counter bloom size in GB. <0 = auto (size by k-mer instances);
  // 0 = bloom disabled (exact count + prune, no false positives); >0 = that
  // many GB of bloom bits.
  double bloom_gb = -1.0;
  // store per-position k-mer ids on each pile (faster Jaccard, more memory)
  // instead of recomputing them on demand (default: recompute)
  bool store_kmer_ids = false;
  // merge collinear overlap fragments per pair (outermost-anchor span) instead
  // of keeping only the longest fragment and extrapolating it (default: longest)
  bool merge_fragments = false;
  // align correction overlaps with edlib infix (HW) on a padded target window
  // and recover exact target coords, instead of global (NW) on fixed coords
  bool infix_align = false;
  // dump pre-alignment backbone correction overlaps (with features) and stop,
  // skipping the alignment step (debugging)
  bool dump_overlaps_only = false;
  // dump EVERY raw ram overlap fragment (no dedup, no merge, no features, no
  // alignment) to all_overlaps.paf and stop -- for comparing raw ram
  // coordinates against minimap2 (debugging)
  bool dump_all_overlaps = false;
  // random forest exported by scripts/export_rf_model.py; when set, correction
  // overlaps are computed, scored and selected in-process (top rf_top_k per
  // rf_window_size target window) instead of loading a selected PAF
  std::string rf_model_path = "";
  std::uint32_t rf_window_size = 4096;
  std::uint32_t rf_top_k = 30;
  // overlap must span the whole window to be eligible (else: touch it)
  bool rf_require_full_window = true;
  // write every candidate's probability and selection to rf_selection.tsv
  bool dump_rf_scores = false;
  // keep ram's chain anchors ~this many bp apart and align correction
  // overlaps piecewise between them (0: one global alignment per overlap)
  std::uint32_t anchor_spacing = 500;
  // give up aligning a correction overlap once its edit distance exceeds this
  // fraction of its length (0: no limit)
  double max_edit_rate = 0.3;
  // path to the precomputed FastK counts (used when --internal-count is off)
  std::string fastk_path = "/mnt/share1_Jabba/ftomas/fastk_counts/badread_ONT_k21/chr19";

  std::string error_corrected_reads = "";
  std::string herro_snps_path = "";
  std::string load_paf = "";

  std::int32_t num_polishing_rounds = 0;

  std::string ul_read_path;

  double disagreement = 0.01;
  std::string gfa_path {"final.gfa"};
  std::string input_gfa_path = "";
  bool skip_contruction = false;
  bool skip_loading_fasta = false;
  bool resume = false;
  bool checkpoints = true;

  std::string root_path{"out/raven"};

  bool load_input_gfa = false;
  bool load_input_paf = false;
  std::string input_paf_path = "";
  std::string sequence_path;
  std::string gt_overlaps;

  std::string gfa_after_overlap_graph_construction_filename {"after_overlap_graph_constructor.gfa"};
  std::string gfa_after_construction_filename {"after_construction.gfa"};
  std::string gfa_after_transitive_filename {"after_transitive.gfa"};
  std::string gfa_after_bubble_filename {"after_bubble.gfa"};
  std::string gfa_after_force_filename {"after_force.gfa"};
  std::string gfa_post_construction_filename {"post_construction.gfa"};
  std::string gfa_post_cleaning_filename {"post_cleaning.gfa"};
  std::string gfa_unitig_graph_filename {"unitig_graph.gfa"};

  std::string paf_initial_overlaps_filename {"initialOverlaps.paf"};
  std::string paf_after_snp_filename {"afterSNP.paf"};
  std::string paf_after_contained_filename {"afterContained.paf"};
  std::string paf_after_chimeric_filename {"afterChimeric.paf"};
  std::string paf_before_parsing_edges_filename {"beforeEdges.paf"};

  std::string cereal_filename {"backup.cereal"};

  std::uint32_t num_threads = 1;

  std::size_t max_overlaps = 16;
  std::uint32_t min_unitig_size = 9999;
  std::uint16_t valid_region_length_threshold = 1260;
  std::uint16_t valid_region_coverage_threshold = 4;
  bool paf = false;
  bool print_gfa_seq = false;
  bool correction_overlaps = false;
};

int ProcessParameters(int argc, char **argv, Program_Parameters& param);


#endif //RAVEN_OPTION_MANAGER_H
