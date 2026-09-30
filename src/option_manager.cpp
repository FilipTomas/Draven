
#include <getopt.h>
#include <iostream>
#include <filesystem>
#include "option_manager.h"

enum program_opt {
  opt_error_corrected_reads,
  opt_herro_snps,
  opt_split,
  opt_graphical_fragment_assembly,
  opt_resume,
  opt_disable_checkpoints,
  opt_version,
  opt_help,
  opt_output,
  opt_input_gfa,
  opt_input_paf,
  opt_paf,
  opt_load_paf,
  opt_load_gt_overlaps,
  opt_print_gfa_seq,
  opt_correction_overlaps,

  opt_kmer_len,
  opt_window_len,
  opt_frequency,
  opt_hpc,
  opt_bandwidth,
  opt_chain_n,
  opt_match_n,
  opt_gap_size,
  opt_fraction,
  opt_coverage,
  opt_minimizers,
  opt_no_fuse,
  opt_internal_count,
  opt_bloom_gb,
  opt_store_kmer_ids,
  opt_merge_fragments,
  opt_infix_align,
  opt_dump_overlaps_only,
  opt_dump_all_overlaps,
  opt_rf_model,
  opt_rf_window_size,
  opt_rf_top_k,
  opt_rf_any_overlap,
  opt_dump_rf_scores,
  opt_fastk_counts,
  opt_threads,

  opt_max_overlaps,
  opt_disagreement,
  opt_valid_region_length,
  opt_valid_region_coverage,
  opt_cleaning_rounds,
  opt_min_unitig_size,
  opt_ploidy,

  opt_ultralong_phasing,
  opt_polishing_rounds
};

static struct option options[] = {
    {"error-corrected-reads", required_argument, nullptr, opt_error_corrected_reads},
    {"herro-snps", required_argument, nullptr, opt_herro_snps},
    {"split", required_argument, nullptr, opt_split},
    {"graphical-fragment-assembly", required_argument, nullptr, opt_graphical_fragment_assembly},
    {"resume", no_argument, nullptr, opt_resume},
    {"disable-checkpoints", no_argument, nullptr, opt_disable_checkpoints},
    {"version", no_argument, nullptr, opt_version},
    {"help", no_argument, nullptr, opt_help},
    {"output", required_argument, nullptr, opt_output},
    {"input-gfa", required_argument, nullptr, opt_input_gfa},
    {"input-paf", required_argument, nullptr, opt_input_paf},
    {"paf", no_argument, nullptr, opt_paf},
    {"load-paf", required_argument, nullptr, opt_load_paf},
    {"load-gt-overlaps", required_argument, nullptr, opt_load_gt_overlaps},
    {"print-seq", no_argument, nullptr, opt_print_gfa_seq},
    {"correction-overlaps", no_argument, nullptr, opt_correction_overlaps},

    {"kmer-len", required_argument, nullptr, opt_kmer_len},
    {"window-len", required_argument, nullptr, opt_window_len},
    {"frequency", required_argument, nullptr, opt_frequency},
    {"hpc", no_argument, nullptr, opt_hpc},
    {"bandwidth", required_argument, nullptr, opt_bandwidth},
    {"chain_n", required_argument, nullptr, opt_chain_n},
    {"match_n", required_argument, nullptr, opt_match_n},
    {"gap_size", required_argument, nullptr, opt_gap_size},
    {"fraction", required_argument, nullptr, opt_fraction},
    {"coverage", required_argument, nullptr, opt_coverage},
    {"minimizers", no_argument, nullptr, opt_minimizers},
    {"no-fuse", no_argument, nullptr, opt_no_fuse},
    {"internal-count", no_argument, nullptr, opt_internal_count},
    {"bloom-gb", required_argument, nullptr, opt_bloom_gb},
    {"store-kmer-ids", no_argument, nullptr, opt_store_kmer_ids},
    {"merge-fragments", no_argument, nullptr, opt_merge_fragments},
    {"infix-align", no_argument, nullptr, opt_infix_align},
    {"dump-overlaps-only", no_argument, nullptr, opt_dump_overlaps_only},
    {"dump-all-overlaps", no_argument, nullptr, opt_dump_all_overlaps},
    {"rf-model", required_argument, nullptr, opt_rf_model},
    {"rf-window-size", required_argument, nullptr, opt_rf_window_size},
    {"rf-top-k", required_argument, nullptr, opt_rf_top_k},
    {"rf-any-overlap", no_argument, nullptr, opt_rf_any_overlap},
    {"dump-rf-scores", no_argument, nullptr, opt_dump_rf_scores},
    {"fastk-counts", required_argument, nullptr, opt_fastk_counts},
    {"threads", required_argument, nullptr, opt_threads},

    {"max-overlaps", required_argument, nullptr, opt_max_overlaps},
    {"disagreement", required_argument, nullptr, opt_disagreement},
    {"valid-region-length", required_argument, nullptr, opt_valid_region_length},
    {"valid-region-coverage", required_argument, nullptr, opt_valid_region_coverage},
    {"cleaning-rounds", required_argument, nullptr, opt_cleaning_rounds},
    {"min-unitig-size", required_argument, nullptr, opt_min_unitig_size},
    {"ploidy", required_argument, nullptr, opt_ploidy},

    {"ultralong-phasing", required_argument, nullptr, opt_ultralong_phasing},
    {"polishing-rounds", required_argument, nullptr, opt_polishing_rounds},

    {nullptr, 0, nullptr, 0}
};

std::string optstr = "K:W:F:HB:C:M:G:t:X:D:R:";

void Help() {
  std::cout <<
            "Usage: raven [options ...] <sequences>\n"
            "\n"
            "  # default output is to stdout in FASTA format\n"
            "  <sequences>\n"
            "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
            "\n"
            "Options:\n"
            "  Input/Output:\n"
            "    --error-corrected-reads <string>\n"
            "      default: \"\"\n"
            "      path to error corrected reads\n"
            "    --herro-snps <string>\n"
            "      default: \"\"\n"
            "      path to herro snps\n"
            "    --split <int>\n"
            "      default: 0\n"
            "      graph coloring\n"
            "    --graphical-fragment-assembly <string>\n"
            "      prints the assembly graph in GFA format\n"
            "    --resume\n"
            "      resume previous run from last checkpoint\n"
            "    --disable-checkpoints\n"
            "      disable checkpoint file creation\n"
            "    --version\n"
            "      prints the version number\n"
            "    --help\n"
            "      prints this help text\n"
            "    --output <string>\n"
            "      output file name, if it is not set output is written to stdout\n"
            "      for diploid assembly, outputs will be written in 2 files with suffixes -1, -2\n"
            "    --input-gfa <string>\n"
            "      input GFA file name, if it is set raven will skip construction phase\n"
            "    --input-paf <string>\n"
            "      input PAF file name, if it is set raven will load overlaps from PAF file\n"
            "    --paf\n"
            "      overlaps are stored to files in PAF format\n"
            "    --load-paf <string>\n"
            "      load overlaps from PAF files\n"
            "    --print-seq\n"
            "      print sequences in GFA file\n"
            "\n"
            "  Overlap:\n"
            "    -K, --kmer-len <int>\n"
            "      default: 51\n"
            "      length of minimizers used to find overlaps\n"
            "    -W, --window-len <int>\n"
            "      default: 30\n"
            "      length of sliding window from which minimizers are sampled\n"
            "    -F, --frequency <double>\n"
            "      default: 0.001\n"
            "      threshold for ignoring most frequent minimizers\n"
            "    -H, --hpc\n"
            "     default: false\n"
            "     use HPC overlaps\n"
            "    -B, --bandwidth <int>\n"
            "      default: 500\n"
            "      bandwidth\n"
            "    -C, --chain_n <int>\n"
            "      default: 4\n"
            "      chain_n\n"
            "    -M, --match_n <int>\n"
            "      default: 100\n"
            "      match_n\n"
            "    -G, --gap_size <int>\n"
            "      default: 10000\n"
            "      gap_size\n"
            "    --fraction <double>\n"
            "      default: 0.2\n"
            "      fraction of sequences to use for kmer counting\n"
            "    --coverage <int>\n"
            "      default: 25\n"
            "      haploid sequencing coverage\n"
            "    --no-fuse\n"
            "      default: false\n"
            "      use the un-fused 3-pass FastK path (sketch, index, map done\n"
            "      separately) instead of the fused single-pass path; FastK counts\n"
            "      are still loaded (for apples-to-apples benchmarking)\n"
            "    --internal-count\n"
            "      default: false\n"
            "      count k-mers in-process instead of loading a precomputed FastK\n"
            "      table (singletons are dropped); k-mer counts come from the\n"
            "      internal counter\n"
            "    --bloom-gb <double>\n"
            "      default: -1 (auto: size by k-mer instances)\n"
            "      internal-counter bloom size in GB; 0 disables the bloom (exact\n"
            "      count + prune, no false positives, uses more memory)\n"
            "    --store-kmer-ids\n"
            "      default: false\n"
            "      store per-position k-mer ids on each pile (faster Jaccard, more\n"
            "      memory) instead of recomputing them on demand\n"
            "    --merge-fragments\n"
            "      default: false\n"
            "      merge collinear overlap fragments per read pair (span the\n"
            "      outermost anchors) instead of keeping only the longest fragment\n"
            "      and extrapolating it -- more accurate coordinates for error-rich\n"
            "      reads\n"
            "    --infix-align\n"
            "      default: false\n"
            "      align correction overlaps with edlib infix (HW) on a padded\n"
            "      target window and recover exact target coordinates from the\n"
            "      alignment, instead of global (NW) on fixed seed coordinates\n"
            "    --dump-overlaps-only\n"
            "      default: false\n"
            "      dump pre-alignment backbone correction overlaps with features\n"
            "      and stop, skipping the alignment step (debugging)\n"
            "    --dump-all-overlaps\n"
            "      default: false\n"
            "      dump EVERY raw ram overlap fragment (no dedup/merge/features/\n"
            "      alignment) to all_overlaps.paf and stop -- for comparing raw\n"
            "      ram coordinates against minimap2 (debugging)\n"
            "    --rf-model <string>\n"
            "      random forest exported with scripts/export_rf_model.py; compute,\n"
            "      score and select correction overlaps in-process (top-K per\n"
            "      target window) instead of loading a selected PAF\n"
            "    --rf-window-size <int>\n"
            "      default: 4096\n"
            "      target window size for --rf-model selection\n"
            "    --rf-top-k <int>\n"
            "      default: 30\n"
            "      overlaps kept per target window for --rf-model selection\n"
            "    --rf-any-overlap\n"
            "      default: false\n"
            "      any overlap touching a window is eligible, instead of only\n"
            "      those spanning the whole window\n"
            "    --dump-rf-scores\n"
            "      default: false\n"
            "      write every candidate's probability and selection to\n"
            "      rf_selection.tsv\n"
            "    --fastk-counts <string>\n"
            "      path to the precomputed FastK counts (used when --internal-count\n"
            "      is off); build FastK on the same reads for an apples-to-apples\n"
            "      comparison against --internal-count\n"
            "    -t, --threads <int>\n"
            "      default: 1\n"
            "      number of threads\n"
            "\n"
            "  Layout:\n"
            "    -X, --max-overlaps <long unsigned int>\n"
            "      default: 32\n"
            "      maximum number of overlaps that will be taken during find overlaps and create piles stage\n"
            "    -D, --disagreement <double>\n"
            "      default: 0.1\n"
            "      maximal percentage of different anntoated bases in overlaps\n"
            "    -R --valid-region <int>\n"
            "      default: 4\n"
            "      valid region minimal length\n"
            "    --cleaning-rounds <int>\n"
            "      default: 0\n"
            "      number of cleaning rounds\n"
            "    --min-unitig-size <int>\n"
            "      default: 9999\n"
            "      minimal unitig length\n"
            "    --ploidy <int>\n"
            "      default: 2\n"
            "      number of haplotypes\n"
            "\n"
            "  Phasing:\n"
            "    --ultralong-phasing <string>\n"
            "       path to ul reads used for phasing\n"
            "\n"
            "  Consensus:\n"
            "    --polishing-rounds <int>\n"
            "      default: 0\n"
            "      number of times racon is invoked\n";

}

int ProcessParameters(int argc, char **argv, Program_Parameters& param) {
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case opt_error_corrected_reads:
        param.error_corrected_reads = optarg;
        break;
      case opt_herro_snps:
        param.herro_snps_path = optarg;
        break;
      case opt_split:
        param.split = std::atoi(optarg);
        break;
      case opt_graphical_fragment_assembly:
        param.gfa_path = optarg;
        break;
      case opt_resume:
        param.resume = true;
        break;
      case opt_disable_checkpoints:
        param.checkpoints = false;
        break;
      case opt_version:
        std::cout << RAVEN_VERSION << std::endl;
        return 0;
      case opt_help:
        Help();
        return 0;
      case opt_output:
        param.root_path = optarg;
        {
          std::filesystem::path noext{""};
          std::filesystem::path path_prefix {param.root_path};
          path_prefix.replace_extension(noext);
          if (!path_prefix.has_filename())
            path_prefix /= "raven";
          param.root_path = path_prefix.string();

          param.gfa_path = param.root_path + "_" + param.gfa_path;
          param.gfa_after_overlap_graph_construction_filename = param.root_path + "_" + param.gfa_after_overlap_graph_construction_filename;
          param.gfa_after_construction_filename = param.root_path + "_" + param.gfa_after_construction_filename;
          param.gfa_after_transitive_filename = param.root_path + "_" + param.gfa_after_transitive_filename;
          param.gfa_after_bubble_filename = param.root_path + "_" + param.gfa_after_bubble_filename;
          param.gfa_after_force_filename = param.root_path + "_" + param.gfa_after_force_filename;
          param.gfa_post_construction_filename = param.root_path + "_" + param.gfa_post_construction_filename;
          param.gfa_post_cleaning_filename = param.root_path + "_" + param.gfa_post_cleaning_filename;
          param.gfa_unitig_graph_filename = param.root_path + "_" + param.gfa_unitig_graph_filename;

          param.paf_initial_overlaps_filename = param.root_path + "_" + param.paf_initial_overlaps_filename;
          param.paf_after_snp_filename = param.root_path + "_" + param.paf_after_snp_filename;
          param.paf_after_contained_filename = param.root_path + "_" + param.paf_after_contained_filename;
          param.paf_after_chimeric_filename = param.root_path + "_" + param.paf_after_chimeric_filename;
          param.paf_before_parsing_edges_filename = param.root_path + "_" + param.paf_before_parsing_edges_filename;

          param.cereal_filename = param.root_path + "_" + param.cereal_filename;
        }
        break;
      case opt_input_gfa:
        param.load_input_gfa = true;
        param.input_gfa_path = optarg;
        param.skip_contruction = true;
        param.skip_loading_fasta = true;
        break;
      case opt_input_paf:
        param.load_input_paf = true;
        param.input_paf_path = optarg;
        param.skip_contruction = true;
        break;
      case opt_paf:
        param.paf = true;
        break;
      case opt_load_paf:
        param.load_paf = optarg;
        break;
      case opt_load_gt_overlaps:
        param.gt_overlaps = optarg;
        break;
      case opt_print_gfa_seq:
        param.print_gfa_seq = true;
        break;
      case opt_kmer_len:
      case 'K':
        param.kmer_len = std::atoi(optarg);
        break;
      case opt_window_len:
      case 'W':
        param.window_len = std::atoi(optarg);
        break;
      case opt_frequency:
      case 'F':
        param.freq = std::atof(optarg);
        break;
      case opt_hpc:
      case 'H':
        param.hpc = true;
        break;
      case opt_bandwidth:
      case 'B':
        param.bandwidth = std::atoi(optarg);
        break;
      case opt_chain_n:
      case 'C':
        param.chain_n = std::atoi(optarg);
        break;
      case opt_match_n:
      case 'M':
        param.match_n = std::atoi(optarg);
        break;
      case opt_gap_size:
      case 'G':
        param.gap_size = std::atoi(optarg);
        break;
      case opt_fraction:
        param.freq = std::atof(optarg);
        break;
      case opt_coverage:
        param.coverage = std::atoi(optarg);
        break;
      case opt_minimizers:
        param.minimizers = true;
        break;
      case opt_no_fuse:
        param.no_fuse = true;
        break;
      case opt_internal_count:
        param.internal_count = true;
        break;
      case opt_bloom_gb:
        param.bloom_gb = std::atof(optarg);
        break;
      case opt_store_kmer_ids:
        param.store_kmer_ids = true;
        break;
      case opt_merge_fragments:
        param.merge_fragments = true;
        break;
      case opt_infix_align:
        param.infix_align = true;
        break;
      case opt_dump_overlaps_only:
        param.dump_overlaps_only = true;
        break;
      case opt_dump_all_overlaps:
        param.dump_all_overlaps = true;
        break;
      case opt_rf_model:
        param.rf_model_path = optarg;
        break;
      case opt_rf_window_size:
        param.rf_window_size = std::atoi(optarg);
        break;
      case opt_rf_top_k:
        param.rf_top_k = std::atoi(optarg);
        break;
      case opt_rf_any_overlap:
        param.rf_require_full_window = false;
        break;
      case opt_dump_rf_scores:
        param.dump_rf_scores = true;
        break;
      case opt_fastk_counts:
        param.fastk_path = optarg;
        break;
      case opt_threads:
      case 't':
        param.num_threads = atoi(optarg);
        break;
      case opt_max_overlaps:
      case 'X':
        param.max_overlaps = std::atof(optarg);
        break;
      case opt_disagreement:
      case 'D':
        param.disagreement = std::atof(optarg);
        break;
      case opt_valid_region_length:
      case 'R':
        param.valid_region_length_threshold = std::atoi(optarg);
        break;
      case opt_valid_region_coverage:
        param.valid_region_coverage_threshold = std::atoi(optarg);
        break;
      case opt_cleaning_rounds:
        break;
      case opt_min_unitig_size:
        param.min_unitig_size = std::atoi(optarg);
        break;
      case opt_ploidy:
        param.ploidy = std::atoi(optarg);
        break;
      case opt_ultralong_phasing:
        param.ul_read_path = optarg;
        break;
      case opt_polishing_rounds:
        param.num_polishing_rounds = atoi(optarg);
        break;
      case opt_correction_overlaps:
        param.correction_overlaps = true;
        break;
      default:
        return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc && !param.skip_contruction) {
    std::cerr << "[raven::] error: missing input file!" << std::endl;
    return 0;
  }

  param.sequence_path = argv[optind];

  return 1;
}
