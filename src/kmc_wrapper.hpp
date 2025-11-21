#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

#include "biosoup/nucleic_acid.hpp"

// KMC API headers (ensure they are on your include path and you link kmc_api)
#include "kmc_file.h"   // CKMCFile
#include "kmer_api.h"   // CKmerAPI

// Tiny helper around KMC counting + read-only API
class KmcCounter {
public:
  // Construct with K, resources and basic knobs.
  // kmc_bin: absolute path to 'kmc' executable (e.g. <conda-env>/bin/kmc).
  // workdir: working directory for temporary files + DB when building.
  KmcCounter(std::string kmc_bin,
             std::string workdir,
             int k,
             int threads   = 4,
             int mem_gb    = 4,      // KMC requires >= 2
             bool canonical = true,  // -b
             uint32_t min_count = 1) // -ci
  : kmc_bin_(std::move(kmc_bin)),
    workdir_(std::move(workdir)),
    k_(k),
    threads_(threads),
    mem_gb_(mem_gb),
    canonical_(canonical),
    min_count_(min_count) {
    if (k_ <= 0) throw std::runtime_error("K must be > 0");
  }

  // Build DB from in-memory sequences (writes a temporary FASTA).
  // fraction in (0,1] selects only the first ceil(N*fraction) sequences.
  void build_from_sequences(const std::vector<biosoup::NucleicAcid>& sequences,
                            double fraction = 1.0) {
    namespace fs = std::filesystem;
    fs::create_directories(workdir_);
    tmp_fasta_ = (fs::path(workdir_) / "kmc_input.fa").string();
    db_prefix_ = (fs::path(workdir_) / "kmc_db").string();
    tmp_dir_   = (fs::path(workdir_) / "tmp").string();
    fs::create_directories(tmp_dir_);

    write_fasta_(sequences, fraction, tmp_fasta_);
    run_kmc_fasta_(tmp_fasta_, db_prefix_, tmp_dir_);
    open_db_();
  }

  // Build DB from existing files using KMC list-file (@file).
  // Pass FASTA/FASTQ paths in 'inputs'. Set is_fasta=false for FASTQ.
  void build_from_filelist(const std::vector<std::string>& inputs,
                           bool is_fasta = true) {
    namespace fs = std::filesystem;
    fs::create_directories(workdir_);
    db_prefix_ = (fs::path(workdir_) / "kmc_db").string();
    tmp_dir_   = (fs::path(workdir_) / "tmp").string();
    fs::create_directories(tmp_dir_);

    const std::string list_path = (fs::path(workdir_) / "inputs.lst").string();
    {
      std::ofstream out(list_path);
      if (!out) throw std::runtime_error("Cannot write list file: " + list_path);
      for (const auto& p : inputs) out << p << "\n";
    }
    run_kmc_list_(list_path, db_prefix_, tmp_dir_, is_fasta);
    open_db_();
  }

  // Open an already-built KMC DB (prefix path without .kmc_pre/.kmc_suf).
  void open_existing(const std::string& db_prefix) {
    close(); // close previous if open
    db_prefix_ = db_prefix;
    if (!db_.OpenForRA(db_prefix_)) {
      throw std::runtime_error("Failed to open KMC DB: " + db_prefix_);
    }
    if (db_.KmerLength() != k_) {
      throw std::runtime_error("K in DB (" + std::to_string(db_.KmerLength()) +
                               ") != requested K (" + std::to_string(k_) + ")");
    }
    opened_ = true;
    // Optionally grab params if needed:
    // uint32 min_cnt, uint32 max_cnt; db_.GetMinCount(min_cnt); db_.GetMaxCount(max_cnt);
  }

  // Histogram of counts: count -> number of distinct k-mers at that count.
  // Iterates once over all distinct k-mers (fast from SSD).
  std::unordered_map<uint32_t, uint64_t> histogram() {
    ensure_open_();
    std::unordered_map<uint32_t, uint64_t> hist;
    CKmerAPI kmer(k_);
    uint32_t cnt = 0;
    db_.RestartListing();
    while (db_.ReadNextKmer(kmer, cnt)) {
      ++hist[cnt];
    }
    db_.RestartListing();
    return hist;
  }

  // Count for an arbitrary k-mer string (length must equal k).
  // If DB was built with -b (canonical), pass canonical strings or let KMC handle it internally.
  uint32_t count_kmer(const std::string& s) {
    ensure_open_();
    if ((int)s.size() != k_) throw std::runtime_error("count_kmer: wrong k-mer length");
    CKmerAPI q(k_);
    q.from_string(s.c_str());
    uint32_t cnt = 0;
    db_.CheckKmer(q, cnt);
    return cnt;
  }

  // Accessors
  int k() const { return k_; }
  const std::string& db_prefix() const { return db_prefix_; }
  const std::string& workdir() const { return workdir_; }

  // Close DB explicitly (also called in destructor).
  void close() {
    if (opened_) {
      try { db_.Close(); } catch (...) {}
      opened_ = false;
    }
  }

  ~KmcCounter() { close(); }

private:
  // config
  std::string kmc_bin_;
  std::string workdir_;
  int k_;
  int threads_;
  int mem_gb_;
  bool canonical_;
  uint32_t min_count_;

  // state
  std::string tmp_fasta_;
  std::string db_prefix_;
  std::string tmp_dir_;
  CKMCFile db_;
  bool opened_ = false;

  void ensure_open_() const {
    if (!opened_) throw std::runtime_error("KMC DB is not open");
  }

  // ---------- building helpers ----------
  void write_fasta_(const std::vector<biosoup::NucleicAcid>& seqs,
                    double fraction,
                    const std::string& path) {
    if (fraction <= 0.0 || fraction > 1.0)
      throw std::runtime_error("fraction must be in (0,1]");
    size_t n = static_cast<size_t>(seqs.size() * fraction);
    if (n == 0 && !seqs.empty()) n = 1;

    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write FASTA: " + path);

    for (size_t i = 0; i < n; ++i) {
      out << '>' << seqs[i].name << '\n';
      out << seqs[i].InflateData() << '\n';  // biosoup returns std::string
    }
  }

  void run_kmc_fasta_(const std::string& fasta,
                      const std::string& db_prefix,
                      const std::string& tmp_dir) {
    std::ostringstream cmd;
    cmd << shell_escape_(kmc_bin_)
        << " -k"  << k_
        << " -t"  << threads_
        << " -m"  << mem_gb_         // must be >= 2
        << " -ci" << min_count_
        << (canonical_ ? " -b" : "")
        << " -fa " << shell_escape_(fasta)
        << " "     << shell_escape_(db_prefix)
        << " "     << shell_escape_(tmp_dir);
    run_or_throw_(cmd.str());
  }

  void run_kmc_list_(const std::string& list_file,
                     const std::string& db_prefix,
                     const std::string& tmp_dir,
                     bool is_fasta) {
    std::ostringstream cmd;
    cmd << shell_escape_(kmc_bin_)
        << " -k"  << k_
        << " -t"  << threads_
        << " -m"  << mem_gb_
        << " -ci" << min_count_
        << (canonical_ ? " -b" : "")
        << (is_fasta ? " -fa " : " -fq ")
        << "@" << shell_escape_(list_file)
        << " "  << shell_escape_(db_prefix)
        << " "  << shell_escape_(tmp_dir);
    run_or_throw_(cmd.str());
  }

  void open_db_() {
    if (!db_.OpenForRA(db_prefix_)) {
      throw std::runtime_error("Failed to open KMC DB: " + db_prefix_);
    }
    if (db_.KmerLength() != k_) {
      throw std::runtime_error("K in DB (" + std::to_string(db_.KmerLength()) +
                               ") != requested K (" + std::to_string(k_) + ")");
    }
    opened_ = true;
  }

  // ---------- shell helpers ----------
  static std::string shell_escape_(const std::string& s) {
    // simple single-quote escape for POSIX shells
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('\'');
    for (char c : s) {
      if (c == '\'') { out += "'\\''"; } else { out.push_back(c); }
    }
    out.push_back('\'');
    return out;
  }

  static void run_or_throw_(const std::string& cmd) {
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
      throw std::runtime_error("KMC command failed (" + std::to_string(ret) + "): " + cmd);
    }
  }
};
