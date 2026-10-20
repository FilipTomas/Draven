
#ifndef RAVEN_GRAPH_CONSTRUCTOR_H
#define RAVEN_GRAPH_CONSTRUCTOR_H

#include <functional>
#include <memory>
#include <string>
#include "biosoup/nucleic_acid.hpp"
#include "pile.hpp"
#include "option_manager.h"
#include "biosoup/timer.hpp"
#include "ram/kmer_counter.hpp"
#include "correction_candidate.h"
#include <numeric>

namespace raven {

class Graph_Constructor {
public:
  Graph_Constructor() = default;

  Graph_Constructor(Graph &graph,
                    std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

  // break chimeric sequences, remove contained sequences and overlaps not
  // spanning bridged repeats at sequence ends
  void Construct(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
    Program_Parameters &param);

  // Stream reads from `sequence_path` while counting their k-mers in-process,
  // so parsing (single-threaded I/O) overlaps counting (pool). Fills
  // `sequences` and populates the internal counter that Construct then reuses
  // (used on the --internal-count path instead of a separate Parse + count).
  void LoadAndCount(
    const std::string &sequence_path,
    Program_Parameters &param,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences);
  
  // Writes correction.paf (or the --dump-* debug output) and ends the process
  // with exit status 0 without running destructors.
  [[noreturn]] void ConstructCorrectionOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
    Program_Parameters &param);

  void LoadFromGfa(const std::string &gfa_path);
  void LoadFromPaf(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, const std::string &paf_path);
private:
  void ConstructOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                         std::vector<std::vector<extended_overlap>> &extended_overlaps,
                         biosoup::Timer &timer,
                         Program_Parameters &param);
  void ConstructOverlapsForCorrection(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                         std::vector<std::vector<extended_overlap>> &extended_overlaps,
                         biosoup::Timer &timer,
                         Program_Parameters &param);
                         
  void AlignBackboneReads(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                         std::vector<std::vector<extended_overlap>> &extended_overlaps,
                         biosoup::Timer &timer,
                         Program_Parameters &param);

  void ConstructOverlapsFromGT(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                               std::vector<std::vector<extended_overlap>> &extended_overlaps,
                               biosoup::Timer &timer,
                               Program_Parameters &param);
  void ConstructOverlapGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                              std::vector<std::vector<extended_overlap>> &overlaps,
                              biosoup::Timer &timer,
                              Program_Parameters &param);
                              
  void ConstructAssemblyGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                              std::vector<std::vector<extended_overlap>> &overlaps,
                              biosoup::Timer &timer,
                              Program_Parameters &param);

  void ConstructOnlyBackboneGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::vector<std::vector<extended_overlap>>& overlaps,
    biosoup::Timer& timer,
    Program_Parameters& param);
                        
  void ConstructBackboneGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    std::vector<std::vector<extended_overlap>>& overlaps,
    biosoup::Timer& timer,
    Program_Parameters& param);
  
  void ConstructAssemblyGraphInPhases(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                               std::vector<std::vector<extended_overlap>> &overlaps,
                                               biosoup::Timer &timer,
                                               Program_Parameters &param);

  void LoadOverlapsFromPaf(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                           std::vector<std::vector<extended_overlap>> &extended_overlaps,
                           bool load_cigar,
                           Program_Parameters &param);

  // pile_aligned(i), if set, runs on the worker thread as soon as pile i
  // is aligned (e.g. to write it out and free it).
  void LoadOverlapsAndAlign(Program_Parameters &param, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps,
                            const std::function<void(std::uint32_t)> &pile_aligned = nullptr);
  void AlignCorrectionOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                               std::vector<std::vector<extended_overlap>> &extended_overlaps,
                               std::uint32_t kmer_len,
                               double max_edit_rate,
                               const std::function<void(std::uint32_t)> &pile_aligned = nullptr);
  void SelectCorrectionOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                Program_Parameters &param);
  void MapSequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                    std::vector<std::vector<extended_overlap>> &extended_overlaps,
                    biosoup::Timer &timer,
                    Program_Parameters &param);

  void MapSequencesFast(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                        std::vector<std::vector<extended_overlap>> &extended_overlaps,
                        biosoup::Timer &timer,
                        Program_Parameters &param);
                        
  void TrimAndAnnotatePiles(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                            std::vector<std::vector<extended_overlap>> &extended_overlaps,
                            biosoup::Timer &timer,
                            Program_Parameters &param);
  void ResolveSnps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                   std::vector<std::vector<extended_overlap>> &extended_overlaps,
                   biosoup::Timer &timer,
                   Program_Parameters &param);
  void ResolveOverlapType(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                           std::vector<std::vector<extended_overlap>> &extended_overlaps,
                           biosoup::Timer &timer,
                           Program_Parameters &param);
void ResolveContainedReads(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
                            std::vector<std::vector<extended_overlap>>& extended_overlaps,
                            biosoup::Timer& timer,
                            std::uint32_t containment_threshold);
  // PASS 1 of ResolveContainedReads for one overlap copy (overlap type,
  // containment votes, backbone marks); returns the read that got a vote, or
  // -1. Candidate: extended_overlap or correction_candidate.
  template <typename Candidate>
  std::uint32_t ApplyContainmentPass1(
      Candidate& ext_ovlp,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
      std::uint32_t* containment_votes);
  // Correction mapping: apply PASS 1 to a copy as it is stored in
  // correction_candidates_, keeping nothing in the piles of reads that have
  // reached the containment threshold.
  void StoreCorrectionCandidate(
      correction_candidate&& candidate,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);
  // Correction path: what ResolveContainedReads does after PASS 1, which ran
  // as mapping stored the candidates (its PASS 2 categories go unused here).
  void ResolveContainedCandidates(biosoup::Timer& timer);
  void ResolveContainedReadsGT(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                              std::vector<std::vector<extended_overlap>> &extended_overlaps,
                              biosoup::Timer &timer);
  void ResolveChimericSequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                std::vector<std::vector<extended_overlap>> &extended_overlaps,
                                biosoup::Timer &timer);

  void PrintPiles(const std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences);

  void LoadHerroSNPs(const std::string &herro_snps_path,
                     std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences);
  void LoadAnnotations(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                       std::vector<std::vector<extended_overlap>> &extended_overlaps,
                       Program_Parameters &param);
  void LoadOverlaps(const std::string &overlaps_path,
                    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                    std::vector<std::vector<extended_overlap>> &extended_overlaps,
                    bool load_cigar);
  void LoadGTOverlaps(const std::string &overlaps_path,
                    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                    std::vector<std::vector<extended_overlap>> &extended_overlaps,
                    bool load_cigar);
  
  void FindHistogramPeaks();

  Graph &graph_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  // Populated by LoadAndCount (overlapped load+count); reused by MapSequencesFast.
  std::unique_ptr<ram::KmerCounter> internal_counter_;
  // Correction path: containment votes applied while mapping stores overlaps
  // (empty when not tracking), and the threshold for strong containment.
  std::vector<std::uint32_t> containment_votes_;
  std::uint32_t containment_threshold_ = 3;
  // Correction path: candidate overlaps per read, from mapping to selection
  std::vector<std::vector<correction_candidate>> correction_candidates_;

  double disagreement_;
};

} // raven

#endif //RAVEN_GRAPH_CONSTRUCTOR_H
