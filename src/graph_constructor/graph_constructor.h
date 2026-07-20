
#ifndef RAVEN_GRAPH_CONSTRUCTOR_H
#define RAVEN_GRAPH_CONSTRUCTOR_H

#include <memory>
#include <string>
#include "biosoup/nucleic_acid.hpp"
#include "pile.hpp"
#include "option_manager.h"
#include "biosoup/timer.hpp"
#include "ram/kmer_counter.hpp"
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
  
  void ConstructCorrectionOverlaps(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
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

  void LoadOverlapsAndAlign(Program_Parameters &param, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps);
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

  double disagreement_;
};

} // raven

#endif //RAVEN_GRAPH_CONSTRUCTOR_H
