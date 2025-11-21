// Copyright (c) 2020 Robert Vaser

#ifndef RAVEN_GRAPH_HPP_
#define RAVEN_GRAPH_HPP_

#include <cstdint>
#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <utility>
#include <sstream>

#include <iostream>

#include "biosoup/nucleic_acid.hpp"
#include "cereal/access.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/unordered_set.hpp"
#include "cereal/types/vector.hpp"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

#include "pile.hpp"
#include "node.h"
#include "edge.h"
#include "marked_edge.h"
#include "graph_constructor/extended_overlap.h"
#include "state_manager.h"

namespace raven {
  constexpr std::uint8_t use_frequencies = 0;
  constexpr std::uint8_t variant_call_th = 3;
  constexpr double freq_low_th = 0.333;
  constexpr double freq_high_th = 0.667;
  constexpr std::uint8_t print_snp_data = 1;

  class Graph {
  public:
    Graph() = default;  // needed for cereal

    Graph(bool checkpoints, std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

    Graph(const Graph &) = delete;

    Graph &operator=(const Graph &) = default;

    Graph(Graph &&) = default;

    Graph &operator=(Graph &&) = default;

    ~Graph() = default;

    bool use_checkpoints() {
      return checkpoints_;
    }

    void DuplicateGraph();

    // ignore nodes that are less than epsilon away from any junction node
    std::uint32_t CreateUnitigs(std::uint32_t epsilon = 0);

    std::uint32_t CreateUnitigsAlternate(std::uint32_t epsilon = 0);

    void CreateUnitigGraph(const std::string &gfa_path);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetUnitigPairs(bool drop_unpolished = false);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetAssembledData(bool primary = true);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetUnitigs(bool drop_unpolished = false);

    // draw with misc/plotter.py
    void PrintJson(const std::string &path) const;

    // draw with Cytoscape
    void PrintCsv(const std::string &path) const;

    // draw with Bandage
    void PrintGfa(const std::string &path, const bool print_seq) const;

    // print overlaps in PAF file format
    void PrintOverlaps(std::vector<std::vector<extended_overlap>> overlaps,
                       std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, bool print_cigar,
                       const std::string &path) const;

    // draw unitig graph with Bandage
    void PrintUnitigGfa(const std::string &path, const bool print_seq) const;

    // cereal load wrapper
    void Load(const std::string &cereal_path);

    //
    //void LoadFromGfa(const std::string &nput_gfa_path);

    // cereal store wrapper
    void Store(const std::string &cereal_path) const;

    std::unordered_set<std::uint32_t> FindRemovableEdges(const std::vector<Node *> &path);

    void RemoveEdges(const std::unordered_set<std::uint32_t> &indices, bool remove_nodes = false);

    void RemoveUnitigEdges(const std::unordered_set<std::uint32_t> &indices, bool remove_nodes = false);

    void RemoveDiploidEdges(const std::vector<MarkedEdge> &indices, bool remove_nodes = false);

    void RemoveAlternateEdges(const std::vector<MarkedEdge> &indices, bool remove_nodes = false);

    // label small circular contigs as unitigs
    std::uint32_t SalvagePlasmids();

    void SalvageHaplotypes();

    void SalvageHaplotypesPrimary();

    void SalvageHaplotypesAlternative();

    std::vector<std::shared_ptr<Node>> nodes_;
    std::vector<std::shared_ptr<Edge>> edges_;

    std::vector<std::unordered_set<std::uint32_t>> annotations_;
    std::vector<std::unordered_set<std::uint32_t>> annotations_compressed_;
    std::vector<std::vector<std::uint32_t>> anno_;
    std::vector<std::unique_ptr<Pile>> piles_;
    //std::vector<std::vector<std::pair<std::unique_ptr<std::uint64_t>, std::uint16_t>>> minimizers_;

    std::vector<std::shared_ptr<Node>> unitig_nodes_;
    std::vector<std::shared_ptr<Edge>> unitig_edges_;

    std::vector<std::shared_ptr<Node>> nodes_alternate_;
    std::vector<std::shared_ptr<Edge>> edges_alternate_;

    StateManager state_manager_;

  private:
    friend cereal::access;

    template<class Archive>
    void save(Archive &archive) const {  // NOLINT
      std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;
      for (const auto &it: edges_) {
        if (it && !it->is_rc()) {
          connections.emplace_back(it->tail->id, it->head->id);
        } else {
          connections.emplace_back();  // dummy
        }
      }

      archive(state_manager_.state(), annotations_, piles_, nodes_, edges_, connections);
    }

    template<class Archive>
    void load(Archive &archive) {  // NOLINT
      std::vector<std::pair<std::uint32_t, std::uint32_t>> connections;

      GraphState loadedState;
      archive(loadedState, annotations_, piles_, nodes_, edges_, connections);

      state_manager_.set_state(loadedState);

      for (std::uint32_t i = 0; i < nodes_.size(); i += 2) {
        if (nodes_[i]) {
          nodes_[i]->pair = nodes_[i + 1].get();
          nodes_[i + 1]->pair = nodes_[i].get();
        }
      }
      for (std::uint32_t i = 0; i < edges_.size(); i += 2) {
        if (edges_[i]) {
          edges_[i]->pair = edges_[i + 1].get();
          edges_[i + 1]->pair = edges_[i].get();

          edges_[i]->tail = nodes_[connections[i].first].get();
          edges_[i]->head = nodes_[connections[i].second].get();
          edges_[i + 1]->head = edges_[i]->tail->pair;
          edges_[i + 1]->tail = edges_[i]->head->pair;

          edges_[i]->tail->outedges.emplace_back(edges_[i].get());
          edges_[i]->head->inedges.emplace_back(edges_[i].get());
          edges_[i + 1]->tail->outedges.emplace_back(edges_[i + 1].get());
          edges_[i + 1]->head->inedges.emplace_back(edges_[i + 1].get());
        }
      }

      Node::num_objects = nodes_.size();
      Edge::num_objects = edges_.size();
    }

    bool checkpoints_;
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  };

}  // namespace raven

namespace biosoup {

  template<class Archive>
  void serialize(Archive &archive, NucleicAcid &sequence) {  // NOLINT
    archive(
        sequence.id,
        sequence.name,
        sequence.deflated_data,
        sequence.block_quality,
        sequence.inflated_len,
        sequence.is_reverse_complement);
  }

}  // namespace biosoup

#endif  // RAVEN_GRAPH_HPP_
