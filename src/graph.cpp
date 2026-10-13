// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <limits>

#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"
#include "biosoup/timer.hpp"
#include "marked_edge.h"
#include "graph_constructor/extended_overlap.h"

namespace raven {

  Graph::Graph(bool checkpoints, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
      : nodes_(),
        edges_(),
        annotations_(),
        annotations_compressed_(),
        piles_(),
       // minimizers_(),
        state_manager_(GraphState::Construct_Graph),
        checkpoints_(checkpoints),
        thread_pool_(thread_pool ?
                     thread_pool :
                     std::make_shared<thread_pool::ThreadPool>(1)) {}


  void Graph::DuplicateGraph() {
    for (const auto &it: nodes_) {
      if (it == nullptr)
        continue;

      auto node = std::make_shared<Node>(it->sequence, Node::num_objects_alternate++);
      nodes_alternate_.emplace_back(node);
      if (it->pair != nullptr)
        node->pair = it->pair->alternate;
      if (node->pair != nullptr)
        node->pair->pair = node.get();
      node->color = it->color;
      it->alternate = node.get();
      node->alternate = it.get();
      node->is_primary = false;
    }

    for (const auto &it: edges_) {
      if (it == nullptr)
        continue;

      auto edge = std::make_shared<Edge>(it->tail->alternate, it->head->alternate, it->length,
                                         Edge::num_objects_alternate++);
      it->alternate = edge.get();
      edge->alternate = it.get();

      edge->pair = it->pair->alternate;
      if (edge->pair != nullptr)
        edge->pair->pair = edge.get();

      edges_alternate_.emplace_back(edge);
    }
  }

  std::uint32_t Graph::CreateUnitigs(std::uint32_t epsilon) {
    std::unordered_set<std::uint32_t> marked_edges;
    std::vector<std::shared_ptr<Node>> unitigs;
    std::vector<std::shared_ptr<Edge>> unitig_edges;
    std::vector<std::uint32_t> node_updates(nodes_.size(), 0);
    std::vector<char> is_visited(nodes_.size(), 0);

    for (const auto &it: nodes_) {
      if (it == nullptr || is_visited[it->id] || it->is_junction()) {
        continue;
      }

      std::uint32_t extension = 1;

      bool is_circular = false;
      auto begin = it.get();
      while (!begin->is_junction()) {  // extend left
        is_visited[begin->id] = 1;
        is_visited[begin->pair->id] = 1;
        if (begin->indegree() == 0 ||
            begin->inedges.front()->tail->is_junction()) {
          break;
        }
        begin = begin->inedges.front()->tail;
        ++extension;
        if (begin == it.get()) {
          is_circular = true;
          break;
        }
      }

      auto end = it.get();
      while (!end->is_junction()) {  // extend right
        is_visited[end->id] = 1;
        is_visited[end->pair->id] = 1;
        if (end->outdegree() == 0 ||
            end->outedges.front()->head->is_junction()) {
          break;
        }
        end = end->outedges.front()->head;
        ++extension;
        if (end == it.get()) {
          is_circular = true;
          break;
        }
      }

      if (!is_circular && begin == end) {
        continue;
      }
      if (!is_circular && extension < 2 * epsilon + 2) {
        continue;
      }

      if (begin != end) {  // remove nodes near junctions
        for (std::uint32_t i = 0; i < epsilon; ++i) {
          begin = begin->outedges.front()->head;
        }
        for (std::uint32_t i = 0; i < epsilon; ++i) {
          end = end->inedges.front()->tail;
        }
      }

      auto unitig = std::make_shared<Node>(begin, end);
      unitigs.emplace_back(unitig);
      unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair));
      unitig->pair = unitigs.back().get();
      unitig->pair->pair = unitig.get();

      if (begin != end) {  // connect unitig to graph
        if (begin->indegree()) {
          marked_edges.emplace(begin->inedges.front()->id);
          marked_edges.emplace(begin->inedges.front()->pair->id);

          auto edge = std::make_shared<Edge>(
              begin->inedges.front()->tail,
              unitig.get(),
              begin->inedges.front()->length);
          unitig_edges.emplace_back(edge);
          unitig_edges.emplace_back(std::make_shared<Edge>(
              unitig->pair,
              begin->inedges.front()->pair->head,
              begin->inedges.front()->pair->length + unitig->pair->sequence.inflated_len -
              begin->pair->sequence.inflated_len));  // NOLINT
          edge->pair = unitig_edges.back().get();
          edge->pair->pair = edge.get();
        }
        if (end->outdegree()) {
          marked_edges.emplace(end->outedges.front()->id);
          marked_edges.emplace(end->outedges.front()->pair->id);

          auto edge = std::make_shared<Edge>(
              unitig.get(),
              end->outedges.front()->head,
              end->outedges.front()->length + unitig->sequence.inflated_len -
              end->sequence.inflated_len);  // NOLINT
          unitig_edges.emplace_back(edge);
          unitig_edges.emplace_back(std::make_shared<Edge>(
              end->outedges.front()->pair->tail,
              unitig->pair,
              end->outedges.front()->pair->length));
          edge->pair = unitig_edges.back().get();
          edge->pair->pair = edge.get();
        }
      }

      auto jt = begin;
      while (true) {
        marked_edges.emplace(jt->outedges.front()->id);
        marked_edges.emplace(jt->outedges.front()->pair->id);

        // update transitive edges
        node_updates[jt->id & ~1UL] = unitig->id;
        unitig->transitive.insert(
            nodes_[jt->id & ~1UL]->transitive.begin(),
            nodes_[jt->id & ~1UL]->transitive.end());

        if ((jt = jt->outedges.front()->head) == end) {
          break;
        }
      }
    }

    nodes_.insert(nodes_.end(), unitigs.begin(), unitigs.end());
    edges_.insert(edges_.end(), unitig_edges.begin(), unitig_edges.end());
    RemoveEdges(marked_edges, true);

    for (const auto &it: nodes_) {  // update transitive edges
      if (it) {
        std::unordered_set<std::uint32_t> valid;
        for (auto jt: it->transitive) {
          valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
        }
        it->transitive.swap(valid);
      }
    }

    return unitigs.size() / 2;
  }

  std::uint32_t Graph::CreateUnitigsAlternate(std::uint32_t epsilon) {
    //std::unordered_set<std::uint32_t> marked_edges;
    std::vector<MarkedEdge> marked_edges;
    std::vector<std::shared_ptr<Node>> unitigs;
    std::vector<std::shared_ptr<Edge>> unitig_edges;
    std::vector<std::uint32_t> node_updates(nodes_alternate_.size(), 0);
    std::vector<char> is_visited(nodes_alternate_.size(), 0);

    for (const auto &it: nodes_alternate_) {
      if (it == nullptr || is_visited[it->id] || it->is_junction()) {
        continue;
      }

      std::uint32_t extension = 1;

      bool is_circular = false;
      auto begin = it.get();
      while (!begin->is_junction()) {  // extend left
        is_visited[begin->id] = 1;
        is_visited[begin->pair->id] = 1;
        if (begin->indegree() == 0 ||
            begin->inedges.front()->tail->is_junction()) {
          break;
        }
        begin = begin->inedges.front()->tail;
        ++extension;
        if (begin == it.get()) {
          is_circular = true;
          break;
        }
      }

      auto end = it.get();
      while (!end->is_junction()) {  // extend right
        is_visited[end->id] = 1;
        is_visited[end->pair->id] = 1;
        if (end->outdegree() == 0 ||
            end->outedges.front()->head->is_junction()) {
          break;
        }
        end = end->outedges.front()->head;
        ++extension;
        if (end == it.get()) {
          is_circular = true;
          break;
        }
      }

      if (!is_circular && begin == end) {
        continue;
      }
      if (!is_circular && extension < 2 * epsilon + 2) {
        continue;
      }

      if (begin != end) {  // remove nodes near junctions
        for (std::uint32_t i = 0; i < epsilon; ++i) {
          begin = begin->outedges.front()->head;
        }
        for (std::uint32_t i = 0; i < epsilon; ++i) {
          end = end->inedges.front()->tail;
        }
      }

      auto unitig = std::make_shared<Node>(begin, end, Node::num_objects_alternate++);
      unitigs.emplace_back(unitig);
      unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair, Node::num_objects_alternate++));
      unitig->pair = unitigs.back().get();
      unitig->pair->pair = unitig.get();

      if (begin != end) {  // connect unitig to graph
        if (begin->indegree()) {
          marked_edges.emplace_back(begin->inedges.front()->id, 2);
          marked_edges.emplace_back(begin->inedges.front()->pair->id, 2);

          auto edge = std::make_shared<Edge>(
              begin->inedges.front()->tail,
              unitig.get(),
              begin->inedges.front()->length, Edge::num_objects_alternate++);
          unitig_edges.emplace_back(edge);
          unitig_edges.emplace_back(std::make_shared<Edge>(
              unitig->pair,
              begin->inedges.front()->pair->head,
              begin->inedges.front()->pair->length + unitig->pair->sequence.inflated_len -
              begin->pair->sequence.inflated_len,
              Edge::num_objects_alternate++));  // NOLINT
          edge->pair = unitig_edges.back().get();
          edge->pair->pair = edge.get();
        }
        if (end->outdegree()) {
          marked_edges.emplace_back(end->outedges.front()->id, 2);
          marked_edges.emplace_back(end->outedges.front()->pair->id, 2);

          auto edge = std::make_shared<Edge>(
              unitig.get(),
              end->outedges.front()->head,
              end->outedges.front()->length + unitig->sequence.inflated_len - end->sequence.inflated_len,
              Edge::num_objects_alternate++);  // NOLINT
          unitig_edges.emplace_back(edge);
          unitig_edges.emplace_back(std::make_shared<Edge>(
              end->outedges.front()->pair->tail,
              unitig->pair,
              end->outedges.front()->pair->length,
              Edge::num_objects_alternate++));
          edge->pair = unitig_edges.back().get();
          edge->pair->pair = edge.get();
        }
      }

      auto jt = begin;
      while (true) {
        marked_edges.emplace_back(jt->outedges.front()->id, 2);
        marked_edges.emplace_back(jt->outedges.front()->pair->id, 2);

        // update transitive edges
        node_updates[jt->id & ~1UL] = unitig->id;
        unitig->transitive.insert(
            nodes_alternate_[jt->id & ~1UL]->transitive.begin(),
            nodes_alternate_[jt->id & ~1UL]->transitive.end());

        if ((jt = jt->outedges.front()->head) == end) {
          break;
        }
      }
    }

    nodes_alternate_.insert(nodes_alternate_.end(), unitigs.begin(), unitigs.end());
    edges_alternate_.insert(edges_alternate_.end(), unitig_edges.begin(), unitig_edges.end());
    RemoveAlternateEdges(marked_edges, true);

    for (const auto &it: nodes_alternate_) {  // update transitive edges
      if (it) {
        std::unordered_set<std::uint32_t> valid;
        for (auto jt: it->transitive) {
          valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
        }
        it->transitive.swap(valid);
      }
    }

    return unitigs.size() / 2;
  }

  void Graph::CreateUnitigGraph(const std::string &gfa_path) {
    std::unordered_set<std::uint32_t> marked_edges;
    std::vector<std::shared_ptr<Node>> unitigs;
    std::vector<std::shared_ptr<Edge>> unitig_edges;
    std::vector<std::uint32_t> node_updates(nodes_.size(), 0);
    std::vector<char> is_visited(nodes_.size(), 0);

    for (const auto &it: nodes_) {

      if (it == nullptr || is_visited[it->id]) {
        continue;
      }

      std::uint32_t extension = 1;

      auto begin = it.get();

      while (true) {  // extend left
        is_visited[begin->id] = 1;
        is_visited[begin->pair->id] = 1;

        if (begin->indegree() > 1 || begin->indegree() == 0) {
          break;
        }

        auto next = begin->inedges.front()->tail;
        if (next->outdegree() != 1) {
          break;
        }

        begin = next;
        ++extension;
        if (begin == it.get()) {
          break;
        }
      }


      auto end = it.get();
      while (true) {  // extend right
        is_visited[end->id] = 1;
        is_visited[end->pair->id] = 1;

        if (end->outdegree() > 1 || end->outdegree() == 0) {
          break;
        }

        auto next = end->outedges.front()->head;

        if (next->indegree() != 1) {
          break;
        }
        end = next;
        ++extension;
        if (end == it.get()) {
          break;
        }
      }

      auto unitig = std::make_shared<Node>(begin, end, true);
      unitigs.emplace_back(unitig);
      unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair, true));
      unitig->pair = unitigs.back().get();
      unitig->pair->pair = unitig.get();
    }

    auto GetUnitigFromNode = [&](Node *query_node) -> Node * {
      for (auto &unitig: unitigs) {
        for (auto unitig_node: unitig->unitig_nodes) {
          if (unitig_node->id == query_node->id) {
            return unitig.get();
          };
        };
      };
      return nullptr;
    };

    std::uint32_t counter = 0;
    for (auto &unitig: unitigs) {
      if (!unitig->is_circular) {  // connect unitig to graph
        counter++;
        for (auto &inedge: unitig->front_inedges) {
          marked_edges.emplace(inedge->id);
          marked_edges.emplace(inedge->pair->id);

          auto tail_unitig = GetUnitigFromNode(inedge->tail);

          auto edge = std::make_shared<Edge>(
              tail_unitig,
              unitig.get(),
              tail_unitig->sequence.inflated_len - inedge->tail->sequence.inflated_len + inedge->length
          );
          unitig_edges.emplace_back(edge);

          auto rc_head_unitig = GetUnitigFromNode(inedge->pair->head);
          auto rc_edge = std::make_shared<Edge>(
              unitig->pair,
              rc_head_unitig,
              unitig->pair->sequence.inflated_len - inedge->pair->head->sequence.inflated_len +
              inedge->length
          );
          unitig_edges.emplace_back(rc_edge);
          edge->pair = unitig_edges.back().get();
          edge->pair->pair = edge.get();
        };

        for (auto &outedge: unitig->back_outedges) {
          marked_edges.emplace(outedge->id);
          marked_edges.emplace(outedge->pair->id);

          auto head_unitig = GetUnitigFromNode(outedge->head);
          auto edge = std::make_shared<Edge>(
              unitig.get(),
              head_unitig,
              unitig->sequence.inflated_len - outedge->tail->sequence.inflated_len + outedge->length
          );
          unitig_edges.emplace_back(edge);

          auto rc_tail_unitig = GetUnitigFromNode(outedge->pair->head);
          auto rc_edge = std::make_shared<Edge>(
              rc_tail_unitig,
              unitig->pair,
              rc_tail_unitig->sequence.inflated_len - outedge->pair->head->sequence.inflated_len +
              outedge->length
          );
          unitig_edges.emplace_back(rc_edge);
          edge->pair = unitig_edges.back().get();
          edge->pair->pair = edge.get();
        };
      };

    }
    unitig_nodes_.insert(unitig_nodes_.end(), unitigs.begin(), unitigs.end());
    unitig_edges_.insert(unitig_edges_.end(), unitig_edges.begin(), unitig_edges.end());

    PrintUnitigGfa(gfa_path, true);
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> Graph::GetUnitigPairs(bool drop_unpolished) {
    //CreateUnitigs();

    biosoup::NucleicAcid::num_objects = 0;

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
    for (const auto &it: nodes_alternate_) {
      if (it == nullptr || it->is_rc() || !it->is_unitig) {
        continue;
      }
      if (drop_unpolished && !it->is_polished) {
        continue;
      }

      std::string name = it->sequence.name +
                         " LN:i:" + std::to_string(it->sequence.inflated_len) +
                         " RC:i:" + std::to_string(it->count) +
                         " XO:i:" + std::to_string(it->is_circular);

      dst.emplace_back(new biosoup::NucleicAcid(
          name,
          it->sequence.InflateData()));
    }

    return dst;
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> Graph::GetAssembledData(bool primary) {

    biosoup::NucleicAcid::num_objects = 0;
    auto node_list = nodes_;
    if (!primary)
      node_list = nodes_alternate_;

    if (node_list.size() == 0)
      return std::vector<std::unique_ptr<biosoup::NucleicAcid>>{};

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
    for (const auto &it: node_list) {
      if (it == nullptr || it->is_rc()) {
        continue;
      }

      std::string name = it->sequence.name +
                         " LN:i:" + std::to_string(it->sequence.inflated_len) +
                         " RC:i:" + std::to_string(it->count) +
                         " XO:i:" + std::to_string(it->is_circular);

      dst.emplace_back(new biosoup::NucleicAcid(
          name,
          it->sequence.InflateData()));
    }

    return dst;
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> Graph::GetUnitigs(bool drop_unpolished) {

    CreateUnitigs();

    biosoup::NucleicAcid::num_objects = 0;

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
    for (const auto &it: nodes_) {
      if (it == nullptr || it->is_rc() || !it->is_unitig) {
        continue;
      }
      if (drop_unpolished && !it->is_polished) {
        continue;
      }

      std::string name = it->sequence.name +
                         " LN:i:" + std::to_string(it->sequence.inflated_len) +
                         " RC:i:" + std::to_string(it->count) +
                         " XO:i:" + std::to_string(it->is_circular);

      dst.emplace_back(new biosoup::NucleicAcid(
          name,
          it->sequence.InflateData()));
    }

    return dst;
  }

  void Graph::PrintJson(const std::string &path) const {
    if (path.empty()) {
      return;
    }

    std::ofstream os(path);
    cereal::JSONOutputArchive archive(os);
    for (const auto &it: piles_) {
      if (it->is_invalid()) {
        continue;
      }
      archive(cereal::make_nvp(std::to_string(it->id()), *(it.get())));
    }
  }

  void Graph::PrintCsv(const std::string &path) const {
    if (path.empty()) {
      return;
    }

    std::ofstream os(path);
    for (const auto &it: nodes_) {
      if ((it == nullptr) || it->is_rc() ||
          (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
        continue;
      }
      os << it->id << " [" << it->sequence.id << "]"
         << " LN:i:" << it->sequence.inflated_len
         << " RC:i:" << it->count
         << ","
         << it->pair->id << " [" << it->pair->sequence.id << "]"
         << " LN:i:" << it->pair->sequence.inflated_len
         << " RC:i:" << it->pair->count
         << ",0,-"
         << std::endl;
    }
    for (const auto &it: edges_) {
      if (it == nullptr) {
        continue;
      }
      os << it->tail->id << " [" << it->tail->sequence.id << "]"
         << " LN:i:" << it->tail->sequence.inflated_len
         << " RC:i:" << it->tail->count
         << ","
         << it->head->id << " [" << it->head->sequence.id << "]"
         << " LN:i:" << it->head->sequence.inflated_len
         << " RC:i:" << it->head->count
         << ",1,"
         << it->id << " " << it->length << " " << it->weight
         << std::endl;
    }
    for (const auto &it: nodes_) {  // circular edges TODO(rvaser): check
      if (it == nullptr || !it->is_circular) {
        continue;
      }
      os << it->id << " [" << it->sequence.id << "]"
         << " LN:i:" << it->sequence.inflated_len
         << " RC:i:" << it->count
         << ","
         << it->id << " [" << it->sequence.id << "]"
         << " LN:i:" << it->sequence.inflated_len
         << " RC:i:" << it->count
         << ",1,-"
         << std::endl;
    }
    os.close();
  }

  void Graph::PrintGfa(const std::string &path, const bool print_seq) const {
    if (path.empty()) {
      return;
    }

    std::ofstream os(path);
    for (const auto &it: nodes_) {
      // if ((it == nullptr) || it->is_rc()){
      //   continue;
      //   }
      // if ((it == nullptr)){
      //   continue;
      // }
      if ((it == nullptr) || it->is_rc() ||
          (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
        continue;
      }
      os << "S\t" << it->sequence.name
         << "\t" << (print_seq ? it->sequence.InflateData() : "*")
         << "\tLN:i:" << it->sequence.inflated_len
         << "\tRC:i:" << it->count
         << "\tCL:z:" << (it->color ? "blue" : "orange")
         << std::endl;
      if (it->is_circular) {
        os << "L\t" << it->sequence.name << "\t" << '+'
           << "\t" << it->sequence.name << "\t" << '+'
           << "\t0M"
           << std::endl;
      }
    }
    for (const auto &it: edges_) {
   if (it == nullptr || it->is_rc()) {
     continue;
    }
      // if (it == nullptr) {
      //   continue;
      // }
      os << "L\t" << it->tail->sequence.name << "\t" << (it->tail->is_rc() ? '-' : '+')  // NOLINT
         << "\t" << it->head->sequence.name << "\t" << (it->head->is_rc() ? '-' : '+')  // NOLINT
         << "\t" << it->tail->sequence.inflated_len - it->length << 'M'
         << std::endl;
    }
    os.close();
  }

  void Graph::PrintOverlaps(std::vector<std::vector<extended_overlap>> overlaps,
                            std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, bool print_cigar,
                            const std::string &path) const {
    if (path.empty()) {
      return;
    }
    
    std::ofstream os(path);
    os << std::defaultfloat
    << std::setprecision(std::numeric_limits<double>::max_digits10);

    for (const auto &it: overlaps) {
      for (const auto &jt: it) {
        os << sequences[jt.overlap.lhs_id]->name
           << "\t" << sequences[jt.overlap.lhs_id]->inflated_len  // length
           << "\t" << jt.overlap.lhs_begin
           << "\t" << jt.overlap.lhs_end
           << "\t" << (jt.overlap.strand ? "+" : "-")
           << "\t" << sequences[jt.overlap.rhs_id]->name
           << "\t" << sequences[jt.overlap.rhs_id]->inflated_len  // length
           << "\t" << jt.overlap.rhs_begin
           << "\t" << jt.overlap.rhs_end
           << "\t" << jt.edlib_alignment.matches // residue matches
           << "\t" << jt.edlib_alignment.block_length // alignment block length
           << "\t" << jt.overlap.score
           << "\t" << jt.graph_overlap_type
           << "\t" << jt.lhs_begin_original
           << "\t" << jt.lhs_end_original
           << "\t" << jt.rhs_begin_original
           << "\t" << jt.rhs_end_original
           << "\t" << jt.extended_length
           << "\t" << jt.found_length
           << "\t" << jt.found_to_extended_length
           << "\t" << jt.score_to_length
           << "\t" << jt.lhs_hap
           << "\t" << jt.rhs_hap
           << "\t" << jt.lhs_err
           << "\t" << jt.rhs_err
           << "\t" << jt.lhs_rep
           << "\t" << jt.rhs_rep
           << "\t" << jt.hap_ratio
           << "\t" << (jt.q_hor ? "1" : "0")
           << "\t" << (jt.t_hor ? "1" : "0")
           << "\t" << jt.jaccard_index
           << "\t" << jt.diploid_jaccard_full_overlap
           << "\t" << jt.diploid_jaccard_longer
           << "\t" << jt.diploid_jaccard_shorter
           << "\t" << jt.longer_overhang_len
           << "\t" << jt.shorter_overhang_len
           << "\t" << jt.jaccard_non_extended_overlap
          //  << "\t" << jt.hap_rate_shorter_overhang_lhs
          //  << "\t" << jt.hap_rate_longer_overhang_lhs
          //  << "\t" << jt.hap_rate_shorter_overhang_rhs
          //  << "\t" << jt.hap_rate_longer_overhang_rhs
          //  << "\t" << jt.dip_rate_shorter_overhang_lhs
          //  << "\t" << jt.dip_rate_longer_overhang_lhs
          //  << "\t" << jt.dip_rate_shorter_overhang_rhs
          //  << "\t" << jt.dip_rate_longer_overhang_rhs
          //  << "\t" << jt.err_rate_shorter_overhang_lhs
          //  << "\t" << jt.err_rate_longer_overhang_lhs
          //  << "\t" << jt.err_rate_shorter_overhang_rhs
          //  << "\t" << jt.err_rate_longer_overhang_rhs
//            << "\t" << (piles_[jt.overlap.lhs_id]->is_strong_contained() ? "1" : "0")
//            << "\t" << (piles_[jt.overlap.rhs_id]->is_strong_contained() ? "1" : "0")
//            << "\t" << (piles_[jt.overlap.lhs_id]->is_backbone() ? "1" : "0")
//            << "\t" << (piles_[jt.overlap.rhs_id]->is_backbone() ? "1" : "0")
//  //          << "\t" << static_cast<std::uint8_t>(jt.cat)   
//            << "\t" << (jt.candidate_overlap ? "1" : "0")
           << "\t" << jt.classification_label
           << std::endl;
      }
    }

    os.close();
  }

CorrectionPafText Graph::FormatCorrectionOverlaps(
    const std::vector<extended_overlap>& pile,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) const {

  // ========================================================================
  // INTERNAL EDLIB ORIENTATION
  //
  //   query  = lhs = backbone
  //   target = rhs = supporting read
  //
  // OUTPUT PAF ORIENTATION
  //
  //   query  = rhs = supporting read
  //   target = lhs = backbone
  //
  // Therefore:
  //
  //   Edlib I -> PAF D
  //   Edlib D -> PAF I
  //   M       -> M
  //
  // Terminal I/D operations should already have been removed inside
  // LoadOverlapsAndAlign(), with the corresponding overlap coordinates
  // adjusted there.
  // ========================================================================

  auto transpose_cigar =
      [](const std::string& cigar) -> std::string {

    std::string result;
    result.reserve(cigar.size());

    for (char c : cigar) {

      if (c == 'I') {
        result += 'D';

      } else if (c == 'D') {
        result += 'I';

      } else {
        result += c;
      }
    }

    return result;
  };


  auto first_cigar_op =
      [](const std::string& cigar) -> char {

    for (char c : cigar) {

      if (!std::isdigit(
              static_cast<unsigned char>(c))) {

        return c;
      }
    }

    return '\0';
  };


  auto last_cigar_op =
      [](const std::string& cigar) -> char {

    for (auto it = cigar.rbegin();
         it != cigar.rend();
         ++it) {

      if (!std::isdigit(
              static_cast<unsigned char>(*it))) {

        return *it;
      }
    }

    return '\0';
  };


  CorrectionPafText result;
  std::ostringstream os;
  std::ostringstream log;
  std::uint64_t& written = result.written;
  std::uint64_t& skipped = result.skipped;

  for (const auto& jt : pile) {


    // ====================================================================
    // Only output overlaps whose alignment + terminal trimming succeeded.
    // ====================================================================

    if (!jt.correction_overlap) {

      ++skipped;
      continue;
    }


    // ====================================================================
    // PAF coordinates
    //
    // query  = rhs/support
    // target = lhs/backbone
    // ====================================================================

    const std::uint32_t qs =
        jt.overlap.rhs_begin;

    const std::uint32_t qe =
        jt.overlap.rhs_end;

    const std::uint32_t ts =
        jt.overlap.lhs_begin;

    const std::uint32_t te =
        jt.overlap.lhs_end;


    const std::uint32_t qlen =
        sequences[
            jt.overlap.rhs_id
        ]->inflated_len;

    const std::uint32_t tlen =
        sequences[
            jt.overlap.lhs_id
        ]->inflated_len;


    // ====================================================================
    // Coordinate validation
    // ====================================================================

    if (qe <= qs ||
        te <= ts ||
        qe > qlen ||
        te > tlen) {

      log
          << "Invalid correction overlap coordinates:"
          << " qs=" << qs
          << " qe=" << qe
          << " qlen=" << qlen
          << " ts=" << ts
          << " te=" << te
          << " tlen=" << tlen
          << std::endl;

      ++skipped;
      continue;
    }


    // ====================================================================
    // Transpose CIGAR into output PAF orientation
    // ====================================================================

    const std::string fixed_cigar =
        transpose_cigar(
            jt.edlib_alignment.cigar);


    if (fixed_cigar.empty()) {

      ++skipped;
      continue;
    }


    // ====================================================================
    // HERRO safety:
    //
    // terminal I/D should already have been removed.
    //
    // We require the final PAF CIGAR to begin and end with an aligned
    // operation.
    // ====================================================================

    const char first_op =
        first_cigar_op(
            fixed_cigar);

    const char last_op =
        last_cigar_op(
            fixed_cigar);


    if ((first_op != 'M' &&
         first_op != '=' &&
         first_op != 'X') ||
        (last_op != 'M' &&
         last_op != '=' &&
         last_op != 'X')) {

      log
          << "Terminal gap survived trimming:"
          << "\n  query  = "
          << sequences[jt.overlap.rhs_id]->name
          << "\n  target = "
          << sequences[jt.overlap.lhs_id]->name
          << "\n  strand = "
          << (jt.overlap.strand ? '+' : '-')
          << "\n  cigar  = "
          << fixed_cigar.substr(
                 0,
                 std::min<std::size_t>(
                     fixed_cigar.size(),
                     200))
          << std::endl;

      ++skipped;
      continue;
    }


    // ====================================================================
    // Verify that CIGAR consumption exactly matches PAF coordinates.
    //
    // PAF semantics:
    //
    //   M / = / X : consume query + target
    //   I         : consume query
    //   D         : consume target
    // ====================================================================

    std::uint64_t cigar_query_consumed = 0;
    std::uint64_t cigar_target_consumed = 0;

    std::uint64_t num = 0;


    for (char c : fixed_cigar) {

      if (std::isdigit(
              static_cast<unsigned char>(c))) {

        num =
            num * 10 +
            static_cast<std::uint64_t>(
                c - '0');

        continue;
      }


      if (c == 'M' ||
          c == '=' ||
          c == 'X') {

        cigar_query_consumed += num;
        cigar_target_consumed += num;

      } else if (c == 'I') {

        cigar_query_consumed += num;

      } else if (c == 'D') {

        cigar_target_consumed += num;

      } else {

        log
            << "Unexpected CIGAR operation: "
            << c
            << std::endl;

        num = 0;
        continue;
      }


      num = 0;
    }


    const std::uint64_t query_span =
        static_cast<std::uint64_t>(
            qe - qs);

    const std::uint64_t target_span =
        static_cast<std::uint64_t>(
            te - ts);


    if (cigar_query_consumed != query_span ||
        cigar_target_consumed != target_span) {

      log
          << "CIGAR/PAF span mismatch:"
          << "\n  query span      = "
          << query_span
          << "\n  cigar query     = "
          << cigar_query_consumed
          << "\n  target span     = "
          << target_span
          << "\n  cigar target    = "
          << cigar_target_consumed
          << "\n  strand          = "
          << (jt.overlap.strand ? '+' : '-')
          << "\n  query           = "
          << sequences[jt.overlap.rhs_id]->name
          << "\n  target          = "
          << sequences[jt.overlap.lhs_id]->name
          << "\n  cigar           = "
          << fixed_cigar.substr(
                 0,
                 std::min<std::size_t>(
                     fixed_cigar.size(),
                     200))
          << std::endl;

      ++skipped;
      continue;
    }


    // ====================================================================
    // OUTPUT PAF
    // ====================================================================

    os
        // Query = supporting read
        << sequences[
               jt.overlap.rhs_id
           ]->name
        << '\t'

        << qlen
        << '\t'

        << qs
        << '\t'

        << qe
        << '\t'

        // Strand relative to backbone
        << (jt.overlap.strand
                ? '+'
                : '-')
        << '\t'

        // Target = backbone
        << sequences[
               jt.overlap.lhs_id
           ]->name
        << '\t'

        << tlen
        << '\t'

        << ts
        << '\t'

        << te
        << '\t'

        // Alignment statistics
        << jt.edlib_alignment.matches
        << '\t'

        << jt.edlib_alignment.block_length
        << '\t'

        << 255
        << '\t'

        // HERRO expects cg:Z as final field
        << "cg:Z:"
        << fixed_cigar
        << '\n';


    ++written;
  }

  result.paf = os.str();
  result.log = log.str();
  return result;

}

void Graph::PrintCorrectionOverlaps(
    const std::vector<std::vector<extended_overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::string& path,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) const {

  if (path.empty()) {
    return;
  }


  biosoup::Timer timer{};
  timer.Start();

  std::ofstream out(path);

  if (!out) {

    std::cerr
        << "Could not open correction overlap file: "
        << path
        << std::endl;

    return;
  }


  auto format_pile = [&](std::size_t pile_id) -> CorrectionPafText {
    return FormatCorrectionOverlaps(overlaps[pile_id], sequences);
  };


  std::uint64_t written = 0;
  std::uint64_t skipped = 0;

  auto emit = [&](CorrectionPafText&& text) {
    out << text.paf;
    std::cerr << text.log;
    written += text.written;
    skipped += text.skipped;
  };

  if (thread_pool == nullptr) {
    for (std::size_t i = 0; i < overlaps.size(); ++i) {
      emit(format_pile(i));
    }
  } else {
    // Formatted in parallel but written in pile order, so the file is the
    // same as a serial write; a bounded window of piles is held at once.
    const std::size_t window = 4 * thread_pool->num_threads();
    std::deque<std::future<CorrectionPafText>> pending;
    std::size_t next = 0;
    while (next < overlaps.size() || !pending.empty()) {
      while (next < overlaps.size() && pending.size() < window) {
        pending.emplace_back(thread_pool->Submit(format_pile, next++));
      }
      emit(pending.front().get());
      pending.pop_front();
    }
  }


  out.close();


  std::cerr
      << "[raven::Graph::PrintCorrectionOverlaps] printed "
      << written
      << " overlaps, skipped "
      << skipped
      << " "
      << std::fixed << timer.Stop() << "s"
      << std::endl;
}
/*
void Graph::PrintCorrectionOverlaps(
    const std::vector<std::vector<extended_overlap>>& overlaps,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::string& path) const {

  if (path.empty()) {
    return;
  }


  // ========================================================================
  // Convert the Edlib CIGAR orientation to the PAF orientation.
  //
  // Edlib was calculated as:
  //
  //     query  = lhs = backbone
  //     target = rhs = support read
  //
  // PAF is written as:
  //
  //     query  = rhs = support read
  //     target = lhs = backbone
  //
  // Therefore:
  //
  //     Edlib I -> PAF D
  //     Edlib D -> PAF I
  //
  // M remains M.
  //
  // No coordinate trimming is needed here because Edlib NW aligned exactly
  // the complete intervals represented by overlap.{lhs,rhs}_{begin,end}.
  // ========================================================================

  auto transpose_cigar =
      [](const std::string& cigar) -> std::string {

    std::string result;
    result.reserve(cigar.size());

    for (char c : cigar) {

      if (c == 'I') {
        result += 'D';

      } else if (c == 'D') {
        result += 'I';

      } else {
        result += c;
      }
    }

    return result;
  };


  std::ofstream os(path);

  if (!os) {
    std::cerr
        << "Could not open correction overlap file: "
        << path
        << std::endl;

    return;
  }


  std::uint64_t written = 0;
  std::uint64_t skipped = 0;


  for (const auto& it : overlaps) {

    for (const auto& jt : it) {

      // ====================================================================
      // Output orientation:
      //
      // query  = rhs = supporting read
      // target = lhs = backbone
      // ====================================================================

      const std::uint32_t qs =
          jt.overlap.rhs_begin;

      const std::uint32_t qe =
          jt.overlap.rhs_end;

      const std::uint32_t ts =
          jt.overlap.lhs_begin;

      const std::uint32_t te =
          jt.overlap.lhs_end;


      // ====================================================================
      // Basic coordinate validation
      // ====================================================================

      if (qe <= qs ||
          te <= ts) {

        ++skipped;
        continue;
      }


      const std::uint32_t qlen =
          sequences[jt.overlap.rhs_id]->inflated_len;

      const std::uint32_t tlen =
          sequences[jt.overlap.lhs_id]->inflated_len;


      if (qe > qlen ||
          te > tlen) {

        std::cerr
            << "Invalid correction overlap coordinates:"
            << " qs=" << qs
            << " qe=" << qe
            << " qlen=" << qlen
            << " ts=" << ts
            << " te=" << te
            << " tlen=" << tlen
            << std::endl;

        ++skipped;
        continue;
      }


      // ====================================================================
      // Transpose Edlib CIGAR into PAF query/target orientation
      // ====================================================================

      const std::string fixed_cigar =
          transpose_cigar(
              jt.edlib_alignment.cigar);


      if (fixed_cigar.empty()) {

        ++skipped;
        continue;
      }


      // ====================================================================
      // DEBUG / SANITY CHECK
      //
      // Verify that the CIGAR consumes exactly the PAF intervals.
      //
      // PAF:
      //     M/I consume query
      //     M/D consume target
      // ====================================================================

      std::uint64_t cigar_query_consumed = 0;
      std::uint64_t cigar_target_consumed = 0;

      std::uint64_t num = 0;

      for (char c : fixed_cigar) {

        if (std::isdigit(
                static_cast<unsigned char>(c))) {

          num =
              num * 10 +
              static_cast<std::uint64_t>(
                  c - '0');

          continue;
        }


        if (c == 'M' ||
            c == '=' ||
            c == 'X') {

          cigar_query_consumed += num;
          cigar_target_consumed += num;

        } else if (c == 'I') {

          cigar_query_consumed += num;

        } else if (c == 'D') {

          cigar_target_consumed += num;
        }


        num = 0;
      }


      const std::uint64_t query_span =
          static_cast<std::uint64_t>(qe) - qs;

      const std::uint64_t target_span =
          static_cast<std::uint64_t>(te) - ts;


      if (cigar_query_consumed != query_span ||
          cigar_target_consumed != target_span) {

        std::cerr
            << "CIGAR/PAF span mismatch:"
            << "\n  query span  = "
            << query_span
            << ", cigar consumes = "
            << cigar_query_consumed
            << "\n  target span = "
            << target_span
            << ", cigar consumes = "
            << cigar_target_consumed
            << "\n  strand = "
            << (jt.overlap.strand ? '+' : '-')
            << "\n  cigar = "
            << fixed_cigar.substr(
                   0,
                   std::min<std::size_t>(
                       fixed_cigar.size(),
                       200))
            << std::endl;

        ++skipped;
        continue;
      }


      // ====================================================================
      // Output PAF
      // ====================================================================

      os
          // Query = supporting read
          << sequences[jt.overlap.rhs_id]->name
          << '\t'
          << qlen
          << '\t'
          << qs
          << '\t'
          << qe
          << '\t'

          // Relative strand
          << (jt.overlap.strand ? '+' : '-')
          << '\t'

          // Target = backbone
          << sequences[jt.overlap.lhs_id]->name
          << '\t'
          << tlen
          << '\t'
          << ts
          << '\t'
          << te
          << '\t'

          // Alignment statistics
          << jt.edlib_alignment.matches
          << '\t'
          << jt.edlib_alignment.block_length
          << '\t'
          << 255
          << '\t'

          // Keep cg:Z as LAST field for HERRO
          << "cg:Z:"
          << fixed_cigar
          << '\n';


      ++written;
    }
  }


  os.close();


  std::cerr
      << "Printed "
      << written
      << " correction overlaps, skipped "
      << skipped
      << std::endl;
}*/

 /*void Graph::PrintCorrectionOverlaps(
      const std::vector<std::vector<extended_overlap>>& overlaps,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
      const std::string& path) const {

    if (path.empty()) return;

    auto fix_cigar = [](std::string cigar,
                        std::uint32_t& qstart, std::uint32_t& qend,
                        std::uint32_t& tstart, std::uint32_t& tend,
                        std::uint32_t& trimmed_len) -> std::string {
        // parse + swap I<->D
        std::vector<std::pair<int,char>> ops;
        int num = 0;
        for (char c : cigar) {
            if (std::isdigit(c)) { num = num*10 + (c-'0'); }
            else {
                char op = c;
                if (op == 'I') op = 'D';
                else if (op == 'D') op = 'I';
                ops.emplace_back(num, op);
                num = 0;
            }
        }
        trimmed_len = 0;
        // trim leading non-M
        size_t lo = 0;
        while (lo < ops.size() && ops[lo].second != 'M') {
            if (ops[lo].second == 'I') qstart += ops[lo].first;
            else if (ops[lo].second == 'D') tstart += ops[lo].first;
            trimmed_len += ops[lo].first;
            lo++;
        }
        // trim trailing non-M
        size_t hi = ops.size();
        while (hi > lo && ops[hi-1].second != 'M') {
            if (ops[hi-1].second == 'I') qend -= ops[hi-1].first;
            else if (ops[hi-1].second == 'D') tend -= ops[hi-1].first;
            trimmed_len += ops[hi-1].first;
            hi--;
        }
        std::string result;
        for (size_t k = lo; k < hi; k++) {
            result += std::to_string(ops[k].first) + ops[k].second;
        }
        return result;
    };

    std::ofstream os(path);

    for (const auto& it : overlaps) {
      for (const auto& jt : it) {

        if (!jt.correction_overlap) continue;
        std::uint32_t qs = jt.rhs_begin_original, qe = jt.rhs_end_original;
        std::uint32_t ts = jt.lhs_begin_original, te = jt.lhs_end_original;
        
        if(jt.classification_label == 1){
          qs = jt.overlap.rhs_begin;
          qe = jt.overlap.rhs_end;
          ts = jt.overlap.lhs_begin;
          te = jt.overlap.lhs_end;
        }
        std::uint32_t trimmed_len = 0;
        std::string fixed_cigar = fix_cigar(jt.edlib_alignment.cigar, qs, qe, ts, te, trimmed_len);

        os << sequences[jt.overlap.rhs_id]->name << '\t'
          << sequences[jt.overlap.rhs_id]->inflated_len << '\t'
          << qs << '\t' << qe << '\t'
          << (jt.overlap.strand ? '+' : '-') << '\t'
          << sequences[jt.overlap.lhs_id]->name << '\t'
          << sequences[jt.overlap.lhs_id]->inflated_len << '\t'
          << ts << '\t' << te << '\t'
          << jt.edlib_alignment.matches << '\t'
          << (jt.edlib_alignment.block_length - trimmed_len) << '\t'
          << 255 << '\t'
          << "cg:Z:" << fixed_cigar << '\n';
      }
    }

    os.close();
  }

  void Graph::PrintCorrectionOverlapsWithFeatures(
      const std::vector<std::vector<extended_overlap>>& overlaps,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
      const std::string& path) const {

    std::ofstream os(path);

    for (const auto& it : overlaps) {
      for (const auto& jt : it) {

        if (!jt.correction_overlap) continue;

        os << sequences[jt.overlap.rhs_id]->name << '\t'
          << sequences[jt.overlap.rhs_id]->inflated_len << '\t'
          << jt.overlap.rhs_begin << '\t' << jt.overlap.rhs_end << '\t'
          << (jt.overlap.strand ? '+' : '-') << '\t'
          << sequences[jt.overlap.lhs_id]->name << '\t'
          << sequences[jt.overlap.lhs_id]->inflated_len << '\t'
          << jt.overlap.lhs_begin << '\t' << jt.overlap.lhs_end << '\t'
          << jt.found_matches << '\t'
          << jt.extended_length << '\t'
          << 255 << '\t';
          //<< jt.;
      }
    }

    os.close();
  }*/ 

  void Graph::PrintBackboneOverlapsWithFeatures(
      const std::vector<std::vector<extended_overlap>>& overlaps,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
      const std::string& path) const {
    if (path.empty()) return;

    std::ofstream os(path);
    os << std::defaultfloat
      << std::setprecision(std::numeric_limits<double>::max_digits10);

    for (std::size_t i = 0; i < overlaps.size(); ++i) {
      // Only reads that will be corrected:
      // backbone and not strong-contained.
      //
      // These are lhs reads internally, but they are printed as
      // PAF TARGETS because they are the reads to be corrected.
      if (i >= piles_.size() ||
          !piles_[i]->is_backbone() || piles_[i]->is_strong_contained()) {
        continue;
      }

      for (const auto& jt : overlaps[i]) {
        // PAF:
        //   query  = rhs (supporting read)
        //   target = lhs (read to correct)
        os << sequences[jt.overlap.rhs_id]->name
          << "\t" << sequences[jt.overlap.rhs_id]->inflated_len
          << "\t" << jt.overlap.rhs_begin
          << "\t" << jt.overlap.rhs_end
          << "\t" << (jt.overlap.strand ? "+" : "-")
          << "\t" << sequences[jt.overlap.lhs_id]->name
          << "\t" << sequences[jt.overlap.lhs_id]->inflated_len
          << "\t" << jt.overlap.lhs_begin
          << "\t" << jt.overlap.lhs_end
          << "\t" << jt.edlib_alignment.matches
          << "\t" << jt.edlib_alignment.block_length
          << "\t" << jt.overlap.score
          << "\t" << jt.graph_overlap_type

          // Swap lhs/rhs-specific feature fields as well so that
          // they follow query/target order in the printed record.
          << "\t" << jt.rhs_begin_original
          << "\t" << jt.rhs_end_original
          << "\t" << jt.lhs_begin_original
          << "\t" << jt.lhs_end_original

          << "\t" << jt.extended_length
          << "\t" << jt.found_length
          << "\t" << jt.found_to_extended_length
          << "\t" << jt.score_to_length

          << "\t" << jt.rhs_hap
          << "\t" << jt.lhs_hap

          << "\t" << jt.rhs_err
          << "\t" << jt.lhs_err

          << "\t" << jt.rhs_rep
          << "\t" << jt.lhs_rep

          << "\t" << jt.hap_ratio

          // q = rhs, t = lhs now
          << "\t" << (jt.t_hor ? "1" : "0")
          << "\t" << (jt.q_hor ? "1" : "0")

          << "\t" << jt.jaccard_index
          << "\t" << jt.diploid_jaccard_full_overlap
          << "\t" << jt.diploid_jaccard_longer
          << "\t" << jt.diploid_jaccard_shorter
          << "\t" << jt.longer_overhang_len
          << "\t" << jt.shorter_overhang_len
          << "\t" << jt.jaccard_non_extended_overlap
          << "\t" << jt.classification_label
          << '\n';
      }
    }

    os.close();
  }

  void Graph::PrintAllRawOverlaps(
      const std::vector<std::vector<extended_overlap>>& overlaps,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
      const std::string& path) const {
    if (path.empty()) return;

    std::ofstream os(path);

    // Standard PAF columns using the raw ram anchor coordinates (no extension,
    // no dedup, no merge). Every fragment ram found is printed once, oriented
    // with read i (lhs) as the query. Columns 10/11 carry the anchor count
    // (score) and the raw span length as stand-ins for matches/block_length;
    // only columns 1-9 (the coordinates) are meaningful for the minimap2
    // comparison.
    for (const auto& it : overlaps) {
      for (const auto& jt : it) {
        const std::uint32_t q_id = jt.overlap.lhs_id;
        const std::uint32_t t_id = jt.overlap.rhs_id;
        const std::uint32_t qs = jt.lhs_begin_original;
        const std::uint32_t qe = jt.lhs_end_original;
        const std::uint32_t ts = jt.rhs_begin_original;
        const std::uint32_t te = jt.rhs_end_original;

        os << sequences[q_id]->name << '\t'
           << sequences[q_id]->inflated_len << '\t'
           << qs << '\t' << qe << '\t'
           << (jt.overlap.strand ? '+' : '-') << '\t'
           << sequences[t_id]->name << '\t'
           << sequences[t_id]->inflated_len << '\t'
           << ts << '\t' << te << '\t'
           << jt.found_matches << '\t'
           << jt.found_length << '\t'
           << 255 << '\n';
      }
    }

    os.close();
  }
  // void Graph::PrintCorrectionOverlaps(
  //     const std::vector<std::vector<extended_overlap>>& overlaps,
  //     const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
  //     const std::string& path) const {

  //   if (path.empty()) return;

  //   std::ofstream os(path);

  //   for (const auto& it : overlaps) {
  //     for (const auto& jt : it) {

  //       if (!jt.correction_overlap) continue;

  //       os << sequences[jt.overlap.lhs_id]->name         << '\t'  // query name
  //         << sequences[jt.overlap.lhs_id]->inflated_len << '\t'  // query length
  //         << jt.overlap.lhs_begin                       << '\t'  // query start
  //         << jt.overlap.lhs_end                         << '\t'  // query end
  //         << (jt.overlap.strand ? '+' : '-')            << '\t'  // strand
  //         << sequences[jt.overlap.rhs_id]->name         << '\t'  // target name
  //         << sequences[jt.overlap.rhs_id]->inflated_len << '\t'  // target length
  //         << jt.overlap.rhs_begin                       << '\t'  // target start
  //         << jt.overlap.rhs_end                         << '\t'  // target end
  //         << jt.edlib_alignment.matches                 << '\t'  // matches
  //         << jt.edlib_alignment.block_length            << '\t'  // block length
  //         << 255                                        << '\t'  // mapq
  //         << "cg:Z:" << jt.edlib_alignment.cigar        << '\n';
  //     }
  //   }

  //   os.close();
  // }
  void Graph::PrintPiles(){
// for (const auto &it : graph_.piles_) {
    //   std::ofstream sketch_out(sequences[it->id()]->name + ".csv");
    //   std::cerr << "Sketching: " << sequences[it->id()]->name << std::endl;
    //   if(sequences[it->id()]->inflated_len < 500){
    //     std::cerr << "Skipped sketching: " << sequences[it->id()]->name << std::endl;
    //     continue;
    //   }

    //    // auto sketch_results = it->get_kmer_types();
    //     // auto kmer_ids = it->get_k_kmer_ids();
    //     // auto multiplicity_data = it->get_sketch_data();
    //     // for(int i = 0; i < multiplicity_data.size(); i++){
    //     //   sketch_out << kmer_ids[i] << "\t" << multiplicity_data[i] << std::endl;
    //     // };
    //     auto sketch_results = it->get_kmer_types();
    //     for(auto rez : sketch_results){
    //       switch (rez)
    //       {
    //       case kMerType::None:
    //         sketch_out << "N" << std::endl;
    //         break;
    //       case kMerType::Haploid:
    //         sketch_out << "H" << std::endl;
    //         break;
    //       case kMerType::Diploid:
    //         sketch_out << "D"  << std::endl;
    //         break;
    //       case kMerType::Repetitive:
    //         sketch_out << "R"  << std::endl;
    //         break;
    //       case kMerType::Error:
    //         sketch_out << "E"  << std::endl;
    //         break;
    //       default:
    //         break;
    //       }
    //     }
    //     // auto avg_quality_data = it->get_avg_base_qualities();
    //     // auto min_quality_data = it->get_min_base_qualities();
    //     // int safe_count = std::min({
    //     //     (int)multiplicity_data.size(),
    //     //     (int)sketch_results.size() / 10,
    //     //     (int)kmer_ids.size() / 10
    //     // });
    //     // for(int i = 0; i < safe_count; i++) {
    //     //   switch (sketch_results[i*10])
    //     //   {
    //     //   case kMerType::None:
    //     //     sketch_out << "N" << "\t" << multiplicity_data[i] << "\t" << avg_quality_data[i] << "\t" << min_quality_data[i] << "\t" << kmer_ids[i*10] << std::endl;
    //     //     break;
    //     //   case kMerType::Haploid:
    //     //     sketch_out << "H" << "\t" << multiplicity_data[i] << "\t" << avg_quality_data[i] << "\t" << min_quality_data[i] << "\t" << kmer_ids[i*10] << std::endl;
    //     //     break;
    //     //   case kMerType::Diploid:
    //     //     sketch_out << "D" << "\t" << multiplicity_data[i] << "\t" << avg_quality_data[i] << "\t" << min_quality_data[i] << "\t" << kmer_ids[i*10] << std::endl;
    //     //     break;
    //     //   case kMerType::Repetitive:
    //     //     sketch_out << "R" << "\t" << multiplicity_data[i] << "\t" << avg_quality_data[i] << "\t" << min_quality_data[i] << "\t" << kmer_ids[i*10] << std::endl;
    //     //     break;
    //     //   case kMerType::Error:
    //     //     sketch_out << "E" << "\t" << multiplicity_data[i] << "\t" << avg_quality_data[i] << "\t" << min_quality_data[i] << "\t" << kmer_ids[i*10] << std::endl;
    //     //     break;
    //     //   default:
    //     //     break;
    //     //   }
    //     // }
    //   sketch_out.close();
    // }
  };

  void Graph::PrintUnitigGfa(const std::string &path, const bool print_seq) const {
    if (path.empty()) {
      return;
    }
    std::ofstream os(path);
    for (const auto &it: unitig_nodes_) {
      if ((it == nullptr) || it->is_rc() ||
          (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
        continue;
      }

      os << "S\t" << it->sequence.name
         << "\t" << (print_seq ? it->sequence.InflateData() : "*")
         << "\tLN:i:" << it->sequence.inflated_len
         << "\tRC:i:" << it->count
         << "\tCL:z:" << (it->color ? "blue" : "orange")
         << std::endl;
      if (it->is_circular) {
        os << "A\t" << it->sequence.name << "\t"
           << it->sequence.name << "\t"
           << std::endl;
        continue;
      }
      for (const auto &unitig_node: it->unitig_nodes) {
        os << "A\t" << it->sequence.name << "\t"
           << unitig_node->sequence.name << "\t"
           << std::endl;
      };
    }
    for (const auto &it: unitig_edges_) {
      if (it == nullptr || it->is_rc()) {
        continue;
      }
      os << "L\t" << it->tail->sequence.name << "\t" << (it->tail->is_rc() ? '-' : '+')  // NOLINT
         << "\t" << it->head->sequence.name << "\t" << (it->head->is_rc() ? '-' : '+')  // NOLINT
         << "\t" << it->tail->sequence.inflated_len - it->length << 'M'
         << std::endl;
    }
    os.close();
  }

  void Graph::Store(const std::string &cereal_path) const {
    std::ofstream os(cereal_path);
    try {
      cereal::BinaryOutputArchive archive(os);
      archive(*this);
    } catch (std::exception &) {
      throw std::logic_error("[raven::Graph::Store] error: unable to store archive");
    }
  }

  void Graph::Load(const std::string &cereal_path) {
    std::ifstream is(cereal_path);
    try {
      cereal::BinaryInputArchive archive(is);
      archive(*this);
    } catch (std::exception &) {
      throw std::logic_error("[raven::Graph::Load] error: unable to load archive");
    }
  }

  std::unordered_set<std::uint32_t> Graph::FindRemovableEdges(const std::vector<Node *> &path) {
    if (path.empty()) {
      return std::unordered_set<std::uint32_t>{};
    }

    auto find_edge = [](Node *tail, Node *head) -> Edge * {
      for (auto it: tail->outedges) {
        if (it->head == head) {
          return it;
        }
      }
      return nullptr;
    };

    // find first node with multiple in edges
    std::int32_t pref = -1;
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
      if (path[i]->indegree() > 1) {
        pref = i;
        break;
      }
    }

    // find last node with multiple out edges
    std::int32_t suff = -1;
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
      if (path[i]->outdegree() > 1) {
        suff = i;
      }
    }

    std::unordered_set<std::uint32_t> dst;
    if (pref == -1 && suff == -1) {  // remove whole path
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        auto it = find_edge(path[i], path[i + 1]);
        dst.emplace(it->id);
        dst.emplace(it->pair->id);
      }
      return dst;
    }

    if (pref != -1 && path[pref]->outdegree() > 1) {  // complex path
      return dst;  // empty
    }
    if (suff != -1 && path[suff]->indegree() > 1) {  // complex path
      return dst;  // empty
    }

    if (pref == -1) {  // remove everything after last suffix node
      for (std::uint32_t i = suff; i < path.size() - 1; ++i) {
        auto it = find_edge(path[i], path[i + 1]);
        dst.emplace(it->id);
        dst.emplace(it->pair->id);
      }
    } else if (suff == -1) {  // remove everything before first prefix node
      for (std::int32_t i = 0; i < pref; ++i) {
        auto it = find_edge(path[i], path[i + 1]);
        dst.emplace(it->id);
        dst.emplace(it->pair->id);
      }
    } else if (suff < pref) {  // remove everything in between
      for (std::int32_t i = suff; i < pref; ++i) {
        auto it = find_edge(path[i], path[i + 1]);
        dst.emplace(it->id);
        dst.emplace(it->pair->id);
      }
    }
    return dst;  // empty
  }

  void Graph::RemoveEdges(const std::unordered_set<std::uint32_t> &indices, bool remove_nodes) {

    auto erase_remove = [](std::vector<Edge *> &edges, Edge *marked) -> void {
      edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
    };

    std::unordered_set<std::uint32_t> node_indices;
    for (auto i: indices) {
      if (remove_nodes) {
        node_indices.emplace(edges_[i]->tail->id);
        node_indices.emplace(edges_[i]->head->id);
      }
      erase_remove(edges_[i]->tail->outedges, edges_[i].get());
      erase_remove(edges_[i]->head->inedges, edges_[i].get());
    }
    if (remove_nodes) {
      for (auto i: node_indices) {
        if (nodes_[i]->outdegree() == 0 && nodes_[i]->indegree() == 0) {
          nodes_[i].reset();
        }
      }
    }
    for (auto i: indices) {
      edges_[i].reset();
    }
  }

  void Graph::RemoveUnitigEdges(const std::unordered_set<std::uint32_t> &indices, bool remove_nodes) {

    auto erase_remove = [](std::vector<Edge *> &edges, Edge *marked) -> void {
      edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
    };

    std::unordered_set<std::uint32_t> node_indices;
    for (auto i: indices) {
      if (remove_nodes) {
        node_indices.emplace(unitig_edges_[i]->tail->id);
        node_indices.emplace(unitig_edges_[i]->head->id);
      }
      erase_remove(unitig_edges_[i]->tail->outedges, unitig_edges_[i].get());
      erase_remove(unitig_edges_[i]->head->inedges, unitig_edges_[i].get());
    }
    if (remove_nodes) {
      for (auto i: node_indices) {
        if (unitig_nodes_[i]->outdegree() == 0 && unitig_nodes_[i]->indegree() == 0) {
          unitig_nodes_[i].reset();
        }
      }
    }
    for (auto i: indices) {
      unitig_edges_[i].reset();
    }
  }

  void Graph::RemoveDiploidEdges(const std::vector<MarkedEdge> &indices, bool remove_nodes) {
    auto erase_remove = [](std::vector<Edge *> &edges, Edge *marked) -> void {
      edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
    };

    std::unordered_set<std::uint32_t> node_indices;
    std::unordered_set<std::uint32_t> node_indices_alternate;

    for (auto i: indices) {
      Edge *edge = edges_[i.id].get();
      Edge *edge_alternate = edge->alternate;

      if (remove_nodes && (i.where == 0 || i.where == 1)) {
        node_indices.emplace(edge->tail->id);
        node_indices.emplace(edge->head->id);
      }
      if (remove_nodes && (i.where == 0 || i.where == 2)) {
        if (edge_alternate != nullptr) {
          node_indices_alternate.emplace(edge_alternate->tail->id);
          node_indices_alternate.emplace(edge_alternate->head->id);
        }
      }

      if (i.where == 0 || i.where == 1) {
        erase_remove(edges_[i.id]->tail->outedges, edge);
        erase_remove(edges_[i.id]->head->inedges, edge);
      }
      if (i.where == 0 || i.where == 2) {
        if (edge_alternate != nullptr) {
          erase_remove(edge_alternate->tail->outedges, edge_alternate);
          erase_remove(edge_alternate->head->inedges, edge_alternate);
        }
      }
    }

    if (remove_nodes) {
      for (auto i: node_indices) {
        if (nodes_[i]->outdegree() == 0 && nodes_[i]->indegree() == 0) {
          nodes_[i].reset();
        }
      }
      for (auto i: node_indices_alternate) {
        if (i < nodes_alternate_.size() && nodes_alternate_[i] != nullptr &&
            nodes_alternate_[i]->outdegree() == 0 &&
            nodes_alternate_[i]->indegree() == 0) {
          nodes_alternate_[i].reset();
        }
      }
    }

    for (auto i: indices) {
      Edge *edge = edges_[i.id].get();
      Edge *edge_alternate = edge->alternate;

      if (i.where == 0 || i.where == 1) {
        if (edges_.size() > i.id)
          edges_[i.id].reset();
      }
      if (i.where == 0 || i.where == 2) {
        if (edge_alternate != nullptr && edges_alternate_.size() > edge_alternate->id) {
          edges_alternate_[edge_alternate->id].reset();
        }
      }
    }
  }

  void Graph::RemoveAlternateEdges(const std::vector<MarkedEdge> &indices, bool remove_nodes) {
    auto erase_remove = [](std::vector<Edge *> &edges, Edge *marked) -> void {
      edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
    };

    std::unordered_set<std::uint32_t> node_indices_alternate;

    for (auto i: indices) {
      Edge *edge = edges_alternate_[i.id].get();

      if (remove_nodes && (i.where == 0 || i.where == 2)) {
        node_indices_alternate.emplace(edge->tail->id);
        node_indices_alternate.emplace(edge->head->id);
      }

      if (i.where == 0 || i.where == 2) {
        erase_remove(edges_alternate_[i.id]->tail->outedges, edges_alternate_[i.id].get());
        erase_remove(edges_alternate_[i.id]->head->inedges, edges_alternate_[i.id].get());
      }
    }

    if (remove_nodes) {
      for (auto i: node_indices_alternate) {
        if (nodes_alternate_[i]->outdegree() == 0 && nodes_alternate_[i]->indegree() == 0) {
          nodes_alternate_[i].reset();
        }
      }
    }

    for (auto i: indices) {
      edges_alternate_[i.id].reset();
    }
  }

  std::uint32_t Graph::SalvagePlasmids() {
    CreateUnitigs();

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> plasmids;
    for (const auto &it: nodes_) {
      if (!it || it->is_rc() || it->is_unitig || !it->is_circular) {
        continue;
      }
      plasmids.emplace_back(new biosoup::NucleicAcid(it->sequence));
    }
    if (plasmids.empty()) {
      return 0;
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
    for (const auto &it: nodes_) {
      if (!it || it->is_rc() || !it->is_unitig) {
        continue;
      }
      unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
    }
    if (unitigs.empty()) {
      return 0;
    }

    // remove duplicates within plasmids
    std::sort(plasmids.begin(), plasmids.end(),
              [](const std::unique_ptr<biosoup::NucleicAcid> &lhs,
                 const std::unique_ptr<biosoup::NucleicAcid> &rhs) -> bool {
                return lhs->inflated_len < rhs->inflated_len;
              });
    for (std::uint32_t i = 0; i < plasmids.size(); ++i) {
      plasmids[i]->id = i;
    }

    ram::MinimizerEngine minimizer_engine{thread_pool_};
    minimizer_engine.Minimize(plasmids.begin(), plasmids.end());
    minimizer_engine.Filter(0.001);
    for (auto &it: plasmids) {
      if (!minimizer_engine.Map(it, true, true).empty()) {
        it.reset();
      }
    }
    plasmids.erase(
        std::remove(plasmids.begin(), plasmids.end(), nullptr),
        plasmids.end());

    // remove duplicates between plasmids and unitigs
    minimizer_engine.Minimize(unitigs.begin(), unitigs.end(), true);
    minimizer_engine.Filter(0.001);
    for (auto &it: plasmids) {
      if (!minimizer_engine.Map(it, false, false).empty()) {
        it.reset();
      }
    }
    plasmids.erase(
        std::remove(plasmids.begin(), plasmids.end(), nullptr),
        plasmids.end());

    // update nodes
    for (const auto &it: plasmids) {
      const auto &node = nodes_[std::atoi(&it->name[3])];
      node->is_unitig = node->pair->is_unitig = true;
      node->sequence.name[0] = node->pair->sequence.name[0] = 'U';
    }

    return plasmids.size();
  }

  void Graph::SalvageHaplotypes() {
    ram::MinimizerEngine minimizer_engine{thread_pool_};

    while (true) {
      auto num_nodes = nodes_.size();

      CreateUnitigs();
      if (num_nodes == nodes_.size()) {
        break;
      }

      // extend primary
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
      for (const auto &it: nodes_) {
        if (!it || it->is_rc() || !it->is_unitig) {
          continue;
        }
        unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
        unitigs.back()->id = unitigs.size() - 1;
      }
      if (unitigs.empty()) {
        return;
      }

      minimizer_engine.Minimize(unitigs.begin(), unitigs.end());
      minimizer_engine.Filter(0.001);

      std::vector<std::vector<biosoup::Overlap>> overlaps(unitigs.size());
      for (const auto &it: unitigs) {
        for (const auto &jt: minimizer_engine.Map(it, true, true)) {
          std::uint32_t lhs_len = unitigs[jt.lhs_id]->inflated_len;
          std::uint32_t lhs_begin = jt.lhs_begin;
          std::uint32_t lhs_end = jt.lhs_end;

          std::uint32_t rhs_len = unitigs[jt.rhs_id]->inflated_len;
          std::uint32_t rhs_begin = jt.strand ? jt.rhs_begin : rhs_len - jt.rhs_end;  // NOLINT
          std::uint32_t rhs_end = jt.strand ? jt.rhs_end : rhs_len - jt.rhs_begin;

          std::uint32_t overhang =
              std::min(lhs_begin, rhs_begin) +
              std::min(lhs_len - lhs_end, rhs_len - rhs_end);

          if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
              rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
            continue;
          }
          if (lhs_end - lhs_begin < lhs_len * 0.9 &&
              rhs_end - rhs_begin < rhs_len * 0.9) {
            continue;
          }
          if ((lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) ||  // NOLINT
              (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end)) {  // NOLINT
            continue;
          }
          if (lhs_len > rhs_len) {
            overlaps[jt.lhs_id].emplace_back(jt);
          } else {
            overlaps[jt.rhs_id].emplace_back(
                jt.rhs_id, jt.rhs_begin, jt.rhs_end,
                jt.lhs_id, jt.lhs_begin, jt.lhs_end,
                jt.score,
                jt.strand);
          }
        }
      }
      for (auto &it: overlaps) {
        if (it.empty() || it.size() > 2) {
          continue;
        }

        const auto &unitig = unitigs[it.front().lhs_id];
        auto data = unitig->InflateData();

        for (auto &jt: it) {
          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
            auto tmp = jt.rhs_begin;
            jt.rhs_begin = unitigs[jt.rhs_id]->inflated_len - jt.rhs_end;
            jt.rhs_end = unitigs[jt.rhs_id]->inflated_len - tmp;
          }

          if (jt.lhs_begin > jt.rhs_begin) {
            data += unitigs[jt.rhs_id]->InflateData(jt.rhs_end);
          } else {
            data = unitigs[jt.rhs_id]->InflateData(0, jt.rhs_begin) + data;
          }

          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
          }
        }

        const auto &node = nodes_[std::atoi(&unitig->name[3])];

        auto na = biosoup::NucleicAcid(node->sequence.name, data);
        node->sequence = na;

        na.name = node->pair->sequence.name;
        na.ReverseAndComplement();
        node->pair->sequence = na;
      }

      // reconstruct alternative
      overlaps.clear();
      unitigs.clear();
      for (const auto &it: nodes_) {
        if (!it || it->is_rc() || !it->is_unitig) {
          continue;
        }
        unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
        unitigs.back()->id = unitigs.size() - 1;
      }
      overlaps.resize(unitigs.size());

      minimizer_engine.Minimize(unitigs.begin(), unitigs.end());
      minimizer_engine.Filter(0.001);

      for (const auto &it: unitigs) {
        for (const auto &jt: minimizer_engine.Map(it, true, true)) {
          std::uint32_t lhs_len = unitigs[jt.lhs_id]->inflated_len;
          std::uint32_t lhs_begin = jt.lhs_begin;
          std::uint32_t lhs_end = jt.lhs_end;

          std::uint32_t rhs_len = unitigs[jt.rhs_id]->inflated_len;
          std::uint32_t rhs_begin = jt.strand ? jt.rhs_begin : rhs_len - jt.rhs_end;  // NOLINT
          std::uint32_t rhs_end = jt.strand ? jt.rhs_end : rhs_len - jt.rhs_begin;

          std::uint32_t overhang =
              std::min(lhs_begin, rhs_begin) +
              std::min(lhs_len - lhs_end, rhs_len - rhs_end);

          if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
              rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
            continue;
          }
          if (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end) {
            overlaps[jt.lhs_id].emplace_back(jt);
          } else if (lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) {  // NOLINT
            overlaps[jt.rhs_id].emplace_back(
                jt.rhs_id, jt.rhs_begin, jt.rhs_end,
                jt.lhs_id, jt.lhs_begin, jt.lhs_end,
                jt.score,
                jt.strand);
          }
        }
      }
      for (auto &it: overlaps) {
        if (it.empty()) {
          continue;
        }

        const auto &unitig = unitigs[it.front().lhs_id];
        if (nodes_[std::atoi(&unitig->name[3])] == nullptr) {
          continue;
        }

        std::sort(it.begin(), it.end(),
                  [](const biosoup::Overlap &lhs,
                     const biosoup::Overlap &rhs) -> bool {
                    return lhs.lhs_begin < rhs.lhs_begin;
                  });
        it.emplace_back(unitig->id, -1, -1, unitig->id, -1, -1, 0);  // dummy

        std::string data = unitig->InflateData(0, it.front().lhs_begin);

        std::unordered_set<std::uint32_t> marked_edges;
        for (std::uint32_t j = 0; j < it.size() - 1; ++j) {
          auto &jt = it[j];

          const auto &n = nodes_[std::atoi(&unitigs[jt.rhs_id]->name[3])];
          if (n == nullptr) {
            continue;
          }

          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
            auto tmp = jt.rhs_begin;
            jt.rhs_begin = unitigs[jt.rhs_id]->inflated_len - jt.rhs_end;
            jt.rhs_end = unitigs[jt.rhs_id]->inflated_len - tmp;
          }

          data += unitigs[jt.rhs_id]->InflateData(
              jt.rhs_begin,
              jt.rhs_end - jt.rhs_begin);
          if (jt.lhs_end < it[j + 1].lhs_begin) {
            data += unitig->InflateData(
                jt.lhs_end,
                it[j + 1].lhs_begin - jt.lhs_end);
          }

          for (const auto &kt: n->inedges) {
            marked_edges.emplace(kt->id);
            marked_edges.emplace(kt->pair->id);
          }
          for (const auto &kt: n->outedges) {
            marked_edges.emplace(kt->id);
            marked_edges.emplace(kt->pair->id);
          }

          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
          }
        }
        RemoveEdges(marked_edges);

        it.pop_back();
        for (const auto &jt: it) {
          auto id = std::atoi(&unitigs[jt.rhs_id]->name[3]);
          nodes_[id ^ 1].reset();
          nodes_[id].reset();
        }

        auto na = biosoup::NucleicAcid("", data);

        auto node = std::make_shared<Node>(na);
        node->sequence.name = "Utg" + std::to_string(node->id & (~1UL));
        node->is_unitig = true;
        nodes_.emplace_back(node);

        na.ReverseAndComplement();
        nodes_.emplace_back(std::make_shared<Node>(na));

        node->pair = nodes_.back().get();
        node->pair->pair = node.get();
        node->pair->sequence.name = node->sequence.name;
        node->pair->is_unitig = true;
      }
    }
  }

  void Graph::SalvageHaplotypesPrimary() {

    std::cerr << "[raven::Graph::Assemble] SalvageHaplotypesPrimary "
              << std::endl;

    ram::MinimizerEngine minimizer_engine{thread_pool_};

    while (true) {
      auto num_nodes = nodes_.size();

      CreateUnitigs();
      if (num_nodes == nodes_.size()) {
        break;
      }

      // extend primary
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
      for (const auto &it: nodes_) {
        if (!it || it->is_rc() || !it->is_unitig) {
          continue;
        }
        unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
        unitigs.back()->id = unitigs.size() - 1;
      }
      if (unitigs.empty()) {
        return;
      }

      minimizer_engine.Minimize(unitigs.begin(), unitigs.end());
      minimizer_engine.Filter(0.001);

      std::vector<std::vector<biosoup::Overlap>> overlaps(unitigs.size());

      for (const auto &it: unitigs) {
        for (const auto &jt: minimizer_engine.Map(it, true, true)) {
          std::uint32_t lhs_len = unitigs[jt.lhs_id]->inflated_len;
          std::uint32_t lhs_begin = jt.lhs_begin;
          std::uint32_t lhs_end = jt.lhs_end;

          std::uint32_t rhs_len = unitigs[jt.rhs_id]->inflated_len;
          std::uint32_t rhs_begin = jt.strand ? jt.rhs_begin : rhs_len - jt.rhs_end;  // NOLINT
          std::uint32_t rhs_end = jt.strand ? jt.rhs_end : rhs_len - jt.rhs_begin;

          std::uint32_t overhang =
              std::min(lhs_begin, rhs_begin) +
              std::min(lhs_len - lhs_end, rhs_len - rhs_end);

          if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
              rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
            continue;
          }
          if (lhs_end - lhs_begin < lhs_len * 0.9 &&
              rhs_end - rhs_begin < rhs_len * 0.9) {
            continue;
          }
          if ((lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) ||  // NOLINT
              (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end)) {  // NOLINT
            continue;
          }
          if (lhs_len > rhs_len) {
            overlaps[jt.lhs_id].emplace_back(jt);
          } else {
            overlaps[jt.rhs_id].emplace_back(
                jt.rhs_id, jt.rhs_begin, jt.rhs_end,
                jt.lhs_id, jt.lhs_begin, jt.lhs_end,
                jt.score,
                jt.strand);
          }
        }
      }
      for (auto &it: overlaps) {
        if (it.empty() || it.size() > 2) {
          continue;
        }

        const auto &unitig = unitigs[it.front().lhs_id];
        auto data = unitig->InflateData();

        for (auto &jt: it) {
          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
            auto tmp = jt.rhs_begin;
            jt.rhs_begin = unitigs[jt.rhs_id]->inflated_len - jt.rhs_end;
            jt.rhs_end = unitigs[jt.rhs_id]->inflated_len - tmp;
          }

          if (jt.lhs_begin > jt.rhs_begin) {
            data += unitigs[jt.rhs_id]->InflateData(jt.rhs_end);
          } else {
            data = unitigs[jt.rhs_id]->InflateData(0, jt.rhs_begin) + data;
          }

          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
          }
        }

        int ind = std::atoi(&unitig->name[3]);
        const auto &node = nodes_[ind];

        auto na = biosoup::NucleicAcid(node->sequence.name, data);
        node->sequence = na;

        na.name = node->pair->sequence.name;
        na.ReverseAndComplement();
        node->pair->sequence = na;
      }
    }
  }

  void Graph::SalvageHaplotypesAlternative() {

    std::cerr << "[raven::Graph::Assemble] SalvageHaplotypesAlternative "
              << std::endl;

    ram::MinimizerEngine minimizer_engine{thread_pool_};

    while (true) {
      auto num_nodes = nodes_alternate_.size();

      CreateUnitigsAlternate();
      if (num_nodes == nodes_alternate_.size()) {
        break;
      }

      // reconstruct alternative
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> unitigs;
      for (const auto &it: nodes_alternate_) {
        if (!it || it->is_rc() || !it->is_unitig) {
          continue;
        }
        unitigs.emplace_back(new biosoup::NucleicAcid(it->sequence));
        unitigs.back()->id = unitigs.size() - 1;
      }
      if (unitigs.empty()) {
        return;
      }

      minimizer_engine.Minimize(unitigs.begin(), unitigs.end());
      minimizer_engine.Filter(0.001);

      std::vector<std::vector<biosoup::Overlap>> overlaps(unitigs.size());

      for (const auto &it: unitigs) {
        for (const auto &jt: minimizer_engine.Map(it, true, true)) {
          std::uint32_t lhs_len = unitigs[jt.lhs_id]->inflated_len;
          std::uint32_t lhs_begin = jt.lhs_begin;
          std::uint32_t lhs_end = jt.lhs_end;

          std::uint32_t rhs_len = unitigs[jt.rhs_id]->inflated_len;
          std::uint32_t rhs_begin = jt.strand ? jt.rhs_begin : rhs_len - jt.rhs_end;  // NOLINT
          std::uint32_t rhs_end = jt.strand ? jt.rhs_end : rhs_len - jt.rhs_begin;

          std::uint32_t overhang =
              std::min(lhs_begin, rhs_begin) +
              std::min(lhs_len - lhs_end, rhs_len - rhs_end);

          if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
              rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
            continue;
          }
          if (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end) {
            overlaps[jt.lhs_id].emplace_back(jt);
          } else if (lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) {  // NOLINT
            overlaps[jt.rhs_id].emplace_back(
                jt.rhs_id, jt.rhs_begin, jt.rhs_end,
                jt.lhs_id, jt.lhs_begin, jt.lhs_end,
                jt.score,
                jt.strand);
          }
        }
      }
      for (auto &it: overlaps) {
        if (it.empty()) {
          continue;
        }

        const auto &unitig = unitigs[it.front().lhs_id];
        if (nodes_alternate_[std::atoi(&unitig->name[3])] == nullptr) {
          continue;
        }

        std::sort(it.begin(), it.end(),
                  [](const biosoup::Overlap &lhs,
                     const biosoup::Overlap &rhs) -> bool {
                    return lhs.lhs_begin < rhs.lhs_begin;
                  });
        it.emplace_back(unitig->id, -1, -1, unitig->id, -1, -1, 0);  // dummy

        std::string data = unitig->InflateData(0, it.front().lhs_begin);

        //std::unordered_set<std::uint32_t> marked_edges;
        std::vector<MarkedEdge> marked_edges;
        for (std::uint32_t j = 0; j < it.size() - 1; ++j) {
          auto &jt = it[j];

          const auto &n = nodes_alternate_[std::atoi(&unitigs[jt.rhs_id]->name[3])];
          if (n == nullptr) {
            continue;
          }

          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
            auto tmp = jt.rhs_begin;
            jt.rhs_begin = unitigs[jt.rhs_id]->inflated_len - jt.rhs_end;
            jt.rhs_end = unitigs[jt.rhs_id]->inflated_len - tmp;
          }

          data += unitigs[jt.rhs_id]->InflateData(
              jt.rhs_begin,
              jt.rhs_end - jt.rhs_begin);
          if (jt.lhs_end < it[j + 1].lhs_begin) {
            data += unitig->InflateData(
                jt.lhs_end,
                it[j + 1].lhs_begin - jt.lhs_end);
          }

          for (const auto &kt: n->inedges) {
            marked_edges.emplace_back(kt->id, 2);
            marked_edges.emplace_back(kt->pair->id, 2);
          }
          for (const auto &kt: n->outedges) {
            marked_edges.emplace_back(kt->id, 2);
            marked_edges.emplace_back(kt->pair->id, 2);
          }

          if (!jt.strand) {
            unitigs[jt.rhs_id]->ReverseAndComplement();
          }
        }
        RemoveAlternateEdges(marked_edges);

        it.pop_back();
        for (const auto &jt: it) {
          auto id = std::atoi(&unitigs[jt.rhs_id]->name[3]);
          nodes_alternate_[id ^ 1].reset();
          nodes_alternate_[id].reset();
        }

        auto na = biosoup::NucleicAcid("", data);

        auto node = std::make_shared<Node>(na, Node::num_objects_alternate++);
        node->sequence.name = "Utg" + std::to_string(node->id & (~1UL));
        node->is_unitig = true;
        nodes_alternate_.emplace_back(node);

        na.ReverseAndComplement();
        nodes_alternate_.emplace_back(std::make_shared<Node>(na, Node::num_objects_alternate++));

        node->pair = nodes_alternate_.back().get();
        node->pair->pair = node.get();
        node->pair->sequence.name = node->sequence.name;
        node->pair->is_unitig = true;
      }
    }
  }

}  // namespace raven
