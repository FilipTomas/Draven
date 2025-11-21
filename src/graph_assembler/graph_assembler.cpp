
#include <fstream>
#include <random>
#include "graph_assembler.h"
#include "biosoup/timer.hpp"
#include "edlib.h"

namespace raven {
  Graph_Assembler::Graph_Assembler(Graph &graph, Program_Parameters &param, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
      : graph_(graph), param_(param), thread_pool_(thread_pool ?
                                    thread_pool :
                                    std::make_shared<thread_pool::ThreadPool>(1)) {
  }

  void Graph_Assembler::Assemble() {
    if (!graph_.state_manager_.assemble_any_state()) {
      return;
    }

    biosoup::Timer timer{};

    // remove transitive edges
    if (graph_.state_manager_.state() == GraphState::Assemble_Transitive_Edges) {
      timer.Start();

      RemoveTransitiveEdges();

      std::cerr << "[raven::Graph::Assemble] removed transitive edges "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_transitive_filename, false);
      //exit(0);
    }

    if (graph_.state_manager_.state() == GraphState::Assemble_Transitive_Edges) {
      timer.Start();

      //CreateUnitigGraph();

      std::cerr << "[raven::Graph::Assemble] created bubble chain "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      // PrintGfa("after_bubble_chain.gfa");
    }

    // checkpoint
    if (graph_.state_manager_.state() == GraphState::Assemble_Transitive_Edges) {
      graph_.state_manager_.advance_state();
      if (graph_.use_checkpoints()) {
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    // remove tips and bubbles
    if (graph_.state_manager_.state() == GraphState::Assemble_Tips_Bubbles) {
      timer.Start();

      while (true) {
        std::uint32_t num_changes = RemoveTips();
        num_changes += RemoveBubbles();
        if (num_changes == 0) {
          break;
        }
      }

      std::cerr << "[raven::Graph::Assemble] removed tips and bubbles "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_bubble_filename, false);
    }

    // checkpoint
    if (graph_.state_manager_.state() == GraphState::Assemble_Tips_Bubbles) {
      graph_.state_manager_.advance_state();
      if (graph_.use_checkpoints()) {
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    // remove long edges
    if (graph_.state_manager_.state() == GraphState::Assemble_Force_Directed) {
      timer.Start();

      if (graph_.annotations_.empty()) {
        graph_.CreateUnitigs(42);  // speed up force directed layout
      }
      RemoveLongEdges(16);

      std::cerr << "[raven::Graph::Assemble] removed long edges "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_force_filename, false);
    }

    // checkpoint
    if (graph_.state_manager_.state() == GraphState::Assemble_Force_Directed) {
      graph_.state_manager_.advance_state();
      if (graph_.use_checkpoints()) {
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    std::cerr << "[raven::Graph::Assemble] "
              << std::fixed << timer.elapsed_time() << "s"
              << std::endl;
  }

  void Graph_Assembler::AssembleDiploids() {
    if (graph_.state_manager_.state() != GraphState::Assemble_Diploids) {
      return;
    }
    biosoup::Timer timer{};

    timer.Start();

    graph_.DuplicateGraph();

    std::cerr << "[raven::Graph::AssembleDiploids] Duplicated graph "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    while (true) {
      std::uint32_t num_changes = RemoveSnpBubbles();

      if (num_changes == 0) {
        break;
      }
    }

    std::cerr << "[raven::Graph::AssembleDiploids] removed SNP bubbles "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    //graph_.SalvageHaplotypesPrimary();
    //graph_.SalvageHaplotypesAlternative();
    graph_.SalvageHaplotypes();

    std::cerr << "[raven::Graph::AssembleDiploids] Created haplotypes "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    std::cerr << "[raven::Graph::AssembleDiploids] "
              << std::fixed << timer.elapsed_time() << "s"
              << std::endl;
  }

  void Graph_Assembler::AssembleHaploids(){
      if (graph_.state_manager_.state() != GraphState::Assemble_Diploids) {
      return;
    }
    biosoup::Timer timer{};

    timer.Start();

    graph_.SalvageHaplotypes();

    std::cerr << "[raven::Graph::AssembleHaploids] Created haplotypes "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  void Graph_Assembler::UlAssemble(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_sequences) {
    if (!graph_.state_manager_.assemble_any_state()) {
      return;
    }

    biosoup::Timer timer{};
    graph_.PrintGfa(param_.gfa_after_construction_filename, false);

    // remove transitive edges
    if (graph_.state_manager_.state() == GraphState::Assemble_Transitive_Edges) {
      timer.Start();

      RemoveTransitiveEdges();

      std::cerr << "[raven::Graph::Assemble] removed transitive edges "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_transitive_filename, false);
    }

    if (graph_.state_manager_.state() == GraphState::Assemble_Transitive_Edges) {
      timer.Start();

      graph_.CreateUnitigGraph(param_.gfa_unitig_graph_filename);

      std::cerr << "[raven::Graph::Assemble] created bubble chain "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      // PrintGfa("after_bubble_chain.gfa");
    }

    if (graph_.state_manager_.state() == GraphState::Assemble_Transitive_Edges) {
      timer.Start();

      ResolveGraphWithUl(ul_sequences);

      std::cerr << "[raven::Graph::Assemble] created bubble chain "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

    }

  }

  // std::ofstream transitive_pairs("transitive_pairs.txt");
  
  std::uint32_t Graph_Assembler::RemoveTransitiveEdges() {
    auto is_comparable = [](double a, double b) -> bool {
      double eps = 0.12;
      return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
             (b >= a * (1 - eps) && b <= a * (1 + eps));
    };

    std::vector<Edge *> candidate(graph_.nodes_.size(), nullptr);
    std::unordered_set<std::uint32_t> marked_edges;
    for (const auto &it: graph_.nodes_) {
      if (it == nullptr) {
        continue;
      }
      for (auto jt: it->outedges) {
        candidate[jt->head->id] = jt;
      }
      for (auto jt: it->outedges) {
        for (auto kt: jt->head->outedges) {
          if (candidate[kt->head->id] &&
              is_comparable(jt->length + kt->length, candidate[kt->head->id]->length)) {  // NOLINT
            marked_edges.emplace(candidate[kt->head->id]->id);
            marked_edges.emplace(candidate[kt->head->id]->pair->id);
            // transitive_pairs << candidate[kt->head->id]->tail->sequence.name << "->" << candidate[kt->head->id]->head->sequence.name 
            //                  << "\t" << graph_.edges_[jt->id]->tail->sequence.name << "->" << graph_.edges_[jt->id]->head->sequence.name 
            //                  << "\t" << graph_.edges_[kt->id]->tail->sequence.name << "->" << graph_.edges_[kt->id]->head->sequence.name << std::endl;
          }
        }
      }
      for (auto jt: it->outedges) {
        candidate[jt->head->id] = nullptr;
      }
    }

    for (auto i: marked_edges) {  // store for force directed layout
      if (i & 1) {
        auto lhs = graph_.edges_[i]->tail->id & ~1UL;
        auto rhs = graph_.edges_[i]->head->id & ~1UL;
        graph_.nodes_[lhs]->transitive.emplace(rhs);
        graph_.nodes_[rhs]->transitive.emplace(lhs);
      }
    }

    graph_.RemoveEdges(marked_edges);
    return marked_edges.size() / 2;
  }

  std::uint32_t Graph_Assembler::RemoveTips() {
    std::uint32_t num_tips = 0;
    std::vector<char> is_visited(graph_.nodes_.size(), 0);

    for (const auto &it: graph_.nodes_) {
      if (it == nullptr || is_visited[it->id] || !it->is_tip()) {
        continue;
      }
      bool is_circular = false;
      std::uint32_t num_sequences = 0;

      auto end = it.get();
      while (!end->is_junction()) {
        num_sequences += end->count;
        is_visited[end->id] = 1;
        is_visited[end->pair->id] = 1;
        if (end->outdegree() == 0 ||
            end->outedges.front()->head->is_junction()) {
          break;
        }
        end = end->outedges.front()->head;
        if (end == it.get()) {
          is_circular = true;
          break;
        }
      }

      if (is_circular || end->outdegree() == 0 || num_sequences > 5) {
        continue;
      }

      std::unordered_set<std::uint32_t> marked_edges;
      for (auto jt: end->outedges) {
        if (jt->head->indegree() > 1) {
          marked_edges.emplace(jt->id);
          marked_edges.emplace(jt->pair->id);
        }
      }
      if (marked_edges.size() / 2 == end->outedges.size()) {  // delete whole
        auto begin = it.get();
        while (begin != end) {
          marked_edges.emplace(begin->outedges.front()->id);
          marked_edges.emplace(begin->outedges.front()->pair->id);
          begin = begin->outedges.front()->head;
        }
        ++num_tips;
      }
      graph_.RemoveEdges(marked_edges, true);
    }

    return num_tips;
  }

  std::uint32_t Graph_Assembler::RemoveBubbles() {
    std::vector<std::uint32_t> distance(graph_.nodes_.size(), 0);
    std::vector<Node *> predecessor(graph_.nodes_.size(), nullptr);

    //region annotations_ helper functions
    auto annotation_extract = [&](
        std::uint32_t i,
        std::uint32_t begin,
        std::uint32_t end,
        std::uint32_t len,
        bool strand) -> std::unordered_set<std::uint32_t> {
      std::unordered_set<std::uint32_t> dst;
      if (graph_.annotations_[i].empty()) {
        return dst;
      }
      for (const auto &it: graph_.annotations_[i]) {
        if (begin <= it && it <= end) {
          dst.emplace(strand ? it : len - 1 - it);
        }
      }
      return dst;
    };
    //endregion annotations_ helper functions

    //region path helper functions
    auto path_extract = [&](Node *begin, Node *end) -> std::vector<Node *> {
      std::vector<Node *> dst;
      while (end != begin) {
        dst.emplace_back(end);
        end = predecessor[end->id];
      }
      dst.emplace_back(begin);
      std::reverse(dst.begin(), dst.end());
      return dst;
    };
    auto path_is_simple = [](const std::vector<Node *> &path) -> bool {
      if (path.empty()) {
        return false;
      }
      for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
        if (path[i]->is_junction()) {
          return false;  // complex
        }
      }
      return true;  // without branches
    };
    auto path_sequence = [](const std::vector<Node *> &path) -> std::string {
      std::string data{};
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        for (auto it: path[i]->outedges) {
          if (it->head == path[i + 1]) {
            data += it->Label();
            break;
          }
        }
      }
      data += path.back()->sequence.InflateData();
      return data;
    };
    auto path_annotation = [&](const std::vector<Node *> &path)
        -> std::unordered_set<std::uint32_t> {
      std::unordered_set<std::uint32_t> dst;
      std::uint32_t offset = 0;
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        for (auto it: path[i]->outedges) {
          if (it->head == path[i + 1]) {
            const auto &annotations = annotation_extract(
                it->tail->sequence.id,
                0,
                it->length,
                it->tail->sequence.inflated_len,
                !it->tail->is_rc());
            for (const auto &jt: annotations) {
              dst.emplace(offset + jt);
            }
            offset += it->length;
          }
        }
      }
      const auto &annotations = annotation_extract(
          path.back()->sequence.id,
          0,
          path.back()->sequence.inflated_len,
          path.back()->sequence.inflated_len,
          !path.back()->is_rc());
      for (const auto &jt: annotations) {
        dst.emplace(offset + jt);
      }
      return dst;
    };
    auto bubble_pop = [&](
        const std::vector<Node *> &lhs,
        const std::vector<Node *> &rhs) -> std::unordered_set<std::uint32_t> {
      if (lhs.empty() || rhs.empty()) {
        return std::unordered_set<std::uint32_t>{};
      }

      // check BFS result
      std::unordered_set<Node *> bubble;
      bubble.insert(lhs.begin(), lhs.end());
      bubble.insert(rhs.begin(), rhs.end());
      if (lhs.size() + rhs.size() - 2 != bubble.size()) {
        return std::unordered_set<std::uint32_t>{};
      }
      for (const auto &it: lhs) {
        if (bubble.count(it->pair) != 0) {
          return std::unordered_set<std::uint32_t>{};
        }
      }

      if (!path_is_simple(lhs) || !path_is_simple(rhs)) {  // complex path(s)
        // check poppability
        if (graph_.FindRemovableEdges(lhs).empty() && graph_.FindRemovableEdges(rhs).empty()) {
          return std::unordered_set<std::uint32_t>{};
        }

        // check similarity
        auto l = path_sequence(lhs);
        auto r = path_sequence(rhs);
        if (std::min(l.size(), r.size()) < std::max(l.size(), r.size()) * 0.8) {
          return std::unordered_set<std::uint32_t>{};
        }

        EdlibAlignResult result = edlibAlign(
            l.c_str(), l.size(),
            r.c_str(), r.size(),
            edlibDefaultAlignConfig());
        double score = 0;
        if (result.status == EDLIB_STATUS_OK) {
          score = 1 - result.editDistance /
                      static_cast<double>(std::max(l.size(), r.size()));
          edlibFreeAlignResult(result);
        }
        if (score < 0.8) {
          return std::unordered_set<std::uint32_t>{};
        }
        if (!graph_.annotations_.empty()) {
          std::unordered_set<std::uint32_t> marked_edges;
          if (!path_is_simple(lhs)) {
            marked_edges = graph_.FindRemovableEdges(lhs);
            if (marked_edges.size() > 2) {
              marked_edges.clear();
            }
          }
          if (marked_edges.empty() && !path_is_simple(rhs)) {
            marked_edges = graph_.FindRemovableEdges(rhs);
            if (marked_edges.size() > 2) {
              marked_edges.clear();
            }
          }
          return marked_edges;
        }
      }

      if (!graph_.annotations_.empty()) {
        auto la = path_annotation(lhs);
        auto ra = path_annotation(rhs);

        if (!la.empty() && !ra.empty()) {
          auto l = path_sequence(lhs);
          auto r = path_sequence(rhs);

          EdlibAlignResult result = edlibAlign(
              l.c_str(), l.size(),
              r.c_str(), r.size(),
              edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

          if (result.status == EDLIB_STATUS_OK) {
            std::uint32_t lhs_pos = 0;
            std::uint32_t rhs_pos = 0;

            std::uint32_t mismatches = 0;
            std::uint32_t snps = 0;

            for (int a = 0; a < result.alignmentLength; ++a) {
              if (la.find(lhs_pos) != la.end() ||
                  ra.find(rhs_pos) != ra.end()) {
                ++snps;
                if (result.alignment[a] == 3) {
                  ++mismatches;
                }
              }
              switch (result.alignment[a]) {
                case 0:
                case 3: {
                  ++lhs_pos;
                  ++rhs_pos;
                  break;
                }
                case 1: {
                  ++lhs_pos;
                  break;
                }
                case 2: {
                  ++rhs_pos;
                  break;
                }
                default:
                  break;
              }
            }

            edlibFreeAlignResult(result);

            if (mismatches / static_cast<double>(snps) > 0.1) {  // disagreement_) {
              return std::unordered_set<std::uint32_t>{};
            }
          }
        }
      }

      std::uint32_t lhs_count = 0;
      for (auto jt: lhs) {
        lhs_count += jt->count;
      }
      std::uint32_t rhs_count = 0;
      for (auto jt: rhs) {
        rhs_count += jt->count;
      }
      auto marked_edges = graph_.FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
      if (marked_edges.empty()) {
        marked_edges = graph_.FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
      }
      return marked_edges;
    };
    //endregion path helper functions

    std::uint32_t num_bubbles = 0;
    for (const auto &it: graph_.nodes_) {
      if (it == nullptr || it->outdegree() < 2) {
        continue;
      }

      // BFS
      Node *begin = it.get();
      Node *end = nullptr;
      Node *other_end = nullptr;
      std::deque<Node *> que{begin};
      std::vector<Node *> visited{1, begin};
      while (!que.empty() && !end) {
        auto jt = que.front();
        que.pop_front();

        for (auto kt: jt->outedges) {
          if (kt->head == begin) {  // cycle
            continue;
          }
          if (distance[jt->id] + kt->length > 5000000) {  // out of reach
            continue;
          }
          if (kt->head->id > distance.size())
            continue;
          distance[kt->head->id] = distance[jt->id] + kt->length;
          visited.emplace_back(kt->head);
          que.emplace_back(kt->head);

          if (predecessor[kt->head->id]) {  // found bubble
            end = kt->head;
            other_end = jt;
            break;
          }

          predecessor[kt->head->id] = jt;
        }
      }

      std::unordered_set<std::uint32_t> marked_edges;
      if (end) {
        auto lhs = path_extract(begin, end);
        auto rhs = path_extract(begin, other_end);
        rhs.emplace_back(end);
        marked_edges = bubble_pop(lhs, rhs);
      }

      for (auto jt: visited) {
        distance[jt->id] = 0;
        predecessor[jt->id] = nullptr;
      }

      graph_.RemoveEdges(marked_edges, true);
      num_bubbles += 1 - marked_edges.empty();
    }

    return num_bubbles;
  }

  std::uint32_t Graph_Assembler::RemoveSnpBubbles() {
    std::vector<std::uint32_t> distance(graph_.nodes_.size(), 0);
    std::vector<Node *> predecessor(graph_.nodes_.size(), nullptr);

    //region annotations_ helper functions
    auto annotation_extract = [&](
        std::uint32_t i,
        std::uint32_t begin,
        std::uint32_t end,
        std::uint32_t len,
        bool strand) -> std::unordered_set<std::uint32_t> {
      std::unordered_set<std::uint32_t> dst;
      if (graph_.annotations_[i].empty()) {
        return dst;
      }
      for (const auto &it: graph_.annotations_[i]) {
        if (begin <= it && it <= end) {
          dst.emplace(strand ? it : len - 1 - it);
        }
      }
      return dst;
    };
    //endregion annotations_ helper functions

    //region path helper functions
    auto path_extract = [&](Node *begin, Node *end) -> std::vector<Node *> {
      std::vector<Node *> dst;
      while (end != begin) {
        dst.emplace_back(end);
        end = predecessor[end->id];
      }
      dst.emplace_back(begin);
      std::reverse(dst.begin(), dst.end());
      return dst;
    };
    auto path_is_simple = [](const std::vector<Node *> &path) -> bool {
      if (path.empty()) {
        return false;
      }
      for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
        if (path[i]->is_junction()) {
          return false;  // complex
        }
      }
      return true;  // without branches
    };
    auto path_sequence = [](const std::vector<Node *> &path) -> std::string {
      std::string data{};
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        for (auto it: path[i]->outedges) {
          if (it->head == path[i + 1]) {
            data += it->Label();
            break;
          }
        }
      }
      data += path.back()->sequence.InflateData();
      return data;
    };
    auto path_annotation = [&](const std::vector<Node *> &path)
        -> std::unordered_set<std::uint32_t> {
      std::unordered_set<std::uint32_t> dst;
      std::uint32_t offset = 0;
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        for (auto it: path[i]->outedges) {
          if (it->head == path[i + 1]) {
            const auto &annotations = annotation_extract(
                it->tail->sequence.id,
                0,
                it->length,
                it->tail->sequence.inflated_len,
                !it->tail->is_rc());
            for (const auto &jt: annotations) {
              dst.emplace(offset + jt);
            }
            offset += it->length;
          }
        }
      }
      const auto &annotations = annotation_extract(
          path.back()->sequence.id,
          0,
          path.back()->sequence.inflated_len,
          path.back()->sequence.inflated_len,
          !path.back()->is_rc());
      for (const auto &jt: annotations) {
        dst.emplace(offset + jt);
      }
      return dst;
    };
    auto bubble_pop = [&](
        const std::vector<Node *> &lhs,
        const std::vector<Node *> &rhs) -> std::vector<MarkedEdge> {
      if (lhs.empty() || rhs.empty()) {
        return std::vector<MarkedEdge>{};
      }

      // check BFS result
      std::unordered_set<Node *> bubble;
      bubble.insert(lhs.begin(), lhs.end());
      bubble.insert(rhs.begin(), rhs.end());
      if (lhs.size() + rhs.size() - 2 != bubble.size()) {
        return std::vector<MarkedEdge>{};
      }
      for (const auto &it: lhs) {
        if (bubble.count(it->pair) != 0) {
          return std::vector<MarkedEdge>{};
        }
      }

      if (!path_is_simple(lhs) || !path_is_simple(rhs)) {  // complex path(s)
        // check poppability
        if (graph_.FindRemovableEdges(lhs).empty() && graph_.FindRemovableEdges(rhs).empty()) {
          return std::vector<MarkedEdge>{};
        }

        // check similarity
        auto l = path_sequence(lhs);
        auto r = path_sequence(rhs);
        if (std::min(l.size(), r.size()) < std::max(l.size(), r.size()) * 0.8) {
          return std::vector<MarkedEdge>{};
        }

        EdlibAlignResult result = edlibAlign(
            l.c_str(), l.size(),
            r.c_str(), r.size(),
            edlibDefaultAlignConfig());
        double score = 0;
        if (result.status == EDLIB_STATUS_OK) {
          score = 1 - result.editDistance /
                      static_cast<double>(std::max(l.size(), r.size()));
          edlibFreeAlignResult(result);
        }
        if (score < 0.8) {
          return std::vector<MarkedEdge>{};
        }
        if (!graph_.annotations_.empty()) {
          std::unordered_set<std::uint32_t> marked_edges;
          if (!path_is_simple(lhs)) {
            marked_edges = graph_.FindRemovableEdges(lhs);
            if (marked_edges.size() > 2) {
              marked_edges.clear();
            }
          }
          if (marked_edges.empty() && !path_is_simple(rhs)) {
            marked_edges = graph_.FindRemovableEdges(rhs);
            if (marked_edges.size() > 2) {
              marked_edges.clear();
            }
          }
          std::vector<MarkedEdge> marked_edge_list;
          for (const auto &it: marked_edges)
            marked_edge_list.emplace_back(it);
          return marked_edge_list;
        }
      }

      if (!graph_.annotations_.empty()) {
        auto la = path_annotation(lhs);
        auto ra = path_annotation(rhs);

        if (!la.empty() && !ra.empty()) {
          auto l = path_sequence(lhs);
          auto r = path_sequence(rhs);

          EdlibAlignResult result = edlibAlign(
              l.c_str(), l.size(),
              r.c_str(), r.size(),
              edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

          if (result.status == EDLIB_STATUS_OK) {
            std::uint32_t lhs_pos = 0;
            std::uint32_t rhs_pos = 0;

            std::uint32_t mismatches = 0;
            std::uint32_t snps = 0;

            for (int a = 0; a < result.alignmentLength; ++a) {
              if (la.find(lhs_pos) != la.end() ||
                  ra.find(rhs_pos) != ra.end()) {
                ++snps;
                if (result.alignment[a] == 3) {
                  ++mismatches;
                }
              }
              switch (result.alignment[a]) {
                case 0:
                case 3: {
                  ++lhs_pos;
                  ++rhs_pos;
                  break;
                }
                case 1: {
                  ++lhs_pos;
                  break;
                }
                case 2: {
                  ++rhs_pos;
                  break;
                }
                default:
                  break;
              }
            }

            edlibFreeAlignResult(result);

            if (mismatches / static_cast<double>(snps) > 0.1) {  // disagreement_) {
              std::vector<MarkedEdge> marked_edge_list;

              auto marked_edges = graph_.FindRemovableEdges(lhs);
              for (const auto &it: marked_edges)
                marked_edge_list.emplace_back(it, 1);

              marked_edges = graph_.FindRemovableEdges(rhs);
              for (const auto &it: marked_edges)
                marked_edge_list.emplace_back(it, 2);

              return marked_edge_list;
            }
          }
        }
      }

      std::uint32_t lhs_count = 0;
      for (auto jt: lhs) {
        lhs_count += jt->count;
      }
      std::uint32_t rhs_count = 0;
      for (auto jt: rhs) {
        rhs_count += jt->count;
      }
      auto marked_edges = graph_.FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
      if (marked_edges.empty()) {
        marked_edges = graph_.FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
      }
      std::vector<MarkedEdge> marked_edge_list;
      for (const auto &it: marked_edges)
        marked_edge_list.emplace_back(it);
      return marked_edge_list;
    };
    //endregion path helper functions

    std::uint32_t num_bubbles = 0;
    for (const auto &it: graph_.nodes_) {
      if (it == nullptr || it->outdegree() < 2) {
        continue;
      }

      // BFS
      Node *begin = it.get();
      Node *end = nullptr;
      Node *other_end = nullptr;
      std::deque<Node *> que{begin};
      std::vector<Node *> visited{1, begin};
      while (!que.empty() && !end) {
        auto jt = que.front();
        que.pop_front();

        for (auto kt: jt->outedges) {
          if (kt->head == begin) {  // cycle
            continue;
          }
          if (distance[jt->id] + kt->length > 5000000) {  // out of reach
            continue;
          }
          if (kt->head->id > distance.size())
            continue;
          distance[kt->head->id] = distance[jt->id] + kt->length;
          visited.emplace_back(kt->head);
          que.emplace_back(kt->head);

          if (predecessor[kt->head->id]) {  // found bubble
            end = kt->head;
            other_end = jt;
            break;
          }

          predecessor[kt->head->id] = jt;
        }
      }

      std::vector<MarkedEdge> marked_edges;
      if (end) {
        auto lhs = path_extract(begin, end);
        auto rhs = path_extract(begin, other_end);
        rhs.emplace_back(end);
        marked_edges = bubble_pop(lhs, rhs);
      }

      for (auto jt: visited) {
        distance[jt->id] = 0;
        predecessor[jt->id] = nullptr;
      }

      graph_.RemoveDiploidEdges(marked_edges, true);
      num_bubbles += 1 - marked_edges.empty();
      std::cerr << "[raven::Graph::Assemble] RemoveBubbles: num_bubbles=" << num_bubbles << std::endl;
    }

    return num_bubbles;
  }

  std::uint32_t Graph_Assembler::RemoveLongEdges(std::uint32_t num_rounds) {
    std::uint32_t num_long_edges = 0;

    for (std::uint32_t i = 0; i < num_rounds; ++i) {
      CreateForceDirectedLayout();

      std::unordered_set<std::uint32_t> marked_edges;
      for (const auto &it: graph_.nodes_) {
        if (it == nullptr || it->outdegree() < 2) {
          continue;
        }
        for (auto jt: it->outedges) {
          for (auto kt: it->outedges) {
            if (jt != kt && jt->weight * 2.0 < kt->weight) {  // TODO(rvaser)
              marked_edges.emplace(kt->id);
              marked_edges.emplace(kt->pair->id);
            }
          }
        }
      }
      graph_.RemoveEdges(marked_edges);
      num_long_edges += marked_edges.size() / 2;

      RemoveTips();
    }

    return num_long_edges;
  }

  void Graph_Assembler::CreateForceDirectedLayout(const std::string &path) {
    std::ofstream os;
    bool is_first = true;
    if (!path.empty()) {
      os.open(path);
      os << "{" << std::endl;
    }

    std::vector<std::unordered_set<std::uint32_t>> components;
    std::vector<char> is_visited(graph_.nodes_.size(), 0);
    for (std::uint32_t i = 0; i < graph_.nodes_.size(); ++i) {
      if (graph_.nodes_[i] == nullptr || is_visited[i]) {
        continue;
      }

      components.resize(components.size() + 1);

      std::deque<std::uint32_t> que = {i};
      while (!que.empty()) {
        std::uint32_t j = que.front();
        que.pop_front();

        if (is_visited[j]) {
          continue;
        }
        const auto &node = graph_.nodes_[j];
        is_visited[node->id] = 1;
        is_visited[node->pair->id] = 1;
        components.back().emplace((node->id >> 1) << 1);

        for (auto it: node->inedges) {
          que.emplace_back(it->tail->id);
        }
        for (auto it: node->outedges) {
          que.emplace_back(it->head->id);
        }
      }
    }
    std::vector<char>().swap(is_visited);

    std::sort(components.begin(), components.end(),
              [](const std::unordered_set<std::uint32_t> &lhs,
                 const std::unordered_set<std::uint32_t> &rhs) {
                return lhs.size() > rhs.size();
              });

    static std::uint64_t seed = 21;
    seed <<= 1;

    std::mt19937 generator(seed);
    std::uniform_real_distribution<> distribution(0., 1.);

    struct Point {
      Point() = default;

      Point(double x, double y)
          : x(x),
            y(y) {}

      bool operator==(const Point &other) const {
        return x == other.x && y == other.y;
      }

      Point operator+(const Point &other) const {
        return Point(x + other.x, y + other.y);
      }

      Point &operator+=(const Point &other) {
        x += other.x;
        y += other.y;
        return *this;
      }

      Point operator-(const Point &other) const {
        return Point(x - other.x, y - other.y);
      }

      Point operator*(double c) const {
        return Point(x * c, y * c);
      }

      Point &operator/=(double c) {
        x /= c;
        y /= c;
        return *this;
      }

      double Norm() const {
        return sqrt(x * x + y * y);
      }

      double x;
      double y;
    };

    struct Quadtree {
      Quadtree(Point nucleus, double width)
          : nucleus(nucleus),
            width(width),
            center(0, 0),
            mass(0),
            subtrees() {
      }

      bool Add(const Point &p) {
        if (nucleus.x - width > p.x || p.x > nucleus.x + width ||
            nucleus.y - width > p.y || p.y > nucleus.y + width) {
          return false;
        }
        ++mass;
        if (mass == 1) {
          center = p;
        } else if (subtrees.empty()) {
          if (center == p) {
            return true;
          }
          double w = width / 2;
          subtrees.emplace_back(Point(nucleus.x + w, nucleus.y + w), w);
          subtrees.emplace_back(Point(nucleus.x - w, nucleus.y + w), w);
          subtrees.emplace_back(Point(nucleus.x - w, nucleus.y - w), w);
          subtrees.emplace_back(Point(nucleus.x + w, nucleus.y - w), w);
          for (auto &it: subtrees) {
            if (it.Add(center)) {
              break;
            }
          }
        }
        for (auto &it: subtrees) {
          if (it.Add(p)) {
            break;
          }
        }
        return true;
      }

      void Centre() {
        if (subtrees.empty()) {
          return;
        }
        center = Point(0, 0);
        for (auto &it: subtrees) {
          it.Centre();
          center += it.center * it.mass;
        }
        center /= mass;
      }

      Point Force(const Point &p, double k) const {
        auto delta = p - center;
        auto distance = delta.Norm();
        if (width * 2 / distance < 1) {
          return delta * (mass * (k * k) / (distance * distance));
        }
        delta = Point(0, 0);
        for (const auto &it: subtrees) {
          delta += it.Force(p, k);
        }
        return delta;
      }

      Point nucleus;
      double width;
      Point center;
      std::uint32_t mass;
      std::vector<Quadtree> subtrees;
    };

    std::uint32_t c = 0;
    for (const auto &component: components) {
      if (component.size() < 6) {
        continue;
      }

      bool has_junctions = false;
      for (const auto &it: component) {
        if (graph_.nodes_[it]->is_junction()) {
          has_junctions = true;
          break;
        }
      }
      if (has_junctions == false) {
        continue;
      }

      // update transitive edges
      for (const auto &n: component) {
        std::unordered_set<std::uint32_t> valid;
        for (const auto &m: graph_.nodes_[n]->transitive) {
          if (component.find(m) != component.end()) {
            valid.emplace(m);
          }
        }
        graph_.nodes_[n]->transitive.swap(valid);
      }

      std::uint32_t num_iterations = 100;
      double k = sqrt(1. / static_cast<double>(component.size()));
      double t = 0.1;
      double dt = t / static_cast<double>(num_iterations + 1);

      std::vector<Point> points(graph_.nodes_.size());
      for (const auto &it: component) {
        points[it].x = distribution(generator);
        points[it].y = distribution(generator);
      }

      for (std::uint32_t i = 0; i < num_iterations; ++i) {
        Point x = {0, 0}, y = {0, 0};
        for (const auto &n: component) {
          x.x = std::min(x.x, points[n].x);
          x.y = std::max(x.y, points[n].x);
          y.x = std::min(y.x, points[n].y);
          y.y = std::max(y.y, points[n].y);
        }
        double w = (x.y - x.x) / 2, h = (y.y - y.x) / 2;

        Quadtree tree(Point(x.x + w, y.x + h), std::max(w, h) + 0.01);
        for (const auto &n: component) {
          tree.Add(points[n]);
        }
        tree.Centre();

        std::vector<std::future<void>> thread_futures;
        std::vector<Point> displacements(graph_.nodes_.size(), Point(0, 0));

        auto thread_task = [&](std::uint32_t n) -> void {
          auto displacement = tree.Force(points[n], k);
          for (auto e: graph_.nodes_[n]->inedges) {
            auto m = (e->tail->id >> 1) << 1;
            auto delta = points[n] - points[m];
            auto distance = delta.Norm();
            if (distance < 0.01) {
              distance = 0.01;
            }
            displacement += delta * (-1. * distance / k);
          }
          for (auto e: graph_.nodes_[n]->outedges) {
            auto m = (e->head->id >> 1) << 1;
            auto delta = points[n] - points[m];
            auto distance = delta.Norm();
            if (distance < 0.01) {
              distance = 0.01;
            }
            displacement += delta * (-1. * distance / k);
          }
          for (const auto &m: graph_.nodes_[n]->transitive) {
            auto delta = points[n] - points[m];
            auto distance = delta.Norm();
            if (distance < 0.01) {
              distance = 0.01;
            }
            displacement += delta * (-1. * distance / k);
          }
          auto length = displacement.Norm();
          if (length < 0.01) {
            length = 0.1;
          }
          displacements[n] = displacement * (t / length);
          return;
        };

        for (const auto &n: component) {
          thread_futures.emplace_back(thread_pool_->Submit(thread_task, n));
        }
        for (const auto &it: thread_futures) {
          it.wait();
        }
        for (const auto &n: component) {
          points[n] += displacements[n];
        }

        t -= dt;
      }

      for (const auto &it: graph_.edges_) {
        if (it == nullptr || it->id & 1) {
          continue;
        }
        auto n = (it->tail->id >> 1) << 1;
        auto m = (it->head->id >> 1) << 1;

        if (component.find(n) != component.end() &&
            component.find(m) != component.end()) {
          it->weight = (points[n] - points[m]).Norm();
          it->pair->weight = it->weight;
        }
      }

      if (!path.empty()) {
        if (!is_first) {
          os << "," << std::endl;
        }
        is_first = false;

        os << "    \"component_" << c++ << "\": {" << std::endl;

        bool is_first_node = true;
        os << "      \"nodes\": {" << std::endl;
        for (const auto &it: component) {
          if (!is_first_node) {
            os << "," << std::endl;
          }
          is_first_node = false;
          os << "        \"" << it << "\": [";
          os << points[it].x << ", ";
          os << points[it].y << ", ";
          os << (graph_.nodes_[it]->is_junction() ? 1 : 0) << ", ";
          os << graph_.nodes_[it]->count << "]";
        }
        os << std::endl << "      }," << std::endl;

        bool is_first_edge = true;
        os << "      \"edges\": [" << std::endl;
        for (const auto &it: component) {
          for (auto e: graph_.nodes_[it]->inedges) {
            auto o = (e->tail->id >> 1) << 1;
            if (it < o) {
              continue;
            }
            if (!is_first_edge) {
              os << "," << std::endl;
            }
            is_first_edge = false;
            os << "        [\"" << it << "\", \"" << o << "\", 0]";
          }
          for (auto e: graph_.nodes_[it]->outedges) {
            auto o = (e->head->id >> 1) << 1;
            if (it < o) {
              continue;
            }
            if (!is_first_edge) {
              os << "," << std::endl;
            }
            is_first_edge = false;
            os << "        [\"" << it << "\", \"" << o << "\", 0]";
          }
          for (const auto &o: graph_.nodes_[it]->transitive) {
            if (it < o) {
              continue;
            }
            if (!is_first_edge) {
              os << "," << std::endl;
            }
            is_first_edge = false;
            os << "        [\"" << it << "\", \"" << o << "\", 1]";
          }
        }
        os << std::endl << "      ]" << std::endl;
        os << "    }";
      }
    }

    if (!path.empty()) {
      os << std::endl << "}";
      os << std::endl;
      os.close();
    }
  }

  void Graph_Assembler::ResolveGraphWithUl(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_reads) {
    //get unitig sequences
    //
  }
} // raven