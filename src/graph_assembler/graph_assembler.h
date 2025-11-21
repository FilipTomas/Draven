
#ifndef RAVEN_GRAPH_ASSEMBLER_H
#define RAVEN_GRAPH_ASSEMBLER_H

#include <cstdint>
#include <memory>
#include "biosoup/nucleic_acid.hpp"
#include "graph.hpp"
#include "option_manager.h"
#include "marked_edge.h"

namespace raven {

  class Graph_Assembler {
  public:
    Graph_Assembler() = default;

    Graph_Assembler(Graph &graph, Program_Parameters &param, std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

    // simplify with transitive reduction, tip prunning and bubble popping
    void Assemble();
    
    void AssembleHaploids();
    
    void AssembleDiploids();

    void UlAssemble(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_sequences);

    // inspired by (Myers 1995) & (Myers 2005)
    std::uint32_t RemoveTransitiveEdges();

    std::uint32_t RemoveTips();

    std::uint32_t RemoveBubbles();

    std::uint32_t RemoveSnpBubbles();

    // remove long edges in force directed layout
    std::uint32_t RemoveLongEdges(std::uint32_t num_round);

  private:
    Graph &graph_;
    Program_Parameters &param_;
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

    // use (Fruchterman & Reingold 1991) with (Barnes & Hut 1986) approximation
    // (draw with misc/plotter.py)
    void CreateForceDirectedLayout(const std::string &path = "");

    void ResolveGraphWithUl(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_reads);
  };

} // raven

#endif //RAVEN_GRAPH_ASSEMBLER_H
