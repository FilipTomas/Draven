#include "node.h"

namespace raven {

  Node::Node(const biosoup::NucleicAcid &sequence)
      : id(num_objects++),
        sequence(sequence),
        count(1),
        is_unitig(),
        is_circular(),
        is_polished(),
        transitive(),
        inedges(),
        outedges(),
        pair(),
        alternate() {}

  Node::Node(const biosoup::NucleicAcid &sequence, std::uint32_t id)
      : id(id),
        sequence(sequence),
        count(1),
        is_unitig(),
        is_circular(),
        is_polished(),
        transitive(),
        inedges(),
        outedges(),
        pair(),
        alternate() {}

  Node::Node(Node *begin, Node *end)
      : id(num_objects++),
        sequence(),
        count(),
        is_unitig(),
        is_circular(begin == end),
        is_polished(),
        transitive(),
        inedges(),
        outedges(),
        pair(),
        alternate() {
    std::string data{};

    auto it = begin;
    while (true) {
      data += it->outedges.front()->Label();
      count += it->count;
      if ((it = it->outedges.front()->head) == end) {
        break;
      }
    }
    if (begin != end) {
      data += end->sequence.InflateData();
      count += end->count;
    }

    is_unitig = count > 5 && data.size() > min_unitig_size;

    sequence = biosoup::NucleicAcid(
        (is_unitig ? "Utg" : "Ctg") + std::to_string(id & (~1UL)),
        data);
  }

  Node::Node(Node *begin, Node *end, std::uint32_t id)
      : id(id),
        sequence(),
        count(),
        is_unitig(),
        is_circular(begin == end),
        is_polished(),
        transitive(),
        inedges(),
        outedges(),
        pair(),
        alternate() {

    std::string data{};

    auto it = begin;
    while (true) {
      data += it->outedges.front()->Label();
      count += it->count;
      if ((it = it->outedges.front()->head) == end) {
        break;
      }
    }
    if (begin != end) {
      data += end->sequence.InflateData();
      count += end->count;
    }

    is_unitig = count > 5 && data.size() > 9999;

    sequence = biosoup::NucleicAcid(
        (is_unitig ? "Utg" : "Ctg") + std::to_string(id & (~1UL)),
        data);
  }

  Node::Node(Node *begin, Node *end, bool is_unitig)
      : id(num_objects++),
        sequence(),
        count(),
        is_unitig(is_unitig),
        is_circular(false),
        is_polished(),
        transitive(),
        front_inedges(),
        front_outedges(),
        back_inedges(),
        back_outedges(),
        pair(),
        unitig_nodes(),
        alternate() {

    std::string data{};

    auto it = begin;
    if (begin == end) {
      data += it->sequence.InflateData();
      count += it->count;
      unitig_nodes.emplace_back(it);
    } else {
      while (true) {
        data += it->outedges.front()->Label();
        count += it->count;
        unitig_nodes.emplace_back(it);
        if ((it = it->outedges.front()->head) == end) {
          unitig_nodes.emplace_back(it);
          break;
        }
      }

      if (begin != end) {
        data += end->sequence.InflateData();
        count += end->count;
      }
    }


    std::uint32_t front_inedge_count{0};
    std::uint32_t back_inedge_count{0};
    std::uint32_t front_outedge_count{0};
    std::uint32_t back_outedge_count{0};

    for (auto &inedge: begin->inedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), inedge->tail) == unitig_nodes.end()) {
        front_inedges.emplace_back(inedge);
        front_inedge_count++;
      }
    }

    for (auto &outedge: begin->outedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), outedge->head) == unitig_nodes.end()) {
        front_outedges.emplace_back(outedge);
        front_outedge_count++;
      }
    }

    for (auto &inedge: end->inedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), inedge->tail) == unitig_nodes.end()) {
        back_inedges.emplace_back(inedge);
        back_inedge_count++;
      }
    }


    for (auto &outedge: end->outedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), outedge->head) == unitig_nodes.end()) {
        back_outedges.emplace_back(outedge);
        back_outedge_count++;
      }
    }
    sequence = biosoup::NucleicAcid(
        "Utg" + std::to_string(id & (~1UL)),
        data);
  }

  Node::Node(Node *begin, Node *end, bool is_unitig, std::uint32_t id)
      : id(id),
        sequence(),
        count(),
        is_unitig(is_unitig),
        is_circular(begin == end),
        is_polished(),
        transitive(),
        inedges(),
        outedges(),
        pair(),
        unitig_nodes(),
        alternate() {

    std::string data{};


    auto it = begin;
    if (begin == end) {
      data += it->sequence.InflateData();
      count += it->count;
      unitig_nodes.emplace_back(it);
    } else {
      while (true) {
        data += it->outedges.front()->Label();
        count += it->count;
        unitig_nodes.emplace_back(it);
        if ((it = it->outedges.front()->head) == end) {
          unitig_nodes.emplace_back(it);
          break;
        }
      }

      if (begin != end) {
        data += end->sequence.InflateData();
        count += end->count;
      }
    }


    std::uint32_t front_inedge_count{0};
    std::uint32_t back_inedge_count{0};
    std::uint32_t front_outedge_count{0};
    std::uint32_t back_outedge_count{0};

    for (auto &inedge: begin->inedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), inedge->tail) == unitig_nodes.end()) {
        front_inedges.emplace_back(inedge);
        front_inedge_count++;
      }
    }

    for (auto &outedge: begin->outedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), outedge->head) == unitig_nodes.end()) {
        front_outedges.emplace_back(outedge);
        front_outedge_count++;
      }
    }

    for (auto &inedge: end->inedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), inedge->tail) == unitig_nodes.end()) {
        back_inedges.emplace_back(inedge);
        back_inedge_count++;
      }
    }


    for (auto &outedge: end->outedges) {
      if (std::find(unitig_nodes.begin(), unitig_nodes.end(), outedge->head) == unitig_nodes.end()) {
        back_outedges.emplace_back(outedge);
        back_outedge_count++;
      }
    }
    sequence = biosoup::NucleicAcid(
        "Utg" + std::to_string(id & (~1UL)),
        data);
  }

  std::atomic<std::uint32_t> Node::num_objects{0};
  std::atomic<std::uint32_t> Node::num_objects_alternate{0};
  std::uint32_t min_unitig_size = 9999;
} // raven