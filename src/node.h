#ifndef RAVEN_NODE_H
#define RAVEN_NODE_H

#include <unordered_set>
#include "biosoup/nucleic_acid.hpp"
#include "edge.h"

namespace raven {

  extern std::uint32_t min_unitig_size;

  struct Edge;

  struct Node {
  public:
    Node() = default;  // needed for cereal

    explicit Node(const biosoup::NucleicAcid &sequence);

    explicit Node(const biosoup::NucleicAcid &sequence, std::uint32_t id);

    Node(Node *begin, Node *end);

    Node(Node *begin, Node *end, std::uint32_t id);

    Node(Node *begin, Node *end, bool is_unitig);

    Node(Node *begin, Node *end, bool is_unitig, std::uint32_t id);

    Node(const Node &) = delete;

    Node &operator=(const Node &) = delete;

    Node(Node &&) = default;

    Node &operator=(Node &&) = default;

    bool operator==(const Node &other) {
      return other.sequence.id == this->sequence.id;
    };

    ~Node() = default;

    std::uint32_t indegree() const {
      return inedges.size();
    }

    std::uint32_t outdegree() const {
      return outedges.size();
    }

    bool is_rc() const {
      return id & 1;
    }

    bool is_junction() const {
      return outdegree() > 1 || indegree() > 1;
    }

    bool is_tip() const {
      return outdegree() > 0 && indegree() == 0 && count < 6;
    }

    template<class Archive>
    void serialize(Archive &archive) {  // NOLINT
      archive(
          id,
          sequence,
          count,
          is_unitig,
          is_circular,
          is_polished,
          transitive,
          color);
    }

    static std::atomic<std::uint32_t> num_objects;
    static std::atomic<std::uint32_t> num_objects_alternate;

    std::uint32_t id;
    biosoup::NucleicAcid sequence;
    std::uint32_t count;
    bool is_unitig;
    bool is_circular;
    bool is_polished;
    std::unordered_set<std::uint32_t> transitive;
    unsigned color = 0;

    std::vector<Edge *> inedges;
    std::vector<Edge *> outedges;
    std::vector<Edge *> front_inedges;
    std::vector<Edge *> front_outedges;
    std::vector<Edge *> back_inedges;
    std::vector<Edge *> back_outedges;
    Node *pair;

    std::vector<Node *> unitig_nodes;

    Node *alternate;
    bool is_primary = true;
  };

} // raven

#endif //RAVEN_NODE_H
