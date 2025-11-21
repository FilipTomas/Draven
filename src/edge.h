#ifndef RAVEN_EDGE_H
#define RAVEN_EDGE_H

#include <cstdint>
#include "node.h"

namespace raven {
  struct Node;

  struct Edge {
  public:
    Edge() = default;  // needed for cereal

    Edge(Node *tail, Node *head, std::uint32_t length);

    Edge(Node *tail, Node *head, std::uint32_t length, std::uint32_t id);

    Edge(const Edge &) = delete;

    Edge &operator=(const Edge &) = delete;

    ~Edge() = default;

    bool is_rc() const {
      return id & 1;
    }

    std::string Label() const;

    template<class Archive>
    void serialize(Archive &archive) {  // NOLINT
      archive(id, length, weight);
    }

    static std::atomic<std::uint32_t> num_objects;
    static std::atomic<std::uint32_t> num_objects_alternate;

    std::uint32_t id;
    std::uint32_t length;
    double weight;
    Node *tail;
    Node *head;
    Edge *pair;
    Edge *alternate;
  };

} // raven

#endif //RAVEN_EDGE_H
