#include "edge.h"

namespace raven {

  Edge::Edge(Node *tail, Node *head, std::uint32_t length)
      : id(num_objects++),
        length(length),
        weight(0),
        tail(tail),
        head(head),
        pair(),
        alternate() {
    tail->outedges.emplace_back(this);
    head->inedges.emplace_back(this);
  }

  Edge::Edge(Node *tail, Node *head, std::uint32_t length, std::uint32_t id)
      : id(id),
        length(length),
        weight(0),
        tail(tail),
        head(head),
        pair(),
        alternate() {
    tail->outedges.emplace_back(this);
    head->inedges.emplace_back(this);
  }

  std::string Edge::Label() const {
    return tail->sequence.InflateData(0, length);
  }

  std::atomic<std::uint32_t> Edge::num_objects{0};
  std::atomic<std::uint32_t> Edge::num_objects_alternate{0};
} // raven