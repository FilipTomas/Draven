
#ifndef RAVEN_MARKED_EDGE_H
#define RAVEN_MARKED_EDGE_H

#include <cstdint>

namespace raven {

  struct MarkedEdge {
  public:
    MarkedEdge() = default;

    ~MarkedEdge() = default;

    MarkedEdge(std::uint32_t id)
        : id(id),
          where(0) {}

    MarkedEdge(std::uint32_t id, int where)
        : id(id),
          where(where) {}

    std::uint32_t id;
    int where = 0;
  };

}

#endif //RAVEN_MARKED_EDGE_H
