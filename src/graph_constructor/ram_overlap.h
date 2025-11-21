#ifndef RAVEN_RAM_OVERLAP_H
#define RAVEN_RAM_OVERLAP_H

#include <cstdint>
#include <vector>

#include "biosoup/overlap.hpp"
#include "extended_overlap.h"


class ram_overlap{
public:
    std::vector<extended_overlap> overlap_fragments;
    std::uint32_t lhs_id;
    std::uint32_t rhs_id;

    std::uint32_t lhs_begin;
    std::uint32_t lhs_end;
    std::uint32_t rhs_begin;
    std::uint32_t rhs_end;

    std::uint16_t n_fragments;
    std::vector<std::pair<std::uint16_t, std::uint16_t>> fragment_gaps;


    extended_overlap resolve();
};


#endif // RAVEN_RAM_OVERLAP_H