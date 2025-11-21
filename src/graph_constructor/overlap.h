
#ifndef RAVEN_OVERLAP_H
#define RAVEN_OVERLAP_H

#include <string>
#include <vector>

namespace raven {
  struct OverlapDescriptor {
    std::string name;
    std::string type;
    std::string data;
  };

  struct Overlap {
  public:
    Overlap() = default;

    Overlap(
        const char *q_name, std::uint32_t q_name_len,
        std::uint32_t q_len,
        std::uint32_t q_begin,
        std::uint32_t q_end,
        char orientation,
        const char *t_name, std::uint32_t t_name_len,
        std::uint32_t t_len,
        std::uint32_t t_begin,
        std::uint32_t t_end,
        std::uint32_t score,
        std::uint32_t overlap_len,
        std::uint32_t quality,
        std::vector<OverlapDescriptor> descriptors)
        :
        lhs_id(0),
        lhs_begin(q_begin),
        lhs_end(q_end),
        rhs_id(0),
        rhs_begin(t_begin),
        rhs_end(t_end),
        score(score),
        strand(orientation == '+'),
        alignment(),
        q_name(q_name, q_name_len),
        q_len(q_len),
        t_name(t_name, t_name_len),
        t_len(t_len),
        overlap_len(overlap_len),
        paf_quality(quality),
        descriptors(descriptors) {}

    Overlap(
        const char *q_name, std::uint32_t q_name_len,
        std::uint32_t flag,
        const char *t_name, std::uint32_t t_name_len,
        std::uint32_t t_begin,
        std::uint32_t map_quality,
        const char *cigar, std::uint32_t cigar_len,
        const char *t_next_name, std::uint32_t t_next_name_len,
        std::uint32_t t_next_begin,
        std::uint32_t template_len,
        const char *data, std::uint32_t data_len,
        const char *quality, std::uint32_t quality_len,
        std::vector<OverlapDescriptor> descriptors)
        :
        lhs_id(0),
        lhs_begin(0),
        lhs_end(0),
        rhs_id(0),
        rhs_begin(t_begin),
        rhs_end(0),
        score(0),
        alignment(cigar, cigar_len),
        q_name(q_name, q_name_len),
        t_name(t_name, t_name_len),
        sam_flag(flag),
        map_quality(map_quality),
        t_next_name(t_next_name, t_next_name_len),
        t_next_begin(t_next_begin),
        template_len(template_len),
        data(data, data_len),
        sam_quality(quality, quality_len),
        descriptors(descriptors) {}

    Overlap(const Overlap &) = default;

    Overlap &operator=(const Overlap &) = default;

    Overlap(Overlap &&) = default;

    Overlap &operator=(Overlap &&) = default;

    ~Overlap() = default;

    std::uint32_t lhs_id;
    std::uint32_t lhs_begin;
    std::uint32_t lhs_end;
    std::uint32_t rhs_id;
    std::uint32_t rhs_begin;
    std::uint32_t rhs_end;
    std::uint32_t score;  // based on k-mer matches or alignment score
    bool strand;  // (optional) Watson-Crick strand
    std::string alignment;  // (optional) cigar string

    //PAF
    std::string q_name;
    std::uint32_t q_len;
    std::string t_name;
    std::uint32_t t_len;
    std::uint32_t overlap_len;
    std::uint32_t paf_quality;

    //SAM
    std::uint32_t sam_flag;
    std::uint32_t map_quality;
    std::string t_next_name;
    std::uint32_t t_next_begin;
    std::uint32_t template_len;
    std::string data;
    std::string sam_quality;

    //descriptors
    std::vector<OverlapDescriptor> descriptors;
  };
}

#endif //RAVEN_OVERLAP_H
