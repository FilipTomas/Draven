// Copyright (c) 2020 Robert Vaser

#ifndef RAVEN_PILE_HPP_
#define RAVEN_PILE_HPP_

#include <algorithm>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>
#include <set>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "cereal/cereal.hpp"
#include "cereal/access.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"
#include "graph_constructor/extended_overlap.h"

namespace raven {

  constexpr std::uint32_t kPSS = 4;  // shrink 2 ^ kPSS times
  enum class kMerType : std::uint8_t { Haploid = 0, Diploid = 1, Error = 2, Repetitive = 3, None = 4 };
  class Pile {
  public:
    // coverage: allocate the per-base coverage (2 B per base) that the
    // assembly path's layers, median and region finding use; the correction
    // path uses only k-mer types and flags.
    Pile(std::uint32_t id, std::uint32_t len, bool coverage = true);

    Pile(const Pile &) = delete;

    Pile &operator=(const Pile &) = delete;

    Pile(Pile &&) = default;

    Pile &operator=(Pile &&) = default;

    ~Pile() = default;

    std::vector<std::uint16_t> get_data() const {
      return data_;
    }

    std::vector<std::uint64_t> get_k_kmer_ids() const {
      return k_mer_ids_;
    }

    // per k-mer position (expanded from the stored windows)
    std::vector<kMerType> get_kmer_types() const {
      std::vector<kMerType> types(num_kmer_types_);
      for (std::uint32_t i = 0; i < num_kmer_types_; ++i) {
        types[i] = kmer_type_at(i);
      }
      return types;
    }

    float return_repetitve(std::uint32_t start, std::uint32_t end){
      return kmer_type_fraction(start, end, kMerType::Repetitive);
    }

    float return_haploid(std::uint32_t start, std::uint32_t end){
      return kmer_type_fraction(start, end, kMerType::Haploid);
    }

    float return_diploid(std::uint32_t start, std::uint32_t end){
      return kmer_type_fraction(start, end, kMerType::Diploid);
    }

    float return_erroneous(std::uint32_t start, std::uint32_t end){
      return kmer_type_fraction(start, end, kMerType::Error);
    }

    float return_ratio(std::uint32_t start, uint32_t end, kMerType kmer_type){
      return kmer_type_fraction(start, end, kmer_type);
    }

    std::uint32_t id() const {
      return id_;
    }

    std::uint32_t begin() const {
      return begin_ << kPSS;
    }

    std::uint32_t end() const {
      return end_ << kPSS;
    }

    std::uint32_t length() const {
      return end() - begin();
    }

    std::uint16_t median() const {
      return median_;
    }

    bool is_invalid() const {
      return is_invalid_;
    }

    void set_is_invalid() {
      is_invalid_ = true;
    }

    bool is_contained() const {
      return is_contained_;
    }

    void set_is_contained() {
      is_contained_ = true;
    }

    bool is_chimeric() const {
      return is_chimeric_;
    }

    bool is_maybe_chimeric() const {
      return !chimeric_regions_.empty();
    }

    void set_is_chimeric() {
      is_chimeric_ = true;
    }

    bool is_repetitive() const {
      return is_repetitive_;
    }

    void set_is_repetitive() {
      is_repetitive_ = true;
    }

    bool is_hor() const {
      return is_hor_;
    }

    void set_is_hor() {
      is_hor_ = true;
    }

    float get_hor_percent() const {
      return hor_percent_;
    }

    void set_hor_percent(float p) {
      hor_percent_ = p;
    }

    // one type per k-mer position
    void set_k_mer_types(std::vector<kMerType> tmp){
      num_kmer_types_ = tmp.size();
      kmer_types_ = std::move(tmp);
      kmer_type_window_ = 1;
    }

    // One type per `window` k-mer positions (the CNN's per-window classes),
    // for num_positions positions; positions past the last window take its
    // type. 1/window of the memory of set_k_mer_types, same lookups.
    void set_k_mer_type_windows(std::vector<kMerType> windows,
                                std::uint32_t window,
                                std::uint32_t num_positions){
      kmer_types_ = std::move(windows);
      kmer_type_window_ = window;
      num_kmer_types_ = kmer_types_.empty() ? 0 : num_positions;
    }

    void clear_k_mer_types(){
      std::vector<kMerType>().swap(kmer_types_);
      num_kmer_types_ = 0;
    }

    void set_backbone(){
      backbone_pile_ = true;
    }

    bool is_backbone() const {
      return backbone_pile_;
    }

    void set_strong_contained(){
      strong_contained_ = true;
    }

    bool is_strong_contained() const {
      return strong_contained_;
    }

    // add coverage
    void AddLayers(
        std::vector<biosoup::Overlap>::const_iterator begin,
        std::vector<biosoup::Overlap>::const_iterator end);

    void AddExtendedLayers(
      std::vector<extended_overlap>::const_iterator begin,
      std::vector<extended_overlap>::const_iterator end);

    // mark repetitive k-mers
    void AddKmers(
        const std::vector<std::uint32_t> &kmers,
        std::uint32_t kmer_len,
        const std::unique_ptr<biosoup::NucleicAcid> &sequence);

    // store longest region with values greater or equal than given coverage
    void FindValidRegion(std::uint16_t coverage, std::uint16_t length);

    // fill valid region with zeroes
    void ClearValidRegion();

    // fill everything outside valid region with zeroes
    void ClearInvalidRegion();

    // store median of valid region
    void FindMedian();

    // store coverage drops
    void FindChimericRegions();

    // update valid region to longest non-chimeric given the component median
    void ClearChimericRegions(std::uint16_t median);

    // store coverage spikes given component median, and
    //   tightly packed groups of repetitive k-mers
    void FindRepetitiveRegions(std::uint16_t median);

    // increase confidence in repetitive regions given an overlap
    void UpdateRepetitiveRegions(const biosoup::Overlap &o);

    // define relationship between repetitive regions and given overlap
    bool CheckRepetitiveRegions(const biosoup::Overlap &o);

    // remove all repetitive regions
    void ClearRepetitiveRegions();

    void set_k_kmer_ids(std::vector<std::uint64_t> k_mer_ids);

    std::vector<std::pair<std::uint32_t, std::uint32_t>> find_region_positions(kMerType type, std::uint32_t span_start, std::uint32_t span_end);

    // Append the k-mer ids in [start, end] to `out` (no dedup/sort here;
    // callers concatenate across regions then sort+unique once).
    void append_kmers_in_region(std::uint32_t start, std::uint32_t end,
                                std::vector<std::uint64_t>& out) const;

    
  private:
    Pile() = default;

    template<class Archive>
    void serialize(Archive &archive) {  // NOLINT
      archive(
          CEREAL_NVP(id_),
          CEREAL_NVP(begin_),
          CEREAL_NVP(end_),
          CEREAL_NVP(median_),
          CEREAL_NVP(is_invalid_),
          CEREAL_NVP(is_contained_),
          CEREAL_NVP(is_chimeric_),
          CEREAL_NVP(is_repetitive_),
          CEREAL_NVP(data_),
          CEREAL_NVP(kmers_),
          CEREAL_NVP(chimeric_regions_),
          CEREAL_NVP(repetitive_regions_));
    }

    friend cereal::access;

    using Region = std::pair<std::uint32_t, std::uint32_t>;

    // clear invalid region after update
    void UpdateValidRegion(std::uint32_t begin, std::uint32_t end);

    // merge overlapping regions
    static std::vector<Region> MergeRegions(const std::vector<Region> &regions);

    // find drop and spike regions
    std::vector<Region> FindSlopes(double q);

    std::uint32_t id_;
    std::uint32_t begin_;
    std::uint32_t end_;
    std::uint16_t median_;
    bool is_invalid_;
    bool is_contained_;
    bool is_chimeric_;
    bool is_repetitive_;
    bool is_hor_;
    float hor_percent_ = 0.f;
    std::vector<std::uint16_t> data_;
    std::vector<std::uint64_t> k_mer_ids_;
    std::vector<bool> kmers_;
    // k-mer types, one per kmer_type_window_ positions, for num_kmer_types_
    // positions
    std::vector<kMerType> kmer_types_;
    std::uint32_t kmer_type_window_ = 1;
    std::uint32_t num_kmer_types_ = 0;

    kMerType kmer_type_at(std::uint32_t i) const {
      const std::size_t w = std::min<std::size_t>(i / kmer_type_window_,
                                                  kmer_types_.size() - 1);
      return kmer_types_[w];
    }

    // f(first, last, type) for consecutive runs of positions in [first, last]
    // (last < num_kmer_types_) that share a window, hence a type
    template <typename F>
    void for_each_kmer_type_segment(std::uint32_t first, std::uint32_t last,
                                    F f) const {
      const std::size_t last_window = kmer_types_.size() - 1;
      for (std::uint32_t i = first; i <= last;) {
        const std::size_t w = std::min<std::size_t>(i / kmer_type_window_, last_window);
        const std::uint32_t segment_last = w == last_window ?
            last :
            std::min<std::uint32_t>(last, (w + 1) * kmer_type_window_ - 1);
        f(i, segment_last, kmer_types_[w]);
        if (segment_last == last) {
          break;
        }
        i = segment_last + 1;
      }
    }

    // fraction of positions in [start, end) of the given type
    float kmer_type_fraction(std::uint32_t start, std::uint32_t end,
                             kMerType type) const {
      if (num_kmer_types_ == 0 || start >= end || end > num_kmer_types_) {
        return 0.0;
      }
      std::uint32_t total = 0;
      for_each_kmer_type_segment(start, end - 1,
          [&](std::uint32_t a, std::uint32_t b, kMerType t) {
            if (t == type) {
              total += b - a + 1;
            }
          });
      return static_cast<float>(total) / (end - start);
    }
    std::vector<Region> chimeric_regions_;
    std::vector<Region> repetitive_regions_;
    bool backbone_pile_;
    bool strong_contained_;
  };

}  // namespace raven

#endif  // RAVEN_PILE_HPP_
