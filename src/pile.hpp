// Copyright (c) 2020 Robert Vaser

#ifndef RAVEN_PILE_HPP_
#define RAVEN_PILE_HPP_

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "cereal/cereal.hpp"
#include "cereal/access.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"
#include "graph_constructor/extended_overlap.h"

namespace raven {

  constexpr std::uint32_t kPSS = 4;  // shrink 2 ^ kPSS times
  enum class kMerType : std::uint8_t { Haploid = 0, Diploid = 1, Repetitive = 2, Error = 3, None = 4 };
  class Pile {
  public:
    Pile(std::uint32_t id, std::uint32_t len);

    Pile(const Pile &) = delete;

    Pile &operator=(const Pile &) = delete;

    Pile(Pile &&) = default;

    Pile &operator=(Pile &&) = default;

    ~Pile() = default;

    std::vector<std::uint16_t> get_data() const {
      return data_;
    }

    std::vector<std::uint16_t> get_sketch_data() const {
      return sketch_data_;
    }

    std::vector<std::uint64_t> get_k_kmer_ids() const {
      return k_mer_ids_;
    }

    std::vector<kMerType> get_kmer_types() const {
      return kmer_types_;
    }

    float return_repetitve(std::uint32_t start, std::uint32_t end){
      if(kmer_types_.empty() || start >= end || end > kmer_types_.size()){
        return 0.0;
    }

      float total = 0.0;
      for(std::uint32_t i = start; i < end; ++i){
        if(kmer_types_[i] == kMerType::Repetitive){
          total += 1.0;
        }
      }
      return total / (end - start);
    }

    float return_haploid(std::uint32_t start, std::uint32_t end){
      if(kmer_types_.empty() || start >= end || end > kmer_types_.size()){
        return 0.0;
      }

      float total = 0.0;
      for(std::uint32_t i = start; i < end; ++i){
        if(kmer_types_[i] == kMerType::Haploid){
          total += 1.0;
        }
      }
      return total / (end - start);
    }
    
    float return_diploid(std::uint32_t start, std::uint32_t end){
      if(kmer_types_.empty() || start >= end || end > kmer_types_.size()){
        return 0.0;
      }

      float total = 0.0;
      for(std::uint32_t i = start; i < end; ++i){
        if(kmer_types_[i] == kMerType::Diploid){
          total += 1.0;
        }
      }
      return total / (end - start);

    }

    float return_erroneous(std::uint32_t start, std::uint32_t end){

      if(kmer_types_.empty() || start >= end || end > kmer_types_.size()){
        return 0.0;
      }

      float total = 0.0;
      for(std::uint32_t i = start; i < end; ++i){
        if(kmer_types_[i] == kMerType::Error){
          total += 1.0;
        }
      }
      return total / (end - start);
    }


    void check_HOR(std::uint32_t hom_peak);

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

    // classify k-mers in sketch data
    void classify_sketch_kmers(std::uint32_t hom_peak, std::uint32_t window_size, std::uint32_t kmer_len, double downsample_factor);

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

    void set_sketch(std::vector<std::uint16_t> sketch);

    void set_k_kmer_ids(std::vector<std::uint64_t> k_mer_ids);

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
    std::vector<std::uint16_t> data_;
    std::vector<std::uint16_t> sketch_data_;
    std::vector<std::uint64_t> k_mer_ids_;
    std::vector<bool> kmers_;
    std::vector<kMerType> kmer_types_;
    std::vector<Region> chimeric_regions_;
    std::vector<Region> repetitive_regions_;
  };

}  // namespace raven

#endif  // RAVEN_PILE_HPP_
