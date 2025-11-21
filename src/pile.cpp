// Copyright (c) 2020 Robert Vaser

#include "pile.hpp"

#include <algorithm>
#include <deque>
#include <limits>
#include <string>

namespace raven {

  template<typename T>
  T clamp(T v) {
    return (v < std::numeric_limits<std::uint16_t>::max()) ?
           v : std::numeric_limits<std::uint16_t>::max();
  }

  //enum class kMerType : std::uint8_t { Haploid = 0, Diploid = 1, Repetitive = 2, Error = 3 };

  Pile::Pile(std::uint32_t id, std::uint32_t len)
      : id_(id),
        begin_(0),
        end_(len),
        median_(0),
        is_invalid_(0),
        is_contained_(0),
        is_chimeric_(0),
        is_repetitive_(0),
        is_hor_(0),
        data_(end_, 0),
        sketch_data_(),
        k_mer_ids_(),
        kmers_(),
        chimeric_regions_(),
        repetitive_regions_() {}

  void Pile::AddLayers(
      std::vector<biosoup::Overlap>::const_iterator begin,
      std::vector<biosoup::Overlap>::const_iterator end) {
    if (begin >= end) {
      return;
    }

    std::vector<std::uint32_t> boundaries;
    for (auto it = begin; it != end; ++it) {
      if (it->lhs_id == id_) {
        boundaries.emplace_back(((it->lhs_begin >> kPSS) + 1) << 1);
        boundaries.emplace_back(((it->lhs_end >> kPSS) - 1) << 1 | 1);
      } else if (it->rhs_id == id_) {
        boundaries.emplace_back(((it->rhs_begin >> kPSS) + 1) << 1);
        boundaries.emplace_back(((it->rhs_end >> kPSS) - 1) << 1 | 1);
      }
    }
    std::sort(boundaries.begin(), boundaries.end());

    std::uint32_t coverage = 0;
    std::uint32_t last_boundary = 0;
    for (const auto &it: boundaries) {
      if (coverage > 0) {
        for (std::uint32_t i = last_boundary; i < (it >> 1); ++i) {
          data_[i] = clamp(data_[i] + coverage);
        }
      }
      last_boundary = it >> 1;
      coverage += it & 1 ? -1 : 1;
    }
  }

  void Pile::set_sketch(std::vector<std::uint16_t> sketch) {
    if (sketch.empty()) {
      return;
    }
    sketch_data_ = sketch;

  }

  void Pile::set_k_kmer_ids(std::vector<std::uint64_t> k_mer_ids) {
    if (k_mer_ids.empty()) {
      return;
    }
    k_mer_ids_ = k_mer_ids;

  }

  void Pile::check_HOR(std::uint32_t hom_peak) {
    if(sketch_data_.empty()) {
      return;
    }
    const auto cnt = std::count_if(sketch_data_.begin(), sketch_data_.end(), [hom_peak](std::uint16_t v) {
      return v >= hom_peak * 1.5;
    });
    if(cnt > sketch_data_.size() * 0.7) {
      set_is_hor();
    }
  }

  void Pile::classify_sketch_kmers(std::uint32_t hom_peak, std::uint32_t window_size, std::uint32_t kmer_len, double downsample_factor) {
    if(sketch_data_.empty() || k_mer_ids_.empty()) {
      return;
    }

    auto window_median_without_repetitve = [&](std::vector<std::uint16_t>::const_iterator begin,
                                         std::vector<std::uint16_t>::const_iterator end) -> std::uint16_t {
      std::vector<std::uint16_t> window(begin, end);
      window.erase(
          std::remove_if(
              window.begin(),
              window.end(),
              [hom_peak, downsample_factor](std::uint16_t v) { return (v / downsample_factor) > (hom_peak * 1.5); }
          ),
          window.end()
      );
      std::nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
      if (window.empty()) {
        return 0;
      } else {
        return window[window.size() / 2];
      }
    };

    auto classify_point = [&](std::uint16_t value, std::uint16_t median) -> kMerType {
      if (median == 0){
        return kMerType::Repetitive;
      }
      if (value) {
        if (median * 0.3 < value && value < median * 0.75) {
          return kMerType::Haploid;
        } else if (value <= 0.3 * median || (value / downsample_factor) < 0.1 * hom_peak) {
          return kMerType::Error;
        } else if ((value / downsample_factor) > hom_peak * 1.5) {
          return kMerType::Repetitive;
        } else {
          return kMerType::Diploid;
        }
      } else {
        return kMerType::Diploid;
      }
    };

    const std::size_t n = sketch_data_.size();
    const std::size_t W = window_size;

    std::size_t begin = 0;
    for (; begin + W <= n; begin += W) {
        auto first = sketch_data_.begin() + begin;
        auto last  = first + W;
        std::uint16_t median = window_median_without_repetitve(first, last);
        for (auto it = first; it != last; ++it)
            kmer_types_.emplace_back(classify_point(*it, median));
    };

    // final partial window [begin, n)
    if (begin < n) {
        auto first = sketch_data_.begin() + begin;
        auto last  = sketch_data_.begin() + n;
        auto first_for_median = last - W;
        std::uint16_t median = window_median_without_repetitve(first_for_median, last);
        for (auto it = first; it != last; ++it)
            kmer_types_.emplace_back(classify_point(*it, median));
    };
    kmer_types_.insert(kmer_types_.begin(), kmer_len - 1, kMerType::None);
  }

  // void Pile::PrintKmerTypes() {
  //   std::cout 
  //   for (std::size_t i = 0; i < kmer_types_.size(); ++i) {
  //     switch (kmer_types_[i]) {
  //       case kMerType::Haploid:
  //         std::cout << "H";
  //         break;
  //       case kMerType::Diploid:
  //         std::cout << "D";
  //         break;
  //       case kMerType::Repetitive:
  //         std::cout << "R";
  //         break;
  //       case kMerType::Error:
  //         std::cout << "E";
  //         break;
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  
  void Pile::AddExtendedLayers(
      std::vector<extended_overlap>::const_iterator begin,
      std::vector<extended_overlap>::const_iterator end) {
    if (begin >= end) {
      return;
    }

    std::vector<std::uint32_t> boundaries;
    std::uint32_t lhs_id_tmp;
    std::uint32_t rhs_id_tmp;

    if(id_ == 27 || id_ == 13){
      std::cout << "[raven::Pile::AddExtendedLayers] id: " << id_ << std::endl;
    }
    for (auto it = begin; it != end; ++it) {
      if(it->overlap.lhs_id == id_) {
        // boundaries.emplace_back(((it->overlap.lhs_begin >> kPSS) + 1) << 1);
        // boundaries.emplace_back(((it->overlap.lhs_end >> kPSS) - 1) << 1 | 1);
        boundaries.emplace_back(((it->overlap.lhs_begin >> kPSS)) << 1);
        boundaries.emplace_back(((it->overlap.lhs_end >> kPSS)) << 1 | 1);        
      } else if (it->overlap.rhs_id == id_) {
        // boundaries.emplace_back(((it->overlap.rhs_begin >> kPSS) + 1) << 1);
        // boundaries.emplace_back(((it->overlap.rhs_end >> kPSS) - 1) << 1 | 1);
        boundaries.emplace_back(((it->overlap.rhs_begin >> kPSS)) << 1);
        boundaries.emplace_back(((it->overlap.rhs_end >> kPSS)) << 1 | 1);
      }
      lhs_id_tmp = it->overlap.lhs_id;
      rhs_id_tmp = it->overlap.rhs_id;
    }
    std::sort(boundaries.begin(), boundaries.end());
    auto pile_id_debug = id_;
    std::uint32_t coverage = 0;
    std::uint32_t last_boundary = 0;
    auto data_size = data_.size();
    for (const auto &it2: boundaries) {
      if (coverage > 0) {
        for (std::uint32_t i = last_boundary; i < (it2 >> 1); ++i) {
          if(data_size <= i){
            // std::cerr << "[raven::Pile::AddExtendedLayers] data size: " << data_size
            //           << ", i: " << i << ", id: " << pile_id_debug << std::endl;
            continue;
          }
          data_[i] = clamp(data_[i] + coverage);
        }
      }
      last_boundary = it2 >> 1;
      coverage += it2 & 1 ? -1 : 1;
    }
  }  

  void Pile::AddKmers(
      const std::vector<std::uint32_t> &kmers,
      std::uint32_t kmer_len,
      const std::unique_ptr<biosoup::NucleicAcid> &sequence) {
    if (kmers.empty()) {
      return;
    }
    if (kmers_.empty()) {
      kmers_.resize(data_.size() + 1, 0);
    }
    for (const auto &it: kmers) {
      // filter low-complexity regions
      std::string kmer = sequence->InflateData(it, kmer_len);

      std::vector<std::string> polymers;
      for (const auto &it: kmer) {
        polymers.emplace_back(1, it);
      }
      polymers.erase(std::unique(polymers.begin(), polymers.end()), polymers.end());  // NOLINT
      kmer = std::accumulate(polymers.begin(), polymers.end(), std::string());
      if (kmer.size() < kmer_len / 2 + 1) {
        continue;
      }

      polymers.clear();
      for (auto it = kmer.begin(); it != kmer.end(); ++it) {
        if ((it - kmer.begin()) % 2 == 1) {
          polymers.back() += *it;
        } else {
          polymers.emplace_back(1, *it);
        }
      }
      polymers.erase(std::unique(polymers.begin(), polymers.end()), polymers.end());  // NOLINT
      kmer = std::accumulate(polymers.begin(), polymers.end(), std::string());
      if (kmer.size() < kmer_len / 2 + 1) {
        continue;
      }

      polymers.clear();
      for (auto it = kmer.begin(); it != kmer.end(); ++it) {
        if (!polymers.empty() && (it - kmer.begin()) % 2 == 0) {
          polymers.back() += *it;
        } else {
          polymers.emplace_back(1, *it);
        }
      }
      polymers.erase(std::unique(polymers.begin(), polymers.end()), polymers.end());  // NOLINT
      kmer = std::accumulate(polymers.begin(), polymers.end(), std::string());
      if (kmer.size() < kmer_len / 2 + 1) {
        continue;
      }

      kmers_[it >> kPSS] = 1;
    }
  }

  void Pile::FindValidRegion(std::uint16_t coverage_threshold, std::uint16_t length_threshold){
    std::uint32_t begin = 0;
    std::uint32_t end = 0;
    bool any_above_th_flag = false;
    for (std::uint32_t i = begin_; i < end_; ++i) {
      if (data_[i] < coverage_threshold) {
        continue;
      }
      for (std::uint32_t j = i + 1; j < end_; ++j) {
        if (data_[j] >= coverage_threshold) {
          any_above_th_flag = true;
          continue;
        }
        if (end - begin < j - i) {
          begin = i;
          end = j;
        }
        i = j;
        break;
      }
    }
    if(any_above_th_flag && end == 0){
      end = end_;
    }
    UpdateValidRegion(begin, end);
  }

  void Pile::UpdateValidRegion(std::uint32_t begin, std::uint32_t end) {
    if (begin >= end || end - begin < 32 >> kPSS){
      set_is_invalid();
      return;
    }
    for (std::uint32_t i = begin_; i < begin; ++i) {
      data_[i] = 0;
    }
    for (std::uint32_t i = end; i < end_; ++i) {
      data_[i] = 0;
    }
    begin_ = begin;
    end_ = end;
  }

  void Pile::ClearValidRegion() {
    std::fill(data_.begin() + begin_, data_.begin() + end_, 0);
  }

  void Pile::ClearInvalidRegion() {
    std::fill(data_.begin(), data_.begin() + begin_, 0);
    std::fill(data_.begin() + end_, data_.end(), 0);
  }

  void Pile::FindMedian() {
    std::vector<std::uint16_t> tmp(data_.begin() + begin_, data_.begin() + end_);
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
    median_ = tmp[tmp.size() / 2];
  }

  void Pile::FindChimericRegions() {
    auto slopes = FindSlopes(1.82);
    if (slopes.empty()) {
      return;
    }

    for (std::uint32_t i = 0; i < slopes.size() - 1; ++i) {
      if (!(slopes[i].first & 1) && (slopes[i + 1].first & 1)) {
        chimeric_regions_.emplace_back(slopes[i].first >> 1, slopes[i + 1].second);  // NOLINT
      }
    }
    chimeric_regions_ = MergeRegions(chimeric_regions_);
  }

  void Pile::ClearChimericRegions(std::uint16_t median) {
    auto is_chimeric_region = [&](const Region &r) -> bool {
      for (std::uint32_t i = r.first; i <= r.second; ++i) {
        if (clamp(data_[i] * 1.82) <= median) {
          return true;
        }
      }
      return false;
    };

    std::uint32_t begin = 0;
    std::uint32_t end = 0;
    std::uint32_t last = begin_;
    std::vector<Region> unresolved_regions;
    for (const auto &it: chimeric_regions_) {
      if (begin_ > it.first || end_ < it.second) {
        continue;
      }
      if (is_chimeric_region(it)) {
        if (it.first - last > end - begin) {
          begin = last;
          end = it.first;
        }
        last = it.second;
      } else {
        unresolved_regions.emplace_back(it);
      }
    }
    if (end_ - last > end - begin) {
      begin = last;
      end = end_;
    }

    if (begin != begin_ || end != end_) {
      set_is_chimeric();
    }
    chimeric_regions_.swap(unresolved_regions);

    UpdateValidRegion(begin, end);
  }

  void Pile::FindRepetitiveRegions(std::uint16_t median) {
    // check kmers_
    if (!kmers_.empty()) {
      std::uint32_t w = 479 >> kPSS;
      std::uint32_t group = 12;

      Region region;
      std::size_t count = 0;
      for (std::uint32_t i = 0; i < kmers_.size(); ++i) {
        if (kmers_[i] == 0) {
          continue;
        }
        if (count && i - region.second <= w) {
          region.second = i;
          ++count;
          continue;
        }
        if (count > group) {
          repetitive_regions_.emplace_back(region);
          set_is_repetitive();
        }
        region = {i, i};
        count = 1;
      }
      if (count > group) {
        repetitive_regions_.emplace_back(region);
        set_is_repetitive();
      }
    }

    // check data_
    auto slopes = FindSlopes(1.42);
    if (!slopes.empty()) {
      auto is_repetitive_region = [&](const Region &begin, const Region &end) -> bool {  // NOLINT
        if (((end.first >> 1) + end.second) / 2 -
            ((begin.first >> 1) + begin.second) / 2 > 0.84 * (end_ - begin_)) {
          return false;
        }
        bool found_peak = false;
        std::uint16_t peak_value =
            clamp(1.42 * std::max(data_[begin.second], data_[end.first >> 1]));
        std::uint16_t min_value = clamp(1.42 * median);
        std::uint32_t num_valid = 0;

        for (std::uint32_t i = begin.second + 1; i < (end.first >> 1); ++i) {
          if (data_[i] > min_value) {
            ++num_valid;
          }
          if (data_[i] > peak_value) {
            found_peak = true;
          }
        }

        if (!found_peak || num_valid < 0.9 * ((end.first >> 1) - begin.second)) {
          return false;
        }
        return true;
      };

      for (std::uint32_t i = 0; i < slopes.size() - 1; ++i) {
        if (!(slopes[i].first & 1)) {
          continue;
        }
        for (std::uint32_t j = i + 1; j < slopes.size(); ++j) {
          if (slopes[j].first & 1) {
            continue;
          }
          if (is_repetitive_region(slopes[i], slopes[j])) {
            repetitive_regions_.emplace_back(
                (slopes[i].second) - 0.336 * (slopes[i].second - (slopes[i].first >> 1)),  // NOLINT
                (slopes[j].first >> 1) + 0.336 * (slopes[j].second - (slopes[j].first >> 1)));  // NOLINT
            set_is_repetitive();
          }
        }
      }
    }

    repetitive_regions_ = MergeRegions(repetitive_regions_);
    for (auto &it: repetitive_regions_) {
      it.first = std::max(begin_, it.first) << 1;
      it.second = std::min(end_, it.second);
    }
  }

  void Pile::UpdateRepetitiveRegions(const biosoup::Overlap &o) {
    if (repetitive_regions_.empty() || (id_ != o.lhs_id && id_ != o.rhs_id)) {
      return;
    }

    std::uint32_t begin = (id_ == o.lhs_id ? o.lhs_begin : o.rhs_begin) >> kPSS;
    std::uint32_t end = (id_ == o.lhs_id ? o.lhs_end : o.rhs_end) >> kPSS;
    std::uint32_t fuzz = 420 >> kPSS;
    std::uint32_t offset = 0.1 * (end_ - begin_);

    for (auto &it: repetitive_regions_) {
      if (begin < it.second && (it.first >> 1) < end) {
        if ((it.first >> 1) < begin_ + offset && begin - begin_ < end_ - end) {
          if (end >= it.second + fuzz) {
            it.first |= 1;
          }
        } else if (it.second > end_ - offset && begin - begin_ > end_ - end) {
          if (begin + fuzz <= (it.first >> 1)) {
            it.first |= 1;
          }
        }
      }
    }
  }

  bool Pile::CheckRepetitiveRegions(const biosoup::Overlap &o) {
    if (repetitive_regions_.empty() || (id_ != o.lhs_id && id_ != o.rhs_id)) {
      return false;
    }

    std::uint32_t begin = (id_ == o.lhs_id ? o.lhs_begin : o.rhs_begin) >> kPSS;
    std::uint32_t end = (id_ == o.lhs_id ? o.lhs_end : o.rhs_end) >> kPSS;
    std::uint32_t fuzz = 420 >> kPSS;
    std::uint32_t offset = 0.1 * (end_ - begin_);

    for (const auto &it: repetitive_regions_) {
      if (begin < it.second && (it.first >> 1) < end) {
        if ((it.first >> 1) < begin_ + offset) {
          if (end < it.second + fuzz && (it.first & 1)) {
            return true;
          }
        } else if (it.second > end_ - offset) {
          if (begin + fuzz > (it.first >> 1) && (it.first & 1)) {
            return true;
          }
        }
      }
    }

    return false;
  }

  void Pile::ClearRepetitiveRegions() {
    repetitive_regions_.clear();
  }

  std::vector<Pile::Region> Pile::MergeRegions(const std::vector<Region> &src) {
    std::vector<Region> dst;
    std::vector<bool> is_merged(src.size(), 0);
    for (std::uint32_t i = 0; i < src.size(); ++i) {
      if (is_merged[i]) {
        continue;
      }
      Region r = src[i];
      while (true) {
        is_merged[i] = false;
        for (std::uint32_t j = i + 1; j < src.size(); ++j) {
          if (is_merged[j]) {
            continue;
          }
          if (r.first < src[j].second && r.second > src[j].first) {
            is_merged[i] = true;
            is_merged[j] = true;
            r.first = std::min(r.first, src[j].first);
            r.second = std::max(r.second, src[j].second);
          }
        }
        if (!is_merged[i]) {
          break;
        }
      }
      dst.emplace_back(r);
    }
    return dst;
  }

  std::vector<Pile::Region> Pile::FindSlopes(double q) {
    using Subpile = std::deque<std::pair<std::int32_t, std::uint16_t>>;
    auto subpile_add = [](Subpile &s, std::uint16_t value, std::int32_t position) -> void {  // NOLINT
      while (!s.empty() && s.back().second <= value) {
        s.pop_back();
      }
      s.emplace_back(position, value);
    };
    auto subpile_update = [](Subpile &s, std::int32_t position) {
      while (!s.empty() && s.front().first <= position) {
        s.pop_front();
      }
    };

    // find slopes
    std::vector<Region> dst;

    std::int32_t w = 847 >> kPSS;
    std::int32_t data_size = data_.size();

    Subpile left_subpile;
    std::uint32_t first_down = 0, last_down = 0;
    bool found_down = false;

    Subpile right_subpile;
    std::uint32_t first_up = 0, last_up = 0;
    bool found_up = false;

    // find slope regions
    for (std::int32_t i = 0; i < w; ++i) {
      subpile_add(right_subpile, data_[i], i);
    }
    for (std::int32_t i = 0; i < data_size; ++i) {
      if (i > 0) {
        subpile_add(left_subpile, data_[i - 1], i - 1);
      }
      subpile_update(left_subpile, i - 1 - w);

      if (i < data_size - w) {
        subpile_add(right_subpile, data_[i + w], i + w);
      }
      subpile_update(right_subpile, i);

      std::uint16_t d = clamp(data_[i] * q);
      if (i != 0 && left_subpile.front().second > d) {
        if (found_down) {
          if (i - last_down > 1) {
            dst.emplace_back(first_down << 1 | 0, last_down);
            first_down = i;
          }
        } else {
          found_down = true;
          first_down = i;
        }
        last_down = i;
      }
      if (i != (data_size - 1) && right_subpile.front().second > d) {
        if (found_up) {
          if (i - last_up > 1) {
            dst.emplace_back(first_up << 1 | 1, last_up);
            first_up = i;
          }
        } else {
          found_up = true;
          first_up = i;
        }
        last_up = i;
      }
    }
    if (found_down) {
      dst.emplace_back(first_down << 1 | 0, last_down);
    }
    if (found_up) {
      dst.emplace_back(first_up << 1 | 1, last_up);
    }
    if (dst.empty()) {
      return dst;
    }

    // separate overlaping slopes
    while (true) {
      std::sort(dst.begin(), dst.end());

      bool is_changed = false;
      for (std::uint32_t i = 0; i < dst.size() - 1; ++i) {
        if (dst[i].second < (dst[i + 1].first >> 1)) {
          continue;
        }

        if (dst[i].first & 1) {
          right_subpile.clear();
          found_up = false;
          std::uint32_t subpile_begin = dst[i].first >> 1;
          std::uint32_t subpile_end = std::min(dst[i].second, dst[i + 1].second);

          for (std::uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
            subpile_add(right_subpile, data_[j], j);
          }
          for (std::uint32_t j = subpile_begin; j < subpile_end; ++j) {
            subpile_update(right_subpile, j);
            if (clamp(data_[j] * q) < right_subpile.front().second) {
              if (found_up) {
                if (j - last_up > 1) {
                  dst.emplace_back(first_up << 1 | 1, last_up);
                  first_up = j;
                }
              } else {
                found_up = true;
                first_up = j;
              }
              last_up = j;
            }
          }
          if (found_up) {
            dst.emplace_back(first_up << 1 | 1, last_up);
          }
          dst[i].first = subpile_end << 1 | 1;

        } else {
          if (dst[i].second == (dst[i + 1].first >> 1)) {
            continue;
          }

          left_subpile.clear();
          found_down = false;

          std::uint32_t subpile_begin =
              std::max(dst[i].first >> 1, dst[i + 1].first >> 1);
          std::uint32_t subpile_end = dst[i].second;

          for (std::uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
            if (left_subpile.empty() == false &&
                clamp(data_[j] * q) < left_subpile.front().second) {
              if (found_down) {
                if (j - last_down > 1) {
                  dst.emplace_back(first_down << 1, last_down);
                  first_down = j;
                }
              } else {
                found_down = true;
                first_down = j;
              }
              last_down = j;
            }
            subpile_add(left_subpile, data_[j], j);
          }
          if (found_down) {
            dst.emplace_back(first_down << 1, last_down);
          }
          dst[i].second = subpile_begin;
        }

        is_changed = true;
        break;
      }

      if (!is_changed) {
        break;
      }
    }

    // narrow slopes
    for (std::uint32_t i = 0; i < dst.size() - 1; ++i) {
      if ((dst[i].first & 1) && !(dst[i + 1].first & 1)) {
        std::uint32_t subpile_begin = dst[i].second;
        std::uint32_t subpile_end = dst[i + 1].first >> 1;

        if (subpile_end - subpile_begin > static_cast<std::uint32_t>(w)) {
          continue;
        }

        std::uint16_t max_coverage = 0;
        for (std::uint32_t j = subpile_begin + 1; j < subpile_end; ++j) {
          max_coverage = std::max(max_coverage, data_[j]);
        }

        std::uint32_t valid_point = dst[i].first >> 1;
        for (std::uint32_t j = dst[i].first >> 1; j <= subpile_begin; ++j) {
          if (max_coverage > clamp(data_[j] * q)) {
            valid_point = j;
          }
        }
        dst[i].second = valid_point;

        valid_point = dst[i + 1].second;
        for (uint32_t j = subpile_end; j <= dst[i + 1].second; ++j) {
          if (max_coverage > clamp(data_[j] * q)) {
            valid_point = j;
            break;
          }
        }
        dst[i + 1].first = valid_point << 1 | 0;
      }
    }

    return dst;
  }

}  // namespace raven
