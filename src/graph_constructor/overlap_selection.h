#ifndef RAVEN_OVERLAP_SELECTION_H
#define RAVEN_OVERLAP_SELECTION_H

#include <cstdint>
#include <string>
#include <vector>

#include "extended_overlap.h"

// In-process port of the correction-overlap selection from
// Herro_assembly_tools/jupyter_notebooks/mm2_draven_match.ipynb:
//
//   1. per-overlap features (base + derived + per-target ranks)
//   2. random forest trained in sklearn, exported with
//      scripts/export_rf_model.py
//   3. inference -> P(useful overlap)
//   4. per-target windowing: keep the top-K overlaps of every window,
//      final set = union across windows
//
// A "pile" below is extended_overlaps[i] of one target (lhs) read, the unit
// the notebook groups by (tn).

namespace overlap_selection {

// Feature names in model column order (notebook's `feature_cols`).
const std::vector<std::string>& FeatureNames();

// Row-major [pile.size() x FeatureNames().size()] feature matrix, with the
// notebook's fillna(0) applied. query_lengths[k] is the inflated length of
// pile[k].overlap.rhs_id.
std::vector<float> ComputeFeatures(
    const std::vector<extended_overlap>& pile,
    std::uint32_t target_length,
    const std::vector<std::uint32_t>& query_lengths);

class RandomForest {
 public:
  // Throws std::runtime_error on a malformed file or when the feature names
  // stored in the file differ from FeatureNames().
  void Load(const std::string& path);

  // sklearn predict_proba(x)[:, 1]
  double PredictProba(const float* x) const;

  bool empty() const { return trees_.empty(); }
  std::size_t num_trees() const { return trees_.size(); }

 private:
  struct Node {
    std::int32_t left;
    std::int32_t right;
    std::int32_t feature;
    double threshold;
    double p_positive;
  };

  std::vector<std::vector<Node>> trees_;
};

// Notebook's select_window_topk for one target: split [0, target_length)
// into window_size windows, keep the top_k highest-probability overlaps of
// each window. With require_full_window an overlap is eligible only if it
// spans the whole window, otherwise any overlap touching it is.
std::vector<bool> SelectWindowTopK(
    const std::vector<extended_overlap>& pile,
    const std::vector<double>& probabilities,
    std::uint32_t target_length,
    std::uint32_t window_size,
    std::uint32_t top_k,
    bool require_full_window);

}  // namespace overlap_selection

#endif  // RAVEN_OVERLAP_SELECTION_H
