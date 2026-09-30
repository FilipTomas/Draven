#include "overlap_selection.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace overlap_selection {

namespace {

constexpr double kEps = 1e-8;

// Everything is computed in double, like pandas, and only rounded to float32
// at the end, like sklearn does on predict.
enum Column : std::size_t {
  // base
  kScore, kOlType, kOlLength, kOgLength, kOgToOlLengthRatio,
  kScoreToOgLengthRatio, kQHap, kTHap, kQErr, kTErr, kQRep, kTRep,
  kHapRatio, kQHor, kTHor, kJaccardIndex, kDiploidJaccardFull,
  kDiploidJaccardLonger, kDiploidJaccardShorter, kLongerOverhangLen,
  kShorterOverhangLen, kJaccradNonExtended,
  // derived
  kQErrVsTErr, kErrDifference, kQHapVsTHap, kHapDifference, kQRepVsTRep,
  kRepDifference, kOverlapFractionTarget, kOverlapFractionQuery,
  kTargetSpanFraction, kQuerySpanFraction, kTargetStartFraction,
  kTargetEndFraction, kQueryStartFraction, kQueryEndFraction,
  kTotalOverhang, kOverhangDifference,
  kNumValueColumns
};

const char* const kValueColumnNames[kNumValueColumns] = {
  "score", "ol_type", "ol_length", "og_length", "og_to_ol_length_ratio",
  "score_to_og_length_ratio", "q_hap", "t_hap", "q_err", "t_err", "q_rep",
  "t_rep", "hap_ratio", "q_hor", "t_hor", "jaccard_index",
  "diploid_jaccard_full", "diploid_jaccard_longer", "diploid_jaccard_shorter",
  "longer_overhang_len", "shorter_overhang_len", "jaccrad_non_extended",
  "q_err_vs_t_err", "err_difference", "q_hap_vs_t_hap", "hap_difference",
  "q_rep_vs_t_rep", "rep_difference", "overlap_fraction_target",
  "overlap_fraction_query", "target_span_fraction", "query_span_fraction",
  "target_start_fraction", "target_end_fraction", "query_start_fraction",
  "query_end_fraction", "total_overhang", "overhang_difference",
};

// Columns ranked within each target (notebook's `rank_features`).
const Column kRankColumns[] = {
  kScore, kOlLength, kOgLength, kOgToOlLengthRatio, kScoreToOgLengthRatio,
  kQHap, kQErr, kQRep, kHapRatio, kJaccardIndex, kDiploidJaccardFull,
  kDiploidJaccardLonger, kDiploidJaccardShorter, kJaccradNonExtended,
  kLongerOverhangLen, kShorterOverhangLen,
  kQErrVsTErr, kErrDifference, kQHapVsTHap, kHapDifference, kQRepVsTRep,
  kRepDifference, kOverlapFractionTarget, kOverlapFractionQuery,
  kTargetSpanFraction, kQuerySpanFraction, kTotalOverhang,
  kOverhangDifference,
};
constexpr std::size_t kNumRankColumns =
    sizeof(kRankColumns) / sizeof(kRankColumns[0]);
constexpr std::size_t kNumFeatures = kNumValueColumns + kNumRankColumns;

// Same values as the correction_overlaps_features.paf columns printed by
// Graph::PrintBackboneOverlapsWithFeatures (query = rhs, target = lhs).
void FillValues(const extended_overlap& o, double tl, double ql, double* v) {
  v[kScore] = o.overlap.score;
  v[kOlType] = o.graph_overlap_type;
  v[kOlLength] = o.extended_length;
  v[kOgLength] = o.found_length;
  v[kOgToOlLengthRatio] = o.found_to_extended_length;
  v[kScoreToOgLengthRatio] = o.score_to_length;
  v[kQHap] = o.rhs_hap;
  v[kTHap] = o.lhs_hap;
  v[kQErr] = o.rhs_err;
  v[kTErr] = o.lhs_err;
  v[kQRep] = o.rhs_rep;
  v[kTRep] = o.lhs_rep;
  v[kHapRatio] = o.hap_ratio;
  // the PAF writer prints t_hor in the q_hor column and vice versa
  v[kQHor] = o.t_hor ? 1 : 0;
  v[kTHor] = o.q_hor ? 1 : 0;
  v[kJaccardIndex] = o.jaccard_index;
  v[kDiploidJaccardFull] = o.diploid_jaccard_full_overlap;
  v[kDiploidJaccardLonger] = o.diploid_jaccard_longer;
  v[kDiploidJaccardShorter] = o.diploid_jaccard_shorter;
  v[kLongerOverhangLen] = o.longer_overhang_len;
  v[kShorterOverhangLen] = o.shorter_overhang_len;
  v[kJaccradNonExtended] = o.jaccard_non_extended_overlap;

  const double qs = o.overlap.rhs_begin, qe = o.overlap.rhs_end;
  const double ts = o.overlap.lhs_begin, te = o.overlap.lhs_end;

  v[kQErrVsTErr] = (v[kQErr] + kEps) / (v[kTErr] + kEps);
  v[kErrDifference] = std::abs(v[kQErr] - v[kTErr]);
  v[kQHapVsTHap] = (v[kQHap] + kEps) / (v[kTHap] + kEps);
  v[kHapDifference] = std::abs(v[kQHap] - v[kTHap]);
  v[kQRepVsTRep] = (v[kQRep] + kEps) / (v[kTRep] + kEps);
  v[kRepDifference] = std::abs(v[kQRep] - v[kTRep]);

  v[kOverlapFractionTarget] = v[kOlLength] / tl;
  v[kOverlapFractionQuery] = v[kOlLength] / ql;
  v[kTargetSpanFraction] = (te - ts) / tl;
  v[kQuerySpanFraction] = (qe - qs) / ql;

  v[kTargetStartFraction] = ts / tl;
  v[kTargetEndFraction] = te / tl;
  v[kQueryStartFraction] = qs / ql;
  v[kQueryEndFraction] = qe / ql;

  v[kTotalOverhang] = v[kLongerOverhangLen] + v[kShorterOverhangLen];
  v[kOverhangDifference] =
      std::abs(v[kLongerOverhangLen] - v[kShorterOverhangLen]);
}

// pandas Series.rank(method="average", pct=True): NaNs stay NaN and are not
// counted in the denominator.
void PercentRank(const std::vector<double>& x, std::vector<double>& out) {
  out.assign(x.size(), std::numeric_limits<double>::quiet_NaN());

  std::vector<std::size_t> order;
  order.reserve(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!std::isnan(x[i])) order.push_back(i);
  }
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });

  const double n = static_cast<double>(order.size());
  for (std::size_t i = 0; i < order.size();) {
    std::size_t j = i + 1;
    while (j < order.size() && x[order[j]] == x[order[i]]) ++j;
    // 1-based ranks i+1 .. j share their average
    const double rank = (static_cast<double>(i + 1 + j) / 2.0) / n;
    for (std::size_t k = i; k < j; ++k) out[order[k]] = rank;
    i = j;
  }
}

float ToFeature(double v) {
  return std::isnan(v) ? 0.0f : static_cast<float>(v);
}

}  // namespace

const std::vector<std::string>& FeatureNames() {
  static const std::vector<std::string> names = [] {
    std::vector<std::string> n(kValueColumnNames,
                               kValueColumnNames + kNumValueColumns);
    for (Column c : kRankColumns) {
      n.emplace_back(std::string(kValueColumnNames[c]) + "_rank_t");
    }
    return n;
  }();
  return names;
}

std::vector<float> ComputeFeatures(
    const std::vector<extended_overlap>& pile,
    std::uint32_t target_length,
    const std::vector<std::uint32_t>& query_lengths) {
  const std::size_t n = pile.size();

  std::vector<double> values(n * kNumValueColumns);
  for (std::size_t i = 0; i < n; ++i) {
    FillValues(pile[i], target_length, query_lengths[i],
               &values[i * kNumValueColumns]);
  }

  std::vector<float> features(n * kNumFeatures);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t c = 0; c < kNumValueColumns; ++c) {
      features[i * kNumFeatures + c] = ToFeature(values[i * kNumValueColumns + c]);
    }
  }

  std::vector<double> column(n), ranks;
  for (std::size_t r = 0; r < kNumRankColumns; ++r) {
    for (std::size_t i = 0; i < n; ++i) {
      column[i] = values[i * kNumValueColumns + kRankColumns[r]];
    }
    PercentRank(column, ranks);
    for (std::size_t i = 0; i < n; ++i) {
      features[i * kNumFeatures + kNumValueColumns + r] = ToFeature(ranks[i]);
    }
  }

  return features;
}

void RandomForest::Load(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open()) {
    throw std::runtime_error("[raven::RandomForest] cannot open " + path);
  }
  auto fail = [&](const std::string& what) {
    throw std::runtime_error("[raven::RandomForest] " + path + ": " + what);
  };

  std::string tag;
  int version = 0;
  if (!(in >> tag >> version) || tag != "DRAVEN_RF" || version != 1) {
    fail("not a DRAVEN_RF 1 file (export with scripts/export_rf_model.py)");
  }

  std::size_t n_features = 0;
  if (!(in >> tag >> n_features) || tag != "n_features") fail("bad header");

  const auto& expected = FeatureNames();
  if (n_features != expected.size()) {
    fail("model has " + std::to_string(n_features) + " features, draven "
         "computes " + std::to_string(expected.size()));
  }
  for (std::size_t f = 0; f < n_features; ++f) {
    std::string name;
    if (!(in >> name)) fail("truncated feature list");
    if (name != expected[f]) {
      fail("feature " + std::to_string(f) + " is '" + name +
           "', draven expects '" + expected[f] + "'");
    }
  }

  std::size_t n_trees = 0;
  if (!(in >> tag >> n_trees) || tag != "n_trees") fail("bad tree count");

  trees_.assign(n_trees, {});
  for (auto& tree : trees_) {
    std::size_t n_nodes = 0;
    if (!(in >> tag >> n_nodes) || tag != "tree" || n_nodes == 0) {
      fail("bad tree header");
    }
    tree.resize(n_nodes);
    for (auto& node : tree) {
      if (!(in >> node.left >> node.right >> node.feature >> node.threshold >>
            node.p_positive)) {
        fail("truncated tree");
      }
      const bool leaf = node.left == -1;
      if (!leaf && (node.left < 0 || node.right < 0 ||
                    static_cast<std::size_t>(node.left) >= n_nodes ||
                    static_cast<std::size_t>(node.right) >= n_nodes ||
                    node.feature < 0 ||
                    static_cast<std::size_t>(node.feature) >= n_features)) {
        fail("node index out of range");
      }
    }
  }
}

double RandomForest::PredictProba(const float* x) const {
  double sum = 0.0;
  for (const auto& tree : trees_) {
    const Node* node = &tree[0];
    while (node->left != -1) {
      // sklearn compares the float32 feature against a double threshold
      node = &tree[static_cast<double>(x[node->feature]) <= node->threshold
                       ? node->left : node->right];
    }
    sum += node->p_positive;
  }
  return sum / static_cast<double>(trees_.size());
}

std::vector<bool> SelectWindowTopK(
    const std::vector<extended_overlap>& pile,
    const std::vector<double>& probabilities,
    std::uint32_t target_length,
    std::uint32_t window_size,
    std::uint32_t top_k,
    bool require_full_window) {
  std::vector<bool> selected(pile.size(), false);
  if (window_size == 0 || top_k == 0) return selected;

  const std::uint64_t n_windows =
      (static_cast<std::uint64_t>(target_length) + window_size - 1) / window_size;

  std::vector<std::size_t> eligible;
  for (std::uint64_t w = 0; w < n_windows; ++w) {
    const std::uint64_t ws = w * window_size;
    const std::uint64_t we =
        std::min<std::uint64_t>((w + 1) * window_size, target_length);

    eligible.clear();
    for (std::size_t i = 0; i < pile.size(); ++i) {
      const std::uint64_t ts = pile[i].overlap.lhs_begin;
      const std::uint64_t te = pile[i].overlap.lhs_end;
      const bool ok = require_full_window ? (ts <= ws && te >= we)
                                          : (ts < we && te > ws);
      if (ok) eligible.push_back(i);
    }

    if (eligible.size() > top_k) {
      std::nth_element(eligible.begin(), eligible.begin() + (top_k - 1),
                       eligible.end(), [&](std::size_t a, std::size_t b) {
                         return probabilities[a] > probabilities[b];
                       });
      eligible.resize(top_k);
    }
    for (std::size_t i : eligible) selected[i] = true;
  }

  return selected;
}

}  // namespace overlap_selection
