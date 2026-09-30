#!/usr/bin/env python3
"""Export a trained sklearn RandomForestClassifier to draven's text format.

Usage from the notebook (after `model.fit(...)`):

    import sys; sys.path.append("/home/ftomas/github/draven/scripts")
    from export_rf_model import export_rf
    export_rf(model, feature_cols, "correction_rf.txt")

or from the shell with a joblib/pickle dump:

    python export_rf_model.py model.joblib correction_rf.txt

draven checks that the feature names in the file match (in order) the
features it computes, see src/graph_constructor/overlap_selection.cpp.

Format (whitespace separated):

    DRAVEN_RF 1
    n_features <F>
    <feature name>            x F
    n_trees <T>
    tree <n_nodes>            x T, each followed by n_nodes lines:
    <left> <right> <feature> <threshold> <p_positive>

Leaves have left == right == -1 and p_positive is the leaf's
(class-weighted) fraction of the positive class, which is exactly what
sklearn's predict_proba averages over the trees. Internal nodes go left
when float32(x[feature]) <= threshold, as in sklearn.
"""

import sys

import numpy as np


def export_rf(model, feature_cols, path):
    feature_cols = list(feature_cols)

    if model.n_features_in_ != len(feature_cols):
        raise ValueError(
            f"model has {model.n_features_in_} features, "
            f"feature_cols has {len(feature_cols)}")

    if hasattr(model, "feature_names_in_") and \
            list(model.feature_names_in_) != feature_cols:
        raise ValueError("feature_cols order differs from the training order")

    classes = list(model.classes_)
    if len(classes) != 2:
        raise ValueError(f"expected a binary classifier, got classes {classes}")
    pos = classes.index(True) if True in classes else 1

    with open(path, "w") as out:
        out.write("DRAVEN_RF 1\n")
        out.write(f"n_features {len(feature_cols)}\n")
        for name in feature_cols:
            out.write(f"{name}\n")

        out.write(f"n_trees {len(model.estimators_)}\n")
        for est in model.estimators_:
            tree = est.tree_
            value = tree.value[:, 0, :].astype(np.float64)
            total = value.sum(axis=1)
            total[total == 0] = 1.0
            p_pos = value[:, pos] / total

            out.write(f"tree {tree.node_count}\n")
            for n in range(tree.node_count):
                left = int(tree.children_left[n])
                right = int(tree.children_right[n])
                if left == -1:
                    out.write(f"-1 -1 -1 0 {p_pos[n]!r}\n")
                else:
                    out.write(f"{left} {right} {int(tree.feature[n])} "
                              f"{float(tree.threshold[n])!r} 0\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: export_rf_model.py <model.joblib> <out.txt>")

    import joblib

    model = joblib.load(sys.argv[1])
    if not hasattr(model, "feature_names_in_"):
        sys.exit("model was not fit on a DataFrame; call export_rf() with "
                 "the feature_cols list from the notebook instead")
    export_rf(model, model.feature_names_in_, sys.argv[2])
    print(f"wrote {sys.argv[2]}")
