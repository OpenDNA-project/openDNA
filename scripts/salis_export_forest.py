"""Turn the Synthesis Success Calculator's trained forest into app/salis_forest.json.gz.

The authors' model is a pickled scikit-learn 0.20.4 RandomForestClassifier (SSC_forest.pkl in
https://github.com/hsalis/SalisLabCode, Synthesis_Success_Calculator/, GPLv3). A pickle runs code when it is
loaded, so this script never lets it: only the few classes the file needs are rebuilt, as inert stand-ins,
and anything else raises. It needs numpy; scikit-learn and pandas are NOT needed.

    python scripts/salis_export_forest.py path/to/SSC_forest.pkl app/salis_forest.json.gz

The output keeps, per tree and node, what prediction needs: left, right, feature, threshold and, for each class,
the share of the node (p0 "Canceled", p1 "Synthesized"), value / sum(value) as scikit-learn's predict_proba computes it.
"""
import gzip
import json
import pickle
import sys



class Stub:
    def __init__(self, *a, **k):
        self.args = a

    def __setstate__(self, st):
        self.state = st


ALLOWED = {
    ("sklearn.ensemble.forest", "RandomForestClassifier"): Stub,
    ("sklearn.tree.tree", "DecisionTreeClassifier"): Stub,
    ("sklearn.tree._tree", "Tree"): Stub,
    ("pandas.core.indexes.base", "_new_Index"): lambda cls, d: d,
    ("pandas.core.indexes.base", "Index"): Stub,
}


class Reader(pickle.Unpickler):
    def find_class(self, mod, name):
        if (mod, name) in ALLOWED:
            return ALLOWED[(mod, name)]
        if mod.split(".")[0] == "numpy" and name in ("_reconstruct", "ndarray", "dtype", "scalar"):
            return super().find_class(mod, name)
        if (mod, name) in (("collections", "OrderedDict"), ("builtins", "object"), ("__builtin__", "object"), ("copy_reg", "_reconstructor")):
            return super().find_class(mod, name)
        raise pickle.UnpicklingError(f"blocked global {mod}.{name}")


def main(src, dst):
    with open(src, "rb") as f:
        model = Reader(f, encoding="latin1").load()
    state = model["forest"].state
    classes = [str(c) for c in state["classes_"]]
    assert classes == ["Canceled", "Synthesized"], classes
    trees = []
    for est in state["estimators_"]:
        t = est.state["tree_"].state
        nodes, values = t["nodes"], t["values"][:, 0, :]
        total = values.sum(axis=1)
        total[total == 0.0] = 1.0
        trees.append({
            "left": nodes["left_child"].tolist(), "right": nodes["right_child"].tolist(),
            "feature": nodes["feature"].tolist(), "threshold": nodes["threshold"].tolist(),
            "p0": (values[:, 0] / total).tolist(), "p1": (values[:, 1] / total).tolist(),
        })
    out = {"sklearn": str(state["_sklearn_version"]), "features": [str(x) for x in model["features"]["data"]], "trees": trees}
    with gzip.open(dst, "wt", encoding="ascii", compresslevel=9) as f:
        json.dump(out, f, separators=(",", ":"))
    print(len(trees), "trees,", sum(len(t["left"]) for t in trees), "nodes ->", dst)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
