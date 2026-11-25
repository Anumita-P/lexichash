"""merge_prediction_files in LexicHash_1/LexicHash/run_module.py: the vectorised version must give the
same rows as the old row-by-row loop (resolve-list item 10)."""
import ast
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RM = os.path.join(ROOT, "LexicHash_1", "LexicHash", "run_module.py")

OLD = '''
def merge_prediction_files(paths):
    merged = {}
    for path in paths:
        df = pd.read_csv(path, sep="\\t", header=None, names=["i1", "i2", "score", "direction"])
        for _, row in df.iterrows():
            key = (int(row.i1), int(row.i2), row.direction)
            score = int(row.score)
            if key not in merged or score > merged[key]:
                merged[key] = score
    merged_df = pd.DataFrame([(i1, i2, score, d) for (i1, i2, d), score in merged.items()],
                             columns=["i1", "i2", "score", "direction"])
    return merged_df
'''


def get(source):
    ns = {"pd": pd}
    for node in ast.parse(source).body:
        if isinstance(node, ast.FunctionDef) and node.name == "merge_prediction_files":
            exec(compile(ast.Module([node], []), "x", "exec"), ns)
    return ns["merge_prediction_files"]


def main():
    old = get(OLD)
    new = get(open(RM, encoding="utf-8").read())
    rng = np.random.default_rng(0)
    with tempfile.TemporaryDirectory() as d:
        paths = []
        for s in range(4):
            n = 3000
            df = pd.DataFrame({"i1": rng.integers(1, 60, n), "i2": rng.integers(60, 120, n),
                               "score": rng.integers(13, 129, n), "d": rng.choice(["+", "-"], n)})
            p = os.path.join(d, "s%d.tsv" % s)
            df.to_csv(p, sep="\t", header=False, index=False)
            paths.append(p)
        a, b = old(paths), new(paths)
    key = ["i1", "i2", "direction"]
    a = a.sort_values(key).reset_index(drop=True)
    b = b.sort_values(key).reset_index(drop=True)
    print("rows: old %d  new %d" % (len(a), len(b)))
    assert list(a.columns) == list(b.columns)
    assert a.equals(b), "merged tables differ"
    print("OK")


if __name__ == "__main__":
    main()
