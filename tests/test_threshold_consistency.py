#!/usr/bin/env python3
"""Check that every evaluation script uses the single ground-truth threshold from
evaluation/lh_eval.py (THETA = 0.2, comparison ">=").

Before this check the repository mixed 0 (gengraphs.py), 0.3 with ">" (default of get_gt_df in
metric_utils.py) and 0.2 with ">=" (the plot scripts), so ROC/PRC numbers from different scripts
were not comparable.

Run from the repository root:  python3 tests/test_threshold_consistency.py
"""
import glob
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "evaluation"))
from lh_eval import THETA  # noqa: E402

problems = []


def read(path):
    with open(os.path.join(REPO, path), encoding="utf-8", errors="replace") as fh:
        return fh.read()


def check(cond, msg):
    if not cond:
        problems.append(msg)


check(THETA == 0.2, "THETA in evaluation/lh_eval.py is %r, expected 0.2" % THETA)

# 1. metric_utils.py copies: default 0.2 and ">="
for path in glob.glob(os.path.join(REPO, "**", "metric_utils.py"), recursive=True):
    rel = os.path.relpath(path, REPO)
    s = read(rel)
    defaults = re.findall(r"def get_gt_df\w*\([^)]*min_theta=([0-9.]+)", s)
    check(defaults and all(float(d) == THETA for d in defaults), "%s: min_theta default is %s" % (rel, defaults))
    check("thetas > min_theta" not in s, "%s: strict '>' comparison with min_theta" % rel)
    check(s.count("thetas >= min_theta") == 2, "%s: expected two 'thetas >= min_theta'" % rel)

# 2. no call passes a different threshold explicitly
call = re.compile(r"get_gt_df(?:_nctc)?\(\s*[^,()]+,\s*[^,()]+,\s*([^)]+)\)")
for path in glob.glob(os.path.join(REPO, "**", "*.py"), recursive=True):
    rel = os.path.relpath(path, REPO)
    for line in read(rel).splitlines():
        if line.lstrip().startswith(("def ", "#")):  # definitions and comments are not calls
            continue
        for m in call.finditer(line):
            check(False, "%s: get_gt_df called with an explicit threshold (%s)" % (rel, m.group(1).strip()))

# 3. plot scripts take the threshold from lh_eval and have no 0.2 literal of their own
for rel in ("test accuracy/roc_prc/plot.py", "test accuracy/roc_prc/plot copy.py",
            "comparisons/expwala/plot copy.py", "comparisons/multipleruns/plot copy.py"):
    s = read(rel)
    check("from lh_eval import THETA" in s, "%s: does not import THETA" % rel)
    check(not re.search(r"jaccard\s*>=\s*0\.\d|jaccard_threshold\s*=\s*0\.\d", s), "%s: has its own threshold literal" % rel)

if problems:
    print("FAILED:")
    for p in problems:
        print("  - " + p)
    sys.exit(1)
print("PASS: one threshold (theta >= %g) in all evaluation scripts" % THETA)
