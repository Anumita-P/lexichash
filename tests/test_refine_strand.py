#!/usr/bin/env python3
"""refine_match_pairs must score a '-' pair on the read / reverse-complement rows (batch 3 fix).

Builds a tiny sketch matrix by hand (3 reads + their reverse complements, 4 masks) in which read 0
and the reverse complement of read 1 carry identical min-hash k-mers, i.e. a perfect '-' overlap
between reads 0 and 1. The score of the pair (0, 1, '-') must then be the maximum, and must not depend
on the unrelated forward rows 0 and 1. Run from the repository root:
    python3 tests/test_refine_strand.py
"""
import os
import sys

sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "experiment with lexichash", "LexicHash-main"))
import numpy as np  # noqa: E402
import LexicHash.LexicHash as LH  # noqa: E402

LH.MAX_K, LH.MIN_K = 32, 14
n, H = 3, 4
rng = np.random.RandomState(0)
S = rng.randint(0, 2 ** 63, size=(2 * n, H), dtype=np.int64).astype(np.uint64)
S[0] = S[n + 1]                      # read 0 == revcomp of read 1, in every mask -> '-' overlap
S[1] = S[0] ^ np.uint64(1 << 63)     # forward rows 0 and 1 differ in the first base only

pairs = {(0, 1, "-"): 32, (0, 1, "+"): 0}
df = LH.refine_match_pairs(pairs, S, decay=0.1, top_k=100, n_seq=n)
score = {(r.id1, r.id2, r.is_rc): r.score for r in df.itertuples()}
full = sum(32 * 0.9 ** i for i in range(H))
ok_minus = abs(score[(0, 1, "-")] - full) < 1e-9
ok_order = score[(0, 1, "-")] > score[(0, 1, "+")]
old = LH.refine_match_pairs(pairs, S, decay=0.1, top_k=100)  # n_seq not given: old behaviour
print("'-' pair score with strand-aware rows: %.3f (perfect = %.3f)" % (score[(0, 1, "-")], full))
print("'-' pair score on forward rows only  : %.3f" % old[(old.id2 == 1) & (old.is_rc == "-")].score.iloc[0])
if ok_minus and ok_order:
    print("PASS")
else:
    print("FAIL")
    sys.exit(1)
