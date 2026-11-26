"""Check evaluation/multiscale_eval.py against brute force (resolve-list item 10).

1. lexicographic_first() returns the substring that is smallest when every base b at position j is
   ranked by (b - mask[j]) mod 4 (brute force over all start positions, end of read padded with the
   lexicographically last base).
2. y_all_pairs() equals the common-prefix length computed base by base.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "evaluation"))
import multiscale_eval as M


def brute_first(seq, mask, K, MIN_K):
    best, best_key = None, None
    for i in range(len(seq) - MIN_K):
        w = seq[i:i + K]
        digits = [(int(w[j]) - int(mask[j])) % 4 for j in range(len(w))]
        digits += [3] * (K - len(w))
        key = tuple(digits)
        if best_key is None or key < best_key:
            best_key, best = key, i
    w = seq[best:best + K]
    out = list(w) + [int(mask[j]) ^ 3 for j in range(len(w), K)]
    return np.array(out, dtype=np.uint8)


def main():
    rng = np.random.default_rng(0)
    bad = 0
    n_tests = 0
    for K, MIN_K in [(16, 13), (32, 14), (64, 31), (128, 63)]:
        for _ in range(40):
            L = int(rng.integers(150, 700))
            seq = rng.integers(0, 4, size=L).astype(np.uint8)
            mask = rng.integers(0, 4, size=K).astype(np.uint8)
            if rng.random() < 0.5:                       # make the mask match part of the read
                st = int(rng.integers(0, L - K))
                m0 = min(20, K); mask[:m0] = seq[st:st + m0]
            got = M.lexicographic_first(seq, mask, K, MIN_K)
            want = brute_first(seq, mask, K, MIN_K)
            n_tests += 1
            if not np.array_equal(got, want):
                bad += 1
    print("lexicographic_first: %d of %d cases differ from brute force" % (bad, n_tests))
    assert bad == 0

    for K in (16, 32, 64, 128):
        n, H = 12, 7
        S = rng.integers(0, 4, size=(n, H, K)).astype(np.uint8)
        for i in range(n):                                 # create shared prefixes of random length
            for j in range(i):
                if rng.random() < 0.6:
                    m = int(rng.integers(0, H)); L = int(rng.integers(1, K + 1))
                    S[i, m, :L] = S[j, m, :L]
        y, iu, ju = M.y_all_pairs(M.pack(S), K, chunk=5)
        ref = np.zeros_like(y)
        for r, (a, b) in enumerate(zip(iu, ju)):
            for h in range(H):
                d = np.nonzero(S[a, h] != S[b, h])[0]
                ref[r, h] = K if len(d) == 0 else d[0]
        print("y_all_pairs K=%d: %d of %d values differ from base-by-base prefix" % (K, int((y != ref).sum()), y.size))
        assert np.array_equal(y, ref)
    print("OK")


if __name__ == "__main__":
    main()
