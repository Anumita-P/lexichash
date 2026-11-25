#!/usr/bin/env python3
"""Scorer ablation: does the exponential-decay score beat simpler scores?  (resolve-list item 3)

For every read pair we have, per mask m, the shared-prefix length y_m of the two min-hash
k-mers. A "scorer" turns the vector (y_1 .. y_H) into one number:

  maxy            original LexicHash: max_m y_m
  expdecay        sort y descending, sum_i y_(i) * (1 - decay)^i over the top_k values   (ours)
  meantop         mean of the top_k values of y
  meanall         mean of all H values
  count_ge        number of masks with y_m >= t
  lik             sum_m 4^(max(y_m - min_k, 0))   (each extra shared base is 4x less likely by chance)

Two candidate pools are evaluated:
  all_pairs       every pair of reads is scored (upper bound for a scorer, not a pipeline)
  pipeline_aN     only the alpha*n pairs the prefix-tree step returns (alpha = N); other pairs
                  get score 0. This is how the real pipeline uses the scorer (as a re-ranker).

Three settings for reverse complements and for choosing which sketches go into y:
  no_rc               C++ run with --no_rc; y from the forward rows (a, b) only
  rc_strand_aware     reverse complements on; '+' pair (a, b): rows (a, b) and (a', b'); '-' pair:
                      rows (a, b') and (a', b) (a' = reverse complement of a); per mask the larger y
  rc_fwd_rows_only    reverse complements on, but y from the forward rows (a, b) for every pair, also
                      for '-' pairs. This is what refine_match_pairs() in "experiment with lexichash"
                      does today.
NCTC1080 is effectively single-strand (all ground-truth lines are '+', and for 97.6% of the true
pairs the same-strand sketches give the longer match), so the reverse-complement settings mainly
add chance matches there. A dataset with mixed strands is needed to judge them properly.

Sketches and candidate pairs come from the C++ program (identical to Python, see
tests/equivalence_py_vs_cpp.py). Metrics use evaluation/lh_eval.py (theta >= 0.2).

Usage (from the repository root):
  g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native -lz
  python3 evaluation/scorer_ablation.py --cpp-bin ./lexichash \
      --fasta LexicHash-main/data/NCTC1080/NCTC1080_reads.fasta.gz \
      --gt LexicHash-main/data/NCTC1080/NCTC1080_daligner_ground_truth.txt \
      --seeds 1 2 3 4 5 --out evaluation/results/scorer_ablation_NCTC1080
"""
import argparse
import itertools
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lh_eval as E  # noqa: E402

MAX_K = 32
MIN_K = 14
DECAYS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.9]
TOP_KS = [10, 25, 50, 100, None]  # None = all masks
MEAN_KS = [10, 50, 100, None]
COUNT_TS = [14, 16, 18, 20, 22, 24]


# ----------------------------------------------------------------------------- sketches
def run_cpp(binary, fasta, workdir, n_hash, seed, alpha, mask_file, tag, rc=True):
    out = os.path.join(workdir, tag)
    os.makedirs(out, exist_ok=True)
    cmd = [binary, "--out", out, "--fasta", fasta, "--n_hash", str(n_hash), "--alpha", str(alpha),
           "--n_cpu", str(min(8, os.cpu_count() or 1)), "--max_k", str(MAX_K), "--min_k", str(MIN_K),
           "--seed", str(seed), "--mask_file", mask_file]
    if not rc:
        cmd.append("--no_rc")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("C++ run failed:\n" + r.stdout[-400:] + r.stderr[-400:])
    return out


def read_sketches(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            rows.append(np.array(line.split(":", 1)[1].split(), dtype=np.uint64))
    return np.vstack(rows)


def read_pairs(path):
    """pairwise_results.txt -> list of (a, b, strand) with 0-based a < b."""
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 4:
                a, b = int(p[0]) - 1, int(p[1]) - 1
                out.append((a, b, p[3]))
    return out


# ----------------------------------------------------------------------------- y vectors
def prefix_len(x):
    """Shared prefix length (in bases, 0..32) from the XOR of two packed 32-mers."""
    hi = (x >> np.uint64(32)).astype(np.float64)
    lo = (x & np.uint64(0xFFFFFFFF)).astype(np.float64)
    e_hi = np.frexp(hi)[1]
    e_lo = np.frexp(lo)[1]
    bit_len = np.where(hi > 0, 32 + e_hi, e_lo)
    return ((64 - bit_len) // 2).astype(np.uint8)


def y_matrix(S, n, a, b, plus, mode):
    """y vectors (len(a), H) for pairs (a[i], b[i]); plus[i] True for '+' strand."""
    if mode in ("no_rc", "rc_fwd_rows_only"):
        return prefix_len(S[a] ^ S[b])
    y_same = np.maximum(prefix_len(S[a] ^ S[b]), prefix_len(S[a + n] ^ S[b + n]))
    y_opp = np.maximum(prefix_len(S[a] ^ S[b + n]), prefix_len(S[a + n] ^ S[b]))
    return np.where(plus[:, None], y_same, y_opp)


# ----------------------------------------------------------------------------- scorers
def scorer_names(H):
    names = [("maxy", "")]
    for d, k in itertools.product(DECAYS, TOP_KS):
        names.append(("expdecay", "decay=%g,top_k=%s" % (d, k or "all")))
    for k in MEAN_KS:
        names.append(("meantop", "top_k=%s" % (k or "all")))
    for t in COUNT_TS:
        names.append(("count_ge", "t=%d" % t))
    names.append(("lik", ""))
    return names


def apply_scorers(y, H):
    """dict {(scorer, params): scores array} for a y matrix (m, H)."""
    ys = np.sort(y, axis=1)[:, ::-1].astype(np.float64)  # descending
    res = {("maxy", ""): ys[:, 0]}
    cum = np.cumsum(ys, axis=1)
    for d, k in itertools.product(DECAYS, TOP_KS):
        kk = min(k or H, H)
        w = (1.0 - d) ** np.arange(kk)
        res[("expdecay", "decay=%g,top_k=%s" % (d, k or "all"))] = ys[:, :kk] @ w
    for k in MEAN_KS:
        kk = min(k or H, H)
        res[("meantop", "top_k=%s" % (k or "all"))] = cum[:, kk - 1] / kk
    for t in COUNT_TS:
        res[("count_ge", "t=%d" % t)] = (y >= t).sum(axis=1).astype(np.float64)
    res[("lik", "")] = np.power(4.0, np.clip(y.astype(np.int16) - MIN_K, 0, None)).sum(axis=1)
    return res


# ----------------------------------------------------------------------------- pools
def score_all_pairs(S, n, mode, H, chunk=20000):
    iu, ju = np.triu_indices(n, k=1)
    out = {}
    for start in range(0, len(iu), chunk):
        a, b = iu[start:start + chunk], ju[start:start + chunk]
        flags = (True,) if mode in ("no_rc", "rc_fwd_rows_only") else (True, False)
        per_strand = [apply_scorers(y_matrix(S, n, a, b, np.full(len(a), f), mode), H) for f in flags]
        for key in per_strand[0]:
            v = per_strand[0][key]
            for other in per_strand[1:]:
                v = np.maximum(v, other[key])
            out.setdefault(key, []).append(v.astype(np.float32))
    scores = {k: np.concatenate(v) for k, v in out.items()}
    return scores, iu, ju


def score_candidates(S, n, cands, mode, H):
    a = np.array([c[0] for c in cands])
    b = np.array([c[1] for c in cands])
    plus = np.array([c[2] == "+" for c in cands])
    scored = apply_scorers(y_matrix(S, n, a, b, plus, mode), H)
    return scored, a, b


def to_full_array(scores_arr, a, b, n):
    """Scores of candidate pairs -> array over all pairs (np.triu_indices order); 0 elsewhere.
    A pair that appears with both strands keeps the larger score."""
    full = np.zeros(n * (n - 1) // 2, dtype=np.float64)
    np.maximum.at(full, E.pair_index(a, b, n), scores_arr)
    return full


def load_gt_top_idx(path, n, k):
    """Pair indices of the k longest ground-truth alignments (no threshold)."""
    best = {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 3:
                i, j, al = int(p[0]), int(p[1]), float(p[2])
                key = (i, j) if i < j else (j, i)
                best[key] = max(best.get(key, 0.0), al)
    top = sorted(best.items(), key=lambda x: -x[1])[:k]
    return {E.pair_index(i, j, n) for (i, j), _ in top}


def evaluate_scores(full, y_true, gt_top, k=1000):
    m = E.evaluate_arrays(full, y_true)
    top_pred = set(E.top_indices(full, k).tolist())
    m["top1000_overlap"] = 100.0 * len(top_pred & gt_top) / k
    return m


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpp-bin", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--n-hash", type=int, default=500)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--alphas", type=float, nargs="+", default=[1, 2, 5, 10])
    ap.add_argument("--rc-seeds", type=int, nargs="*", default=[1],
                    help="seeds for which the (about 6x slower) reverse-complement settings are also run")
    ap.add_argument("--rc-alphas", type=float, nargs="*", default=None,
                    help="alphas used in the reverse-complement runs (default: same as --alphas)")
    ap.add_argument("--out", required=True, help="output prefix (writes <out>.csv and <out>_summary.csv)")
    ap.add_argument("--theta", type=float, default=E.THETA)
    ap.add_argument("--no-all-pairs", action="store_true", help="skip the all-pairs pool (faster)")
    args = ap.parse_args()
    binary = os.path.abspath(args.cpp_bin)

    lens = E.read_lengths(args.fasta)
    n = len(lens)
    true_pairs = E.load_true_pairs(args.gt, lens, args.theta)
    y_true = E.true_pair_labels(true_pairs, n)
    gt_top = load_gt_top_idx(args.gt, n, 1000)
    print("reads: %d   true pairs (theta >= %g): %d of %d" % (n, args.theta, len(true_pairs), n * (n - 1) // 2))

    rows = []

    def save():
        df = pd.DataFrame(rows)
        df.to_csv(args.out + ".csv", index=False)
        return df

    with tempfile.TemporaryDirectory() as wd:
        for seed in args.seeds:
            t0 = time.time()
            for rc in ((False, True) if seed in args.rc_seeds else (False,)):
                mask_file = os.path.join(wd, "masks_%d.txt" % seed)
                outs = {}
                alphas = args.alphas if (not rc or args.rc_alphas is None) else args.rc_alphas
                for alpha in alphas:
                    outs[alpha] = run_cpp(binary, args.fasta, wd, args.n_hash, seed, alpha, mask_file,
                                          "s%d_a%g_rc%d" % (seed, alpha, rc), rc=rc)
                S = read_sketches(os.path.join(outs[alphas[0]], "sketches.txt"))
                assert S.shape == ((2 if rc else 1) * n, args.n_hash), S.shape
                H = S.shape[1]
                modes = ("rc_strand_aware", "rc_fwd_rows_only") if rc else ("no_rc",)
                for mode in modes:
                    for alpha in alphas:
                        cands = read_pairs(os.path.join(outs[alpha], "pairwise_results.txt"))
                        sc, a, b = score_candidates(S, n, cands, mode, H)
                        for (name, par), arr in sc.items():
                            met = evaluate_scores(to_full_array(arr, a, b, n), y_true, gt_top)
                            rows.append(dict(seed=seed, pool="pipeline_a%g" % alpha, setting=mode, scorer=name, params=par, **met))
                    if not args.no_all_pairs:
                        sc, iu, ju = score_all_pairs(S, n, mode, H)
                        for (name, par), arr in sc.items():
                            met = evaluate_scores(arr.astype(np.float64), y_true, gt_top)
                            rows.append(dict(seed=seed, pool="all_pairs", setting=mode, scorer=name, params=par, **met))
            save()
            print("seed %d done in %.0f s" % (seed, time.time() - t0), flush=True)

    df = save()
    keys = ["pool", "setting", "scorer", "params"]
    metrics = ["auroc", "auprc", "precision_at_n_true", "recall_scored", "top1000_overlap"]
    agg = df.groupby(keys)[metrics].agg(["mean", "std"])
    agg.columns = ["%s_%s" % c for c in agg.columns]
    agg = agg.reset_index()
    agg.to_csv(args.out + "_summary.csv", index=False)
    print("wrote %s.csv and %s_summary.csv" % (args.out, args.out))


if __name__ == "__main__":
    main()
