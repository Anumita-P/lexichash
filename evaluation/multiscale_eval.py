"""Multi-scale LexicHash evaluation (resolve-list item 10).

The repo's multi-scale code (LexicHash_1/LexicHash/run_module.py) runs LexicHash four times with
different (min_k, max_k, n_hash, alpha) and keeps, for every pair, the largest reported match length
(merge_prediction_files). This script re-creates that pipeline with exact, vectorised y-values and
evaluates it with the shared evaluator, next to controls that spend the same number of masks at a
single scale. It answers: is the multi-scale gain a scale effect or just "more masks"?

Sketches follow LexicHash.lexicographic_first (copied below): for each read and mask, the
lexicographically first substring of length <= max_k under the mask's per-position base order.
The match length y of a read pair for one mask is the length of the common prefix of the two
sketches (capped at max_k). A pair is reported by a stage if y >= min_k of that stage and the
pair is within that stage's alpha*n best pairs (ranked by max y, ties broken at random).

Usage:  python3 multiscale_eval.py --fasta reads.fasta.gz --gt gt.txt --out results/multiscale_NCTC1080 --seeds 1 2 3
"""
import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lh_eval as E
import mask_strategies as MS
import scorer_ablation as SA

# (min_k, max_k, n_hash, alpha) of run_module.py in LexicHash_1 (the code, not the report text)
STAGES_CODE = [(63, 128, 700, 10.0), (31, 64, 600, 100.0), (14, 32, 600, 1000.0), (13, 16, 400, 10.0)]

_G = {}


def _init(enc_reads, masks, max_k, min_k):
    _G.update(reads=enc_reads, masks=masks, K=max_k, MIN_K=min_k)


def lexicographic_first(seq, mask, K, MIN_K):
    """Same algorithm as LexicHash.lexicographic_first; returns the K base codes of the minimum."""
    def get_candidate_locs():
        j = 0
        cur = np.arange(len(seq) - MIN_K)
        while j < K:
            nxt = np.where(seq[cur] == mask[j])[0]
            if len(nxt) == 0:
                return cur - j, j
            elif len(nxt) == 1:
                return [cur[nxt[0]] - j], j + 1
            cur = cur[nxt] + 1
            if cur[-1] == len(seq):
                cur = cur[:-1]
            j += 1
        return cur - K, K

    def extend(cands, j):
        while len(cands) > 1 and j < K:
            best, nxt_c = 4, []
            for loc in cands:
                if loc + j < len(seq):
                    bm = (int(seq[loc + j]) - int(mask[j])) % 4
                    if bm < best:
                        best, nxt_c = bm, [loc]
                    elif bm == best:
                        nxt_c.append(loc)
            if not nxt_c:
                break
            cands = nxt_c
            j += 1
        return cands[0]

    cands, n_match = get_candidate_locs()
    idx = extend(cands, n_match)
    out = np.empty(K, dtype=np.uint8)
    window = seq[idx:idx + K]
    out[:len(window)] = window
    if len(window) < K:                      # substring at the end of the read: extend with the
        out[len(window):] = mask[-(K - len(window)):] ^ 3   # lexicographically last bases
    return out


def _sketch_read(i):
    seq = _G["reads"][i]
    return np.array([lexicographic_first(seq, m, _G["K"], _G["MIN_K"]) for m in _G["masks"]])


def sketch_all(enc_reads, masks, K, MIN_K, procs):
    with Pool(procs, initializer=_init, initargs=(enc_reads, masks, K, MIN_K)) as p:
        out = p.map(_sketch_read, range(len(enc_reads)), chunksize=8)
    return np.stack(out)                     # (n, H, K) uint8 base codes


def pack(S):
    """(n, H, K) base codes -> (n, H, C) uint64, 32 bases per word, first base in the top bits."""
    n, H, K = S.shape
    C = (K + 31) // 32
    pad = np.zeros((n, H, C * 32), dtype=np.uint8)
    pad[:, :, :K] = S
    w = (np.uint64(4) ** np.arange(31, -1, -1, dtype=np.uint64))
    return (pad.reshape(n, H, C, 32).astype(np.uint64) * w).sum(axis=3).astype(np.uint64)


def y_all_pairs(P, K, chunk=2000):
    """y (shared prefix length, capped at K) for every pair i<j: uint8 array (n_pairs, H)."""
    n, H, C = P.shape
    iu, ju = np.triu_indices(n, k=1)
    out = np.empty((len(iu), H), dtype=np.uint8)
    for s in range(0, len(iu), chunk):
        a, b = iu[s:s + chunk], ju[s:s + chunk]
        pl = SA.prefix_len(P[a] ^ P[b]).astype(np.int16)          # (m, H, C), 32 where equal
        full = np.cumprod(pl == 32, axis=2)
        before = np.concatenate([np.ones_like(full[:, :, :1]), full[:, :, :-1]], axis=2)
        out[s:s + chunk] = np.minimum((pl * before).sum(axis=2), K)
    return out, iu, ju


def reported_pool(maxy, min_k, alpha, n, rng):
    """Boolean mask of the pairs a stage returns: y >= min_k and within the alpha*n best pairs."""
    ok = maxy >= min_k
    n_pairs = int(n * alpha)
    cand = np.where(ok)[0]
    if len(cand) > n_pairs:
        order = np.lexsort((rng.random(len(cand)), -maxy[cand]))
        cand = cand[order[:n_pairs]]
    sel = np.zeros(len(maxy), dtype=bool)
    sel[cand] = True
    return sel


def expdecay(y, decay=0.1, top_k=100):
    ys = np.sort(y, axis=1)[:, ::-1][:, :top_k].astype(np.float64)
    return ys @ ((1.0 - decay) ** np.arange(ys.shape[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--baseline-hash", type=int, default=500)
    args = ap.parse_args()

    reads = MS.read_fasta(args.fasta)
    enc = [MS.encode(s) for s in reads]
    lens = E.read_lengths(args.fasta)
    n = len(lens)
    true_pairs = E.load_true_pairs(args.gt, lens, E.THETA)
    y_true = E.true_pair_labels(true_pairs, n)
    gt_top = SA.load_gt_top_idx(args.gt, n, 1000)
    print("reads %d  true pairs %d  shortest read %d" % (n, len(true_pairs), lens.min()), flush=True)
    total_masks = sum(s[2] for s in STAGES_CODE)
    rows = []

    def add(seed, name, score, extra=""):
        met = SA.evaluate_scores(score.astype(np.float64), y_true, gt_top)
        rows.append(dict(seed=seed, method=name, extra=extra, **met))
        print("  %-52s auprc %.3f  auroc %.4f  p@n %.3f  top1000 %.1f" % (
            name, met["auprc"], met["auroc"], met["precision_at_n_true"], met["top1000_overlap"]), flush=True)

    for seed in args.seeds:
        t0 = time.time()
        rng = np.random.default_rng(seed)
        print("seed %d" % seed, flush=True)
        stage_y, stage_max = [], []
        for si, (min_k, K, H, alpha) in enumerate(STAGES_CODE):
            masks = rng.integers(0, 4, size=(H, K)).astype(np.uint8)
            S = sketch_all(enc, masks, K, min_k, args.procs)
            y, iu, ju = y_all_pairs(pack(S), K)
            stage_y.append(y)
            stage_max.append(y.max(axis=1))
            print("  stage %d (k %d-%d, %d masks) sketched+compared in %.0f s" % (si + 1, min_k, K, H, time.time() - t0), flush=True)
            del S

        # --- the repo's pipeline: per stage the pairs it returns, then the max reported length
        merged = np.zeros(len(iu))
        for (min_k, K, H, alpha), mx in zip(STAGES_CODE, stage_max):
            sel = reported_pool(mx, min_k, alpha, n, rng)
            merged = np.maximum(merged, np.where(sel, mx, 0))
        add(seed, "multiscale as coded (max of reported lengths)", merged, "%d masks" % total_masks)
        for si, (min_k, K, H, alpha) in enumerate(STAGES_CODE):
            sel = reported_pool(stage_max[si], min_k, alpha, n, rng)
            add(seed, "  stage %d alone (k %d-%d, %d masks, alpha %g)" % (si + 1, min_k, K, H, alpha),
                np.where(sel, stage_max[si], 0), "%d masks" % H)
        # distinct score values: how many ties does the merged score have?
        print("  merged score: %d distinct values over %d pairs" % (len(np.unique(merged)), len(merged)), flush=True)

        # --- controls with a single scale k=32 (stage 3 masks only, then all masks of the k=32 and
        #     shorter stages truncated to 32 are NOT comparable; use fresh random k=32 masks instead)
        for H in (args.baseline_hash, total_masks):
            masks = rng.integers(0, 4, size=(H, 32)).astype(np.uint8)
            S = sketch_all(enc, masks, 32, 14, args.procs)
            y, _, _ = y_all_pairs(pack(S), 32)
            add(seed, "single scale k=32, %d masks, max y" % H, y.max(axis=1), "%d masks" % H)
            add(seed, "single scale k=32, %d masks, exp-decay(0.1,100)" % H, expdecay(y), "%d masks" % H)
            del S, y
        # --- exp-decay on the 4-stage y pooled (each stage's y rescaled to [0,1] by its max_k)
        pooled = np.empty(len(iu))
        for s in range(0, len(iu), 20000):
            blk = np.concatenate([yy[s:s + 20000].astype(np.float32) / K
                                  for yy, (_, K, _, _) in zip(stage_y, STAGES_CODE)], axis=1)
            pooled[s:s + 20000] = expdecay(blk)
        add(seed, "multiscale pooled, y/max_k, exp-decay(0.1,100)", pooled, "%d masks" % total_masks)
        pd.DataFrame(rows).to_csv(args.out + ".csv", index=False)
        print("seed %d done in %.0f s" % (seed, time.time() - t0), flush=True)

    df = pd.DataFrame(rows)
    agg = df.groupby("method")[["auroc", "auprc", "precision_at_n_true", "top1000_overlap"]].agg(["mean", "std"])
    agg.columns = ["%s_%s" % c for c in agg.columns]
    agg.reset_index().to_csv(args.out + "_summary.csv", index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
