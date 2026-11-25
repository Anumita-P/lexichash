"""Compare mask-generation strategies with the shared evaluator (resolve-list item 9).

For every strategy and seed: build the masks (mask_strategies.py), write them to a file, run the
C++ binary with --mask_file (so the sketches use exactly these masks), then score the pairs with the
scorers of scorer_ablation.py and evaluate against the ground truth (lh_eval.py).

Pools: all_pairs (every pair scored) and pipeline_a<alpha> (the alpha*n pairs the C++ prefix-tree
step returns). Setting: no_rc (single-strand data such as NCTC1080).

Example:
  python3 mask_ablation.py --cpp-bin build/lexichash --fasta NCTC1080_reads.fasta.gz \
      --gt NCTC1080_daligner_ground_truth.txt --n-hash 500 --seeds 1 2 3 4 5 --out results/mask_ablation_NCTC1080_h500
"""
import argparse
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lh_eval as E
import mask_strategies as MS
import scorer_ablation as SA

# the few scorers that matter for this comparison (the full grid is in scorer_ablation.py)
SA.DECAYS = [0.1]
SA.TOP_KS = [100]
SA.MEAN_KS = [100]
SA.COUNT_TS = []


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpp-bin", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--n-hash", type=int, default=500)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--strategies", nargs="+", default=list(MS.STRATEGIES))
    ap.add_argument("--alphas", type=float, nargs="+", default=[5])
    ap.add_argument("--out", required=True)
    ap.add_argument("--rc", action="store_true", help="also use reverse complements (strand-aware scoring)")
    args = ap.parse_args()
    binary = os.path.abspath(args.cpp_bin)
    K = SA.MAX_K

    reads = MS.read_fasta(args.fasta)
    lens = E.read_lengths(args.fasta)
    n = len(lens)
    assert n == len(reads)
    true_pairs = E.load_true_pairs(args.gt, lens, E.THETA)
    y_true = E.true_pair_labels(true_pairs, n)
    gt_top = SA.load_gt_top_idx(args.gt, n, 1000)
    print("reads %d, true pairs %d" % (n, len(true_pairs)), flush=True)

    rows, infos = [], []
    mode = "rc_strand_aware" if args.rc else "no_rc"
    with tempfile.TemporaryDirectory() as wd:
        for seed in args.seeds:
            for name in args.strategies:
                t0 = time.time()
                masks, info = MS.make_masks(name, reads, args.n_hash, K, seed)
                mf = os.path.join(wd, "m_%s_%d.txt" % (name, seed))
                MS.write_mask_file(mf, masks)
                infos.append(dict(seed=seed, strategy=name, **info))
                outs = {a: SA.run_cpp(binary, args.fasta, wd, args.n_hash, seed, a, mf,
                                      "%s_%d_a%g" % (name, seed, a), rc=args.rc) for a in args.alphas}
                S = SA.read_sketches(os.path.join(outs[args.alphas[0]], "sketches.txt"))
                H = S.shape[1]
                for a in args.alphas:
                    cands = SA.read_pairs(os.path.join(outs[a], "pairwise_results.txt"))
                    sc, aa, bb = SA.score_candidates(S, n, cands, mode, H)
                    for (sn, par), arr in sc.items():
                        met = SA.evaluate_scores(SA.to_full_array(arr, aa, bb, n), y_true, gt_top)
                        rows.append(dict(seed=seed, strategy=name, pool="pipeline_a%g" % a, scorer=sn, params=par, **met))
                sc, iu, ju = SA.score_all_pairs(S, n, mode, H)
                for (sn, par), arr in sc.items():
                    met = SA.evaluate_scores(arr.astype(np.float64), y_true, gt_top)
                    rows.append(dict(seed=seed, strategy=name, pool="all_pairs", scorer=sn, params=par, **met))
                print("seed %d %-14s %.0f s  %s" % (seed, name, time.time() - t0, info), flush=True)
                pd.DataFrame(rows).to_csv(args.out + ".csv", index=False)
                remove_run_dirs(wd, "%s_%d_a" % (name, seed))
    df = pd.DataFrame(rows)
    df.to_csv(args.out + ".csv", index=False)
    pd.DataFrame(infos).to_csv(args.out + "_masks_info.csv", index=False)
    keys = ["strategy", "pool", "scorer", "params"]
    metrics = ["auroc", "auprc", "precision_at_n_true", "recall_scored", "top1000_overlap"]
    agg = df.groupby(keys)[metrics].agg(["mean", "std"])
    agg.columns = ["%s_%s" % c for c in agg.columns]
    agg.reset_index().to_csv(args.out + "_summary.csv", index=False)
    print("wrote", args.out)


def remove_run_dirs(wd, prefix):
    """Delete the per-run output folders (they hold large sketch files)."""
    import glob
    import shutil
    for d in glob.glob(os.path.join(wd, prefix + "*")):
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    main()
