"""Markdown tables for the mask-strategy and multi-scale results (resolve-list items 9 and 10).

  python3 summarize_batch4.py masks  results/mask_ablation_NCTC1080_h20.csv  > results/mask_ablation_NCTC1080_h20.md
  python3 summarize_batch4.py multiscale results/multiscale_NCTC1080.csv     > results/multiscale_NCTC1080.md
"""
import sys

import numpy as np
import pandas as pd


def welch_ci(a, b):
    """Difference of means a - b and its 95% interval (normal approximation of Welch's test)."""
    a, b = np.asarray(a), np.asarray(b)
    d = a.mean() - b.mean()
    se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return d, d - 1.96 * se, d + 1.96 * se


def masks(path):
    df = pd.read_csv(path)
    info_path = path.replace(".csv", "_masks_info.csv")
    seeds = df.seed.nunique()
    print("Seeds: %d. Entries are AUPRC, mean +- std over seeds; the last column is the difference to `random` "
          "(mean, 95%% interval; an interval containing 0 means no detectable difference).\n" % seeds)
    for pool in ["all_pairs", "pipeline_a5"]:
        for scorer, par, label in [("maxy", "", "max y (original LexicHash score)"),
                                   ("expdecay", "decay=0.1,top_k=100", "exp-decay (0.1, top 100)")]:
            d = df[(df.pool == pool) & (df.scorer == scorer) & (df.params.fillna("") == par)]
            if d.empty:
                continue
            print("### %s, %s\n" % (pool, label))
            print("| strategy | AUPRC | AUROC | top-1000 overlap (%) | AUPRC minus random |")
            print("|---|---|---|---|---|")
            ref = d[d.strategy == "random"].auprc.values
            order = d.groupby("strategy").auprc.mean().sort_values(ascending=False).index
            for st in order:
                x = d[d.strategy == st]
                if st == "random":
                    diff = "(reference)"
                else:
                    m, lo, hi = welch_ci(x.auprc.values, ref)
                    diff = "%+.3f [%+.3f, %+.3f]" % (m, lo, hi)
                print("| %s | %.3f +- %.3f | %.4f +- %.4f | %.1f +- %.1f | %s |" % (
                    st, x.auprc.mean(), x.auprc.std(), x.auroc.mean(), x.auroc.std(),
                    x.top1000_overlap.mean(), x.top1000_overlap.std(), diff))
            print()
    try:
        inf = pd.read_csv(info_path)
        cols = [c for c in ["mean_count_of_picked", "reads_high", "reads_low", "source_reads", "unique_masks"] if c in inf]
        print("### Mask-set properties (mean over seeds)\n")
        print(inf.groupby("strategy")[cols].mean().round(2).to_markdown())
    except FileNotFoundError:
        pass


def multiscale(path):
    df = pd.read_csv(path)
    print("Seeds: %d. NCTC1080, no reverse complements, all pairs scored, ground truth theta >= 0.2.\n" % df.seed.nunique())
    print("| method | masks | AUPRC | AUROC | precision@n_true | top-1000 overlap (%) |")
    print("|---|---|---|---|---|---|")
    for m, x in df.groupby("method", sort=False):
        print("| %s | %s | %.3f +- %.3f | %.4f | %.3f | %.1f |" % (
            m.strip(), x.extra.iloc[0].split()[0], x.auprc.mean(), x.auprc.std(), x.auroc.mean(),
            x.precision_at_n_true.mean(), x.top1000_overlap.mean()))


if __name__ == "__main__":
    {"masks": masks, "multiscale": multiscale}[sys.argv[1]](sys.argv[2])
