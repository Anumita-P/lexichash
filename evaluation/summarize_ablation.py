#!/usr/bin/env python3
"""Turn the raw ablation CSV(s) written by scorer_ablation.py into a Markdown report.

Usage:
  python3 evaluation/summarize_ablation.py --csv part1.csv part2.csv --title "NCTC1080" --out report.md

For each (setting, pool) it lists: original max-y, the exp-decay score as it is in the code today
(decay 0.1, top_k 100), the best exp-decay setting of the grid, the best top-k mean, the best
"number of masks with y >= t" and the likelihood-style score. Numbers are mean +- standard
deviation over mask seeds. "Best" is chosen by mean AUPRC on the same data, so it is optimistic;
it shows what the score family can reach, not an unbiased estimate.
"""
import argparse

import numpy as np
import pandas as pd


def fmt(m, s):
    return "%.3f +- %.3f" % (m, 0.0 if np.isnan(s) else s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="+", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    d = pd.concat([pd.read_csv(f) for f in args.csv], ignore_index=True)
    d["params"] = d["params"].fillna("")
    lines = ["# Scorer ablation: %s" % args.title, ""]
    lines.append("Seeds per row: %s. Metrics: AUROC and AUPRC over all pairs of reads (unscored pairs = 0), "
                 "ground truth theta >= 0.2. Mean +- std over mask seeds." % ", ".join(
                     "%s=%d" % (k, v) for k, v in d.groupby("setting")["seed"].nunique().items()))
    lines.append("")

    pool_order = ["pipeline_a1", "pipeline_a2", "pipeline_a5", "pipeline_a10", "all_pairs"]
    for setting in [s for s in ["no_rc", "rc_strand_aware", "rc_fwd_rows_only"] if s in set(d.setting)]:
        lines.append("## Setting: %s" % setting)
        lines.append("")
        for pool in [p for p in pool_order if p in set(d[d.setting == setting].pool)]:
            x = d[(d.setting == setting) & (d.pool == pool)]
            g = x.groupby(["scorer", "params"])
            agg = g[["auroc", "auprc", "precision_at_n_true", "recall_scored", "top1000_overlap"]].agg(["mean", "std"])
            agg.columns = ["%s_%s" % c for c in agg.columns]
            agg = agg.reset_index()
            seeds = x.seed.nunique()

            def best(scorer):
                a = agg[agg.scorer == scorer]
                return a.loc[a.auprc_mean.idxmax()]

            picks = [("max y (original LexicHash)", agg[agg.scorer == "maxy"].iloc[0]),
                     ("exp-decay as in code (decay 0.1, top_k 100)",
                      agg[(agg.scorer == "expdecay") & (agg.params == "decay=0.1,top_k=100")].iloc[0]),
                     ("exp-decay, best of grid", best("expdecay")),
                     ("mean of top-k y, best", best("meantop")),
                     ("count of masks with y >= t, best", best("count_ge")),
                     ("likelihood-style sum 4^(y-14)", agg[agg.scorer == "lik"].iloc[0])]
            lines.append("### Pool: %s  (%d seed%s; true pairs inside the pool: %.1f%%)" % (
                pool, seeds, "" if seeds == 1 else "s", 100 * agg.recall_scored_mean.iloc[0]))
            lines.append("")
            lines.append("| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |")
            lines.append("|---|---|---|---|---|---|")
            for name, r in picks:
                lines.append("| %s | %s | %s | %s | %s | %.1f |" % (
                    name, r.params or "-", fmt(r.auprc_mean, r.auprc_std), fmt(r.auroc_mean, r.auroc_std),
                    fmt(r.precision_at_n_true_mean, r.precision_at_n_true_std), r.top1000_overlap_mean))
            # paired difference to max y
            base = x[x.scorer == "maxy"].set_index("seed").auprc
            for name, r in picks[1:4]:
                sel = x[(x.scorer == r.scorer) & (x.params == r.params)].set_index("seed").auprc
                diff = (sel - base).dropna()
                if len(diff):
                    lines.append("")
                    lines.append("AUPRC gain of %s over max y: %+.3f (better in %d of %d seeds)." % (
                        name, diff.mean(), int((diff > 0).sum()), len(diff)))
            lines.append("")

        # exp-decay grid for the pool with the most room (all pairs)
        for pool in ["pipeline_a5", "all_pairs"]:
            x = d[(d.setting == setting) & (d.pool == pool) & (d.scorer == "expdecay")]
            if x.empty:
                continue
            x = x.copy()
            x["decay"] = x.params.str.extract(r"decay=([0-9.]+)")[0].astype(float)
            x["top_k"] = x.params.str.extract(r"top_k=(\w+)")[0]
            piv = x.groupby(["decay", "top_k"]).auprc.mean().unstack()
            order = [c for c in ["10", "25", "50", "100", "all"] if c in piv.columns]
            piv = piv[order]
            lines.append("### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool %s" % pool)
            lines.append("")
            lines.append("| decay | " + " | ".join("top_k=" + c for c in order) + " |")
            lines.append("|---|" + "---|" * len(order))
            for dec, row in piv.iterrows():
                lines.append("| %g | " % dec + " | ".join("%.3f" % v for v in row.values) + " |")
            lines.append("")

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
