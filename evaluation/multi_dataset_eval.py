#!/usr/bin/env python3
"""Scorer comparison on many NCTC data sets (resolve item 11: do the results hold beyond NCTC1080?).

The LexicHash paper evaluates 40 NCTC data sets (Figure 6). The Spectral Jaccard Similarity
repository (github.com/TavorB/spectral_jaccard_similarity) contains, for 41 NCTC data sets, about 1,000
filtered reads each (folder filtered_fasta/) and the DALIGNER ground truth (folder groundTruths/).
This script runs evaluation/scorer_ablation.py on each of them and summarises the result over the data sets.

Ground-truth files there have 5 columns  i j a len_i len_j  with 1-BASED read numbers
(your NCTC1080 file in LexicHash-main/data is the same data, 0-based with a '+' column; checked line by line).
They are converted to 0-based  i j a +  files in <out-dir>/gt/.

What is compared (strand ignored, as everywhere in evaluation/lh_eval.py; theta >= 0.2):
  maxy                the original LexicHash score (largest shared prefix over all masks)
  expdecay(0.1,100)   the exponential-decay score of the report
  meantop100          plain mean of the 100 largest shared-prefix lengths
on two candidate pools: all read pairs, and the alpha * n pairs the pipeline returns.

Usage (from the repository root):
  git clone --depth 1 --filter=blob:none --sparse https://github.com/TavorB/spectral_jaccard_similarity.git sjs
  (cd sjs && git sparse-checkout set filtered_fasta groundTruths)
  python3 evaluation/multi_dataset_eval.py --cpp-bin ./lexichash --sjs-repo sjs \
      --out-dir evaluation/results/multi_nctc --seeds 1 2 3 --alphas 5
Re-running skips data sets that already have a result file; add --summarize-only to only rebuild the tables.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lh_eval as E  # noqa: E402

CANDIDATES = [("expdecay(0.1,100)", "expdecay", "decay=0.1,top_k=100"),
              ("meantop100", "meantop", "top_k=100")]


def list_datasets(sjs):
    names = sorted(f[:-len("_reads.fasta")] for f in os.listdir(os.path.join(sjs, "filtered_fasta"))
                   if f.endswith("_reads.fasta"))
    return names


def convert_gt(src, dst):
    """SJS ground truth (1-based, 5 columns) -> LexicHash format (0-based, i j a +)."""
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            p = line.split()
            if len(p) < 3:
                continue
            fo.write("%d\t%d\t%s\t+\n" % (int(p[0]) - 1, int(p[1]) - 1, p[2]))
            n += 1
    return n


def fasta_stats(path):
    n_reads, total, gc = 0, 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n_reads += 1
                continue
            s = line.strip().upper()
            total += len(s)
            gc += s.count("G") + s.count("C")
    return n_reads, total / max(n_reads, 1), gc / max(total, 1)


def run_dataset(name, a, gt_dir, runs_dir):
    fasta = os.path.join(a.sjs_repo, "filtered_fasta", name + "_reads.fasta")
    src = os.path.join(a.sjs_repo, "groundTruths", name + "_daligner_ground_truth.txt")
    gt = os.path.join(gt_dir, name + ".txt")
    convert_gt(src, gt)
    out = os.path.join(runs_dir, name)
    if os.path.exists(out + ".csv") and os.path.exists(out + "_summary.csv"):
        print("%s: result exists, skipped" % name, flush=True)
        return
    cmd = [sys.executable, os.path.join(HERE, "scorer_ablation.py"), "--cpp-bin", a.cpp_bin, "--fasta", fasta,
           "--gt", gt, "--n-hash", str(a.n_hash), "--seeds"] + [str(s) for s in a.seeds] + \
          ["--alphas"] + [str(x) for x in a.alphas] + ["--rc-seeds", "--out", out]
    print("%s: running scorer_ablation ..." % name, flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("%s FAILED:\n%s\n%s" % (name, res.stdout[-600:], res.stderr[-600:]), flush=True)
        for ext in (".csv", "_summary.csv"):
            if os.path.exists(out + ext):
                os.remove(out + ext)


def collect(a, gt_dir, runs_dir):
    rows = []
    for name in a.datasets:
        f = os.path.join(runs_dir, name + ".csv")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d = d[d.setting == "no_rc"]
        fasta = os.path.join(a.sjs_repo, "filtered_fasta", name + "_reads.fasta")
        n_reads, mean_len, gc = fasta_stats(fasta)
        lens = E.read_lengths(fasta)
        n_true = len(E.load_true_pairs(os.path.join(gt_dir, name + ".txt"), lens))
        pools = ["all_pairs"] + ["pipeline_a%g" % x for x in a.alphas]
        for pool in pools:
            dp = d[d.pool == pool]
            sets = [("maxy", dp[dp.scorer == "maxy"])] + \
                   [(lab, dp[(dp.scorer == sc) & (dp.params == par)]) for lab, sc, par in CANDIDATES]
            for lab, x in sets:
                if len(x) == 0:
                    continue
                rows.append(dict(dataset=name, n_reads=n_reads, mean_len=mean_len, gc=gc, n_true=n_true, pool=pool,
                                 scorer=lab, auroc=x.auroc.mean(), auprc=x.auprc.mean(),
                                 precision_at_n_true=x.precision_at_n_true.mean(), n_seeds=len(x)))
    return pd.DataFrame(rows)


def boot_ci(x, n=4000, seed=0):
    x = np.asarray(x, dtype=float)
    rng = np.random.RandomState(seed)
    m = [rng.choice(x, len(x)).mean() for _ in range(n)]
    return np.percentile(m, 2.5), np.percentile(m, 97.5)


def summarise(df, out_dir, alphas):
    df.to_csv(os.path.join(out_dir, "multi_dataset_long.csv"), index=False)
    lines = ["# Scorer comparison over NCTC data sets", "",
             "%d data sets; no reverse complements; theta >= 0.2; mean over seeds per data set. "
             "Delta = candidate minus the original max-y score, paired by data set; the interval is a bootstrap "
             "over data sets (95 %%)." % df.dataset.nunique(), ""]
    for pool in ["all_pairs"] + ["pipeline_a%g" % x for x in alphas]:
        p = df[df.pool == pool]
        if p.empty:
            continue
        wide = p.pivot(index="dataset", columns="scorer", values=["auroc", "auprc"])
        lines += ["## Pool: %s" % pool, "",
                  "| scorer | mean AUROC | mean AUPRC | median AUPRC | mean dAUPRC vs max-y | 95 % CI | data sets better | dAUROC (mean) |",
                  "|---|---|---|---|---|---|---|---|"]
        for sc in ["maxy"] + [c[0] for c in CANDIDATES]:
            if ("auprc", sc) not in wide:
                continue
            u, r = wide[("auprc", sc)], wide[("auroc", sc)]
            if sc == "maxy":
                lines.append("| maxy | %.3f | %.3f | %.3f | - | - | - | - |" % (r.mean(), u.mean(), u.median()))
            else:
                du, dr = u - wide[("auprc", "maxy")], r - wide[("auroc", "maxy")]
                lo, hi = boot_ci(du.values)
                lines.append("| %s | %.3f | %.3f | %.3f | %+.3f | %+.3f to %+.3f | %d / %d | %+.3f |"
                             % (sc, r.mean(), u.mean(), u.median(), du.mean(), lo, hi, (du > 0).sum(), len(du), dr.mean()))
        lines.append("")
    # per data set table for the first pipeline pool and all pairs
    ds = df[df.scorer == "maxy"].groupby("dataset")[["n_reads", "mean_len", "gc", "n_true"]].first()
    lines += ["## Per data set (AUPRC)", "", "| data set | reads | mean len | GC | true pairs | " +
              " | ".join("%s maxy | %s exp | %s dExp" % (pl, pl, pl) for pl in ["all_pairs", "pipeline_a%g" % alphas[0]]) + " |",
              "|---|---|---|---|---|" + "---|---|---|" * 2]
    for name, r in ds.iterrows():
        cells = []
        for pl in ["all_pairs", "pipeline_a%g" % alphas[0]]:
            q = df[(df.dataset == name) & (df.pool == pl)].set_index("scorer")
            if "maxy" in q.index and "expdecay(0.1,100)" in q.index:
                m, e = q.loc["maxy", "auprc"], q.loc["expdecay(0.1,100)", "auprc"]
                cells += ["%.3f" % m, "%.3f" % e, "%+.3f" % (e - m)]
            else:
                cells += ["-", "-", "-"]
        lines.append("| %s | %d | %.0f | %.2f | %d | %s |" % (name, r.n_reads, r.mean_len, r.gc, r.n_true, " | ".join(cells)))
    lines.append("")
    with open(os.path.join(out_dir, "multi_dataset_summary.md"), "w") as fh:
        fh.write("\n".join(lines))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pools = ["all_pairs", "pipeline_a%g" % alphas[0]]
        fig, ax = plt.subplots(1, len(pools), figsize=(5.2 * len(pools), 4.6))
        for k, pl in enumerate(pools):
            w = df[df.pool == pl].pivot(index="dataset", columns="scorer", values="auprc")
            if "expdecay(0.1,100)" not in w:
                continue
            ax[k].scatter(w["maxy"], w["expdecay(0.1,100)"], s=22)
            ax[k].plot([0, 1], [0, 1], "k--", lw=0.8)
            ax[k].set_xlabel("AUPRC, max-y (original)")
            ax[k].set_ylabel("AUPRC, exponential decay (0.1, top 100)")
            ax[k].set_title(pl)
            ax[k].set_xlim(0, 1)
            ax[k].set_ylim(0, 1)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "multi_dataset_auprc.png"), dpi=140)
    except Exception as ex:  # plotting is optional
        print("figure skipped: %s" % ex)
    print("wrote %s/multi_dataset_summary.md, multi_dataset_long.csv" % out_dir)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpp-bin", required=True)
    ap.add_argument("--sjs-repo", required=True, help="clone of TavorB/spectral_jaccard_similarity (filtered_fasta, groundTruths)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--datasets", nargs="*", default=None, help="default: all found in filtered_fasta/")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--alphas", type=float, nargs="+", default=[5])
    ap.add_argument("--n-hash", type=int, default=500)
    ap.add_argument("--summarize-only", action="store_true")
    a = ap.parse_args()
    a.cpp_bin = os.path.abspath(a.cpp_bin)
    a.datasets = a.datasets or list_datasets(a.sjs_repo)
    gt_dir, runs_dir = os.path.join(a.out_dir, "gt"), os.path.join(a.out_dir, "runs")
    os.makedirs(gt_dir, exist_ok=True)
    os.makedirs(runs_dir, exist_ok=True)
    if not a.summarize_only:
        for name in a.datasets:
            run_dataset(name, a, gt_dir, runs_dir)
    df = collect(a, gt_dir, runs_dir)
    if df.empty:
        sys.exit("no results to summarise")
    summarise(df, a.out_dir, a.alphas)


if __name__ == "__main__":
    main()
