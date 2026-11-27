#!/usr/bin/env python3
"""Runtime / memory benchmark on the E. coli PacBio data set (resolve item 2; paper Section 3.2, Figure 7, Figure S7).

Paper set-up (main text and supplement): E. coli SRR11434956 (PacBio RS II), 144,238 reads; reads with an average
error rate of 30 % or more are removed with fastq-filter, leaving 113,378. Data sets of 10, 30, 50, 70 and 90
thousand reads are drawn from it by downsampling. 100 hash functions, T = 5n pairs returned, LexicHash k_max = 32
(k_min = 16), MinHash k = 12. The full-set numbers quoted for ONE cpu: LexicHash 100.1 CPU-min and 515.7 MB,
MinHash 339.3 CPU-min and 1.95 GB. No ground truth is needed: this is a pure runtime/memory experiment.

Two steps.

1) prepare: stream the .fastq.gz once, drop reads whose average error rate (mean of 10^(-Q/10) over the bases,
   Phred+33) is >= the threshold, write the rest to <out-dir>/filtered.fasta. The kept number is compared with
   113,378 (the paper); counts for a few other thresholds are printed so a mismatch can be explained.

     python3 evaluation/ecoli_runtime_benchmark.py prepare --fastq SRR11434956_subreads.fastq.gz --out-dir ~/ecoli/work

2) run: for each size, draw a seeded random subset of filtered.fasta (nested: the 10k set is inside the 30k set, ...),
   and time each method on ONE cpu. CPU time (user + system) and peak resident memory come from the kernel's
   accounting for the child process (os.wait4), wall time is measured too.
     cpp           the C++ program in mycode/ (our re-implementation), writes sketches.txt and pairwise_results.txt
     py_lexichash  the original Python package in LexicHash-main (needs numpy, scipy, tqdm ...)
     py_minhash    the same package with --method minhash --k 12  (slow: the paper needed 339 CPU-min on the full set)

     python3 evaluation/ecoli_runtime_benchmark.py run --work-dir ~/ecoli/work --cpp-bin ./lexichash \
         --sizes 10000 30000 50000 70000 90000 0 --methods cpp py_lexichash
   (size 0 = the whole filtered set). Results are appended to <work-dir>/runtime_results.csv; already finished
   (method, size) pairs are skipped, so the run can be resumed. `summary` prints the table and draws the figure.
"""
import argparse
import csv
import gzip
import os
import random
import subprocess
import sys
import time

PAPER_N_FILTERED = 113378
PAPER_N_ALL = 144238
PAPER = {"cpp": None, "py_lexichash": (100.1, 515.7), "py_minhash": (339.3, 1950.0)}  # CPU-min, MB (full set, 1 cpu)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


# ---------------------------------------------------------------- prepare
def prepare(a):
    import numpy as np
    os.makedirs(a.out_dir, exist_ok=True)
    lut = np.full(256, np.nan)
    for q in range(0, 94):
        lut[q + 33] = 10.0 ** (-q / 10.0)
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]
    n_below = {t: 0 for t in thresholds}
    n = kept = 0
    lens_kept = 0
    opener = gzip.open if a.fastq.endswith(".gz") else open
    out_path = os.path.join(a.out_dir, "filtered.fasta")
    t0 = time.time()
    with opener(a.fastq, "rb") as fh, open(out_path, "w") as out:
        while True:
            h = fh.readline()
            if not h:
                break
            s, p, q = fh.readline().rstrip(b"\r\n"), fh.readline(), fh.readline().rstrip(b"\r\n")
            if not h.startswith(b"@") or not p.startswith(b"+") or len(s) != len(q):
                sys.exit("record %d is not a plain 4-line FASTQ record (header %r); wrapped FASTQ is not supported here"
                         % (n + 1, h[:40]))
            n += 1
            if len(s) == 0:
                continue
            e = lut[np.frombuffer(q, dtype=np.uint8)]
            if np.isnan(e).any():
                sys.exit("record %d has a quality character outside Phred+33 range" % n)
            err = float(e.mean())
            for t in thresholds:
                if err < t:
                    n_below[t] += 1
            if err < a.max_error_rate:
                kept += 1
                lens_kept += len(s)
                out.write(">r%d\n%s\n" % (n - 1, s.decode()))
            if n % 20000 == 0:
                print("  %d reads, %.0f s" % (n, time.time() - t0), flush=True)
    print("reads in file: %d (paper: %d)" % (n, PAPER_N_ALL))
    print("kept with average error rate < %.2f: %d (paper: %d); mean kept length %.0f" %
          (a.max_error_rate, kept, PAPER_N_FILTERED, lens_kept / max(kept, 1)))
    print("reads below other thresholds: " + ", ".join("<%.1f: %d" % (t, n_below[t]) for t in thresholds))
    if n == PAPER_N_ALL and kept != PAPER_N_FILTERED:
        print("NOTE: the kept number differs from the paper's 113,378. fastq-filter may compute the average error "
              "rate differently (e.g. from the mean quality score instead of the mean error probability); "
              "report the number you got, it does not change the runtime scaling.")
    with open(os.path.join(a.out_dir, "prepare.txt"), "w") as fh:
        fh.write("reads_in_file\t%d\nkept\t%d\nthreshold\t%g\n" % (n, kept, a.max_error_rate))
    print("wrote", out_path)


# ---------------------------------------------------------------- run
def fasta_records(path):
    """Byte offsets and sequence lengths of the records of a FASTA file written by prepare (one line per sequence)."""
    offs, lens = [], []
    with open(path, "rb") as fh:
        pos = 0
        for line in fh:
            if line.startswith(b">"):
                offs.append(pos)
            else:
                lens.append(len(line.rstrip(b"\r\n")))
            pos += len(line)
    return offs, lens


def make_subset(src, offs, lens, size, seed, min_length, dst):
    """Seeded nested random subset (size 0 = all) of the reads that are at least min_length long.

    The length cut is needed because the original Python code crashes on reads shorter than min_k
    (IndexError, reproduced with a 10-base read) and the C++ program skips reads shorter than max_k, so both
    would otherwise see different read sets."""
    idx = [k for k in range(len(offs)) if lens[k] >= min_length]
    random.Random(seed).shuffle(idx)            # nested subsets: prefixes of one shuffled order
    chosen = sorted(idx[:size]) if size else sorted(idx)
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        for k in chosen:
            fi.seek(offs[k])
            fo.write(fi.readline())
            fo.write(fi.readline())
    return len(chosen)


def timed(cmd, cwd=None):
    t0 = time.time()
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    _, status, ru = os.wait4(p.pid, 0)
    wall = time.time() - t0
    err = p.stderr.read().decode(errors="replace")
    p.stderr.close()
    p.returncode = os.waitstatus_to_exitcode(status)
    return p.returncode, wall, ru.ru_utime + ru.ru_stime, ru.ru_maxrss / 1024.0, err


def command(method, a, fasta, out_dir):
    if method == "cpp":
        return [a.cpp_bin, "--fasta", fasta, "--out", out_dir, "--n_hash", str(a.n_hash), "--min_k", "16",
                "--max_k", "32", "--alpha", str(a.alpha), "--n_cpu", "1", "--seed", "1"], None
    run_module = os.path.join(a.py_dir, "run_module.py")
    base = [sys.executable, run_module, "--fasta", fasta, "--out", out_dir, "--n_hash", str(a.n_hash),
            "--alpha", str(a.alpha), "--n_cpu", "1"]
    if method == "py_lexichash":
        return base + ["--method", "lexichash", "--min_k", "16", "--max_k", "32"], a.py_dir
    return base + ["--method", "minhash", "--k", "12"], a.py_dir


def run(a):
    src = os.path.join(a.work_dir, "filtered.fasta")
    if not os.path.exists(src):
        sys.exit("run `prepare` first (missing %s)" % src)
    res_path = os.path.join(a.work_dir, "runtime_results.csv")
    done = set()
    if os.path.exists(res_path):
        with open(res_path) as fh:
            done = {(r["method"], int(r["requested_size"])) for r in csv.DictReader(fh) if r["status"] == "ok"}
    offs, lens = fasta_records(src)
    n_ok = sum(1 for x in lens if x >= a.min_length)
    print("filtered.fasta: %d reads, %d of them at least %d bases long (the others are not used)"
          % (len(offs), n_ok, a.min_length))
    new_file = not os.path.exists(res_path)
    with open(res_path, "a", newline="") as fh:
        w = csv.writer(fh)
        if new_file:
            w.writerow(["method", "requested_size", "n_reads", "wall_s", "cpu_s", "cpu_min", "peak_mem_mb", "n_hash",
                        "alpha", "status"])
        for size in a.sizes:
            sub = os.path.join(a.work_dir, "subset_%d.fasta" % size)
            n = None
            for method in a.methods:
                if (method, size) in done:
                    print("%s size %d: done before, skipped" % (method, size))
                    continue
                if n is None:
                    n = make_subset(src, offs, lens, size, a.seed, a.min_length, sub)
                out_dir = os.path.join(a.work_dir, "out_%s_%d" % (method, size))
                os.makedirs(os.path.join(out_dir, "tmp"), exist_ok=True)
                cmd, cwd = command(method, a, sub, out_dir)
                print("%s, %d reads ..." % (method, n), flush=True)
                rc, wall, cpu, mem, err = timed(cmd, cwd)
                status = "ok" if rc == 0 else "failed(%d)" % rc
                if rc != 0:
                    print("   FAILED (not recorded, will be retried on the next run):\n" + err[-800:])
                    continue
                w.writerow([method, size, n, "%.1f" % wall, "%.1f" % cpu, "%.2f" % (cpu / 60), "%.1f" % mem,
                            a.n_hash, a.alpha, status])
                fh.flush()
                print("   wall %.0f s, cpu %.1f min, peak %.0f MB, %s" % (wall, cpu / 60, mem, status), flush=True)
                if not a.keep_outputs:
                    for f in os.listdir(out_dir):
                        p = os.path.join(out_dir, f)
                        if os.path.isfile(p):
                            os.remove(p)
            if n is not None and not a.keep_subsets:
                os.remove(sub)


# ---------------------------------------------------------------- summary
def summary(a):
    res_path = os.path.join(a.work_dir, "runtime_results.csv")
    rows = [r for r in csv.DictReader(open(res_path)) if r["status"] == "ok"]
    methods = sorted({r["method"] for r in rows})
    lines = ["| method | reads | CPU-min | wall-min | peak memory (MB) |", "|---|---|---|---|---|"]
    for m in methods:
        for r in sorted((r for r in rows if r["method"] == m), key=lambda r: int(r["n_reads"])):
            lines.append("| %s | %s | %.1f | %.1f | %.0f |" % (m, r["n_reads"], float(r["cpu_min"]),
                                                              float(r["wall_s"]) / 60, float(r["peak_mem_mb"])))
    lines += ["", "Paper, full set, one CPU: LexicHash (Python) 100.1 CPU-min / 515.7 MB; MinHash 339.3 CPU-min / 1.95 GB."]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(a.work_dir, "runtime_summary.md"), "w") as fh:
        fh.write(text + "\n")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for m in methods:
            rr = sorted((r for r in rows if r["method"] == m), key=lambda r: int(r["n_reads"]))
            x = [int(r["n_reads"]) / 1000 for r in rr]
            ax[0].plot(x, [float(r["cpu_min"]) for r in rr], "o-", label=m)
            ax[1].plot(x, [float(r["peak_mem_mb"]) for r in rr], "o-", label=m)
        ax[0].set_ylabel("CPU time (min, 1 cpu)")
        ax[1].set_ylabel("peak memory (MB)")
        for z in ax:
            z.set_xlabel("reads (thousands)")
            z.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(a.work_dir, "runtime_scaling.png"), dpi=140)
    except Exception as ex:
        print("figure skipped:", ex)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--fastq", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--max-error-rate", type=float, default=0.3)
    p.set_defaults(f=prepare)
    r = sub.add_parser("run")
    r.add_argument("--work-dir", required=True)
    r.add_argument("--cpp-bin", default=os.path.join(REPO, "lexichash"))
    r.add_argument("--py-dir", default=os.path.join(REPO, "LexicHash-main"))
    r.add_argument("--methods", nargs="+", default=["cpp"], choices=["cpp", "py_lexichash", "py_minhash"])
    r.add_argument("--sizes", type=int, nargs="+", default=[10000, 30000, 50000, 70000, 90000, 0])
    r.add_argument("--n-hash", type=int, default=100)
    r.add_argument("--alpha", type=float, default=5)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--min-length", type=int, default=32,
                   help="ignore reads shorter than this (default 32 = max_k; see make_subset)")
    r.add_argument("--keep-subsets", action="store_true")
    r.add_argument("--keep-outputs", action="store_true")
    r.set_defaults(f=run)
    s = sub.add_parser("summary")
    s.add_argument("--work-dir", required=True)
    s.set_defaults(f=summary)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
