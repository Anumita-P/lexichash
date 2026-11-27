#!/usr/bin/env python3
"""Turn DALIGNER output into a LexicHash ground-truth file, and check it against a reference.

Why this exists (resolve item 2): the paper builds its ground truth from DALIGNER alignments
(NCTC datasets) with a threshold of 0.2 on the alignment fraction a / (l1 + l2 - a). The paper does
not give the DALIGNER parameters or say how the alignment length `a` is measured, so this script
(1) converts `LAshow` output to the file format LexicHash reads, and (2) tries the plausible
definitions of `a` and tells you which one best reproduces the shipped
NCTC1080_daligner_ground_truth.txt.

Input  : output of   LAshow reads.dam reads.reads.las > lashow.txt
         lines such as
           1.1  23.1 n   <     0.. 3,395] x [ 3,862.. 7,885>  ~  26.58%   ( 8,066 x 7,885 bps, ...)
         Read numbers are 1-based (DAZZ_DB) and become 0-based in the output.
Output : tab-separated  i  j  a  strand   (0-based read numbers in FASTA order, a = bases)

Usage
  python3 daligner_ground_truth.py lashow.txt --fasta reads.fasta --reference shipped_gt.txt --calibrate
  python3 daligner_ground_truth.py lashow.txt --strand n --a-def A --agg max --out my_gt.txt

Definitions tried by --calibrate
  strand : n    only same-strand alignments (the shipped file has strand '+' only)
           both n and c alignments
  a-def  : A    length of the interval on read A      (aepos - abpos)
           B    length of the interval on read B
           mean (A + B) / 2
           max  max(A, B)
  agg    : max  longest single alignment of the ordered pair (default)
           sum  sum over all alignments of the pair (DALIGNER often reports several pieces)
"""
import argparse
import itertools
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lh_eval import THETA, load_true_pairs, read_lengths  # noqa: E402

LINE = re.compile(
    r"^\s*(\d+)\.\d+\s+(\d+)\.\d+\s+([nc])\s+"
    r"[\[<]\s*([\d,]+)\s*\.\.\s*([\d,]+)\s*[\]>]\s+x\s+"
    r"[\[<]\s*([\d,]+)\s*\.\.\s*([\d,]+)\s*[\]>]")


def _num(s):
    return int(s.replace(",", ""))


UNPARSED = []   # alignment-looking lines that the pattern did not recognise (first few kept)


def parse_lashow(path):
    """List of (a_read0, b_read0, strand, a_len, b_len); read numbers 0-based.

    For complement ('c') alignments LAshow prints the B interval backwards (start > end),
    so lengths are absolute differences."""
    recs = []
    with open(path) as fh:
        for line in fh:
            m = LINE.match(line)
            if not m:
                if line.strip() and ".." in line:
                    UNPARSED.append(line.rstrip())
                continue
            a, b, st = int(m.group(1)) - 1, int(m.group(2)) - 1, m.group(3)
            alen = abs(_num(m.group(5)) - _num(m.group(4)))
            blen = abs(_num(m.group(7)) - _num(m.group(6)))
            recs.append((a, b, st, alen, blen))
    return recs


def measure(alen, blen, a_def):
    if a_def == "A":
        return alen
    if a_def == "B":
        return blen
    if a_def == "mean":
        return (alen + blen) / 2.0
    if a_def == "max":
        return float(max(alen, blen))
    raise ValueError(a_def)


def build(recs, strand, a_def, agg):
    """dict (i, j) -> (a, strand_char) for ordered pairs, i != j."""
    out = {}
    for a, b, st, alen, blen in recs:
        if a == b or (strand == "n" and st != "n"):
            continue
        v = measure(alen, blen, a_def)
        key = (a, b)
        if agg == "sum":
            out[key] = out.get(key, 0.0) + v
        elif v > out.get(key, 0.0):
            out[key] = v
    return {k: (v, "+") for k, v in out.items()}


def write_gt(pairs, path):
    with open(path, "w") as fh:
        for (i, j), (a, s) in sorted(pairs.items()):
            fh.write("%d\t%d\t%.1f\t%s\n" % (i, j, a, s))


def write_tmp(pairs, path):
    write_gt(pairs, path)
    return path


def read_reference(path):
    ref = {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 3:
                ref[(int(p[0]), int(p[1]))] = float(p[2])
    return ref


def calibrate(recs, fasta, reference, workdir):
    lens = read_lengths(fasta)
    ref_true = load_true_pairs(reference, lens)
    ref_pairs = read_reference(reference)
    ref_unordered = {(min(i, j), max(i, j)) for i, j in ref_pairs}
    print("reference: %d lines, %d unordered pairs, %d true pairs at theta >= %.2f"
          % (len(ref_pairs), len(ref_unordered), len(ref_true), THETA))
    print("DALIGNER : %d alignment records parsed (%d with n, %d with c)\n"
          % (len(recs), sum(r[2] == "n" for r in recs), sum(r[2] == "c" for r in recs)))
    hdr = ("strand a-def agg | pairs(unord) true | common true | only ref | only new | "
           "Jaccard | a equal to ref (%) | median a/ref")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for strand, a_def, agg in itertools.product(("n", "both"), ("A", "B", "mean", "max"), ("max", "sum")):
        pairs = build(recs, strand, a_def, agg)
        tmp = write_tmp(pairs, os.path.join(workdir, "_cal_tmp.txt"))
        true = load_true_pairs(tmp, lens)
        common = len(true & ref_true)
        only_ref, only_new = len(ref_true - true), len(true - ref_true)
        union = len(true | ref_true)
        jac = common / union if union else 0.0
        both = [(k, pairs[k][0], ref_pairs[k]) for k in pairs if k in ref_pairs]
        eq = 100.0 * sum(abs(a - r) < 0.5 for _, a, r in both) / len(both) if both else 0.0
        ratios = [a / r for _, a, r in both if r > 0]
        med = statistics.median(ratios) if ratios else float("nan")
        unord = {(min(i, j), max(i, j)) for i, j in pairs}
        rows.append((jac, eq, strand, a_def, agg))
        print("%-6s %-5s %-4s | %9d %5d | %11d | %8d | %8d | %7.3f | %17.1f | %.3f"
              % (strand, a_def, agg, len(unord), len(true), common, only_ref, only_new, jac, eq, med))
    os.remove(os.path.join(workdir, "_cal_tmp.txt"))
    best = max(rows)   # highest Jaccard of the true-pair sets, ties broken by exact-length agreement
    print("\nBest match to the reference: strand=%s a-def=%s agg=%s (Jaccard of true-pair sets %.3f, "
          "a equal to the reference for %.1f %% of common pairs)" % (best[2], best[3], best[4], best[0], best[1]))
    print("Jaccard close to 1.0 and 'a equal' close to 100 means the DALIGNER settings and the length "
          "definition reproduce the shipped ground truth. If both stay low for every row, the "
          "shipped file was made with other DALIGNER parameters (try -e, -l, -k, -h, -s).")


def diagnose(recs, fasta, reference, strand, a_def, agg, n_show=12):
    """Why are some true pairs of the reference missing from the new ground truth?"""
    lens = read_lengths(fasta)
    ref_true = load_true_pairs(reference, lens)
    pairs = build(recs, strand, a_def, agg)
    new_best = {}
    for (i, j), (a, _) in pairs.items():
        k = (min(i, j), max(i, j))
        new_best[k] = max(new_best.get(k, 0.0), a)
    ref_a = {}
    for (i, j), a in read_reference(reference).items():
        k = (min(i, j), max(i, j))
        ref_a[k] = max(ref_a.get(k, 0.0), a)
    missing = sorted(ref_true - {k for k in new_best if k in ref_true and
                                 new_best[k] / (lens[k[0]] + lens[k[1]] - new_best[k]) >= THETA})
    no_aln = [k for k in missing if k not in new_best]
    print("\nTrue pairs of the reference that are not true in the new file: %d" % len(missing))
    print("  of these, no alignment at all in the DALIGNER output: %d; alignment found but shorter: %d"
          % (len(no_aln), len(missing) - len(no_aln)))
    print("  pair (0-based)   len1   len2   ref a   new a   ref theta  new theta")
    for k in missing[:n_show]:
        l1, l2 = lens[k[0]], lens[k[1]]
        r = ref_a[k]
        n = new_best.get(k, 0.0)
        print("  %5d %5d   %6d %6d %7.0f %7.0f   %7.3f   %7.3f"
              % (k[0], k[1], l1, l2, r, n, r / (l1 + l2 - r), n / (l1 + l2 - n)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lashow", help="output of LAshow")
    ap.add_argument("--fasta", help="reads FASTA (needed for --calibrate)")
    ap.add_argument("--reference", help="existing ground-truth file to compare with")
    ap.add_argument("--calibrate", action="store_true", help="try all definitions against --reference")
    ap.add_argument("--strand", choices=("n", "both"), default="n")
    ap.add_argument("--a-def", choices=("A", "B", "mean", "max"), default="A")
    ap.add_argument("--agg", choices=("max", "sum"), default="max")
    ap.add_argument("--diagnose", action="store_true",
                    help="list reference true pairs missing from the new file (uses --strand/--a-def/--agg)")
    ap.add_argument("--out", help="write the ground-truth file here")
    a = ap.parse_args()
    recs = parse_lashow(a.lashow)
    if not recs:
        sys.exit("no alignment lines recognised in %s" % a.lashow)
    if UNPARSED:
        print("WARNING: %d lines contain '..' but were not recognised as alignments; first ones:" % len(UNPARSED))
        for u in UNPARSED[:5]:
            print("   " + u)
    if a.diagnose:
        if not (a.fasta and a.reference):
            sys.exit("--diagnose needs --fasta and --reference")
        diagnose(recs, a.fasta, a.reference, a.strand, a.a_def, a.agg)
    if a.calibrate:
        if not (a.fasta and a.reference):
            sys.exit("--calibrate needs --fasta and --reference")
        calibrate(recs, a.fasta, a.reference, os.path.dirname(os.path.abspath(a.lashow)))
    if a.out:
        pairs = build(recs, a.strand, a.a_def, a.agg)
        write_gt(pairs, a.out)
        print("wrote %d lines to %s" % (len(pairs), a.out))


if __name__ == "__main__":
    main()
