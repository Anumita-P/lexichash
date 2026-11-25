#!/usr/bin/env python3
"""Tests for how the C++ program reads its input (batch 2 fixes).

Each check builds a small synthetic FASTA file, runs the C++ binary on it with the SAME
mask file, and compares the outputs (sketches.txt and pairwise_results.txt).

  1. gzip      : a .fasta.gz file gives exactly the same output as the plain .fasta
  2. CRLF      : Windows line endings and a missing final newline change nothing
  3. N / case  : a read containing N, IUPAC codes and lower case gives the same output as the
                 same read with those characters removed and upper-cased (the Python
                 clean_seqs behaviour). Before the fix N was silently read as A.
  4. short read: a read shorter than max_k keeps its place in the numbering. Output for a
                 file with the short read at record 3 equals the output for the file without
                 it, with every index from 3 onwards shifted by one. Before the fix the short
                 read was dropped and all later reads were renumbered.
  5. threads   : --n_cpu 1 and --n_cpu 8 give identical output on an N-heavy file

Usage (from the repository root):
  g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native -lz
  python3 tests/test_input_handling.py --cpp-bin ./lexichash
"""
import argparse
import gzip
import os
import random
import subprocess
import sys
import tempfile

import numpy as np

MAX_K, MIN_K, N_HASH = 20, 10, 60


def make_reads(seed=1, n_reads=40, genome_len=2500):
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(genome_len))
    reads = []
    for _ in range(n_reads):
        length = rng.randint(250, 450)
        start = rng.randint(0, genome_len - length)
        reads.append(genome[start:start + length])
    return reads


def write_fasta(path, reads, eol="\n", final_newline=True, gz=False):
    text = eol.join(">r%d%s%s" % (i, eol, s) for i, s in enumerate(reads))
    if final_newline:
        text += eol
    data = text.encode()
    if gz:
        with gzip.open(path, "wb") as fh:
            fh.write(data)
    else:
        with open(path, "wb") as fh:
            fh.write(data)


def write_masks(path, seed=7):
    masks = np.random.RandomState(seed).randint(4, size=(N_HASH, MAX_K))
    with open(path, "w") as fh:
        for m in masks:
            fh.write("".join(str(x) for x in m) + "\n")


def run(binary, fasta, mask_file, workdir, name, n_cpu=2):
    out = os.path.join(workdir, name)
    os.makedirs(out, exist_ok=True)
    cmd = [binary, "--out", out, "--fasta", fasta, "--n_hash", str(N_HASH), "--alpha", "0",
           "--n_cpu", str(n_cpu), "--max_k", str(MAX_K), "--min_k", str(MIN_K),
           "--mask_file", mask_file]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError("run %s failed (%d):\n%s\n%s" % (name, res.returncode, res.stdout[-500:], res.stderr[-500:]))
    with open(os.path.join(out, "pairwise_results.txt")) as fh:
        pairs = sorted(l.strip() for l in fh if l.strip())
    with open(os.path.join(out, "sketches.txt")) as fh:
        sketches = [l.rstrip("\n") for l in fh]
    return pairs, sketches


def shift_pairs(pairs, removed_1based):
    """Map pair lines from a file WITHOUT record `removed_1based` back to the file with it."""
    shifted = []
    for line in pairs:
        a, b, ml, strand = line.split("\t")
        a, b = int(a), int(b)
        a += a >= removed_1based
        b += b >= removed_1based
        shifted.append("%d\t%d\t%s\t%s" % (a, b, ml, strand))
    return sorted(shifted)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpp-bin", required=True)
    args = ap.parse_args()
    binary = os.path.abspath(args.cpp_bin)

    failures = []

    def check(name, ok, detail=""):
        print("%-8s %s %s" % ("PASS" if ok else "FAIL", name, detail))
        if not ok:
            failures.append(name)

    reads = make_reads()
    with tempfile.TemporaryDirectory() as wd:
        masks = os.path.join(wd, "masks.txt")
        write_masks(masks)

        base = os.path.join(wd, "base.fasta")
        write_fasta(base, reads)
        ref_pairs, ref_sketches = run(binary, base, masks, wd, "base")
        check("baseline produces pairs", len(ref_pairs) > 0, "(%d pairs)" % len(ref_pairs))

        # 1. gzip
        write_fasta(os.path.join(wd, "base.fasta.gz"), reads, gz=True)
        p, s = run(binary, os.path.join(wd, "base.fasta.gz"), masks, wd, "gz")
        check("1 gzip input == plain input", p == ref_pairs and s == ref_sketches)

        # 2. CRLF and no final newline
        write_fasta(os.path.join(wd, "crlf.fasta"), reads, eol="\r\n", final_newline=False)
        p, s = run(binary, os.path.join(wd, "crlf.fasta"), masks, wd, "crlf")
        check("2 CRLF + no final newline == plain input", p == ref_pairs and s == ref_sketches)

        # 3. N, IUPAC codes and lower case in one read
        dirty = list(reads)
        r = reads[5]
        d = list(r)
        for pos in (10, 11, 80, 150, 151, 152, 200):
            d[pos] = "N"
        d[40] = "R"
        d[100] = "y"
        for pos in range(120, 135):
            d[pos] = d[pos].lower()
        dirty[5] = "".join(d)
        cleaned = list(reads)
        cleaned[5] = "".join(c for c in dirty[5].upper() if c in "ACGT")
        write_fasta(os.path.join(wd, "dirty.fasta"), dirty)
        write_fasta(os.path.join(wd, "cleaned.fasta"), cleaned)
        p1, s1 = run(binary, os.path.join(wd, "dirty.fasta"), masks, wd, "dirty")
        p2, s2 = run(binary, os.path.join(wd, "cleaned.fasta"), masks, wd, "cleaned")
        check("3 N / IUPAC / lower case removed like the Python code", p1 == p2 and s1 == s2)

        # 4. a short read keeps its index
        with_short = list(reads)
        with_short.insert(2, "ACGTACGTAC")  # record 3 (1-based), shorter than max_k
        write_fasta(os.path.join(wd, "short.fasta"), with_short)
        p, s = run(binary, os.path.join(wd, "short.fasta"), masks, wd, "short")
        check("4 short read keeps numbering", p == shift_pairs(ref_pairs, 3))
        check("4 short read is never in a pair", all(3 not in (int(l.split("\t")[0]), int(l.split("\t")[1])) for l in p))

        # 5. thread count does not change the result (N-heavy input, many OpenMP threads)
        heavy = []
        rng = random.Random(3)
        for r in reads:
            heavy.append("".join("N" if rng.random() < 0.02 else c for c in r))
        write_fasta(os.path.join(wd, "heavy.fasta"), heavy)
        p1, s1 = run(binary, os.path.join(wd, "heavy.fasta"), masks, wd, "heavy1", n_cpu=1)
        p8, s8 = run(binary, os.path.join(wd, "heavy.fasta"), masks, wd, "heavy8", n_cpu=8)
        check("5 --n_cpu 1 == --n_cpu 8 on N-heavy input", p1 == p8 and s1 == s8)

    print()
    if failures:
        print("FAILED: " + ", ".join(failures))
        sys.exit(1)
    print("ALL PASSED")


if __name__ == "__main__":
    main()
