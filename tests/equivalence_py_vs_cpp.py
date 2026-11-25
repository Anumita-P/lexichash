#!/usr/bin/env python3
"""Check that the C++ port and the Python reference give identical results.

Both programs are fed the SAME masks (written to a mask file), then we compare:
  1. the sketch matrix (every min-hash, for reads and reverse complements), and
  2. the pairs found by the prefix-tree step, with their match lengths.

With identical masks the two implementations should agree exactly, including
which pairs survive the alpha * n_seq cap, because both visit masks, depths
and groups in the same order.

Usage (from the repository root):
  g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native
  python3 tests/equivalence_py_vs_cpp.py --cpp-bin ./lexichash \
      --fasta LexicHash-main/data/NCTC1080/NCTC1080_reads.fasta.gz \
      --n-hash 100 --alpha 1

Use --alpha 0 for "all pairs" (Python alpha=None, C++ alpha=0).
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cpp-bin", required=True, help="compiled C++ binary")
    p.add_argument("--fasta", required=True, help="reads (.fasta or .fasta.gz)")
    p.add_argument("--n-hash", type=int, default=100)
    p.add_argument("--max-k", type=int, default=32)
    p.add_argument("--min-k", type=int, default=14)
    p.add_argument("--alpha", type=float, default=1.0, help="0 means all pairs")
    p.add_argument("--seed", type=int, default=1, help="seed for the shared masks")
    p.add_argument("--n-cpu", type=int, default=2)
    p.add_argument("--no-rc", action="store_true", help="skip reverse complements")
    p.add_argument("--py-root", default=os.path.join(REPO, "LexicHash-main"),
                   help="folder that contains the Python package LexicHash/")
    return p.parse_args()


def write_plain_fasta(src, dst):
    opener = gzip.open if src.endswith(".gz") else open
    with opener(src, "rt") as fin, open(dst, "w") as fout:
        shutil.copyfileobj(fin, fout)


def read_cpp_sketches(path, n_hash):
    rows = []
    with open(path) as fh:
        for line in fh:
            _, _, rest = line.partition(":")
            vals = rest.split()
            rows.append([int(v) for v in vals] if vals else [0] * n_hash)
    return np.array(rows, dtype=np.uint64)


def read_cpp_pairs(path):
    pairs = {}
    with open(path) as fh:
        for line in fh:
            i1, i2, k, strand = line.split()
            pairs[(int(i1) - 1, int(i2) - 1, strand)] = int(k)
    return pairs


def main():
    a = parse_args()
    rc = not a.no_rc
    work = tempfile.mkdtemp(prefix="lh_equiv_")
    print(f"work folder: {work}")

    # 1. shared masks
    masks = np.random.RandomState(a.seed).randint(4, size=(a.n_hash, a.max_k))
    mask_path = os.path.join(work, "masks.txt")
    with open(mask_path, "w") as fh:
        for m in masks:
            fh.write("".join(str(int(x)) for x in m) + "\n")

    # 2. Python reference
    sys.path.insert(0, a.py_root)
    import pyfastx
    from LexicHash import LexicHash as LH

    seqs = pyfastx.Fastx(a.fasta)
    py_sketches, n_seq = LH.sketching(seqs, masks, n_hash=a.n_hash, n_cpu=a.n_cpu, rc=rc,
                                      max_k=a.max_k, min_k=a.min_k)
    alpha_py = a.alpha if a.alpha > 0 else None
    py_pairs = LH.pairwise_comparison(py_sketches, n_seq, n_hash=a.n_hash, n_cpu=a.n_cpu,
                                      alpha=alpha_py, max_k=a.max_k, min_k=a.min_k)
    print(f"python: {n_seq} reads, sketch matrix {py_sketches.shape}, {len(py_pairs)} pairs")

    # 3. C++
    fasta = os.path.join(work, "reads.fasta")
    write_plain_fasta(a.fasta, fasta)
    cmd = [a.cpp_bin, "--out", work, "--fasta", fasta, "--n_hash", str(a.n_hash),
           "--alpha", str(a.alpha), "--n_cpu", str(a.n_cpu), "--min_k", str(a.min_k),
           "--max_k", str(a.max_k), "--mask_file", mask_path]
    if not rc:
        cmd.append("--no_rc")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-2000:], res.stderr[-2000:])
        sys.exit("C++ run failed")
    cpp_sketches = read_cpp_sketches(os.path.join(work, "sketches.txt"), a.n_hash)
    cpp_pairs = read_cpp_pairs(os.path.join(work, "pairwise_results.txt"))
    print(f"c++:    sketch matrix {cpp_sketches.shape}, {len(cpp_pairs)} pairs")

    # 4. compare
    ok = True
    if py_sketches.shape != cpp_sketches.shape:
        print(f"FAIL sketch shapes differ: {py_sketches.shape} vs {cpp_sketches.shape}")
        ok = False
    else:
        diff = py_sketches != cpp_sketches
        n_bad = int(diff.sum())
        print(f"sketch cells that differ: {n_bad} of {diff.size}"
              f" ({int(diff.any(axis=1).sum())} reads affected)")
        ok &= n_bad == 0

    only_py = set(py_pairs) - set(cpp_pairs)
    only_cpp = set(cpp_pairs) - set(py_pairs)
    both = set(py_pairs) & set(cpp_pairs)
    len_diff = [k for k in both if py_pairs[k] != cpp_pairs[k]]
    print(f"pairs only in python: {len(only_py)}, only in c++: {len(only_cpp)}, "
          f"in both: {len(both)}, match length differs: {len(len_diff)}")
    ok &= not only_py and not only_cpp and not len_diff

    print("RESULT:", "IDENTICAL" if ok else "DIFFERENT")
    shutil.rmtree(work, ignore_errors=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
