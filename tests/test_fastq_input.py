#!/usr/bin/env python3
"""Tests for FASTQ input of the C++ program (batch 6).

The C++ reader used to understand FASTA only. It now detects the format from the first
non-empty line ('>' = FASTA, '@' = FASTQ). Every check builds a FASTQ file, runs the binary with
the SAME mask file as on the equivalent FASTA, and compares sketches.txt and pairwise_results.txt.

  1. plain FASTQ            == FASTA with the same sequences
  2. gzip FASTQ             == FASTA
  3. CRLF + no final newline== FASTA
  4. tricky quality lines   : quality strings that start with '@' or '+', and that are as long as the read,
                              must not be taken for headers / separators
  5. N, IUPAC, lower case   : cleaned exactly as for FASTA (output equals the run on the cleaned read)
  6. multi-line sequence    : a read split over several lines before the '+' line
  7. empty read             : a record with no bases keeps its place in the numbering
                              (output equals the file without it, indices from it onwards shifted by one)
  8. '+' line repeats the header: "+r5" instead of "+" gives the same output

Usage (from the repository root):
  g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native -lz
  python3 tests/test_fastq_input.py --cpp-bin ./lexichash
"""
import argparse
import gzip
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_input_handling import make_reads, run, shift_pairs, write_fasta, write_masks  # noqa: E402


def quality(n, kind):
    if kind == "at":      # starts with '@' (a valid Phred+33 character, value 31)
        return "@" + "I" * (n - 1)
    if kind == "plus":    # starts with '+' (Phred 10)
        return "+" + "5" * (n - 1)
    return "I" * n


def write_fastq(path, reads, eol="\n", final_newline=True, gz=False, qkind="plain", wrap=0, plus_name=False):
    recs = []
    for i, s in enumerate(reads):
        seq_lines = [s] if not wrap else [s[j:j + wrap] for j in range(0, len(s), wrap)] or [""]
        rec = ["@r%d" % i] + seq_lines + [("+r%d" % i) if plus_name else "+"]
        rec.append(quality(len(s), qkind) if s else "")
        recs.append(eol.join(rec))
    text = eol.join(recs)
    if final_newline:
        text += eol
    data = text.encode()
    opener = gzip.open if gz else open
    with opener(path, "wb") as fh:
        fh.write(data)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpp-bin", required=True)
    args = ap.parse_args()
    binary = os.path.abspath(args.cpp_bin)
    failures = []

    def check(name, ok, detail=""):
        print("%-5s %s %s" % ("PASS" if ok else "FAIL", name, detail))
        if not ok:
            failures.append(name)

    reads = make_reads()
    with tempfile.TemporaryDirectory() as wd:
        masks = os.path.join(wd, "masks.txt")
        write_masks(masks)
        fa = os.path.join(wd, "base.fasta")
        write_fasta(fa, reads)
        ref_pairs, ref_sketches = run(binary, fa, masks, wd, "base")
        check("baseline (FASTA) produces pairs", len(ref_pairs) > 0, "(%d pairs)" % len(ref_pairs))

        def same(name, path, expect_pairs=ref_pairs, expect_sketches=ref_sketches):
            p, s = run(binary, path, masks, wd, name.split()[0])
            check(name, p == expect_pairs and s == expect_sketches)

        fq = os.path.join(wd, "a.fastq")
        write_fastq(fq, reads)
        same("1 plain FASTQ == FASTA", fq)

        fqz = os.path.join(wd, "a.fastq.gz")
        write_fastq(fqz, reads, gz=True)
        same("2 gzip FASTQ == FASTA", fqz)

        fqc = os.path.join(wd, "crlf.fastq")
        write_fastq(fqc, reads, eol="\r\n", final_newline=False)
        same("3 CRLF + no final newline", fqc)

        for kind in ("at", "plus"):
            f = os.path.join(wd, "q_%s.fastq" % kind)
            write_fastq(f, reads, qkind=kind)
            same("4 quality line starting with '%s'" % ("@" if kind == "at" else "+"), f)

        # 5. dirty read, same as the FASTA test: N, IUPAC and lower case
        dirty = list(reads)
        d = list(reads[5])
        for pos in (10, 11, 80, 150, 151, 152, 200):
            d[pos] = "N"
        d[40], d[100] = "R", "y"
        for pos in range(120, 135):
            d[pos] = d[pos].lower()
        dirty[5] = "".join(d)
        cleaned = list(reads)
        cleaned[5] = "".join(c for c in dirty[5].upper() if c in "ACGT")
        cf = os.path.join(wd, "clean.fastq")
        write_fastq(cf, cleaned)
        cp, cs = run(binary, cf, masks, wd, "cleanref")
        df = os.path.join(wd, "dirty.fastq")
        write_fastq(df, dirty)
        same("5 N / IUPAC / lower case == cleaned read", df, cp, cs)

        # 6. multi-line sequence
        ml = os.path.join(wd, "wrap.fastq")
        write_fastq(ml, reads, wrap=60)
        same("6 sequence wrapped over several lines", ml)

        # 7. empty read keeps its place
        with_empty = list(reads[:2]) + [""] + list(reads[2:])
        ef = os.path.join(wd, "empty.fastq")
        write_fastq(ef, with_empty)
        p, s = run(binary, ef, masks, wd, "empty")
        shifted_ref = shift_pairs(ref_pairs, 3)
        check("7 empty read keeps its index", p == shifted_ref, "(%d pairs)" % len(p))

        # 8. '+' line repeating the header
        pf = os.path.join(wd, "plusname.fastq")
        write_fastq(pf, reads, plus_name=True)
        same("8 '+name' separator line", pf)

    print("\n%s" % ("ALL PASSED" if not failures else "FAILED: " + ", ".join(failures)))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
