"""Tests for the mask strategies (resolve-list item 9).

A. The inverse-frequency sampler in LexicHash.py: a toy read set with known k-mer counts shows that
   the old code (kept in this file as OLD_SOURCE) is NOT inverse-frequency and the fixed code is.
B. The GC sampler splits at the median GC, so both pools are populated.
C. evaluation/mask_strategies.py: every strategy is reproducible from its seed, returns valid masks,
   and the file written for `lexic_hash --mask_file` has the format the C++ reader expects.
D. The port invfreq_fixed in mask_strategies.py agrees with the fixed function in LexicHash.py.
"""
import ast
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "evaluation"))
import mask_strategies as MS

LH = "experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py"
BASE_TO_INT = {'A': 0, 'T': 1, 'G': 2, 'C': 3}


def load_funcs(source, names):
    tree = ast.parse(source)
    ns = {"np": np, "BASE_TO_INT": BASE_TO_INT}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(compile(ast.Module([node], []), LH, "exec"), ns)
    return {n: ns[n] for n in names}


def current_source():
    return open(os.path.join(ROOT, LH), encoding="utf-8").read()


OLD_SOURCE = """
def get_masks_freq_weighted(n_masks, max_k, sequences):
    from collections import Counter
    seq_strings = [s[1] for s in sequences]
    kmers = []
    for seq in seq_strings:
        if len(seq) < max_k:
            continue
        for i in range(len(seq) - max_k + 1):
            kmers.append(seq[i:i+max_k])
    kmer_counts = Counter(kmers)
    weights = np.array([1 / (kmer_counts[k] + 1) for k in kmers])
    probs = weights / weights.sum()
    sampled = np.random.choice(kmers, size=n_masks, p=probs)
    masks = np.array([[BASE_TO_INT[c] for c in k] for k in sampled])
    return masks

def get_masks_cluster_gc(n_masks, max_k, sequences):
    low_gc, high_gc = [], []
    for _, seq in sequences:
        if not seq:
            continue
        gc_content = (seq.count('G') + seq.count('C')) / len(seq)
        if gc_content > 0.5:
            high_gc.append(seq)
        else:
            low_gc.append(seq)
    masks = []
    while len(masks) < n_masks:
        pool = low_gc if np.random.rand() < 0.5 else high_gc
        if not pool:
            continue
        seq = np.random.choice(pool)
        if len(seq) < max_k:
            continue
        start = np.random.randint(0, len(seq) - max_k + 1)
        mask_str = seq[start:start + max_k]
        mask = [BASE_TO_INT[b] for b in mask_str]
        masks.append(mask)
    return np.array(masks)
"""  # the two functions exactly as they were before resolve-batch4 (logic unchanged, comments removed)


def old_source():
    return OLD_SOURCE


def fraction_from_unique_group(masks, group_codes, k):
    inv = {v: k_ for k_, v in BASE_TO_INT.items()}
    n = 0
    for m in masks:
        if "".join(inv[int(x)] for x in m) in group_codes:
            n += 1
    return n / len(masks)


def toy_reads(k):
    rng = np.random.default_rng(3)
    r1 = "".join(rng.choice(list("ACGT"), 40))
    r3 = "".join(rng.choice(list("ACGT"), 40))
    kmers = lambda s: {s[i:i + k] for i in range(len(s) - k + 1)}
    only1, only3 = kmers(r1) - kmers(r3), kmers(r3) - kmers(r1)
    assert len(only1) == 33 and len(only3) == 33
    return [("a", r1), ("b", r1), ("c", r3)], only3


def test_invfreq():
    k = 8
    seqs, group_single = toy_reads(k)       # r1's k-mers have count 2, r3's have count 1
    # P(k-mer drawn is a count-1 k-mer): inverse-frequency  = 33*1 / (33*1 + 33*0.5) = 2/3
    #                                    old (c/(c+1))      = 33*(1/2) / (33*(1/2) + 33*2*(1/3)) = 3/7
    old = load_funcs(old_source(), ["get_masks_freq_weighted"])["get_masks_freq_weighted"]
    new = load_funcs(current_source(), ["get_masks_freq_weighted"])["get_masks_freq_weighted"]
    np.random.seed(0)
    f_old = fraction_from_unique_group(old(6000, k, seqs), group_single, k)
    np.random.seed(0)
    f_new = fraction_from_unique_group(new(6000, k, seqs), group_single, k)
    print("share of masks that are count-1 k-mers:  old %.3f (3/7 = 0.429)   fixed %.3f (2/3 = 0.667)" % (f_old, f_new))
    assert abs(f_old - 3 / 7) < 0.02
    assert abs(f_new - 2 / 3) < 0.02
    return new, seqs, group_single, k


def test_gc():
    rng = np.random.default_rng(5)
    reads = []
    for i in range(200):        # 194 reads with GC about 0.3, 6 reads with GC 0.58 (NCTC1080: 6 of 1,003 above 0.5)
        p = 0.58 if i < 6 else rng.uniform(0.25, 0.38)
        reads.append((str(i), "".join(rng.choice(list("ATGC"), 300, p=[(1 - p) / 2, (1 - p) / 2, p / 2, p / 2]))))
    old = load_funcs(old_source(), ["get_masks_cluster_gc"])["get_masks_cluster_gc"]
    new = load_funcs(current_source(), ["get_masks_cluster_gc"])["get_masks_cluster_gc"]
    inv = {v: k_ for k_, v in BASE_TO_INT.items()}
    rich_text = "\n".join(s for _, s in reads[:6])

    def share_from_rich(masks):
        return float(np.mean(["".join(inv[int(x)] for x in m) in rich_text for m in masks]))

    np.random.seed(1)
    s_old = share_from_rich(old(600, 32, reads))
    np.random.seed(1)
    s_new = share_from_rich(new(600, 32, reads))
    print("share of masks taken from the 6 GC-rich reads (3%% of reads):  old %.2f   fixed %.2f" % (s_old, s_new))
    assert s_old > 0.4          # about half of all masks come from 6 reads
    assert s_new < 0.15         # spread over the reads of the upper half


def test_ports(new_fn, seqs, group_single, k):
    # D: the port in mask_strategies.py gives the same share of count-1 k-mers as the fixed function
    reads = [s for _, s in seqs]
    masks, info = MS.make_masks("invfreq_fixed", reads, 6000, k, seed=0)
    f = fraction_from_unique_group(masks, group_single, k)
    print("port invfreq_fixed: share of count-1 k-mers %.3f (2/3)" % f)
    assert abs(f - 2 / 3) < 0.02
    # C: reproducible, valid, file format
    rng = np.random.default_rng(1)
    big = ["".join(rng.choice(list("ACGT"), 400)) for _ in range(60)]
    big = big + big[:30]            # repeated reads, so that solid (count >= 2) k-mers exist
    import tempfile
    for name in MS.STRATEGIES:
        a, _ = MS.make_masks(name, big, 50, 32, seed=7)
        b, _ = MS.make_masks(name, big, 50, 32, seed=7)
        c, _ = MS.make_masks(name, big, 50, 32, seed=8)
        assert a.shape == (50, 32) and a.max() <= 3
        assert np.array_equal(a, b), name + " is not reproducible"
        assert not np.array_equal(a, c), name + " ignores the seed"
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.txt")
            MS.write_mask_file(p, a)
            lines = open(p).read().split("\n")[:-1]
            assert len(lines) == 50 and all(len(l) == 32 and set(l) <= set("0123") for l in lines)
    print("all %d strategies: reproducible, valid, file format ok" % len(MS.STRATEGIES))


if __name__ == "__main__":
    new_fn, seqs, grp, k = test_invfreq()
    test_gc()
    test_ports(new_fn, seqs, grp, k)
    print("OK")
