"""Mask-generation strategies, one function per strategy, reproducible and testable.

Resolve-list item 9. The original strategies live in
`experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py` (get_masks_uniform_blend,
get_masks_freq_weighted, get_masks_cluster_gc, get_masks_mutated_data,
get_masks_progressive_refresh). Nothing in the repo calls them, and they use the global numpy
random state, so the results in `othermetric/` cannot be re-created from the code. This module
makes each strategy a function of (reads, n_masks, max_k, rng) so that every run can be repeated
from a seed and the masks written to the file format read by `lexic_hash --mask_file`.

Strategies
----------
random          iid uniform mask (the original LexicHash).
blend           half random, half a random 32-mer of a random read (original code, ported as is).
invfreq_orig    the original "inverse frequency" sampler, ported as is. Every k-mer OCCURRENCE gets
                weight 1/(count+1), so a k-mer seen c times is drawn with probability proportional
                to c/(c+1): almost uniform over occurrences, NOT inverse-frequency (bug, see log).
invfreq_fixed   probability of a distinct k-mer proportional to 1/count (what the report describes).
solid           uniform over distinct k-mers that occur at least twice (control: k-mers that exist
                in more than one read, so they can support an overlap).
gc_orig         original GC-aware sampler: reads split at GC 0.5, pools chosen with probability 0.5.
gc_median       same, but split at the median read GC so both pools are populated.
mutated         random 32-mer of a random read, every base mutated with probability 0.1 (original).
buffered        70% random masks, 30% drawn (with replacement) from one 32-mer of each of the first
                300 reads (original code; the report calls it "rolling buffer", it is not).

Base codes: A=0, T=1, G=2, C=3 (same as BASE_TO_INT in LexicHash.py and as the C++ code).
"""
import gzip
import numpy as np

BASE_TO_INT = {"A": 0, "T": 1, "G": 2, "C": 3, "a": 0, "t": 1, "g": 2, "c": 3}
LUT = np.full(256, 255, dtype=np.uint8)
for _b, _v in BASE_TO_INT.items():
    LUT[ord(_b)] = _v


def read_fasta(path):
    """List of sequences (str) of a FASTA file, plain or .gz (multi-line records supported)."""
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, cur = [], []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            elif line.strip():
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def encode(seq):
    """str -> uint8 codes 0..3 (255 for a non-ACGT character)."""
    return LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def window(codes, start, k):
    return codes[start:start + k]


def _usable(reads, max_k):
    """Indices of reads that are long enough and contain only ACGT."""
    out = []
    for i, s in enumerate(reads):
        if len(s) >= max_k and (encode(s) < 4).all():
            out.append(i)
    return out


def _random_window(rng, enc, ids, max_k):
    i = ids[rng.integers(len(ids))]
    start = rng.integers(0, len(enc[i]) - max_k + 1)
    return enc[i][start:start + max_k], i


def _kmer_codes(enc_read, k):
    """uint64 code of every k-mer of one encoded read (k <= 32), first base in the top bits."""
    n = len(enc_read) - k + 1
    v = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        v = (v << np.uint64(2)) | enc_read[j:j + n].astype(np.uint64)
    return v


def _decode_kmer(code, k):
    out = np.empty(k, dtype=np.uint8)
    for j in range(k):
        out[k - 1 - j] = (int(code) >> (2 * j)) & 3
    return out


def _kmer_table(enc, ids, k):
    """All k-mer occurrences of the usable reads -> (distinct codes, counts, occurrence codes)."""
    allk = np.concatenate([_kmer_codes(enc[i], k) for i in ids])
    uniq, inv, counts = np.unique(allk, return_inverse=True, return_counts=True)
    return uniq, counts, allk, inv


# ----------------------------------------------------------------------------- strategies
def random_masks(reads, n, k, rng):
    return rng.integers(0, 4, size=(n, k)).astype(np.uint8), {}


def blend(reads, n, k, rng):
    enc = [encode(s) for s in reads]
    ids = _usable(reads, k)
    n_rand = n // 2
    rand = rng.integers(0, 4, size=(n_rand, k)).astype(np.uint8)
    data, src = [], []
    while len(data) < n - n_rand:
        w, i = _random_window(rng, enc, ids, k)
        data.append(w)
        src.append(i)
    return np.vstack([rand, np.array(data)]), {"source_reads": len(set(src))}


def invfreq_orig(reads, n, k, rng):
    enc = [encode(s) for s in reads]
    ids = _usable(reads, k)
    uniq, counts, allk, inv = _kmer_table(enc, ids, k)
    w = 1.0 / (counts[inv] + 1.0)            # one weight per OCCURRENCE (as in the original)
    p = w / w.sum()
    pick = rng.choice(len(allk), size=n, p=p)
    masks = np.array([_decode_kmer(allk[j], k) for j in pick])
    return masks, {"mean_count_of_picked": float(counts[inv[pick]].mean())}


def invfreq_fixed(reads, n, k, rng):
    enc = [encode(s) for s in reads]
    ids = _usable(reads, k)
    uniq, counts, allk, inv = _kmer_table(enc, ids, k)
    w = 1.0 / counts.astype(np.float64)       # one weight per DISTINCT k-mer
    p = w / w.sum()
    pick = rng.choice(len(uniq), size=n, p=p)
    masks = np.array([_decode_kmer(uniq[j], k) for j in pick])
    return masks, {"mean_count_of_picked": float(counts[pick].mean())}


def solid(reads, n, k, rng):
    enc = [encode(s) for s in reads]
    ids = _usable(reads, k)
    uniq, counts, allk, inv = _kmer_table(enc, ids, k)
    keep = np.where(counts >= 2)[0]
    pick = keep[rng.integers(len(keep), size=n)]
    masks = np.array([_decode_kmer(uniq[j], k) for j in pick])
    return masks, {"mean_count_of_picked": float(counts[pick].mean()), "n_solid": int(len(keep))}


def _gc_sampler(reads, n, k, rng, split):
    enc = [encode(s) for s in reads]
    ids = _usable(reads, k)
    gc = {i: float(((enc[i] == 2) | (enc[i] == 3)).mean()) for i in ids}
    thr = 0.5 if split is None else float(np.median(list(gc.values())))
    high = [i for i in ids if gc[i] > thr]
    low = [i for i in ids if gc[i] <= thr]
    masks, src = [], []
    while len(masks) < n:
        pool = low if rng.random() < 0.5 else high
        if not pool:
            continue
        w, i = _random_window(rng, enc, pool, k)
        masks.append(w)
        src.append(i)
    return np.array(masks), {"gc_threshold": thr, "reads_high": len(high), "reads_low": len(low),
                             "source_reads": len(set(src))}


def gc_orig(reads, n, k, rng):
    return _gc_sampler(reads, n, k, rng, None)


def gc_median(reads, n, k, rng):
    return _gc_sampler(reads, n, k, rng, "median")


def mutated(reads, n, k, rng, rate=0.1):
    enc = [encode(s) for s in reads]
    ids = _usable(reads, k)
    out = []
    while len(out) < n:
        w, _ = _random_window(rng, enc, ids, k)
        w = w.copy()
        for j in range(k):
            if rng.random() < rate:
                w[j] = rng.choice([b for b in range(4) if b != w[j]])
        out.append(w)
    return np.array(out), {}


def buffered(reads, n, k, rng, buffer_size=300):
    enc = [encode(s) for s in reads]
    buf = []
    for i in range(min(buffer_size, len(reads))):
        if len(reads[i]) < k or not (enc[i] < 4).all():
            continue
        start = rng.integers(0, len(enc[i]) - k + 1)
        buf.append(enc[i][start:start + k])
    n_buf = int(n * 0.3)
    rand = rng.integers(0, 4, size=(n - n_buf, k)).astype(np.uint8)
    pick = rng.integers(len(buf), size=n_buf)
    return np.vstack([rand, np.array(buf)[pick]]), {"buffer_size": len(buf),
                                                    "distinct_from_buffer": int(len(set(pick.tolist())))}


STRATEGIES = {
    "random": random_masks, "blend": blend, "invfreq_orig": invfreq_orig, "invfreq_fixed": invfreq_fixed,
    "solid": solid, "gc_orig": gc_orig, "gc_median": gc_median, "mutated": mutated, "buffered": buffered,
}


def make_masks(name, reads, n, k, seed):
    rng = np.random.default_rng(seed)
    masks, info = STRATEGIES[name](reads, n, k, rng)
    masks = np.asarray(masks, dtype=np.uint8)
    assert masks.shape == (n, k) and masks.max() <= 3
    info["unique_masks"] = int(len({m.tobytes() for m in masks}))
    return masks, info


def write_mask_file(path, masks):
    """One mask per line, characters 0-3: the format of `lexic_hash --mask_file`."""
    with open(path, "w", newline="\n") as fh:
        for m in masks:
            fh.write("".join(str(int(x)) for x in m) + "\n")
