"""Shared evaluation code for LexicHash experiments.

One place for the ground-truth threshold, so every script and every table uses the same
definition (resolve-list item 4).

Ground truth
------------
A read pair is a true overlap when its alignment fraction

    theta = a / (l1 + l2 - a)        (a = alignment length, l1 and l2 = read lengths)

is at least THETA = 0.2. This is the definition used in the LexicHash paper (Bioinformatics
2023): "a threshold of 0.2 for the alignment fraction". The comparison is ">=".

Scoring
-------
A scorer gives each read pair a number; higher means "more likely to overlap". Pairs that were
not scored (not in the candidate list) count as score 0, i.e. they are ranked last but still
count as negatives or missed positives. All C(n, 2) read pairs are therefore always included,
so ROC/PRC are not inflated by looking only at the pairs the method returned (this is also what
get_overlaps() in the original utils/metric_utils.py does by padding with zeros).
Pairs are unordered and the strand is ignored (the best score over both strands is used).
"""
import gzip
import numpy as np

THETA = 0.2  # alignment-fraction threshold used for every evaluation


def read_lengths(fasta_path):
    """Lengths of all records of a FASTA/FASTQ file (plain or .gz), in file order."""
    import pyfastx
    return np.array([len(rec[1]) for rec in pyfastx.Fastx(fasta_path)], dtype=np.int64)


def load_true_pairs(gt_path, seq_lens, theta=THETA):
    """Set of (i, j), i < j, 0-based, whose alignment fraction is >= theta.

    gt_path: DALIGNER-style file, columns: i, j, alignment_length[, strand], tab/space separated,
    0-based read numbers. If a pair appears several times the longest alignment is used.
    Self-alignments (i == j, a read aligned to its own reverse complement) are not read pairs and
    are ignored; the P. falciparum ground truth contains 911 of them above the threshold.
    """
    best = {}
    with open(gt_path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 3:
                continue
            i, j, a = int(p[0]), int(p[1]), float(p[2])
            if i == j:
                continue
            key = (i, j) if i < j else (j, i)
            if a > best.get(key, 0.0):
                best[key] = a
    true_pairs = set()
    for (i, j), a in best.items():
        denom = seq_lens[i] + seq_lens[j] - a
        if denom > 0 and a / denom >= theta:
            true_pairs.add((i, j))
    return true_pairs


def pair_index(i, j, n):
    """Index of the unordered pair (i < j) in the flattened upper triangle of an n x n matrix."""
    return i * n - i * (i + 1) // 2 + (j - i - 1)


def true_pair_labels(true_pairs, n_seq):
    """uint8 label array over all unordered pairs, in np.triu_indices(n_seq, k=1) order."""
    y = np.zeros(n_seq * (n_seq - 1) // 2, dtype=np.uint8)
    for i, j in true_pairs:
        y[pair_index(i, j, n_seq)] = 1
    return y


def top_indices(s, k):
    """Indices of the k best scores (ties broken by pair order), working on the non-zero part first.

    Same result as np.argsort(-s, kind="stable")[:k] but much faster when most scores are 0."""
    nz = np.flatnonzero(s > 0)
    if k <= len(nz):
        return nz[np.argsort(-s[nz], kind="stable")[:k]]
    return np.argsort(-s, kind="stable")[:k]


def evaluate_arrays(s, y):
    """Metrics for score array s and label array y over ALL pairs (score 0 = not scored).

    Returns auroc, auprc (average precision), n_true, n_scored (score > 0), recall_scored (share
    of true pairs with score > 0) and precision_at_n_true (precision among the n_true best pairs;
    ties are broken by pair order).

    The many pairs with score 0 are one tie block, so they are passed to scikit-learn as two
    weighted rows (positives and negatives) instead of millions of identical rows. This gives the
    same numbers as using the full arrays, only faster.
    """
    from sklearn.metrics import roc_auc_score, average_precision_score
    s = np.asarray(s, dtype=np.float64)
    n_true = int(y.sum())
    nz = s > 0
    pos0 = int(y[~nz].sum())
    neg0 = int((~nz).sum()) - pos0
    labels = np.concatenate([y[nz], [1, 0]])
    vals = np.concatenate([s[nz], [0.0, 0.0]])
    w = np.concatenate([np.ones(int(nz.sum())), [pos0, neg0]])
    top = top_indices(s, n_true)
    return {
        "auroc": float(roc_auc_score(labels, vals, sample_weight=w)),
        "auprc": float(average_precision_score(labels, vals, sample_weight=w)),
        "n_true": n_true,
        "n_scored": int(nz.sum()),
        "recall_scored": float((nz & (y == 1)).sum() / n_true) if n_true else float("nan"),
        "precision_at_n_true": float(y[top].sum() / max(1, n_true)),
    }


def evaluate(scores, true_pairs, n_seq):
    """Same as evaluate_arrays, for a dict {(i, j): score} with i < j (0-based); missing pairs = 0."""
    s = np.zeros(n_seq * (n_seq - 1) // 2)
    for (i, j), v in scores.items():
        s[pair_index(i, j, n_seq)] = v
    return evaluate_arrays(s, true_pair_labels(true_pairs, n_seq))
