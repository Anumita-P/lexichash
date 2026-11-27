# Resolve log

Running record of every code change made while resolving the open items before the arXiv paper.
One section per batch (= one branch / pull request). For each change: the file, the exact lines,
whether it was added, edited or removed, what it does, and why.

**How to read the line numbers.** They refer to the file as it is at the end of that batch's branch
(the "new" numbering). Where a line was replaced, the old line number is given as well
(`old L242 -> new L242`). Later batches can shift these numbers, so use the batch's branch or commit to look them up.

**Change types:** ADD = new lines, EDIT = existing lines changed, DEL = lines removed, NEW FILE = whole file is new.

---

## Batch 1 - branch `resolve-batch1`
Commit `0289628` (code). Resolve-list items covered: 1 (same-mask equivalence test), 5 (exp-decay loop bug), 6 (alpha cast bug).

### 1.1 Python: `alpha` cap must be an integer (5 copies)
Function: `pairwise_comparison`. The line `n_pairs = n_seq * alpha if alpha is not None else ...` became
`n_pairs = int(n_seq * alpha) if alpha is not None else ...`.

**Why:** with a fractional alpha (for example alpha = 2.5 and an odd number of sequences) `n_seq * alpha` is a float.
That float is later used as a cap on the number of pairs, so the cap did not behave correctly for fractional products.
Casting to `int` makes the cap an exact integer count in every case.

| File | Lines | Type |
|---|---|---|
| `LexicHash-main/LexicHash/LexicHash.py` | L242 (old L242) | EDIT |
| `LexicHash_1/LexicHash/LexicHash/LexicHash.py` | L240 (old L240) | EDIT |
| `LexicHash_1orig/LexicHash/LexicHash/LexicHash.py` | L240 (old L240) | EDIT |
| `experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py` | L405 (old L405) | EDIT |
| `mycode/le.py` | L231 (old L231) | EDIT |

### 1.2 Python: exp-decay score looped over the wrong range
File: `experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py`

| Lines | Type | Change | Why |
|---|---|---|---|
| L574 (old L574) | EDIT | `def expdecayscore(y_score, decay)` -> `def expdecayscore(y_score, decay, top_k=100)` | The 100 was hard-coded in the loop. It is now a parameter so the scorer can be varied in the ablation. |
| L575-L579 | ADD | Docstring explaining the score: sort the match lengths in descending order and weight the i-th largest by `(1 - decay) ** i`. | The score was undocumented. |
| old L577-L582 -> new L581-L583 | EDIT / DEL | Removed the debug comment, `n = len(...)` and blank lines, and replaced `for i in range(max(n,100))` with `for i in range(min(len(revsorted_y_score), top_k))`. | `max(n,100)` runs at least 100 iterations, so it crashed with an `IndexError` whenever there were fewer than 100 hashes. It also never took only the top-k as intended. `min` is what was meant. |
| L586 | ADD | Blank line before the next import. | Formatting only. |
| L589 (old L587) | EDIT | `def refine_match_pairs(pair_match_lens, sketches)` -> `def refine_match_pairs(pair_match_lens, sketches, decay=0.1, top_k=100)` | `decay` and `top_k` were hard-coded inside; the new defaults keep the old behaviour (decay 0.1). |
| L607 (old L605) | EDIT | `expdecayscore(y_score, 0.1)` -> `expdecayscore(y_score, decay, top_k)` | Passes the new parameters through. |

### 1.3 C++: seeded masks and mask files, so C++ and Python can use identical masks
File: `mycode/lexic_hash.cpp` (CRLF line endings kept as in the original)

| Lines | Type | Change | Why |
|---|---|---|---|
| L48-L50 | ADD | Globals `long long mask_seed = -1;` and `string mask_file_path;` | Hold the new `--seed` and `--mask_file` options. Seed -1 means "random as before". |
| old L51-L52 -> new L54-L61 | EDIT | In `get_masks`: `mt19937 gen(rd())` became `mt19937 gen;` followed by an if/else. If `mask_seed >= 0` the generator is seeded with it, otherwise with `random_device` as before. | Makes runs reproducible. Without a seed the behaviour is unchanged. |
| L73-L119 | ADD | New functions `load_masks(path, max_k, masks)` and `save_masks(path, masks)` and the comment describing the format (one mask per line, `max_k` characters of 0-3). `load_masks` strips trailing CR/space, skips blank lines, and exits with an error for a wrong line length or a character outside 0-3. | The Python and C++ programs draw random numbers differently, so the only way to compare them is to give both the same masks. |
| L756-L763 | ADD | Parsing for `--seed N` and `--mask_file PATH` in `main`. | Exposes the above on the command line. |
| L786-L798 | ADD | Validation `1 <= max_k <= 32` and `1 <= min_k <= max_k`, then `MIN_K = min_k; MAX_K = max_k;` | **Bug fix.** The sketching code reads the upper-case globals `MIN_K`/`MAX_K` (14 and 32), so `--min_k` and `--max_k` had no effect on hashing and gave silently inconsistent results. `max_k > 32` cannot work because k-mers are packed into 64 bits, so it is now rejected. |
| old L753 -> new L830-L844 | EDIT | `masks = get_masks(n_hash, max_k)` became: load masks from `--mask_file` if the file exists (and set `n_hash` to the number loaded); otherwise generate them and, if a path was given, save them to that file. | Lets one run write the masks and the other read the same file. |

### 1.4 New test: same masks, same results?
| File | Lines | Type | What / why |
|---|---|---|---|
| `tests/equivalence_py_vs_cpp.py` | L1-L142 | NEW FILE | Writes one shared mask file, runs the Python sketching and pairwise comparison and the C++ program on it, then compares the sketch matrices and the selected pair sets. Prints IDENTICAL or DIFFERENT. This is the evidence that the C++ port is a faithful reimplementation, which the paper needs. |
| `tests/README.md` | L1-L43 | NEW FILE | How to run the test, the new flags, and the results table (identical on NCTC1080 in every configuration tried). |

### 1.5 Not changed in this batch
- The other copies of `LexicHash.py` do not contain `expdecayscore`, so 1.2 applies only to the experiment copy.
- No change to the scoring, the pair selection, or the output format.

---

## Batch 2 - branch `resolve-batch2` (starts from `resolve-batch1`)
Resolve-list item covered: 7 (C++ input handling: gzip, N bases, short reads, shared-map race).
Line numbers below are for `mycode/lexic_hash.cpp` on this branch. Batch 1 additions earlier in the file moved down (for example `get_masks` is now at L98).

### 2.1 C++: read gzip-compressed FASTA directly
| Lines | Type | Change | Why |
|---|---|---|---|
| L3 | ADD | `#include <zlib.h>` | `gzopen` / `gzgets` read both plain and gzip files with one code path. |
| old L18-L44 -> new L33-L93 | EDIT | `read_fasta` rewritten: `ifstream` replaced by `gzFile fp = gzopen(...)` (L41), a 1 MiB read buffer (`gzbuffer`, L46) and a line reader built on `gzgets` (L71-L82) that also copes with long lines, CRLF endings (L50-L51) and a last line with no newline (L81-L82). | The Python code reads `.fasta.gz` directly (`pyfastx`), but the C++ program could only read plain text, so a `.gz` file gave garbage and users had to decompress by hand. Tested: the gzip output is identical to the plain output. |
| `README.md` L9 | EDIT | Added `-lz` to the build command. | zlib must be linked now. |

### 2.2 C++: N and other non-ACGT characters
| Lines | Type | Change | Why |
|---|---|---|---|
| L37-L38, L64-L68 | EDIT (inside the rewritten `read_fasta`) | Each character goes through `BASE_CODE`; lower case is upper-cased and anything that is not A, C, G or T is dropped. | **Bug fix.** Before, an `N` (or any IUPAC code) was looked up with `BASE_TO_INT[c]`, which silently inserted the key with value 0, so every N became an A and created false matches. The Python code (`clean_seqs` in `seq_utils.py`) deletes these characters, so now both agree. |
| L59-L62 | ADD | A sequence line that appears before any `>` header starts a record. | Keeps the old behaviour for header-less files. |

### 2.3 C++: data race on the shared base table
| Lines | Type | Change | Why |
|---|---|---|---|
| L19-L31 | ADD | New read-only `BASE_CODE` array (256 entries, -1 for non-bases). The old `BASE_TO_INT` map (L16-L17) is kept but no longer used for sketching. | **Bug fix.** `to_int_seq` ran inside the OpenMP loop in `sketching()` and called `BASE_TO_INT[c]`. `operator[]` on an `unordered_map` inserts when the key is missing, so several threads could write to the same global map at once when they met an unknown character (undefined behaviour, can crash or corrupt results). A plain array is only read, so it is thread-safe. |
| old L154 -> new L200-L209 | EDIT | `to_int_seq` uses `BASE_CODE` and calls `abort()` with a message if it ever sees a non-base. | After 2.2 this cannot happen, and failing loudly is better than silently treating it as A. |

### 2.4 C++: short reads no longer shift the read numbers
| Lines | Type | Change | Why |
|---|---|---|---|
| old L31 and L41 -> removed | DEL (inside `read_fasta`) | Removed the `if (current_seq.size() >= max_k)` filters that dropped short reads before numbering. | **Bug fix.** The output file lists reads by 1-based number. After a short read was dropped, every later read had the wrong number, so the pairs no longer matched the FASTA record numbers or the ground-truth files used for ROC/PRC evaluation. |
| L85-L91 | ADD | Counts reads shorter than `max_k` (after removing non-ACGT) and prints a warning with the count. | The user should know how many reads are left out. |
| L515 | unchanged | `if (seqs[i].size() < max_k) continue;` in `sketching()` already existed. Short reads now reach this line, get an empty sketch and are skipped; `get_matching_sets` already ignores empty sketches (L684). | A short read keeps its number but never appears in a pair. Note: the Python code instead pads short reads, so it can pair them; this is a known difference for reads shorter than `max_k`. |

### 2.5 New test
| File | Lines | Type | What / why |
|---|---|---|---|
| `tests/test_input_handling.py` | L1-L172 | NEW FILE | Five checks on small synthetic FASTA files: gzip = plain; CRLF and no final newline; N / IUPAC / lower case equal the cleaned read; short read keeps its number (output equals the run without it, shifted by one); `--n_cpu 1` equals `--n_cpu 8` on N-heavy input. |
| `tests/README.md` | L8 (edit), L44-L58 (add) | EDIT / ADD | Added `-lz` to the build line and a section describing the new test and the input rules. |

**Test results.** New binary: all checks pass. The equivalence test against Python (`tests/equivalence_py_vs_cpp.py`, 50 masks, alpha 1) is still IDENTICAL: 0 of 100,300 sketch cells differ, and all 1,003 pairs match. The same input checks on the batch-1 binary: the N check, the short-read check (320 pairs differ) and the gzip check all fail, so the tests catch the bugs.

### 2.6 Not changed / still open
- FASTQ input (the Python code accepts it through `pyfastx`): not supported in C++ yet.
- Short reads are skipped in C++ but padded in Python (see 2.4).

---

## Batch 3 - branch `resolve-batch3` (starts from `main` after batch 2)
Resolve-list items covered: 3 (scorer ablation) and 4 (one 0.2 threshold everywhere). One extra bug found on the way (3.3).
No C++ code was changed in this batch; `mycode/lexic_hash.cpp` is as in batch 2.

### 3.1 One threshold for every evaluation (item 4)
**What was inconsistent.** Four different ways of deciding which read pairs are true overlaps were in use:
`gengraphs.py` passed threshold 0 (every ground-truth alignment counts as true: 5,260 pairs on NCTC1080);
the default of `get_gt_df` / `get_gt_df_nctc` in `metric_utils.py` was 0.3 with a strict `>` (3,335 pairs; this is what
`trymetrics.py` in `LexicHash_1/` and `LexicHash_1orig/` used); the ROC/PRC plot scripts used 0.2 with `>=` (3,857 pairs);
and `test.py` scripts do not threshold at all (top-k overlap with the longest ground-truth alignments).
The LexicHash paper states 0.2 for the alignment fraction `a / (l1 + l2 - a)`, so 0.2 is now used everywhere.
Results from the different scripts could not be compared before.

| File | Lines | Type | Change |
|---|---|---|---|
| `evaluation/lh_eval.py` | L27 | NEW FILE (123 lines) | `THETA = 0.2`, the single definition. The file also holds `load_true_pairs` (L36-L61), `evaluate_arrays` (L87-L115) and the helpers; see 3.2. |
| `LexicHash-main/LexicHash/utils/metric_utils.py`, `LexicHash_1/LexicHash/LexicHash/utils/metric_utils.py`, `LexicHash_1orig/LexicHash/LexicHash/utils/metric_utils.py`, `experiment with lexichash/LexicHash-main/LexicHash/utils/metric_utils.py` | L12, L16, L19, L26 (same lines in all four) | EDIT | L12 and L19: `min_theta=0.3` -> `min_theta=0.2`. L16 and L26: `thetas > min_theta` -> `thetas >= min_theta`. These are copies of the original LexicHash code, so they now differ from upstream here on purpose. |
| `LexicHash-main/LexicHash/gengraphs.py` | L19 | EDIT | `get_gt_df(gt_path, seq_lens, 0)` -> `get_gt_df(gt_path, seq_lens)` (default 0.2). |
| `test accuracy/roc_prc/plot.py` | L25-L28 | EDIT | Replaced the literal `jaccard_threshold = 0.2` by an import of `THETA` (L25-L27) and `jaccard_threshold = THETA` (L28). The value is the same; the source is now shared. |
| `test accuracy/roc_prc/plot copy.py` | L38-L41 | EDIT | Same change. |
| `comparisons/expwala/plot copy.py` | L9-L11 (add), L36 (edit) | EDIT | Import of `THETA`; `if jaccard >= 0.2` -> `if jaccard >= THETA`. |
| `comparisons/multipleruns/plot copy.py` | L9-L11 (add), L39 (edit) | EDIT | Same change. |
| `tests/test_threshold_consistency.py` | L1-L66 | NEW FILE | Fails if any `metric_utils.py` default differs from `THETA`, if a strict `>` is used, if a script passes an explicit threshold to `get_gt_df`, or if a plot script has its own literal. It fails on the old code and passes now. |

**Not changed.** The `test.py` scripts (`test accuracy/`, `comparisons/*/`, `othermetric/`) report "top-k predicted pairs that are among the
top-k longest ground-truth alignments". That is a different metric and has no threshold; keep it separate from ROC/PRC in the report.
`plot.py` computes ROC/PRC only over the predicted pairs, so its curves stop early and the AUC is not comparable with `evaluation/lh_eval.py`
(which counts all pairs). `pr_vs_k(..., thresh=300)` in `metric_utils.py` is the alignment-score threshold from the paper for a different plot and is left as is.

### 3.2 Shared evaluator and scorer ablation (item 3)
| File | Lines | Type | What / why |
|---|---|---|---|
| `evaluation/lh_eval.py` | L1-L123 | NEW FILE | Ground truth (alignment fraction >= 0.2) and metrics: ROC-AUC, PR-AUC, precision among the n_true best pairs. All C(n,2) pairs are included; pairs not returned count as score 0. Self-alignments (a read aligned to itself) are ignored; the P. falciparum ground truth has 911 of them above the threshold. |
| `evaluation/scorer_ablation.py` | L1-L275 | NEW FILE | Compares max y, exp-decay (8 decays x 5 top_k), mean of top-k, count of masks with y >= t, and a likelihood-style sum, on the pairs the pipeline returns (alpha = 1, 2, 5, 10) and on all pairs, with and without reverse complements. Sketches and candidate pairs come from the C++ program. |
| `evaluation/summarize_ablation.py` | L1-L106 | NEW FILE | Turns the CSV into the Markdown report. |
| `evaluation/README.md` | L1-L59 | NEW FILE | How to run it and how to read it. |
| `evaluation/results/scorer_ablation_NCTC1080.{csv,md}` | - | NEW FILES | 5 seeds without reverse complements; reverse-complement settings for seed 1 only (6 times slower). 500 masks. |
| `evaluation/results/scorer_ablation_Pfalciparum.{csv,md}` | - | NEW FILES | 2 seeds, reverse complements on and off, 200 masks, alpha 2, 5, 10 (all-pairs pool skipped for time). |
| `tests/README.md` | L8 and L60-L71 | EDIT / ADD | Documents the new tests. |

**Checks of the ablation code.** Vectorised match lengths: 0 mismatches against `yscore()` in the Python code (6,000 values). Vectorised exp-decay:
max difference 6e-14 against `expdecayscore()`. Evaluator: perfect scores give AUROC and AUPRC 1.0, random scores give 0.51 and the prevalence.
The sparse evaluator gives the same numbers as the dense one on synthetic data with many ties.

### 3.3 Bug: exp-decay re-ranker used the wrong sketch rows for reverse-strand pairs
File: `experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py`

| Lines | Type | Change | Why |
|---|---|---|---|
| L59 | EDIT | `refine_match_pairs(pair_match_lens, sketches)` -> `refine_match_pairs(pair_match_lens, sketches, n_seq=n_seq)` | The function needs the number of reads to find the reverse-complement rows. |
| L589 | EDIT | Added parameter `n_seq=None`. | Same. |
| L595-L598 | ADD | Docstring lines for `sketches` (2*n_seq rows when reverse complements are on) and `n_seq`. | |
| L609-L622 | EDIT | New branch: with reverse complements, a '+' pair uses rows (a, b) and (a', b'), a '-' pair uses rows (a, b') and (a', b) (a' = row of the reverse complement), per mask the larger match. The old single call `yscore(ind1, ind2, sketches)` stays as the else branch (L622). | **Bug fix.** The pair key holds read numbers below n_seq, so the old code always compared the forward rows of the two reads. For a '-' pair (read on one strand, mate on the other) that is the wrong comparison. This is the same rule the pair search uses (`bottom_up_pair_aggregation`). |
| `tests/test_refine_strand.py` L1-L39 | NEW FILE | | Hand-built sketches with a perfect '-' overlap: the pair scores the maximum with the fix and fails on the old code (TypeError, no `n_seq`). The fixed Python function also agrees with the vectorised scorer to 1e-13 on 600 real P. falciparum candidate pairs. |

### 3.4 Results (what the ablation shows)
Numbers are AUPRC, mean +- std over mask seeds, ground truth theta >= 0.2. Full tables: `evaluation/results/*.md`.

**NCTC1080, no reverse complements, 5 seeds, 500 masks.** The pairs the pipeline returns contain 23.7% (alpha 1), 46% (alpha 2), 95.3% (alpha 5) and 98.6% (alpha 10) of the true pairs.

| Pool | max y | exp-decay, as in code (decay 0.1, top_k 100) | exp-decay, best of grid (decay 0.01; top_k all / 100 / 50) | mean of top-k y, best k (100 / 100 / 50) |
|---|---|---|---|---|
| alpha 1 | 0.221 +- 0.002 | 0.240 +- 0.002 | 0.241 +- 0.001 | 0.241 +- 0.001 |
| alpha 5 | 0.825 +- 0.007 | 0.915 +- 0.005 | 0.930 +- 0.004 | 0.932 +- 0.004 |
| all pairs | 0.846 +- 0.006 | 0.948 +- 0.002 | 0.960 +- 0.001 (top_k 50) | 0.961 +- 0.001 (top_k 50) |

**P. falciparum, 2 seeds, 200 masks, alpha 10** (the pool contains 39% to 58% of the true pairs depending on the setting).

| Setting | max y | exp-decay (0.1, 100) | best exp-decay | mean of top-k |
|---|---|---|---|---|
| reverse complements, strand-aware rows (fixed) | 0.124 +- 0.059 | 0.461 +- 0.056 | 0.473 +- 0.059 | 0.470 +- 0.060 |
| reverse complements, forward rows only (old code) | 0.120 +- 0.057 | 0.298 +- 0.057 | 0.317 +- 0.055 | 0.313 +- 0.055 |
| no reverse complements | 0.083 +- 0.045 | 0.270 +- 0.036 | 0.277 +- 0.035 | 0.272 +- 0.033 |

**Conclusions that the data support.**
1. A score that uses many masks beats the original max-of-masks score on both datasets, in every seed: +0.09 AUPRC at alpha 5 on NCTC1080 and +0.34 (strand-aware, alpha 10) on P. falciparum.
2. The exponential weighting itself is not what helps. A plain mean of the top-k matches does as well (within 0.005 AUPRC in every setting), and on NCTC1080 AUPRC rises monotonically as the decay gets smaller (best at the smallest value tried, 0.01, where the weights are nearly flat). The paper should say "aggregating evidence over many masks", not claim a benefit of the exponential form, unless a further experiment shows one. The default `decay = 0.1` is slightly worse than 0.01 on NCTC1080 (-0.015 AUPRC at alpha 5).
3. The gain depends on the pool. At alpha 1 only 24% of the true pairs are inside the pool on NCTC1080, so AUROC cannot move and AUPRC gains only +0.02. The scorer matters most when alpha is 5 or more.
4. The strand fix (3.3) matters on mixed-strand data: +0.16 AUPRC at alpha 10 on P. falciparum. On NCTC1080 it makes no difference because the data are effectively single-strand (all ground-truth lines are '+', and for 97.6% of the true pairs the same-strand sketches give the longer match). Running with reverse complements on a single-strand dataset costs a lot there (seed 1, alpha 5, max y: 0.83 without vs 0.38 with strand-aware rows), because the reverse strand adds chance candidates.
5. The count-of-masks and likelihood-style scores are not better than exp-decay.

**Limits.** "Best of grid" is chosen on the same data. P. falciparum has only 2 seeds with a large spread between them (std about 0.05), the NCTC1080 reverse-complement settings only 1 seed, and 200 masks were used on P. falciparum (500 on NCTC1080). The E. coli dataset for the runtime benchmark would be a third check.

---

## Batch 4 - branch `resolve-batch4` (starts from `main` after batch 3)
Resolve-list items covered: 9 (mask strategies) and 10 (multi-scale).

### 4.0 What was found (before any change)

Mask strategies (`experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py`, L226-L392):
1. **None of the five strategy functions is called anywhere.** `find_overlaps` (L45) always calls `get_masks`
   (random). The results in `othermetric/` were made with hand-edited copies, use the global numpy random state and
   no seed, so they cannot be reproduced from the repo. `othermetric/test.py` also ranks pairs by an integer match
   length with thousands of ties (the `kmax32` files hold 17,000-27,000 pairs, most with length 14-16), so its
   "top-k" numbers depend on the order of lines in the file.
2. **`get_masks_freq_weighted` is not inverse-frequency.** It gives every k-mer *occurrence* the weight
   `1/(count+1)`; a k-mer seen `c` times is therefore drawn with probability proportional to `c/(c+1)`, which is
   close to uniform over occurrences (it even prefers repeated k-mers slightly). On NCTC1080 the mean count of a
   drawn k-mer is 1.00-1.01: with noisy long reads almost every 32-mer occurs once, so even a correct
   inverse-frequency sampler has almost no frequency signal to use there.
3. **`get_masks_cluster_gc` splits reads at a fixed GC of 0.5.** On NCTC1080 only 6 of 1,003 reads are above 0.5
   (mean GC 0.39), so half of all masks came from 6 reads. On P. falciparum (mean read GC 0.24) only 1.2% of the 2,000 reads are above 0.5.
4. **`get_masks_progressive_refresh` is not a rolling buffer** (the report says it is) and not "70% from the full
   dataset": it takes one 32-mer from each of the *first 300 reads*, draws 30% of the masks from them *with
   replacement* (duplicates: on average 468 unique masks of 500) and fills the rest with purely random masks.
5. The report says masks were "implemented in Python" and tested "with nhash = 20 on 2,000 pairs"; the 20-mask
   comparison has a standard deviation of about 0.02 AUPRC between seeds, so single runs cannot rank strategies.

Multi-scale (`LexicHash_1/LexicHash/run_module.py`):
6. **The stages in the code are not the stages in the report.** Report: (127,256), (63,128), (31,64), (14,32).
   Code: (63,128) 700 masks, (31,64) 600, (14,32) 600, (13,16) 400. There is no (127,256) stage, and stage 4
   (13-16) lies inside the range of stage 3 (14-32).
7. **The merge keeps the maximum raw match length over stages**, so scores from different scales (up to 128 vs up
   to 16) are compared directly, and the merged score has only 56-58 distinct values over 502,503 pairs.
8. **Total masks are 2,300 against 500 for the naive method**, so any comparison to naive LexicHash mixes "more
   scales" with "more masks".
9. `merge_prediction_files` loops over every row with `iterrows()` (slow for million-row stage files).

### 4.1 Python: inverse-frequency weights per distinct k-mer
File `experiment with lexichash/LexicHash-main/LexicHash/LexicHash.py`, function `get_masks_freq_weighted`.
EDIT old L277-L278 -> new L277-L281 (comment, `distinct = list(kmer_counts.keys())`, `weights = np.array([1.0 / kmer_counts[k] for k in distinct])`),
EDIT old L281 -> new L284 (`np.random.choice(distinct, ...)`).

**Why:** the weights are now one per distinct k-mer and proportional to `1/count`, which is what the report describes.
Checked by `tests/test_mask_strategies.py` (toy reads with known counts: share of count-1 k-mers 0.422 with the old code against the
expected 3/7, 0.657 with the new code against the expected 2/3).

### 4.2 Python: GC split at the median read GC
Same file, function `get_masks_cluster_gc`. EDIT old L298 -> new L301-L305 (median threshold `threshold = float(np.median(gcs))`) and
EDIT old L303 -> new L310 (`gc_content > threshold`).

**Why:** both pools are populated. Test: with 6 GC-rich reads out of 200, the old code took 51% of its masks from those
6 reads, the new code 3%.

### 4.3 Python: vectorised merge of the stage files
File `LexicHash_1/LexicHash/run_module.py`, function `merge_prediction_files`. EDIT old L47-L56 -> new L47-L55 (CRLF file).
The row loop became `pd.concat` + `groupby(["i1","i2","direction"]).score.max()`. **Why:** same rows, much faster.
`tests/test_merge_stages.py` compares old and new on 4 random stage files (5,785 merged rows, identical).

### 4.4 NEW FILE `evaluation/mask_strategies.py` (227 lines)
All nine strategies (random, blend, invfreq_orig, invfreq_fixed, solid, gc_orig, gc_median, mutated, buffered) as functions
of `(reads, n, max_k, rng)`, reproducible from a seed, plus `write_mask_file` for `lexic_hash --mask_file`.
`invfreq_orig` and `gc_orig` are ports of the old code, kept so the old behaviour can be compared with the fix.
`solid` (uniform over distinct k-mers seen in at least 2 reads) is a control.

### 4.5 NEW FILE `evaluation/mask_ablation.py` (106 lines)
For each strategy and seed it writes the masks, runs the C++ binary with them, scores all pairs and the alpha*n candidate pairs
with max y and exp-decay (0.1, 100) and evaluates with `lh_eval.py`. Results: `evaluation/results/mask_ablation_NCTC1080_h{20,500}.{csv,md}`.

### 4.6 NEW FILE `evaluation/multiscale_eval.py` (218 lines)
Re-creates the four stages of `run_module.py` with exact, vectorised y-values (sketches follow `lexicographic_first`),
including the per-stage `min_k` and alpha*n pool, and the max merge; plus controls with the same number of masks at one scale,
and a pooled exp-decay merge. Results: `evaluation/results/multiscale_NCTC1080.{csv,md}`.

### 4.7 NEW FILES `evaluation/summarize_batch4.py` (69 lines) and tests
`tests/test_mask_strategies.py` (175 lines), `tests/test_multiscale_sketch.py` (71 lines; `lexicographic_first` equals brute force in 160 of 160 cases for
k 16-128, `y_all_pairs` equals the base-by-base prefix in all values), `tests/test_merge_stages.py` (62 lines).

### 4.8 Results

NCTC1080, no reverse complements, ground truth theta >= 0.2, AUPRC, mean +- std over seeds. Full tables with intervals: `evaluation/results/*.md`.

**Mask strategies, 500 masks, 5 seeds (the real operating point).** Every strategy is within 0.01 AUPRC of random masks with max y
and within 0.005 with exp-decay (all intervals for the difference to random contain 0 except `gc_orig` with max y, -0.006 [-0.011, -0.002],
and `buffered` with exp-decay, -0.003 [-0.006, -0.000]).

| score, pool | random | invfreq_orig | invfreq_fixed | gc_orig | gc_median | blend | mutated | buffered | solid |
|---|---|---|---|---|---|---|---|---|---|
| max y, all pairs | 0.852 | 0.851 | 0.845 | 0.845 | 0.846 | 0.843 | 0.848 | 0.847 | 0.832 |
| exp-decay, all pairs | 0.951 | 0.953 | 0.953 | 0.952 | 0.953 | 0.951 | 0.951 | 0.948 | 0.951 |
| exp-decay, alpha 5 pool | 0.922 | 0.920 | 0.920 | 0.920 | 0.917 | 0.920 | 0.922 | 0.917 | 0.903 |

**Mask strategies, 20 masks (the setting of the report), 20 seeds.** Std between seeds is about 0.02, so the ranking is noisy.
With all pairs scored, no strategy beats random (0.455 max y, 0.491 exp-decay); `gc_orig` is worse by 0.025 [0.007, 0.042] (max y) and
0.022 [0.008, 0.036] (exp-decay), `invfreq_orig` by 0.017 / 0.018, and `mutated` by 0.012 [-0.001, 0.024] / 0.013. Only inside the alpha-5 candidate pool do masks made of
k-mers that occur in at least two reads help (`solid`: +0.016 / +0.022 over random, intervals above 0).

**Multi-scale, 3 seeds, 500 masks single-scale baseline against 2,300 masks.**

| method | masks | AUPRC |
|---|---|---|
| multi-scale as coded (max of reported lengths) | 2300 | 0.858 +- 0.004 |
| stage 1 alone (k 63-128) | 700 | 0.010 (returns no pair at all: noisy reads almost never share 63 bases) |
| stage 2 alone (k 31-64) | 600 | 0.395 |
| stage 3 alone (k 14-32) | 600 | 0.851 +- 0.013 |
| stage 4 alone (k 13-16) | 400 | 0.642 |
| single scale k=32, max y | 500 / 2300 | 0.849 / 0.836 |
| single scale k=32, exp-decay (0.1, 100) | 500 / 2300 | 0.949 / 0.941 |
| multi-scale, y/max_k pooled over all masks, exp-decay (0.1, 100) | 2300 | 0.960 +- 0.002 |

**Conclusions that the data support.**
1. The mask strategies do not matter at 500 masks, and the report's statement that the mutated strategy is slightly better than random is
   not supported (-0.003 [-0.008, 0.002] with max y; -0.012 with 20 masks). The honest statement for the paper is "no tested strategy differs from random masks".
2. The idea behind inverse-frequency and GC sampling does not apply to this data: almost every 32-mer occurs once, and one dataset has almost no GC spread.
   Both original versions were also implemented incorrectly (4.0, items 2 and 3); after the fix they are still not better than random.
3. Multi-scale as coded gives 0.858 against 0.849 for naive max-y with 500 masks, but its own k=14-32 stage (600 masks) gives 0.851, i.e. the gain is the extra
   masks of one stage. Stage 1 is empty on this data, stage 2 is poor, the max merge compares match lengths on different scales and creates ties.
   Merging the same information with exp-decay after scaling by max_k gives 0.960 (+0.011 over the best single-scale exp-decay, 0.949), a real but small gain
   that needs the y-values of every mask, which `run_module.py` does not keep.

**Limits.** NCTC1080 only (no reverse complements). P. falciparum was not run: with reverse complements the 9 strategies x seeds would take about 6x longer
and its GC spread is even narrower (1.2% of reads above 0.5). Multi-scale has 3 seeds. The multi-scale pool emulation breaks ties at random, which is how
the original prefix-tree step behaves only on average. The best-of-grid style choice was not used here (one fixed scorer, exp-decay 0.1 / 100).
**Not changed:** the multi-scale max merge in `run_module.py`, the strategy functions that are still not wired into `find_overlaps`
(no `--mask_strategy` option yet); both are decisions for item 11 (what the paper claims).

<!-- Add each new batch below this line using the same layout. -->

---

## Batch 5 - branch `daligner-ground-truth` (resolve item 2, first part: ground truth)
The E. coli runtime benchmark needs a DALIGNER ground truth. The paper gives only the tool (DALIGNER for the NCTC data sets,
minimap2 for P. falciparum) and the threshold (alignment fraction >= 0.2), not the DALIGNER parameters or how the alignment length
`a` is measured.

| File | Lines | Type | What / why |
|---|---|---|---|
| `evaluation/daligner_ground_truth.py` | L1-L180 | NEW FILE | Converts `LAshow` output to the tab-separated `i j a strand` file LexicHash reads (0-based read numbers, `a` in bases). `--calibrate` tries 16 definitions (same-strand only or both; length of the A interval, the B interval, their mean or maximum; longest alignment or sum of the pieces per pair) against the shipped `NCTC1080_daligner_ground_truth.txt` and prints, for each, the number of true pairs (theta >= 0.2), the overlap with the shipped set (Jaccard) and the share of pairs with an identical length. This shows whether the DALIGNER settings reproduce the paper's ground truth before the same recipe is used on E. coli. Tested on 6 pasted `LAshow` lines only (parsing; the first line reproduces `0 22 3395.0 +`). |

Finding while testing: with default DALIGNER settings the pair (reads 1, 307) has a single alignment of 2,402 bases, whereas the shipped file has 3,969.0, so the shipped file was probably made with other settings or another length definition; the calibration run on the full NCTC1080 output decides this.

Result of the first calibration run on the full NCTC1080 DALIGNER output (default `daligner -T4`; 19,032 alignment lines parsed, 7,424 of them reverse-complement):
best definition = same-strand alignments only, `a` = length of the interval on read A, sum over the pieces of a pair: 3,694 of the shipped file's 3,857 true pairs are reproduced
(Jaccard 0.956; 163 missing, 6 extra), median length ratio 1.000, identical length for 52.8 % of the common pairs.
Added in the same branch: `--diagnose` (lists true pairs of the reference that are missing from the new file and says whether DALIGNER found no alignment or a shorter one)
and a warning that lists lines of the `LAshow` output the pattern did not recognise (the first run parsed 19,032 of 19,250 records).
Second calibration step (same branch): `--diagnose` showed that the 218 unparsed lines are those with a read number >= 1,000, which LAshow prints with a thousands separator (`1,002.1`); the pattern now accepts `[\d,]+` for both read numbers (about 1 % of the records, all involving the last 3 reads, so it explains some of the 58 pairs "with no alignment at all"). It also showed that for 105 true pairs the shipped length is clearly longer than the sum of the DALIGNER pieces (for example reads 30 and 459: 9,516 against 2,388), so a third way of combining the pieces was added: `--agg span` (from the start of the first piece to the end of the last piece on each read). The calibration table now has 24 rows.

Calibration result after the parameter sweep (NCTC1080, `daligner -T4`, strand n, a = A interval, best of sum/span):
default 0.966 (3,738 true pairs against 3,857 in the shipped file); `-l500` 0.972 (3,871 true pairs, 53 missing, 67 extra); `-h25` 0.967.
`-e` below 0.70 is rejected by DALIGNER ("Average correlation must be in [.7,1.)"), so those sweep runs were invalid.
No setting reproduces the shipped lengths exactly (about 50 % of common pairs have an identical length), so tuning was stopped: the remaining difference is small.

What the supplement says (Supplementary Section B, "Experiment Details"): NCTC data sets use about 1,000 reads picked so that each read has at least 5 significant DALIGNER overlaps, ground truth by DALIGNER at theta = 0.2 (no DALIGNER parameters given); P. falciparum ground truth by minimap2 2.17 (mapping quality 60; 2,000 reads with degree >= 40); the E. coli data set (SRR11434956, PacBio RS II, 144,238 reads; 113,378 after removing reads with more than 30 % error with fastq-filter) is used only for the runtime and memory comparison, so the runtime benchmark (resolve item 2) needs no ground truth. Paper section 3.2: 100 hash functions, k = 12 for MinHash, kmax = 32 for LexicHash, top T = 5n, 80 CPUs, data sets of 10, 30, 50, 70 and 90 thousand reads (down-sampled) plus the full set; reported full-set cost: LexicHash 100.1 CPU-min and 515.7 MB, MinHash 339.3 CPU-min and 1.95 GB (one CPU).

---

## Batch 6 - branch `fastq-support` (starts from `daligner-ground-truth`; resolve item 2, second part: input needed for the E. coli runtime benchmark)
The E. coli data set of the paper (SRR11434956) is distributed as FASTQ and the C++ program read only FASTA.

| File | Lines | Type | What / why |
|---|---|---|---|
| `mycode/lexic_hash.cpp` | L3 | EDIT | Comment on the zlib include: FASTA or FASTQ. |
| `mycode/lexic_hash.cpp` | L33-L126 (`read_fasta`, comment L33-L42) | EDIT | The reader detects the format from the first non-empty line ('>' = FASTA, '@' = FASTQ). FASTQ records are `@header`, one or more sequence lines, a `+` line, then quality lines; the quality lines are skipped by length (as many characters as the sequence had) because a quality string may start with `@` or `+`. Sequences are cleaned exactly as before (upper case, only A, C, G, T kept), every record keeps its index (an empty read stays in the numbering), the function name and the FASTA path are unchanged. CRLF line endings of the file are kept. |
| `tests/test_fastq_input.py` | L1-L140 | NEW FILE | 8 checks (plain, gzip, CRLF without final newline, quality lines starting with `@` or `+`, N/IUPAC/lower case, wrapped sequence, empty read, `+name` line); every FASTQ run is compared with the FASTA run on the same mask file. All pass; the same test fails on the binary of batch 5. |
| `tests/README.md` | end of file | ADD | Section "FASTQ input (batch 6)". |

Regression checks with the new binary: `tests/test_input_handling.py` (all pass) and `tests/equivalence_py_vs_cpp.py` on NCTC1080 with 100 masks (0 of 200,600 sketch cells differ, 1,003 pairs identical).
Not changed: `--method` is still ignored, reads shorter than max_k are still skipped by the C++ code, there is still no option for mask strategies.
