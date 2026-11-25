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

<!-- Add each new batch below this line using the same layout. -->
