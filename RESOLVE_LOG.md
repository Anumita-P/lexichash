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

<!-- Add each new batch below this line using the same layout. -->
