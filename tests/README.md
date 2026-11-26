# Python vs C++ equivalence test

`equivalence_py_vs_cpp.py` feeds the same masks to the Python reference
(`LexicHash-main/`) and to the C++ port (`mycode/lexic_hash.cpp`) and compares
the sketch matrix and the pairs found by the prefix-tree step.

```bash
g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native -lz
pip install numpy pandas pyfastx tqdm matplotlib scikit-learn
python3 tests/equivalence_py_vs_cpp.py --cpp-bin ./lexichash \
    --fasta LexicHash-main/data/NCTC1080/NCTC1080_reads.fasta.gz \
    --n-hash 100 --alpha 1
```

`--alpha 0` means "all pairs". Other options: `--max-k`, `--min-k`, `--seed`, `--no-rc`.

## New C++ flags

* `--seed N` seeds the mask generator (default: a fresh random seed per run).
* `--mask_file PATH` loads masks from PATH if it exists, otherwise writes the
  generated masks there. One mask per line, `max_k` characters from `0123`
  (0=A, 1=T, 2=G, 3=C, the same encoding as the Python code).

## Results on NCTC1080 (1,003 reads), 50 to 100 masks

| run | sketch cells that differ | pairs only in one program | match lengths that differ |
|---|---|---|---|
| alpha 1, max_k 32, min_k 14 | 0 of 200,600 | 0 | 0 |
| alpha 0 (all 12,761 pairs) | 0 of 100,300 | 0 | 0 |
| alpha 0.1 (cap of 100 pairs) | 0 of 100,300 | 0 | 0 |
| alpha 3, max_k 24, min_k 10 | 0 of 100,300 | 0 | 0 |
| alpha 1, no reverse complements | 0 of 50,150 | 0 | 0 |

## Two bugs this exposed in the code before the fixes

1. Python `n_pairs = n_seq * alpha` was a float. With alpha 0.1 and 1,003 reads
   it is 100.3, the `len(...) == n_pairs` stop test never fires, and the run
   returns all 12,761 pairs instead of 100. Fixed with `int(...)`.
2. C++ `--max_k` and `--min_k` changed the prefix-tree depth but not the
   sketching code, which read separate globals fixed at 32 and 14. With
   `--max_k 24` it hashed 32-mers (sketch values up to about 1.7e19) and then
   read the wrong bits in the tree. The globals now follow the flags, and
   max_k above 32 is rejected.

## Input handling test (batch 2)

`test_input_handling.py` runs the C++ binary on small synthetic FASTA files and checks that
gzip input, CRLF line endings, N / IUPAC / lower-case characters, short reads and the thread
count do not change the result in unexpected ways (details in the file header).

```bash
python3 tests/test_input_handling.py --cpp-bin ./lexichash
```

The C++ program now needs zlib (`-lz`) because it reads `.fasta.gz` directly; no need to
decompress first. Characters other than A, C, G, T are removed and lower case is upper-cased,
as in the Python code. Reads shorter than `max_k` keep their position in the numbering and are
simply not sketched, so output indices always equal the record number in the FASTA file.

## Evaluation tests (batch 3)

* `test_threshold_consistency.py` checks that every evaluation script uses the one ground-truth
  threshold `THETA = 0.2` from `evaluation/lh_eval.py` and that no script passes another value.
* `test_refine_strand.py` checks that `refine_match_pairs` scores a reverse-strand ('-') pair on the
  read / reverse-complement sketch rows (it needs the Python packages numpy and pandas).

```bash
python3 tests/test_threshold_consistency.py
python3 tests/test_refine_strand.py
```

## Mask-strategy and multi-scale tests (batch 4)

* `test_mask_strategies.py` checks the inverse-frequency and GC samplers on toy reads with known answers (the old code is kept in the
  file for comparison), and that every strategy in `evaluation/mask_strategies.py` is reproducible from its seed and writes a valid mask file.
* `test_multiscale_sketch.py` checks the sketch function and the vectorised match lengths of `evaluation/multiscale_eval.py` against brute force.
* `test_merge_stages.py` checks that the vectorised `merge_prediction_files` returns the same rows as the old loop.

```bash
python3 tests/test_mask_strategies.py
python3 tests/test_multiscale_sketch.py
python3 tests/test_merge_stages.py
```
