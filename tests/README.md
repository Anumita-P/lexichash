# Python vs C++ equivalence test

`equivalence_py_vs_cpp.py` feeds the same masks to the Python reference
(`LexicHash-main/`) and to the C++ port (`mycode/lexic_hash.cpp`) and compares
the sketch matrix and the pairs found by the prefix-tree step.

```bash
g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native
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
