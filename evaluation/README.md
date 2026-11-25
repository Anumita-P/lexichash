# Evaluation

Shared evaluation code and the scorer ablation (resolve-list items 3 and 4).

| File | What it is |
|---|---|
| `lh_eval.py` | The one place that defines the ground truth and the metrics. `THETA = 0.2`: a pair is a true overlap when `a / (l1 + l2 - a) >= 0.2` (alignment length `a`, read lengths `l1`, `l2`), as in the LexicHash paper. ROC-AUC and PR-AUC are computed over **all** read pairs; pairs a method did not return count as score 0. |
| `scorer_ablation.py` | Compares ways of turning the per-mask match lengths of a read pair into one score. |
| `summarize_ablation.py` | Turns the raw CSV into a Markdown report. |
| `results/` | `scorer_ablation_<dataset>.csv` (raw, one row per seed, pool, setting and scorer) and `scorer_ablation_<dataset>.md` (the report). |

`tests/test_threshold_consistency.py` checks that every evaluation script uses `THETA` and that no
script passes a different threshold.

## Running the ablation

```bash
g++ -std=c++17 -fopenmp mycode/lexic_hash.cpp -o lexichash -O3 -march=native -lz
python3 evaluation/scorer_ablation.py --cpp-bin ./lexichash \
    --fasta LexicHash-main/data/NCTC1080/NCTC1080_reads.fasta.gz \
    --gt LexicHash-main/data/NCTC1080/NCTC1080_daligner_ground_truth.txt \
    --seeds 1 2 3 4 5 --rc-seeds 1 --out evaluation/results/scorer_ablation_NCTC1080
python3 evaluation/summarize_ablation.py --csv evaluation/results/scorer_ablation_NCTC1080.csv \
    --title "NCTC1080" --out evaluation/results/scorer_ablation_NCTC1080.md
```

The reverse-complement settings are about 6 times slower on NCTC1080, so `--rc-seeds` limits them to
some seeds. One seed of the NCTC1080 run takes about 2 minutes without reverse complements.

## What is compared

For a pair of reads, each mask `m` gives a shared-prefix length `y_m` of the two min-hash k-mers.

| Scorer | Definition |
|---|---|
| `maxy` | `max_m y_m` (original LexicHash) |
| `expdecay` | sort `y` descending, `sum_i y_(i) * (1 - decay)^i` over the `top_k` largest (our score) |
| `meantop` | mean of the `top_k` largest `y` |
| `count_ge` | number of masks with `y_m >= t` |
| `lik` | `sum_m 4^max(y_m - 14, 0)` |

**Pools.** `all_pairs` scores every pair (an upper bound for a scorer). `pipeline_aN` scores only the
`N * n_reads` pairs the prefix-tree step returns, which is how the scorer is used in the real pipeline.

**Reverse-complement settings.** `no_rc`: run without reverse complements. `rc_strand_aware`: reverse
complements on and the sketch rows chosen to match the strand of the pair. `rc_fwd_rows_only`: reverse
complements on, but `y` taken from the forward rows only, which is what `refine_match_pairs()` in
`experiment with lexichash/` does.

## Things to keep in mind

* "Best of grid" is picked by AUPRC on the same data, so it is optimistic.
* ROC/PRC use the unordered pair and ignore the strand.
* The scorer only reorders the pairs the prefix-tree step returned. How many true pairs that step
  returns (`alpha`) limits what any scorer can reach; the report lists the share of true pairs inside
  each pool.
* Match lengths are read from the sketches produced by the C++ program, which is identical to the Python
  code on the same masks (`tests/equivalence_py_vs_cpp.py`). The vectorised `y` and exp-decay scores were
  checked against `yscore()` and `expdecayscore()` in the Python code (0 mismatches in 6,000 values).
