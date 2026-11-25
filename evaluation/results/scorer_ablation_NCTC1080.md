# Scorer ablation: NCTC1080 (1,003 reads, 500 masks, max_k 32, min_k 14)

Seeds per row: no_rc=5, rc_fwd_rows_only=1, rc_strand_aware=1. Metrics: AUROC and AUPRC over all pairs of reads (unscored pairs = 0), ground truth theta >= 0.2. Mean +- std over mask seeds.

## Setting: no_rc

### Pool: pipeline_a1  (5 seeds; true pairs inside the pool: 23.7%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.221 +- 0.002 | 0.618 +- 0.001 | 0.241 +- 0.001 | 29.7 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.240 +- 0.002 | 0.618 +- 0.001 | 0.241 +- 0.001 | 29.8 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.241 +- 0.001 | 0.618 +- 0.001 | 0.241 +- 0.001 | 29.8 |
| mean of top-k y, best | top_k=100 | 0.241 +- 0.001 | 0.618 +- 0.001 | 0.241 +- 0.001 | 29.8 |
| count of masks with y >= t, best | t=14 | 0.241 +- 0.001 | 0.618 +- 0.001 | 0.241 +- 0.001 | 29.8 |
| likelihood-style sum 4^(y-14) | - | 0.233 +- 0.003 | 0.618 +- 0.001 | 0.241 +- 0.001 | 29.8 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.019 (better in 5 of 5 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.020 (better in 5 of 5 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.020 (better in 5 of 5 seeds).

### Pool: pipeline_a2  (5 seeds; true pairs inside the pool: 46.0%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.416 +- 0.006 | 0.730 +- 0.002 | 0.462 +- 0.003 | 25.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.457 +- 0.004 | 0.730 +- 0.002 | 0.462 +- 0.003 | 36.1 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.461 +- 0.003 | 0.730 +- 0.002 | 0.462 +- 0.003 | 40.5 |
| mean of top-k y, best | top_k=100 | 0.461 +- 0.003 | 0.730 +- 0.002 | 0.462 +- 0.003 | 40.2 |
| count of masks with y >= t, best | t=14 | 0.460 +- 0.003 | 0.730 +- 0.002 | 0.462 +- 0.003 | 38.4 |
| likelihood-style sum 4^(y-14) | - | 0.435 +- 0.005 | 0.730 +- 0.002 | 0.462 +- 0.003 | 31.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.041 (better in 5 of 5 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.045 (better in 5 of 5 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.045 (better in 5 of 5 seeds).

### Pool: pipeline_a5  (5 seeds; true pairs inside the pool: 95.3%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.825 +- 0.007 | 0.976 +- 0.002 | 0.825 +- 0.008 | 25.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.915 +- 0.005 | 0.976 +- 0.002 | 0.885 +- 0.003 | 37.2 |
| exp-decay, best of grid | decay=0.01,top_k=100 | 0.930 +- 0.004 | 0.976 +- 0.002 | 0.899 +- 0.002 | 41.3 |
| mean of top-k y, best | top_k=100 | 0.932 +- 0.004 | 0.976 +- 0.002 | 0.900 +- 0.003 | 42.3 |
| count of masks with y >= t, best | t=14 | 0.922 +- 0.004 | 0.976 +- 0.002 | 0.889 +- 0.004 | 40.2 |
| likelihood-style sum 4^(y-14) | - | 0.849 +- 0.007 | 0.976 +- 0.002 | 0.828 +- 0.007 | 31.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.090 (better in 5 of 5 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.105 (better in 5 of 5 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.106 (better in 5 of 5 seeds).

### Pool: pipeline_a10  (5 seeds; true pairs inside the pool: 98.6%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.844 +- 0.006 | 0.992 +- 0.001 | 0.825 +- 0.008 | 25.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.942 +- 0.003 | 0.993 +- 0.001 | 0.885 +- 0.004 | 37.2 |
| exp-decay, best of grid | decay=0.01,top_k=100 | 0.957 +- 0.002 | 0.993 +- 0.001 | 0.903 +- 0.002 | 41.3 |
| mean of top-k y, best | top_k=100 | 0.958 +- 0.002 | 0.993 +- 0.001 | 0.905 +- 0.003 | 42.3 |
| count of masks with y >= t, best | t=14 | 0.943 +- 0.002 | 0.993 +- 0.001 | 0.890 +- 0.003 | 40.2 |
| likelihood-style sum 4^(y-14) | - | 0.870 +- 0.005 | 0.992 +- 0.001 | 0.828 +- 0.007 | 31.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.098 (better in 5 of 5 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.113 (better in 5 of 5 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.114 (better in 5 of 5 seeds).

### Pool: all_pairs  (5 seeds; true pairs inside the pool: 99.1%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.846 +- 0.006 | 0.997 +- 0.001 | 0.825 +- 0.008 | 25.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.948 +- 0.002 | 0.999 +- 0.000 | 0.885 +- 0.004 | 37.2 |
| exp-decay, best of grid | decay=0.01,top_k=50 | 0.960 +- 0.001 | 0.997 +- 0.000 | 0.901 +- 0.002 | 40.2 |
| mean of top-k y, best | top_k=50 | 0.961 +- 0.001 | 0.997 +- 0.000 | 0.903 +- 0.002 | 40.6 |
| count of masks with y >= t, best | t=14 | 0.943 +- 0.002 | 0.995 +- 0.000 | 0.890 +- 0.003 | 40.2 |
| likelihood-style sum 4^(y-14) | - | 0.869 +- 0.005 | 0.992 +- 0.001 | 0.828 +- 0.007 | 31.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.102 (better in 5 of 5 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.115 (better in 5 of 5 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.115 (better in 5 of 5 seeds).

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool pipeline_a5

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.911 | 0.925 | 0.929 | 0.930 | 0.930 |
| 0.02 | 0.910 | 0.924 | 0.927 | 0.929 | 0.929 |
| 0.05 | 0.909 | 0.920 | 0.923 | 0.923 | 0.923 |
| 0.1 | 0.906 | 0.914 | 0.915 | 0.915 | 0.915 |
| 0.2 | 0.899 | 0.903 | 0.903 | 0.903 | 0.903 |
| 0.3 | 0.892 | 0.893 | 0.893 | 0.893 | 0.893 |
| 0.5 | 0.877 | 0.877 | 0.877 | 0.877 | 0.877 |
| 0.9 | 0.852 | 0.852 | 0.852 | 0.852 | 0.852 |

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool all_pairs

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.943 | 0.957 | 0.960 | 0.957 | 0.949 |
| 0.02 | 0.942 | 0.956 | 0.959 | 0.959 | 0.958 |
| 0.05 | 0.941 | 0.953 | 0.955 | 0.956 | 0.956 |
| 0.1 | 0.938 | 0.947 | 0.948 | 0.948 | 0.948 |
| 0.2 | 0.930 | 0.934 | 0.935 | 0.935 | 0.935 |
| 0.3 | 0.922 | 0.923 | 0.923 | 0.923 | 0.923 |
| 0.5 | 0.905 | 0.905 | 0.905 | 0.905 | 0.905 |
| 0.9 | 0.877 | 0.877 | 0.877 | 0.877 | 0.877 |

## Setting: rc_strand_aware

### Pool: pipeline_a1  (1 seed; true pairs inside the pool: 15.6%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.100 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.118 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.123 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| mean of top-k y, best | top_k=all | 0.125 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| count of masks with y >= t, best | t=14 | 0.123 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| likelihood-style sum 4^(y-14) | - | 0.114 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.019 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.023 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.025 (better in 1 of 1 seeds).

### Pool: pipeline_a2  (1 seed; true pairs inside the pool: 29.9%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.178 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 13.5 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.212 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 24.5 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.222 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 27.1 |
| mean of top-k y, best | top_k=all | 0.225 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 29.3 |
| count of masks with y >= t, best | t=14 | 0.222 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 27.1 |
| likelihood-style sum 4^(y-14) | - | 0.202 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 21.8 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.034 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.043 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.046 (better in 1 of 1 seeds).

### Pool: pipeline_a5  (1 seed; true pairs inside the pool: 68.5%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.383 +- 0.000 | 0.840 +- 0.000 | 0.556 +- 0.000 | 10.3 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.463 +- 0.000 | 0.841 +- 0.000 | 0.610 +- 0.000 | 28.0 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.486 +- 0.000 | 0.841 +- 0.000 | 0.623 +- 0.000 | 32.5 |
| mean of top-k y, best | top_k=all | 0.491 +- 0.000 | 0.841 +- 0.000 | 0.625 +- 0.000 | 37.9 |
| count of masks with y >= t, best | t=14 | 0.486 +- 0.000 | 0.841 +- 0.000 | 0.623 +- 0.000 | 33.4 |
| likelihood-style sum 4^(y-14) | - | 0.424 +- 0.000 | 0.840 +- 0.000 | 0.557 +- 0.000 | 22.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.080 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.104 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.109 (better in 1 of 1 seeds).

### Pool: pipeline_a10  (1 seed; true pairs inside the pool: 98.7%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.531 +- 0.000 | 0.990 +- 0.000 | 0.556 +- 0.000 | 10.3 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.635 +- 0.000 | 0.991 +- 0.000 | 0.624 +- 0.000 | 28.0 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.670 +- 0.000 | 0.991 +- 0.000 | 0.650 +- 0.000 | 32.8 |
| mean of top-k y, best | top_k=all | 0.674 +- 0.000 | 0.991 +- 0.000 | 0.655 +- 0.000 | 38.2 |
| count of masks with y >= t, best | t=14 | 0.668 +- 0.000 | 0.991 +- 0.000 | 0.646 +- 0.000 | 33.1 |
| likelihood-style sum 4^(y-14) | - | 0.575 +- 0.000 | 0.990 +- 0.000 | 0.557 +- 0.000 | 22.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.104 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.139 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.143 (better in 1 of 1 seeds).

### Pool: all_pairs  (1 seed; true pairs inside the pool: 99.7%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.533 +- 0.000 | 0.996 +- 0.000 | 0.556 +- 0.000 | 10.3 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.640 +- 0.000 | 0.997 +- 0.000 | 0.624 +- 0.000 | 28.0 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.674 +- 0.000 | 0.995 +- 0.000 | 0.650 +- 0.000 | 32.8 |
| mean of top-k y, best | top_k=100 | 0.675 +- 0.000 | 0.995 +- 0.000 | 0.650 +- 0.000 | 32.3 |
| count of masks with y >= t, best | t=14 | 0.671 +- 0.000 | 0.996 +- 0.000 | 0.646 +- 0.000 | 33.1 |
| likelihood-style sum 4^(y-14) | - | 0.578 +- 0.000 | 0.995 +- 0.000 | 0.557 +- 0.000 | 22.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.107 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.140 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.141 (better in 1 of 1 seeds).

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool pipeline_a5

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.453 | 0.468 | 0.477 | 0.483 | 0.486 |
| 0.02 | 0.453 | 0.467 | 0.475 | 0.480 | 0.481 |
| 0.05 | 0.452 | 0.465 | 0.470 | 0.472 | 0.472 |
| 0.1 | 0.451 | 0.461 | 0.463 | 0.463 | 0.463 |
| 0.2 | 0.448 | 0.452 | 0.453 | 0.453 | 0.453 |
| 0.3 | 0.444 | 0.446 | 0.446 | 0.446 | 0.446 |
| 0.5 | 0.438 | 0.438 | 0.438 | 0.438 | 0.438 |
| 0.9 | 0.426 | 0.426 | 0.426 | 0.426 | 0.426 |

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool all_pairs

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.627 | 0.650 | 0.662 | 0.670 | 0.674 |
| 0.02 | 0.627 | 0.649 | 0.659 | 0.665 | 0.667 |
| 0.05 | 0.625 | 0.644 | 0.652 | 0.653 | 0.653 |
| 0.1 | 0.623 | 0.637 | 0.640 | 0.640 | 0.640 |
| 0.2 | 0.617 | 0.623 | 0.623 | 0.623 | 0.623 |
| 0.3 | 0.610 | 0.613 | 0.613 | 0.613 | 0.613 |
| 0.5 | 0.598 | 0.599 | 0.599 | 0.599 | 0.599 |
| 0.9 | 0.581 | 0.581 | 0.581 | 0.581 | 0.581 |

## Setting: rc_fwd_rows_only

### Pool: pipeline_a1  (1 seed; true pairs inside the pool: 15.6%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.154 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.160 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.161 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| mean of top-k y, best | top_k=100 | 0.161 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |
| count of masks with y >= t, best | t=14 | 0.161 +- 0.000 | 0.578 +- 0.000 | 0.162 +- 0.000 | 18.6 |
| likelihood-style sum 4^(y-14) | - | 0.156 +- 0.000 | 0.577 +- 0.000 | 0.161 +- 0.000 | 18.6 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.007 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.008 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.008 (better in 1 of 1 seeds).

### Pool: pipeline_a2  (1 seed; true pairs inside the pool: 29.9%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.283 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 30.2 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.301 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 33.0 |
| exp-decay, best of grid | decay=0.01,top_k=all | 0.303 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 33.7 |
| mean of top-k y, best | top_k=100 | 0.303 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 33.8 |
| count of masks with y >= t, best | t=14 | 0.303 +- 0.000 | 0.650 +- 0.000 | 0.303 +- 0.000 | 33.1 |
| likelihood-style sum 4^(y-14) | - | 0.292 +- 0.000 | 0.649 +- 0.000 | 0.302 +- 0.000 | 30.4 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.019 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.021 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.021 (better in 1 of 1 seeds).

### Pool: pipeline_a5  (1 seed; true pairs inside the pool: 68.5%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.609 +- 0.000 | 0.842 +- 0.000 | 0.685 +- 0.000 | 26.6 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.670 +- 0.000 | 0.842 +- 0.000 | 0.685 +- 0.000 | 37.9 |
| exp-decay, best of grid | decay=0.01,top_k=100 | 0.679 +- 0.000 | 0.842 +- 0.000 | 0.685 +- 0.000 | 41.6 |
| mean of top-k y, best | top_k=100 | 0.680 +- 0.000 | 0.842 +- 0.000 | 0.685 +- 0.000 | 42.8 |
| count of masks with y >= t, best | t=14 | 0.676 +- 0.000 | 0.842 +- 0.000 | 0.686 +- 0.000 | 40.5 |
| likelihood-style sum 4^(y-14) | - | 0.629 +- 0.000 | 0.842 +- 0.000 | 0.684 +- 0.000 | 32.0 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.062 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.070 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.071 (better in 1 of 1 seeds).

### Pool: pipeline_a10  (1 seed; true pairs inside the pool: 98.4%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.847 +- 0.000 | 0.993 +- 0.000 | 0.827 +- 0.000 | 26.6 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.943 +- 0.000 | 0.993 +- 0.000 | 0.886 +- 0.000 | 38.0 |
| exp-decay, best of grid | decay=0.01,top_k=100 | 0.959 +- 0.000 | 0.993 +- 0.000 | 0.904 +- 0.000 | 41.7 |
| mean of top-k y, best | top_k=100 | 0.961 +- 0.000 | 0.993 +- 0.000 | 0.907 +- 0.000 | 43.2 |
| count of masks with y >= t, best | t=14 | 0.946 +- 0.000 | 0.992 +- 0.000 | 0.892 +- 0.000 | 40.9 |
| likelihood-style sum 4^(y-14) | - | 0.870 +- 0.000 | 0.993 +- 0.000 | 0.829 +- 0.000 | 32.0 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.096 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.112 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.113 (better in 1 of 1 seeds).

### Pool: all_pairs  (1 seed; true pairs inside the pool: 99.2%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.846 +- 0.000 | 0.998 +- 0.000 | 0.827 +- 0.000 | 26.6 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.949 +- 0.000 | 0.999 +- 0.000 | 0.886 +- 0.000 | 38.0 |
| exp-decay, best of grid | decay=0.01,top_k=50 | 0.961 +- 0.000 | 0.997 +- 0.000 | 0.903 +- 0.000 | 40.1 |
| mean of top-k y, best | top_k=50 | 0.962 +- 0.000 | 0.997 +- 0.000 | 0.904 +- 0.000 | 41.1 |
| count of masks with y >= t, best | t=14 | 0.946 +- 0.000 | 0.996 +- 0.000 | 0.891 +- 0.000 | 40.9 |
| likelihood-style sum 4^(y-14) | - | 0.870 +- 0.000 | 0.993 +- 0.000 | 0.829 +- 0.000 | 32.0 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.103 (better in 1 of 1 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.116 (better in 1 of 1 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.116 (better in 1 of 1 seeds).

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool pipeline_a5

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.667 | 0.675 | 0.678 | 0.679 | 0.679 |
| 0.02 | 0.667 | 0.675 | 0.677 | 0.678 | 0.678 |
| 0.05 | 0.666 | 0.673 | 0.675 | 0.675 | 0.675 |
| 0.1 | 0.664 | 0.670 | 0.670 | 0.670 | 0.670 |
| 0.2 | 0.660 | 0.663 | 0.663 | 0.663 | 0.663 |
| 0.3 | 0.656 | 0.657 | 0.657 | 0.657 | 0.657 |
| 0.5 | 0.647 | 0.647 | 0.647 | 0.647 | 0.647 |
| 0.9 | 0.631 | 0.631 | 0.631 | 0.631 | 0.631 |

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool all_pairs

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.944 | 0.958 | 0.961 | 0.960 | 0.952 |
| 0.02 | 0.944 | 0.957 | 0.960 | 0.960 | 0.960 |
| 0.05 | 0.942 | 0.954 | 0.956 | 0.957 | 0.957 |
| 0.1 | 0.939 | 0.948 | 0.949 | 0.949 | 0.949 |
| 0.2 | 0.931 | 0.935 | 0.936 | 0.936 | 0.936 |
| 0.3 | 0.923 | 0.924 | 0.924 | 0.924 | 0.924 |
| 0.5 | 0.906 | 0.906 | 0.906 | 0.906 | 0.906 |
| 0.9 | 0.877 | 0.877 | 0.877 | 0.877 | 0.877 |

