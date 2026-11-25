# Scorer ablation: P. falciparum (2,000 reads, 200 masks, max_k 32, min_k 14)

Seeds per row: no_rc=2, rc_fwd_rows_only=2, rc_strand_aware=2. Metrics: AUROC and AUPRC over all pairs of reads (unscored pairs = 0), ground truth theta >= 0.2. Mean +- std over mask seeds.

## Setting: no_rc

### Pool: pipeline_a2  (2 seeds; true pairs inside the pool: 21.5%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.065 +- 0.042 | 0.607 +- 0.033 | 0.191 +- 0.064 | 8.4 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.171 +- 0.052 | 0.607 +- 0.033 | 0.214 +- 0.065 | 18.2 |
| exp-decay, best of grid | decay=0.01,top_k=25 | 0.174 +- 0.053 | 0.607 +- 0.033 | 0.213 +- 0.063 | 18.5 |
| mean of top-k y, best | top_k=10 | 0.171 +- 0.051 | 0.607 +- 0.033 | 0.214 +- 0.065 | 18.2 |
| count of masks with y >= t, best | t=14 | 0.157 +- 0.043 | 0.607 +- 0.033 | 0.211 +- 0.061 | 17.2 |
| likelihood-style sum 4^(y-14) | - | 0.083 +- 0.040 | 0.607 +- 0.033 | 0.199 +- 0.059 | 8.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.107 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.109 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.107 (better in 2 of 2 seeds).

### Pool: pipeline_a5  (2 seeds; true pairs inside the pool: 32.5%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.078 +- 0.046 | 0.660 +- 0.033 | 0.191 +- 0.064 | 8.4 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.238 +- 0.049 | 0.661 +- 0.033 | 0.300 +- 0.063 | 22.2 |
| exp-decay, best of grid | decay=0.01,top_k=25 | 0.243 +- 0.050 | 0.661 +- 0.033 | 0.299 +- 0.063 | 22.8 |
| mean of top-k y, best | top_k=10 | 0.239 +- 0.048 | 0.661 +- 0.033 | 0.302 +- 0.062 | 21.9 |
| count of masks with y >= t, best | t=14 | 0.212 +- 0.034 | 0.661 +- 0.033 | 0.289 +- 0.050 | 21.1 |
| likelihood-style sum 4^(y-14) | - | 0.099 +- 0.045 | 0.660 +- 0.033 | 0.199 +- 0.059 | 8.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.159 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.164 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.160 (better in 2 of 2 seeds).

### Pool: pipeline_a10  (2 seeds; true pairs inside the pool: 39.4%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.083 +- 0.045 | 0.693 +- 0.022 | 0.191 +- 0.064 | 8.4 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.270 +- 0.036 | 0.694 +- 0.021 | 0.338 +- 0.044 | 23.2 |
| exp-decay, best of grid | decay=0.02,top_k=25 | 0.277 +- 0.035 | 0.694 +- 0.021 | 0.340 +- 0.041 | 24.4 |
| mean of top-k y, best | top_k=10 | 0.272 +- 0.033 | 0.694 +- 0.021 | 0.342 +- 0.041 | 23.2 |
| count of masks with y >= t, best | t=14 | 0.232 +- 0.016 | 0.694 +- 0.021 | 0.324 +- 0.033 | 22.1 |
| likelihood-style sum 4^(y-14) | - | 0.104 +- 0.043 | 0.693 +- 0.022 | 0.199 +- 0.059 | 8.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.187 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.194 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.189 (better in 2 of 2 seeds).

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool pipeline_a5

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.240 | 0.243 | 0.238 | 0.228 | 0.220 |
| 0.02 | 0.239 | 0.243 | 0.240 | 0.236 | 0.234 |
| 0.05 | 0.237 | 0.242 | 0.241 | 0.241 | 0.241 |
| 0.1 | 0.233 | 0.237 | 0.238 | 0.238 | 0.238 |
| 0.2 | 0.222 | 0.224 | 0.224 | 0.224 | 0.224 |
| 0.3 | 0.207 | 0.208 | 0.208 | 0.208 | 0.208 |
| 0.5 | 0.169 | 0.169 | 0.169 | 0.169 | 0.169 |
| 0.9 | 0.106 | 0.106 | 0.106 | 0.106 | 0.106 |

## Setting: rc_strand_aware

### Pool: pipeline_a2  (2 seeds; true pairs inside the pool: 26.4%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.079 +- 0.045 | 0.632 +- 0.035 | 0.235 +- 0.067 | 9.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.227 +- 0.057 | 0.632 +- 0.035 | 0.264 +- 0.069 | 22.1 |
| exp-decay, best of grid | decay=0.01,top_k=25 | 0.230 +- 0.058 | 0.632 +- 0.035 | 0.264 +- 0.069 | 22.3 |
| mean of top-k y, best | top_k=50 | 0.230 +- 0.059 | 0.632 +- 0.035 | 0.264 +- 0.070 | 22.5 |
| count of masks with y >= t, best | t=14 | 0.223 +- 0.054 | 0.632 +- 0.035 | 0.263 +- 0.069 | 21.9 |
| likelihood-style sum 4^(y-14) | - | 0.145 +- 0.031 | 0.632 +- 0.035 | 0.240 +- 0.072 | 13.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.148 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.151 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.151 (better in 2 of 2 seeds).

### Pool: pipeline_a5  (2 seeds; true pairs inside the pool: 43.5%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.108 +- 0.056 | 0.716 +- 0.044 | 0.235 +- 0.067 | 9.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.360 +- 0.062 | 0.717 +- 0.044 | 0.423 +- 0.081 | 26.6 |
| exp-decay, best of grid | decay=0.01,top_k=25 | 0.367 +- 0.065 | 0.717 +- 0.044 | 0.420 +- 0.082 | 27.9 |
| mean of top-k y, best | top_k=50 | 0.366 +- 0.066 | 0.717 +- 0.044 | 0.415 +- 0.078 | 29.2 |
| count of masks with y >= t, best | t=14 | 0.353 +- 0.058 | 0.717 +- 0.044 | 0.420 +- 0.079 | 27.2 |
| likelihood-style sum 4^(y-14) | - | 0.179 +- 0.043 | 0.716 +- 0.044 | 0.240 +- 0.072 | 13.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.252 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.259 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.258 (better in 2 of 2 seeds).

### Pool: pipeline_a10  (2 seeds; true pairs inside the pool: 58.0%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.124 +- 0.059 | 0.787 +- 0.044 | 0.235 +- 0.067 | 9.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.461 +- 0.056 | 0.788 +- 0.043 | 0.535 +- 0.068 | 28.0 |
| exp-decay, best of grid | decay=0.01,top_k=25 | 0.473 +- 0.059 | 0.788 +- 0.043 | 0.536 +- 0.068 | 29.2 |
| mean of top-k y, best | top_k=50 | 0.470 +- 0.060 | 0.788 +- 0.043 | 0.526 +- 0.068 | 30.5 |
| count of masks with y >= t, best | t=14 | 0.451 +- 0.051 | 0.788 +- 0.043 | 0.532 +- 0.067 | 28.5 |
| likelihood-style sum 4^(y-14) | - | 0.198 +- 0.048 | 0.787 +- 0.044 | 0.240 +- 0.072 | 13.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.337 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.349 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.345 (better in 2 of 2 seeds).

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool pipeline_a5

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.359 | 0.367 | 0.367 | 0.365 | 0.362 |
| 0.02 | 0.359 | 0.367 | 0.367 | 0.366 | 0.366 |
| 0.05 | 0.357 | 0.364 | 0.365 | 0.365 | 0.365 |
| 0.1 | 0.354 | 0.360 | 0.360 | 0.360 | 0.360 |
| 0.2 | 0.345 | 0.348 | 0.348 | 0.348 | 0.348 |
| 0.3 | 0.334 | 0.335 | 0.335 | 0.335 | 0.335 |
| 0.5 | 0.299 | 0.299 | 0.299 | 0.299 | 0.299 |
| 0.9 | 0.191 | 0.191 | 0.191 | 0.191 | 0.191 |

## Setting: rc_fwd_rows_only

### Pool: pipeline_a2  (2 seeds; true pairs inside the pool: 20.3%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.095 +- 0.052 | 0.632 +- 0.035 | 0.250 +- 0.065 | 8.8 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.171 +- 0.056 | 0.632 +- 0.035 | 0.248 +- 0.064 | 17.6 |
| exp-decay, best of grid | decay=0.01,top_k=50 | 0.180 +- 0.056 | 0.632 +- 0.035 | 0.248 +- 0.068 | 19.3 |
| mean of top-k y, best | top_k=50 | 0.179 +- 0.056 | 0.632 +- 0.035 | 0.248 +- 0.068 | 19.5 |
| count of masks with y >= t, best | t=14 | 0.151 +- 0.036 | 0.601 +- 0.023 | 0.204 +- 0.046 | 18.3 |
| likelihood-style sum 4^(y-14) | - | 0.106 +- 0.046 | 0.632 +- 0.035 | 0.239 +- 0.060 | 8.9 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.076 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.085 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.084 (better in 2 of 2 seeds).

### Pool: pipeline_a5  (2 seeds; true pairs inside the pool: 32.6%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.117 +- 0.059 | 0.716 +- 0.044 | 0.269 +- 0.081 | 8.4 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.252 +- 0.064 | 0.716 +- 0.044 | 0.305 +- 0.069 | 21.9 |
| exp-decay, best of grid | decay=0.01,top_k=50 | 0.266 +- 0.064 | 0.716 +- 0.043 | 0.311 +- 0.063 | 24.5 |
| mean of top-k y, best | top_k=50 | 0.265 +- 0.064 | 0.716 +- 0.043 | 0.310 +- 0.062 | 24.6 |
| count of masks with y >= t, best | t=14 | 0.223 +- 0.039 | 0.662 +- 0.026 | 0.305 +- 0.056 | 22.3 |
| likelihood-style sum 4^(y-14) | - | 0.128 +- 0.052 | 0.716 +- 0.043 | 0.270 +- 0.081 | 8.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.135 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.149 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.148 (better in 2 of 2 seeds).

### Pool: pipeline_a10  (2 seeds; true pairs inside the pool: 41.9%)

| scorer | parameters | AUPRC | AUROC | precision@n_true | top-1000 overlap % |
|---|---|---|---|---|---|
| max y (original LexicHash) | - | 0.120 +- 0.057 | 0.787 +- 0.043 | 0.191 +- 0.064 | 8.4 |
| exp-decay as in code (decay 0.1, top_k 100) | decay=0.1,top_k=100 | 0.298 +- 0.057 | 0.787 +- 0.043 | 0.343 +- 0.063 | 23.4 |
| exp-decay, best of grid | decay=0.01,top_k=25 | 0.317 +- 0.055 | 0.787 +- 0.043 | 0.361 +- 0.055 | 24.9 |
| mean of top-k y, best | top_k=50 | 0.313 +- 0.055 | 0.787 +- 0.043 | 0.356 +- 0.055 | 26.1 |
| count of masks with y >= t, best | t=14 | 0.261 +- 0.026 | 0.708 +- 0.021 | 0.349 +- 0.034 | 23.8 |
| likelihood-style sum 4^(y-14) | - | 0.132 +- 0.051 | 0.787 +- 0.043 | 0.199 +- 0.059 | 8.7 |

AUPRC gain of exp-decay as in code (decay 0.1, top_k 100) over max y: +0.179 (better in 2 of 2 seeds).

AUPRC gain of exp-decay, best of grid over max y: +0.198 (better in 2 of 2 seeds).

AUPRC gain of mean of top-k y, best over max y: +0.194 (better in 2 of 2 seeds).

### Exp-decay grid, mean AUPRC (rows: decay, columns: top_k), pool pipeline_a5

| decay | top_k=10 | top_k=25 | top_k=50 | top_k=100 | top_k=all |
|---|---|---|---|---|---|
| 0.01 | 0.255 | 0.265 | 0.266 | 0.261 | 0.255 |
| 0.02 | 0.254 | 0.264 | 0.265 | 0.264 | 0.263 |
| 0.05 | 0.251 | 0.259 | 0.261 | 0.261 | 0.261 |
| 0.1 | 0.245 | 0.251 | 0.252 | 0.252 | 0.252 |
| 0.2 | 0.232 | 0.234 | 0.234 | 0.234 | 0.234 |
| 0.3 | 0.218 | 0.219 | 0.219 | 0.219 | 0.219 |
| 0.5 | 0.189 | 0.190 | 0.190 | 0.190 | 0.190 |
| 0.9 | 0.141 | 0.141 | 0.141 | 0.141 | 0.141 |

