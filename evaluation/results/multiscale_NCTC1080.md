Seeds: 3. NCTC1080, no reverse complements, all pairs scored, ground truth theta >= 0.2.

| method | masks | AUPRC | AUROC | precision@n_true | top-1000 overlap (%) |
|---|---|---|---|---|---|
| multiscale as coded (max of reported lengths) | 2300 | 0.858 +- 0.004 | 0.9972 | 0.829 | 28.3 |
| stage 1 alone (k 63-128, 700 masks, alpha 10) | 700 | 0.010 +- 0.001 | 0.5013 | 0.012 | 0.2 |
| stage 2 alone (k 31-64, 600 masks, alpha 100) | 600 | 0.395 +- 0.010 | 0.7142 | 0.432 | 28.3 |
| stage 3 alone (k 14-32, 600 masks, alpha 1000) | 600 | 0.851 +- 0.013 | 0.9956 | 0.829 | 23.7 |
| stage 4 alone (k 13-16, 400 masks, alpha 10) | 400 | 0.642 +- 0.004 | 0.9882 | 0.656 | 11.8 |
| single scale k=32, 500 masks, max y | 500 | 0.849 +- 0.004 | 0.9981 | 0.825 | 25.4 |
| single scale k=32, 500 masks, exp-decay(0.1,100) | 500 | 0.949 +- 0.001 | 0.9988 | 0.887 | 36.8 |
| single scale k=32, 2300 masks, max y | 2300 | 0.836 +- 0.003 | 0.9987 | 0.818 | 16.6 |
| single scale k=32, 2300 masks, exp-decay(0.1,100) | 2300 | 0.941 +- 0.000 | 0.9994 | 0.874 | 35.2 |
| multiscale pooled, y/max_k, exp-decay(0.1,100) | 2300 | 0.960 +- 0.002 | 0.9983 | 0.896 | 40.4 |
