# Study 2 — Empirical Application (LEVANTE)

These scripts reproduce the empirical results in the paper (Tables 4–7, Figure 5,
Appendix D, and the cross-task bootstrap CIs). They operate on **child-level data
from the LEVANTE project**, which is **not redistributed in this repository**.

## Required data (not included)

Set two paths (via environment variables or by editing the `CONFIG` block at the
top of `01_rescoring_analysis.R`) to point at your local LEVANTE copies:

| Variable | Must contain |
|---|---|
| `LEVANTE_PILOTS` | `01_fetched_data/task_data_nested.rds`, `01_fetched_data/run_data.rds` |
| `LEVANTE_RUNS`   | `all_scores.csv`, `item parameter comparison.csv` |

```r
Sys.setenv(LEVANTE_PILOTS = "/path/to/levante-pilots")
Sys.setenv(LEVANTE_RUNS   = "/path/to/levante-runs")
```

LEVANTE data access: the LEVANTE pilot data analyzed in the paper are publicly
available through the LEVANTE Data Repository 2026.1 release
(https://stanford.redivis.com/datasets/68kn-csrddrz5x). See also Kachergis et al.
(2025) and the LEVANTE project (https://github.com/levante-framework). The raw item-level trial data and
per-child ability estimates are governed by the project's human-subjects data
agreements.

## Run order

```r
# 1. Main empirical pipeline — writes tables, figures, and (locally) comparison.csv
source("study2_empirical/01_rescoring_analysis.R")

# 2. Appendix D — REUSES objects from script 01; run in the SAME R session
source("study2_empirical/02_item_stopping_analysis.R")

# 3. Bootstrap CIs on cross-task validity changes (reads comparison.csv)
source("study2_empirical/03_cross_task_bootstrap.R")
```

## Outputs

- `output/data/comparison.csv` — per-child, per-task ability estimates
  (old vs new parameters, age, site). **Child-level → git-ignored**; generated
  locally when you run script 01.
- `output/tables/results_*.csv` — aggregate result tables reproduced in the paper
  (these are committed: they contain only group-level statistics).
- `figures/paper/fig_scatter.png`, `fig_density.png` — score-shift figures (Fig 5).
- `figures/paper/fig_item_level_shift.png`, `fig_stopping_by_ability.png` — App. D.

## Scoring method

To isolate the *parameter* effect from the *scoring model*, both old and new
parameters are scored with the **same** 1PL EAP estimator (prior N(0, 1.5); grid
[−6, 6], step 0.05). Multigroup IRT scores are also computed to quantify the
residual model effect (Table 6). Children with fewer than 5 responses per task are
excluded; the latest run is retained when multiple administrations exist.
