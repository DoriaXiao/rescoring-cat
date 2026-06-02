# Archive — superseded code (kept for provenance)

These files are **not** part of the reproduction pipeline. They are earlier
versions or exploratory extensions, retained so the project history is auditable.
Use the scripts in `study1_simulation/` and `study2_empirical/` instead.

| File | Status / superseded by |
|---|---|
| `run_full_simulation.R` | Original Study 1 simulation. Superseded by `study1_simulation/01_run_simulation.R` (shared item pool across families + per-rep output enabling SE bands). |
| `run_main_simulation.R` | Old driver that sourced the `R_legacy/` combined-sim engine. Superseded by the standalone v2 scripts. |
| `run_poster_figures.R` | NCME **poster** figures (different styling/scope from the paper). Not used for paper figures. |
| `run_paper_analyses.R` | Combined bootstrap-CI + bias-decomposition script. The bootstrap-CI half is current and lives in `study2_empirical/03_cross_task_bootstrap.R`; the bias-decomposition half assumed the **old** `res.rds` schema (`design`/`select_bank`/`score_with`) and is superseded by `study1_simulation/04_analysis_figures.R`, which reads the current per-rep schema. |
| `archive_levante_rescoring_revised.R` | Earlier empirical analysis. Superseded by `study2_empirical/01_rescoring_analysis.R`. |
| `R_legacy/` | Original combined IRT/CAT engine (`simulation_source.R`, `plotting_source.R`, `plot_traj.R`) used only by the old `run_main_simulation.R` / `run_poster_figures.R`. The current Study 1 scripts are self-contained and do **not** source these. |
| `R_ml/` | Exploratory machine-learning stopping-rule extension (`simulation_ml_source.R`, `comparisons_ml.R`). Not part of the paper. |

> **Schema note.** The current `output/data/res.rds` is per-replication:
> `condition, rep, rmse, bias, cor_r, mse, n, mode, rho, sigma, delta`.
> The archived scripts that reference `design` / `select_bank` / `score_with`
> expect the older long-format output and will not run against the current file.
