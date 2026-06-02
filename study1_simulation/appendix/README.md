# Study 1 — appendix / robustness scripts

Each script is self-contained and run from the repo root. They do **not** modify
the main `output/data/res.rds`. Run order: `*_sim.R` / generators first, then the
matching `*_figure.R`.

| Script | Purpose | Paper |
|---|---|---|
| `pool_sensitivity.R` → `pool_sensitivity_figure.R` | CAT-O/CAT-R curves with 10th–90th percentile bands across 10 item pools | Appendix B (pool robustness) → `fig_sim_sensitivity.png` |
| `perturbed_rescoring_sim.R` → `perturbed_rescoring_figure.R` | Rescoring with an *estimated* bank (CAT-R̂), sweeping N_rescore ∈ {100, 300, 1000} | Appendix C → `fig_sim_perturbed.png`, `results_sim_perturbed.csv` |
| `widespread_spread.R` | Cross-population family with wider difficulty spread, b ~ N(0, 1.5²) | robustness (single pool) → `sensitivity_widespread_*.csv` |
| `widespread_spread_multipool.R` | Same, across 10 pools for stable estimates | robustness (multi-pool) → `sensitivity_widespread_multipool.csv` |
| `diagnose_pools.R` | Per-pool structure/error diagnostics for the multi-pool run | diagnostic (no figure) |
| `random20_check.R` | Random-20 vs CAT-B vs Full-100 reference | reviewer check (no figure) |
