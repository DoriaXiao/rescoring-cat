# Recovering Accuracy from Miscalibrated Item Banks: Scoring Error Is Fixable, Selection Error Is Not

Code and results for:

> Xiao, X., Domingue, B. W., Ram, N., & Frank, M. C. (in press).
> Recovering accuracy from miscalibrated item banks: Scoring error is fixable,
> selection error is not. *Educational and Psychological Measurement*.

A preprint of the originally submitted version, under the earlier title *Good Kitty,
Bad Bank? Rescoring Miscalibrated CATs Improves Accuracy*, is on PsyArXiv:
<https://doi.org/10.31234/osf.io/gk368_v1>.

**Interactive companion app:** [CAT Rescoring Explorer](https://huggingface.co/spaces/doriax/cat-rescoring-explorer)
(source in [`app/`](app/))

---

## What the paper does

Computerized adaptive tests (CATs) depend on accurately calibrated item banks,
but banks are often built from small pilots, transported across populations, or
seeded with placeholder / AI-generated values. **Rescoring** recomputes ability
estimates from the *already-collected* responses using improved item parameters.
The paper decomposes miscalibration error into **scoring error** (fixable by
rescoring) and **selection error** (not recoverable), and shows:

- **Study 1 (simulation):** scoring error accounts for **86–96%** of
  miscalibration-attributable error across three error families; rescoring drives
  post-hoc bias below 0.5% of MSE.
- **Study 2 (empirical, LEVANTE):** rescoring 2,657 children's scores across four
  tasks shifts and reorders scores, and external validation (age validity,
  cross-task convergent validity) confirms the rescored estimates are *better*,
  not just different.

---

## Repository layout

```
rescoring-cat/
├── study1_simulation/        # Study 1 — fully reproducible, no private data
│   ├── 01_run_simulation.R       # main sim → output/data/res.rds   (~2–3 h)
│   ├── 02_stepwise_highrep.R     # 200-rep stepwise → res_stepwise_v2.rds
│   ├── 03_n_to_rho.R             # N→ρ recovery mapping (Sec 3.2)
│   ├── 04_analysis_figures.R     # Figures 1–4 + simulation tables (fast)
│   ├── 05_figure4_stepwise.R     # re-knit Figure 4 from high-rep data (fast)
│   └── appendix/                 # robustness checks (see appendix/ below)
│
├── study2_empirical/         # Study 2 — requires LEVANTE data (NOT bundled)
│   ├── README.md                 # how to obtain the data + path config
│   ├── 01_rescoring_analysis.R   # Tables 4–7, Fig 5, comparison.csv
│   ├── 02_item_stopping_analysis.R # Appendix D (run after 01, same session)
│   └── 03_cross_task_bootstrap.R # bootstrap CIs on validity changes
│
├── app/                      # Gradio companion app (Hugging Face Space)
├── output/
│   ├── data/                     # simulation outputs (.rds/.csv) — committed
│   │                             #   (comparison.csv is child-level → gitignored)
│   └── tables/                   # result tables (CSV) reproduced in the paper
├── figures/
│   └── paper/                    # publication figures
├── archive/                  # superseded / earlier code, kept for provenance
├── rescoring-cat.Rproj
└── README.md
```

> **Working directory.** All scripts resolve paths with the
> [`here`](https://here.r-lib.org/) package relative to the repo root (marked by
> `.here` / `rescoring-cat.Rproj`). Open `rescoring-cat.Rproj` in RStudio, or run
> scripts with the repo root as the working directory.

---

## Quick start

### Reproduce the figures from bundled results (seconds)

The simulation outputs are committed, so you can regenerate every Study 1 figure
without re-running the multi-hour simulation:

```r
source("study1_simulation/04_analysis_figures.R")  # Figures 1–3 + tables
source("study1_simulation/05_figure4_stepwise.R")  # Figure 4 (high-rep)
```

### Full simulation reproduction (~2–3 h)

```r
source("study1_simulation/01_run_simulation.R")    # → output/data/res.rds
source("study1_simulation/02_stepwise_highrep.R")  # → res_stepwise_v2.rds
source("study1_simulation/04_analysis_figures.R")
source("study1_simulation/05_figure4_stepwise.R")
```

### Study 2 (empirical)

Requires access to LEVANTE child-level data, which is **not** included here. See
[`study2_empirical/README.md`](study2_empirical/README.md).

---

## Figure & table → script map

| Paper item | Description | Produced by | File |
|---|---|---|---|
| **Figure 1** | RMSE vs cross-population alignment ρ | `study1/04` | `figures/paper/fig_sim_rmse_rho.png` |
| **Figure 2** | RMSE under random noise σ | `study1/04` | `figures/paper/fig_sim_rmse_noise.png` |
| **Figure 3** | RMSE under structured shifts Δb | `study1/04` | `figures/paper/fig_sim_rmse_shift.png` |
| **Figure 4** | Stepwise RMSE by test length (ρ=.75, 200 reps) | `study1/05` | `figures/paper/fig_sim_stepwise.png` |
| **Figure 5** | Score shifts by site & task | `study2/01` | `figures/paper/fig_density.png` |
| **Table 1** | Selection/scoring conditions | (conceptual) | — |
| **Table 2** | Bias as % of MSE | `study1/04` | `output/tables/results_sim_bias_decomposition.csv` |
| **Table 3** | Item-parameter alignment (old vs new) | `study2/01` | `output/tables/results_by_site.csv` (+ params) |
| **Table 4** | Parameter effect (1PL EAP, old vs new) | `study2/01` | `output/tables/results_table_a_param_effect.csv` |
| **Table 5** | Score differences by task & site | `study2/01` | `output/tables/results_by_site.csv` |
| **Table 6** | Model effect (1PL vs multigroup IRT) | `study2/01` | `output/tables/results_table_c_model_effect.csv` |
| **Table 7** | Correlation with age | `study2/01` | `output/tables/results_age_validity.csv` |
| **Table 8** | Cross-task convergent validity (bootstrap CIs) | `study2/03` | `output/tables/results_cross_task_bootstrap.csv` |
| **Supp. A** | Bias and correlation by ρ | `study1/04` | `figures/paper/fig_sim_bias_rho.png`, `fig_sim_cor_rho.png` |
| **Supp. B** | Pool-to-pool robustness | `study1/appendix/pool_sensitivity*` | `figures/paper/fig_sim_sensitivity.png` |
| **Supp. C** | Estimated rescoring bank (CAT-R̂) | `study1/appendix/perturbed_rescoring*` | `figures/paper/fig_sim_perturbed.png` |
| **Supp. D** | Item-level shift & stopping behavior | `study2/02` | `figures/paper/fig_item_level_shift.png`, `fig_stopping_by_ability.png` |
| Discussion | N→ρ recovery benchmark (cited in text; figure not in the paper) | `study1/03` | `figures/paper/fig_n_to_rho.png` |

Additional robustness scripts in `study1_simulation/appendix/`:
`widespread_spread*.R` (wider difficulty spread, b ~ N(0,1.5²)),
`diagnose_pools.R` (per-pool diagnostics), `random20_check.R` (random-20 reference).

---

## Key design choices (Study 1)

- **Item pool:** 100 items, 2PL. b ~ N(0, 0.75²); a ~ lognormal(0.1, 0.3), lower
  bound 0.5. One fixed pool reused across all conditions.
- **Examinees:** θ ~ N(0, 1), 100 per cell, 50 replications.
- **Scoring:** EAP with N(0, 1.5) prior on a [−6, 6] grid (step 0.05).
- **Five conditions** (Table 1): CAT-O (operational), CAT-R (rescored),
  CAT-B (benchmark), Fixed-20, Full-100.
- **Three error families:** random noise `b' = b + ε`; structured shift
  `b' = b + Δb + ε`; cross-population misalignment
  `b' = μ_b + ρ(b − μ_b) + ε`, with ε scaled so var(b') = var(b).

---

## Dependencies

```r
install.packages(c(
  "tidyverse", "here", "patchwork", "ggrepel", "withr", "tibble",
  "mirt"   # needed for the N→ρ simulation (study1/03) and item calibration
))
```

R ≥ 4.1 recommended. See [`DEPENDENCIES.R`](DEPENDENCIES.R) for a one-line installer.

---

## Data availability

Study 1 is fully self-contained. Study 2 uses child-level data from the LEVANTE
project (Kachergis et al., 2025), publicly available through the LEVANTE Data
Repository 2026.1 release at <https://stanford.redivis.com/datasets/68kn-csrddrz5x>.
The data are not redistributed here; the empirical scripts and all aggregate result
tables are included so the methods are fully transparent. See [`study2_empirical/README.md`](study2_empirical/README.md).

## Citation

If you use this code, please cite the article (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@article{xiao2026rescoring,
  author  = {Xiao, Xingyao and Domingue, Benjamin W. and Ram, Nil{\`a}m and Frank, Michael C.},
  title   = {Recovering Accuracy from Miscalibrated Item Banks: Scoring Error Is Fixable, Selection Error Is Not},
  journal = {Educational and Psychological Measurement},
  year    = {2026},
  note    = {In press}
}
```

## License

Code released under the [MIT License](LICENSE).
