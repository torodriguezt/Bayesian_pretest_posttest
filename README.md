# Bayesian analysis for pretest-posttest binary outcomes with adaptive significance levels

Support code for *"Bayesian hypothesis testing in the bivariate Beta-Binomial model with FBST: adaptive cutoff calibration"*.

C++ engine (via Rcpp) for FBST e-values and adaptive cutoffs k\* without MCMC, validated against Stan.

---

## Structure

```text
BivBetaBinom.cpp          # Core engine: posterior, e-value, simulations, k*
BBpost3.stan              # Stan model (validation reference)
Article.tex               # Manuscript
build_article_tables.R    # Main pipeline: generates Tables A and B
build_kstar_table.R       # k* table builder + helpers
replicate_tab2.R          # Reproduces paper's tab2 (posterior-based)
validate_ev_quad.R        # Validates quadrature e-value vs MCMC (Stan)
make_figures.R            # Regenerates all figures
THKS_run.R                # Full THKS experiment pipeline
output/                   # Generated LaTeX tables and .rds files
Figures/                  # Generated figures
original_poster_code/     # Original poster code (with documented bugs)
```

---

## Quick start

Requirements: R ≥ 4.2, packages `Rcpp`, `dplyr`, `tidyr`, `ggplot2`, `ALA`, `rstan` (validation only).

```r
setwd("path/to/repo")

source("validate_ev_quad.R")       # validate C++ engine vs MCMC (~1 min)
source("build_article_tables.R")   # generate Tables A and B (~3 h, M=2000)
source("make_figures.R")           # regenerate all figures (~10 min, M=1000)
```

**Outputs:**

- `output/kstar_prior_table.tex` — Table A: generic k\*(n₁, n₂) lookup
- `output/kstar_posterior_table.tex` — Table B: per-treatment k\* for THKS

---

## Engine API ([BivBetaBinom.cpp](BivBetaBinom.cpp))

| Function | Description |
| --- | --- |
| `bb_constants(...)` | Precomputes log normalizing constant and posterior exponents |
| `densBB_cpp / vec / grid` | Bivariate posterior density |
| `densBB_H_cpp / vec` | Density restricted to H (θ₁=θ₂) |
| `find_sup_H(consts)` | Supremum under H via grid search |
| `ev_quad(consts, sup_H)` | E-value by 2D quadrature |
| `simulate_evs_H / A` | E-value distributions under prior (eqs 10-11) |
| `simulate_evs_H_post / A_post` | E-value distributions under posterior (eqs 16-17) |
| `find_kstar(ev_H, ev_A, a, b)` | k\* = argmin(a·α + b·β) |

---

## References

- Pereira & Stern (1999). Evidence and credibility: full Bayesian significance test. *Entropy*.
- Olkin & Liu (2003). A bivariate beta distribution. *Statistics & Probability Letters*.
- THKS data: R package `ALA` (`tvsfp`).
