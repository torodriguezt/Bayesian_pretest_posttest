# BivBetaBinomial — Adaptive cutoff k\* for the FBST in the bivariate Beta-Binomial model

Código de soporte del artículo *"Bayesian hypothesis testing in the bivariate Beta-Binomial model with FBST: adaptive cutoff calibration"*.

Implementa un motor en C++ (vía Rcpp) para calcular el e-valor FBST y los cutoffs adaptativos k\* sin recurrir a MCMC, validado contra una implementación de referencia con Stan.

---

## Estructura

```
.
├── BivBetaBinom.cpp              # Motor C++: posterior cerrada, ev, simulaciones, k*
├── BBpost3.stan                  # Modelo Stan (referencia para validación)
├── Article.tex                   # Manuscrito
│
├── build_article_tables.R        # Pipeline principal: genera Tablas A y B
├── build_kstar_table.R           # Helpers + demos para la Tabla A
├── replicate_tab2.R              # Reproduce los k* de tab2 del paper (posterior-based)
├── validate_ev_quad.R            # Valida ev por cuadratura vs MCMC (Stan)
├── THKS_run.R                    # Pipeline completo del experimento THKS
│
├── Figures/                      # Figuras generadas
├── output/                       # Tablas LaTeX y .rds resultantes
└── original_poster_code/         # Código del poster original (con bugs documentados)
```

---

## Cómo correr

Pre-requisitos: R ≥ 4.2, paquetes `Rcpp`, `dplyr`, `tidyr`, `ggplot2`, `ALA` (datos THKS), `rstan` (sólo para validación).

```r
setwd("ruta/al/repo")

# 1) Validar que el motor C++ coincide con MCMC (~1 min)
source("validate_ev_quad.R")

# 2) Generar las tablas del artículo (~3 hs con M = 2000)
source("build_article_tables.R")
#   → output/kstar_prior_table.tex     (Tabla A: genérica k*(n1, n2))
#   → output/kstar_posterior_table.tex (Tabla B: per-tratamiento THKS)

# 3) (Opcional) Reproducir tab2 del paper bajo posterior-based
source("replicate_tab2.R")
```

---

## Hallazgos metodológicos

### 1. Bugs en el código del poster original (`original_poster_code/`)

Al revisar `alpha_final_THKS.R` y `beta_final_THKS.R` se detectaron **dos errores fundamentales**:

- **Estimador MC mal construido**: el indicador `I(ev ≤ k)` se multiplica por la verosimilitud `f(x|θ)` antes de promediar. Esto computa `E_X[f(X|θ)·𝟙]` en lugar de `E_X[𝟙]`, sesgando α y β.
- **θ ~ Uniform(0,1) en lugar de θ ~ f_H** para el cálculo de α. Ignora el peso de la "integral de línea" derivada en el apéndice del artículo.

Ambos están corregidos en este motor.

### 2. La prior KL-óptima restringida a H es impropia

Con (α₀, α₁, α₂) = (0.8374, 0.8411, 0.8053):

> f_H(t) ∝ t^{α₁+α₂−2} · (1−t)^{α₀−2} · (1+t)^{−α}

Como α₀ < 1, el exponente de (1−t) es **−1.16** y `∫₀¹ (1−t)^{−1.16} dt = ∞`. La densidad de la prior restringida a la línea θ₁=θ₂ **diverge en t→1**.

Consecuencia práctica visible en [Figures/error_curves_n25.png](Figures/error_curves_n25.png): para n=25, el α salta abruptamente cerca de 0 y el mínimo de α+β queda en ≈0.97. La formulación posterior-based "se salva" porque la verosimilitud Bin(x|n,θ)→0 cerca de t=1 cuando x<n, dominando la singularidad de la prior.

### 3. ev validado contra MCMC

`ev_quad` (cuadratura de Simpson 2D sobre la posterior cerrada) reproduce `ev_FBST` (10 000 muestras MCMC con Stan) hasta ~3 decimales. Speedup: **~900×** (8 ms vs 7 s).

### 4. Tablas resultantes

- **Tabla A — `output/kstar_prior_table.tex`**: lookup genérico k\*(n₁, n₂) en grilla {10, 20, 30, 40, 50, 75, 100, 150, 200}, M=2000.
- **Tabla B — `output/kstar_posterior_table.tex`**: k\* per-tratamiento del THKS, posterior-based. Reemplaza tab2 del manuscrito.

Para los 4 tratamientos del THKS, la decisión final coincide con el paper original (rechazar H en todos), pero los k\* difieren — los nuestros son los corregidos.

---

## API del motor

Funciones principales en [BivBetaBinom.cpp](BivBetaBinom.cpp):

| Función | Descripción |
|---|---|
| `bb_constants(n1, n2, x1, x2, a0, a1, a2)` | Precomputa `log_C` (incluye `−log ₃F₂(1)`) y exponentes de la posterior |
| `densBB_cpp / vec / grid` | Densidad bivariada en θ |
| `densBB_H_cpp / vec` | Densidad restringida a H |
| `find_sup_H(consts)` | Supremo bajo H (grid search determinístico) |
| `ev_quad(consts, sup_H)` | e-valor por cuadratura 2D |
| `ev_quad_from_data(...)` | Wrapper conveniente |
| `simulate_evs_H / A` | Distribuciones de ev bajo prior (formulación 10-11) |
| `simulate_evs_H_post / A_post` | Distribuciones de ev bajo posterior (formulación 16-17) |
| `find_kstar(ev_H, ev_A, a, b)` | k\* = argmin(a·α + b·β) por enumeración exacta |

---

## Referencias

- Pereira, C. A. B., & Stern, J. M. (1999). Evidence and credibility: full Bayesian significance test for precise hypotheses. *Entropy*.
- Olkin, I., & Liu, R. (2003). A bivariate beta distribution. *Statistics & Probability Letters*.
- Datos THKS: paquete R `ALA` (`tvsfp`).
