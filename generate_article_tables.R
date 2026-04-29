# generate_article_tables.R
# Genera las tablas corregidas para el artículo (sin McNemar).
#
#   output/tab2_corrected.tex  — prior KL-óptima (α≈0.84), posterior-based
#   output/tab3_corrected.tex  — prior informativa (α=10),  posterior-based
#
# Reemplaza directamente tab2/tab3 del manuscrito.

library(ALA)
library(dplyr)
library(Rcpp)

# setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")
dir.create("output", showWarnings = FALSE)

source("priors_config.R")
a0_kl   <- prior_NI["a0"];   a1_kl   <- prior_NI["a1"];   a2_kl   <- prior_NI["a2"]
a0_inf  <- prior_INF["a0"];  a1_inf  <- prior_INF["a1"];  a2_inf  <- prior_INF["a2"]
a0_conf <- prior_CONF["a0"]; a1_conf <- prior_CONF["a1"]; a2_conf <- prior_CONF["a2"]
M      <- 2000

# ---------------------------------------------------------------------------
# Datos THKS
# ---------------------------------------------------------------------------
datos1 <- tvsfp

prep_school <- function(d, sid, sb, tv) {
  x <- d %>% mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
        filter(school == sid, school.based == sb, tv.based == tv)
  list(X = (x %>% group_by(stage) %>% summarise(B = sum(binTHKS)))$B,
       n = (x %>% group_by(stage) %>% summarise(n = n()))$n)
}
prep_global <- function(d, sb, tv) {
  x <- d %>% mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
        filter(school.based == sb, tv.based == tv)
  list(X = (x %>% group_by(stage) %>% summarise(B = sum(binTHKS)))$B,
       n = (x %>% group_by(stage) %>% summarise(n = n()))$n)
}

grupos <- list(
  yy = list(g = prep_school(datos1, "404", "yes", "yes"), lab = "CC + TV"),
  yn = list(g = prep_school(datos1, "408", "yes", "no"),  lab = "CC, no TV"),
  ny = list(g = prep_global(datos1,         "no",  "yes"), lab = "no CC, TV"),
  nn = list(g = prep_school(datos1, "409", "no",  "no"),  lab = "no CC, no TV")
)

# ---------------------------------------------------------------------------
# Función principal: corre un grupo con cualquier prior
# ---------------------------------------------------------------------------
# P(theta1 <= theta2 | x) via 2D trapezoidal quadrature on the posterior grid
compute_ph_x <- function(n1, n2, x1, x2, a0, a1, a2, ngrid = 200) {
  consts <- bb_constants(n1, n2, x1, x2, a0, a1, a2)
  xs <- seq(0.001, 0.999, length.out = ngrid)
  ys <- seq(0.001, 0.999, length.out = ngrid)
  z  <- densBB_grid(xs, ys, consts)
  mask <- outer(seq_len(ngrid), seq_len(ngrid), function(i, j) xs[i] <= ys[j])
  round(sum(z[mask]) / sum(z), 4)
}

run_group <- function(g_entry, a0, a1, a2, M, seed = 7) {
  g <- g_entry$g; lab <- g_entry$lab
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  set.seed(seed)
  ph_x   <- compute_ph_x(n1, n2, x1, x2, a0, a1, a2)
  ev_obs <- ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(n1, n2, x1, x2, a0, a1, a2, M, 401)
  ev_A   <- simulate_evs_A_post(n1, n2, x1, x2, a0, a1, a2, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  data.frame(group  = lab,
             n1 = n1, n2 = n2, x1 = x1, x2 = x2,
             ph_x   = ph_x,
             ev_obs = round(ev_obs,    4),
             k_star = round(opt$k_star,4),
             alpha  = round(opt$alpha, 4),
             beta   = round(opt$beta,  4),
             decision = ifelse(ev_obs <= opt$k_star, "Reject $H$", "Do not reject"))
}

# ---------------------------------------------------------------------------
# Función LaTeX
# ---------------------------------------------------------------------------
to_latex <- function(tab, label, caption) {
  lines <- c(
    "\\begin{table}[!h]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lcccccccccl}",
    "\\toprule",
    paste("Treatment & $P(H|\\mathbf{X})$ & $n_1$ & $n_2$ & $x_1$ & $x_2$",
          "& $\\mathrm{ev}_{\\mathrm{obs}}$ & $k^{*}$",
          "& $\\hat{\\alpha}$ & $\\hat{\\beta}$ & Decision \\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    lines <- c(lines, sprintf(
      "%s & %.4f & %d & %d & %d & %d & %.4f & %.4f & %.4f & %.4f & %s \\\\",
      r$group, r$ph_x, r$n1, r$n2, r$x1, r$x2,
      r$ev_obs, r$k_star, r$alpha, r$beta, r$decision))
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

# ---------------------------------------------------------------------------
# Tab2 — prior KL-óptima (α≈0.84), posterior-based
# ---------------------------------------------------------------------------
cat("\n=== Tab2: prior KL-óptima (M =", M, ") ===\n")
tab2 <- do.call(rbind, lapply(names(grupos), function(k) {
  cat(" -", grupos[[k]]$lab, "\n")
  run_group(grupos[[k]], a0_kl, a1_kl, a2_kl, M)
}))
print(tab2, row.names = FALSE)

tex2 <- to_latex(
  tab2,
  label   = "tab2",
  caption = sprintf(paste0(
    "TVSFP Project: Bayesian hypothesis testing with non-informative prior ",
    "(KL-optimal: $\\alpha_0=%.3f,\\,\\alpha_1=%.3f,\\,\\alpha_2=%.3f$), ",
    "posterior-based formulation, $a=b=1$, $M=%d$."),
    a0_kl, a1_kl, a2_kl, M)
)

# ---------------------------------------------------------------------------
# Tab3 — prior informativa simétrica (KL a Beta(N·0.5, N·0.5)^2, N=50)
# ---------------------------------------------------------------------------
cat("\n=== Tab3: prior informativa simétrica (μ=0.5, N=50) (M =", M, ") ===\n")
tab3 <- do.call(rbind, lapply(names(grupos), function(k) {
  cat(" -", grupos[[k]]$lab, "\n")
  run_group(grupos[[k]], a0_inf, a1_inf, a2_inf, M)
}))
print(tab3, row.names = FALSE)

tex3 <- to_latex(
  tab3,
  label   = "tab3",
  caption = sprintf(paste0(
    "TVSFP Project: Bayesian hypothesis testing with informative prior ",
    "($\\alpha_0=%.2f,\\,\\alpha_1=%.2f,\\,\\alpha_2=%.2f$, ",
    "KL fit to Beta($N\\mu, N(1-\\mu)$)$^2$ with $\\mu=0.5,\\,N=50$), ",
    "posterior-based formulation, $a=b=1$, $M=%d$."),
    a0_inf, a1_inf, a2_inf, M)
)

# ---------------------------------------------------------------------------
# Tab4 — prior en conflicto (KL a Beta(N·0.1, N·0.9)^2, N=50  →  E[θ]≈0.1)
# La prior espera θ≈0.1 pero los datos THKS muestran θ≈0.3-0.6 → conflicto.
# ---------------------------------------------------------------------------
cat("\n=== Tab4: prior en conflicto (μ=0.1, N=50, E[θ]≈0.1) (M =", M, ") ===\n")
tab4 <- do.call(rbind, lapply(names(grupos), function(k) {
  cat(" -", grupos[[k]]$lab, "\n")
  run_group(grupos[[k]], a0_conf, a1_conf, a2_conf, M)
}))
print(tab4, row.names = FALSE)

cat("\n=== P(H|X) para tab4 — filas .tex (reemplazar '---') ===\n")
group_labels <- c("\\textbf{CC - TV}", "\\textbf{CC - No TV}",
                  "\\textbf{No CC - TV}", "\\textbf{No CC - No TV}")
for (i in seq_len(nrow(tab4))) {
  r <- tab4[i, ]
  cat(sprintf("%s & %.4f & %.4f & %.4f & %s & %d  \\\\\n",
              group_labels[i], r$ph_x, r$ev_obs, r$k_star,
              r$decision, r$n1))
}

tex4 <- to_latex(
  tab4,
  label   = "tab4",
  caption = sprintf(paste0(
    "TVSFP Project: Bayesian hypothesis testing with informative prior ",
    "in conflict with data ",
    "($\\alpha_0=%.2f,\\,\\alpha_1=%.2f,\\,\\alpha_2=%.2f$, ",
    "KL fit to Beta($N\\mu, N(1-\\mu)$)$^2$ with $\\mu=0.1,\\,N=50$, ",
    "prior mean $E[\\theta]\\approx 0.1$), ",
    "posterior-based formulation, $a=b=1$, $M=%d$."),
    a0_conf, a1_conf, a2_conf, M)
)

# ---------------------------------------------------------------------------
# Escribir las tres tablas en un único archivo combinado
# ---------------------------------------------------------------------------
all_tables <- c(
  tex2, "",
  tex3, "",
  tex4
)
writeLines(all_tables, "output/thks_tables.tex")
cat("→ output/thks_tables.tex  (tab2 + tab3 + tab4 combinadas)\n")

# ---------------------------------------------------------------------------
# Resumen lado a lado (tres priors)
# ---------------------------------------------------------------------------
cat("\n=== Comparación: KL vs informativa vs conflicto ===\n")
comp <- data.frame(
  group      = tab2$group,
  n          = tab2$n1,
  ev_obs     = tab2$ev_obs,
  kstar_kl   = tab2$k_star,  dec_kl   = tab2$decision,
  kstar_inf  = tab3$k_star,  dec_inf  = tab3$decision,
  kstar_conf = tab4$k_star,  dec_conf = tab4$decision
)
print(comp, row.names = FALSE)
cat("\n→ En el artículo usar: \\input{output/thks_tables.tex}\n")
