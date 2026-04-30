# thks_fbst_mcnemar.R
# FBST (3 priors) + hipotesis no precisa + McNemar — datos reales TVSFP.
# Grupos definidos SOLO por condicion de tratamiento (school.based x tv.based),
# sin filtrar por escuela individual.
#
# Salida: output/thks_fbst_mcnemar.tex

library(ALA)
library(dplyr)
library(tidyr)
library(Rcpp)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")
source("priors_config.R")

M_post <- 2000
set.seed(7)

priors <- list(
  list(label = "Non-informative", a0 = prior_NI["a0"],   a1 = prior_NI["a1"],   a2 = prior_NI["a2"]),
  list(label = "Informative",     a0 = prior_INF["a0"],  a1 = prior_INF["a1"],  a2 = prior_INF["a2"]),
  list(label = "Conflict",        a0 = prior_CONF["a0"], a1 = prior_CONF["a1"], a2 = prior_CONF["a2"])
)

# ============================================================================
# 1. EXTRACCION DE PARES COMPLETOS (sin filtro por escuela)
# ============================================================================

extract_pairs <- function(d, label, sb, tv) {
  wide <- d %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school.based == sb, tv.based == tv) %>%
    select(id, stage, binTHKS) %>%
    pivot_wider(names_from = stage, values_from = binTHKS) %>%
    filter(!is.na(pre), !is.na(post))

  list(
    label = label,
    n     = nrow(wide),
    x1    = sum(wide$pre),
    x2    = sum(wide$post),
    n00   = sum(wide$pre == 0 & wide$post == 0),
    n01   = sum(wide$pre == 0 & wide$post == 1),
    n10   = sum(wide$pre == 1 & wide$post == 0),
    n11   = sum(wide$pre == 1 & wide$post == 1)
  )
}

datos1 <- tvsfp

grupos_raw <- list(
  extract_pairs(datos1, "CC + TV",      sb = "yes", tv = "yes"),
  extract_pairs(datos1, "CC, no TV",    sb = "yes", tv = "no"),
  extract_pairs(datos1, "no CC, TV",    sb = "no",  tv = "yes"),
  extract_pairs(datos1, "no CC, no TV", sb = "no",  tv = "no")
)

# ============================================================================
# 2. McNEMAR (no depende de la prior)
# ============================================================================

mcnemar_pval <- function(n01, n10) {
  bc <- n01 + n10
  if (bc == 0L) return(1.0)
  if (bc < 25L) {
    2 * min(pbinom(n01, bc, 0.5), pbinom(n10, bc, 0.5))
  } else {
    pchisq((n01 - n10)^2 / bc, df = 1, lower.tail = FALSE)
  }
}

# ============================================================================
# 3. FBST + HIPOTESIS NO PRECISA por prior
# ============================================================================

run_fbst <- function(g, prior) {
  a0 <- prior$a0; a1 <- prior$a1; a2 <- prior$a2

  # Hipotesis precisa: FBST
  ev_obs <- ev_quad_from_data(g$n, g$n, g$x1, g$x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(g$n, g$n, g$x1, g$x2, a0, a1, a2, M_post)
  ev_A   <- simulate_evs_A_post(g$n, g$n, g$x1, g$x2, a0, a1, a2, M_post)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)

  # Hipotesis no precisa: P(theta1 <= theta2 | X) via muestras SIR
  samp      <- sample_posterior(M_post, g$n, g$n, g$x1, g$x2, a0, a1, a2)
  prob_H_np <- mean(samp[, 1] <= samp[, 2])

  list(
    ev_obs   = round(ev_obs, 4),
    k_star   = round(opt$k_star, 4),
    reject   = ev_obs <= opt$k_star,
    prob_H   = round(prob_H_np, 4)
  )
}

# ============================================================================
# 4. LOOP PRINCIPAL
# ============================================================================

cat("\n========== FBST (3 priors) + McNemar — TVSFP ==========\n")

rows <- list()

for (g in grupos_raw) {
  pval    <- mcnemar_pval(g$n01, g$n10)
  mcn_dec <- ifelse(pval < 0.05, "Reject", "No reject")

  cat(sprintf("\n--- %s (n=%d, x1=%d, x2=%d | n01=%d, n10=%d) ---\n",
              g$label, g$n, g$x1, g$x2, g$n01, g$n10))
  cat(sprintf("  McNemar: p=%.4f → %s\n", pval, mcn_dec))

  row <- list(
    group   = g$label,
    n       = g$n,
    x1      = g$x1,
    x2      = g$x2,
    n01     = g$n01,
    n10     = g$n10,
    p_mcn   = round(pval, 4),
    mcn_dec = mcn_dec
  )

  for (pr in priors) {
    res <- run_fbst(g, pr)
    dec <- ifelse(res$reject, "Reject", "No reject")
    cat(sprintf("  FBST %-18s P(H_np)=%.4f | ev=%.4f k*=%.4f → %s\n",
                paste0("[", pr$label, "]:"),
                res$prob_H, res$ev_obs, res$k_star, dec))
    key <- gsub("[ -]", "_", tolower(pr$label))
    row[[paste0("prob_H_", key)]] <- res$prob_H
    row[[paste0("ev_",     key)]] <- res$ev_obs
    row[[paste0("ks_",     key)]] <- res$k_star
    row[[paste0("dec_",    key)]] <- dec
  }

  rows[[length(rows) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)

cat("\n========== RESUMEN ==========\n")
print(tab, row.names = FALSE)

# ============================================================================
# 5. TABLAS LATEX
# ============================================================================

dir.create("output", showWarnings = FALSE)

dec_sym <- function(dec) ifelse(dec == "Reject", "Reject $H$", "Do not reject")

# --- Helper: tabla FBST para una prior ------------------------------------
write_fbst_table <- function(tab, key, prior_label, prior_spec,
                              label, filename) {
  ev_col  <- paste0("ev_",     key)
  ks_col  <- paste0("ks_",     key)
  dec_col <- paste0("dec_",    key)
  ph_col  <- paste0("prob_H_", key)

  tex <- c(
    "\\begin{table}[!h]",
    "\\centering",
    paste0("\\caption{TVSFP: FBST results under the ", prior_label,
           " prior (", prior_spec, "), posterior-based formulation, $a=b=1$, $M=",
           M_post, "$. ",
           "$P(H_{\\mathrm{np}}\\mid\\mathbf{X}) = P(\\theta_1 \\leq \\theta_2 \\mid \\mathbf{X})$.}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lccccccc}",
    "\\toprule",
    paste0("\\textbf{Group} & $n$ & $x_1$ & $x_2$ & ",
           "$P(H_{\\mathrm{np}}\\mid\\mathbf{X})$ & ",
           "$ev(\\mathbf{H};\\mathbf{X})$ & $k^*$ & Decision \\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    tex <- c(tex, sprintf(
      "%s & %d & %d & %d & %.4f & %.4f & %.4f & %s \\\\",
      r$group, r$n, r$x1, r$x2,
      r[[ph_col]], r[[ev_col]], r[[ks_col]],
      dec_sym(r[[dec_col]])
    ))
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, filename)
  cat(sprintf("→ %s\n", filename))
}

# Tabla 1: prior no informativa
write_fbst_table(tab,
  key         = "non_informative",
  prior_label = "non-informative KL-optimal",
  prior_spec  = "$\\alpha_0=0.76,\\,\\alpha_1=0.76,\\,\\alpha_2=0.76$",
  label       = "tab:thks_ni",
  filename    = "output/thks_fbst_ni.tex")

# Tabla 2: prior informativa
write_fbst_table(tab,
  key         = "informative",
  prior_label = "informative",
  prior_spec  = "$\\alpha_0=24.99,\\,\\alpha_1=24.80,\\,\\alpha_2=24.88$",
  label       = "tab:thks_inf",
  filename    = "output/thks_fbst_inf.tex")

# Tabla 3: prior en conflicto
write_fbst_table(tab,
  key         = "conflict",
  prior_label = "informative (conflict)",
  prior_spec  = "$\\alpha_0=45.87,\\,\\alpha_1=5.11,\\,\\alpha_2=5.07$",
  label       = "tab:thks_conf",
  filename    = "output/thks_fbst_conf.tex")

# Tabla 4: McNemar
tex_mcn <- c(
  "\\begin{table}[!h]",
  "\\centering",
  paste0("\\caption{TVSFP: McNemar's test for each treatment group. ",
         "$n_{01}$: pre$=0\\to$post$=1$ (improvement); ",
         "$n_{10}$: pre$=1\\to$post$=0$ (deterioration).}"),
  "\\label{tab:thks_mcnemar}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "\\textbf{Group} & $n$ & $n_{01}$ & $n_{10}$ & $p$-value & Decision \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  tex_mcn <- c(tex_mcn, sprintf(
    "%s & %d & %d & %d & %.4f & %s \\\\",
    r$group, r$n, r$n01, r$n10, r$p_mcn, dec_sym(r$mcn_dec)
  ))
}
tex_mcn <- c(tex_mcn, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_mcn, "output/thks_mcnemar.tex")
cat("→ output/thks_mcnemar.tex\n")
