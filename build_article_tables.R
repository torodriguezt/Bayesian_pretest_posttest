# build_article_tables.R
# Genera las DOS tablas para el artículo y las guarda como .tex listos para \input{}.
#
#   Tabla A — k*(n1, n2)  PRIOR-based (genérica, no depende de los datos)
#       Para cualquier futuro experimento con tamaños (n1, n2), da el corte
#       óptimo que minimiza α(k) + β(k), con α y β promediados bajo la prior
#       restringida a H y bajo la prior completa, respectivamente.
#       → guarda en: kstar_prior_table.tex
#
#   Tabla B — k* por tratamiento POSTERIOR-based (específica del THKS)
#       Para cada uno de los 4 grupos del estudio, da el corte óptimo cuando
#       α y β se promedian bajo la POSTERIOR (dados los datos observados).
#       Reemplaza tab2 del artículo (cuyos valores parecen previos a la
#       corrección de la integral de línea / posterior).
#       → guarda en: kstar_posterior_table.tex

library(ALA)
library(dplyr)
library(Rcpp)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")

source("priors_config.R")
a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]

source("build_kstar_table.R", echo = FALSE, max.deparse.length = Inf,
       local = new.env(parent = .GlobalEnv))   # define build_kstar_table, table_to_latex
# (re-cargamos defs sin re-correr los demos)
build_kstar_table <- function(n_grid, M = 500, a_w = 1, b_w = 1, seed = 42,
                               ngrid_quad = 401, verbose = TRUE) {
  N <- length(n_grid)
  K <- matrix(NA_real_, N, N, dimnames = list(n_grid, n_grid))
  Alpha <- K; Beta <- K
  set.seed(seed)
  total <- N * N; done <- 0; t_start <- Sys.time()
  for (i in seq_along(n_grid)) for (j in seq_along(n_grid)) {
    n1 <- n_grid[i]; n2 <- n_grid[j]
    ev_H <- simulate_evs_H(n1, n2, a0, a1, a2, M, ngrid_quad)
    ev_A <- simulate_evs_A(n1, n2, a0, a1, a2, M, ngrid_quad)
    opt  <- find_kstar(ev_H, ev_A, a_w, b_w)
    K[i, j] <- opt$k_star; Alpha[i, j] <- opt$alpha; Beta[i, j] <- opt$beta
    done <- done + 1
    if (verbose) {
      el <- as.numeric(Sys.time() - t_start, units = "secs")
      cat(sprintf("[%3d/%3d] n1=%d n2=%d  k*=%.4f  α=%.4f  β=%.4f  | ETA %.0fs\n",
                  done, total, n1, n2, opt$k_star, opt$alpha, opt$beta,
                  el * (total - done) / done))
    }
  }
  list(k_star = K, alpha = Alpha, beta = Beta, n_grid = n_grid,
       M = M, a = a_w, b = b_w)
}

table_to_latex <- function(tab, label = "tab:kstar",
                           caption = "Cutoff óptimo $k^{*}(n_1, n_2)$ para $a=b=1$.",
                           digits = 4) {
  K <- tab$k_star; ng <- tab$n_grid; N <- length(ng)
  col_spec <- paste(rep("c", N + 1), collapse = "")
  out <- c("\\begin{table}[!h]", "\\centering",
           paste0("\\caption{", caption, "}"),
           paste0("\\label{", label, "}"),
           paste0("\\begin{tabular}{", col_spec, "}"),
           "\\toprule",
           paste0("\\diagbox{$n_1$}{$n_2$} & ",
                  paste(ng, collapse = " & "), " \\\\"),
           "\\midrule")
  for (i in seq_len(N)) {
    row_vals <- formatC(K[i, ], digits = digits, format = "f")
    out <- c(out, paste0(ng[i], " & ", paste(row_vals, collapse = " & "), " \\\\"))
  }
  c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

# ===========================================================================
# Tabla A — prior-based, genérica
# ===========================================================================
cat("\n========== Tabla A: k*(n1, n2) prior-based ==========\n")
n_grid_full <- c(10, 20, 30, 40, 50, 75, 100, 150, 200)
M_prior <- 2000  # con M=500 el k* es muy ruidoso (lo decide la cola de ev_A);
                 # M=2000 reduce el error MC en ~√4 = 2× respecto a M=500.

t0 <- Sys.time()
tabA <- build_kstar_table(n_grid_full, M = M_prior, seed = 123, verbose = TRUE)
cat(sprintf("\nTabla A completada en %.1f min.\n",
            as.numeric(Sys.time() - t0, units = "mins")))

cat("\nMatriz k* (prior-based):\n"); print(round(tabA$k_star, 4))

texA <- table_to_latex(
  tabA, label = "tab:kstar_prior",
  caption = sprintf(
    paste0("Cutoff óptimo $k^{*}(n_1, n_2)$ para la formulación con prior ",
           "(ec. 10-11): $\\alpha$ y $\\beta$ promediadas bajo la prior ",
           "restringida a $H$ y la prior completa, $a=b=1$, $M=%d$ datasets."),
    M_prior))
dir.create("output", showWarnings = FALSE)
writeLines(texA, "output/kstar_prior_table.tex")
saveRDS(tabA, "output/kstar_prior_table.rds")
cat("\n→ Tabla A guardada en output/kstar_prior_table.tex\n")

# ===========================================================================
# Tabla B — posterior-based, específica del THKS
# ===========================================================================
cat("\n========== Tabla B: k* por tratamiento posterior-based ==========\n")

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
  yy = prep_school(datos1, "404", "yes", "yes"),
  yn = prep_school(datos1, "408", "yes", "no"),
  ny = prep_global(datos1,         "no",  "yes"),
  nn = prep_school(datos1, "409", "no",  "no")
)
labels <- c(yy = "CC + TV", yn = "CC, no TV",
            ny = "no CC, TV", nn = "no CC, no TV")

M_post <- 2000

post_row <- function(g, label, M = M_post, seed = 7, ngrid_quad = 401) {
  set.seed(seed)
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  ev_obs <- ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2, ngrid_quad)
  ess    <- sir_ess(n1, n2, x1, x2, a0, a1, a2, N_prop = max(50 * M, 10000))
  ev_H   <- simulate_evs_H_post(n1, n2, x1, x2, a0, a1, a2, M, ngrid_quad)
  ev_A   <- simulate_evs_A_post(n1, n2, x1, x2, a0, a1, a2, M, ngrid_quad)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  data.frame(group = label, n1 = n1, n2 = n2, x1 = x1, x2 = x2,
             ev_obs = ev_obs, k_star = opt$k_star,
             alpha = opt$alpha, beta = opt$beta, ess = round(ess, 0))
}

cat("Procesando 4 grupos del THKS (M =", M_post, ")...\n")
t0 <- Sys.time()
tabB <- do.call(rbind, lapply(names(grupos), function(k) {
  cat("  -", labels[[k]], "\n"); post_row(grupos[[k]], labels[[k]])
}))
cat(sprintf("Tabla B completada en %.1f min.\n",
            as.numeric(Sys.time() - t0, units = "mins")))

print(tabB, row.names = FALSE, digits = 4)

# Decisión: rechaza H si ev_obs ≤ k*
tabB$decision <- ifelse(tabB$ev_obs <= tabB$k_star, "rechaza H", "no rechaza")

texB <- c(
  "\\begin{table}[!h]",
  "\\centering",
  paste0("\\caption{Cutoff óptimo $k^{*}$ por tratamiento (formulación ",
         "posterior, ec. 16-17), $a=b=1$, $M=", M_post, "$ datasets ",
         "simulados desde la posterior dada $(x_1, x_2)$.}"),
  "\\label{tab:kstar_posterior}",
  "\\begin{tabular}{lccccccc}",
  "\\toprule",
  "Tratamiento & $n_1$ & $n_2$ & $x_1$ & $x_2$ & $\\mathrm{ev}_{\\mathrm{obs}}$ & $k^{*}$ & Decisión \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(tabB))) {
  r <- tabB[i, ]
  texB <- c(texB, sprintf("%s & %d & %d & %d & %d & %.4f & %.4f & %s \\\\",
                          r$group, r$n1, r$n2, r$x1, r$x2,
                          r$ev_obs, r$k_star, r$decision))
}
texB <- c(texB, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(texB, "output/kstar_posterior_table.tex")
saveRDS(tabB, "output/kstar_posterior_table.rds")
cat("\n→ Tabla B guardada en output/kstar_posterior_table.tex\n")

cat("\n========== Resumen ==========\n")
cat("kstar_prior_table.tex     →  Tabla genérica k*(n1,n2), prior-based\n")
cat("kstar_posterior_table.tex →  Tabla por tratamiento, posterior-based\n")
cat("\nNota: las asignaciones de escuela en los 4 grupos vienen de THKS_run.R.\n")
cat("Si el artículo usa otra agrupación, hay que ajustar `grupos` y re-correr\n")
cat("sólo la Tabla B (la A es genérica y no depende de los datos).\n")
