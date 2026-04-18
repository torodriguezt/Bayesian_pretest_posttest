# build_kstar_table.R
# Pipeline completo para construir la tabla k*(n1, n2) usando el motor C++.
# Reproduce los pasos del artículo (sec. 2.4):
#   - α(k) = ∮_H π_φe(θ) f_H(θ) dθ   (prior bajo H, restringida a θ1=θ2)
#   - β(k) = ∫∫_{Θ × Ω} f(x|θ) 𝟙(ev>k) f(θ) dx dθ   (prior completa)
#   - k*(n1, n2) = argmin_k [a·α(k) + b·β(k)]

library(Rcpp)
library(ggplot2)
library(tidyr)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")

# --- 1) Hiperparámetros (prior no informativa por KL del artículo) ---------
alphas_opt <- c(0.8373879, 0.8410984, 0.8053298)
a0 <- alphas_opt[1]; a1 <- alphas_opt[2]; a2 <- alphas_opt[3]

# --- 2) Curvas α(k), β(k) y k* para un (n1, n2) dado ----------------------
analyze_pair <- function(n1, n2, M = 500, a_w = 1, b_w = 1, seed = 42,
                         ngrid_quad = 401) {
  set.seed(seed)
  ev_H <- simulate_evs_H(n1, n2, a0, a1, a2, M, ngrid_quad)
  ev_A <- simulate_evs_A(n1, n2, a0, a1, a2, M, ngrid_quad)
  k_grid <- seq(0, 1, length.out = 401)
  curves <- error_curves(ev_H, ev_A, k_grid, a_w, b_w)
  opt    <- find_kstar(ev_H, ev_A, a_w, b_w)
  list(curves = curves, opt = opt, ev_H = ev_H, ev_A = ev_A,
       n1 = n1, n2 = n2)
}

plot_errors <- function(res, save_path = NULL) {
  df <- res$curves
  long <- pivot_longer(df, c("alpha", "beta", "sum"),
                       names_to = "metrica", values_to = "valor")
  long$metrica <- factor(long$metrica,
                         levels = c("alpha", "beta", "sum"),
                         labels = c(expression(alpha[varphi[e]]),
                                    expression(beta[varphi[e]]),
                                    expression(alpha + beta)))
  p <- ggplot(long, aes(k, valor, color = metrica, linetype = metrica)) +
    geom_line(size = 0.8) +
    geom_vline(xintercept = res$opt$k_star, linetype = "dashed", color = "grey40") +
    annotate("text", x = res$opt$k_star, y = 0.05,
             label = sprintf("k* = %.3f", res$opt$k_star),
             hjust = -0.1, size = 3.5) +
    scale_color_manual(values = c("steelblue", "tomato", "black"),
                       labels = c(expression(alpha), expression(beta),
                                  expression(alpha+beta))) +
    scale_linetype_manual(values = c("solid", "solid", "dashed"),
                          labels = c(expression(alpha), expression(beta),
                                     expression(alpha+beta))) +
    labs(x = "k", y = "Error promediado",
         title = sprintf("n1 = %d, n2 = %d", res$n1, res$n2),
         color = NULL, linetype = NULL) +
    theme_bw() + theme(legend.position = "top")
  if (!is.null(save_path)) ggsave(save_path, p, width = 6, height = 4, dpi = 150)
  p
}

# --- 3) Tabla k*(n1, n2) ---------------------------------------------------
build_kstar_table <- function(n_grid, M = 500, a_w = 1, b_w = 1, seed = 42,
                               ngrid_quad = 401, verbose = TRUE) {
  N <- length(n_grid)
  K <- matrix(NA_real_, N, N, dimnames = list(n_grid, n_grid))
  Alpha <- K; Beta <- K
  set.seed(seed)
  total <- N * N; done <- 0
  t_start <- Sys.time()
  for (i in seq_along(n_grid)) {
    for (j in seq_along(n_grid)) {
      n1 <- n_grid[i]; n2 <- n_grid[j]
      ev_H <- simulate_evs_H(n1, n2, a0, a1, a2, M, ngrid_quad)
      ev_A <- simulate_evs_A(n1, n2, a0, a1, a2, M, ngrid_quad)
      opt  <- find_kstar(ev_H, ev_A, a_w, b_w)
      K[i, j]     <- opt$k_star
      Alpha[i, j] <- opt$alpha
      Beta[i, j]  <- opt$beta
      done <- done + 1
      if (verbose) {
        elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
        eta <- elapsed * (total - done) / done
        cat(sprintf("[%3d/%3d] n1=%d n2=%d  k*=%.4f  α=%.4f  β=%.4f  | ETA %.0fs\n",
                    done, total, n1, n2, opt$k_star, opt$alpha, opt$beta, eta))
      }
    }
  }
  list(k_star = K, alpha = Alpha, beta = Beta, n_grid = n_grid,
       M = M, a = a_w, b = b_w)
}

# --- 4) Exportar tabla a LaTeX ---------------------------------------------
table_to_latex <- function(tab, label = "tab:kstar",
                           caption = "Cutoff óptimo $k^{*}(n_1, n_2)$ para $a=b=1$.",
                           digits = 4) {
  K <- tab$k_star
  ng <- tab$n_grid
  N <- length(ng)
  col_spec <- paste(rep("c", N + 1), collapse = "")
  out <- c(
    "\\begin{table}[!h]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0("\\diagbox{$n_1$}{$n_2$} & ",
           paste(ng, collapse = " & "), " \\\\"),
    "\\midrule"
  )
  for (i in seq_len(N)) {
    row_vals <- formatC(K[i, ], digits = digits, format = "f")
    out <- c(out, paste0(ng[i], " & ", paste(row_vals, collapse = " & "), " \\\\"))
  }
  out <- c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  paste(out, collapse = "\n")
}

# ===========================================================================
# DEMO 1: curvas α(k), β(k) para un par (n1, n2) (replica fig. 9-12 del paper)
# ===========================================================================
cat("\n--- Demo 1: curvas de error para n1=n2=25 ---\n")
res25 <- analyze_pair(n1 = 25, n2 = 25, M = 500)
cat(sprintf("k* = %.4f   α = %.4f   β = %.4f\n",
            res25$opt$k_star, res25$opt$alpha, res25$opt$beta))
print(plot_errors(res25, "Figures/error_curves_n25.png"))

# ===========================================================================
# DEMO 2: tabla k*(n1, n2) en una grilla chica para inspección rápida
# ===========================================================================
cat("\n--- Demo 2: tabla k*(n1, n2) en grilla {10, 20, 40, 80} ---\n")
tab_demo <- build_kstar_table(n_grid = c(10, 20, 40, 80), M = 300, seed = 1)
cat("\nMatriz k*:\n")
print(round(tab_demo$k_star, 4))

cat("\nLaTeX:\n")
cat(table_to_latex(tab_demo,
    caption = "Cutoff óptimo $k^{*}(n_1, n_2)$ — demo, M=300."), "\n")

# ===========================================================================
# DEMO 3 (lento, descomentar): tabla completa para el artículo
# ===========================================================================
# cat("\n--- Demo 3: tabla completa (puede tardar) ---\n")
# n_grid_full <- c(10, 20, 30, 40, 50, 75, 100, 150, 200)
# tab_full <- build_kstar_table(n_grid_full, M = 1000, seed = 123)
# saveRDS(tab_full, "kstar_table_full.rds")
# writeLines(
#   table_to_latex(tab_full,
#     caption = "Cutoff óptimo $k^{*}(n_1, n_2)$ obtenido por simulación con $M=1000$ datasets bajo $H$ y bajo la prior completa, $a=b=1$.",
#     label = "tab:kstar_full"),
#   "kstar_table_full.tex"
# )
# cat("Tabla guardada en kstar_table_full.tex\n")

# ===========================================================================
# DEMO 4 (opcional): aplicar a los 4 grupos del THKS y comparar k* con el artículo
# ===========================================================================
# library(ALA); library(dplyr)
# datos1 <- tvsfp
# prep_school <- function(d, sid, sb, tv) {
#   x <- d %>% mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
#         filter(school == sid, school.based == sb, tv.based == tv)
#   list(X = (x %>% group_by(stage) %>% summarise(B=sum(binTHKS)))$B,
#        n = (x %>% group_by(stage) %>% summarise(n=n()))$n)
# }
# grupos <- list(
#   yy = prep_school(datos1, "404", "yes", "yes"),
#   yn = prep_school(datos1, "408", "yes", "no"),
#   nn = prep_school(datos1, "409", "no",  "no")
# )
# for (k in names(grupos)) {
#   g <- grupos[[k]]; n1 <- g$n[1]; n2 <- g$n[2]
#   res <- analyze_pair(n1, n2, M = 1000, seed = 7)
#   ev_obs <- ev_quad_from_data(n1, n2, g$X[1], g$X[2], a0, a1, a2)
#   cat(sprintf("[%s] n1=%d n2=%d  k*=%.4f  ev_obs=%.4f  decision=%s\n",
#               k, n1, n2, res$opt$k_star, ev_obs,
#               ifelse(ev_obs <= res$opt$k_star, "rechaza H", "no rechaza")))
# }
