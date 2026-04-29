# simulation_comparison_priors.R
# Estudio de simulación: efecto de tamaño de muestra (n) y tamaño de efecto (d)
# Compara TRES prioris (todas obtenidas via KL minimization, ver priors_config.R):
#   - No informativa (KL a la uniforme)
#   - Informativa simétrica (KL a Beta(N·0.5, N·0.5)^2, N=50)
#   - Informativa con conflicto (KL a Beta(N·0.1, N·0.9)^2, N=50, E[θ]≈0.1)
#
# Salida:
#   - Tabla resumen (n, d, P(H), ev_obs, k*, decisión para los 3 priors)
#   - Gráficos: curvas de error, heatmaps de decisión

library(ALA)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(tidyr)
library(gridExtra)

sourceCpp("BivBetaBinom.cpp")
source("priors_config.R")

dir.create("output",  showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# ============================================================================
# 1. PRIORIS A COMPARAR
# ============================================================================

priors <- list(
  "KL-optimal"            = prior_NI,
  "Informative (mu=0.5)"  = prior_INF,
  "Conflict (mu=0.1)"     = prior_CONF
)
prior_levels <- names(priors)

cat("Prioris (KL fit):\n"); print(priors_summary); cat("\n")

# ============================================================================
# 2. ESCENARIOS
# ============================================================================

scenarios <- expand.grid(
  n             = c(30, 50, 75, 100, 150),
  scenario_type = c("null", "small", "medium", "large")
) %>%
  mutate(
    theta1  = 0.5,
    theta2  = case_when(
      scenario_type == "null"   ~ 0.50,
      scenario_type == "small"  ~ 0.52,
      scenario_type == "medium" ~ 0.60,
      scenario_type == "large"  ~ 0.80
    ),
    d       = theta2 - theta1,
    d_label = case_when(
      scenario_type == "null"   ~ "0.00 (H true)",
      scenario_type == "small"  ~ "0.02 (small)",
      scenario_type == "medium" ~ "0.10 (medium)",
      scenario_type == "large"  ~ "0.30 (large)"
    )
  )

M_sim <- 1000
set.seed(42)

# ============================================================================
# 3. FUNCION: analizar un escenario con un prior
# ============================================================================

analyze_scenario <- function(n, x1, x2, prior_alpha, M = 1000) {
  a0 <- prior_alpha["a0"]; a1 <- prior_alpha["a1"]; a2 <- prior_alpha["a2"]
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  sup_H  <- find_sup_H(consts)$sup_H
  ev_obs <- ev_quad(consts, sup_H)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  post   <- sample_posterior(5000, n, n, x1, x2, a0, a1, a2)
  list(
    ev_obs   = ev_obs,
    k_star   = opt$k_star,
    alpha    = opt$alpha,
    beta     = opt$beta,
    decision = ifelse(ev_obs <= opt$k_star, "Reject H", "Do not reject"),
    prob_H   = mean(post[, 1] <= post[, 2]),
    ev_H     = ev_H,
    ev_A     = ev_A
  )
}

# ============================================================================
# 4. LOOP PRINCIPAL: tabla resumen
# ============================================================================

cat("========== SIMULATION: PRIOR COMPARISON ==========\n\n")

results <- list()
results_full <- list()  # incluye ev_H, ev_A para los plots

for (i in seq_len(nrow(scenarios))) {
  row <- scenarios[i, ]
  cat(sprintf("n=%d, theta1=%.2f, theta2=%.2f (%s)... ",
              row$n, row$theta1, row$theta2, row$d_label))

  # Una sola realización de datos por escenario (compartida entre los 3 priors)
  x1 <- rbinom(1, row$n, row$theta1)
  x2 <- rbinom(1, row$n, row$theta2)

  for (pname in prior_levels) {
    res <- analyze_scenario(row$n, x1, x2, priors[[pname]], M = M_sim)
    results[[length(results) + 1]] <- data.frame(
      n        = row$n,
      theta1   = row$theta1,
      theta2   = row$theta2,
      d        = row$d,
      d_label  = row$d_label,
      x1 = x1, x2 = x2,
      prior    = pname,
      prob_H   = round(res$prob_H, 3),
      ev_obs   = round(res$ev_obs, 4),
      k_star   = round(res$k_star, 4),
      alpha    = round(res$alpha,  4),
      beta     = round(res$beta,   4),
      decision = res$decision
    )
    results_full[[paste(row$n, row$d, pname, sep = "|")]] <- res
  }
  cat("OK\n")
}
results <- do.call(rbind, results)
results$prior <- factor(results$prior, levels = prior_levels)

cat("\n========== SUMMARY ==========\n")
print(results[, c("n", "d", "prior", "prob_H", "ev_obs", "k_star", "decision")])
write.csv(results, "output/simulation_comparison.csv", row.names = FALSE)

# ============================================================================
# 5. TABLA LATEX  (lado a lado, 3 prioris)
# ============================================================================

results_wide <- results %>%
  select(n, d, prior, prob_H, ev_obs, k_star, decision) %>%
  pivot_wider(names_from = prior,
              values_from = c(prob_H, ev_obs, k_star, decision))

tex <- c(
  "\\begin{table}[!h]", "\\centering", "\\footnotesize",
  paste0("\\caption{Simulation study: effect of sample size ($n$) and effect ",
         "size ($\\Delta = \\theta_2 - \\theta_1$) on FBST decision. ",
         "Three priors fit by KL: non-informative (KL$\\to$U), informative ",
         "($\\mu=0.5$, $N=50$), conflict ($\\mu=0.1$, $N=50$). $M=1000$.}"),
  "\\label{tab:simulation}",
  "\\begin{tabular}{cc|cccc|cccc|cccc}",
  "\\toprule",
  "\\multicolumn{2}{c}{} & \\multicolumn{4}{c}{\\textbf{Non-informative}} & \\multicolumn{4}{c}{\\textbf{Informative}} & \\multicolumn{4}{c}{\\textbf{Conflict}} \\\\",
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}\\cmidrule(lr){11-14}",
  "$n$ & $\\Delta$ & $P(H)$ & $ev$ & $k^*$ & Dec & $P(H)$ & $ev$ & $k^*$ & Dec & $P(H)$ & $ev$ & $k^*$ & Dec \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(results_wide))) {
  r <- results_wide[i, ]
  short <- function(x) ifelse(x == "Reject H", "R", "NR")
  tex <- c(tex, sprintf(
    paste0("%d & %.2f & %.3f & %.4f & %.4f & %s ",
           "& %.3f & %.4f & %.4f & %s ",
           "& %.3f & %.4f & %.4f & %s \\\\"),
    r$n, r$d,
    r$`prob_H_KL-optimal`,            r$`ev_obs_KL-optimal`,            r$`k_star_KL-optimal`,            short(r$`decision_KL-optimal`),
    r$`prob_H_Informative (mu=0.5)`,  r$`ev_obs_Informative (mu=0.5)`,  r$`k_star_Informative (mu=0.5)`,  short(r$`decision_Informative (mu=0.5)`),
    r$`prob_H_Conflict (mu=0.1)`,     r$`ev_obs_Conflict (mu=0.1)`,     r$`k_star_Conflict (mu=0.1)`,     short(r$`decision_Conflict (mu=0.1)`)
  ))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}",
         "\\\\\\smallskip\\footnotesize R = Reject H, NR = Do not reject.",
         "\\end{table}")
writeLines(tex, "output/simulation_table.tex")
cat("\n--> output/simulation_table.tex\n")

# ============================================================================
# 6. CURVAS DE ERROR para escenarios representativos
# ============================================================================

cat("\n========== PLOTS ==========\n")

plot_error_curves_3priors <- function(n_val, d_val) {
  k_grid <- seq(0, 1, length.out = 401)
  df_all <- list()
  k_stars <- numeric()
  for (pname in prior_levels) {
    key <- paste(n_val, d_val, pname, sep = "|")
    if (is.null(results_full[[key]])) next
    res <- results_full[[key]]
    cur <- error_curves(res$ev_H, res$ev_A, k_grid, 1, 1)
    cur$prior <- pname
    df_all[[pname]] <- cur
    k_stars[pname] <- res$k_star
  }
  if (length(df_all) == 0) return(NULL)
  curves <- do.call(rbind, df_all)
  curves$prior <- factor(curves$prior, levels = prior_levels)
  long <- pivot_longer(curves, c("alpha","beta","sum"),
                       names_to = "metric", values_to = "value")
  long$metric <- factor(long$metric,
                        levels = c("alpha","beta","sum"),
                        labels = c("alpha (Type I)", "beta (Type II)", "alpha+beta"))
  vlines <- data.frame(prior = factor(names(k_stars), levels = prior_levels),
                       k_star = unname(k_stars))
  ggplot(long, aes(k, value, color = metric, linetype = metric)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~prior, ncol = 3) +
    geom_vline(data = vlines, aes(xintercept = k_star),
               linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("steelblue","tomato","black")) +
    scale_linetype_manual(values = c("solid","solid","dashed")) +
    labs(x = "k", y = "Average error",
         title = sprintf("Error curves: n=%d, Delta=%.2f", n_val, d_val),
         color = NULL, linetype = NULL) +
    theme_bw() + theme(legend.position = "top")
}

scenarios_to_plot <- expand.grid(n = c(50, 100, 150),
                                 d = c(0.10, 0.30))
for (i in seq_len(nrow(scenarios_to_plot))) {
  n_val <- scenarios_to_plot[i, "n"]
  d_val <- scenarios_to_plot[i, "d"]
  cat(sprintf("  curves: n=%d, d=%.2f\n", n_val, d_val))
  p <- plot_error_curves_3priors(n_val, d_val)
  if (!is.null(p)) {
    ggsave(sprintf("Figures/error_curves_comparison_n%d_d%02d.png",
                   n_val, round(d_val * 100)),
           p, width = 12, height = 4, dpi = 150)
  }
}

# ============================================================================
# 7. HEATMAPS de decisión (uno por prior)
# ============================================================================

cat("  decision heatmap\n")
heatmaps <- lapply(prior_levels, function(pname) {
  ggplot(results %>% filter(prior == pname),
         aes(x = factor(d), y = factor(n), fill = decision)) +
    geom_tile(color = "black", size = 0.6) +
    scale_fill_manual(values = c("Reject H" = "tomato",
                                  "Do not reject" = "steelblue")) +
    labs(x = "Effect size (Delta)", y = "Sample size (n)",
         title = pname, fill = "Decision") +
    theme_bw()
})
p_heat <- do.call(gridExtra::grid.arrange, c(heatmaps, ncol = 3))
ggsave("Figures/decision_heatmap_comparison.png", p_heat,
       width = 14, height = 5, dpi = 150)

# ============================================================================
# 8. ev_obs vs k* (scatter, una faceta por prior)
# ============================================================================

cat("  scatter ev_obs vs k*\n")
p_scatter <- ggplot(results, aes(x = ev_obs, y = k_star,
                                 color = prior, shape = decision)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_grid(~prior) +
  scale_color_manual(values = c("KL-optimal" = "steelblue",
                                "Informative (mu=0.5)" = "tomato",
                                "Conflict (mu=0.1)" = "darkgreen")) +
  scale_shape_manual(values = c("Reject H" = 17, "Do not reject" = 16)) +
  labs(x = "Observed e-value", y = "k* (adaptive cutoff)",
       title = "ev_obs vs k* across priors",
       color = "Prior", shape = "Decision") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Figures/ev_vs_kstar_comparison.png", p_scatter,
       width = 12, height = 5, dpi = 150)

# ============================================================================
# 9. RESUMEN FINAL: acuerdo entre prioris
# ============================================================================

cat("\n========== AGREEMENT MATRIX ==========\n")
agree <- results %>%
  select(n, d, prior, decision) %>%
  pivot_wider(names_from = prior, values_from = decision) %>%
  mutate(
    NI_vs_INF  = `KL-optimal` == `Informative (mu=0.5)`,
    NI_vs_CONF = `KL-optimal` == `Conflict (mu=0.1)`,
    INF_vs_CONF = `Informative (mu=0.5)` == `Conflict (mu=0.1)`
  )
cat(sprintf("  Non-inf vs Informative:  %.1f%% agreement\n",
            mean(agree$NI_vs_INF) * 100))
cat(sprintf("  Non-inf vs Conflict:     %.1f%% agreement\n",
            mean(agree$NI_vs_CONF) * 100))
cat(sprintf("  Informative vs Conflict: %.1f%% agreement\n",
            mean(agree$INF_vs_CONF) * 100))

cat("\nGenerated:\n")
cat("  output/simulation_comparison.csv\n")
cat("  output/simulation_table.tex\n")
cat("  Figures/error_curves_comparison_*.png\n")
cat("  Figures/decision_heatmap_comparison.png\n")
cat("  Figures/ev_vs_kstar_comparison.png\n")
