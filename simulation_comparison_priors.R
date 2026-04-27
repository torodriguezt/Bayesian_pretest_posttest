# simulation_comparison_priors.R
# Estudio de simulación: efecto de tamaño de muestra (n) y tamaño de efecto (d)
# Compara: Priori no informativa (KL-óptima) vs Priori informativa (α=10)
#
# Salida:
#   - Tabla resumen (n, d, P(H), ev_obs, k*, decisión para ambos priors)
#   - Gráficos: curvas de error, superficies posteriores, heatmap de decisiones

library(ALA)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(tidyr)
library(gridExtra)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")

dir.create("output", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# ============================================================================
# 1. PARÁMETROS DE SIMULACIÓN
# ============================================================================

# Hiperparámetros
a0_kl <- 0.8373879
a1_kl <- 0.8410984
a2_kl <- 0.8053298

a0_inf <- 10
a1_inf <- 10
a2_inf <- 10

# Grilla de simulación REALISTA (con casos de NO rechazo)
# Estructura: (theta1, theta2) pairs that represent different scenarios
scenarios <- expand.grid(
  n = c(30, 50, 75, 100, 150),
  scenario_type = c("null", "small", "medium", "large")
) %>%
  mutate(
    # Genera pares (theta1, theta2) según tipo
    theta1 = case_when(
      scenario_type == "null"   ~ 0.5,
      scenario_type == "small"  ~ 0.5,
      scenario_type == "medium" ~ 0.5,
      scenario_type == "large"  ~ 0.5
    ),
    theta2 = case_when(
      scenario_type == "null"   ~ 0.5,      # H₀ verdadera: efecto = 0
      scenario_type == "small"  ~ 0.52,     # efecto pequeño: d = 0.02
      scenario_type == "medium" ~ 0.60,     # efecto mediano: d = 0.10
      scenario_type == "large"  ~ 0.80      # efecto grande: d = 0.30
    ),
    d = theta2 - theta1,
    d_label = case_when(
      scenario_type == "null"   ~ "0.00 (H₀ true)",
      scenario_type == "small"  ~ "0.02 (small effect)",
      scenario_type == "medium" ~ "0.10 (medium effect)",
      scenario_type == "large"  ~ "0.30 (large effect)"
    )
  )

M_sim    <- 1000                         # simulaciones para curvas de error
set.seed(42)

# ============================================================================
# 2. FUNCIÓN AUXILIAR: analizar un escenario (n, d, prior)
# ============================================================================

analyze_scenario <- function(n, x1, x2, a0, a1, a2, M = 1000, label = "scenario") {
  # Receives x1, x2 as PARAMETERS (no longer generates them internally)

  # Calculate constants and supremum
  consts  <- bb_constants(n, n, x1, x2, a0, a1, a2)
  sup_info <- find_sup_H(consts)
  sup_H   <- sup_info$sup_H

  # Observed e-value
  ev_obs <- ev_quad(consts, sup_H)

  # Simulate e-values under H and A
  ev_H <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M)
  ev_A <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M)

  # Calculate k* and errors
  opt <- find_kstar(ev_H, ev_A, 1, 1)
  k_star <- opt$k_star
  alpha  <- opt$alpha
  beta   <- opt$beta

  # Decision
  decision <- ifelse(ev_obs <= k_star, "Reject H₀", "Do not reject")

  # P(H | data) = P(θ1 ≤ θ2 | data)
  posterior_samples <- sample_posterior(5000, n, n, x1, x2, a0, a1, a2)
  prob_H <- mean(posterior_samples[, 1] <= posterior_samples[, 2])

  # Return
  list(
    n = n, label = label,
    x1 = x1, x2 = x2,
    ev_obs = ev_obs,
    k_star = k_star,
    alpha = alpha,
    beta = beta,
    decision = decision,
    prob_H = prob_H,
    ev_H = ev_H,
    ev_A = ev_A,
    consts = consts,
    sup_H = sup_H
  )
}

# ============================================================================
# 3. TABLA RESUMEN PRINCIPAL
# ============================================================================

cat("\n========== SIMULATION: PRIOR COMPARISON ==========\n\n")

results <- data.frame()

for (i in 1:nrow(scenarios)) {
  row_scenario <- scenarios[i, ]
  n_val   <- row_scenario$n
  theta1  <- row_scenario$theta1
  theta2  <- row_scenario$theta2
  d_val   <- row_scenario$d
  d_label <- row_scenario$d_label

  cat(sprintf("n=%d, θ1=%.2f, θ2=%.2f (%s)... ", n_val, theta1, theta2, d_label))

  # Generate data ONCE (use SAME data for both priors)
  x1 <- rbinom(1, n_val, theta1)
  x2 <- rbinom(1, n_val, theta2)

  # Analysis with KL-optimal prior (same x1, x2)
  res_kl <- analyze_scenario(n_val, x1, x2, a0_kl, a1_kl, a2_kl, M=M_sim,
                             label = "KL-optimal")

  # Analysis with informative prior (SAME x1, x2)
  res_inf <- analyze_scenario(n_val, x1, x2, a0_inf, a1_inf, a2_inf, M=M_sim,
                              label = "Informative (α=10)")

  # Add to table (SAME DATA for both priors)
  results <- rbind(results,
    data.frame(
      n = n_val,
      theta1 = theta1,
      theta2 = theta2,
      d = d_val,
      d_label = d_label,
      x1 = x1, x2 = x2,
      prior = "KL-optimal",
      prob_H = round(res_kl$prob_H, 3),
      ev_obs = round(res_kl$ev_obs, 4),
      k_star = round(res_kl$k_star, 4),
      alpha = round(res_kl$alpha, 4),
      beta = round(res_kl$beta, 4),
      decision = res_kl$decision
    ),
    data.frame(
      n = n_val,
      theta1 = theta1,
      theta2 = theta2,
      d = d_val,
      d_label = d_label,
      x1 = x1, x2 = x2,
      prior = "Informative (α=10)",
      prob_H = round(res_inf$prob_H, 3),
      ev_obs = round(res_inf$ev_obs, 4),
      k_star = round(res_inf$k_star, 4),
      alpha = round(res_inf$alpha, 4),
      beta = round(res_inf$beta, 4),
      decision = res_inf$decision
    )
  )

  cat("OK\n")
}

cat("\n========== SUMMARY TABLE ==========\n")
print(results[, c("n", "d", "prior", "prob_H", "ev_obs", "k_star", "decision")])

# Guarda como CSV y LaTeX
write.csv(results, "output/simulation_comparison.csv", row.names = FALSE)

# ============================================================================
# 4. TABLA LATEX
# ============================================================================

# Pivotea para formato lado a lado (similar al ejemplo del usuario)
results_wide <- results %>%
  select(n, d, prior, prob_H, ev_obs, k_star, decision) %>%
  pivot_wider(
    names_from = prior,
    values_from = c(prob_H, ev_obs, k_star, decision)
  )

tex_table <- c(
  "\\begin{table}[!h]",
  "\\centering",
  "\\small",
  "\\caption{Simulation study: effect of sample size ($n$) and effect size ($\\Delta = \\theta_2 - \\theta_1$) on FBST decision. Comparison between KL-optimal prior and informative prior ($\\alpha_0=\\alpha_1=\\alpha_2=10$), with $M=1000$ simulated e-values.}",
  "\\label{tab:simulation}",
  "\\begin{tabular}{cc|cccc|cccc}",
  "\\toprule",
  "\\multicolumn{2}{c}{} & \\multicolumn{4}{c}{\\textbf{KL-optimal prior}} & \\multicolumn{4}{c}{\\textbf{Informative prior}} \\\\",
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
  "$n$ & $\\Delta$ & $P(H)$ & $ev$ & $k^*$ & Decision & $P(H)$ & $ev$ & $k^*$ & Decision \\\\",
  "\\midrule"
)

for (i in seq(1, nrow(results), by = 2)) {
  r_kl  <- results[i, ]
  r_inf <- results[i+1, ]

  tex_table <- c(tex_table,
    sprintf("%d & %.1f & %.3f & %.4f & %.4f & %s & %.3f & %.4f & %.4f & %s \\\\",
            r_kl$n, r_kl$d,
            r_kl$prob_H, r_kl$ev_obs, r_kl$k_star,
            ifelse(r_kl$decision == "Reject H₀", "Reject", "Do not reject"),
            r_inf$prob_H, r_inf$ev_obs, r_inf$k_star,
            ifelse(r_inf$decision == "Reject H₀", "Reject", "Do not reject"))
  )
}

tex_table <- c(tex_table, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_table, "output/simulation_table.tex")

cat("\n→ output/simulation_table.tex\n")

# ============================================================================
# 5. GRÁFICOS: CURVAS DE ERROR (comparación de priors)
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")

plot_error_curves_comparison <- function(n_val, d_val, res_kl, res_inf) {

  k_grid <- seq(0, 1, length.out = 401)

  # Calcula curvas para KL
  curves_kl <- error_curves(res_kl$ev_H, res_kl$ev_A, k_grid, 1, 1)
  curves_kl$prior <- "KL-optimal"

  # Calculate curves for informative
  curves_inf <- error_curves(res_inf$ev_H, res_inf$ev_A, k_grid, 1, 1)
  curves_inf$prior <- "Informative"

  curves_all <- rbind(curves_kl, curves_inf)
  long <- pivot_longer(curves_all, c("alpha", "beta", "sum"),
                       names_to = "metric", values_to = "value")
  long$metric <- factor(long$metric,
                        levels = c("alpha", "beta", "sum"),
                        labels = c("α (Type I error)", "β (Type II error)", "α+β"))

  # Plot
  p <- ggplot(long, aes(k, value, color = metric, linetype = metric)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~prior, labeller = label_parsed) +
    geom_vline(data = data.frame(prior = c("KL-optimal", "Informative"),
                                  k_star = c(res_kl$k_star, res_inf$k_star)),
               aes(xintercept = k_star), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("steelblue", "tomato", "black")) +
    scale_linetype_manual(values = c("solid", "solid", "dashed")) +
    labs(x = "k", y = "Average error",
         title = sprintf("Error curves: n=%d, Δ=%.2f", n_val, d_val),
         color = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.position = "top")

  return(p)
}

# Selecciona algunos escenarios representativos
scenarios_to_plot <- expand.grid(
  n = c(50, 100, 150),
  d = c(0.2, 0.8)
)

for (i in 1:nrow(scenarios_to_plot)) {
  n_val <- scenarios_to_plot[i, "n"]
  d_val <- scenarios_to_plot[i, "d"]

  cat(sprintf("  Gráfico curvas: n=%d, d=%.1f\n", n_val, d_val))

  # Recupera resultados para este escenario
  idx_kl  <- which(results$n == n_val & results$d == d_val &
                   results$prior == "KL-óptima")
  idx_inf <- which(results$n == n_val & results$d == d_val &
                   results$prior == "Informativa (α=10)")

  if (length(idx_kl) > 0 && length(idx_inf) > 0) {
    res_kl <- analyze_scenario(n_val, d_val, a0_kl, a1_kl, a2_kl, M=M_sim)
    res_inf <- analyze_scenario(n_val, d_val, a0_inf, a1_inf, a2_inf, M=M_sim)

    p <- plot_error_curves_comparison(n_val, d_val, res_kl, res_inf)
    ggsave(sprintf("Figures/error_curves_comparison_n%d_d%.1f.png", n_val, d_val*10),
           p, width = 10, height = 4, dpi = 150)
  }
}

# ============================================================================
# 6. HEATMAP: DECISIÓN vs (n, d)
# ============================================================================

cat("  Decision heatmap plot\n")

results_decision <- results %>%
  mutate(decision_num = ifelse(decision == "Reject H₀", 1, 0)) %>%
  select(n, d, prior, decision_num) %>%
  pivot_wider(names_from = prior, values_from = decision_num)

# Heatmap KL
p_kl <- ggplot(results %>% filter(prior == "KL-optimal"),
               aes(x = factor(d), y = factor(n), fill = decision)) +
  geom_tile(color = "black", size = 0.8) +
  scale_fill_manual(values = c("Reject H₀" = "tomato", "Do not reject" = "steelblue")) +
  labs(x = "Effect size (Δ)", y = "Sample size (n)",
       title = "KL-optimal prior",
       fill = "Decision") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

# Heatmap informative
p_inf <- ggplot(results %>% filter(prior == "Informative (α=10)"),
                aes(x = factor(d), y = factor(n), fill = decision)) +
  geom_tile(color = "black", size = 0.8) +
  scale_fill_manual(values = c("Reject H₀" = "tomato", "Do not reject" = "steelblue")) +
  labs(x = "Effect size (Δ)", y = "Sample size (n)",
       title = "Informative prior (α=10)",
       fill = "Decision") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

p_combined <- gridExtra::grid.arrange(p_kl, p_inf, ncol = 2)
ggsave("Figures/decision_heatmap_comparison.png", p_combined,
       width = 10, height = 6, dpi = 150)

# ============================================================================
# 7. GRÁFICO: ev_obs vs k* (scatter)
# ============================================================================

cat("  Scatter plot: ev_obs vs k*\n")

p_scatter <- ggplot(results, aes(x = ev_obs, y = k_star, color = prior, shape = decision)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_grid(~prior) +
  scale_color_manual(values = c("KL-optimal" = "steelblue",
                                "Informative (α=10)" = "tomato")) +
  scale_shape_manual(values = c("Reject H₀" = 17, "Do not reject" = 16)) +
  labs(x = "Observed e-value", y = "k* (adaptive cutoff)",
       title = "Relationship between ev_obs and k*: Prior comparison",
       color = "Prior", shape = "Decision") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("Figures/ev_vs_kstar_comparison.png", p_scatter, width = 10, height = 5, dpi = 150)

# ============================================================================
# 8. FINAL SUMMARY
# ============================================================================

cat("\n========== FINAL SUMMARY ==========\n")
cat("Generated files:\n")
cat("  - output/simulation_comparison.csv       (data table)\n")
cat("  - output/simulation_table.tex            (LaTeX table)\n")
cat("  - Figures/error_curves_comparison_*.png  (error curves)\n")
cat("  - Figures/decision_heatmap_comparison.png (decision heatmap)\n")
cat("  - Figures/ev_vs_kstar_comparison.png     (scatter ev vs k*)\n")

cat("\nConclusions?\n")
cat("Comparing KL-optimal vs Informative:\n")

diffs <- results %>%
  select(n, d, prior, decision) %>%
  pivot_wider(names_from = prior, values_from = decision) %>%
  mutate(agreement = `KL-optimal` == `Informative (α=10)`)

cat(sprintf("  - Agreement in decisions: %.1f%%\n",
            mean(diffs$agreement) * 100))

cat("\n")
