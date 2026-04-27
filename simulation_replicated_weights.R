# simulation_replicated_weights.R
# Replicated simulation study with multiple (a, b) weight configurations
# For each scenario, computes k* under different penalties for Type I vs Type II errors
#
# KEY OPTIMIZATION: ev_H and ev_A are simulated ONCE per (n, d, replicate, prior),
# then find_kstar() is called multiple times with different (a, b) - very cheap.

library(Rcpp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")

dir.create("output", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# ============================================================================
# 1. PARAMETERS
# ============================================================================

# Hyperparameters
a0_kl  <- 0.8373879;  a1_kl  <- 0.8410984;  a2_kl  <- 0.8053298
a0_inf <- 10;          a1_inf <- 10;         a2_inf <- 10

# Weight configurations for k*
# Each row: (label, a, b)
# a = penalty for Type I error (rejecting H when H is true)
# b = penalty for Type II error (failing to reject H when A is true)
weight_configs <- data.frame(
  label = c("a=1, b=1 (balanced)",
            "a=5, b=1 (moderate)",
            "a=20, b=1 (strict)"),
  a_w   = c(1,  5, 20),
  b_w   = c(1,  1, 1),
  stringsAsFactors = FALSE
)

# Scenarios
n_grid <- c(30, 50, 75, 100, 150)
scenarios <- expand.grid(
  n = n_grid,
  scenario_type = c("null", "small", "medium", "large")
) %>%
  mutate(
    theta1 = 0.5,
    theta2 = case_when(
      scenario_type == "null"   ~ 0.50,
      scenario_type == "small"  ~ 0.52,
      scenario_type == "medium" ~ 0.60,
      scenario_type == "large"  ~ 0.80
    ),
    d = theta2 - theta1,
    d_label = case_when(
      scenario_type == "null"   ~ "Δ=0.00 (H₀ true)",
      scenario_type == "small"  ~ "Δ=0.02 (small)",
      scenario_type == "medium" ~ "Δ=0.10 (medium)",
      scenario_type == "large"  ~ "Δ=0.30 (large)"
    )
  )

R_reps <- 30
M_sim  <- 300
set.seed(42)

# ============================================================================
# 2. ANALYSIS FUNCTION
# ============================================================================
# Returns ev_obs and the simulated ev_H, ev_A vectors
# (k* is computed later for each (a, b) configuration)

analyze_one <- function(n, x1, x2, a0, a1, a2, M = 300) {
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  sup_H  <- find_sup_H(consts)$sup_H
  ev_obs <- ev_quad(consts, sup_H)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M)

  list(ev_obs = ev_obs, ev_H = ev_H, ev_A = ev_A)
}

# ============================================================================
# 3. MAIN SIMULATION LOOP
# ============================================================================

cat("\n========== REPLICATED SIMULATION (multiple weights) ==========\n")
cat(sprintf("R = %d replicates per scenario\n", R_reps))
cat(sprintf("M = %d e-value simulations per analysis\n", M_sim))
cat(sprintf("Weight configurations: %d\n", nrow(weight_configs)))
cat(sprintf("Total scenarios: %d\n", nrow(scenarios)))
cat(sprintf("Total analyses: %d × 2 priors = %d\n",
            nrow(scenarios) * R_reps, nrow(scenarios) * R_reps * 2))
cat("(k* will be computed for each weight config without re-simulating)\n\n")

# Storage for all results
all_results <- data.frame()

t_start <- Sys.time()

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  cat(sprintf("[%2d/%d] n=%3d, %s ", i, nrow(scenarios), s$n, s$d_label))

  # Storage for this scenario: matrix [R x weights] of rejections
  rej_kl  <- array(FALSE, dim = c(R_reps, nrow(weight_configs)))
  rej_inf <- array(FALSE, dim = c(R_reps, nrow(weight_configs)))
  ev_kl_vec  <- numeric(R_reps); ev_inf_vec <- numeric(R_reps)

  # Per-weight storage of k* values
  k_kl_mat  <- array(NA, dim = c(R_reps, nrow(weight_configs)))
  k_inf_mat <- array(NA, dim = c(R_reps, nrow(weight_configs)))

  for (r in seq_len(R_reps)) {
    # Generate data ONCE per replicate
    x1 <- rbinom(1, s$n, s$theta1)
    x2 <- rbinom(1, s$n, s$theta2)

    # Simulate ev distributions ONCE per prior
    res_kl  <- analyze_one(s$n, x1, x2, a0_kl,  a1_kl,  a2_kl,  M = M_sim)
    res_inf <- analyze_one(s$n, x1, x2, a0_inf, a1_inf, a2_inf, M = M_sim)

    ev_kl_vec[r]  <- res_kl$ev_obs
    ev_inf_vec[r] <- res_inf$ev_obs

    # For each weight config, compute k* (cheap operation)
    for (w in seq_len(nrow(weight_configs))) {
      a_w <- weight_configs$a_w[w]
      b_w <- weight_configs$b_w[w]

      opt_kl  <- find_kstar(res_kl$ev_H,  res_kl$ev_A,  a_w, b_w)
      opt_inf <- find_kstar(res_inf$ev_H, res_inf$ev_A, a_w, b_w)

      k_kl_mat[r, w]  <- opt_kl$k_star
      k_inf_mat[r, w] <- opt_inf$k_star
      rej_kl[r, w]    <- (res_kl$ev_obs  <= opt_kl$k_star)
      rej_inf[r, w]   <- (res_inf$ev_obs <= opt_inf$k_star)
    }

    if (r %% 10 == 0) cat(".")
  }

  # Aggregate per weight configuration
  for (w in seq_len(nrow(weight_configs))) {
    all_results <- rbind(all_results,
      data.frame(
        n = s$n, d = s$d, d_label = s$d_label,
        weight_label = weight_configs$label[w],
        a_w = weight_configs$a_w[w], b_w = weight_configs$b_w[w],
        prior = "KL-optimal",
        reject_rate = mean(rej_kl[, w]),
        mean_ev = mean(ev_kl_vec),
        mean_k_star = mean(k_kl_mat[, w]),
        sd_k_star   = sd(k_kl_mat[, w])
      ),
      data.frame(
        n = s$n, d = s$d, d_label = s$d_label,
        weight_label = weight_configs$label[w],
        a_w = weight_configs$a_w[w], b_w = weight_configs$b_w[w],
        prior = "Informative (α=10)",
        reject_rate = mean(rej_inf[, w]),
        mean_ev = mean(ev_inf_vec),
        mean_k_star = mean(k_inf_mat[, w]),
        sd_k_star   = sd(k_inf_mat[, w])
      )
    )
  }

  cat(sprintf(" KL[1,5,20]=[%.0f%%, %.0f%%, %.0f%%]\n",
              mean(rej_kl[, 1]) * 100,
              mean(rej_kl[, 2]) * 100,
              mean(rej_kl[, 3]) * 100))
}

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\nTotal time: %.1f minutes\n", elapsed))

# ============================================================================
# 4. SUMMARY: TYPE I ERROR
# ============================================================================

cat("\n========== TYPE I ERROR (when H₀ is true, d=0) ==========\n")

type1 <- all_results %>%
  filter(d == 0) %>%
  group_by(weight_label, prior) %>%
  summarise(mean_type1 = mean(reject_rate), .groups = "drop") %>%
  arrange(weight_label, prior)

print(as.data.frame(type1), row.names = FALSE, digits = 3)

cat("\n========== POWER (when H₁ is true) ==========\n")

power_summary <- all_results %>%
  filter(d > 0) %>%
  group_by(d_label, weight_label, prior) %>%
  summarise(mean_power = mean(reject_rate), .groups = "drop") %>%
  arrange(d_label, weight_label, prior)

print(as.data.frame(power_summary), row.names = FALSE, digits = 3)

# ============================================================================
# 5. SAVE DATA
# ============================================================================

write.csv(all_results, "output/simulation_replicated_weights.csv", row.names = FALSE)

# ============================================================================
# 6. PLOT 1: Rejection rate vs n, by weight configuration
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")

p_power_weights <- ggplot(all_results,
                          aes(x = n, y = reject_rate * 100,
                              color = weight_label, linetype = prior)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = c("a=1, b=1 (balanced)" = "steelblue",
                                 "a=5, b=1 (moderate)" = "darkorange",
                                 "a=20, b=1 (strict)"  = "tomato")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Sample size (n)", y = "Rejection rate (%)",
       title = "FBST: rejection rate vs sample size",
       subtitle = "Different weight configurations and priors",
       color = "Weight (a, b)", linetype = "Prior") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        strip.text = element_text(size = 11, face = "bold"))

ggsave("Figures/power_curves_weights.png", p_power_weights,
       width = 11, height = 8, dpi = 150)
cat("  → Figures/power_curves_weights.png\n")

# ============================================================================
# 7. PLOT 2: Type I error vs n (only d=0)
# ============================================================================

p_type1 <- all_results %>%
  filter(d == 0) %>%
  ggplot(aes(x = n, y = reject_rate * 100,
             color = weight_label, linetype = prior)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 5, linetype = "dotted", color = "grey40") +
  annotate("text", x = 30, y = 7, label = "α = 5% (traditional)",
           color = "grey40", size = 3) +
  scale_color_manual(values = c("a=1, b=1 (balanced)" = "steelblue",
                                 "a=5, b=1 (moderate)" = "darkorange",
                                 "a=20, b=1 (strict)"  = "tomato")) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(x = "Sample size (n)", y = "Type I error rate (%)",
       title = "Type I error rate when H₀ is true (Δ = 0)",
       subtitle = "How much does the prior + weight choice over-reject?",
       color = "Weight (a, b)", linetype = "Prior") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical")

ggsave("Figures/type1_error_weights.png", p_type1,
       width = 10, height = 6, dpi = 150)
cat("  → Figures/type1_error_weights.png\n")

# ============================================================================
# 8. PLOT 3: Heatmap of rejection rates per weight config
# ============================================================================

p_heatmap_weights <- ggplot(all_results,
                            aes(x = factor(n), y = d_label,
                                fill = reject_rate * 100)) +
  geom_tile(color = "black", size = 0.4) +
  geom_text(aes(label = sprintf("%.0f", reject_rate * 100)),
            color = "white", fontface = "bold", size = 3.5) +
  facet_grid(prior ~ weight_label) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 50, limits = c(0, 100),
                       name = "Rej.\nrate (%)") +
  labs(x = "Sample size (n)", y = "Effect size",
       title = "Rejection rates: priors × weight configurations") +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"))

ggsave("Figures/rejection_heatmap_weights.png", p_heatmap_weights,
       width = 14, height = 7, dpi = 150)
cat("  → Figures/rejection_heatmap_weights.png\n")

# ============================================================================
# 9. PLOT 4: k* values vs n for each weight (KL-optimal only)
# ============================================================================

p_kstar <- all_results %>%
  filter(prior == "KL-optimal") %>%
  ggplot(aes(x = n, y = mean_k_star, color = weight_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = pmax(0, mean_k_star - sd_k_star),
                  ymax = pmin(1, mean_k_star + sd_k_star),
                  fill = weight_label), alpha = 0.15, color = NA) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = c("a=1, b=1 (balanced)" = "steelblue",
                                 "a=5, b=1 (moderate)" = "darkorange",
                                 "a=20, b=1 (strict)"  = "tomato")) +
  scale_fill_manual(values = c("a=1, b=1 (balanced)" = "steelblue",
                                "a=5, b=1 (moderate)" = "darkorange",
                                "a=20, b=1 (strict)"  = "tomato")) +
  labs(x = "Sample size (n)", y = "Mean k* ± SD",
       title = "Adaptive cutoff k* by weight configuration (KL-optimal prior)",
       color = "Weight (a, b)", fill = "Weight (a, b)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))

ggsave("Figures/kstar_by_weights.png", p_kstar, width = 11, height = 7, dpi = 150)
cat("  → Figures/kstar_by_weights.png\n")

# ============================================================================
# 10. LATEX TABLE
# ============================================================================

# Type I error table
tex_t1 <- c(
  "\\begin{table}[!h]", "\\centering", "\\small",
  paste0("\\caption{Type I error rates (\\%) under different weight configurations ",
         "and priors. Generated from $R = ", R_reps, "$ replicates per scenario.",
         " H$_0$ is true ($\\Delta = 0$).}"),
  "\\label{tab:type1_weights}",
  "\\begin{tabular}{cc|ccc|ccc}", "\\toprule",
  "\\multicolumn{2}{c}{} & \\multicolumn{3}{c}{\\textbf{KL-optimal}} & \\multicolumn{3}{c}{\\textbf{Informative}} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "$n$ & $\\Delta$ & $a{=}1$ & $a{=}5$ & $a{=}20$ & $a{=}1$ & $a{=}5$ & $a{=}20$ \\\\",
  "\\midrule"
)

# Rebuild table with one row per (n, d) and columns per weight × prior
type1_table <- all_results %>%
  filter(d == 0) %>%
  select(n, d, weight_label, prior, reject_rate) %>%
  pivot_wider(names_from = c(prior, weight_label),
              values_from = reject_rate)

for (i in seq_len(nrow(type1_table))) {
  r <- type1_table[i, ]
  tex_t1 <- c(tex_t1, sprintf(
    "%d & %.2f & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% \\\\",
    r$n, r$d,
    r[["KL-optimal_a=1, b=1 (balanced)"]] * 100,
    r[["KL-optimal_a=5, b=1 (moderate)"]] * 100,
    r[["KL-optimal_a=20, b=1 (strict)"]] * 100,
    r[["Informative (α=10)_a=1, b=1 (balanced)"]] * 100,
    r[["Informative (α=10)_a=5, b=1 (moderate)"]] * 100,
    r[["Informative (α=10)_a=20, b=1 (strict)"]] * 100
  ))
}
tex_t1 <- c(tex_t1, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_t1, "output/type1_weights_table.tex")
cat("  → output/type1_weights_table.tex\n")

# Power table for medium effect
power_table <- all_results %>%
  filter(d == 0.10) %>%
  select(n, d, weight_label, prior, reject_rate) %>%
  pivot_wider(names_from = c(prior, weight_label),
              values_from = reject_rate)

tex_pow <- c(
  "\\begin{table}[!h]", "\\centering", "\\small",
  paste0("\\caption{Power (\\%) under different weight configurations and priors. ",
         "Generated from $R = ", R_reps, "$ replicates per scenario. ",
         "$\\Delta = 0.10$ (medium effect).}"),
  "\\label{tab:power_weights}",
  "\\begin{tabular}{cc|ccc|ccc}", "\\toprule",
  "\\multicolumn{2}{c}{} & \\multicolumn{3}{c}{\\textbf{KL-optimal}} & \\multicolumn{3}{c}{\\textbf{Informative}} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "$n$ & $\\Delta$ & $a{=}1$ & $a{=}5$ & $a{=}20$ & $a{=}1$ & $a{=}5$ & $a{=}20$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(power_table))) {
  r <- power_table[i, ]
  tex_pow <- c(tex_pow, sprintf(
    "%d & %.2f & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% \\\\",
    r$n, r$d,
    r[["KL-optimal_a=1, b=1 (balanced)"]] * 100,
    r[["KL-optimal_a=5, b=1 (moderate)"]] * 100,
    r[["KL-optimal_a=20, b=1 (strict)"]] * 100,
    r[["Informative (α=10)_a=1, b=1 (balanced)"]] * 100,
    r[["Informative (α=10)_a=5, b=1 (moderate)"]] * 100,
    r[["Informative (α=10)_a=20, b=1 (strict)"]] * 100
  ))
}
tex_pow <- c(tex_pow, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_pow, "output/power_weights_table.tex")
cat("  → output/power_weights_table.tex\n")

# ============================================================================
# 11. FINAL SUMMARY
# ============================================================================

cat("\n========== FINAL SUMMARY ==========\n")
cat(sprintf("Total time: %.1f minutes\n", elapsed))
cat(sprintf("R = %d replicates × M = %d e-values × %d weight configs\n",
            R_reps, M_sim, nrow(weight_configs)))

cat("\nGenerated files:\n")
cat("  - output/simulation_replicated_weights.csv  (raw data)\n")
cat("  - output/type1_weights_table.tex             (LaTeX Type I table)\n")
cat("  - output/power_weights_table.tex             (LaTeX Power table)\n")
cat("  - Figures/power_curves_weights.png           (rejection vs n)\n")
cat("  - Figures/type1_error_weights.png            (Type I error focus)\n")
cat("  - Figures/rejection_heatmap_weights.png      (rate heatmap)\n")
cat("  - Figures/kstar_by_weights.png               (k* by weight)\n")

cat("\nKey insight: with a=20, b=1, Type I error should drop to ~5%\n")
cat("  (matching traditional frequentist convention)\n\n")
