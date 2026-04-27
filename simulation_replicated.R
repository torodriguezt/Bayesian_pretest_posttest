# simulation_replicated.R
# Replicated simulation study: rejection rates as a function of (n, d)
# For each (n, d) combination, generates R replicates of data and reports:
#   - Rejection rate (% of times H₀ is rejected)
#   - Mean ev_obs, Mean k*
#   - Comparison between KL-optimal and Informative (α=10) priors

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

# Simulation grid
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

R_reps <- 30      # Replicates per scenario (tradeoff: precision vs time)
M_sim  <- 300     # E-value simulations per analysis (reduced for speed)
set.seed(42)

# Time estimation:
#   20 scenarios × R reps × 2 priors = 1200 analyses (with R=30)
#   Each analysis: ~3-8 seconds (depends on n)
#   Total: ~1-3 hours
#
# To speed up further, reduce R_reps or M_sim
# To increase precision, increase R_reps (recommended: 50-100)

# ============================================================================
# 2. ANALYSIS FUNCTION (single replicate)
# ============================================================================

analyze_one <- function(n, x1, x2, a0, a1, a2, M = 500) {
  consts   <- bb_constants(n, n, x1, x2, a0, a1, a2)
  sup_H    <- find_sup_H(consts)$sup_H
  ev_obs   <- ev_quad(consts, sup_H)
  ev_H     <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M)
  ev_A     <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M)
  opt      <- find_kstar(ev_H, ev_A, 1, 1)

  list(
    ev_obs = ev_obs,
    k_star = opt$k_star,
    alpha  = opt$alpha,
    beta   = opt$beta,
    reject = (ev_obs <= opt$k_star)
  )
}

# ============================================================================
# 3. MAIN SIMULATION LOOP (with replications)
# ============================================================================

cat("\n========== REPLICATED SIMULATION STUDY ==========\n")
cat(sprintf("R = %d replicates per scenario\n", R_reps))
cat(sprintf("M = %d e-value simulations per analysis\n", M_sim))
cat(sprintf("Total scenarios: %d\n", nrow(scenarios)))
cat(sprintf("Total analyses: %d × 2 priors = %d\n\n",
            nrow(scenarios) * R_reps, nrow(scenarios) * R_reps * 2))

# Store all individual results
all_results <- data.frame()

t_start <- Sys.time()

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]

  cat(sprintf("[%2d/%d] n=%3d, %s ", i, nrow(scenarios), s$n, s$d_label))

  # Storage for this scenario
  ev_obs_kl  <- numeric(R_reps); k_kl  <- numeric(R_reps); rej_kl  <- logical(R_reps)
  ev_obs_inf <- numeric(R_reps); k_inf <- numeric(R_reps); rej_inf <- logical(R_reps)

  for (r in seq_len(R_reps)) {
    # Generate data ONCE per replicate (same for both priors)
    x1 <- rbinom(1, s$n, s$theta1)
    x2 <- rbinom(1, s$n, s$theta2)

    # KL-optimal prior
    res_kl <- analyze_one(s$n, x1, x2, a0_kl, a1_kl, a2_kl, M=M_sim)
    ev_obs_kl[r] <- res_kl$ev_obs
    k_kl[r]      <- res_kl$k_star
    rej_kl[r]    <- res_kl$reject

    # Informative prior
    res_inf <- analyze_one(s$n, x1, x2, a0_inf, a1_inf, a2_inf, M=M_sim)
    ev_obs_inf[r] <- res_inf$ev_obs
    k_inf[r]      <- res_inf$k_star
    rej_inf[r]    <- res_inf$reject

    if (r %% 10 == 0) cat(".")
  }

  # Aggregate this scenario
  all_results <- rbind(all_results,
    data.frame(
      n = s$n, d = s$d, d_label = s$d_label,
      prior = "KL-optimal",
      reject_rate = mean(rej_kl),
      mean_ev = mean(ev_obs_kl),
      mean_k_star = mean(k_kl),
      sd_ev = sd(ev_obs_kl),
      sd_k_star = sd(k_kl)
    ),
    data.frame(
      n = s$n, d = s$d, d_label = s$d_label,
      prior = "Informative (α=10)",
      reject_rate = mean(rej_inf),
      mean_ev = mean(ev_obs_inf),
      mean_k_star = mean(k_inf),
      sd_ev = sd(ev_obs_inf),
      sd_k_star = sd(k_inf)
    )
  )

  cat(sprintf(" KL=%.0f%% Inf=%.0f%%\n", mean(rej_kl)*100, mean(rej_inf)*100))
}

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\nTotal time: %.1f minutes\n", elapsed))

# ============================================================================
# 4. SUMMARY TABLE
# ============================================================================

cat("\n========== SUMMARY: REJECTION RATES ==========\n")

# Pivot for side-by-side display
summary_wide <- all_results %>%
  select(n, d, d_label, prior, reject_rate, mean_ev, mean_k_star) %>%
  pivot_wider(names_from = prior,
              values_from = c(reject_rate, mean_ev, mean_k_star))

print(summary_wide, row.names = FALSE, digits = 3)

write.csv(all_results, "output/simulation_replicated.csv", row.names = FALSE)

# ============================================================================
# 5. LATEX TABLE
# ============================================================================

tex_table <- c(
  "\\begin{table}[!h]",
  "\\centering",
  "\\small",
  paste0("\\caption{Replicated simulation study: rejection rates of FBST under ",
         "different sample sizes ($n$) and effect sizes ($\\Delta = \\theta_2 - \\theta_1$). ",
         sprintf("$R = %d$ replicates per scenario, $M = %d$ e-value simulations.}", R_reps, M_sim)),
  "\\label{tab:replicated_simulation}",
  "\\begin{tabular}{cc|ccc|ccc}",
  "\\toprule",
  "\\multicolumn{2}{c}{} & \\multicolumn{3}{c}{\\textbf{KL-optimal}} & \\multicolumn{3}{c}{\\textbf{Informative ($\\alpha=10$)}} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "$n$ & $\\Delta$ & Rej. rate & $\\overline{ev}$ & $\\overline{k^*}$ & Rej. rate & $\\overline{ev}$ & $\\overline{k^*}$ \\\\",
  "\\midrule"
)

for (i in seq(1, nrow(all_results), by = 2)) {
  r_kl  <- all_results[i, ]
  r_inf <- all_results[i+1, ]

  tex_table <- c(tex_table,
    sprintf("%d & %.2f & %.0f\\%% & %.3f & %.3f & %.0f\\%% & %.3f & %.3f \\\\",
            r_kl$n, r_kl$d,
            r_kl$reject_rate * 100, r_kl$mean_ev, r_kl$mean_k_star,
            r_inf$reject_rate * 100, r_inf$mean_ev, r_inf$mean_k_star)
  )
}

tex_table <- c(tex_table, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_table, "output/simulation_replicated_table.tex")
cat("\n→ output/simulation_replicated_table.tex\n")

# ============================================================================
# 6. PLOT 1: Power curve (rejection rate vs n, faceted by effect size)
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")

p_power <- ggplot(all_results, aes(x = n, y = reject_rate * 100,
                                    color = prior, linetype = prior)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = c("KL-optimal" = "steelblue",
                                "Informative (α=10)" = "tomato")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Sample size (n)", y = "Rejection rate (%)",
       title = "FBST power: rejection rate vs sample size",
       subtitle = sprintf("R=%d replicates per scenario", R_reps),
       color = "Prior", linetype = "Prior") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))

ggsave("Figures/power_curves.png", p_power, width = 10, height = 7, dpi = 150)
cat("  → Figures/power_curves.png\n")

# ============================================================================
# 7. PLOT 2: Heatmap of rejection rates
# ============================================================================

p_heatmap <- ggplot(all_results, aes(x = factor(n), y = d_label,
                                      fill = reject_rate * 100)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", reject_rate * 100)),
            color = "white", fontface = "bold", size = 4) +
  facet_wrap(~prior) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 50, limits = c(0, 100),
                       name = "Rejection\nrate (%)") +
  labs(x = "Sample size (n)", y = "Effect size",
       title = "Rejection rates: KL-optimal vs Informative prior") +
  theme_bw() +
  theme(strip.text = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 11))

ggsave("Figures/rejection_heatmap.png", p_heatmap, width = 12, height = 6, dpi = 150)
cat("  → Figures/rejection_heatmap.png\n")

# ============================================================================
# 8. PLOT 3: Mean ev_obs and k* per scenario
# ============================================================================

results_long <- all_results %>%
  select(n, d_label, prior, mean_ev, mean_k_star) %>%
  pivot_longer(c("mean_ev", "mean_k_star"),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("mean_ev", "mean_k_star"),
                         labels = c("Mean ev_obs", "Mean k*")))

p_means <- ggplot(results_long, aes(x = n, y = value, color = prior, linetype = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = c("KL-optimal" = "steelblue",
                                "Informative (α=10)" = "tomato")) +
  labs(x = "Sample size (n)", y = "Value",
       title = "Mean observed e-value and adaptive cutoff k*",
       color = "Prior", linetype = "Metric") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))

ggsave("Figures/ev_kstar_means.png", p_means, width = 10, height = 7, dpi = 150)
cat("  → Figures/ev_kstar_means.png\n")

# ============================================================================
# 9. PLOT 4: Decision agreement between priors
# ============================================================================

agreement_data <- all_results %>%
  select(n, d, d_label, prior, reject_rate) %>%
  pivot_wider(names_from = prior, values_from = reject_rate) %>%
  mutate(
    diff = `KL-optimal` - `Informative (α=10)`,
    abs_diff = abs(diff)
  )

p_agreement <- ggplot(agreement_data, aes(x = factor(n), y = d_label,
                                           fill = diff * 100)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = sprintf("%+.0f%%", diff * 100)),
            color = "black", fontface = "bold", size = 4) +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "steelblue",
                       midpoint = 0, name = "Difference\n(KL - Inf, %)") +
  labs(x = "Sample size (n)", y = "Effect size",
       title = "Difference in rejection rates: KL-optimal minus Informative",
       subtitle = "Positive (blue): KL rejects more | Negative (red): Inf rejects more") +
  theme_bw() +
  theme(axis.text = element_text(size = 11))

ggsave("Figures/decision_agreement.png", p_agreement, width = 9, height = 5, dpi = 150)
cat("  → Figures/decision_agreement.png\n")

# ============================================================================
# 10. FINAL SUMMARY
# ============================================================================

cat("\n========== FINAL SUMMARY ==========\n")
cat(sprintf("Total time: %.1f minutes\n", elapsed))
cat(sprintf("R = %d replicates × M = %d e-values\n", R_reps, M_sim))

cat("\nGenerated files:\n")
cat("  - output/simulation_replicated.csv         (raw data)\n")
cat("  - output/simulation_replicated_table.tex   (LaTeX table)\n")
cat("  - Figures/power_curves.png                 (rejection rate vs n)\n")
cat("  - Figures/rejection_heatmap.png            (rate heatmap)\n")
cat("  - Figures/ev_kstar_means.png               (mean ev and k*)\n")
cat("  - Figures/decision_agreement.png           (KL vs Inf difference)\n")

# Type I error and power summary
cat("\n--- Type I error (when H₀ is true, d=0) ---\n")
type1 <- all_results %>% filter(d == 0)
print(type1[, c("n", "prior", "reject_rate")], row.names = FALSE, digits = 3)

cat("\n--- Power (when H₁ is true) ---\n")
power_summary <- all_results %>%
  filter(d > 0) %>%
  group_by(d_label, prior) %>%
  summarise(mean_power = mean(reject_rate), .groups = "drop")
print(power_summary, row.names = FALSE, digits = 3)

cat("\n")
