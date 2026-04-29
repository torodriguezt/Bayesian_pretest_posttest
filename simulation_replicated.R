# simulation_replicated.R
# Replicated simulation study: rejection rates as a function of (n, d)
# Compares THREE priors (all KL-fitted, see priors_config.R):
#   - Non-informative (KL to U)
#   - Informative symmetric (mu=0.5, N=50)
#   - Conflict (mu=0.1, N=50)

library(Rcpp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

sourceCpp("BivBetaBinom.cpp")
source("priors_config.R")

dir.create("output",  showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# ============================================================================
# 1. PARAMETERS
# ============================================================================

priors <- list(
  "Non-informative"      = prior_NI,
  "Informative (mu=0.5)" = prior_INF,
  "Conflict (mu=0.1)"    = prior_CONF
)
prior_levels <- names(priors)

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
      scenario_type == "null"   ~ "Delta=0.00 (H true)",
      scenario_type == "small"  ~ "Delta=0.02 (small)",
      scenario_type == "medium" ~ "Delta=0.10 (medium)",
      scenario_type == "large"  ~ "Delta=0.30 (large)"
    )
  )

R_reps <- 30
M_sim  <- 300
set.seed(42)

# ============================================================================
# 2. ANALYSIS FUNCTION
# ============================================================================

analyze_one <- function(n, x1, x2, prior_alpha, M = 500) {
  a0 <- prior_alpha["a0"]; a1 <- prior_alpha["a1"]; a2 <- prior_alpha["a2"]
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  sup_H  <- find_sup_H(consts)$sup_H
  ev_obs <- ev_quad(consts, sup_H)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  list(ev_obs = ev_obs, k_star = opt$k_star,
       alpha = opt$alpha, beta = opt$beta,
       reject = (ev_obs <= opt$k_star))
}

# ============================================================================
# 3. MAIN SIMULATION LOOP
# ============================================================================

cat("\n========== REPLICATED SIMULATION STUDY ==========\n")
cat(sprintf("R = %d replicates per scenario\n", R_reps))
cat(sprintf("M = %d e-value simulations per analysis\n", M_sim))
cat(sprintf("Priors compared: %s\n", paste(prior_levels, collapse = ", ")))
cat(sprintf("Total scenarios: %d\n", nrow(scenarios)))
cat(sprintf("Total analyses: %d x %d priors = %d\n\n",
            nrow(scenarios) * R_reps, length(priors),
            nrow(scenarios) * R_reps * length(priors)))

all_results <- list()
t_start <- Sys.time()

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  cat(sprintf("[%2d/%d] n=%3d, %s ", i, nrow(scenarios), s$n, s$d_label))

  # Storage per prior
  store <- setNames(
    lapply(prior_levels, function(p) {
      list(ev = numeric(R_reps), k = numeric(R_reps), rej = logical(R_reps))
    }),
    prior_levels
  )

  for (r in seq_len(R_reps)) {
    # Same data for all priors within a replicate
    x1 <- rbinom(1, s$n, s$theta1)
    x2 <- rbinom(1, s$n, s$theta2)
    for (pname in prior_levels) {
      res <- analyze_one(s$n, x1, x2, priors[[pname]], M = M_sim)
      store[[pname]]$ev[r]  <- res$ev_obs
      store[[pname]]$k[r]   <- res$k_star
      store[[pname]]$rej[r] <- res$reject
    }
    if (r %% 10 == 0) cat(".")
  }

  for (pname in prior_levels) {
    all_results[[length(all_results) + 1]] <- data.frame(
      n           = s$n,
      d           = s$d,
      d_label     = s$d_label,
      prior       = pname,
      reject_rate = mean(store[[pname]]$rej),
      mean_ev     = mean(store[[pname]]$ev),
      mean_k_star = mean(store[[pname]]$k),
      sd_ev       = sd(store[[pname]]$ev),
      sd_k_star   = sd(store[[pname]]$k)
    )
  }

  cat(sprintf(" %s\n",
    paste(sapply(prior_levels, function(p)
      sprintf("%s=%.0f%%", substr(p, 1, 3), mean(store[[p]]$rej) * 100)),
      collapse = " ")))
}

all_results <- do.call(rbind, all_results)
all_results$prior <- factor(all_results$prior, levels = prior_levels)

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\nTotal time: %.1f minutes\n", elapsed))

# ============================================================================
# 4. SUMMARY TABLE
# ============================================================================

cat("\n========== SUMMARY: REJECTION RATES ==========\n")
summary_wide <- all_results %>%
  select(n, d, d_label, prior, reject_rate, mean_ev, mean_k_star) %>%
  pivot_wider(names_from  = prior,
              values_from = c(reject_rate, mean_ev, mean_k_star))
print(summary_wide, row.names = FALSE, digits = 3)
write.csv(all_results, "output/simulation_replicated.csv", row.names = FALSE)

# ============================================================================
# 5. LATEX TABLE
# ============================================================================

tex_table <- c(
  "\\begin{table}[!h]", "\\centering", "\\footnotesize",
  paste0("\\caption{Replicated simulation study: rejection rates for three ",
         "priors (KL-fitted) under different sample sizes ($n$) and effect ",
         sprintf("sizes ($\\Delta$). $R=%d$ replicates, $M=%d$.}", R_reps, M_sim)),
  "\\label{tab:replicated_simulation}",
  "\\begin{tabular}{cc|ccc|ccc|ccc}",
  "\\toprule",
  "\\multicolumn{2}{c}{} & \\multicolumn{3}{c}{\\textbf{Non-informative}} & \\multicolumn{3}{c}{\\textbf{Informative}} & \\multicolumn{3}{c}{\\textbf{Conflict}} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
  "$n$ & $\\Delta$ & Rej. & $\\overline{ev}$ & $\\overline{k^*}$ & Rej. & $\\overline{ev}$ & $\\overline{k^*}$ & Rej. & $\\overline{ev}$ & $\\overline{k^*}$ \\\\",
  "\\midrule"
)

n_per_block <- length(prior_levels)
for (i in seq(1, nrow(all_results), by = n_per_block)) {
  rows <- lapply(0:(n_per_block - 1), function(k) all_results[i + k, ])
  r_ni <- rows[[1]]; r_in <- rows[[2]]; r_co <- rows[[3]]
  tex_table <- c(tex_table,
    sprintf(paste0("%d & %.2f & %.0f\\%% & %.3f & %.3f ",
                   "& %.0f\\%% & %.3f & %.3f ",
                   "& %.0f\\%% & %.3f & %.3f \\\\"),
            r_ni$n, r_ni$d,
            r_ni$reject_rate*100, r_ni$mean_ev, r_ni$mean_k_star,
            r_in$reject_rate*100, r_in$mean_ev, r_in$mean_k_star,
            r_co$reject_rate*100, r_co$mean_ev, r_co$mean_k_star))
}

tex_table <- c(tex_table, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_table, "output/simulation_replicated_table.tex")
cat("\n--> output/simulation_replicated_table.tex\n")

# ============================================================================
# 6. PLOTS
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")

prior_colors <- c("Non-informative" = "steelblue",
                  "Informative (mu=0.5)" = "tomato",
                  "Conflict (mu=0.1)" = "darkgreen")

# Power curves
p_power <- ggplot(all_results,
                  aes(x = n, y = reject_rate * 100, color = prior, linetype = prior)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = prior_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Sample size (n)", y = "Rejection rate (%)",
       title = "FBST power: rejection rate vs sample size",
       subtitle = sprintf("R=%d replicates per scenario", R_reps),
       color = "Prior", linetype = "Prior") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Figures/power_curves.png", p_power, width = 11, height = 7, dpi = 150)
cat("  --> Figures/power_curves.png\n")

# Heatmap of rejection rates
p_heatmap <- ggplot(all_results,
                    aes(x = factor(n), y = d_label, fill = reject_rate * 100)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", reject_rate * 100)),
            color = "white", fontface = "bold", size = 3.5) +
  facet_wrap(~prior, ncol = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 50, limits = c(0, 100),
                       name = "Rej. rate (%)") +
  labs(x = "Sample size (n)", y = "Effect size",
       title = "Rejection rates across priors") +
  theme_bw()
ggsave("Figures/rejection_heatmap.png", p_heatmap, width = 14, height = 5, dpi = 150)
cat("  --> Figures/rejection_heatmap.png\n")

# Pairwise difference heatmap
agreement_data <- all_results %>%
  select(n, d, d_label, prior, reject_rate) %>%
  pivot_wider(names_from = prior, values_from = reject_rate) %>%
  mutate(
    `NI minus INF`  = `Non-informative`     - `Informative (mu=0.5)`,
    `NI minus CONF` = `Non-informative`     - `Conflict (mu=0.1)`,
    `INF minus CONF`= `Informative (mu=0.5)` - `Conflict (mu=0.1)`
  )

diff_long <- agreement_data %>%
  select(n, d_label, `NI minus INF`, `NI minus CONF`, `INF minus CONF`) %>%
  pivot_longer(cols = c("NI minus INF", "NI minus CONF", "INF minus CONF"),
               names_to = "comparison", values_to = "diff")

p_agreement <- ggplot(diff_long, aes(x = factor(n), y = d_label,
                                     fill = diff * 100)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = sprintf("%+.0f", diff * 100)),
            color = "black", fontface = "bold", size = 3.5) +
  facet_wrap(~comparison, ncol = 3) +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "steelblue",
                       midpoint = 0, name = "Diff (pp)") +
  labs(x = "Sample size (n)", y = "Effect size",
       title = "Pairwise differences in rejection rates (percentage points)") +
  theme_bw()
ggsave("Figures/decision_agreement.png", p_agreement, width = 14, height = 5, dpi = 150)
cat("  --> Figures/decision_agreement.png\n")

# Mean ev and k* per prior
results_long <- all_results %>%
  select(n, d_label, prior, mean_ev, mean_k_star) %>%
  pivot_longer(c("mean_ev", "mean_k_star"),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("mean_ev", "mean_k_star"),
                         labels = c("Mean ev_obs", "Mean k*")))

p_means <- ggplot(results_long,
                  aes(x = n, y = value, color = prior, linetype = metric)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = prior_colors) +
  labs(x = "Sample size (n)", y = "Value",
       title = "Mean observed e-value and adaptive cutoff k*",
       color = "Prior", linetype = "Metric") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Figures/ev_kstar_means.png", p_means, width = 11, height = 7, dpi = 150)
cat("  --> Figures/ev_kstar_means.png\n")

# ============================================================================
# 7. TYPE I / POWER SUMMARY
# ============================================================================

cat("\n--- Type I error (d = 0, H true) ---\n")
type1 <- all_results %>% filter(d == 0)
print(type1[, c("n", "prior", "reject_rate")], row.names = FALSE, digits = 3)

cat("\n--- Mean power (d > 0) by effect size ---\n")
power_summary <- all_results %>%
  filter(d > 0) %>%
  group_by(d_label, prior) %>%
  summarise(mean_power = mean(reject_rate), .groups = "drop")
print(power_summary, row.names = FALSE, digits = 3)

cat(sprintf("\nTotal time: %.1f minutes\n", elapsed))
cat("\nGenerated:\n")
cat("  output/simulation_replicated.csv\n")
cat("  output/simulation_replicated_table.tex\n")
cat("  Figures/power_curves.png\n")
cat("  Figures/rejection_heatmap.png\n")
cat("  Figures/decision_agreement.png\n")
cat("  Figures/ev_kstar_means.png\n")
