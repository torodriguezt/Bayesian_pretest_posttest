# simulation_replicated_weights.R
# Replicated simulation study with multiple (a, b) weight configurations
# Compares THREE priors (KL-fitted, see priors_config.R):
#   - Non-informative (KL to U)
#   - Informative symmetric (mu=0.5, N=50)
#   - Conflict (mu=0.1, N=50)
#
# Optimization: ev_H, ev_A are simulated ONCE per (n, d, replicate, prior),
# then find_kstar() is called for each (a, b) - that step is cheap.

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

weight_configs <- data.frame(
  label = c("a=1, b=1 (balanced)",
            "a=5, b=1 (moderate)",
            "a=20, b=1 (strict)"),
  a_w   = c(1, 5, 20),
  b_w   = c(1, 1, 1),
  stringsAsFactors = FALSE
)

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

R_reps <- 200
M_sim  <- 1000
set.seed(42)

# ============================================================================
# 2. ANALYSIS FUNCTION
# ============================================================================

analyze_one <- function(n, x1, x2, prior_alpha, M = 300) {
  a0 <- prior_alpha["a0"]; a1 <- prior_alpha["a1"]; a2 <- prior_alpha["a2"]
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  sup_H  <- find_sup_H(consts)$sup_H
  ev_obs <- ev_quad(consts, sup_H)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M)
  list(ev_obs = ev_obs, ev_H = ev_H, ev_A = ev_A)
}

# ============================================================================
# 3. MAIN LOOP
# ============================================================================

cat("\n========== REPLICATED SIM (multiple weights x 3 priors) ==========\n")
cat(sprintf("R = %d replicates per scenario\n", R_reps))
cat(sprintf("M = %d e-value simulations per analysis\n", M_sim))
cat(sprintf("Priors: %s\n", paste(prior_levels, collapse = ", ")))
cat(sprintf("Weights: %d configs\n", nrow(weight_configs)))
cat(sprintf("Total analyses: %d x %d priors = %d\n\n",
            nrow(scenarios) * R_reps, length(priors),
            nrow(scenarios) * R_reps * length(priors)))

all_results <- list()
t_start <- Sys.time()

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  cat(sprintf("[%2d/%d] n=%3d, %s ", i, nrow(scenarios), s$n, s$d_label))

  # rej[r, w, prior]: reject indicator;  k_mat[r, w, prior]: k* values
  rej   <- array(FALSE, dim = c(R_reps, nrow(weight_configs), length(priors)),
                 dimnames = list(NULL, weight_configs$label, prior_levels))
  k_mat <- array(NA,    dim = dim(rej), dimnames = dimnames(rej))
  ev_obs_mat <- matrix(NA_real_, R_reps, length(priors),
                       dimnames = list(NULL, prior_levels))

  for (r in seq_len(R_reps)) {
    x1 <- rbinom(1, s$n, s$theta1)
    x2 <- rbinom(1, s$n, s$theta2)

    for (pname in prior_levels) {
      res <- analyze_one(s$n, x1, x2, priors[[pname]], M = M_sim)
      ev_obs_mat[r, pname] <- res$ev_obs

      for (w in seq_len(nrow(weight_configs))) {
        opt <- find_kstar(res$ev_H, res$ev_A,
                          weight_configs$a_w[w], weight_configs$b_w[w])
        k_mat[r, w, pname] <- opt$k_star
        rej[r, w, pname]   <- (res$ev_obs <= opt$k_star)
      }
    }
    if (r %% 10 == 0) cat(".")
  }

  # Aggregate per (prior, weight)
  for (pname in prior_levels) {
    for (w in seq_len(nrow(weight_configs))) {
      all_results[[length(all_results) + 1]] <- data.frame(
        n            = s$n,
        d            = s$d,
        d_label      = s$d_label,
        weight_label = weight_configs$label[w],
        a_w          = weight_configs$a_w[w],
        b_w          = weight_configs$b_w[w],
        prior        = pname,
        reject_rate  = mean(rej[, w, pname]),
        mean_ev      = mean(ev_obs_mat[, pname]),
        mean_k_star  = mean(k_mat[, w, pname]),
        sd_k_star    = sd(k_mat[, w, pname])
      )
    }
  }

  # Compact line of progress (a=1 rates per prior)
  cat(sprintf(" a=1: %s\n",
    paste(sapply(prior_levels, function(p)
      sprintf("%s=%.0f%%", substr(p, 1, 3), mean(rej[, 1, p]) * 100)),
      collapse = " ")))
}

all_results <- do.call(rbind, all_results)
all_results$prior <- factor(all_results$prior, levels = prior_levels)

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\nTotal time: %.1f minutes\n", elapsed))

# ============================================================================
# 4. SUMMARIES
# ============================================================================

cat("\n========== TYPE I ERROR (d = 0) ==========\n")
type1 <- all_results %>%
  filter(d == 0) %>%
  group_by(weight_label, prior) %>%
  summarise(mean_type1 = mean(reject_rate), .groups = "drop")
print(as.data.frame(type1), row.names = FALSE, digits = 3)

cat("\n========== POWER (d > 0) ==========\n")
power_summary <- all_results %>%
  filter(d > 0) %>%
  group_by(d_label, weight_label, prior) %>%
  summarise(mean_power = mean(reject_rate), .groups = "drop")
print(as.data.frame(power_summary), row.names = FALSE, digits = 3)

write.csv(all_results, "output/simulation_replicated_weights.csv",
          row.names = FALSE)

# ============================================================================
# 5. PLOTS
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")

weight_colors <- c("a=1, b=1 (balanced)" = "steelblue",
                   "a=5, b=1 (moderate)" = "darkorange",
                   "a=20, b=1 (strict)"  = "tomato")

# Power curves: rows = prior, cols = effect size
p_power <- ggplot(all_results,
                  aes(x = n, y = reject_rate * 100,
                      color = weight_label, linetype = weight_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_grid(prior ~ d_label) +
  scale_color_manual(values = weight_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Sample size (n)", y = "Rejection rate (%)",
       title = "FBST: rejection rate vs n, by prior x weight",
       subtitle = sprintf("R=%d replicates", R_reps),
       color = "Weight (a, b)", linetype = "Weight (a, b)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Figures/power_curves_weights.png", p_power,
       width = 14, height = 8, dpi = 150)
cat("  --> Figures/power_curves_weights.png\n")

# Type I error focus
p_type1 <- all_results %>%
  filter(d == 0) %>%
  ggplot(aes(x = n, y = reject_rate * 100,
             color = weight_label, linetype = prior)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_hline(yintercept = 5, linetype = "dotted", color = "grey40") +
  annotate("text", x = 30, y = 7, label = "5% (frequentist)",
           color = "grey40", size = 3) +
  scale_color_manual(values = weight_colors) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  labs(x = "Sample size (n)", y = "Type I error rate (%)",
       title = "Type I error when H is true (Delta = 0)",
       color = "Weight (a, b)", linetype = "Prior") +
  theme_bw() + theme(legend.position = "bottom", legend.box = "vertical")
ggsave("Figures/type1_error_weights.png", p_type1, width = 11, height = 6, dpi = 150)
cat("  --> Figures/type1_error_weights.png\n")

# Heatmap: prior x weight, n vs effect
p_heatmap <- ggplot(all_results,
                    aes(x = factor(n), y = d_label,
                        fill = reject_rate * 100)) +
  geom_tile(color = "black", size = 0.4) +
  geom_text(aes(label = sprintf("%.0f", reject_rate * 100)),
            color = "white", fontface = "bold", size = 3.0) +
  facet_grid(prior ~ weight_label) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 50, limits = c(0, 100),
                       name = "Rej.\nrate (%)") +
  labs(x = "Sample size (n)", y = "Effect size",
       title = "Rejection rates: priors x weight configurations") +
  theme_bw()
ggsave("Figures/rejection_heatmap_weights.png", p_heatmap,
       width = 14, height = 9, dpi = 150)
cat("  --> Figures/rejection_heatmap_weights.png\n")

# k* by weight (one panel per prior)
p_kstar <- ggplot(all_results,
                  aes(x = n, y = mean_k_star, color = weight_label)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  geom_ribbon(aes(ymin = pmax(0, mean_k_star - sd_k_star),
                  ymax = pmin(1, mean_k_star + sd_k_star),
                  fill = weight_label), alpha = 0.15, color = NA) +
  facet_grid(prior ~ d_label) +
  scale_color_manual(values = weight_colors) +
  scale_fill_manual(values = weight_colors) +
  labs(x = "Sample size (n)", y = "Mean k* (+/- SD)",
       title = "Adaptive cutoff k* by prior x weight",
       color = "Weight (a, b)", fill = "Weight (a, b)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Figures/kstar_by_weights.png", p_kstar, width = 14, height = 8, dpi = 150)
cat("  --> Figures/kstar_by_weights.png\n")

# ============================================================================
# 6. LATEX TABLES (one block per prior, stacked vertically)
# ============================================================================

build_table_for_d <- function(d_target, label, caption) {
  tab <- all_results %>%
    filter(d == d_target) %>%
    select(n, prior, weight_label, reject_rate) %>%
    pivot_wider(names_from = c(prior, weight_label), values_from = reject_rate)

  hd <- c(
    "\\begin{table}[!h]", "\\centering", "\\footnotesize",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{c|ccc|ccc|ccc}", "\\toprule",
    "\\multicolumn{1}{c}{} & \\multicolumn{3}{c}{\\textbf{Non-inf.}} & \\multicolumn{3}{c}{\\textbf{Informative}} & \\multicolumn{3}{c}{\\textbf{Conflict}} \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
    "$n$ & $a{=}1$ & $a{=}5$ & $a{=}20$ & $a{=}1$ & $a{=}5$ & $a{=}20$ & $a{=}1$ & $a{=}5$ & $a{=}20$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    hd <- c(hd, sprintf(
      paste0("%d & %.0f\\%% & %.0f\\%% & %.0f\\%% ",
             "& %.0f\\%% & %.0f\\%% & %.0f\\%% ",
             "& %.0f\\%% & %.0f\\%% & %.0f\\%% \\\\"),
      r$n,
      r[["Non-informative_a=1, b=1 (balanced)"]] * 100,
      r[["Non-informative_a=5, b=1 (moderate)"]] * 100,
      r[["Non-informative_a=20, b=1 (strict)"]]  * 100,
      r[["Informative (mu=0.5)_a=1, b=1 (balanced)"]] * 100,
      r[["Informative (mu=0.5)_a=5, b=1 (moderate)"]] * 100,
      r[["Informative (mu=0.5)_a=20, b=1 (strict)"]]  * 100,
      r[["Conflict (mu=0.1)_a=1, b=1 (balanced)"]] * 100,
      r[["Conflict (mu=0.1)_a=5, b=1 (moderate)"]] * 100,
      r[["Conflict (mu=0.1)_a=20, b=1 (strict)"]]  * 100
    ))
  }
  c(hd, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

tex_t1 <- build_table_for_d(
  d_target = 0,
  label    = "tab:type1_weights",
  caption  = sprintf(paste0(
    "Type I error rates (\\%%) under different weight configurations and ",
    "priors. $R = %d$ replicates, $\\Delta = 0$."), R_reps))
writeLines(tex_t1, "output/type1_weights_table.tex")
cat("  --> output/type1_weights_table.tex\n")

tex_small <- build_table_for_d(
  d_target = 0.02,
  label    = "tab:power_small_weights",
  caption  = sprintf(paste0(
    "Power (\\%%) under different weight configurations and priors. ",
    "$R = %d$ replicates, $\\Delta = 0.02$ (small effect)."), R_reps))
writeLines(tex_small, "output/power_small_weights_table.tex")
cat("  --> output/power_small_weights_table.tex\n")

tex_pow <- build_table_for_d(
  d_target = 0.10,
  label    = "tab:power_medium_weights",
  caption  = sprintf(paste0(
    "Power (\\%%) under different weight configurations and priors. ",
    "$R = %d$ replicates, $\\Delta = 0.10$ (medium effect)."), R_reps))
writeLines(tex_pow, "output/power_medium_weights_table.tex")
cat("  --> output/power_medium_weights_table.tex\n")

tex_large <- build_table_for_d(
  d_target = 0.30,
  label    = "tab:power_large_weights",
  caption  = sprintf(paste0(
    "Power (\\%%) under different weight configurations and priors. ",
    "$R = %d$ replicates, $\\Delta = 0.30$ (large effect)."), R_reps))
writeLines(tex_large, "output/power_large_weights_table.tex")
cat("  --> output/power_large_weights_table.tex\n")

# ============================================================================
# 7. FINAL
# ============================================================================

cat(sprintf("\nTotal time: %.1f minutes\n", elapsed))
cat("\nGenerated files:\n")
cat("  output/simulation_replicated_weights.csv\n")
cat("  output/type1_weights_table.tex\n")
cat("  output/power_small_weights_table.tex\n")
cat("  output/power_medium_weights_table.tex\n")
cat("  output/power_large_weights_table.tex\n")
cat("  Figures/power_curves_weights.png\n")
cat("  Figures/type1_error_weights.png\n")
cat("  Figures/rejection_heatmap_weights.png\n")
cat("  Figures/kstar_by_weights.png\n")
