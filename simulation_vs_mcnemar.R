# simulation_vs_mcnemar.R
# -----------------------------------------------------------------------------
# Replicated simulation comparison: FBST (adaptive k*) vs McNemar's test.
#
# Genera datos pareados a nivel individual (y_i1, y_i2) para poder calcular
# tanto el FBST (usa X1, X2 agregados) como McNemar (usa pares discordantes).
# Prior usada para FBST: no informativa (KL-optima).
#
# Salida:
#   output/simulation_vs_mcnemar.csv
#   output/comparison_mcnemar_d*.tex  (una tabla por tamanio de efecto)
#   Figures/fbst_vs_mcnemar_power.png
# -----------------------------------------------------------------------------

library(Rcpp)
library(ggplot2)
library(dplyr)
library(tidyr)

sourceCpp("BivBetaBinom.cpp")
source("priors_config.R")

dir.create("output",  showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# ============================================================================
# 1. PARAMETROS
# ============================================================================

weight_configs <- data.frame(
  label = c("FBST (a=1)", "FBST (a=5)", "FBST (a=20)"),
  a_w   = c(1, 5, 20),
  b_w   = c(1, 1,  1),
  stringsAsFactors = FALSE
)

# Se incluyen n pequenos (15, 25) para capturar la ventaja del FBST ahi
n_grid <- c(15, 25, 30, 50, 75, 100, 150, 200, 300, 400, 450)

scenarios <- expand.grid(
  n             = n_grid,
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
# RESUME: load existing results and skip already-computed n values
# ============================================================================
csv_path <- "output/simulation_vs_mcnemar.csv"
if (file.exists(csv_path)) {
  existing <- read.csv(csv_path, stringsAsFactors = FALSE)
  existing_n <- sort(unique(existing$n))
  cat(sprintf("Resuming: found existing results for n in {%s}\n",
              paste(existing_n, collapse = ", ")))
  scenarios_to_run <- scenarios[!scenarios$n %in% existing_n, ]
  if (nrow(scenarios_to_run) == 0) {
    cat("All scenarios already computed. Regenerating plots only.\n")
  } else {
    cat(sprintf("New n values to run: {%s}\n\n",
                paste(sort(unique(scenarios_to_run$n)), collapse = ", ")))
  }
} else {
  existing        <- NULL
  scenarios_to_run <- scenarios
}

# ============================================================================
# 2. FUNCIONES
# ============================================================================

# Genera n pares (y_i1, y_i2) independientes dado (theta1, theta2).
# Devuelve los agregados que necesita cada test.
generate_paired <- function(n, theta1, theta2) {
  y1 <- rbinom(n, 1, theta1)
  y2 <- rbinom(n, 1, theta2)
  list(
    X1 = sum(y1),
    X2 = sum(y2),
    b  = sum(y1 == 0L & y2 == 1L),  # (0->1) pares discordantes
    c_ = sum(y1 == 1L & y2 == 0L)   # (1->0) pares discordantes
  )
}

# McNemar exacto cuando b+c < 25, chi-cuadrado de lo contrario.
# Retorna TRUE si rechaza H: theta1 = theta2 con alpha = 0.05.
mcnemar_reject <- function(b, c_, alpha = 0.05) {
  bc <- b + c_
  if (bc == 0L) return(FALSE)
  if (bc < 25L) {
    # Test exacto binomial: bajo H, b ~ Bin(b+c, 0.5)
    pval <- 2 * min(pbinom(b, bc, 0.5), pbinom(c_, bc, 0.5))
  } else {
    pval <- pchisq((b - c_)^2 / bc, df = 1L, lower.tail = FALSE)
  }
  pval < alpha
}

# ============================================================================
# 3. LOOP PRINCIPAL
# ============================================================================

cat("\n========== FBST vs McNemar (simulacion replicada) ==========\n")
cat(sprintf("R = %d replicates | M = %d e-value sims | Prior: KL no-informativa\n",
            R_reps, M_sim))
cat(sprintf("n grid: %s\n", paste(n_grid, collapse = ", ")))
cat(sprintf("Escenarios nuevos: %d\n\n", nrow(scenarios_to_run)))

all_results <- list()
t_start     <- Sys.time()

for (i in seq_len(nrow(scenarios_to_run))) {
  s <- scenarios_to_run[i, ]
  cat(sprintf("[%2d/%d] n=%3d, %s ", i, nrow(scenarios_to_run), s$n, s$d_label))

  mcn_rej  <- logical(R_reps)
  fbst_rej <- matrix(FALSE, R_reps, nrow(weight_configs))

  a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]

  for (r in seq_len(R_reps)) {
    dat <- generate_paired(s$n, s$theta1, s$theta2)

    # McNemar
    mcn_rej[r] <- mcnemar_reject(dat$b, dat$c_)

    # FBST: calcular ev_H y ev_A una sola vez, aplicar cada peso
    ev_obs <- ev_quad_from_data(s$n, s$n, dat$X1, dat$X2, a0, a1, a2)
    ev_H   <- simulate_evs_H_post(s$n, s$n, dat$X1, dat$X2, a0, a1, a2, M_sim)
    ev_A   <- simulate_evs_A_post(s$n, s$n, dat$X1, dat$X2, a0, a1, a2, M_sim)

    for (w in seq_len(nrow(weight_configs))) {
      opt           <- find_kstar(ev_H, ev_A, weight_configs$a_w[w],
                                  weight_configs$b_w[w])
      fbst_rej[r, w] <- (ev_obs <= opt$k_star)
    }
    if (r %% 50 == 0) cat(".")
  }

  # Guardar resultados
  all_results[[length(all_results) + 1]] <- data.frame(
    n = s$n, d = s$d, d_label = s$d_label,
    label = "McNemar", reject_rate = mean(mcn_rej)
  )
  for (w in seq_len(nrow(weight_configs))) {
    all_results[[length(all_results) + 1]] <- data.frame(
      n = s$n, d = s$d, d_label = s$d_label,
      label = weight_configs$label[w],
      reject_rate = mean(fbst_rej[, w])
    )
  }

  cat(sprintf(" McNemar=%.0f%% | FBST(a=1)=%.0f%% FBST(a=5)=%.0f%%\n",
              mean(mcn_rej) * 100,
              mean(fbst_rej[, 1]) * 100,
              mean(fbst_rej[, 2]) * 100))
}

if (length(all_results) > 0) {
  new_results <- do.call(rbind, all_results)
  if (!is.null(existing)) {
    all_results <- dplyr::bind_rows(existing, new_results)
  } else {
    all_results <- new_results
  }
} else {
  all_results <- existing
}
method_order      <- c("McNemar", "FBST (a=1)", "FBST (a=5)", "FBST (a=20)")
all_results$label <- factor(all_results$label, levels = method_order)

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\nTiempo total: %.1f minutos\n", elapsed))

# ============================================================================
# 4. RESUMEN
# ============================================================================

cat("\n========== ERROR TIPO I (d = 0) ==========\n")
type1 <- all_results %>% filter(d == 0) %>%
  group_by(label) %>%
  summarise(type1 = mean(reject_rate), .groups = "drop")
print(as.data.frame(type1), row.names = FALSE, digits = 3)

cat("\n========== PODER (d > 0) ==========\n")
power_tab <- all_results %>% filter(d > 0) %>%
  group_by(d_label, label) %>%
  summarise(power = mean(reject_rate), .groups = "drop")
print(as.data.frame(power_tab), row.names = FALSE, digits = 3)

write.csv(all_results, "output/simulation_vs_mcnemar.csv", row.names = FALSE)

# ============================================================================
# 5. FIGURA: curvas de poder
# ============================================================================

method_colors <- c(
  "McNemar"     = "black",
  "FBST (a=1)"  = "steelblue",
  "FBST (a=5)"  = "darkorange",
  "FBST (a=20)" = "tomato"
)
method_lines <- c(
  "McNemar"     = "dashed",
  "FBST (a=1)"  = "solid",
  "FBST (a=5)"  = "solid",
  "FBST (a=20)" = "solid"
)

p_power <- ggplot(all_results,
                  aes(x = n, y = reject_rate * 100,
                      color = label, linetype = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_hline(yintercept = 5, linetype = "dotted", color = "grey40") +
  annotate("text", x = 30, y = 8, label = "5%", color = "grey40", size = 3) +
  facet_wrap(~d_label, ncol = 2) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_lines) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Sample size (n)", y = "Rejection rate (%)",
       title = "FBST vs McNemar: tasas de rechazo",
       subtitle = sprintf("R=%d replicates | Prior KL no-informativa", R_reps),
       color = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("Figures/fbst_vs_mcnemar_power.png", p_power,
       width = 12, height = 7, dpi = 150)
cat("  --> Figures/fbst_vs_mcnemar_power.png\n")

# ============================================================================
# 6. TABLAS LATEX (una por Delta)
# ============================================================================

build_comparison_table <- function(d_target, tab_label, caption) {
  tab <- all_results %>%
    filter(d == d_target) %>%
    select(n, label, reject_rate) %>%
    pivot_wider(names_from = label, values_from = reject_rate) %>%
    arrange(n)

  hd <- c(
    "\\begin{table}[!h]", "\\centering", "\\footnotesize",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", tab_label, "}"),
    "\\begin{tabular}{c|c|ccc}", "\\toprule",
    paste0("$n$ & \\textbf{McNemar} & \\textbf{FBST ($a=1$)} ",
           "& \\textbf{FBST ($a=5$)} & \\textbf{FBST ($a=20$)} \\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    hd <- c(hd, sprintf(
      "%d & %.0f\\%% & %.0f\\%% & %.0f\\%% & %.0f\\%% \\\\",
      r$n,
      r[["McNemar"]]     * 100,
      r[["FBST (a=1)"]]  * 100,
      r[["FBST (a=5)"]]  * 100,
      r[["FBST (a=20)"]] * 100
    ))
  }
  c(hd, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

d_scenarios <- list(
  list(d = 0.00, lbl = "tab:comp_type1",
       cap = sprintf("Tasa de Error Tipo I (\\%%): FBST vs McNemar. $R=%d$, $\\Delta=0$.", R_reps)),
  list(d = 0.02, lbl = "tab:comp_small",
       cap = sprintf("Poder (\\%%): FBST vs McNemar. $R=%d$, $\\Delta=0.02$ (efecto peque\\~no).", R_reps)),
  list(d = 0.10, lbl = "tab:comp_medium",
       cap = sprintf("Poder (\\%%): FBST vs McNemar. $R=%d$, $\\Delta=0.10$ (efecto mediano).", R_reps)),
  list(d = 0.30, lbl = "tab:comp_large",
       cap = sprintf("Poder (\\%%): FBST vs McNemar. $R=%d$, $\\Delta=0.30$ (efecto grande).", R_reps))
)

for (sc in d_scenarios) {
  d_str <- gsub("\\.", "", sprintf("%.2f", sc$d))
  fname <- sprintf("output/comparison_mcnemar_d%s.tex", d_str)
  writeLines(build_comparison_table(sc$d, sc$lbl, sc$cap), fname)
  cat(sprintf("  --> %s\n", fname))
}

cat(sprintf("\nTiempo total: %.1f minutos\n", elapsed))
cat("\nArchivos generados:\n")
cat("  output/simulation_vs_mcnemar.csv\n")
cat("  output/comparison_mcnemar_d000.tex  (Tipo I)\n")
cat("  output/comparison_mcnemar_d002.tex  (poder Delta=0.02)\n")
cat("  output/comparison_mcnemar_d010.tex  (poder Delta=0.10)\n")
cat("  output/comparison_mcnemar_d030.tex  (poder Delta=0.30)\n")
cat("  Figures/fbst_vs_mcnemar_power.png\n")
