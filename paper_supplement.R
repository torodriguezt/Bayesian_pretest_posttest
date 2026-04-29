# paper_supplement.R
# Suplemento completo para el artículo. Genera:
#
#   output/kstar_informative_table.tex  — Tab3/4 corregidas (prior α=10)
#   Figures/mcnemar_comparison.png      — FBST vs McNemar en escenarios sintéticos
#   Figures/decision_vs_n_clean.png     — k* adaptativo vs tamaño muestral
#   Figures/sensitivity_hyperparams.png — sensibilidad a hiperparámetros
#   Figures/nonrejection_example.png    — ejemplo donde el test NO rechaza

library(Rcpp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ALA)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")
dir.create("output",  showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# Hiperparámetros (KL-fitted, ver priors_config.R)
source("priors_config.R")
a0_kl   <- prior_NI["a0"];   a1_kl   <- prior_NI["a1"];   a2_kl   <- prior_NI["a2"]
a0_inf  <- prior_INF["a0"];  a1_inf  <- prior_INF["a1"];  a2_inf  <- prior_INF["a2"]
a0_conf <- prior_CONF["a0"]; a1_conf <- prior_CONF["a1"]; a2_conf <- prior_CONF["a2"]

M      <- 1500
set.seed(42)

# ===========================================================================
# 1. Tab3/Tab4 CORREGIDAS — prior informativa α=10, posterior-based
#    (mismos 4 grupos THKS que Tab2 pero con prior distinta)
# ===========================================================================
cat("\n=== 1. Tablas corregidas con prior informativa (α=10) ===\n")

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
  yy = list(g = prep_school(datos1, "404", "yes", "yes"), lab = "CC + TV"),
  yn = list(g = prep_school(datos1, "408", "yes", "no"),  lab = "CC, no TV"),
  ny = list(g = prep_global(datos1,         "no",  "yes"), lab = "no CC, TV"),
  nn = list(g = prep_school(datos1, "409", "no",  "no"),  lab = "no CC, no TV")
)

run_group <- function(g_entry, a0, a1, a2, M, seed = 7) {
  g  <- g_entry$g; lab <- g_entry$lab
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  set.seed(seed)
  ev_obs <- ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(n1, n2, x1, x2, a0, a1, a2, M, 401)
  ev_A   <- simulate_evs_A_post(n1, n2, x1, x2, a0, a1, a2, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  data.frame(group = lab, n1 = n1, n2 = n2, x1 = x1, x2 = x2,
             ev_obs = ev_obs, k_star = opt$k_star,
             alpha = opt$alpha, beta = opt$beta,
             decision = ifelse(ev_obs <= opt$k_star, "reject H", "do not reject"))
}

cat("Procesando grupos con prior informativa (α=10)...\n")
tab_inf <- do.call(rbind, lapply(names(grupos), function(k) {
  cat(" -", grupos[[k]]$lab, "\n")
  run_group(grupos[[k]], a0_inf, a1_inf, a2_inf, M)
}))
print(tab_inf, digits = 4, row.names = FALSE)

tex_inf <- c(
  "\\begin{table}[!h]", "\\centering",
  "\\caption{Cutoff óptimo $k^{*}$ por tratamiento, prior informativa",
  "($\\alpha_0=\\alpha_1=\\alpha_2=10$), formulación posterior.}",
  "\\label{tab:kstar_informative}",
  "\\begin{tabular}{lccccccl}", "\\toprule",
  "Tratamiento & $n_1$ & $n_2$ & $x_1$ & $x_2$ & $\\mathrm{ev}_{\\mathrm{obs}}$ & $k^{*}$ & Decisión \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(tab_inf))) {
  r <- tab_inf[i, ]
  tex_inf <- c(tex_inf,
    sprintf("%s & %d & %d & %d & %d & %.4f & %.4f & %s \\\\",
            r$group, r$n1, r$n2, r$x1, r$x2, r$ev_obs, r$k_star, r$decision))
}
tex_inf <- c(tex_inf, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_inf, "output/kstar_informative_table.tex")
cat("→ output/kstar_informative_table.tex\n")

# ===========================================================================
# 1b. McNemar sobre los datos reales del THKS
#     tvsfp tiene id con 2 filas por sujeto (pre/post) → emparejamiento exacto
# ===========================================================================
cat("\n=== 1b. McNemar en datos reales THKS ===\n")

mcnemar_pval <- function(b, c) {
  if (b + c == 0) return(1)
  pchisq((b - c)^2 / (b + c), df = 1, lower.tail = FALSE)
}

# Extrae tabla 2×2 pareada para un subgrupo
paired_2x2 <- function(d) {
  wide <- d %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    select(id, stage, binTHKS) %>%
    pivot_wider(names_from = stage, values_from = binTHKS) %>%
    filter(!is.na(pre), !is.na(post))
  list(
    n  = nrow(wide),
    a  = sum(wide$pre == 1 & wide$post == 1),  # éxito → éxito
    b  = sum(wide$pre == 1 & wide$post == 0),  # éxito → fallo  (discordante)
    c  = sum(wide$pre == 0 & wide$post == 1),  # fallo → éxito  (discordante)
    d  = sum(wide$pre == 0 & wide$post == 0)   # fallo → fallo
  )
}

subsets <- list(
  "CC + TV"      = datos1 %>% filter(school == "404", school.based == "yes", tv.based == "yes"),
  "CC, no TV"    = datos1 %>% filter(school == "408", school.based == "yes", tv.based == "no"),
  "no CC, TV"    = datos1 %>% filter(school.based == "no",  tv.based == "yes"),
  "no CC, no TV" = datos1 %>% filter(school == "409", school.based == "no",  tv.based == "no")
)

mc_real <- do.call(rbind, lapply(names(subsets), function(lab) {
  tbl   <- paired_2x2(subsets[[lab]])
  p_val <- mcnemar_pval(tbl$b, tbl$c)
  cat(sprintf("  %-14s  n=%d  b=%d  c=%d  p=%.4f  → %s\n",
              lab, tbl$n, tbl$b, tbl$c, p_val,
              ifelse(p_val <= 0.05, "rechaza H", "no rechaza")))
  data.frame(group = lab, n = tbl$n,
             b = tbl$b, c = tbl$c,
             concordantes = tbl$a + tbl$d,
             p_mcnemar = p_val,
             mcnemar_dec = ifelse(p_val <= 0.05, "rechaza H", "no rechaza"))
}))

# Tabla combinada FBST (prior KL) + FBST (prior inf) + McNemar
tab_kl_real <- do.call(rbind, lapply(names(grupos), function(k) {
  r <- run_group(grupos[[k]], a0_kl, a1_kl, a2_kl, M)
  r$prior <- "KL-óptima"
  r
}))

combined <- tab_kl_real %>%
  left_join(tab_inf %>% select(group, ev_obs, k_star, decision) %>%
              rename(ev_inf = ev_obs, kstar_inf = k_star, dec_inf = decision),
            by = "group") %>%
  left_join(mc_real %>% select(group, b, c, p_mcnemar, mcnemar_dec),
            by = "group")

cat("\n--- Tabla comparativa completa ---\n")
print(combined[, c("group", "n1",
                   "ev_obs",  "k_star",  "decision",
                   "ev_inf",  "kstar_inf","dec_inf",
                   "b", "c",  "p_mcnemar","mcnemar_dec")],
      digits = 3, row.names = FALSE)

# LaTeX con las tres columnas de decisión
tex_cmp <- c(
  "\\begin{table}[!h]", "\\centering",
  "\\caption{Comparación de decisiones: FBST (prior KL-óptima), FBST (prior informativa $\\alpha=10$) y McNemar, datos THKS.}",
  "\\label{tab:comparison}",
  "\\begin{tabular}{lcccccc}", "\\toprule",
  "Grupo & $n$ & $\\mathrm{ev}_{\\mathrm{obs}}$ (KL) & Dec.~FBST-KL & Dec.~FBST-Inf & $p_{\\mathrm{McN}}$ & Dec.~McNemar \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(combined))) {
  r <- combined[i, ]
  tex_cmp <- c(tex_cmp,
    sprintf("%s & %d & %.4f & %s & %s & %.4f & %s \\\\",
            r$group, r$n1,
            r$ev_obs,  r$decision,
            r$dec_inf, r$p_mcnemar, r$mcnemar_dec))
}
tex_cmp <- c(tex_cmp, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_cmp, "output/comparison_table.tex")
cat("→ output/comparison_table.tex\n")

# ===========================================================================
# 2. Comparación FBST vs McNemar
#    McNemar requiere pares (éxito pretest, éxito posttest) por individuo.
#    Aquí simulamos datos pareados sintéticos con distintos efectos.
# ===========================================================================
cat("\n=== 2. Comparación FBST vs McNemar (escenarios sintéticos) ===\n")

# McNemar: dado tabla 2×2 de pares (s1, s2):
#   b = # pares (1→0),  c = # pares (0→1)
#   estadístico = (b-c)²/(b+c),  p-valor = pchisq(stat, 1, lower.tail=FALSE)
mcnemar_pval <- function(b, c) {
  if (b + c == 0) return(1)
  pchisq((b - c)^2 / (b + c), df = 1, lower.tail = FALSE)
}

# Simula n pares con probabilidades (theta1, theta2) y calcula ambos tests
compare_tests <- function(n, theta1, theta2, a0, a1, a2, M = 1000,
                          alpha_level = 0.05, seed = 1) {
  set.seed(seed)
  s1 <- rbinom(n, 1, theta1)
  s2 <- rbinom(n, 1, theta2)
  x1 <- sum(s1); x2 <- sum(s2)
  b  <- sum(s1 == 1 & s2 == 0)  # 1→0
  c_ <- sum(s1 == 0 & s2 == 1)  # 0→1

  ev_obs <- ev_quad_from_data(n, n, x1, x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M, 401)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)

  p_mc   <- mcnemar_pval(b, c_)

  data.frame(
    n = n, theta1 = theta1, theta2 = theta2,
    x1 = x1, x2 = x2, b = b, c = c_,
    ev_obs = ev_obs, k_star = opt$k_star,
    fbst   = ifelse(ev_obs <= opt$k_star, "reject", "do not reject"),
    mcnemar_p = p_mc,
    mcnemar   = ifelse(p_mc <= alpha_level, "reject", "do not reject")
  )
}

# Escenarios: efecto nulo, pequeño, moderado, grande
escenarios <- list(
  list(n=50,  t1=0.50, t2=0.50, lab="Efecto nulo (θ1=θ2=0.50)"),
  list(n=50,  t1=0.55, t2=0.45, lab="Efecto pequeño (Δ=0.10, n=50)"),
  list(n=25,  t1=0.60, t2=0.40, lab="Efecto moderado (Δ=0.20, n=25)"),
  list(n=25,  t1=0.70, t2=0.30, lab="Efecto grande (Δ=0.40, n=25)"),
  list(n=15,  t1=0.60, t2=0.40, lab="Efecto moderado n pequeño (n=15)"),
  list(n=100, t1=0.52, t2=0.48, lab="Efecto mínimo n grande (Δ=0.04, n=100)"),
  list(n=400, t1=0.55, t2=0.45, lab="Efecto pequeño n grande (Δ=0.10, n=400)")
)

cmp <- do.call(rbind, lapply(seq_along(escenarios), function(i) {
  e <- escenarios[[i]]
  cat(sprintf("  %s...\n", e$lab))
  r <- compare_tests(e$n, e$t1, e$t2, a0_kl, a1_kl, a2_kl, M, seed = i * 10)
  r$escenario <- e$lab
  r
}))

cat("\nResultados FBST vs McNemar:\n")
print(cmp[, c("escenario", "ev_obs", "k_star", "fbst", "mcnemar_p", "mcnemar")],
      digits = 3, row.names = FALSE)

# Figura: ev_obs vs k* con color=decisión FBST, forma=decisión McNemar
cmp$acuerdo <- ifelse(cmp$fbst == cmp$mcnemar, "acuerdo", "desacuerdo")
p_mc <- ggplot(cmp, aes(x = ev_obs, y = k_star)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = fbst, shape = mcnemar), size = 4) +
  geom_text(aes(label = sprintf("n=%d\nΔ=%.2f", n, abs(theta1 - theta2))),
            hjust = -0.15, size = 2.8) +
  scale_color_manual(values = c("reject" = "tomato", "do not reject" = "steelblue"),
                     name = "FBST") +
  scale_shape_manual(values = c("reject" = 17, "do not reject" = 16),
                     name = "McNemar") +
  annotate("text", x = 0.25, y = 0.72, label = "ev < k*: reject →",
           color = "tomato", size = 3) +
  annotate("text", x = 0.65, y = 0.25, label = "← ev > k*: do not reject",
           color = "steelblue", size = 3) +
  labs(x = "observed ev", y = "optimal k* (posterior-based)",
       title = "FBST vs McNemar: when do they agree?",
       subtitle = "Points above the diagonal: reject. Below: do not reject.") +
  theme_bw() + theme(legend.position = "right")
ggsave("Figures/mcnemar_comparison.png", p_mc, width = 8, height = 5, dpi = 150)
cat("→ Figures/mcnemar_comparison.png\n")

# ===========================================================================
# 3. k* adaptativo vs tamaño muestral (figura limpia para el paper)
#    Usa proporciones exactas sin redondeo: varía n y calcula x=round(p*n),
#    pero reporta p̂ real para que el lector vea las fluctuaciones.
# ===========================================================================
cat("\n=== 3. k* adaptativo vs n ===\n")

run_n_series <- function(p1, p2, n_vec, a0, a1, a2, M, seed = 42, label) {
  lapply(n_vec, function(n) {
    x1 <- round(p1 * n); x2 <- round(p2 * n)
    set.seed(seed)
    ev_obs <- ev_quad_from_data(n, n, x1, x2, a0, a1, a2)
    ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M, 401)
    ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M, 401)
    opt    <- find_kstar(ev_H, ev_A, 1, 1)
    data.frame(n = n, ev_obs = ev_obs, k_star = opt$k_star,
               ab = opt$alpha + opt$beta,
               decision = ifelse(ev_obs <= opt$k_star, "reject H", "do not reject"),
               label = label)
  }) |> do.call(what = rbind)
}

cat("  Serie 1: Δp = 0.15...\n")
s1 <- run_n_series(0.55, 0.40, c(10,20,30,40,50,75,100), a0_kl, a1_kl, a2_kl, M,
                   label = "Δp = 0.15 (p1=0.55, p2=0.40)")
cat("  Serie 2: Δp = 0.08...\n")
s2 <- run_n_series(0.50, 0.42, c(20,30,50,75,100,150,200), a0_kl, a1_kl, a2_kl, M,
                   label = "Δp = 0.08 (p1=0.50, p2=0.42)")
df_n <- rbind(s1, s2)

p_n <- ggplot(df_n, aes(x = n)) +
  geom_ribbon(aes(ymin = pmin(ev_obs, k_star), ymax = pmax(ev_obs, k_star),
                  fill = decision), alpha = 0.15) +
  geom_line(aes(y = ev_obs, color = "observed ev"), linewidth = 1) +
  geom_line(aes(y = k_star, color = "optimal k*"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = ev_obs, color = "observed ev", shape = decision), size = 3) +
  scale_color_manual(values = c("observed ev" = "steelblue", "optimal k*" = "tomato")) +
  scale_fill_manual(values  = c("reject H" = "tomato", "do not reject" = "steelblue")) +
  scale_shape_manual(values = c("reject H" = 17, "do not reject" = 16)) +
  facet_wrap(~label, scales = "free_x") +
  labs(x = "n  (n1 = n2)", y = "value",
       title = "k* adaptive to sample size (posterior formulation)",
       subtitle = "The test rejects H when ev_obs falls below k*",
       color = NULL, shape = "Decision", fill = "Decision") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Figures/decision_vs_n_clean.png", p_n, width = 10, height = 5, dpi = 150)
cat("→ Figures/decision_vs_n_clean.png\n")

# ===========================================================================
# 4. Sensibilidad a hiperparámetros: KL-óptimo vs informativo (α=10)
#    Para los 4 grupos del THKS muestra cómo cambian ev_obs, k*, decisión.
# ===========================================================================
cat("\n=== 4. Sensibilidad a hiperparámetros ===\n")

cat("  Prior no informativa (KL→U)...\n")
tab_kl <- do.call(rbind, lapply(names(grupos), function(k) {
  r <- run_group(grupos[[k]], a0_kl, a1_kl, a2_kl, M)
  r$prior <- "Non-informative (KL)"
  r
}))

cat("  Prior informativa simétrica (μ=0.5, N=50)...\n")
tab_i2 <- do.call(rbind, lapply(names(grupos), function(k) {
  r <- run_group(grupos[[k]], a0_inf, a1_inf, a2_inf, M)
  r$prior <- "Informative (μ=0.5)"
  r
}))

cat("  Prior con conflicto (μ=0.1, N=50)...\n")
tab_i3 <- do.call(rbind, lapply(names(grupos), function(k) {
  r <- run_group(grupos[[k]], a0_conf, a1_conf, a2_conf, M)
  r$prior <- "Conflict (μ=0.1)"
  r
}))

sens <- rbind(tab_kl, tab_i2, tab_i3)
sens$prior <- factor(sens$prior,
                     levels = c("Non-informative (KL)",
                                "Informative (μ=0.5)",
                                "Conflict (μ=0.1)"))
cat("\nComparación k* y ev_obs por prior:\n")
print(sens[, c("group", "prior", "ev_obs", "k_star", "decision")],
      digits = 4, row.names = FALSE)

p_sens <- ggplot(sens, aes(x = prior, y = k_star, group = group, color = group)) +
  geom_line(linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(shape = decision), size = 4) +
  geom_hline(data = sens[sens$prior == "Non-informative (KL)", ],
             aes(yintercept = ev_obs, color = group),
             linetype = "dotted", linewidth = 0.5) +
  scale_shape_manual(values = c("reject H" = 17, "do not reject" = 16)) +
  labs(x = "Prior", y = "optimal k*",
       title = "Sensitivity of k* to the hyperparameters",
       subtitle = "Dotted lines = ev_obs (does not change with the prior)",
       color = "Group", shape = "Decision") +
  theme_bw() + theme(legend.position = "right")
ggsave("Figures/sensitivity_hyperparams.png", p_sens, width = 7, height = 5, dpi = 150)
cat("→ Figures/sensitivity_hyperparams.png\n")

# ===========================================================================
# 5. Ejemplo donde el test NO rechaza H (efecto nulo o mínimo)
#    Muestra las curvas de error y la posición de ev_obs respecto a k*.
# ===========================================================================
cat("\n=== 5. Ejemplo de no rechazo ===\n")

# Caso 1: efecto real = 0
n_nr <- 50; x1_nr <- 25; x2_nr <- 25  # θ1=θ2=0.5 exacto
# Caso 2: efecto mínimo con n chico
n_nr2 <- 20; x1_nr2 <- 11; x2_nr2 <- 10

nr_cases <- list(
  list(n=n_nr,  x1=x1_nr,  x2=x2_nr,
       lab = sprintf("No effect: n=%d, x1=%d, x2=%d (p̂1=p̂2=0.50)", n_nr, x1_nr, x2_nr)),
  list(n=n_nr2, x1=x1_nr2, x2=x2_nr2,
       lab = sprintf("Minimal effect: n=%d, x1=%d, x2=%d", n_nr2, x1_nr2, x2_nr2))
)

plots_nr <- lapply(nr_cases, function(case) {
  n <- case$n; x1 <- case$x1; x2 <- case$x2
  set.seed(99)
  ev_obs <- ev_quad_from_data(n, n, x1, x2, a0_kl, a1_kl, a2_kl)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0_kl, a1_kl, a2_kl, M, 401)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0_kl, a1_kl, a2_kl, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  cat(sprintf("  %s\n  ev_obs=%.4f  k*=%.4f  → %s\n\n",
              case$lab, ev_obs, opt$k_star,
              ifelse(ev_obs <= opt$k_star, "RECHAZA H", "NO RECHAZA H")))

  k_grid <- seq(0, 1, length.out = 401)
  curves <- error_curves(ev_H, ev_A, k_grid, 1, 1)
  long   <- pivot_longer(curves, c("alpha", "beta", "sum"),
                         names_to = "m", values_to = "v")
  long$m <- factor(long$m, levels = c("alpha","beta","sum"),
                   labels = c("α (type I error)", "β (type II error)", "α+β"))
  ggplot(long, aes(k, v, color = m, linetype = m)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = opt$k_star, linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = ev_obs,     linetype = "solid",  color = "purple",
               linewidth = 1) +
    annotate("text", x = opt$k_star, y = 0.92,
             label = sprintf("k*=%.3f", opt$k_star), hjust = -0.1, size = 3) +
    annotate("text", x = ev_obs, y = 0.82,
             label = sprintf("ev=%.3f", ev_obs), hjust = 1.1, size = 3,
             color = "purple") +
    annotate("label", x = 0.5, y = 0.5,
             label = ifelse(ev_obs <= opt$k_star, "REJECT H", "DO NOT REJECT H"),
             color = ifelse(ev_obs <= opt$k_star, "tomato", "steelblue"),
             size = 4, fontface = "bold") +
    scale_color_manual(values = c("steelblue","tomato","black")) +
    scale_linetype_manual(values = c("solid","solid","dashed")) +
    labs(x = "k", y = "averaged error", title = case$lab,
         color = NULL, linetype = NULL) +
    theme_bw() + theme(legend.position = "top")
})

# Combina los dos paneles
library(gridExtra)
g <- arrangeGrob(plots_nr[[1]], plots_nr[[2]], ncol = 2)
ggsave("Figures/nonrejection_example.png", g, width = 12, height = 5, dpi = 150)
cat("→ Figures/nonrejection_example.png\n")

# ===========================================================================
# Resumen final
# ===========================================================================
cat("\n========== Resumen ==========\n")
cat("output/kstar_informative_table.tex  — Tab3/4 corregidas (prior α=10)\n")
cat("Figures/mcnemar_comparison.png      — FBST vs McNemar\n")
cat("Figures/decision_vs_n_clean.png     — k* adaptativo vs n\n")
cat("Figures/sensitivity_hyperparams.png — sensibilidad a hiperparámetros\n")
cat("Figures/nonrejection_example.png    — ejemplo de no rechazo\n")
