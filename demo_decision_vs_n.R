# demo_decision_vs_n.R
# Muestra cómo la decisión (rechaza/no rechaza H) depende del tamaño muestral
# cuando el efecto real es pequeño.
#
# Idea: fijar proporciones observadas p1, p2 (efecto pequeño pero real),
# escalar n, y ver dónde ev_obs cruza k*.
#
# Para n chico: ev_obs > k*  →  no rechaza  (poca evidencia)
# Para n grande: ev_obs < k* →  rechaza     (efecto acumulado)

library(Rcpp)
library(ggplot2)
library(dplyr)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")

source("priors_config.R")
a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]

# ---------------------------------------------------------------------------
# Escenario 1: efecto moderado-chico — p1=0.55, p2=0.40 (15 pp de diferencia)
# ---------------------------------------------------------------------------
p1 <- 0.55; p2 <- 0.40
n_grid <- c(10, 15, 20, 25, 30, 40, 50, 75, 100)
M      <- 1500
seed   <- 42

cat(sprintf("\n=== p1=%.2f, p2=%.2f (Δ=%.2f) ===\n", p1, p2, p1 - p2))
cat(sprintf("%-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-12s\n",
            "n", "x1", "x2", "ev_obs", "k*", "α+β", "decisión"))

rows1 <- lapply(n_grid, function(n) {
  x1 <- round(p1 * n); x2 <- round(p2 * n)
  set.seed(seed)
  ev_obs <- ev_quad_from_data(n, n, x1, x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M, 401)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  dec    <- ifelse(ev_obs <= opt$k_star, "REJECT H", "do not reject")
  cat(sprintf("%-6d  %-6d  %-6d  %-6.4f  %-6.4f  %-6.4f  %s\n",
              n, x1, x2, ev_obs, opt$k_star, opt$alpha + opt$beta, dec))
  data.frame(n=n, x1=x1, x2=x2, ev_obs=ev_obs,
             k_star=opt$k_star, ab=opt$alpha+opt$beta,
             decision=dec, scenario="p1=0.55, p2=0.40")
})
df1 <- do.call(rbind, rows1)

# ---------------------------------------------------------------------------
# Escenario 2: efecto muy chico — p1=0.50, p2=0.42 (8 pp de diferencia)
# ---------------------------------------------------------------------------
p1 <- 0.50; p2 <- 0.42
n_grid2 <- c(20, 30, 50, 75, 100, 150, 200)

cat(sprintf("\n=== p1=%.2f, p2=%.2f (Δ=%.2f) ===\n", p1, p2, p1 - p2))
cat(sprintf("%-6s  %-6s  %-6s  %-6s  %-6s  %-6s  %-12s\n",
            "n", "x1", "x2", "ev_obs", "k*", "α+β", "decisión"))

rows2 <- lapply(n_grid2, function(n) {
  x1 <- round(p1 * n); x2 <- round(p2 * n)
  set.seed(seed)
  ev_obs <- ev_quad_from_data(n, n, x1, x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M, 401)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  dec    <- ifelse(ev_obs <= opt$k_star, "REJECT H", "do not reject")
  cat(sprintf("%-6d  %-6d  %-6d  %-6.4f  %-6.4f  %-6.4f  %s\n",
              n, x1, x2, ev_obs, opt$k_star, opt$alpha + opt$beta, dec))
  data.frame(n=n, x1=x1, x2=x2, ev_obs=ev_obs,
             k_star=opt$k_star, ab=opt$alpha+opt$beta,
             decision=dec, scenario="p1=0.50, p2=0.42")
})
df2 <- do.call(rbind, rows2)

# ---------------------------------------------------------------------------
# Gráfico: ev_obs y k* vs n para ambos escenarios
# ---------------------------------------------------------------------------
df_all <- rbind(df1, df2)

p <- ggplot(df_all, aes(x = n)) +
  geom_line(aes(y = ev_obs, color = "observed ev"), linewidth = 1) +
  geom_point(aes(y = ev_obs, color = "observed ev",
                 shape = decision), size = 3) +
  geom_line(aes(y = k_star, color = "optimal k*"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = k_star, color = "optimal k*"), size = 2) +
  scale_color_manual(values = c("observed ev" = "steelblue",
                                "optimal k*"  = "tomato")) +
  scale_shape_manual(values = c("REJECT H" = 17, "do not reject" = 16)) +
  facet_wrap(~scenario, scales = "free_x") +
  labs(x = "n (n1 = n2 = n)", y = "value",
       title = "Decision as a function of sample size (posterior formulation)",
       subtitle = "Reject H when ev_obs crosses below k*",
       color = NULL, shape = "Decision") +
  theme_bw() + theme(legend.position = "bottom")

ggsave("Figures/decision_vs_n.png", p, width = 10, height = 5, dpi = 150)
cat("\n→ Figures/decision_vs_n.png\n")

# ---------------------------------------------------------------------------
# El caso real: grupo yn del THKS (frontera natural del experimento)
# ---------------------------------------------------------------------------
cat("\n=== Caso real: grupo yn (CC sin TV) ===\n")
cat("n1=27, n2=27 — ev_obs=0.397, k*=0.412 → rechaza por margen 0.015\n")
cat("Con n=20 (mismas proporciones), probablemente no rechazaría.\n")
cat("Esto ilustra que el yn es un caso genuinamente borderline.\n")
