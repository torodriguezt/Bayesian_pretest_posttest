# make_figures.R
# Regenera todas las figuras del artículo en Figures/.
#
#   - error_curves_prior_n{N}.png   → α(k), β(k), α+β bajo prior, n1=n2=N
#                                      (réplica de figs 9-12 del paper)
#   - error_curves_post_{grupo}.png → α(k), β(k), α+β bajo posterior,
#                                      grupos del THKS
#   - posterior_surface_{grupo}.png → superficie posterior bivariada
#
# Tiempo aprox: ~5-10 min con M = 1000.

library(ALA)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(tidyr)
library(lattice)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")

dir.create("Figures", showWarnings = FALSE)

source("priors_config.R")
a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
plot_error_curves <- function(ev_H, ev_A, title, save_path,
                              k_grid = seq(0, 1, length.out = 401)) {
  curves <- error_curves(ev_H, ev_A, k_grid, 1, 1)
  opt <- find_kstar(ev_H, ev_A, 1, 1)
  long <- pivot_longer(curves, c("alpha", "beta", "sum"),
                       names_to = "metrica", values_to = "valor")
  long$metrica <- factor(long$metrica,
                         levels = c("alpha", "beta", "sum"),
                         labels = c("alpha", "beta", "alpha+beta"))
  p <- ggplot(long, aes(k, valor, color = metrica, linetype = metrica)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = opt$k_star, linetype = "dashed",
               color = "grey40") +
    annotate("text", x = opt$k_star, y = 0.05,
             label = sprintf("k* = %.3f", opt$k_star),
             hjust = -0.1, size = 3.5) +
    scale_color_manual(values = c("steelblue", "tomato", "black")) +
    scale_linetype_manual(values = c("solid", "solid", "dashed")) +
    labs(x = "k", y = "Averaged error", title = title,
         color = NULL, linetype = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave(save_path, p, width = 6, height = 4, dpi = 150)
  cat("  →", save_path, sprintf("(k*=%.4f, α=%.4f, β=%.4f)\n",
                                opt$k_star, opt$alpha, opt$beta))
  invisible(opt)
}

# ---------------------------------------------------------------------------
# 1) Curvas de error PRIOR-based — réplica de figs 9-12 del paper
#     n1 = n2 ∈ {25, 50, 75, 100}
# ---------------------------------------------------------------------------
cat("\n=== Curvas prior-based (réplica figs 9-12) ===\n")
ns_prior <- c(25, 50, 75, 100)
M_prior  <- 1000

for (N in ns_prior) {
  cat(sprintf("n1 = n2 = %d (M = %d)...\n", N, M_prior))
  set.seed(123 + N)
  ev_H <- simulate_evs_H(N, N, a0, a1, a2, M_prior, 401)
  ev_A <- simulate_evs_A(N, N, a0, a1, a2, M_prior, 401)
  plot_error_curves(
    ev_H, ev_A,
    title = sprintf("Prior-based: n1 = n2 = %d (M = %d)", N, M_prior),
    save_path = sprintf("Figures/error_curves_prior_n%d.png", N))
}

# ---------------------------------------------------------------------------
# 2) Curvas de error POSTERIOR-based — los 4 grupos del THKS
# ---------------------------------------------------------------------------
cat("\n=== Curvas posterior-based (THKS) ===\n")

datos1 <- tvsfp

# Grupos definidos solo por condicion de tratamiento (sin filtro por escuela)
prep_group <- function(d, sb, tv) {
  wide <- d %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school.based == sb, tv.based == tv) %>%
    select(id, stage, binTHKS) %>%
    pivot_wider(names_from = stage, values_from = binTHKS) %>%
    filter(!is.na(pre), !is.na(post))
  list(X = c(sum(wide$pre), sum(wide$post)),
       n = c(nrow(wide), nrow(wide)))
}

grupos <- list(
  yy = list(g = prep_group(datos1, "yes", "yes"), lab = "CC + TV"),
  yn = list(g = prep_group(datos1, "yes", "no"),  lab = "CC, no TV"),
  ny = list(g = prep_group(datos1, "no",  "yes"), lab = "no CC, TV"),
  nn = list(g = prep_group(datos1, "no",  "no"),  lab = "no CC, no TV")
)

M_post <- 1000

for (k in names(grupos)) {
  g <- grupos[[k]]$g; lab <- grupos[[k]]$lab
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  cat(sprintf("%s — n1=%d n2=%d x1=%d x2=%d (M = %d)...\n",
              lab, n1, n2, x1, x2, M_post))
  set.seed(7)
  ev_H <- simulate_evs_H_post(n1, n2, x1, x2, a0, a1, a2, M_post, 401)
  ev_A <- simulate_evs_A_post(n1, n2, x1, x2, a0, a1, a2, M_post, 401)
  plot_error_curves(
    ev_H, ev_A,
    title = sprintf("Posterior-based: %s (n1=%d, n2=%d)", lab, n1, n2),
    save_path = sprintf("Figures/error_curves_post_%s.png", k))
}

# ---------------------------------------------------------------------------
# 3) Superficies posteriores bivariadas — los 4 grupos del THKS
# ---------------------------------------------------------------------------
cat("\n=== Posterior surfaces (THKS) ===\n")
for (k in names(grupos)) {
  g <- grupos[[k]]$g; lab <- grupos[[k]]$lab
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  consts <- bb_constants(n1, n2, x1, x2, a0, a1, a2)
  xs <- seq(0.01, 0.99, length.out = 80)
  ys <- seq(0.01, 0.99, length.out = 80)
  z  <- densBB_grid(xs, ys, consts)
  df <- expand.grid(theta1 = xs, theta2 = ys)
  df$density <- as.vector(z)
  p <- wireframe(
    density ~ theta1 * theta2, data = df,
    xlab = expression(theta[1]),
    ylab = expression(theta[2]),
    zlab = NULL,
    main = sprintf("%s (n1=%d, n2=%d, x1=%d, x2=%d)",
                   lab, n1, n2, x1, x2),
    scales = list(arrows = FALSE, cex = 0.7),
    drape = TRUE,
    col.regions = "gold",
    colorkey = FALSE,
    screen = list(z = -30, x = -60),
    par.settings = list(axis.line = list(col = "transparent"))
  )
  png(sprintf("Figures/posterior_surface_%s.png", k),
      width = 700, height = 600, res = 100)
  print(p)
  dev.off()
  cat(sprintf("  -> Figures/posterior_surface_%s.png\n", k))
}

cat("\n=== Listo ===\n")
cat("Figuras generadas en Figures/:\n")
print(list.files("Figures"))
