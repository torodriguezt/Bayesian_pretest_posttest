# make_simulation_figures.R
# Generates all simulation study figures (article Figs 3-8).
# Replaces Stan-based figures with the C++ engine.
#
# Outputs (in Figures/):
#   prior_noinf.png, post_noinf_n20_t0.1.png
#   prior_inf_Confl.png, post_inf_n20_t0.5.png, post_inf_n20_t0.1_Confl.png
#   mean_noinf_100_stan.png,      mode_noinf_100_stan.png
#   mean_inf_100_NoConfl_stan.png, mode_inf_100_NoConfl_stan.png
#   mean_inf_100_Confl_stan.png,  mode_inf_100_Confl_stan.png
#
# Runtime: ~5-10 min (R_SIM=100, n=2:60, grid=60x60).

library(Rcpp)
library(ggplot2)
library(dplyr)
library(lattice)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")
dir.create("Figures", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Hyperparameters
# ---------------------------------------------------------------------------
a0_kl   <- 0.8373879; a1_kl   <- 0.8410984; a2_kl   <- 0.8053298
a0_inf  <- 10;        a1_inf  <- 10;        a2_inf  <- 10
a0_conf <- 90;        a1_conf <- 10;        a2_conf <- 10

NGRID <- 60    # grid for density/moments (60x60 = fast)
R_SIM <- 100   # Monte Carlo repetitions per n
N_MAX <- 60    # maximum sample size

# ---------------------------------------------------------------------------
# Helper: marginal posterior mean and mode from grid
# ---------------------------------------------------------------------------
post_moments <- function(n, x1, x2, a0, a1, a2, ngrid = NGRID) {
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  xs <- seq(0.001, 0.999, length.out = ngrid)
  z  <- densBB_grid(xs, xs, consts)
  z[!is.finite(z)] <- 0
  f1 <- rowSums(z); f2 <- colSums(z)
  list(
    mean1 = sum(xs * f1) / sum(f1),
    mean2 = sum(xs * f2) / sum(f2),
    mode1 = xs[which.max(f1)],
    mode2 = xs[which.max(f2)]
  )
}

# ---------------------------------------------------------------------------
# Helper: simulation — mean and mode of posterior vs n
# ---------------------------------------------------------------------------
run_simulation <- function(a0, a1, a2, theta_true,
                           n_max = N_MAX, R = R_SIM, seed = 42) {
  set.seed(seed)
  n_vals <- 2:n_max
  cat(sprintf("  n = 2:%d, R=%d reps each...\n", n_max, R))
  do.call(rbind, lapply(n_vals, function(n) {
    res <- replicate(R, {
      x1 <- rbinom(1, n, theta_true)
      x2 <- rbinom(1, n, theta_true)
      m  <- post_moments(n, x1, x2, a0, a1, a2)
      c(m$mean1, m$mean2, m$mode1, m$mode2)
    })
    data.frame(
      n      = n,
      mE1    = mean(res[1, ]),  loE1 = quantile(res[1, ], .025),  hiE1 = quantile(res[1, ], .975),
      mE2    = mean(res[2, ]),  loE2 = quantile(res[2, ], .025),  hiE2 = quantile(res[2, ], .975),
      mMo1   = mean(res[3, ]),  loMo1= quantile(res[3, ], .025),  hiMo1= quantile(res[3, ], .975),
      mMo2   = mean(res[4, ]),  loMo2= quantile(res[4, ], .025),  hiMo2= quantile(res[4, ], .975)
    )
  }))
}

# ---------------------------------------------------------------------------
# Helper: ggplot2 estimation figure (mean or mode vs n)
# ---------------------------------------------------------------------------
plot_estim <- function(df, theta_true, type = c("Mean", "Mode"),
                       title, path) {
  type <- match.arg(type)
  if (type == "Mean") {
    d <- rbind(
      data.frame(n=df$n, est=df$mE1,  lo=df$loE1,  hi=df$hiE1,
                 param = paste0("E(\u03b8\u2081)")),
      data.frame(n=df$n, est=df$mE2,  lo=df$loE2,  hi=df$hiE2,
                 param = paste0("E(\u03b8\u2082)"))
    )
  } else {
    d <- rbind(
      data.frame(n=df$n, est=df$mMo1, lo=df$loMo1, hi=df$hiMo1,
                 param = paste0("Mode(\u03b8\u2081)")),
      data.frame(n=df$n, est=df$mMo2, lo=df$loMo2, hi=df$hiMo2,
                 param = paste0("Mode(\u03b8\u2082)"))
    )
  }
  true_lbl <- sprintf("\u03b8\u2081=\u03b8\u2082=%.1f", theta_true)
  p <- ggplot(d, aes(x = n)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = param), alpha = 0.2) +
    geom_line(aes(y = est, color = param), linewidth = 0.8) +
    geom_hline(aes(yintercept = theta_true, linetype = true_lbl),
               color = "black", linewidth = 0.6) +
    scale_color_manual(name = "Estimated",
                       values = c("#888888", "#00bcd4")) +
    scale_fill_manual(name  = "Estimated",
                      values = c("#888888", "#00bcd4")) +
    scale_linetype_manual(name = "Real", values = "dashed") +
    labs(x = "n", y = type, title = title) +
    theme_bw() +
    theme(legend.position = "right")
  ggsave(path, p, width = 7, height = 5, dpi = 120)
  cat("  ->", path, "\n")
}

# ---------------------------------------------------------------------------
# Helper: bivariate density persp plot → PNG
# ---------------------------------------------------------------------------
plot_density_png <- function(a0, a1, a2, n, x1, x2, path, main,
                              ngrid = 80) {
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  xs <- seq(0.001, 0.999, length.out = ngrid)
  z  <- densBB_grid(xs, xs, consts)
  z[!is.finite(z)] <- 0
  cap <- quantile(z, 0.99)
  z[z > cap] <- cap
  df <- expand.grid(theta1 = xs, theta2 = xs)
  df$density <- as.vector(z)
  p <- wireframe(
    density ~ theta1 * theta2, data = df,
    xlab = expression(theta[1]),
    ylab = expression(theta[2]),
    zlab = NULL,
    main = main,
    scales = list(arrows = FALSE, cex = 0.7),
    drape = TRUE,
    col.regions = "lightblue",
    colorkey = FALSE,
    screen = list(z = -30, x = -60),
    par.settings = list(axis.line = list(col = "transparent"))
  )
  png(path, width = 700, height = 600, res = 100)
  print(p)
  dev.off()
  cat("  ->", path, "\n")
}

# ---------------------------------------------------------------------------
# 1) Prior and posterior density plots
# ---------------------------------------------------------------------------
cat("\n=== Prior densities ===\n")

plot_density_png(a0_kl,   a1_kl,   a2_kl,   0, 0, 0,
  "Figures/prior_noinf.png",
  sprintf("Non-informative prior (\u03b1=%.2f)", a0_kl))

# Informative prior (α=10) — reused for both conflict / no-conflict figures
plot_density_png(a0_inf, a1_inf, a2_inf, 0, 0, 0,
  "Figures/prior_inf_Confl.png",
  "Informative prior (\u03b1\u2080=\u03b1\u2081=\u03b1\u2082=10)")

cat("\n=== Posterior densities (n=20) ===\n")

x01 <- round(0.1 * 20)   # = 2
x05 <- round(0.5 * 20)   # = 10

plot_density_png(a0_kl,   a1_kl,   a2_kl,   20, x01, x01,
  "Figures/post_noinf_n20_t0.1.png",
  sprintf("Posterior: non-inf., n=20, x\u2081=x\u2082=%d", x01))

plot_density_png(a0_inf, a1_inf, a2_inf, 20, x05, x05,
  "Figures/post_inf_n20_t0.5.png",
  sprintf("Posterior: inf., n=20, x\u2081=x\u2082=%d (no conflict)", x05))

plot_density_png(a0_conf, a1_conf, a2_conf, 20, x01, x01,
  "Figures/post_inf_n20_t0.1_Confl.png",
  sprintf("Posterior: conflict prior, n=20, x\u2081=x\u2082=%d", x01))

# ---------------------------------------------------------------------------
# 2) Simulation: mode and mean vs n — three scenarios
# ---------------------------------------------------------------------------
scenarios <- list(
  list(a0=a0_kl,   a1=a1_kl,   a2=a2_kl,   theta=0.1, seed=42,
       label="Non-informative (KL-optimal), \u03b8=0.1",
       mean_out = "Figures/mean_noinf_100_stan.png",
       mode_out = "Figures/mode_noinf_100_stan.png"),
  list(a0=a0_inf,  a1=a1_inf,  a2=a2_inf,  theta=0.5, seed=43,
       label="Informative (no conflict), \u03b8=0.5",
       mean_out = "Figures/mean_inf_100_NoConfl_stan.png",
       mode_out = "Figures/mode_inf_100_NoConfl_stan.png"),
  list(a0=a0_conf, a1=a1_conf, a2=a2_conf, theta=0.1, seed=44,
       label="Informative (conflict, \u03b1\u2080=90), \u03b8=0.1",
       mean_out = "Figures/mean_inf_100_Confl_stan.png",
       mode_out = "Figures/mode_inf_100_Confl_stan.png")
)

for (sc in scenarios) {
  cat(sprintf("\n=== %s ===\n", sc$label))
  df <- run_simulation(sc$a0, sc$a1, sc$a2, sc$theta, seed = sc$seed)
  plot_estim(df, sc$theta, "Mean", sc$label, sc$mean_out)
  plot_estim(df, sc$theta, "Mode", sc$label, sc$mode_out)
}

cat("\n=== Done — all simulation figures in Figures/ ===\n")
print(list.files("Figures", pattern = "\\.png$"))
