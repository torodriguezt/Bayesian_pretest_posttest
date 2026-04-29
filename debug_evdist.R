# debug_evdist.R
# Visualizes the distribution of e-values under H and A for two priors
# This explains why k* changes so much between priors even with large n

library(Rcpp)
library(ggplot2)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")

# Hyperparameters (from priors_config.R, fitted by KL minimization)
source("priors_config.R")
a0_kl   <- prior_NI["a0"];   a1_kl   <- prior_NI["a1"];   a2_kl   <- prior_NI["a2"]
a0_inf  <- prior_INF["a0"];  a1_inf  <- prior_INF["a1"];  a2_inf  <- prior_INF["a2"]
a0_conf <- prior_CONF["a0"]; a1_conf <- prior_CONF["a1"]; a2_conf <- prior_CONF["a2"]

# Scenario: n=100, large effect (θ1=0.5, θ2=0.8)
n <- 100; theta1 <- 0.5; theta2 <- 0.8; d <- 0.3

# Simulate data
set.seed(123)
x1 <- rbinom(1, n, theta1)
x2 <- rbinom(1, n, theta2)

cat(sprintf("\nData: n=%d, x1=%d, x2=%d (theta1=%.1f, theta2=%.1f)\n\n", n, x1, x2, theta1, theta2))

# ====== KL-OPTIMAL ======
cat("KL-optimal prior:\n")
ev_H_kl <- simulate_evs_H_post(n, n, x1, x2, a0_kl, a1_kl, a2_kl, M=1000)
ev_A_kl <- simulate_evs_A_post(n, n, x1, x2, a0_kl, a1_kl, a2_kl, M=1000)
opt_kl <- find_kstar(ev_H_kl, ev_A_kl, 1, 1)
cat(sprintf("  k* = %.4f, alpha = %.4f, beta = %.4f\n", opt_kl$k_star, opt_kl$alpha, opt_kl$beta))
cat(sprintf("  mean(ev_H) = %.4f, mean(ev_A) = %.4f\n\n", mean(ev_H_kl), mean(ev_A_kl)))

# ====== INFORMATIVE (symmetric) ======
cat("Informative symmetric prior (mu=0.5, N=50):\n")
ev_H_inf <- simulate_evs_H_post(n, n, x1, x2, a0_inf, a1_inf, a2_inf, M=1000)
ev_A_inf <- simulate_evs_A_post(n, n, x1, x2, a0_inf, a1_inf, a2_inf, M=1000)
opt_inf <- find_kstar(ev_H_inf, ev_A_inf, 1, 1)
cat(sprintf("  k* = %.4f, alpha = %.4f, beta = %.4f\n", opt_inf$k_star, opt_inf$alpha, opt_inf$beta))
cat(sprintf("  mean(ev_H) = %.4f, mean(ev_A) = %.4f\n\n", mean(ev_H_inf), mean(ev_A_inf)))

# ====== INFORMATIVE WITH CONFLICT ======
cat("Informative prior with conflict (mu=0.1, N=50):\n")
ev_H_conf <- simulate_evs_H_post(n, n, x1, x2, a0_conf, a1_conf, a2_conf, M=1000)
ev_A_conf <- simulate_evs_A_post(n, n, x1, x2, a0_conf, a1_conf, a2_conf, M=1000)
opt_conf <- find_kstar(ev_H_conf, ev_A_conf, 1, 1)
cat(sprintf("  k* = %.4f, alpha = %.4f, beta = %.4f\n", opt_conf$k_star, opt_conf$alpha, opt_conf$beta))
cat(sprintf("  mean(ev_H) = %.4f, mean(ev_A) = %.4f\n\n", mean(ev_H_conf), mean(ev_A_conf)))

# Plot
df <- rbind(
  data.frame(ev=ev_H_kl,   hyp="Under H0", prior="KL-optimal"),
  data.frame(ev=ev_A_kl,   hyp="Under A",  prior="KL-optimal"),
  data.frame(ev=ev_H_inf,  hyp="Under H0", prior="Informative"),
  data.frame(ev=ev_A_inf,  hyp="Under A",  prior="Informative"),
  data.frame(ev=ev_H_conf, hyp="Under H0", prior="Conflict"),
  data.frame(ev=ev_A_conf, hyp="Under A",  prior="Conflict")
)
df$prior <- factor(df$prior, levels=c("KL-optimal","Informative","Conflict"))

vlines <- data.frame(
  prior  = factor(c("KL-optimal","Informative","Conflict"),
                  levels=c("KL-optimal","Informative","Conflict")),
  k_star = c(opt_kl$k_star, opt_inf$k_star, opt_conf$k_star)
)

p <- ggplot(df, aes(x=ev, fill=hyp)) +
  geom_histogram(bins=40, alpha=0.6, position="identity") +
  geom_vline(data=vlines, aes(xintercept=k_star),
             linetype="dashed", color="red", size=1) +
  facet_wrap(~prior) +
  labs(x="E-value", y="Frequency", title="E-value distributions (red dashed = k*)", fill="Hypothesis") +
  theme_bw()

ggsave("Figures/ev_distribution_comparison.png", p, width=10, height=5, dpi=150)
cat("Plot saved: Figures/ev_distribution_comparison.png\n")
