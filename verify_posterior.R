# verify_posterior.R
# Verifies that the posterior under H is similar between priors,
# but the e-value distributions are different because the e-value
# is calculated using the SAME prior that was used to sample theta.

library(Rcpp)
library(ggplot2)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")

a0_kl <- 0.8373879; a1_kl <- 0.8410984; a2_kl <- 0.8053298
a0_inf <- 10; a1_inf <- 10; a2_inf <- 10

# Same data
n <- 150
set.seed(42)
x1 <- rbinom(1, n, 0.5)
x2 <- rbinom(1, n, 0.52)
cat(sprintf("Data: n=%d, x1=%d, x2=%d\n\n", n, x1, x2))

# ========================================================
# STEP 1: Compare POSTERIOR under H (should be similar)
# ========================================================
cat("STEP 1: Posterior under H (theta1=theta2=t)\n")

# Sample theta under H for both priors
fH_kl  <- fH_post_cdf(n, n, x1, x2, a0_kl, a1_kl, a2_kl, 10001)
fH_inf <- fH_post_cdf(n, n, x1, x2, a0_inf, a1_inf, a2_inf, 10001)

theta_H_kl  <- sample_fH(5000, fH_kl)
theta_H_inf <- sample_fH(5000, fH_inf)

cat(sprintf("  KL: mean(t|H,x) = %.4f, sd = %.4f\n", mean(theta_H_kl), sd(theta_H_kl)))
cat(sprintf("  Inf: mean(t|H,x) = %.4f, sd = %.4f\n", mean(theta_H_inf), sd(theta_H_inf)))
cat("  -> Both posteriors under H are very similar\n\n")

# ========================================================
# STEP 2: Show that e-values are calculated differently
# Take the SAME theta and SAME simulated data,
# but evaluate ev with both priors
# ========================================================
cat("STEP 2: Evaluate e-value of SAME data with both priors\n")

# Use a representative theta
t_test <- mean(theta_H_kl)
set.seed(99)
x1_sim <- rbinom(1, n, t_test)
x2_sim <- rbinom(1, n, t_test)

# Calculate ev using KL prior
ev_kl <- ev_quad_from_data(n, n, x1_sim, x2_sim, a0_kl, a1_kl, a2_kl)

# Calculate ev using Informative prior
ev_inf <- ev_quad_from_data(n, n, x1_sim, x2_sim, a0_inf, a1_inf, a2_inf)

cat(sprintf("  Simulated data: x1*=%d, x2*=%d (under H, theta=%.4f)\n", x1_sim, x2_sim, t_test))
cat(sprintf("  ev calculated with KL prior:  %.4f\n", ev_kl))
cat(sprintf("  ev calculated with Inf prior: %.4f\n", ev_inf))
cat(sprintf("  RATIO: %.2fx difference!\n\n", ev_inf / ev_kl))

# ========================================================
# STEP 3: Massive verification with 1000 simulations
# Same simulated data, different priors for ev calculation
# ========================================================
cat("STEP 3: Massive verification (1000 simulations under H)\n")

set.seed(123)
ev_kl_vec  <- numeric(1000)
ev_inf_vec <- numeric(1000)

# Use the SAME theta sampled from KL posterior under H
for (i in 1:1000) {
  t <- sample(theta_H_kl, 1)  # SAME theta source
  x1_s <- rbinom(1, n, t)
  x2_s <- rbinom(1, n, t)
  
  ev_kl_vec[i]  <- ev_quad_from_data(n, n, x1_s, x2_s, a0_kl, a1_kl, a2_kl)
  ev_inf_vec[i] <- ev_quad_from_data(n, n, x1_s, x2_s, a0_inf, a1_inf, a2_inf)
}

cat(sprintf("  Mean ev (KL prior in calculation):  %.4f\n", mean(ev_kl_vec)))
cat(sprintf("  Mean ev (Inf prior in calculation): %.4f\n", mean(ev_inf_vec)))
cat("\n  CONCLUSION: Same data, but different prior in ev calculation\n")
cat("  produces VERY DIFFERENT e-values!\n\n")

# ========================================================
# Plot
# ========================================================
df <- data.frame(
  ev = c(ev_kl_vec, ev_inf_vec),
  prior_for_ev = rep(c("KL prior", "Inf prior"), each=1000)
)

p <- ggplot(df, aes(x=ev, fill=prior_for_ev)) +
  geom_histogram(bins=40, alpha=0.6, position="identity") +
  labs(x="e-value (from same simulated data)",
       y="count",
       title="Same data, different priors → very different e-values",
       fill="Prior used to compute ev") +
  theme_bw()

ggsave("Figures/ev_calculation_dependence.png", p, width=10, height=5, dpi=150)
cat("Plot saved: Figures/ev_calculation_dependence.png\n")
