# Deep debug: What's happening with simulated e-values?

library(Rcpp)
library(ggplot2)

setwd("c:/Users/Tomas/Bayesian_pretest_posttest")
sourceCpp("BivBetaBinom.cpp")

a0_kl <- 0.8373879; a1_kl <- 0.8410984; a2_kl <- 0.8053298
a0_inf <- 10; a1_inf <- 10; a2_inf <- 10

# Reproduce row 25-26: n=150, d=0.2, theta1=0.5, theta2=0.7
n <- 150
theta1 <- 0.5
theta2 <- 0.52  # small effect

set.seed(42)
x1 <- rbinom(1, n, theta1)
x2 <- rbinom(1, n, theta2)

cat(sprintf("SCENARIO: n=%d, theta1=%.2f, theta2=%.2f\n", n, theta1, theta2))
cat(sprintf("DATA: x1=%d, x2=%d\n\n", x1, x2))

# ===== KL-OPTIMAL =====
cat("="*60 %+% "\n")
cat("KL-OPTIMAL PRIOR:\n")
cat("="*60 %+% "\n")

ev_H_kl <- simulate_evs_H_post(n, n, x1, x2, a0_kl, a1_kl, a2_kl, M=1000)
ev_A_kl <- simulate_evs_A_post(n, n, x1, x2, a0_kl, a1_kl, a2_kl, M=1000)

cat(sprintf("ev_H: mean=%.4f, sd=%.4f, min=%.6f, max=%.6f\n", 
            mean(ev_H_kl), sd(ev_H_kl), min(ev_H_kl), max(ev_H_kl)))
cat(sprintf("ev_A: mean=%.4f, sd=%.4f, min=%.6f, max=%.6f\n", 
            mean(ev_A_kl), sd(ev_A_kl), min(ev_A_kl), max(ev_A_kl)))

opt_kl <- find_kstar(ev_H_kl, ev_A_kl, 1, 1)
cat(sprintf("\nk* = %.4f\n", opt_kl$k_star))
cat(sprintf("alpha(k*) = %.4f (proportion of ev_H <= k*)\n", opt_kl$alpha))
cat(sprintf("beta(k*) = %.4f (proportion of ev_A > k*)\n", opt_kl$beta))

# ===== INFORMATIVE =====
cat("\n" %+% "="*60 %+% "\n")
cat("INFORMATIVE PRIOR (alpha=10):\n")
cat("="*60 %+% "\n")

ev_H_inf <- simulate_evs_H_post(n, n, x1, x2, a0_inf, a1_inf, a2_inf, M=1000)
ev_A_inf <- simulate_evs_A_post(n, n, x1, x2, a0_inf, a1_inf, a2_inf, M=1000)

cat(sprintf("ev_H: mean=%.4f, sd=%.4f, min=%.6f, max=%.6f\n", 
            mean(ev_H_inf), sd(ev_H_inf), min(ev_H_inf), max(ev_H_inf)))
cat(sprintf("ev_A: mean=%.4f, sd=%.4f, min=%.6f, max=%.6f\n", 
            mean(ev_A_inf), sd(ev_A_inf), min(ev_A_inf), max(ev_A_inf)))

opt_inf <- find_kstar(ev_H_inf, ev_A_inf, 1, 1)
cat(sprintf("\nk* = %.4f\n", opt_inf$k_star))
cat(sprintf("alpha(k*) = %.4f (proportion of ev_H <= k*)\n", opt_inf$alpha))
cat(sprintf("beta(k*) = %.4f (proportion of ev_A > k*)\n", opt_inf$beta))

# ===== COMPARISON =====
cat("\n" %+% "="*60 %+% "\n")
cat("COMPARISON:\n")
cat("="*60 %+% "\n")
cat(sprintf("k* ratio (inf/kl): %.2f\n", opt_inf$k_star / opt_kl$k_star))
cat(sprintf("alpha ratio: %.2f\n", opt_inf$alpha / opt_kl$alpha))
cat(sprintf("beta ratio: %.2f\n", opt_inf$beta / opt_kl$beta))

cat("\nKEY QUESTION: Why is k* so much higher with informative prior?\n")
cat("Answer: Check the DISTRIBUTIONS of ev_H and ev_A\n\n")

# Plot the distributions
df <- rbind(
  data.frame(ev=ev_H_kl, type="ev_H", prior="KL-optimal"),
  data.frame(ev=ev_A_kl, type="ev_A", prior="KL-optimal"),
  data.frame(ev=ev_H_inf, type="ev_H", prior="Informative"),
  data.frame(ev=ev_A_inf, type="ev_A", prior="Informative")
)

p <- ggplot(df, aes(x=ev, fill=type)) +
  geom_histogram(bins=50, alpha=0.6, position="identity") +
  facet_wrap(~prior) +
  labs(title="Distribution of simulated e-values", x="e-value", y="count") +
  theme_bw()

ggsave("Figures/deep_debug_evdist.png", p, width=10, height=5, dpi=150)
cat("Plot saved: Figures/deep_debug_evdist.png\n")
