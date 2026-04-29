# fit_priors_kl.R
# -----------------------------------------------------------------------------
# Encuentra dos prioris dentro de la familia beta bivariada de Olkin-Liu
# (parametrizada por alpha0, alpha1, alpha2) por minimizacion de divergencia
# Kullback-Leibler:
#
#   1) Priori NO INFORMATIVA  = argmin_alpha  D_KL( U || p_alpha )
#      donde U es la uniforme bivariada en [0,1]^2 (referencia "ignorante").
#      Esto reproduce y formaliza el procedimiento del articulo (sec. 2.4).
#
#   2) Priori INFORMATIVA     = argmin_alpha  D_KL( p_alpha || pi_target )
#      donde pi_target es una distribucion concentrada que codifica la creencia
#      a priori (por defecto: producto de Beta(N*mu_j, N*(1-mu_j)) con N grande).
#
# Salida: output/kl_priors.rds  con ambas configuraciones (alphas, KL, ESS).
# -----------------------------------------------------------------------------

library(Rcpp)
sourceCpp("BivBetaBinom.cpp")  # provee sample_prior(n, a0, a1, a2)

# ============================================================================
# 1. Densidad de Olkin-Liu (forma cerrada)
# ============================================================================
# Derivacion: V0~Gamma(a0), V1~Gamma(a1), V2~Gamma(a2) indep.,
# Y_j = V_j / (V_j + V0). Cambio de variable + integracion en V0 da:
#
#   p(y1,y2) = Gamma(a)/[Gamma(a0)Gamma(a1)Gamma(a2)] *
#              y1^(a1-1) y2^(a2-1) (1-y1)^(a0+a2-1) (1-y2)^(a0+a1-1)
#              (1 - y1*y2)^(-a),       a = a0 + a1 + a2.

log_p_olkin_liu <- function(y1, y2, a0, a1, a2) {
  a <- a0 + a1 + a2
  lgamma(a) - lgamma(a0) - lgamma(a1) - lgamma(a2) +
    (a1 - 1) * log(y1) + (a2 - 1) * log(y2) +
    (a0 + a2 - 1) * log1p(-y1) + (a0 + a1 - 1) * log1p(-y2) -
    a * log1p(-y1 * y2)
}

# ============================================================================
# 2. Estimadores Monte Carlo de KL
# ============================================================================

# D_KL(U || p_alpha) = -E_U[log p_alpha]   (U uniforme => log U = 0)
# Aproximacion: -mean(log p_alpha(theta_j)),  theta_j ~ U[0,1]^2.
# Es la "M-projection" de la uniforme sobre la familia Olkin-Liu:
# encuentra el alpha que mejor "explica" muestras uniformes.
kl_U_to_p <- function(par, theta_unif) {
  a <- exp(par)  # exp() para forzar a > 0
  -mean(log_p_olkin_liu(theta_unif[, 1], theta_unif[, 2], a[1], a[2], a[3]))
}

# D_KL(p_alpha || pi_target) = E_{p_alpha}[log p_alpha - log pi_target]
# Muestras de p_alpha por construccion gamma (sample_prior del C++).
kl_p_to_target <- function(par, log_target_fun, M, eps = 1e-10) {
  a <- exp(par)
  th <- sample_prior(M, a[1], a[2], a[3])
  y1 <- pmin(pmax(th[, 1], eps), 1 - eps)
  y2 <- pmin(pmax(th[, 2], eps), 1 - eps)
  mean(log_p_olkin_liu(y1, y2, a[1], a[2], a[3]) - log_target_fun(y1, y2))
}

# ============================================================================
# 3. Priori NO INFORMATIVA: minimizar D_KL(U || p_alpha)
# ============================================================================

fit_noninformative_kl <- function(M = 100000, seed = 42, init = c(1, 1, 1)) {
  set.seed(seed)
  theta_unif <- cbind(runif(M), runif(M))
  res <- optim(log(init), kl_U_to_p, theta_unif = theta_unif,
               method = "Nelder-Mead",
               control = list(reltol = 1e-10, maxit = 5000))
  alpha_hat <- exp(res$par)
  list(alpha = alpha_hat,
       kl    = res$value,
       N0    = sum(alpha_hat),
       conv  = res$convergence,
       method = "argmin D_KL(U || p_alpha), MC con M = ",
       M = M)
}

# ============================================================================
# 4. Priori INFORMATIVA: minimizar D_KL(p_alpha || pi_target)
# ============================================================================
# pi_target = Beta(N*mu_1, N*(1-mu_1)) x Beta(N*mu_2, N*(1-mu_2))  (independientes)
# - mu controla DONDE se concentra la creencia previa (modo/media).
# - N controla CUAN concentrada es (ESS de la referencia).

make_log_target_indep_beta <- function(mu = c(0.5, 0.5), N = 50) {
  a_t <- N * mu
  b_t <- N * (1 - mu)
  function(y1, y2) {
    dbeta(y1, a_t[1], b_t[1], log = TRUE) +
      dbeta(y2, a_t[2], b_t[2], log = TRUE)
  }
}

fit_informative_kl <- function(mu = c(0.5, 0.5), N_target = 50,
                               M = 5000, n_replicates = 10,
                               seed = 42, init = c(10, 10, 10)) {
  log_target <- make_log_target_indep_beta(mu, N_target)
  # Promedio sobre replicas para suavizar ruido MC en el objetivo
  obj <- function(par) {
    set.seed(seed)
    mean(replicate(n_replicates, kl_p_to_target(par, log_target, M)))
  }
  res <- optim(log(init), obj, method = "Nelder-Mead",
               control = list(reltol = 1e-6, maxit = 3000))
  alpha_hat <- exp(res$par)
  list(alpha   = alpha_hat,
       kl      = res$value,
       N0      = sum(alpha_hat),
       conv    = res$convergence,
       target  = list(mu = mu, N = N_target,
                      family = "indep Beta(N*mu, N*(1-mu))"),
       M = M, n_replicates = n_replicates)
}

# ============================================================================
# 5. Ejecucion
# ============================================================================

cat("========================================================\n")
cat("  Ajuste de prioris (Olkin-Liu) por divergencia KL\n")
cat("========================================================\n\n")

dir.create("output", showWarnings = FALSE)

# ---- No informativa ----
cat("[1] Priori NO INFORMATIVA  -- min D_KL(U || p_alpha)\n")
fit_NI <- fit_noninformative_kl(M = 100000, seed = 42)
cat(sprintf("    alpha = (%.6f, %.6f, %.6f)\n",
            fit_NI$alpha[1], fit_NI$alpha[2], fit_NI$alpha[3]))
cat(sprintf("    ESS N0 = %.3f   |   KL minima = %.6f\n",
            fit_NI$N0, fit_NI$kl))
cat(sprintf("    Hardcode actual:  (0.8373879, 0.8410984, 0.8053298)\n\n"))

# ---- Informativa (simetrica, mu = 0.5) ----
cat("[2a] Priori INFORMATIVA simetrica  -- min D_KL(p_alpha || Beta x Beta)\n")
cat("     Target: mu = (0.5, 0.5), N = 50\n")
fit_INF_sim <- fit_informative_kl(mu = c(0.5, 0.5), N_target = 50,
                                  M = 5000, n_replicates = 10)
cat(sprintf("     alpha = (%.4f, %.4f, %.4f)\n",
            fit_INF_sim$alpha[1], fit_INF_sim$alpha[2], fit_INF_sim$alpha[3]))
cat(sprintf("     ESS N0 = %.3f   |   KL = %.6f\n\n",
            fit_INF_sim$N0, fit_INF_sim$kl))

# ---- Informativa con conflicto (mu = 0.1, mimica Tab4 del paper) ----
cat("[2b] Priori INFORMATIVA con conflicto  -- target mu = (0.1, 0.1), N = 50\n")
fit_INF_conf <- fit_informative_kl(mu = c(0.1, 0.1), N_target = 50,
                                   M = 5000, n_replicates = 10,
                                   init = c(40, 5, 5))
cat(sprintf("     alpha = (%.4f, %.4f, %.4f)\n",
            fit_INF_conf$alpha[1], fit_INF_conf$alpha[2], fit_INF_conf$alpha[3]))
cat(sprintf("     ESS N0 = %.3f   |   KL = %.6f\n\n",
            fit_INF_conf$N0, fit_INF_conf$kl))

# ---- Tabla resumen ----
resumen <- data.frame(
  prior = c("No informativa (KL a U)",
            "Informativa simetrica (mu=0.5, N=50)",
            "Informativa conflicto (mu=0.1, N=50)"),
  alpha0 = c(fit_NI$alpha[1], fit_INF_sim$alpha[1], fit_INF_conf$alpha[1]),
  alpha1 = c(fit_NI$alpha[2], fit_INF_sim$alpha[2], fit_INF_conf$alpha[2]),
  alpha2 = c(fit_NI$alpha[3], fit_INF_sim$alpha[3], fit_INF_conf$alpha[3]),
  ESS_N0 = c(fit_NI$N0, fit_INF_sim$N0, fit_INF_conf$N0),
  KL_min = c(fit_NI$kl, fit_INF_sim$kl, fit_INF_conf$kl)
)
cat("Resumen:\n"); print(resumen, row.names = FALSE)

saveRDS(list(noninformative      = fit_NI,
             informative_sym     = fit_INF_sim,
             informative_conflict = fit_INF_conf,
             resumen             = resumen),
        "output/kl_priors.rds")
cat("\n--> Guardado en output/kl_priors.rds\n")
