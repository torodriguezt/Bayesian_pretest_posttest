# priors_config.R
# -----------------------------------------------------------------------------
# Fuente unica de verdad para los hiperparametros de la beta bivariada
# Olkin-Liu, parametrizada por (alpha0, alpha1, alpha2). Todos los scripts
# del proyecto deben hacer:   source("priors_config.R")
#
# Los tres conjuntos fueron obtenidos por minimizacion de divergencia
# Kullback-Leibler en fit_priors_kl.R (semilla = 42, M = 1e5 / 5e3):
#
#   prior_NI   : argmin D_KL( U[0,1]^2 || p_alpha )            (no informativa)
#   prior_INF  : argmin D_KL( p_alpha || Beta(N*0.5, N*0.5)^2 ) con N=50
#   prior_CONF : argmin D_KL( p_alpha || Beta(N*0.1, N*0.9)^2 ) con N=50
# -----------------------------------------------------------------------------

# --- Priori NO INFORMATIVA (KL-optima a la uniforme bivariada) ---------------
prior_NI <- c(a0 = 0.760595, a1 = 0.762204, a2 = 0.758178)

# --- Priori INFORMATIVA simetrica (mu = 0.5, ESS-target = 50) ----------------
prior_INF <- c(a0 = 24.9938, a1 = 24.8026, a2 = 24.8843)

# --- Priori INFORMATIVA con CONFLICTO (mu = 0.1, ESS-target = 50) ------------
# E[theta_j] = a_j / (a_j + a_0) ~ 0.10, contradice datos THKS (theta ~ 0.3-0.6)
prior_CONF <- c(a0 = 45.8654, a1 = 5.1072, a2 = 5.0664)

# --- Tabla resumen (informativa para reportes) -------------------------------
priors_summary <- data.frame(
  prior  = c("No informativa (KL a U)",
             "Informativa simetrica (mu=0.5, N=50)",
             "Conflicto (mu=0.1, N=50)"),
  alpha0 = c(prior_NI["a0"],  prior_INF["a0"],  prior_CONF["a0"]),
  alpha1 = c(prior_NI["a1"],  prior_INF["a1"],  prior_CONF["a1"]),
  alpha2 = c(prior_NI["a2"],  prior_INF["a2"],  prior_CONF["a2"]),
  ESS_N0 = c(sum(prior_NI),   sum(prior_INF),   sum(prior_CONF)),
  row.names = NULL
)

# --- Compatibilidad hacia atras (alias usados en scripts viejos) -------------
# Los scripts originales esperan (a0, a1, a2) sueltos. Quien quiera el viejo
# nombre lo expone con: a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; ...
