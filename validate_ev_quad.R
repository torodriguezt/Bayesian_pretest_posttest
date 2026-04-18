# validate_ev_quad.R
# Compara el ev calculado por cuadratura 2D (ev_quad, sin Stan)
# contra el ev calculado por MCMC (ev_FBST, con Stan).
# Si coinciden hasta ~3-4 decimales, la cuadratura está validada
# y la podemos usar como motor para construir la tabla k*(n1, n2).

library(ALA)
library(dplyr)
library(Rcpp)
library(GA)
library(rstan)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")

# --- Datos (igual que THKS_run.R) ----------------------------------------
datos1 <- tvsfp

prep_school <- function(data, school_id, sb, tv) {
  d <- data %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school == school_id, school.based == sb, tv.based == tv)
  list(
    X = (d %>% group_by(stage) %>% summarise(Bin = sum(binTHKS)))$Bin,
    n = (d %>% group_by(stage) %>% summarise(n   = n()))$n
  )
}
prep_global <- function(data, sb, tv) {
  d <- data %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school.based == sb, tv.based == tv)
  list(
    X = (d %>% group_by(stage) %>% summarise(Bin = sum(binTHKS)))$Bin,
    n = (d %>% group_by(stage) %>% summarise(n   = n()))$n
  )
}

grupos <- list(
  yy = prep_school(datos1, "404", "yes", "yes"),
  yn = prep_school(datos1, "408", "yes", "no"),
  ny = prep_global(datos1,         "no",  "yes"),
  nn = prep_school(datos1, "409", "no",  "no")
)

alphas_opt <- c(0.8373879, 0.8410984, 0.8053298)
a0 <- alphas_opt[1]; a1 <- alphas_opt[2]; a2 <- alphas_opt[3]

# --- Compilar Stan una sola vez ------------------------------------------
stan_model_obj <- stan_model("BBpost3.stan")

# --- Validación grupo por grupo ------------------------------------------
compare_one <- function(g, label) {
  n1 <- g$n[1]; n2 <- g$n[2]
  x1 <- g$X[1]; x2 <- g$X[2]

  consts <- bb_constants(n1, n2, x1, x2, a0, a1, a2)
  sup_info <- find_sup_H(consts)
  sup_H <- sup_info$sup_H

  # Ruta A: cuadratura (sin Stan)
  t0 <- Sys.time()
  ev_q  <- ev_quad(consts, sup_H, ngrid = 401)
  t_q   <- as.numeric(Sys.time() - t0, units = "secs")

  # Ruta B: MCMC con Stan
  stan_data <- list(P = length(g$X), X = g$X, n = as.integer(g$n),
                    alpha1 = a0, alpha2 = a1, alpha3 = a2)
  t0 <- Sys.time()
  fit <- sampling(stan_model_obj, data = stan_data,
                  iter = 10000, chains = 4, cores = 4,
                  refresh = 0, show_messages = FALSE)
  params <- rstan::extract(fit)
  ev_m <- ev_FBST(params$Theta[, 1], params$Theta[, 2], sup_H, consts)
  t_m  <- as.numeric(Sys.time() - t0, units = "secs")

  data.frame(grupo = label,
             n1 = n1, n2 = n2, x1 = x1, x2 = x2,
             sup_H = sup_H,
             theta_argmax = sup_info$theta,
             ev_quad   = ev_q,
             ev_mcmc   = ev_m,
             diff_abs  = abs(ev_q - ev_m),
             tiempo_quad_s = round(t_q, 4),
             tiempo_mcmc_s = round(t_m, 1))
}

labels <- c(yy = "school + tv", yn = "school - no tv",
            ny = "no school - tv", nn = "no school - no tv")

resultados <- do.call(rbind, lapply(names(grupos), function(k) {
  cat("Procesando", labels[[k]], "...\n")
  compare_one(grupos[[k]], labels[[k]])
}))

cat("\n========== Comparación ev_quad vs ev_mcmc ==========\n")
print(resultados, row.names = FALSE, digits = 5)

cat("\nMáxima diferencia absoluta:", max(resultados$diff_abs), "\n")
cat("Speedup medio (mcmc / quad):",
    round(mean(resultados$tiempo_mcmc_s / resultados$tiempo_quad_s), 1), "x\n")
