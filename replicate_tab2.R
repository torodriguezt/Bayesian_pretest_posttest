# replicate_tab2.R
# Reproduce los k* del Cuadro 2 del artículo usando la formulación
# POSTERIOR-based (eqs 16-17 del artículo).
#
# Diferencia clave con build_kstar_table.R:
#   - prior-based  (10-11): θ ~ prior, da una tabla genérica k*(n1, n2)
#   - posterior-based (16-17): θ ~ posterior dado (x1_obs, x2_obs),
#     da un k* específico por experimento  ← lo que reporta tab2
#
# Si los k* aquí calculados se acercan a tab2, el motor (cuadratura + SIR)
# está validado y podemos usar build_kstar_table con confianza.

library(ALA)
library(dplyr)
library(Rcpp)

setwd("c:/Users/Tomas/BivBetaBinomial_Tomás/BivBetaBinomial_Tomás")
sourceCpp("BivBetaBinom.cpp")

source("priors_config.R")
a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]

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

# Mismas asignaciones que THKS_run.R (atención: yn=27, nn=80, ny=416 difieren
# de lo reportado en tab2 — ver nota al final).
grupos <- list(
  yy = prep_school(datos1, "404", "yes", "yes"),
  yn = prep_school(datos1, "408", "yes", "no"),
  ny = prep_global(datos1,         "no",  "yes"),
  nn = prep_school(datos1, "409", "no",  "no")
)

labels <- c(yy = "CC - TV", yn = "CC - No TV",
            ny = "No CC - TV", nn = "No CC - No TV")

# Valores reportados en tab2 del artículo
expected <- data.frame(
  group    = c("CC - TV", "CC - No TV", "No CC - TV", "No CC - No TV"),
  n_paper  = c(25,        38,           82,           33),
  ev_paper = c(0.0005,    0.0184,       1.0000,       0.9999),
  k_paper  = c(0.5455,    0.4424,       0.2432,       0.4605)
)

analyze_post <- function(g, M = 500, seed = 42, ngrid_quad = 401) {
  set.seed(seed)
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  ev_obs <- ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2, ngrid_quad)
  ess    <- sir_ess(n1, n2, x1, x2, a0, a1, a2, N_prop = max(50 * M, 10000))
  ev_H   <- simulate_evs_H_post(n1, n2, x1, x2, a0, a1, a2, M, ngrid_quad)
  ev_A   <- simulate_evs_A_post(n1, n2, x1, x2, a0, a1, a2, M, ngrid_quad)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)
  list(n1 = n1, n2 = n2, x1 = x1, x2 = x2,
       ev_obs = ev_obs, k_star = opt$k_star,
       alpha = opt$alpha, beta = opt$beta, ess = ess)
}

cat("\n=== Replicación de Tab. 2 (formulación posterior, M = 500) ===\n")
res <- lapply(names(grupos), function(k) {
  cat("Procesando", labels[[k]], "...\n")
  r <- analyze_post(grupos[[k]], M = 500, seed = 7)
  data.frame(group = labels[[k]],
             n1 = r$n1, n2 = r$n2, x1 = r$x1, x2 = r$x2,
             ev_obs = r$ev_obs, k_star = r$k_star,
             alpha = r$alpha, beta = r$beta, ess = round(r$ess, 0))
})
res <- do.call(rbind, res)
print(res, row.names = FALSE, digits = 4)

cat("\n=== Comparación contra valores del artículo ===\n")
cmp <- merge(res, expected, by = "group", sort = FALSE)
cmp$diff_ev <- cmp$ev_obs - cmp$ev_paper
cmp$diff_k  <- cmp$k_star - cmp$k_paper
print(cmp[, c("group", "n1", "n_paper", "ev_obs", "ev_paper", "diff_ev",
              "k_star", "k_paper", "diff_k")],
      row.names = FALSE, digits = 4)

cat("\n--- Notas ---\n")
cat("1) ESS de la SIR: si es < ~200, la posterior es muy concentrada para esa\n")
cat("   propuesta y conviene aumentar N_prop (parámetro de simulate_evs_A_post).\n")
cat("2) Sample-size mismatch: las asignaciones de escuela en este script vienen\n")
cat("   de THKS_run.R y no coinciden con tab2: yn=27 vs 38, nn=80 vs 33,\n")
cat("   ny=416 vs 82.  Sólo yy=25 coincide.  Si los k_star de yn/nn/ny no se\n")
cat("   acercan a tab2, hay que verificar qué subconjunto del THKS usa el\n")
cat("   artículo (¿otra escuela?, ¿agrupación distinta?, ¿filtro adicional?).\n")
cat("3) Si yy reproduce bien (ev_obs ≈ 0.0005, k* ≈ 0.5455), el motor está OK\n")
cat("   y se puede correr build_kstar_table.R con confianza.\n")

