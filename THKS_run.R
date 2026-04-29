# THKS_run.R — Pipeline acelerado con C++ (Rcpp).
# Usa BivBetaBinom.cpp para todas las evaluaciones de la posterior.

library(ALA)
library(dplyr)
library(Rcpp)
library(GA)
library(pracma)

# 1) Compila el módulo C++ (una sola vez por sesión)
sourceCpp("BivBetaBinom.cpp")

# 2) Datos --------------------------------------------------------------
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
  ny = prep_global(datos1,         "no",  "yes"),  # original usa datos1 directo
  nn = prep_school(datos1, "409", "no",  "no")
)

# 3) Hiperparámetros ----------------------------------------------------
source("priors_config.R")
a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]

# 4) Pipeline para un grupo (sin MCMC) ---------------------------------
analyze_group <- function(g, label, do_plot = TRUE) {
  n1 <- g$n[1]; n2 <- g$n[2]
  x1 <- g$X[1]; x2 <- g$X[2]
  a  <- a0 + a1 + a2
  stopifnot((n1 + a) + (n2 + a) > a + (x1 + a1) + (x2 + a2))

  consts <- bb_constants(n1, n2, x1, x2, a0, a1, a2)

  # 4a) Superficie a posteriori
  if (do_plot) {
    xs <- seq(0.01, 0.8, 0.01); ys <- seq(0.01, 0.8, 0.01)
    z  <- densBB_grid(xs, ys, consts)
    persp(xs, ys, z, theta = -30, phi = 25, shade = 0.75,
          col = "gold", expand = 0.5, r = 2, ltheta = 25,
          ticktype = "detailed",
          xlab = "theta_1", ylab = "theta_2", zlab = "",
          main = label)
  }

  # 4b) Moda vía GA (la fitness es una llamada C++ pura)
  GA_mod <- ga(type = "real-valued",
               fitness = function(v) densBB_cpp(v[1], v[2], consts),
               lower = c(1e-4, 1e-4), upper = c(0.99, 0.99),
               popSize = 50, maxiter = 1000, run = 100, monitor = FALSE)
  mode_th <- as.numeric(GA_mod@solution[1, ])

  # 4c) Medias a posteriori vía integral2 (acepta matrices)
  esp1m <- function(t1, t2) {
    matrix(as.vector(t1) * densBB_vec(as.vector(t1), as.vector(t2), consts),
           nrow = nrow(t1))
  }
  esp2m <- function(t1, t2) {
    matrix(as.vector(t2) * densBB_vec(as.vector(t1), as.vector(t2), consts),
           nrow = nrow(t1))
  }
  E_th1 <- integral2(esp1m, 0, 1, 0, 1)$Q
  E_th2 <- integral2(esp2m, 0, 1, 0, 1)$Q

  # 4d) Supremo bajo H (theta1 = theta2)
  GA_H <- ga(type = "real-valued",
             fitness = function(v) densBB_H_cpp(v[1], consts),
             lower = 1e-4, upper = 0.99,
             popSize = 50, maxiter = 1000, run = 100, monitor = FALSE)
  th_sup <- as.numeric(GA_H@solution[1, 1])
  sup_H  <- densBB_H_cpp(th_sup, consts)

  list(consts = consts,
       mode   = mode_th,
       mean   = c(theta1 = E_th1, theta2 = E_th2),
       sup_H  = sup_H,
       th_sup = th_sup)
}

# 5) MCMC + e-valor FBST para un grupo ---------------------------------
mcmc_and_ev <- function(g, res, iter = 10000, cores = 4) {
  library(rstan)
  stan_data <- list(P = length(g$X), X = g$X, n = as.integer(g$n),
                    alpha1 = a0, alpha2 = a1, alpha3 = a2)
  fit    <- stan(file = "BBpost3.stan", data = stan_data,
                 iter = iter, cores = cores)
  params <- rstan::extract(fit)
  th1    <- params$Theta[, 1]
  th2    <- params$Theta[, 2]

  prob_gt <- mean(th1 >  th2)
  prob_le <- mean(th1 <= th2)
  ev10    <- ev_FBST(th1, th2, res$sup_H, res$consts)

  list(fit = fit, prob_gt = prob_gt, prob_le = prob_le, ev = ev10)
}

# 6) Ejecución sobre los 4 grupos --------------------------------------
labels <- c(yy = "school + tv", yn = "school - no tv",
            ny = "no school - tv", nn = "no school - no tv")

resultados <- list()
for (k in names(grupos)) {
  cat("\n===", labels[[k]], "===\n")
  resultados[[k]] <- analyze_group(grupos[[k]], labels[[k]])
  cat("Moda:  theta1 =", resultados[[k]]$mode[1],
      "theta2 =", resultados[[k]]$mode[2], "\n")
  cat("Media: theta1 =", resultados[[k]]$mean[1],
      "theta2 =", resultados[[k]]$mean[2], "\n")
  cat("Sup H:", resultados[[k]]$sup_H,
      "en theta* =", resultados[[k]]$th_sup, "\n")
}

# 7) MCMC + FBST (descomentar para ejecutar; es lo más caro) -----------
mc_yy <- mcmc_and_ev(grupos$yy, resultados$yy)
cat("P(theta1>theta2) =", mc_yy$prob_gt, "\n")
cat("FBST e-value     =", mc_yy$ev,      "\n")

# 8) Benchmark rápido (opcional) ---------------------------------------
library(microbenchmark)
g <- grupos$yy; n1<-g$n[1]; n2<-g$n[2]; x1<-g$X[1]; x2<-g$X[2]
a <- a0+a1+a2; u <- c(a, x1+a1, x2+a2); l <- c(n1+a, n2+a)
consts <- bb_constants(n1, n2, x1, x2, a0, a1, a2)
library(hypergeo)
densBB_R <- function(t1, t2) {
exp(lgamma(n1+a)+lgamma(n2+a)-lgamma(x1+a1)-lgamma(n1-x1+a-a1)
-lgamma(x2+a2)-lgamma(n2-x2+a-a2)
-log(genhypergeo(U=u, L=l, check_mod=TRUE, z=1))
+(a1+x1-1)*log(t1)+(a2+a0+(n1-x1)-1)*log(1-t1)
+(a2+x2-1)*log(t2)+(a1+a0+(n2-x2)-1)*log(1-t2)
-a*log(1-t1*t2))
}
microbenchmark(R   = densBB_R(0.3, 0.4),
                cpp = densBB_cpp(0.3, 0.4, consts), times = 200)
