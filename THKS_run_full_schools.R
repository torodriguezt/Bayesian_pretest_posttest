# THKS_run_full_schools.R
# -----------------------------------------------------------------------------
# Re-corre el analisis FBST sobre los datos THKS usando *todos* los colegios
# por grupo (CC x TV), para los 3 priors KL-fitted.
#
# Diferencia con THKS_run.R: el original filtra por un colegio por grupo
# (404, 408, 409) excepto para "no CC - TV" que ya usa todos. Aqui usamos
# TODOS los colegios para los 4 grupos -> tamaños comparables (~380-421).
#
# Salida: output/thks_full_schools_results.rds + tabla LaTeX 3-prior.
# -----------------------------------------------------------------------------

library(Rcpp)
library(dplyr)

sourceCpp("BivBetaBinom.cpp")
source("priors_config.R")

dir.create("output", showWarnings = FALSE)

# ---- Cargar datos crudos (Harvard SPH) ----
local_path <- "tvsfp_raw.txt"
if (!file.exists(local_path)) {
  options(timeout = 60)
  download.file("https://content.sph.harvard.edu/fitzmaur/ala2e/tvsfp-data.txt",
                local_path, quiet = TRUE)
}
d <- read.table(local_path,
                col.names = c("school", "id", "cc", "tv", "prethks", "postthks"))
d$pre_bin  <- as.integer(d$prethks  >= 3)
d$post_bin <- as.integer(d$postthks >= 3)
d$group <- with(d, ifelse(cc == 1 & tv == 1, "yy",
                   ifelse(cc == 1 & tv == 0, "yn",
                   ifelse(cc == 0 & tv == 1, "ny", "nn"))))

# ---- Resumir por grupo (todos los colegios) ----
group_data <- function(g) {
  s <- subset(d, group == g)
  list(n = nrow(s), x1 = sum(s$pre_bin), x2 = sum(s$post_bin),
       p1 = mean(s$pre_bin), p2 = mean(s$post_bin))
}
groups <- list(
  yy = list(d = group_data("yy"), lab = "CC + TV"),
  yn = list(d = group_data("yn"), lab = "CC, no TV"),
  ny = list(d = group_data("ny"), lab = "no CC, TV"),
  nn = list(d = group_data("nn"), lab = "no CC, no TV")
)

cat("\n=== Tamaños (todos los colegios, criterio consistente) ===\n")
size_tab <- do.call(rbind, lapply(groups, function(g) {
  data.frame(group = g$lab, n = g$d$n, x1_pre = g$d$x1, x2_post = g$d$x2,
             p1 = round(g$d$p1, 3), p2 = round(g$d$p2, 3))
}))
print(size_tab, row.names = FALSE); cat("\n")

# ---- Funcion principal: un grupo, un prior ----
run_group_prior <- function(g_entry, prior_alpha, M = 2000, seed = 7) {
  d_ <- g_entry$d; lab <- g_entry$lab
  n <- d_$n; x1 <- d_$x1; x2 <- d_$x2
  a0 <- prior_alpha["a0"]; a1 <- prior_alpha["a1"]; a2 <- prior_alpha["a2"]
  set.seed(seed)

  # P(H|x) = P(theta1 <= theta2 | x) por cuadratura sobre la grilla
  consts <- bb_constants(n, n, x1, x2, a0, a1, a2)
  ngrid <- 200
  xs <- seq(0.001, 0.999, length.out = ngrid)
  ys <- xs
  z  <- densBB_grid(xs, ys, consts)
  mask <- outer(seq_len(ngrid), seq_len(ngrid), function(i, j) xs[i] <= ys[j])
  ph_x <- sum(z[mask]) / sum(z)

  ev_obs <- ev_quad_from_data(n, n, x1, x2, a0, a1, a2)
  ev_H   <- simulate_evs_H_post(n, n, x1, x2, a0, a1, a2, M, 401)
  ev_A   <- simulate_evs_A_post(n, n, x1, x2, a0, a1, a2, M, 401)
  opt    <- find_kstar(ev_H, ev_A, 1, 1)

  data.frame(group = lab, n = n, x1 = x1, x2 = x2,
             ph_x = round(ph_x, 4),
             ev_obs = round(ev_obs, 4),
             k_star = round(opt$k_star, 4),
             alpha  = round(opt$alpha, 4),
             beta   = round(opt$beta, 4),
             decision = ifelse(ev_obs <= opt$k_star,
                               "Reject H", "Do not reject"))
}

# ---- Correr los 3 priors ----
priors <- list(
  "Non-informative"      = prior_NI,
  "Informative (mu=0.5)" = prior_INF,
  "Conflict (mu=0.1)"    = prior_CONF
)
M <- 2000

results <- list()
for (pname in names(priors)) {
  cat("--- Prior:", pname, "---\n")
  res_p <- do.call(rbind, lapply(groups, function(g) {
    cat(" ", g$lab, "\n")
    run_group_prior(g, priors[[pname]], M = M)
  }))
  res_p$prior <- pname
  results[[pname]] <- res_p
}
results <- do.call(rbind, results)
results$prior <- factor(results$prior, levels = names(priors))
rownames(results) <- NULL

cat("\n=== RESULTADOS (todos los colegios, 3 priors) ===\n")
print(results, row.names = FALSE)

saveRDS(list(results = results, sizes = size_tab),
        "output/thks_full_schools_results.rds")

# ---- LaTeX: una tabla por prior ----
build_table <- function(sub, label, caption) {
  hd <- c(
    "\\begin{table*}[!h]", "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lccccc@{}}",
    "\\toprule",
    "& \\multicolumn{1}{c}{\\textbf{Non-precise}} & \\multicolumn{3}{c}{\\textbf{Precise}} & \\\\",
    "\\cmidrule(lr){2-2}\\cmidrule(lr){3-5}",
    "\\textbf{Study conditions} & $P(\\textbf{H}|\\bX)$ & $ev(\\textbf{H};\\bX)$ & $k^*$ & Decision & $n_1=n_2$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    hd <- c(hd, sprintf(
      "\\textbf{%s} & %.4f & %.4f & %.4f & %s & %d \\\\",
      r$group, r$ph_x, r$ev_obs, r$k_star, r$decision, r$n))
  }
  c(hd, "\\bottomrule", "\\end{tabular*}", "\\end{table*}")
}

tex2 <- build_table(
  results[results$prior == "Non-informative", ],
  label   = "tab2",
  caption = sprintf(paste0(
    "TVSFP Project (all schools per group): Bayesian hypothesis testing with ",
    "non-informative prior (KL-optimal: $\\alpha_0=%.3f,\\,\\alpha_1=%.3f,\\,",
    "\\alpha_2=%.3f$), posterior-based formulation, $a=b=1$, $M=%d$."),
    prior_NI["a0"], prior_NI["a1"], prior_NI["a2"], M))

tex3 <- build_table(
  results[results$prior == "Informative (mu=0.5)", ],
  label   = "tab3",
  caption = sprintf(paste0(
    "TVSFP Project (all schools per group): Bayesian hypothesis testing with ",
    "informative prior ($\\alpha_0=%.2f,\\,\\alpha_1=%.2f,\\,\\alpha_2=%.2f$, ",
    "KL fit to Beta($N\\mu, N(1-\\mu)$)$^2$ with $\\mu=0.5,\\,N=50$), $M=%d$."),
    prior_INF["a0"], prior_INF["a1"], prior_INF["a2"], M))

tex4 <- build_table(
  results[results$prior == "Conflict (mu=0.1)", ],
  label   = "tab4",
  caption = sprintf(paste0(
    "TVSFP Project (all schools per group): Bayesian hypothesis testing with ",
    "informative prior in conflict with data ",
    "($\\alpha_0=%.2f,\\,\\alpha_1=%.2f,\\,\\alpha_2=%.2f$, ",
    "KL fit to Beta($N\\mu, N(1-\\mu)$)$^2$ with $\\mu=0.1,\\,N=50$, ",
    "$E[\\theta]\\approx 0.1$), $M=%d$."),
    prior_CONF["a0"], prior_CONF["a1"], prior_CONF["a2"], M))

writeLines(c(tex2, "", tex3, "", tex4), "output/thks_full_schools_tables.tex")
cat("\n--> output/thks_full_schools_tables.tex (3 tablas, todos los colegios)\n")
cat("--> output/thks_full_schools_results.rds\n")
