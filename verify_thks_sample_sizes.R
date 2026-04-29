# verify_thks_sample_sizes.R
# Verifica los tamanos de muestra usados en el analisis THKS
# y los compara con lo que se obtiene del dataset completo.
#
# El codigo original (original_poster_code/THKS_data.R) filtra por colegio
# de forma INCONSISTENTE:
#   - yy (CC + TV)     -> filter(school == "404")           -> n chico
#   - yn (CC, no TV)   -> filter(school == "408")           -> n chico
#   - ny (no CC, TV)   -> NO filtra colegio (usa datos1)    -> n grande (~416)
#   - nn (no CC, no TV)-> filter(school == "409")           -> n chico
#
# Verificamos: que da si usamos TODOS los colegios para los 4 grupos.

options(timeout = 60)

url <- "https://content.sph.harvard.edu/fitzmaur/ala2e/tvsfp-data.txt"
local_path <- "tvsfp_raw.txt"
if (!file.exists(local_path)) {
  cat("Descargando dataset desde Harvard SPH...\n")
  download.file(url, local_path, quiet = TRUE)
}

# Estructura (segun documentacion ALA):
#   1: school id
#   2: student id
#   3: cc        (curriculum / school-based: 1 = yes, 0 = no)
#   4: tv        (TV-based: 1 = yes, 0 = no)
#   5: prethks   (pretest THKS, 0-7)
#   6: postthks  (posttest THKS, 0-7)
d <- read.table(local_path,
                col.names = c("school", "id", "cc", "tv", "prethks", "postthks"))
cat("Filas:", nrow(d), " | Colegios unicos:", length(unique(d$school)), "\n\n")

# Binarizamos: high THKS si >= 3 (mismo criterio del codigo original)
d$pre_bin  <- as.integer(d$prethks  >= 3)
d$post_bin <- as.integer(d$postthks >= 3)

# Etiqueta de grupo
d$group <- with(d, ifelse(cc == 1 & tv == 1, "yy",
                   ifelse(cc == 1 & tv == 0, "yn",
                   ifelse(cc == 0 & tv == 1, "ny", "nn"))))

# ---------------------------------------------------------------------------
# (1) Tamanos por colegio dentro de cada grupo (para entender el filtrado)
# ---------------------------------------------------------------------------
cat("=== Distribucion de colegios por grupo (n estudiantes por colegio) ===\n")
tab_school_group <- table(d$school, d$group)
print(tab_school_group)
cat("\n")

# ---------------------------------------------------------------------------
# (2) Replicar el filtrado del codigo original
# ---------------------------------------------------------------------------
cat("=== Filtrado ORIGINAL (un colegio por grupo, salvo ny) ===\n")
orig <- list(
  yy = subset(d, school == 404 & group == "yy"),
  yn = subset(d, school == 408 & group == "yn"),
  ny = subset(d, group == "ny"),                 # <- TODOS los colegios
  nn = subset(d, school == 409 & group == "nn")
)
res_orig <- do.call(rbind, lapply(names(orig), function(g) {
  s <- orig[[g]]
  data.frame(group = g, fuente = "original",
             n = nrow(s),
             x1_pre  = sum(s$pre_bin),
             x2_post = sum(s$post_bin),
             p1_pre  = round(mean(s$pre_bin), 3),
             p2_post = round(mean(s$post_bin), 3))
}))
print(res_orig); cat("\n")

# ---------------------------------------------------------------------------
# (3) Filtrado CONSISTENTE: todos los colegios para los 4 grupos
# ---------------------------------------------------------------------------
cat("=== Filtrado CONSISTENTE (todos los colegios, por grupo) ===\n")
res_full <- do.call(rbind, lapply(c("yy","yn","ny","nn"), function(g) {
  s <- subset(d, group == g)
  data.frame(group = g, fuente = "todos colegios",
             n = nrow(s),
             x1_pre  = sum(s$pre_bin),
             x2_post = sum(s$post_bin),
             p1_pre  = round(mean(s$pre_bin), 3),
             p2_post = round(mean(s$post_bin), 3))
}))
print(res_full); cat("\n")

# ---------------------------------------------------------------------------
# (4) Comparacion lado a lado
# ---------------------------------------------------------------------------
cat("=== Comparacion (n original vs n consistente) ===\n")
cmp <- merge(
  res_orig[, c("group", "n", "x1_pre", "x2_post")],
  res_full[, c("group", "n", "x1_pre", "x2_post")],
  by = "group", suffixes = c("_orig", "_full")
)
cmp$ratio_n <- round(cmp$n_full / cmp$n_orig, 1)
print(cmp); cat("\n")

# ---------------------------------------------------------------------------
# (5) Notas finales
# ---------------------------------------------------------------------------
cat("=== Notas ===\n")
cat("- Colegios con cada combinacion CC x TV:\n")
xt <- with(d, table(school, cc, tv))
present <- apply(xt > 0, c(2,3), function(x) paste(names(x)[x], collapse=","))
cat("  CC=1 TV=1 (yy):", present["1","1"], "\n")
cat("  CC=1 TV=0 (yn):", present["1","0"], "\n")
cat("  CC=0 TV=1 (ny):", present["0","1"], "\n")
cat("  CC=0 TV=0 (nn):", present["0","0"], "\n")

saveRDS(list(per_school = tab_school_group,
             original = res_orig,
             full     = res_full,
             comparison = cmp),
        "output/thks_sample_size_check.rds")
cat("\n--> Guardado en output/thks_sample_size_check.rds\n")
