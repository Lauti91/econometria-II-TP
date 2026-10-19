# =============================================================================
# TP Econometria II - Candidata 5: Indice Li Keqiang (China)
# Construye en R, desde las fuentes originales, la base trimestral/mensual:
#   - PBI oficial trimestral de China (NBS, via chinadata.live)
#   - Generacion de electricidad mensual (OCDE via FRED  +  NBS via chinadata.live)
#   - Carga ferroviaria mensual (NBS, via chinadata.live)
# Salida: series/ (CSV tal cual se bajan) y, en output/,
#         base_candidata5_li_keqiang_china.xlsx (Notas, Trimestral, Mensual) + CSV
# Paquetes: solo R base + openxlsx (opcional: si no esta, escribe solo CSV)
# =============================================================================

# Ejecutar desde la raiz del repositorio: Rscript scripts/01_construir_base.R
dir_raw <- "series"
dir_out <- "output"
dir.create(dir_raw, showWarnings = FALSE)
dir.create(dir_out, showWarnings = FALSE)

# ---- 1. Descarga de las fuentes originales ---------------------------------
cd_api <- "https://chinadata.live/api/v2/data/"
fuentes <- c(
  gdp_nominal = paste0(cd_api, "china-gdp-quarterly?format=csv"),   # 1e8 yuan, trimestral
  gdp_indice  = paste0(cd_api, "china-gdp-index?format=csv"),       # ano previo = 100
  gdp_real    = paste0(cd_api, "china-gdp-real?format=csv"),        # real base 2011
  elec_cd     = paste0(cd_api, "china-electricity-generation?format=csv"),  # 1e8 kWh
  rail_cd     = paste0(cd_api, "china-rail-freight?format=csv"),            # 1e4 t
  elec_oecd   = "https://fred.stlouisfed.org/graph/fredgraph.csv?id=CHNPRENEL01MLM" # GWh
)
archivos <- file.path(dir_raw, paste0(names(fuentes), ".csv"))
names(archivos) <- names(fuentes)

for (k in names(fuentes)) {
  ok <- tryCatch({
    download.file(fuentes[[k]], archivos[[k]], mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) { message("FALLO descarga ", k, ": ", conditionMessage(e)); FALSE })
  if (!ok && !file.exists(archivos[[k]])) stop("No se pudo obtener ", k)
}

# ---- 2. Lectura y utilidades -----------------------------------------------
rd <- function(k) read.csv(archivos[[k]], stringsAsFactors = FALSE)

# "1992-Q1" -> clave numerica anio*4 + (q-1)
qkey <- function(s) as.integer(substr(s, 1, 4)) * 4L + as.integer(substr(s, 7, 7)) - 1L
qlab <- function(k) sprintf("%dQ%d", k %/% 4L, k %% 4L + 1L)
# "2007-02" -> clave numerica anio*12 + (m-1)
mkey <- function(s) as.integer(substr(s, 1, 4)) * 12L + as.integer(substr(s, 6, 7)) - 1L
mlab <- function(k) sprintf("%d-%02d", k %/% 12L, k %% 12L + 1L)

# ---- 3. Electricidad mensual -----------------------------------------------
# 3a. OCDE/FRED (GWh), 1999-01 .. 2023-11
o <- rd("elec_oecd")
oecd <- setNames(o$CHNPRENEL01MLM, mkey(substr(o$observation_date, 1, 7)))

# Correccion: desde 2012 en varios anios el "febrero" de la OCDE trae el
# acumulado ene-feb (= 2 x enero). Se reparte en partes iguales (Feb := Ene).
anios_corregidos <- c()
for (y in unique(as.integer(names(oecd)) %/% 12L)) {
  kj <- as.character(y * 12L); kf <- as.character(y * 12L + 1L)
  if (kj %in% names(oecd) && kf %in% names(oecd) &&
      abs(oecd[kf] / oecd[kj] - 2) < 0.005) {
    oecd[kf] <- oecd[kj]
    anios_corregidos <- c(anios_corregidos, y)
  }
}
message("Anios con febrero corregido: ", paste(anios_corregidos, collapse = ", "))

# 3b. chinadata.live (1e8 kWh -> x100 = GWh), se usa desde donde termina la OCDE
e <- rd("elec_cd")
elec_cd <- setNames(e$value * 100, mkey(e$date))
ext <- elec_cd[as.integer(names(elec_cd)) > max(as.integer(names(oecd)))]

# Control: en el solape mar-dic ambas fuentes deben coincidir (ratio ~ 1)
sol <- intersect(names(oecd), names(elec_cd))
sol <- sol[(as.integer(sol) %% 12L) >= 2L]
stopifnot(all(abs(oecd[sol] / elec_cd[sol] - 1) < 0.001))

elec <- c(oecd, ext)
elec_fuente <- c(rep("OCDE/FRED (Feb corregido)", length(oecd)),
                 rep("chinadata.live (NBS) x100", length(ext)))
names(elec_fuente) <- names(elec)

# ---- 4. Carga ferroviaria mensual (1e4 t) ----------------------------------
r <- rd("rail_cd")
rail <- setNames(r$value, mkey(r$date))

# ---- 5. Panel mensual ------------------------------------------------------
k_min <- mkey("1999-01")
k_max <- max(as.integer(names(elec)), as.integer(names(rail)))
kk <- k_min:k_max
M <- data.frame(
  mes = mlab(kk),
  electricidad_GWh = unname(elec[as.character(kk)]),
  electricidad_fuente = unname(elec_fuente[as.character(kk)]),
  carga_ferroviaria_10k_t = unname(rail[as.character(kk)]),
  stringsAsFactors = FALSE
)

# ---- 6. Panel trimestral ---------------------------------------------------
# Suma de 3 meses; si falta alguno queda NA (no se imputa).
sum_trim <- function(x, kk) {
  tq <- (kk %/% 12L) * 4L + ((kk %% 12L) %/% 3L)
  res <- tapply(x, tq, function(v) if (length(v) == 3L && !anyNA(v)) sum(v) else NA_real_)
  setNames(as.numeric(res), names(res))
}
eq <- sum_trim(M$electricidad_GWh, kk)
rq <- sum_trim(M$carga_ferroviaria_10k_t, kk)

g  <- rd("gdp_nominal"); gi <- rd("gdp_indice"); gr <- rd("gdp_real")
kq <- qkey("1992-Q1"):qkey("2026-Q2")
pick <- function(df, kq) unname(setNames(df$value, qkey(df$date))[as.character(kq)])
rez4 <- function(x) c(rep(NA_real_, 4), head(x, -4))

Q <- data.frame(trimestre = qlab(kq), stringsAsFactors = FALSE)
Q$pbi_nominal_1e8_yuan            <- pick(g, kq)
Q$pbi_crec_interanual_oficial_pct <- pick(gi, kq) - 100
Q$pbi_real_1e8_yuan_2011          <- pick(gr, kq)
Q$electricidad_GWh                <- unname(eq[as.character(kq)])
Q$carga_ferroviaria_10k_t         <- unname(rq[as.character(kq)])
Q$elec_crec_interanual_pct        <- 100 * (Q$electricidad_GWh / rez4(Q$electricidad_GWh) - 1)
Q$rail_crec_interanual_pct        <- 100 * (Q$carga_ferroviaria_10k_t / rez4(Q$carga_ferroviaria_10k_t) - 1)
Q$indice_2comp_pct                <- (Q$elec_crec_interanual_pct + Q$rail_crec_interanual_pct) / 2
Q$dln_elec_trim                   <- c(NA, 100 * diff(log(Q$electricidad_GWh)))
Q$dln_pbi_nom_trim                <- c(NA, 100 * diff(log(Q$pbi_nominal_1e8_yuan)))

# ---- 7. Controles y conteos ------------------------------------------------
cat("\n--- Conteos de observaciones ---\n")
for (v in c("pbi_nominal_1e8_yuan", "pbi_crec_interanual_oficial_pct",
            "electricidad_GWh", "elec_crec_interanual_pct",
            "carga_ferroviaria_10k_t", "rail_crec_interanual_pct", "indice_2comp_pct")) {
  ok <- which(!is.na(Q[[v]]))
  cat(sprintf("%-34s n=%3d  %s .. %s\n", v, length(ok), Q$trimestre[min(ok)], Q$trimestre[max(ok)]))
}
cc <- Q[Q$trimestre >= "2000Q1" & Q$trimestre <= "2023Q4", ]
cat("\nPar PBI x electricidad, 2000Q1-2023Q4 contiguo:",
    sum(complete.cases(cc[, c("pbi_crec_interanual_oficial_pct", "elec_crec_interanual_pct")])), "obs\n")
cr <- Q[Q$trimestre >= "2006Q1" & Q$trimestre <= "2023Q4", ]
cat("Par con carga ferroviaria / indice 2comp, 2006Q1-2023Q4 contiguo:",
    sum(complete.cases(cr[, c("pbi_crec_interanual_oficial_pct", "rail_crec_interanual_pct")])), "obs\n")
cat("Correlacion PBI x electricidad (interanual):",
    round(cor(Q$pbi_crec_interanual_oficial_pct, Q$elec_crec_interanual_pct, use = "complete.obs"), 3), "\n")

# ---- 8. Salida -------------------------------------------------------------
Notas <- data.frame(Nota = c(
  "Base candidata 5 (Indice Li Keqiang). Construida en R desde las fuentes originales.",
  "PBI nominal trimestral (1e8 yuan, flujo trimestral) 1992Q1-2026Q2 y tasa interanual oficial (indice ano previo=100, menos 100) 1993Q1-2026Q2: chinadata.live /api/v2/data/china-gdp-quarterly y china-gdp-index (espejo NBS).",
  "PBI real (base 2011) solo desde 2011Q1 (china-gdp-real).",
  paste0("Electricidad: OCDE/FRED CHNPRENEL01MLM (GWh) 1999-01..2023-11, empalmada con chinadata.live china-electricity-generation (x100 -> GWh) desde 2023-12. ",
         "Desde 2012 en varios anios el febrero OCDE trae el acumulado ene-feb (=2x enero); se corrige repartiendo en partes iguales (anios: ",
         paste(anios_corregidos, collapse = ", "), "). En chinadata.live faltan Ene y casi siempre Feb, por eso Q1 de 2024-2026 queda vacio."),
  "Carga ferroviaria (1e4 t): chinadata.live china-rail-freight, 2005-01..2026-07; faltan 2024-01/02, 2025-01/02 y 2026-01 (acumulados ene-feb del NBS), por lo que Q1 de 2024-2026 queda vacio.",
  "Prestamos bancarios (tercer componente del Li Keqiang): sin serie larga gratuita accesible. El indice es de 2 componentes (promedio simple de crecimiento interanual de electricidad y carga ferroviaria).",
  "Trimestral: suma de 3 meses; si falta alguno, NA (no se imputa)."
))

write.csv(Q, file.path(dir_out, "trimestral_li_keqiang.csv"), row.names = FALSE, na = "")
write.csv(M, file.path(dir_out, "mensual_li_keqiang.csv"), row.names = FALSE, na = "")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(list(Notas = Notas, Trimestral = Q, Mensual = M),
                       file.path(dir_out, "base_candidata5_li_keqiang_china.xlsx"))
  message("Listo: ", file.path(dir_out, "base_candidata5_li_keqiang_china.xlsx"))
} else message("openxlsx no instalado: se escribieron solo los CSV")

