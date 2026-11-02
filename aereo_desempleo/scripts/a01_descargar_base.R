# =============================================================================
# a01_descargar_base.R -- Par alternativo 4: trafico aereo de pasajeros y desempleo (EE.UU.)
# Descarga desde FRED (CSV directo, sin registro), guarda los crudos y arma la base mensual.
#   RPM        millas-pasajero pagas, vuelos regulares domesticos e internacionales de aerolineas de
#              EE.UU. (BTS), miles, mensual, SIN desestacionalizar
#   UNRATENSA  tasa de desempleo (BLS), %, mensual, SIN desestacionalizar
#   UNRATE     idem desestacionalizada (solo para comparar, regla 8 de la clase)
#   USREC      recesiones NBER, para sombrear los graficos
# Dato faltante: el BLS no publico la tasa de desempleo de octubre de 2025 (cierre del gobierno
# federal, oct-nov 2025). Se completa con el promedio de septiembre y noviembre y queda marcado;
# no afecta la muestra principal del analisis (2000-2019).
# Ejecutar desde aereo_desempleo/:  Rscript scripts/a01_descargar_base.R
# =============================================================================
ids <- c("RPM", "UNRATENSA", "UNRATE", "USREC")
fred <- function(id) {
  d <- read.csv(sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id), na.strings = ".")
  names(d) <- c("fecha", "valor"); d$fecha <- as.Date(d$fecha); d
}
crudos <- lapply(setNames(ids, ids), fred)
for (id in ids) write.csv(crudos[[id]], sprintf("datos/crudos/%s.csv", id), row.names = FALSE)
writeLines(sprintf("Descarga FRED: %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "datos/crudos/fecha_descarga.txt")

mensual <- function(id, nombre) { d <- crudos[[id]]; setNames(data.frame(format(d$fecha, "%Y-%m"), d$valor), c("mes", nombre)) }
B <- Reduce(function(a, b) merge(a, b, by = "mes", all.x = TRUE),
            list(mensual("RPM", "rpm_miles"), mensual("UNRATENSA", "desempleo_nsa"), mensual("UNRATE", "desempleo_sa"), mensual("USREC", "recesion")))
B <- B[!is.na(B$rpm_miles), ]                                # 2000-01 a la ultima publicacion del BTS

# Imputacion del faltante del BLS (y solo de ese)
falta <- which(is.na(B$desempleo_nsa))
stopifnot(identical(B$mes[falta], "2025-10"))
B$desempleo_imputado <- 0L; B$desempleo_imputado[falta] <- 1L
for (v in c("desempleo_nsa", "desempleo_sa")) B[[v]][falta] <- (B[[v]][falta - 1] + B[[v]][falta + 1]) / 2

# ---- Controles -------------------------------------------------------------------------------
meses <- format(seq(as.Date(paste0(B$mes[1], "-01")), as.Date(paste0(B$mes[nrow(B)], "-01")), by = "month"), "%Y-%m")
stopifnot(identical(B$mes, meses), !anyNA(B), all(B$rpm_miles > 0), all(B$desempleo_nsa > 0))
# Control NSA vs SA del desempleo: el promedio de 12 meses debe coincidir (factores estacionales ~ suman cero)
k <- 12:nrow(B); dif <- sapply(k, function(i) mean(B$desempleo_nsa[(i - 11):i]) - mean(B$desempleo_sa[(i - 11):i]))
cat(sprintf("Base mensual: %s a %s, %d meses (muestra principal 2000-2019: %d meses)\n", B$mes[1], B$mes[nrow(B)], nrow(B),
            sum(B$mes <= "2019-12")))
cat(sprintf("Control desempleo NSA vs SA, promedios de 12 meses: dif. media %.3f pp, max |dif| %.3f pp\n", mean(dif), max(abs(dif))))
cat(sprintf("Imputado: %s (desempleo NSA = %.2f)\n", B$mes[falta], B$desempleo_nsa[falta]))
write.csv(B, "datos/base_mensual.csv", row.names = FALSE)
print(head(B, 3)); print(B[B$mes %in% c("2001-08", "2001-09", "2020-03", "2020-04"), ]); print(tail(B, 2))
