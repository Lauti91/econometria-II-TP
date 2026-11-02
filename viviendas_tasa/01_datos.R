# =============================================================================
# 01_datos.R -- Arma la base mensual a partir de los CSV de FRED guardados en datos/crudos/
#   HOUSTNSA      inicios de viviendas privadas, miles de unidades, mensual, SIN desestacionalizar
#   HOUST         idem desestacionalizada (tasa anual, miles): solo para comparar (regla 8)
#   MORTGAGE30US  tasa hipotecaria fija a 30 anios (Freddie Mac, PMMS), %, SEMANAL
#   GS10          bono del Tesoro a 10 anios, %, mensual (contexto)
#   USREC         recesiones NBER (1 = mes en recesion)
# La tasa semanal se pasa a mensual con el promedio de las semanas fechadas en cada mes (como FRED).
# Por defecto NO vuelve a descargar: usa los crudos del 19-oct-2026, para que los resultados del
# informe sean reproducibles. Con DESCARGAR <- TRUE baja todo de nuevo desde FRED.
# Ejecutar desde viviendas_tasa/:  Rscript 01_datos.R
# =============================================================================
DESCARGAR <- FALSE
ids <- c("HOUSTNSA", "HOUST", "MORTGAGE30US", "GS10", "USREC")
dir.create("datos/crudos", showWarnings = FALSE, recursive = TRUE)
if (DESCARGAR) {
  for (id in ids) {
    d <- read.csv(sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id), na.strings = ".")
    names(d) <- c("fecha", "valor"); write.csv(d[!is.na(d$valor), ], sprintf("datos/crudos/%s.csv", id), row.names = FALSE)
  }
  writeLines(sprintf("Descarga FRED: %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "datos/crudos/fecha_descarga.txt")
}
crudos <- lapply(setNames(ids, ids), function(id) { d <- read.csv(sprintf("datos/crudos/%s.csv", id)); d$fecha <- as.Date(d$fecha); d })

# ---- Tasa hipotecaria: semanal -> mensual ---------------------------------------------------------
H <- crudos$MORTGAGE30US; H$mes <- format(H$fecha, "%Y-%m")
tasa <- aggregate(valor ~ mes, H, mean); tasa$semanas <- aggregate(valor ~ mes, H, length)$valor
names(tasa)[2] <- "tasa_hip_30"
mensual <- function(id, nombre) setNames(data.frame(format(crudos[[id]]$fecha, "%Y-%m"), crudos[[id]]$valor), c("mes", nombre))
B <- Reduce(function(a, b) merge(a, b, by = "mes", all.x = TRUE),
            list(mensual("HOUSTNSA", "inicios_nsa"), mensual("HOUST", "inicios_sa_anual"), tasa,
                 mensual("GS10", "bono10"), mensual("USREC", "recesion")))
# Muestra comun: desde el primer mes completo de la tasa (1971-04) hasta el ultimo dato de inicios,
# sin el mes en curso de la descarga (la tasa de ese mes estaria incompleta)
descarga <- as.Date(sub("Descarga FRED: ([0-9-]+).*", "\\1", readLines("datos/crudos/fecha_descarga.txt")))
fin <- max(B$mes[!is.na(B$inicios_nsa)])
B <- B[B$mes >= "1971-04" & B$mes <= fin & B$mes < format(descarga, "%Y-%m"), ]

# ---- Controles ------------------------------------------------------------------------------------
esperados <- format(seq(as.Date(paste0(B$mes[1], "-01")), as.Date(paste0(B$mes[nrow(B)], "-01")), by = "month"), "%Y-%m")
stopifnot(identical(B$mes, esperados), !anyNA(B[, -1]), all(B$semanas %in% 4:5), all(B$inicios_nsa > 0), all(B$tasa_hip_30 > 0))
# La suma de 12 meses de la serie NSA y el promedio de 12 meses de la SA (tasa anual) deben ser parecidos
k <- 12:nrow(B)
dif <- 100 * (sapply(k, function(i) sum(B$inicios_nsa[(i - 11):i])) / sapply(k, function(i) mean(B$inicios_sa_anual[(i - 11):i])) - 1)
cat(sprintf("Base mensual: %s a %s (%d meses). Control NSA vs SA: dif. media %.2f%%, max |dif| %.2f%%\n",
            B$mes[1], B$mes[nrow(B)], nrow(B), mean(dif), max(abs(dif))))
write.csv(B, "datos/base_mensual.csv", row.names = FALSE)
