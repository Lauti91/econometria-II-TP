# =============================================================================
# datos_extra.R -- Series adicionales de FRED para las extensiones (scripts 15 a 19)
#   Se bajan UNA sola vez a datos/crudos/ (si el archivo ya existe no se vuelve a bajar, para no cambiar
#   el vintage) y la fecha queda en datos/crudos/fecha_descarga_bloqueC.txt.
#     MSACSR         meses de oferta de viviendas nuevas en venta (SA)
#     WPUSI012011    indice de precios al productor: materiales de construccion (NSA)
#     WPU081         indice de precios al productor: madera y derivados (NSA)
#     CES2000000003  salario horario medio de la construccion (desde 2006)
#     CSUSHPISA      indice Case-Shiller nacional de precios de viviendas (SA, desde 1987)
#     HOUST1F, HOUST2F, HOUST5F   inicios: unifamiliares, 2-4 unidades, 5 o mas (SA, tasa anual, miles)
#     PERMIT1, PERMIT5            permisos: unifamiliares, 5 o mas (SA, tasa anual, miles)
# Se carga con source("R/datos_extra.R") despues de R/funciones.R.
# =============================================================================
IDS_EXTRA <- c("MSACSR", "WPUSI012011", "WPU081", "CES2000000003", "CSUSHPISA", "HOUST1F", "HOUST2F", "HOUST5F", "PERMIT1", "PERMIT5")
bajar_extra <- function(ids = IDS_EXTRA) {
  falta <- ids[!file.exists(sprintf("datos/crudos/%s.csv", ids))]
  if (!length(falta)) return(invisible(NULL))
  for (id in falta) { d <- read.csv(sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id), na.strings = ".")
    names(d) <- c("fecha", "valor"); write.csv(d[!is.na(d$valor), ], sprintf("datos/crudos/%s.csv", id), row.names = FALSE) }
  cat(sprintf("Descarga FRED (%s): %s\n", paste(falta, collapse = ", "), format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
      file = "datos/crudos/fecha_descarga_bloqueC.txt", append = TRUE)
}
# Serie mensual alineada con los meses de una serie de referencia (ts). Los meses sin dato quedan en NA:
# cada script controla que no haya NA dentro de su ventana de estimacion.
leer_extra <- function(id, ref) {
  d <- read.csv(sprintf("datos/crudos/%s.csv", id)); v <- setNames(d$valor, substr(d$fecha, 1, 7))[etiqueta_mes(ref)]
  ts(as.numeric(v), start = start(ref), frequency = frequency(ref))
}
