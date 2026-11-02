# =============================================================================
# utils_a.R -- Herramientas del par 4 (se usa con source(), desde aereo_desempleo/)
# Reutiliza las herramientas mensuales de viviendas_tasa (version anterior, utils_m.R) (rejilla, ADF con dummies, recesiones, etc.) y
# solo redefine la carga de datos y los colores.
# =============================================================================
source("scripts/utils_m.R")
col_rpm <- "#2A62B8"; col_des <- "#D9701F"

cargar_base <- function() {
  B <- read.csv("datos/base_mensual.csv", stringsAsFactors = FALSE)
  ini <- c(as.integer(substr(B$mes[1], 1, 4)), as.integer(substr(B$mes[1], 6, 7)))
  list(B = B,
       rpm = ts(100 * log(B$rpm_miles), start = ini, frequency = 12),        # ln millas-pasajero (x100), NSA
       des = ts(B$desempleo_nsa, start = ini, frequency = 12),              # % de la fuerza laboral, NSA
       des_sa = ts(B$desempleo_sa, start = ini, frequency = 12),
       rec = ts(B$recesion, start = ini, frequency = 12))
}
# Efecto calendario: febrero de anio bisiesto tiene 29 dias (+3.6%). En una serie de flujo mensual
# (millas voladas en el mes) eso sube el dato; con la diferencia estacional el efecto no desaparece,
# porque febrero bisiesto se compara con un febrero de 28 dias. Regresor: 1 en febrero bisiesto.
bisiesto <- function(x) {
  t <- as.numeric(time(x)) + 1e-6; a <- floor(t); m <- floor((t - a) * 12) + 1
  ts(as.numeric(m == 2 & ((a %% 4 == 0 & a %% 100 != 0) | a %% 400 == 0)), start = start(x), frequency = 12)
}
