# =============================================================================
# 02_exploracion.R  --  Exploracion de las series (paso 0 y 1 del ciclo Box-Jenkins)
# TP Econometria II, candidata 5: Indice Li Keqiang (China)
#
# Series (muestra principal 2000Q1-2023Q4, 96 trimestres):
#   Y = PBI oficial       X = Electricidad
# Dos tratamientos en paralelo (regla acordada, ver claude/reglas_TP.md):
#   (a) crecimiento interanual (%)           -> ya sin estacionalidad (diferencia 1-B^4)
#   (b) logaritmo de los niveles             -> conserva tendencia y estacionalidad
#
# Respaldo teorico (Notas propias):
#   - FAC/FACP, bandas de Bartlett y Quenouille, Box-Pierce y Ljung-Box: Nota 1; Nota 6, Sec. 9
#   - Dickey-Fuller (tau, tau_mu, tau_tau), ADF, Said-Dickey p* <= n^(1/3): Nota 6, Sec. 4
#   - ET vs ED, diferenciacion, sobrediferenciacion: Nota 6, Sec. 2 y 3
#   - Transformacion logaritmica estabilizadora de varianza: Nota 6, Sec. 5
#   - Idea rectora de la clase (Mirta, 14/10): el test no es una receta; el correlograma
#     "grueso" y los cambios estructurales deben leerse junto con el test.
#
# Ejecutar desde la raiz del repositorio:  Rscript scripts/02_exploracion.R
# Paquetes: urca  (install.packages("urca"))
# =============================================================================


install.packages("urca")
suppressPackageStartupMessages(library(urca))
dir.create("output/figuras", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tablas",  recursive = TRUE, showWarnings = FALSE)

# ---- 1. Datos ----------------------------------------------------------------
Q <- read.csv("output/trimestral_li_keqiang.csv", stringsAsFactors = FALSE)
Q$anio <- as.integer(substr(Q$trimestre, 1, 4))
Q$trim <- as.integer(substr(Q$trimestre, 6, 6))

# PBI real en nivel (indice): el PBI de chinadata.live es nominal (incluye precios) y el
# PBI real solo existe desde 2011. Se encadena la tasa interanual oficial (real) sobre el
# nivel de 1992 (anio base, nominal = real). Es una construccion propia y aproximada: la
# tasa oficial esta redondeada a 0,1 pp, por lo que el error de encadenado se acumula.
pbi_idx <- rep(NA_real_, nrow(Q))
pbi_idx[1:4] <- Q$pbi_nominal_1e8_yuan[1:4]
for (i in 5:nrow(Q)) pbi_idx[i] <- pbi_idx[i - 4] * (1 + Q$pbi_crec_interanual_oficial_pct[i] / 100)
Q$pbi_real_idx <- pbi_idx

# Control de la construccion. El PBI real publicado (chinadata.live, desde 2011) esta
# empalmado con cambios de anio base (saltos de nivel al pasar de 2015 a 2016, de 2020 a 2021
# y de 2025 a 2026), por lo que NO sirve comparar niveles. Lo que si debe coincidir es la tasa
# interanual: la tasa implicita en el PBI real publicado vs. la tasa oficial con la que se encadeno.
# Solo pueden diferir los trimestres cuyo "t" y "t-4" caen en bloques de base distinta
# (los primeros cuatro trimestres posteriores a cada cambio de base: 2016, 2021 y 2026).
k  <- which(!is.na(Q$pbi_real_1e8_yuan_2011))
g4 <- 100 * (Q$pbi_real_1e8_yuan_2011[k][-(1:4)] / Q$pbi_real_1e8_yuan_2011[k][1:(length(k) - 4)] - 1)
of <- Q$pbi_crec_interanual_oficial_pct[k][-(1:4)]
tr <- Q$trimestre[k][-(1:4)]
disc <- abs(g4 - of) > 0.15
cat(sprintf("Control PBI real: tasa interanual implicita en el publicado vs. oficial -> coinciden (<=0.15pp) en %d de %d trimestres.\n",
            sum(!disc), length(disc)))
cat("  Trimestres que no coinciden (cambios de base del publicado):", paste(tr[disc], collapse = ", "), "\n")
stopifnot(all(substr(tr[disc], 1, 4) %in% c("2016", "2021", "2026")))
# Limitacion conocida: el encadenado hereda la forma estacional del anio semilla (1992),
# fija en el tiempo; afecta solo al tratamiento en log-niveles.

# Muestra principal
S <- Q[Q$trimestre >= "2000Q1" & Q$trimestre <= "2023Q4", ]
stopifnot(nrow(S) == 96, !anyNA(S$pbi_crec_interanual_oficial_pct), !anyNA(S$elec_crec_interanual_pct))

mk <- function(x, start = c(2000, 1)) ts(x, start = start, frequency = 4)
series <- list(
  "PBI oficial, crec. interanual (%)"      = mk(S$pbi_crec_interanual_oficial_pct),
  "Electricidad, crec. interanual (%)"     = mk(S$elec_crec_interanual_pct),
  "ln PBI real (indice encadenado)"        = mk(log(S$pbi_real_idx)),
  "ln Electricidad (GWh, trimestral)"      = mk(log(S$electricidad_GWh))
)
tipo <- c("interanual", "interanual", "log-nivel", "log-nivel")
names(tipo) <- names(series)

# Eventos para leer los graficos (contexto, regla 5 de la clase)
eventos <- data.frame(t = c(2008.75, 2020.0), etiqueta = c("Crisis global 2008-09", "COVID 2020"))

# ---- 2. Graficos de evolucion --------------------------------------------------
png("output/figuras/02_evolucion.png", width = 1400, height = 1000, res = 130)
par(mfrow = c(2, 2), mar = c(3, 4, 3, 1))
for (nm in names(series)) {
  plot(series[[nm]], main = nm, ylab = "", xlab = "", col = "steelblue", lwd = 1.6)
  abline(v = eventos$t, lty = 2, col = "firebrick")
  if (tipo[nm] == "interanual") abline(h = 0, lty = 3, col = "grey50")
  mtext(paste(eventos$etiqueta, collapse = "  |  "), side = 1, line = 2, cex = 0.6, col = "firebrick")
}
dev.off()

# ---- 3. Correlogramas con bandas de Bartlett/Quenouille (+-1.96/sqrt(T)) --------
LAGS <- 24
png("output/figuras/02_correlogramas.png", width = 1600, height = 1700, res = 130)
par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))
for (nm in names(series)) {
  x <- as.numeric(series[[nm]])   # sin atributos ts: el eje queda en rezagos (trimestres), no en anios
  acf(x,  lag.max = LAGS, main = paste("FAC -",  nm), xlab = "Rezago (trimestres)")
  pacf(x, lag.max = LAGS, main = paste("FACP -", nm), xlab = "Rezago (trimestres)")
}
dev.off()

# ---- 4. Resumen de correlogramas y Ljung-Box -----------------------------------
resumen <- do.call(rbind, lapply(names(series), function(nm) {
  x <- series[[nm]]; n <- length(x); banda <- 1.96 / sqrt(n)
  r <- acf(x, lag.max = LAGS, plot = FALSE)$acf[-1]
  p <- pacf(x, lag.max = LAGS, plot = FALSE)$acf
  dentro <- which(abs(r) < banda)
  data.frame(
    serie = nm, tratamiento = tipo[nm], n = n, banda_95 = round(banda, 3),
    rho1 = round(r[1], 3), rho4 = round(r[4], 3), rho8 = round(r[8], 3), rho12 = round(r[12], 3),
    n_signif_FAC_de_24 = sum(abs(r) > banda),
    primer_rezago_FAC_en_banda = if (length(dentro)) dentro[1] else NA,
    n_signif_FACP_de_24 = sum(abs(p) > banda),
    FACP_signif_rezagos = paste(which(abs(p) > banda), collapse = ","),
    Q_LB_4  = round(Box.test(x, 4,  "Ljung-Box")$statistic, 1),
    p_LB_4  = signif(Box.test(x, 4,  "Ljung-Box")$p.value, 2),
    Q_LB_12 = round(Box.test(x, 12, "Ljung-Box")$statistic, 1),
    p_LB_12 = signif(Box.test(x, 12, "Ljung-Box")$p.value, 2),
    Q_LB_24 = round(Box.test(x, 24, "Ljung-Box")$statistic, 1),
    p_LB_24 = signif(Box.test(x, 24, "Ljung-Box")$p.value, 2),
    row.names = NULL)
}))
write.csv(resumen, "output/tablas/02_correlogramas_resumen.csv", row.names = FALSE)
cat("\n===== CORRELOGRAMAS y LJUNG-BOX (banda 95% = 1.96/sqrt(T)) =====\n")
print(resumen[, c("serie", "banda_95", "rho1", "rho4", "rho8", "n_signif_FAC_de_24",
                  "primer_rezago_FAC_en_banda", "FACP_signif_rezagos", "p_LB_4", "p_LB_12", "p_LB_24")],
      row.names = FALSE)

# ---- 5. Dickey-Fuller aumentado (3 variantes) -------------------------------------
# Rezago maximo por Said-Dickey: p* <= n^(1/3) (Nota 6, Sec. 4); el optimo por AIC.
# Regla de la Nota 6 / Dolado-Jenkinson-Sosvilla: partir de la especificacion mas general y
# descender; igual se reportan las tres para poder discutirlas.
adf_filas <- function(x, nombre, muestra) {
  n <- length(x); kmax <- floor(n^(1 / 3))
  do.call(rbind, lapply(c("none", "drift", "trend"), function(tp) {
    m  <- ur.df(x, type = tp, lags = kmax, selectlags = "AIC")
    cv <- m@cval[1, ]; tau <- m@teststat[1, 1]
    data.frame(serie = nombre, muestra = muestra, n = n, spec = c(none = "tau", drift = "tau_mu", trend = "tau_tau")[tp],
               rezagos_AIC = m@lags, estadistico = round(tau, 3),
               cv1 = cv["1pct"], cv5 = cv["5pct"], cv10 = cv["10pct"],
               rechaza_H0_raiz_unitaria_5pct = tau < cv["5pct"], row.names = NULL)
  }))
}

adf_todas <- list()
for (nm in names(series)) {
  adf_todas[[length(adf_todas) + 1]] <- adf_filas(series[[nm]], nm, "2000Q1-2023Q4")
}
# Transformaciones de los log-niveles: diferencia (ln PBI, ln Elec) y diferencia estacional
for (nm in names(series)[3:4]) {
  x <- series[[nm]]
  adf_todas[[length(adf_todas) + 1]] <- adf_filas(diff(x),       paste0("D ", nm),   "2000Q2-2023Q4")
  adf_todas[[length(adf_todas) + 1]] <- adf_filas(diff(x, 4),    paste0("D4 ", nm),  "2001Q1-2023Q4")
  adf_todas[[length(adf_todas) + 1]] <- adf_filas(diff(diff(x, 4)), paste0("D D4 ", nm), "2001Q2-2023Q4")
}

# Sensibilidad a los cambios estructurales (regla 7 y 9 de la clase): mismo test, otras muestras
for (nm in names(series)) {
  x <- series[[nm]]
  adf_todas[[length(adf_todas) + 1]] <- adf_filas(window(x, end = c(2019, 4)),   nm, "2000Q1-2019Q4 (sin COVID)")
  adf_todas[[length(adf_todas) + 1]] <- adf_filas(window(x, start = c(2010, 1)), nm, "2010Q1-2023Q4")
}
ADF <- do.call(rbind, adf_todas)
write.csv(ADF, "output/tablas/02_adf.csv", row.names = FALSE)

cat("\n===== DICKEY-FULLER AUMENTADO (H0: raiz unitaria; rechazo si estadistico < cv) =====\n")
print(ADF[, c("serie", "muestra", "n", "spec", "rezagos_AIC", "estadistico", "cv5", "rechaza_H0_raiz_unitaria_5pct")],
      row.names = FALSE)

# ---- 6. Estacionalidad (regla 8 de la clase) --------------------------------------
# Autocorrelaciones en rezagos estacionales (4, 8, 12) de la primera diferencia de los
# log-niveles: si son grandes y significativas, hay componente estacional (Mirta: "el 7 se
# repite en el 14 ... te da la pauta de que hay una componente estacional").
cat("\n===== ESTACIONALIDAD: FAC de la primera diferencia de los log-niveles =====\n")
for (nm in names(series)[3:4]) {
  d <- diff(series[[nm]]); r <- acf(d, lag.max = 12, plot = FALSE)$acf[-1]
  cat(sprintf("%-36s rho1=%.2f rho2=%.2f rho3=%.2f rho4=%.2f rho8=%.2f rho12=%.2f | banda=%.2f\n",
              paste("D", nm), r[1], r[2], r[3], r[4], r[8], r[12], 1.96 / sqrt(length(d))))
}
# Medias por trimestre del crecimiento trimestral (patron estacional de la electricidad)
dl <- 100 * diff(series[["ln Electricidad (GWh, trimestral)"]])
cat("Crecimiento trimestral medio de la electricidad (%) por trimestre del anio:\n")
print(round(tapply(dl, cycle(dl), mean), 2))

cat("\nListo. Figuras en output/figuras/, tablas en output/tablas/.\n")