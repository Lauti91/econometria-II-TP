# =============================================================================
# 06_pronostico.R -- Evaluacion fuera de muestra (Nota 5, Sec. 9) y pronostico hacia adelante
#   Estimacion: 1985-01 a 2017-12. Evaluacion: 2018-01 a 2019-12 (antes del COVID), origen movil:
#   en cada mes se reestiman los parametros y se pronostica a h = 1, 3, 6 y 12. El ORDEN de cada
#   modelo se elige solo con la muestra de estimacion (Nota 5, Sec. 9: si no, la comparacion esta
#   contaminada) y queda fijo.
#   Inicios: lineas aereas (F3), el rival F2 con su constante, su promedio (Nota 6, Sec. 11) y dos
#   referencias ingenuas: ingenuo estacional x_{T+h-12} y paseo aleatorio estacional.
#   Tasa: ARIMA(0,1,1) contra el paseo aleatorio.
#   Metricas: RMSE y MAE en % (100 ln), MAPE y U de Theil sobre los niveles (miles de viviendas, %).
#   Al final: pronostico de 12 meses desde el ultimo dato (2026-08) con bandas del 95%.
# Ejecutar desde viviendas_tasa/:  Rscript 06_pronostico.R
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); x_all <- L$viv; y_all <- L$tasa
FIN_EST <- fecha_num(2017, 12); FIN_EV <- fecha_num(2019, 12); H <- c(1, 3, 6, 12)

# ---- Orden de cada modelo con la muestra de estimacion ---------------------------------------------------
x_est <- window(x_all, end = FIN_EST); d0 <- time(x_est)[14]
R3 <- rejilla(window(diff(diff(x_est, 12)), start = d0), NULL, FALSE); o3 <- unlist(R3[R3$admisible, ][1, c("p", "q", "P", "Q")])
R2 <- rejilla(window(diff(x_est, 12), start = d0), NULL, TRUE); o2 <- unlist(R2[R2$admisible, ][1, c("p", "q", "P", "Q")])
Rt <- rejilla(diff(window(y_all, end = FIN_EST)), NULL, FALSE, 3, 2, 0, 0); ot <- unlist(Rt[Rt$admisible, ][1, c("p", "q")])
cat("Orden elegido con 1985-2017: F3", o3, " F2", o2, " tasa", ot, "\n")
guardar_valor("OrdenFtres", R3[R3$admisible, ][1, "modelo"]); guardar_valor("OrdenFdos", R2[R2$admisible, ][1, "modelo"])
guardar_valor("OrdenTasa", Rt[Rt$admisible, ][1, "modelo"])

ttrend <- function(z) cbind(c12 = (seq_along(z)) / 12)          # (1 - B^12)(t/12) = 1: la constante de F2
ajus_niv <- function(z, tipo) suppressWarnings(switch(tipo,
  F3 = arima(z, order = c(o3[1], 1, o3[2]), seasonal = list(order = c(o3[3], 1, o3[4]), period = 12), method = "ML"),
  F2 = arima(z, order = c(o2[1], 0, o2[2]), seasonal = list(order = c(o2[3], 1, o2[4]), period = 12), xreg = ttrend(z), include.mean = FALSE, method = "ML"),
  T  = arima(z, order = c(ot[1], 1, ot[2]), method = "ML")))
pron <- function(m, z, tipo, h) if (tipo == "F2") predict(m, h, newxreg = cbind(c12 = (length(z) + 1:h) / 12)) else predict(m, h)

# ---- Origen movil -------------------------------------------------------------------------------------
tiempos <- as.numeric(time(x_all)); origenes <- tiempos[tiempos >= FIN_EST - 1e-6 & tiempos < FIN_EV - 1e-6]
filas <- list()
for (T0 in origenes) {
  x <- window(x_all, end = T0); y <- window(y_all, end = T0); n <- length(x); hm <- min(12, round((FIN_EV - T0) * 12))
  real <- as.numeric(window(x_all, start = T0 + 1 / 12, end = T0 + hm / 12)); real_t <- as.numeric(window(y_all, start = T0 + 1 / 12, end = T0 + hm / 12))
  p3 <- as.numeric(pron(ajus_niv(x, "F3"), x, "F3", hm)$pred); p2 <- as.numeric(pron(ajus_niv(x, "F2"), x, "F2", hm)$pred)
  P <- list("L\\'ineas a\\'ereas (F3)" = p3, "Rival F2" = p2, "Promedio F3 y F2" = (p3 + p2) / 2,
            "Ingenuo estacional" = as.numeric(x)[n - 12 + 1:hm],
            "Paseo aleatorio estacional" = as.numeric(x)[n] + as.numeric(x)[n - 12 + 1:hm] - as.numeric(x)[n - 12])
  for (nm in names(P)) filas[[length(filas) + 1]] <- data.frame(serie = "Inicios", modelo = nm, origen = T0, h = 1:hm, real = real, pron = P[[nm]])
  pt <- as.numeric(pron(ajus_niv(y, "T"), y, "T", hm)$pred)
  filas[[length(filas) + 1]] <- data.frame(serie = "Tasa", modelo = "ARIMA(0,1,1)", origen = T0, h = 1:hm, real = real_t, pron = pt)
  filas[[length(filas) + 1]] <- data.frame(serie = "Tasa", modelo = "Paseo aleatorio", origen = T0, h = 1:hm, real = real_t, pron = as.numeric(y)[n])
}
PR <- do.call(rbind, filas); PR$e <- PR$real - PR$pron
PR$Yr <- exp(PR$real / 100); PR$Yp <- exp(PR$pron / 100)                # niveles (miles de viviendas / % de tasa)
write.csv(PR, "salidas/tablas/06_pronosticos.csv", row.names = FALSE)
MET <- do.call(rbind, lapply(split(PR[PR$h %in% H, ], list(PR$serie[PR$h %in% H], PR$modelo[PR$h %in% H], PR$h[PR$h %in% H]), drop = TRUE), function(d)
  data.frame(serie = d$serie[1], modelo = d$modelo[1], h = d$h[1], N = nrow(d), RMSE = sqrt(mean(d$e^2)), MAE = mean(abs(d$e)),
             MAPE = 100 * mean(abs((d$Yr - d$Yp) / d$Yr)), U = sqrt(mean((d$Yr - d$Yp)^2)) / (sqrt(mean(d$Yr^2)) + sqrt(mean(d$Yp^2))))))
ref <- ifelse(MET$serie == "Inicios", "Paseo aleatorio estacional", "Paseo aleatorio")
MET$RMSE_rel <- MET$RMSE / MET$RMSE[match(paste(MET$serie, ref, MET$h), paste(MET$serie, MET$modelo, MET$h))]
write.csv(MET, "salidas/tablas/06_metricas.csv", row.names = FALSE)
orden_m <- c("L\\'ineas a\\'ereas (F3)", "Rival F2", "Promedio F3 y F2", "Ingenuo estacional", "Paseo aleatorio estacional", "ARIMA(0,1,1)", "Paseo aleatorio")
MI <- MET[MET$serie == "Inicios", ]; MI <- MI[order(match(MI$modelo, orden_m), MI$h), ]
ancho <- function(D, col, d) { W <- reshape(D[, c("modelo", "h", col)], idvar = "modelo", timevar = "h", direction = "wide"); W[, -1] <- lapply(W[, -1], num, d); W }
TI <- merge(ancho(MI, "RMSE", 1), ancho(MI, "MAPE", 1), by = "modelo", sort = FALSE); TI$U <- num(MI$U[MI$h == 1], 3)[match(TI$modelo, MI$modelo[MI$h == 1])]
TI <- TI[order(match(TI$modelo, orden_m)), ]
tabla_tex(TI, "06_metricas_inicios.tex", rojo = celdas(1, 1:10), "lccccccccc",
          c("Modelo", "RMSE $h{=}1$", "$h{=}3$", "$h{=}6$", "$h{=}12$", "MAPE $h{=}1$", "$h{=}3$", "$h{=}6$", "$h{=}12$", "$U$ ($h{=}1$)"))
MT <- MET[MET$serie == "Tasa", ]; MT <- MT[order(match(MT$modelo, orden_m), MT$h), ]
TT <- merge(ancho(MT, "RMSE", 2), ancho(MT, "MAPE", 2), by = "modelo", sort = FALSE)
tabla_tex(TT, "06_metricas_tasa.tex", rojo = celdas(1, 1:9), "lcccccccc", c("Modelo", "RMSE $h{=}1$", "$h{=}3$", "$h{=}6$", "$h{=}12$", "MAPE $h{=}1$", "$h{=}3$", "$h{=}6$", "$h{=}12$"))
print(MET[order(MET$serie, MET$h, MET$RMSE), ], digits = 3)
rel <- function(mod, h) MET$RMSE_rel[MET$modelo == mod & MET$h == h]
guardar_valor("MejoraUno", num(100 * (1 - rel(orden_m[1], 1)), 0)); guardar_valor("MejoraTres", num(100 * (1 - rel(orden_m[1], 3)), 0))
guardar_valor("MejoraSeis", num(100 * (1 - rel(orden_m[1], 6)), 0)); guardar_valor("MejoraDoce", num(100 * (1 - rel(orden_m[1], 12)), 0))
guardar_valor("MejoraTasaUno", num(100 * (1 - rel("ARIMA(0,1,1)", 1)), 0)); guardar_valor("MejoraTasaDoce", num(100 * (1 - rel("ARIMA(0,1,1)", 12)), 0))
guardar_valor("MapeFtresUno", num(MET$MAPE[MET$modelo == orden_m[1] & MET$h == 1], 1))
guardar_valor("MapeFtresDoce", num(MET$MAPE[MET$modelo == orden_m[1] & MET$h == 12], 1))
guardar_valor("RmseIngDoce", num(MET$RMSE[MET$modelo == "Ingenuo estacional" & MET$h == 12], 1))
guardar_valor("RmseFtresDoce", num(MET$RMSE[MET$modelo == orden_m[1] & MET$h == 12], 1))
guardar_valor("NEval", length(origenes))

# ---- Figuras: abanico desde 2017-12 (24 meses) y pronostico desde el ultimo dato -----------------------------
abanico <- function(z, h, tipo, xr_ao = NULL) {
  m <- if (is.null(xr_ao)) ajus_niv(z, tipo) else suppressWarnings(arima(z, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12), xreg = xr_ao, method = "ML"))
  P <- if (is.null(xr_ao)) pron(m, z, tipo, h) else predict(m, h, newxreg = matrix(0, h, 1))
  data.frame(t = as.numeric(time(P$pred)), m = exp(as.numeric(P$pred) / 100), lo = exp((as.numeric(P$pred) - 1.96 * as.numeric(P$se)) / 100),
             hi = exp((as.numeric(P$pred) + 1.96 * as.numeric(P$se)) / 100))
}
A1 <- abanico(x_est, 24, "F3"); obs1 <- window(x_all, start = 2015, end = FIN_EV)
A2 <- abanico(x_all, 12, "F3", xr_ao = indicadora(x_all, 2020, 4)); obs2 <- window(x_all, start = 2021)
panel <- function(A, obs, tit, sub) ggplot() + geom_ribbon(data = A, aes(t, ymin = lo, ymax = hi), fill = rojo, alpha = 0.10) +
  geom_line(aes(as.numeric(time(obs)), exp(as.numeric(obs) / 100)), colour = tinta, linewidth = 0.6) +
  geom_line(data = A, aes(t, m), colour = rojo, linewidth = 0.8, linetype = "22") + labs(title = tit, subtitle = sub, x = NULL, y = NULL) + tema
g <- panel(A1, obs1, "Pronostico desde 2017-12 (24 meses) contra lo observado", "negro: observado (miles de inicios); azul: pronostico y banda del 95%") /
  panel(A2, obs2, "Pronostico desde el ultimo dato (2026-08), 12 meses", "modelo estimado con toda la muestra (incluye el AO de abril de 2020)")
guardar_fig(g, "06_pronostico.png", 10, 7.5)
cob <- mean(window(x_all, start = FIN_EST + 1 / 12, end = FIN_EV) >= 100 * log(A1$lo) & window(x_all, start = FIN_EST + 1 / 12, end = FIN_EV) <= 100 * log(A1$hi))
guardar_valor("CoberturaAbanico", num(100 * cob, 0))
escribir_valores("valores_06.tex")
