# =============================================================================
# 12_episodio_2022.R -- El episodio 2022-2024: cuanto debian caer los inicios con la suba de tasas?
#   Entre diciembre de 2021 y octubre de 2023 la tasa hipotecaria a 30 anios paso de ~3,1 % a ~7,6 %.
#   Ejercicio fuera de la muestra: se estima la funcion de transferencia con datos HASTA 2021-12 y se
#   pronostican los 36 meses de 2022-01 a 2024-12 de tres maneras:
#     a) Pronostico condicional: con la trayectoria OBSERVADA de la tasa como input.
#     b) Contrafactico: con la tasa congelada en su nivel de 2021-12 (X_t = 0 desde 2022-01).
#     c) Univariado: el mismo modelo de ruido, sin input.
#   a) - b) es la caida de los inicios atribuible a la suba de tasas segun el modelo. Se compara con lo observado.
#   Herramientas y fuentes:
#     - Funcion de transferencia e intervencion: apunte parte 2, pp. 1-21. Forma FIJA, la elegida en 09 con la
#       muestra completa: b = 2, s = 0, r = 1, phi = 0,75; ruido ARMA(1,2)(1,1)12; COVID como TC (delta de 08).
#       (La forma se eligio mirando toda la muestra; los parametros se reestiman solo con datos hasta 2021-12.
#       Como control, phi tambien se vuelve a elegir por perfil de verosimilitud con la submuestra.)
#     - Pronostico con funcion de transferencia y input conocido, y varianza del error de pronostico con los
#       pesos psi del ruido: Box, G., Jenkins, G., Reinsel, G. y Ljung, G. (2015), Time Series Analysis:
#       Forecasting and Control, 5.a ed., Wiley, cap. 12 (sec. 12.5); pesos psi y pronostico ARIMA: apunte parte 1.
#       Las bandas son condicionales al input y a los parametros estimados (no incluyen su incertidumbre).
#     - Modelo de lineas aereas ARIMA(0,1,1)(0,1,1)12 de la serie sin desestacionalizar (script 04; Box et al.
#       2015, cap. 9), como segunda referencia univariada: se compara por totales anuales.
#   Los niveles se obtienen con exp() del pronostico en logaritmos (mediana, sin correccion lognormal).
#   Salidas: tablas 12_episodio.tex, 12_totales.tex, 12_parametros.tex, 12_trayectorias.csv;
#   figura 12_episodio.png; valores_12.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 12_episodio_2022.R   (despues de 08)
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
B <- 2; PHI <- 0.75; RUIDO <- c(1, 2, 1, 1); INI <- 7; H <- 36           # forma fija (09) y horizonte
Y <- diff(L$viv_sa); X <- diff(L$tasa); nY <- length(Y); mes <- etiqueta_mes(Y)
TC <- as.numeric(diff(ts(as.numeric(stats::filter(indicadora(L$viv_sa, 2020, 4), I8$delta, method = "recursive")), start = start(L$viv_sa), frequency = S)))
fin <- which(mes == "2021-12"); fut <- fin + 1:H                          # origen del pronostico y meses pronosticados
stopifnot(mes[fut[H]] == "2024-12")
# Regresor del numerador: X_{t-b} filtrado por 1/(1 - phi B)
regresor <- function(Xv, phi) { xf <- as.numeric(stats::filter(Xv, phi, method = "recursive")); c(rep(0, B), xf[seq_len(nY - B)]) }
serie <- function(i) ts(as.numeric(Y)[i], start = time(Y)[i[1]], frequency = S)
estimar <- function(hasta, phi, input = TRUE) {
  xr <- cbind(D2 = regresor(X, phi), COVID = TC)[INI:hasta, if (input) 1:2 else 2, drop = FALSE]
  ajustar(serie(INI:hasta), RUIDO[1], RUIDO[2], RUIDO[3], RUIDO[4], xr, TRUE)
}

# ---- 1. Estimacion hasta 2021-12 y comparacion con la muestra completa ---------------------------------------
m_sub <- estimar(fin, PHI); m_tot <- estimar(nY, PHI); m_uni <- estimar(fin, PHI, input = FALSE)
grilla <- seq(0.30, 0.95, by = 0.05)
perfil <- sapply(grilla, function(ph) { m <- estimar(fin, ph); if (is.null(m)) NA else as.numeric(logLik(m)) })
phi_sub <- grilla[which.max(perfil)]; m_phi <- estimar(fin, phi_sub)
LR_phi <- 2 * (max(perfil, na.rm = TRUE) - as.numeric(logLik(m_sub)))        # H0: phi = 0,75; chi2(1), VC 5 % = 3,84
cat(sprintf("phi por perfil con datos hasta 2021-12: %.2f (muestra completa: %.2f); LR de phi = 0,75: %.2f (p = %.2f)\n", phi_sub, PHI, LR_phi, pchisq(LR_phi, 1, lower.tail = FALSE)))
par_tab <- function(m, phi) { ee <- sqrt(diag(m$var.coef)); cf <- m$coef
  c(delta = cf[["D2"]], delta_ee = ee[["D2"]], gan = cf[["D2"]] / (1 - phi), gan_ee = ee[["D2"]] / (1 - phi), cf[c("ar1", "ma1", "ma2", "sar1", "sma1", "COVID")], sigma = sqrt(m$sigma2), n = length(residuals(m))) }
P <- rbind(completa = par_tab(m_tot, PHI), sub = par_tab(m_sub, PHI), sub_phi = par_tab(m_phi, phi_sub))
print(round(P, 3))
ult <- L$tasa_pp[length(L$tasa_pp)]; pp <- 100 * log((ult + 1) / ult)       # +1 pp sobre la tasa de 2026-08, como en 09
filaP <- function(i, etq, phi) data.frame(etq, P[i, "n"], num(phi, 2), sprintf("%s (%s)", num(P[i, "delta"], 3), num(P[i, "delta_ee"], 3)),
                                          sprintf("%s (%s)", num(P[i, "gan"], 2), num(P[i, "gan_ee"], 2)), num(P[i, "gan"] * pp, 1),
                                          num(P[i, "ar1"], 2), num(P[i, "sar1"], 2), num(P[i, "COVID"], 1), num(P[i, "sigma"], 2))
tabla_tex(rbind(filaP(1, "1985--2026 (script 09)", PHI), filaP(2, "1985--2021, $\\phi$ fijo", PHI), filaP(3, "1985--2021, $\\phi$ por perfil", phi_sub)),
          "12_parametros.tex", rojo = celdas(2, 1:10), "lccccccccc",
          c("Muestra", "$n$", "$\\phi$", "$\\hat\\delta_0$ (ee)", "Ganancia $g$ (ee cond.)", "\\% por +1 pp", "AR(1)", "SAR(12)", "$\\hat\\omega$ COVID", "$\\hat\\sigma$"))
D <- diagnostico(m_sub, sum(RUIDO))

# ---- 2. Pronosticos desde 2021-12 ------------------------------------------------------------------------------
X_cf <- X; X_cf[fut[1]:nY] <- 0                                            # tasa congelada desde 2022-01
pronostico <- function(m, phi, Xv, input = TRUE) {
  nx <- cbind(D2 = regresor(Xv, phi), COVID = TC)[fut, if (input) 1:2 else 2, drop = FALSE]
  xreg <- nx                                                               # predict.Arima busca "xreg" en este entorno (solo mira sus columnas)
  as.numeric(predict(m, n.ahead = H, newxreg = nx)$pred)
}
nivel0 <- as.numeric(L$viv_sa)[fin + 1]                                    # 100 ln(inicios) en 2021-12 (Y_t = nivel_t - nivel_{t-1})
stopifnot(etiqueta_mes(L$viv_sa)[fin + 1] == "2021-12")
camino <- function(d) nivel0 + cumsum(d)
obs <- as.numeric(L$viv_sa)[fin + 1 + 1:H]
cond <- camino(pronostico(m_sub, PHI, X)); contra <- camino(pronostico(m_sub, PHI, X_cf)); uni <- camino(pronostico(m_uni, PHI, X, FALSE))
cond_phi <- camino(pronostico(m_phi, phi_sub, X)); contra_phi <- camino(pronostico(m_phi, phi_sub, X_cf))
# Bandas del nivel: Var[e(h)] = sigma^2 * sum_{j<h} (psi*_j)^2, con psi* los pesos psi acumulados (el nivel es la suma de Y)
psi <- cumsum(c(1, ARMAtoMA(ar = m_sub$model$phi, ma = m_sub$model$theta, lag.max = H - 1)))
ee_h <- sqrt(m_sub$sigma2 * cumsum(psi^2))
# Efecto de la tasa = condicional - contrafactico = delta_0 * (suma de los regresores); su ee, condicional a phi
suma_reg <- cumsum(regresor(X, PHI)[fut] - regresor(X_cf, PHI)[fut])
efecto <- cond - contra; efecto_ee <- abs(suma_reg) * sqrt(m_sub$var.coef["D2", "D2"])
stopifnot(max(abs(efecto - m_sub$coef[["D2"]] * suma_reg)) < 1e-6)
# Lineas aereas sobre la serie sin desestacionalizar, estimado hasta 2021-12 (como en 04/08, con el TC del COVID)
x <- L$viv; nx <- length(x); z <- ts(as.numeric(stats::filter(indicadora(x, 2020, 4), I8$delta, method = "recursive")), start = start(x), frequency = S)
Wn <- diff(diff(x, 12)); Zn <- as.numeric(diff(diff(z, 12))); i_fin <- which(etiqueta_mes(Wn) == "2021-12")
m_air <- ajustar(ts(as.numeric(Wn)[1:i_fin], start = start(Wn), frequency = S), 0, 1, 0, 1, cbind(COVID = Zn[1:i_fin]), FALSE)
xreg <- cbind(COVID = Zn[i_fin + 1:H]); w_hat <- as.numeric(predict(m_air, n.ahead = H, newxreg = xreg)$pred)
xa <- as.numeric(x)[1:(i_fin + 13)]                                        # nivel hasta 2021-12; despues, x_t = w_t + x_{t-1} + x_{t-12} - x_{t-13}
for (h in 1:H) { t <- i_fin + 13 + h; xa[t] <- w_hat[h] + xa[t - 1] + xa[t - 12] - xa[t - 13] }
aire <- xa[i_fin + 13 + 1:H]; obs_nsa <- as.numeric(x)[i_fin + 13 + 1:H]
stopifnot(etiqueta_mes(x)[i_fin + 13] == "2021-12")

TR <- data.frame(mes = mes[fut], h = 1:H, tasa_pp = as.numeric(L$tasa_pp)[fin + 1 + 1:H], observado = obs, condicional = cond, cond_lo = cond - 1.96 * ee_h, cond_hi = cond + 1.96 * ee_h,
                 contrafactico = contra, univariado = uni, efecto_tasa = efecto, efecto_ee = efecto_ee, observado_nsa = obs_nsa, lineas_aereas = aire)
write.csv(TR, "salidas/tablas/12_trayectorias.csv", row.names = FALSE)

# ---- 3. Cuadros: promedios anuales respecto de 2021 y totales anuales -----------------------------------------
anio <- as.integer(substr(TR$mes, 1, 4)); base21 <- mean(as.numeric(L$viv_sa)[etiqueta_mes(L$viv_sa) >= "2021-01" & etiqueta_mes(L$viv_sa) <= "2021-12"])
pct <- function(d) 100 * (exp(d / 100) - 1)                               # de puntos logaritmicos a %
AN <- do.call(rbind, lapply(2022:2024, function(a) { i <- anio == a
  data.frame(anio = a, tasa = mean(TR$tasa_pp[i]), obs = pct(mean(obs[i]) - base21), cond = pct(mean(cond[i]) - base21), contra = pct(mean(contra[i]) - base21),
             uni = pct(mean(uni[i]) - base21), efecto = pct(mean(efecto[i])), ef_lo = pct(mean(efecto[i]) - 1.96 * mean(efecto_ee[i])), ef_hi = pct(mean(efecto[i]) + 1.96 * mean(efecto_ee[i])),
             efecto_phi = pct(mean((cond_phi - contra_phi)[i])), cond_phi = pct(mean(cond_phi[i]) - base21),
             T_obs = sum(exp(obs[i] / 100)), T_cond = sum(exp(cond[i] / 100)), T_contra = sum(exp(contra[i] / 100)), T_uni = sum(exp(uni[i] / 100)),
             T_obs_nsa = sum(exp(obs_nsa[i] / 100)), T_aire = sum(exp(aire[i] / 100))) }))
print(round(AN, 1), row.names = FALSE)
tasa21 <- mean(as.numeric(L$tasa_pp)[etiqueta_mes(L$tasa_pp) >= "2021-01" & etiqueta_mes(L$tasa_pp) <= "2021-12"])
tabla_tex(data.frame(AN$anio, num(AN$tasa, 2), num(AN$obs, 1), num(AN$cond, 1), num(AN$contra, 1), num(AN$uni, 1), sprintf("%s (%s; %s)", num(AN$efecto, 1), num(AN$ef_lo, 1), num(AN$ef_hi, 1))),
          "12_episodio.tex", rojo = celdas(1:3, c(3, 7)), "ccccccc",
          c("A\\~no", "Tasa media (\\%)", "Observado", "FT con la tasa observada", "FT con la tasa congelada", "Univariado", "Efecto de la tasa (IC 95\\,\\%)"))
tabla_tex(data.frame(AN$anio, num(AN$T_obs, 0), num(AN$T_cond, 0), num(AN$T_contra, 0), num(AN$T_uni, 0), num(AN$T_obs_nsa, 0), num(AN$T_aire, 0)),
          "12_totales.tex", rojo = celdas(1:3, 2:3), "ccccccc",
          c("A\\~no", "Observado (SA)", "FT, tasa observada", "FT, tasa congelada", "Univariado (SA)", "Observado (NSA)", "L\\'ineas a\\'ereas (NSA)"))
# Error de pronostico en los 36 meses (en 100 ln)
rmse <- function(p) sqrt(mean((obs - p)^2)); err <- c(cond = rmse(cond), contra = rmse(contra), uni = rmse(uni), aire = sqrt(mean((obs_nsa - aire)^2)))
sesgo <- c(cond = mean(obs - cond), contra = mean(obs - contra), uni = mean(obs - uni))
cat("RMSE 36 meses (100 ln):", paste(names(err), round(err, 1), collapse = " | "), "\nSesgo medio (obs - pron):", paste(names(sesgo), round(sesgo, 1), collapse = " | "), "\n")
cat("Meses observados dentro de la banda del 95 %:", sum(obs >= TR$cond_lo & obs <= TR$cond_hi), "de", H, "\n")

# ---- Valores para el texto ---------------------------------------------------------------------------------
for (a in 2022:2024) { r <- AN[AN$anio == a, ]; s <- c(`2022` = "Uno", `2023` = "Dos", `2024` = "Tres")[as.character(a)]
  guardar_valor(paste0("EpObs", s), num(r$obs, 1)); guardar_valor(paste0("EpCond", s), num(r$cond, 1)); guardar_valor(paste0("EpContra", s), num(r$contra, 1))
  guardar_valor(paste0("EpUni", s), num(r$uni, 1)); guardar_valor(paste0("EpEfecto", s), num(r$efecto, 1)); guardar_valor(paste0("EpEfectoLo", s), num(r$ef_lo, 1)); guardar_valor(paste0("EpEfectoHi", s), num(r$ef_hi, 1))
  guardar_valor(paste0("EpEfectoPhi", s), num(r$efecto_phi, 1)); guardar_valor(paste0("EpTasa", s), num(r$tasa, 2))
  guardar_valor(paste0("EpTotObs", s), num(r$T_obs, 0)); guardar_valor(paste0("EpTotCond", s), num(r$T_cond, 0)); guardar_valor(paste0("EpTotContra", s), num(r$T_contra, 0)) }
for (h in c(12, 24, 36)) { s <- c(`12` = "Doce", `24` = "VeinteCuatro", `36` = "TreintaSeis")[as.character(h)]
  guardar_valor(paste0("EpEfectoMes", s), num(pct(efecto[h]), 1)) }
tri <- anio %in% 2022:2024
guardar_valor("EpObsTrienio", num(pct(mean(obs) - base21), 1)); guardar_valor("EpCondTrienio", num(pct(mean(cond) - base21), 1))
guardar_valor("EpContraTrienio", num(pct(mean(contra) - base21), 1)); guardar_valor("EpEfectoTrienio", num(pct(mean(efecto)), 1))
guardar_valor("EpEfectoMax", num(pct(min(efecto)), 1)); guardar_valor("EpEfectoMaxMes", TR$mes[which.min(efecto)])
guardar_valor("EpTasaIni", num(as.numeric(L$tasa_pp)[fin + 1], 2)); guardar_valor("EpTasaMax", num(max(TR$tasa_pp), 2)); guardar_valor("EpTasaMaxMes", TR$mes[which.max(TR$tasa_pp)])
guardar_valor("EpTasaBase", num(tasa21, 2))
guardar_valor("EpGanSub", num(P["sub", "gan"], 2)); guardar_valor("EpGanSubEE", num(P["sub", "gan_ee"], 2)); guardar_valor("EpGanSubPP", num(P["sub", "gan"] * pp, 1))
guardar_valor("EpGanTot", num(P["completa", "gan"], 2)); guardar_valor("EpPhiSub", num(phi_sub, 2)); guardar_valor("EpGanPhiSub", num(P["sub_phi", "gan"], 2))
# Que fraccion de la caida que el modelo atribuye a la tasa se observo (2023-2024, respecto del contrafactico)
i23 <- anio >= 2023; fraccion <- mean(obs[i23] - contra[i23]) / mean(efecto[i23])
cat(sprintf("Caida observada respecto del contrafactico / caida atribuida a la tasa (2023-24): %.2f\n", fraccion))
guardar_valor("EpFraccion", num(100 * fraccion, 0))
guardar_valor("EpPhiLR", num(LR_phi, 2)); guardar_valor("EpPhiLRp", num(pchisq(LR_phi, 1, lower.tail = FALSE), 2))
guardar_valor("EpRMSECond", num(err["cond"], 1)); guardar_valor("EpRMSEContra", num(err["contra"], 1)); guardar_valor("EpRMSEUni", num(err["uni"], 1))
guardar_valor("EpSesgoCond", num(sesgo["cond"], 1)); guardar_valor("EpSesgoContra", num(sesgo["contra"], 1))
guardar_valor("EpDentroBanda", sum(obs >= TR$cond_lo & obs <= TR$cond_hi)); guardar_valor("EpLBSub", num(D$LB24, 2)); guardar_valor("EpNSub", D$n)

# ---- Figura ------------------------------------------------------------------------------------------------
# Arriba: la tasa observada y la congelada. Abajo: inicios (miles, tasa anual desestacionalizada) observados y pronosticados.
anual <- function(v) 12 * exp(v / 100); t_f <- as.numeric(time(Y))[fut]
hist <- data.frame(t = as.numeric(time(L$viv_sa)), v = anual(as.numeric(L$viv_sa)), tasa = as.numeric(L$tasa_pp)); hist <- hist[hist$t >= 2019 & hist$t < 2025, ]
PR <- rbind(data.frame(t = t_f, v = anual(cond), s = "FT con la tasa observada"), data.frame(t = t_f, v = anual(contra), s = "FT con la tasa congelada"),
            data.frame(t = t_f, v = anual(uni), s = "Univariado"))
colores <- c("FT con la tasa observada" = rojo, "FT con la tasa congelada" = col_tasa, "Univariado" = "#7B8794")
t0 <- as.numeric(time(L$viv_sa))[fin + 1]
g1 <- ggplot(hist, aes(t, tasa)) + geom_vline(xintercept = t0, colour = tinta2, linetype = "22", linewidth = 0.3) +
  geom_line(colour = col_tasa, linewidth = 0.7) + annotate("segment", x = t0, xend = max(t_f), y = hist$tasa[which.min(abs(hist$t - t0))], yend = hist$tasa[which.min(abs(hist$t - t0))], colour = col_tasa, linetype = "22", linewidth = 0.5) +
  annotate("text", x = max(t_f), y = hist$tasa[which.min(abs(hist$t - t0))] + 0.35, label = "congelada en dic-2021 (contrafactico)", hjust = 1, size = 3, colour = tinta2) +
  labs(title = "Tasa hipotecaria a 30 anios", subtitle = "% anual; la linea vertical marca el origen del pronostico (diciembre de 2021)", x = NULL, y = NULL) + tema
g2 <- ggplot() + geom_vline(xintercept = t0, colour = tinta2, linetype = "22", linewidth = 0.3) +
  geom_ribbon(data = data.frame(t = t_f, lo = anual(TR$cond_lo), hi = anual(TR$cond_hi)), aes(t, ymin = lo, ymax = hi), fill = rojo, alpha = 0.10) +
  geom_line(data = hist, aes(t, v), colour = tinta, linewidth = 0.55) +
  geom_line(data = PR, aes(t, v, colour = s, linetype = s), linewidth = 0.8) +
  scale_colour_manual(values = colores, name = NULL) + scale_linetype_manual(values = c("FT con la tasa observada" = "solid", "FT con la tasa congelada" = "solid", "Univariado" = "22"), name = NULL) +
  annotate("text", x = 2019.1, y = max(hist$v) * 0.99, label = "observado", hjust = 0, size = 3, colour = tinta) +
  labs(title = "Inicios de viviendas: observado y pronosticos desde diciembre de 2021", subtitle = "miles de unidades, tasa anual desestacionalizada; banda: 95 % del pronostico condicional (input y parametros dados)", x = NULL, y = NULL) +
  tema + theme(legend.position = "bottom")
guardar_fig(g1 / g2 + plot_layout(heights = c(1, 2.1)), "12_episodio.png", 10, 7.2)
escribir_valores("valores_12.tex")
