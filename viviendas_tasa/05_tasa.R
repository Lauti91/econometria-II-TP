# =============================================================================
# 05_tasa.R -- Ciclo Box-Jenkins de y = 100 ln(tasa hipotecaria), 1985-01 a 2026-08
#   Filtro: (1-B) (03: ADF y FAC del nivel). Sin estacionalidad (02: F de dummies no rechaza).
#   Identificacion: rejilla ARMA(p<=3, q<=2), con y sin constante; admisibilidad y SBC (Nota 6, Sec. 8-9)
#   Diagnostico: Ljung-Box, Jarque-Bera, cuadrados, sobreajuste (Sec. 10); estabilidad por subperiodos
# Ejecutar desde viviendas_tasa/:  Rscript 05_tasa.R
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); y <- L$tasa; W <- diff(y)

# ---- Identificacion ---------------------------------------------------------------------------------
R <- rbind(transform(rejilla(W, NULL, TRUE, 3, 2, 0, 0), constante = "s\\'i"), transform(rejilla(W, NULL, FALSE, 3, 2, 0, 0), constante = "no"))
R <- R[order(R$SBC), ]; R$dSBC <- R$SBC - min(R$SBC)
top <- head(R, 8)
tabla_tex(data.frame(top$modelo, top$constante, num(top$sigma, 2), num(top$LB12, 2), num(top$LB24, 2), num(top$LB36, 2),
                     ifelse(top$invertible & top$estacionario, "s\\'i", "no"), num(top$AICn, 4), num(top$SCn, 4), num(top$dSBC, 1)),
          "05_rejilla.tex", rojo = celdas(1, 1:10), "lccccccccc", c("Modelo", "Const.", "$\\hat\\sigma$", "LB12", "LB24", "LB36", "Estac./invert.", "AIC", "SC", "$\\Delta$SBC"))
g <- R[R$admisible, ][1, ]; print(head(R[, c("modelo", "constante", "SBC", "dSBC", "LB12", "LB24", "admisible")], 8))
stopifnot(g$p == 0, g$q == 1)
media <- g$constante == "s\\'i"

# ---- Estimacion y diagnostico -------------------------------------------------------------------------
m <- ajustar(W, 0, 1, 0, 0, NULL, media); k <- 1
guardar_salida(salida_arma(m, W, "D(LMORT)", "Tasa: ARIMA(0,1,1)"), "05_modelo_tasa.png")
saveRDS(list(m = m, W = W), "salidas/tablas/05_modelo.rds")
cf <- m$coef; se <- sqrt(diag(m$var.coef)); D <- diagnostico(m, k)
guardar_valor("ThetaT", num(-cf["ma1"], 3)); guardar_valor("ThetaTee", num(se["ma1"], 3)); guardar_valor("ConstT", if (media) "con constante" else "sin constante")
guardar_valor("SigmaT", num(sqrt(m$sigma2), 2)); guardar_valor("NT", length(W))
guardar_valor("LBTdoce", num(D$LB12, 2)); guardar_valor("LBTveinticuatro", num(D$LB24, 2)); guardar_valor("LBTtreintayseis", num(D$LB36, 2))
guardar_valor("JBT", num(D$JB, 0)); guardar_valor("JBTp", pval(D$JB_p)); guardar_valor("LBcuadT", pval(D$LB2_24))
guardar_valor("RhoTeoT", num(cf["ma1"] / (1 + cf["ma1"]^2), 2))
e <- as.numeric(residuals(m)); es <- e / sqrt(m$sigma2); at <- which(abs(es) > 3.5)
guardar_valor("AtipicosT", paste(sprintf("%s (%s)", etiqueta_mes(W)[at], num(es[at], 1)), collapse = "; "))
guardar_valor("CurtosisT", num(jarque_bera(e, k)["curtosis"], 1))
guardar_salida(salida_correlograma(residuals(m), "Correlograma de los residuos: tasa", 36, k, rojo_k = c(12, 24, 36)), "05_corr_residuos.png")
SO <- sobreajuste(W, c(0, 1, 0, 0), NULL, media)
SO <- SO[SO$agregado %in% c("AR(1)", "MA(2)"), ]
tabla_tex(data.frame(SO$agregado, num(SO$coef, 3), num(SO$t, 2), num(SO$dSBC, 1)), "05_sobreajuste.tex", rojo = celdas(seq_len(nrow(SO)), 3), "lccc",
          c("Par\\'ametro agregado", "Estimaci\\'on", "$t$", "$\\Delta$SBC"))
# Residuos por subperiodo: la volatilidad cambia (heterocedasticidad, Nota 6 Sec. 10.2)
tt <- as.numeric(time(W))
fila <- function(nombre, i) { mm <- ajustar(ts(as.numeric(W)[i], start = tt[which(i)[1]], frequency = 12), 0, 1, 0, 0, NULL, media)
  data.frame(nombre, sprintf("%s (%s)", num(-mm$coef["ma1"], 2), num(sqrt(mm$var.coef["ma1", "ma1"]), 2)), num(sqrt(mm$sigma2), 2), sum(i)) }
EST <- rbind(fila("1985--2026 (modelo elegido)", rep(TRUE, length(W))), fila("1985--2007", tt < 2008), fila("2008--2026", tt >= 2008))
tabla_tex(EST, "05_estabilidad.tex", rojo = celdas(2:3, 3), "lccc", c("Muestra", "$\\hat\\theta$ (ee)", "$\\hat\\sigma$", "$n$"))
print(EST)
# Figura de residuos
C <- correlograma(e, 36)
g1 <- ggplot(data.frame(t = tt, e = es), aes(t, e)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3, fill = gris, alpha = 0.3) +
  geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_line(colour = col_tasa, linewidth = 0.35) +
  geom_point(data = data.frame(t = tt, e = es)[abs(es) > 3.5, ], colour = rojo, size = 2.2) +
  scale_x_continuous(breaks = seq(1985, 2025, 10)) + labs(title = "Residuos estandarizados", subtitle = "franja: +-3", x = NULL, y = NULL) + tema
g2 <- ggplot(C, aes(k, AC)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -C$banda[1], ymax = C$banda[1], fill = gris, alpha = 0.45) +
  geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_segment(aes(xend = k, yend = 0, colour = abs(AC) > banda), linewidth = 1.1) +
  scale_colour_manual(values = c(`TRUE` = col_tasa, `FALSE` = "#9AA5B1"), guide = "none") + scale_x_continuous(breaks = c(1, 12, 24, 36)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) + labs(title = "FAC de los residuos", subtitle = sprintf("Ljung-Box p: %.2f / %.2f / %.2f (12/24/36)", D$LB12, D$LB24, D$LB36), x = "Rezago", y = NULL) + tema
g3 <- ggplot(data.frame(teo = qnorm(ppoints(length(es)))[rank(es)], mue = es), aes(teo, mue)) +
  geom_abline(slope = 1, intercept = 0, colour = tinta2, linetype = "22") + geom_point(colour = col_tasa, size = 0.7) +
  labs(title = "Q-Q normal", subtitle = sprintf("Jarque-Bera = %.0f (p < 0,001)", D$JB), x = "Cuantil teorico", y = NULL) + tema
guardar_fig(g1 | g2 | g3, "05_residuos.png", 12, 3.6)
escribir_valores("valores_05.tex")
