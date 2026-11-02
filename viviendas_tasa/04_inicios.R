# =============================================================================
# 04_inicios.R -- Ciclo Box-Jenkins de x = 100 ln(inicios), 1985-01 a 2026-08 (Nota 6, Sec. 6-10)
#   1. Filtros rivales, comparados por SBC sobre las mismas fechas (desde 1986-02):
#        F1: (1-B)x + dummies mensuales      estacionalidad determinista
#        F2: (1-B^12)x + constante            estacionalidad estocastica, sin raiz regular
#        F2t: (1-B^12)x + tendencia
#        F3: (1-B)(1-B^12)x                   modelo de lineas aereas si sale MA(1) x SMA(1)
#      Como todo filtro F(B) tiene F(0) = 1, el jacobiano es unitario y los SBC son comparables.
#   2. Identificacion sobre F3: rejilla ARMA(p<=3, q<=2)(P<=1, Q<=1)_12, admisibilidad y SBC (Sec. 8-9)
#   3. Estimacion y salida; FAC teorica del modelo contra la empirica (Nota 4)
#   4. Diagnostico: residuos, correlograma con Q, Jarque-Bera, cuadrados, sobreajuste (Sec. 10)
#   5. Estabilidad por subperiodos (regla 9) y serie desestacionalizada (regla 8)
# Atipico aditivo (AO) en 2020-04: cierre por COVID (contexto en 02). Se modela como regresor filtrado.
# Ejecutar desde viviendas_tasa/:  Rscript 04_inicios.R
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); x <- L$viv
filtrar <- function(z, f) switch(f, F1 = diff(z), F2 = diff(z, 12), F2t = diff(z, 12), F3 = diff(diff(z, 12)))
desde <- time(x)[14]                                         # primer mes disponible con F3: 1986-02
construir <- function(z, f, ao = TRUE) {
  W <- window(filtrar(z, f), start = desde)
  xr <- switch(f, F1 = dummies_est(W), F2t = cbind(TREND = seq_along(W)), NULL)
  if (ao) { a <- matrix(as.numeric(window(filtrar(indicadora(z, 2020, 4), f), start = desde)), ncol = 1, dimnames = list(NULL, "AO2020M04"))
            xr <- if (is.null(xr)) a else cbind(as.matrix(xr), a) }
  list(W = W, xreg = xr, media = f != "F3")
}

# ---- 1. Filtros rivales -----------------------------------------------------------------------------
desc <- c(F1 = "$(1-B)x$ + dummies", F2 = "$(1-B^{12})x$", F2t = "$(1-B^{12})x$ + tendencia", F3 = "$(1-B)(1-B^{12})x$")
REJ <- lapply(setNames(names(desc), names(desc)), function(f) { b <- construir(x, f); rejilla(b$W, b$xreg, b$media) })
FIL <- do.call(rbind, lapply(names(REJ), function(f) { R <- REJ[[f]]; a <- R[R$admisible, ]; m <- if (nrow(a)) a[1, ] else R[1, ]
  data.frame(filtro = f, desc = desc[f], modelo = m$modelo, k = m$p + m$q + m$P + m$Q, admisible = m$admisible, n_adm = nrow(a), SBC = m$SBC,
             sigma = m$sigma, min_mod_AR = m$min_mod_AR, LB12 = m$LB12, LB24 = m$LB24, LB36 = m$LB36, row.names = NULL) }))
FIL$dSBC <- FIL$SBC - min(FIL$SBC)
write.csv(FIL, "salidas/tablas/04_filtros.csv", row.names = FALSE); print(FIL[, -2])
tabla_tex(data.frame(FIL$filtro, FIL$desc, FIL$modelo, num(FIL$sigma, 2), sprintf("%s / %s / %s", num(FIL$LB12, 2), num(FIL$LB24, 2), num(FIL$LB36, 2)),
                     num(FIL$SBC, 1), num(FIL$dSBC, 1)),
          "04_filtros.tex", rojo = celdas(which(FIL$filtro == "F3"), 1:7), "lllcccc", c("", "Filtro", "Mejor modelo admisible", "$\\hat\\sigma$ (\\%)", "Ljung--Box $p$ (12/24/36)", "SBC", "$\\Delta$SBC"))
gan_f <- FIL$filtro[which.min(FIL$SBC)]; stopifnot(gan_f == "F3")
guardar_valor("FiltroDSBCdos", num(FIL$dSBC[FIL$filtro == "F2"], 1)); guardar_valor("FiltroDSBCuno", num(FIL$dSBC[FIL$filtro == "F1"], 1))
guardar_valor("FiltroDosRaizAR", num(FIL$min_mod_AR[FIL$filtro == "F2"], 3)); guardar_valor("FiltroDSBCdost", num(FIL$dSBC[FIL$filtro == "F2t"], 1)); guardar_valor("FiltroUnoModelo", FIL$modelo[FIL$filtro == "F1"])

# ---- 2. Identificacion sobre F3 -----------------------------------------------------------------------
b <- construir(x, "F3"); R3 <- REJ$F3
top <- head(R3, 8)
tabla_tex(data.frame(top$modelo, top$p + top$q + top$P + top$Q, num(top$sigma, 2), num(top$LB12, 2), num(top$LB24, 2), num(top$LB36, 2),
                     ifelse(top$invertible & top$estacionario, "s\\'i", "no"), num(top$AICn, 4), num(top$SCn, 4), num(top$dSBC, 1)),
          "04_rejilla.tex", rojo = celdas(1, 1:10), "lccccccccc",
          c("Modelo", "$k$", "$\\hat\\sigma$", "LB12", "LB24", "LB36", "Estac./invert.", "AIC", "SC", "$\\Delta$SBC"))
guardar_valor("RejillaN", nrow(R3)); guardar_valor("RejillaAdm", sum(R3$admisible))
stopifnot(R3$p[1] == 0, R3$q[1] == 1, R3$P[1] == 0, R3$Q[1] == 1)

# ---- 3. Estimacion ---------------------------------------------------------------------------------------
m <- ajustar(b$W, 0, 1, 0, 1, b$xreg, FALSE); k <- 2
guardar_salida(salida_arma(m, b$W, "D(LHOUST,1,12)", "Inicios: ARIMA(0,1,1)(0,1,1)12 con AO 2020M04"), "04_modelo_inicios.png")
saveRDS(list(m = m, W = b$W, xreg = b$xreg), "salidas/tablas/04_modelo.rds")
cf <- m$coef; se <- sqrt(diag(m$var.coef))
th <- -cf["ma1"]; TH <- -cf["sma1"]                         # signo de las Notas: Theta(B) = 1 - theta B
guardar_valor("ThetaV", num(th, 3)); guardar_valor("ThetaVee", num(se["ma1"], 3))
guardar_valor("SThetaV", num(TH, 3)); guardar_valor("SThetaVee", num(se["sma1"], 3))
guardar_valor("AOV", num(cf["AO2020M04"], 1)); guardar_valor("AOVee", num(se["AO2020M04"], 1))
guardar_valor("AOVpct", num(100 * (exp(cf["AO2020M04"] / 100) - 1), 1))
guardar_valor("SigmaV", num(sqrt(m$sigma2), 2)); guardar_valor("NV", length(b$W))
guardar_valor("RaizMAV", num(1 / th, 2)); guardar_valor("RaizSMAV", num((1 / TH)^(1 / 12), 3))
# FAC teorica del modelo de lineas aereas (Nota 4: MA(1) x SMA(1)_12) contra la empirica de F3
ma <- c(cf["ma1"], rep(0, 10), cf["sma1"], cf["ma1"] * cf["sma1"])
rho_teo <- ARMAacf(ma = ma, lag.max = 24)[-1]; rho_emp <- correlograma(b$W, 24)$AC
sel <- c(1, 2, 11, 12, 13, 24)
tabla_tex(data.frame(sel, num(rho_emp[sel], 2), num(rho_teo[sel], 2), c("$-\\theta/(1+\\theta^2)$", "0", "$\\theta\\Theta/[(1+\\theta^2)(1+\\Theta^2)]$",
                                                                           "$-\\Theta/(1+\\Theta^2)$", "$\\theta\\Theta/[(1+\\theta^2)(1+\\Theta^2)]$", "0")),
          "04_fac_teorica.tex", rojo = celdas(c(1, 3, 4, 5), 2:3), "cccl", c("Rezago", "FAC emp\\'irica", "FAC te\\'orica", "F\\'ormula (Nota 4)"))

# ---- 4. Diagnostico ---------------------------------------------------------------------------------------
D <- diagnostico(m, k); e <- as.numeric(residuals(m)); es <- e / sqrt(m$sigma2)
guardar_valor("LBVdoce", num(D$LB12, 2)); guardar_valor("LBVveinticuatro", num(D$LB24, 2)); guardar_valor("LBVtreintayseis", num(D$LB36, 2))
guardar_valor("JBV", num(D$JB, 1)); guardar_valor("JBVp", pval(D$JB_p)); guardar_valor("LBcuadV", pval(D$LB2_24))
at <- which(abs(es) > 3.5)
guardar_valor("AtipicosV", if (length(at)) paste(sprintf("%s (%s)", etiqueta_mes(b$W)[at], num(es[at], 1)), collapse = "; ") else "ninguno")
guardar_salida(salida_correlograma(residuals(m), "Correlograma de los residuos: inicios", 36, k, rojo_k = c(12, 24, 36)), "04_corr_residuos.png")
# Sin el AO, abril de 2020 queda como un residuo de varias desviaciones estandar
m_sin <- ajustar(b$W, 0, 1, 0, 1, NULL, FALSE); e_sin <- residuals(m_sin) / sqrt(m_sin$sigma2)
guardar_valor("ResidAbrSinAO", num(e_sin[etiqueta_mes(b$W) == "2020-04"], 1))
guardar_valor("JBVsinAO", num(diagnostico(m_sin, 2)$JB, 0))
# Figura de residuos
tt <- as.numeric(time(b$W)); C <- correlograma(e, 36)
g1 <- ggplot(data.frame(t = tt, e = es), aes(t, e)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3, fill = gris, alpha = 0.3) +
  geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_line(colour = col_viv, linewidth = 0.35) +
  scale_x_continuous(breaks = seq(1985, 2025, 10)) + labs(title = "Residuos estandarizados", subtitle = "franja: +-3", x = NULL, y = NULL) + tema
g2 <- ggplot(C, aes(k, AC)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -C$banda[1], ymax = C$banda[1], fill = gris, alpha = 0.45) +
  geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_segment(aes(xend = k, yend = 0, colour = abs(AC) > banda), linewidth = 1.1) +
  scale_colour_manual(values = c(`TRUE` = col_viv, `FALSE` = "#9AA5B1"), guide = "none") + scale_x_continuous(breaks = c(1, 12, 24, 36)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) + labs(title = "FAC de los residuos", subtitle = sprintf("Ljung-Box p: %.2f / %.2f / %.2f (12/24/36)", D$LB12, D$LB24, D$LB36), x = "Rezago", y = NULL) + tema
g3 <- ggplot(data.frame(teo = qnorm(ppoints(length(es)))[rank(es)], mue = es), aes(teo, mue)) +
  geom_abline(slope = 1, intercept = 0, colour = tinta2, linetype = "22") + geom_point(colour = col_viv, size = 0.7) +
  labs(title = "Q-Q normal", subtitle = sprintf("Jarque-Bera = %.1f (p = %.3f)", D$JB, D$JB_p), x = "Cuantil teorico", y = NULL) + tema
guardar_fig(g1 | g2 | g3, "04_residuos.png", 12, 3.6)
# Sobreajuste
SO <- sobreajuste(b$W, c(0, 1, 0, 1), b$xreg, FALSE)
tabla_tex(data.frame(SO$agregado, num(SO$coef, 3), num(SO$t, 2), num(SO$dSBC, 1)), "04_sobreajuste.tex", rojo = celdas(seq_len(nrow(SO)), 3), "lccc",
          c("Par\\'ametro agregado", "Estimaci\\'on", "$t$", "$\\Delta$SBC"))
guardar_valor("SobreMaxT", num(max(abs(SO$t)), 2))

# ---- 5. Estabilidad y robustez ---------------------------------------------------------------------------
fila <- function(nombre, mm) { c_ <- mm$coef; s_ <- sqrt(diag(mm$var.coef))
  data.frame(nombre, sprintf("%s (%s)", num(-c_["ma1"], 2), num(s_["ma1"], 2)), sprintf("%s (%s)", num(-c_["sma1"], 2), num(s_["sma1"], 2)),
             num(sqrt(mm$sigma2), 2), length(residuals(mm))) }
submuestra <- function(a) { i <- if (a == 1) as.numeric(time(b$W)) < 2008 else as.numeric(time(b$W)) >= 2008
  W <- ts(as.numeric(b$W)[i], start = time(b$W)[which(i)[1]], frequency = 12); xr <- if (a == 2) b$xreg[i, , drop = FALSE] else NULL
  ajustar(W, 0, 1, 0, 1, xr, FALSE) }
m1 <- submuestra(1); m2 <- submuestra(2)
# Serie desestacionalizada (regla 8)
Wsa <- window(diff(L$viv_sa), start = desde)
aosa <- matrix(as.numeric(window(diff(indicadora(L$viv_sa, 2020, 4)), start = desde)), ncol = 1, dimnames = list(NULL, "AO"))
Rsa <- rejilla(Wsa, aosa, TRUE, 3, 2, 1, 1); Rsa <- Rsa[Rsa$admisible, ]       # el mejor admisible por SBC
msa <- ajustar(Wsa, Rsa$p[1], Rsa$q[1], Rsa$P[1], Rsa$Q[1], aosa, TRUE)
Csa <- correlograma(Wsa, 24); guardar_valor("RhoSAdoce", num(Csa$AC[12], 2)); guardar_valor("RhoSAveinticuatro", num(Csa$AC[24], 2))
EST <- rbind(fila("1986--2026 (modelo elegido)", m), fila("1986--2007", m1), fila("2008--2026", m2))
tabla_tex(EST, "04_estabilidad.tex", rojo = celdas(2:3, 4), "lcccc", c("Muestra", "$\\hat\\theta$ (ee)", "$\\hat\\Theta$ (ee)", "$\\hat\\sigma$", "$n$"))
guardar_valor("ThetaSubUno", num(-m1$coef["ma1"], 2)); guardar_valor("ThetaSubDos", num(-m2$coef["ma1"], 2))
guardar_valor("SThetaSubUno", num(-m1$coef["sma1"], 2)); guardar_valor("SThetaSubDos", num(-m2$coef["sma1"], 2))
guardar_valor("SigmaSubUno", num(sqrt(m1$sigma2), 2)); guardar_valor("SigmaSubDos", num(sqrt(m2$sigma2), 2))
guardar_valor("SAmodelo", Rsa$modelo[1]); guardar_valor("SAsigma", num(sqrt(msa$sigma2), 2))
guardar_valor("SAlbdoce", num(Rsa$LB12[1], 2)); guardar_valor("SAlbveinticuatro", num(Rsa$LB24[1], 2)); guardar_valor("SAlbtreintayseis", num(Rsa$LB36[1], 2))
guardar_salida(salida_arma(msa, Wsa, "D(LHOUST_SA)", sprintf("Inicios desestacionalizados: (1-B)x, %s", sub("\\\\ensuremath\\{_\\{12\\}\\}", "_12", Rsa$modelo[1])), rojo_tab = grepl("^s", names(msa$coef))), "04_modelo_sa.png")
guardar_salida(salida_correlograma(Wsa, "Correlograma de D(LHOUST_SA)", 24, rojo_k = c(12, 24)), "04_corr_sa.png")
print(EST); cat("SA:", Rsa$modelo[1], "\n")
escribir_valores("valores_04.tex")
