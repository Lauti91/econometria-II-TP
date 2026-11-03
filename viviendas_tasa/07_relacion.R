# =============================================================================
# 07_relacion.R -- La relacion entre las series, recien despues de tener cada una estacionaria (regla 4)
#   1. Regresion en niveles y la alarma de regresion espuria R^2 > d_W (Granger-Newbold; Nota 6, Sec. 3.3)
#   2. Correlacion cruzada de las innovaciones de los modelos univariados (04 y 05) y test de Haugh:
#      S = n sum_{k=-M}^{M} r_k^2 ~ chi2(2M+1) bajo independencia. k > 0: la tasa anticipa a los inicios;
#      k < 0: los inicios anticipan a la tasa (retroalimentacion via la Fed). Por direccion se excluye k = 0
#      (contemporaneo, sin direccion): M grados de libertad de cada lado.
#   3. Rezagos distribuidos con el MISMO filtro (1-B)(1-B^12) de ambos lados: si y = v(B)x + n, entonces
#      F(B)y = v(B)F(B)x + F(B)n y la respuesta v(B) no cambia. Es un anticipo de la funcion de
#      transferencia (segunda parte): la inferencia por MCO es orientativa (el ruido F(B)n esta autocorrelacionado).
#   Las dos series en 100 ln: los coeficientes son elasticidades (% de inicios por 1% de tasa).
# Ejecutar desde viviendas_tasa/:  Rscript 07_relacion.R   (despues de 04 y 05)
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); x <- L$viv; y <- L$tasa
M4 <- readRDS("salidas/tablas/04_modelo.rds"); M5 <- readRDS("salidas/tablas/05_modelo.rds")

# ---- 1. Alarma de regresion espuria -------------------------------------------------------------------
alarma <- function(m, nombre) { s <- summary(m); e <- residuals(m)
  data.frame(nombre, num(coef(m)[2], 3), num(s$coefficients[2, 3], 2), num(s$r.squared, 3), num(dw(e), 3), ifelse(s$r.squared > dw(e), "S\\'i", "No")) }
m_niv <- lm(LHOUST ~ LMORT, data.frame(LHOUST = as.numeric(x), LMORT = as.numeric(y)))
m_dum <- lm(as.numeric(x) ~ as.numeric(y) + dummies_est(x))
F <- function(z) diff(diff(z, 12))
m_fil <- lm(as.numeric(F(x)) ~ as.numeric(F(y)) + as.numeric(F(indicadora(x, 2020, 4))))
AL <- rbind(alarma(m_niv, "Niveles"), alarma(m_dum, "Niveles + dummies mensuales"), alarma(m_fil, "$(1-B)(1-B^{12})$ en ambas (+ AO)"))
tabla_tex(AL, "07_alarma.tex", rojo = rbind(celdas(1, 2), celdas(2, 6), celdas(3, 4)), "lccccc", c("Regresi\\'on de $\\ln$ inicios sobre $\\ln$ tasa", "Pendiente", "$t$", "$R^2$", "$d_W$", "¿$R^2>d_W$?"))
guardar_salida(salida_mco(m_niv, "LHOUST", "Regresion en niveles: LHOUST sobre LMORT", paste(range(em(x)), collapse = " "), rojo_var = "LMORT", rojo_stats = c("R-squared", "Durbin-Watson stat")), "07_reg_niveles.png")
sn <- summary(m_niv); sd_ <- summary(m_dum)
guardar_valor("NivR", num(sn$r.squared, 3)); guardar_valor("NivDW", num(dw(residuals(m_niv)), 3)); guardar_valor("NivT", num(sn$coefficients[2, 3], 1))
guardar_valor("NivPend", num(coef(m_niv)[2], 2))
guardar_valor("DumR", num(sd_$r.squared, 3)); guardar_valor("DumDW", num(dw(residuals(m_dum)), 3))
guardar_valor("FilR", num(summary(m_fil)$r.squared, 3)); guardar_valor("FilDW", num(dw(residuals(m_fil)), 2))

# ---- 2. Correlacion cruzada de innovaciones y Haugh --------------------------------------------------------
ev <- ts(as.numeric(residuals(M4$m)), start = start(M4$W), frequency = 12); et <- ts(as.numeric(residuals(M5$m)), start = start(M5$W), frequency = 12)
Z <- ts.intersect(ev, et)
ccf_tab <- function(desde, hasta, nombre) {
  z <- window(Z, start = desde, end = hasta); n <- nrow(z)
  r <- ccf(as.numeric(z[, "ev"]), as.numeric(z[, "et"]), lag.max = 24, plot = FALSE)   # ccf(a, b) da corr(a_{t+k}, b_t): r(k) = corr(et_t, ev_{t+k})
  data.frame(periodo = nombre, k = as.numeric(r$lag), r = as.numeric(r$acf), banda = 1.96 / sqrt(n), n = n)
}
CC <- rbind(ccf_tab(start(Z), end(Z), "1986--2026"), ccf_tab(start(Z), c(2007, 12), "1986--2007"), ccf_tab(c(2008, 1), end(Z), "2008--2026"))
CC$sig <- abs(CC$r) > CC$banda
haugh <- function(d, M, lado = "ambos") { d <- d[abs(d$k) <= M & (lado == "ambos" | (lado == "pos" & d$k > 0) | (lado == "neg" & d$k < 0)), ]
  S <- d$n[1] * sum(d$r^2); c(S = S, gl = nrow(d), p = pchisq(S, nrow(d), lower.tail = FALSE)) }
HA <- do.call(rbind, lapply(split(CC, factor(CC$periodo, levels = unique(CC$periodo))), function(d) data.frame(periodo = d$periodo[1], n = d$n[1],
  M6 = haugh(d, 6)["S"], gl6 = haugh(d, 6)["gl"], p6 = haugh(d, 6)["p"], M12 = haugh(d, 12)["S"], gl12 = haugh(d, 12)["gl"], p12 = haugh(d, 12)["p"],
  pos = haugh(d, 12, "pos")["p"], neg = haugh(d, 12, "neg")["p"])))
tabla_tex(data.frame(HA$periodo, HA$n, sprintf("%s (%d)", num(HA$M6, 1), HA$gl6), sapply(HA$p6, pval), sprintf("%s (%d)", num(HA$M12, 1), HA$gl12), sapply(HA$p12, pval),
                     sapply(HA$pos, pval), sapply(HA$neg, pval)),
          "07_haugh.tex", rojo = which(cbind(NA, NA, NA, HA$p6, NA, HA$p12, HA$pos, HA$neg) < 0.05, arr.ind = TRUE), "lccccccc", c("Per\\'iodo", "$n$", "$S$, $M{=}6$ (gl)", "$p$", "$S$, $M{=}12$ (gl)", "$p$", "$p$, solo $k>0$", "$p$, solo $k<0$"))
sig <- CC[CC$sig & CC$periodo == "1986--2026", ]
guardar_valor("CCsig", paste(sprintf("$k=%d$ (%s)", sig$k, num(sig$r, 2)), collapse = ", "))
guardar_valor("HaughSeisP", pval(HA$p6[1])); guardar_valor("HaughDocePos", pval(HA$pos[1])); guardar_valor("HaughDoceNeg", pval(HA$neg[1]))
guardar_valor("BandaCC", num(CC$banda[1], 3))
print(HA); print(sig)

# ---- 3. Rezagos distribuidos con el mismo filtro ----------------------------------------------------------
K <- 12
rd <- function(desde, hasta) {
  fy <- F(x); fx <- F(y); ao <- F(indicadora(x, 2020, 4))
  X <- do.call(cbind, lapply(0:K, function(j) stats::lag(fx, -j))); A <- window(ts.intersect(fy, X, ao), start = desde, end = hasta)
  D <- data.frame(A[, 1], A[, 2:(K + 2)]); names(D) <- c("Y", sprintf("DLMORT_%d", 0:K))
  if (any(A[, K + 3] != 0)) D$AO2020M04 <- A[, K + 3]
  m <- lm(Y ~ ., D); V <- vcov(m); b <- coef(m)[sprintf("DLMORT_%d", 0:K)]
  ac <- sapply(0:K, function(j) { i <- sprintf("DLMORT_%d", 0:j); c(sum(b[1:(j + 1)]), sqrt(sum(V[i, i]))) })
  list(m = m, D = D, acum = data.frame(k = 0:K, v = b, acum = ac[1, ], ee = ac[2, ]))
}
R0 <- rd(c(1987, 2), c(2026, 8)); R1 <- rd(c(1987, 2), c(2007, 12)); R2 <- rd(c(2008, 1), c(2026, 8))
nombres <- c(setNames(c("D(LMORT,1,12)", sprintf("D(LMORT(-%d),1,12)", 1:K)), sprintf("DLMORT_%d", 0:K)), AO2020M04 = "D(AO2020M04,1,12)")
guardar_salida(salida_mco(R0$m, "D(LHOUST,1,12)", "Rezagos distribuidos 0-12, filtro (1-B)(1-B^12) en ambas",
                          sprintf("%s %s", em(window(F(x), start = c(1987, 2)))[1], em(F(x))[length(F(x))]), nombres = nombres,
                          rojo_var = unname(nombres[setdiff(names(which(summary(R0$m)$coefficients[, 4] < 0.05)), "(Intercept)")])), "07_reg_rezagos.png")
ACU <- rbind(transform(R0$acum, periodo = "1987--2026"), transform(R1$acum, periodo = "1987--2007"), transform(R2$acum, periodo = "2008--2026"))
write.csv(ACU, "salidas/tablas/07_respuesta.csv", row.names = FALSE)
pp <- 100 * log((exp(y[length(y)] / 100) + 1) / exp(y[length(y)] / 100))      # +1 pp sobre la tasa de 2026-08, en % de la tasa
sel <- c(0, 3, 6, 12)
TA <- do.call(rbind, lapply(sel, function(j) { f <- function(R) R$acum[R$acum$k == j, ]
  data.frame(j, sprintf("%s (%s)", num(f(R0)$acum, 2), num(f(R0)$ee, 2)), num(f(R0)$acum * pp, 1), sprintf("%s (%s)", num(f(R1)$acum, 2), num(f(R1)$ee, 2)),
             sprintf("%s (%s)", num(f(R2)$acum, 2), num(f(R2)$ee, 2))) }))
tabla_tex(TA, "07_respuesta.tex", rojo = celdas(4, 1:5), "ccccc", c("Meses", "Elasticidad acumulada (ee)", "\\% de inicios por +1 pp", "1987--2007", "2008--2026"))
a12 <- R0$acum[R0$acum$k == 12, ]; a3 <- R0$acum[R0$acum$k == 3, ]
guardar_valor("ElasDoce", num(a12$acum, 2)); guardar_valor("ElasDoceee", num(a12$ee, 2)); guardar_valor("ElasTres", num(a3$acum, 2)); guardar_valor("ElasTresee", num(a3$ee, 2))
guardar_valor("PPequiv", num(pp, 1)); guardar_valor("PPDoce", num(a12$acum * pp, 1)); guardar_valor("PPDoceLo", num((a12$acum - 1.96 * a12$ee) * pp, 0))
guardar_valor("PPDoceHi", num((a12$acum + 1.96 * a12$ee) * pp, 0)); guardar_valor("TasaUltima", num(exp(y[length(y)] / 100), 2))
guardar_valor("RDdw", num(dw(residuals(R0$m)), 2)); guardar_valor("RDrdos", num(summary(R0$m)$r.squared, 3))

# ---- Figura -----------------------------------------------------------------------------------------------
CCp <- CC; CCp$periodo <- factor(gsub("--", "-", CCp$periodo), levels = gsub("--", "-", unique(CC$periodo)))
g1 <- ggplot(CCp, aes(k, r)) +
  geom_rect(data = unique(CCp[, c("periodo", "banda")]), aes(xmin = -24.5, xmax = 24.5, ymin = -banda, ymax = banda), inherit.aes = FALSE, fill = gris, alpha = 0.45) +
  geom_vline(xintercept = 0, colour = tinta2, linewidth = 0.3) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.1) + scale_colour_manual(values = c(`TRUE` = rojo, `FALSE` = "#9AA5B1"), guide = "none") +
  facet_wrap(~periodo, ncol = 1) + scale_x_continuous(breaks = seq(-24, 24, 6)) +
  labs(title = "Correlacion cruzada de las innovaciones: corr(tasa_t, inicios_t+k)", subtitle = "k > 0: la tasa anticipa a los inicios; k < 0: al reves. Banda +-1,96/raiz(n)",
       x = "k (meses)", y = NULL) + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
ACp <- ACU; ACp$periodo <- factor(ACp$periodo, levels = c("1987--2026", "1987--2007", "2008--2026"))
g2 <- ggplot(ACp, aes(k, acum, colour = periodo)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_ribbon(data = ACp[ACp$periodo == "1987--2026", ], aes(ymin = acum - 1.96 * ee, ymax = acum + 1.96 * ee), fill = tinta, alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) + scale_colour_manual(values = c(rojo, "#9AA5B1", col_viv), name = NULL,
                                                                             labels = c("1987-2026", "1987-2007", "2008-2026")) +
  scale_x_continuous(breaks = 0:12) +
  labs(title = "Respuesta acumulada de los inicios a la tasa (elasticidad)", subtitle = "Suma de v_0..v_k; sombra: +-1,96 ee, 1987-2026 (MCO, orientativo)",
       x = "Meses despues del cambio de tasa", y = "% de inicios por 1% de tasa") + tema + theme(legend.position = "top")
guardar_fig((g1 | g2) + plot_layout(widths = c(1.1, 1)), "07_relacion.png", 12.5, 6.5)
escribir_valores("valores_07.tex")
