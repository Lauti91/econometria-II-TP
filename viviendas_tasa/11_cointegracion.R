# =============================================================================
# 11_cointegracion.R -- Hay una relacion de largo plazo entre los NIVELES de inicios y tasa?
#   Pregunta: si 100 ln(inicios desestacionalizados) y 100 ln(tasa hipotecaria) cointegran, lo correcto
#   seria un modelo de correccion de errores (VEC); si no, se justifica trabajar en diferencias (scripts
#   09 y 10) y se descarta el VEC. Las dos series son I(1) (ADF de 03).
#   Tres pruebas, todas con H0: NO hay cointegracion:
#     1. Engle-Granger (apunte parte 2, pp. 49-55): MCO en niveles y ADF sin constante sobre los residuos.
#        Los residuos son estimados, asi que NO valen las tablas de Dickey-Fuller: valores criticos de
#        MacKinnon, J. (2010), "Critical Values for Cointegration Tests", Queen's Economics Department
#        Working Paper 1227, cuadro 2, N = 2, con constante (superficie de respuesta en 1/T).
#        Se repite invirtiendo el orden de las variables (apunte p. 55: el resultado puede cambiar).
#     2. Johansen (apunte pp. 56-60): traza y maximo autovalor, urca::ca.jo, constante dentro de la
#        relacion de cointegracion (ecdet = "const"). Valores criticos: Osterwald-Lenum, M. (1992),
#        Oxford Bulletin of Economics and Statistics, 54(3), 461-472 (los que trae urca).
#        K (rezagos del VAR en niveles) por AIC, HQ y SBC, y sensibilidad con K = 2, 4, 10.
#     3. Bandas de Pesaran, Shin y Smith (apunte p. 55): Pesaran, M. H., Shin, Y. y Smith, R. J. (2001),
#        "Bounds testing approaches to the analysis of level relationships", Journal of Applied
#        Econometrics, 16(3), 289-326. F conjunto de los niveles rezagados en el modelo de correccion de
#        errores condicional; cotas del cuadro CI(iii) (caso III: constante sin restringir, sin tendencia),
#        k = 1. HQ: Hannan y Quinn (1979), JRSS B, 41(2), 190-195.
#   COVID: el cierre de abril de 2020 entra como en 08 (cambio transitorio, delta de 08). En un modelo con
#   niveles rezagados, el efecto omega*delta^j sobre el nivel se captura con el pulso de abril y el TC
#   rezagado un mes; los valores criticos tabulados no contemplan estas variables (se informa la sensibilidad).
#   Nota economica: la tasa tiene una tendencia descendente de 40 anios (de ~13 % en 1985 a ~3 % en 2021)
#   que los inicios no comparten: es la razon de fondo por la que no se espera una relacion estable en niveles.
#   Salidas: tablas 11_eg.tex, 11_johansen.tex, 11_pss.tex; eviews 11_largo_plazo.png, 11_adf_residuos.png;
#   figura 11_cointegracion.png; valores_11.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 11_cointegracion.R   (despues de 08)
# =============================================================================
source("R/funciones.R"); suppressPackageStartupMessages(library(urca))
L <- cargar_base(); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
y <- L$viv_sa; x <- L$tasa; n <- length(y)
pulso <- indicadora(y, 2020, 4)
TC <- ts(as.numeric(stats::filter(pulso, I8$delta, method = "recursive")), start = start(y), frequency = S)
TC1 <- c(0, TC[-n])                                                       # TC rezagado un mes
DUM <- cbind(COVID = as.numeric(pulso), TC1 = TC1)

# ---- 1. Engle-Granger -------------------------------------------------------------------------------------
vc_mackinnon <- function(T) c(`1%` = -3.89644 - 10.9519 / T - 22.527 / T^2, `5%` = -3.33613 - 6.1101 / T - 6.823 / T^2,
                              `10%` = -3.04445 - 4.2412 / T - 2.720 / T^2)   # MacKinnon (2010), N = 2, constante
PMAX <- 24
rezagos_sbc <- function(u) {                                              # SBC del ADF sobre fechas comunes
  du <- diff(u); t0 <- (PMAX + 1):length(du)
  sbc <- sapply(0:PMAX, function(p) { Z <- data.frame(d = du[t0], u1 = u[t0]); for (j in seq_len(p)) Z[[paste0("d", j)]] <- du[t0 - j]
    BIC(lm(d ~ . - 1, data = Z)) })
  which.min(sbc) - 1
}
engle_granger <- function(D, dep, reg) {
  m <- lm(reformulate(reg, dep), data = D); u <- ts(as.numeric(residuals(m)), start = start(y), frequency = S)
  p_sbc <- rezagos_sbc(as.numeric(u)); lags <- c(0, 4, 12, p_sbc)
  A <- lapply(lags, function(p) adf(u, p, "nc", v = "RESID"))
  list(m = m, u = u, p_sbc = p_sbc, A = A,
       tab = data.frame(dep = dep, lags = lags, sbc = c(FALSE, FALSE, FALSE, TRUE), tau = sapply(A, `[[`, "tau"), n = sapply(A, `[[`, "n")))
}
D <- data.frame(LHOUST_SA = as.numeric(y), LMORT = as.numeric(x))
EG <- list(dir = engle_granger(D, "LHOUST_SA", "LMORT"), inv = engle_granger(D, "LMORT", "LHOUST_SA"))
TEG <- do.call(rbind, lapply(EG, function(e) { s <- summary(e$m); v <- t(sapply(e$tab$n, vc_mackinnon))
  cbind(e$tab, beta = coef(e$m)[2], R2 = s$r.squared, DW = dw(residuals(e$m)), vc5 = v[, "5%"], vc10 = v[, "10%"]) }))
TEG$rechaza <- TEG$tau < TEG$vc5
print(TEG, row.names = FALSE)
etq_dep <- c(LHOUST_SA = "Inicios sobre tasa", LMORT = "Tasa sobre inicios")
tabla_tex(data.frame(ifelse(duplicated(TEG$dep), "", etq_dep[TEG$dep]), ifelse(duplicated(TEG$dep), "", num(TEG$beta, 2)),
                     ifelse(duplicated(TEG$dep), "", num(TEG$R2, 2)), ifelse(duplicated(TEG$dep), "", num(TEG$DW, 2)),
                     ifelse(TEG$sbc, sprintf("%d (SBC)", TEG$lags), TEG$lags), num(TEG$tau, 2), num(TEG$vc5, 2), num(TEG$vc10, 2),
                     ifelse(TEG$rechaza, "Cointegran", "No rechaza")),
          "11_eg.tex", rojo = celdas(which(TEG$sbc), 5:9), "lcccccccc",
          c("Regresi\\'on en niveles", "Pendiente", "$R^2$", "DW", "Rezagos del ADF", "$\\tau$ de los residuos", "VC 5\\,\\%", "VC 10\\,\\%", "Decisi\\'on (5\\,\\%)"))
e <- EG$dir; i_s <- which(e$tab$sbc)
guardar_valor("EGBeta", num(coef(e$m)[2], 2)); guardar_valor("EGRdos", num(summary(e$m)$r.squared, 2)); guardar_valor("EGDW", num(dw(residuals(e$m)), 2))
guardar_valor("EGTauCero", num(e$tab$tau[1], 2)); guardar_valor("EGTauCuatro", num(e$tab$tau[2], 2)); guardar_valor("EGTauDoce", num(e$tab$tau[3], 2))
guardar_valor("EGTauSBC", num(e$tab$tau[i_s], 2)); guardar_valor("EGLagSBC", e$p_sbc)
guardar_valor("EGVCcinco", num(vc_mackinnon(e$tab$n[i_s])["5%"], 2)); guardar_valor("EGVCdiez", num(vc_mackinnon(e$tab$n[i_s])["10%"], 2))
guardar_valor("EGInvBeta", num(coef(EG$inv$m)[2], 2)); guardar_valor("EGInvTauSBC", num(EG$inv$tab$tau[which(EG$inv$tab$sbc)], 2)); guardar_valor("EGInvLagSBC", EG$inv$p_sbc)
guardar_valor("EGTauMin", num(min(TEG$tau), 2)); guardar_valor("EGTauMax", num(max(TEG$tau), 2))
guardar_valor("EGRechazos", sum(TEG$rechaza))

# Salidas EViews: regresion de largo plazo y ADF de sus residuos (con los rezagos del SBC)
muestra <- sprintf("%s %s", em(y)[1], em(y)[n])
guardar_salida(salida_mco(e$m, "LHOUST_SA", "Regresion de largo plazo (Engle-Granger, paso 1)", muestra, rojo_var = "LMORT",
                          rojo_stats = c("R-squared", "Durbin-Watson stat")), "11_largo_plazo.png")
A <- e$A[[i_s]]; vc <- vc_mackinnon(A$n); cf <- summary(A$m)$coefficients
nombres <- c(X_1 = "RESID(-1)", setNames(sprintf("D(RESID(-%d))", seq_len(A$lags)), sprintf("DX_%d", seq_len(A$lags))))
previo <- list(filas = list(c("Engle-Granger tau-statistic", fmt(A$tau), ""), c("Critical values (MacKinnon 2010, N = 2, constant):"),
                            c("      1% level", fmt(vc["1%"])), c("      5% level", fmt(vc["5%"])), c("     10% level", fmt(vc["10%"]))),
               post = "*Residuos estimados: no valen los valores p de Dickey-Fuller.")
guardar_salida(salida_mco(A$m, "D(RESID)", "Engle-Granger, paso 2: ADF de los residuos", sprintf("%s %s", A$fechas[1], A$fechas[length(A$fechas)]),
                          previo = previo, cab_previa = c("Null Hypothesis: RESID has a unit root (no cointegration)", "Exogenous: None",
                                                          sprintf("Lag Length: %d (SBC, maxlag = %d)", A$lags, PMAX)),
                          nombres = nombres, rojo_var = "RESID(-1)", rojo_previo = TRUE), "11_adf_residuos.png")

# ---- 2. Johansen ---------------------------------------------------------------------------------------------
Z <- cbind(LHOUST_SA = as.numeric(y), LMORT = as.numeric(x))
# Criterios de informacion del VAR en niveles (constante y COVID), sobre fechas comunes
KMAX <- 12
criterios <- do.call(rbind, lapply(1:KMAX, function(K) { t0 <- (KMAX + 1):n
  X <- cbind(1, DUM[t0, ], do.call(cbind, lapply(1:K, function(j) Z[t0 - j, ])))
  E <- Z[t0, ] - X %*% solve(crossprod(X), crossprod(X, Z[t0, ])); Tn <- length(t0); m <- 2 * ncol(X); ld <- log(det(crossprod(E) / Tn))
  data.frame(K = K, AIC = ld + 2 * m / Tn, HQ = ld + 2 * m * log(log(Tn)) / Tn, SBC = ld + m * log(Tn) / Tn) }))
K_ic <- sapply(c("AIC", "HQ", "SBC"), function(cn) criterios$K[which.min(criterios[[cn]])])
cat("K por criterio:", paste(names(K_ic), K_ic, collapse = " | "), "\n")
johansen <- function(K, dum = DUM) {
  dv <- if (is.null(dum)) NULL else as.matrix(dum)
  tr <- ca.jo(Z, type = "trace", ecdet = "const", K = K, spec = "transitory", dumvar = dv)
  ei <- ca.jo(Z, type = "eigen", ecdet = "const", K = K, spec = "transitory", dumvar = dv)
  f <- function(o, h) { i <- which(trimws(rownames(o@cval)) == h); c(est = o@teststat[i], vc5 = unname(o@cval[i, "5pct"])) }
  data.frame(K = K, tr0 = f(tr, "r = 0  |")[["est"]], tr0_vc = f(tr, "r = 0  |")[["vc5"]], tr1 = f(tr, "r <= 1 |")[["est"]], tr1_vc = f(tr, "r <= 1 |")[["vc5"]],
             ei0 = f(ei, "r = 0  |")[["est"]], ei0_vc = f(ei, "r = 0  |")[["vc5"]])
}
Ks <- sort(unique(pmax(2, c(K_ic, 2, 4, 10))))                           # ca.jo exige K >= 2
JO <- do.call(rbind, lapply(Ks, johansen))
JO$criterio <- sapply(JO$K, function(K) { c_ <- names(K_ic)[pmax(2, K_ic) == K]; if (length(c_)) paste(c_, collapse = ", ") else "--" })
JO$dec <- ifelse(JO$tr0 > JO$tr0_vc & JO$ei0 > JO$ei0_vc, "Rechaza", ifelse(JO$tr0 > JO$tr0_vc | JO$ei0 > JO$ei0_vc, "Solo una", "No rechaza"))
print(JO, row.names = FALSE)
tabla_tex(data.frame(JO$K, JO$criterio, num(JO$tr0, 2), num(JO$tr0_vc, 2), num(JO$ei0, 2), num(JO$ei0_vc, 2), num(JO$tr1, 2), num(JO$tr1_vc, 2), JO$dec),
          "11_johansen.tex", rojo = celdas(which(JO$criterio != "--"), 1:9), "ccccccccc",
          c("$K$", "Elegido por", "Traza, $r=0$", "VC 5\\,\\%", "$\\lambda_{\\max}$, $r=0$", "VC 5\\,\\%", "Traza, $r\\leq1$", "VC 5\\,\\%", "$H_0$: $r=0$ (5\\,\\%)"))
# Sensibilidad al tratamiento del COVID (mismos K): sin variables y solo con el pulso
JS <- do.call(rbind, lapply(Ks, function(K) data.frame(K = K, completo = JO$tr0[JO$K == K], pulso = johansen(K, DUM[, "COVID", drop = FALSE])$tr0, sin = johansen(K, NULL)$tr0)))
cat("Traza r = 0 segun tratamiento del COVID:\n"); print(JS, row.names = FALSE)
# Vector que se obtendria si se aceptara r = 1 (K del AIC): beta normalizado en inicios y velocidades de ajuste alfa.
# Sirve para ver si la relacion tendria sentido economico (inicios = constante + b * tasa, con b < 0 esperado).
vj <- ca.jo(Z, type = "trace", ecdet = "const", K = max(2, K_ic[["AIC"]]), spec = "transitory", dumvar = as.matrix(DUM))
bv <- vj@V[, 1] / vj@V[1, 1]; av <- vj@W[, 1] * vj@V[1, 1]
cat(sprintf("Vector con r = 1 (K = %d): inicios = %.1f %+.2f tasa | alfa inicios %.4f, alfa tasa %.4f\n", max(2, K_ic[["AIC"]]), -bv[3], -bv[2], av[1], av[2]))
guardar_valor("JohBetaTasa", num(-bv[2], 2)); guardar_valor("JohAlfaInicios", num(av[1], 3)); guardar_valor("JohAlfaTasa", num(av[2], 3))
for (cn in names(K_ic)) guardar_valor(paste0("JohK", cn), K_ic[cn])
k_s <- JO$K == max(2, K_ic["SBC"]); k_a <- JO$K == max(2, K_ic["AIC"])
guardar_valor("JohTrazaSBC", num(JO$tr0[k_s], 2)); guardar_valor("JohTrazaAIC", num(JO$tr0[k_a], 2)); guardar_valor("JohTrazaVC", num(JO$tr0_vc[1], 2))
guardar_valor("JohEigSBC", num(JO$ei0[k_s], 2)); guardar_valor("JohEigAIC", num(JO$ei0[k_a], 2)); guardar_valor("JohEigVC", num(JO$ei0_vc[1], 2))
guardar_valor("JohTrazaMin", num(min(JO$tr0), 2)); guardar_valor("JohTrazaMax", num(max(JO$tr0), 2))
guardar_valor("JohRechazos", sum(JO$dec == "Rechaza")); guardar_valor("JohCasos", nrow(JO))
guardar_valor("JohKRechaza", if (any(JO$dec != "No rechaza")) paste(JO$K[JO$dec != "No rechaza"], collapse = ", ") else "ninguno")
guardar_valor("JohSinCovidMin", num(min(JS$sin), 2)); guardar_valor("JohSinCovidMax", num(max(JS$sin), 2))

# ---- 3. Bandas de Pesaran, Shin y Smith ---------------------------------------------------------------------
# MCE condicional: D y_t = c + rho y_{t-1} + theta x_{t-1} + sum_{i=1..p} g_i D y_{t-i} + sum_{j=0..q} b_j D x_{t-j} + COVID + e_t
# H0 (no hay relacion en niveles): rho = theta = 0. F contra las cotas I(0) e I(1): por debajo de la inferior no se
# rechaza; por encima de la superior se rechaza; en el medio, no concluye. Cotas: PSS (2001), cuadro CI(iii), k = 1.
COTAS_F <- data.frame(nivel = c("10\\,\\%", "5\\,\\%", "1\\,\\%"), I0 = c(4.04, 4.94, 6.84), I1 = c(4.78, 5.73, 7.84))
COTAS_T <- data.frame(I0 = c(-2.57, -2.86, -3.43), I1 = c(-2.91, -3.22, -3.82))    # t de rho, cuadro CII(iii), k = 1
dy <- c(NA, diff(as.numeric(y))); dx <- c(NA, diff(as.numeric(x))); QMAX <- 12; t0 <- (QMAX + 2):n
mce <- function(p, q) {
  Zr <- data.frame(dy = dy[t0], y1 = as.numeric(y)[t0 - 1], x1 = as.numeric(x)[t0 - 1], COVID = DUM[t0, "COVID"], TC1 = DUM[t0, "TC1"])
  for (i in seq_len(p)) Zr[[sprintf("dy%d", i)]] <- dy[t0 - i]
  for (j in 0:q) Zr[[sprintf("dx%d", j)]] <- dx[t0 - j]
  m <- lm(dy ~ ., data = Zr); m0 <- lm(dy ~ . - y1 - x1, data = Zr); cf <- summary(m)$coefficients
  data.frame(p = p, q = q, SBC = BIC(m), AIC = AIC(m), F = anova(m0, m)$F[2], t = cf["y1", 3], rho = cf["y1", 1],
             largo = -cf["x1", 1] / cf["y1", 1], LB12 = Box.test(residuals(m), 12, "Ljung-Box")$p.value, LB24 = Box.test(residuals(m), 24, "Ljung-Box")$p.value)
}
GR <- do.call(rbind, lapply(0:QMAX, function(p) do.call(rbind, lapply(0:QMAX, function(q) mce(p, q)))))
PS <- rbind(cbind(crit = "SBC", GR[which.min(GR$SBC), ]), cbind(crit = "AIC", GR[which.min(GR$AIC), ]), cbind(crit = "$p=q=12$", GR[GR$p == 12 & GR$q == 12, ]))
zona <- function(F, i = 2) if (F < COTAS_F$I0[i]) "No rechaza" else if (F > COTAS_F$I1[i]) "Rechaza" else "No concluye"
PS$dec <- sapply(PS$F, zona)
print(PS, row.names = FALSE)
tabla_tex(data.frame(PS$crit, sprintf("(%d, %d)", PS$p, PS$q), num(PS$F, 2), num(COTAS_F$I0[2], 2), num(COTAS_F$I1[2], 2), num(PS$t, 2),
                     num(COTAS_T$I0[2], 2), num(COTAS_T$I1[2], 2), num(PS$LB24, 2), PS$dec),
          "11_pss.tex", rojo = celdas(1, 1:10), "lccccccccc",
          c("Rezagos por", "$(p,q)$", "$F$ de los niveles", "Cota $I(0)$ 5\\,\\%", "Cota $I(1)$ 5\\,\\%", "$t$ de $\\rho$", "Cota $I(0)$", "Cota $I(1)$", "LB24 $p$", "Decisi\\'on (5\\,\\%)"))
guardar_valor("PSSF", num(PS$F[1], 2)); guardar_valor("PSSFdoce", num(PS$F[3], 2)); guardar_valor("PSSp", PS$p[1]); guardar_valor("PSSq", PS$q[1])
guardar_valor("PSSCotaInf", num(COTAS_F$I0[2], 2)); guardar_valor("PSSCotaSup", num(COTAS_F$I1[2], 2))
guardar_valor("PSSCotaInfDiez", num(COTAS_F$I0[1], 2)); guardar_valor("PSSCotaSupDiez", num(COTAS_F$I1[1], 2))
guardar_valor("PSSt", num(PS$t[1], 2)); guardar_valor("PSSDecision", PS$dec[1]); guardar_valor("PSSFMax", num(max(PS$F), 2))

# ---- Figura: la relacion en niveles y los residuos de Engle-Granger -------------------------------------------
G <- data.frame(t = as.numeric(time(y)), y = as.numeric(y), x = as.numeric(x), u = as.numeric(e$u))
G$decada <- paste0(floor(G$t / 10) * 10, "s")
g1 <- ggplot(G, aes(x, y)) + geom_point(aes(colour = decada), size = 1, alpha = 0.75) +
  geom_abline(intercept = coef(e$m)[1], slope = coef(e$m)[2], colour = rojo, linewidth = 0.7) +
  scale_colour_manual(values = c("1980s" = "#1B3A6B", "1990s" = "#2A62B8", "2000s" = "#6FA3E0", "2010s" = "#E8A25B", "2020s" = "#B8520F"), name = NULL) +
  labs(title = "Inicios y tasa en niveles", subtitle = sprintf("recta MCO: pendiente %+.2f (signo contrario al esperado); R2 = %.2f, DW = %.2f", coef(e$m)[2], summary(e$m)$r.squared, dw(residuals(e$m))),
       x = "100 ln(tasa hipotecaria)", y = "100 ln(inicios desestacionalizados)") + tema + theme(legend.position = "bottom")
g2 <- ggplot(G, aes(t, u)) + capa_recesion(recesiones(L$rec)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_line(colour = col_viv, linewidth = 0.5) +
  labs(title = "Residuos de la regresion de largo plazo", subtitle = sprintf("no vuelven a cero con rapidez: tau = %.2f con %d rezagos (VC 5 %%: %.2f); en gris, recesiones NBER", e$tab$tau[i_s], e$p_sbc, vc_mackinnon(e$tab$n[i_s])["5%"]),
       x = NULL, y = "desvio (100 ln)") + tema
guardar_fig((g1 | g2) + plot_layout(widths = c(1, 1.25)), "11_cointegracion.png", 12, 4.6)
escribir_valores("valores_11.tex")
