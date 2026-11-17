# =============================================================================
# 16_asimetria_lockin.R -- La respuesta de los inicios a la tasa es asimetrica o depende del "lock-in"?
#   Motivo: en 2022-24 la funcion de transferencia lineal predijo el doble de la caida observada (script 12).
#   Hipotesis escritas ANTES de estimar (plan, seccion 3b):
#     A. Asimetria. Para que explique la sobreprediccion, las subas tendrian que pesar MENOS que las bajas:
#        |g+| < |g-|.
#     B. Lock-in. Con la tasa de mercado muy por encima de la de los creditos ya tomados, los duenios no venden,
#        hay poca oferta de usadas y la construccion nueva resiste: el efecto de una suba se atenua, delta_1 > 0
#        (con delta_0 < 0).
#   Modelos (forma de 09 fija: b = 2, s = 0, r = 1, phi = 0,75; ruido ARMA(1,2)(1,1)12; COVID como TC):
#     A. Y_t = c + [B^2/(1 - phi B)] (delta+ X+_t + delta- X-_t) + omega TC_t + N_t, X+ = max(X, 0), X- = min(X, 0).
#        H0: delta+ = delta-, con Wald y con LR contra el modelo lineal. Sensibilidad: phi en {0,60; 0,75; 0,90}.
#        Variante por tamanio: movimientos grandes (|X| > desvio de X) contra chicos.
#     B. Y_t = c + [B^2/(1 - phi B)] (delta_0 + delta_1 gap_t) X_t + omega TC_t + N_t: el lock-in entra como
#        MODULADOR de la ganancia. No entra como regresor aditivo porque (1-B) gap = (1-B) tasa - (1-B) m y m es
#        suave: seria casi colineal con la propia tasa. gap_t centrada en su media de 1985-2021.
#        H0: delta_1 = 0 (t y LR). Variante por regimen: ganancias distintas con gap > 0 y con gap <= 0.
#   La brecha (gap): tasa vigente menos la tasa media de las hipotecas vigentes. Ese dato (FHFA, National Mortgage
#   Database) no esta en FRED y el sitio de FHFA no es accesible desde este entorno (un intento, 13-nov-2026).
#   Se usa un SUSTITUTO construido con la propia MORTGAGE30US, solo con tasas pasadas:
#     m60, m84: promedio de los 60 y 84 meses anteriores;  mexp: pesos 0,5^(s/60) sobre hasta 360 meses anteriores
#     (renormalizados con la historia disponible desde 1971).  gap = tasa - m, en puntos porcentuales.
#   Limitaciones del sustituto: no refleja las olas de refinanciacion (2020-21), asi que SOBREESTIMA la tasa media
#   de la cartera y SUBESTIMA la brecha de 2022-24; y supone que la cartera se renueva de manera uniforme.
#   Episodio 2022-24 (procedimiento de 12): cada variante se estima con datos hasta 2021-12 y se pronostica con la
#   tasa y la brecha observadas (la brecha usa solo tasas pasadas y la del mes: no hay look-ahead).
#   ADVERTENCIA: la brecha de 2022-24 queda fuera del rango observado en 1985-2021; es una extrapolacion.
#   Fuentes (ninguna esta en el apunte). Verificadas el 13-nov-2026 contra RePEc/Springer, salvo donde se indica:
#     - Asimetria por sumas parciales positivas y negativas: Granger, C. y Lee, T. (1989), "Investigation of production,
#       sales and inventory relationships using multicointegration and non-symmetric error correction models", Journal
#       of Applied Econometrics, 4(S1), S145-S159; Shin, Y., Yu, B. y Greenwood-Nimmo, M. (2014), "Modelling asymmetric
#       cointegration and dynamic multipliers in a nonlinear ARDL framework", en Festschrift in Honor of Peter Schmidt,
#       Springer, cap. 9, pp. 281-314 (DOI 10.1007/978-1-4899-8008-3_9).
#     - Lock-in hipotecario: Quigley, J. (1987), "Interest rate variations, mortgage prepayments and household
#       mobility", Review of Economics and Statistics, 69(4), 636-643; Ferreira, F., Gyourko, J. y Tracy, J. (2010),
#       "Housing busts and household mobility", Journal of Urban Economics, 68(1), 34-45.
#     - Modelos por regimen (umbral): Tong, H. (1990), Non-linear Time Series, Oxford University Press.
#   Preregistro: se corren e informan TODAS las variantes listadas; el script cuenta los modelos estimados.
#   Salidas: 16_asimetria.tex, 16_lockin.tex, 16_episodio_variantes.tex, 16_episodio_variantes.csv;
#   eviews 16_modelo_asimetria.png, 16_modelo_lockin.png; figura 16_lockin.png; valores_16.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 16_asimetria_lockin.R   (despues de 08)
# =============================================================================
source("R/funciones.R"); source("R/tf.R")
L <- cargar_base(); LC <- cargar_base(completa = TRUE); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
tf_prueba(L, I8$delta)
Y <- diff(L$viv_sa); X <- as.numeric(diff(L$tasa)); nY <- length(Y); mes <- etiqueta_mes(Y); TCm <- cbind(COVID = tf_tc(L$viv_sa, I8$delta))
fin <- which(mes == TF_ORIGEN); Xcf <- X; Xcf[(fin + 1):nY] <- 0; pre <- seq_len(nY) <= fin
ult <- L$tasa_pp[length(L$tasa_pp)]; pp <- 100 * log((ult + 1) / ult)
W <- tf_serie(Y, TF_INI:nY); n_esp <- 0
d2 <- function(v, phi = TF_PHI) tf_input(v, TF_B, 0, phi)[, 1]              # v_{t-2} / (1 - phi B)
base_xr <- function(Xv = X, phi = TF_PHI) cbind(D2 = d2(Xv, phi), TCm)
m0 <- tf_ajustar(Y, base_xr()); m0_sub <- tf_ajustar(Y, base_xr(), hasta = fin)
E0 <- tf_episodio(Y, L$viv_sa, base_xr(X), base_xr(Xcf))
# Contraste de igualdad de dos coeficientes (Wald) y LR contra el modelo lineal con el mismo phi y la misma muestra
dos_coef <- function(m, a, b, m_lin) { d <- m$coef[c(a, b)]; V <- m$var.coef[c(a, b), c(a, b)]; z <- (d[1] - d[2]) / sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
  LR <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(m_lin)))
  c(d1 = d[[1]], ee1 = sqrt(V[1, 1]), d2 = d[[2]], ee2 = sqrt(V[2, 2]), z = unname(z), p_wald = 2 * pnorm(-abs(unname(z))), LR = LR, p_lr = pchisq(max(LR, 0), 1, lower.tail = FALSE), dSBC = BIC(m) - BIC(m_lin)) }

# ---- A. Asimetria ----------------------------------------------------------------------------------------
xr_signo <- function(Xv, phi = TF_PHI) cbind(DPOS = d2(pmax(Xv, 0), phi), DNEG = d2(pmin(Xv, 0), phi), TCm)
sx <- sd(X[pre]); grande <- abs(X) > sx                                    # umbral fijado con datos hasta 2021
xr_tam <- function(Xv, phi = TF_PHI) cbind(DGRA = d2(Xv * grande, phi), DCHI = d2(Xv * !grande, phi), TCm)
asim <- function(xr, a, b, phi, hasta = nY) { n_esp <<- n_esp + 1; m <- tf_ajustar(Y, xr(X, phi), hasta = hasta); ml <- tf_ajustar(Y, base_xr(X, phi), hasta = hasta)
  c(phi = phi, hasta = hasta, dos_coef(m, a, b, ml)) }
AS <- as.data.frame(rbind(asim(xr_signo, "DPOS", "DNEG", 0.60), asim(xr_signo, "DPOS", "DNEG", 0.75), asim(xr_signo, "DPOS", "DNEG", 0.90),
                          asim(xr_signo, "DPOS", "DNEG", 0.75, fin), asim(xr_tam, "DGRA", "DCHI", 0.75), asim(xr_tam, "DGRA", "DCHI", 0.75, fin)))
AS$caso <- c("Subas / bajas, $\\phi=0{,}60$", "Subas / bajas, $\\phi=0{,}75$", "Subas / bajas, $\\phi=0{,}90$", "Subas / bajas, hasta 2021-12", "Grandes / chicos", "Grandes / chicos, hasta 2021-12")
AS$g1 <- AS$d1 / (1 - AS$phi); AS$g2 <- AS$d2 / (1 - AS$phi); AS$g1_ee <- AS$ee1 / (1 - AS$phi); AS$g2_ee <- AS$ee2 / (1 - AS$phi)
cat("\nA. Asimetria (g1 = subas o grandes; g2 = bajas o chicos)\n"); print(AS[, c("caso", "g1", "g1_ee", "g2", "g2_ee", "p_wald", "p_lr", "dSBC")], row.names = FALSE, digits = 3)
cat(sprintf("Meses con suba de tasa: %d de %d; con movimiento grande (|X| > %.2f): %d\n", sum(X > 0), nY, sx, sum(grande)))
tabla_tex(data.frame(AS$caso, sprintf("%s (%s)", num(AS$g1, 2), num(AS$g1_ee, 2)), sprintf("%s (%s)", num(AS$g2, 2), num(AS$g2_ee, 2)), num(AS$z, 2), sapply(AS$p_wald, pval), sapply(AS$p_lr, pval), num(AS$dSBC, 1)),
          "16_asimetria.tex", rojo = celdas(2, 1:7), "lcccccc",
          c("Caso", "$g$ de subas o grandes (ee cond.)", "$g$ de bajas o chicos (ee cond.)", "$z$ de la diferencia", "Wald: $p$", "LR: $p$", "$\\Delta$SBC vs lineal"))
m_as <- tf_ajustar(Y, xr_signo(X))
stopifnot(max(abs(xr_signo(X)[, "DPOS"] + xr_signo(X)[, "DNEG"] - base_xr()[, "D2"])) < 1e-9)   # con delta+ = delta- se vuelve al modelo lineal
# CUIDADO con la lectura (efecto "trinquete"): la separacion es por el signo del cambio MENSUAL. Un zigzag de la tasa
# que termina donde empezo suma subas y bajas que no se cancelan si delta+ y delta- difieren, y el promedio de ese efecto
# queda absorbido en la constante. Por eso en este modelo NO se interpreta "condicional menos tasa congelada" como efecto
# de la tasa: congelar la tasa tambien apaga el zigzag normal. Si se interpreta el pronostico condicional del nivel.
# La variante por regimen de la brecha (parte B) separa fases de suba y de baja sin este problema.
suba_ep <- sum(pmax(X[fin + 1:TF_H], 0)); baja_ep <- sum(pmin(X[fin + 1:TF_H], 0))
cat(sprintf("Episodio 2022-24: subas mensuales acumuladas %+.1f, bajas acumuladas %+.1f (neto %+.1f, en 100 ln)\n", suba_ep, baja_ep, suba_ep + baja_ep))

# ---- B. Lock-in: brecha entre la tasa vigente y la de la cartera (sustituto) --------------------------------
r <- as.numeric(LC$tasa_pp); nr <- length(r); mr <- etiqueta_mes(LC$tasa_pp)
m_mov <- function(Wn) sapply(seq_len(nr), function(t) if (t - Wn < 1) NA else mean(r[(t - Wn):(t - 1)]))
m_exp <- sapply(seq_len(nr), function(t) { s <- seq_len(min(360, t - 1)); if (length(s) < 60) return(NA); w <- 0.5^(s / 60); sum(w * r[t - s]) / sum(w) })
M <- list(m60 = m_mov(60), m84 = m_mov(84), mexp = m_exp)
al <- match(mes, mr); stopifnot(!anyNA(al), all(abs(r[al] - as.numeric(L$tasa_pp)[-1]) < 1e-9))   # mismas fechas que Y y X
GAP <- lapply(M, function(m) { g <- (r - m)[al]; stopifnot(!anyNA(g)); g })
GAPc <- lapply(GAP, function(g) g - mean(g[pre]))                          # centrada con datos hasta 2021
etq_m <- c(m60 = "Promedio de 60 meses", m84 = "Promedio de 84 meses", mexp = "Pesos exponenciales")
RG <- do.call(rbind, lapply(names(GAP), function(k) { g <- GAP[[k]]; data.frame(k = k, media_pre = mean(g[pre]), min_pre = min(g[pre]), max_pre = max(g[pre]), max_post = max(g[!pre]), mes_max = mes[!pre][which.max(g[!pre])],
                                                                            pos_pre = sum(g[pre] > 0), fuera = sum(g[E0$fut] > max(g[pre]))) }))
cat("\nBrecha (pp): rango hasta 2021-12 y maximo posterior\n"); print(RG, row.names = FALSE, digits = 3)
xr_lock <- function(Xv, k, phi = TF_PHI) cbind(D2 = d2(Xv, phi), DGAP = d2(Xv * GAPc[[k]], phi), TCm)
xr_reg <- function(Xv, k, phi = TF_PHI) cbind(DALTA = d2(Xv * (GAP[[k]] > 0), phi), DBAJA = d2(Xv * (GAP[[k]] <= 0), phi), TCm)
lock <- function(k, hasta = nY) { n_esp <<- n_esp + 1; m <- tf_ajustar(Y, xr_lock(X, k), hasta = hasta); ml <- if (hasta == nY) m0 else m0_sub
  d <- m$coef[c("D2", "DGAP")]; ee <- sqrt(diag(m$var.coef)[c("D2", "DGAP")]); LR <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(ml)))
  data.frame(k = k, hasta = hasta, d0 = d[[1]], d0_ee = ee[[1]], d1 = d[[2]], d1_ee = ee[[2]], t = d[[2]] / ee[[2]], p_t = 2 * pnorm(-abs(d[[2]] / ee[[2]])),
             LR = LR, p_lr = pchisq(max(LR, 0), 1, lower.tail = FALSE), dSBC = BIC(m) - BIC(ml), LB24 = diagnostico(m, sum(TF_RUIDO))$LB24,
             g_media = d[[1]] / (1 - TF_PHI), g_2023 = (d[[1]] + d[[2]] * mean(GAPc[[k]][substr(mes, 1, 4) == "2023"])) / (1 - TF_PHI)) }
LK <- do.call(rbind, c(lapply(names(GAP), lock), lapply(names(GAP), lock, hasta = fin)))
regimen <- function(k, hasta = nY) { n_esp <<- n_esp + 1; m <- tf_ajustar(Y, xr_reg(X, k), hasta = hasta); c(hasta = hasta, dos_coef(m, "DALTA", "DBAJA", if (hasta == nY) m0 else m0_sub)) }
RGM <- as.data.frame(rbind(regimen("m60"), regimen("m60", fin)))
cat("\nB. Lock-in continuo: (delta_0 + delta_1 gap) con cada sustituto\n"); print(LK, row.names = FALSE, digits = 3)
cat("\nB. Lock-in por regimen (m60): ganancia con gap > 0 y con gap <= 0\n"); print(round(cbind(RGM, g_alta = RGM$d1 / (1 - TF_PHI), g_baja = RGM$d2 / (1 - TF_PHI)), 3), row.names = FALSE)
muestra_txt <- function(h) ifelse(h == nY, "1985--2026", "1985--2021")
tabla_tex(rbind(data.frame(a = unname(etq_m[LK$k]), b = muestra_txt(LK$hasta), c = sprintf("%s (%s)", num(LK$d0, 3), num(LK$d0_ee, 3)), d = sprintf("%s (%s)", num(LK$d1, 3), num(LK$d1_ee, 3)),
                           e = num(LK$t, 2), f = sapply(LK$p_lr, pval), g = num(LK$dSBC, 1), h = num(LK$g_media, 2), i = num(LK$g_2023, 2)),
                data.frame(a = "R\\'egimen ($m_{60}$): gap $>0$ / gap $\\leq 0$", b = muestra_txt(RGM$hasta), c = sprintf("%s (%s)", num(RGM$d1, 3), num(RGM$ee1, 3)), d = sprintf("%s (%s)", num(RGM$d2, 3), num(RGM$ee2, 3)),
                           e = num(RGM$z, 2), f = sapply(RGM$p_lr, pval), g = num(RGM$dSBC, 1), h = num(RGM$d2 / (1 - TF_PHI), 2), i = num(RGM$d1 / (1 - TF_PHI), 2))),
          "16_lockin.tex", rojo = celdas(1, 1:9), "llccccccc",
          c("Sustituto de la tasa de la cartera", "Muestra", "$\\hat\\delta_0$ (ee)", "$\\hat\\delta_1$ (ee)", "$t$ de $\\delta_1$", "LR: $p$", "$\\Delta$SBC vs lineal", "$g$ con brecha media", "$g$ con la brecha de 2023"))
m_lk <- tf_ajustar(Y, xr_lock(X, "m60"))

# ---- C. Episodio 2022-24 con cada variante (estimacion hasta 2021-12) ----------------------------------------
sig <- function(p) !is.na(p) && p < ALFA
as_sig <- sig(AS$p_wald[2]); lk_sig <- sig(LK$p_t[LK$k == "m60" & LK$hasta == nY])
VAR_EP <- list(asim = list(etq = "Asimetria (subas / bajas)", xr = function(Xv) xr_signo(Xv), k = 1, dSBC = AS$dSBC[2]),
               lock60 = list(etq = "Lock-in, promedio de 60 meses", xr = function(Xv) xr_lock(Xv, "m60"), k = 1, dSBC = LK$dSBC[1]),
               lock84 = list(etq = "Lock-in, promedio de 84 meses", xr = function(Xv) xr_lock(Xv, "m84"), k = 1, dSBC = LK$dSBC[2]),
               lockexp = list(etq = "Lock-in, pesos exponenciales", xr = function(Xv) xr_lock(Xv, "mexp"), k = 1, dSBC = LK$dSBC[3]),
               regimen = list(etq = "Lock-in por regimen (m60)", xr = function(Xv) xr_reg(Xv, "m60"), k = 1, dSBC = RGM$dSBC[1]))
if (as_sig && lk_sig) VAR_EP$ambos <- list(etq = "Asimetria + lock-in (m60)", k = 2, dSBC = NA,
  xr = function(Xv) cbind(DPOS = d2(pmax(Xv, 0)), DNEG = d2(pmin(Xv, 0)), DGAP = d2(Xv * GAPc$m60), TCm))
cat(sprintf("\nCombinacion asimetria + lock-in: %s (asimetria p = %.3f; lock-in p = %.3f; regla: ambas < 0,05 en la muestra completa)\n",
            if (as_sig && lk_sig) "se estima" else "NO se estima", AS$p_wald[2], LK$p_t[1]))
EP <- lapply(VAR_EP, function(v) { n_esp <<- n_esp + 1; tf_episodio(Y, L$viv_sa, v$xr(X), v$xr(Xcf)) })
res_ep <- function(E, etq) { a <- tf_anual(E); e <- tf_error(E); data.frame(variante = etq, p22 = a[1], p23 = a[2], p24 = a[3], ef_tasa24 = tf_efecto_anual(E)[3], rmse = e[["rmse"]], sesgo = e[["sesgo"]]) }
RE <- rbind(res_ep(E0, "Base (script 12)"), do.call(rbind, lapply(names(EP), function(k) res_ep(EP[[k]], VAR_EP[[k]]$etq))))
RE$dSBC <- c(0, sapply(VAR_EP, `[[`, "dSBC")); obs_a <- tf_anual(E0, E0$obs)
ef_asim_ni <- RE$ef_tasa24[2]; RE$ef_tasa24[2] <- NA                      # no interpretable en el modelo por signo mensual (ver nota en la parte A)
cat("\nC. Episodio 2022-24. Observado:", round(obs_a, 1), "\n"); print(RE, row.names = FALSE, digits = 3)
tabla_tex(rbind(data.frame(v = "Observado", a = num(obs_a[1], 1), b = num(obs_a[2], 1), c = num(obs_a[3], 1), d = "--", e = "--", f = "--", g = "--"),
                data.frame(v = RE$variante, a = num(RE$p22, 1), b = num(RE$p23, 1), c = num(RE$p24, 1), d = ifelse(is.na(RE$ef_tasa24), "n.i.", num(RE$ef_tasa24, 1)), e = num(RE$rmse, 1), f = num(RE$sesgo, 1), g = ifelse(is.na(RE$dSBC), "--", num(RE$dSBC, 1)))),
          "16_episodio_variantes.tex", rojo = celdas(1:2, 2:4), "lccccccc",
          c("Modelo (datos hasta 2021-12)", "2022", "2023", "2024", "Efecto de la tasa, 2024", "RMSE 36 meses", "Sesgo medio", "$\\Delta$SBC (1985--2026)"))
write.csv(do.call(rbind, lapply(names(EP), function(k) tf_variante(6, k, VAR_EP[[k]]$etq, EP[[k]], k_extra = VAR_EP[[k]]$k, dSBC = VAR_EP[[k]]$dSBC,
                                                                 nota = if (grepl("lock|regimen|ambos", k)) "brecha de 2022-24 fuera del rango de estimacion (extrapolacion)" else ""))),
          "salidas/tablas/16_episodio_variantes.csv", row.names = FALSE)

# ---- Valores para el texto -------------------------------------------------------------------------------------
a <- AS[2, ]; guardar_valor("AsGanPos", num(a$g1, 2)); guardar_valor("AsGanPosEE", num(a$g1_ee, 2)); guardar_valor("AsGanNeg", num(a$g2, 2)); guardar_valor("AsGanNegEE", num(a$g2_ee, 2))
guardar_valor("AsWald", pval(a$p_wald)); guardar_valor("AsLR", pval(a$p_lr)); guardar_valor("AsDSBC", num(a$dSBC, 1)); guardar_valor("AsWaldSub", pval(AS$p_wald[4]))
guardar_valor("AsWaldMin", pval(min(AS$p_wald[1:3]))); guardar_valor("AsWaldMax", pval(max(AS$p_wald[1:3]))); guardar_valor("AsTamWald", pval(AS$p_wald[5]))
guardar_valor("AsGanGrande", num(AS$g1[5], 2)); guardar_valor("AsGanChico", num(AS$g2[5], 2)); guardar_valor("AsMesesSuba", sum(X > 0)); guardar_valor("AsMeses", nY)
for (k in names(GAP)) { s <- c(m60 = "Sesenta", m84 = "Ochenta", mexp = "Exp")[k]; f <- LK[LK$k == k & LK$hasta == nY, ]; fs <- LK[LK$k == k & LK$hasta == fin, ]; g <- RG[RG$k == k, ]
  guardar_valor(paste0("Lk", s, "DeltaUno"), num(f$d1, 3)); guardar_valor(paste0("Lk", s, "DeltaUnoEE"), num(f$d1_ee, 3)); guardar_valor(paste0("Lk", s, "T"), num(f$t, 2)); guardar_valor(paste0("Lk", s, "P"), pval(f$p_t))
  guardar_valor(paste0("Lk", s, "DSBC"), num(f$dSBC, 1)); guardar_valor(paste0("Lk", s, "TSub"), num(fs$t, 2)); guardar_valor(paste0("Lk", s, "PSub"), pval(fs$p_t))
  guardar_valor(paste0("Lk", s, "GanMedia"), num(f$g_media, 2)); guardar_valor(paste0("Lk", s, "GanDosTres"), num(f$g_2023, 2))
  guardar_valor(paste0("Lk", s, "MaxPre"), num(g$max_pre, 2)); guardar_valor(paste0("Lk", s, "MaxPost"), num(g$max_post, 2)); guardar_valor(paste0("Lk", s, "Fuera"), g$fuera) }
guardar_valor("LkRegGanAlta", num(RGM$d1[1] / (1 - TF_PHI), 2)); guardar_valor("LkRegGanBaja", num(RGM$d2[1] / (1 - TF_PHI), 2)); guardar_valor("LkRegWald", pval(RGM$p_wald[1])); guardar_valor("LkRegMesesAlta", RG$pos_pre[1])
guardar_valor("LkCombinada", if (as_sig && lk_sig) "s\\'i" else "no"); guardar_valor("AsLkModelos", n_esp)
for (i in seq_len(nrow(RE))) { s <- c("Base", "Asim", "LockSesenta", "LockOchenta", "LockExp", "Regimen", "Ambos")[i]
  guardar_valor(paste0("EpV", s, "Tres"), num(RE$p24[i], 1)); guardar_valor(paste0("EpV", s, "RMSE"), num(RE$rmse[i], 1)); guardar_valor(paste0("EpV", s, "Sesgo"), num(RE$sesgo[i], 1)) }
guardar_valor("EpVObsTres", num(obs_a[3], 1)); guardar_valor("AsSubasEp", num(suba_ep, 1)); guardar_valor("AsBajasEp", num(baja_ep, 1)); guardar_valor("AsEfectoNoInterp", num(ef_asim_ni, 1))

# ---- Salidas EViews y figura --------------------------------------------------------------------------------------
guardar_salida(salida_arma(m_as, W, "D(LHOUST_SA)", "Funcion de transferencia con subas (DPOS) y bajas (DNEG) de la tasa"), "16_modelo_asimetria.png")
guardar_salida(salida_arma(m_lk, W, "D(LHOUST_SA)", "Funcion de transferencia con la ganancia modulada por la brecha (DGAP)"), "16_modelo_lockin.png")
t_y <- as.numeric(time(Y))
TS <- rbind(data.frame(t = t_y, y = as.numeric(L$tasa_pp)[-1], s = "Tasa vigente"), data.frame(t = t_y, y = M$m60[al], s = "Cartera: promedio de 60 meses"),
            data.frame(t = t_y, y = M$m84[al], s = "Cartera: promedio de 84 meses"), data.frame(t = t_y, y = M$mexp[al], s = "Cartera: pesos exponenciales"))
col_s <- c("Tasa vigente" = col_tasa, "Cartera: promedio de 60 meses" = col_viv, "Cartera: promedio de 84 meses" = "#6FA3E0", "Cartera: pesos exponenciales" = "#6B4EA0")
g1 <- ggplot(TS, aes(t, y, colour = s)) + geom_line(linewidth = 0.6) + scale_colour_manual(values = col_s, name = NULL) +
  labs(title = "Tasa vigente y sustitutos de la tasa de la cartera", subtitle = "% anual; los sustitutos usan solo tasas pasadas", x = NULL, y = NULL) + tema + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
GP <- do.call(rbind, lapply(names(GAP), function(k) data.frame(t = t_y, y = GAP[[k]], s = sub("Cartera: ", "", names(col_s)[match(k, c("x", "m60", "m84", "mexp"))]))))
g2 <- ggplot(GP, aes(t, y, colour = s)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_hline(yintercept = RG$max_pre[1], colour = rojo, linetype = "22", linewidth = 0.4) +
  geom_vline(xintercept = 2022, colour = tinta2, linetype = "22", linewidth = 0.3) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = setNames(unname(col_s[2:4]), sub("Cartera: ", "", names(col_s)[2:4])), name = NULL) +
  labs(title = "Brecha: tasa vigente menos tasa de la cartera", subtitle = sprintf("puntos porcentuales; en rojo, el maximo de 1985-2021 con 60 meses (%.2f)", RG$max_pre[1]), x = NULL, y = NULL) + tema + theme(legend.position = "bottom")
anual <- function(x) 12 * exp(x / 100); t_n <- as.numeric(time(L$viv_sa)); t_f <- t_y[E0$fut]
hist <- data.frame(t = t_n, v = anual(as.numeric(L$viv_sa))); hist <- hist[hist$t >= 2019 & hist$t < 2025, ]
PR <- rbind(data.frame(t = t_f, v = anual(E0$cond), s = "Base (lineal)"), data.frame(t = t_f, v = anual(EP$asim$cond), s = "Asimetria"),
            data.frame(t = t_f, v = anual(EP$lock60$cond), s = "Lock-in (60 meses)"), data.frame(t = t_f, v = anual(EP$regimen$cond), s = "Lock-in por regimen"))
g3 <- ggplot() + geom_vline(xintercept = t_n[fin + 1], colour = tinta2, linetype = "22", linewidth = 0.3) + geom_line(data = hist, aes(t, v), colour = tinta, linewidth = 0.55) +
  geom_line(data = PR, aes(t, v, colour = s), linewidth = 0.8) +
  scale_colour_manual(values = c("Base (lineal)" = rojo, "Asimetria" = col_viv, "Lock-in (60 meses)" = col_tasa, "Lock-in por regimen" = "#6B4EA0"), name = NULL) +
  annotate("text", x = 2019.1, y = max(hist$v) * 0.99, label = "observado", hjust = 0, size = 3, colour = tinta) +
  labs(title = "Inicios: pronostico condicional desde diciembre de 2021 con cada variante", subtitle = "miles de unidades, tasa anual desestacionalizada; parametros estimados con datos hasta 2021-12", x = NULL, y = NULL) + tema + theme(legend.position = "bottom")
guardar_fig((g1 | g2) / g3 + plot_layout(heights = c(1, 1.3)), "16_lockin.png", 12, 8.6)
cat("Modelos estimados en este script:", n_esp, "\n")
escribir_valores("valores_16.tex")
