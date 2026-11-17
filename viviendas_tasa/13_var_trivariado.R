# =============================================================================
# 13_var_trivariado.R -- VAR trivariado: bono del Tesoro a 10 anios -> tasa hipotecaria -> inicios
#   Pregunta: el efecto de la tasa hipotecaria sobre los inicios (09: funcion de transferencia; 10: VAR
#   bivariado) es en realidad el del costo del dinero de largo plazo (el bono)? Agrega informacion el bono?
#   La tasa hipotecaria sigue al bono a 10 anios (se fija como bono + margen), asi que el bono es la fuente
#   mas exogena respecto del mercado de viviendas.
#   Variables (estacionarias, 03): (1-B) 100 ln(bono), (1-B) 100 ln(tasa), (1-B) 100 ln(inicios SA).
#   Exogenas: constante y el COVID con la forma de 08 (cambio transitorio).
#   Herramientas (las mismas funciones de 10_var.R, generalizadas a 3 variables):
#     - VAR en forma reducida por MCO, estabilidad, Cholesky, impulso-respuesta por la VMA, descomposicion
#       de varianza y causalidad de Granger: apunte parte 2, pp. 22-48.
#     - Rezagos: AIC, SBC (apunte); HQ: Hannan y Quinn (1979), JRSS B, 41(2), 190-195; LR secuencial con la
#       correccion de Sims: Sims (1980), Econometrica, 48(1), 1-48, como lo presenta Enders (2015), Applied
#       Econometric Time Series, 4.a ed., cap. 5. Mismo criterio que en 10: se parte del LR y se agregan
#       rezagos hasta que los residuos de cada ecuacion sean ruido blanco.
#     - Portmanteau multivariado: Hosking (1980), JASA, 75(371), 602-608.
#     - Bandas: wild bootstrap (Rademacher), Goncalves y Kilian (2004), Journal of Econometrics, 123(1), 89-120.
#     - Exogeneidad por bloques (los inicios no ayudan a predecir al bloque financiero): LR del sistema con
#       la correccion de Sims (Enders 2015, cap. 5).
#   Identificacion: Cholesky con el orden bono -> tasa -> inicios (dentro del mes, el bono no responde a la
#   tasa hipotecaria ni a los inicios; la tasa responde al bono pero no a los inicios). El shock de "tasa"
#   es entonces el movimiento de la tasa hipotecaria que NO viene del bono (el margen). Sensibilidad: orden
#   tasa -> bono -> inicios. El VAR estructural se limita a esta comparacion de ordenes.
#   Salidas: tablas 13_rezagos.tex, 13_sigma.tex, 13_irf.tex, 13_fevd.tex, 13_granger.tex, 13_sensibilidad.tex,
#   13_comparacion.tex, 13_irf.csv; figura 13_irf.png; eviews 13_ecuacion_inicios.png; valores_13.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 13_var_trivariado.R   (despues de 08, 09 y 10)
# =============================================================================
source("R/funciones.R")
set.seed(2026)
L <- cargar_base(); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
X <- cbind(DLGS10 = diff(100 * log(L$bono)), DLMORT = diff(L$tasa), DLHOUST_SA = diff(L$viv_sa))   # orden: bono, tasa, inicios
TC <- diff(ts(as.numeric(stats::filter(indicadora(L$viv_sa, 2020, 4), I8$delta, method = "recursive")), start = start(L$viv_sa), frequency = S))
Ymat <- as.matrix(X); Ex <- cbind(C = 1, COVID = as.numeric(TC)); nT <- nrow(Ymat); nk <- ncol(Ymat); nomb <- colnames(Ymat)
H <- 36; PMAX <- 13; B_BOOT <- 1000
iB <- 1; iT <- 2; iV <- 3                                                   # posiciones: bono, tasa, inicios
pp <- 100 * log((L$tasa_pp[length(L$tasa_pp)] + 1) / L$tasa_pp[length(L$tasa_pp)])   # +1 pp sobre la tasa de 2026-08, en % de la tasa
fechas <- em(X[, 1])

# ---- Herramientas (de 10_var.R) ------------------------------------------------------------------------
disenio <- function(Y, p, ini = p) { idx <- (ini + 1):nrow(Y)
  Z <- cbind(Ex[idx, , drop = FALSE], do.call(cbind, lapply(1:p, function(i) { z <- Y[idx - i, , drop = FALSE]; colnames(z) <- paste0(nomb, "_", i); z })))
  list(Y = Y[idx, , drop = FALSE], Z = Z, idx = idx) }
mco <- function(Y, p, ini = p) { d <- disenio(Y, p, ini); Bh <- solve(crossprod(d$Z), crossprod(d$Z, d$Y)); U <- d$Y - d$Z %*% Bh
  list(B = Bh, U = U, Z = d$Z, Y = d$Y, idx = d$idx, p = p, Sigma = crossprod(U) / (nrow(U) - ncol(d$Z))) }
coefA <- function(Bh, p) lapply(1:p, function(i) t(Bh[paste0(nomb, "_", i), , drop = FALSE]))   # A_i: fila = ecuacion
compania <- function(A) { p <- length(A); M <- matrix(0, nk * p, nk * p); M[1:nk, ] <- do.call(cbind, A)
  if (p > 1) M[(nk + 1):(nk * p), 1:(nk * (p - 1))] <- diag(nk * (p - 1)); M }
vma <- function(A, H) { p <- length(A); Ph <- list(diag(nk)); for (h in 1:H) Ph[[h + 1]] <- Reduce(`+`, lapply(1:min(h, p), function(i) A[[i]] %*% Ph[[h + 1 - i]])); Ph }
irf <- function(A, Sigma, H, orden = 1:nk) {
  P <- matrix(0, nk, nk); P[orden, orden] <- t(chol(Sigma[orden, orden]))  # Cholesky en el orden pedido
  lapply(vma(A, H), function(F) F %*% P) }                                  # [[h+1]][i, j]: respuesta de i al shock j
lb <- function(u, K, k = 0) Box.test(u, K, "Ljung-Box", fitdf = k)$p.value
portmanteau <- function(U, h, p) { T <- nrow(U); U <- scale(U, scale = FALSE); C0i <- solve(crossprod(U) / T)
  Q <- T^2 * sum(sapply(1:h, function(j) { Cj <- crossprod(U[(j + 1):T, ], U[1:(T - j), ]) / T; sum(diag(t(Cj) %*% C0i %*% Cj %*% C0i)) / (T - j) }))
  c(Q = Q, gl = nk^2 * (h - p), p = pchisq(Q, nk^2 * (h - p), lower.tail = FALSE)) }
resp <- function(A, Sg, orden = 1:nk) { R <- irf(A, Sg, H, orden); arr <- simplify2array(R)   # nk x nk x (H+1)
  list(dif = arr, niv = aperm(apply(arr, c(1, 2), cumsum), c(2, 3, 1))) }   # niv: respuesta acumulada = nivel (en % por ser 100 ln)
ORD <- c(1, 2, 3); ORD_ALT <- c(2, 1, 3)                                    # bono-tasa-inicios y tasa-bono-inicios

# ---- 1. Eleccion de rezagos (misma muestra para todos los p) -------------------------------------------
SEL <- do.call(rbind, lapply(1:PMAX, function(p) { m <- mco(Ymat, p, PMAX); U <- m$U; T <- nrow(U); S_ <- crossprod(U) / T; k <- length(m$B)
  data.frame(p = p, T = T, c = ncol(m$Z), ld = log(det(S_)), AIC = log(det(S_)) + 2 * k / T, SBC = log(det(S_)) + k * log(T) / T, HQ = log(det(S_)) + 2 * k * log(log(T)) / T) }))
SEL$LR <- c(NA, (SEL$T[-1] - SEL$c[-1]) * (SEL$ld[-PMAX] - SEL$ld[-1])); SEL$pLR <- pchisq(SEL$LR, nk^2, lower.tail = FALSE)
p_ic <- c(AIC = SEL$p[which.min(SEL$AIC)], SBC = SEL$p[which.min(SEL$SBC)], HQ = SEL$p[which.min(SEL$HQ)])
p_lr <- { s <- which(SEL$pLR < ALFA); if (length(s)) max(s) else 1 }
blanco <- function(p) { m <- mco(Ymat, p); all(sapply(1:nk, function(i) all(sapply(K_LB, function(K) lb(m$U[, i], K)) > ALFA))) }
p_sel <- p_lr; while (!blanco(p_sel) && p_sel < PMAX) p_sel <- p_sel + 1
p_biv <- as.integer(gsub("\\D", "", grep("VARp}", readLines("salidas/tablas/valores_10.tex"), value = TRUE)))   # p del VAR bivariado (10)
cat("Rezagos -> AIC:", p_ic["AIC"], " SBC:", p_ic["SBC"], " HQ:", p_ic["HQ"], " LR:", p_lr, " elegido (residuos blancos):", p_sel, " | bivariado:", p_biv, "\n")
marcados <- unique(c(p_ic, p_lr, p_sel))
tabla_tex(data.frame(SEL$p, num(SEL$AIC, 4), num(SEL$SBC, 4), num(SEL$HQ, 4), ifelse(is.na(SEL$LR), "--", num(SEL$LR, 2)), sapply(SEL$pLR, pval),
                     sapply(SEL$p, function(p) if (p %in% marcados) (if (blanco(p)) "s\\'i" else "no") else "")),
          "13_rezagos.tex", rojo = rbind(celdas(which(SEL$p == p_ic["AIC"]), 2), celdas(which(SEL$p == p_ic["SBC"]), 3), celdas(which(SEL$p == p_ic["HQ"]), 4),
                                         celdas(which(SEL$p == p_lr), 5:6), celdas(which(SEL$p == p_sel), c(1, 7))), "ccccccc",
          c("$p$", "AIC", "SBC", "HQ", "LR (Sims)", "$p$ del LR", "Residuos blancos"))
guardar_valor("TRIpAIC", p_ic["AIC"]); guardar_valor("TRIpSBC", p_ic["SBC"]); guardar_valor("TRIpHQ", p_ic["HQ"]); guardar_valor("TRIpLR", p_lr); guardar_valor("TRIp", p_sel)

# ---- 2. Forma reducida, Sigma, estabilidad y diagnostico ------------------------------------------------
p <- p_sel; V <- mco(Ymat, p); A <- coefA(V$B, p); U <- V$U; Te <- nrow(U); Sg <- V$Sigma
muestra <- paste(fechas[V$idx[1]], fechas[nT])
Co <- cov2cor(Sg); ev <- eigen(compania(A))$values; modmax <- max(Mod(ev))
etq_v <- c("Bono", "Tasa", "Inicios")
tabla_tex(data.frame(etq_v, num(sqrt(diag(Sg)), 2), num(Co[, 1], 3), num(Co[, 2], 3), num(Co[, 3], 3)), "13_sigma.tex", rojo = rbind(celdas(2, 3), celdas(3, 3:4)), "lcccc",
          c("", "Desv\\'io del residuo", "Correl.: bono", "Correl.: tasa", "Correl.: inicios"))
DG <- do.call(rbind, lapply(1:nk, function(i) data.frame(eq = nomb[i], LB12 = lb(U[, i], 12), LB24 = lb(U[, i], 24), LB36 = lb(U[, i], 36), JB = jarque_bera(U[, i])["JB"], row.names = NULL)))
PM <- portmanteau(U, 24, p)
guardar_valor("TRIrhoBonoTasa", num(Co[iB, iT], 2)); guardar_valor("TRIrhoBonoViv", num(Co[iB, iV], 3)); guardar_valor("TRIrhoTasaViv", num(Co[iT, iV], 3))
guardar_valor("TRImodmax", num(modmax, 3)); guardar_valor("TRIT", Te); guardar_valor("TRImuestra", gsub("M", "-", muestra)); guardar_valor("TRIpm", pval(PM["p"]))
guardar_valor("TRIlbViv", num(DG$LB24[iV], 2)); guardar_valor("TRIlbTasa", num(DG$LB24[iT], 2)); guardar_valor("TRIlbBono", num(DG$LB24[iB], 2))
guardar_valor("TRIsigBono", num(sqrt(Sg[iB, iB]), 2)); guardar_valor("TRIsigTasa", num(sqrt(Sg[iT, iT]), 2))
# Ecuacion de los inicios, con el formato de EViews
{ D <- data.frame(y = Ymat[V$idx, iV], V$Z[, -1]); m <- lm(y ~ ., D)
  nm <- c(COVID = "COVID_TC", setNames(sprintf("%s(-%s)", sub("_[0-9]+$", "", colnames(V$Z)[-(1:2)]), sub("^.*_", "", colnames(V$Z)[-(1:2)])), colnames(V$Z)[-(1:2)]))
  guardar_salida(salida_mco(m, nomb[iV], sprintf("VAR(%d) trivariado, ecuacion de %s", p, nomb[iV]), muestra, nombres = nm,
                            rojo_var = nm[grepl("^DLGS10_|^DLMORT_", names(nm))], rojo_stats = c("Durbin-Watson stat")), "13_ecuacion_inicios.png", 7.4) }

# ---- 3. Impulso-respuesta con bandas (wild bootstrap) ----------------------------------------------------
R0 <- resp(A, Sg, ORD); R0a <- resp(A, Sg, ORD_ALT)
nd <- nk * nk * (H + 1)
bootstrap <- function(V, p, B) { U <- V$U; Te <- nrow(U)
  replicate(B, { Yb <- Ymat; Ub <- U * sample(c(-1, 1), Te, replace = TRUE)
    for (t in seq_along(V$idx)) { s <- V$idx[t]; z <- c(Ex[s, ], unlist(lapply(1:p, function(i) Yb[s - i, ]))); Yb[s, ] <- drop(z %*% V$B) + Ub[t, ] }
    Vb <- mco(Yb, p); Ab <- coefA(Vb$B, p); c(as.vector(resp(Ab, Vb$Sigma, ORD)$niv), as.vector(resp(Ab, Vb$Sigma, ORD_ALT)$niv)) }) }
BOOT <- bootstrap(V, p, B_BOOT)
pos <- function(i, j, h, alt = FALSE) (if (alt) nd else 0) + h * nk * nk + (j - 1) * nk + i     # respuesta en nivel de i al shock j en h
ic <- function(i, j, h, alt = FALSE) quantile(BOOT[pos(i, j, h, alt), ], c(0.025, 0.975))
# Elasticidad implicita: respuesta de los inicios / respuesta de la tasa hipotecaria ante el mismo shock
elas <- function(R, j, h) R$niv[iV, j, h + 1] / R$niv[iT, j, h + 1]
elas_ic <- function(j, h, alt = FALSE) quantile(BOOT[pos(iV, j, h, alt), ] / BOOT[pos(iT, j, h, alt), ], c(0.025, 0.975))
hs <- c(0, 1, 3, 6, 12, 24, 36)
IR <- do.call(rbind, lapply(hs, function(h) data.frame(h, tasa_bono = R0$niv[iT, iB, h + 1],
  viv_bono = R0$niv[iV, iB, h + 1], viv_bono_lo = ic(iV, iB, h)[1], viv_bono_hi = ic(iV, iB, h)[2],
  viv_tasa = R0$niv[iV, iT, h + 1], viv_tasa_lo = ic(iV, iT, h)[1], viv_tasa_hi = ic(iV, iT, h)[2],
  viv_bono_alt = R0a$niv[iV, iB, h + 1], viv_tasa_alt = R0a$niv[iV, iT, h + 1], viv_tasa_alt_lo = ic(iV, iT, h, TRUE)[1], viv_tasa_alt_hi = ic(iV, iT, h, TRUE)[2], row.names = NULL)))
sig_txt <- function(x, lo, hi) ifelse(lo > 0 | hi < 0, rojo_tex(sprintf("%s [%s; %s]", num(x, 2), num(lo, 2), num(hi, 2))), sprintf("%s [%s; %s]", num(x, 2), num(lo, 2), num(hi, 2)))
tabla_tex(data.frame(IR$h, num(IR$tasa_bono, 2), sig_txt(IR$viv_bono, IR$viv_bono_lo, IR$viv_bono_hi), sig_txt(IR$viv_tasa, IR$viv_tasa_lo, IR$viv_tasa_hi),
                     sig_txt(IR$viv_tasa_alt, IR$viv_tasa_alt_lo, IR$viv_tasa_alt_hi), num(IR$viv_bono_alt, 2)),
          "13_irf.tex", alineacion = "cccccc",
          encabezado = c("Meses", "Tasa ante shock de bono", "Inicios ante shock de bono [IC 95\\%]", "Inicios ante shock de tasa (margen) [IC 95\\%]",
                         "Orden inverso: inicios ante shock de tasa [IC 95\\%]", "Orden inverso: inicios ante shock de bono"))
write.csv(IR, "salidas/tablas/13_irf.csv", row.names = FALSE)
for (h in c(12, 24)) { s <- c(`12` = "Doce", `24` = "Veinticuatro")[as.character(h)]; r <- IR[IR$h == h, ]
  guardar_valor(paste0("TRIvivBono", s), num(r$viv_bono, 2)); guardar_valor(paste0("TRIvivBonoLo", s), num(r$viv_bono_lo, 2)); guardar_valor(paste0("TRIvivBonoHi", s), num(r$viv_bono_hi, 2))
  guardar_valor(paste0("TRIvivTasa", s), num(r$viv_tasa, 2)); guardar_valor(paste0("TRIvivTasaLo", s), num(r$viv_tasa_lo, 2)); guardar_valor(paste0("TRIvivTasaHi", s), num(r$viv_tasa_hi, 2))
  guardar_valor(paste0("TRItasaBono", s), num(r$tasa_bono, 2))
  guardar_valor(paste0("TRIelasBono", s), num(elas(R0, iB, h), 2)); guardar_valor(paste0("TRIelasBonoPP", s), num(elas(R0, iB, h) * pp, 1))
  e <- elas_ic(iB, h); guardar_valor(paste0("TRIelasBonoLo", s), num(e[1] * pp, 1)); guardar_valor(paste0("TRIelasBonoHi", s), num(e[2] * pp, 1))
  guardar_valor(paste0("TRIelasAltPP", s), num(elas(R0a, iT, h) * pp, 1)); e <- elas_ic(iT, h, TRUE)
  guardar_valor(paste0("TRIelasAltLo", s), num(e[1] * pp, 1)); guardar_valor(paste0("TRIelasAltHi", s), num(e[2] * pp, 1)) }
guardar_valor("TRIshockBono", num(R0$niv[iB, iB, 1], 2)); guardar_valor("TRIshockMargen", num(R0$niv[iT, iT, 1], 2)); guardar_valor("TRItasaBonoCero", num(R0$niv[iT, iB, 1], 2))
guardar_valor("TRIboot", B_BOOT)

# ---- 4. Descomposicion de varianza (los dos ordenes) -------------------------------------------------------
fevd <- function(R, h) { Th <- R$dif[, , 1:h, drop = FALSE]; n_ <- apply(Th^2, c(1, 2), sum); 100 * n_ / rowSums(n_) }
FV <- do.call(rbind, lapply(c(1, 6, 12, 24, 36), function(h) { f <- fevd(R0, h); fa <- fevd(R0a, h)
  data.frame(h, viv_bono = f[iV, iB], viv_tasa = f[iV, iT], viv_propia = f[iV, iV], tasa_bono = f[iT, iB], tasa_viv = f[iT, iV], bono_viv = f[iB, iV],
             viv_bono_alt = fa[iV, iB], viv_tasa_alt = fa[iV, iT], tasa_bono_alt = fa[iT, iB]) }))
tabla_tex(data.frame(FV$h, num(FV$viv_bono, 1), num(FV$viv_tasa, 1), num(FV$viv_bono + FV$viv_tasa, 1), num(FV$viv_tasa_alt, 1), num(FV$viv_bono_alt, 1), num(FV$tasa_bono, 1), num(FV$tasa_viv, 1)),
          "13_fevd.tex", rojo = celdas(which(FV$h == 24), 2:8), "cccccccc",
          c("Horizonte", "Inicios: bono", "Inicios: tasa (margen)", "Inicios: bloque financiero", "Orden inverso: tasa", "Orden inverso: bono", "Tasa: bono", "Tasa: inicios"))
f24 <- FV[FV$h == 24, ]
guardar_valor("TRIfevdVivBono", num(f24$viv_bono, 1)); guardar_valor("TRIfevdVivTasa", num(f24$viv_tasa, 1)); guardar_valor("TRIfevdVivFin", num(f24$viv_bono + f24$viv_tasa, 1))
guardar_valor("TRIfevdTasaBono", num(f24$tasa_bono, 1)); guardar_valor("TRIfevdTasaViv", num(f24$tasa_viv, 1)); guardar_valor("TRIfevdBonoViv", num(f24$bono_viv, 1))
guardar_valor("TRIfevdVivTasaAlt", num(f24$viv_tasa_alt, 1)); guardar_valor("TRIfevdVivBonoAlt", num(f24$viv_bono_alt, 1))

# ---- 5. Causalidad de Granger, variable por variable y por bloques -------------------------------------------
# F de exclusion de los p rezagos de "causa" en la ecuacion de "efecto" (los demas rezagos quedan en la ecuacion)
granger <- function(efecto, causa) { fuera <- unlist(lapply(causa, function(j) paste0(nomb[j], "_", 1:p))); Zr <- V$Z[, !colnames(V$Z) %in% fuera]
  u_r <- V$Y[, efecto] - Zr %*% solve(crossprod(Zr), crossprod(Zr, V$Y[, efecto])); u_u <- U[, efecto]; q <- length(fuera)
  F <- ((sum(u_r^2) - sum(u_u^2)) / q) / (sum(u_u^2) / (Te - ncol(V$Z))); c(F = F, gl1 = q, gl2 = Te - ncol(V$Z), p = pf(F, q, Te - ncol(V$Z), lower.tail = FALSE)) }
# Exogeneidad por bloques: los rezagos de los inicios no entran en las ecuaciones del bono y de la tasa (LR del sistema, Sims)
bloque <- function() { Zr <- V$Z[, !colnames(V$Z) %in% paste0(nomb[iV], "_", 1:p)]; Yb <- V$Y[, c(iB, iT)]
  Ur <- Yb - Zr %*% solve(crossprod(Zr), crossprod(Zr, Yb)); Uu <- U[, c(iB, iT)]
  LR <- (Te - ncol(V$Z)) * (log(det(crossprod(Ur) / Te)) - log(det(crossprod(Uu) / Te))); c(F = LR, gl1 = 2 * p, gl2 = NA, p = pchisq(LR, 2 * p, lower.tail = FALSE)) }
GR <- rbind(granger(iV, iB), granger(iV, iT), granger(iV, c(iB, iT)), granger(iT, iB), granger(iB, iT), granger(iB, iV), granger(iT, iV), bloque())
h0 <- c("El bono no causa a los inicios (dada la tasa)", "La tasa no causa a los inicios (dado el bono)", "Ni el bono ni la tasa causan a los inicios",
        "El bono no causa a la tasa", "La tasa no causa al bono", "Los inicios no causan al bono", "Los inicios no causan a la tasa",
        "Los inicios no causan al bloque bono--tasa (LR, $\\chi^2$)")
tabla_tex(data.frame(h0, num(GR[, "F"], 2), ifelse(is.na(GR[, "gl2"]), sprintf("%d", GR[, "gl1"]), sprintf("(%d, %d)", GR[, "gl1"], GR[, "gl2"])), sapply(GR[, "p"], pval)),
          "13_granger.tex", rojo = celdas(which(GR[, "p"] < ALFA), 4), "lccc", c("$H_0$", "Estad\\'istico", "gl", "$p$"))
for (i in 1:8) guardar_valor(paste0("TRIgr", c("BonoViv", "TasaViv", "FinViv", "BonoTasa", "TasaBono", "VivBono", "VivTasa", "VivBloque")[i]), pval(GR[i, "p"]))

# ---- 6. Sensibilidad al numero de rezagos y comparacion con 09 y 10 ---------------------------------------
p_sens <- sort(unique(c(p_ic, p_lr, p_sel, p_biv, 12)))
SENS <- do.call(rbind, lapply(p_sens, function(q) { Vq <- mco(Ymat, q); Aq <- coefA(Vq$B, q); Rq <- resp(Aq, Vq$Sigma, ORD); Ra <- resp(Aq, Vq$Sigma, ORD_ALT)
  data.frame(p = q, elas_bono = elas(Rq, iB, 12) * pp, elas_alt = elas(Ra, iT, 12) * pp, viv_margen = Rq$niv[iV, iT, 13], lb24 = lb(Vq$U[, iV], 24), pm24 = portmanteau(Vq$U, 24, q)["p"], row.names = NULL) }))
crit <- sapply(SENS$p, function(q) paste(c(names(p_ic)[p_ic == q], if (q == p_lr) "LR", if (q == p_sel) "elegido", if (q == p_biv) "bivariado", if (q == 12) "un a\\~no"), collapse = ", "))
tabla_tex(data.frame(SENS$p, crit, num(SENS$elas_bono, 1), num(SENS$elas_alt, 1), num(SENS$viv_margen, 2), num(SENS$lb24, 3), sapply(SENS$pm24, pval)),
          "13_sensibilidad.tex", rojo = celdas(which(SENS$p == p), 1:7), "clccccc",
          c("$p$", "Criterio", "\\% inicios por +1 pp (shock de bono)", "Id., orden inverso (shock de tasa)", "Inicios ante shock de margen, 12 meses", "LB24 $p$ (ec.\\ inicios)", "Portmanteau 24 $p$"))
print(SENS)
# Comparacion: % de inicios por +1 pp de tasa hipotecaria segun cada modelo (los de 09 y 10 se leen de sus salidas)
valor_tex <- function(archivo, nombre) sub(sprintf("^\\\\newcommand\\{\\\\%s\\}\\{(.*)\\}$", nombre), "\\1", grep(sprintf("\\{\\\\%s\\}", nombre), readLines(file.path("salidas/tablas", archivo)), value = TRUE))
eB <- elas_ic(iB, 12) * pp; eA <- elas_ic(iT, 12, TRUE) * pp
tabla_tex(data.frame(c("Funci\\'on de transferencia (09), ganancia de largo plazo", "VAR bivariado (10), 12 meses", sprintf("VAR trivariado, shock de bono, 12 meses ($p=%d$)", p),
                       "VAR trivariado, orden inverso (shock de tasa), 12 meses"),
                     c(valor_tex("valores_09.tex", "TFTasaPP"), valor_tex("valores_10.tex", "VARelasPPDoce"), num(elas(R0, iB, 12) * pp, 1), num(elas(R0a, iT, 12) * pp, 1)),
                     c("--", "--", sprintf("[%s; %s]", num(eB[1], 1), num(eB[2], 1)), sprintf("[%s; %s]", num(eA[1], 1), num(eA[2], 1)))),
          "13_comparacion.tex", rojo = celdas(3, 1:3), "lcc", c("Modelo", "\\% de inicios por +1 pp de tasa hipotecaria", "IC 95\\,\\% (bootstrap)"))
print(SEL[, c("p", "AIC", "SBC", "HQ", "LR", "pLR")]); print(round(Co, 3)); print(DG); print(PM); print(IR); print(FV); print(cbind(h0, round(GR, 4))); cat("max|raiz|:", modmax, "\n")
cat(sprintf("%% inicios por +1 pp: shock de bono %.1f [%.1f; %.1f] | orden inverso, shock de tasa %.1f [%.1f; %.1f]\n", elas(R0, iB, 12) * pp, eB[1], eB[2], elas(R0a, iT, 12) * pp, eA[1], eA[2]))

# ---- Figura -----------------------------------------------------------------------------------------------
banda <- function(i, j, alt = FALSE) t(sapply(0:H, function(h) ic(i, j, h, alt)))
panel <- function(i, j, tit) { b <- banda(i, j); data.frame(panel = tit, h = 0:H, y = R0$niv[i, j, ], lo = b[, 1], hi = b[, 2], alt = R0a$niv[i, j, ]) }
FI <- rbind(panel(iT, iB, "Tasa hipotecaria ante shock de bono"), panel(iV, iB, "Inicios ante shock de bono"), panel(iV, iT, "Inicios ante shock de tasa (margen)"))
FI$panel <- factor(FI$panel, levels = unique(FI$panel)); FI$clave <- FI$panel != levels(FI$panel)[1]
g <- ggplot(FI, aes(h)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = clave), alpha = 0.18) + geom_line(aes(y = alt), colour = tinta2, linetype = "22", linewidth = 0.5) +
  geom_line(aes(y = y, colour = clave), linewidth = 0.9) +
  scale_colour_manual(values = c(`TRUE` = rojo, `FALSE` = col_tasa), guide = "none") + scale_fill_manual(values = c(`TRUE` = rojo, `FALSE` = col_tasa), guide = "none") +
  facet_wrap(~panel, scales = "free_y", nrow = 1) + scale_x_continuous(breaks = seq(0, H, 6)) +
  labs(title = sprintf("Impulso-respuesta del VAR(%d) trivariado, en niveles (respuesta acumulada, %%)", p),
       subtitle = sprintf("Shocks de 1 desvio, Cholesky bono - tasa - inicios; banda: wild bootstrap 95%% (%d replicas); punteado: orden tasa - bono - inicios", B_BOOT),
       x = "Meses", y = "%") + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
guardar_fig(g, "13_irf.png", 12, 4.2)
escribir_valores("valores_13.tex")
