# =============================================================================
# 10_var.R -- VAR bivariado tasa hipotecaria - inicios de viviendas
#   Apunte parte 2, pp. 22-35; Enders cap. 5, secs. 4-8. Motivo: hay indicios de retroalimentacion
#   (07: los inicios anticipan a la tasa en un mes), asi que no conviene suponer exogena a ninguna: se
#   tratan simetricamente (apunte p. 22).
#   Variables (estacionarias, 03): x1 = (1-B) 100 ln(tasa), x2 = (1-B) 100 ln(inicios desestacionalizados).
#   Exogenas: constante y el COVID con la forma de 08 (cambio transitorio).
#   1. Rezagos: AIC, SBC, HQ y LR secuencial con la correccion de Sims, (T - c)(ln|S_r| - ln|S_u|) ~ chi2(4)
#      (Enders sec. 8). Se toma el LR y, si los residuos no son ruido blanco, se agregan rezagos.
#   2. Forma reducida por MCO ecuacion por ecuacion, sin recortar rezagos (Sims; apunte pp. 30-31)
#   3. Sigma y correlacion de los residuos (p. 24); estabilidad: raices de la matriz compania (p. 25)
#   4. Identificacion por Cholesky con la tasa primero (b21 = 0 del apunte p. 33: los inicios no mueven la
#      tasa dentro del mes); se recuperan b12 y las varianzas de los shocks estructurales (p. 34)
#   5. Impulso-respuesta por la representacion VMA (p. 35), con bandas bootstrap al 95%, y en niveles
#      (respuesta acumulada) para leerla como % de inicios. Sensibilidad al orden de Cholesky y al numero
#      de rezagos (el efecto entra a los 3-6 meses: con p chico no se capta).
#   6. ADELANTO (todavia no visto en clase): causalidad de Granger y descomposicion de varianza.
#   Se programa a mano (sin el paquete vars) para que cada paso se corresponda con el apunte.
# Ejecutar desde viviendas_tasa/:  Rscript 10_var.R   (despues de 08)
# =============================================================================
source("R/funciones.R")
set.seed(2026)
L <- cargar_base(); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
X <- cbind(DLMORT = diff(L$tasa), DLHOUST_SA = diff(L$viv_sa))          # orden: tasa primero
TC <- diff(ts(as.numeric(stats::filter(indicadora(L$viv_sa, 2020, 4), I8$delta, method = "recursive")), start = start(L$viv_sa), frequency = S))
Ymat <- as.matrix(X); Ex <- cbind(C = 1, COVID = as.numeric(TC)); nT <- nrow(Ymat); nk <- ncol(Ymat); nomb <- colnames(Ymat)
H <- 36; PMAX <- 13; B_BOOT <- 1000
pp <- 100 * log((L$tasa_pp[length(L$tasa_pp)] + 1) / L$tasa_pp[length(L$tasa_pp)])   # +1 pp sobre la tasa de 2026-08, en % de la tasa
fechas <- em(X[, 1])

# ---- Herramientas ------------------------------------------------------------------------------------
# Regresores de la forma reducida X_t = A0 + A1 X_{t-1} + ... + Ap X_{t-p} + e_t, desde la observacion ini+1
disenio <- function(Y, p, ini = p) { idx <- (ini + 1):nrow(Y)
  Z <- cbind(Ex[idx, , drop = FALSE], do.call(cbind, lapply(1:p, function(i) { z <- Y[idx - i, , drop = FALSE]; colnames(z) <- paste0(nomb, "_", i); z })))
  list(Y = Y[idx, , drop = FALSE], Z = Z, idx = idx) }
mco <- function(Y, p, ini = p) { d <- disenio(Y, p, ini); Bh <- solve(crossprod(d$Z), crossprod(d$Z, d$Y)); U <- d$Y - d$Z %*% Bh
  list(B = Bh, U = U, Z = d$Z, idx = d$idx, p = p, Sigma = crossprod(U) / (nrow(U) - ncol(d$Z))) }
coefA <- function(Bh, p) lapply(1:p, function(i) t(Bh[paste0(nomb, "_", i), , drop = FALSE]))   # A_i: fila = ecuacion
compania <- function(A) { p <- length(A); M <- matrix(0, nk * p, nk * p); M[1:nk, ] <- do.call(cbind, A)
  if (p > 1) M[(nk + 1):(nk * p), 1:(nk * (p - 1))] <- diag(nk * (p - 1)); M }
# Matrices de la VMA: Phi_0 = I, Phi_h = sum_i A_i Phi_{h-i}; respuesta ortogonalizada Theta_h = Phi_h P, P = chol(Sigma)
vma <- function(A, H) { p <- length(A); Ph <- list(diag(nk)); for (h in 1:H) Ph[[h + 1]] <- Reduce(`+`, lapply(1:min(h, p), function(i) A[[i]] %*% Ph[[h + 1 - i]])); Ph }
irf <- function(A, Sigma, H, orden = 1:nk) {
  P <- matrix(0, nk, nk); P[orden, orden] <- t(chol(Sigma[orden, orden]))  # Cholesky en el orden pedido
  lapply(vma(A, H), function(F) F %*% P) }                                  # [[h+1]][i, j]: respuesta de i al shock j
lb <- function(u, K, k = 0) Box.test(u, K, "Ljung-Box", fitdf = k)$p.value
portmanteau <- function(U, h, p) { T <- nrow(U); U <- scale(U, scale = FALSE); C0i <- solve(crossprod(U) / T)
  Q <- T^2 * sum(sapply(1:h, function(j) { Cj <- crossprod(U[(j + 1):T, ], U[1:(T - j), ]) / T; sum(diag(t(Cj) %*% C0i %*% Cj %*% C0i)) / (T - j) }))
  c(Q = Q, gl = nk^2 * (h - p), p = pchisq(Q, nk^2 * (h - p), lower.tail = FALSE)) }

# ---- 1. Eleccion de rezagos (misma muestra para todos los p) -------------------------------------------
SEL <- do.call(rbind, lapply(1:PMAX, function(p) { m <- mco(Ymat, p, PMAX); U <- m$U; T <- nrow(U); S <- crossprod(U) / T; k <- length(m$B)
  data.frame(p = p, T = T, c = ncol(m$Z), ld = log(det(S)), AIC = log(det(S)) + 2 * k / T, SBC = log(det(S)) + k * log(T) / T, HQ = log(det(S)) + 2 * k * log(log(T)) / T) }))
SEL$LR <- c(NA, (SEL$T[-1] - SEL$c[-1]) * (SEL$ld[-PMAX] - SEL$ld[-1])); SEL$pLR <- pchisq(SEL$LR, nk^2, lower.tail = FALSE)
p_ic <- c(AIC = SEL$p[which.min(SEL$AIC)], SBC = SEL$p[which.min(SEL$SBC)], HQ = SEL$p[which.min(SEL$HQ)])
p_lr <- { s <- which(SEL$pLR < ALFA); if (length(s)) max(s) else 1 }   # de lo general a lo particular: el p mas largo que todavia aporta
blanco <- function(p) { m <- mco(Ymat, p); all(sapply(1:nk, function(i) all(sapply(K_LB, function(K) lb(m$U[, i], K)) > ALFA))) }
p_sel <- p_lr; while (!blanco(p_sel) && p_sel < PMAX) p_sel <- p_sel + 1
cat("Rezagos -> AIC:", p_ic["AIC"], " SBC:", p_ic["SBC"], " HQ:", p_ic["HQ"], " LR:", p_lr, " elegido (residuos blancos):", p_sel, "\n")
tabla_tex(data.frame(SEL$p, num(SEL$AIC, 4), num(SEL$SBC, 4), num(SEL$HQ, 4), ifelse(is.na(SEL$LR), "--", num(SEL$LR, 2)), sapply(SEL$pLR, pval),
                     ifelse(sapply(SEL$p, function(p) p %in% c(p_ic, p_lr, p_sel)), ifelse(sapply(SEL$p, function(p) if (p %in% c(p_ic, p_lr, p_sel)) blanco(p) else NA), "s\\'i", "no"), "")),
          "10_rezagos.tex", rojo = rbind(celdas(which(SEL$p == p_ic["AIC"]), 2), celdas(which(SEL$p == p_ic["SBC"]), 3), celdas(which(SEL$p == p_ic["HQ"]), 4),
                                         celdas(which(SEL$p == p_lr), 5:6), celdas(which(SEL$p == p_sel), c(1, 7))), "ccccccc",
          c("$p$", "AIC", "SBC", "HQ", "LR (Sims)", "$p$ del LR", "Residuos blancos"))
guardar_valor("VARpAIC", p_ic["AIC"]); guardar_valor("VARpSBC", p_ic["SBC"]); guardar_valor("VARpHQ", p_ic["HQ"]); guardar_valor("VARpLR", p_lr); guardar_valor("VARp", p_sel)

# ---- 2. Forma reducida por MCO --------------------------------------------------------------------------
p <- p_sel; V <- mco(Ymat, p); A <- coefA(V$B, p); U <- V$U; Te <- nrow(U)
muestra <- paste(fechas[V$idx[1]], fechas[nT])
for (i in 1:nk) { D <- data.frame(y = Ymat[V$idx, i], V$Z[, -1]); m <- lm(y ~ ., D)
  nm <- c(COVID = "COVID_TC", setNames(sprintf("%s(-%s)", sub("_[0-9]+$", "", colnames(V$Z)[-(1:2)]), sub("^.*_", "", colnames(V$Z)[-(1:2)])), colnames(V$Z)[-(1:2)]))
  otra <- nomb[-i]; guardar_salida(salida_mco(m, nomb[i], sprintf("VAR(%d), ecuacion de %s", p, nomb[i]), muestra, nombres = nm,
                                               rojo_var = nm[grepl(paste0("^", otra, "_"), names(nm))], rojo_stats = c("Durbin-Watson stat")),
                                   sprintf("10_ecuacion_%s.png", tolower(nomb[i])), 7.4) }

# ---- 3. Sigma, correlacion de residuos y estabilidad ----------------------------------------------------
Sg <- V$Sigma; rho_e <- Sg[1, 2] / sqrt(Sg[1, 1] * Sg[2, 2]); t_rho <- rho_e * sqrt(Te - 2) / sqrt(1 - rho_e^2)
ev <- eigen(compania(A))$values; modmax <- max(Mod(ev))
guardar_valor("VARsigTasa", num(sqrt(Sg[1, 1]), 2)); guardar_valor("VARsigViv", num(sqrt(Sg[2, 2]), 2)); guardar_valor("VARcov", num(Sg[1, 2], 2))
guardar_valor("VARrho", num(rho_e, 3)); guardar_valor("VARrhoT", num(t_rho, 2)); guardar_valor("VARmodmax", num(modmax, 3)); guardar_valor("VART", Te)
guardar_valor("VARmuestra", gsub("M", "-", muestra))
# Recuperacion de la forma estructural con b21 = 0 (apunte p. 34): e_tasa = eps_tasa; e_viv = eps_viv - b12 eps_tasa
b12 <- -Sg[1, 2] / Sg[1, 1]; s_eZ <- sqrt(Sg[1, 1]); s_eY <- sqrt(Sg[2, 2] - b12^2 * Sg[1, 1])
guardar_valor("VARbdoce", num(b12, 3)); guardar_valor("VARsigEpsTasa", num(s_eZ, 2)); guardar_valor("VARsigEpsViv", num(s_eY, 2))
# Diagnostico de residuos
DG <- do.call(rbind, lapply(1:nk, function(i) { jb <- jarque_bera(U[, i]); data.frame(eq = nomb[i], LB12 = lb(U[, i], 12), LB24 = lb(U[, i], 24), LB36 = lb(U[, i], 36), JB = jb["JB"], JBp = jb["p"], row.names = NULL) }))
PM <- rbind(portmanteau(U, 12, p + 0), portmanteau(U, 24, p))
if (p >= 12) PM[1, ] <- NA                                                  # con p >= h el estadistico no tiene grados de libertad
tabla_tex(data.frame(c(sprintf("Ecuaci\\'on de %s", c("la tasa", "los inicios")), "Sistema (portmanteau multivariado)"),
                     c(num(DG$LB12, 3), if (is.na(PM[1, "p"])) "--" else pval(PM[1, "p"])), c(num(DG$LB24, 3), pval(PM[2, "p"])), c(num(DG$LB36, 3), "--"),
                     c(sprintf("%s (%s)", num(DG$JB, 1), sapply(DG$JBp, pval)), "--")),
          "10_diagnostico.tex", rojo = celdas(1:3, 3), "lcccc", c("", "LB12 $p$", "LB24 $p$", "LB36 $p$", "Jarque--Bera ($p$)"))
guardar_valor("VARpm", pval(PM[2, "p"])); guardar_valor("VARlbVivVeinticuatro", num(DG$LB24[2], 2)); guardar_valor("VARlbTasaVeinticuatro", num(DG$LB24[1], 2))
guardar_valor("VARjbTasa", num(DG$JB[1], 0)); guardar_valor("VARjbViv", num(DG$JB[2], 0))
SIG <- data.frame(c("Tasa", "Inicios"), num(Sg[, 1], 3), num(Sg[, 2], 3), c(num(1, 3), num(rho_e, 3)), c(num(rho_e, 3), num(1, 3)))
tabla_tex(SIG, "10_sigma.tex", rojo = rbind(celdas(1, 5), celdas(2, 4)), "lcccc", c("", "$\\Sigma$: tasa", "$\\Sigma$: inicios", "Correl.: tasa", "Correl.: inicios"))

# ---- 4-5. Impulso-respuesta con bandas bootstrap ---------------------------------------------------------
# Wild bootstrap (Rademacher): u*_t = w_t u_t con w_t = +-1 comun a las dos ecuaciones. Conserva la covarianza
# contemporanea y la heterocedasticidad de cada mes (la tasa tiene colas pesadas y episodios de volatilidad),
# a diferencia del bootstrap de filas iid (Runkle 1987). Se simula el VAR con las condiciones iniciales
# observadas, se re-estima y se recalcula la IRF. Bandas: percentiles 2,5 y 97,5.
resp <- function(A, Sg, orden = 1:nk) { R <- irf(A, Sg, H, orden); arr <- simplify2array(R)   # nk x nk x (H+1)
  list(dif = arr, niv = aperm(apply(arr, c(1, 2), cumsum), c(2, 3, 1))) }   # niv: respuesta acumulada = nivel (en % por ser 100 ln)
R0 <- resp(A, Sg); R0r <- resp(A, Sg, c(2, 1))
bootstrap <- function(V, p, B) { U <- V$U; Te <- nrow(U)
  replicate(B, { Yb <- Ymat; Ub <- U * sample(c(-1, 1), Te, replace = TRUE)
    for (t in seq_along(V$idx)) { s <- V$idx[t]; z <- c(Ex[s, ], unlist(lapply(1:p, function(i) Yb[s - i, ]))); Yb[s, ] <- drop(z %*% V$B) + Ub[t, ] }
    Vb <- mco(Yb, p); r <- resp(coefA(Vb$B, p), Vb$Sigma); c(as.vector(r$dif), as.vector(r$niv)) }) }
BOOT <- bootstrap(V, p, B_BOOT)
nd <- nk * nk * (H + 1)
banda_b <- function(i0) { q <- apply(BOOT[i0 + seq_len(nd), , drop = FALSE], 1, quantile, c(0.025, 0.975)); list(lo = array(q[1, ], c(nk, nk, H + 1)), hi = array(q[2, ], c(nk, nk, H + 1))) }
Bn <- banda_b(nd)                                                          # bandas de la respuesta en niveles
# Elasticidad implicita a 12 y 24 meses ante el shock de tasa: nivel de inicios / nivel de tasa
elas <- function(h) R0$niv[2, 1, h + 1] / R0$niv[1, 1, h + 1]
elas_b <- sapply(c(12, 24), function(h) { iv <- nd + (h * nk * nk) + (1 - 1) * nk + 2; it <- nd + (h * nk * nk) + 1; quantile(BOOT[iv, ] / BOOT[it, ], c(0.025, 0.975)) })
IR <- do.call(rbind, lapply(c(0, 1, 3, 6, 12, 24, 36), function(h) data.frame(h,
  tasa_tasa = R0$niv[1, 1, h + 1], viv_tasa = R0$niv[2, 1, h + 1], viv_tasa_lo = Bn$lo[2, 1, h + 1], viv_tasa_hi = Bn$hi[2, 1, h + 1],
  tasa_viv = R0$niv[1, 2, h + 1], tasa_viv_lo = Bn$lo[1, 2, h + 1], tasa_viv_hi = Bn$hi[1, 2, h + 1], viv_tasa_rev = R0r$niv[2, 1, h + 1], tasa_viv_rev = R0r$niv[1, 2, h + 1])))
sig_txt <- function(x, lo, hi) ifelse(lo > 0 | hi < 0, rojo_tex(sprintf("%s [%s; %s]", num(x, 2), num(lo, 2), num(hi, 2))), sprintf("%s [%s; %s]", num(x, 2), num(lo, 2), num(hi, 2)))
tabla_tex(data.frame(IR$h, num(IR$tasa_tasa, 2), sig_txt(IR$viv_tasa, IR$viv_tasa_lo, IR$viv_tasa_hi), num(IR$viv_tasa_rev, 2),
                     sig_txt(IR$tasa_viv, IR$tasa_viv_lo, IR$tasa_viv_hi), num(IR$tasa_viv_rev, 2)),
          "10_irf.tex", alineacion = "cccccc", encabezado = c("Meses", "Tasa ante shock de tasa", "Inicios ante shock de tasa [IC 95\\%]", "Id., orden inverso",
                                                             "Tasa ante shock de inicios [IC 95\\%]", "Id., orden inverso"))
write.csv(IR, "salidas/tablas/10_irf.csv", row.names = FALSE)
guardar_valor("VARshockTasa", num(R0$niv[1, 1, 1], 2)); guardar_valor("VARshockTasaPP", num(L$tasa_pp[length(L$tasa_pp)] * (exp(R0$niv[1, 1, 1] / 100) - 1), 2))
for (h in c(12, 24)) { s <- c(`12` = "Doce", `24` = "Veinticuatro")[as.character(h)]; r <- IR[IR$h == h, ]
  guardar_valor(paste0("VARvivTasa", s), num(r$viv_tasa, 2)); guardar_valor(paste0("VARvivTasaLo", s), num(r$viv_tasa_lo, 2)); guardar_valor(paste0("VARvivTasaHi", s), num(r$viv_tasa_hi, 2))
  guardar_valor(paste0("VARtasaViv", s), num(r$tasa_viv, 2)); guardar_valor(paste0("VARtasaVivLo", s), num(r$tasa_viv_lo, 2)); guardar_valor(paste0("VARtasaVivHi", s), num(r$tasa_viv_hi, 2))
  guardar_valor(paste0("VARelas", s), num(elas(h), 2)); guardar_valor(paste0("VARelasPP", s), num(elas(h) * pp, 1))
  eb <- elas_b[, as.character(h) == c("12", "24")]; guardar_valor(paste0("VARelasLo", s), num(eb[1], 2)); guardar_valor(paste0("VARelasHi", s), num(eb[2], 2)) }
guardar_valor("VARboot", B_BOOT)

# ---- Sensibilidad al numero de rezagos -----------------------------------------------------------------
# Respuesta acumulada a 12 meses con cada p que sugiere algun criterio, mas p = 12 (un anio). Con p chico los
# residuos de la ecuacion de inicios quedan autocorrelacionados y el efecto (que entra a los 3-6 meses) se pierde.
p_sens <- sort(unique(c(p_ic, p_lr, p_sel, 12)))
SENS <- do.call(rbind, lapply(p_sens, function(q) { Vq <- mco(Ymat, q); Rq <- resp(coefA(Vq$B, q), Vq$Sigma)
  Bq <- if (q == p) BOOT else bootstrap(Vq, q, 500); i12 <- nd + 12 * nk * nk
  data.frame(p = q, viv = Rq$niv[2, 1, 13], viv_lo = quantile(Bq[i12 + 2, ], 0.025), viv_hi = quantile(Bq[i12 + 2, ], 0.975),
             tasa = Rq$niv[1, 2, 13], tasa_lo = quantile(Bq[i12 + 3, ], 0.025), tasa_hi = quantile(Bq[i12 + 3, ], 0.975),
             elas = Rq$niv[2, 1, 13] / Rq$niv[1, 1, 13], lb24 = lb(Vq$U[, 2], 24), pm24 = portmanteau(Vq$U, 24, q)["p"], row.names = NULL) }))
crit <- sapply(SENS$p, function(q) paste(c(names(p_ic)[p_ic == q], if (q == p_lr) "LR", if (q == p_sel) "elegido", if (q == 12) "un a\\~no"), collapse = ", "))
tabla_tex(data.frame(SENS$p, crit, sig_txt(SENS$viv, SENS$viv_lo, SENS$viv_hi), num(SENS$elas * pp, 1), sig_txt(SENS$tasa, SENS$tasa_lo, SENS$tasa_hi),
                     num(SENS$lb24, 3), sapply(SENS$pm24, pval)),
          "10_sensibilidad_p.tex", rojo = celdas(which(SENS$p == p), 1:2), "clccccc",
          c("$p$", "Criterio", "Inicios ante shock de tasa, 12 meses [IC 95\\%]", "\\% inicios por +1 pp", "Tasa ante shock de inicios, 12 meses [IC 95\\%]",
            "LB24 $p$ (ec.\\ inicios)", "Portmanteau 24 $p$"))
write.csv(SENS, "salidas/tablas/10_sensibilidad_p.csv", row.names = FALSE)
for (q in SENS$p) guardar_valor(paste0("VARsensPP", c("Uno", "Dos", "Tres", "Cuatro", "Cinco", "Seis", "Siete", "Ocho", "Nueve", "Diez", "Once", "Doce", "Trece")[q]),
                                num(SENS$elas[SENS$p == q] * pp, 1))
print(SENS)

# ---- 6. ADELANTO: Granger y descomposicion de varianza ----------------------------------------------------
granger <- function(i) { j <- setdiff(1:nk, i); Zr <- V$Z[, !colnames(V$Z) %in% paste0(nomb[j], "_", 1:p)]
  u_r <- V$Y[, i] - Zr %*% solve(crossprod(Zr), crossprod(Zr, V$Y[, i])); u_u <- U[, i]
  F <- ((sum(u_r^2) - sum(u_u^2)) / p) / (sum(u_u^2) / (Te - ncol(V$Z))); c(F = F, p = pf(F, p, Te - ncol(V$Z), lower.tail = FALSE)) }
V$Y <- Ymat[V$idx, ]
GR <- rbind(granger(2), granger(1))                                         # fila 1: la tasa causa a los inicios; fila 2: al reves
tabla_tex(data.frame(c("La tasa no causa (Granger) a los inicios", "Los inicios no causan (Granger) a la tasa"), num(GR[, "F"], 2), sprintf("(%d, %d)", p, Te - ncol(V$Z)), sapply(GR[, "p"], pval)),
          "10_granger.tex", rojo = celdas(which(GR[, "p"] < ALFA), 4), "lccc", c("$H_0$", "$F$", "gl", "$p$"))
guardar_valor("VARgrTasaViv", pval(GR[1, "p"])); guardar_valor("VARgrVivTasa", pval(GR[2, "p"]))
fevd <- function(h) { Th <- R0$dif[, , 1:h, drop = FALSE]; num_ <- apply(Th^2, c(1, 2), sum); num_ / rowSums(num_) }
FV <- do.call(rbind, lapply(c(1, 6, 12, 24, 36), function(h) { f <- fevd(h); data.frame(h, tasa_por_viv = 100 * f[1, 2], viv_por_tasa = 100 * f[2, 1]) }))
tabla_tex(data.frame(FV$h, num(FV$viv_por_tasa, 1), num(FV$tasa_por_viv, 1)), "10_fevd.tex", rojo = celdas(which(FV$h == 24), 2:3), "ccc",
          c("Horizonte (meses)", "\\% de la varianza de los inicios explicado por la tasa", "\\% de la varianza de la tasa explicado por los inicios"))
guardar_valor("VARfevdViv", num(FV$viv_por_tasa[FV$h == 24], 1)); guardar_valor("VARfevdTasa", num(FV$tasa_por_viv[FV$h == 24], 1))
print(SEL[, c("p", "AIC", "SBC", "HQ", "LR", "pLR")]); print(DG); print(PM); print(IR); print(GR); print(FV); cat("rho e:", rho_e, " max|raiz|:", modmax, "\n")

# ---- Figuras ---------------------------------------------------------------------------------------------
etq_v <- c("Tasa", "Inicios")
FI <- do.call(rbind, lapply(1:nk, function(i) do.call(rbind, lapply(1:nk, function(j) data.frame(
  panel = sprintf("Respuesta de %s ante shock de %s", tolower(etq_v[i]), tolower(etq_v[j])), clave = i != j,
  h = 0:H, y = R0$niv[i, j, ], lo = Bn$lo[i, j, ], hi = Bn$hi[i, j, ], rev = R0r$niv[i, j, ])))))
FI$panel <- factor(FI$panel, levels = unique(FI$panel))
g1 <- ggplot(FI, aes(h)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = clave), alpha = 0.18) + geom_line(aes(y = rev), colour = tinta2, linetype = "22", linewidth = 0.5) +
  geom_line(aes(y = y, colour = clave), linewidth = 0.9) +
  scale_colour_manual(values = c(`TRUE` = rojo, `FALSE` = col_viv), guide = "none") + scale_fill_manual(values = c(`TRUE` = rojo, `FALSE` = col_viv), guide = "none") +
  facet_wrap(~panel, scales = "free_y") + scale_x_continuous(breaks = seq(0, H, 6)) +
  labs(title = sprintf("Impulso-respuesta del VAR(%d), en niveles (respuesta acumulada, %%)", p),
       subtitle = sprintf("Shocks de 1 desvio, Cholesky con la tasa primero; banda bootstrap 95%% (%d replicas); punteado: orden inverso", B_BOOT),
       x = "Meses", y = "%") + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
guardar_fig(g1, "10_irf.png", 11, 7)
RZ <- data.frame(re = Re(ev), im = Im(ev)); circ <- data.frame(a = seq(0, 2 * pi, length.out = 200))
g2 <- ggplot() + geom_path(data = circ, aes(cos(a), sin(a)), colour = tinta2, linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = gris) + geom_vline(xintercept = 0, colour = gris) +
  geom_point(data = RZ, aes(re, im), colour = col_viv, size = 1.8) + coord_equal() +
  labs(title = "Raices de la matriz compania", subtitle = sprintf("modulo maximo = %.3f < 1: VAR estable", modmax), x = "Real", y = "Imaginaria") + tema
guardar_fig(g2, "10_raices.png", 4.6, 4.6)
escribir_valores("valores_10.tex")
