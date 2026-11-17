# =============================================================================
# tf.R -- Funciones comunes de la funcion de transferencia para las extensiones (scripts 15 a 19)
#   Son las mismas de 09_transferencia.R, 12_episodio_2022.R y 14_permisos.R, sin variables globales,
#   para no repetir el codigo en cada script. Los scripts 01 a 14 NO usan este archivo (no se tocaron).
#   Apunte parte 2, pp. 8-21 (funcion de transferencia en 6 pasos) y pp. 1-7 (intervencion).
#   Pronostico con input conocido y varianza del error con los pesos psi: Box, Jenkins, Reinsel y Ljung (2015),
#   Time Series Analysis: Forecasting and Control, 5.a ed., Wiley, caps. 11 y 12.
#   tf_prueba() es la prueba de regresion: tiene que reproducir la ganancia de 09 (-0,45; ee 0,15) y el
#   efecto de la tasa en 2024 de 12 (-31,9 %). Cada script la corre antes de usar estas funciones.
# Se carga con source("R/tf.R") despues de R/funciones.R.
# =============================================================================
TF_B <- 2; TF_PHI <- 0.75; TF_RUIDO <- c(1, 2, 1, 1); TF_INI <- 7; TF_H <- 36; TF_ORIGEN <- "2021-12"   # forma elegida en 09 y episodio de 12

# ---- Regresores -----------------------------------------------------------------------------------------
# Cambio transitorio del COVID (abril de 2020) en la serie diferenciada, con el delta dado
tf_tc <- function(nivel, delta) as.numeric(diff(ts(as.numeric(stats::filter(indicadora(nivel, 2020, 4), delta, method = "recursive")), start = start(nivel), frequency = S)))
tf_rezago <- function(v, j) c(rep(0, j), v[seq_len(length(v) - j)])       # v_{t-j}, con ceros al inicio
tf_filtro <- function(v, phi) as.numeric(if (phi != 0) stats::filter(as.numeric(v), phi, method = "recursive") else as.numeric(v))   # v / (1 - phi B)
# Columnas X_{t-b}, ..., X_{t-b-s}, cada una filtrada por 1/(1 - phi B)
tf_input <- function(X, b, s, phi = 0, pref = "D") { xf <- tf_filtro(X, phi)
  matrix(sapply(b:(b + s), function(j) tf_rezago(xf, j)), length(xf), dimnames = list(NULL, sprintf("%s%d", pref, b:(b + s)))) }
tf_serie <- function(Y, i) ts(as.numeric(Y)[i], start = time(Y)[i[1]], frequency = S)

# ---- Estimacion y pronostico ----------------------------------------------------------------------------
# Maxima verosimilitud exacta, como ajustar(). Si el optimizador falla desde el arranque por defecto (pasa con
# algunos regresores adicionales), se reintenta arrancando desde la estimacion por minimos cuadrados condicionales
# ("CSS-ML"): el criterio que se maximiza es el mismo, solo cambia el punto de partida.
tf_ajustar <- function(Y, XR, ruido = TF_RUIDO, ini = TF_INI, hasta = length(Y)) {
  W <- tf_serie(Y, ini:hasta); xr <- XR[ini:hasta, , drop = FALSE]
  m <- ajustar(W, ruido[1], ruido[2], ruido[3], ruido[4], xr, TRUE)
  if (is.null(m)) m <- tryCatch(suppressWarnings(arima(W, order = c(ruido[1], 0, ruido[2]), seasonal = list(order = c(ruido[3], 0, ruido[4]), period = S), xreg = xr,
                                                       include.mean = TRUE, method = "CSS-ML", optim.control = list(maxit = 1000))), error = function(e) NULL)
  m
}
tf_sbc <- function(m, extra = 0) -2 * as.numeric(logLik(m)) + (length(m$coef) + extra) * log(length(residuals(m)))
tf_predecir <- function(m, XRfut) { xreg <- xr <- XRfut                     # predict.Arima busca el objeto xreg de la llamada original en este entorno (solo mira sus columnas)
  as.numeric(predict(m, n.ahead = nrow(XRfut), newxreg = XRfut)$pred) }
# ee del pronostico del NIVEL a h pasos: sigma * raiz(suma de los psi acumulados al cuadrado). Condicional a inputs y parametros.
tf_ee_nivel <- function(m, H) sqrt(m$sigma2 * cumsum(cumsum(c(1, ARMAtoMA(ar = m$model$phi, ma = m$model$theta, lag.max = H - 1)))^2))
pct <- function(d) 100 * (exp(d / 100) - 1)                               # de puntos logaritmicos (100 ln) a %

# Episodio: estima con datos hasta 'origen' y pronostica H meses con los regresores observados (XR) y, si se da,
# con los del contrafactico (XR_cf). Y = diff(nivel). Devuelve los caminos del nivel (100 ln) y lo observado.
tf_episodio <- function(Y, nivel, XR, XR_cf = NULL, ruido = TF_RUIDO, ini = TF_INI, origen = TF_ORIGEN, H = TF_H) {
  mes <- etiqueta_mes(Y); fin <- which(mes == origen); fut <- fin + 1:H; mn <- etiqueta_mes(nivel)
  stopifnot(length(fin) == 1, mn[fin + 1] == origen, !anyNA(XR[ini:(fin + H), ]))
  m <- tf_ajustar(Y, XR, ruido, ini, fin); n0 <- as.numeric(nivel)[fin + 1]; camino <- function(R) n0 + cumsum(tf_predecir(m, R[fut, , drop = FALSE]))
  anio0 <- substr(origen, 1, 4)
  list(m = m, fin = fin, fut = fut, mes = mes[fut], anio = as.integer(substr(mes[fut], 1, 4)), obs = as.numeric(nivel)[fin + 1 + 1:H],
       cond = camino(XR), contra = if (is.null(XR_cf)) NULL else camino(XR_cf), ee = tf_ee_nivel(m, H),
       base = mean(as.numeric(nivel)[substr(mn, 1, 4) == anio0]))          # promedio del anio del origen (2021)
}
# Promedios anuales respecto del promedio de 2021 (%), error de pronostico y efecto de la tasa (condicional - contrafactico)
tf_anual <- function(E, v = E$cond) sapply(sort(unique(E$anio)), function(a) pct(mean(v[E$anio == a]) - E$base))
tf_efecto_anual <- function(E) sapply(sort(unique(E$anio)), function(a) pct(mean((E$cond - E$contra)[E$anio == a])))
tf_error <- function(E, v = E$cond) c(rmse = sqrt(mean((E$obs - v)^2)), sesgo = mean(E$obs - v))
# Fila(s) del archivo de variantes del episodio, con formato comun para la sintesis (script 19)
tf_variante <- function(etapa, clave, etiqueta, E, v = E$cond, k_extra = 0, dSBC = NA, nota = "")
  data.frame(etapa = etapa, clave = clave, etiqueta = etiqueta, mes = E$mes, h = seq_along(v), observado = E$obs, prediccion = v, base21 = E$base,
             k_extra = k_extra, dSBC = dSBC, nota = nota, stringsAsFactors = FALSE)

# ---- Preblanqueo y FCC (pasos 1 a 3 del apunte, p. 19) ------------------------------------------------------
# r(k) = corr(alfa_t, beta_{t+k}): k > 0, el input anticipa al output. Pmax = Qmax = 0 para inputs sin estacionalidad.
tf_preblanqueo <- function(X, Yl, K = 24, Pmax = 0, Qmax = 0) {
  R <- rbind(transform(rejilla(X, NULL, TRUE, 3, 2, Pmax, Qmax), media = TRUE), transform(rejilla(X, NULL, FALSE, 3, 2, Pmax, Qmax), media = FALSE))
  Ra <- R[R$admisible, ]; if (!nrow(Ra)) Ra <- R; Ra <- Ra[order(Ra$SBC), ]; g <- Ra[1, ]
  mX <- ajustar(X, g$p, g$q, g$P, g$Q, NULL, g$media); banda <- 2 / sqrt(length(X))
  filtro <- function(z) { z <- z - mean(z); ar <- mX$model$phi; ma <- mX$model$theta
    if (length(ar)) z <- stats::filter(z, c(1, -ar), sides = 1)
    if (length(ma)) z <- stats::filter(ifelse(is.na(z), 0, z), -ma, method = "recursive")
    as.numeric(z) }
  a <- filtro(as.numeric(X)); b <- filtro(as.numeric(Yl)); quemar <- max(12, length(mX$model$phi) + 1)
  ok <- !is.na(a) & !is.na(b) & seq_along(a) > quemar; stopifnot(all(diff(which(ok)) == 1))
  a <- ts(a[ok], start = time(X)[which(ok)[1]], frequency = S); b <- b[ok]
  r <- ccf(b, as.numeric(a), lag.max = K, plot = FALSE)                    # ccf(b, a) de R: corr(b_{t+k}, a_t)
  C <- data.frame(k = as.numeric(r$lag), r = as.numeric(r$acf)); C$sig <- abs(C$r) > banda; C$v <- ifelse(C$k >= 0, sd(b) / sd(a) * C$r, NA)
  stopifnot(abs(cor(as.numeric(a)[-length(a)], b[-1]) - C$r[C$k == 1]) < 0.02)   # control de direccion, a mano, para k = +1
  Qk <- function(ks) { d <- C[C$k %in% ks, ]; Q <- length(a) * (length(a) + 2) * sum(d$r^2 / (length(a) - abs(d$k))); c(Q = Q, gl = nrow(d), p = pchisq(Q, nrow(d), lower.tail = FALSE)) }
  list(modelo = g$modelo, mX = mX, alfa = a, C = C, banda = banda, Qpos = Qk(0:12), Qneg = Qk(-(1:12)), Qpos1 = Qk(1:12))
}
# Correlacion cruzada de los residuos de un modelo con el input preblanqueado (paso 6): p del Q para k >= 0 y k < 0
tf_cc_residuos <- function(m, W, alfa, gl_menos = 0, K = 24) {
  Z <- na.omit(ts.intersect(e = ts(as.numeric(residuals(m)), start = start(W), frequency = S), a = alfa))
  ce <- ccf(as.numeric(Z[, "e"]), as.numeric(Z[, "a"]), lag.max = K, plot = FALSE); CE <- data.frame(k = as.numeric(ce$lag), r = as.numeric(ce$acf)); n <- nrow(Z)
  Q <- function(ks) { d <- CE[CE$k %in% ks, ]; n * (n + 2) * sum(d$r^2 / (n - abs(d$k))) }
  c(pos = pchisq(Q(0:12), 13 - gl_menos, lower.tail = FALSE), neg = pchisq(Q(-(1:12)), 12, lower.tail = FALSE))
}

# ---- Ganancia con la incertidumbre de phi (como en 09) ----------------------------------------------------
# phi se elige por perfil, asi que var.coef solo da Var(delta | phi). Var(phi) = -1 / l''(phi) (curvatura del perfil);
# la pendiente de delta_hat(phi) da Cov(delta, phi); despues, metodo delta sobre g = sum(delta) / (1 - phi).
tf_ee_ganancia <- function(ajuste_phi, phi, nd, h = 0.02) {
  M <- lapply(c(phi - h, phi, phi + h), ajuste_phi)
  if (any(sapply(M, is.null))) return(c(ee = NA, ee_cond = NA, ee_phi = NA))
  ll <- sapply(M, function(m) as.numeric(logLik(m))); d <- sapply(M, function(m) sum(m$coef[nd]))
  curv <- (ll[3] - 2 * ll[2] + ll[1]) / h^2; v_cond <- sum(M[[2]]$var.coef[nd, nd])
  if (curv >= 0) return(c(ee = NA, ee_cond = sqrt(v_cond) / (1 - phi), ee_phi = NA))
  v_phi <- -1 / curv; pend <- (d[3] - d[1]) / (2 * h)
  V <- matrix(c(v_cond + pend^2 * v_phi, pend * v_phi, pend * v_phi, v_phi), 2); gr <- c(1 / (1 - phi), d[2] / (1 - phi)^2)
  c(ee = sqrt(drop(t(gr) %*% V %*% gr)), ee_cond = sqrt(v_cond) / (1 - phi), ee_phi = sqrt(v_phi))
}

# ---- Seleccion de (b, s, r) y estimacion completa de una forma (pasos 4 a 6; como en 09 y 14) --------------
tf_seleccion <- function(Y, X, TCm, ruido, ini = TF_INI, grilla_phi = seq(-0.8, 0.8, by = 0.2)) { out <- list()
  for (b in 0:3) for (s in 0:3) for (r in 0:1) { phis <- if (r == 0) 0 else grilla_phi
    ms <- lapply(phis, function(ph) tf_ajustar(Y, cbind(tf_input(X, b, s, ph), TCm), ruido, ini)); ll <- sapply(ms, function(m) if (is.null(m)) -Inf else as.numeric(logLik(m)))
    i <- which.max(ll); m <- ms[[i]]; if (is.null(m)) next
    out[[length(out) + 1]] <- data.frame(b = b, s = s, r = r, phi = phis[i], SBC = tf_sbc(m, r), logLik = ll[i], row.names = NULL) }
  G <- do.call(rbind, out); G$dSBC <- G$SBC - min(G$SBC); G[order(G$SBC), ] }
tf_forma <- function(Y, X, TCm, b, s, r, phi0, ruido_ini, alfa = NULL, ruido_fijo = NULL, ini = TF_INI) {
  xr <- function(ph) cbind(tf_input(X, b, s, ph), TCm); phi <- 0
  if (r == 1) { fina <- seq(max(-0.95, phi0 - 0.15), min(0.95, phi0 + 0.15), by = 0.05)
    ll <- sapply(fina, function(ph) { m <- tf_ajustar(Y, xr(ph), ruido_ini, ini); if (is.null(m)) -Inf else as.numeric(logLik(m)) }); phi <- fina[which.max(ll)] }
  W <- tf_serie(Y, ini:length(Y))
  if (is.null(ruido_fijo)) { R <- rejilla(W, xr(phi)[ini:length(Y), , drop = FALSE], TRUE, 3, 2, 1, 1); Ra <- R[R$admisible, ]; if (!nrow(Ra)) Ra <- R; ru <- unlist(Ra[1, c("p", "q", "P", "Q")]); et <- Ra$modelo[1] } else {
    ru <- ruido_fijo; et <- etiqueta_modelo(ru[1], ru[2], ru[3], ru[4]) }
  m <- tf_ajustar(Y, xr(phi), ru, ini)
  nd <- sprintf("D%d", b:(b + s)); dl <- m$coef[nd]; gan <- sum(dl) / (1 - phi); gan_ee <- sqrt(sum(m$var.coef[nd, nd])) / (1 - phi); gan_cond <- gan_ee; phi_ee <- NA
  if (r == 1) { E <- tf_ee_ganancia(function(ph) tf_ajustar(Y, xr(ph), ru, ini), phi, nd); if (!is.na(E["ee"])) { gan_ee <- unname(E["ee"]); phi_ee <- unname(E["ee_phi"]) } }
  vj <- sapply(0:240, function(j) { if (j < b) return(0); sum(sapply(b:min(j, b + s), function(i) dl[sprintf("D%d", i)] * (if (r == 1) phi^(j - i) else as.numeric(j == i)))) })
  acum <- cumsum(vj) / sum(vj)
  list(b = b, s = s, r = r, phi = phi, phi_ee = phi_ee, ruido = et, ru = ru, m = m, W = W, D = diagnostico(m, sum(ru)), vj = unname(vj), gan = gan, gan_ee = gan_ee, gan_cond = gan_cond,
       medio = sum((0:240) * vj) / sum(vj), mitad = which(acum >= 0.5)[1] - 1, p_cc = if (is.null(alfa)) c(pos = NA, neg = NA) else tf_cc_residuos(m, W, alfa, s + 1 + r), n_par = length(m$coef) + r)
}

# ---- Prueba de regresion ------------------------------------------------------------------------------------
# Reproduce los dos numeros de referencia: ganancia de 09 (-0,45; ee 0,15) y efecto de la tasa en 2024 de 12 (-31,9 %).
tf_prueba <- function(L, delta) {
  Y <- diff(L$viv_sa); X <- diff(L$tasa); TC <- cbind(COVID = tf_tc(L$viv_sa, delta))
  f <- tf_forma(Y, X, TC, TF_B, 0, 1, TF_PHI, TF_RUIDO, ruido_fijo = TF_RUIDO)
  Xcf <- X; Xcf[(which(etiqueta_mes(Y) == TF_ORIGEN) + 1):length(X)] <- 0
  E <- tf_episodio(Y, L$viv_sa, cbind(tf_input(X, TF_B, 0, TF_PHI), TC), cbind(tf_input(Xcf, TF_B, 0, TF_PHI), TC)); ef <- tf_efecto_anual(E)[3]
  cat(sprintf("Prueba de regresion de R/tf.R: ganancia %.3f (ee %.3f) [09: -0,45; 0,15] | efecto de la tasa en 2024 %.1f %% [12: -31,9]\n", f$gan, f$gan_ee, ef))
  stopifnot(f$phi == TF_PHI, abs(f$gan + 0.449) < 0.005, abs(f$gan_ee - 0.15) < 0.01, abs(ef + 31.9) < 0.1)
  invisible(list(forma = f, episodio = E))
}
