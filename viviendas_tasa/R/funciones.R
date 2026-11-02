# =============================================================================
# funciones.R -- Herramientas comunes del par "inicios de viviendas y tasa hipotecaria" (EE.UU.)
# Se carga con source("R/funciones.R") desde la carpeta viviendas_tasa/.
#
# Contenido
#   1. Parametros del trabajo (muestra, colores, estilo de graficos)
#   2. Datos
#   3. Herramientas de series: dummies, Box-Cox, ADF, correlograma con Q de Ljung-Box, ajuste ARMA
#   4. Salidas con formato EViews (ARMA, MCO, test ADF, correlograma)
#   5. Exportar a LaTeX (valores sueltos y tablas)
#
# Convencion de signos: las Notas escriben Theta(B) = 1 - theta B. R (stats::arima) y EViews usan
# 1 + theta B. Las salidas EViews se dejan con el signo de R/EViews; en el informe se aclara.
# =============================================================================

for (pk in c("ggplot2", "patchwork", "ragg", "urca")) {
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ---- 1. Parametros -------------------------------------------------------------------------------
INI <- c(1985, 1)          # muestra principal: desde 1985 (se deja afuera el regimen Volcker, ver 02)
FIN <- c(2026, 8)          # ultimo dato de HOUSTNSA en la descarga del 19-oct-2026
S   <- 12                  # periodo estacional
K_LB <- c(12, 24, 36)      # rezagos del Ljung-Box
ALFA <- 0.05

col_viv  <- "#2A62B8"      # inicios: azul
col_tasa <- "#D9701F"      # tasa: naranja
col_bono <- "#9AA5B1"
rojo <- "#C8102E"                                                   # resaltado de lo importante
tinta <- "#1F2933"; tinta2 <- "#52606D"; gris <- "#C5CCD3"; fondo <- "#FFFFFF"

tema <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = fondo, colour = NA),
        panel.background = element_rect(fill = fondo, colour = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#E4E7EB", linewidth = 0.3),
        axis.text = element_text(colour = tinta2, size = 9), axis.title = element_text(colour = tinta2, size = 9),
        plot.title = element_text(colour = tinta, face = "bold", size = 11, margin = margin(b = 2)),
        plot.subtitle = element_text(colour = tinta2, size = 9, margin = margin(b = 6)),
        plot.margin = margin(8, 12, 6, 8))
guardar_fig <- function(g, archivo, ancho, alto)
  ggsave(file.path("salidas/figuras", archivo), g, width = ancho, height = alto, dpi = 200, device = ragg::agg_png, bg = "white")

# ---- 2. Datos ------------------------------------------------------------------------------------
# base_mensual.csv la arma 01_datos.R. Series:
#   viv  = 100 ln(HOUSTNSA)   inicios de viviendas sin desestacionalizar (miles de unidades)
#   tasa = 100 ln(MORTGAGE30US promedio mensual)
#   viv_sa = 100 ln(HOUST / 12) desestacionalizada (solo para comparar, regla 8)
cargar_base <- function(completa = FALSE) {
  B <- read.csv("datos/base_mensual.csv", stringsAsFactors = FALSE)
  ini <- c(as.integer(substr(B$mes[1], 1, 4)), as.integer(substr(B$mes[1], 6, 7)))
  mk <- function(v) ts(v, start = ini, frequency = S)
  L <- list(inicios = mk(B$inicios_nsa), inicios_sa = mk(B$inicios_sa_anual / 12), tasa_pp = mk(B$tasa_hip_30),
            bono = mk(B$bono10), rec = mk(B$recesion))
  L$viv <- 100 * log(L$inicios); L$viv_sa <- 100 * log(L$inicios_sa); L$tasa <- 100 * log(L$tasa_pp)
  if (!completa) L <- lapply(L, window, start = INI, end = FIN)
  L
}
etiqueta_mes <- function(x) { t <- as.numeric(time(x)) + 1e-6; sprintf("%d-%02d", floor(t), floor((t - floor(t)) * 12) + 1) }
em <- function(x) sub("-", "M", etiqueta_mes(x))                       # 1986M02, como EViews
fecha_num <- function(a, m) a + (m - 1) / 12
indicadora <- function(x, a, m) ts(as.numeric(abs(time(x) - fecha_num(a, m)) < 1e-6), start = start(x), frequency = frequency(x))

# Recesiones NBER para sombrear
recesiones <- function(rec) {
  r <- as.numeric(rec); t <- as.numeric(time(rec)); d <- diff(c(0, r, 0))
  data.frame(ini = t[which(d == 1)], fin = t[which(d == -1) - 1] + 1 / 12)
}
capa_recesion <- function(R) annotate("rect", xmin = R$ini, xmax = R$fin, ymin = -Inf, ymax = Inf, fill = "#E4E7EB", alpha = 0.7)

# ---- 3. Herramientas de series -------------------------------------------------------------------
# Dummies estacionales centradas (mes j menos diciembre): suman cero en el anio
dummies_est <- function(x) {
  q <- cycle(x); D <- sapply(1:11, function(j) (q == j) - (q == 12)); colnames(D) <- month.abb[1:11]; D
}

# Relacion rango-mediana (Nota 6, Sec. 5.4): ln R_i = a + b ln M_i por anio; lambda = 1 - b
rango_mediana <- function(y) {
  a <- floor(as.numeric(time(y)) + 1e-6)
  M <- tapply(as.numeric(y), a, median); R <- tapply(as.numeric(y), a, function(v) diff(range(v)))
  n <- tapply(as.numeric(y), a, length); M <- M[n == 12]; R <- R[n == 12]
  m <- lm(log(R) ~ log(M)); b <- coef(m)[2]; ee <- sqrt(vcov(m)[2, 2])
  list(datos = data.frame(anio = as.integer(names(M)), mediana = as.numeric(M), rango = as.numeric(R)),
       b = unname(b), ee = ee, lambda = unname(1 - b), ic = sort(unname(1 - b + c(1, -1) * 1.96 * ee)))
}
# Transformacion de Box-Cox escalada por 100 (lambda = 0: 100 ln y)
box_cox <- function(y, l) if (abs(l) < 1e-9) 100 * log(y) else 100 * (y^l - 1) / l

# Regresion de Dickey-Fuller aumentada (Nota 6, Sec. 4.3-4.4)
#   Dx_t = [a] + [b t] + g x_{t-1} + sum_{j=1}^p c_j Dx_{t-j} + [dummies] + e_t ;  tau = t(g)
# det: "nc" (tau), "c" (tau_mu), "ct" (tau_tau). Las dummies estacionales centradas no cambian la
# distribucion limite de tau (Dickey, Bell y Miller, 1986). p-valores: MacKinnon (1996), urca::punitroot.
adf <- function(x, p, det = c("ct", "c", "nc"), dummies = FALSE, v = "X") {
  det <- match.arg(det); dx <- diff(x); t0 <- (p + 1):length(dx)
  Z <- data.frame(DX = as.numeric(dx)[t0], X_1 = as.numeric(x)[t0])
  for (j in seq_len(p)) Z[[sprintf("DX_%d", j)]] <- as.numeric(dx)[t0 - j]
  if (det == "ct") Z$TREND <- t0
  if (dummies) Z <- cbind(Z, dummies_est(dx)[t0, , drop = FALSE])
  m <- if (det == "nc") lm(DX ~ . - 1, data = Z) else lm(DX ~ ., data = Z)
  cf <- summary(m)$coefficients; tau <- cf["X_1", 3]; n <- length(t0)
  vc <- sapply(c(0.01, 0.05, 0.10), function(a) urca::qunitroot(a, N = n, trend = det, statistic = "t"))
  list(m = m, tau = tau, p = urca::punitroot(tau, N = n, trend = det, statistic = "t"), vc = setNames(vc, c("1%", "5%", "10%")),
       t_trend = if (det == "ct") cf["TREND", 3] else NA_real_, n = n, det = det, lags = p, dummies = dummies, v = v,
       fechas = em(dx)[t0], LB_p = Box.test(residuals(m), 24, "Ljung-Box")$p.value)
}

# Correlograma con Q de Ljung-Box (Nota 6, Sec. 9.1-9.2). gl = K - k (k: parametros ARMA, para residuos)
correlograma <- function(z, K = 36, k = 0) {
  z <- as.numeric(z); n <- length(z)
  ac <- acf(z, K, plot = FALSE, demean = TRUE)$acf[-1]; pac <- pacf(z, K, plot = FALSE)$acf[, 1, 1]
  Q <- n * (n + 2) * cumsum(ac^2 / (n - 1:K)); gl <- (1:K) - k
  data.frame(k = 1:K, AC = ac, PAC = pac, Q = Q, gl = gl, p = ifelse(gl > 0, pchisq(Q, pmax(gl, 1), lower.tail = FALSE), NA), banda = 1.96 / sqrt(n), n = n)
}

jarque_bera <- function(e, k = 0) {
  e <- e - mean(e); n <- length(e); s <- sqrt(mean(e^2)); As <- mean(e^3) / s^3; Ku <- mean(e^4) / s^4
  JB <- (n - k) / 6 * (As^2 + (Ku - 3)^2 / 4)
  c(asimetria = As, curtosis = Ku, JB = JB, p = pchisq(JB, 2, lower.tail = FALSE))
}
# Raices de los polinomios AR y MA estimados (incluida la parte estacional), en la convencion z
raices <- function(m) {
  ar <- m$model$phi; ma <- m$model$theta
  list(AR = if (length(ar)) polyroot(c(1, -ar)) else complex(0), MA = if (length(ma)) polyroot(c(1, ma)) else complex(0))
}
# ARMA(p,q)x(P,Q)_12 por maxima verosimilitud exacta sobre la serie ya filtrada W
ajustar <- function(W, p, q, P = 0, Q = 0, xreg = NULL, media = TRUE) {
  tryCatch(suppressWarnings(arima(W, order = c(p, 0, q), seasonal = list(order = c(P, 0, Q), period = S), xreg = xreg,
                                  include.mean = media, method = "ML", optim.control = list(maxit = 1000))),
           error = function(e) NULL)
}
# Diagnostico de un modelo ajustado (Nota 6, Sec. 10)
diagnostico <- function(m, k) {
  e <- as.numeric(residuals(m)); rt <- raices(m)
  lb <- sapply(K_LB, function(K) Box.test(e, K, "Ljung-Box", fitdf = k)$p.value)
  jb <- jarque_bera(e, k)
  l <- as.numeric(logLik(m)); n <- length(e); kc <- length(m$coef)          # criterios por observacion, como EViews
  data.frame(n = n, sigma = sqrt(m$sigma2), logLik = l, AIC = AIC(m), SBC = BIC(m), AICn = (-2 * l + 2 * kc) / n, SCn = (-2 * l + kc * log(n)) / n,
             min_mod_AR = if (length(rt$AR)) min(Mod(rt$AR)) else NA, min_mod_MA = if (length(rt$MA)) min(Mod(rt$MA)) else NA,
             LB12 = lb[1], LB24 = lb[2], LB36 = lb[3], JB = jb["JB"], JB_p = jb["p"],
             LB2_24 = Box.test(e^2, 24, "Ljung-Box")$p.value, row.names = NULL)
}
# Rejilla de identificacion: todos los ARMA(p,q)x(P,Q) chicos sobre W, con diagnostico y admisibilidad
rejilla <- function(W, xreg = NULL, media = TRUE, pmax = 3, qmax = 2, Pmax = 1, Qmax = 1) {
  out <- list()
  for (p in 0:pmax) for (q in 0:qmax) for (P in 0:Pmax) for (Q in 0:Qmax) {
    m <- ajustar(W, p, q, P, Q, xreg, media); if (is.null(m) || m$code != 0) next
    d <- diagnostico(m, p + q + P + Q)
    out[[length(out) + 1]] <- data.frame(modelo = etiqueta_modelo(p, q, P, Q), p = p, q = q, P = P, Q = Q, d)
  }
  R <- do.call(rbind, out)
  R$estacionario <- is.na(R$min_mod_AR) | R$min_mod_AR > 1.01
  R$invertible <- is.na(R$min_mod_MA) | R$min_mod_MA > 1.01
  R$ruido_blanco <- R$LB12 > ALFA & R$LB24 > ALFA & R$LB36 > ALFA
  R$admisible <- R$estacionario & R$invertible & R$ruido_blanco
  R$dSBC <- R$SBC - min(R$SBC)
  R[order(R$SBC), ]
}
etiqueta_modelo <- function(p, q, P = 0, Q = 0) {
  r <- sprintf("ARMA(%d,%d)", p, q)
  if (P > 0 || Q > 0) r <- sprintf("%s(%d,%d)\\ensuremath{_{12}}", r, P, Q)
  r
}
# Sobreajuste (Nota 6, Sec. 10.4): agrega un parametro por vez y mira su t y el cambio en SBC
sobreajuste <- function(W, base, xreg = NULL, media = TRUE) {
  m0 <- ajustar(W, base[1], base[2], base[3], base[4], xreg, media)
  extras <- list(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
  nombre <- c("AR", "MA", "SAR", "SMA")
  do.call(rbind, lapply(1:4, function(i) {
    o <- base + extras[[i]]; m <- ajustar(W, o[1], o[2], o[3], o[4], xreg, media)
    if (is.null(m) || m$code != 0) return(NULL)
    cn <- c(sprintf("ar%d", o[1]), sprintf("ma%d", o[2]), sprintf("sar%d", o[3]), sprintf("sma%d", o[4]))[i]
    t <- m$coef[cn] / sqrt(m$var.coef[cn, cn])
    data.frame(agregado = sprintf("%s(%d)", nombre[i], c(o[1], o[2], 12 * o[3], 12 * o[4])[i]), coef = unname(m$coef[cn]),
               t = unname(t), dSBC = BIC(m) - BIC(m0), row.names = NULL)
  }))
}

# ---- 4. Salidas con formato EViews -------------------------------------------------------------
# Criterios de informacion POR OBSERVACION (como EViews y las filminas de Corfield, Nota 6 Sec. 8.2):
#   AIC = -2l/n + 2k/n ; SC = -2l/n + k ln(n)/n ; HQ = -2l/n + 2k ln(ln n)/n
# Covarianza de stats::arima: Hessiano numerico (EViews usa por defecto OPG); pueden diferir en ee.
fmt <- function(x, d = 6) sub("^-(0\\.0+)$", "\\1", formatC(x, format = "f", digits = d))
fmt_raiz <- function(z) {
  z <- z[order(-Mod(z), -Im(z))]; out <- character(0); usado <- logical(length(z))
  for (i in seq_along(z)) {
    if (usado[i]) next
    if (abs(Im(z[i])) < 1e-8) { out <- c(out, fmt(Re(z[i]), 2)); usado[i] <- TRUE
    } else {
      j <- which(!usado & seq_along(z) != i & abs(z - Conj(z[i])) < 1e-6)[1]; usado[c(i, j)] <- TRUE
      out <- c(out, sprintf("%s+/-%si", fmt(Re(z[i]), 2), fmt(abs(Im(z[i])), 2)))
    }
  }
  out
}
dw <- function(e) sum(diff(e)^2) / sum(e^2)

# Motor de dibujo: tabla de texto con lineas, al estilo de la ventana de EViews
dibujar <- function(titulo, cabecera, tab, izq, der, extra = NULL, previo = NULL, rojo_tab = NULL, rojo_stats = NULL, rojo_previo = FALSE) {
  if (is.null(rojo_tab)) rojo_tab <- rep(FALSE, nrow(tab))
  fuente <- "Liberation Sans"; serif <- "Latin Modern Roman"; tam <- 3.15
  el <- list(); seg <- list(); y <- 0
  add <- function(x, y, lab, h = 0.5, face = "plain", fam = fuente, s = tam, col = tinta)
    el[[length(el) + 1]] <<- data.frame(x = x, y = y, lab = lab, h = h, face = face, fam = fam, s = s, col = col)
  linea <- function(y) seg[[length(seg) + 1]] <<- data.frame(x = 0.5, xend = 99.5, y = y, yend = y)
  add(50, y, titulo, 0.5, "bold", serif, 5.0); y <- y - 2.1; y_top <- y + 0.65
  for (ln in cabecera) { add(2, y, ln, 0); y <- y - 1 }
  if (!is.null(previo)) {
    y <- y - 0.2; linea(y + 0.55)
    add(78, y - 0.1, "t-Statistic", 1, "bold", s = tam - 0.1); add(97, y - 0.1, "Prob.*", 1, "bold", s = tam - 0.1)
    y <- y - 1.15; linea(y + 0.5)
    for (i in seq_along(previo$filas)) { r <- previo$filas[[i]]; cc <- if (rojo_previo && i == 1) rojo else tinta; ff <- if (rojo_previo && i == 1) "bold" else "plain"
      add(2, y, r[1], 0, ff, col = cc); if (length(r) > 1) add(78, y, r[2], 1, ff, col = cc); if (length(r) > 2) add(97, y, r[3], 1, ff, col = cc); y <- y - 1 }
    y <- y - 0.1; linea(y + 0.5)
    for (ln in previo$post) { add(2, y, ln, 0); y <- y - 1 }
  }
  y <- y - 0.2; linea(y + 0.55); xs <- c(17, 43, 60, 78, 97)
  add(xs[1], y - 0.4, "Variable", 0.5, "bold", s = tam - 0.1)
  for (j in 2:5) add(xs[j], y - 0.4, c("Coefficient", "Std. Error", "t-Statistic", "Prob.")[j - 1], 1, "bold", s = tam - 0.1)
  y <- y - 1.35; linea(y + 0.5)
  for (i in seq_len(nrow(tab))) { cc <- if (rojo_tab[i]) rojo else tinta; ff <- if (rojo_tab[i]) "bold" else "plain"
    add(xs[1], y, tab$var[i], 0.5, ff, col = cc); for (j in 2:5) add(xs[j], y, tab[[j]][i], 1, ff, col = cc); y <- y - 1 }
  y <- y - 0.1; linea(y + 0.5)
  for (i in seq_along(izq$lab)) {
    ci <- if (izq$lab[i] %in% rojo_stats) rojo else tinta; cd <- if (der$lab[i] %in% rojo_stats) rojo else tinta
    fi <- if (izq$lab[i] %in% rojo_stats) "bold" else "plain"; fd <- if (der$lab[i] %in% rojo_stats) "bold" else "plain"
    add(2, y, izq$lab[i], 0, fi, col = ci); add(47, y, izq$val[i], 1, fi, col = ci); add(52, y, der$lab[i], 0, fd, col = cd); add(97, y, der$val[i], 1, fd, col = cd); y <- y - 1 }
  if (!is.null(extra)) { y <- y - 0.1; linea(y + 0.5); for (ln in extra) { add(2, y, ln$lab, 0); add(30, y, ln$val, 0); y <- y - 1 } }
  y_bot <- y + 0.45
  E <- do.call(rbind, el); Sg <- do.call(rbind, seg)
  p <- ggplot() +
    annotate("rect", xmin = 0.5, xmax = 99.5, ymin = y_bot, ymax = y_top, fill = "white", colour = "#6B7480", linewidth = 0.4) +
    geom_segment(data = Sg, aes(x = x, xend = xend, y = y, yend = yend), colour = "#6B7480", linewidth = 0.3) +
    geom_text(data = E, aes(x = x, y = y, label = lab, hjust = h, fontface = face, family = fam), size = E$s, colour = E$col) +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + scale_y_continuous(limits = c(y_bot - 0.5, 1.2), expand = c(0, 0)) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  attr(p, "alto") <- (1.2 - y_bot + 0.5) * 0.205
  p
}
nombre_coef <- function(v) {
  v <- sub("^intercept$", "C", v); v <- sub("^ar([0-9]+)$", "AR(\\1)", v); v <- sub("^ma([0-9]+)$", "MA(\\1)", v)
  v <- sub("^sar([0-9]+)$", "SAR(12)", v); sub("^sma([0-9]+)$", "SMA(12)", v)
}
# Modelo ARMA estimado con stats::arima sobre W (serie ya filtrada)
salida_arma <- function(m, W, dep, titulo, rojo_tab = TRUE, rojo_stats = c("Schwarz criterion")) {
  e <- as.numeric(residuals(m)); n <- length(e); b <- m$coef; k <- length(b); l <- as.numeric(logLik(m))
  se <- sqrt(diag(m$var.coef))[names(b)]; tt <- b / se; p <- 2 * pt(-abs(tt), df = n - k)
  tab <- data.frame(var = nombre_coef(names(b)), coef = fmt(b), se = fmt(se), t = fmt(tt), p = fmt(p, 4))
  ssr <- sum(e^2); sst <- sum((as.numeric(W) - mean(W))^2); r2 <- 1 - ssr / sst; r2a <- 1 - (1 - r2) * (n - 1) / (n - k)
  ic <- c(-2 * l / n + 2 * k / n, -2 * l / n + k * log(n) / n, -2 * l / n + 2 * k * log(log(n)) / n)
  f <- em(W)
  cab <- c(paste("Dependent Variable:", dep), "Method: ARMA Maximum Likelihood (exact, BFGS)",
           sprintf("Sample: %s %s", f[1], f[n]), sprintf("Included observations: %d", n),
           if (m$code == 0) "Convergence achieved" else "WARNING: convergence not achieved",
           "Coefficient covariance computed using the observed Hessian")
  izq <- list(lab = c("R-squared", "Adjusted R-squared", "S.E. of regression", "Sum squared resid", "Log likelihood", "Jarque-Bera (Prob.)"),
              val = c(fmt(r2), fmt(r2a), fmt(sqrt(ssr / (n - k))), fmt(ssr), fmt(l), sprintf("%s (%s)", fmt(jarque_bera(e, k)["JB"], 3), fmt(jarque_bera(e, k)["p"], 4))))
  der <- list(lab = c("Mean dependent var", "S.D. dependent var", "Akaike info criterion", "Schwarz criterion", "Hannan-Quinn criter.", "Durbin-Watson stat"),
              val = c(fmt(mean(W)), fmt(sd(W)), fmt(ic[1]), fmt(ic[2]), fmt(ic[3]), fmt(dw(e))))
  rt <- raices(m)
  reng <- function(lab, z) {
    v <- if (length(z)) fmt_raiz(1 / z) else "(none)"; g <- split(v, ceiling(seq_along(v) / 4))
    lapply(seq_along(g), function(i) list(lab = if (i == 1) lab else "", val = paste(g[[i]], collapse = "     ")))
  }
  dibujar(titulo, cab, tab, izq, der, c(reng("Inverted AR Roots", rt$AR), reng("Inverted MA Roots", rt$MA)),
          rojo_tab = rep_len(rojo_tab, nrow(tab)), rojo_stats = rojo_stats)
}
# Regresion lm()
salida_mco <- function(m, dep, titulo, muestra, previo = NULL, cab_previa = NULL, nombres = NULL, rojo_var = NULL, rojo_stats = NULL, rojo_previo = FALSE) {
  s <- summary(m); cf <- s$coefficients; e <- residuals(m); n <- length(e); k <- nrow(cf)
  nm <- sub("^\\(Intercept\\)$", "C", rownames(cf)); if (!is.null(nombres)) nm <- unname(ifelse(is.na(nombres[nm]), nm, nombres[nm]))
  tab <- data.frame(var = nm, coef = fmt(cf[, 1]), se = fmt(cf[, 2]), t = fmt(cf[, 3]), p = fmt(cf[, 4], 4))
  y <- model.response(model.frame(m)); l <- as.numeric(logLik(m)) ; kk <- length(coef(m))
  ic <- c(-2 * l / n + 2 * kk / n, -2 * l / n + kk * log(n) / n, -2 * l / n + 2 * kk * log(log(n)) / n)
  fs <- s$fstatistic
  cab <- c(cab_previa, paste("Dependent Variable:", dep), "Method: Least Squares", paste("Sample:", muestra), sprintf("Included observations: %d", n))
  izq <- list(lab = c("R-squared", "Adjusted R-squared", "S.E. of regression", "Sum squared resid", "Log likelihood", "F-statistic", "Prob(F-statistic)"),
              val = c(fmt(s$r.squared), fmt(s$adj.r.squared), fmt(s$sigma), fmt(sum(e^2)), fmt(l),
                      if (is.null(fs)) "NA" else fmt(fs[1]), if (is.null(fs)) "NA" else fmt(pf(fs[1], fs[2], fs[3], lower.tail = FALSE))))
  der <- list(lab = c("Mean dependent var", "S.D. dependent var", "Akaike info criterion", "Schwarz criterion", "Hannan-Quinn criter.", "Durbin-Watson stat", ""),
              val = c(fmt(mean(y)), fmt(sd(y)), fmt(ic[1]), fmt(ic[2]), fmt(ic[3]), fmt(dw(e)), ""))
  dibujar(titulo, cab, tab, izq, der, previo = previo, rojo_tab = tab$var %in% rojo_var, rojo_stats = rojo_stats, rojo_previo = rojo_previo)
}
# Test ADF completo (cabecera + ecuacion del test)
salida_adf <- function(A, titulo) {
  cf <- summary(A$m)$coefficients; v <- A$v; p <- A$lags
  nombres <- c(X_1 = sprintf("%s(-1)", v), setNames(sprintf("D(%s(-%d))", v, seq_len(p)), sprintf("DX_%d", seq_len(p))), TREND = "@TREND")
  exo <- c(nc = "None", c = "Constant", ct = "Constant, Linear Trend")[A$det]
  if (A$dummies) exo <- paste0(exo, ", Seasonal dummies (Jan-Nov)")
  previo <- list(filas = list(c("Augmented Dickey-Fuller test statistic", fmt(A$tau), fmt(A$p, 4)), c("Test critical values:"),
                              c("      1% level", fmt(A$vc["1%"])), c("      5% level", fmt(A$vc["5%"])), c("     10% level", fmt(A$vc["10%"]))),
                 post = "*MacKinnon (1996) one-sided p-values.")
  cab <- c(sprintf("Null Hypothesis: %s has a unit root", v), paste("Exogenous:", exo), sprintf("Lag Length: %d (fixed)", p))
  salida_mco(A$m, sprintf("D(%s)", v), titulo, sprintf("%s %s", A$fechas[1], A$fechas[length(A$fechas)]),
             previo = previo, cab_previa = cab, nombres = nombres, rojo_var = sprintf("%s(-1)", v), rojo_previo = TRUE)
}
# Correlograma con el formato de EViews: barras de FAC y FACP, AC, PAC, Q-Stat y Prob
salida_correlograma <- function(z, titulo, K = 24, k = 0, sub = NULL, rojo_k = NULL) {
  C <- correlograma(z, K, k); b <- C$banda[1]; n <- C$n[1]
  fuente <- "Liberation Mono"; tam <- 2.9
  y <- -(1:K)                                                     # una fila por rezago
  x_ac <- 12; x_pac <- 33; w <- 10.5                              # centro y semiancho de cada panel de barras
  bar <- function(x0, r) data.frame(xmin = pmin(x0, x0 + r * w), xmax = pmax(x0, x0 + r * w), ymin = y - 0.32, ymax = y + 0.32,
                                    sig = ifelse(C$k %in% rojo_k, "rojo", ifelse(abs(r) > b, "si", "no")))
  Bc <- rbind(bar(x_ac, C$AC), bar(x_pac, C$PAC))
  txt <- data.frame(y = y, k = C$k, AC = fmt(C$AC, 3), PAC = fmt(C$PAC, 3), Q = fmt(C$Q, 3),
                    p = ifelse(is.na(C$p), "", fmt(C$p, 3)))
  txt$col <- ifelse(C$k %in% rojo_k, rojo, tinta); txt$ff <- ifelse(C$k %in% rojo_k, "bold", "plain")
  p <- ggplot() +
    annotate("rect", xmin = 0.5, xmax = 99.5, ymin = -K - 0.9, ymax = 3.6, fill = "white", colour = "#6B7480", linewidth = 0.4) +
    annotate("text", x = 50, y = 5.2, label = titulo, family = "Latin Modern Roman", fontface = "bold", size = 5.0, colour = tinta) +
    annotate("text", x = 2, y = 2.6, hjust = 0, size = 3.1, family = "Liberation Sans", colour = tinta,
             label = sprintf("Sample: %s   Included observations: %d%s", paste(range(em(z)), collapse = " "), n,
                             if (k > 0) sprintf("   Q-statistic probabilities adjusted for %d ARMA term(s)", k) else "")) +
    annotate("segment", x = 0.5, xend = 99.5, y = 1.9, yend = 1.9, colour = "#6B7480", linewidth = 0.3) +
    annotate("segment", x = 0.5, xend = 99.5, y = 0.45, yend = 0.45, colour = "#6B7480", linewidth = 0.3) +
    annotate("text", x = c(x_ac, x_pac), y = 1.15, label = c("Autocorrelation", "Partial Correlation"), size = 3.0, family = "Liberation Sans", fontface = "bold", colour = tinta) +
    annotate("text", x = c(47, 60, 73, 86, 97), y = 1.15, label = c("k", "AC", "PAC", "Q-Stat", "Prob"), hjust = 1, size = 3.0,
             family = "Liberation Sans", fontface = "bold", colour = tinta) +
    annotate("segment", x = c(x_ac, x_pac), xend = c(x_ac, x_pac), y = 0.3, yend = -K - 0.6, colour = tinta2, linewidth = 0.3) +
    annotate("segment", x = c(x_ac - b * w, x_ac + b * w, x_pac - b * w, x_pac + b * w), xend = c(x_ac - b * w, x_ac + b * w, x_pac - b * w, x_pac + b * w),
             y = 0.3, yend = -K - 0.6, colour = "#9AA5B1", linewidth = 0.3, linetype = "22") +
    geom_rect(data = Bc, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = sig)) +
    scale_fill_manual(values = c(si = col_viv, no = "#B8C2CC", rojo = rojo), guide = "none") +
    geom_text(data = txt, aes(x = 47, y = y, label = k, fontface = ff), hjust = 1, size = tam, family = fuente, colour = txt$col) +
    geom_text(data = txt, aes(x = 60, y = y, label = AC, fontface = ff), hjust = 1, size = tam, family = fuente, colour = txt$col) +
    geom_text(data = txt, aes(x = 73, y = y, label = PAC, fontface = ff), hjust = 1, size = tam, family = fuente, colour = txt$col) +
    geom_text(data = txt, aes(x = 86, y = y, label = Q, fontface = ff), hjust = 1, size = tam, family = fuente, colour = txt$col) +
    geom_text(data = txt, aes(x = 97, y = y, label = p, fontface = ff), hjust = 1, size = tam, family = fuente, colour = txt$col) +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + scale_y_continuous(limits = c(-K - 1.4, 6.3), expand = c(0, 0)) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  if (!is.null(sub)) p <- p + annotate("text", x = 2, y = -K - 1.15, hjust = 0, vjust = 1, size = 2.7, family = "Liberation Sans", colour = tinta2, label = sub)
  attr(p, "alto") <- (K + 7.7) * 0.205
  p
}
guardar_salida <- function(p, archivo, ancho = 7.4)
  ggsave(file.path("salidas/eviews", archivo), p, width = ancho, height = attr(p, "alto"), dpi = 220, device = ragg::agg_png)

# ---- 5. Exportar a LaTeX ---------------------------------------------------------------------------
# Cada script escribe sus valores en salidas/tablas/valores_XX.tex como \newcommand{\Nombre}{valor};
# el informe los lee con \input, asi los numeros del texto salen siempre del codigo.
num <- function(x, d = 2) {                 # formato espanol: coma decimal, signo menos matematico
  x <- round(as.numeric(x), d) + 0         # evita "-0"
  s <- formatC(x, format = "f", digits = d, big.mark = ""); s <- sub("^-", "\\\\ensuremath{-}", s); gsub("\\.", "{,}", s)
}
pval <- function(p) if (is.na(p)) "--" else if (p < 0.001) "\\ensuremath{<0{,}001}" else num(p, 3)
rojo_tex <- function(x) sprintf("\\textcolor{red}{%s}", x)        # resalta una celda de tabla LaTeX
marcar <- function(D, filas = NULL, cols = NULL) {                    # resalta filas enteras o celdas [filas, cols]
  D[] <- lapply(D, as.character)
  if (is.null(cols)) cols <- seq_along(D)
  for (i in filas) for (j in cols) D[i, j] <- rojo_tex(D[i, j])
  D
}
valores <- new.env()
guardar_valor <- function(nombre, valor) assign(nombre, valor, envir = valores)
escribir_valores <- function(archivo) {
  nm <- sort(ls(valores))
  writeLines(sprintf("\\newcommand{\\%s}{%s}", nm, sapply(nm, function(n) get(n, envir = valores))), file.path("salidas/tablas", archivo))
  rm(list = nm, envir = valores)
}
# Tabla booktabs a partir de un data.frame ya formateado (columnas de texto)
celdas <- function(filas, cols) as.matrix(expand.grid(filas, cols))     # celdas a resaltar en rojo
tabla_tex <- function(D, archivo, alineacion, encabezado = names(D), rojo = NULL) {
  D[] <- lapply(D, as.character)
  if (!is.null(rojo)) for (i in seq_len(nrow(rojo))) D[rojo[i, 1], rojo[i, 2]] <- rojo_tex(D[rojo[i, 1], rojo[i, 2]])
  filas <- apply(D, 1, function(r) paste(r, collapse = " & "))
  writeLines(c(sprintf("\\begin{tabular}{%s}", alineacion), "\\toprule", paste(encabezado, collapse = " & "), "\\\\ \\midrule",
               paste0(filas, " \\\\"), "\\bottomrule", "\\end{tabular}"), file.path("salidas/tablas", archivo))
}
