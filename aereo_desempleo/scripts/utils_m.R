# =============================================================================
# utils_m.R -- Version mensual de las herramientas de scripts/utils.R (se usa con source())
# Reutiliza el estilo, jarque_bera() y modulos() del TP principal, y redefine lo que dependia de la
# frecuencia trimestral: dummies, ADF, ajuste con MA estacional de rezago s = frequency(x), rejilla y
# Ljung-Box en 12, 24 y 36 rezagos.
# =============================================================================
source("scripts/utils.R"); source("scripts/eviews.R")
col_viv <- "#2A62B8"; col_tasa <- "#D9701F"; col_bono <- "#9AA5B1"       # inicios azul, tasa naranja

cargar_base <- function() {
  B <- read.csv("datos/base_mensual.csv", stringsAsFactors = FALSE)
  ini <- c(as.integer(substr(B$mes[1], 1, 4)), as.integer(substr(B$mes[1], 6, 7)))
  list(B = B,
       viv  = ts(100 * log(B$inicios_nsa), start = ini, frequency = 12),     # ln inicios (x100), NSA
       viv_sa = ts(100 * log(B$inicios_sa_anual / 12), start = ini, frequency = 12),
       tasa = ts(B$tasa_hip_30, start = ini, frequency = 12),               # % anual
       bono = ts(B$bono10, start = ini, frequency = 12),
       rec  = ts(B$recesion, start = ini, frequency = 12))
}
etiqueta_mes <- function(x) { t <- as.numeric(time(x)) + 1e-6; sprintf("%d-%02d", floor(t), floor((t - floor(t)) * 12) + 1) }
fecha_num <- function(a, m) a + (m - 1) / 12

# Dummies estacionales centradas (mes j menos el ultimo mes): suman cero en el anio
dummies_est <- function(x) {
  s <- frequency(x); q <- cycle(x)
  D <- sapply(1:(s - 1), function(j) (q == j) - (q == s)); colnames(D) <- if (s == 12) month.abb[1:11] else paste0("Q", 1:(s - 1)); D
}
adf_est <- function(x, p, tendencia = TRUE, dummies = TRUE) {
  dx <- diff(x); t0 <- (p + 1):length(dx)
  X <- cbind(nivel_rez = as.numeric(x)[t0])
  for (j in seq_len(p)) X <- cbind(X, dx[t0 - j])
  if (tendencia) X <- cbind(X, tend = t0)
  if (dummies) X <- cbind(X, dummies_est(dx)[t0, ])
  m <- lm(dx[t0] ~ X); s <- summary(m)$coefficients
  c(tau = s["Xnivel_rez", 3], t_tend = if (tendencia) s["Xtend", 3] else NA_real_,
    LB_p = Box.test(residuals(m), 36, "Ljung-Box")$p.value, n = length(t0))
}
# Valores criticos de Dickey-Fuller para n = 500 (Fuller 1976; los de urca::ur.df para esa fila)
VC_ADF <- rbind(tau_mu = c(`1%` = -3.44, `5%` = -2.87, `10%` = -2.57),
                tau_tau = c(`1%` = -3.98, `5%` = -3.42, `10%` = -3.13))

P_MAX <- 4; Q_MAX <- 2; ALFA <- 0.05; UMBRAL_CANCEL <- 0.10
K_LB <- c(12, 24, 36)

ajustar <- function(x, p, q, Qs, xreg = NULL, media = TRUE) {
  tryCatch(suppressWarnings(
    arima(x, order = c(p, 0, q), seasonal = list(order = c(0, 0, Qs), period = frequency(x)), xreg = xreg,
          include.mean = media, method = "ML", optim.control = list(maxit = 500))),
    error = function(e) NULL)
}
evaluar <- function(m, x, p, q, Qs) {
  k <- p + q + Qs; e <- as.numeric(residuals(m)); rt <- modulos(m)
  Hreg <- if (q > 0) polyroot(c(1, m$coef[paste0("ma", 1:q)])) else complex(0)
  cancel <- if (length(rt$G) && length(Hreg)) min(abs(outer(rt$G, Hreg, "-"))) else Inf
  lb <- sapply(K_LB, function(K) Box.test(e, K, "Ljung-Box", fitdf = k)$p.value)
  data.frame(p = p, q = q, Q_est = Qs, k = k, logLik = as.numeric(logLik(m)), AIC = AIC(m), SBC = BIC(m),
             sigma2 = m$sigma2, min_mod_G = if (length(rt$G)) min(Mod(rt$G)) else Inf,
             min_mod_H = if (length(rt$H)) min(Mod(rt$H)) else Inf, dist_AR_MA = cancel,
             LB12 = lb[1], LB24 = lb[2], LB36 = lb[3], JB_p = jarque_bera(e, k)["p"], row.names = NULL)
}
rejilla <- function(x, xreg = NULL, media = TRUE, Qs_max = 1) {
  cand <- list()
  for (p in 0:P_MAX) for (q in 0:Q_MAX) for (Qs in 0:Qs_max) {
    m <- ajustar(x, p, q, Qs, xreg, media); if (is.null(m) || m$code != 0) next
    cand[[length(cand) + 1]] <- evaluar(m, x, p, q, Qs)
  }
  C <- do.call(rbind, cand)
  C$estacionario <- C$min_mod_G > 1.001; C$invertible <- C$min_mod_H > 1.001
  C$sin_factor_comun <- C$dist_AR_MA > UMBRAL_CANCEL
  C$LB_ok <- with(C, LB12 > ALFA & LB24 > ALFA & LB36 > ALFA)
  C$admisible <- with(C, estacionario & invertible & sin_factor_comun & LB_ok)
  C
}
etiqueta_modelo <- function(p, q, Qs, s = 12) sprintf("ARMA(%d,%d)%s", p, q, ifelse(Qs == 1, sprintf("+SMA(%d)", s), ""))

# Sombreado de recesiones NBER para ggplot
recesiones <- function(rec) {
  r <- as.numeric(rec); t <- as.numeric(time(rec)); d <- diff(c(0, r, 0))
  data.frame(ini = t[which(d == 1)], fin = t[which(d == -1) - 1] + 1 / 12)
}
capa_recesion <- function(R) annotate("rect", xmin = R$ini, xmax = R$fin, ymin = -Inf, ymax = Inf, fill = "#E4E7EB", alpha = 0.7)
