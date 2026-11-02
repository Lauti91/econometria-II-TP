# =============================================================================
# utils.R  --  Estilo de graficos y carga de datos compartidos (se usa con source())
# =============================================================================

for (pk in c("ggplot2", "patchwork")) {
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ---- Estilo ---------------------------------------------------------------------
# Un color por variable (PBI azul, electricidad naranja), el mismo en todos los graficos.
col_pbi <- "#2A62B8"; col_elec <- "#D9701F"; col_rail <- "#14946A"   # validados (CVD y contraste)
tinta <- "#1F2933"; tinta2 <- "#52606D"; gris <- "#C5CCD3"; fondo <- "#FCFCFB"

tema <- theme_minimal(base_size = 11) +
  theme(
    plot.background    = element_rect(fill = fondo, colour = NA),
    panel.background   = element_rect(fill = fondo, colour = NA),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E4E7EB", linewidth = 0.3),
    axis.text  = element_text(colour = tinta2, size = 9),
    axis.title = element_text(colour = tinta2, size = 9),
    plot.title = element_text(colour = tinta, face = "bold", size = 11, margin = margin(b = 2)),
    plot.subtitle = element_text(colour = tinta2, size = 9, margin = margin(b = 6)),
    plot.margin = margin(8, 12, 6, 8)
  )

tema_global <- theme(
  plot.background = element_rect(fill = fondo, colour = NA),
  plot.title   = element_text(colour = tinta, face = "bold", size = 14),
  plot.caption = element_text(colour = tinta2, size = 8, hjust = 0))

# ---- Datos ----------------------------------------------------------------------
# Muestra principal 2000Q1-2023Q4 (96 trimestres): crecimiento interanual (%) del PBI oficial
# y de la electricidad, construidos por 01_construir_base.R.
cargar_yoy <- function() {
  Q <- read.csv("output/trimestral_li_keqiang.csv", stringsAsFactors = FALSE)
  S <- Q[Q$trimestre >= "2000Q1" & Q$trimestre <= "2023Q4", ]
  stopifnot(nrow(S) == 96, !anyNA(S$pbi_crec_interanual_oficial_pct), !anyNA(S$elec_crec_interanual_pct))
  list(
    "PBI oficial, crec. interanual (%)"  = ts(S$pbi_crec_interanual_oficial_pct,  start = c(2000, 1), frequency = 4),
    "Electricidad, crec. interanual (%)" = ts(S$elec_crec_interanual_pct,         start = c(2000, 1), frequency = 4))
}

etiqueta_trim <- function(x) {
  t <- as.numeric(time(x)); sprintf("%dQ%d", floor(t + 1e-6), round((t - floor(t + 1e-6)) * 4) + 1)
}

# Log-niveles (x 100, para que las diferencias se lean en %): PBI real y electricidad, 2000Q1-2023Q4;
# carga ferroviaria (10.000 t, NBS), 2005Q1-2023Q4 (la serie mensual arranca en 2005-01 y desde 2024
# el NBS publica enero y febrero acumulados, por eso el primer trimestre de 2024 en adelante falta).
# PBI real: el publicado (base 2011) solo existe desde 2011 y cambia de base en 2016 y 2021, asi que
# se encadena la tasa interanual oficial g_t por trimestre, R_t = R_{t-4} (1 + g_t/100), hacia
# adelante y hacia atras desde los cuatro niveles reales publicados de 2011 (semilla). Control: en
# 2011-2015 (misma base que la semilla) el encadenado debe reproducir el publicado.
cargar_niveles <- function() {
  Q <- read.csv("output/trimestral_li_keqiang.csv", stringsAsFactors = FALSE)
  g <- Q$pbi_crec_interanual_oficial_pct; n <- nrow(Q)
  i11 <- match(paste0("2011Q", 1:4), Q$trimestre)
  R <- rep(NA_real_, n); R[i11] <- Q$pbi_real_1e8_yuan_2011[i11]
  for (i in (max(i11) + 1):n) R[i] <- R[i - 4] * (1 + g[i] / 100)
  for (i in (min(i11) - 1):1) if (!is.na(g[i + 4])) R[i] <- R[i + 4] / (1 + g[i + 4] / 100)
  k <- Q$trimestre >= "2011Q1" & Q$trimestre <= "2015Q4"
  dif <- 100 * (R[k] / Q$pbi_real_1e8_yuan_2011[k] - 1)
  stopifnot(max(abs(dif)) < 0.25)          # el redondeo a 0,1 pp de la tasa oficial explica el resto
  S <- Q$trimestre >= "2000Q1" & Q$trimestre <= "2023Q4"
  F <- Q$trimestre >= "2005Q1" & Q$trimestre <= "2023Q4"
  stopifnot(sum(S) == 96, !anyNA(R[S]), !anyNA(Q$electricidad_GWh[S]), sum(F) == 76, !anyNA(Q$carga_ferroviaria_10k_t[F]))
  list(niveles = data.frame(trimestre = Q$trimestre[S], pbi_real = R[S], electricidad_GWh = Q$electricidad_GWh[S],
                            carga_ferroviaria_10k_t = Q$carga_ferroviaria_10k_t[S]),
       control_max_dif_pct = max(abs(dif)),
       "ln PBI real (x100)"     = ts(100 * log(R[S]), start = c(2000, 1), frequency = 4),
       "ln Electricidad (x100)" = ts(100 * log(Q$electricidad_GWh[S]), start = c(2000, 1), frequency = 4),
       "ln Carga ferroviaria (x100)" = ts(100 * log(Q$carga_ferroviaria_10k_t[F]), start = c(2005, 1), frequency = 4))
}

# Dummies estacionales centradas (Q1..Q3 menos Q4): suman cero en el anio, por eso conviven con la
# constante y no alteran la distribucion limite del Dickey-Fuller (Dickey, Bell y Miller, 1986).
dummies_est <- function(x) {
  q <- cycle(x); D <- sapply(1:3, function(j) (q == j) - (q == 4)); colnames(D) <- paste0("Q", 1:3); D
}

# ADF (Nota 6, Sec. 4.4) con dummies estacionales centradas y, opcionalmente, tendencia.
# Regresion: Dx_t = a + b t + g x_{t-1} + sum_j c_j Dx_{t-j} + dummies + e_t; tau = t de g.
adf_est <- function(x, p, tendencia = TRUE, dummies = TRUE) {
  dx <- diff(x); t0 <- (p + 1):length(dx)
  X <- cbind(nivel_rez = as.numeric(x)[t0])          # x[t0] es x_{t-1} respecto de dx[t0]
  for (j in seq_len(p)) X <- cbind(X, dx[t0 - j])
  if (tendencia) X <- cbind(X, tend = t0)
  if (dummies) X <- cbind(X, dummies_est(dx)[t0, ])
  m <- lm(dx[t0] ~ X); s <- summary(m)$coefficients
  c(tau = s["Xnivel_rez", 3], t_tend = if (tendencia) s["Xtend", 3] else NA_real_,
    LB8_p = Box.test(residuals(m), 8, "Ljung-Box")$p.value, n = length(t0))
}
# Valores criticos (MacKinnon, n ~ 100), los mismos que reporta urca::ur.df
VC_ADF <- rbind(tau_mu = c(`1%` = -3.51, `5%` = -2.89, `10%` = -2.58),
                tau_tau = c(`1%` = -4.04, `5%` = -3.45, `10%` = -3.15))

# ---- Herramientas Box-Jenkins (usadas por 03 y 04) -----------------------------------------
# Convencion de signos: Theta(B) = 1 - theta_1 B - ... (Notas y BJRL). R usa 1 + theta B.
# Parametros del ejercicio
P_MAX <- 4          # rezago maximo AR: cota de Said-Dickey n^(1/3) ~ 4 (Nota 6, Sec. 4.6)
Q_MAX <- 2          # rezago maximo MA regular
ALFA  <- 0.05
K_LB  <- c(8, 12, 16)   # rezagos del Ljung-Box
UMBRAL_CANCEL <- 0.10   # criterio propio: raiz AR y raiz MA a menos de 0.10 => posible factor comun

ajustar <- function(x, p, q, Qs, xreg = NULL, media = TRUE) {
  tryCatch(suppressWarnings(
    arima(x, order = c(p, 0, q), seasonal = list(order = c(0, 0, Qs), period = 4), xreg = xreg,
          include.mean = media, method = "ML")),
    error = function(e) NULL)
}

modulos <- function(m) {
  ar <- m$model$phi; ma <- m$model$theta     # polinomios ya expandidos (incluida la parte estacional)
  list(G = if (length(ar)) polyroot(c(1, -ar)) else complex(0),
       H = if (length(ma)) polyroot(c(1,  ma)) else complex(0))
}

jarque_bera <- function(e, k) {
  e <- e - mean(e); n <- length(e); s <- sqrt(mean(e^2))
  As <- mean(e^3) / s^3; Ku <- mean(e^4) / s^4
  JB <- (n - k) / 6 * (As^2 + (Ku - 3)^2 / 4)     # Nota 3, Sec. 6.2 (n-p) y Nota 4 (n-p-q)
  c(As = As, K = Ku, JB = JB, p = pchisq(JB, 2, lower.tail = FALSE))
}

evaluar <- function(m, x, p, q, Qs) {
  k <- p + q + Qs                      # parametros ARMA estimados (sin la media)
  e <- as.numeric(residuals(m)); rt <- modulos(m)
  # Factor comun (Nota 4, Sec. 3.2): se compara la parte AR con la parte MA *regular*. El MA
  # estacional no puede cancelarse con la parte AR porque no hay AR estacional en la rejilla.
  Hreg   <- if (q > 0) polyroot(c(1, m$coef[paste0("ma", 1:q)])) else complex(0)
  cancel <- if (length(rt$G) && length(Hreg)) min(abs(outer(rt$G, Hreg, "-"))) else Inf
  lb <- sapply(K_LB, function(K) if (K > k) Box.test(e, K, "Ljung-Box", fitdf = k)$p.value else NA_real_)
  jb <- jarque_bera(e, k)
  data.frame(
    p = p, q = q, Q_est = Qs, k = k, logLik = as.numeric(logLik(m)), AIC = AIC(m), SBC = BIC(m),
    sigma2 = m$sigma2,
    min_mod_G = if (length(rt$G)) min(Mod(rt$G)) else Inf,       # estacionariedad: > 1
    min_mod_H = if (length(rt$H)) min(Mod(rt$H)) else Inf,       # invertibilidad: > 1
    dist_AR_MA = cancel,
    LB8 = lb[1], LB12 = lb[2], LB16 = lb[3], JB_p = jb["p"], row.names = NULL)
}


# Ajusta la rejilla ARMA(p<=P_MAX, q<=Q_MAX) con y sin MA estacional (rezago 4) sobre x y marca
# los modelos admisibles (estacionario, invertible, sin factor comun, Ljung-Box > ALFA en K_LB).
rejilla <- function(x, xreg = NULL, media = TRUE) {
  cand <- list()
  for (p in 0:P_MAX) for (q in 0:Q_MAX) for (Qs in 0:1) {
    m <- ajustar(x, p, q, Qs, xreg, media); if (is.null(m)) next
    cand[[length(cand) + 1]] <- evaluar(m, x, p, q, Qs)
  }
  C <- do.call(rbind, cand)
  C$estacionario <- C$min_mod_G > 1
  C$invertible   <- C$min_mod_H > 1
  C$sin_factor_comun <- C$dist_AR_MA > UMBRAL_CANCEL
  C$LB_ok <- with(C, !is.na(LB8) & !is.na(LB12) & !is.na(LB16) & LB8 > ALFA & LB12 > ALFA & LB16 > ALFA)
  C$admisible <- with(C, estacionario & invertible & sin_factor_comun & LB_ok)
  C
}
