# =============================================================================
# 08_salidas_eviews.R -- Salidas "tipo EViews" de los modelos ganadores de 05 (muestra principal:
# completa 2000Q1-2023Q4, con atipicos aditivos 2008Q4 y 2020Q1) y de la ecuacion ADF.
# Ejecutar desde china_li_keqiang/: Rscript scripts/08_salidas_eviews.R
# Convencion: como en EViews, el MA entra con signo +theta (en las Notas es Theta(B) = 1 - theta B,
# de modo que el theta de las Notas es el opuesto del que se lee aca).
# =============================================================================
source("scripts/utils.R"); source("scripts/eviews.R")
dir.create("output/salidas_eviews", showWarnings = FALSE, recursive = TRUE)

L <- cargar_niveles()
X <- L[c("ln PBI real (x100)", "ln Electricidad (x100)", "ln Carga ferroviaria (x100)")]
corto <- c("PBI real", "Electricidad", "Ferrocarril"); names(corto) <- names(X)
dep_x <- c("PBI_REAL", "ELEC", "RAIL"); names(dep_x) <- corto      # log-niveles (x100)
inicio <- function(x) time(x)[6]
dummies_est <- function(x) { q <- cycle(x); D <- sapply(1:3, function(j) (q == j) - (q == 4)); colnames(D) <- paste0("Q", 1:3); D }
filtrar <- function(z, f) switch(substr(f, 1, 3), "F1:" = diff(z), "F2:" = diff(z, 4), "F2t" = diff(z, 4), "F3:" = diff(diff(z, 4)))
nombre_dep <- function(v, f) switch(substr(f, 1, 3), "F1:" = sprintf("D(L%s)", v), "F2:" = sprintf("D(L%s,0,4)", v),
                                    "F2t" = sprintf("D(L%s,0,4)", v), "F3:" = sprintf("D(L%s,1,4)", v))
construir <- function(x, f, ao = NULL) {          # idem 05_log_niveles.R
  W <- window(filtrar(x, f), start = inicio(x))
  xr <- switch(substr(f, 1, 3), "F1:" = dummies_est(W), "F2t" = cbind(tend = seq_along(W)), NULL)
  for (a in ao) {
    I <- ts(as.numeric(abs(time(x) - a) < 1e-6), start = start(x), frequency = 4)
    z <- matrix(as.numeric(window(filtrar(I, f), start = inicio(x))), ncol = 1,
                dimnames = list(NULL, sprintf("AO%dQ%d", floor(a), round((a %% 1) * 4) + 1)))
    xr <- if (is.null(xr)) z else cbind(as.matrix(xr), z)
  }
  list(W = W, xreg = if (is.null(xr)) NULL else as.matrix(xr))
}

FIL <- read.csv("output/tablas/05_filtros.csv", stringsAsFactors = FALSE)
FIL <- FIL[FIL$muestra == "completa + AO 2008Q4 y 2020Q1", ]
ao <- c(2008.75, 2020)

# Ganador por serie (menor SBC) + el rival empatado de la electricidad (F2t, dSBC < 2)
pedidos <- do.call(rbind, lapply(split(FIL, FIL$serie), function(d) d[d$dSBC < 2 | d$dSBC == 0, ]))
pedidos <- pedidos[order(match(pedidos$serie, corto), pedidos$dSBC), ]
n_salida <- 0
for (i in seq_len(nrow(pedidos))) {
  g <- pedidos[i, ]; nm <- names(corto)[corto == g$serie]
  x <- window(X[[nm]], start = start(X[[nm]]), end = c(2023, 4)); b <- construir(x, g$filtro, ao)
  m <- ajustar(b$W, g$p, g$q, g$Qs, b$xreg); stopifnot(!is.null(m))
  n_salida <- n_salida + 1
  tit <- sprintf("Salida %d. %s: %s, %s", n_salida, g$serie, sub("^F[0-9t]+: ", "", g$filtro), g$modelo)
  S <- salida_arma(m, b$W, nombre_dep(dep_x[g$serie], g$filtro), tit)
  guardar_salida(S, sprintf("output/salidas_eviews/08_%d_%s_%s.png", n_salida, c("PBI real" = "pbi", "Electricidad" = "elec", "Ferrocarril" = "rail")[g$serie],
                            sub(":.*", "", g$filtro)))
}

# ---- Ecuacion del test ADF (PBI real, completa; p = 4, con tendencia y dummies) -----------------------
# Misma regresion que adf_est(): dx_t = a + b t + rho x_{t-1} + sum_j g_j dx_{t-j} + dummies + e_t
adf_lm <- function(x, p, tendencia = TRUE, dummies = TRUE) {
  dx <- diff(x); t0 <- (p + 1):length(dx); D <- data.frame(DX = as.numeric(dx)[t0], X_1 = as.numeric(x)[t0])
  for (j in seq_len(p)) D[[sprintf("DX_%d", j)]] <- as.numeric(dx)[t0 - j]
  if (tendencia) D$TREND <- t0
  if (dummies) D <- cbind(D, dummies_est(dx)[t0, ])
  list(m = lm(DX ~ ., data = D), fechas = etiqueta_trim(dx)[t0])
}
for (nm in names(corto)) {
  x <- X[[nm]]; A <- adf_lm(x, 4); v <- dep_x[corto[nm]]; cf <- summary(A$m)$coefficients
  nombres <- c(X_1 = sprintf("L%s(-1)", v), setNames(sprintf("D(L%s(-%d))", v, 1:4), sprintf("DX_%d", 1:4)), TREND = "@TREND")
  vc <- VC_ADF["tau_tau", ]
  previo <- list(filas = list(c("Augmented Dickey-Fuller test statistic", fmt(cf["X_1", 3], 6)), c("Test critical values:"),
                              c("      1% level", fmt(vc["1%"], 6)), c("      5% level", fmt(vc["5%"], 6)), c("     10% level", fmt(vc["10%"], 6))),
                 post = c("Critical values: Dickey-Fuller tau_tau (n ~ 100), as in urca::ur.df.", "Test Equation: the Prob. of the lagged level does not apply (non-standard distribution)."))
  cab <- c(sprintf("Null Hypothesis: L%s has a unit root", v), "Exogenous: Constant, Linear Trend, Seasonal dummies (Q1-Q3)", "Lag Length: 4 (fixed)")
  S <- salida_mco(A$m, sprintf("D(L%s)", v), sprintf("Salida %d. Test ADF: %s", n_salida + match(nm, names(corto)), corto[nm]),
                  sprintf("%s %s", A$fechas[1], A$fechas[length(A$fechas)]), "Least Squares", previo = previo, cab_previa = cab, nombres = nombres)
  guardar_salida(S, sprintf("output/salidas_eviews/08_adf_%s.png", c("PBI real" = "pbi", "Electricidad" = "elec", "Ferrocarril" = "rail")[corto[nm]]))
}
cat("Listo: output/salidas_eviews/\n"); print(list.files("output/salidas_eviews"))
