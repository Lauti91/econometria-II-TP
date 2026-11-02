# =============================================================================
# a06_salidas_eviews.R -- Salidas con formato EViews del par 4 (usa ../scripts/eviews.R)
# Variables: LRPM = 100 ln RPM; UNR = desempleo NSA (pp); LUNR = 100 ln UNR.
# El MA aparece con el signo de EViews (+theta); en las Notas es el opuesto.
# =============================================================================
source("scripts/utils_a.R")
D <- cargar_base(); G <- readRDS("output/tablas/a03_ganadores.rds")
dir <- "output/salidas_eviews"; n_sal <- 0
em <- function(z) sub("-", "M", etiqueta_mes(z))
nombre_dep <- function(serie, trans, filtro) {
  v <- if (serie == "RPM") "LRPM" else if (trans == "100 ln") "LUNR" else "UNR"
  switch(filtro, F0 = v, F1 = sprintf("D(%s)", v), F2 = sprintf("D(%s,0,12)", v), F2t = sprintf("D(%s,0,12)", v), F3 = sprintf("D(%s,1,12)", v))
}
for (nm in names(G)) {
  g <- G[[nm]]$g; n_sal <- n_sal + 1
  S <- salida_arma(G[[nm]]$m, G[[nm]]$b$W, nombre_dep(g$serie, g$trans, g$filtro), sprintf("Salida %d. %s, %s", n_sal, g$serie, g$muestra))
  guardar_salida(S, sprintf("%s/a06_%02d_%s.png", dir, n_sal, gsub("[^0-9A-Za-z]+", "_", tolower(paste(g$serie, substr(g$muestra, 1, 14))))))
}

# ---- ADF --------------------------------------------------------------------------------------------------
adf_lm <- function(x, p, tendencia, dummies) {
  dx <- diff(x); t0 <- (p + 1):length(dx); Z <- data.frame(DX = as.numeric(dx)[t0], X_1 = as.numeric(x)[t0])
  for (j in seq_len(p)) Z[[sprintf("DX_%d", j)]] <- as.numeric(dx)[t0 - j]
  if (tendencia) Z$TREND <- t0
  if (dummies) Z <- cbind(Z, dummies_est(dx)[t0, ])
  list(m = lm(DX ~ ., data = Z), fechas = em(dx)[t0])
}
adf_salida <- function(x, v, p, tendencia, dummies, muestra_txt) {
  A <- adf_lm(x, p, tendencia, dummies); cf <- summary(A$m)$coefficients
  nombres <- c(X_1 = sprintf("%s(-1)", v), setNames(sprintf("D(%s(-%d))", v, 1:p), sprintf("DX_%d", 1:p)), TREND = "@TREND")
  vc <- VC_ADF[if (tendencia) "tau_tau" else "tau_mu", ]
  previo <- list(filas = list(c("Augmented Dickey-Fuller test statistic", fmt(cf["X_1", 3])), c("Test critical values:"),
                              c("      1% level", fmt(vc["1%"])), c("      5% level", fmt(vc["5%"])), c("     10% level", fmt(vc["10%"]))),
                 post = c(sprintf("Critical values: Dickey-Fuller %s, n = 500 (Fuller 1976); with n ~ 230 they are slightly more negative.", if (tendencia) "tau_tau" else "tau_mu"),
                          "Test Equation: the Prob. of the lagged level does not apply (non-standard distribution)."))
  exo <- paste(c("Constant", if (tendencia) "Linear Trend", if (dummies) "Seasonal dummies (Jan-Nov)"), collapse = ", ")
  cab <- c(sprintf("Null Hypothesis: %s has a unit root", v), paste("Exogenous:", exo), sprintf("Lag Length: %d (fixed)", p))
  n_sal <<- n_sal + 1
  salida_mco(A$m, sprintf("D(%s)", v), sprintf("Salida %d. Test ADF: %s, %s", n_sal, v, muestra_txt),
             sprintf("%s %s", A$fechas[1], A$fechas[length(A$fechas)]), "Least Squares", previo = previo, cab_previa = cab, nombres = nombres)
}
r19 <- window(D$rpm, end = c(2019, 12)); u19 <- window(D$des, end = c(2019, 12))
guardar_salida(adf_salida(r19, "LRPM", 12, TRUE, TRUE, "2000-2019"), sprintf("%s/a06_adf_lrpm_2000_2019.png", dir))
guardar_salida(adf_salida(D$rpm, "LRPM", 12, TRUE, TRUE, "2000-2026 (con COVID)"), sprintf("%s/a06_adf_lrpm_con_covid.png", dir))
guardar_salida(adf_salida(u19, "UNR", 12, TRUE, TRUE, "2000-2019"), sprintf("%s/a06_adf_unr.png", dir))
guardar_salida(adf_salida(diff(u19), "DUNR", 12, FALSE, TRUE, "2000-2019 (primera diferencia)"), sprintf("%s/a06_adf_dunr.png", dir))

# ---- Regresiones de la relacion (a04) -------------------------------------------------------------------------
pulso <- function(x, a, b) ts(as.numeric(abs(time(x) - fecha_num(a, b)) < 1e-6), start = start(x), frequency = 12)
IND <- cbind(AO2001M09 = pulso(r19, 2001, 9), AO2001M10 = pulso(r19, 2001, 10), AO2001M11 = pulso(r19, 2001, 11),
             AO2001M12 = pulso(r19, 2001, 12), BISIESTO = bisiesto(r19))
Zn <- data.frame(LRPM = as.numeric(r19), UNR = as.numeric(u19), dummies_est(r19), as.data.frame(IND))
n_sal <- n_sal + 1
guardar_salida(salida_mco(lm(LRPM ~ ., data = Zn), "LRPM", sprintf("Salida %d. Regresion en niveles con meses (alarma: R2 > DW)", n_sal),
                          sprintf("%s %s", em(r19)[1], em(r19)[length(r19)])), sprintf("%s/a06_reg_niveles.png", dir))
F <- function(z) diff(diff(z, 12))
Z <- ts.intersect(F(r19), do.call(cbind, lapply(0:6, function(k) stats::lag(F(u19), -k))), F(IND))
colnames(Z) <- c("Y", sprintf("DUNR_%d", 0:6), colnames(IND))
nombres <- c(setNames(c("D(UNR,1,12)", sprintf("D(UNR(-%d),1,12)", 1:6)), sprintf("DUNR_%d", 0:6)),
             setNames(sprintf("D(%s,1,12)", colnames(IND)), colnames(IND)))
n_sal <- n_sal + 1
guardar_salida(salida_mco(lm(Y ~ ., data = as.data.frame(Z)), "D(LRPM,1,12)", sprintf("Salida %d. Rezagos distribuidos 0-6, filtro (1-B)(1-B^12)", n_sal),
                          sprintf("%s %s", em(Z)[1], em(Z)[nrow(Z)]), nombres = nombres), sprintf("%s/a06_reg_rezagos.png", dir))
cat("Salidas:\n"); print(list.files(dir))
