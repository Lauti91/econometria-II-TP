# =============================================================================
# a05_pronostico.R -- Evaluacion fuera de muestra (Nota 5, Sec. 9; Nota 6, Sec. 11)
# Estimacion 2000-01 a 2016-12; evaluacion 2017-01 a 2019-12 (antes del COVID), origen movil,
# horizontes h = 1, 3, 6, 12. Orden elegido solo con la muestra de estimacion; parametros reestimados
# en cada origen. Los filtros se escriben como modelos para x, con los regresores sin filtrar (arima
# los diferencia con el mismo operador):
#   F0: (p,0,q)(0,0,Q) + constante + dummies          F1: (p,1,q) + deriva t + dummies acumuladas
#   F2: (p,0,q)(0,1,Q) + t/12                         F2t: idem + (t^2/24 + t/2)      F3: (p,1,q)(0,1,Q)
# Intervenciones en todos: pulsos sep-dic 2001 y febrero bisiesto.
# Referencias: ingenuo estacional x_{T+h-12} y paseo aleatorio estacional x_T + x_{T+h-12} - x_{T-12}.
# Desempleo: el ganador de a03 (100 ln, (1-B)(1-B^12)) con orden elegido en 2000-2016, contra las
# mismas referencias; errores en pp.
# =============================================================================
source("scripts/utils_a.R")
D <- cargar_base()
x_all <- window(D$rpm, end = c(2019, 12)); d_all <- window(100 * log(D$des), end = c(2019, 12))
FIN_EST <- fecha_num(2016, 12); FIN_EVAL <- fecha_num(2019, 12); H <- c(1, 3, 6, 12)
FILTROS <- c("F0", "F1", "F2", "F2t", "F3")

# Regresores sobre una grilla que llega 12 meses mas alla de la evaluacion
grilla <- ts(numeric(length(x_all) + 12), start = start(x_all), frequency = 12)
tg <- as.numeric(time(grilla)); t_idx <- seq_along(tg); mes <- cycle(grilla)
pul <- sapply(9:12, function(m) as.numeric(abs(tg - fecha_num(2001, m)) < 1e-6)); colnames(pul) <- sprintf("AO2001M%02d", 9:12)
INT <- cbind(pul, BISIESTO = as.numeric(bisiesto(grilla)))
Dm <- sapply(1:11, function(j) (mes == j) - (mes == 12)); colnames(Dm) <- month.abb[1:11]
REG <- list(F0 = cbind(Dm, INT), F1 = cbind(deriva = t_idx, apply(Dm, 2, cumsum), INT), F2 = cbind(c12 = t_idx / 12, INT),
            F2t = cbind(c12 = t_idx / 12, tend12 = t_idx^2 / 24 + t_idx / 2, INT), F3 = INT)
REG <- lapply(REG, function(r) ts(r, start = start(grilla), frequency = 12))
modelo_x <- function(x, f, p, q, Qs, reg = REG) {
  xr <- window(reg[[f]], start = start(x), end = end(x)); xr <- xr[, colSums(abs(xr)) > 0, drop = FALSE]
  d <- if (f %in% c("F1", "F3")) 1 else 0; Dd <- if (f %in% c("F2", "F2t", "F3")) 1 else 0
  m <- tryCatch(suppressWarnings(arima(x, order = c(p, d, q), seasonal = list(order = c(0, Dd, Qs), period = 12), xreg = xr,
                                       include.mean = f == "F0", method = "ML", optim.control = list(maxit = 500))), error = function(e) NULL)
  if (!is.null(m)) attr(m, "cols") <- colnames(xr)
  m
}
pronosticar <- function(m, x, f, h, reg = REG) {
  xr <- window(reg[[f]], start = time(x)[length(x)] + 1 / 12, end = time(x)[length(x)] + h / 12)[, attr(m, "cols"), drop = FALSE]
  predict(m, n.ahead = h, newxreg = xr)
}

# ---- Orden con la muestra de estimacion (rejilla sobre el filtro, como en a03) ----------------------------
filtrar <- function(z, f) switch(f, F0 = z, F1 = diff(z), F2 = diff(z, 12), F2t = diff(z, 12), F3 = diff(diff(z, 12)))
x_est <- window(x_all, end = FIN_EST); desde <- time(x_est)[14]
ORD <- do.call(rbind, lapply(FILTROS, function(f) {
  W <- window(filtrar(x_est, f), start = desde)
  I <- window(ts(INT, start = start(grilla), frequency = 12), end = FIN_EST)
  Ii <- do.call(cbind, lapply(seq_len(ncol(I)), function(j) as.numeric(window(filtrar(I[, j], f), start = desde)))); colnames(Ii) <- colnames(I)
  xr <- cbind(switch(f, F0 = dummies_est(W), F1 = dummies_est(W), F2t = cbind(tend = seq_along(W)), NULL), Ii)
  C <- rejilla(W, xr, media = f != "F3"); A <- C[C$admisible, ]; s <- if (nrow(A)) A[which.min(A$SBC), ] else C[which.min(C$SBC), ]
  data.frame(filtro = f, p = s$p, q = s$q, Qs = s$Q_est, SBC = s$SBC, admisible = s$admisible)
}))
ORD$dSBC <- ORD$SBC - min(ORD$SBC); ORD <- ORD[order(ORD$SBC), ]
cat("\n=== RPM: orden elegido con 2000-2016 ===\n"); print(transform(ORD, SBC = round(SBC, 2), dSBC = round(dSBC, 2)), row.names = FALSE)
mejores2 <- ORD$filtro[1:2]
dW <- window(diff(diff(window(d_all, end = FIN_EST), 12)), start = desde)
Cd <- rejilla(dW, media = FALSE); Ad <- Cd[Cd$admisible, ]; od <- if (nrow(Ad)) Ad[which.min(Ad$SBC), ] else Cd[which.min(Cd$SBC), ]
cat(sprintf("\nDesempleo: orden elegido con 2000-2016 en (1-B)(1-B^12) 100 ln: %s (admisible: %s)\n", etiqueta_modelo(od$p, od$q, od$Q_est), od$admisible))
ORD <- rbind(ORD, data.frame(filtro = "Desempleo F3", p = od$p, q = od$q, Qs = od$Q_est, SBC = od$SBC, admisible = od$admisible, dSBC = NA))
write.csv(ORD, "output/tablas/a05_orden_estimacion.csv", row.names = FALSE)

# ---- Origen movil -------------------------------------------------------------------------------------------
origenes <- as.numeric(time(x_all))[as.numeric(time(x_all)) >= FIN_EST - 1e-6 & as.numeric(time(x_all)) < FIN_EVAL - 1e-6]
filas <- list()
ingenuos <- function(x, hmax) { n <- length(x); list(`Ingenuo estacional` = as.numeric(x)[n - 12 + 1:hmax],
  `Paseo aleatorio estacional` = as.numeric(x)[n] + as.numeric(x)[n - 12 + 1:hmax] - as.numeric(x)[n - 12]) }
for (T0 in origenes) {
  x <- window(x_all, end = T0); hmax <- min(12, round((FIN_EVAL - T0) * 12)); etq <- etiqueta_mes(window(x_all, start = T0, end = T0))
  real <- as.numeric(window(x_all, start = T0 + 1 / 12, end = T0 + hmax / 12)); pr <- list()
  for (f in FILTROS) {
    o <- ORD[ORD$filtro == f, ]; m <- modelo_x(x, f, o$p, o$q, o$Qs); if (is.null(m)) next
    P <- pronosticar(m, x, f, hmax); pr[[f]] <- list(m = as.numeric(P$pred), s = as.numeric(P$se))
  }
  pr[["Comb. 2 mejores"]] <- list(m = (pr[[mejores2[1]]]$m + pr[[mejores2[2]]]$m) / 2, s = rep(NA, hmax))
  for (nm in names(ingenuos(x, hmax))) pr[[nm]] <- list(m = ingenuos(x, hmax)[[nm]], s = rep(NA, hmax))
  for (nm in names(pr)) filas[[length(filas) + 1]] <- data.frame(serie = "RPM", origen = etq, modelo = nm, h = 1:hmax, real = real, pron = pr[[nm]]$m, se = pr[[nm]]$s)
  # Desempleo (en pp)
  d <- window(d_all, end = T0); real_d <- exp(as.numeric(window(d_all, start = T0 + 1 / 12, end = T0 + hmax / 12)) / 100)
  md <- tryCatch(arima(d, order = c(od$p, 1, od$q), seasonal = list(order = c(0, 1, od$Q_est), period = 12), method = "ML"), error = function(e) NULL)
  if (!is.null(md)) filas[[length(filas) + 1]] <- data.frame(serie = "Desempleo", origen = etq, modelo = "F3 (ganador a03)", h = 1:hmax,
      real = real_d, pron = exp(as.numeric(predict(md, n.ahead = hmax)$pred) / 100), se = NA)
  for (nm in names(ingenuos(d, hmax))) filas[[length(filas) + 1]] <- data.frame(serie = "Desempleo", origen = etq, modelo = nm, h = 1:hmax,
      real = real_d, pron = exp(ingenuos(d, hmax)[[nm]] / 100), se = NA)
}
PR <- do.call(rbind, filas); PR$error <- PR$real - PR$pron; PR$cubre <- ifelse(is.na(PR$se), NA, abs(PR$error) <= 1.96 * PR$se)
write.csv(PR, "output/tablas/a05_pronosticos.csv", row.names = FALSE)
sel <- PR$h %in% H
met <- do.call(rbind, lapply(split(PR[sel, ], list(PR$serie[sel], PR$modelo[sel], PR$h[sel]), drop = TRUE), function(d)
  data.frame(serie = d$serie[1], modelo = d$modelo[1], h = d$h[1], n = nrow(d), RMSE = sqrt(mean(d$error^2)), MAE = mean(abs(d$error)),
             cobertura95 = if (all(is.na(d$cubre))) NA else mean(d$cubre))))
met$RMSE_rel_PAE <- met$RMSE / met$RMSE[match(paste(met$serie, "Paseo aleatorio estacional", met$h), paste(met$serie, met$modelo, met$h))]
met <- met[order(met$serie, met$h, met$RMSE), ]
write.csv(met, "output/tablas/a05_metricas.csv", row.names = FALSE)
cat("\n=== Metricas 2017-2019 (RPM: % de 100 ln; desempleo: pp). RMSE_rel_PAE < 1: mejor que el paseo aleatorio estacional ===\n")
num <- sapply(met, is.numeric); P <- met; P[num] <- lapply(met[num], round, 3); print(P, row.names = FALSE)

# ---- Figura ---------------------------------------------------------------------------------------------------
Mi <- met[met$serie == "RPM", ]
g1 <- ggplot(Mi, aes(factor(h), RMSE, fill = modelo)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  scale_fill_manual(values = c(F0 = "#E4E7EB", F1 = "#C5CCD3", F2 = "#9AA5B1", F2t = "#7B8794", F3 = col_rpm, `Comb. 2 mejores` = "#14946A",
                               `Ingenuo estacional` = "#F2C29B", `Paseo aleatorio estacional` = col_des), name = NULL) +
  labs(title = "RPM: RMSE fuera de muestra (2017-2019) por horizonte", subtitle = "error en % (100 ln); menor es mejor", x = "Horizonte (meses)", y = NULL) +
  tema + theme(legend.position = "right")
f1 <- mejores2[1]; o <- ORD[ORD$filtro == f1, ]; m <- modelo_x(x_est, f1, o$p, o$q, o$Qs); P <- pronosticar(m, x_est, f1, 36)
tt <- as.numeric(time(P$pred)); obs <- window(x_all, start = 2014)
g2 <- ggplot() + annotate("rect", xmin = FIN_EST, xmax = FIN_EVAL + 1 / 12, ymin = -Inf, ymax = Inf, fill = "#F5F7FA") +
  geom_ribbon(aes(x = tt, ymin = exp((P$pred - 1.96 * P$se) / 100) / 1e6, ymax = exp((P$pred + 1.96 * P$se) / 100) / 1e6), fill = col_rpm, alpha = 0.15) +
  geom_line(aes(as.numeric(time(obs)), exp(as.numeric(obs) / 100) / 1e6), colour = tinta, linewidth = 0.6) +
  geom_line(aes(tt, exp(as.numeric(P$pred) / 100) / 1e6), colour = col_rpm, linewidth = 0.8, linetype = "22") +
  labs(title = sprintf("Pronostico desde 2016-12 con %s (%s), 36 meses", f1, etiqueta_modelo(o$p, o$q, o$Qs)),
       subtitle = "negro: observado (miles de millones de millas); azul: pronostico y banda del 95%", x = NULL, y = NULL) + tema
ggsave("output/figuras/a05_pronostico.png", (g1 / g2) + plot_annotation(theme = tema_global), width = 11, height = 9, dpi = 200)
