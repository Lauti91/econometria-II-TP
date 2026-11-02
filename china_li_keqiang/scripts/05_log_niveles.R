# =============================================================================
# 05_log_niveles.R -- Quinta pasada: PBI real, electricidad y carga ferroviaria en log-niveles
# Ciclo Box-Jenkins desde el nivel (Nota 6, Sec. 12, pasos 1 a 6):
#   1. Transformacion para estabilizar la varianza: relacion rango-mediana (Nota 6, Sec. 5.4)
#   2. Graficos y FAC en niveles (Nota 6, Sec. 9.3: fila "no estacionario")
#   3. Estacionalidad: deterministica (dummies) o estocastica (raices estacionales, (1-B^4)),
#      Nota 6, Sec. 2.3. El HEGY queda fuera del programa: se usa evidencia indirecta.
#   4. Raiz unitaria en z = 1: ADF con dummies (Sec. 4.4), esquema descendente de Dickey-Pantula
#      (Sec. 4.5), ET vs ED (Sec. 3.3), por subperiodos (regla 6 de la clase: quiebres y atipicos)
#   5. Filtros rivales comparados en una muestra comun por SBC (Sec. 8) y diagnostico (Sec. 10)
# Ejecutar desde china_li_keqiang/: Rscript scripts/05_log_niveles.R
# =============================================================================
source("scripts/utils.R")
dir.create("output/tablas", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figuras", showWarnings = FALSE, recursive = TRUE)

L <- cargar_niveles()
X <- L[c("ln PBI real (x100)", "ln Electricidad (x100)", "ln Carga ferroviaria (x100)")]
cat(sprintf("Control del PBI real encadenado vs publicado 2011-2015: diferencia maxima %.2f%%\n", L$control_max_dif_pct))
write.csv(L$niveles, "output/tablas/05_niveles.csv", row.names = FALSE)
colores <- c(col_pbi, col_elec, col_rail); names(colores) <- names(X)
corto <- c("PBI real", "Electricidad", "Ferrocarril"); names(corto) <- names(X)
archivo <- c("pbi", "elec", "rail"); names(archivo) <- corto
inicio <- function(x) time(x)[6]          # primer dato de (1-B)(1-B^4)x: 2001Q2 (PBI, elec.), 2006Q2 (ferrocarril)
rec <- function(x, a, b) window(x, start = a, end = b)

# ---- 1. Transformacion de varianza: rango-mediana (Nota 6, Sec. 5.4) ----------------------
# Grupos = anios (4 trimestres). Pendiente de ln(rango) sobre ln(mediana) = 1 - lambda.
bc <- list(); bc_pts <- list()
for (j in 1:3) {
  v <- L$niveles[[c("pbi_real", "electricidad_GWh", "carga_ferroviaria_10k_t")[j]]]
  a <- rep(2000:2023, each = 4); a <- a[!is.na(v)]; v <- v[!is.na(v)]
  M <- tapply(v, a, median); R <- tapply(v, a, function(z) diff(range(z)))
  m <- lm(log(R) ~ log(M)); s <- summary(m)$coefficients
  bc[[j]] <- data.frame(serie = unname(corto[j]), pendiente = s[2, 1], ee = s[2, 2],
                        lambda = 1 - s[2, 1], IC95_lambda = sprintf("[%.2f, %.2f]", 1 - s[2, 1] - 1.96 * s[2, 2], 1 - s[2, 1] + 1.96 * s[2, 2]))
  bc_pts[[j]] <- data.frame(serie = unname(corto[j]), lM = as.numeric(log(M)), lR = as.numeric(log(R)))
}
BC <- do.call(rbind, bc); write.csv(BC, "output/tablas/05_box_cox.csv", row.names = FALSE)
cat("\n=== 1. Relacion rango-mediana (lambda = 1 - pendiente; lambda = 0 => logaritmo) ===\n")
print(transform(BC, pendiente = round(pendiente, 3), ee = round(ee, 3), lambda = round(lambda, 2)), row.names = FALSE)

# ---- 2. Niveles y FAC en niveles ---------------------------------------------------------------
ev_crisis <- c(2008.5, 2009.5); ev_covid <- c(2020, 2021)
p_niv <- function(nm) {
  x <- X[[nm]]; d <- data.frame(t = as.numeric(time(x)), y = as.numeric(x))
  ggplot(d, aes(t, y)) +
    annotate("rect", xmin = ev_crisis[1], xmax = ev_crisis[2], ymin = -Inf, ymax = Inf, fill = gris, alpha = 0.45) +
    annotate("rect", xmin = ev_covid[1], xmax = ev_covid[2], ymin = -Inf, ymax = Inf, fill = gris, alpha = 0.45) +
    geom_line(colour = colores[nm], linewidth = 0.8) +
    scale_x_continuous(breaks = seq(2000, 2025, 5)) +
    labs(title = corto[nm], subtitle = "100 x ln(nivel); franjas: crisis 2008-09 y COVID 2020", x = NULL, y = NULL) + tema
}
p_fac <- function(nm, lag.max = 24) {
  x <- as.numeric(X[[nm]]); b <- 1.96 / sqrt(length(x))
  d <- data.frame(k = 1:lag.max, r = acf(x, lag.max = lag.max, plot = FALSE)$acf[-1])
  ggplot(d, aes(k, r)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -b, ymax = b, fill = gris, alpha = 0.45) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_segment(aes(xend = k, yend = 0), colour = colores[nm], linewidth = 1.2, lineend = "round") +
    scale_x_continuous(breaks = c(1, seq(4, lag.max, 4))) + scale_y_continuous(limits = c(-0.6, 1)) +
    labs(title = paste("FAC en niveles -", corto[nm]), subtitle = "decaimiento casi lineal: raiz unitaria (Nota 6, Sec. 9.3)",
         x = "Rezago (trimestres)", y = NULL) + tema
}
f2 <- (p_niv(names(X)[1]) | p_niv(names(X)[2]) | p_niv(names(X)[3])) /
      (p_fac(names(X)[1]) | p_fac(names(X)[2]) | p_fac(names(X)[3])) +
  plot_annotation(title = "PBI real, electricidad y carga ferroviaria en log-niveles (ferrocarril desde 2005)",
                  caption = "PBI real: tasa interanual oficial encadenada desde los niveles publicados de 2011. Fuente: NBS (chinadata.live), OCDE/FRED.",
                  theme = tema_global)
ggsave("output/figuras/05_niveles.png", f2, width = 13, height = 7.5, dpi = 150, bg = fondo)

# ---- 3. Estacionalidad: deterministica o estocastica? -------------------------------------
# (a) Medias estacionales de Dx por subperiodo: si la estacionalidad es deterministica, el patron
#     se mantiene; si hay raices estacionales, el patron deriva (Nota 6, Sec. 2.3).
# (b) F (indicativo) de igualdad de las dummies entre subperiodos.
# (c) Firma de sobrediferenciacion estacional: en (1-B)(1-B^4)x = (1-theta B)(1-THETA B^4)e,
#     THETA cercano a 1 indica que (1-B^4) sobraba (analogo estacional de la Nota 6, Sec. 2.5).
sub3 <- list("hasta 2007" = c(2000, 2007), "2008-2015" = c(2008, 2015), "2016-2023" = c(2016, 2023))
est <- list(); bb <- list()
for (nm in names(X)) {
  dx <- diff(X[[nm]]); tt <- as.numeric(floor(time(dx) + 1e-6)); q <- as.numeric(cycle(dx))
  per <- cut(tt, c(1999, 2007, 2015, 2023), labels = names(sub3))
  medias <- tapply(dx, list(per, q), mean)
  D <- dummies_est(dx)
  m0 <- lm(dx ~ per + D); m1 <- lm(dx ~ per * D)
  Ft <- anova(m0, m1)
  ma <- arima(X[[nm]], order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 4), method = "ML")
  TH <- -ma$coef["sma1"]; seTH <- sqrt(ma$var.coef["sma1", "sma1"])
  rr <- acf(residuals(m0), lag.max = 12, plot = FALSE)$acf[-1]
  est[[nm]] <- data.frame(
    serie = unname(corto[nm]),
    patron_hasta_07 = paste(sprintf("%.1f", medias[1, ]), collapse = " / "),
    patron_2008_15 = paste(sprintf("%.1f", medias[2, ]), collapse = " / "),
    patron_2016_23 = paste(sprintf("%.1f", medias[3, ]), collapse = " / "),
    F_cambio_patron = Ft$F[2], p_F = Ft$`Pr(>F)`[2],
    THETA_estacional = TH, ee_THETA = seTH, IC95_sup = TH + 1.96 * seTH,
    FAC_resid_dummies_r4 = rr[4], r8 = rr[8], r12 = rr[12], row.names = NULL)
  bb[[nm]] <- data.frame(serie = unname(corto[nm]), anio = tt, trimestre = paste0("T", q), dx = as.numeric(dx),
                         media = medias[cbind(as.integer(per), q)])
}
EST <- do.call(rbind, est); write.csv(EST, "output/tablas/05_estacionalidad.csv", row.names = FALSE)
cat("\n=== 3. Estacionalidad (patron = medias de Dx en T1 / T2 / T3 / T4, en %) ===\n")
print(EST[, 1:4], row.names = FALSE)
print(transform(EST[, c(1, 5:12)], F_cambio_patron = round(F_cambio_patron, 2), p_F = round(p_F, 4),
                THETA_estacional = round(THETA_estacional, 3), ee_THETA = round(ee_THETA, 3), IC95_sup = round(IC95_sup, 3),
                FAC_resid_dummies_r4 = round(FAC_resid_dummies_r4, 2), r8 = round(r8, 2), r12 = round(r12, 2)), row.names = FALSE)

BB <- do.call(rbind, bb)
g_bb <- ggplot(BB, aes(anio, dx)) +
  geom_step(aes(y = media), colour = tinta2, linewidth = 0.5, linetype = "22", direction = "mid") +
  geom_line(aes(colour = serie), linewidth = 0.6) + geom_point(aes(colour = serie), size = 0.9) +
  scale_colour_manual(values = c(`PBI real` = col_pbi, Electricidad = col_elec, Ferrocarril = col_rail), guide = "none") +
  facet_grid(serie ~ trimestre, scales = "free_y") +
  scale_x_continuous(breaks = c(2005, 2015)) +
  labs(title = "Crecimiento trimestral (Dx, %) por trimestre del anio",
       subtitle = "Linea punteada: media del trimestre hasta 2007 (ferrocarril: desde 2005), 2008-15 y 2016-23. Si el patron estacional fuera fijo, la punteada seria horizontal.",
       x = NULL, y = NULL) +
  tema + theme(strip.text = element_text(colour = tinta, face = "bold"), panel.spacing = unit(10, "pt"))
ggsave("output/figuras/05_estacionalidad.png", g_bb, width = 11, height = 8, dpi = 150, bg = fondo)

# ---- 4. Raiz unitaria en z = 1 (Nota 6, Sec. 4) ------------------------------------------
# Con p >= 3 rezagos de Dx, el ADF absorbe las raices estacionales: 1 - B^4 = (1-B)(1+B+B^2+B^3),
# asi que el factor (1+B+B^2+B^3) entra en los rezagos aumentados y tau sigue contrastando z = 1.
# Esquema descendente (Sec. 4.5): primero si Dx tiene raiz unitaria (y seria I(2)), despues x.
muestras <- list("completa (hasta 2023)" = c(2000, 2023), "sin COVID (hasta 2019)" = c(2000, 2019),
                 "inicio-2011" = c(2000, 2011), "2012-2023" = c(2012, 2023))
adf_t <- list()
for (nm in names(X)) for (mu in names(muestras)) for (p in c(3, 4, 8)) {
  w <- muestras[[mu]]; x <- rec(X[[nm]], c(max(w[1], start(X[[nm]])[1]), 1), c(w[2], 4))
  a1 <- adf_est(diff(x), p, tendencia = FALSE); a2 <- adf_est(diff(x), p, tendencia = TRUE)
  a3 <- adf_est(x, p, tendencia = TRUE)
  adf_t[[length(adf_t) + 1]] <- data.frame(
    serie = unname(corto[nm]), muestra = mu, p = p, n = a3["n"],
    Dx_tau_mu = a1["tau"], Dx_tau_tau = a2["tau"], Dx_t_tend = a2["t_tend"],
    x_tau_tau = a3["tau"], x_t_tend = a3["t_tend"], LB8_p_x = a3["LB8_p"], row.names = NULL)
}
ADF <- do.call(rbind, adf_t)
ADF$Dx_rechaza_5 <- ADF$Dx_tau_mu < VC_ADF["tau_mu", "5%"]
ADF$Dx_tend_rechaza_5 <- ADF$Dx_tau_tau < VC_ADF["tau_tau", "5%"]
ADF$x_rechaza_5 <- ADF$x_tau_tau < VC_ADF["tau_tau", "5%"]
write.csv(ADF, "output/tablas/05_adf.csv", row.names = FALSE)
cat("\n=== 4. ADF con dummies estacionales (VC 5%: tau_mu -2.89, tau_tau -3.45; con n chico son algo mas exigentes) ===\n")
print(transform(ADF[, c("serie", "muestra", "p", "n", "Dx_tau_mu", "Dx_tau_tau", "x_tau_tau", "Dx_rechaza_5", "Dx_tend_rechaza_5", "x_rechaza_5")],
                Dx_tau_mu = round(Dx_tau_mu, 2), Dx_tau_tau = round(Dx_tau_tau, 2), x_tau_tau = round(x_tau_tau, 2)), row.names = FALSE)

# ---- 5. Filtros rivales en una muestra comun ------------------------------------------------
# Todos son modelos para el mismo x_t (log-nivel): cada filtro W = F(B) x con F(0) = 1 tiene
# jacobiano unitario, de modo que la verosimilitud de W en las mismas fechas es la de x_t dado los
# valores iniciales. Por eso, con la misma muestra (desde el 6to trimestre de cada serie: 2001Q2, o
# 2006Q2 para el ferrocarril), los SBC SI son comparables
# entre filtros, y sigma es el error de pronostico a un paso del log-nivel en todos los casos.
#   F1: (1-B) x       + constante + dummies      estacionalidad deterministica, x ~ I(1)
#   F2: (1-B^4) x     + constante                estacionalidad estocastica, yoy estacionario (=03)
#   F2t:(1-B^4) x     + constante + tendencia    yoy estacionario alrededor de una tendencia (desaceleracion)
#   F3: (1-B)(1-B^4) x (+ constante)             yoy con raiz unitaria (=04)
# Intervenciones (regla 9 de la clase): atipicos aditivos (AO) en el log-nivel en fechas con evento
# conocido, 2008Q4 (crisis global) y 2020Q1 (COVID). El indicador se define sobre x y se filtra con
# el mismo F(B) que la serie, para que el modelo de x sea el mismo en todos los filtros (en yoy, un AO
# en 2020Q1 aparece como -w en 2020Q1 y +w en 2021Q1). Se informan ambas versiones, con y sin AO.
filtros <- c("F1: (1-B)x + dummies", "F2: (1-B^4)x", "F2t: (1-B^4)x + tendencia", "F3: (1-B)(1-B^4)x")
muestras_bj <- list(
  "completa (hasta 2023Q4)"                = list(fin = c(2023, 4), ao = NULL),
  "completa + AO 2008Q4 y 2020Q1"          = list(fin = c(2023, 4), ao = c(2008.75, 2020)),
  "sin COVID (hasta 2019Q4)"               = list(fin = c(2019, 4), ao = NULL),
  "sin COVID + AO 2008Q4"                  = list(fin = c(2019, 4), ao = 2008.75))
filtrar <- function(z, f) switch(substr(f, 1, 3), "F1:" = diff(z), "F2:" = diff(z, 4), "F2t" = diff(z, 4), "F3:" = diff(diff(z, 4)))
construir <- function(x, f, ao = NULL) {
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
res_f <- list(); mejores <- list()
for (nm in names(X)) for (mu in names(muestras_bj)) for (f in filtros) {
  x <- rec(X[[nm]], start(X[[nm]]), muestras_bj[[mu]]$fin)
  b <- construir(x, f, muestras_bj[[mu]]$ao); C <- rejilla(b$W, b$xreg); D <- C[C$admisible, ]
  fila <- data.frame(serie = unname(corto[nm]), muestra = mu, filtro = f, n = length(b$W),
                     candidatos = nrow(C), admisibles = nrow(D))
  if (nrow(D)) {
    s <- D[which.min(D$SBC), ]
    fila <- cbind(fila, modelo = sprintf("ARMA(%d,%d)%s", s$p, s$q, ifelse(s$Q_est == 1, "+SMA(4)", "")),
                  p = s$p, q = s$q, Qs = s$Q_est, SBC = s$SBC, sigma = sqrt(s$sigma2),
                  min_mod_AR = s$min_mod_G, min_mod_MA = s$min_mod_H, LB12_p = s$LB12, JB_p = s$JB_p)
  } else {
    s <- C[which.min(C$SBC), ]   # sin admisibles: se informa el de menor SBC y por que falla
    fila <- cbind(fila, modelo = sprintf("ninguno (mejor SBC: ARMA(%d,%d)%s, LB12 p=%.3f)", s$p, s$q,
                                         ifelse(s$Q_est == 1, "+SMA(4)", ""), s$LB12),
                  p = NA, q = NA, Qs = NA, SBC = NA, sigma = NA, min_mod_AR = NA, min_mod_MA = NA, LB12_p = NA, JB_p = NA)
  }
  res_f[[length(res_f) + 1]] <- fila
}
FIL <- do.call(rbind, res_f)
FIL$dSBC <- ave(FIL$SBC, FIL$serie, FIL$muestra, FUN = function(v) v - min(v, na.rm = TRUE))
write.csv(FIL, "output/tablas/05_filtros.csv", row.names = FALSE)
cat("\n=== 5. Mejor modelo admisible por filtro (muestra comun por serie; SBC comparable entre filtros de una misma serie) ===\n")
num <- sapply(FIL, is.numeric); FILp <- FIL; FILp[num] <- lapply(FIL[num], round, 3)
print(FILp[, c("serie", "muestra", "filtro", "admisibles", "modelo", "SBC", "dSBC", "sigma", "min_mod_AR", "min_mod_MA", "LB12_p", "JB_p")], row.names = FALSE)

# ---- 6. Modelo ganador por serie y muestra: coeficientes, raices, diagnostico, sobreajuste ---
GAN <- do.call(rbind, lapply(split(FIL, list(FIL$serie, FIL$muestra), drop = TRUE),
                             function(d) d[which.min(d$SBC), ]))
coef_t <- list(); diag_t <- list(); sob_t <- list(); paneles <- list()
for (i in seq_len(nrow(GAN))) {
  g <- GAN[i, ]; nm <- names(corto)[corto == g$serie]
  x <- rec(X[[nm]], start(X[[nm]]), muestras_bj[[g$muestra]]$fin); b <- construir(x, g$filtro, muestras_bj[[g$muestra]]$ao)
  m <- ajustar(b$W, g$p, g$q, g$Qs, b$xreg); k <- g$p + g$q + g$Qs; n <- length(b$W)
  etiqueta <- sprintf("%s | %s", g$filtro, g$modelo)
  est_ <- m$coef; se <- sqrt(diag(m$var.coef))[names(est_)]
  signo <- ifelse(grepl("^(ma|sma)", names(est_)), -1, 1)          # MA en la convencion de las Notas
  coef_t[[i]] <- data.frame(serie = g$serie, muestra = g$muestra, modelo = etiqueta, parametro = names(est_),
                            estimacion = signo * est_, ee = se, t = signo * est_ / se, row.names = NULL)
  rt <- modulos(m); e <- as.numeric(residuals(m)); es <- e / sqrt(m$sigma2)
  jb <- jarque_bera(e, k); lb2 <- Box.test(e^2, 12, "Ljung-Box"); fechas <- etiqueta_trim(b$W)
  atip <- which(abs(es) > 3)
  diag_t[[i]] <- data.frame(serie = g$serie, muestra = g$muestra, modelo = etiqueta, n = n, sigma = sqrt(m$sigma2),
    raices_AR = paste(round(sort(Mod(rt$G)), 3), collapse = " "), raices_MA_min = if (length(rt$H)) min(Mod(rt$H)) else NA,
    LB8_p = Box.test(e, 8, "Ljung-Box", fitdf = k)$p.value, LB12_p = Box.test(e, 12, "Ljung-Box", fitdf = k)$p.value,
    LB16_p = Box.test(e, 16, "Ljung-Box", fitdf = k)$p.value, JB_p = jb["p"], LB_cuadrados12_p = lb2$p.value,
    atipicos = paste(fechas[atip], collapse = " "), row.names = NULL)
  extras <- list(c(g$p + 1, g$q, g$Qs), c(g$p, g$q + 1, g$Qs), if (g$Qs == 0) c(g$p, g$q, 1))
  for (ex in extras) {
    if (is.null(ex)) next
    m2 <- ajustar(b$W, ex[1], ex[2], ex[3], b$xreg); if (is.null(m2)) next
    nuevo <- if (ex[1] > g$p) paste0("ar", ex[1]) else if (ex[2] > g$q) paste0("ma", ex[2]) else "sma1"
    if (!nuevo %in% names(m2$coef)) next
    tn <- m2$coef[nuevo] / sqrt(diag(m2$var.coef))[nuevo]
    sob_t[[length(sob_t) + 1]] <- data.frame(serie = g$serie, muestra = g$muestra, base = etiqueta, extra = nuevo,
      t_extra = as.numeric(tn), signif_5pct = abs(tn) > 1.96, dSBC = BIC(m2) - BIC(m), row.names = NULL)
  }
  col <- colores[nm]; banda <- 1.96 / sqrt(n)
  g1 <- ggplot(data.frame(t = as.numeric(time(b$W)), e = es), aes(t, e)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3, fill = gris, alpha = 0.35) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_line(colour = col, linewidth = 0.6) + geom_point(colour = col, size = 1.1) +
    scale_x_continuous(breaks = seq(2000, 2025, 5)) +
    labs(title = "Residuos estandarizados", subtitle = "franja: +-3", x = NULL, y = NULL) + tema
  ra <- acf(e, lag.max = 16, plot = FALSE)$acf[-1]
  g2 <- ggplot(data.frame(k = 1:16, r = ra, sig = abs(ra) > banda), aes(k, r)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -banda, ymax = banda, fill = gris, alpha = 0.45) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.3, lineend = "round") +
    scale_colour_manual(values = c(`TRUE` = unname(col), `FALSE` = "#9AA5B1"), guide = "none") +
    scale_x_continuous(breaks = c(1, 4, 8, 12, 16)) + scale_y_continuous(limits = c(-0.6, 1)) +
    labs(title = "FAC de los residuos", subtitle = sprintf("banda +-1.96/raiz(T); Ljung-Box(12) p=%.2f", diag_t[[i]]$LB12_p),
         x = "Rezago (trimestres)", y = NULL) + tema
  g3 <- ggplot(data.frame(teo = qnorm(ppoints(n))[order(order(es))], mue = es), aes(teo, mue)) +
    geom_abline(slope = 1, intercept = 0, colour = tinta2, linewidth = 0.4, linetype = "22") +
    geom_point(colour = col, size = 1.3) +
    labs(title = "Q-Q normal", subtitle = sprintf("Jarque-Bera p=%.3f", jb["p"]), x = "Cuantil teorico", y = NULL) + tema
  paneles[[i]] <- wrap_elements(full = (g1 | g2 | g3) +
    plot_annotation(title = sprintf("%s - %s", g$serie, etiqueta), subtitle = g$muestra, theme = tema_global))
}
COEF <- do.call(rbind, coef_t); DIAG <- do.call(rbind, diag_t); SOB <- do.call(rbind, sob_t)
write.csv(COEF, "output/tablas/05_coeficientes.csv", row.names = FALSE)
write.csv(DIAG, "output/tablas/05_diagnostico.csv", row.names = FALSE)
write.csv(SOB,  "output/tablas/05_sobreajuste.csv", row.names = FALSE)
cat("\n=== 6. Modelos ganadores: coeficientes (MA con el signo de las Notas) ===\n")
print(transform(COEF, estimacion = round(estimacion, 3), ee = round(ee, 3), t = round(t, 2)), row.names = FALSE)
cat("\n--- Diagnostico ---\n"); print(DIAG, row.names = FALSE, digits = 3)
cat("\n--- Sobreajuste (Nota 6, Sec. 10.4) ---\n"); print(SOB, row.names = FALSE, digits = 3)
for (s in unique(GAN$serie)) {
  idx <- which(GAN$serie == s)
  idx <- idx[grepl("AO", GAN$muestra[idx])]          # figura: especificacion con intervenciones
  ggsave(sprintf("output/figuras/05_residuos_%s.png", archivo[s]),
         wrap_plots(paneles[idx], ncol = 1), width = 11, height = 4.2 * length(idx), dpi = 150, bg = fondo)
}
cat("\nFiguras: output/figuras/05_niveles.png, 05_estacionalidad.png, 05_residuos_{pbi,elec,rail}.png\n")
