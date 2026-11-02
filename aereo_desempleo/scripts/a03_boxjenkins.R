# =============================================================================
# a03_boxjenkins.R -- Ciclo Box-Jenkins de cada serie, muestra principal 2000-2019 (Nota 6, Sec. 6-12)
# Filtros rivales para el mismo x en las mismas fechas (jacobiano unitario => SBC comparable):
#   F0: x + dummies (estacionaria en nivel)    F1: (1-B) x + dummies       F2: (1-B^12) x
#   F2t: (1-B^12) x + tendencia                 F3: (1-B)(1-B^12) x (sin constante)
# RPM (x = 100 ln millas-pasajero). Intervenciones con explicacion historica (regla 5):
#   - 11-S: pulsos en sep, oct, nov y dic de 2001 (la caida se va disipando mes a mes: es la forma de
#     una intervencion con decaimiento gradual, que se modela bien en la segunda parte)
#   - febrero de anio bisiesto (29 dias en vez de 28): efecto calendario esperado de +3.6%
# Desempleo (pp y 100 ln): pp vs. log se comparan sumando el jacobiano, l_pp = l_log + sum ln(100/y).
# Ejecutar desde aereo_desempleo/:  Rscript scripts/a03_boxjenkins.R
# =============================================================================
source("scripts/utils_a.R")
D <- cargar_base()
rpm <- window(D$rpm, end = c(2019, 12)); des <- window(D$des, end = c(2019, 12))
des_log <- ts(100 * log(des), start = start(des), frequency = 12)

filtrar <- function(z, f) switch(f, F0 = z, F1 = diff(z), F2 = diff(z, 12), F2t = diff(z, 12), F3 = diff(diff(z, 12)))
desc_filtro <- c(F0 = "x + dummies", F1 = "(1-B)x + dummies", F2 = "(1-B^12)x", F2t = "(1-B^12)x + tendencia", F3 = "(1-B)(1-B^12)x")
pulso <- function(x, a, b) ts(as.numeric(abs(time(x) - fecha_num(a, b)) < 1e-6), start = start(x), frequency = 12)
intervenciones <- function(x, nombres) {
  L <- list(AO2001M09 = pulso(x, 2001, 9), AO2001M10 = pulso(x, 2001, 10), AO2001M11 = pulso(x, 2001, 11), AO2001M12 = pulso(x, 2001, 12),
            BISIESTO = bisiesto(x))
  L[nombres]
}
construir <- function(x, f, ind, desde) {
  W <- window(filtrar(x, f), start = desde)
  xr <- switch(f, F0 = dummies_est(W), F1 = dummies_est(W), F2t = cbind(tend = seq_along(W)), NULL)
  for (nm in names(ind)) {
    z <- matrix(as.numeric(window(filtrar(ind[[nm]], f), start = desde)), ncol = 1, dimnames = list(NULL, nm))
    xr <- if (is.null(xr)) z else cbind(as.matrix(xr), z)
  }
  list(W = W, xreg = if (is.null(xr)) NULL else as.matrix(xr), media = f != "F3")
}
INT_RPM <- c("AO2001M09", "AO2001M10", "AO2001M11", "AO2001M12", "BISIESTO")
FILTROS <- c("F0", "F1", "F2", "F2t", "F3")
specs <- list(
  list(serie = "RPM", trans = "100 ln", x = rpm, muestra = "2000-2019", ind = character(0)),
  list(serie = "RPM", trans = "100 ln", x = rpm, muestra = "2000-2019 + 11-S + bisiesto", ind = INT_RPM),
  list(serie = "Desempleo", trans = "pp", x = des, muestra = "2000-2019", ind = character(0)),
  list(serie = "Desempleo", trans = "100 ln", x = des_log, muestra = "2000-2019", ind = character(0)))

res <- list(); modelos <- list()
for (sp in specs) {
  x <- sp$x; desde <- time(x)[14]
  jac <- if (sp$trans == "100 ln") sum(log(100 / exp(as.numeric(window(x, start = desde)) / 100))) else 0
  for (f in FILTROS) {
    b <- construir(x, f, intervenciones(x, sp$ind), desde)
    C <- rejilla(b$W, b$xreg, media = b$media); A <- C[C$admisible, ]
    s <- if (nrow(A)) A[which.min(A$SBC), ] else C[which.min(C$SBC), ]
    clave <- paste(sp$serie, sp$trans, sp$muestra, f, sep = " | ")
    modelos[[clave]] <- list(sp = sp, f = f, b = b, s = s)
    res[[clave]] <- data.frame(serie = sp$serie, trans = sp$trans, muestra = sp$muestra, filtro = f, desc = desc_filtro[f],
      n = length(b$W), admisibles = nrow(A), modelo = etiqueta_modelo(s$p, s$q, s$Q_est), p = s$p, q = s$q, Qs = s$Q_est,
      admisible = s$admisible, SBC = s$SBC, SBC_escala_pp = s$SBC - 2 * jac, sigma = sqrt(s$sigma2),
      min_mod_AR = s$min_mod_G, min_mod_MA = s$min_mod_H, LB12 = s$LB12, LB24 = s$LB24, LB36 = s$LB36, JB_p = s$JB_p, row.names = NULL)
  }
}
FIL <- do.call(rbind, res)
FIL$dSBC <- ave(FIL$SBC_escala_pp, FIL$serie, FIL$muestra, FUN = function(v) v - min(v))
write.csv(FIL, "output/tablas/a03_filtros.csv", row.names = FALSE)
cat("\n=== Mejor modelo por filtro (admisible si lo hay); dSBC dentro de serie y muestra (pp y log comparables) ===\n")
num <- sapply(FIL, is.numeric); P <- FIL; P[num] <- lapply(FIL[num], round, 3)
print(P[, c("serie", "trans", "muestra", "filtro", "admisibles", "modelo", "admisible", "SBC_escala_pp", "dSBC", "sigma", "min_mod_AR", "LB12", "LB24", "LB36", "JB_p")], row.names = FALSE)

# ---- Ganadores (admisible primero) -----------------------------------------------------------------------
GAN <- do.call(rbind, lapply(split(FIL, list(FIL$serie, FIL$muestra), drop = TRUE), function(d) {
  a <- d[d$admisible, ]; if (nrow(a)) a[which.min(a$SBC_escala_pp), ] else d[which.min(d$SBC_escala_pp), ] }))
GAN <- GAN[order(GAN$serie, GAN$muestra), ]
coef_t <- list(); diag_t <- list(); sob_t <- list(); ganadores <- list()
for (i in seq_len(nrow(GAN))) {
  g <- GAN[i, ]; b <- modelos[[paste(g$serie, g$trans, g$muestra, g$filtro, sep = " | ")]]$b
  m <- ajustar(b$W, g$p, g$q, g$Qs, b$xreg, b$media)
  ganadores[[paste(g$serie, g$muestra, sep = " | ")]] <- list(m = m, g = g, b = b)
  est <- m$coef; se <- sqrt(diag(m$var.coef))[names(est)]; signo <- ifelse(grepl("^(ma|sma)", names(est)), -1, 1)
  coef_t[[i]] <- data.frame(serie = g$serie, muestra = g$muestra, modelo = sprintf("%s %s: %s", g$trans, g$desc, g$modelo),
                            parametro = names(est), estimacion = signo * est, ee = se, t = signo * est / se, row.names = NULL)
  e <- as.numeric(residuals(m)); es <- e / sqrt(m$sigma2); rt <- modulos(m)
  diag_t[[i]] <- data.frame(serie = g$serie, muestra = g$muestra, modelo = sprintf("%s %s: %s", g$trans, g$desc, g$modelo), n = length(e),
    sigma = sqrt(m$sigma2), raices_AR = paste(round(sort(Mod(rt$G)), 3), collapse = " "), LB12 = g$LB12, LB24 = g$LB24, LB36 = g$LB36,
    JB_p = g$JB_p, LB_cuadrados24_p = Box.test(e^2, 24, "Ljung-Box")$p.value,
    atipicos_3.5 = paste(etiqueta_mes(b$W)[abs(es) > 3.5], collapse = " "), row.names = NULL)
  for (ex in list(c(g$p + 1, g$q, g$Qs), c(g$p, g$q + 1, g$Qs), if (g$Qs == 0) c(g$p, g$q, 1))) {
    if (is.null(ex)) next
    m2 <- ajustar(b$W, ex[1], ex[2], ex[3], b$xreg, b$media); if (is.null(m2) || m2$code != 0) next
    nuevo <- if (ex[1] > g$p) paste0("ar", ex[1]) else if (ex[2] > g$q) paste0("ma", ex[2]) else "sma1"
    tn <- m2$coef[nuevo] / sqrt(diag(m2$var.coef))[nuevo]
    sob_t[[length(sob_t) + 1]] <- data.frame(serie = g$serie, muestra = g$muestra, base = g$modelo, extra = nuevo,
      t_extra = as.numeric(tn), signif_5 = abs(tn) > 1.96, dSBC = BIC(m2) - BIC(m), row.names = NULL)
  }
}
COEF <- do.call(rbind, coef_t); DIAG <- do.call(rbind, diag_t); SOB <- do.call(rbind, sob_t)
write.csv(COEF, "output/tablas/a03_coeficientes.csv", row.names = FALSE)
write.csv(DIAG, "output/tablas/a03_diagnostico.csv", row.names = FALSE)
write.csv(SOB, "output/tablas/a03_sobreajuste.csv", row.names = FALSE)
cat("\n=== Ganadores: coeficientes (MA con el signo de las Notas) ===\n")
COEFp <- COEF; COEFp[5:7] <- lapply(COEF[5:7], round, 3); print(COEFp, row.names = FALSE)
cat("\n=== Diagnostico ===\n"); DIAGp <- DIAG; num <- sapply(DIAG, is.numeric); DIAGp[num] <- lapply(DIAG[num], round, 3); print(DIAGp, row.names = FALSE)
cat("\n=== Sobreajuste ===\n"); SOBp <- SOB; SOBp[c(5, 7)] <- lapply(SOB[c(5, 7)], round, 2); print(SOBp, row.names = FALSE)
saveRDS(ganadores, "output/tablas/a03_ganadores.rds")

# ---- Figura de residuos ------------------------------------------------------------------------------------
panel_res <- function(G, col, titulo) {
  m <- G$m; e <- as.numeric(residuals(m)); es <- e / sqrt(m$sigma2); n <- length(e); banda <- 1.96 / sqrt(n); tt <- as.numeric(time(G$b$W))
  g1 <- ggplot(data.frame(t = tt, e = es), aes(t, e)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3, fill = gris, alpha = 0.35) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_line(colour = col, linewidth = 0.45) +
    scale_x_continuous(breaks = seq(2000, 2020, 5)) + labs(title = "Residuos estandarizados", subtitle = "franja: +-3", x = NULL, y = NULL) + tema
  ra <- acf(e, 36, plot = FALSE)$acf[-1]
  g2 <- ggplot(data.frame(k = 1:36, r = ra, sig = abs(ra) > banda), aes(k, r)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -banda, ymax = banda, fill = gris, alpha = 0.45) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.1, lineend = "round") +
    scale_colour_manual(values = c(`TRUE` = col, `FALSE` = "#9AA5B1"), guide = "none") + scale_x_continuous(breaks = c(1, 12, 24, 36)) +
    scale_y_continuous(limits = c(-0.35, 0.35)) +
    labs(title = "FAC de los residuos", subtitle = sprintf("Ljung-Box p: 12 -> %.2f; 24 -> %.2f; 36 -> %.2f", G$g$LB12, G$g$LB24, G$g$LB36), x = "Rezago (meses)", y = NULL) + tema
  g3 <- ggplot(data.frame(teo = qnorm(ppoints(n))[order(order(es))], mue = es), aes(teo, mue)) +
    geom_abline(slope = 1, intercept = 0, colour = tinta2, linewidth = 0.4, linetype = "22") + geom_point(colour = col, size = 0.9) +
    labs(title = "Q-Q normal", subtitle = sprintf("Jarque-Bera p = %.2g", G$g$JB_p), x = "Cuantil teorico", y = NULL) + tema
  wrap_elements(full = (g1 | g2 | g3) + plot_annotation(title = titulo, theme = tema_global))
}
G1 <- ganadores[["RPM | 2000-2019 + 11-S + bisiesto"]]; G2 <- ganadores[["Desempleo | 2000-2019"]]
fig <- panel_res(G1, col_rpm, sprintf("RPM, 2000-2019 + 11-S + bisiesto: %s %s, %s", G1$g$trans, G1$g$desc, G1$g$modelo)) /
       panel_res(G2, col_des, sprintf("Desempleo, 2000-2019: %s %s, %s", G2$g$trans, G2$g$desc, G2$g$modelo))
ggsave("output/figuras/a03_residuos.png", fig, width = 13, height = 8, dpi = 200)
