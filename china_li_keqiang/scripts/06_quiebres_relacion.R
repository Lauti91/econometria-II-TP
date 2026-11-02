# =============================================================================
# 06_quiebres_relacion.R -- Cambios estructurales y relacion del PBI oficial con dos indicadores
# fisicos: electricidad (desde 2000) y carga ferroviaria (desde 2005).
# Parte de los modelos elegidos en 05 (muestra completa con intervenciones AO 2008Q4 y AO 2020Q1):
#   PBI real y electricidad: (1-B)(1-B^4) x = (1 - THETA B^4) e
#   Ferrocarril:             (1-B) x = constante + dummies estacionales + e   (paseo al azar con deriva)
# Cada modelo es un modelo para el log-nivel x_t, asi que sus residuos son, en los tres casos, el
# error de pronostico a un paso de 100 ln(nivel): las "sorpresas" son comparables entre series
# aunque el filtro no sea el mismo (ver la nota de jacobiano unitario en 05).
#   A. Quiebre en la varianza de las innovaciones (no estacionariedad de 2do orden, Nota 6 Sec. 5;
#      LB sobre residuos al cuadrado, Nota 3 Sec. 6.3) y barrido de la fecha (regla 9: subperiodos).
#   B. Relacion: regresion estatica en el filtro comun (1-B)(1-B^4) con la alarma R2 > d_W
#      (Nota 6, Sec. 3.3) y correlacion contemporanea de las innovaciones, por subperiodo.
#      La FCC con rezagos y la funcion de transferencia (preblanqueo) son de la 2da parte.
# Ejecutar desde china_li_keqiang/: Rscript scripts/06_quiebres_relacion.R
# =============================================================================
source("scripts/utils.R")
dir.create("output/tablas", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figuras", showWarnings = FALSE, recursive = TRUE)

L <- cargar_niveles()
X <- list("PBI real" = L[["ln PBI real (x100)"]], "Electricidad" = L[["ln Electricidad (x100)"]],
          "Ferrocarril" = L[["ln Carga ferroviaria (x100)"]])
colores <- c("PBI real" = col_pbi, "Electricidad" = col_elec, "Ferrocarril" = col_rail)
indicadores <- c("Electricidad", "Ferrocarril")
AO <- function(x, a) ts(as.numeric(abs(time(x) - a) < 1e-6), start = start(x), frequency = 4)

# ---- Modelos univariados de 05 e innovaciones ----------------------------------------------
modelo <- function(s) {
  x <- X[[s]]; ini <- time(x)[6]
  if (s == "Ferrocarril") {
    W <- window(diff(x), start = ini)
    xr <- cbind(dummies_est(W), AO2008Q4 = window(diff(AO(x, 2008.75)), start = ini),
                AO2020Q1 = window(diff(AO(x, 2020)), start = ini))
    m <- arima(W, order = c(0, 0, 0), xreg = xr, method = "ML")
  } else {
    W <- window(diff(diff(x, 4)), start = ini)
    xr <- cbind(AO2008Q4 = window(diff(diff(AO(x, 2008.75), 4)), start = ini),
                AO2020Q1 = window(diff(diff(AO(x, 2020), 4)), start = ini))
    m <- arima(W, order = c(0, 0, 0), seasonal = list(order = c(0, 0, 1), period = 4), xreg = xr, method = "ML")
  }
  list(m = m, e = ts(as.numeric(residuals(m)), start = ini, frequency = 4))
}
MOD <- lapply(setNames(names(X), names(X)), modelo)
E <- lapply(MOD, `[[`, "e")
cat("Desvio de las innovaciones, muestra completa (pp):", sprintf("%s %.2f", names(E), sapply(E, sd)), "\n")

# Subperiodos. Las innovaciones se recortan con window(); para pares de series se usa la interseccion.
periodos <- list("inicio-2012Q4" = c(2000, 2012.75), "2006Q2-2012Q4" = c(2006.25, 2012.75),
                 "2013Q1-2019Q4" = c(2013, 2019.75), "2020Q1-2023Q4" = c(2020, 2023.75),
                 "sin COVID (hasta 2019Q4)" = c(2000, 2019.75), "completa (hasta 2023Q4)" = c(2000, 2023.75))
corte <- function(e, per) { a <- max(per[1], time(e)[1]); window(e, start = a, end = per[2]) }

# ---- A. Varianza de las innovaciones por subperiodo ----------------------------------------
va <- list()
for (s in names(X)) for (pn in c("inicio-2012Q4", "2013Q1-2019Q4", "2020Q1-2023Q4")) {
  e <- corte(E[[s]], periodos[[pn]])
  va[[length(va) + 1]] <- data.frame(serie = s, periodo = pn, desde = etiqueta_trim(e)[1], n = length(e),
                                     desvio_innov = sd(e), LB_cuadrados8_p = Box.test(e^2, 8, "Ljung-Box")$p.value)
}
VA <- do.call(rbind, va)
FV <- do.call(rbind, lapply(names(X), function(s) {
  e1 <- corte(E[[s]], periodos[["inicio-2012Q4"]]); e2 <- corte(E[[s]], periodos[["2013Q1-2019Q4"]])
  f <- var(e1) / var(e2)          # supone normalidad e independencia: indicativo
  data.frame(serie = s, razon_varianzas = f, gl1 = length(e1) - 1, gl2 = length(e2) - 1,
             p_valor = pf(f, length(e1) - 1, length(e2) - 1, lower.tail = FALSE))
}))
write.csv(VA, "output/tablas/06_varianza_subperiodos.csv", row.names = FALSE)
write.csv(FV, "output/tablas/06_test_varianza.csv", row.names = FALSE)
cat("\n=== A. Desvio de las innovaciones (pp de log-nivel) por subperiodo ===\n")
print(transform(VA, desvio_innov = round(desvio_innov, 2), LB_cuadrados8_p = round(LB_cuadrados8_p, 3)), row.names = FALSE)
cat("\nF = var(inicio-2012) / var(2013-2019):\n")
print(transform(FV, razon_varianzas = round(razon_varianzas, 1), p_valor = signif(p_valor, 3)), row.names = FALSE)

# Barrido de la fecha de quiebre (hasta 2019, para que el COVID no domine), con >= 12 trimestres por lado
barrido <- do.call(rbind, lapply(names(X), function(s) {
  e <- as.numeric(corte(E[[s]], periodos[["sin COVID (hasta 2019Q4)"]]))
  t0 <- as.numeric(time(corte(E[[s]], periodos[["sin COVID (hasta 2019Q4)"]])))
  do.call(rbind, lapply(13:(length(e) - 12), function(j)
    data.frame(serie = s, fecha_quiebre = t0[j], razon = var(e[1:(j - 1)]) / var(e[j:length(e)]))))
}))
write.csv(barrido, "output/tablas/06_barrido_quiebre.csv", row.names = FALSE)
cat("\nFecha de quiebre que maximiza la razon de varianzas antes/despues (hasta 2019):\n")
print(do.call(rbind, lapply(split(barrido, barrido$serie), function(d) {
  j <- which.max(d$razon); data.frame(serie = d$serie[1], fecha = etiqueta_trim(ts(0, start = d$fecha_quiebre[j], frequency = 4)),
                                      razon = round(d$razon[j], 1))
})), row.names = FALSE)

# ---- B. Relacion PBI - indicador ----------------------------------------------------------
dw <- function(e) sum(diff(e)^2) / sum(e^2)
F3 <- function(z) diff(diff(z, 4))
rel <- list()
for (ind in indicadores) for (pn in names(periodos)) {
  per <- periodos[[pn]]
  a <- max(per[1], time(E[[ind]])[1], time(E[["PBI real"]])[1])
  if (pn == "inicio-2012Q4" && ind == "Ferrocarril") next      # coincide con 2006Q2-2012Q4
  # (1) regresion estatica en el filtro comun (1-B)(1-B^4), con los AO filtrados igual
  y <- window(F3(X[["PBI real"]]), start = a, end = per[2]); x <- window(F3(X[[ind]]), start = a, end = per[2])
  ao <- cbind(window(F3(AO(X[["PBI real"]], 2008.75)), start = a, end = per[2]),
              window(F3(AO(X[["PBI real"]], 2020)), start = a, end = per[2]))
  ao <- ao[, colSums(abs(ao)) > 0, drop = FALSE]
  m <- if (ncol(ao)) lm(as.numeric(y) ~ as.numeric(x) + ao) else lm(as.numeric(y) ~ as.numeric(x))
  sm <- summary(m); e <- residuals(m)
  # (2) innovaciones univariadas
  ei <- window(E[["PBI real"]], start = a, end = per[2]); ex <- window(E[[ind]], start = a, end = per[2])
  rel[[length(rel) + 1]] <- data.frame(
    indicador = ind, periodo = pn, desde = etiqueta_trim(y)[1], n = length(y),
    pendiente = coef(m)[2], t_ingenuo = sm$coefficients[2, 3], R2 = sm$r.squared, DW = dw(e),
    alarma_R2_mayor_DW = sm$r.squared > dw(e), LB8_p = Box.test(e, 8, "Ljung-Box")$p.value,
    cor_innov = cor(ei, ex), resp_innov = unname(coef(lm(as.numeric(ei) ~ as.numeric(ex)))[2]), row.names = NULL)
}
REL <- do.call(rbind, rel)
write.csv(REL, "output/tablas/06_relacion.csv", row.names = FALSE)
cat("\n=== B. PBI vs indicadores: regresion en (1-B)(1-B^4) con AO e innovaciones, por subperiodo ===\n")
num <- sapply(REL, is.numeric); RELp <- REL; RELp[num] <- lapply(REL[num], round, 3)
print(RELp, row.names = FALSE)
cat("resp_innov: pp que se mueve la sorpresa del PBI por cada pp de sorpresa del indicador.\n")

# ---- C. Figura ----------------------------------------------------------------------------
rol <- do.call(rbind, lapply(names(X), function(s) {
  e <- as.numeric(E[[s]]); t0 <- as.numeric(time(E[[s]])); j <- 8:length(e)
  data.frame(serie = s, t = t0[j], sd = sapply(j, function(i) sd(e[(i - 7):i])))
}))
etq <- do.call(rbind, lapply(split(rol, rol$serie), function(d) d[which.max(d$t), ]))
etq$sd_lab <- etq$sd + c("PBI real" = 0, "Electricidad" = -0.25, "Ferrocarril" = 0.25)[etq$serie]
g1 <- ggplot(rol, aes(t, sd, colour = serie)) +
  annotate("rect", xmin = 2020, xmax = 2021, ymin = -Inf, ymax = Inf, fill = gris, alpha = 0.45) +
  geom_vline(xintercept = 2013, colour = tinta2, linewidth = 0.3, linetype = "22") +
  geom_line(linewidth = 0.9) +
  geom_text(data = etq, aes(y = sd_lab, label = serie), hjust = 0, nudge_x = 0.4, size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = colores, guide = "none") +
  scale_x_continuous(breaks = seq(2005, 2025, 5), expand = expansion(mult = c(0.01, 0.16))) +
  labs(title = "Volatilidad de las sorpresas (innovaciones)",
       subtitle = "desvio movil de 8 trimestres, pp del log-nivel; linea punteada: 2013", x = NULL, y = NULL) + tema
sc <- do.call(rbind, lapply(indicadores, function(ind) do.call(rbind, lapply(c("2006Q2-2012Q4", "2013Q1-2019Q4", "2020Q1-2023Q4"), function(pn) {
  per <- periodos[[pn]]
  data.frame(indicador = ind, periodo = pn, x = as.numeric(window(E[[ind]], start = per[1], end = per[2])),
             y = as.numeric(window(E[["PBI real"]], start = per[1], end = per[2])))
}))))
txt <- REL[REL$periodo %in% unique(sc$periodo), c("indicador", "periodo", "cor_innov", "resp_innov")]
txt$lab <- sprintf("r = %.2f\nrespuesta = %.2f", txt$cor_innov, txt$resp_innov)
g2 <- ggplot(sc, aes(x, y)) +
  geom_hline(yintercept = 0, colour = gris, linewidth = 0.3) + geom_vline(xintercept = 0, colour = gris, linewidth = 0.3) +
  geom_point(aes(colour = indicador), size = 1.6, alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = tinta, linewidth = 0.5) +
  geom_text(data = txt, aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.2, size = 3, colour = tinta2, lineheight = 0.95) +
  scale_colour_manual(values = colores, guide = "none") +
  facet_grid(indicador ~ periodo, scales = "free_x") +
  labs(title = "Sorpresa del PBI vs sorpresa de cada indicador (muestra comun desde 2006Q2)",
       subtitle = "residuos de los modelos univariados (pp); respuesta = pendiente de MCO",
       x = "Sorpresa del indicador", y = "Sorpresa del PBI") +
  tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
f <- g1 / g2 + plot_layout(heights = c(1, 1.6)) +
  plot_annotation(title = "Despues de 2013 el PBI oficial casi no tiene sorpresas; la electricidad y el ferrocarril si",
                  caption = "Modelos del script 05 (con atipicos aditivos en 2008Q4 y 2020Q1). Fuente: NBS (chinadata.live), OCDE/FRED.",
                  theme = tema_global)
ggsave("output/figuras/06_quiebre_relacion.png", f, width = 11, height = 11, dpi = 150, bg = fondo)
cat("\nFigura: output/figuras/06_quiebre_relacion.png\n")
