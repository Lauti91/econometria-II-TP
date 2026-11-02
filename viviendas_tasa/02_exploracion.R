# =============================================================================
# 02_exploracion.R -- Cada serie antes de cualquier modelo
#   1. Graficos con contexto historico (regla 5)
#   2. Por que la muestra empieza en 1985: la volatilidad de la tasa en el regimen Volcker
#   3. Transformacion estabilizadora de varianza (Nota 6, Sec. 5.4): rango-mediana y perfil de Box-Cox
#   4. Estacionalidad (regla 8): F de dummies fijas y F de cambio del patron entre subperiodos
# Ejecutar desde viviendas_tasa/:  Rscript 02_exploracion.R
# =============================================================================
source("R/funciones.R")
LC <- cargar_base(completa = TRUE); L <- cargar_base()
R <- recesiones(LC$rec); t <- as.numeric(time(LC$inicios))

# ---- 1. Niveles con contexto (muestra completa 1971-2026; la principal empieza en 1985) -------------
hito <- function(a, m, lab, y_lab, serie) { t0 <- fecha_num(a, m); data.frame(t = t0, lab = lab, y_lab = y_lab, y = approx(t, serie, t0)$y) }
capa_hitos <- function(H) list(
  geom_segment(data = H, aes(x = t, xend = t, y = y, yend = y_lab), colour = tinta2, linewidth = 0.3, linetype = "22"),
  geom_point(data = H, aes(t, y), colour = ifelse(grepl("COVID|Volcker", H$lab), rojo, tinta), size = 1.8),
  geom_label(data = H, aes(t, y_lab, label = lab), colour = ifelse(grepl("COVID|Volcker", H$lab), rojo, tinta), fill = "white", size = 2.7, label.size = 0, lineheight = 0.9,
             label.padding = unit(0.12, "lines")))
sa <- as.numeric(LC$inicios_sa)
Hv <- rbind(hito(1972, 12, "Boom 1972", 245, sa), hito(1982, 1, "Volcker: minimo 1982", 22, sa),
            hito(2006, 1, "Pico de la burbuja\nene-2006", 238, sa), hito(2009, 4, "Minimo abr-2009", 12, sa),
            hito(2020, 4, "COVID abr-2020", 40, sa), hito(2022, 4, "Suba de tasas 2022", 205, sa))
inicio_muestra <- annotate("rect", xmin = -Inf, xmax = 1985, ymin = -Inf, ymax = Inf, fill = "#F5F7FA", alpha = 0.9)
g1 <- ggplot(data.frame(t = t, nsa = as.numeric(LC$inicios), sa = sa)) + inicio_muestra + capa_recesion(R) +
  geom_vline(xintercept = 1985, colour = rojo, linewidth = 0.5, linetype = "22") +
  geom_line(aes(t, nsa), colour = col_viv, linewidth = 0.35, alpha = 0.85) + geom_line(aes(t, sa), colour = tinta, linewidth = 0.6) +
  capa_hitos(Hv) + scale_x_continuous(breaks = seq(1970, 2025, 5)) + scale_y_continuous(limits = c(0, 260)) +
  labs(title = "Inicios de viviendas privadas, EE.UU. (miles de unidades por mes)",
       subtitle = "azul: sin desestacionalizar (HOUSTNSA, la que se modela); negro: desestacionalizada (HOUST/12); gris: recesiones NBER; la muestra principal empieza en 1985",
       x = NULL, y = NULL) + tema
tp <- as.numeric(LC$tasa_pp)
Ht <- rbind(hito(1981, 10, "Volcker: 18,6% oct-1981", 20.2, tp), hito(2008, 11, "Fed compra MBS\nnov-2008", 11.5, tp),
            hito(2021, 1, "Minimo 2,7%\nene-2021", 0.8, tp), hito(2023, 10, "Maximo 7,8%\noct-2023", 15.5, tp))
g2 <- ggplot(data.frame(t = t, tasa = tp, bono = as.numeric(LC$bono))) + inicio_muestra + capa_recesion(R) +
  geom_vline(xintercept = 1985, colour = rojo, linewidth = 0.5, linetype = "22") +
  geom_line(aes(t, bono), colour = col_bono, linewidth = 0.5) + geom_line(aes(t, tasa), colour = col_tasa, linewidth = 0.7) +
  capa_hitos(Ht) + scale_x_continuous(breaks = seq(1970, 2025, 5)) + scale_y_continuous(limits = c(0, 21)) +
  labs(title = "Tasa hipotecaria fija a 30 anios (% anual, promedio mensual de la serie semanal)",
       subtitle = "naranja: MORTGAGE30US; gris: bono del Tesoro a 10 anios (GS10), al que la tasa sigue de cerca", x = NULL, y = NULL) + tema
guardar_fig(g1 / g2, "02_niveles.png", 12, 8.5)

# ---- 2. Volatilidad de la tasa por subperiodo -------------------------------------------------------
dtp <- diff(LC$tasa_pp); tt <- as.numeric(time(dtp))
tramos <- list("1971--1978" = c(1971, 1979), "1979--1984 (Volcker)" = c(1979, 1985), "1985--2007" = c(1985, 2008),
               "2008--2019" = c(2008, 2020), "2020--2026" = c(2020, 2027))
VT <- data.frame(periodo = names(tramos), desvio = sapply(tramos, function(r) sd(dtp[tt >= r[1] & tt < r[2]])))
tabla_tex(data.frame(VT$periodo, num(VT$desvio, 2)), "02_volatilidad_tasa.tex", rojo = celdas(2, 1:2), "lc", c("Per\\'iodo", "Desv\\'io de $\\Delta$ tasa (pp)"))
guardar_valor("VolVolcker", num(VT$desvio[2], 2)); guardar_valor("VolNormal", num(VT$desvio[3], 2))
guardar_valor("VolRatio", num(VT$desvio[2] / VT$desvio[3], 1))

# ---- 3. Transformacion: rango-mediana y perfil de verosimilitud de Box-Cox ----------------------------
RMv <- rango_mediana(L$inicios); RMt <- rango_mediana(L$tasa_pp)
D <- rbind(data.frame(serie = "Inicios (miles)", RMv$datos), data.frame(serie = "Tasa (%)", RMt$datos))
etiq <- c(`Inicios (miles)` = sprintf("pendiente %.2f  ->  lambda = %.2f", RMv$b, RMv$lambda),
          `Tasa (%)` = sprintf("pendiente %.2f  ->  lambda = %.2f", RMt$b, RMt$lambda))
g3 <- ggplot(D, aes(log(mediana), log(rango), colour = serie)) + geom_point(size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, formula = y ~ x) +
  facet_wrap(~serie, scales = "free", labeller = labeller(serie = function(s) paste0(s, ": ", etiq[s]))) +
  scale_colour_manual(values = c(col_viv, col_tasa), guide = "none") +
  labs(title = "Relacion rango-mediana por anio (1985-2026)", subtitle = "ln(rango) contra ln(mediana); lambda = 1 - pendiente (0: logaritmo, 1/2: raiz, 1: sin transformar)",
       x = "ln(mediana del anio)", y = "ln(rango del anio)") + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
guardar_fig(g3, "02_rango_mediana.png", 10, 4.2)

# Perfil de verosimilitud (Nota 6, Sec. 5.4, opcion ii): l(lambda) + (lambda - 1) * sum ln y (jacobiano),
# con el modelo que despues se identifica para cada serie (04 y 05)
perfil <- function(y, orden, estacional, xreg = NULL) sapply(c(0, 0.25, 0.5, 0.75, 1), function(l) {
  z <- ts(box_cox(as.numeric(y), l), start = start(y), frequency = S)
  m <- arima(z, order = orden, seasonal = list(order = estacional, period = S), xreg = xreg, method = "ML")
  d <- orden[2] + S * estacional[2]; m$loglik + (l - 1) * sum(log(as.numeric(y)[(d + 1):length(y)]))
})
PV <- data.frame(lambda = c(0, 0.25, 0.5, 0.75, 1),
                 inicios = perfil(L$inicios, c(0, 1, 1), c(0, 1, 1), indicadora(L$inicios, 2020, 4)),
                 tasa = perfil(L$tasa_pp, c(0, 1, 1), c(0, 0, 0)))
PV$inicios <- PV$inicios - max(PV$inicios); PV$tasa <- PV$tasa - max(PV$tasa)
tabla_tex(data.frame(c("0 (logaritmo)", "0,25", "0,5 (ra\\'iz)", "0,75", "1 (sin transformar)"), num(PV$inicios, 1), num(PV$tasa, 1)),
          "02_perfil_boxcox.tex", "lcc", c("$\\lambda$", "Inicios", "Tasa"))
guardar_valor("LamInicios", num(RMv$lambda, 2)); guardar_valor("LamTasa", num(RMt$lambda, 2))
guardar_valor("LamInicioslo", num(RMv$ic[1], 2)); guardar_valor("LamIniciosup", num(RMv$ic[2], 2))
guardar_valor("LamTasalo", num(RMt$ic[1], 2)); guardar_valor("LamTasaup", num(RMt$ic[2], 2))
guardar_valor("PerfilIniciosLog", num(PV$inicios[1], 1)); guardar_valor("PerfilTasaLog", num(PV$tasa[1], 1))
guardar_valor("PerfilIniciosUno", num(PV$inicios[5], 1)); guardar_valor("PerfilTasaUno", num(PV$tasa[5], 1))
cat("Rango-mediana: inicios lambda =", round(RMv$lambda, 2), " tasa lambda =", round(RMt$lambda, 2), "\n"); print(round(PV, 1))

# ---- 4. Estacionalidad --------------------------------------------------------------------------------
# Regresion de la primera diferencia sobre 11 dummies centradas: F de estacionalidad fija; y F de
# cambio del patron entre 1985-2007 y 2008-2026 (interaccion mes x subperiodo). Si el patron migra,
# la estacionalidad es estocastica y se trata con (1 - B^12) (Nota 6, Sec. 2.3).
estac <- function(x) {
  dx <- diff(x); mes <- factor(cycle(dx)); sp <- factor(as.numeric(time(dx)) >= 2008)
  m0 <- lm(dx ~ sp); m1 <- lm(dx ~ sp + mes); m2 <- lm(dx ~ sp * mes); a <- anova(m0, m1); b <- anova(m1, m2)
  c(F_fija = a$F[2], p_fija = a$`Pr(>F)`[2], gl1 = a$Df[2], gl2 = a$Res.Df[2], F_cambio = b$F[2], p_cambio = b$`Pr(>F)`[2],
    gl1c = b$Df[2], gl2c = b$Res.Df[2], R2 = summary(m1)$r.squared)
}
Ev <- estac(L$viv); Et <- estac(L$tasa)
tabla_tex(data.frame(c("$\\Delta\\ln$ inicios", "$\\Delta\\ln$ tasa"),
                     c(sprintf("$F_{%d,%d}=%s$", Ev["gl1"], Ev["gl2"], num(Ev["F_fija"], 1)), sprintf("$F_{%d,%d}=%s$", Et["gl1"], Et["gl2"], num(Et["F_fija"], 2))),
                     c(pval(Ev["p_fija"]), pval(Et["p_fija"])), c(num(Ev["R2"], 2), num(Et["R2"], 3)),
                     c(sprintf("$F_{%d,%d}=%s$", Ev["gl1c"], Ev["gl2c"], num(Ev["F_cambio"], 2)), sprintf("$F_{%d,%d}=%s$", Et["gl1c"], Et["gl2c"], num(Et["F_cambio"], 2))),
                     c(pval(Ev["p_cambio"]), pval(Et["p_cambio"]))),
          "02_estacionalidad.tex", rojo = rbind(celdas(1, 5:6), celdas(2, 2:3)), "lccccc", c("Serie", "Estac.\\ fija", "$p$", "$R^2$ dummies", "Cambio de patr\\'on", "$p$"))
guardar_valor("EstFv", num(Ev["F_fija"], 1)); guardar_valor("EstRv", num(Ev["R2"], 2)); guardar_valor("EstFcv", num(Ev["F_cambio"], 2))
guardar_valor("EstPcv", pval(Ev["p_cambio"])); guardar_valor("EstFt", num(Et["F_fija"], 2)); guardar_valor("EstPt", pval(Et["p_fija"]))
perfil_m <- aggregate(list(dx = as.numeric(diff(L$viv))), list(mes = cycle(diff(L$viv)), sub = ifelse(as.numeric(time(diff(L$viv))) >= 2008, "2008-2026", "1985-2007")), mean)
g4 <- ggplot(perfil_m, aes(mes, dx, colour = sub)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) + scale_colour_manual(values = c("#9AA5B1", col_viv), name = NULL) +
  scale_x_continuous(breaks = 1:12, labels = c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")) +
  labs(title = "Patron estacional de los inicios: variacion mensual promedio de 100 ln(inicios), por subperiodo",
       subtitle = "La amplitud se achica despues de 2008: el patron estacional no es fijo", x = NULL, y = "% mensual") +
  tema + theme(legend.position = "top")
guardar_fig(g4, "02_estacionalidad.png", 9, 4.2)
escribir_valores("valores_02.tex")
