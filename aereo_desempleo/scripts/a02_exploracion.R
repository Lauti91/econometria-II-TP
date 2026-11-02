# =============================================================================
# a02_exploracion.R -- Cada serie por separado (reglas 4 a 8 de la clase)
#   1. Graficos con contexto historico                                                regla 5
#   2. Rango-media por anio (Nota 6, Sec. 5.4)
#   3. FAC y FACP por filtro, muestra principal 2000-2019 (Nota 6, Sec. 9.3)           regla 6
#   4. Estacionalidad: fija vs. cambiante entre subperiodos (Nota 6, Sec. 2.3)        regla 8
#   5. ADF con dummies, Dickey-Pantula, por subperiodo, incluida la muestra con COVID  regla 7
# Muestra principal: 2000-01 a 2019-12 (240 meses). El COVID (abr-2020: -95% de millas) no es un
# atipico sino un cambio transitorio de varios anios: se trata como subperiodo aparte y queda para
# el analisis de intervencion (segunda parte). Ejecutar desde aereo_desempleo/.
# =============================================================================
source("scripts/utils_a.R")
D <- cargar_base(); B <- D$B; R <- recesiones(D$rec); t <- as.numeric(time(D$rpm))
rpm <- window(D$rpm, end = c(2019, 12)); des <- window(D$des, end = c(2019, 12))

# ---- 1. Niveles con contexto -----------------------------------------------------------------------
hito <- function(t0, lab, y_lab, serie) data.frame(t = t0, lab = lab, y_lab = y_lab, y = approx(t, serie, t0)$y)
capa_hitos <- function(H) list(
  geom_segment(data = H, aes(x = t, xend = t, y = y, yend = y_lab), colour = tinta2, linewidth = 0.3, linetype = "22"),
  geom_point(data = H, aes(t, y), colour = tinta, size = 1.8),
  geom_label(data = H, aes(t, y_lab, label = lab), colour = tinta, fill = fondo, size = 2.8, label.size = 0, lineheight = 0.9, label.padding = unit(0.12, "lines")))
rpm_mm <- B$rpm_miles / 1e6                                   # miles de millones de millas
H1 <- rbind(hito(fecha_num(2001, 9), "11-S\nsep-2001", 20, rpm_mm), hito(fecha_num(2009, 2), "Gran Recesion\nfeb-2009", 25, rpm_mm),
            hito(fecha_num(2020, 4), "COVID abr-2020\n(-95%)", 30, rpm_mm), hito(fecha_num(2023, 7), "Supera 2019\njul-2023", 118, rpm_mm))
g1 <- ggplot(data.frame(t = t, y = rpm_mm)) + capa_recesion(R) +
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F5E6DA", alpha = 0.5) +
  geom_line(aes(t, y), colour = col_rpm, linewidth = 0.55) + capa_hitos(H1) +
  scale_x_continuous(breaks = seq(2000, 2026, 2)) + scale_y_continuous(limits = c(0, 125)) +
  labs(title = "Millas-pasajero pagas, aerolineas de EE.UU. (miles de millones por mes, sin desestacionalizar)",
       subtitle = "gris: recesiones NBER; rosado: 2020-2026, fuera de la muestra principal (COVID)", x = NULL, y = NULL) + tema
H2 <- rbind(hito(fecha_num(2003, 6), "Pos burbuja .com\njun-2003", 9, B$desempleo_nsa), hito(fecha_num(2010, 1), "Maximo\nene-2010", 13, B$desempleo_nsa),
            hito(fecha_num(2020, 4), "COVID\nabr-2020", 16.5, B$desempleo_nsa), hito(fecha_num(2025, 10), "Sin dato oct-2025\n(cierre del gobierno)", 9.5, B$desempleo_nsa))
g2 <- ggplot(data.frame(t = t, nsa = B$desempleo_nsa, sa = B$desempleo_sa)) + capa_recesion(R) +
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F5E6DA", alpha = 0.5) +
  geom_line(aes(t, nsa), colour = col_des, linewidth = 0.55) + geom_line(aes(t, sa), colour = tinta, linewidth = 0.5) + capa_hitos(H2) +
  scale_x_continuous(breaks = seq(2000, 2026, 2)) + scale_y_continuous(limits = c(0, 18)) +
  labs(title = "Tasa de desempleo, EE.UU. (%)", subtitle = "naranja: sin desestacionalizar (UNRATENSA, la que se modela); negro: desestacionalizada (UNRATE)",
       x = NULL, y = NULL) + tema
ggsave("output/figuras/a02_niveles.png", (g1 / g2) + plot_annotation(theme = tema_global), width = 12, height = 8.5, dpi = 200)

# ---- 2. Rango-media (muestra principal) ---------------------------------------------------------------
rango_media <- function(y, anio) {
  M <- do.call(rbind, tapply(y, anio, function(v) c(m = mean(v), r = diff(range(v)))))
  b <- coef(lm(log(M[, "r"]) ~ log(M[, "m"])))[2]; c(pendiente = unname(b), lambda = unname(1 - b))
}
anio <- floor(as.numeric(time(rpm)) + 1e-6)
BC <- rbind(rpm = rango_media(exp(as.numeric(rpm) / 100), anio), desempleo = rango_media(as.numeric(des), anio))
cat("\n=== 2. Rango-media anual 2000-2019: lambda = 1 - b (0 = logaritmo, 1 = niveles) ===\n"); print(round(BC, 2))
write.csv(data.frame(serie = rownames(BC), BC), "output/tablas/a02_box_cox.csv", row.names = FALSE)

# ---- 3. FAC y FACP por filtro ----------------------------------------------------------------------------
corr_df <- function(x, nombre, L = 48) {
  a <- acf(x, L, plot = FALSE)$acf[-1]; p <- pacf(x, L, plot = FALSE)$acf
  rbind(data.frame(serie = nombre, tipo = "FAC", k = 1:L, r = a), data.frame(serie = nombre, tipo = "FACP", k = 1:L, r = p))
}
filtros <- list("ln RPM: nivel" = rpm, "ln RPM: (1-B)" = diff(rpm), "ln RPM: (1-B^12)" = diff(rpm, 12), "ln RPM: (1-B)(1-B^12)" = diff(diff(rpm, 12)),
                "desempleo: nivel" = des, "desempleo: (1-B)" = diff(des), "desempleo: (1-B^12)" = diff(des, 12), "desempleo: (1-B)(1-B^12)" = diff(diff(des, 12)))
CR <- do.call(rbind, Map(corr_df, filtros, names(filtros))); CR$serie <- factor(CR$serie, levels = names(filtros))
CR$banda <- 1.96 / sqrt(sapply(filtros, length)[as.character(CR$serie)]); CR$sig <- abs(CR$r) > CR$banda
CR$col <- ifelse(grepl("^desempleo", CR$serie), "des", "rpm")
bandas <- unique(CR[, c("serie", "banda")])
g3 <- ggplot(CR, aes(k, r)) +
  geom_rect(data = bandas, aes(xmin = 0, xmax = 49, ymin = -banda, ymax = banda), inherit.aes = FALSE, fill = gris, alpha = 0.45) +
  geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_vline(xintercept = c(12, 24, 36, 48), colour = "#E4E7EB", linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = interaction(col, sig)), linewidth = 0.9, lineend = "round") +
  scale_colour_manual(values = c(rpm.TRUE = col_rpm, des.TRUE = col_des, rpm.FALSE = "#9AA5B1", des.FALSE = "#9AA5B1"), guide = "none") +
  facet_grid(serie ~ tipo, switch = "y") + scale_x_continuous(breaks = c(1, 12, 24, 36, 48)) + scale_y_continuous(limits = c(-1, 1), breaks = c(-0.5, 0, 0.5)) +
  labs(title = "Correlogramas simple y parcial por filtro, 2000-2019 (48 rezagos; banda +-1.96/raiz(T))", x = "Rezago (meses)", y = NULL) +
  tema + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, colour = tinta, size = 8.5), strip.placement = "outside",
               strip.text.x = element_text(colour = tinta, face = "bold"))
ggsave("output/figuras/a02_correlogramas.png", g3 + plot_annotation(theme = tema_global), width = 12, height = 14, dpi = 200)

# ---- 4. Estacionalidad ------------------------------------------------------------------------------------
tm <- as.numeric(time(rpm))
subper <- cut(tm, c(-Inf, 2007, 2014, Inf), labels = c("2000-2006", "2007-2013", "2014-2019"), right = FALSE)
estac <- function(x, nombre, excluir = NULL) {
  dx <- as.numeric(diff(x)); mes <- factor(cycle(x)[-1]); sp <- subper[-1]
  if (length(excluir)) { ok <- !(round(as.numeric(time(x))[-1], 4) %in% round(excluir, 4)); dx[!ok] <- NA }
  m0 <- lm(dx ~ sp); m1 <- lm(dx ~ sp + mes); m2 <- lm(dx ~ sp * mes); a <- anova(m0, m1); b <- anova(m1, m2)
  data.frame(serie = nombre, F_fija = a$F[2], p_fija = a$`Pr(>F)`[2], F_cambio = b$F[2], p_cambio = b$`Pr(>F)`[2], R2_dummies = summary(m1)$r.squared)
}
# El 11-S cae dentro de 2000-2006 y puede fabricar un "cambio de patron": se repite el test sin las
# variaciones de sep-2001 a ene-2002 (las afectadas por la caida y la recuperacion).
ES <- rbind(estac(rpm, "ln RPM (x100)"), estac(rpm, "ln RPM sin sep-2001 a ene-2002", fecha_num(2001, 9:12) |> c(fecha_num(2002, 1))),
            estac(des, "desempleo (pp)"))
cat("\n=== 4. Estacionalidad en la primera diferencia, 2000-2019: fija y cambio de patron entre 2000-06, 2007-13, 2014-19 ===\n")
print(transform(ES, F_fija = round(F_fija, 2), p_fija = signif(p_fija, 3), F_cambio = round(F_cambio, 2), p_cambio = signif(p_cambio, 3), R2_dummies = round(R2_dummies, 3)), row.names = FALSE)
write.csv(ES, "output/tablas/a02_estacionalidad.csv", row.names = FALSE)
perfil <- rbind(data.frame(serie = "1. ln RPM: variacion mensual promedio (%)", aggregate(list(dx = as.numeric(diff(rpm))), list(mes = cycle(rpm)[-1], subperiodo = subper[-1]), mean)),
                data.frame(serie = "2. Desempleo: variacion mensual promedio (pp)", aggregate(list(dx = as.numeric(diff(des))), list(mes = cycle(des)[-1], subperiodo = subper[-1]), mean)))
g4 <- ggplot(perfil, aes(mes, dx, colour = subperiodo)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) + facet_wrap(~serie, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("#9AA5B1", "#52606D", col_rpm), name = NULL) +
  scale_x_continuous(breaks = 1:12, labels = c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")) +
  labs(title = "Patron estacional por subperiodo (2000-2019)",
       subtitle = sprintf("F cambio de patron: RPM %.2f (p = %.2g); RPM sin 11-S %.2f (p = %.2g); desempleo %.2f (p = %.2g)",
                          ES$F_cambio[1], ES$p_cambio[1], ES$F_cambio[2], ES$p_cambio[2], ES$F_cambio[3], ES$p_cambio[3]), x = NULL, y = NULL) +
  tema + theme(legend.position = "top", strip.text = element_text(colour = tinta, face = "bold"))
ggsave("output/figuras/a02_estacionalidad.png", g4 + plot_annotation(theme = tema_global), width = 10, height = 7.5, dpi = 200)

# ---- 5. ADF por subperiodo --------------------------------------------------------------------------------
muestras <- list("2000-2019" = c(2000, 2019.95), "2000-2009" = c(2000, 2009.95), "2010-2019" = c(2010, 2019.95),
                 "completa 2000-2026 (con COVID)" = c(2000, 2026.5))
adf_tabla <- function(x, nombre) {
  do.call(rbind, lapply(names(muestras), function(mu) {
    z <- window(x, start = muestras[[mu]][1], end = min(muestras[[mu]][2], max(time(x))))
    do.call(rbind, lapply(c(6, 12), function(p) {
      a1 <- adf_est(diff(z), p, tendencia = FALSE); a2 <- adf_est(z, p, tendencia = FALSE); a3 <- adf_est(z, p, tendencia = TRUE)
      data.frame(serie = nombre, muestra = mu, p = p, n = a3["n"], Dx_tau_mu = a1["tau"], x_tau_mu = a2["tau"], x_tau_tau = a3["tau"],
                 t_tend = a3["t_tend"], LB36_p = a3["LB_p"], row.names = NULL)
    }))
  }))
}
ADF <- rbind(adf_tabla(D$rpm, "ln RPM"), adf_tabla(D$des, "desempleo"))
ADF$Dx_rechaza <- ADF$Dx_tau_mu < VC_ADF["tau_mu", "5%"]
ADF$x_rechaza_mu <- ADF$x_tau_mu < VC_ADF["tau_mu", "5%"]; ADF$x_rechaza_tau <- ADF$x_tau_tau < VC_ADF["tau_tau", "5%"]
cat("\n=== 5. ADF con dummies mensuales centradas (VC 5%, n grande: tau_mu -2.87, tau_tau -3.42; con n = 120 son algo mas exigentes) ===\n")
num <- sapply(ADF, is.numeric); P <- ADF; P[num] <- lapply(ADF[num], round, 2); print(P, row.names = FALSE)
write.csv(ADF, "output/tablas/a02_adf.csv", row.names = FALSE)
