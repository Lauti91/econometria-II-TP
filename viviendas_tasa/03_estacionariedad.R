# =============================================================================
# 03_estacionariedad.R -- Cuantas diferencias hacen falta (regla 6: el correlograma manda; regla 7: el
# Dickey-Fuller no es una receta). Muestra principal 1985-01 a 2026-08, en logaritmos (02).
#   1. Correlogramas por filtro (FAC y FACP con banda +-1,96/raiz(T), Nota 6 Sec. 9.1)
#      inicios: nivel, (1-B), (1-B^12), (1-B)(1-B^12); tasa: nivel, (1-B)
#   2. ADF con el esquema descendente de Dickey-Pantula (Nota 6, Sec. 4.5): primero se contrasta una
#      segunda raiz (ADF sobre la diferencia) y recien despues la raiz del nivel, con tau_tau y, si la
#      tendencia no es significativa, tau_mu (Sec. 4.3). Rezagos: 12 (Said-Dickey n^(1/3) ~ 8, y al
#      menos un anio con datos mensuales); sensibilidad con 24 y con el p que minimiza el SBC.
#   3. ADF sobre (1-B^12)x: si no rechaza, a la diferencia estacional le falta un (1-B)  -> filtro F3
# Ejecutar desde viviendas_tasa/:  Rscript 03_estacionariedad.R
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); viv <- L$viv; tasa <- L$tasa

# ---- 1. Correlogramas -------------------------------------------------------------------------------
filtros <- list("ln inicios: nivel" = viv, "ln inicios: (1-B)" = diff(viv), "ln inicios: (1-B^12)" = diff(viv, 12),
                "ln inicios: (1-B)(1-B^12)" = diff(diff(viv, 12)), "ln tasa: nivel" = tasa, "ln tasa: (1-B)" = diff(tasa))
K <- 36
CR <- do.call(rbind, Map(function(z, nm) { C <- correlograma(z, K)
  rbind(data.frame(serie = nm, tipo = "FAC", k = C$k, r = C$AC, banda = C$banda), data.frame(serie = nm, tipo = "FACP", k = C$k, r = C$PAC, banda = C$banda)) },
  filtros, names(filtros)))
CR$serie <- factor(CR$serie, levels = names(filtros)); CR$sig <- abs(CR$r) > CR$banda
CR$col <- ifelse(grepl("tasa", CR$serie), "tasa", "viv")
CR$col[CR$serie %in% c("ln inicios: (1-B)(1-B^12)", "ln tasa: (1-B)")] <- "final"     # filtros elegidos: en rojo
g <- ggplot(CR, aes(k, r)) +
  geom_rect(aes(xmin = 0, xmax = K + 1, ymin = -banda, ymax = banda), fill = gris, alpha = 0.05) +
  geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_vline(xintercept = c(12, 24, 36), colour = "#E4E7EB", linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = interaction(col, sig)), linewidth = 1, lineend = "round") +
  scale_colour_manual(values = c(viv.TRUE = col_viv, tasa.TRUE = col_tasa, final.TRUE = rojo, viv.FALSE = "#9AA5B1", tasa.FALSE = "#9AA5B1", final.FALSE = "#9AA5B1"), guide = "none") +
  facet_grid(serie ~ tipo, switch = "y") + scale_x_continuous(breaks = c(1, 12, 24, 36)) + scale_y_continuous(limits = c(-1, 1), breaks = c(-0.5, 0, 0.5)) +
  labs(title = "FAC y FACP por filtro, 1985-2026 (36 rezagos)", subtitle = "En color: fuera de la banda +-1,96/raiz(T). En rojo: los filtros elegidos. Lineas verticales en 12, 24 y 36.",
       x = "Rezago (meses)", y = NULL) +
  tema + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, colour = tinta, size = 8.5), strip.placement = "outside",
               strip.text.x = element_text(colour = tinta, face = "bold"))
guardar_fig(g, "03_correlogramas.png", 11, 10)

# Salidas tipo EViews (24 rezagos) de las cinco series clave
corr_ev <- list(
  list(viv, "Correlograma de LHOUST (nivel)", "03_corr_lhoust.png", c(1, 12, 24)),
  list(diff(viv), "Correlograma de D(LHOUST)", "03_corr_dlhoust.png", c(12, 24)),
  list(diff(diff(viv, 12)), "Correlograma de D(LHOUST,1,12)", "03_corr_d112lhoust.png", c(1, 11, 12, 13)),
  list(tasa, "Correlograma de LMORT (nivel)", "03_corr_lmort.png", c(1, 24)),
  list(diff(tasa), "Correlograma de D(LMORT)", "03_corr_dlmort.png", 1))
for (c_ in corr_ev) guardar_salida(salida_correlograma(c_[[1]], c_[[2]], 24, rojo_k = c_[[4]]), c_[[3]])
# Valores que se citan en el texto
Cv0 <- correlograma(viv, 24); Cv1 <- correlograma(diff(viv), 36); Cv3 <- correlograma(diff(diff(viv, 12)), 36)
Ct0 <- correlograma(tasa, 24); Ct1 <- correlograma(diff(tasa), 36)
guardar_valor("RhoNivUno", num(Cv0$AC[1], 2)); guardar_valor("RhoNivVeinticuatro", num(Cv0$AC[24], 2))
guardar_valor("RhoDifDoce", num(Cv1$AC[12], 2)); guardar_valor("RhoDifVeinticuatro", num(Cv1$AC[24], 2)); guardar_valor("RhoDifTreintayseis", num(Cv1$AC[36], 2))
Cv2 <- correlograma(diff(viv, 12), 24); guardar_valor("RhoFdosUno", num(Cv2$AC[1], 2)); guardar_valor("RhoFdosDoce", num(Cv2$AC[12], 2)); guardar_valor("RhoFdosVeinticuatro", num(Cv2$AC[24], 2))
guardar_valor("RhoFtresUno", num(Cv3$AC[1], 2)); guardar_valor("RhoFtresDoce", num(Cv3$AC[12], 2))
guardar_valor("RhoFtresOnce", num(Cv3$AC[11], 2)); guardar_valor("RhoFtresTrece", num(Cv3$AC[13], 2))
guardar_valor("RhoFtresVeinticuatro", num(Cv3$AC[24], 2))
guardar_valor("RhoTasaNivUno", num(Ct0$AC[1], 2)); guardar_valor("RhoTasaNivVeinticuatro", num(Ct0$AC[24], 2))
guardar_valor("RhoTasaDifUno", num(Ct1$AC[1], 2)); guardar_valor("RhoTasaDifDos", num(Ct1$AC[2], 2))
guardar_valor("PacTasaDifDos", num(Ct1$PAC[2], 2))

# ---- 2 y 3. ADF ---------------------------------------------------------------------------------------
# p elegido por SBC entre 0 y 24 sobre una muestra comun (la de p = 24)
p_sbc <- function(x, det, dummies) {
  sbc <- sapply(0:24, function(p) { A <- adf(x, p, det, dummies); e <- residuals(A$m)[(24 - p + 1):A$n]
    n <- length(e); k <- length(coef(A$m)); n * log(sum(e^2) / n) + k * log(n) })
  (0:24)[which.min(sbc)]
}
pruebas <- list(
  list(serie = "$\\ln$ inicios", paso = "1: $\\Delta x$ (2da ra\\'iz)", x = diff(viv), det = "c", dum = TRUE, v = "D(LHOUST)"),
  list(serie = "$\\ln$ inicios", paso = "2: nivel, $\\tau_\\tau$", x = viv, det = "ct", dum = TRUE, v = "LHOUST"),
  list(serie = "$\\ln$ inicios", paso = "2: nivel, $\\tau_\\mu$", x = viv, det = "c", dum = TRUE, v = "LHOUST"),
  list(serie = "$\\ln$ inicios", paso = "$(1-B^{12})x$", x = diff(viv, 12), det = "c", dum = FALSE, v = "D(LHOUST,0,12)"),
  list(serie = "$\\ln$ tasa", paso = "1: $\\Delta y$ (2da ra\\'iz)", x = diff(tasa), det = "c", dum = FALSE, v = "D(LMORT)"),
  list(serie = "$\\ln$ tasa", paso = "2: nivel, $\\tau_\\tau$", x = tasa, det = "ct", dum = FALSE, v = "LMORT"),
  list(serie = "$\\ln$ tasa", paso = "2: nivel, $\\tau_\\mu$", x = tasa, det = "c", dum = FALSE, v = "LMORT"))
AD <- do.call(rbind, lapply(pruebas, function(P) {
  ps <- p_sbc(P$x, P$det, P$dum)
  r <- lapply(c(12, 24, ps), function(p) adf(P$x, p, P$det, P$dum, P$v))
  data.frame(serie = P$serie, paso = P$paso, det = P$det, n = r[[1]]$n, tau12 = r[[1]]$tau, p12 = r[[1]]$p, vc5 = r[[1]]$vc["5%"],
             t_tend = r[[1]]$t_trend, tau24 = r[[2]]$tau, p24 = r[[2]]$p, p_sbc = ps, tau_sbc = r[[3]]$tau, psbc = r[[3]]$p, row.names = NULL)
}))
write.csv(AD, "salidas/tablas/03_adf.csv", row.names = FALSE)
print(transform(AD[, -1], tau12 = round(tau12, 2), p12 = round(p12, 3), tau24 = round(tau24, 2), p24 = round(p24, 3), tau_sbc = round(tau_sbc, 2), psbc = round(psbc, 3), t_tend = round(t_tend, 2)))
dec <- function(p) ifelse(p < 0.05, "Rechaza", "No rechaza")
tabla_tex(data.frame(AD$serie, AD$paso, AD$n, num(AD$tau12, 2), num(AD$vc5, 2), sapply(AD$p12, pval), dec(AD$p12),
                     num(AD$tau24, 2), sprintf("%s (%d)", num(AD$tau_sbc, 2), AD$p_sbc)),
          "03_adf.tex", rojo = rbind(celdas(1:7, 7), celdas(4, 9)), "llccccccc",
          c("Serie", "Contraste", "$n$", "$\\tau$ (12 rez.)", "VC 5\\%", "$p$", "Decisi\\'on", "$\\tau$ (24 rez.)", "$\\tau$ (SBC, $p$)"))
# Valores que se citan
v <- function(i, col, d = 2) num(AD[i, col], d)
guardar_valor("AdfVivDif", v(1, "tau12")); guardar_valor("AdfVivNivTT", v(2, "tau12")); guardar_valor("AdfVivNivTTp", pval(AD$p12[2]))
guardar_valor("AdfVivTend", v(2, "t_tend")); guardar_valor("AdfVivNivTM", v(3, "tau12")); guardar_valor("AdfVivNivTMp", pval(AD$p12[3]))
guardar_valor("AdfVivEst", v(4, "tau12")); guardar_valor("AdfVivEstp", pval(AD$p12[4])); guardar_valor("AdfVivEstVC", v(4, "vc5"))
guardar_valor("AdfVivEstSBC", v(4, "tau_sbc")); guardar_valor("AdfVivEstSBCp", AD$p_sbc[4]); guardar_valor("AdfVivEstSBCpv", pval(AD$psbc[4]))
guardar_valor("AdfTasaDif", v(5, "tau12")); guardar_valor("AdfTasaNivTT", v(6, "tau12")); guardar_valor("AdfTasaNivTTp", pval(AD$p12[6]))
guardar_valor("AdfTasaTend", v(6, "t_tend")); guardar_valor("AdfTasaNivTM", v(7, "tau12")); guardar_valor("AdfTasaNivTMp", pval(AD$p12[7]))
guardar_valor("AdfVCtt", v(2, "vc5")); guardar_valor("AdfVCtm", v(3, "vc5"))

# Salidas EViews de los ADF
ev <- list(c(2, "03_adf_lhoust.png", "Test ADF: LHOUST en nivel (tau_tau)"), c(1, "03_adf_dlhoust.png", "Test ADF: D(LHOUST)"),
           c(4, "03_adf_d12lhoust.png", "Test ADF: D(LHOUST,0,12)"), c(6, "03_adf_lmort.png", "Test ADF: LMORT en nivel (tau_tau)"),
           c(5, "03_adf_dlmort.png", "Test ADF: D(LMORT)"))
for (e_ in ev) { P <- pruebas[[as.integer(e_[1])]]; guardar_salida(salida_adf(adf(P$x, 12, P$det, P$dum, P$v), e_[3]), e_[2]) }
escribir_valores("valores_03.tex")
