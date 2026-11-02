# =============================================================================
# a04_relacion.R -- Relacion entre las series, muestra 2000-2019 (mismas herramientas que viviendas_tasa)
#   1. Regresion en niveles y alarma de regresion espuria R^2 > d_W (Nota 6, Sec. 3.3)
#   2. Correlacion cruzada de las innovaciones de cada modelo univariado (Haugh 1976).
#      k > 0: el desempleo anticipa al trafico; k < 0: el trafico anticipa al desempleo. Como el
#      desempleo es un indicador rezagado del ciclo, es esperable ver correlacion en k < 0 sin que
#      eso sea causalidad del trafico al desempleo: ambas responden al ciclo.
#   3. Rezagos distribuidos con F(B) = (1-B)(1-B^12) de ambos lados (la respuesta v(B) no cambia),
#      con las intervenciones del 11-S y el bisiesto filtradas igual.
#   4. Subperiodos 2001-2009 vs. 2010-2019.
# Ejecutar desde aereo_desempleo/.
# =============================================================================
source("scripts/utils_a.R")
D <- cargar_base(); G <- readRDS("output/tablas/a03_ganadores.rds")
rpm <- window(D$rpm, end = c(2019, 12)); des <- window(D$des, end = c(2019, 12))
pulso <- function(x, a, b) ts(as.numeric(abs(time(x) - fecha_num(a, b)) < 1e-6), start = start(x), frequency = 12)
IND <- cbind(AO2001M09 = pulso(rpm, 2001, 9), AO2001M10 = pulso(rpm, 2001, 10), AO2001M11 = pulso(rpm, 2001, 11),
             AO2001M12 = pulso(rpm, 2001, 12), BISIESTO = bisiesto(rpm))

# ---- 1. Niveles --------------------------------------------------------------------------------------
reg_alarma <- function(y, X, nombre) {
  m <- lm(y ~ X); s <- summary(m); e <- residuals(m); dw <- sum(diff(e)^2) / sum(e^2)
  data.frame(regresion = nombre, n = length(e), pendiente_desempleo = coef(m)[2], t = s$coefficients[2, 3],
             R2 = s$r.squared, DW = dw, alarma_R2_mayor_DW = s$r.squared > dw, row.names = NULL)
}
F <- function(z) diff(diff(z, 12))
AL <- rbind(reg_alarma(as.numeric(rpm), as.numeric(des), "niveles: ln RPM ~ desempleo"),
            reg_alarma(as.numeric(rpm), cbind(as.numeric(des), dummies_est(rpm), IND), "niveles + meses + 11-S + bisiesto"),
            reg_alarma(as.numeric(F(rpm)), cbind(as.numeric(F(des)), F(IND)), "filtrada: (1-B)(1-B^12) en ambas"))
cat("\n=== 1. Regresion estatica y alarma de regresion espuria (R^2 > DW) ===\n")
num <- sapply(AL, is.numeric); P <- AL; P[num] <- lapply(AL[num], round, 3); print(P, row.names = FALSE)
write.csv(AL, "output/tablas/a04_alarma.csv", row.names = FALSE)

# ---- 2. Correlacion cruzada de innovaciones -------------------------------------------------------------
inno <- function(Gi) ts(as.numeric(residuals(Gi$m)), start = start(Gi$b$W), frequency = 12)
Z <- ts.intersect(er = inno(G[["RPM | 2000-2019 + 11-S + bisiesto"]]), ed = inno(G[["Desempleo | 2000-2019"]]))
ccf_tab <- function(desde, hasta, nombre) {
  z <- window(Z, start = desde, end = hasta); n <- nrow(z)
  r <- ccf(as.numeric(z[, "ed"]), as.numeric(z[, "er"]), lag.max = 24, plot = FALSE)    # r(k) = corr(ed_t, er_{t+k})
  data.frame(periodo = nombre, k = as.numeric(r$lag), r = as.numeric(r$acf), banda = 1.96 / sqrt(n), n = n)
}
CC <- rbind(ccf_tab(start(Z), end(Z), "2001-2019"), ccf_tab(start(Z), c(2009, 12), "2001-2009"), ccf_tab(c(2010, 1), end(Z), "2010-2019"))
CC$sig <- abs(CC$r) > CC$banda
write.csv(CC, "output/tablas/a04_correlacion_cruzada.csv", row.names = FALSE)
cat("\n=== 2. Correlaciones cruzadas significativas: corr(innovacion desempleo_t, innovacion RPM_{t+k}) ===\n")
print(transform(CC[CC$sig, c("periodo", "k", "r", "banda")], r = round(r, 3), banda = round(banda, 3)), row.names = FALSE)
haugh <- function(d, M) { d <- d[abs(d$k) <= M, ]; S <- d$n[1] * sum(d$r^2); c(S = S, gl = 2 * M + 1, p = pchisq(S, 2 * M + 1, lower.tail = FALSE)) }
HA <- do.call(rbind, lapply(split(CC, CC$periodo), function(d) data.frame(periodo = d$periodo[1], rbind(c(M = 6, haugh(d, 6)), c(M = 12, haugh(d, 12))))))
cat("\n=== Test de Haugh de independencia entre las innovaciones ===\n"); print(transform(HA, S = round(S, 1), p = signif(p, 3)), row.names = FALSE)
write.csv(HA, "output/tablas/a04_haugh.csv", row.names = FALSE)

# ---- 3. Rezagos distribuidos con el mismo filtro -------------------------------------------------------------
rezagos_dist <- function(desde, hasta, nombre, K = 12) {
  y <- F(rpm); x <- F(des)
  X <- do.call(cbind, lapply(0:K, function(k) stats::lag(x, -k))); colnames(X) <- paste0("x_", 0:K)
  A <- window(ts.intersect(y, X, F(IND)), start = desde, end = hasta)
  yy <- A[, 1]; XX <- A[, 2:(K + 2)]; II <- A[, -(1:(K + 2)), drop = FALSE]; II <- II[, colSums(II != 0) > 0, drop = FALSE]
  m <- lm(yy ~ XX + II); s <- summary(m)$coefficients; e <- residuals(m); V <- vcov(m)
  data.frame(periodo = nombre, k = 0:K, v = s[2:(K + 2), 1], ee = s[2:(K + 2), 2], acumulada = cumsum(s[2:(K + 2), 1]),
             ee_acumulada = sapply(0:K, function(j) sqrt(sum(V[2:(j + 2), 2:(j + 2)]))),
             DW = sum(diff(e)^2) / sum(e^2), R2 = summary(m)$r.squared, n = length(e), row.names = NULL)
}
RD <- rbind(rezagos_dist(c(2002, 2), c(2019, 12), "2002-2019"), rezagos_dist(c(2002, 2), c(2009, 12), "2002-2009"),
            rezagos_dist(c(2010, 1), c(2019, 12), "2010-2019"))
write.csv(RD, "output/tablas/a04_rezagos_distribuidos.csv", row.names = FALSE)
cat("\n=== 3. Respuesta de 100 ln(RPM) a +1 pp de desempleo: v_k y respuesta acumulada (MCO, ee orientativos) ===\n")
print(transform(RD[RD$k %in% c(0, 1, 2, 3, 6, 9, 12), c("periodo", "k", "v", "ee", "acumulada", "ee_acumulada", "DW", "R2")],
                v = round(v, 2), ee = round(ee, 2), acumulada = round(acumulada, 2), ee_acumulada = round(ee_acumulada, 2), DW = round(DW, 2), R2 = round(R2, 3)), row.names = FALSE)

# ---- Figura -----------------------------------------------------------------------------------------------------
CCp <- CC; CCp$periodo <- factor(CCp$periodo, levels = c("2001-2019", "2001-2009", "2010-2019"))
g1 <- ggplot(CCp, aes(k, r)) +
  geom_rect(data = unique(CCp[, c("periodo", "banda")]), aes(xmin = -24.5, xmax = 24.5, ymin = -banda, ymax = banda), inherit.aes = FALSE, fill = gris, alpha = 0.45) +
  geom_vline(xintercept = 0, colour = tinta2, linewidth = 0.3) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.1, lineend = "round") +
  scale_colour_manual(values = c(`TRUE` = col_des, `FALSE` = "#9AA5B1"), guide = "none") +
  facet_wrap(~periodo, ncol = 1) + scale_x_continuous(breaks = seq(-24, 24, 6)) +
  labs(title = "Correlacion cruzada de las innovaciones: corr(desempleo_t, RPM_t+k)",
       subtitle = "k > 0: el desempleo anticipa al trafico; k < 0: el trafico anticipa al desempleo. Banda +-1.96/raiz(n)",
       x = "k (meses)", y = NULL) + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
RDp <- RD; RDp$periodo <- factor(RDp$periodo, levels = c("2002-2019", "2002-2009", "2010-2019"))
g2 <- ggplot(RDp, aes(k, acumulada, colour = periodo)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_ribbon(data = RDp[RDp$periodo == "2002-2019", ], aes(ymin = acumulada - 1.96 * ee_acumulada, ymax = acumulada + 1.96 * ee_acumulada),
              fill = tinta, alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(tinta, "#9AA5B1", col_rpm), name = NULL) + scale_x_continuous(breaks = 0:12) +
  labs(title = "Respuesta acumulada del trafico a +1 pp de desempleo",
       subtitle = "Suma de v_0..v_k, rezagos distribuidos en (1-B)(1-B^12), en %.\nSombra: +-1.96 ee de la suma, 2002-2019 (MCO; ee orientativos)",
       x = "Meses despues de la suba", y = "%") + tema + theme(legend.position = "top")
ggsave("output/figuras/a04_relacion.png", (g1 | g2) + plot_layout(widths = c(1.1, 1)) + plot_annotation(theme = tema_global), width = 13, height = 7.5, dpi = 200)
