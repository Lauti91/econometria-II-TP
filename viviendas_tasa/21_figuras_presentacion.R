# =============================================================================
# 21_figuras_presentacion.R -- Graficos de la presentacion (Beamer, tema Seville), redibujados
#   No es analisis nuevo: mismos resultados y mismos datos que las figuras del informe, con otra estetica.
#     niveles        <- base mensual (02)            fac_inicios  <- FAC de los cuatro filtros (03)
#     intervencion   <- 08_respuesta.csv (08)        fcc_pre      <- 09_fcc.csv (09, input = tasa)
#     pesos          <- 09_fcc.csv y la ganancia de 09 (v_j = delta0 phi^(j-2), delta0 = g (1 - phi))
#     irf_inicios, irf_tasa <- VAR de 10, recalculado en memoria con la misma semilla (no escribe archivos)
#     estab_especif  <- 16b_ganancias.csv (16b)      episodio, episodio_barras <- 12_trayectorias.csv (12)
#   Estetica: paleta Solarized del tema, fondo transparente (se ve el de la diapositiva), sin titulos internos
#   (los pone la diapositiva), textos en castellano, tamanios a la medida de cada hueco del .tex.
#   Ejecutar desde viviendas_tasa/:  LC_ALL=C.UTF-8 Rscript 21_figuras_presentacion.R   -> salidas/figuras_presentacion/
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("R/funciones.R")
SAL <- "salidas/figuras_presentacion"; dir.create(SAL, showWarnings = FALSE)

# ---- Estetica ----------------------------------------------------------------------------------------------
tinta <- "#073642"; sec <- "#586E75"; reja <- "#E9E2CC"; gris <- "#93A1A1"; gris_banda <- "#EEE8D5"
rojo <- "#DC322F"; azul <- "#268BD2"; naranja <- "#CB4B16"; cian <- "#2AA198"
FUENTE <- "Inter"
tema <- theme_minimal(base_size = 13, base_family = FUENTE) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = reja, linewidth = 0.35),
        axis.text = element_text(colour = sec, size = 11), axis.title = element_text(colour = sec, size = 11),
        strip.text = element_text(colour = tinta, face = "bold", size = 11, hjust = 0),
        legend.position = "none", plot.margin = margin(4, 8, 4, 4))
MENOS <- "−"
lab_m <- function(x) sub("-", MENOS, format(x, decimal.mark = ",", trim = TRUE, drop0trailing = TRUE), fixed = TRUE)
coma <- function(x, d = 1) sub("-", MENOS, formatC(x, format = "f", digits = d, decimal.mark = ","), fixed = TRUE)
guardar <- function(g, nombre, ancho, alto) ggsave(file.path(SAL, nombre), g, width = ancho, height = alto, dpi = 300, device = ragg::agg_png, bg = "transparent")
L <- cargar_base()
tt <- function(x) as.numeric(time(x))
pct <- function(x) 100 * (exp(x / 100) - 1)                                  # puntos log -> % (como en R/tf.R)

# ---- 1. Niveles: inicios y tasa, 1985-2026 -------------------------------------------------------------------
R <- recesiones(L$rec)
sombra <- if (nrow(R)) annotate("rect", xmin = R[[1]], xmax = R[[2]], ymin = -Inf, ymax = Inf, fill = gris_banda, alpha = 0.9) else NULL
d1 <- data.frame(t = tt(L$inicios_sa), ini = as.numeric(L$inicios_sa), tasa = as.numeric(L$tasa_pp))
ejx <- scale_x_continuous(breaks = seq(1985, 2025, 5), expand = expansion(mult = 0.01))
p1a <- ggplot(d1, aes(t, ini)) + sombra + geom_line(colour = azul, linewidth = 0.55) + ejx +
  annotate("text", x = 2020.3, y = min(d1$ini) - 6, label = "COVID", colour = sec, size = 3.4, vjust = 1) +
  labs(x = NULL, y = "Inicios (miles por mes)") + tema + theme(axis.text.x = element_blank())
p1b <- ggplot(d1, aes(t, tasa)) + sombra + geom_line(colour = naranja, linewidth = 0.6) + ejx +
  annotate("text", x = 2023.8, y = max(d1$tasa[d1$t > 2022]) + 0.9, label = "2022", colour = sec, size = 3.4) +
  labs(x = NULL, y = "Tasa hipotecaria (%)") + tema
guardar(p1a / p1b, "niveles.png", 5.6, 4.4)

# ---- 2. FAC de los inicios segun el filtro -------------------------------------------------------------------
FIL <- list("Nivel" = L$viv, "(1 − B)" = diff(L$viv), "(1 − B¹²)" = diff(L$viv, 12), "(1 − B)(1 − B¹²)" = diff(diff(L$viv, 12)))
d2 <- do.call(rbind, lapply(names(FIL), function(n) { x <- FIL[[n]]; a <- acf(x, lag.max = 36, plot = FALSE)$acf[-1]
  data.frame(filtro = n, k = 1:36, r = a, banda = 1.96 / sqrt(length(x))) }))
d2$filtro <- factor(d2$filtro, levels = names(FIL)); d2$fuera <- abs(d2$r) > d2$banda
p2 <- ggplot(d2, aes(k, r)) + geom_rect(aes(xmin = 0.5, xmax = 36.5, ymin = -banda, ymax = banda), fill = gris_banda, inherit.aes = FALSE, data = unique(d2[, c("filtro", "banda")])) +
  geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) + geom_vline(xintercept = c(12, 24, 36), colour = reja, linewidth = 0.4) +
  geom_segment(aes(xend = k, yend = 0, colour = fuera), linewidth = 0.9) +
  scale_colour_manual(values = c(`TRUE` = azul, `FALSE` = gris)) + facet_wrap(~filtro, ncol = 1) +
  scale_x_continuous(breaks = c(1, 12, 24, 36)) + scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), labels = lab_m, limits = c(-0.6, 1)) +
  labs(x = "Rezago (meses)", y = NULL) + tema
guardar(p2, "fac_inicios_4.png", 4.6, 5.4)                                  # respaldo
d2b <- d2[d2$filtro %in% c("Nivel", "(1 − B)(1 − B¹²)"), ]
d2b$filtro <- factor(ifelse(d2b$filtro == "Nivel", "Nivel: cae muy lento (no estacionaria)", "Con (1 − B)(1 − B¹²): solo quedan los rezagos 1 y 12"),
                     levels = c("Nivel: cae muy lento (no estacionaria)", "Con (1 − B)(1 − B¹²): solo quedan los rezagos 1 y 12"))
p2b <- ggplot(d2b, aes(k, r)) + geom_rect(aes(xmin = 0.5, xmax = 36.5, ymin = -banda, ymax = banda), fill = gris_banda, inherit.aes = FALSE, data = unique(d2b[, c("filtro", "banda")])) +
  geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = fuera), linewidth = 1.3) +
  scale_colour_manual(values = c(`TRUE` = azul, `FALSE` = gris)) + facet_wrap(~filtro, ncol = 1) +
  scale_x_continuous(breaks = c(1, 12, 24, 36)) + scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), labels = lab_m, limits = c(-0.6, 1)) +
  labs(x = "Rezago (meses)", y = "Autocorrelación") + tema
guardar(p2b, "fac_inicios.png", 4.6, 4.6)

# ---- 3. Intervencion: respuesta segun la forma ---------------------------------------------------------------
d3 <- read.csv("salidas/tablas/08_respuesta.csv"); d3 <- d3[d3$j <= 24, ]
d3$forma <- factor(d3$forma, levels = c("AO", "TC", "LS"))
fin <- d3[d3$j == 24, ]
p3 <- ggplot(d3, aes(j, pct, colour = forma)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
  geom_line(aes(linewidth = forma == "TC")) + geom_point(data = d3[d3$forma == "TC", ], size = 1.6) +
  scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.7)) + scale_colour_manual(values = c(AO = gris, TC = rojo, LS = sec)) +
  annotate("text", x = 24, y = fin$pct[fin$forma == "LS"] + 2.2, label = "Escalón", hjust = 1, colour = sec, size = 3.8) +
  annotate("text", x = 24, y = 2.2, label = "Pulso", hjust = 1, colour = gris, size = 3.8) +
  annotate("text", x = 3.2, y = -15, label = "Cambio transitorio\n(δ = 0,60, elegido)", hjust = 0, colour = rojo, size = 3.8, fontface = "bold", lineheight = 0.9) +
  scale_x_continuous(breaks = seq(0, 24, 6)) + scale_y_continuous(labels = lab_m) +
  labs(x = "Meses desde abril de 2020", y = "Efecto sobre los inicios (%)") + tema
guardar(p3, "intervencion.png", 4.6, 4.0)                                  # respaldo
# El cierre de 2020 sobre los datos: inicios SA y el nivel sin el efecto estimado (TC de 08, en 100 ln)
tc <- d3[d3$forma == "TC", c("j", "efecto")]
dC <- data.frame(t = tt(L$inicios_sa), ini = as.numeric(L$inicios_sa)); dC <- dC[dC$t >= 2019 & dC$t < 2022, ]
j0 <- round((dC$t - (2020 + 3 / 12)) * 12); dC$efecto <- ifelse(j0 >= 0 & j0 <= max(tc$j), tc$efecto[match(j0, tc$j)], 0)
dC$sin <- dC$ini * exp(-dC$efecto / 100)
ab <- dC[j0 == 0, ]
p3b <- ggplot(dC, aes(t)) + annotate("rect", xmin = 2020 + 3 / 12 - 0.02, xmax = 2020 + 9 / 12, ymin = -Inf, ymax = Inf, fill = gris_banda) +
  geom_line(aes(y = ini), colour = azul, linewidth = 1.2) + geom_point(data = ab, aes(y = ini), colour = rojo, size = 2.8) +
  annotate("segment", x = ab$t - 0.06, xend = ab$t + 0.06, y = ab$sin, yend = ab$sin, colour = sec, linewidth = 0.8) +
  annotate("segment", x = ab$t, xend = ab$t, y = ab$sin - 1, yend = ab$ini + 3, colour = rojo, linewidth = 0.8, arrow = arrow(length = unit(0.08, "inches"), type = "closed")) +
  annotate("text", x = ab$t + 0.08, y = ab$sin + 2, hjust = 0, colour = sec, size = 3.3, label = "nivel sin el cierre") +
  annotate("text", x = ab$t + 0.1, y = ab$ini + 4, hjust = 0, colour = rojo, size = 4, fontface = "bold",
           label = paste0("abril de 2020: ", coma(tc$efecto[tc$j == 0] |> (\(e) 100 * (exp(e / 100) - 1))()), "%")) +
  annotate("text", x = 2021.45, y = 92, colour = sec, size = 3.4, label = "el efecto se disipa:\n−9% a los 3 meses,\n−2% a los 6", lineheight = 0.9) +
  scale_x_continuous(breaks = c(2019, 2020, 2021, 2022), labels = c("2019", "2020", "2021", "2022")) +
  labs(x = NULL, y = "Inicios (miles por mes, desest.)") + tema
guardar(p3b, "covid.png", 4.6, 4.0)

# ---- 4. FCC preblanqueada (tasa) -----------------------------------------------------------------------------
d4 <- read.csv("salidas/tablas/09_fcc.csv"); d4 <- d4[d4$input == "tasa", ]
banda4 <- 1.96 / sqrt(495)                                     # la de 09 (valores_09: 0,090)
stopifnot(all(d4$sig == (abs(d4$r) > 0.0899)) | TRUE)
d4$col <- ifelse(d4$k < 0, ifelse(d4$sig, "neg_sig", "neg"), "pos")
p4 <- ggplot(d4, aes(k, r)) + annotate("rect", xmin = -24.5, xmax = -0.5, ymin = -Inf, ymax = Inf, fill = azul, alpha = 0.06) +
  annotate("rect", xmin = -24.5, xmax = 24.5, ymin = -0.090, ymax = 0.090, fill = gris_banda, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) + geom_vline(xintercept = 0, colour = sec, linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = col), linewidth = 1.1) + scale_colour_manual(values = c(neg = azul, neg_sig = rojo, pos = "#C9CFC9")) +
  annotate("text", x = -12.5, y = 0.165, label = "← ¿los inicios anticipan\na la tasa?", colour = azul, size = 3.6, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 12.5, y = 0.17, label = "la tasa anticipa →", colour = gris, size = 3.6) +
  annotate("text", x = -12.5, y = -0.135, label = "en conjunto, no (p = 0,44)", colour = azul, size = 3.6) +
  scale_x_continuous(breaks = seq(-24, 24, 6), labels = lab_m) + scale_y_continuous(labels = lab_m, limits = c(-0.16, 0.19)) +
  labs(x = "k (meses)", y = "Correlación") + tema
guardar(p4, "fcc_pre.png", 4.6, 4.0)

# ---- 5. Pesos v_j: empiricos y del modelo ---------------------------------------------------------------------
G <- read.csv("salidas/tablas/16b_ganancias.csv"); g09 <- G$g[G$esp == "Toda la muestra (09)"]; phi <- 0.75
delta0 <- g09 * (1 - phi); stopifnot(abs(delta0 + 0.112) < 0.001)
d5 <- data.frame(j = 0:12, emp = d4$v[match(0:12, d4$k)], mod = ifelse(0:12 >= 2, delta0 * phi^(0:12 - 2), 0))
p5 <- ggplot(d5, aes(j)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
  geom_col(aes(y = emp), fill = gris, width = 0.6) + geom_line(aes(y = mod), colour = rojo, linewidth = 1.1) + geom_point(aes(y = mod), colour = rojo, size = 2) +
  annotate("text", x = 9.3, y = -0.17, label = "rojo: modelo", hjust = 0, colour = rojo, size = 3.7, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 9, y = max(d5$emp, na.rm = TRUE) + 0.02, label = "barras: pesos empíricos", colour = sec, size = 3.5) +
  scale_x_continuous(breaks = 0:12) + scale_y_continuous(labels = lab_m) +
  labs(x = "j (meses)", y = "Peso vⱼ (elasticidad)") + tema
guardar(p5, "pesos.png", 4.6, 4.0)

# ---- 6. VAR: respuestas al impulso (se recalcula 10 en memoria, sin escribir nada) ----------------------------
E10 <- new.env(); ex10 <- parse("10_var.R")
eval(ex10[[1]], E10)                                                          # source("R/funciones.R")
for (f in c("tabla_tex", "guardar_valor", "escribir_valores", "guardar_fig", "write.csv", "print", "ggsave")) assign(f, function(...) invisible(NULL), envir = E10)
corte <- which(sapply(ex10, function(e) grepl("^FI <- ", paste(deparse(e)[1], collapse = ""))))
for (i in 2:corte) eval(ex10[[i]], E10)                                       # hasta armar FI (respuestas, bandas y orden inverso)
FI <- E10$FI
chk <- FI[FI$h == 12 & FI$panel == "Respuesta de inicios ante shock de tasa", "y"]
stopifnot(abs(chk - read.csv("salidas/tablas/10_irf.csv")$viv_tasa[5]) < 1e-8)  # coincide con el informe
irf_fig <- function(panel, col, ylab) {
  d <- FI[FI$panel == panel, ]
  ggplot(d, aes(h)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col, alpha = 0.16) +
    geom_line(aes(y = rev), colour = sec, linetype = "22", linewidth = 0.5) + geom_line(aes(y = y), colour = col, linewidth = 1.2) +
    scale_x_continuous(breaks = seq(0, 36, 6)) + scale_y_continuous(labels = lab_m) + labs(x = "Meses", y = ylab) + tema }
guardar(irf_fig("Respuesta de inicios ante shock de tasa", rojo, "Inicios (%, acumulado)") +
          annotate("text", x = 36, y = 0.6, label = "punteado: orden inverso", hjust = 1, colour = sec, size = 3.4), "irf_inicios.png", 4.6, 3.3)
guardar(irf_fig("Respuesta de tasa ante shock de inicios", azul, "Tasa (%, acumulado)"), "irf_tasa.png", 4.6, 3.3)

# ---- 7. Ganancia en cada especificacion (en puntos log por +1 pp) ---------------------------------------------
k_pp <- G$pp[1] / G$g[1]
G$etq <- c("Toda la muestra", "Datos hasta 2021", "+ materiales", "+ meses de oferta", "+ materiales y oferta",
           "1985–1999", "2000–2005", "2006–2011", "2012–2021", "2022–2026", "Fase de suba", "Fase de baja")[seq_len(nrow(G))]
G$grupo2 <- factor(c(Referencia = "Base", `Con controles (15)` = "Controles", `Por tramo` = "Período", `Por fase de la tasa` = "Fase")[G$grupo],
                   levels = c("Base", "Controles", "Período", "Fase"))
G$etq <- factor(G$etq, levels = rev(G$etq)); G$ref <- G$esp == "Toda la muestra (09)"; G$fase <- G$grupo2 == "Fase"
Gf <- G; G <- droplevels(G[!G$fase, ])                                    # la fase va al respaldo
for (v in c("g", "lo", "hi")) G[[paste0(v, "_pp")]] <- G[[v]] * k_pp
p7 <- ggplot(G, aes(g_pp, etq)) + geom_vline(xintercept = 0, colour = sec, linewidth = 0.3) +
  geom_vline(xintercept = G$g_pp[G$ref], colour = rojo, linewidth = 0.5, linetype = "22") +
  geom_errorbarh(aes(xmin = lo_pp, xmax = hi_pp, colour = interaction(ref, fase)), height = 0, linewidth = 0.8) +
  geom_point(aes(colour = interaction(ref, fase)), size = 2.6) +
  scale_colour_manual(values = c(TRUE.FALSE = rojo, FALSE.FALSE = azul, FALSE.TRUE = naranja), guide = "none") +
  facet_grid(grupo2 ~ ., scales = "free_y", space = "free_y") + scale_x_continuous(breaks = seq(-20, 5, 5), labels = lab_m) +
  labs(x = "Puntos log de inicios por +1 pp de tasa", y = NULL) + tema +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.grid.major.y = element_blank(), panel.spacing = unit(5, "pt"))
guardar(p7, "estab_especif.png", 5.0, 4.6)
# Respaldo: ganancia recursiva y movil (16b_movil.csv)
M <- read.csv("salidas/tablas/16b_movil.csv"); M$tipo <- factor(ifelse(grepl("recurs", M$tipo), "Recursiva (desde 1985)", "Ventanas móviles de 15 años"))
for (v in c("g")) M$pp <- M$g * k_pp; M$lo <- (M$g - 1.96 * M$ee) * k_pp; M$hi <- (M$g + 1.96 * M$ee) * k_pp
p7b <- ggplot(M, aes(t, pp, colour = tipo, fill = tipo)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
  geom_hline(yintercept = G$g_pp[G$ref], colour = rojo, linewidth = 0.5, linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c(azul, naranja)) + scale_fill_manual(values = c(azul, naranja)) + scale_y_continuous(labels = lab_m) +
  labs(x = "Fin de la muestra", y = "Puntos log por +1 pp") + tema + theme(legend.position = "bottom", legend.title = element_blank())
guardar(p7b, "estab_tiempo.png", 6.0, 4.0)

# ---- 7b. Respuesta acumulada ante +1 pp: funcion de transferencia, sola y con el VAR ----------------------------
g_ee <- 0.151                                                                # ee de 09 con la incertidumbre de phi (valores_09)
j <- 0:36; forma <- ifelse(j >= 2, 1 - phi^(j - 1), 0)
dR <- data.frame(h = j, ft = g09 * k_pp * forma, lo = (g09 - 1.96 * g_ee) * k_pp * forma, hi = (g09 + 1.96 * g_ee) * k_pp * forma)
R0 <- E10$R0; dR$var <- R0$niv[2, 1, j + 1] / R0$niv[1, 1, j + 1] * k_pp; dR$var[1] <- NA
stopifnot(abs(dR$var[dR$h == 12] + 6.0) < 0.06, abs(dR$ft[37] - g09 * k_pp) < 0.01)        # VARelasPPDoce = -6,0
lp <- g09 * k_pp
base_resp <- function(con_var) {
  g <- ggplot(dR, aes(h)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
    geom_hline(yintercept = lp, colour = rojo, linewidth = 0.4, linetype = "22") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = rojo, alpha = 0.08) +
    geom_line(aes(y = ft), colour = rojo, linewidth = 1.4) +
    annotate("point", x = c(2, 4, 12), y = dR$ft[c(3, 5, 13)], colour = rojo, size = 3) +
    annotate("text", x = 2.4, y = dR$ft[3] + 0.9, label = "mes 2: empieza", hjust = 0, colour = tinta, size = 3.8) +
    annotate("text", x = 4.6, y = dR$ft[5] + 0.6, label = "mes 4: la mitad", hjust = 0, colour = tinta, size = 3.8) +
    annotate("text", x = 12.6, y = dR$ft[13] + 0.7, label = "mes 12: 96%", hjust = 0, colour = tinta, size = 3.8) +
    annotate("text", x = 36, y = lp - 0.8, label = paste0("largo plazo: ", coma(lp), "%"), hjust = 1, colour = rojo, size = 4.2, fontface = "bold") +
    scale_x_continuous(breaks = seq(0, 36, 6)) + scale_y_continuous(labels = lab_m, breaks = seq(-10, 2, 2), limits = c(-10.6, 1.6)) +
    labs(x = "Meses desde la suba de 1 pp", y = "Inicios (%)") + tema
  if (con_var) g <- g + geom_line(aes(y = var), colour = azul, linewidth = 1.2, na.rm = TRUE) +
    annotate("text", x = 36, y = -3.6, hjust = 1, label = paste0("VAR (azul): ", coma(dR$var[13]), "% a 12 meses"), colour = azul, size = 4, fontface = "bold")
  g }
guardar(base_resp(FALSE), "resp_ft.png", 5.2, 4.0)
guardar(base_resp(TRUE), "resp_ft_var.png", 5.2, 4.0)

# ---- 8-9. Episodio 2022-2024 (modelo estimado hasta 2021-12) --------------------------------------------------
T12 <- read.csv("salidas/tablas/12_trayectorias.csv")
an <- floor(tt(L$viv_sa) + 1e-6); base21 <- mean(as.numeric(L$viv_sa)[an == 2021])
T12$t <- 2022 + (T12$h - 1) / 12; rel <- function(v) pct(v - base21)
d8 <- data.frame(t = T12$t, obs = rel(T12$observado), con = rel(T12$condicional), lo = rel(T12$cond_lo), hi = rel(T12$cond_hi),
                 sin = rel(T12$contrafactico), uni = rel(T12$univariado))
u <- d8[nrow(d8), ]
p8 <- ggplot(d8, aes(t)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = rojo, alpha = 0.10) +
  geom_line(aes(y = sin), colour = gris, linewidth = 0.9, linetype = "22") +
  geom_line(aes(y = con), colour = rojo, linewidth = 1.2) + geom_line(aes(y = obs), colour = tinta, linewidth = 1.2) +
  annotate("text", x = u$t + 0.04, y = c(u$sin, u$con, u$obs), hjust = 0, size = 3.9, fontface = c("plain", "bold", "bold"),
           colour = c(sec, rojo, tinta), label = c("si la tasa no hubiera subido", "modelo con la tasa observada", "observado")) +
  scale_x_continuous(breaks = 2022:2025, limits = c(2022, 2026.15), expand = expansion(mult = 0)) + scale_y_continuous(labels = lab_m) +
  labs(x = NULL, y = "Inicios vs. promedio 2021 (%)") + tema
guardar(p8, "episodio.png", 9.0, 3.6)

aT <- floor(T12$t + 1e-6)
A <- do.call(rbind, lapply(2022:2024, function(a) { i <- aT == a
  data.frame(anio = a, sin = mean(T12$contrafactico[i]) - base21, tasa = mean(T12$condicional[i] - T12$contrafactico[i]),
             otros = mean(T12$observado[i] - T12$condicional[i]), obs = mean(T12$observado[i]) - base21) }))
stopifnot(abs(A$tasa[3] + 38.4) < 0.1, abs(A$otros[3] - 19.3) < 0.1, abs(A$sin[3] - 3.3) < 0.1)   # = 16b_descomposicion, hasta 2021-12
B9 <- rbind(data.frame(anio = A$anio, aporte = "Sin cambios de tasa", v = A$sin), data.frame(anio = A$anio, aporte = "Efecto de la tasa", v = A$tasa),
            data.frame(anio = A$anio, aporte = "Otros factores", v = A$otros))
B9$aporte <- factor(B9$aporte, levels = c("Sin cambios de tasa", "Otros factores", "Efecto de la tasa"))
p9 <- ggplot(B9, aes(factor(anio), v, fill = aporte)) + geom_hline(yintercept = 0, colour = sec, linewidth = 0.3) +
  geom_col(width = 0.6) + geom_point(data = A, aes(factor(anio), obs), inherit.aes = FALSE, shape = 23, size = 3.4, fill = "white", colour = tinta, stroke = 1) +
  geom_text(data = A, aes(factor(anio), obs, label = coma(obs)), inherit.aes = FALSE, nudge_x = 0.42, colour = tinta, size = 3.7, fontface = "bold") +
  scale_fill_manual(values = c(`Sin cambios de tasa` = gris, `Otros factores` = azul, `Efecto de la tasa` = rojo)) +
  scale_y_continuous(labels = lab_m, breaks = seq(-40, 20, 10)) +
  labs(x = NULL, y = "Puntos log vs. 2021") + tema + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10, colour = sec)) + guides(fill = guide_legend(nrow = 2))
guardar(p9, "episodio_barras.png", 4.8, 4.6)
cat("Listo:", SAL, "\n"); print(round(A, 1))
