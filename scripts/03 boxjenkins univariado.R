# =============================================================================
# 03_boxjenkins_univariado.R  --  Ciclo Box-Jenkins sobre cada serie (identificacion,
# estimacion, diagnostico)  |  TP Econometria II, candidata 5: Indice Li Keqiang
#
# Series: PBI oficial y electricidad, crecimiento interanual (%), trimestral.
# Dos muestras, para leer el efecto del COVID (regla 7 y 9 de la clase del 14/10):
#   - completa  2000Q1-2023Q4 (96 obs)
#   - sin COVID 2000Q1-2019Q4 (80 obs)
# (Todavia sin variables de intervencion: eso es de la segunda mitad del curso.)
#
# Respaldo teorico (Notas propias):
#   - Ciclo identificacion-estimacion-diagnostico-pronostico ........ Nota 6, Sec. 6
#   - Parsimonia, factor comun ...................................... Nota 6, Sec. 7; Nota 4, Sec. 3.2
#   - AIC / SBC, y preferencia por SBC .............................. Nota 6, Sec. 8
#   - FAC/FACP, bandas de Bartlett y Quenouille; Ljung-Box .......... Nota 6, Sec. 9
#   - Diagnostico: residuos, correlograma de residuos, sobreajuste ... Nota 6, Sec. 10
#   - Condiciones de estacionariedad e invertibilidad (raices) ....... Nota 3, Sec. 2; Nota 4, Sec. 2.4 y 3.2
#   - Ljung-Box con gl = K - p - q (K - parametros ARMA) ............. Nota 4, Sec. 3.6
#   - Jarque-Bera con (n - p - q)/6 y grafico Q-Q .................... Nota 3, Sec. 6.2; Nota 4, Sec. 3.6
#   - Residuos al cuadrado (dependencia no lineal) ................... Nota 3, Sec. 6.3
#   - Estimacion por maxima verosimilitud (no lineal si hay MA) ...... Nota 4, Sec. 4
#   - Si dos modelos son comparables (dif. de SBC < 2): no hay ganador contundente ... Nota 6, Sec. 11.5
# Clase de Mirta (14/10): los rezagos estacionales se modelan con terminos MA; un correlograma
# grueso o un AR de orden 3-4 obliga a dudar de la estacionariedad y a mostrar las raices.
#
# Termino estacional: el crecimiento interanual es (aprox.) la diferencia estacional (1-B^4) de
# los log-niveles (Nota 6, Sec. 2.3), asi que deja una estructura de MA en el rezago 4. Se la
# modela con un MA estacional multiplicativo (1 - Theta B^4), es decir, un MA con coeficientes
# en los rezagos 4 (y 5 por el producto con la parte regular). Las Notas 1-6 no desarrollan los
# modelos estacionales multiplicativos como tales: se justifica por la regla 8 de Mirta y se
# estima con el mismo aparato de maxima verosimilitud de la Nota 4.
#
# Convencion de signos: Theta(B) = 1 - theta_1 B - ... (la de las Notas y de BJRL). R usa
# 1 + theta B, por lo que los coeficientes MA que se reportan aqui son los de R con el signo cambiado.
#
# Ejecutar desde la raiz del repositorio:  Rscript scripts/03_boxjenkins_univariado.R
# Paquetes: ggplot2, patchwork (se instalan solos si faltan)
# =============================================================================

source("scripts/utils.R")
dir.create("output/figuras", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tablas",  recursive = TRUE, showWarnings = FALSE)

series  <- cargar_yoy()
colores <- c(col_pbi, col_elec); names(colores) <- names(series)
muestras <- list("completa 2000Q1-2023Q4" = c(2000, 1, 2023, 4),
                 "sin COVID 2000Q1-2019Q4" = c(2000, 1, 2019, 4))

# Parametros del ejercicio
P_MAX <- 4          # rezago maximo AR: cota de Said-Dickey n^(1/3) ~ 4 (Nota 6, Sec. 4.5)
Q_MAX <- 2          # rezago maximo MA regular
ALFA  <- 0.05
K_LB  <- c(8, 12, 16)   # rezagos del Ljung-Box
UMBRAL_CANCEL <- 0.10   # criterio propio: raiz AR y raiz MA a menos de 0.10 => posible factor comun

# ---- 1. Utilidades ------------------------------------------------------------------
ajustar <- function(x, p, q, Qs) {
  tryCatch(suppressWarnings(
    arima(x, order = c(p, 0, q), seasonal = list(order = c(0, 0, Qs), period = 4), method = "ML")),
    error = function(e) NULL)
}

modulos <- function(m) {
  ar <- m$model$phi; ma <- m$model$theta     # polinomios ya expandidos (incluida la parte estacional)
  list(G = if (length(ar)) polyroot(c(1, -ar)) else complex(0),
       H = if (length(ma)) polyroot(c(1,  ma)) else complex(0))
}

jarque_bera <- function(e, k) {
  e <- e - mean(e); n <- length(e); s <- sqrt(mean(e^2))
  As <- mean(e^3) / s^3; Ku <- mean(e^4) / s^4
  JB <- (n - k) / 6 * (As^2 + (Ku - 3)^2 / 4)     # Nota 3, Sec. 6.2 (n-p) y Nota 4 (n-p-q)
  c(As = As, K = Ku, JB = JB, p = pchisq(JB, 2, lower.tail = FALSE))
}

evaluar <- function(m, x, p, q, Qs) {
  k <- p + q + Qs                      # parametros ARMA estimados (sin la media)
  e <- as.numeric(residuals(m)); rt <- modulos(m)
  # Factor comun (Nota 4, Sec. 3.2): se compara la parte AR con la parte MA *regular*. El MA
  # estacional no puede cancelarse con la parte AR porque no hay AR estacional en la rejilla.
  Hreg   <- if (q > 0) polyroot(c(1, m$coef[paste0("ma", 1:q)])) else complex(0)
  cancel <- if (length(rt$G) && length(Hreg)) min(abs(outer(rt$G, Hreg, "-"))) else Inf
  lb <- sapply(K_LB, function(K) if (K > k) Box.test(e, K, "Ljung-Box", fitdf = k)$p.value else NA_real_)
  jb <- jarque_bera(e, k)
  data.frame(
    p = p, q = q, Q_est = Qs, k = k, logLik = as.numeric(logLik(m)), AIC = AIC(m), SBC = BIC(m),
    sigma2 = m$sigma2,
    min_mod_G = if (length(rt$G)) min(Mod(rt$G)) else Inf,       # estacionariedad: > 1
    min_mod_H = if (length(rt$H)) min(Mod(rt$H)) else Inf,       # invertibilidad: > 1
    dist_AR_MA = cancel,
    LB8 = lb[1], LB12 = lb[2], LB16 = lb[3], JB_p = jb["p"], row.names = NULL)
}

# ---- 2. Rejilla de candidatos -----------------------------------------------------------
cand <- list()
for (nm in names(series)) for (mu in names(muestras)) {
  w <- muestras[[mu]]; x <- window(series[[nm]], start = w[1:2], end = w[3:4])
  for (p in 0:P_MAX) for (q in 0:Q_MAX) for (Qs in 0:1) {
    m <- ajustar(x, p, q, Qs); if (is.null(m)) next
    cand[[length(cand) + 1]] <- cbind(serie = nm, muestra = mu, n = length(x), evaluar(m, x, p, q, Qs))
  }
}
CAND <- do.call(rbind, cand)
CAND$estacionario <- CAND$min_mod_G > 1
CAND$invertible   <- CAND$min_mod_H > 1
CAND$sin_factor_comun <- CAND$dist_AR_MA > UMBRAL_CANCEL
CAND$LB_ok <- with(CAND, !is.na(LB8) & !is.na(LB12) & !is.na(LB16) & LB8 > ALFA & LB12 > ALFA & LB16 > ALFA)
CAND$admisible <- with(CAND, estacionario & invertible & sin_factor_comun & LB_ok)
CAND <- CAND[order(CAND$serie, CAND$muestra, CAND$SBC), ]
write.csv(CAND, "output/tablas/03_candidatos.csv", row.names = FALSE)

cat("\n=== Candidatos evaluados:", nrow(CAND), "| admisibles (estacionario, invertible, sin factor comun, Ljung-Box en 8/12/16 > 5%):",
    sum(CAND$admisible), "===\n")

# ---- 3. Seleccion por SBC entre los admisibles (Nota 6, Sec. 8.3) -----------------------
sel <- list(); top <- list()
for (nm in names(series)) for (mu in names(muestras)) {
  d <- CAND[CAND$serie == nm & CAND$muestra == mu & CAND$admisible, ]
  if (!nrow(d)) next
  d <- d[order(d$SBC), ]
  d$dSBC <- d$SBC - min(d$SBC)
  top[[length(top) + 1]] <- head(d, 5)
  sel[[length(sel) + 1]] <- d[1, ]
}
TOP <- do.call(rbind, top); SEL <- do.call(rbind, sel)
write.csv(TOP, "output/tablas/03_top_admisibles.csv", row.names = FALSE)

cat("\n=== Mejores admisibles por SBC (dSBC < 2 => evidencia debil a favor de un unico modelo, Nota 6 Sec. 11.5) ===\n")
print(TOP[, c("serie", "muestra", "p", "q", "Q_est", "AIC", "SBC", "dSBC", "LB12", "JB_p")], row.names = FALSE, digits = 4)

sin_modelo <- setdiff(paste(rep(names(series), each = length(muestras)), names(muestras)),
                      paste(SEL$serie, SEL$muestra))
if (length(sin_modelo)) cat("\nSIN modelo admisible en:", paste(sin_modelo, collapse = " | "), "\n")

# ---- 4. Detalle de los modelos elegidos: coeficientes, raices, diagnostico ---------------
coef_tab <- list(); diag_tab <- list(); sobre <- list()
png_paneles <- list()

for (i in seq_len(nrow(SEL))) {
  nm <- SEL$serie[i]; mu <- SEL$muestra[i]; p <- SEL$p[i]; q <- SEL$q[i]; Qs <- SEL$Q_est[i]
  w <- muestras[[mu]]; x <- window(series[[nm]], start = w[1:2], end = w[3:4])
  m <- ajustar(x, p, q, Qs); k <- p + q + Qs; n <- length(x)
  nombre_mod <- sprintf("ARMA(%d,%d)%s", p, q, if (Qs) " + MA estacional(4)" else "")
  
  # Coeficientes (convencion de las Notas: MA con signo cambiado)
  est <- m$coef; se <- sqrt(diag(m$var.coef))[names(est)]
  signo <- ifelse(grepl("^(ma|sma)", names(est)), -1, 1)
  coef_tab[[i]] <- data.frame(serie = nm, muestra = mu, modelo = nombre_mod, parametro = names(est),
                              estimacion = signo * est, ee = se, t = signo * est / se, row.names = NULL)
  
  # Raices (Mirta: mostrar las raices estimadas; Nota 3 Sec. 2)
  rt <- modulos(m)
  cat(sprintf("\n--- %s | %s | %s (n=%d) ---\n", nm, mu, nombre_mod, n))
  print(transform(coef_tab[[i]][, c("parametro", "estimacion", "ee", "t")],
                  estimacion = round(estimacion, 3), ee = round(ee, 3), t = round(t, 2)), row.names = FALSE)
  if (length(rt$G)) cat("Raices AR (G): modulos", paste(round(sort(Mod(rt$G)), 3), collapse = ", "),
                        "| menor modulo:", round(min(Mod(rt$G)), 3), "(estacionariedad exige > 1)\n")
  if (length(rt$H)) cat("Raices MA (H): modulos", paste(round(sort(Mod(rt$H)), 3), collapse = ", "),
                        "| menor modulo:", round(min(Mod(rt$H)), 3), "(invertibilidad exige > 1)\n")
  
  # Diagnostico de residuos (Nota 6 Sec. 10; Nota 3 Sec. 6)
  e  <- as.numeric(residuals(m)); ee_std <- e / sqrt(m$sigma2)
  jb <- jarque_bera(e, k); lb2 <- Box.test(e^2, 12, "Ljung-Box")
  atip <- which(abs(ee_std) > 3)
  fechas <- etiqueta_trim(x)
  diag_tab[[i]] <- data.frame(
    serie = nm, muestra = mu, modelo = nombre_mod, n = n, sigma = sqrt(m$sigma2),
    LB8_p = SEL$LB8[i], LB12_p = SEL$LB12[i], LB16_p = SEL$LB16[i],
    asimetria = jb["As"], curtosis = jb["K"], JB = jb["JB"], JB_p = jb["p"],
    LB_cuadrados12_p = lb2$p.value,
    n_atipicos_3sd = length(atip), atipicos = paste(fechas[atip], collapse = " "), row.names = NULL)
  cat(sprintf("Ljung-Box (gl=K-%d): p(8)=%.3f p(12)=%.3f p(16)=%.3f | JB=%.2f (p=%.3f; As=%.2f, K=%.2f) | LB residuos^2 (12): p=%.3f | |e|>3: %s\n",
              k, SEL$LB8[i], SEL$LB12[i], SEL$LB16[i], jb["JB"], jb["p"], jb["As"], jb["K"], lb2$p.value,
              if (length(atip)) paste(fechas[atip], collapse = " ") else "ninguno"))
  
  # Sobreajuste (Nota 6 Sec. 10.4): un coeficiente de mas, significativo?
  extras <- list(c(p + 1, q, Qs), c(p, q + 1, Qs), if (Qs == 0) c(p, q, 1))
  for (ex in extras) {
    if (is.null(ex) || ex[1] > P_MAX + 1) next
    m2 <- ajustar(x, ex[1], ex[2], ex[3]); if (is.null(m2)) next
    nuevo <- if (ex[1] > p) paste0("ar", ex[1]) else if (ex[2] > q) paste0("ma", ex[2]) else "sma1"
    if (!nuevo %in% names(m2$coef)) next
    t_nuevo <- m2$coef[nuevo] / sqrt(diag(m2$var.coef))[nuevo]
    sobre[[length(sobre) + 1]] <- data.frame(
      serie = nm, muestra = mu, base = nombre_mod, extra = nuevo, t_extra = as.numeric(t_nuevo),
      signif_5pct = abs(t_nuevo) > 1.96, dSBC = BIC(m2) - BIC(m), row.names = NULL)
  }
  
  # Grafico de residuos
  col <- colores[nm]; banda <- 1.96 / sqrt(n)
  d1 <- data.frame(t = as.numeric(time(x)), e = ee_std)
  g1 <- ggplot(d1, aes(t, e)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3, fill = gris, alpha = 0.35) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_line(colour = col, linewidth = 0.6) + geom_point(colour = col, size = 1.1) +
    scale_x_continuous(breaks = seq(2000, 2025, 5)) +
    labs(title = "Residuos estandarizados", subtitle = "franja: +-3", x = NULL, y = NULL) + tema
  ra <- acf(e, lag.max = 16, plot = FALSE)$acf[-1]
  d2 <- data.frame(k = 1:16, r = as.numeric(ra), sig = abs(ra) > banda)
  g2 <- ggplot(d2, aes(k, r)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -banda, ymax = banda, fill = gris, alpha = 0.45) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.3, lineend = "round") +
    scale_colour_manual(values = c(`TRUE` = unname(col), `FALSE` = "#9AA5B1"), guide = "none") +
    scale_x_continuous(breaks = c(1, 4, 8, 12, 16)) + scale_y_continuous(limits = c(-0.6, 1)) +
    labs(title = "FAC de los residuos", subtitle = sprintf("banda +-1.96/sqrt(T); Ljung-Box(12) p=%.2f", SEL$LB12[i]),
         x = "Rezago (trimestres)", y = NULL) + tema
  qq <- data.frame(teo = qnorm(ppoints(n))[order(order(ee_std))], mue = ee_std)
  g3 <- ggplot(qq, aes(teo, mue)) +
    geom_abline(slope = 1, intercept = 0, colour = tinta2, linewidth = 0.4, linetype = "22") +
    geom_point(colour = col, size = 1.3) +
    labs(title = "Q-Q normal", subtitle = sprintf("Jarque-Bera p=%.3f", jb["p"]), x = "Cuantil teorico", y = NULL) + tema
  png_paneles[[i]] <- (g1 | g2 | g3) +
    plot_annotation(title = sprintf("%s - %s", sub(" \\(.*", "", nm), nombre_mod),
                    subtitle = mu, theme = tema_global)
}

COEF <- do.call(rbind, coef_tab); DIAG <- do.call(rbind, diag_tab); SOBRE <- do.call(rbind, sobre)
write.csv(COEF,  "output/tablas/03_coeficientes.csv", row.names = FALSE)
write.csv(DIAG,  "output/tablas/03_diagnostico.csv",  row.names = FALSE)
write.csv(SOBRE, "output/tablas/03_sobreajuste.csv",  row.names = FALSE)

cat("\n=== Sobreajuste: coeficiente adicional (Nota 6, Sec. 10.4) ===\n")
print(SOBRE, row.names = FALSE, digits = 3)

# Una figura por serie, con las dos muestras apiladas
for (nm in names(series)) {
  idx <- which(SEL$serie == nm); if (!length(idx)) next
  fig <- wrap_plots(png_paneles[idx], ncol = 1)
  fn  <- paste0("output/figuras/03_residuos_", if (grepl("PBI", nm)) "pbi" else "elec", ".png")
  ggsave(fn, fig, width = 11, height = 4.2 * length(idx), dpi = 150, bg = fondo)
}

# ---- 5. Estabilidad frente al COVID (regla 7 y 9 de la clase) --------------------------------
cat("\n=== Coeficientes de los modelos elegidos, con y sin COVID ===\n")
print(transform(COEF, estimacion = round(estimacion, 3), ee = round(ee, 3), t = round(t, 2)),
      row.names = FALSE)

cat("\nListo. Tablas en output/tablas/03_*.csv, figuras en output/figuras/03_residuos_*.png.\n")