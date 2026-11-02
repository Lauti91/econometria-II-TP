# =============================================================================
# eviews.R  --  Salidas con el formato de EViews (se usa con source(), despues de utils.R)
#   salida_arma(m, W, ...)  ->  modelo ARMA/ARIMA estimado con stats::arima (ML exacto)
#   salida_mco(m, ...)      ->  regresion lm() (incluye la ecuacion del test ADF)
#   guardar_salida(S, archivo)
# Los estadisticos replican las definiciones de EViews (criterios de informacion POR OBSERVACION):
#   AIC = -2 l/n + 2k/n ;  SC = -2 l/n + k ln(n)/n ;  HQ = -2 l/n + 2k ln(ln n)/n
#   R2 = 1 - SSR/SST (SST sobre la variable dependiente tal como entra al modelo)
#   S.E. regresion = raiz(SSR/(n-k)) ;  F = [R2/(k-1)] / [(1-R2)/(n-k)]
# k = cantidad de coeficientes (incluidas constante y regresores, sin la varianza).
# Diferencias con EViews: stats::arima calcula la covarianza con el Hessiano numerico (EViews, por
# defecto, con el producto externo de gradientes) y la convencion del MA es la de EViews/R: +theta.
# =============================================================================
suppressPackageStartupMessages(library(ragg))

fmt <- function(x, d = 6) sub("^-(0\\.0+)$", "\\1", formatC(x, format = "f", digits = d))
fmt_raiz <- function(z) {                    # raices invertidas como en EViews: reales o a +- bi
  z <- z[order(-Mod(z))]; out <- character(0); usados <- logical(length(z))
  for (i in seq_along(z)) {
    if (usados[i]) next
    if (abs(Im(z[i])) < 1e-8) { out <- c(out, fmt(Re(z[i]), 2)); usados[i] <- TRUE
    } else {
      j <- which(!usados & seq_along(z) != i & abs(z - Conj(z[i])) < 1e-6)[1]
      usados[c(i, j)] <- TRUE
      out <- c(out, sprintf("%s +/- %si", fmt(Re(z[i]), 2), fmt(abs(Im(z[i])), 2)))
    }
  }
  out
}
nombre_coef <- function(v, s = 4) {
  v <- sub("^intercept$", "C", v); v <- sub("^ar([0-9]+)$", "AR(\\1)", v)
  v <- sub("^ma([0-9]+)$", "MA(\\1)", v); sub("^tend$", "@TREND", sub("^sma([0-9]+)$", sprintf("SMA(%d)", s), v))
}
dw <- function(e) sum(diff(e)^2) / sum(e^2)

# ---- Bloque comun: arma la lista de "estadisticos" y dibuja ------------------------------------
dibujar <- function(titulo, cabecera, tab, izq, der, extra = NULL, previo = NULL) {
  fuente <- "Liberation Sans"; serif <- "Latin Modern Roman"; tam <- 3.15
  el <- list(); seg <- list(); y <- 0
  add <- function(x, y, lab, h = 0.5, face = "plain", fam = fuente, s = tam, col = tinta)
    el[[length(el) + 1]] <<- data.frame(x = x, y = y, lab = lab, h = h, face = face, fam = fam, s = s, col = col)
  linea <- function(y, x0 = 0.5, x1 = 99.5, lw = 0.35, col = "#6B7480")
    seg[[length(seg) + 1]] <<- data.frame(x = x0, xend = x1, y = y, yend = y, lw = lw, col = col)

  add(50, y, titulo, 0.5, "bold", serif, 5.2); y <- y - 2.1
  y_top <- y + 0.65
  for (ln in cabecera) { add(2, y, ln, 0); y <- y - 1 }
  if (!is.null(previo)) {                                           # bloque tipo "ADF test statistic"
    y <- y - 0.2; linea(y + 0.55)
    add(97, y - 0.1, "t-Statistic", 1, "bold", s = tam - 0.1); y <- y - 1.15; linea(y + 0.5)
    for (r in previo$filas) { add(2, y, r[1], 0); if (length(r) > 1) add(97, y, r[2], 1); y <- y - 1 }
    y <- y - 0.1; linea(y + 0.5)
    for (ln in previo$post) { add(2, y, ln, 0); y <- y - 1 }
  }
  y <- y - 0.2; linea(y + 0.55)
  xs <- c(17, 43, 60, 78, 97)                                        # variable | coef | ee | t | prob
  add(xs[1], y - 0.4, "Variable", 0.5, "bold", s = tam - 0.1)
  for (j in 2:5) add(xs[j], y - 0.4, c("Coefficient", "Std. Error", "t-Statistic", "Prob.")[j - 1], 1, "bold", s = tam - 0.1)
  y <- y - 1.35; linea(y + 0.5)
  for (i in seq_len(nrow(tab))) {
    add(xs[1], y, tab$var[i], 0.5)
    for (j in 2:5) add(xs[j], y, tab[[j]][i], 1)
    y <- y - 1
  }
  y <- y - 0.1; linea(y + 0.5)
  for (i in seq_along(izq$lab)) {
    add(2, y, izq$lab[i], 0); add(47, y, izq$val[i], 1)
    add(52, y, der$lab[i], 0); add(97, y, der$val[i], 1); y <- y - 1
  }
  if (!is.null(extra)) {
    y <- y - 0.1; linea(y + 0.5)
    for (ln in extra) { add(2, y, ln$lab, 0, "plain"); add(47, y, ln$val, 0, "plain", s = tam); y <- y - 1 }
  }
  y_bot <- y + 0.45
  rect <- data.frame(xmin = 0.5, xmax = 99.5, ymin = y_bot, ymax = y_top)
  E <- do.call(rbind, el); S <- do.call(rbind, seg)
  p <- ggplot() +
    geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "#6B7480", linewidth = 0.4) +
    geom_segment(data = S, aes(x = x, xend = xend, y = y, yend = yend), colour = S$col, linewidth = 0.3) +
    geom_text(data = E, aes(x = x, y = y, label = lab, hjust = h, fontface = face, family = fam), size = E$s, colour = E$col) +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_bot - 0.5, 1.2), expand = c(0, 0)) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  attr(p, "alto") <- (1.2 - y_bot + 0.5) * 0.205
  p
}

# ---- ARMA / ARIMA (stats::arima) ----------------------------------------------------------------
# W: serie tal como entra al modelo (ya diferenciada); dep: nombre de la variable dependiente.
salida_arma <- function(m, W, dep, titulo, metodo = "ARMA Maximum Likelihood (BFGS, exact)") {
  e <- as.numeric(residuals(m)); n <- length(e); b <- m$coef; k <- length(b); l <- as.numeric(logLik(m))
  est <- b; se <- sqrt(diag(m$var.coef))[names(b)]; tt <- est / se
  p <- 2 * pt(-abs(tt), df = n - k)
  tab <- data.frame(var = nombre_coef(names(b), frequency(W)), coef = fmt(est), se = fmt(se), t = fmt(tt), p = fmt(p, 4))
  ssr <- sum(e^2); sst <- sum((as.numeric(W) - mean(W))^2); r2 <- 1 - ssr / sst
  r2a <- 1 - (1 - r2) * (n - 1) / (n - k)
  tiene_c <- "intercept" %in% names(b); kF <- k - tiene_c
  Fst <- if (kF > 0) (r2 / kF) / ((1 - r2) / (n - k)) else NA
  pF <- if (kF > 0) pf(Fst, kF, n - k, lower.tail = FALSE) else NA
  ic <- c(-2 * l / n + 2 * k / n, -2 * l / n + k * log(n) / n, -2 * l / n + 2 * k * log(log(n)) / n)
  f <- if (frequency(W) == 4) etiqueta_trim(W) else { t <- as.numeric(time(W)) + 1e-6; sprintf("%dM%02d", floor(t), floor((t - floor(t)) * frequency(W)) + 1) }
  cab <- c(paste("Dependent Variable:", dep), paste("Method:", metodo),
           sprintf("Sample: %s %s", f[1], f[length(f)]), sprintf("Included observations: %d", n),
           if (m$code == 0) "Convergence achieved" else "WARNING: convergence not achieved",
           "Coefficient covariance computed using the observed Hessian")
  izq <- list(lab = c("R-squared", "Adjusted R-squared", "S.E. of regression", "Sum squared resid", "Log likelihood", "F-statistic", "Prob(F-statistic)"),
              val = c(fmt(r2), fmt(r2a), fmt(sqrt(ssr / (n - k))), fmt(ssr), fmt(l), if (is.na(Fst)) "NA" else fmt(Fst, 6), if (is.na(pF)) "NA" else fmt(pF)))
  der <- list(lab = c("Mean dependent var", "S.D. dependent var", "Akaike info criterion", "Schwarz criterion", "Hannan-Quinn criter.", "Durbin-Watson stat", ""),
              val = c(fmt(mean(W)), fmt(sd(W)), fmt(ic[1]), fmt(ic[2]), fmt(ic[3]), fmt(dw(e)), ""))
  rt <- modulos(m)
  # Raices en renglones de a 4 (con MA estacional de rezago 12 hay 12 o mas raices)
  renglones <- function(lab, z) {
    v <- if (length(z)) fmt_raiz(1 / z) else "(none)"
    g <- split(v, ceiling(seq_along(v) / 4))
    lapply(seq_along(g), function(i) list(lab = if (i == 1) lab else "", val = paste(g[[i]], collapse = "    ")))
  }
  ex <- c(renglones("Inverted AR Roots", rt$G), renglones("Inverted MA Roots", rt$H))
  dibujar(titulo, cab, tab, izq, der, ex)
}

# ---- Minimos cuadrados (lm) -----------------------------------------------------------------------
salida_mco <- function(m, dep, titulo, muestra, metodo = "Least Squares", previo = NULL, cab_previa = NULL, nombres = NULL) {
  s <- summary(m); cf <- s$coefficients; e <- residuals(m); n <- length(e); k <- nrow(cf)
  nm <- sub("^\\(Intercept\\)$", "C", rownames(cf)); if (!is.null(nombres)) nm <- unname(ifelse(is.na(nombres[nm]), nm, nombres[nm]))
  tab <- data.frame(var = nm, coef = fmt(cf[, 1]), se = fmt(cf[, 2]), t = fmt(cf[, 3]), p = fmt(cf[, 4], 4))
  y <- model.response(model.frame(m)); l <- as.numeric(logLik(m)); ssr <- sum(e^2)
  ic <- c(-2 * l / n + 2 * k / n, -2 * l / n + k * log(n) / n, -2 * l / n + 2 * k * log(log(n)) / n)
  fs <- s$fstatistic
  cab <- c(cab_previa, paste("Dependent Variable:", dep), paste("Method:", metodo), paste("Sample:", muestra),
           sprintf("Included observations: %d", n))
  izq <- list(lab = c("R-squared", "Adjusted R-squared", "S.E. of regression", "Sum squared resid", "Log likelihood", "F-statistic", "Prob(F-statistic)"),
              val = c(fmt(s$r.squared), fmt(s$adj.r.squared), fmt(s$sigma), fmt(ssr), fmt(l), fmt(fs[1]),
                      fmt(pf(fs[1], fs[2], fs[3], lower.tail = FALSE))))
  der <- list(lab = c("Mean dependent var", "S.D. dependent var", "Akaike info criterion", "Schwarz criterion", "Hannan-Quinn criter.", "Durbin-Watson stat", ""),
              val = c(fmt(mean(y)), fmt(sd(y)), fmt(ic[1]), fmt(ic[2]), fmt(ic[3]), fmt(dw(e)), ""))
  dibujar(titulo, cab, tab, izq, der, previo = previo)
}

guardar_salida <- function(p, archivo, ancho = 7.4)
  ggsave(archivo, p, width = ancho, height = attr(p, "alto"), dpi = 220, device = ragg::agg_png)
