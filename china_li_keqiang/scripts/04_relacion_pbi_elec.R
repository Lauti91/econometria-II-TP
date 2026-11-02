# =============================================================================
# 04_relacion_pbi_elec.R -- Relacion PBI oficial vs electricidad: yoy y Delta(yoy)
# Objetivo (regla 1 de Mirta): ¿el crecimiento oficial del PBI es consistente con un indicador fisico?
# Se compara dos tratamientos de las series y se aplica la alarma de regresion espuria.
#   A) Box-Jenkins univariado sobre Delta(yoy): firmas de sobrediferenciacion (Nota 6, Sec. 2.5)
#   B) Regresion estatica en yoy y en Delta(yoy), ambos sentidos, con R2 > d_W (Nota 6, Sec. 3.3)
# Fuera de alcance (2da parte de la materia): funcion de transferencia / CCF con prewhitening,
# cointegracion (Engle-Granger).
# =============================================================================
source("scripts/utils.R")
dir.create("output/tablas", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figuras", showWarnings = FALSE, recursive = TRUE)

S <- cargar_yoy(); nm_pbi <- names(S)[1]; nm_el <- names(S)[2]
muestras <- list("completa 2000Q1-2023Q4" = c(2000,1,2023,4), "sin COVID 2000Q1-2019Q4" = c(2000,1,2019,4))
recorte <- function(x, w) window(x, start = w[1:2], end = w[3:4])
d1 <- function(x) diff(x)                      # Delta = (1 - B)

# ---- A. Box-Jenkins sobre Delta(yoy) vs sobre yoy ---------------------------------------------
# Si yoy ya es estacionaria (I(0)), diferenciar de mas introduce una raiz MA unitaria (Nota 6, Sec. 2.5):
# Theta(1) = 0, min|H| ~ 1, y la varianza del error sube. Se compara con el modelo elegido en yoy.
mejor <- function(x) {
  C <- rejilla(x); D <- C[C$admisible, ]
  if (!nrow(D)) return(NULL)
  D[which.min(D$SBC), ]
}
A <- list()
for (nm in names(S)) for (mu in names(muestras)) for (tr in c("yoy", "Delta(yoy)")) {
  x <- recorte(S[[nm]], muestras[[mu]]); if (tr != "yoy") x <- d1(x)
  b <- mejor(x)
  if (is.null(b)) { A[[length(A)+1]] <- data.frame(serie = nm, muestra = mu, tratamiento = tr, n = length(x), modelo = "ninguno admisible"); next }
  A[[length(A)+1]] <- data.frame(serie = nm, muestra = mu, tratamiento = tr, n = length(x),
    modelo = sprintf("ARMA(%d,%d)%s", b$p, b$q, ifelse(b$Q_est == 1, "+SMA(4)", "")),
    SBC = b$SBC, sigma = sqrt(b$sigma2), min_mod_AR = b$min_mod_G, min_mod_MA = b$min_mod_H,
    LB12_p = b$LB12)
}
TA <- do.call(rbind, A)
TA[sapply(TA, is.numeric)] <- lapply(TA[sapply(TA, is.numeric)], round, 3)
write.csv(TA, "output/tablas/04_bj_delta.csv", row.names = FALSE)
cat("\n=== A. Modelo elegido (min SBC entre admisibles) en yoy vs Delta(yoy) ===\n"); print(TA, row.names = FALSE)
cat("Lectura: min_mod_MA cercano a 1 en Delta(yoy) y sigma mayor que en yoy => sobrediferenciacion (Nota 6, Sec. 2.5).\n")

# ---- B. Regresion estatica -------------------------------------------------------------------
dw  <- function(e) sum(diff(e)^2) / sum(e^2)
reg <- function(y, x, etiqueta) {
  m <- lm(y ~ x); e <- residuals(m); s <- summary(m)
  r1 <- unname(acf(e, plot = FALSE, lag.max = 1)$acf[2])
  lb <- Box.test(e, 8, "Ljung-Box")$p.value
  data.frame(spec = etiqueta, n = length(e), pendiente = coef(m)[2], t_ingenuo = s$coefficients[2, 3],
             R2 = s$r.squared, DW = dw(e), alarma_R2_mayor_DW = s$r.squared > dw(e),
             rho1_resid = r1, LB8_p = lb, row.names = NULL)
}
B <- list(); RES <- list()
for (mu in names(muestras)) {
  w <- muestras[[mu]]
  p <- recorte(S[[nm_pbi]], w); el <- recorte(S[[nm_el]], w)
  for (tr in c("yoy", "Delta(yoy)")) {
    py <- if (tr == "yoy") p else d1(p); ey <- if (tr == "yoy") el else d1(el)
    B[[length(B)+1]] <- cbind(muestra = mu, tratamiento = tr, direccion = "PBI ~ Electricidad", reg(as.numeric(py), as.numeric(ey), ""))
    B[[length(B)+1]] <- cbind(muestra = mu, tratamiento = tr, direccion = "Electricidad ~ PBI", reg(as.numeric(ey), as.numeric(py), ""))
    if (mu == names(muestras)[1])
      RES[[tr]] <- list(y = py, x = ey, e = residuals(lm(as.numeric(py) ~ as.numeric(ey))))
  }
}
TB <- do.call(rbind, B); TB$spec <- NULL
TB[sapply(TB, is.numeric)] <- lapply(TB[sapply(TB, is.numeric)], round, 3)
write.csv(TB, "output/tablas/04_regresion_cruzada.csv", row.names = FALSE)
cat("\n=== B. Regresion estatica y alarma de Granger-Newbold (R2 > d_W) ===\n"); print(TB, row.names = FALSE)
cat("Nota: t_ingenuo y R2 suponen residuos ruido blanco; con rho1 alto o DW bajo no son validos (Nota 6, Sec. 3.3).\n")
cat("Un valor de DW cercano a 2 y alarma FALSE es condicion necesaria, no suficiente: ver LB8_p.\n")

# ---- C. Figura -------------------------------------------------------------------------------
cov <- function(x) { t <- as.numeric(time(x)); t >= 2020 & t < 2022 }
sc <- function(r, tr, xl, yl) {
  df <- data.frame(x = as.numeric(r$x), y = as.numeric(r$y), cov = cov(r$y))
  ggplot(df, aes(x, y)) +
    geom_hline(yintercept = 0, colour = gris, linewidth = 0.3) + geom_vline(xintercept = 0, colour = gris, linewidth = 0.3) +
    geom_point(aes(colour = cov), size = 1.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, colour = tinta, linewidth = 0.5, formula = y ~ x) +
    scale_colour_manual(values = c(`FALSE` = col_pbi, `TRUE` = "#B42318"), labels = c("2000-2019 y 2022-23", "COVID 2020-21"), name = NULL) +
    labs(title = tr, x = xl, y = yl) + tema + theme(legend.position = "bottom")
}
ac <- function(e, tit) {
  a <- acf(e, plot = FALSE, lag.max = 16); b <- 1.96 / sqrt(length(e))
  df <- data.frame(k = a$lag[-1], r = a$acf[-1])
  ggplot(df, aes(k, r)) + annotate("rect", xmin = 0.5, xmax = 16.5, ymin = -b, ymax = b, fill = gris, alpha = 0.4) +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_segment(aes(xend = k, yend = 0), colour = col_pbi, linewidth = 0.8) +
    scale_y_continuous(limits = c(-0.6, 1)) + scale_x_continuous(breaks = c(1, seq(4, 16, 4))) +
    labs(title = tit, x = "Rezago", y = "FAC de residuos") + tema
}
f <- (sc(RES[["yoy"]], "En crecimiento interanual", "Electricidad, yoy (%)", "PBI, yoy (%)") |
      sc(RES[["Delta(yoy)"]], "En primeras diferencias de yoy", "Var. de electricidad yoy (p.p.)", "Var. del PBI yoy (p.p.)")) /
     (ac(RES[["yoy"]]$e, "Residuos de la regresion en yoy") | ac(RES[["Delta(yoy)"]]$e, "Residuos de la regresion en dif. de yoy")) +
  plot_annotation(title = "PBI oficial vs electricidad: la relacion aguanta al cambiar el tratamiento?",
                  subtitle = "Muestra completa 2000Q1-2023Q4. Banda gris: +/-1.96/raiz(T) (Bartlett bajo ruido blanco).",
                  theme = tema_global)
ggsave("output/figuras/04_relacion.png", f, width = 11, height = 8, dpi = 150)
cat("\nFigura: output/figuras/04_relacion.png\n")
