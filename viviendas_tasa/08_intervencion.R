# =============================================================================
# 08_intervencion.R -- Analisis de intervencion: cierre por COVID (abril de 2020) sobre los inicios
#   Apunte parte 2, pp. 1-7; Enders cap. 5, sec. 1.
#   El modelo es x_t = omega * xi_t(B) Z_t + N_t, con N_t el ruido ARIMA(0,1,1)(0,1,1)12 de 04 y
#   Z_t un pulso en 2020-04. Se compara la FORMA de la intervencion (apunte p. 2):
#     AO: xi = 1                  efecto solo en abril (lo que usa 04)
#     TC: xi = 1/(1 - delta B)    efecto que se desvanece a la tasa delta (0 < delta < 1)
#     LS: xi = 1/(1 - B)          escalon: cambio permanente de nivel desde abril
#   Como el ruido es I(1)xI(1)_12, se estima sobre W = (1-B)(1-B^12)x y el regresor se filtra igual
#   (apunte p. 7: con raiz unitaria, un pulso en la serie diferenciada es permanente en el nivel).
#   delta del TC se elige por perfil de verosimilitud; la forma, por SBC (apunte p. 7).
# Ejecutar desde viviendas_tasa/:  Rscript 08_intervencion.R
# =============================================================================
source("R/funciones.R")
L <- cargar_base(); x <- L$viv
desde <- time(x)[14]                                              # mismas fechas que 04: desde 1986-02
W <- window(diff(diff(x, 12)), start = desde)
pulso <- indicadora(x, 2020, 4)
# Regresor de intervencion en nivel segun la forma, y su version filtrada para la ecuacion en W
regresor <- function(forma, delta = NA) {
  z <- switch(forma, AO = pulso, LS = cumsum(pulso), TC = stats::filter(pulso, delta, method = "recursive"))
  z <- ts(as.numeric(z), start = start(x), frequency = S)
  matrix(as.numeric(window(diff(diff(z, 12)), start = desde)), ncol = 1, dimnames = list(NULL, "COVID"))
}

# ---- 1. TC: perfil de verosimilitud en delta ----------------------------------------------------------
grilla <- seq(0.05, 0.95, by = 0.05)
perfil <- sapply(grilla, function(d) { m <- ajustar(W, 0, 1, 0, 1, regresor("TC", d), FALSE); if (is.null(m)) NA else as.numeric(logLik(m)) })
delta <- grilla[which.max(perfil)]
cat("delta del TC por perfil:", delta, "\n")

# ---- 2. Estimacion de las formas y comparacion ---------------------------------------------------------
formas <- c(sin = "Sin intervenci\\'on", AO = "Pulso (AO)", TC = sprintf("Cambio transitorio (TC, $\\delta=%s$)", num(delta, 2)), LS = "Escal\\'on (LS)")
M <- list(sin = ajustar(W, 0, 1, 0, 1, NULL, FALSE), AO = ajustar(W, 0, 1, 0, 1, regresor("AO"), FALSE),
          TC = ajustar(W, 0, 1, 0, 1, regresor("TC", delta), FALSE), LS = ajustar(W, 0, 1, 0, 1, regresor("LS"), FALSE))
extra <- c(sin = 0, AO = 0, TC = 1, LS = 0)                         # delta cuenta como parametro estimado
COMP <- do.call(rbind, lapply(names(M), function(f) { m <- M[[f]]; d <- diagnostico(m, 2); n <- d$n
  kk <- length(m$coef) + extra[f]; sbc <- -2 * d$logLik + kk * log(n)
  om <- if (f == "sin") NA else m$coef["COVID"]; se <- if (f == "sin") NA else sqrt(m$var.coef["COVID", "COVID"])
  data.frame(forma = f, omega = om, ee = se, t = om / se, theta = -m$coef["ma1"], Theta = -m$coef["sma1"], sigma = d$sigma,
             logLik = d$logLik, SBC = sbc, LB24 = d$LB24, JB = d$JB, row.names = NULL) }))
COMP$dSBC <- COMP$SBC - min(COMP$SBC)
print(COMP)
gana <- COMP$forma[which.min(COMP$SBC)]
write.csv(COMP, "salidas/tablas/08_formas.csv", row.names = FALSE)
f_omega <- function(i) if (is.na(COMP$omega[i])) "--" else sprintf("%s (%s)", num(COMP$omega[i], 1), num(COMP$ee[i], 1))
tabla_tex(data.frame(formas[COMP$forma], sapply(seq_len(nrow(COMP)), f_omega), num(COMP$theta, 3), num(COMP$Theta, 3), num(COMP$sigma, 2),
                     num(COMP$LB24, 2), num(COMP$JB, 1), num(COMP$SBC, 1), num(COMP$dSBC, 1)),
          "08_formas.tex", rojo = celdas(which(COMP$forma == gana), 1:9), "lcccccccc",
          c("Forma", "$\\hat\\omega$ (ee)", "$\\hat\\theta$", "$\\hat\\Theta$", "$\\hat\\sigma$", "LB24 $p$", "JB", "SBC", "$\\Delta$SBC"))
guardar_valor("IntDelta", num(delta, 2)); guardar_valor("IntGana", gana)
for (f in c("AO", "TC", "LS")) { i <- COMP$forma == f
  guardar_valor(paste0("IntOmega", f), num(COMP$omega[i], 1)); guardar_valor(paste0("IntT", f), num(COMP$t[i], 1))
  guardar_valor(paste0("IntDSBC", f), num(COMP$dSBC[i], 1)) }

# ---- 3. Respuesta al impulso de la intervencion (apunte p. 4-5) ---------------------------------------
# Efecto sobre el nivel j meses despues: omega * xi_j, con xi_j = 1 (j=0) en AO, delta^j en TC, 1 en LS.
# En % de inicios: 100 (exp(efecto/100) - 1), porque x = 100 ln(inicios).
j <- 0:24
resp <- do.call(rbind, lapply(c("AO", "TC", "LS"), function(f) { om <- COMP$omega[COMP$forma == f]
  xi <- switch(f, AO = as.numeric(j == 0), TC = delta^j, LS = rep(1, length(j)))
  data.frame(forma = f, j = j, efecto = om * xi, pct = 100 * (exp(om * xi / 100) - 1)) }))
write.csv(resp, "salidas/tablas/08_respuesta.csv", row.names = FALSE)
r_g <- resp[resp$forma == gana, ]
guardar_valor("IntImpactoPct", num(r_g$pct[1], 1))
guardar_valor("IntTresPct", num(r_g$pct[r_g$j == 3], 1)); guardar_valor("IntSeisPct", num(r_g$pct[r_g$j == 6], 1))
if (gana == "TC") guardar_valor("IntVidaMedia", num(log(0.5) / log(delta), 1))   # meses hasta que el efecto cae a la mitad
etq <- c(AO = "Pulso (AO)", TC = sprintf("Transitorio (TC, delta = %.2f)", delta), LS = "Escalon (LS)")
g1 <- ggplot(resp, aes(j, pct, colour = forma == gana, group = forma)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  geom_text(data = resp[(resp$forma == "AO" & resp$j == 1) | (resp$forma == "TC" & resp$j == 3) | (resp$forma == "LS" & resp$j == 24), ],
            aes(label = etq[forma]), hjust = c(-0.05, -0.08, 1), vjust = c(-0.6, 1.2, -0.6), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = rojo, `FALSE` = "#9AA5B1"), guide = "none") + scale_x_continuous(breaks = seq(0, 24, 6)) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "Respuesta de los inicios a la intervencion de abril de 2020", subtitle = sprintf("en rojo, la forma elegida por SBC (%s)", gana),
       x = "Meses desde abril de 2020", y = "% de los inicios") + tema
# Perfil de verosimilitud del delta
g2 <- ggplot(data.frame(d = grilla, l = perfil), aes(d, l)) + geom_line(colour = col_viv, linewidth = 0.6) + geom_point(colour = col_viv, size = 1) +
  geom_vline(xintercept = delta, colour = rojo, linetype = "22") +
  labs(title = "Perfil de verosimilitud del TC", subtitle = sprintf("maximo en delta = %.2f", delta), x = "delta", y = "log-verosimilitud") + tema
guardar_fig(g1 | g2, "08_intervencion.png", 11, 3.8)

# ---- 4. Salida EViews del modelo elegido y diagnostico -------------------------------------------------
m <- M[[gana]]
guardar_salida(salida_arma(m, W, "D(LHOUST,1,12)", sprintf("Inicios: ARIMA(0,1,1)(0,1,1)12 con intervencion %s 2020M04", gana)), "08_modelo_intervencion.png")
guardar_salida(salida_correlograma(residuals(m), "Correlograma de los residuos: modelo con intervencion", 36, 2, rojo_k = c(12, 24, 36)), "08_corr_residuos.png")
es <- residuals(m) / sqrt(m$sigma2)
guardar_valor("IntResidAbr", num(es[etiqueta_mes(W) == "2020-04"], 1)); guardar_valor("IntResidMay", num(es[etiqueta_mes(W) == "2020-05"], 1))
saveRDS(list(delta = delta, forma = gana), "salidas/tablas/08_intervencion.rds")   # lo usan 09 y 10
escribir_valores("valores_08.tex")
