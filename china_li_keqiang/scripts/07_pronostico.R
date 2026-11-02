# =============================================================================
# 07_pronostico.R -- Evaluacion de pronosticos fuera de muestra (Nota 5, Sec. 9; Nota 6, Sec. 11)
# Variable pronosticada: crecimiento interanual (yoy, %) = x_{t} - x_{t-4}, con x = 100 ln(nivel).
# Diseno:
#   - Estimacion hasta 2015Q4; evaluacion 2016Q1-2019Q4 (antes del COVID, para no evaluar con un
#     shock que ningun modelo univariado puede anticipar). PBI y electricidad desde 2001Q2,
#     ferrocarril desde 2006Q2 (primer dato de (1-B)(1-B^4)x en cada serie).
#   - Los cuatro filtros de 05 se estiman como modelos del log-nivel x (arima con d y D), con las
#     partes deterministicas como regresores y el AO 2008Q4. Cada uno arranca x tantos trimestres
#     antes como diferencias tenga, de modo que la verosimilitud de todos cubre las mismas fechas:
#     asi el SBC es comparable entre filtros (argumento del jacobiano unitario de 05).
#       F1: (1-B) x  = deriva + dummies + ARMA        F2: (1-B^4) x = media + ARMA
#       F2t: (1-B^4) x = media + tendencia + ARMA     F3: (1-B)(1-B^4) x = ARMA
#   - El orden (p,q,Q) se elige SOLO con la muestra de estimacion (Nota 5, Sec. 9).
#   - Origen movil con ventana creciente: en cada origen se reestiman los parametros y se
#     pronostica x a 1-4 trimestres; yoy_{T+h} = x_{T+h} - x_{T+h-4}, y como T+h-4 <= T para h <= 4,
#     el error del yoy es el error de x y su error estandar es el de predict().
#   - Combinaciones (Nota 6, Sec. 11): pesos iguales de los dos mejores filtros por SBC en la
#     estimacion, y pesos SBC sobre los cuatro (Sec. 11.3). Ingenuo: yoy_{T+h} = yoy_T.
#   - Metricas: RMSE, MAE, U de Theil acotada (Nota 5), RMSE relativo al ingenuo, cobertura 95%.
# Ejecutar desde china_li_keqiang/: Rscript scripts/07_pronostico.R
# =============================================================================
source("scripts/utils.R")
dir.create("output/tablas", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figuras", showWarnings = FALSE, recursive = TRUE)

L <- cargar_niveles()
X <- list("PBI real" = L[["ln PBI real (x100)"]], "Electricidad" = L[["ln Electricidad (x100)"]],
          "Ferrocarril" = L[["ln Carga ferroviaria (x100)"]])
colores <- c("PBI real" = col_pbi, "Electricidad" = col_elec, "Ferrocarril" = col_rail)
H <- 4; FIN_EST <- 2015.75; ULT_OBS <- 2019.75
FILTROS <- c("F1", "F2", "F2t", "F3")
retardo <- c(F1 = 1, F2 = 4, F2t = 4, F3 = 5)        # observaciones de x que consume cada filtro

# Regresores en el nivel x (arima los diferencia con el mismo d y D que la serie):
#   F1: cumsum de dummies centradas -> dummies; t -> deriva.   F2: t/4 -> media (la (1-B^4) de t/4 es 1).
#   F2t: ademas t^2/8 + t/2, cuya (1-B^4) es t -> tendencia en yoy.   Todos: AO 2008Q4 en el nivel.
regresores <- function(tiempos, f) {
  t <- seq_along(tiempos); q <- round((tiempos %% 1) * 4) + 1
  D <- sapply(1:3, function(j) cumsum((q == j) - (q == 4))); colnames(D) <- paste0("S", 1:3)
  ao <- as.numeric(abs(tiempos - 2008.75) < 1e-6)
  switch(f, F1 = cbind(D, deriva = t, AO = ao), F2 = cbind(media = t / 4, AO = ao),
         F2t = cbind(media = t / 4, tend = t^2 / 8 + t / 2, AO = ao), F3 = cbind(AO = ao))
}
ajuste <- function(x, t0, f, ord) {
  xx <- window(x, start = t0 - retardo[f] / 4)
  tiempos <- as.numeric(time(xx)); n <- length(xx)
  todo <- regresores(c(tiempos, tail(tiempos, 1) + (1:H) / 4), f)
  xr <- todo[1:n, , drop = FALSE]; nx <- todo[n + 1:H, , drop = FALSE]
  d <- if (f %in% c("F1", "F3")) 1 else 0; D <- if (f == "F1") 0 else 1
  m <- tryCatch(suppressWarnings(arima(xx, order = c(ord[1], d, ord[2]), seasonal = list(order = c(0, D, ord[3]), period = 4),
                                       xreg = xr, method = "ML")), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  # predict() se llama aca adentro: predict.Arima reevalua el xreg de la llamada en este entorno
  list(m = m, pr = predict(m, n.ahead = H, newxreg = nx), n = n - retardo[f])
}

# ---- 1. Eleccion del orden con la muestra de estimacion (sobre W = F(B) x, como en 05) ----------
W_de <- function(x, f) switch(f, F1 = diff(x), F2 = diff(x, 4), F2t = diff(x, 4), F3 = diff(diff(x, 4)))
orden <- list(); sel <- list()
for (s in names(X)) {
  x <- X[[s]]; t0 <- time(x)[6]; xe <- window(x, end = FIN_EST)
  for (f in FILTROS) {
    W <- window(W_de(xe, f), start = t0)
    xr <- regresores(as.numeric(time(xe)), f)
    xr <- switch(f, F1 = diff(xr), F2 = diff(xr, 4), F2t = diff(xr, 4), F3 = diff(diff(xr, 4)))
    xr <- tail(xr, length(W)); xr <- xr[, apply(xr, 2, function(v) any(v != v[1])) | colnames(xr) == "AO", drop = FALSE]
    xr <- xr[, colnames(xr) != "media", drop = FALSE]          # la media la pone rejilla()
    # F1 (deriva) y F2/F2t (media) llevan constante; F3 no, igual que en el modelo de pronostico
    C <- rejilla(W, xr, media = f != "F3"); Dm <- C[C$admisible, ]; b <- Dm[which.min(Dm$SBC), ]
    orden[[s]][[f]] <- c(b$p, b$q, b$Q_est)
    sel[[length(sel) + 1]] <- data.frame(serie = s, filtro = f, n = length(W), desde = etiqueta_trim(W)[1],
      modelo = sprintf("ARMA(%d,%d)%s", b$p, b$q, ifelse(b$Q_est == 1, "+SMA(4)", "")), SBC = b$SBC, admisibles = nrow(Dm))
  }
}
SEL <- do.call(rbind, sel)
SEL$dSBC <- ave(SEL$SBC, SEL$serie, FUN = function(v) v - min(v))
write.csv(SEL, "output/tablas/07_orden_estimacion.csv", row.names = FALSE)
cat("=== 1. Orden elegido con la muestra de estimacion (hasta 2015Q4) ===\n")
print(transform(SEL, SBC = round(SBC, 2), dSBC = round(dSBC, 2)), row.names = FALSE)
top2 <- lapply(split(SEL, SEL$serie), function(d) d$filtro[order(d$SBC)][1:2])

# ---- 2. Pronostico con origen movil -----------------------------------------------------------
res <- list()
for (s in names(X)) {
  x <- X[[s]]; t0 <- time(x)[6]; tt <- as.numeric(time(x))
  origenes <- tt[tt >= FIN_EST - 1e-6 & tt <= ULT_OBS - 0.25 + 1e-6]
  for (o in origenes) {
    xo <- window(x, end = o); fc <- list(); bic <- c()
    for (f in FILTROS) {
      aj <- ajuste(xo, t0, f, orden[[s]][[f]]); if (is.null(aj)) next
      obj <- o + (1:H) / 4
      base <- sapply(obj - 1, function(b) as.numeric(window(x, start = b, end = b)))      # x_{T+h-4}
      fc[[f]] <- data.frame(pron = as.numeric(aj$pr$pred) - base, ee = as.numeric(aj$pr$se))
      bic[f] <- -2 * aj$m$loglik + log(aj$n) * (length(aj$m$coef) + 1)
    }
    w <- exp(-0.5 * (bic - min(bic))); w <- w / sum(w)
    fc[["Comb. 2 mejores"]] <- data.frame(pron = (fc[[top2[[s]][1]]]$pron + fc[[top2[[s]][2]]]$pron) / 2, ee = NA)
    fc[["Comb. pesos SBC"]] <- data.frame(pron = Reduce(`+`, lapply(names(w), function(f) w[f] * fc[[f]]$pron)), ee = NA)
    yoy_T <- as.numeric(window(x, start = o, end = o)) - as.numeric(window(x, start = o - 1, end = o - 1))
    fc[["Ingenuo"]] <- data.frame(pron = rep(yoy_T, H), ee = NA)
    for (mo in names(fc)) for (h in 1:H) {
      t_obj <- o + h / 4; if (t_obj > ULT_OBS + 1e-6) next
      obs <- as.numeric(window(x, start = t_obj, end = t_obj)) - as.numeric(window(x, start = t_obj - 1, end = t_obj - 1))
      res[[length(res) + 1]] <- data.frame(serie = s, modelo = mo, origen = o, h = h, objetivo = t_obj,
        observado = obs, pronostico = fc[[mo]]$pron[h], ee = fc[[mo]]$ee[h],
        peso_SBC = if (mo %in% names(w)) unname(w[mo]) else NA)
    }
  }
}
R <- do.call(rbind, res); R$error <- R$observado - R$pronostico
R$cubre95 <- ifelse(is.na(R$ee), NA, abs(R$error) <= 1.96 * R$ee)
write.csv(R, "output/tablas/07_pronosticos.csv", row.names = FALSE)

# ---- 3. Metricas ---------------------------------------------------------------------------
met <- do.call(rbind, lapply(split(R, list(R$serie, R$modelo, R$h), drop = TRUE), function(d) {
  ing <- R[R$serie == d$serie[1] & R$modelo == "Ingenuo" & R$h == d$h[1], ]
  data.frame(serie = d$serie[1], modelo = d$modelo[1], h = d$h[1], N = nrow(d),
             RMSE = sqrt(mean(d$error^2)), MAE = mean(abs(d$error)),
             U_Theil = sqrt(mean(d$error^2)) / (sqrt(mean(d$observado^2)) + sqrt(mean(d$pronostico^2))),
             RMSE_rel_ingenuo = sqrt(mean(d$error^2)) / sqrt(mean(ing$error^2)),
             cobertura95 = if (all(is.na(d$cubre95))) NA else mean(d$cubre95),
             peso_SBC_medio = mean(d$peso_SBC))
}))
met <- met[order(met$serie, met$h, met$RMSE), ]
write.csv(met, "output/tablas/07_metricas.csv", row.names = FALSE)
cat("\n=== 3. Evaluacion fuera de muestra 2016Q1-2019Q4 (yoy, pp) ===\n")
num <- sapply(met, is.numeric); mp <- met; mp[num] <- lapply(met[num], round, 3)
print(mp[mp$h %in% c(1, 4), ], row.names = FALSE)
cat("\nDos mejores filtros por SBC en la estimacion (combinacion con pesos iguales):\n")
print(sapply(top2, paste, collapse = " + "))

# ---- 4. Figura: pronosticos a 1 trimestre del mejor filtro de cada serie --------------------
P1 <- do.call(rbind, lapply(names(X), function(s) R[R$serie == s & R$h == 1 & R$modelo == top2[[s]][1], ]))
P1$etiqueta <- paste0("Mejor filtro en la estimacion: ", P1$modelo)
hist_ <- do.call(rbind, lapply(names(X), function(s) {
  y <- window(diff(X[[s]], 4), start = c(2012, 1), end = c(2019, 4))
  data.frame(serie = s, t = as.numeric(time(y)), y = as.numeric(y))
}))
P1$serie <- factor(P1$serie, levels = names(X)); hist_$serie <- factor(hist_$serie, levels = names(X))
lab <- do.call(rbind, lapply(names(X), function(s) {
  m1 <- met[met$serie == s & met$h == 1 & met$modelo == top2[[s]][1], ]
  data.frame(serie = factor(s, levels = names(X)), lab = sprintf("%s | RMSE %.2f pp | relativo al ingenuo %.2f | cobertura 95%%: %.0f%%",
                                                                 top2[[s]][1], m1$RMSE, m1$RMSE_rel_ingenuo, 100 * m1$cobertura95))
}))
g <- ggplot() +
  annotate("rect", xmin = 2015.875, xmax = 2019.875, ymin = -Inf, ymax = Inf, fill = gris, alpha = 0.3) +
  geom_ribbon(data = P1, aes(objetivo, ymin = pronostico - 1.96 * ee, ymax = pronostico + 1.96 * ee), fill = gris, alpha = 0.7) +
  geom_line(data = hist_, aes(t, y, colour = serie), linewidth = 0.9) +
  geom_point(data = P1, aes(objetivo, pronostico), colour = tinta, size = 1.6, shape = 4, stroke = 0.9) +
  geom_text(data = lab, aes(x = 2012, y = Inf, label = lab), hjust = 0, vjust = 1.4, size = 3, colour = tinta2) +
  scale_colour_manual(values = colores, guide = "none") +
  facet_wrap(~serie, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(title = "Pronosticos fuera de muestra a un trimestre, 2016-2019",
       subtitle = "Linea: crecimiento interanual observado (%). Cruces: pronostico del mejor filtro; sombra: intervalo del 95%; franja: tramo de evaluacion.",
       x = NULL, y = NULL) +
  tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
ggsave("output/figuras/07_pronostico.png", g, width = 11, height = 9, dpi = 150, bg = fondo)
cat("\nFigura: output/figuras/07_pronostico.png\n")
