# =============================================================================
# 16b_estabilidad_ganancia.R -- La ganancia de la tasa es estable? Por que se sobrepredijo 2022-24?
#   SEGUNDA RONDA (exploratoria): estos analisis se definieron DESPUES de ver los resultados de 15 y 16; no
#   estaban preregistrados. Se informan todos los que se corrieron y el script cuenta los modelos estimados.
#   Responde a la regla de cambios estructurales del TP ("trabajar por subperiodos o modelar cambios por tramos").
#   Modelo de 09, con la forma fija (b = 2, s = 0, r = 1, phi = 0,75; ruido ARMA(1,2)(1,1)12; COVID como TC).
#   1. Ganancia por subperiodo: un delta por tramo (1985-99, 2000-05, 2006-11, 2012-21, 2022-26, segun la fecha del
#      cambio de la tasa), con el ruido comun. H0: los cinco iguales (LR, chi2 con 4 gl). Quiebre en 2022: t del
#      termino adicional desde 2022 (contraste tipo Chow sobre un coeficiente; Chow, G. (1960), Econometrica,
#      28(3), 591-605).
#   2. Estimacion recursiva (muestra que crece desde 1985) y por ventanas moviles de 15 anios: trayectoria de la
#      ganancia. Idea de la estimacion recursiva: Brown, R., Durbin, J. y Evans, J. (1975), Journal of the Royal
#      Statistical Society B, 37(2), 149-192; aca solo se grafica la ganancia, sin el CUSUM.
#   3. Fase del ciclo de tasas: ganancia cuando la tasa esta por encima de la de 12 meses antes (fase de suba) y
#      cuando esta por debajo (fase de baja). Es la version SIN efecto trinquete de la asimetria de 16: no separa
#      por el signo del cambio mensual sino por la fase, y cada cambio mensual entra entero en una fase.
#   4. Descomposicion del episodio 2022-24. Identidad contable del modelo, sin contraste nuevo:
#        observado = [pronostico sin cambios de tasa] + [efecto de la tasa] + [otros factores]
#      donde el efecto de la tasa es (condicional - tasa congelada) y "otros factores" es el error del pronostico
#      condicional (shocks posteriores a 2021-12 que no son la tasa). Dos juegos de parametros: los de 12 (datos
#      hasta 2021-12) y los de toda la muestra con el estado del modelo en 2021-12 (arima con parametros fijos).
#      El error se compara con su desvio teorico (pesos psi, como en 12) y se contrasta si los residuos de un paso
#      de 2022-24 tienen media cero (z = media / (sigma / raiz(36))).
#   5. Resumen de la ganancia en todas las especificaciones (09, controles de 15, subperiodos, fases).
#   Salidas: 16b_subperiodos.tex, 16b_fase.tex, 16b_descomposicion.tex, 16b_ganancias.csv, 16b_movil.csv,
#   16b_episodio_variantes.csv; figuras 16b_estabilidad.png y 16b_descomposicion.png; valores_16b.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 16b_estabilidad_ganancia.R   (despues de 08 y 15)
# =============================================================================
source("R/funciones.R"); source("R/datos_extra.R"); source("R/tf.R")
L <- cargar_base(); LC <- cargar_base(completa = TRUE); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
tf_prueba(L, I8$delta)
Y <- diff(L$viv_sa); X <- as.numeric(diff(L$tasa)); nY <- length(Y); mes <- etiqueta_mes(Y); an <- as.integer(substr(mes, 1, 4))
TCm <- cbind(COVID = tf_tc(L$viv_sa, I8$delta)); fin <- which(mes == TF_ORIGEN); Xcf <- X; Xcf[(fin + 1):nY] <- 0
ult <- L$tasa_pp[length(L$tasa_pp)]; pp <- 100 * log((ult + 1) / ult); den <- 1 - TF_PHI; n_esp <- 0
d2 <- function(v) tf_input(v, TF_B, 0, TF_PHI)[, 1]
base_xr <- function(Xv = X) cbind(D2 = d2(Xv), TCm)
m0 <- tf_ajustar(Y, base_xr()); g0 <- m0$coef[["D2"]] / den; g0_ee <- sqrt(m0$var.coef["D2", "D2"]) / den
lr <- function(m, ref, gl) { s <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(ref))); c(LR = s, p = pchisq(max(s, 0), gl, lower.tail = FALSE)) }
gan <- function(m, nm) data.frame(g = unname(m$coef[nm]) / den, ee = unname(sqrt(diag(m$var.coef)[nm])) / den)   # ee condicional a phi
dif2 <- function(m, a, b) { d <- m$coef[c(a, b)]; V <- m$var.coef[c(a, b), c(a, b)]; z <- (d[[1]] - d[[2]]) / sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2]); c(z = z, p = 2 * pnorm(-abs(z))) }

# ---- 1. Ganancia por subperiodo y quiebre en 2022 -----------------------------------------------------------
cortes <- list(PA = c(1985, 1999), PB = c(2000, 2005), PC = c(2006, 2011), PD = c(2012, 2021), PE = c(2022, 2026))
etq_p <- c(PA = "1985--1999", PB = "2000--2005", PC = "2006--2011 (colapso inmobiliario)", PD = "2012--2021", PE = "2022--2026")
xr_per <- function(Xv) cbind(sapply(cortes, function(c) d2(Xv * (an >= c[1] & an <= c[2]))), TCm)
m_per <- tf_ajustar(Y, xr_per(X)); n_esp <- n_esp + 1; GP <- cbind(tramo = names(cortes), gan(m_per, names(cortes)), meses = sapply(cortes, function(c) sum(an >= c[1] & an <= c[2])))
LRp <- lr(m_per, m0, 4)
m_q <- tf_ajustar(Y, cbind(D2 = d2(X), DPOST = d2(X * (an >= 2022)), TCm)); n_esp <- n_esp + 1
t_q <- m_q$coef[["DPOST"]] / sqrt(m_q$var.coef["DPOST", "DPOST"]); p_q <- 2 * pnorm(-abs(t_q))
cat("\n1. Ganancia por subperiodo (ee condicional a phi)\n"); print(GP, row.names = FALSE, digits = 3)
cat(sprintf("H0: cinco ganancias iguales -> LR = %.2f, p = %.3f | quiebre en 2022: t = %.2f, p = %.3f | SBC por tramos - SBC comun = %.1f\n", LRp["LR"], LRp["p"], t_q, p_q, BIC(m_per) - BIC(m0)))
tabla_tex(rbind(data.frame(a = unname(etq_p[GP$tramo]), b = GP$meses, c = sprintf("%s (%s)", num(GP$g, 2), num(GP$ee, 2)), d = num(GP$g * pp, 1)),
                data.frame(a = "Toda la muestra (09)", b = nY, c = sprintf("%s (%s)", num(g0, 2), num(g0_ee, 2)), d = num(g0 * pp, 1))),
          "16b_subperiodos.tex", rojo = celdas(c(5, 6), 3:4), "lccc", c("Tramo (fecha del cambio de la tasa)", "Meses", "Ganancia $g$ (ee cond.)", "\\% por +1 pp"))

# ---- 2. Estimacion recursiva y ventanas moviles de 15 anios -----------------------------------------------------
i_cov <- which(mes == "2020-04")
ajuste_ventana <- function(ini, hasta) { XR <- if (ini <= i_cov - 6 && hasta >= i_cov + 6) base_xr() else base_xr()[, "D2", drop = FALSE]   # el COVID solo si la ventana lo contiene
  m <- tryCatch(tf_ajustar(Y, XR, ini = ini, hasta = hasta), error = function(e) NULL); n_esp <<- n_esp + 1
  if (is.null(m)) return(NULL); data.frame(desde = mes[ini], hasta = mes[hasta], t = as.numeric(time(Y))[hasta], gan(m, "D2"), n = hasta - ini + 1) }
cierres <- c(which(substr(mes, 6, 7) == "12" & an >= 1999), nY)
REC <- do.call(rbind, lapply(cierres, function(h) ajuste_ventana(TF_INI, h)))
MOV <- do.call(rbind, lapply(cierres[cierres - 179 >= TF_INI], function(h) ajuste_ventana(h - 179, h)))
write.csv(rbind(cbind(tipo = "recursiva", REC), cbind(tipo = "movil 15 anios", MOV)), "salidas/tablas/16b_movil.csv", row.names = FALSE)
cat(sprintf("\n2. Recursiva: ganancia entre %.2f y %.2f (%d muestras). Ventanas moviles de 15 anios: entre %.2f y %.2f (%d ventanas)\n",
            min(REC$g), max(REC$g), nrow(REC), min(MOV$g), max(MOV$g), nrow(MOV)))

# ---- 3. Fase del ciclo de tasas (sin efecto trinquete) ------------------------------------------------------------
r <- as.numeric(LC$tasa_pp); al <- match(mes, etiqueta_mes(LC$tasa_pp)); stopifnot(!anyNA(al), all(al > 12))
sube <- r[al] > r[al - 12]                                                 # la tasa esta por encima de la de un anio antes
xr_fase <- function(Xv) cbind(DSUBE = d2(Xv * sube), DBAJA = d2(Xv * !sube), TCm)
stopifnot(max(abs(xr_fase(X)[, "DSUBE"] + xr_fase(X)[, "DBAJA"] - base_xr()[, "D2"])) < 1e-9)
fase <- function(hasta) { n_esp <<- n_esp + 1; m <- tf_ajustar(Y, xr_fase(X), hasta = hasta); ref <- tf_ajustar(Y, base_xr(), hasta = hasta); G <- gan(m, c("DSUBE", "DBAJA")); d <- dif2(m, "DSUBE", "DBAJA")
  data.frame(hasta = mes[hasta], g_sube = G$g[1], ee_sube = G$ee[1], g_baja = G$g[2], ee_baja = G$ee[2], z = d[["z"]], p = d[["p"]], p_lr = lr(m, ref, 1)[["p"]], dSBC = BIC(m) - BIC(ref), meses_sube = sum(sube[TF_INI:hasta])) }
FA <- rbind(fase(nY), fase(fin)); m_fase <- tf_ajustar(Y, xr_fase(X))
E0 <- tf_episodio(Y, L$viv_sa, base_xr(X), base_xr(Xcf)); EF <- tf_episodio(Y, L$viv_sa, xr_fase(X), xr_fase(Xcf)); n_esp <- n_esp + 1
cat("\n3. Fase del ciclo de tasas\n"); print(FA, row.names = FALSE, digits = 3)
cat(sprintf("Episodio con el modelo por fases (datos hasta 2021-12): %s | RMSE %.1f, sesgo %.1f (base: RMSE %.1f, sesgo %.1f)\n", paste(round(tf_anual(EF), 1), collapse = " / "),
            tf_error(EF)[["rmse"]], tf_error(EF)[["sesgo"]], tf_error(E0)[["rmse"]], tf_error(E0)[["sesgo"]]))
tabla_tex(data.frame(ifelse(FA$hasta == mes[nY], "1985--2026", "1985--2021"), FA$meses_sube, sprintf("%s (%s)", num(FA$g_sube, 2), num(FA$ee_sube, 2)), sprintf("%s (%s)", num(FA$g_baja, 2), num(FA$ee_baja, 2)),
                     num(FA$z, 2), sapply(FA$p, pval), sapply(FA$p_lr, pval), num(FA$dSBC, 1)),
          "16b_fase.tex", rojo = celdas(1, 3:6), "lccccccc", c("Muestra", "Meses en fase de suba", "$g$ en fase de suba (ee cond.)", "$g$ en fase de baja (ee cond.)", "$z$", "Wald: $p$", "LR: $p$", "$\\Delta$SBC vs lineal"))
write.csv(tf_variante("6b", "fase", "Ganancia por fase del ciclo de tasas", EF, k_extra = 1, dSBC = FA$dSBC[1], nota = "segunda ronda, exploratoria"), "salidas/tablas/16b_episodio_variantes.csv", row.names = FALSE)

# ---- 4. Descomposicion del episodio 2022-24 -----------------------------------------------------------------------
# Parametros de toda la muestra con el estado del modelo en 2021-12: arima con todos los coeficientes fijos
xr <- base_xr()[TF_INI:fin, , drop = FALSE]
m_fijo <- arima(tf_serie(Y, TF_INI:fin), order = c(TF_RUIDO[1], 0, TF_RUIDO[2]), seasonal = list(order = c(TF_RUIDO[3], 0, TF_RUIDO[4]), period = S), xreg = xr, include.mean = TRUE,
                fixed = m0$coef, transform.pars = FALSE, method = "ML")
stopifnot(max(abs(m_fijo$coef - m0$coef)) < 1e-12)
n0 <- as.numeric(L$viv_sa)[fin + 1]; camino <- function(m, R) n0 + cumsum(tf_predecir(m, R[E0$fut, , drop = FALSE]))
DES <- list(ante = list(etq = "Par\\'ametros con datos hasta 2021-12", cond = E0$cond, contra = E0$contra, ee = E0$ee, g = E0$m$coef[["D2"]] / den),
            post = list(etq = "Par\\'ametros de toda la muestra", cond = camino(m_fijo, base_xr(X)), contra = camino(m_fijo, base_xr(Xcf)), ee = tf_ee_nivel(m0, TF_H), g = g0))
obs <- E0$obs; b21 <- E0$base; anio <- E0$anio
TD <- do.call(rbind, lapply(names(DES), function(k) { d <- DES[[k]]; do.call(rbind, lapply(2022:2024, function(a) { i <- anio == a
  data.frame(par = k, anio = a, sin_tasa = mean(d$contra[i]) - b21, tasa = mean(d$cond[i] - d$contra[i]), otros = mean(obs[i] - d$cond[i]), observado = mean(obs[i]) - b21, ee = mean(d$ee[i])) })) }))
TD$otros_ee <- TD$otros / TD$ee                                            # error del pronostico condicional en desvios teoricos
stopifnot(max(abs(TD$sin_tasa + TD$tasa + TD$otros - TD$observado)) < 1e-9)
e36 <- sapply(DES, function(d) (obs - d$cond)[c(12, 24, 36)] / d$ee[c(12, 24, 36)])
res <- as.numeric(residuals(m0))[(fin + 1:TF_H) - TF_INI + 1]; z_res <- mean(res) / (sqrt(m0$sigma2) / sqrt(TF_H))   # residuos de un paso de 2022-01 a 2024-12
stopifnot(etiqueta_mes(tf_serie(Y, TF_INI:nY))[(fin + 1) - TF_INI + 1] == "2022-01")
cat("\n4. Descomposicion del episodio (promedios anuales, en 100 ln respecto del promedio de 2021)\n"); print(TD, row.names = FALSE, digits = 3)
cat("Error del pronostico condicional en desvios teoricos a 12, 24 y 36 meses:\n"); print(round(e36, 2))
cat(sprintf("Residuos de un paso en 2022-24 (modelo de toda la muestra): media %+.2f (sigma %.2f), z = %.2f, p = %.3f; positivos %d de %d\n", mean(res), sqrt(m0$sigma2), z_res, 2 * pnorm(-abs(z_res)), sum(res > 0), TF_H))
tabla_tex(data.frame(unlist(lapply(DES, function(d) c(d$etq, "", ""))), TD$anio, num(TD$sin_tasa, 1), num(TD$tasa, 1), num(TD$otros, 1), num(TD$observado, 1), num(pct(TD$observado), 1), num(TD$otros_ee, 2)),
          "16b_descomposicion.tex", rojo = celdas(c(3, 6), 3:6), "lccccccc",     # en puntos logaritmicos (100 ln) para que los aportes sumen exacto
          c("", "A\\~no", "Sin cambios de tasa", "Efecto de la tasa", "Otros factores", "Observado", "Observado, en \\%", "Otros factores, en desv\\'ios del pron\\'ostico"))

# ---- 5. La ganancia en todas las especificaciones ------------------------------------------------------------------
FZ <- read.csv("salidas/tablas/15_formas.csv"); forma1 <- function(k) FZ[FZ$control == k, ][which.min(FZ$SBC[FZ$control == k]), ]
zc <- list(mat = as.numeric(diff(100 * log(leer_extra("WPUSI012011", L$viv_sa)))), ofe = as.numeric(diff(100 * log(leer_extra("MSACSR", L$viv_sa)))))
ZR <- lapply(setNames(names(zc), names(zc)), function(k) tf_input(zc[[k]], forma1(k)$b, forma1(k)$s, 0, c(mat = "M", ofe = "O")[k]))
con <- function(Z) { n_esp <<- n_esp + 1; gan(tf_ajustar(Y, cbind(base_xr(), Z)), "D2") }
GA <- rbind(cbind(grupo = "Referencia", esp = "Toda la muestra (09)", data.frame(g = g0, ee = g0_ee)), cbind(grupo = "Referencia", esp = "Datos hasta 2021-12 (12)", gan(E0$m, "D2")),
            cbind(grupo = "Con controles (15)", esp = "+ materiales", con(ZR$mat)), cbind(grupo = "Con controles (15)", esp = "+ meses de oferta", con(ZR$ofe)),
            cbind(grupo = "Con controles (15)", esp = "+ materiales y oferta", con(cbind(ZR$mat, ZR$ofe))),
            cbind(grupo = "Por tramo", esp = c("1985-1999", "2000-2005", "2006-2011", "2012-2021", "2022-2026"), GP[, c("g", "ee")]),
            cbind(grupo = "Por fase de la tasa", esp = c("Fase de suba", "Fase de baja"), data.frame(g = c(FA$g_sube[1], FA$g_baja[1]), ee = c(FA$ee_sube[1], FA$ee_baja[1]))))
GA$lo <- GA$g - 1.96 * GA$ee; GA$hi <- GA$g + 1.96 * GA$ee; GA$pp <- GA$g * pp; rownames(GA) <- NULL
write.csv(GA, "salidas/tablas/16b_ganancias.csv", row.names = FALSE); cat("\n5. Ganancia en todas las especificaciones\n"); print(GA[, c("grupo", "esp", "g", "ee", "pp")], row.names = FALSE, digits = 3)

# ---- Valores para el texto -----------------------------------------------------------------------------------------
for (k in names(cortes)) { s <- c(PA = "Uno", PB = "Dos", PC = "Tres", PD = "Cuatro", PE = "Cinco")[k]; f <- GP[GP$tramo == k, ]
  guardar_valor(paste0("EstGan", s), num(f$g, 2)); guardar_valor(paste0("EstGanEE", s), num(f$ee, 2)); guardar_valor(paste0("EstPP", s), num(f$g * pp, 1)) }
guardar_valor("EstLR", num(LRp[["LR"]], 2)); guardar_valor("EstLRp", pval(LRp[["p"]])); guardar_valor("EstQuiebreT", num(t_q, 2)); guardar_valor("EstQuiebreP", pval(p_q))
guardar_valor("EstRecMin", num(min(REC$g), 2)); guardar_valor("EstRecMax", num(max(REC$g), 2)); guardar_valor("EstMovMin", num(min(MOV$g), 2)); guardar_valor("EstMovMax", num(max(MOV$g), 2))
guardar_valor("EstFaseSube", num(FA$g_sube[1], 2)); guardar_valor("EstFaseSubeEE", num(FA$ee_sube[1], 2)); guardar_valor("EstFaseBaja", num(FA$g_baja[1], 2)); guardar_valor("EstFaseBajaEE", num(FA$ee_baja[1], 2))
guardar_valor("EstFaseP", pval(FA$p[1])); guardar_valor("EstFasePSub", pval(FA$p[2])); guardar_valor("EstFaseDSBC", num(FA$dSBC[1], 1)); guardar_valor("EstFaseEpTres", num(tf_anual(EF)[3], 1))
guardar_valor("EstFaseRMSE", num(tf_error(EF)[["rmse"]], 1)); guardar_valor("EstFaseSesgo", num(tf_error(EF)[["sesgo"]], 1))
for (k in names(DES)) for (a in 2022:2024) { s <- paste0(c(ante = "Ante", post = "Post")[k], c(`2022` = "Uno", `2023` = "Dos", `2024` = "Tres")[as.character(a)]); f <- TD[TD$par == k & TD$anio == a, ]
  guardar_valor(paste0("DesSin", s), num(pct(f$sin_tasa), 1)); guardar_valor(paste0("DesTasa", s), num(pct(f$tasa), 1)); guardar_valor(paste0("DesOtros", s), num(pct(f$otros), 1)); guardar_valor(paste0("DesOtrosEE", s), num(f$otros_ee, 2)) }
guardar_valor("DesErrTreintaSeisAnte", num(e36[3, "ante"], 2)); guardar_valor("DesErrTreintaSeisPost", num(e36[3, "post"], 2)); guardar_valor("DesErrMaxAnte", num(max(abs((obs - DES$ante$cond) / DES$ante$ee)), 2))
guardar_valor("DesResMedia", num(mean(res), 2)); guardar_valor("DesResZ", num(z_res, 2)); guardar_valor("DesResP", pval(2 * pnorm(-abs(z_res)))); guardar_valor("DesResPositivos", sum(res > 0))
guardar_valor("EstGanMin", num(min(GA$g[GA$grupo != "Por fase de la tasa"]), 2)); guardar_valor("EstGanMax", num(max(GA$g[GA$grupo != "Por fase de la tasa"]), 2)); guardar_valor("EstModelos", n_esp)

# ---- Figuras ---------------------------------------------------------------------------------------------------------
# A. Estabilidad: trayectoria de la ganancia (recursiva y ventanas moviles) y ganancia por especificacion
g1 <- ggplot() + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) + geom_hline(yintercept = g0, colour = rojo, linetype = "22", linewidth = 0.4) +
  geom_ribbon(data = REC, aes(t, ymin = g - 1.96 * ee, ymax = g + 1.96 * ee), fill = col_viv, alpha = 0.13) + geom_line(data = REC, aes(t, g), colour = col_viv, linewidth = 0.9) +
  geom_errorbar(data = MOV, aes(t, ymin = g - 1.96 * ee, ymax = g + 1.96 * ee), colour = col_tasa, width = 0, linewidth = 0.4, alpha = 0.7) + geom_point(data = MOV, aes(t, g), colour = col_tasa, size = 1.7) +
  annotate("text", x = min(REC$t), y = min(MOV$g - 1.96 * MOV$ee) - 0.02, label = "azul: muestra que crece desde 1985\nnaranja: ventanas moviles de 15 anios", hjust = 0, vjust = 0, size = 3, colour = tinta2, lineheight = 1) +
  labs(title = "La ganancia de la tasa a lo largo del tiempo", subtitle = sprintf("estimada con datos hasta cada fecha, IC 95 %%; punteado rojo: toda la muestra (%.2f)", g0),
       x = "Ultimo dato usado en la estimacion", y = "ganancia g") + tema
GA$esp <- factor(GA$esp, levels = rev(GA$esp)); GA$grupo <- factor(GA$grupo, levels = unique(GA$grupo))
g2 <- ggplot(GA, aes(g, esp)) + geom_vline(xintercept = 0, colour = tinta2, linewidth = 0.3) + geom_vline(xintercept = g0, colour = rojo, linetype = "22", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "#7B8794", linewidth = 0.5) + geom_point(aes(colour = grupo == "Por fase de la tasa"), size = 2.3) +
  scale_colour_manual(values = c(`TRUE` = col_tasa, `FALSE` = col_viv), guide = "none") + facet_grid(grupo ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "La ganancia en cada especificacion", subtitle = "punto: estimacion; barra: IC 95 % (ee cond.); punteado rojo: referencia", x = "ganancia g (% de inicios por 1 % de tasa, largo plazo)", y = NULL) +
  tema + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, colour = tinta, face = "bold", size = 8), panel.spacing = unit(0.5, "lines"))
guardar_fig((g1 | g2) + plot_layout(widths = c(1, 1.15)), "16b_estabilidad.png", 14, 5.4)

# B. Descomposicion del episodio (parametros de toda la muestra): trayectorias y aportes anuales
t_f <- as.numeric(time(Y))[E0$fut]; d <- DES$post; rel <- function(v) pct(v - b21)
suave <- as.numeric(stats::filter(c(as.numeric(L$viv_sa)[fin + 0:1], obs), rep(1 / 3, 3), sides = 1))[-(1:2)]   # promedio de 3 meses del observado
CA <- data.frame(t = t_f, contra = rel(d$contra), cond = rel(d$cond), obs = rel(obs), suave = rel(suave))
g3 <- ggplot(CA, aes(t)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_ribbon(aes(ymin = pmin(suave, contra), ymax = contra), fill = rojo, alpha = 0.16) + geom_ribbon(aes(ymin = cond, ymax = pmax(pmin(suave, contra), cond)), fill = col_viv, alpha = 0.22) +
  geom_line(aes(y = obs), colour = tinta, linewidth = 0.3, alpha = 0.45) + geom_line(aes(y = suave), colour = tinta, linewidth = 0.9) +
  geom_line(aes(y = contra), colour = "#7B8794", linewidth = 0.8, linetype = "22") + geom_line(aes(y = cond), colour = rojo, linewidth = 0.9) +
  annotate("text", x = 2024.95, y = CA$contra[36] + 3, label = "sin cambios de tasa", hjust = 1, size = 3, colour = "#52606D") +
  annotate("text", x = 2024.95, y = CA$cond[36] - 3, label = "con la tasa observada", hjust = 1, size = 3, colour = rojo) +
  annotate("text", x = 2023.45, y = (CA$contra[18] + CA$suave[18]) / 2 + 1, label = "caida observada frente al escenario sin cambios de tasa", size = 3.1, colour = rojo, fontface = "bold") +
  annotate("text", x = 2023.95, y = (CA$suave[24] + CA$cond[24]) / 2 - 1, label = "parte del efecto de la tasa compensada por otros factores", size = 3.1, colour = col_viv, fontface = "bold") +
  annotate("text", x = 2022.05, y = max(CA$obs) + 2, label = "observado (linea gruesa: promedio de 3 meses)", hjust = 0, size = 3, colour = tinta) +
  labs(title = "Inicios 2022-2024: que parte explica la tasa", subtitle = "% respecto del promedio de 2021. Efecto de la tasa segun el modelo = banda roja + banda azul (de la linea punteada a la roja).\nParametros de toda la muestra, estado del modelo en diciembre de 2021", x = NULL, y = NULL) + tema
BA <- do.call(rbind, lapply(seq_len(nrow(TD)), function(i) data.frame(par = ifelse(TD$par[i] == "ante", "Parametros hasta 2021-12", "Parametros de toda la muestra"), anio = TD$anio[i],
                                                                     comp = c("Sin cambios de tasa", "Efecto de la tasa", "Otros factores"), v = c(TD$sin_tasa[i], TD$tasa[i], TD$otros[i]))))
BA$comp <- factor(BA$comp, levels = c("Sin cambios de tasa", "Otros factores", "Efecto de la tasa")); OB <- transform(TD, par = ifelse(par == "ante", "Parametros hasta 2021-12", "Parametros de toda la muestra"))
g4 <- ggplot(BA, aes(factor(anio), v)) + geom_col(aes(fill = comp), width = 0.62, colour = "white", linewidth = 0.6) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_point(data = OB, aes(factor(anio), observado), shape = 23, size = 3, fill = "white", colour = tinta, stroke = 0.9) +
  scale_fill_manual(values = c("Sin cambios de tasa" = "#B8C1CB", "Otros factores" = col_viv, "Efecto de la tasa" = rojo), name = NULL) + facet_wrap(~par) +
  labs(title = "Aportes a la variacion de los inicios respecto de 2021", subtitle = "puntos logaritmicos (100 ln), promedio de cada anio; rombo: observado = suma de los tres aportes", x = NULL, y = NULL) +
  tema + theme(legend.position = "bottom", strip.text = element_text(colour = tinta, face = "bold"))
guardar_fig(g3 / g4 + plot_layout(heights = c(1.15, 1)), "16b_descomposicion.png", 10.5, 9.2)
cat("Modelos estimados en este script:", n_esp, "\n")
escribir_valores("valores_16b.tex")
