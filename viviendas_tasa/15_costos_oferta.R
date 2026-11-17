# =============================================================================
# 15_costos_oferta.R -- El efecto de la tasa sobrevive al controlar por costos y oferta?
#   Preguntas: (1) la ganancia de la funcion de transferencia tasa -> inicios (09: -0,45) cambia al agregar
#   como segundo input el costo de los materiales o la oferta de viviendas nuevas en venta? (2) cuanto de la
#   sobreprediccion del episodio 2022-24 (script 12) se corrige con esos controles?
#   Modelo: Y_t = c + [delta_0 B^2 / (1 - phi B)] X_t + sum_j gamma_j Z_{t-j} + omega TC_t + N_t
#     Y = (1-B) 100 ln(inicios SA); X = (1-B) 100 ln(tasa); Z = control. La parte de la tasa queda FIJA en la
#     forma de 09 (b = 2, s = 0, r = 1, phi = 0,75; ruido ARMA(1,2)(1,1)12; COVID como TC con el delta de 08).
#   Herramientas y fuentes:
#     - Funcion de transferencia en 6 pasos (preblanqueo, FCC, forma por SBC, diagnostico): apunte parte 2,
#       pp. 8-21. Aca se aplica a un segundo input; el caso de varios inputs esta en Box, Jenkins, Reinsel y
#       Ljung (2015), Time Series Analysis: Forecasting and Control, 5.a ed., Wiley, cap. 12.
#     - ADF: parte 1 del apunte; valores p de MacKinnon (1996), Journal of Applied Econometrics, 11(6), 601-618.
#     - Contraste de Wald de los gamma (chi2) y ee de la ganancia por metodo delta con el perfil de phi (como en 09).
#     - Episodio: mismo procedimiento que 12 (estimacion hasta 2021-12 y pronostico condicional de 36 meses).
#   Especificaciones (preregistradas en el plan; se corren e informan TODAS):
#     1. materiales: Z = (1-B) 100 ln(WPUSI012011)      2. madera: Z = (1-B) 100 ln(WPU081)
#     3. oferta: Z = MSACSR (nivel o diferencia segun el ADF)      4. materiales + oferta
#     5. salarios de la construccion (CES2000000003), solo desde 2006 y solo dentro de la muestra
#   Forma de cada control: rezagos b_Z a b_Z + s_Z, con b_Z en {1, 2, 3} (minimo 1, para evitar simultaneidad),
#   s_Z en {0, 1, 2}, sin denominador, por SBC. Para el episodio se usa la forma elegida con la muestra completa
#   y se reestiman los parametros con datos hasta 2021-12 (igual que en 12 con la forma de la tasa).
#   CUIDADO (mal control): los meses de oferta son CONSECUENCIA de la tasa (tasa mas alta -> menos ventas -> mas
#   meses de oferta). Controlar por una variable del canal sesga el coeficiente de la tasa hacia cero: si |g| cae
#   al agregarla se lee como canal, no como confusion. Los precios de los materiales tambien responden a la demanda
#   de viviendas; por eso el rezago minimo 1 y la lectura prudente.
#   Salidas: 15_controles.tex, 15_episodio_controles.tex, 15_episodio_variantes.csv, 15_formas.csv;
#   eviews 15_modelo_materiales.png; figura 15_costos.png; valores_15.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 15_costos_oferta.R   (despues de 08)
# =============================================================================
source("R/funciones.R"); source("R/datos_extra.R"); source("R/tf.R")
L <- cargar_base(); I8 <- readRDS("salidas/tablas/08_intervencion.rds"); bajar_extra()
tf_prueba(L, I8$delta)
Y <- diff(L$viv_sa); X <- diff(L$tasa); nY <- length(Y); mes <- etiqueta_mes(Y); TCm <- cbind(COVID = tf_tc(L$viv_sa, I8$delta))
ult <- L$tasa_pp[length(L$tasa_pp)]; pp <- 100 * log((ult + 1) / ult)       # +1 pp sobre la tasa de 2026-08, como en 09
Xcf <- X; Xcf[(which(mes == TF_ORIGEN) + 1):nY] <- 0                       # tasa congelada desde 2022-01
base_xr <- function(Xv = X, phi = TF_PHI) cbind(tf_input(Xv, TF_B, 0, phi), TCm)
n_esp <- 0                                                                  # contador de modelos estimados (preregistro)

# ---- 0. Referencia (09) y series de control ---------------------------------------------------------------
m0 <- tf_ajustar(Y, base_xr()); W <- tf_serie(Y, TF_INI:nY)
g0 <- m0$coef[["D2"]] / (1 - TF_PHI); g0_ee <- tf_ee_ganancia(function(ph) tf_ajustar(Y, base_xr(X, ph)), TF_PHI, "D2")[["ee"]]
PBx <- tf_preblanqueo(X, Y - m0$coef[["COVID"]] * TCm[, 1])                 # alfa de la tasa, para el paso 6
Nhat <- Y - drop(base_xr() %*% m0$coef[c("D2", "COVID")])                   # output sin la parte de la tasa ni el COVID
niv <- list(mat = 100 * log(leer_extra("WPUSI012011", L$viv_sa)), mad = 100 * log(leer_extra("WPU081", L$viv_sa)), ofe = 100 * log(leer_extra("MSACSR", L$viv_sa)))
stopifnot(!anyNA(unlist(niv)))
A_ofe <- adf(niv$ofe, 12, "c", v = "LMSACSR"); ofe_dif <- A_ofe$p > ALFA    # con raiz unitaria, la oferta entra en diferencias
cat(sprintf("ADF de 100 ln(MSACSR), constante, 12 rezagos: tau = %.2f (p = %.3f) -> entra en %s\n", A_ofe$tau, A_ofe$p, if (ofe_dif) "diferencias" else "nivel"))
Z <- list(mat = as.numeric(diff(niv$mat)), mad = as.numeric(diff(niv$mad)), ofe = if (ofe_dif) as.numeric(diff(niv$ofe)) else as.numeric(niv$ofe - mean(niv$ofe))[-1])
etq <- c(mat = "Materiales de construcci\\'on", mad = "Madera y derivados", ofe = "Meses de oferta", mat_ofe = "Materiales + oferta", sal = "Salarios (desde 2006)")
etq_plano <- c(base = "Base (script 12)", mat = "+ materiales", mad = "+ madera", ofe = "+ meses de oferta", mat_ofe = "+ materiales y oferta")

# ---- 1. Pasos 1 a 3 para cada control: preblanqueo y FCC contra el output --------------------------------
PBz <- lapply(Z, function(z) tf_preblanqueo(ts(z, start = start(Y), frequency = S), Nhat, Pmax = 1, Qmax = 1))
for (k in names(PBz)) cat(sprintf("%s: filtro %s | FCC k = 1..12: p = %.3f; k < 0: p = %.3f | fuera de banda en k = %s\n", k, gsub("\\\\ensuremath\\{_\\{12\\}\\}", "12", PBz[[k]]$modelo),
                                  PBz[[k]]$Qpos1["p"], PBz[[k]]$Qneg["p"], paste(PBz[[k]]$C$k[PBz[[k]]$C$sig & abs(PBz[[k]]$C$k) <= 12], collapse = " ")))

# ---- 2. Forma de cada control por SBC y estimacion conjunta ------------------------------------------------
elegir_forma <- function(z, pref, ini = TF_INI) { G <- expand.grid(b = 1:3, s = 0:2)
  G$SBC <- sapply(seq_len(nrow(G)), function(i) { m <- tf_ajustar(Y, cbind(base_xr(), tf_input(z, G$b[i], G$s[i], 0, pref)), ini = ini); if (is.null(m)) NA else BIC(m) })
  n_esp <<- n_esp + nrow(G); G[order(G$SBC), ] }
FORMAS <- lapply(setNames(names(Z), names(Z)), function(k) elegir_forma(Z[[k]], c(mat = "M", mad = "L", ofe = "O")[k]))
ZR <- lapply(setNames(names(Z), names(Z)), function(k) tf_input(Z[[k]], FORMAS[[k]]$b[1], FORMAS[[k]]$s[1], 0, c(mat = "M", mad = "L", ofe = "O")[k]))
ZR$mat_ofe <- cbind(ZR$mat, ZR$ofe)
write.csv(do.call(rbind, lapply(names(FORMAS), function(k) cbind(control = k, FORMAS[[k]]))), "salidas/tablas/15_formas.csv", row.names = FALSE)
evaluar <- function(ZRm, ini = TF_INI, ref = m0) {
  m <- tf_ajustar(Y, cbind(base_xr(), ZRm), ini = ini); n_esp <<- n_esp + 1; nz <- colnames(ZRm); gam <- m$coef[nz]; Vg <- m$var.coef[nz, nz, drop = FALSE]
  ee <- tf_ee_ganancia(function(ph) tf_ajustar(Y, cbind(base_xr(X, ph), ZRm), ini = ini), TF_PHI, "D2")
  wald <- drop(t(gam) %*% solve(Vg) %*% gam)
  list(m = m, g = m$coef[["D2"]] / (1 - TF_PHI), g_ee = if (is.na(ee[["ee"]])) ee[["ee_cond"]] else ee[["ee"]], sg = sum(gam), sg_ee = sqrt(sum(Vg)), wald_p = pchisq(wald, length(nz), lower.tail = FALSE),
       dSBC = BIC(m) - BIC(ref), LB24 = diagnostico(m, sum(TF_RUIDO))$LB24, cc = tf_cc_residuos(m, tf_serie(Y, ini:nY), PBx$alfa, 2), k = length(nz),
       forma = paste(sub("^[A-Z]", "", nz), collapse = ", "),
       sg_txt = paste(sapply(unique(substr(nz, 1, 1)), function(l) { i <- nz[substr(nz, 1, 1) == l]; sprintf("%s (%s)", num(sum(gam[i]), 2), num(sqrt(sum(Vg[i, i])), 2)) }), collapse = "; "))
}
R <- lapply(ZR, evaluar)
# 5. Salarios: solo desde 2006, contra la referencia reestimada en la misma submuestra
sal <- 100 * log(leer_extra("CES2000000003", L$viv_sa)); zs <- as.numeric(diff(sal)); ini_s <- which(!is.na(zs))[1] + 6; zs[is.na(zs)] <- 0
m0s <- tf_ajustar(Y, base_xr(), ini = ini_s); Fs <- elegir_forma(zs, "S", ini_s); Rs <- NULL
if (!is.null(m0s)) { ZRs <- tf_input(zs, Fs$b[1], Fs$s[1], 0, "S"); Rs <- tryCatch(evaluar(ZRs, ini_s, m0s), error = function(e) NULL)
  g0s <- m0s$coef[["D2"]] / (1 - TF_PHI); g0s_ee <- sqrt(m0s$var.coef["D2", "D2"]) / (1 - TF_PHI) }

fila <- function(e, r, f = NULL) data.frame(e, if (is.null(f)) "--" else sprintf("$b_Z=%d$, $s_Z=%d$", f$b[1], f$s[1]), sprintf("%s (%s)", num(r$g, 2), num(r$g_ee, 2)), num(r$g * pp, 1),
                                           r$sg_txt, pval(r$wald_p), num(r$dSBC, 1), num(r$LB24, 2), pval(r$cc[["pos"]]))
CU <- rbind(data.frame("Sin controles (09)", "--", sprintf("%s (%s)", num(g0, 2), num(g0_ee, 2)), num(g0 * pp, 1), "--", "--", "0", num(diagnostico(m0, 5)$LB24, 2), pval(tf_cc_residuos(m0, W, PBx$alfa, 2)[["pos"]])) |> setNames(paste0("V", 1:9)),
            setNames(fila(etq["mat"], R$mat, FORMAS$mat), paste0("V", 1:9)), setNames(fila(etq["mad"], R$mad, FORMAS$mad), paste0("V", 1:9)),
            setNames(fila(etq["ofe"], R$ofe, FORMAS$ofe), paste0("V", 1:9)), setNames(fila(etq["mat_ofe"], R$mat_ofe), paste0("V", 1:9)))
if (!is.null(Rs)) CU <- rbind(CU, data.frame("Sin controles, desde 2006", "--", sprintf("%s (%s)", num(g0s, 2), num(g0s_ee, 2)), num(g0s * pp, 1), "--", "--", "0", num(diagnostico(m0s, 5)$LB24, 2), "--") |> setNames(paste0("V", 1:9)),
                              setNames(fila(etq["sal"], Rs, Fs), paste0("V", 1:9)))
tabla_tex(CU, "15_controles.tex", rojo = celdas(1:5, 3:4), "llccccccc",
          c("Control", "Rezagos del control", "Ganancia de la tasa $g$ (ee)", "\\% por +1 pp", "$\\sum\\hat\\gamma$ (ee)", "Wald de $\\gamma$: $p$", "$\\Delta$SBC", "LB24 $p$", "$\\hat\\varepsilon$ vs $\\alpha$, $k\\geq0$"))
cat("\nGanancia de la tasa con cada control (referencia:", round(g0, 3), "ee", round(g0_ee, 3), ")\n")
print(data.frame(esp = names(R), forma = sapply(R, `[[`, "forma"), g = sapply(R, `[[`, "g"), g_ee = sapply(R, `[[`, "g_ee"), sum_gamma = sapply(R, `[[`, "sg"), sg_ee = sapply(R, `[[`, "sg_ee"),
                 wald_p = sapply(R, `[[`, "wald_p"), dSBC = sapply(R, `[[`, "dSBC"), LB24 = sapply(R, `[[`, "LB24")), row.names = FALSE, digits = 3)
if (!is.null(Rs)) cat(sprintf("Salarios (desde %s, n = %d): g = %.3f (ee %.3f) contra %.3f sin control en la misma muestra; sum gamma = %.2f (ee %.2f), Wald p = %.3f, dSBC = %.1f\n",
                              mes[ini_s], length(residuals(Rs$m)), Rs$g, Rs$g_ee, g0s, Rs$sg, Rs$sg_ee, Rs$wald_p, Rs$dSBC))

# ---- 3. Episodio 2022-24 con controles (estimacion hasta 2021-12) -------------------------------------------
E0 <- tf_episodio(Y, L$viv_sa, base_xr(X), base_xr(Xcf))
EP <- lapply(ZR, function(ZRm) { n_esp <<- n_esp + 1; fin <- which(mes == TF_ORIGEN)
  E <- tf_episodio(Y, L$viv_sa, cbind(base_xr(X), ZRm), cbind(base_xr(Xcf), ZRm))
  stopifnot(ofe_dif); Z0 <- ZRm                                                               # controles congelados desde 2022-01: Z_s = 0 para s >= 2022-01 (son diferencias)
  for (j in seq_len(ncol(ZRm))) { lagj <- as.integer(sub("^[A-Z]", "", colnames(ZRm)[j])); Z0[(fin + 1 + lagj):nY, j] <- 0 }
  E$ef_control <- E$cond - (as.numeric(L$viv_sa)[fin + 1] + cumsum(tf_predecir(E$m, cbind(base_xr(X), Z0)[E$fut, , drop = FALSE]))); E })
res_ep <- function(E, nombre) { a <- tf_anual(E); e <- tf_error(E); ef <- tf_efecto_anual(E)
  data.frame(variante = nombre, p22 = a[1], p23 = a[2], p24 = a[3], ef_tasa24 = ef[3], ef_ctrl24 = if (is.null(E$ef_control)) NA else pct(mean(E$ef_control[E$anio == 2024])),
             rmse = e[["rmse"]], sesgo = e[["sesgo"]], g_sub = E$m$coef[["D2"]] / (1 - TF_PHI)) }
RE <- rbind(res_ep(E0, "base"), do.call(rbind, lapply(names(EP), function(k) res_ep(EP[[k]], k))))
obs_a <- tf_anual(E0, E0$obs); print(round(obs_a, 1)); print(RE, row.names = FALSE, digits = 3)
tabla_tex(rbind(data.frame(v = "Observado", a = num(obs_a[1], 1), b = num(obs_a[2], 1), c = num(obs_a[3], 1), d = "--", e = "--", f = "--", g = "--"),
                data.frame(v = unname(sub("\\+", "$+$", etq_plano[RE$variante])), a = num(RE$p22, 1), b = num(RE$p23, 1), c = num(RE$p24, 1), d = num(RE$ef_tasa24, 1),
                           e = ifelse(is.na(RE$ef_ctrl24), "--", num(RE$ef_ctrl24, 1)), f = num(RE$rmse, 1), g = num(RE$sesgo, 1))),
          "15_episodio_controles.tex", rojo = celdas(1:2, 2:4), "lccccccc",
          c("Modelo (datos hasta 2021-12)", "2022", "2023", "2024", "Efecto de la tasa, 2024", "Efecto del control, 2024", "RMSE 36 meses", "Sesgo medio"))
write.csv(do.call(rbind, lapply(names(EP), function(k) tf_variante(5, k, etq_plano[[k]], EP[[k]], k_extra = R[[k]]$k, dSBC = R[[k]]$dSBC,
                                                                 nota = if (grepl("ofe", k)) "la oferta es parte del canal de la tasa" else ""))), "salidas/tablas/15_episodio_variantes.csv", row.names = FALSE)

# ---- Valores para el texto ----------------------------------------------------------------------------------
guardar_valor("CtlGanBase", num(g0, 2)); guardar_valor("CtlGanBaseEE", num(g0_ee, 2)); guardar_valor("CtlOfertaADF", num(A_ofe$tau, 2)); guardar_valor("CtlOfertaADFp", pval(A_ofe$p))
guardar_valor("CtlOfertaForma", if (ofe_dif) "diferencias" else "nivel"); guardar_valor("CtlModelos", n_esp)
for (k in names(R)) { s <- c(mat = "Mat", mad = "Mad", ofe = "Ofe", mat_ofe = "MatOfe")[k]; r <- R[[k]]; e <- RE[RE$variante == k, ]
  guardar_valor(paste0("Ctl", s, "Gan"), num(r$g, 2)); guardar_valor(paste0("Ctl", s, "GanEE"), num(r$g_ee, 2)); guardar_valor(paste0("Ctl", s, "PP"), num(r$g * pp, 1))
  guardar_valor(paste0("Ctl", s, "Gamma"), num(r$sg, 2)); guardar_valor(paste0("Ctl", s, "GammaEE"), num(r$sg_ee, 2)); guardar_valor(paste0("Ctl", s, "Wald"), pval(r$wald_p))
  guardar_valor(paste0("Ctl", s, "DSBC"), num(r$dSBC, 1)); guardar_valor(paste0("Ctl", s, "CambioEE"), num(abs(r$g - g0) / g0_ee, 2))
  guardar_valor(paste0("Ctl", s, "EpTres"), num(e$p24, 1)); guardar_valor(paste0("Ctl", s, "RMSE"), num(e$rmse, 1)); guardar_valor(paste0("Ctl", s, "Sesgo"), num(e$sesgo, 1)) }
guardar_valor("CtlBaseRMSE", num(RE$rmse[1], 1)); guardar_valor("CtlBaseSesgo", num(RE$sesgo[1], 1)); guardar_valor("CtlBaseEpTres", num(RE$p24[1], 1)); guardar_valor("CtlObsTres", num(obs_a[3], 1))
if (!is.null(Rs)) { guardar_valor("CtlSalGan", num(Rs$g, 2)); guardar_valor("CtlSalGanBase", num(g0s, 2)); guardar_valor("CtlSalWald", pval(Rs$wald_p)); guardar_valor("CtlSalDesde", mes[ini_s]) }

# ---- Salida EViews y figura -----------------------------------------------------------------------------------
guardar_salida(salida_arma(R$mat$m, W, "D(LHOUST_SA)", "Funcion de transferencia con materiales de construccion como segundo input"), "15_modelo_materiales.png")
t_niv <- as.numeric(time(L$viv_sa)); v <- t_niv >= 2015
IDX <- rbind(data.frame(t = t_niv[v], y = 100 * exp((niv$mat[v] - niv$mat[which(t_niv >= 2019)[1]]) / 100), s = "Materiales de construccion"),
             data.frame(t = t_niv[v], y = 100 * exp((niv$mad[v] - niv$mad[which(t_niv >= 2019)[1]]) / 100), s = "Madera y derivados"))
g1 <- ggplot(IDX, aes(t, y, colour = s)) + geom_vline(xintercept = 2022, colour = tinta2, linetype = "22", linewidth = 0.3) + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("Materiales de construccion" = col_viv, "Madera y derivados" = col_tasa), name = NULL) +
  scale_x_continuous(breaks = seq(2015, 2025, 2)) +
  labs(title = "Precios al productor", subtitle = "indice, enero de 2019 = 100", x = NULL, y = NULL) + tema + theme(legend.position = "bottom")
g2 <- ggplot(data.frame(t = t_niv[v], y = exp(niv$ofe[v] / 100)), aes(t, y)) + geom_vline(xintercept = 2022, colour = tinta2, linetype = "22", linewidth = 0.3) + geom_line(colour = tinta, linewidth = 0.7) +
  scale_x_continuous(breaks = seq(2015, 2025, 2)) +
  labs(title = "Meses de oferta de viviendas nuevas", subtitle = "stock en venta / ventas del mes", x = NULL, y = NULL) + tema
anual <- function(x) 12 * exp(x / 100); t_f <- as.numeric(time(Y))[E0$fut]
hist <- data.frame(t = t_niv, v = anual(as.numeric(L$viv_sa))); hist <- hist[hist$t >= 2019 & hist$t < 2025, ]
PR <- rbind(data.frame(t = t_f, v = anual(E0$cond), s = "Base (solo tasa)"), data.frame(t = t_f, v = anual(EP$mat$cond), s = "+ materiales"),
            data.frame(t = t_f, v = anual(EP$ofe$cond), s = "+ meses de oferta"), data.frame(t = t_f, v = anual(EP$mat_ofe$cond), s = "+ materiales y oferta"))
g3 <- ggplot() + geom_vline(xintercept = t_niv[E0$fin + 1], colour = tinta2, linetype = "22", linewidth = 0.3) + geom_line(data = hist, aes(t, v), colour = tinta, linewidth = 0.55) +
  geom_line(data = PR, aes(t, v, colour = s), linewidth = 0.8) +
  scale_colour_manual(values = c("Base (solo tasa)" = rojo, "+ materiales" = col_viv, "+ meses de oferta" = col_tasa, "+ materiales y oferta" = "#6B4EA0"), name = NULL) +
  annotate("text", x = 2019.1, y = max(hist$v) * 0.99, label = "observado", hjust = 0, size = 3, colour = tinta) +
  labs(title = "Inicios: pronostico condicional desde diciembre de 2021, con y sin controles", subtitle = "miles de unidades, tasa anual desestacionalizada; tasa y controles observados", x = NULL, y = NULL) + tema + theme(legend.position = "bottom")
guardar_fig((g1 | g2) / g3 + plot_layout(heights = c(1, 1.5)), "15_costos.png", 11, 8)
cat("Modelos estimados en este script (formas de los controles + especificaciones + episodio):", n_esp, "\n")
escribir_valores("valores_15.tex")
