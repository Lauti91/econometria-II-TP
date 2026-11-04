# =============================================================================
# 09_transferencia.R -- Funcion de transferencia: tasa hipotecaria (input) -> inicios (output)
#   Apunte parte 2, pp. 8-21; Enders cap. 5, sec. 2. Procedimiento de 6 pasos del apunte (p. 19):
#     1. ARMA del input X = (1-B) 100 ln(tasa)  ->  innovaciones alfa_t (X preblanqueado)
#     2. Se filtra el output Y = (1-B) 100 ln(inicios desestacionalizados) con el MISMO filtro -> beta_t
#     3. FCC r(k) = corr(alfa_t, beta_{t+k}), bandas +-2/raiz(n); pesos v_k = (sd(beta)/sd(alfa)) r(k)
#     4-5. Se estima Y_t = c + [delta(B) B^b / (1 - phi B)^r] X_t + omega TC_t + N_t, con N_t ARMA (ruido);
#          (b, s, r) por SBC y el ruido se re-identifica sobre el modelo elegido
#     6. Diagnostico: residuos ruido blanco y SIN correlacion cruzada con alfa (si no, revisar delta(B))
#   Se usa la serie desestacionalizada para que ambas variables tengan la misma diferenciacion (1-B):
#   asi v(B) se lee como elasticidad en niveles (apunte p. 21, caso c). La estacionalidad no es parte
#   de la relacion. El COVID entra con la forma elegida en 08 (cambio transitorio, delta de 08).
#   Exogeneidad (apunte p. 8: E(X_t eps_{t-j}) = 0): se contrasta con la FCC en k < 0. Si los inicios
#   anticipan a la tasa hay retroalimentacion y el MFT no es el modelo adecuado -> VAR (apunte p. 22, script 10).
#   Variante: el bono del Tesoro a 10 anios como input, mas defendible como exogeno.
# Ejecutar desde viviendas_tasa/:  Rscript 09_transferencia.R   (despues de 08)
# =============================================================================
source("R/funciones.R")
L <- cargar_base()
I8 <- readRDS("salidas/tablas/08_intervencion.rds")
Y <- diff(L$viv_sa)                                                        # output: % mensual de inicios (SA)
ENT <- list(tasa = diff(L$tasa), bono = diff(100 * log(L$bono)))          # inputs: % mensual de cada tasa
TC <- diff(ts(as.numeric(stats::filter(indicadora(L$viv_sa, 2020, 4), I8$delta, method = "recursive")), start = start(L$viv_sa), frequency = S))
TCm <- cbind(COVID = as.numeric(TC))
K <- 24; n <- length(Y); banda <- 2 / sqrt(n)
ult <- c(tasa = L$tasa_pp[length(L$tasa_pp)], bono = L$bono[length(L$bono)])
pp <- 100 * log((ult + 1) / ult)                                         # +1 pp sobre cada tasa en 2026-08, en % de esa tasa

# ---- Ruido univariado del output (punto de partida del paso 5) --------------------------------------
RY <- rejilla(Y, TCm, TRUE, 3, 2, 1, 1); RY <- RY[RY$admisible, ]
ruido0 <- unlist(RY[1, c("p", "q", "P", "Q")])
mY <- ajustar(Y, ruido0[1], ruido0[2], ruido0[3], ruido0[4], TCm, TRUE)
Ylimpia <- Y - mY$coef["COVID"] * TC                                       # sin el COVID, para preblanquear
cat("Ruido univariado de Y:", RY$modelo[1], "\n")

# ---- Pasos 1-3: preblanqueo y FCC ----------------------------------------------------------------------
preblanqueo <- function(X) {
  R <- rbind(transform(rejilla(X, NULL, TRUE, 3, 2, 0, 0), media = TRUE), transform(rejilla(X, NULL, FALSE, 3, 2, 0, 0), media = FALSE))
  R <- R[R$admisible, ]; R <- R[order(R$SBC), ]; g <- R[1, ]
  mX <- ajustar(X, g$p, g$q, 0, 0, NULL, g$media)
  # El mismo filtro phi_X(B)/theta_X(B) sobre X y sobre Y (la media se resta antes de filtrar)
  filtro <- function(z) { z <- z - mean(z); ar <- mX$model$phi; ma <- mX$model$theta
    if (length(ar)) z <- stats::filter(z, c(1, -ar), sides = 1)                      # phi_X(B) z
    if (length(ma)) z <- stats::filter(ifelse(is.na(z), 0, z), -ma, method = "recursive") # / theta_X(B), signo de R: 1 + theta B
    ts(as.numeric(z), start = start(X), frequency = S) }
  a <- filtro(X); b <- filtro(Ylimpia)
  ok <- !is.na(a) & !is.na(b) & seq_along(a) > 12                       # descarta el arranque del filtro (condicion inicial 0)
  a <- ts(a[ok], start = time(X)[which(ok)[1]], frequency = S); b <- b[ok]   # alfa como ts: conserva las fechas
  stopifnot(all(diff(which(ok)) == 1))                                     # tramo contiguo
  r <- ccf(b, as.numeric(a), lag.max = K, plot = FALSE)                                # r(k) = corr(beta_{t+k}, alfa_t): k > 0, el input anticipa
  C <- data.frame(k = as.numeric(r$lag), r = as.numeric(r$acf)); C$sig <- abs(C$r) > banda
  C$v <- ifelse(C$k >= 0, sd(b) / sd(a) * C$r, NA)                         # pesos de respuesta al impulso (apunte p. 19)
  Qk <- function(lados, M = 12) { d <- C[if (lados == "pos") C$k >= 0 & C$k <= M else C$k < 0 & C$k >= -M, ]
    Q <- length(a) * (length(a) + 2) * sum(d$r^2 / (length(a) - abs(d$k))); c(Q = Q, gl = nrow(d), p = pchisq(Q, nrow(d), lower.tail = FALSE)) }
  list(modelo = g$modelo, media = g$media, mX = mX, alfa = a, C = C, Qpos = Qk("pos"), Qneg = Qk("neg"))
}
PB <- lapply(ENT, preblanqueo)
for (e in names(PB)) { cat(e, "-> filtro", PB[[e]]$modelo, "| sig:", paste(PB[[e]]$C$k[PB[[e]]$C$sig], collapse = " "), "\n") }

# ---- Pasos 4-5: estimacion del MFT -------------------------------------------------------------------
# Regresores del numerador: X_{t-b}, ..., X_{t-b-s}, cada uno filtrado por 1/(1 - phi B) si r = 1
regresores <- function(X, b, s, phi = 0) {
  stopifnot(all(time(X) == time(Y)))
  xf <- as.numeric(if (phi != 0) stats::filter(X, phi, method = "recursive") else X); nY <- length(Y)
  M <- sapply(b:(b + s), function(j) c(rep(0, j), xf[seq_len(nY - j)]))     # columna j: X_{t-j}
  matrix(M, nY, dimnames = list(NULL, sprintf("D%d", b:(b + s))))
}
# Para que todos los modelos usen las mismas fechas, la estimacion empieza despues del rezago maximo
INI <- 7                                                                   # b + s <= 6
estimar <- function(X, b, s, r, phi, ruido) {
  xr <- cbind(regresores(X, b, s, if (r == 1) phi else 0), TCm)[INI:length(Y), , drop = FALSE]
  W <- ts(as.numeric(Y)[INI:length(Y)], start = time(Y)[INI], frequency = S)
  ajustar(W, ruido[1], ruido[2], ruido[3], ruido[4], xr, TRUE)
}
# Error estandar de la ganancia g = sum(delta) / (1 - phi) contando la incertidumbre de phi.
# phi se elige por perfil de verosimilitud (arima no lo estima), asi que var.coef solo da Var(delta | phi).
# Con un perfil aproximadamente cuadratico: Var(phi) = -1 / l''(phi); la pendiente del camino delta_hat(phi)
# da Cov(delta, phi) = pend * Var(phi), y Var(delta) = Var(delta | phi) + pend^2 Var(phi). Despues, metodo delta.
ee_ganancia <- function(ajuste_phi, phi, nd, h = 0.02) {
  M <- lapply(c(phi - h, phi, phi + h), ajuste_phi)
  if (any(sapply(M, is.null))) return(c(ee = NA, ee_cond = NA, ee_phi = NA))
  ll <- sapply(M, function(m) as.numeric(logLik(m))); d <- sapply(M, function(m) sum(m$coef[nd]))
  curv <- (ll[3] - 2 * ll[2] + ll[1]) / h^2; v_cond <- sum(M[[2]]$var.coef[nd, nd])
  if (curv >= 0) return(c(ee = NA, ee_cond = sqrt(v_cond) / (1 - phi), ee_phi = NA))
  v_phi <- -1 / curv; pend <- (d[3] - d[1]) / (2 * h)
  V <- matrix(c(v_cond + pend^2 * v_phi, pend * v_phi, pend * v_phi, v_phi), 2)
  gr <- c(1 / (1 - phi), d[2] / (1 - phi)^2)
  c(ee = sqrt(drop(t(gr) %*% V %*% gr)), ee_cond = sqrt(v_cond) / (1 - phi), ee_phi = sqrt(v_phi))
}
grilla_phi <- seq(-0.8, 0.8, by = 0.2)                                     # grilla gruesa; despues se refina
seleccion <- function(X) {
  out <- list()
  for (b in 0:3) for (s in 0:3) if (b + s <= 6) for (r in 0:1) {
    phis <- if (r == 0) 0 else grilla_phi
    ms <- lapply(phis, function(ph) estimar(X, b, s, r, ph, ruido0)); ll <- sapply(ms, function(m) if (is.null(m)) -Inf else as.numeric(logLik(m)))
    i <- which.max(ll); m <- ms[[i]]; if (is.null(m)) next
    kk <- length(m$coef) + r; nn <- length(residuals(m))
    out[[length(out) + 1]] <- data.frame(b = b, s = s, r = r, phi = phis[i], SBC = -2 * ll[i] + kk * log(nn), logLik = ll[i], row.names = NULL)
  }
  G <- do.call(rbind, out); G$dSBC <- G$SBC - min(G$SBC); G[order(G$SBC), ]
}
# Sin input (solo ruido + COVID): referencia para ver si la tasa aporta
m_nulo <- { xr <- TCm[INI:length(Y), , drop = FALSE]; W <- ts(as.numeric(Y)[INI:length(Y)], start = time(Y)[INI], frequency = S); ajustar(W, ruido0[1], ruido0[2], ruido0[3], ruido0[4], xr, TRUE) }
SBC_nulo <- BIC(m_nulo)

ajuste_final <- function(X, nombre) {
  G <- seleccion(X); g <- G[1, ]
  if (g$r == 1) {                                                          # refina phi alrededor del ganador (paso 0,05)
    fina <- seq(max(-0.95, g$phi - 0.15), min(0.95, g$phi + 0.15), by = 0.05)
    ll <- sapply(fina, function(ph) { m <- estimar(X, g$b, g$s, 1, ph, ruido0); if (is.null(m)) -Inf else as.numeric(logLik(m)) })
    g$phi <- fina[which.max(ll)] }
  xr <- cbind(regresores(X, g$b, g$s, if (g$r == 1) g$phi else 0), TCm)[INI:length(Y), , drop = FALSE]
  W <- ts(as.numeric(Y)[INI:length(Y)], start = time(Y)[INI], frequency = S)
  # Paso 5: re-identificacion del ruido con el input ya incluido
  R <- rejilla(W, xr, TRUE, 3, 2, 1, 1); R <- R[R$admisible, ]; ru <- unlist(R[1, c("p", "q", "P", "Q")])
  m <- ajustar(W, ru[1], ru[2], ru[3], ru[4], xr, TRUE)
  # Paso 6: diagnostico. Residuos vs innovaciones del input (mismas fechas)
  e <- ts(as.numeric(residuals(m)), start = start(W), frequency = S)
  a <- ts(c(rep(NA, length(Y) - length(PB[[nombre]]$alfa)), PB[[nombre]]$alfa), start = start(Y), frequency = S)   # alfa sin las 12 primeras
  Z <- na.omit(ts.intersect(e, a)); ce <- ccf(as.numeric(Z[, "e"]), as.numeric(Z[, "a"]), lag.max = K, plot = FALSE)
  CE <- data.frame(k = as.numeric(ce$lag), r = as.numeric(ce$acf))
  Qe <- function(pos) { d <- CE[if (pos) CE$k >= 0 & CE$k <= 12 else CE$k < 0 & CE$k >= -12, ]; Q <- nrow(Z) * (nrow(Z) + 2) * sum(d$r^2 / (nrow(Z) - abs(d$k)))
    gl <- if (pos) nrow(d) - (g$s + 1 + g$r) else nrow(d)                # Box-Jenkins: se descuentan r+s+1 solo para k >= 0
    pchisq(Q, gl, lower.tail = FALSE) }
  # Pesos v_j implicados (recursion del apunte p. 20) y ganancia de estado estable g = sum v_j (p. 8)
  dl <- m$coef[sprintf("D%d", g$b:(g$b + g$s))]; V <- m$var.coef[names(dl), names(dl)]
  vj <- sapply(0:K, function(j) { if (j < g$b) return(0); sum(sapply(g$b:min(j, g$b + g$s), function(i) dl[sprintf("D%d", i)] * (if (g$r == 1) g$phi^(j - i) else as.numeric(j == i)))) })
  gan <- sum(dl) / (1 - (if (g$r == 1) g$phi else 0)); gan_ee <- sqrt(sum(V)) / (1 - (if (g$r == 1) g$phi else 0))
  gan_ee_cond <- gan_ee; phi_ee <- NA
  if (g$r == 1) {                                                          # ee que incorpora la incertidumbre de phi
    E <- ee_ganancia(function(ph) ajustar(W, ru[1], ru[2], ru[3], ru[4], cbind(regresores(X, g$b, g$s, ph), TCm)[INI:length(Y), , drop = FALSE], TRUE),
                     g$phi, names(dl))
    if (!is.na(E["ee"])) { gan_ee <- unname(E["ee"]); phi_ee <- unname(E["ee_phi"]) } else warning("perfil de phi no concavo: ee condicional") }
  D <- diagnostico(m, sum(ru))
  vj <- unname(vj)
  list(G = G, g = g, ruido = R$modelo[1], m = m, W = W, D = D, vj = vj, gan = gan, gan_ee = gan_ee, gan_ee_cond = gan_ee_cond, phi_ee = phi_ee, CE = CE,
       p_cc_pos = Qe(TRUE), p_cc_neg = Qe(FALSE), dSBC_nulo = BIC(m) - SBC_nulo)
}
TF <- lapply(setNames(names(ENT), names(ENT)), function(e) ajuste_final(ENT[[e]], e))

# ---- Cuadros -------------------------------------------------------------------------------------------
etq_in <- c(tasa = "Tasa hipotecaria 30 a\\~nos", bono = "Bono del Tesoro 10 a\\~nos")
CU <- do.call(rbind, lapply(names(TF), function(e) { t <- TF[[e]]; P <- PB[[e]]
  data.frame(etq_in[e], P$modelo, pval(P$Qpos["p"]), pval(P$Qneg["p"]),
             sprintf("$b=%d$, $s=%d$, $r=%d$%s", t$g$b, t$g$s, t$g$r, if (t$g$r == 1) sprintf(" ($\\phi=%s$)", num(t$g$phi, 2)) else ""),
             t$ruido, sprintf("%s (%s)", num(t$gan, 2), num(t$gan_ee, 2)), num(t$gan * pp[e], 1), num(t$D$LB24, 2), pval(t$p_cc_pos), pval(t$p_cc_neg), num(t$dSBC_nulo, 1)) }))
tabla_tex(CU, "09_transferencia.tex", rojo = celdas(1, 7:8), "lllllllccccc",
          c("Input", "Filtro del input", "FCC $k\\geq0$: $p$", "FCC $k<0$: $p$", "Funci\\'on de transferencia", "Ruido", "Ganancia $g$ (ee)", "\\% inicios por +1 pp del input",
            "LB24 $p$", "$\\hat\\varepsilon$ vs $\\alpha$, $k\\geq0$", "$\\hat\\varepsilon$ vs $\\alpha$, $k<0$", "SBC vs sin input"))
for (e in names(TF)) { t <- TF[[e]]; P <- PB[[e]]; s <- tools::toTitleCase(e)
  guardar_valor(paste0("TF", s, "Filtro"), P$modelo); guardar_valor(paste0("TF", s, "Ruido"), t$ruido)
  guardar_valor(paste0("TF", s, "Pos"), pval(P$Qpos["p"])); guardar_valor(paste0("TF", s, "Neg"), pval(P$Qneg["p"]))
  guardar_valor(paste0("TF", s, "Gan"), num(t$gan, 2)); guardar_valor(paste0("TF", s, "GanEE"), num(t$gan_ee, 2)); guardar_valor(paste0("TF", s, "PP"), num(t$gan * pp[e], 1))
  guardar_valor(paste0("TF", s, "Esp"), sprintf("$b=%d$, $s=%d$, $r=%d$", t$g$b, t$g$s, t$g$r)); guardar_valor(paste0("TF", s, "DSBC"), num(t$dSBC_nulo, 1))
  guardar_valor(paste0("TF", s, "Sig"), { c_ <- P$C[P$C$sig, ]; if (nrow(c_)) paste(sprintf("$k=%d$ (%s)", c_$k, num(c_$r, 2)), collapse = ", ") else "ninguno" })
  guardar_valor(paste0("TF", s, "CCneg"), pval(t$p_cc_neg))
  guardar_valor(paste0("TF", s, "GanEEcond"), num(t$gan_ee_cond, 2)); if (!is.na(t$phi_ee)) guardar_valor(paste0("TF", s, "PhiEE"), num(t$phi_ee, 2)) }
guardar_valor("TFBanda", num(banda, 3)); guardar_valor("TFRuidoUni", RY$modelo[1]); guardar_valor("TFdeltaTC", num(I8$delta, 2))
G <- head(TF$tasa$G, 6)
tabla_tex(data.frame(G$b, G$s, G$r, ifelse(G$r == 1, num(G$phi, 1), "--"), num(G$SBC, 1), num(G$dSBC, 1)), "09_seleccion_tasa.tex", rojo = celdas(1, 1:6), "cccccc",
          c("Retardo $b$", "$s$", "$r$", "$\\phi$ (grilla 0,2)", "SBC (ruido inicial)", "$\\Delta$SBC"))
write.csv(do.call(rbind, lapply(names(PB), function(e) cbind(PB[[e]]$C, input = e))), "salidas/tablas/09_fcc.csv", row.names = FALSE)
print(CU[, -1])

# ---- Robustez: la ganancia de la tasa con otros modelos de ruido -------------------------------------
# Misma funcion de transferencia (b, s, r, phi elegidos), distintos ruidos. Si la ganancia se mantiene, el
# resultado no depende de la especificacion del ruido.
gt <- TF$tasa$g; xr_t <- cbind(regresores(ENT$tasa, gt$b, gt$s, if (gt$r == 1) gt$phi else 0), TCm)[INI:length(Y), , drop = FALSE]
W_t <- ts(as.numeric(Y)[INI:length(Y)], start = time(Y)[INI], frequency = S)
ruidos <- list(c(1, 2, 1, 1), c(0, 1, 1, 1), c(0, 1, 0, 1), c(1, 1, 1, 1), c(0, 2, 1, 1), c(1, 0, 1, 1))
RB <- do.call(rbind, lapply(ruidos, function(o) { m <- ajustar(W_t, o[1], o[2], o[3], o[4], xr_t, TRUE); if (is.null(m)) return(NULL)
  dl <- m$coef[sprintf("D%d", gt$b:(gt$b + gt$s))]; den <- 1 - (if (gt$r == 1) gt$phi else 0)
  ee <- if (gt$r == 1) ee_ganancia(function(ph) ajustar(W_t, o[1], o[2], o[3], o[4], cbind(regresores(ENT$tasa, gt$b, gt$s, ph), TCm)[INI:length(Y), , drop = FALSE], TRUE),
                                   gt$phi, names(dl))["ee"] else NA
  if (is.na(ee)) ee <- sqrt(sum(m$var.coef[names(dl), names(dl)])) / den
  data.frame(ruido = etiqueta_modelo(o[1], o[2], o[3], o[4]), gan = sum(dl) / den, ee = unname(ee),
             SBC = BIC(m), LB24 = diagnostico(m, sum(o))$LB24) }))
RB$dSBC <- RB$SBC - min(RB$SBC)
tabla_tex(data.frame(RB$ruido, sprintf("%s (%s)", num(RB$gan, 2), num(RB$ee, 2)), num(RB$gan * pp["tasa"], 1), num(RB$LB24, 3), num(RB$dSBC, 1)),
          "09_robustez_ruido.tex", rojo = celdas(1, 1:5), "lcccc", c("Ruido", "Ganancia $g$ (ee)", "\\% inicios por +1 pp", "LB24 $p$", "$\\Delta$SBC"))
guardar_valor("TFRobMin", num(min(RB$gan), 2)); guardar_valor("TFRobMax", num(max(RB$gan), 2))
guardar_valor("TFRobPPMin", num(max(RB$gan) * pp["tasa"], 1)); guardar_valor("TFRobPPMax", num(min(RB$gan) * pp["tasa"], 1))
guardar_valor("TFTasaPhi", num(gt$phi, 2)); print(RB)

# ---- Salidas EViews ---------------------------------------------------------------------------------
for (e in names(TF)) { t <- TF[[e]]
  guardar_salida(salida_arma(t$m, t$W, "D(LHOUST_SA)", sprintf("Funcion de transferencia: %s -> inicios", if (e == "tasa") "D(LMORT)" else "D(LGS10)")),
                 sprintf("09_modelo_%s.png", e))
  k_ru <- sum(grepl("^(ar|ma|sar|sma)[0-9]", names(t$m$coef)))             # parametros del ruido, para los gl del Q
  guardar_salida(salida_correlograma(residuals(t$m), sprintf("Residuos de la funcion de transferencia (%s)", e), 36, k_ru, rojo_k = c(12, 24, 36)),
                 sprintf("09_corr_residuos_%s.png", e)) }
# Correlograma del input preblanqueado: tiene que ser ruido blanco (paso 1 del apunte)
for (e in names(PB)) guardar_salida(salida_correlograma(PB[[e]]$alfa, sprintf("Correlograma de alfa: %s preblanqueada", if (e == "tasa") "tasa" else "bono"), 24,
                                                        sum(grepl("^(ar|ma)[0-9]", names(PB[[e]]$mX$coef))), rojo_k = c(12, 24)), sprintf("09_corr_alfa_%s.png", e))

# ---- Figura: FCC de los preblanqueados y pesos v_j ---------------------------------------------------
FC <- do.call(rbind, lapply(names(PB), function(e) transform(PB[[e]]$C, input = unname(c(tasa = "Tasa hipotecaria", bono = "Bono 10 anios")[e]))))
FC$input <- factor(FC$input, levels = c("Tasa hipotecaria", "Bono 10 anios"))
g1 <- ggplot(FC, aes(k, r)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -banda, ymax = banda, fill = gris, alpha = 0.45) +
  geom_vline(xintercept = 0, colour = tinta2, linewidth = 0.3) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.1) + scale_colour_manual(values = c(`TRUE` = rojo, `FALSE` = "#9AA5B1"), guide = "none") +
  facet_wrap(~input, ncol = 1) + scale_x_continuous(breaks = seq(-24, 24, 6)) +
  labs(title = "FCC del input y el output preblanqueados", subtitle = "k > 0: el input anticipa a los inicios; k < 0: los inicios anticipan al input (retroalimentacion)",
       x = "k (meses)", y = NULL) + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
VJ <- do.call(rbind, lapply(names(TF), function(e) data.frame(input = unname(c(tasa = "Tasa hipotecaria", bono = "Bono 10 anios")[e]), j = 0:K,
                                                              emp = PB[[e]]$C$v[PB[[e]]$C$k >= 0], mod = TF[[e]]$vj)))
VJ <- VJ[VJ$j <= 12, ]; VJ$input <- factor(VJ$input, levels = c("Tasa hipotecaria", "Bono 10 anios"))
g2 <- ggplot(VJ, aes(j)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_col(aes(y = emp), fill = "#9AA5B1", width = 0.6) + geom_line(aes(y = mod), colour = rojo, linewidth = 0.8) + geom_point(aes(y = mod), colour = rojo, size = 1.4) +
  facet_wrap(~input, ncol = 1) + scale_x_continuous(breaks = 0:12) +
  labs(title = "Pesos de respuesta al impulso v_j", subtitle = "barras: desde la FCC (v_j = sd(beta)/sd(alfa) r_j); rojo: los del modelo estimado",
       x = "j (meses)", y = "% de inicios por 1% del input") + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
guardar_fig((g1 | g2) + plot_layout(widths = c(1.2, 1)), "09_transferencia.png", 12.5, 6.5)
escribir_valores("valores_09.tex")
