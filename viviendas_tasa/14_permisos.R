# =============================================================================
# 14_permisos.R -- Permisos de construccion como output alternativo: tasa hipotecaria -> permisos
#   El permiso es la decision de construir; el inicio es su ejecucion (hay que conseguir financiamiento,
#   contratistas y buen tiempo). Hipotesis: si la tasa actua sobre la DECISION, los permisos tienen que
#   responder antes que los inicios (retardo menor) y con una ganancia de largo plazo parecida.
#   Datos (FRED, se bajan una sola vez y se guardan en datos/crudos/ con su fecha):
#     PERMIT     permisos de viviendas privadas nuevas, desestacionalizados, tasa anual, miles
#     PERMITNSA  idem sin desestacionalizar, miles (solo para el control de consistencia)
#   Pasos:
#     1. Box-Jenkins rapido de 100 ln(permisos SA): ADF (parte 1; valores p de MacKinnon 1996, Journal of
#        Applied Econometrics, 11(6), 601-618), identificacion por rejilla y SBC, diagnostico. El COVID entra
#        como cambio transitorio (apunte parte 2, pp. 1-7) con su propio delta por perfil de verosimilitud.
#     2. Funcion de transferencia tasa -> permisos con el MISMO procedimiento de 09 (apunte parte 2, pp. 8-21):
#        preblanqueo, FCC, (b, s, r) por SBC, phi por perfil, re-identificacion del ruido, diagnostico.
#        Las funciones son las de 09_transferencia.R, con los permisos como output.
#     3. Comparacion con los inicios: retardo b, retardo medio de la respuesta (suma j v_j / suma v_j),
#        ganancia de largo plazo y % por +1 pp de tasa. Como control, los permisos con la forma de los inicios
#        (b = 2, s = 0, r = 1) y los inicios con la forma de los permisos.
#   No se estima la funcion permisos -> inicios (fuera del alcance).
#   Salidas: tablas 14_univariado.tex, 14_seleccion.tex, 14_transferencia.tex, 14_fcc.csv; eviews
#   14_modelo_permisos.png, 14_corr_residuos.png; figura 14_permisos.png; valores_14.tex.
# Ejecutar desde viviendas_tasa/:  Rscript 14_permisos.R   (despues de 08 y 09)
# =============================================================================
source("R/funciones.R")
# ---- 0. Datos ---------------------------------------------------------------------------------------------
ids <- c("PERMIT", "PERMITNSA"); DESCARGAR <- !all(file.exists(sprintf("datos/crudos/%s.csv", ids)))   # solo si faltan
if (DESCARGAR) {
  for (id in ids) { d <- read.csv(sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id), na.strings = ".")
    names(d) <- c("fecha", "valor"); write.csv(d[!is.na(d$valor), ], sprintf("datos/crudos/%s.csv", id), row.names = FALSE) }
  writeLines(sprintf("Descarga FRED (PERMIT, PERMITNSA): %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "datos/crudos/fecha_descarga_permisos.txt")
}
L <- cargar_base(); I8 <- readRDS("salidas/tablas/08_intervencion.rds")
leer <- function(id) { d <- read.csv(sprintf("datos/crudos/%s.csv", id)); m <- substr(d$fecha, 1, 7); v <- setNames(d$valor, m)[etiqueta_mes(L$viv_sa)]
  stopifnot(!anyNA(v)); ts(as.numeric(v), start = start(L$viv_sa), frequency = S) }
perm <- leer("PERMIT") / 12; perm_nsa <- leer("PERMITNSA")                 # miles por mes
k <- 12:length(perm); dif <- 100 * (sapply(k, function(i) sum(perm_nsa[(i - 11):i])) / sapply(k, function(i) sum(perm[(i - 11):i])) - 1)
cat(sprintf("Permisos %s a %s. Control NSA vs SA (sumas de 12 meses): dif. media %.2f%%, max |dif| %.2f%%\n", etiqueta_mes(perm)[1], etiqueta_mes(perm)[length(perm)], mean(dif), max(abs(dif))))
lp <- 100 * log(perm)                                                      # 100 ln(permisos SA)
Y <- diff(lp); X <- diff(L$tasa); YI <- diff(L$viv_sa)                     # outputs e input, en % mensual
K <- 24; n <- length(Y); banda <- 2 / sqrt(n)
ult <- L$tasa_pp[length(L$tasa_pp)]; pp <- 100 * log((ult + 1) / ult)       # +1 pp sobre la tasa de 2026-08
tc <- function(delta) diff(ts(as.numeric(stats::filter(indicadora(lp, 2020, 4), delta, method = "recursive")), start = start(lp), frequency = S))

# ---- 1. Box-Jenkins de los permisos -----------------------------------------------------------------------
A0 <- adf(lp, 12, "c", v = "LPERMIT"); A1 <- adf(Y, 12, "c", v = "D(LPERMIT)")
TCm <- cbind(COVID = as.numeric(tc(I8$delta)))
R1 <- rejilla(Y, TCm, TRUE, 3, 2, 1, 1); R1a <- R1[R1$admisible, ]; ru_u <- unlist((if (nrow(R1a)) R1a else R1)[1, c("p", "q", "P", "Q")])
grilla_d <- seq(0.05, 0.95, by = 0.05)                                     # delta del TC por perfil, con el ruido elegido
perfil_d <- sapply(grilla_d, function(d) { m <- ajustar(Y, ru_u[1], ru_u[2], ru_u[3], ru_u[4], cbind(COVID = as.numeric(tc(d))), TRUE); if (is.null(m)) NA else as.numeric(logLik(m)) })
delta_p <- grilla_d[which.max(perfil_d)]; TCm <- cbind(COVID = as.numeric(tc(delta_p))); TC <- tc(delta_p)
RY <- rejilla(Y, TCm, TRUE, 3, 2, 1, 1); n_adm <- sum(RY$admisible); n_tot <- nrow(RY); RYa <- RY[RY$admisible, ]
if (!nrow(RYa)) { warning("ningun ruido admisible para los permisos: se toma el de menor SBC"); RYa <- RY }
ruido0 <- unlist(RYa[1, c("p", "q", "P", "Q")]); mY <- ajustar(Y, ruido0[1], ruido0[2], ruido0[3], ruido0[4], TCm, TRUE); DU <- diagnostico(mY, sum(ruido0))
Ylimpia <- Y - mY$coef["COVID"] * TC
cat(sprintf("ADF nivel: tau = %.2f (p = %.3f); diferencia: tau = %.2f (p = %.4f). Ruido: %s (%d de %d admisibles). delta COVID = %.2f, omega = %.1f\n",
            A0$tau, A0$p, A1$tau, A1$p, RYa$modelo[1], n_adm, n_tot, delta_p, mY$coef["COVID"]))
tabla_tex(data.frame(c("ADF de $100\\ln(\\text{permisos})$, con constante, 12 rezagos", "ADF de la primera diferencia", "Modelo de la diferencia (SBC entre los admisibles)",
                       "COVID: cambio transitorio, $\\delta$ por perfil", "Ljung--Box de los residuos (12, 24, 36)", "$\\hat\\sigma$ (\\% mensual)"),
                     c(sprintf("$\\tau=%s$ ($p=%s$)", num(A0$tau, 2), pval(A0$p)), sprintf("$\\tau=%s$ ($p=%s$)", num(A1$tau, 2), pval(A1$p)), RYa$modelo[1],
                       sprintf("$\\delta=%s$, $\\hat\\omega=%s$ (%s\\,\\%% en abril de 2020)", num(delta_p, 2), num(mY$coef["COVID"], 1), num(100 * (exp(mY$coef["COVID"] / 100) - 1), 1)),
                       sprintf("$p=%s$; %s; %s", num(DU$LB12, 2), num(DU$LB24, 2), num(DU$LB36, 2)), num(DU$sigma, 2))),
          "14_univariado.tex", rojo = celdas(3, 2), "ll", c("Paso", "Resultado"))
guardar_valor("PermADFNivel", num(A0$tau, 2)); guardar_valor("PermADFNivelP", pval(A0$p)); guardar_valor("PermADFDif", num(A1$tau, 2)); guardar_valor("PermADFDifP", pval(A1$p))
guardar_valor("PermRuidoUni", RYa$modelo[1]); guardar_valor("PermDelta", num(delta_p, 2)); guardar_valor("PermOmega", num(mY$coef["COVID"], 1))
guardar_valor("PermOmegaPct", num(100 * (exp(mY$coef["COVID"] / 100) - 1), 1)); guardar_valor("PermSigmaUni", num(DU$sigma, 2)); guardar_valor("PermLBUni", num(DU$LB24, 2))
guardar_valor("PermSigmaIni", num(sd(YI), 2)); guardar_valor("PermSD", num(sd(Y), 2)); guardar_valor("PermCorrIni", num(cor(Y, YI), 2))

# ---- 2. Funcion de transferencia tasa -> permisos (funciones de 09) ------------------------------------------
preblanqueo <- function(X, Yl) {
  R <- rbind(transform(rejilla(X, NULL, TRUE, 3, 2, 0, 0), media = TRUE), transform(rejilla(X, NULL, FALSE, 3, 2, 0, 0), media = FALSE))
  R <- R[R$admisible, ]; R <- R[order(R$SBC), ]; g <- R[1, ]
  mX <- ajustar(X, g$p, g$q, 0, 0, NULL, g$media)
  filtro <- function(z) { z <- z - mean(z); ar <- mX$model$phi; ma <- mX$model$theta
    if (length(ar)) z <- stats::filter(z, c(1, -ar), sides = 1)
    if (length(ma)) z <- stats::filter(ifelse(is.na(z), 0, z), -ma, method = "recursive")
    ts(as.numeric(z), start = start(X), frequency = S) }
  a <- filtro(X); b <- filtro(Yl)
  ok <- !is.na(a) & !is.na(b) & seq_along(a) > 12
  a <- ts(a[ok], start = time(X)[which(ok)[1]], frequency = S); b <- b[ok]
  r <- ccf(b, as.numeric(a), lag.max = K, plot = FALSE)                    # r(k) = corr(beta_{t+k}, alfa_t): k > 0, el input anticipa
  C <- data.frame(k = as.numeric(r$lag), r = as.numeric(r$acf)); C$sig <- abs(C$r) > banda
  C$v <- ifelse(C$k >= 0, sd(b) / sd(a) * C$r, NA)
  Qk <- function(lados, M = 12) { d <- C[if (lados == "pos") C$k >= 0 & C$k <= M else C$k < 0 & C$k >= -M, ]
    Q <- length(a) * (length(a) + 2) * sum(d$r^2 / (length(a) - abs(d$k))); c(Q = Q, gl = nrow(d), p = pchisq(Q, nrow(d), lower.tail = FALSE)) }
  list(modelo = g$modelo, mX = mX, alfa = a, C = C, Qpos = Qk("pos"), Qneg = Qk("neg"))
}
PB <- preblanqueo(X, Ylimpia)
# Control de direccion (ccf de R): el mismo calculo a mano para k = +1, corr(alfa_t, beta_{t+1})
{ a <- as.numeric(PB$alfa); f <- function(z) { z <- z - mean(z); th <- PB$mX$model$theta; if (length(th)) z <- stats::filter(z, -th, method = "recursive"); as.numeric(z) }
  b <- f(Ylimpia)[-(1:12)]; stopifnot(length(PB$mX$model$phi) > 0 || abs(cor(a[-length(a)], b[-1]) - PB$C$r[PB$C$k == 1]) < 0.01) }
cat("Filtro del input:", PB$modelo, "| FCC significativa en k =", paste(PB$C$k[PB$C$sig], collapse = " "), "\n")
INI <- 7
regresores <- function(Yv, b, s, phi = 0) { xf <- as.numeric(if (phi != 0) stats::filter(X, phi, method = "recursive") else X); nY <- length(Yv)
  M <- sapply(b:(b + s), function(j) c(rep(0, j), xf[seq_len(nY - j)])); matrix(M, nY, dimnames = list(NULL, sprintf("D%d", b:(b + s)))) }
serie <- function(Yv) ts(as.numeric(Yv)[INI:length(Yv)], start = time(Yv)[INI], frequency = S)
estimar <- function(Yv, TCv, b, s, r, phi, ruido) ajustar(serie(Yv), ruido[1], ruido[2], ruido[3], ruido[4], cbind(regresores(Yv, b, s, if (r == 1) phi else 0), TCv)[INI:length(Yv), , drop = FALSE], TRUE)
ee_ganancia <- function(ajuste_phi, phi, nd, h = 0.02) {                    # ee de g = sum(delta)/(1 - phi) con la incertidumbre de phi (ver 09)
  M <- lapply(c(phi - h, phi, phi + h), ajuste_phi)
  if (any(sapply(M, is.null))) return(c(ee = NA, ee_cond = NA, ee_phi = NA))
  ll <- sapply(M, function(m) as.numeric(logLik(m))); d <- sapply(M, function(m) sum(m$coef[nd]))
  curv <- (ll[3] - 2 * ll[2] + ll[1]) / h^2; v_cond <- sum(M[[2]]$var.coef[nd, nd])
  if (curv >= 0) return(c(ee = NA, ee_cond = sqrt(v_cond) / (1 - phi), ee_phi = NA))
  v_phi <- -1 / curv; pend <- (d[3] - d[1]) / (2 * h)
  V <- matrix(c(v_cond + pend^2 * v_phi, pend * v_phi, pend * v_phi, v_phi), 2); gr <- c(1 / (1 - phi), d[2] / (1 - phi)^2)
  c(ee = sqrt(drop(t(gr) %*% V %*% gr)), ee_cond = sqrt(v_cond) / (1 - phi), ee_phi = sqrt(v_phi))
}
grilla_phi <- seq(-0.8, 0.8, by = 0.2)
seleccion <- function(Yv, TCv, ruido) { out <- list()
  for (b in 0:3) for (s in 0:3) for (r in 0:1) { phis <- if (r == 0) 0 else grilla_phi
    ms <- lapply(phis, function(ph) estimar(Yv, TCv, b, s, r, ph, ruido)); ll <- sapply(ms, function(m) if (is.null(m)) -Inf else as.numeric(logLik(m)))
    i <- which.max(ll); m <- ms[[i]]; if (is.null(m)) next
    out[[length(out) + 1]] <- data.frame(b = b, s = s, r = r, phi = phis[i], SBC = -2 * ll[i] + (length(m$coef) + r) * log(length(residuals(m))), logLik = ll[i], row.names = NULL) }
  G <- do.call(rbind, out); G$dSBC <- G$SBC - min(G$SBC); G[order(G$SBC), ] }
# Estimacion completa de una forma (b, s, r): phi fino por perfil, ruido re-identificado (o fijo), ganancia y diagnostico
forma <- function(Yv, TCv, b, s, r, phi0, ruido_ini, alfa = NULL, ruido_fijo = NULL) {
  phi <- 0
  if (r == 1) { fina <- seq(max(-0.95, phi0 - 0.15), min(0.95, phi0 + 0.15), by = 0.05)
    ll <- sapply(fina, function(ph) { m <- estimar(Yv, TCv, b, s, 1, ph, ruido_ini); if (is.null(m)) -Inf else as.numeric(logLik(m)) }); phi <- fina[which.max(ll)] }
  xr <- cbind(regresores(Yv, b, s, phi), TCv)[INI:length(Yv), , drop = FALSE]; W <- serie(Yv)
  if (is.null(ruido_fijo)) { R <- rejilla(W, xr, TRUE, 3, 2, 1, 1); Ra <- R[R$admisible, ]; if (!nrow(Ra)) Ra <- R; ru <- unlist(Ra[1, c("p", "q", "P", "Q")]); et <- Ra$modelo[1] } else {
    ru <- ruido_fijo; et <- etiqueta_modelo(ru[1], ru[2], ru[3], ru[4]) }
  m <- ajustar(W, ru[1], ru[2], ru[3], ru[4], xr, TRUE)
  nd <- sprintf("D%d", b:(b + s)); dl <- m$coef[nd]; gan <- sum(dl) / (1 - phi); gan_ee <- sqrt(sum(m$var.coef[nd, nd])) / (1 - phi); gan_cond <- gan_ee; phi_ee <- NA
  if (r == 1) { E <- ee_ganancia(function(ph) ajustar(W, ru[1], ru[2], ru[3], ru[4], cbind(regresores(Yv, b, s, ph), TCv)[INI:length(Yv), , drop = FALSE], TRUE), phi, nd)
    if (!is.na(E["ee"])) { gan_ee <- unname(E["ee"]); phi_ee <- unname(E["ee_phi"]) } }
  vj <- sapply(0:240, function(j) { if (j < b) return(0); sum(sapply(b:min(j, b + s), function(i) dl[sprintf("D%d", i)] * (if (r == 1) phi^(j - i) else as.numeric(j == i)))) })
  p_cc <- c(pos = NA, neg = NA)
  if (!is.null(alfa)) { e <- ts(as.numeric(residuals(m)), start = start(W), frequency = S)
    a <- ts(c(rep(NA, length(Yv) - length(alfa)), alfa), start = start(Yv), frequency = S); Z <- na.omit(ts.intersect(e, a))
    ce <- ccf(as.numeric(Z[, "e"]), as.numeric(Z[, "a"]), lag.max = K, plot = FALSE); CE <- data.frame(k = as.numeric(ce$lag), r = as.numeric(ce$acf))
    Qe <- function(pos) { d <- CE[if (pos) CE$k >= 0 & CE$k <= 12 else CE$k < 0 & CE$k >= -12, ]; Q <- nrow(Z) * (nrow(Z) + 2) * sum(d$r^2 / (nrow(Z) - abs(d$k)))
      pchisq(Q, if (pos) nrow(d) - (s + 1 + r) else nrow(d), lower.tail = FALSE) }
    p_cc <- c(pos = Qe(TRUE), neg = Qe(FALSE)) }
  acum <- cumsum(vj) / sum(vj)
  list(b = b, s = s, r = r, phi = phi, phi_ee = phi_ee, ruido = et, ru = ru, m = m, W = W, D = diagnostico(m, sum(ru)), vj = unname(vj), gan = gan, gan_ee = gan_ee, gan_cond = gan_cond,
       medio = sum((0:240) * vj) / sum(vj), mitad = which(acum >= 0.5)[1] - 1, p_cc = p_cc, n_par = length(m$coef) + r)
}
G <- seleccion(Y, TCm, ruido0); print(head(G, 8), row.names = FALSE); g <- G[1, ]
TFp <- forma(Y, TCm, g$b, g$s, g$r, g$phi, ruido0, PB$alfa)
m_nulo <- ajustar(serie(Y), TFp$ru[1], TFp$ru[2], TFp$ru[3], TFp$ru[4], TCm[INI:length(Y), , drop = FALSE], TRUE)   # mismo ruido, sin input
dSBC_nulo <- (-2 * as.numeric(logLik(TFp$m)) + TFp$n_par * log(length(residuals(TFp$m)))) - BIC(m_nulo)
G6 <- head(G, 6)
tabla_tex(data.frame(G6$b, G6$s, G6$r, ifelse(G6$r == 1, num(G6$phi, 1), "--"), num(G6$SBC, 1), num(G6$dSBC, 1)), "14_seleccion.tex", rojo = celdas(1, 1:6), "cccccc",
          c("Retardo $b$", "$s$", "$r$", "$\\phi$ (grilla 0,2)", "SBC (ruido inicial)", "$\\Delta$SBC"))

# ---- 3. Comparacion con los inicios -------------------------------------------------------------------------
TCi <- cbind(COVID = as.numeric(diff(ts(as.numeric(stats::filter(indicadora(L$viv_sa, 2020, 4), I8$delta, method = "recursive")), start = start(L$viv_sa), frequency = S))))
RU_I <- c(1, 2, 1, 1)                                                      # ruido de los inicios elegido en 09
TFi <- forma(YI, TCi, 2, 0, 1, 0.75, RU_I, ruido_fijo = RU_I)              # inicios, forma de 09 (control: reproduce su ganancia)
TFp_i <- forma(Y, TCm, 2, 0, 1, 0.75, ruido0, PB$alfa)                     # permisos con la forma de los inicios
TFi_p <- forma(YI, TCi, g$b, g$s, g$r, g$phi, RU_I, ruido_fijo = RU_I)     # inicios con la forma de los permisos
sbc <- function(t) -2 * as.numeric(logLik(t$m)) + t$n_par * log(length(residuals(t$m)))
fila <- function(etq, t) data.frame(etq, sprintf("$b=%d$, $s=%d$, $r=%d$", t$b, t$s, t$r), if (t$r == 1) num(t$phi, 2) else "--", t$ruido,
                                    sprintf("%s (%s)", num(t$gan, 2), num(t$gan_ee, 2)), num(t$gan * pp, 1), num(t$medio, 1), t$mitad, num(t$D$LB24, 2), num(sbc(t), 1))
CMP <- rbind(fila("Permisos, forma elegida por SBC", TFp), fila("Permisos, forma de los inicios", TFp_i), fila("Inicios, forma de 09", TFi), fila("Inicios, forma de los permisos", TFi_p))
tabla_tex(CMP, "14_transferencia.tex", rojo = rbind(celdas(1, 2:8), celdas(3, 2:8)), "lccccccccc",
          c("Output y forma", "Funci\\'on de transferencia", "$\\phi$", "Ruido", "Ganancia $g$ (ee)", "\\% por +1 pp", "Retardo medio (meses)", "Meses hasta la mitad del efecto", "LB24 $p$", "SBC"))
for (x in list(list("Perm", TFp), list("PermFI", TFp_i), list("Ini", TFi), list("IniFP", TFi_p))) { s <- x[[1]]; t <- x[[2]]
  guardar_valor(paste0("TFP", s, "Gan"), num(t$gan, 2)); guardar_valor(paste0("TFP", s, "GanEE"), num(t$gan_ee, 2)); guardar_valor(paste0("TFP", s, "PP"), num(t$gan * pp, 1))
  guardar_valor(paste0("TFP", s, "Medio"), num(t$medio, 1)); guardar_valor(paste0("TFP", s, "Mitad"), t$mitad); guardar_valor(paste0("TFP", s, "Phi"), num(t$phi, 2))
  guardar_valor(paste0("TFP", s, "Esp"), sprintf("$b=%d$, $s=%d$, $r=%d$", t$b, t$s, t$r)); guardar_valor(paste0("TFP", s, "B"), t$b) }
guardar_valor("TFPPermRuido", TFp$ruido); guardar_valor("TFPPermLB", num(TFp$D$LB24, 2)); guardar_valor("TFPPermFiltro", PB$modelo)
guardar_valor("TFPPermPos", pval(PB$Qpos["p"])); guardar_valor("TFPPermNeg", pval(PB$Qneg["p"])); guardar_valor("TFPPermCCpos", pval(TFp$p_cc["pos"])); guardar_valor("TFPPermCCneg", pval(TFp$p_cc["neg"]))
guardar_valor("TFPPermDSBC", num(dSBC_nulo, 1)); guardar_valor("TFPPermDSBCFormaIni", num(sbc(TFp_i) - sbc(TFp), 1)); guardar_valor("TFPIniDSBCFormaPerm", num(sbc(TFi_p) - sbc(TFi), 1))
guardar_valor("TFPPermSig", { c_ <- PB$C[PB$C$sig, ]; if (nrow(c_)) paste(sprintf("$k=%d$ (%s)", c_$k, num(c_$r, 2)), collapse = ", ") else "ninguno" })
guardar_valor("TFPPermPhiEE", if (is.na(TFp$phi_ee)) "--" else num(TFp$phi_ee, 2))
write.csv(PB$C, "salidas/tablas/14_fcc.csv", row.names = FALSE)
print(CMP[, -1], row.names = FALSE)
cat(sprintf("Permisos: FCC k>=0 p = %.3f, k<0 p = %.3f | residuos vs alfa: k>=0 p = %.3f, k<0 p = %.3f | SBC vs sin input: %.1f\n", PB$Qpos["p"], PB$Qneg["p"], TFp$p_cc["pos"], TFp$p_cc["neg"], dSBC_nulo))

# ---- Salidas EViews -----------------------------------------------------------------------------------------
guardar_salida(salida_arma(TFp$m, TFp$W, "D(LPERMIT)", "Funcion de transferencia: D(LMORT) -> permisos"), "14_modelo_permisos.png")
guardar_salida(salida_correlograma(residuals(TFp$m), "Residuos de la funcion de transferencia (permisos)", 36, sum(TFp$ru), rojo_k = c(12, 24, 36)), "14_corr_residuos.png")

# ---- Figura ---------------------------------------------------------------------------------------------------
# Izquierda: FCC del preblanqueo, permisos e inicios (la de inicios, de 09). Derecha: respuesta acumulada a +1 % de la tasa.
FI9 <- read.csv("salidas/tablas/09_fcc.csv"); FI9 <- FI9[FI9$input == "tasa", c("k", "r", "sig")]
FC <- rbind(cbind(PB$C[, c("k", "r", "sig")], out = "Permisos"), cbind(FI9, out = "Inicios")); FC$out <- factor(FC$out, levels = c("Permisos", "Inicios")); FC <- FC[abs(FC$k) <= 12, ]
g1 <- ggplot(FC, aes(k, r)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -banda, ymax = banda, fill = gris, alpha = 0.45) +
  geom_vline(xintercept = 0, colour = tinta2, linewidth = 0.3) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_segment(aes(xend = k, yend = 0, colour = sig), linewidth = 1.3) + scale_colour_manual(values = c(`TRUE` = rojo, `FALSE` = "#9AA5B1"), guide = "none") +
  facet_wrap(~out, ncol = 1) + scale_x_continuous(breaks = seq(-12, 12, 3)) +
  labs(title = "FCC de la tasa y el output preblanqueados", subtitle = "k > 0: la tasa anticipa al output; en rojo, fuera de la banda", x = "k (meses)", y = NULL) + tema + theme(strip.text = element_text(colour = tinta, face = "bold"))
AC <- rbind(data.frame(j = 0:24, y = cumsum(TFp$vj)[1:25], out = "Permisos"), data.frame(j = 0:24, y = cumsum(TFi$vj)[1:25], out = "Inicios"))
g2 <- ggplot(AC, aes(j, y, colour = out)) + geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
  geom_hline(yintercept = c(TFp$gan, TFi$gan), colour = c(col_tasa, col_viv), linetype = "22", linewidth = 0.35) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) + scale_colour_manual(values = c(Permisos = col_tasa, Inicios = col_viv), name = NULL) + scale_x_continuous(breaks = seq(0, 24, 3)) +
  labs(title = "Respuesta acumulada a +1 % de la tasa", subtitle = sprintf("punteado: ganancia de largo plazo (permisos %.2f, inicios %.2f)", TFp$gan, TFi$gan), x = "Meses desde el cambio de la tasa", y = "% del output") +
  tema + theme(legend.position = "bottom")
guardar_fig((g1 | g2) + plot_layout(widths = c(1, 1.1)), "14_permisos.png", 12, 5.6)
escribir_valores("valores_14.tex")
