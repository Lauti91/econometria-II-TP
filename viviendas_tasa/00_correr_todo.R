# =============================================================================
# 00_correr_todo.R -- Corre el analisis completo del par en orden. Desde viviendas_tasa/:
#   Rscript 00_correr_todo.R
# Deja figuras en salidas/figuras, salidas tipo EViews en salidas/eviews y tablas y valores para el
# informe en salidas/tablas. Despues se compila informe/informe_viviendas.tex.
# =============================================================================
for (s in c("01_datos.R", "02_exploracion.R", "03_estacionariedad.R", "04_inicios.R", "05_tasa.R", "06_pronostico.R", "07_relacion.R", "08_intervencion.R", "09_transferencia.R", "10_var.R")) {
  cat("\n==========", s, "==========\n"); source(s, echo = FALSE, local = new.env())
}
