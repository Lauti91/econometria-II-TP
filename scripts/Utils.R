# =============================================================================
# utils.R  --  Estilo de graficos y carga de datos compartidos (se usa con source())
# =============================================================================

for (pk in c("ggplot2", "patchwork")) {
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ---- Estilo ---------------------------------------------------------------------
# Un color por variable (PBI azul, electricidad naranja), el mismo en todos los graficos.
col_pbi <- "#2A62B8"; col_elec <- "#D9701F"
tinta <- "#1F2933"; tinta2 <- "#52606D"; gris <- "#C5CCD3"; fondo <- "#FCFCFB"

tema <- theme_minimal(base_size = 11) +
  theme(
    plot.background    = element_rect(fill = fondo, colour = NA),
    panel.background   = element_rect(fill = fondo, colour = NA),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E4E7EB", linewidth = 0.3),
    axis.text  = element_text(colour = tinta2, size = 9),
    axis.title = element_text(colour = tinta2, size = 9),
    plot.title = element_text(colour = tinta, face = "bold", size = 11, margin = margin(b = 2)),
    plot.subtitle = element_text(colour = tinta2, size = 9, margin = margin(b = 6)),
    plot.margin = margin(8, 12, 6, 8)
  )

tema_global <- theme(
  plot.background = element_rect(fill = fondo, colour = NA),
  plot.title   = element_text(colour = tinta, face = "bold", size = 14),
  plot.caption = element_text(colour = tinta2, size = 8, hjust = 0))

# ---- Datos ----------------------------------------------------------------------
# Muestra principal 2000Q1-2023Q4 (96 trimestres): crecimiento interanual (%) del PBI oficial
# y de la electricidad, construidos por 01_construir_base.R.
cargar_yoy <- function() {
  Q <- read.csv("output/trimestral_li_keqiang.csv", stringsAsFactors = FALSE)
  S <- Q[Q$trimestre >= "2000Q1" & Q$trimestre <= "2023Q4", ]
  stopifnot(nrow(S) == 96, !anyNA(S$pbi_crec_interanual_oficial_pct), !anyNA(S$elec_crec_interanual_pct))
  list(
    "PBI oficial, crec. interanual (%)"  = ts(S$pbi_crec_interanual_oficial_pct,  start = c(2000, 1), frequency = 4),
    "Electricidad, crec. interanual (%)" = ts(S$elec_crec_interanual_pct,         start = c(2000, 1), frequency = 4))
}

etiqueta_trim <- function(x) {
  t <- as.numeric(time(x)); sprintf("%dQ%d", floor(t + 1e-6), round((t - floor(t + 1e-6)) * 4) + 1)
}