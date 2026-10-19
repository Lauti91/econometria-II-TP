# econometria-II-TP

Trabajo práctico final de Econometría II (cátedra González, FCE-UBA).
Candidata 5: **Índice Li Keqiang** — PBI oficial de China vs. indicadores físicos de actividad
(electricidad y carga ferroviaria).

## Estructura

```
scripts/   código en R (se numera en orden de ejecución)
series/    series originales tal cual se descargan de cada fuente (sin modificar)
output/    bases construidas, tablas y gráficos
```

## Cómo reproducir

Desde la raíz del repositorio:

```
Rscript scripts/01_construir_base.R
```

Necesita R base y, opcionalmente, `openxlsx` (sin él escribe solo CSV).
El script vuelve a descargar las series a `series/` y reconstruye `output/`.

## Series (descargadas el 2026-10-19)

| Archivo en `series/` | Contenido | Unidad / frecuencia | Cobertura | Fuente |
|---|---|---|---|---|
| `gdp_nominal.csv` | PBI nominal | 1e8 yuan, trimestral | 1992Q1–2026Q2 | chinadata.live (espejo NBS), `china-gdp-quarterly` |
| `gdp_indice.csv` | Crecimiento interanual oficial (año previo = 100) | índice, trimestral | 1993Q1–2026Q2 | chinadata.live, `china-gdp-index` |
| `gdp_real.csv` | PBI real (base 2011) | 1e8 yuan, trimestral | 2011Q1–2026Q2 | chinadata.live, `china-gdp-real` |
| `elec_oecd.csv` | Producción de electricidad | GWh, mensual | 1999-01–2023-11 | OCDE vía FRED, `CHNPRENEL01MLM` |
| `elec_cd.csv` | Generación de electricidad | 1e8 kWh, mensual | 2007-02–2026-08 | chinadata.live (NBS), `china-electricity-generation` |
| `rail_cd.csv` | Carga ferroviaria | 1e4 toneladas, mensual | 2005-01–2026-07 | chinadata.live (NBS), `china-rail-freight` |

## Decisiones de construcción (ver hoja `Notas` de la base)

- **Electricidad:** se usa la serie OCDE hasta 2023-11 y chinadata.live desde 2023-12 (×100 para pasar a GWh).
  En el solape (marzo–diciembre) ambas coinciden; el script frena con error si dejan de coincidir.
- **Febrero duplicado:** en la serie OCDE, 2012–13 y 2016–22, el "febrero" trae el acumulado enero–febrero
  (= 2 × enero). Se reparte en partes iguales.
- **Huecos desde 2024:** chinadata.live no trae enero ni (casi siempre) febrero, ni la carga ferroviaria de
  enero–febrero de 2024–25 y enero de 2026. Los primeros trimestres de 2024–2026 quedan vacíos; no se imputa.
- **Índice Li Keqiang:** de 2 componentes (promedio simple del crecimiento interanual de electricidad y
  carga ferroviaria). No hay serie larga gratuita de préstamos bancarios, el tercer componente original.
- **PBI:** los niveles de chinadata.live difieren entre 0,4% y 2,6% de la serie OCDE/FRED
  (discontinuada en 2023) en los 127 trimestres en común; es consistente con revisiones.

## Muestras disponibles (crecimiento interanual, trimestral)

| Par | Observaciones contiguas | Período |
|---|---|---|
| PBI oficial × electricidad | 96 | 2000Q1–2023Q4 |
| PBI oficial × carga ferroviaria / índice 2 comp. | 72 | 2006Q1–2023Q4 |

Correlación PBI × electricidad (interanual): 0,76.
