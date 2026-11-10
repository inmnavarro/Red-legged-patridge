# Red-legged partridge: environmental favourability and spring abundance in Spain

Analysis code associated with the manuscript *From favourability to a spring abundance index: an assessment of red-legged partridge in mainland Spain* (under review).

The study combines environmental favourability modelling with spring transect observations from the **Observatorio Cinegético** (2021–2024). Quantile regressions are used to explore the relationship between favourability and a cell-level abundance index, assess sensitivity to zero-detection surveys and quantile choice, and calculate spatial block-bootstrap intervals for projected scenario sums.

## Repository contents

All scripts are stored in [`sintaxis/`](sintaxis/):

| File | Purpose |
| --- | --- |
| `01_MODELO_AMBIENTAL_20_VARIABLES.R` | Refit the final 20-predictor logistic model, calculate environmental favourability, check saved predictions, and reconstruct classification statistics. |
| `02_VERIFICAR_MODELO_20_VARIABLES.sps` | Independently check the saved logistic coefficients and favourability values in SPSS. |
| `03_HISTORIAL_SPSS_ORIGINAL_RECUPERADO_NO_EJECUTAR.sps` | Document the recovered original SPSS command history. **Archival reference only; do not execute without review.** |
| `04_ABUNDANCIA_Y_CUANTILES.R` | Reconstruct representative transects, retain valid zero detections, aggregate the abundance index by grid cell, fit quantile regressions, compare scenarios, create figures, and carry out exploratory spatial cross-validation. |
| `05_BOOTSTRAP_ESPACIAL_IC95.R` | Refit quantile models under 1,000 previously generated spatial block-bootstrap resamples and calculate percentile intervals. |

## Software and required inputs

The analyses were developed in **R 4.6.0**. The abundance and bootstrap scripts require the R package [`quantreg`](https://cran.r-project.org/package=quantreg). The SPSS verification scripts require IBM SPSS Statistics.

The R scripts expect a `datos/` directory alongside `sintaxis/`, with the following input files:

- `PERDIZROJA_extraido.csv` — archived occurrence and environmental predictors used to reconstruct the favourability model.
- `transectos_abundancia_original.csv` — transect records used in the original abundance analysis.
- `transectos_tendencias_incluye_ceros.csv` — transect records including zero detections.
- `favorabilidad_espana.csv` — grid-level environmental favourability predictions.
- `celdas_modelo_original.csv` — original cell-level abundance-index data used for comparison.
- `datos_celdas_686.csv` — revised cell-level input for quantile regression and bootstrap.
- `bootstrap_bloques_rep.csv` — the saved 1,000 block-resampling selections used to reproduce the bootstrap intervals.

**These input datasets are not included in this public repository.** The field observations originate from the Observatorio Cinegético programme and are not being redistributed here. Requests for access should be addressed to the programme's data custodians and are subject to their data-sharing conditions. Consequently, the repository documents the analytical workflow, but the scripts cannot currently be executed end-to-end from the public files alone.

## Running the scripts with authorised data

Place the required input files in `datos/` and create empty `resultados/` and `figuras/` directories at the repository root:

```text
Red-legged-patridge/
├── README.md
├── sintaxis/
│   ├── 01_MODELO_AMBIENTAL_20_VARIABLES.R
│   ├── 02_VERIFICAR_MODELO_20_VARIABLES.sps
│   ├── 03_HISTORIAL_SPSS_ORIGINAL_RECUPERADO_NO_EJECUTAR.sps
│   ├── 04_ABUNDANCIA_Y_CUANTILES.R
│   └── 05_BOOTSTRAP_ESPACIAL_IC95.R
├── datos/       # restricted input files; not published
├── resultados/  # output tables
└── figuras/     # output figures
```

In R, install `quantreg` if needed:

```r
install.packages("quantreg")
```

Then open and **Source** the R scripts in RStudio, in this order: `01_MODELO_AMBIENTAL_20_VARIABLES.R`, `04_ABUNDANCIA_Y_CUANTILES.R`, and `05_BOOTSTRAP_ESPACIAL_IC95.R`. The scripts write outputs to `resultados/` and `figuras/`. The SPSS scripts serve as verification and documentation; they are not prerequisites for running the R analyses.

The bootstrap script uses the stored block-resampling matrix rather than generating new resamples; this is necessary to reproduce the reported percentile intervals.

## Interpretation

The abundance response is a **transect-derived index expressed as individuals per nominal 100 km² grid cell**, not a complete census or a detectability-corrected estimate. Q50 and Q70 are distinct conditional-quantile projection scenarios, **not confidence bounds** for the number of wild breeding birds. The spatial bootstrap intervals describe variation conditional on the fitted modelling workflow and fixed favourability surface; they do not account for all sources of sampling bias, imperfect detection, environmental-model uncertainty, or releases of farm-reared birds.

The occurrence model's zero-coded category includes non-detections and background/pseudo-absence records; these should not all be interpreted as confirmed biological absences.

## Contact

For questions about the analysis, please contact the corresponding author through the manuscript or the repository's GitHub issue tracker. For access to original survey records, contact the Observatorio Cinegético data custodians.
