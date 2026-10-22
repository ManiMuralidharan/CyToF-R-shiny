# CyTOF Analysis App

An interactive Shiny application and a matching scripted Quarto pipeline for mass and flow cytometry analysis, built on [CATALYST](https://bioconductor.org/packages/CATALYST) and [diffcyt](https://bioconductor.org/packages/diffcyt) and following the Nowicka et al. (F1000Research, 2017/2019) workflow.

Both interfaces share the same functions in `R/`, so an analysis explored in the app can be reproduced exactly from a config file.

## Modules

| Tab | What it does | Key methods |
|---|---|---|
| 1. Data | Load FCS + metadata + panel, validate, arcsinh transform | `flowCore`, `CATALYST::prepData` |
| 2. QC & EDA | Marker densities, cell counts, NRS, pseudobulk MDS, median heatmap | `plotExprs`, `plotNRS`, `pbMDS` |
| 3. Clustering & DR | FlowSOM + ConsensusClusterPlus, delta-area, UMAP / t-SNE / PCA, marker distributions | `CATALYST::cluster`, `runDR` |
| 4. Merging | Editable or uploaded merge/annotation table, downloadable | `mergeClusters` |
| 5. Batch | Harmony, cyCombine or limma; iLISI, cLISI, ASW before vs after | `harmony`, `cyCombine`, `lisi` |
| 6. Diff. abundance | GLMM with random effects (default), edgeR, voom; covariates; volcano, forest, paired plots, p-value histogram | `diffcyt` |
| 7. Diff. state | LMM with random effects (default) or limma | `diffcyt` |
| 8. MEM / ROC | Marker Enrichment Modeling labels, per-marker AUC | native MEM, `pROC` |
| 9. Export & Report | SCE `.rds`, per-cell CSV, abundances, parameter log (YAML), HTML report | `rmarkdown` |

### Design decisions

- **Batch correction never overwrites the raw `exprs` assay.** Corrected values go to a new assay or embedding and are used for visual inspection. For differential testing, add `batch` as a covariate; correcting first and then testing can remove or invent differences.
- **Paired designs are modelled properly.** Pick `patient_id` as the random effect: GLMM (DA) and LMM (DS) use it as a random intercept; edgeR/voom/limma fall back to fixed effects or blocking.
- **Contrasts are explicit.** You choose the reference and comparison levels, so results never depend on alphabetical level order.
- **Every run is logged.** The parameter log is stored in the exported SCE and can be downloaded as YAML.

## Install

R ≥ 4.4 and Bioconductor ≥ 3.20 are recommended.

```bash
git clone https://github.com/ManiMuralidharan/cytof-shiny.git
cd cytof-shiny
Rscript install.R
```

`lisi` and `cyCombine` come from GitHub and are optional; the related features show a message if they are missing.

For a pinned environment after installing:

```r
renv::init(bare = TRUE); renv::snapshot()
```

Commit the resulting `renv.lock`.

## Run the app

```r
shiny::runApp()          # from the repository root
```

Choose **Example** on the Data tab to try the bundled Bodenmiller BCR-XL subset (8 paired samples, reference vs BCR-XL stimulation; an artificial two-level `batch` column is added only to demonstrate the batch tab).

Or with Docker:

```bash
docker build -t cytof-shiny .
docker run -p 3838:3838 cytof-shiny     # open http://localhost:3838
```

## Your own data

Put files under `data/` (git-ignored). Templates are in `inst/templates/`.

**Metadata** (one row per FCS file):

| column | required | notes |
|---|---|---|
| `file_name` | yes | exact FCS file name |
| `sample_id` | yes | unique |
| `condition` | yes | factor of interest |
| `patient_id` | recommended | enables paired/random-effect models |
| `batch` | optional | acquisition batch |
| others | optional | usable as covariates |

**Panel** (one row per channel to keep): `fcs_colname`, `antigen`, `marker_class` (`type` = clustering markers, `state` = functional markers for DS, `none` = ignored).

Pre-gate to live single cells (and debarcode / bead-normalise) before loading; the app starts from clean FCS files.

## Scripted analysis (for papers)

Edit `analysis/config.yml`, then:

```bash
quarto render analysis/analysis.qmd --execute-dir .
```

Outputs go to `results/`: figures as PDF and 300-dpi PNG, DA/DS/MEM CSV tables, `sce.rds`, `abundances.csv` and `sessionInfo.txt`.

## Tests

```bash
Rscript tests/testthat.R
```

Tests run the full pipeline on the example data. GitHub Actions runs them on every push (`.github/workflows/check.yml`).

## Repository layout

```
app.R                 Shiny UI + server
R/utils.R             helpers
R/pipeline.R          analysis functions (load, cluster, merge, batch, DA, DS, MEM, ROC)
R/plots.R             extra plots (volcano, forest, paired, LISI, MEM, ROC)
analysis/             config.yml + analysis.qmd (scripted pipeline)
inst/report/          HTML report template used by the app
inst/templates/       metadata, panel and merge table templates
tests/                testthat suite
install.R, Dockerfile, .github/workflows/
```

## Citation

If you use this tool, please cite the underlying methods:

- Nowicka M. et al. CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research* (2017).
- Crowell H.L. et al. CATALYST (Bioconductor).
- Weber L.M. et al. diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Communications Biology* (2019).
- Van Gassen S. et al. FlowSOM. *Cytometry A* (2015).
- Korsunsky I. et al. Harmony. *Nature Methods* (2019).
- Pedersen C.B. et al. cyCombine. *Nature Communications* (2022).
- Diggins K.E. et al. Marker Enrichment Modeling. *Nature Methods* (2017).

## License

MIT, see `LICENSE`.
Dr. Muralidharan Mani 
University of Wisconsin-Madison