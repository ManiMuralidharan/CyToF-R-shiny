# Install all dependencies. Run once:  Rscript install.R
# For a locked, exact environment afterwards run: renv::init() then renv::snapshot()

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran <- c("shiny", "shinydashboard", "shinyWidgets", "DT", "ggplot2", "cowplot",
          "readxl", "yaml", "uwot", "harmony", "cluster", "pROC", "circlize",
          "rmarkdown", "knitr", "testthat", "lme4")
bioc <- c("CATALYST", "diffcyt", "flowCore", "SingleCellExperiment",
          "SummarizedExperiment", "ComplexHeatmap", "scater", "limma", "edgeR")
github <- c(lisi = "immunogenomics/lisi", cyCombine = "biosurf/cyCombine")

need <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
if (length(x <- need(cran))) install.packages(x)
if (length(x <- need(bioc))) BiocManager::install(x, update = FALSE, ask = FALSE)
for (p in names(github)) if (!requireNamespace(p, quietly = TRUE))
  tryCatch(remotes::install_github(github[[p]], upgrade = "never"),
           error = function(e) message("Optional package ", p, " not installed: ", conditionMessage(e)))

missing <- need(c(cran, bioc))
if (length(missing)) stop("Still missing: ", paste(missing, collapse = ", ")) else message("All required packages installed.")
