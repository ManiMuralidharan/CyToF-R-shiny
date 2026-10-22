# Shared helpers ------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || all(is.na(a)) || (is.character(a) && !any(nzchar(a)))) b else a
}

#' Read a .csv / .tsv / .xlsx table, deciding by the original file name
read_table_any <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  df <- switch(ext,
    csv  = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv  = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    txt  = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    xlsx = as.data.frame(readxl::read_excel(path)),
    xls  = as.data.frame(readxl::read_excel(path)),
    stop("Unsupported table format: .", ext)
  )
  df[] <- lapply(df, function(x) if (is.character(x)) trimws(x) else x)
  df
}

#' Experiment info (one row per sample) from a CATALYST SCE
get_ei <- function(sce) {
  ei <- S4Vectors::metadata(sce)$experiment_info
  if (is.null(ei)) {
    ei <- unique(as.data.frame(SummarizedExperiment::colData(sce)))
  }
  as.data.frame(ei, stringsAsFactors = FALSE)
}

#' Metadata columns usable as factors (excluding ids)
ei_factor_cols <- function(sce) {
  setdiff(colnames(get_ei(sce)), c("n_cells", "file_name"))
}

#' Names of clusterings available (meta2..metaK plus merges)
clustering_names <- function(sce) {
  cc <- tryCatch(CATALYST::cluster_codes(sce), error = function(e) NULL)
  if (is.null(cc)) character(0) else setdiff(colnames(cc), "som100")
}

#' Draw a ComplexHeatmap object (or print ggplot) inside renderPlot
draw_any <- function(p) {
  if (inherits(p, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(p) else print(p)
  invisible(NULL)
}
