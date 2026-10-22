# Core analysis functions. Used by the Shiny app (app.R), the scripted
# Quarto pipeline (analysis/analysis.qmd) and the tests. No Shiny code here.

# 1. Data -------------------------------------------------------------------

#' Load CATALYST's bundled example (subsampled Bodenmiller BCR-XL, PBMC)
load_example <- function() {
  env <- new.env()
  utils::data("PBMC_fs", "PBMC_panel", "PBMC_md", package = "CATALYST", envir = env)
  md <- as.data.frame(env$PBMC_md, stringsAsFactors = FALSE)
  md[] <- lapply(md, as.character)
  # Two artificial acquisition batches so the batch module can be demonstrated
  pts <- unique(md$patient_id)
  md$batch <- ifelse(md$patient_id %in% pts[seq_len(ceiling(length(pts) / 2))], "B1", "B2")
  fs <- env$PBMC_fs
  # Names inside the bundled flowSet differ between CATALYST versions:
  # align them with the metadata so prepData can match files.
  ids <- basename(as.character(flowCore::fsApply(fs, flowCore::identifier)))
  if (all(md$file_name %in% ids)) {
    flowCore::sampleNames(fs) <- ids            # e.g. "V1" -> "PBMC_patient1_BCRXL.fcs"
  } else {
    if (length(fs) != nrow(md)) stop("Example flowSet and metadata have different lengths.")
    flowCore::sampleNames(fs) <- md$file_name
    fs <- sync_identifiers(fs)
  }
  list(fs = fs, panel = as.data.frame(env$PBMC_panel), md = md)
}

#' Read user FCS files + metadata + panel
load_user <- function(fcs_paths, fcs_names = basename(fcs_paths),
                      md_path, panel_path,
                      md_name = md_path, panel_name = panel_path) {
  fs <- flowCore::read.flowSet(fcs_paths, transformation = FALSE,
                               truncate_max_range = FALSE)
  flowCore::sampleNames(fs) <- fcs_names
  fs <- sync_identifiers(fs)
  md <- read_table_any(md_path, md_name)
  panel <- read_table_any(panel_path, panel_name)
  list(fs = fs, panel = panel, md = md)
}

#' All names a flowSet's frames can be matched by (sample names + FCS identifiers)
fs_names <- function(fs) {
  ids <- tryCatch(as.character(flowCore::fsApply(fs, flowCore::identifier)), error = function(e) character(0))
  unique(c(flowCore::sampleNames(fs), ids, basename(ids)))
}

#' Make frame identifiers equal to sample names (prepData matches on identifiers)
sync_identifiers <- function(fs) {
  for (i in seq_along(fs)) {
    ff <- fs[[i]]
    flowCore::identifier(ff) <- flowCore::sampleNames(fs)[i]
    fs[[i]] <- ff
  }
  fs
}

#' Validate metadata/panel against the flowSet; returns character vector of problems
validate_inputs <- function(fs, panel, md) {
  p <- character(0)
  need_md <- c("file_name", "sample_id", "condition")
  miss <- setdiff(need_md, colnames(md))
  if (length(miss)) p <- c(p, paste("Metadata missing column(s):", paste(miss, collapse = ", ")))
  need_pn <- c("fcs_colname", "antigen", "marker_class")
  miss <- setdiff(need_pn, colnames(panel))
  if (length(miss)) p <- c(p, paste("Panel missing column(s):", paste(miss, collapse = ", ")))
  if (length(p)) return(p)
  if (anyDuplicated(md$sample_id)) p <- c(p, "sample_id values must be unique")
  nf <- setdiff(md$file_name, fs_names(fs))
  if (length(nf)) p <- c(p, paste("file_name not found among FCS files:", paste(head(nf, 5), collapse = ", ")))
  chans <- flowCore::colnames(fs)
  nc <- setdiff(panel$fcs_colname, chans)
  if (length(nc)) p <- c(p, paste("Panel channels not in FCS:", paste(head(nc, 5), collapse = ", ")))
  bad <- setdiff(unique(panel$marker_class), c("type", "state", "none"))
  if (length(bad)) p <- c(p, paste("marker_class must be type/state/none; found:", paste(bad, collapse = ", ")))
  p
}

#' Build the SingleCellExperiment (arcsinh transform, factor metadata)
build_sce <- function(fs, panel, md, cofactor = 5, factors = NULL) {
  probs <- validate_inputs(fs, panel, md)
  if (length(probs)) stop(paste(probs, collapse = "\n"), call. = FALSE)
  if (!"patient_id" %in% colnames(md)) md$patient_id <- md$sample_id
  if (!"batch" %in% colnames(md)) md$batch <- "batch1"
  if (is.null(factors)) {
    factors <- setdiff(colnames(md), c("file_name", "sample_id"))
  }
  md <- md[md$file_name %in% fs_names(fs), , drop = FALSE]
  sce <- CATALYST::prepData(
    fs, panel, md,
    md_cols = list(file = "file_name", id = "sample_id", factors = factors),
    panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
    transform = TRUE, cofactor = cofactor, FACS = FALSE
  )
  sce
}

# 2. Clustering & DR --------------------------------------------------------

run_clustering <- function(sce, features = "type", xdim = 10, ydim = 10,
                           maxK = 20, seed = 1234) {
  CATALYST::cluster(sce, features = features, xdim = xdim, ydim = ydim,
                    maxK = maxK, seed = seed, verbose = FALSE)
}

#' Dimensionality reduction on a per-sample subsample (cells = per sample)
run_dr <- function(sce, dr = "UMAP", cells = 1000, features = "type", seed = 1234) {
  set.seed(seed)
  CATALYST::runDR(sce, dr = dr, cells = cells, features = features)
}

#' Indices of cells that carry a (non-NA) embedding
embedded_cells <- function(sce, dr) {
  m <- SingleCellExperiment::reducedDim(sce, dr)
  which(stats::complete.cases(m))
}

#' Merge metaclusters with a two-column table (old_cluster, new_cluster)
apply_merge <- function(sce, k, table, id = "merging1") {
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  colnames(table)[1:2] <- c("old_cluster", "new_cluster")
  table$new_cluster[!nzchar(table$new_cluster)] <- table$old_cluster[!nzchar(table$new_cluster)]
  CATALYST::mergeClusters(sce, k = k, table = table, id = id, overwrite = TRUE)
}

# 3. Batch correction -------------------------------------------------------
# All corrections are NON-destructive: assay "exprs" is never modified.
#  - harmony   -> reducedDim "HARMONY" (+ "UMAP_HARMONY")
#  - cycombine -> assay "exprs_cycombine" (+ "UMAP_CYCOMBINE")
#  - limma     -> assay "exprs_limma"     (+ "UMAP_LIMMA")
# Differential testing always uses "exprs" with batch as a covariate.

marker_pca <- function(sce, assay = "exprs", features = CATALYST::type_markers(sce), npcs = 20) {
  x <- t(SummarizedExperiment::assay(sce, assay)[features, , drop = FALSE])
  npcs <- min(npcs, ncol(x))
  pc <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = npcs)
  pc$x
}

#' UMAP for a subset of cells, returned as full-length matrix with NA rows
umap_subset <- function(mat, idx, seed = 1234) {
  set.seed(seed)
  emb <- uwot::umap(mat[idx, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, verbose = FALSE)
  out <- matrix(NA_real_, nrow(mat), 2, dimnames = list(NULL, c("UMAP1", "UMAP2")))
  out[idx, ] <- emb
  out
}

batch_correct <- function(sce, method = c("harmony", "cycombine", "limma"),
                          batch_var = "batch", npcs = 20, cy_xdim = 8, cy_ydim = 8,
                          cy_norm = "scale", ref_dr = "UMAP", seed = 1234) {
  method <- match.arg(method)
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!batch_var %in% colnames(cd)) stop("Batch variable not in colData: ", batch_var)
  if (length(unique(cd[[batch_var]])) < 2) stop("Batch variable has a single level; nothing to correct.")
  tm <- CATALYST::type_markers(sce)

  # Cells for the after-correction UMAP: same cells as the 'before' UMAP when possible
  idx <- if (ref_dr %in% SingleCellExperiment::reducedDimNames(sce)) embedded_cells(sce, ref_dr) else
    sort(sample(ncol(sce), min(ncol(sce), 2e4)))

  pca <- marker_pca(sce, features = tm, npcs = npcs)
  SingleCellExperiment::reducedDim(sce, "PCA_MARKERS") <- pca

  if (method == "harmony") {
    hm <- harmony::RunHarmony(pca, cd, vars_use = batch_var, verbose = FALSE)
    SingleCellExperiment::reducedDim(sce, "HARMONY") <- hm
    SingleCellExperiment::reducedDim(sce, "UMAP_HARMONY") <- umap_subset(hm, idx, seed)
    S4Vectors::metadata(sce)$batch <- list(method = "harmony", var = batch_var,
      corrected_space = "HARMONY", umap = "UMAP_HARMONY")
  } else if (method == "cycombine") {
    if (!requireNamespace("cyCombine", quietly = TRUE)) stop("Install cyCombine: remotes::install_github('biosurf/cyCombine')")
    ex <- t(SummarizedExperiment::assay(sce, "exprs")[tm, , drop = FALSE])
    df <- data.frame(ex, check.names = FALSE)
    df$batch <- as.character(cd[[batch_var]])
    df$sample <- as.character(cd$sample_id)
    df$id <- seq_len(nrow(df))
    corr <- cyCombine::batch_correct(df, markers = tm, xdim = cy_xdim, ydim = cy_ydim,
                                     norm_method = cy_norm, seed = seed)
    corr <- corr[order(corr$id), ]
    new <- SummarizedExperiment::assay(sce, "exprs")
    new[tm, ] <- t(as.matrix(corr[, tm]))
    SummarizedExperiment::assay(sce, "exprs_cycombine", withDimnames = FALSE) <- new
    space <- marker_pca(sce, assay = "exprs_cycombine", features = tm, npcs = npcs)
    SingleCellExperiment::reducedDim(sce, "UMAP_CYCOMBINE") <- umap_subset(space, idx, seed)
    S4Vectors::metadata(sce)$batch <- list(method = "cycombine", var = batch_var,
      corrected_space = "exprs_cycombine", umap = "UMAP_CYCOMBINE")
  } else {
    new <- limma::removeBatchEffect(SummarizedExperiment::assay(sce, "exprs"),
                                    batch = factor(cd[[batch_var]]))
    SummarizedExperiment::assay(sce, "exprs_limma", withDimnames = FALSE) <- new
    space <- marker_pca(sce, assay = "exprs_limma", features = tm, npcs = npcs)
    SingleCellExperiment::reducedDim(sce, "UMAP_LIMMA") <- umap_subset(space, idx, seed)
    S4Vectors::metadata(sce)$batch <- list(method = "limma", var = batch_var,
      corrected_space = "exprs_limma", umap = "UMAP_LIMMA")
  }
  sce
}

#' Low-dimensional space after correction (for diagnostics)
corrected_space <- function(sce, npcs = 20) {
  b <- S4Vectors::metadata(sce)$batch
  if (is.null(b)) return(NULL)
  if (b$method == "harmony") return(SingleCellExperiment::reducedDim(sce, "HARMONY"))
  marker_pca(sce, assay = b$corrected_space, npcs = npcs)
}

#' Batch diagnostics on a cell subsample: iLISI, cLISI, ASW(batch), ASW(cell type)
batch_diagnostics <- function(sce, batch_var, label_k = NULL, n = 5000, npcs = 20, seed = 1234) {
  after <- corrected_space(sce, npcs)
  if (is.null(after)) stop("Run a batch correction first.")
  before <- if ("PCA_MARKERS" %in% SingleCellExperiment::reducedDimNames(sce))
    SingleCellExperiment::reducedDim(sce, "PCA_MARKERS") else marker_pca(sce, npcs = npcs)
  set.seed(seed)
  idx <- sort(sample(ncol(sce), min(n, ncol(sce))))
  cd <- as.data.frame(SummarizedExperiment::colData(sce))[idx, , drop = FALSE]
  meta <- data.frame(batch = as.character(cd[[batch_var]]))
  if (!is.null(label_k)) meta$celltype <- as.character(CATALYST::cluster_ids(sce, label_k)[idx])
  labs <- intersect(c("batch", "celltype"), colnames(meta))

  one <- function(X, stage) {
    X <- X[idx, , drop = FALSE]
    out <- data.frame(stage = stage, iLISI = NA_real_, cLISI = NA_real_)[rep(1, nrow(X)), ]
    if (requireNamespace("lisi", quietly = TRUE)) {
      l <- lisi::compute_lisi(X, meta, labs)
      out$iLISI <- l$batch
      if ("celltype" %in% labs) out$cLISI <- l$celltype
    }
    d <- stats::dist(X[seq_len(min(2000, nrow(X))), , drop = FALSE])
    sub <- seq_len(min(2000, nrow(X)))
    asw_b <- mean(cluster::silhouette(as.integer(factor(meta$batch[sub])), d)[, 3])
    asw_c <- if ("celltype" %in% labs) mean(cluster::silhouette(as.integer(factor(meta$celltype[sub])), d)[, 3]) else NA
    list(cells = out, summary = data.frame(stage = stage,
      median_iLISI = stats::median(out$iLISI), median_cLISI = stats::median(out$cLISI),
      ASW_batch = asw_b, ASW_celltype = asw_c, n_batches = length(unique(meta$batch))))
  }
  b <- one(before, "Before"); a <- one(after, "After")
  list(cells = rbind(b$cells, a$cells), summary = rbind(b$summary, a$summary))
}

# 4. Differential analysis --------------------------------------------------

#' Build design / formula / contrast for diffcyt
#'   condition_var : factor of interest; ref / comp : its two levels to compare
#'   covariates    : extra fixed effects (e.g. "batch")
#'   random        : random-effect column (e.g. "patient_id") or NULL
make_model <- function(sce, condition_var, ref, comp, covariates = NULL, random = NULL) {
  ei <- get_ei(sce)
  for (v in c(condition_var, covariates, random)) {
    if (!v %in% colnames(ei)) stop("Column not in experiment_info: ", v)
    ei[[v]] <- factor(ei[[v]])
  }
  lv <- levels(ei[[condition_var]])
  if (!all(c(ref, comp) %in% lv)) stop("Levels must be among: ", paste(lv, collapse = ", "))
  ei[[condition_var]] <- stats::relevel(ei[[condition_var]], ref = ref)
  covariates <- setdiff(covariates, c(condition_var, random))
  covariates <- covariates[vapply(covariates, function(v) length(unique(ei[[v]])) > 1, logical(1))]
  fixed <- c(condition_var, covariates)

  mm <- stats::model.matrix(stats::reformulate(fixed), data = ei)
  if (qr(mm)$rank < ncol(mm)) stop("Design is not full rank (a covariate is confounded with the condition). Remove it.")
  target <- paste0(condition_var, comp)
  cvec <- as.numeric(colnames(mm) == target)
  if (sum(cvec) != 1) stop("Could not locate coefficient ", target)

  design <- diffcyt::createDesignMatrix(ei, cols_design = fixed)
  formula <- diffcyt::createFormula(ei, cols_fixed = fixed,
                                    cols_random = if (!is.null(random) && random != "none") random else NULL)
  list(ei = ei, design = design, formula = formula, contrast = diffcyt::createContrast(cvec),
       condition_var = condition_var, ref = ref, comp = comp, fixed = fixed,
       random = if (!is.null(random) && random != "none") random else NULL)
}

#' Differential abundance. method: "GLMM" (recommended, uses random effect), "edgeR", "voom"
run_da <- function(sce, model, k, method = c("GLMM", "edgeR", "voom")) {
  method <- match.arg(method)
  S4Vectors::metadata(sce)$experiment_info <- model$ei
  args <- list(d_input = sce, contrast = model$contrast, analysis_type = "DA",
               clustering_to_use = k, verbose = FALSE)
  if (method == "GLMM") {
    args$formula <- model$formula; args$method_DA <- "diffcyt-DA-GLMM"
  } else {
    args$design <- model$design
    args$method_DA <- paste0("diffcyt-DA-", method)
    if (method == "voom" && !is.null(model$random))
      args$block_id <- model$ei[[model$random]]
  }
  res <- do.call(diffcyt::diffcyt, args)
  tab <- as.data.frame(SummarizedExperiment::rowData(res$res))
  tab <- merge(tab, da_effect_sizes(res, model), by = "cluster_id", sort = FALSE)
  list(res = res, table = tab[order(tab$p_adj), ], model = model, k = k, method = method)
}

#' Mean cluster proportions per group and log2 fold change (comp / ref)
da_effect_sizes <- function(res, model) {
  cnt <- SummarizedExperiment::assay(res$d_counts)
  prop <- sweep(cnt, 2, colSums(cnt), "/")
  grp <- model$ei[[model$condition_var]][match(colnames(prop), model$ei$sample_id)]
  m_ref <- rowMeans(prop[, grp == model$ref, drop = FALSE])
  m_cmp <- rowMeans(prop[, grp == model$comp, drop = FALSE])
  eps <- 1e-6
  data.frame(cluster_id = rownames(prop),
             prop_ref = m_ref, prop_comp = m_cmp,
             log2FC_prop = log2((m_cmp + eps) / (m_ref + eps)))
}

#' Differential state. method: "LMM" (recommended with random effect) or "limma"
run_ds <- function(sce, model, k, method = c("LMM", "limma")) {
  method <- match.arg(method)
  if (!length(CATALYST::state_markers(sce))) stop("No 'state' markers in the panel; DS needs marker_class = 'state'.")
  S4Vectors::metadata(sce)$experiment_info <- model$ei
  args <- list(d_input = sce, contrast = model$contrast, analysis_type = "DS",
               clustering_to_use = k, verbose = FALSE)
  if (method == "LMM") {
    args$formula <- model$formula; args$method_DS <- "diffcyt-DS-LMM"
  } else {
    args$design <- model$design; args$method_DS <- "diffcyt-DS-limma"
    if (!is.null(model$random)) args$block_id <- model$ei[[model$random]]
  }
  res <- do.call(diffcyt::diffcyt, args)
  tab <- as.data.frame(SummarizedExperiment::rowData(res$res))
  tab <- merge(tab, ds_effect_sizes(res, model), by = c("cluster_id", "marker_id"), sort = FALSE)
  list(res = res, table = tab[order(tab$p_adj), ], model = model, k = k, method = method)
}

#' Difference in mean (per-sample median) expression, comp - ref
ds_effect_sizes <- function(res, model) {
  dm <- res$d_medians
  mk <- SummarizedExperiment::assayNames(dm)
  grp <- model$ei[[model$condition_var]][match(colnames(dm), model$ei$sample_id)]
  out <- lapply(mk, function(m) {
    a <- SummarizedExperiment::assay(dm, m)
    data.frame(cluster_id = rownames(a), marker_id = m,
               median_ref = rowMeans(a[, grp == model$ref, drop = FALSE], na.rm = TRUE),
               median_comp = rowMeans(a[, grp == model$comp, drop = FALSE], na.rm = TRUE))
  })
  out <- do.call(rbind, out)
  out$delta_median <- out$median_comp - out$median_ref
  out
}

# 5. MEM & ROC --------------------------------------------------------------

#' Marker Enrichment Modeling (Diggins et al. 2017), computed natively.
#' Population vs all remaining cells; values scaled to -10..+10.
compute_mem <- function(sce, k, features = CATALYST::type_markers(sce), iqr_thresh = NULL) {
  x <- t(SummarizedExperiment::assay(sce, "exprs")[features, , drop = FALSE])
  x <- pmax(x, 0)
  cl <- CATALYST::cluster_ids(sce, k)
  iqr_all <- apply(x, 2, stats::IQR)
  thr <- if (is.null(iqr_thresh)) 0.5 else iqr_thresh
  mem <- t(vapply(levels(cl), function(c) {
    inp <- cl == c
    mag_p <- apply(x[inp, , drop = FALSE], 2, stats::median)
    mag_r <- apply(x[!inp, , drop = FALSE], 2, stats::median)
    iqr_p <- pmax(apply(x[inp, , drop = FALSE], 2, stats::IQR), thr)
    iqr_r <- pmax(apply(x[!inp, , drop = FALSE], 2, stats::IQR), thr)
    raw <- abs(mag_p - mag_r) + iqr_r / iqr_p - 1
    raw[mag_p - mag_r < 0] <- -raw[mag_p - mag_r < 0]
    raw
  }, numeric(length(features))))
  colnames(mem) <- features
  mem <- mem / max(abs(mem), na.rm = TRUE) * 10
  mem
}

mem_labels <- function(mem, n = 4) {
  vapply(rownames(mem), function(r) {
    v <- round(mem[r, ]); v <- v[order(-abs(v))][seq_len(min(n, length(v)))]
    paste0(names(v), ifelse(v >= 0, "+", "-"), abs(v), collapse = " ")
  }, character(1))
}

#' Per-marker AUC for one cluster vs rest, on a balanced subsample
compute_roc <- function(sce, k, cluster, features = CATALYST::type_markers(sce),
                        n_max = 5000, seed = 1234) {
  cl <- CATALYST::cluster_ids(sce, k)
  set.seed(seed)
  pos <- which(cl == cluster); neg <- which(cl != cluster)
  if (!length(pos) || !length(neg)) stop("Cluster has no cells or contains all cells.")
  idx <- c(pos[sample.int(length(pos), min(n_max, length(pos)))],
           neg[sample.int(length(neg), min(n_max, length(neg)))])
  y <- cl[idx] == cluster
  ex <- SummarizedExperiment::assay(sce, "exprs")[features, idx, drop = FALSE]
  curves <- list(); auc <- numeric(0)
  for (m in features) {
    r <- pROC::roc(y, as.numeric(ex[m, ]), quiet = TRUE, direction = "<", levels = c(FALSE, TRUE))
    auc[m] <- as.numeric(pROC::auc(r))
    curves[[m]] <- data.frame(marker = m, FPR = 1 - r$specificities, TPR = r$sensitivities)
  }
  tab <- data.frame(marker = names(auc), AUC = auc)
  tab$AUC_discriminative <- pmax(tab$AUC, 1 - tab$AUC)
  tab$direction <- ifelse(tab$AUC >= 0.5, "high in cluster", "low in cluster")
  list(table = tab[order(-tab$AUC_discriminative), ], curves = do.call(rbind, curves))
}

# 6. Export -----------------------------------------------------------------

cell_table <- function(sce) {
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  for (k in clustering_names(sce)) cd[[k]] <- CATALYST::cluster_ids(sce, k)
  cd
}

abundance_table <- function(sce, k) {
  cl <- CATALYST::cluster_ids(sce, k)
  tb <- as.data.frame(table(sample_id = sce$sample_id, cluster_id = cl))
  tb$prop <- stats::ave(tb$Freq, tb$sample_id, FUN = function(x) x / sum(x))
  tb
}
