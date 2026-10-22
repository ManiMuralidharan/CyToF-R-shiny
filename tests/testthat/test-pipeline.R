suppressPackageStartupMessages({ library(CATALYST); library(SingleCellExperiment) })
root <- normalizePath(file.path("..", ".."))
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)

d <- load_example()

test_that("inputs validate and SCE builds", {
  expect_length(validate_inputs(d$fs, d$panel, d$md), 0)
  bad <- d$md; bad$sample_id <- NULL
  expect_true(length(validate_inputs(d$fs, d$panel, bad)) > 0)
  sce <- build_sce(d$fs, d$panel, d$md)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true(all(c("condition", "patient_id", "batch") %in% colnames(get_ei(sce))))
})

sce <- build_sce(d$fs, d$panel, d$md)
sce <- run_clustering(sce, maxK = 8)
sce <- run_dr(sce, "UMAP", cells = 200)

test_that("clustering, merge and embedding work", {
  expect_true("meta8" %in% clustering_names(sce))
  expect_true(length(embedded_cells(sce, "UMAP")) > 0)
  tab <- data.frame(old_cluster = 1:8, new_cluster = rep(c("A", "B"), 4))
  s2 <- apply_merge(sce, "meta8", tab, "m1")
  expect_setequal(levels(cluster_ids(s2, "m1")), c("A", "B"))
})

test_that("model is built with correct contrast", {
  m <- make_model(sce, "condition", "Ref", "BCRXL", covariates = "batch", random = "patient_id")
  expect_equal(sum(m$contrast), 1)
  expect_equal(levels(m$ei$condition)[1], "Ref")
  expect_error(make_model(sce, "condition", "Ref", "nope"))
})

test_that("DA and DS run and add effect sizes", {
  m <- make_model(sce, "condition", "Ref", "BCRXL", random = "patient_id")
  da <- run_da(sce, m, "meta8", "edgeR")
  expect_true(all(c("p_adj", "log2FC_prop") %in% colnames(da$table)))
  ds <- run_ds(sce, m, "meta8", "limma")
  expect_true("delta_median" %in% colnames(ds$table))
})

test_that("batch correction is non-destructive and diagnostics run", {
  before <- assay(sce, "exprs")
  s2 <- batch_correct(sce, "limma", "batch", npcs = 5)
  expect_identical(assay(s2, "exprs"), before)
  expect_true("exprs_limma" %in% assayNames(s2))
  dg <- batch_diagnostics(s2, "batch", label_k = "meta8", n = 800, npcs = 5)
  expect_equal(nrow(dg$summary), 2)
})

test_that("MEM and ROC", {
  mem <- compute_mem(sce, "meta8")
  expect_equal(dim(mem), c(8L, length(type_markers(sce))))
  expect_lte(max(abs(mem)), 10 + 1e-8)
  r <- compute_roc(sce, "meta8", "1", n_max = 300)
  expect_true(all(r$table$AUC >= 0 & r$table$AUC <= 1))
})
