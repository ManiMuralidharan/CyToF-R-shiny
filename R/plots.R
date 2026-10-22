# Plot helpers not provided by CATALYST --------------------------------------

theme_cytof <- function() ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

plot_volcano <- function(tab, x, fdr = 0.05, label = "cluster_id", xlab = x) {
  tab$sig <- !is.na(tab$p_adj) & tab$p_adj < fdr
  tab$lab <- if (all(c("cluster_id", "marker_id") %in% colnames(tab)))
    paste(tab$cluster_id, tab$marker_id, sep = ":") else tab[[label]]
  top <- utils::head(tab[tab$sig, ][order(tab$p_adj[tab$sig]), ], 15)
  ggplot2::ggplot(tab, ggplot2::aes(.data[[x]], -log10(p_val), colour = sig)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#c0392b"),
                                 name = paste0("FDR < ", fdr)) +
    ggplot2::geom_text(data = top, ggplot2::aes(label = lab), vjust = -0.6, size = 3, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(x = xlab, y = expression(-log[10]~p)) + theme_cytof()
}

plot_pval_hist <- function(tab) {
  ggplot2::ggplot(tab, ggplot2::aes(p_val)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, 0.05), fill = "#3c8dbc", colour = "white") +
    ggplot2::labs(x = "raw p-value", y = "count",
                  subtitle = "Should be flat with a spike near 0 if signal is present") + theme_cytof()
}

plot_da_forest <- function(tab, fdr = 0.05) {
  tab <- tab[order(tab$log2FC_prop), ]
  tab$cluster_id <- factor(tab$cluster_id, levels = tab$cluster_id)
  tab$sig <- !is.na(tab$p_adj) & tab$p_adj < fdr
  ggplot2::ggplot(tab, ggplot2::aes(log2FC_prop, cluster_id, colour = sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#c0392b"), name = paste0("FDR < ", fdr)) +
    ggplot2::labs(x = "log2 fold change of mean proportion", y = NULL) + theme_cytof()
}

#' Paired slope plot of cluster proportions per patient
plot_paired_props <- function(sce, k, model, clusters = NULL) {
  ab <- abundance_table(sce, k)
  ei <- model$ei
  ab$group <- ei[[model$condition_var]][match(ab$sample_id, ei$sample_id)]
  ab <- ab[ab$group %in% c(model$ref, model$comp), ]
  ab$group <- factor(ab$group, levels = c(model$ref, model$comp))
  if (!is.null(clusters)) ab <- ab[ab$cluster_id %in% clusters, ]
  pid <- if (!is.null(model$random)) model$random else "sample_id"
  ab$pair <- ei[[pid]][match(ab$sample_id, ei$sample_id)]
  ggplot2::ggplot(ab, ggplot2::aes(group, prop * 100)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey95") +
    { if (!is.null(model$random)) ggplot2::geom_line(ggplot2::aes(group = pair), colour = "grey50", alpha = 0.6) } +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = 2) +
    ggplot2::facet_wrap(~cluster_id, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "% of sample cells") + theme_cytof() +
    ggplot2::theme(legend.position = "none")
}

plot_cluster_sizes <- function(sce, k) {
  df <- as.data.frame(table(cluster = CATALYST::cluster_ids(sce, k)))
  ggplot2::ggplot(df, ggplot2::aes(stats::reorder(cluster, Freq), Freq)) +
    ggplot2::geom_col(fill = "#3c8dbc") + ggplot2::coord_flip() +
    ggplot2::labs(x = "Cluster", y = "Cells") + theme_cytof()
}

plot_batch_umaps <- function(sce, batch_var, before = "UMAP") {
  b <- S4Vectors::metadata(sce)$batch
  p1 <- CATALYST::plotDR(sce, before, color_by = batch_var) + ggplot2::ggtitle("Before correction")
  if (is.null(b)) return(p1)
  p2 <- CATALYST::plotDR(sce, b$umap, color_by = batch_var) +
    ggplot2::ggtitle(paste("After", b$method))
  cowplot::plot_grid(p1, p2, ncol = 2)
}

plot_lisi <- function(diag) {
  d <- diag$cells
  d$stage <- factor(d$stage, levels = c("Before", "After"))
  long <- rbind(data.frame(stage = d$stage, metric = "iLISI (batch mixing, higher = better)", value = d$iLISI),
                data.frame(stage = d$stage, metric = "cLISI (cell-type purity, ~1 = better)", value = d$cLISI))
  long <- long[!is.na(long$value), ]
  ggplot2::ggplot(long, ggplot2::aes(value, fill = stage)) +
    ggplot2::geom_density(alpha = 0.5) + ggplot2::facet_wrap(~metric, scales = "free") +
    ggplot2::labs(x = "LISI", y = "density") + theme_cytof()
}

plot_mem <- function(mem) {
  ComplexHeatmap::Heatmap(round(mem, 1), name = "MEM",
    col = circlize::colorRamp2(c(-10, 0, 10), c("#2166ac", "white", "#b2182b")),
    row_labels = paste0(rownames(mem), "  ", mem_labels(mem, 3)),
    row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9),
    cell_fun = function(j, i, x, y, w, h, f) grid::grid.text(round(mem[i, j]), x, y, gp = grid::gpar(fontsize = 7)))
}

plot_roc <- function(roc, n = 6) {
  top <- utils::head(roc$table$marker, n)
  d <- roc$curves[roc$curves$marker %in% top, ]
  lab <- roc$table[roc$table$marker %in% top, ]
  d$marker <- factor(d$marker, levels = top,
                     labels = paste0(top, " (AUC ", sprintf("%.2f", lab$AUC[match(top, lab$marker)]), ")"))
  ggplot2::ggplot(d, ggplot2::aes(FPR, TPR, colour = marker)) +
    ggplot2::geom_path() + ggplot2::geom_abline(linetype = "dashed") +
    ggplot2::coord_equal() + ggplot2::labs(colour = NULL) + theme_cytof()
}
