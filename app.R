# ============================================================================
# CyTOF Analysis Shiny App (v3)
# CATALYST + diffcyt workflow: Data > QC > Clustering > Merging > Batch >
# Differential abundance > Differential state > MEM/ROC > Report
# Run:  shiny::runApp()   (from the repository root)
# ============================================================================

suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(shinyWidgets); library(DT)
  library(ggplot2); library(CATALYST); library(diffcyt)
  library(SingleCellExperiment); library(SummarizedExperiment)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f, local = TRUE)
options(shiny.maxRequestSize = 4 * 1024^3)

box_p <- function(..., title, width = 12) box(..., title = title, width = width, status = "primary", solidHeader = TRUE)
dl_btn <- function(id, label = "PDF") downloadButton(id, label, class = "btn-xs")

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CyTOF Analysis v3"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("1. Data", tabName = "data", icon = icon("database")),
      menuItem("2. QC & EDA", tabName = "qc", icon = icon("chart-line")),
      menuItem("3. Clustering & DR", tabName = "cluster", icon = icon("project-diagram")),
      menuItem("4. Merging", tabName = "merge", icon = icon("object-group")),
      menuItem("5. Batch", tabName = "batch", icon = icon("layer-group")),
      menuItem("6. Diff. Abundance", tabName = "da", icon = icon("balance-scale")),
      menuItem("7. Diff. State", tabName = "ds", icon = icon("flask")),
      menuItem("8. MEM / ROC", tabName = "mem", icon = icon("dna")),
      menuItem("9. Export & Report", tabName = "report", icon = icon("file-export"))
    ),
    hr(), div(style = "padding:0 15px;font-size:11px;color:#b8c7ce;", textOutput("status_line"))
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".content-wrapper{background:#f4f6f9}.small-note{color:#666;font-size:12px}"))),
    tabItems(
      # ---- 1. Data ----
      tabItem("data",
        fluidRow(
          box_p(width = 4, title = "Data source",
            radioButtons("data_source", NULL, c("Example (CATALYST PBMC, Bodenmiller BCR-XL subset)" = "example",
                                                "Upload FCS + metadata + panel" = "upload")),
            conditionalPanel("input.data_source=='upload'",
              fileInput("fcs_files", "FCS files", multiple = TRUE, accept = ".fcs"),
              fileInput("meta_file", "Metadata (.csv/.tsv/.xlsx)", accept = c(".csv", ".tsv", ".xlsx")),
              fileInput("panel_file", "Panel (.csv/.tsv/.xlsx)", accept = c(".csv", ".tsv", ".xlsx")),
              p(class = "small-note", "Metadata needs file_name, sample_id, condition (optional: patient_id, batch, other covariates).",
                br(), "Panel needs fcs_colname, antigen, marker_class (type / state / none). Templates in inst/templates/.")
            ),
            numericInput("cofactor", "arcsinh cofactor (5 mass cytometry, 150 flow)", 5, 1, 1000),
            actionBttn("load_data", "Load & build SCE", style = "gradient", color = "primary", block = TRUE, icon = icon("play"))
          ),
          box_p(width = 8, title = "Inputs",
            tabsetPanel(
              tabPanel("Metadata", DTOutput("meta_tbl")),
              tabPanel("Panel", DTOutput("panel_tbl")),
              tabPanel("Experiment info", verbatimTextOutput("exp_info"))
            ))
        ),
        fluidRow(valueBoxOutput("vb_samples", 3), valueBoxOutput("vb_cells", 3),
                 valueBoxOutput("vb_type", 3), valueBoxOutput("vb_state", 3))
      ),
      # ---- 2. QC ----
      tabItem("qc",
        fluidRow(box_p(title = "QC controls",
          fluidRow(
            column(3, selectInput("qc_color_by", "Colour by", NULL)),
            column(9, pickerInput("qc_markers", "Markers for density plot", NULL, multiple = TRUE,
                                  options = list(`actions-box` = TRUE, `live-search` = TRUE)))
          ))),
        fluidRow(
          box_p(width = 7, title = tagList("Marker densities (arcsinh)", dl_btn("dl_qc_exprs")), plotOutput("qc_exprs", height = 520)),
          box_p(width = 5, title = tagList("Cell counts", dl_btn("dl_qc_counts")), plotOutput("qc_counts", height = 240),
                tagList(h5("Non-redundancy score"), plotOutput("qc_nrs", height = 240)))
        ),
        fluidRow(
          box_p(width = 6, title = tagList("Pseudobulk MDS", dl_btn("dl_qc_mds")), plotOutput("qc_mds", height = 380)),
          box_p(width = 6, title = tagList("Median expression (samples x markers)", dl_btn("dl_qc_heat")), plotOutput("qc_heat", height = 380))
        )
      ),
      # ---- 3. Clustering ----
      tabItem("cluster",
        fluidRow(
          box_p(width = 4, title = "FlowSOM + ConsensusClusterPlus",
            fluidRow(column(6, numericInput("xdim", "SOM x", 10, 3, 30)), column(6, numericInput("ydim", "SOM y", 10, 3, 30))),
            numericInput("maxK", "maxK (metaclusters)", 20, 3, 50),
            pickerInput("cluster_markers", "Clustering markers (default: type)", NULL, multiple = TRUE,
                        options = list(`actions-box` = TRUE)),
            selectInput("dr_method", "Embedding", c("UMAP", "TSNE", "PCA")),
            numericInput("dr_cells", "Cells per sample for embedding", 1000, 100, 1e5),
            numericInput("seed", "Seed", 1234),
            actionBttn("run_cluster", "Run clustering & embedding", color = "success", style = "gradient", block = TRUE, icon = icon("play")),
            hr(), selectInput("k_choose", "Show clustering", NULL),
            h5("Delta area (choose k at the elbow)"), plotOutput("delta_plot", height = 220)
          ),
          box_p(width = 8, title = "Cluster results",
            tabsetPanel(
              tabPanel("Heatmap", dl_btn("dl_cl_heat"), plotOutput("cl_heat", height = 620)),
              tabPanel("Embedding", fluidRow(column(4, selectInput("dr_color", "Colour by", NULL))),
                       dl_btn("dl_cl_dr"), plotOutput("cl_dr", height = 560)),
              tabPanel("Abundances", dl_btn("dl_cl_abund"), plotOutput("cl_abund", height = 520)),
              tabPanel("Cluster size", plotOutput("cl_size", height = 460)),
              tabPanel("Marker distributions", dl_btn("dl_cl_exprs"), plotOutput("cl_exprs", height = 700))
            ))
        )
      ),
      # ---- 4. Merging ----
      tabItem("merge",
        fluidRow(
          box_p(width = 4, title = "Manual merging / annotation",
            selectInput("merge_k", "Source clustering", NULL),
            textInput("merge_id", "New clustering name", "merging1"),
            fileInput("merge_upload", "Or upload a merge table (old_cluster,new_cluster)", accept = c(".csv", ".xlsx")),
            actionBttn("apply_merge", "Apply merge", color = "warning", style = "gradient", block = TRUE, icon = icon("object-group")),
            downloadButton("dl_merge_tbl", "Download merge table (CSV)"),
            p(class = "small-note", "Double-click a cell in 'new_cluster' to rename. Clusters sharing a name are merged.")
          ),
          box_p(width = 8, title = "Annotation",
            DTOutput("merge_tbl"), hr(), plotOutput("merge_heat", height = 450), plotOutput("merge_dr", height = 450))
        )
      ),
      # ---- 5. Batch ----
      tabItem("batch",
        fluidRow(
          box_p(width = 4, title = "Batch correction (non-destructive)",
            selectInput("batch_var", "Batch variable", NULL),
            radioButtons("batch_method", "Method", c("Harmony (PCA space)" = "harmony",
              "cyCombine (SOM ComBat)" = "cycombine", "limma removeBatchEffect" = "limma")),
            numericInput("batch_npcs", "PCs", 20, 2, 50),
            conditionalPanel("input.batch_method=='cycombine'",
              selectInput("cy_norm", "Normalisation", c("scale", "rank")),
              fluidRow(column(6, numericInput("cy_xdim", "SOM x", 8, 3, 20)), column(6, numericInput("cy_ydim", "SOM y", 8, 3, 20)))),
            actionBttn("run_batch", "Run correction", color = "warning", style = "gradient", block = TRUE, icon = icon("layer-group")),
            hr(), selectInput("diag_k", "Cell-type labels for cLISI / ASW", NULL),
            actionBttn("run_diag", "Compute diagnostics", color = "primary", style = "simple", block = TRUE),
            p(class = "small-note", "Original 'exprs' is never overwritten. Differential tests use the raw data with batch as a covariate (recommended).")
          ),
          box_p(width = 8, title = "Before / after",
            tabsetPanel(
              tabPanel("Embedding", verbatimTextOutput("batch_status"), dl_btn("dl_batch_umap"), plotOutput("batch_umap", height = 480)),
              tabPanel("LISI", plotOutput("batch_lisi", height = 420)),
              tabPanel("Metric summary", DTOutput("batch_diag_tbl"),
                       p(class = "small-note", "iLISI near number of batches = good mixing; cLISI near 1 and ASW_celltype high = biology preserved; ASW_batch near 0 = no batch structure."))
            ))
        )
      ),
      # ---- 6. DA ----
      tabItem("da",
        fluidRow(
          box_p(width = 4, title = "Differential abundance model",
            selectInput("da_k", "Clustering", NULL),
            selectInput("da_cond", "Condition", NULL),
            fluidRow(column(6, selectInput("da_ref", "Reference", NULL)), column(6, selectInput("da_comp", "Compare", NULL))),
            pickerInput("da_cov", "Fixed covariates (e.g. batch)", NULL, multiple = TRUE),
            selectInput("da_random", "Random effect (paired / repeated)", NULL),
            selectInput("da_method", "Method", c("GLMM (binomial, random effects)" = "GLMM", "edgeR" = "edgeR", "voom" = "voom")),
            numericInput("fdr_da", "FDR", 0.05, 0.001, 0.5, 0.01),
            actionBttn("run_da", "Run DA", color = "success", style = "gradient", block = TRUE)
          ),
          box_p(width = 8, title = "DA results",
            tabsetPanel(
              tabPanel("Table", downloadButton("dl_da_csv", "CSV"), DTOutput("da_tbl")),
              tabPanel("Heatmap", dl_btn("dl_da_heat"), plotOutput("da_heat", height = 520)),
              tabPanel("Volcano", plotOutput("da_volcano", height = 480)),
              tabPanel("Forest", plotOutput("da_forest", height = 480)),
              tabPanel("Proportions", dl_btn("dl_da_props"), plotOutput("da_props", height = 600)),
              tabPanel("p-value histogram", plotOutput("da_phist", height = 380))
            ))
        )
      ),
      # ---- 7. DS ----
      tabItem("ds",
        fluidRow(
          box_p(width = 4, title = "Differential state model",
            selectInput("ds_k", "Clustering", NULL),
            selectInput("ds_cond", "Condition", NULL),
            fluidRow(column(6, selectInput("ds_ref", "Reference", NULL)), column(6, selectInput("ds_comp", "Compare", NULL))),
            pickerInput("ds_cov", "Fixed covariates", NULL, multiple = TRUE),
            selectInput("ds_random", "Random effect", NULL),
            selectInput("ds_method", "Method", c("LMM (random effects)" = "LMM", "limma" = "limma")),
            numericInput("fdr_ds", "FDR", 0.05, 0.001, 0.5, 0.01),
            numericInput("ds_top", "Top hits in heatmap", 50, 5, 500),
            actionBttn("run_ds", "Run DS", color = "success", style = "gradient", block = TRUE)
          ),
          box_p(width = 8, title = "DS results",
            tabsetPanel(
              tabPanel("Table", downloadButton("dl_ds_csv", "CSV"), DTOutput("ds_tbl")),
              tabPanel("Heatmap", dl_btn("dl_ds_heat"), plotOutput("ds_heat", height = 620)),
              tabPanel("Volcano", plotOutput("ds_volcano", height = 480)),
              tabPanel("p-value histogram", plotOutput("ds_phist", height = 380))
            ))
        )
      ),
      # ---- 8. MEM / ROC ----
      tabItem("mem",
        fluidRow(
          box_p(width = 4, title = "Marker enrichment & ROC",
            selectInput("mem_k", "Clustering", NULL),
            actionBttn("run_mem", "Compute MEM", color = "success", style = "gradient", block = TRUE),
            hr(), selectInput("roc_cluster", "Cluster for ROC", NULL),
            actionBttn("run_roc", "Compute ROC / AUC", color = "primary", style = "simple", block = TRUE)
          ),
          box_p(width = 8, title = "Results",
            tabsetPanel(
              tabPanel("MEM heatmap", dl_btn("dl_mem"), plotOutput("mem_heat", height = 620)),
              tabPanel("MEM labels", DTOutput("mem_tbl")),
              tabPanel("ROC curves", plotOutput("roc_plot", height = 480)),
              tabPanel("AUC table", DTOutput("roc_tbl"))
            ))
        )
      ),
      # ---- 9. Export ----
      tabItem("report",
        fluidRow(box_p(title = "Export",
          fluidRow(
            column(3, downloadButton("dl_sce", "SingleCellExperiment (.rds)", style = "width:100%")),
            column(3, downloadButton("dl_cells", "Per-cell metadata + clusters (.csv.gz)", style = "width:100%")),
            column(3, downloadButton("dl_abund", "Cluster abundances (.csv)", style = "width:100%")),
            column(3, downloadButton("dl_report", "HTML report", style = "width:100%"))
          ),
          hr(), h4("Analysis log (parameters used)"), verbatimTextOutput("log_out"),
          downloadButton("dl_log", "Download log (.yml)"),
          hr(), h4("Session info"), verbatimTextOutput("session_info")))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(sce = NULL, md = NULL, panel = NULL, da = NULL, ds = NULL,
                       mem = NULL, roc = NULL, diag = NULL, log = list(),
                       batch_msg = "No correction run yet.", merge_map = NULL)

  # Run an expression with progress + error notification; returns NULL on error
  safely <- function(msg, expr) {
    withProgress(message = msg, value = 0.4, {
      tryCatch(expr, error = function(e) {
        showNotification(paste0(msg, " failed: ", conditionMessage(e)), type = "error", duration = 15)
        NULL
      })
    })
  }
  log_step <- function(name, params) rv$log[[name]] <- c(params, list(time = format(Sys.time())))
  pdf_handler <- function(name, plot_fun, w = 10, h = 7) downloadHandler(
    filename = function() paste0(name, "_", Sys.Date(), ".pdf"),
    content = function(file) { grDevices::pdf(file, width = w, height = h); on.exit(grDevices::dev.off()); draw_any(plot_fun()) })

  output$status_line <- renderText({
    if (is.null(rv$sce)) "No data loaded" else
      sprintf("%s cells | %d samples | %d clusterings", format(ncol(rv$sce), big.mark = ","),
              nrow(get_ei(rv$sce)), length(clustering_names(rv$sce)))
  })

  # ---- 1. Data ----
  observeEvent(input$load_data, {
    out <- safely("Loading data", {
      d <- if (input$data_source == "example") load_example() else {
        validate(need(input$fcs_files, "Upload FCS files"), need(input$meta_file, "Upload metadata"),
                 need(input$panel_file, "Upload panel"))
        load_user(input$fcs_files$datapath, input$fcs_files$name,
                  input$meta_file$datapath, input$panel_file$datapath,
                  input$meta_file$name, input$panel_file$name)
      }
      d$sce <- build_sce(d$fs, d$panel, d$md, cofactor = input$cofactor)
      d
    })
    req(out)
    rv$sce <- out$sce; rv$md <- out$md; rv$panel <- out$panel
    rv$da <- rv$ds <- rv$mem <- rv$roc <- rv$diag <- NULL; rv$log <- list()
    log_step("data", list(source = input$data_source, cofactor = input$cofactor,
                          n_cells = ncol(out$sce), n_samples = nrow(get_ei(out$sce))))
    sce <- out$sce; tm <- type_markers(sce); fac <- ei_factor_cols(sce)
    cdn <- intersect(c("condition", "sample_id", "patient_id", "batch", fac), colnames(colData(sce)))
    updateSelectInput(session, "qc_color_by", choices = unique(cdn), selected = "condition")
    updatePickerInput(session, "qc_markers", choices = rownames(sce), selected = head(tm, 8))
    updatePickerInput(session, "cluster_markers", choices = rownames(sce), selected = tm)
    updateSelectInput(session, "batch_var", choices = fac, selected = if ("batch" %in% fac) "batch" else fac[1])
    for (id in c("da", "ds")) {
      updateSelectInput(session, paste0(id, "_cond"), choices = fac, selected = "condition")
      updatePickerInput(session, paste0(id, "_cov"), choices = fac, selected = character(0))
      updateSelectInput(session, paste0(id, "_random"), choices = c("none", fac),
                        selected = if ("patient_id" %in% fac) "patient_id" else "none")
    }
    showNotification(sprintf("Loaded %s cells from %d samples", format(ncol(sce), big.mark = ","), nrow(get_ei(sce))), type = "message")
  })

  output$meta_tbl <- renderDT(datatable(req(rv$md), options = list(scrollX = TRUE, pageLength = 15)))
  output$panel_tbl <- renderDT(datatable(req(rv$panel), options = list(scrollX = TRUE, pageLength = 25)))
  output$exp_info <- renderPrint(get_ei(req(rv$sce)))
  vb <- function(val, lab, ic, col) valueBox(val, lab, icon = icon(ic), color = col)
  output$vb_samples <- renderValueBox(vb(if (is.null(rv$sce)) 0 else nrow(get_ei(rv$sce)), "Samples", "vial", "aqua"))
  output$vb_cells <- renderValueBox(vb(if (is.null(rv$sce)) 0 else format(ncol(rv$sce), big.mark = ","), "Cells", "users", "green"))
  output$vb_type <- renderValueBox(vb(if (is.null(rv$sce)) 0 else length(type_markers(rv$sce)), "Type markers", "tag", "yellow"))
  output$vb_state <- renderValueBox(vb(if (is.null(rv$sce)) 0 else length(state_markers(rv$sce)), "State markers", "flask", "purple"))

  # factor level selectors for DA / DS
  for (id in c("da", "ds")) local({
    id <- id
    observeEvent(input[[paste0(id, "_cond")]], {
      req(rv$sce); v <- input[[paste0(id, "_cond")]]; req(v)
      lv <- sort(unique(as.character(get_ei(rv$sce)[[v]])))
      updateSelectInput(session, paste0(id, "_ref"), choices = lv, selected = lv[1])
      updateSelectInput(session, paste0(id, "_comp"), choices = lv, selected = lv[min(2, length(lv))])
    })
  })

  # ---- 2. QC ----
  qc_exprs <- reactive(plotExprs(req(rv$sce), features = input$qc_markers %||% type_markers(rv$sce), color_by = req(input$qc_color_by)))
  qc_counts <- reactive(plotCounts(req(rv$sce), group_by = "sample_id", color_by = req(input$qc_color_by)))
  qc_mds <- reactive(pbMDS(req(rv$sce), color_by = req(input$qc_color_by), label_by = "sample_id"))
  qc_heat <- reactive(plotExprHeatmap(req(rv$sce), features = "type", by = "sample_id", scale = "last", bars = TRUE))
  output$qc_exprs <- renderPlot(qc_exprs())
  output$qc_counts <- renderPlot(qc_counts())
  output$qc_nrs <- renderPlot(plotNRS(req(rv$sce), features = "type", color_by = req(input$qc_color_by)))
  output$qc_mds <- renderPlot(qc_mds())
  output$qc_heat <- renderPlot(draw_any(qc_heat()))
  output$dl_qc_exprs <- pdf_handler("qc_densities", qc_exprs, 14, 10)
  output$dl_qc_counts <- pdf_handler("qc_counts", qc_counts)
  output$dl_qc_mds <- pdf_handler("qc_mds", qc_mds)
  output$dl_qc_heat <- pdf_handler("qc_heatmap", qc_heat, 12, 8)

  # ---- 3. Clustering ----
  observeEvent(input$run_cluster, {
    req(rv$sce)
    feats <- input$cluster_markers %||% "type"
    sce <- safely("Clustering", {
      s <- run_clustering(rv$sce, features = feats, xdim = input$xdim, ydim = input$ydim,
                          maxK = input$maxK, seed = input$seed)
      run_dr(s, dr = input$dr_method, cells = input$dr_cells, features = feats, seed = input$seed)
    })
    req(sce)
    rv$sce <- sce
    log_step("clustering", list(features = feats, xdim = input$xdim, ydim = input$ydim, maxK = input$maxK,
                                dr = input$dr_method, dr_cells_per_sample = input$dr_cells, seed = input$seed))
    refresh_k(selected = paste0("meta", input$maxK))
    cdn <- c(clustering_names(sce), intersect(c("condition", "sample_id", "patient_id", "batch"), colnames(colData(sce))), rownames(sce))
    updateSelectInput(session, "dr_color", choices = cdn, selected = paste0("meta", input$maxK))
    showNotification("Clustering and embedding done", type = "message")
  })

  refresh_k <- function(selected = NULL) {
    ks <- clustering_names(rv$sce)
    sel <- selected %||% tail(ks, 1)
    for (id in c("k_choose", "merge_k", "da_k", "ds_k", "mem_k", "diag_k"))
      updateSelectInput(session, id, choices = ks, selected = sel)
  }

  has_dr <- reactive({ req(rv$sce); input$dr_method %in% reducedDimNames(rv$sce) })
  output$delta_plot <- renderPlot({ req(rv$sce); d <- metadata(rv$sce)$delta_area; req(d); print(d) })
  cl_heat <- reactive(plotExprHeatmap(req(rv$sce), features = "type", by = "cluster_id",
                                      k = req(input$k_choose), scale = "first", q = 0.01, perc = TRUE, bars = TRUE))
  cl_dr <- reactive({ req(has_dr(), input$dr_color); plotDR(rv$sce, input$dr_method, color_by = input$dr_color) })
  cl_abund <- reactive(plotAbundances(req(rv$sce), k = req(input$k_choose), by = "sample_id", group_by = "condition"))
  cl_exprs <- reactive(plotClusterExprs(req(rv$sce), k = req(input$k_choose), features = "type"))
  output$cl_heat <- renderPlot(draw_any(cl_heat()))
  output$cl_dr <- renderPlot(cl_dr())
  output$cl_abund <- renderPlot(cl_abund())
  output$cl_size <- renderPlot(plot_cluster_sizes(req(rv$sce), req(input$k_choose)))
  output$cl_exprs <- renderPlot(cl_exprs())
  output$dl_cl_heat <- pdf_handler("cluster_heatmap", cl_heat, 12, 9)
  output$dl_cl_dr <- pdf_handler("embedding", cl_dr, 9, 7)
  output$dl_cl_abund <- pdf_handler("abundances", cl_abund, 12, 7)
  output$dl_cl_exprs <- pdf_handler("cluster_marker_distributions", cl_exprs, 14, 12)

  # ---- 4. Merging ----
  observeEvent(list(input$merge_k, rv$sce), {
    req(rv$sce, input$merge_k)
    cl <- levels(cluster_ids(rv$sce, input$merge_k))
    rv$merge_map <- data.frame(old_cluster = cl, new_cluster = cl, stringsAsFactors = FALSE)
  }, ignoreInit = TRUE)
  observeEvent(input$merge_upload, {
    t <- tryCatch(read_table_any(input$merge_upload$datapath, input$merge_upload$name), error = function(e) NULL)
    if (is.null(t) || ncol(t) < 2) return(showNotification("Merge table needs two columns", type = "error"))
    t <- data.frame(old_cluster = as.character(t[[1]]), new_cluster = as.character(t[[2]]))
    rv$merge_map <- t
  })
  output$merge_tbl <- renderDT({
    req(rv$merge_map)
    datatable(isolate(rv$merge_map), rownames = FALSE, selection = "none",
              editable = list(target = "cell", disable = list(columns = 0)),
              options = list(pageLength = 50, dom = "t"))
  })
  observeEvent(input$merge_tbl_cell_edit, {
    info <- input$merge_tbl_cell_edit
    rv$merge_map[info$row, info$col + 1] <- as.character(info$value)
  })
  observeEvent(input$apply_merge, {
    req(rv$sce, rv$merge_map, input$merge_k)
    id <- make.names(input$merge_id)
    sce <- safely("Merging", apply_merge(rv$sce, input$merge_k, rv$merge_map, id))
    req(sce); rv$sce <- sce
    log_step(paste0("merge_", id), list(source_k = input$merge_k, table = rv$merge_map))
    refresh_k(selected = id)
    updateSelectInput(session, "dr_color", choices = c(clustering_names(sce), rownames(sce)), selected = id)
    showNotification(paste("Created clustering", id), type = "message")
  })
  output$dl_merge_tbl <- downloadHandler(function() paste0("merge_table_", Sys.Date(), ".csv"),
                                         function(f) utils::write.csv(rv$merge_map, f, row.names = FALSE))
  merged_id <- reactive({ req(rv$sce); id <- make.names(input$merge_id); req(id %in% clustering_names(rv$sce)); id })
  output$merge_heat <- renderPlot(draw_any(plotExprHeatmap(rv$sce, features = "type", by = "cluster_id",
                                                           k = merged_id(), scale = "first", bars = TRUE)))
  output$merge_dr <- renderPlot({ req(has_dr()); plotDR(rv$sce, input$dr_method, color_by = merged_id()) })

  # ---- 5. Batch ----
  observeEvent(input$run_batch, {
    req(rv$sce, input$batch_var)
    ref <- if (has_dr() && input$dr_method == "UMAP") "UMAP" else "none"
    sce <- safely("Batch correction", batch_correct(rv$sce, input$batch_method, input$batch_var,
      npcs = input$batch_npcs, cy_xdim = input$cy_xdim, cy_ydim = input$cy_ydim,
      cy_norm = input$cy_norm, ref_dr = ref, seed = input$seed))
    req(sce); rv$sce <- sce; rv$diag <- NULL
    rv$batch_msg <- sprintf("%s correction on '%s' done. Stored as: %s (raw 'exprs' unchanged).",
                            input$batch_method, input$batch_var, metadata(sce)$batch$corrected_space)
    log_step("batch", list(method = input$batch_method, var = input$batch_var, npcs = input$batch_npcs))
  })
  output$batch_status <- renderText(rv$batch_msg)
  batch_umap <- reactive({
    req(rv$sce, input$batch_var)
    validate(need("UMAP" %in% reducedDimNames(rv$sce), "Run clustering with UMAP embedding first."))
    plot_batch_umaps(rv$sce, input$batch_var)
  })
  output$batch_umap <- renderPlot(batch_umap())
  output$dl_batch_umap <- pdf_handler("batch_umap", batch_umap, 14, 6)
  observeEvent(input$run_diag, {
    req(rv$sce)
    rv$diag <- safely("Diagnostics", batch_diagnostics(rv$sce, metadata(rv$sce)$batch$var %||% input$batch_var,
                                                       label_k = input$diag_k %||% NULL, npcs = input$batch_npcs))
  })
  output$batch_lisi <- renderPlot({ validate(need(rv$diag, "Click 'Compute diagnostics'.")); plot_lisi(rv$diag) })
  output$batch_diag_tbl <- renderDT({ req(rv$diag); datatable(rv$diag$summary, rownames = FALSE) |>
      formatRound(2:5, 3) })

  # ---- 6. DA ----
  observeEvent(input$run_da, {
    req(rv$sce, input$da_k, input$da_cond, input$da_ref, input$da_comp)
    if (input$da_ref == input$da_comp) return(showNotification("Reference and comparison must differ", type = "error"))
    res <- safely("Differential abundance", {
      m <- make_model(rv$sce, input$da_cond, input$da_ref, input$da_comp, input$da_cov, input$da_random)
      run_da(rv$sce, m, input$da_k, input$da_method)
    })
    req(res); rv$da <- res
    log_step("DA", list(k = input$da_k, condition = input$da_cond, ref = input$da_ref, comp = input$da_comp,
                        covariates = input$da_cov, random = input$da_random, method = input$da_method))
  })
  da_heat <- reactive({ req(rv$da); plotDiffHeatmap(rv$sce, rowData(rv$da$res$res), k = rv$da$k,
                                                    fdr = input$fdr_da, all = TRUE, sort_by = "padj") })
  da_props <- reactive({ req(rv$da); plot_paired_props(rv$sce, rv$da$k, rv$da$model) })
  output$da_tbl <- renderDT({ req(rv$da); datatable(rv$da$table, rownames = FALSE, options = list(scrollX = TRUE)) |>
      formatSignif(c("p_val", "p_adj"), 3) |> formatRound(intersect(c("logFC", "log2FC_prop", "prop_ref", "prop_comp"), colnames(rv$da$table)), 4) })
  output$da_heat <- renderPlot(draw_any(da_heat()))
  output$da_volcano <- renderPlot({ req(rv$da); plot_volcano(rv$da$table, "log2FC_prop", input$fdr_da,
    xlab = sprintf("log2 FC proportion (%s / %s)", rv$da$model$comp, rv$da$model$ref)) })
  output$da_forest <- renderPlot({ req(rv$da); plot_da_forest(rv$da$table, input$fdr_da) })
  output$da_props <- renderPlot(da_props())
  output$da_phist <- renderPlot({ req(rv$da); plot_pval_hist(rv$da$table) })
  output$dl_da_heat <- pdf_handler("DA_heatmap", da_heat, 14, 9)
  output$dl_da_props <- pdf_handler("DA_proportions", da_props, 12, 10)
  output$dl_da_csv <- downloadHandler(function() paste0("DA_results_", Sys.Date(), ".csv"),
                                      function(f) utils::write.csv(rv$da$table, f, row.names = FALSE))

  # ---- 7. DS ----
  observeEvent(input$run_ds, {
    req(rv$sce, input$ds_k, input$ds_cond, input$ds_ref, input$ds_comp)
    if (input$ds_ref == input$ds_comp) return(showNotification("Reference and comparison must differ", type = "error"))
    res <- safely("Differential state", {
      m <- make_model(rv$sce, input$ds_cond, input$ds_ref, input$ds_comp, input$ds_cov, input$ds_random)
      run_ds(rv$sce, m, input$ds_k, input$ds_method)
    })
    req(res); rv$ds <- res
    log_step("DS", list(k = input$ds_k, condition = input$ds_cond, ref = input$ds_ref, comp = input$ds_comp,
                        covariates = input$ds_cov, random = input$ds_random, method = input$ds_method))
  })
  ds_heat <- reactive({ req(rv$ds); plotDiffHeatmap(rv$sce, rowData(rv$ds$res$res), k = rv$ds$k,
                                                    fdr = input$fdr_ds, top_n = input$ds_top, sort_by = "padj") })
  output$ds_tbl <- renderDT({ req(rv$ds); datatable(rv$ds$table, rownames = FALSE, filter = "top", options = list(scrollX = TRUE)) |>
      formatSignif(c("p_val", "p_adj"), 3) })
  output$ds_heat <- renderPlot(draw_any(ds_heat()))
  output$ds_volcano <- renderPlot({ req(rv$ds); plot_volcano(rv$ds$table, "delta_median", input$fdr_ds,
    xlab = sprintf("Difference in median arcsinh expression (%s - %s)", rv$ds$model$comp, rv$ds$model$ref)) })
  output$ds_phist <- renderPlot({ req(rv$ds); plot_pval_hist(rv$ds$table) })
  output$dl_ds_heat <- pdf_handler("DS_heatmap", ds_heat, 14, 12)
  output$dl_ds_csv <- downloadHandler(function() paste0("DS_results_", Sys.Date(), ".csv"),
                                      function(f) utils::write.csv(rv$ds$table, f, row.names = FALSE))

  # ---- 8. MEM / ROC ----
  observeEvent(input$run_mem, {
    req(rv$sce, input$mem_k)
    rv$mem <- safely("MEM", compute_mem(rv$sce, input$mem_k))
    req(rv$mem)
    updateSelectInput(session, "roc_cluster", choices = rownames(rv$mem))
  })
  observeEvent(input$mem_k, { req(rv$sce, input$mem_k)
    updateSelectInput(session, "roc_cluster", choices = levels(cluster_ids(rv$sce, input$mem_k))) })
  mem_heat <- reactive(plot_mem(req(rv$mem)))
  output$mem_heat <- renderPlot(draw_any(mem_heat()))
  output$dl_mem <- pdf_handler("MEM_heatmap", mem_heat, 12, 9)
  output$mem_tbl <- renderDT({ req(rv$mem); datatable(data.frame(cluster = rownames(rv$mem), MEM_label = mem_labels(rv$mem, 6)), rownames = FALSE) })
  observeEvent(input$run_roc, {
    req(rv$sce, input$mem_k, input$roc_cluster)
    rv$roc <- safely("ROC", compute_roc(rv$sce, input$mem_k, input$roc_cluster, seed = input$seed))
  })
  output$roc_plot <- renderPlot(plot_roc(req(rv$roc)))
  output$roc_tbl <- renderDT({ req(rv$roc); datatable(rv$roc$table, rownames = FALSE) |> formatRound(c("AUC", "AUC_discriminative"), 3) })

  # ---- 9. Export ----
  output$dl_sce <- downloadHandler(function() paste0("cytof_sce_", Sys.Date(), ".rds"),
    function(f) { req(rv$sce); s <- rv$sce; metadata(s)$analysis_log <- rv$log; saveRDS(s, f) })
  output$dl_cells <- downloadHandler(function() paste0("cells_", Sys.Date(), ".csv.gz"),
    function(f) { req(rv$sce); con <- gzfile(f, "w"); on.exit(close(con)); utils::write.csv(cell_table(rv$sce), con, row.names = FALSE) })
  output$dl_abund <- downloadHandler(function() paste0("abundances_", Sys.Date(), ".csv"),
    function(f) { req(rv$sce, input$k_choose); utils::write.csv(abundance_table(rv$sce, input$k_choose), f, row.names = FALSE) })
  output$log_out <- renderPrint(str(rv$log, max.level = 3, give.attr = FALSE))
  output$dl_log <- downloadHandler(function() paste0("analysis_log_", Sys.Date(), ".yml"),
    function(f) yaml::write_yaml(rv$log, f))
  output$dl_report <- downloadHandler(
    filename = function() paste0("cytof_report_", Sys.Date(), ".html"),
    content = function(file) {
      req(rv$sce)
      withProgress(message = "Rendering report", value = 0.5, {
        tmp <- tempfile("report_"); dir.create(tmp)
        file.copy("inst/report/report.Rmd", file.path(tmp, "report.Rmd"))
        st <- file.path(tmp, "state.rds")
        saveRDS(list(sce = rv$sce, da = rv$da, ds = rv$ds, mem = rv$mem, diag = rv$diag,
                     log = rv$log, k = input$k_choose, dr = input$dr_method,
                     fdr_da = input$fdr_da, fdr_ds = input$fdr_ds), st)
        rmarkdown::render(file.path(tmp, "report.Rmd"), output_file = file,
                          params = list(state = st, rdir = normalizePath("R")),
                          envir = new.env(parent = globalenv()), quiet = TRUE)
      })
    })
  output$session_info <- renderPrint(sessionInfo())
}

shinyApp(ui, server)
