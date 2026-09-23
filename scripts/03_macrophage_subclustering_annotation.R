# ==============================================================================
# MYELOID/MACROPHAGE SUBCLUSTERING
# Continuation of 02_integrated_seurat.rds
# Based on Xu et al. (CCR-22-0296) 
# ==============================================================================
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(batchelor)
  library(BiocSingular)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

set.seed(42)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# To select another location before sourcing the script, run for example:
# Sys.setenv(SCRNA_OUTPUT_DIR = "D:/my_project_outputs")
user_profile <- Sys.getenv("USERPROFILE", unset = "")
if (.Platform$OS.type == "windows" && nzchar(user_profile)) {
  default_output_root <- file.path(user_profile, "scRNAseq_project_output")
} else {
  default_output_root <- getwd()
}

output_root <- Sys.getenv(
  "SCRNA_OUTPUT_DIR",
  unset = default_output_root
)
figures_dir <- file.path(output_root, "figures")

if (!dir.exists(figures_dir)) {
  created <- dir.create(
    figures_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  if (!created && !dir.exists(figures_dir)) {
    stop("Could not create the figure directory: ", figures_dir)
  }
}

figures_dir <- normalizePath(
  figures_dir,
  winslash = "/",
  mustWork = TRUE
)

# Test write permission before starting the time-consuming analysis.
write_test_file <- file.path(figures_dir, ".png_write_test")
write_test_ok <- tryCatch({
  writeLines("write test", write_test_file)
  file.exists(write_test_file)
}, error = function(e) FALSE)

if (file.exists(write_test_file)) unlink(write_test_file)
if (!isTRUE(write_test_ok)) {
  stop(
    "The figure directory is not writable: ", figures_dir,
    ". Set a different directory with ",
    "Sys.setenv(SCRNA_OUTPUT_DIR = 'C:/scRNAseq_project_output')."
  )
}

message("Figures will be saved in: ", figures_dir)

save_png_plot <- function(plot, filename, width, height, dpi = 300) {
  if (!inherits(plot, c("gg", "ggplot", "patchwork"))) {
    warning("The supplied object is not a recognized ggplot/patchwork object.")
  }

  if (!is.numeric(width) || !is.numeric(height) ||
      width <= 0 || height <= 0) {
    stop("width and height must be positive numeric values in inches.")
  }

  # All figure calls are redirected to the validated output directory.
  output_file <- file.path(figures_dir, basename(filename))
  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  device_opened <- FALSE
  device_errors <- character(0)

  # Preferred device on Windows and for high-resolution text rendering.
  if (requireNamespace("ragg", quietly = TRUE)) {
    tryCatch({
      ragg::agg_png(
        filename = output_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        background = "white"
      )
      device_opened <- TRUE
    }, error = function(e) {
      device_errors <<- c(device_errors, paste("ragg:", conditionMessage(e)))
    })
  }

  # Cairo is the first fallback when available.
  if (!device_opened && isTRUE(capabilities("cairo"))) {
    tryCatch({
      grDevices::png(
        filename = output_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        bg = "white",
        type = "cairo-png"
      )
      device_opened <- TRUE
    }, error = function(e) {
      device_errors <<- c(device_errors, paste("cairo-png:", conditionMessage(e)))
    })
  }

  # Final fallback: the platform's standard PNG device.
  if (!device_opened) {
    tryCatch({
      grDevices::png(
        filename = output_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        bg = "white"
      )
      device_opened <- TRUE
    }, error = function(e) {
      device_errors <<- c(device_errors, paste("standard png:", conditionMessage(e)))
    })
  }

  if (!device_opened) {
    stop(
      "No PNG device could be opened for: ", output_file,
      if (length(device_errors) > 0) {
        paste0("\n", paste(device_errors, collapse = "\n"))
      } else {
        ""
      }
    )
  }

  opened_device <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && opened_device %in% open_devices) {
      try(grDevices::dev.off(which = opened_device), silent = TRUE)
    }
  }, add = TRUE)

  print(plot)
  grDevices::dev.off()

  if (!file.exists(output_file) || file.info(output_file)$size == 0) {
    stop("The PNG device closed, but no valid file was created: ", output_file)
  }

  message("Saved figure: ", output_file)
  invisible(output_file)
}

run_mnn_reduction <- function(object, reduction_name = "mnn", k = 5, d = 10) {
  DefaultAssay(object) <- "RNA"
  object <- NormalizeData(object, scale.factor = 10000, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst",
                                 nfeatures = 2000, verbose = FALSE)
  regressors <- intersect(c("nCount_RNA", "percent.mt"), colnames(object[[]]))
  object <- ScaleData(object, features = VariableFeatures(object),
                      vars.to.regress = regressors, verbose = FALSE)

  sce <- as.SingleCellExperiment(object, assay = "RNA")
  mnn_out <- batchelor::fastMNN(
    sce,
    batch = sce$orig.ident,
    subset.row = VariableFeatures(object),
    k = k,
    d = d,
    BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE)
  )
  embedding <- SingleCellExperiment::reducedDim(mnn_out, "corrected")
  if (!all(colnames(object) %in% rownames(embedding))) {
    stop("MNN cell names do not match the Seurat object.")
  }
  embedding <- embedding[colnames(object), , drop = FALSE]
  colnames(embedding) <- paste0("MNN_", seq_len(ncol(embedding)))
  object[[reduction_name]] <- CreateDimReducObject(
    embeddings = embedding, key = "MNN_", assay = "RNA"
  )
  object
}

# ------------------------------------------------------------------------------
# 1. Load script-02 output and subset the broad Monocytic compartment
# ------------------------------------------------------------------------------

input_rds <- "results/02_integrated_seurat.rds"
if (!file.exists(input_rds)) stop("Missing: ", input_rds)
combined_seurat <- readRDS(input_rds)
if (!all(c("cell_type", "orig.ident") %in% colnames(combined_seurat[[]]))) {
  stop("cell_type and/or orig.ident metadata is missing.")
}

broad_myeloid <- subset(combined_seurat, subset = cell_type == "Monocytic")
rm(combined_seurat)
message("Broad Monocytic cells: ", ncol(broad_myeloid))

# ------------------------------------------------------------------------------
# 2. Initial MNN and diagnostic clustering used only for lineage cleaning
# ------------------------------------------------------------------------------

broad_myeloid <- run_mnn_reduction(broad_myeloid, "mnn.initial", k = 5, d = 10)
initial_dims <- 1:min(10, ncol(Embeddings(broad_myeloid, "mnn.initial")))
broad_myeloid <- FindNeighbors(
  broad_myeloid, reduction = "mnn.initial", dims = initial_dims,
  graph.name = c("initial_nn", "initial_snn"), verbose = FALSE
)
broad_myeloid <- FindClusters(
  broad_myeloid, graph.name = "initial_snn", resolution = 0.40,
  algorithm = 1, random.seed = 42,
  cluster.name = "initial_myeloid_cluster", verbose = FALSE
)
broad_myeloid$initial_myeloid_cluster <- as.character(
  broad_myeloid$initial_myeloid_cluster
)
Idents(broad_myeloid) <- "initial_myeloid_cluster"

initial_markers <- FindAllMarkers(
  broad_myeloid, assay = "RNA", test.use = "wilcox", only.pos = TRUE,
  logfc.threshold = 0.25, min.pct = 0.10
) %>% filter(p_val < 0.05)

write.csv(initial_markers,
          "results/03_initial_broad_myeloid_markers.csv", row.names = FALSE)
write.csv(
  initial_markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 20, with_ties = FALSE) %>% ungroup(),
  "results/03_initial_top20_markers.csv", row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Multi-gene lineage audit and contaminant removal
# ------------------------------------------------------------------------------

lineage_panels <- list(
  Myeloid_core = c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "CD68"),
  B_cell = c("PAX5", "MS4A1", "CD79A", "CD19", "BANK1"),
  T_cell = c("CD3D", "CD3E", "TRBC1", "TRBC2"),
  cDC = c("CD1C", "FCER1A", "CD1E", "CLEC10A"),
  pDC = c("CLEC4C", "LILRA4", "GZMB", "DNASE1L3"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "WFDC2"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "COL3A1")
)
lineage_panels <- lapply(lineage_panels, function(x) intersect(x, rownames(broad_myeloid)))

rna_data <- GetAssayData(broad_myeloid, assay = "RNA", layer = "data")
cluster_ids <- sort(unique(broad_myeloid$initial_myeloid_cluster))

panel_detection <- function(cluster_id, genes) {
  cells <- colnames(broad_myeloid)[broad_myeloid$initial_myeloid_cluster == cluster_id]
  if (length(genes) == 0 || length(cells) == 0) return(NA_real_)
  mean(Matrix::rowMeans(rna_data[genes, cells, drop = FALSE] > 0))
}

lineage_audit <- bind_rows(lapply(cluster_ids, function(cluster_id) {
  scores <- vapply(lineage_panels,
                   function(genes) panel_detection(cluster_id, genes), numeric(1))
  data.frame(cluster = cluster_id,
             n_cells = sum(broad_myeloid$initial_myeloid_cluster == cluster_id),
             as.list(scores), check.names = FALSE)
}))

# Cycling cells are retained only if they meet Myeloid_core support.
lineage_audit <- lineage_audit %>%
  mutate(
    dominant_immune_contaminant = pmax(
      B_cell,
      T_cell,
      cDC,
      pDC,
      na.rm = TRUE
    ),
    
    # Epithelial and fibroblast signals are retained as warnings rather than
    # automatic exclusion criteria because ambient RNA may generate low-level
    # keratin or stromal-gene expression in genuine myeloid cells.
    epithelial_warning = Epithelial >= 0.25,
    fibroblast_warning = Fibroblast >= 0.20,
    
    keep_as_myeloid =
      Myeloid_core >= 0.45 &
      Myeloid_core > dominant_immune_contaminant &
      B_cell < 0.25 &
      T_cell < 0.25 &
      cDC < 0.40 &
      pDC < 0.20,
    
    decision = ifelse(
      keep_as_myeloid,
      "KEEP",
      "EXCLUDE"
    )
  )

write.csv(lineage_audit,
          "results/03_initial_cluster_lineage_audit.csv", row.names = FALSE)
print(lineage_audit)

lineage_features <- unique(unlist(lineage_panels))
p_lineage <- DotPlot(
  broad_myeloid, features = lineage_features,
  group.by = "initial_myeloid_cluster", assay = "RNA"
) + RotatedAxis() + ggtitle("Initial-cluster lineage audit") +
  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1))
save_png_plot(p_lineage, "figures/03_initial_lineage_audit_dotplot.png", 15, 7)

clusters_to_keep <- lineage_audit$cluster[lineage_audit$keep_as_myeloid]
if (length(clusters_to_keep) < 2) {
  stop("Fewer than two clusters passed lineage validation. Review the audit CSV.")
}
cells_to_keep <- colnames(broad_myeloid)[
  broad_myeloid$initial_myeloid_cluster %in% clusters_to_keep
]
clean_macrophage <- subset(broad_myeloid, cells = cells_to_keep)
message("Cells retained after lineage cleaning: ", ncol(clean_macrophage))
message("Initial clusters retained: ", paste(clusters_to_keep, collapse = ", "))

# ------------------------------------------------------------------------------
# 4. MNN and clustering on the myeloid/macrophage compartment
# ------------------------------------------------------------------------------

clean_macrophage <- run_mnn_reduction(clean_macrophage, "mnn.clean", k = 5, d = 10)
clean_dims <- 1:min(10, ncol(Embeddings(clean_macrophage, "mnn.clean")))
clean_macrophage <- FindNeighbors(
  clean_macrophage, reduction = "mnn.clean", dims = clean_dims,
  graph.name = c("clean_nn", "clean_snn"), verbose = FALSE
)

# Xu et al. reported graph-based clustering with "optimal parameters". Survey 0.50-2.00 coarsely, then perform
# a fine scan only around the coarse solution closest to ten clusters.
count_clusters_at_resolution <- function(resolution_value) {
  tmp <- FindClusters(
    clean_macrophage,
    graph.name = "clean_snn",
    resolution = resolution_value,
    algorithm = 1,
    random.seed = 42,
    verbose = FALSE
  )
  length(unique(Idents(tmp)))
}

# Fast wide scan followed by a focused fine scan.
coarse_grid <- seq(0.50, 2.00, by = 0.10)
coarse_summary <- data.frame(
  resolution = coarse_grid,
  n_clusters = NA_integer_,
  scan = "coarse"
)

for (i in seq_along(coarse_grid)) {
  coarse_summary$n_clusters[i] <- count_clusters_at_resolution(coarse_grid[i])
  message("Coarse resolution ", coarse_grid[i], ": ",
          coarse_summary$n_clusters[i], " clusters")
}

coarse_summary$distance_from_paper <- abs(coarse_summary$n_clusters - 10L)
best_coarse_i <- which.min(coarse_summary$distance_from_paper)
best_coarse_resolution <- coarse_summary$resolution[best_coarse_i]

fine_grid <- round(seq(
  max(0.01, best_coarse_resolution - 0.10),
  best_coarse_resolution + 0.10,
  by = 0.01
), 2)
fine_summary <- data.frame(
  resolution = fine_grid,
  n_clusters = NA_integer_,
  scan = "fine"
)

for (i in seq_along(fine_grid)) {
  fine_summary$n_clusters[i] <- count_clusters_at_resolution(fine_grid[i])
  message("Fine resolution ", fine_grid[i], ": ",
          fine_summary$n_clusters[i], " clusters")
}

fine_summary$distance_from_paper <- abs(fine_summary$n_clusters - 10L)
resolution_summary <- bind_rows(coarse_summary, fine_summary) %>%
  arrange(resolution, scan) %>%
  distinct(resolution, .keep_all = TRUE)

write.csv(resolution_summary,
          "results/03_clean_macrophage_resolution_survey.csv",
          row.names = FALSE)
print(resolution_summary)

exact_ten_i <- which(fine_summary$n_clusters == 10L)

if (length(exact_ten_i) > 0) {
  selected_i <- exact_ten_i[ceiling(length(exact_ten_i) / 2)]
  final_resolution <- fine_summary$resolution[selected_i]
  selected_n_clusters <- fine_summary$n_clusters[selected_i]
  selection_reason <- paste0(
    "middle fine-grid resolution producing exactly 10 clusters: ",
    final_resolution
  )
} else {
  selected_i <- which.min(fine_summary$distance_from_paper)
  final_resolution <- fine_summary$resolution[selected_i]
  selected_n_clusters <- fine_summary$n_clusters[selected_i]
  selection_reason <- paste0(
    "no tested resolution produced 10 clusters; resolution ",
    final_resolution, " produced the closest result (",
    selected_n_clusters, " clusters)"
  )
  warning(selection_reason)
}

message(
  "Selected clean-macrophage resolution: ", final_resolution,
  " — ", selection_reason
)
writeLines(
  c(
    paste0("selected_resolution=", final_resolution),
    paste0("selected_n_clusters=", selected_n_clusters),
    paste0("selection_reason=", selection_reason)
  ),
  "results/03_clean_macrophage_selected_resolution.txt"
)

p_resolution <- ggplot(
  resolution_summary,
  aes(x = resolution, y = n_clusters)
) +
  geom_line(linewidth = 0.7, color = "#2B83BA") +
  geom_point(size = 2, color = "#2B83BA") +
  geom_hline(yintercept = 10, linetype = 2, color = "#D73027") +
  geom_vline(xintercept = final_resolution, linetype = 3) +
  scale_x_continuous(breaks = seq(0.50, 2.00, by = 0.25)) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title = "Clean macrophage resolution survey",
    subtitle = paste("Selected resolution:", final_resolution),
    x = "Seurat resolution", y = "Number of clusters"
  ) +
  theme_classic()
save_png_plot(
  p_resolution,
  "figures/03_clean_macrophage_resolution_survey.png",
  7, 5
)

clean_macrophage <- FindClusters(
  clean_macrophage, graph.name = "clean_snn", resolution = final_resolution,
  algorithm = 1, random.seed = 42,
  cluster.name = "macrophage_cluster_id", verbose = FALSE
)
clean_macrophage$macrophage_cluster_id <- as.character(
  clean_macrophage$macrophage_cluster_id
)
Idents(clean_macrophage) <- "macrophage_cluster_id"
message(
  "Final clean-macrophage cluster count: ",
  length(unique(clean_macrophage$macrophage_cluster_id))
)

clean_macrophage <- RunUMAP(
  clean_macrophage, reduction = "mnn.clean", dims = clean_dims,
  n.neighbors = 30, min.dist = 0.3, metric = "cosine", seed.use = 42,
  reduction.name = "macrophage.umap", reduction.key = "macUMAP_",
  verbose = FALSE
)

final_markers <- FindAllMarkers(
  clean_macrophage, assay = "RNA", test.use = "wilcox", only.pos = TRUE,
  logfc.threshold = 0.25, min.pct = 0.10
) %>% filter(p_val < 0.05)
write.csv(final_markers, "results/03_clean_macrophage_all_markers.csv",
          row.names = FALSE)
write.csv(
  final_markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 20, with_ties = FALSE) %>% ungroup(),
  "results/03_clean_macrophage_top20_markers.csv", row.names = FALSE
)
# ------------------------------------------------------------------------------
# 5. Xu et al. marker-evidence analysis by numeric cluster
# ------------------------------------------------------------------------------

Idents(clean_macrophage) <- "macrophage_cluster_id"

xu_defining_markers <- c(
  "C0-OLFML3"   = "OLFML3",
  "C1-TCOF1"    = "TCOF1",
  "C2-LPAR6"    = "LPAR6",
  "C3-CCL3L3"   = "CCL3L3",
  "C4-CD36"     = "CD36",
  "C5-TMEM176B" = "TMEM176B",
  "C6-SLC2A1"   = "SLC2A1",
  "C7-APOBEC3A" = "APOBEC3A",
  "C8-CD2"      = "CD2",
  "C9-SIGLEC15" = "SIGLEC15"
)

present_xu_markers <- xu_defining_markers[
  xu_defining_markers %in% rownames(clean_macrophage)
]
missing_xu_markers <- setdiff(
  unname(xu_defining_markers),
  unname(present_xu_markers)
)
if (length(missing_xu_markers) > 0) {
  warning(
    "Xu marker genes absent from the RNA assay: ",
    paste(missing_xu_markers, collapse = ", ")
  )
}

p_xu_cluster_check <- DotPlot(
  clean_macrophage,
  features = present_xu_markers,
  group.by = "macrophage_cluster_id",
  assay = "RNA",
  dot.scale = 7,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("Xu et al. marker evidence by numeric cluster") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    plot.margin = margin(10, 30, 40, 10)
  ) +
  coord_cartesian(clip = "off")
save_png_plot(
  p_xu_cluster_check,
  "figures/03_Xu_marker_evidence_by_numeric_cluster.png",
  12, 7
)

top20_for_annotation <- final_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()
write.csv(
  top20_for_annotation,
  "results/03_top20_markers_for_Xu_annotation.csv",
  row.names = FALSE
)

xu_average_expression <- AggregateExpression(
  clean_macrophage,
  assays = "RNA",
  features = unname(present_xu_markers),
  group.by = "macrophage_cluster_id",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA
write.csv(
  as.data.frame(xu_average_expression),
  "results/03_Xu_marker_average_expression.csv",
  row.names = TRUE
)

xu_marker_zscores <- t(scale(t(as.matrix(xu_average_expression))))
xu_marker_zscores[!is.finite(xu_marker_zscores)] <- 0
write.csv(
  as.data.frame(xu_marker_zscores),
  "results/03_Xu_marker_zscores.csv",
  row.names = TRUE
)

# ==============================================================================
# 5B. VALIDATION OF AMBIGUOUS CLUSTERS
# ==============================================================================

validation_genes <- intersect(
  c(
    "LST1", "CSF1R", "FCER1G", "TYROBP",
    "C1QA", "C1QB", "C1QC", "CD68", "CTSD",
    "CLEC9A", "CD1C", "FCER1A", "CLEC10A", "CST3",
    "CD79A", "CD79B", "MS4A1", "CD37", "HLA-DRA",
    "CDKN1C", "OLFML3", "CX3CR1", "IGF2", "NUPR1"
  ),
  rownames(clean_macrophage)
)

p_validation <- DotPlot(
  clean_macrophage,
  features = validation_genes,
  group.by = "macrophage_cluster_id",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  ggtitle("Macrophage, dendritic and mixed-lineage validation") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    plot.margin = margin(10, 35, 45, 10)
  ) +
  coord_cartesian(clip = "off")

save_png_plot(
  p_validation,
  "figures/03_ambiguous_cluster_validation_dotplot.png",
  16, 6.5
)

mixed_features <- intersect(
  c("CD79B", "LST1", "CSF1R", "FCER1G", "CDKN1C"),
  rownames(clean_macrophage)
)

p_mixed <- FeaturePlot(
  clean_macrophage,
  features = mixed_features,
  reduction = "macrophage.umap",
  ncol = 3,
  order = TRUE
)

save_png_plot(
  p_mixed,
  "figures/03_cluster9_mixed_lineage_featureplots.png",
  14, 9
)

# ------------------------------------------------------------------------------
# 6. Data-driven state proposals using multi-gene programmes
# ------------------------------------------------------------------------------

state_programs <- list(
  C1Q_resident = c("C1QA", "C1QB", "C1QC", "FOLR2", "LYVE1", "SELENOP"),
  FCN1_inflammatory_monocyte = c("FCN1", "S100A8", "S100A9", "S100A12", "VCAN"),
  APOBEC3A_IFN = c("APOBEC3A", "CXCL10", "CXCL11", "IFIT3", "ISG20"),
  SPP1_hypoxic = c("SPP1", "SLC2A1", "HK2", "BNIP3", "ADM"),
  CD36_lipid = c("CD36", "LPL", "APOE", "FABP5", "MSR1"),
  TMEM176 = c("TMEM176A", "TMEM176B"),
  Cycling = c("MKI67", "TOP2A", "UBE2C", "ASPM")
)
state_programs <- lapply(state_programs,
                         function(x) intersect(x, rownames(clean_macrophage)))
for (program in names(state_programs)) {
  if (length(state_programs[[program]]) >= 2) {
    clean_macrophage <- AddModuleScore(
      clean_macrophage, features = list(state_programs[[program]]),
      name = paste0("state_", program), seed = 42
    )
  }
}

score_columns <- grep("^state_.*1$", colnames(clean_macrophage[[]]), value = TRUE)
state_summary <- clean_macrophage[[]] %>%
  select(macrophage_cluster_id, all_of(score_columns)) %>%
  group_by(macrophage_cluster_id) %>%
  summarise(across(all_of(score_columns), mean), .groups = "drop")

score_matrix <- as.matrix(state_summary[, score_columns, drop = FALSE])
rownames(score_matrix) <- state_summary$macrophage_cluster_id
score_z <- scale(score_matrix)
score_z[!is.finite(score_z)] <- 0
dominant_program <- gsub(
  "^state_|1$", "", colnames(score_z)[max.col(score_z)], perl = TRUE
)
state_summary$proposed_state <- paste0(
  "M", state_summary$macrophage_cluster_id, "-", dominant_program
)
write.csv(state_summary, "results/03_clean_macrophage_state_proposals.csv",
          row.names = FALSE)

state_map <- setNames(state_summary$proposed_state,
                      state_summary$macrophage_cluster_id)
clean_macrophage$macrophage_state_proposed <- unname(
  state_map[clean_macrophage$macrophage_cluster_id]
)
Idents(clean_macrophage) <- "macrophage_state_proposed"

# ------------------------------------------------------------------------------
# 7. Final evidence-supported macrophage/myeloid annotation
# ------------------------------------------------------------------------------

# Labels integrate cluster markers, multi-gene programmes, Xu-marker evidence,
# and the macrophage/DC/mixed-lineage validation plots. Xu-like nomenclature is
# retained only for clusters with specific supporting evidence.
final_label_map <- c(
  "0" = "Mac-FOLR2-LYVE1-resident",
  "1" = "Mac-Cycling",
  "2" = "Mac-OLFML3-CX3CR1-resident",
  "3" = "Mac-TNF-CCL4-inflammatory",
  "4" = "C7-APOBEC3A-like",
  "5" = "Mono-FCN1-S100A8-S100A9",
  "6" = "C6-SLC2A1-like",
  "7" = "Mac-NUPR1-IGF2-stress",
  "8" = "Myeloid-CLEC10A-APC",
  "9" = "Mono-CDKN1C"
)

observed_clusters <- sort(unique(clean_macrophage$macrophage_cluster_id))
missing_map_clusters <- setdiff(observed_clusters, names(final_label_map))
if (length(missing_map_clusters) > 0) {
  stop(
    "Clusters missing from final_label_map: ",
    paste(missing_map_clusters, collapse = ", ")
  )
}

clean_macrophage$final_annotation <- unname(
  final_label_map[clean_macrophage$macrophage_cluster_id]
)
if (anyNA(clean_macrophage$final_annotation)) {
  stop("At least one numeric cluster lacks a final annotation.")
}
clean_macrophage$final_macrophage_label <- factor(
  clean_macrophage$final_annotation,
  levels = unname(final_label_map)
)
Idents(clean_macrophage) <- "final_macrophage_label"

# Distinct non-grey palette. The two Xu-like states retain their corresponding
# Figure 4A-inspired red and blue; internally named states use distinct colors.
final_macrophage_colors <- c(
  "Mac-FOLR2-LYVE1-resident" = "#8C510A",
  "Mac-Cycling" = "#D01C8B",
  "Mac-OLFML3-CX3CR1-resident" = "#018571",
  "Mac-TNF-CCL4-inflammatory" = "#E66101",
  "C7-APOBEC3A-like" = "#E56E6A",
  "Mono-FCN1-S100A8-S100A9" = "#FDB863",
  "C6-SLC2A1-like" = "#2D8DBD",
  "Mac-NUPR1-IGF2-stress" = "#5E3C99",
  "Myeloid-CLEC10A-APC" = "#1B7837",
  "Mono-CDKN1C" = "#35978F"
)
missing_color_labels <- setdiff(
  levels(clean_macrophage$final_macrophage_label),
  names(final_macrophage_colors)
)
if (length(missing_color_labels) > 0) {
  stop(
    "No plotting color defined for: ",
    paste(missing_color_labels, collapse = ", ")
  )
}

annotation_audit <- data.frame(
  macrophage_cluster_id = observed_clusters,
  final_label = vapply(
    observed_clusters,
    function(cluster_id) {
      as.character(unique(clean_macrophage$final_macrophage_label[
        clean_macrophage$macrophage_cluster_id == cluster_id
      ])[1])
    },
    character(1)
  ),
  label_status = ifelse(
    observed_clusters %in% c("4", "6"),
    "Xu-like; marker-supported",
    "Internal evidence-supported"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  annotation_audit,
  "results/03_final_macrophage_annotation_audit.csv",
  row.names = FALSE
)
print(table(
  clean_macrophage$macrophage_cluster_id,
  clean_macrophage$final_macrophage_label
))

# ------------------------------------------------------------------------------
# ==============================================================================
# 7B. INVESTIGATION OF UNRESOLVED MACROPHAGE CLUSTERS
# ==============================================================================

# Identify numeric clusters currently carrying an unresolved annotation.
unresolved_clusters <- clean_macrophage[[]] %>%
  filter(grepl("^Unresolved", final_macrophage_label)) %>%
  pull(macrophage_cluster_id) %>%
  as.character() %>%
  unique() %>%
  sort()

message(
  "Unresolved numeric clusters: ",
  paste(unresolved_clusters, collapse = ", ")
)

if (length(unresolved_clusters) > 0) {

# ------------------------------------------------------------------------------
# 7B.1 Positive markers for each unresolved cluster
# ------------------------------------------------------------------------------

Idents(clean_macrophage) <- "macrophage_cluster_id"

unresolved_marker_list <- vector(
  "list",
  length(unresolved_clusters)
)
names(unresolved_marker_list) <- unresolved_clusters

for (cluster_id in unresolved_clusters) {
  
  message("Calculating markers for unresolved cluster ", cluster_id)
  
  marker_result <- FindMarkers(
    clean_macrophage,
    ident.1 = cluster_id,
    assay = "RNA",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  if (nrow(marker_result) > 0) {
    marker_result$gene <- rownames(marker_result)
    marker_result$macrophage_cluster_id <- cluster_id
    rownames(marker_result) <- NULL
    unresolved_marker_list[[cluster_id]] <- marker_result
  }
}

unresolved_markers <- bind_rows(unresolved_marker_list)

write.csv(
  unresolved_markers,
  "results/03_unresolved_cluster_markers.csv",
  row.names = FALSE
)

fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(unresolved_markers)
)[1]

if (is.na(fc_column)) {
  stop("No log-fold-change column was found in the unresolved marker table.")
}

unresolved_top20 <- unresolved_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(macrophage_cluster_id) %>%
  slice_max(
    order_by = .data[[fc_column]],
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  unresolved_top20,
  "results/03_unresolved_top20_markers.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7B.2 Multi-gene biological-state evidence
# ------------------------------------------------------------------------------

resolution_signatures <- list(
  
  Myeloid_core = c(
    "LST1", "TYROBP", "FCER1G", "AIF1",
    "CTSS", "CD68", "CSF1R"
  ),
  
  C1Q_resident = c(
    "C1QA", "C1QB", "C1QC", "APOE",
    "TREM2", "GPNMB", "MARCO", "MSR1"
  ),
  
  Lipid_associated = c(
    "CD36", "LPL", "APOC1", "FABP5",
    "PLIN2", "LGALS3", "MSR1"
  ),
  
  FCN1_inflammatory_monocyte = c(
    "FCN1", "VCAN", "S100A8", "S100A9",
    "CCR2", "SELL", "CTSS", "IL1B"
  ),
  
  CX3CR1_FCGR3A = c(
    "CX3CR1", "FCGR3A", "LILRB1",
    "IFITM3", "SAT1", "CDKN1C"
  ),
  
  CCL3_inflammatory = c(
    "CCL3", "CCL3L1", "CCL3L3", "CCL4",
    "TNF", "NFKBIA", "IL1B", "IER3"
  ),
  
  APOBEC3A_IFN = c(
    "APOBEC3A", "ISG15", "IFI6", "IFIT1",
    "IFIT3", "MX1", "STAT1", "CXCL10"
  ),
  
  Hypoxic_glycolytic = c(
    "SLC2A1", "VEGFA", "LDHA", "ENO1",
    "PGK1", "BNIP3", "HILPDA"
  ),
  
  Stress_response = c(
    "FOS", "JUN", "JUNB", "ATF3",
    "HSPA1A", "HSPA1B", "DDIT3"
  ),
  
  Cycling = c(
    "MKI67", "TOP2A", "TYMS",
    "UBE2C", "CENPF", "STMN1"
  )
)

resolution_signatures <- lapply(
  resolution_signatures,
  intersect,
  y = rownames(clean_macrophage)
)

resolution_signatures <- resolution_signatures[
  lengths(resolution_signatures) >= 3
]

# Remove duplicated genes because DotPlot cannot use duplicated feature levels.
resolution_features <- unique(
  unlist(resolution_signatures, use.names = FALSE)
)

unresolved_cells_for_plot <- colnames(clean_macrophage)[
  clean_macrophage$macrophage_cluster_id %in% unresolved_clusters
]

unresolved_object <- subset(
  clean_macrophage,
  cells = unresolved_cells_for_plot
)

p_unresolved_signatures <- DotPlot(
  unresolved_object,
  features = resolution_features,
  group.by = "macrophage_cluster_id",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  ggtitle("Evidence for unresolved macrophage clusters") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    plot.margin = margin(10, 35, 45, 10)
  ) +
  coord_cartesian(clip = "off")

save_png_plot(
  p_unresolved_signatures,
  "figures/03_unresolved_cluster_signature_dotplot.png",
  16,
  7
)

# ------------------------------------------------------------------------------
# 7B.3 Sample representation of unresolved clusters
# ------------------------------------------------------------------------------

unresolved_sample_distribution <- clean_macrophage[[]] %>%
  filter(macrophage_cluster_id %in% unresolved_clusters) %>%
  count(
    macrophage_cluster_id,
    orig.ident,
    name = "n_cells"
  ) %>%
  group_by(macrophage_cluster_id) %>%
  mutate(
    cluster_total = sum(n_cells),
    sample_fraction = n_cells / cluster_total
  ) %>%
  ungroup()

write.csv(
  unresolved_sample_distribution,
  "results/03_unresolved_sample_distribution.csv",
  row.names = FALSE
)

# Identify clusters excessively dominated by one sample.
unresolved_sample_summary <- unresolved_sample_distribution %>%
  group_by(macrophage_cluster_id) %>%
  summarise(
    n_samples = sum(n_cells > 0),
    dominant_sample = orig.ident[which.max(sample_fraction)],
    maximum_sample_fraction = max(sample_fraction),
    sample_dominance_warning =
      maximum_sample_fraction >= 0.80,
    .groups = "drop"
  )

write.csv(
  unresolved_sample_summary,
  "results/03_unresolved_sample_summary.csv",
  row.names = FALSE
)

print(unresolved_sample_summary)

} else {
  message(
    "No unresolved clusters remain after the final evidence-supported map; ",
    "Section 7B was skipped."
  )
}

# Return identity to the final macrophage annotation before Section 8.
Idents(clean_macrophage) <- "final_macrophage_label"
# ------------------------------------------------------------------------------

# 8. Cleaned UMAP and annotation-evidence plots
# ------------------------------------------------------------------------------

# Obtain UMAP coordinates and calculate the centre of each numeric cluster.
umap_coordinates <- as.data.frame(
  Embeddings(clean_macrophage, "macrophage.umap")
)

umap_coordinates$macrophage_cluster_id <- as.character(
  clean_macrophage$macrophage_cluster_id
)

cluster_centres <- umap_coordinates %>%
  group_by(macrophage_cluster_id) %>%
  summarise(
    UMAP_1 = median(.data[[colnames(umap_coordinates)[1]]]),
    UMAP_2 = median(.data[[colnames(umap_coordinates)[2]]]),
    .groups = "drop"
  )

# Cells are coloured by the full biological annotation.
# Only numeric cluster IDs are displayed inside the UMAP.
p_umap <- DimPlot(
  clean_macrophage,
  reduction = "macrophage.umap",
  group.by = "final_macrophage_label",
  cols = final_macrophage_colors,
  label = FALSE,
  pt.size = 0.45
) +
  ggrepel::geom_label_repel(
    data = cluster_centres,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = macrophage_cluster_id
    ),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold",
    color = "black",
    fill = scales::alpha("white", 0.78),
    label.size = 0.20,
    box.padding = 0.35,
    point.padding = 0.20,
    min.segment.length = 0,
    seed = 42
  ) +
  ggtitle("Evidence-supported macrophage subclusters") +
  labs(color = "Subcluster") +
  theme_classic() +
  theme(
    plot.margin = margin(10, 30, 10, 10),
    legend.text = element_text(size = 8)
  )

save_png_plot(
  p_umap,
  "figures/03_clean_macrophage_subcluster_umap_numbered.png",
  9,
  7
)
# ------------------------------------------------------------------------------
# Presentation-ready UMAP with short biological labels
# ------------------------------------------------------------------------------

short_label_map <- c(
  "0" = "Resident",
  "1" = "Cycling",
  "2" = "OLFML3",
  "3" = "Inflammatory",
  "4" = "APOBEC3A",
  "5" = "FCN1",
  "6" = "SLC2A1",
  "7" = "Stress",
  "8" = "CLEC10A-APC",
  "9" = "CDKN1C"
)

# Full biological descriptions are used only as legend labels. The shorter
# labels above remain inside the UMAP to avoid obscuring the cell distribution.
full_biological_label_map <- c(
  "Mac-FOLR2-LYVE1-resident" =
    "FOLR2+ LYVE1+ resident macrophages",
  "Mac-Cycling" =
    "MKI67+ proliferating macrophages",
  "Mac-OLFML3-CX3CR1-resident" =
    "OLFML3+ CX3CR1+ resident/APC macrophages",
  "Mac-TNF-CCL4-inflammatory" =
    "TNF+ CCL4+ inflammatory macrophages",
  "C7-APOBEC3A-like" =
    "APOBEC3A+ IFN-responsive macrophages",
  "Mono-FCN1-S100A8-S100A9" =
    "FCN1+ S100A8/A9+ classical monocytes",
  "C6-SLC2A1-like" =
    "SLC2A1+ hypoxic-glycolytic TAMs",
  "Mac-NUPR1-IGF2-stress" =
    "NUPR1+ IGF2+ stress-adapted macrophages",
  "Myeloid-CLEC10A-APC" =
    "CLEC10A+ dendritic-like APCs",
  "Mono-CDKN1C" =
    "CDKN1C+ non-classical/intermediate monocytes"
)

missing_legend_labels <- setdiff(
  levels(clean_macrophage$final_macrophage_label),
  names(full_biological_label_map)
)
if (length(missing_legend_labels) > 0) {
  stop(
    "Full biological legend labels are missing for: ",
    paste(missing_legend_labels, collapse = ", ")
  )
}

cluster_centres$short_label <- unname(
  short_label_map[cluster_centres$macrophage_cluster_id]
)

if (anyNA(cluster_centres$short_label)) {
  stop("A short biological label is missing for at least one cluster.")
}

p_umap_biological <- DimPlot(
  clean_macrophage,
  reduction = "macrophage.umap",
  group.by = "final_macrophage_label",
  label = FALSE,
  pt.size = 0.45
) +
  scale_color_manual(
    values = final_macrophage_colors,
    breaks = levels(clean_macrophage$final_macrophage_label),
    labels = unname(
      full_biological_label_map[
        levels(clean_macrophage$final_macrophage_label)
      ]
    ),
    drop = FALSE
  ) +
  ggrepel::geom_label_repel(
    data = cluster_centres,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = short_label
    ),
    inherit.aes = FALSE,
    size = 3.3,
    fontface = "bold",
    color = "black",
    fill = scales::alpha("white", 0.82),
    label.size = 0.20,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 42,
    max.overlaps = Inf
  ) +
  ggtitle("Evidence-based myeloid and macrophage subclusters") +
  labs(color = "Myeloid/macrophage state") +
  theme_classic() +
  theme(
    plot.margin = margin(10, 55, 10, 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

save_png_plot(
  p_umap_biological,
  "figures/03_clean_macrophage_subcluster_umap_biological_labels.png",
  12,
  7.5
)


paper_marker_genes <- c(
  "OLFML3", "TCOF1", "LPAR6", "CCL3L3", "CD36",
  "TMEM176B", "SLC2A1", "APOBEC3A", "CD2", "SIGLEC15"
)
paper_marker_genes <- intersect(paper_marker_genes, rownames(clean_macrophage))
p_paper_markers <- DotPlot(
  clean_macrophage, features = paper_marker_genes,
  group.by = "final_macrophage_label", assay = "RNA"
) + RotatedAxis() +
  ggtitle("Xu et al. label-marker evidence; labels not forced") + theme_classic()
save_png_plot(p_paper_markers,
              "figures/03_clean_Xu_marker_evidence_dotplot.png", 12, 7)

# Figure 4E genes from the paper.
fig4e_genes <- c(
  "CCL8", "CXCL9", "CXCL10", "CXCL11", "CCL20", "CXCL2", "CXCL3",
  "IL10", "IL16", "IL1B", "TNF", "TNFSF10", "EREG", "IGF1", "HGF",
  "PDGFB", "VEGFA"
)
fig4e_genes <- intersect(fig4e_genes, rownames(clean_macrophage))
p_4e <- DotPlot(
  clean_macrophage, features = fig4e_genes,
  group.by = "final_macrophage_label", assay = "RNA",
  dot.scale = 6, scale = TRUE
) + RotatedAxis() +
  scale_color_gradient(low = "#F1EEEE", high = "#F03B2C") +
  ggtitle("Figure 4E-style signaling-gene expression") + theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
save_png_plot(p_4e, "figures/03_clean_fig4e_signaling_dotplot.png", 14, 7)

# ------------------------------------------------------------------------------
# 9. Clinical groups from Xu et al. Supplementary Table S1 and Figure 4D
# ------------------------------------------------------------------------------

stage_map <- c(
  "norm1"="Normal", "norm2"="Normal", "norm3"="Normal",
  "norm4"="Normal", "norm5"="Normal",
  "cancer1"="Advanced", "cancer2"="Advanced",
  "cancer3"="Early", "cancer4"="Early",
  "cancer5"="Advanced", "cancer6"="Advanced", "cancer7"="Early"
)
figo_map <- c(
  "norm1"="Normal", "norm2"="Normal", "norm3"="Normal",
  "norm4"="Normal", "norm5"="Normal",
  "cancer1"="IIIB", "cancer2"="IIB", "cancer3"="IC2",
  "cancer4"="IC2", "cancer5"="IIB", "cancer6"="IIIC", "cancer7"="IC2"
)
extract_sample_id <- function(x) {
  y <- tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
  ok <- grepl("norm[1-5]|cancer[1-7]", y)
  out <- rep(NA_character_, length(y))
  out[ok] <- sub(".*(norm[1-5]|cancer[1-7]).*", "\\1", y[ok])
  out
}
clean_macrophage$paper_sample_id <- extract_sample_id(clean_macrophage$orig.ident)
clean_macrophage$clinical_stage_group <- factor(
  unname(stage_map[clean_macrophage$paper_sample_id]),
  levels = c("Normal", "Early", "Advanced")
)
clean_macrophage$figo_stage <- unname(figo_map[clean_macrophage$paper_sample_id])
if (anyNA(clean_macrophage$clinical_stage_group)) {
  stop("Some orig.ident values do not match Supplementary Table S1.")
}

stage_audit <- unique(data.frame(
  orig.ident = clean_macrophage$orig.ident,
  paper_sample_id = clean_macrophage$paper_sample_id,
  clinical_stage_group = clean_macrophage$clinical_stage_group,
  figo_stage = clean_macrophage$figo_stage
))
write.csv(stage_audit, "results/03_clean_sample_stage_audit.csv", row.names = FALSE)

proportion_df <- clean_macrophage[[]] %>%
  count(final_macrophage_label, clinical_stage_group, name = "n") %>%
  group_by(final_macrophage_label) %>%
  mutate(fraction = n / sum(n)) %>% ungroup()

p_4d <- ggplot(
  proportion_df,
  aes(x = fraction, y = final_macrophage_label,
      fill = clinical_stage_group)
) + geom_col(width = 0.72) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2),
                     position = "top", expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Normal"="#2B83BA", "Early"="#F07C34", "Advanced"="#A6A6A6"
  ), drop = FALSE) +
  labs(title = "Macrophage states by disease stage", x = NULL, y = NULL,
       fill = NULL) + theme_classic() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom")
save_png_plot(p_4d, "figures/03_clean_fig4d_stage_fractions.png", 9, 6)
write.csv(proportion_df, "results/03_clean_stage_fractions.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 10. Save the cleaned and evidence-annotated object
# ------------------------------------------------------------------------------

saveRDS(clean_macrophage, "results/03_clean_macrophages_subclustered.rds")
saveRDS(clean_macrophage, "results/03_clean_macrophages_Xu_annotated.rds")
writeLines(capture.output(sessionInfo()),
           "results/03_clean_macrophage_session_info.txt")

message("Cleaned macrophage analysis complete.")
message("Review lineage audit and state-proposal CSVs before final naming.")
