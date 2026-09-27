# ==============================================================================
# 2.1 REFINE GLOBAL MYELOID ANNOTATION AND REANNOTATE B-LIKE CELLS
# ==============================================================================
# Loads results/02_integrated_seurat.rds, preserves the Xu annotation, detects
# B-dominant clusters inside the original Monocytic compartment, compares them
# with the global B cell/plasma population, and creates a refined annotation.
# The original RDS and `cell_type` metadata are never overwritten.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

set.seed(42)
INPUT_RDS <- "results/02_integrated_seurat.rds"
OUTPUT_RDS <- "results/02.1_integrated_seurat_refined.rds"
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

safe_ggsave <- function(filename, plot, width, height, dpi = 300) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width, height = height, units = "in",
                  res = dpi, background = "white")
    print(plot)
    grDevices::dev.off()
  } else {
    ggsave(filename, plot, width = width, height = height, dpi = dpi,
           bg = "white", device = "png")
  }
}

# 1. Load and validate script-02 output -----------------------------------------
if (!file.exists(INPUT_RDS)) stop("Input object was not found: ", INPUT_RDS)
message("Loading: ", INPUT_RDS)
integrated_seurat <- readRDS(INPUT_RDS)

required_metadata <- c("cell_type", "orig.ident")
missing_metadata <- setdiff(required_metadata, colnames(integrated_seurat[[]]))
if (length(missing_metadata) > 0) {
  stop("Required metadata columns are missing: ",
       paste(missing_metadata, collapse = ", "))
}

DefaultAssay(integrated_seurat) <- "RNA"
original_cell_type <- as.character(integrated_seurat$cell_type)
monocytic_labels <- unique(original_cell_type[
  grepl("^monocytic$|^monocyte", original_cell_type, ignore.case = TRUE)
])
if (length(monocytic_labels) == 0) {
  stop("No Monocytic annotation was found in `cell_type`. Available labels: ",
       paste(sort(unique(original_cell_type)), collapse = ", "))
}

monocytic_cells <- colnames(integrated_seurat)[
  original_cell_type %in% monocytic_labels
]
message("Original Monocytic cells: ", length(monocytic_cells))
if (length(monocytic_cells) < 100) stop("Too few Monocytic cells for review.")
broad_myeloid <- subset(integrated_seurat, cells = monocytic_cells)

# 2. Recluster the Monocytic compartment for lineage auditing -----------------
DefaultAssay(broad_myeloid) <- "RNA"
broad_myeloid <- NormalizeData(broad_myeloid, scale.factor = 10000,
                                verbose = FALSE)
broad_myeloid <- FindVariableFeatures(broad_myeloid, selection.method = "vst",
                                      nfeatures = 2000, verbose = FALSE)
regressors <- intersect(c("nCount_RNA", "percent.mt"),
                        colnames(broad_myeloid[[]]))
broad_myeloid <- ScaleData(broad_myeloid,
                           features = VariableFeatures(broad_myeloid),
                           vars.to.regress = regressors, verbose = FALSE)

# Use fastMNN when available; otherwise use PCA. B-like clusters are selected
# from lineage evidence, so no numeric cluster ID is hard-coded.
use_mnn <- all(vapply(
  c("batchelor", "SingleCellExperiment", "BiocSingular"),
  requireNamespace, logical(1), quietly = TRUE
))

if (use_mnn) {
  message("Using fastMNN for the lineage-audit embedding.")
  sce <- as.SingleCellExperiment(broad_myeloid, assay = "RNA")
  mnn_out <- batchelor::fastMNN(
    sce, batch = sce$orig.ident,
    subset.row = VariableFeatures(broad_myeloid), k = 5, d = 10,
    BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE)
  )
  audit_embedding <- SingleCellExperiment::reducedDim(mnn_out, "corrected")
  if (!all(colnames(broad_myeloid) %in% rownames(audit_embedding))) {
    stop("fastMNN output cell names do not match the Seurat object.")
  }
  audit_embedding <- audit_embedding[colnames(broad_myeloid), , drop = FALSE]
  colnames(audit_embedding) <- paste0("MNNaudit_",
                                      seq_len(ncol(audit_embedding)))
  broad_myeloid[["mnn.audit"]] <- CreateDimReducObject(
    embeddings = audit_embedding, key = "MNNaudit_", assay = "RNA"
  )
  audit_reduction <- "mnn.audit"
  audit_dims <- seq_len(min(10, ncol(audit_embedding)))
} else {
  warning("fastMNN packages unavailable; using PCA for the lineage audit. ",
          "Cluster numbers may differ from script 03, but marker-based ",
          "B-like cluster selection remains active.")
  broad_myeloid <- RunPCA(broad_myeloid,
                          features = VariableFeatures(broad_myeloid),
                          npcs = 20, seed.use = 42, verbose = FALSE)
  audit_reduction <- "pca"
  audit_dims <- 1:10
}

broad_myeloid <- FindNeighbors(
  broad_myeloid, reduction = audit_reduction, dims = audit_dims,
  graph.name = c("audit_nn", "audit_snn"), verbose = FALSE
)
broad_myeloid <- FindClusters(
  broad_myeloid, graph.name = "audit_snn", resolution = 0.40,
  algorithm = 1, random.seed = 42, cluster.name = "audit_cluster",
  verbose = FALSE
)
broad_myeloid$audit_cluster <- as.character(broad_myeloid$audit_cluster)

# 3. Multi-gene cluster-level lineage audit -----------------------------------
lineage_panels <- list(
  Myeloid_core = c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "CD68"),
  B_cell = c("PAX5", "MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  Plasma_cell = c("JCHAIN", "MZB1", "SDC1", "XBP1", "DERL3"),
  T_cell = c("CD3D", "CD3E", "TRBC1", "TRBC2"),
  cDC = c("CD1C", "FCER1A", "CD1E", "CLEC10A"),
  pDC = c("CLEC4C", "LILRA4", "GZMB", "DNASE1L3")
)
lineage_panels <- lapply(lineage_panels,
                         function(x) intersect(x, rownames(broad_myeloid)))
empty_panels <- names(lineage_panels)[lengths(lineage_panels) == 0]
if (length(empty_panels) > 0) {
  stop("No detected genes for panels: ", paste(empty_panels, collapse = ", "))
}

rna_data <- GetAssayData(broad_myeloid, assay = "RNA", layer = "data")
audit_cluster_ids <- sort(unique(broad_myeloid$audit_cluster))
panel_detection <- function(cluster_id, genes) {
  cells <- colnames(broad_myeloid)[broad_myeloid$audit_cluster == cluster_id]
  if (length(cells) == 0 || length(genes) == 0) return(NA_real_)
  mean(Matrix::rowMeans(rna_data[genes, cells, drop = FALSE] > 0))
}

lineage_audit <- bind_rows(lapply(audit_cluster_ids, function(cluster_id) {
  scores <- vapply(lineage_panels,
                   function(g) panel_detection(cluster_id, g), numeric(1))
  data.frame(audit_cluster = cluster_id,
             n_cells = sum(broad_myeloid$audit_cluster == cluster_id),
             as.list(scores), check.names = FALSE)
})) %>%
  mutate(
    strongest_non_B = pmax(Myeloid_core, T_cell, cDC, pDC, na.rm = TRUE),
    B_like_cluster = B_cell >= 0.40 & B_cell > Myeloid_core &
      B_cell > strongest_non_B,
    proposed_B_subtype = case_when(
      B_like_cluster & Plasma_cell >= 0.30 & Plasma_cell > B_cell ~
        "Plasma-like",
      B_like_cluster ~ "B-cell-like",
      TRUE ~ "Not B-like"
    ),
    reannotation_decision = ifelse(
      B_like_cluster, "REASSIGN_TO_B_CELL_PLASMA", "KEEP_AS_MYELOID"
    )
  )

write.csv(lineage_audit, "results/02.1_monocytic_lineage_audit.csv",
          row.names = FALSE)
print(lineage_audit)

b_like_clusters <- lineage_audit$audit_cluster[lineage_audit$B_like_cluster]
if (length(b_like_clusters) == 0) {
  warning("No B-dominant clusters passed the conservative criteria. ",
          "No cells will be reassigned automatically.")
  b_like_cells <- character(0)
} else {
  b_like_cells <- colnames(broad_myeloid)[
    broad_myeloid$audit_cluster %in% b_like_clusters
  ]
}
message("B-like audit clusters: ",
        ifelse(length(b_like_clusters) == 0, "none",
               paste(b_like_clusters, collapse = ", ")))
message("B-like cells proposed for reassignment: ", length(b_like_cells))

# Lineage-audit dot plot.
lineage_features <- unique(unlist(lineage_panels))
p_lineage_audit <- DotPlot(
  broad_myeloid, features = lineage_features,
  group.by = "audit_cluster", assay = "RNA", dot.scale = 6
) + RotatedAxis() +
  ggtitle("Lineage audit of the original Monocytic compartment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7))
safe_ggsave("figures/02.1_monocytic_lineage_audit_dotplot.png",
            p_lineage_audit, 15, 7)

# 4. Compare B-like cells with the global B cell/plasma compartment -----------
global_b_plasma_cells <- colnames(integrated_seurat)[
  grepl("B.*plasma|B[ _-]*cell|plasma", original_cell_type,
        ignore.case = TRUE) &
    !(colnames(integrated_seurat) %in% b_like_cells)
]

if (length(b_like_cells) > 0 && length(global_b_plasma_cells) > 0) {
  comparison_object <- subset(
    integrated_seurat,
    cells = unique(c(b_like_cells, global_b_plasma_cells))
  )
  DefaultAssay(comparison_object) <- "RNA"
  comparison_object <- NormalizeData(comparison_object, verbose = FALSE)
  comparison_object$B_lineage_comparison <- ifelse(
    colnames(comparison_object) %in% b_like_cells,
    "B-like cells from Monocytic compartment",
    "Original global B cell/plasma compartment"
  )

  comparison_features <- unique(unlist(
    lineage_panels[c("B_cell", "Plasma_cell", "Myeloid_core")]
  ))
  comparison_features <- intersect(comparison_features,
                                   rownames(comparison_object))
  p_b_comparison <- DotPlot(
    comparison_object, features = comparison_features,
    group.by = "B_lineage_comparison", assay = "RNA",
    dot.scale = 7, scale = TRUE
  ) + RotatedAxis() +
    ggtitle("B-like cells versus the global B cell/plasma compartment") +
    labs(x = NULL, y = NULL) + theme_classic() +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
          plot.margin = margin(10, 45, 45, 10))
  safe_ggsave("figures/02.1_B_like_vs_global_B_plasma_dotplot.png",
              p_b_comparison, 13, 5.5)

  Idents(comparison_object) <- "B_lineage_comparison"
  comparison_markers <- FindMarkers(
    comparison_object,
    ident.1 = "B-like cells from Monocytic compartment",
    ident.2 = "Original global B cell/plasma compartment",
    assay = "RNA", test.use = "wilcox", logfc.threshold = 0,
    min.pct = 0.05
  )
  comparison_markers$gene <- rownames(comparison_markers)
  write.csv(comparison_markers,
            "results/02.1_B_like_vs_global_B_plasma_markers.csv",
            row.names = FALSE)
} else {
  warning("B-like/global B-cell comparison skipped because a group was empty.")
}

# 5. Create refined annotation while preserving the original ------------------
integrated_seurat$cell_type_Xu <- original_cell_type
refined_cell_type <- original_cell_type
refined_cell_type[colnames(integrated_seurat) %in% monocytic_cells] <- "Myeloid"
if (length(b_like_cells) > 0) {
  refined_cell_type[colnames(integrated_seurat) %in% b_like_cells] <-
    "B cell_plasma"
}
integrated_seurat$cell_type_refined <- factor(refined_cell_type)
integrated_seurat$B_like_reannotation <- ifelse(
  colnames(integrated_seurat) %in% b_like_cells,
  "Reassigned from Monocytic to B cell_plasma", "Not reassigned"
)

integrated_seurat$myeloid_audit_cluster <- NA_character_
audit_positions <- match(colnames(broad_myeloid), colnames(integrated_seurat))
integrated_seurat$myeloid_audit_cluster[audit_positions] <-
  broad_myeloid$audit_cluster

annotation_transition <- as.data.frame(table(
  Xu_annotation = integrated_seurat$cell_type_Xu,
  refined_annotation = integrated_seurat$cell_type_refined
)) %>% filter(Freq > 0)
write.csv(annotation_transition,
          "results/02.1_annotation_transition_table.csv", row.names = FALSE)

b_like_cell_audit <- data.frame(
  cell = b_like_cells,
  orig.ident = integrated_seurat$orig.ident[
    match(b_like_cells, colnames(integrated_seurat))
  ],
  original_annotation = original_cell_type[
    match(b_like_cells, colnames(integrated_seurat))
  ],
  refined_annotation = refined_cell_type[
    match(b_like_cells, colnames(integrated_seurat))
  ],
  audit_cluster = broad_myeloid$audit_cluster[
    match(b_like_cells, colnames(broad_myeloid))
  ],
  stringsAsFactors = FALSE
)
write.csv(b_like_cell_audit, "results/02.1_reannotated_B_like_cells.csv",
          row.names = FALSE)

# 6. Refined global UMAP -------------------------------------------------------
available_reductions <- Reductions(integrated_seurat)
if ("umap" %in% available_reductions) {
  global_umap <- "umap"
} else {
  umap_candidates <- grep("umap", available_reductions, value = TRUE,
                          ignore.case = TRUE)
  global_umap <- if (length(umap_candidates) == 0) NA_character_ else
    umap_candidates[1]
}

if (!is.na(global_umap)) {
  p_refined_global <- DimPlot(
    integrated_seurat, reduction = global_umap,
    group.by = "cell_type_refined", label = TRUE, repel = TRUE,
    pt.size = 0.25
  ) + ggtitle("Refined global cell-type annotation") +
    labs(color = "Refined cell type") + theme_classic() +
    theme(plot.margin = margin(10, 35, 10, 10))
  safe_ggsave("figures/02.1_refined_global_cell_type_umap.png",
              p_refined_global, 10, 7.5)

  if (length(b_like_cells) > 0) {
    p_b_highlight <- DimPlot(
      integrated_seurat, reduction = global_umap,
      cells.highlight = b_like_cells, cols = "grey85",
      cols.highlight = "#D01C8B", sizes.highlight = 0.6, pt.size = 0.20
    ) + ggtitle("B-like cells reassigned from the Monocytic compartment") +
      theme_classic()
    safe_ggsave("figures/02.1_reannotated_B_like_cells_global_umap.png",
                p_b_highlight, 9, 7)
  }
} else {
  warning("No UMAP reduction found; UMAP outputs were skipped.")
}

# 7. Save the refined object and audit outputs --------------------------------
saveRDS(integrated_seurat, OUTPUT_RDS)
writeLines(c(
  paste0("input_rds=", INPUT_RDS),
  paste0("output_rds=", OUTPUT_RDS),
  paste0("original_monocytic_cells=", length(monocytic_cells)),
  paste0("B_like_clusters=", paste(b_like_clusters, collapse = ";")),
  paste0("B_like_cells_reassigned=", length(b_like_cells)),
  paste0("remaining_cells_labelled_myeloid=",
         sum(integrated_seurat$cell_type_refined == "Myeloid")),
  paste0("audit_reduction=", audit_reduction)
), "results/02.1_refined_annotation_summary.txt")
writeLines(capture.output(sessionInfo()), "results/02.1_session_info.txt")

message("Refined global annotation complete.")
message("Original annotation retained in: cell_type_Xu")
message("Refined annotation saved in: cell_type_refined")
message("Refined object saved to: ", OUTPUT_RDS)
