
# ==============================================================================
# 02. NORMALIZATION, CCA INTEGRATION, CLUSTERING, AND GLOBAL ANNOTATION
#
# Input:  results/postQC_merged_12_samples.rds
# Output: results/02_integrated_seurat.rds and reproduction of Xu et al. Fig. 1B/C
# Design: integration is performed by sample, never by clinical condition.
# ==============================================================================
#   - seurat_clusters are explicitly used as identities before FindAllMarkers()
#   - cell_type is assigned from the validated biological interpretation
#     of the 8 clusters using the paper's marker panel
#   - UMAP is plotted using cell_type
#
# ==============================================================================
# ------------------------------------------------------------------------------
# 0. Load packages and prepare folders
# ------------------------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(ggrepel)

RANDOM_SEED <- 42L
N_HVG <- 2000L
N_PCS <- 10L
CLUSTER_RESOLUTION <- 0.08
set.seed(RANDOM_SEED)

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)


# ------------------------------------------------------------------------------
# 1. Load Pre-Merged Post-QC Object
# ------------------------------------------------------------------------------

merged_rds_path <- "results/postQC_merged_12_samples.rds"
if (!file.exists(merged_rds_path)) {
  stop("Missing Stage-01 object: ", merged_rds_path,
       ". Run 01_qc_doublet_removal_and_merge.R first.")
}
combined_seurat <- readRDS(merged_rds_path)
if (!inherits(combined_seurat, "Seurat")) stop("Input is not a Seurat object.")
if (!"orig.ident" %in% colnames(combined_seurat[[]])) {
  stop("Required metadata column is missing: orig.ident")
}

DefaultAssay(combined_seurat) <- "RNA"

# ------------------------------------------------------------------------------
# 2. Normalization, HVGs, and Scaling
# ------------------------------------------------------------------------------

# Normalize RNA data
combined_seurat <- NormalizeData(
  combined_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identify highly variable genes
combined_seurat <- FindVariableFeatures(
  combined_seurat,
  selection.method = "vst",
  nfeatures = N_HVG
)

# Save top 10 HVGs
top10_hvgs <- head(VariableFeatures(combined_seurat), 10)

write.csv(
  data.frame(Gene = top10_hvgs),
  file = "results/top10_highly_variable_genes.csv",
  row.names = FALSE
)

# HVG plot
hvg_base <- VariableFeaturePlot(combined_seurat) +
  theme_classic()

hvg_data <- hvg_base$data

hvg_top10 <- hvg_data[top10_hvgs, ]

hvg_plot <- hvg_base +
  geom_label_repel(
    data = hvg_top10,
    aes(label = rownames(hvg_top10)),
    fontface = "bold.italic",
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    fill = alpha("white", 0.85),
    color = "black",
    max.overlaps = Inf
  ) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/01_highly_variable_genes.png",
  hvg_plot,
  width = 9,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# Regress out:
#   - total UMI counts
#   - mitochondrial percentage
#
# This follows the stated Methods approach.
combined_seurat <- ScaleData(
  combined_seurat,
  vars.to.regress = c("nCount_RNA", "percent.mt")
)
# ------------------------------------------------------------------------------
# 3. PCA
# ------------------------------------------------------------------------------

# Top 2000 HVGs -> first 10 PCs
combined_seurat <- RunPCA(
  combined_seurat,
  npcs = N_PCS,
  verbose = FALSE
)

# Elbow plot
elbow_p <- ElbowPlot(
  combined_seurat,
  ndims = N_PCS
) +
  theme_classic()

png(
  "figures/02_pca_elbow_plot.png",
  width = 1800,
  height = 1500,
  res = 300,
  bg = "white"
)

print(elbow_p)

dev.off()

# PCA heatmap
png(
  "figures/03_pca_heatmap_dims_1_10.png",
  width = 2400,
  height = 3000,
  res = 300,
  bg = "white"
)

DimHeatmap(
  combined_seurat,
  dims = seq_len(N_PCS),
  cells = 500,
  balanced = TRUE
)

dev.off()


# ------------------------------------------------------------------------------
# 4. Pre-integration sample diagnostic, then CCA integration
# ------------------------------------------------------------------------------

# The guide requires checking sample mixing before integration. Use a temporary
# object so this diagnostic does not alter the integrated analysis or its saved
# reductions. Inspect this figure and report the evidence supporting integration.
set.seed(RANDOM_SEED)
integration_check <- RunUMAP(
  combined_seurat,
  reduction = "pca",
  dims = seq_len(N_PCS),
  reduction.name = "umap.unintegrated",
  reduction.key = "preUMAP_",
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = RANDOM_SEED,
  verbose = FALSE
)

p_preintegration <- DimPlot(
  integration_check,
  reduction = "umap.unintegrated",
  group.by = "orig.ident",
  pt.size = 0.15,
  raster = FALSE
) +
  labs(
    title = "Pre-integration UMAP by biological sample",
    subtitle = "Use sample separation/mixing to justify the integration decision",
    colour = "Sample"
  ) +
  theme_classic()

ggsave(
  "figures/04_preintegration_UMAP_by_sample.png",
  p_preintegration,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)
rm(integration_check)
invisible(gc())

# Split RNA assay by sample
combined_seurat[["RNA"]] <- split(
  combined_seurat[["RNA"]],
  f = combined_seurat$orig.ident
)

message("Running Seurat v5 CCA Integration...")

combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = seq_len(N_PCS),
  verbose = FALSE
)

# Join RNA layers again
combined_seurat[["RNA"]] <- JoinLayers(
  combined_seurat[["RNA"]]
)
# ------------------------------------------------------------------------------
# 5. UMAP and Graph-Based Clustering
# ------------------------------------------------------------------------------

# UMAP using integrated CCA
combined_seurat <- RunUMAP(
  combined_seurat,
  reduction = "integrated.cca",
  dims = seq_len(N_PCS),
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  spread = 1.0,
  seed.use = RANDOM_SEED
)

# Find neighbors
combined_seurat <- FindNeighbors(
  combined_seurat,
  reduction = "integrated.cca",
  dims = seq_len(N_PCS)
)

# Clustering

combined_seurat <- FindClusters(
  combined_seurat,
  reduction = "integrated.cca",
  dims = seq_len(N_PCS),
  resolution = CLUSTER_RESOLUTION,
  random.seed = RANDOM_SEED
)

# Explicitly set cluster identities.

Idents(combined_seurat) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# 5. Check clustering BEFORE annotation
# ------------------------------------------------------------------------------

message("Cluster sizes:")

print(
  table(
    combined_seurat$seurat_clusters,
    useNA = "ifany"
  )
)

message("Cluster levels:")

print(
  levels(combined_seurat$seurat_clusters)
)

# Save cluster-only UMAP as a diagnostic
p_clusters <- DimPlot(
  combined_seurat,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  ggtitle("UMAP by Seurat Cluster") +
  theme_classic()

print(p_clusters)

ggsave(
  "figures/05_umap_by_seurat_cluster.png",
  p_clusters,
  width = 8,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 6. Marker Discovery
# ------------------------------------------------------------------------------

message(
  "Running FindAllMarkers using Wilcoxon test..."
)

cluster_markers <- FindAllMarkers(
  combined_seurat,
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Keep statistically significant markers
cluster_markers <- cluster_markers %>%
  filter(p_val < 0.05)

write.csv(
  cluster_markers,
  "results/cluster_markers_all.csv",
  row.names = FALSE
)

# Top 10 markers per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()

write.csv(
  top_markers,
  "results/top10_markers_per_cluster.csv",
  row.names = FALSE
)

print(
  top_markers %>%
    select(
      cluster,
      gene,
      avg_log2FC,
      pct.1,
      pct.2
    ) %>%
    arrange(
      cluster,
      desc(avg_log2FC)
    ),
  n = 80
)

# ------------------------------------------------------------------------------
# 7. Figure 1C Marker Panel
# ------------------------------------------------------------------------------

marker_panel <- list(
  
  "T cell" = c(
    "CD3D",
    "CD3E",
    "CD8A",
    "CD4"
  ),
  
  "Epithelia" = c(
    "KRT18",
    "EPCAM",
    "CD24",
    "KRT19"
  ),
  
  "Monocytic" = c(
    "CD14",
    "C1QA"
  ),
  
  "Endothelia" = c(
    "VCAM1",
    "PECAM1",
    "CLDN5"
  ),
  
  "Cell cycle cell" = c(
    "MKI67",
    "TOP2A"
  ),
  
  "Fibroblast" = c(
    "DCN",
    "OGN"
  ),
  
  "B cell_plasma" = c(
    "CD79A",
    "JCHAIN"
  ),
  
  "SMC_myoFibroblast" = c(
    "ACTA2",
    "MYH11",
    "TAGLN"
  )
)

# ------------------------------------------------------------------------------
# 7A. Check that all paper markers exist
# ------------------------------------------------------------------------------

marker_panel_present <- lapply(
  marker_panel,
  function(x) intersect(
    x,
    rownames(combined_seurat)
  )
)

print(marker_panel_present)


# ------------------------------------------------------------------------------
# 7B. Marker DotPlot,Fig.1C
# ------------------------------------------------------------------------------

dp <- DotPlot(
  combined_seurat,
  features = marker_panel_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      size = 7,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title.x = element_text(size = 12),
    panel.spacing.x = grid::unit(0, "pt"),
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 10,
      face = "plain"
    )
  )

print(dp)

ggsave(
  "figures/04_marker_dotplot_by_cluster.png",
  dp,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 7C. Average expression of paper markers
# ------------------------------------------------------------------------------

avg_marker_expr <- AverageExpression(
  combined_seurat,
  features = unique(
    unlist(marker_panel_present)
  ),
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)$RNA

paper_marker_summary <- round(
  avg_marker_expr,
  2
)

print(paper_marker_summary)

write.csv(
  paper_marker_summary,
  "results/paper_marker_average_expression.csv"
)

# ------------------------------------------------------------------------------
# 8. VALIDATED BIOLOGICAL ANNOTATION
# ------------------------------------------------------------------------------
# The mapping below is based on the observed marker expression in this
# dataset and the marker panel used for Figure 1C.


paper_cluster_map <- c(
  
  "0" = "T cell",
  
  "1" = "Epithelia",
  
  "2" = "Fibroblast",
  
  "3" = "Monocytic",
  
  "4" = "Endothelia",
  
  "5" = "Cell cycle cell",
  
  "6" = "B cell_plasma",
  
  "7" = "SMC_myoFibroblast"
)


# ------------------------------------------------------------------------------
# 8A. Apply cluster -> cell-type mapping
# ------------------------------------------------------------------------------

combined_seurat$cell_type <- unname(
  paper_cluster_map[
    as.character(
      combined_seurat$seurat_clusters
    )
  ]
)


# ------------------------------------------------------------------------------
# 8B. Verify annotation
# ------------------------------------------------------------------------------

message("Cell-type counts:")

print(
  table(
    combined_seurat$cell_type,
    useNA = "ifany"
  )
)

message("Cluster-to-cell-type mapping:")

print(
  table(
    combined_seurat$seurat_clusters,
    combined_seurat$cell_type,
    useNA = "ifany"
  )
)


# Check for missing annotations
if (any(is.na(combined_seurat$cell_type))) {
  
  warning(
    "Some cells have NA cell_type annotations. ",
    "Check whether all seurat_clusters are included in paper_cluster_map."
  )
  
} else {
  
  message(
    "SUCCESS: No cells have missing cell_type annotations."
  )
}


# Check that every cluster has a mapping
unmapped_clusters <- setdiff(
  levels(combined_seurat$seurat_clusters),
  names(paper_cluster_map)
)

if (length(unmapped_clusters) > 0) {
  
  warning(
    "Unmapped clusters: ",
    paste(unmapped_clusters, collapse = ", ")
  )
  
} else {
  
  message(
    "SUCCESS: All clusters have a cell-type assignment."
  )
}


# ------------------------------------------------------------------------------
# 8C. Save annotation table
# ------------------------------------------------------------------------------

annotation_table <- data.frame(
  cluster = names(paper_cluster_map),
  cell_type = unname(paper_cluster_map)
)

write.csv(
  annotation_table,
  "results/cluster_to_celltype_annotation.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Figure 1B UMAP
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 9A. Paper-style UMAP
# ------------------------------------------------------------------------------

paper_colors <- c(
  "T cell" = "#E6B8D7",
  "Epithelia" = "#4DAF4A",
  "Fibroblast" = "#56B4E9",
  "Monocytic" = "#E69F00",
  "Endothelia" = "#F4A582",
  "Cell cycle cell" = "#C7D36F",
  "B cell_plasma" = "#355C7D",
  "SMC_myoFibroblast" = "#D55E5E"
)

p_umap <- DimPlot(
  combined_seurat,
  reduction = "umap",
  group.by = "cell_type",
  cols = paper_colors,
  label = TRUE,
  repel = TRUE
)

print(p_umap)

ggsave(
  "figures/UMAP_paper_colors.png",
  p_umap,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
# Check that all actual cell types have colors
missing_colors <- setdiff(
  unique(combined_seurat$cell_type),
  names(paper_colors)
)

if (length(missing_colors) > 0) {
  
  stop(
    "Missing colors for: ",
    paste(missing_colors, collapse = ", ")
  )
  
}


# ------------------------------------------------------------------------------
# 9B. Final Figure 1B UMAP
# ------------------------------------------------------------------------------

p_fig1b <- DimPlot(
  combined_seurat,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3,
  cols = paper_colors
) +
  ggtitle(
    "Reproduced Figure 1B: Global UMAP Landscape"
  ) +
  theme_classic()

print(p_fig1b)


# Save Figure 1B
ggsave(
  "figures/Figure1B_reproduced_umap.png",
  p_fig1b,
  width = 8,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 10. Final Quality-Control Checks
# ------------------------------------------------------------------------------

message("================================================")
message("FINAL QC CHECKS")
message("================================================")

message("Number of cells:")
print(ncol(combined_seurat))

message("Number of genes:")
print(nrow(combined_seurat))

message("Clusters:")
print(table(combined_seurat$seurat_clusters))

message("Cell types:")
print(table(combined_seurat$cell_type, useNA = "ifany"))

message("Cell-type levels:")
print(unique(combined_seurat$cell_type))

message("Default assay:")
print(DefaultAssay(combined_seurat))

message("UMAP available:")
print("umap" %in% Reductions(combined_seurat))

message("Integrated CCA available:")
print("integrated.cca" %in% Reductions(combined_seurat))


# ------------------------------------------------------------------------------
# 11. Save Final Integrated Seurat Object
# ------------------------------------------------------------------------------

saveRDS(
  combined_seurat,
  "results/02_integrated_seurat.rds"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    "results",
    "session_info.txt"
  )
)

write.csv(
  data.frame(
    parameter = c("random_seed", "normalization", "scale_factor", "n_hvg",
                  "regressors", "integration", "n_pcs", "umap_neighbors",
                  "umap_min_dist", "cluster_resolution"),
    value = c(RANDOM_SEED, "LogNormalize", 10000, N_HVG,
              "nCount_RNA; percent.mt", "CCA by orig.ident", N_PCS, 30,
              0.3, CLUSTER_RESOLUTION)
  ),
  "results/02_analysis_parameters.csv",
  row.names = FALSE
)

message("================================================")
message("Module 1 complete.")
message("================================================")

message(
  "Final Figure 1B saved to: ",
  "figures/Figure1B_reproduced_umap.png"
)

message(
  "Final Seurat object saved to: ",
  "results/02_integrated_seurat.rds"
)

