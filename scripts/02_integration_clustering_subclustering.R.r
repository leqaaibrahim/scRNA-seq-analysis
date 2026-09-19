# ==============================================================================
# INTEGRATED scRNA-SEQ ANALYSIS: PREPROCESSING, INTEGRATION, & SUBCLUSTERING
# ==============================================================================

library(Seurat)
library(scran)             # Required for MNN subclustering
library(SingleCellExperiment) # For object conversion
library(SeuratWrappers)    # Enables RunMNN integration
library(ggplot2)
library(dplyr)

set.seed(42)

# Ensure output directories exist
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load Pre-Merged Post-QC Object
# ------------------------------------------------------------------------------
merged_rds_path <- "results/postQC_merged_12_samples.rds"

if (!file.exists(merged_rds_path)) {
  stop("Merged RDS file not found in results/ directory!")
}

combined_seurat <- readRDS(merged_rds_path)

# ------------------------------------------------------------------------------
# 2. Normalization & Feature Selection
# ------------------------------------------------------------------------------
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)

# --- HVG Plotting & Export ---
top10_hvgs <- head(VariableFeatures(combined_seurat), 10)

# Save top 10 HVGs to CSV
write.csv(data.frame(Gene = top10_hvgs), file = "results/top10_highly_variable_genes.csv", row.names = FALSE)

# Generate and save labeled HVG Plot
hvg_plot <- LabelPoints(
  plot = VariableFeaturePlot(combined_seurat),
  points = top10_hvgs,
  repel = TRUE
)
ggsave(filename = "figures/01_highly_variable_genes.png", plot = hvg_plot, width = 8, height = 6, dpi = 300)

# Scale data and regress out UMI counts and % mitochondrial reads
combined_seurat <- ScaleData(
  combined_seurat, 
  vars.to.regress = c("nCount_RNA", "percent.mt")
)

# ------------------------------------------------------------------------------
# 3. PCA & Visualizations
# ------------------------------------------------------------------------------
combined_seurat <- RunPCA(combined_seurat, npcs = 30, verbose = FALSE)

# --- Save Elbow Plot ---
png("figures/02_pca_elbow_plot.png", width = 1800, height = 1500, res = 300)
print(ElbowPlot(combined_seurat, ndims = 30))
dev.off()

# --- Save PCA Heatmap ---
png("figures/03_pca_heatmap_dims_1_10.png", width = 2400, height = 3000, res = 300)
DimHeatmap(
  combined_seurat,
  dims = 1:10,
  cells = 500,
  balanced = TRUE
)
dev.off()

# ------------------------------------------------------------------------------
# 4. Seurat v5 CCA Integration & Primary UMAP Clustering
# ------------------------------------------------------------------------------
# Integrate across sample identity layers
combined_seurat[["RNA"]] <- split(combined_seurat[["RNA"]], f = combined_seurat$orig.ident)

combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# Re-join layers for downstream RNA-level analyses
combined_seurat[["RNA"]] <- JoinLayers(combined_seurat[["RNA"]])

# Construct UMAP and Graph-Based Clusters using top 10 PCs
combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.cca", dims = 1:10)
combined_seurat <- FindNeighbors(combined_seurat, reduction = "integrated.cca", dims = 1:10)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# --- Export Cluster Counts CSV ---
cluster_counts <- as.data.frame(table(Cluster = combined_seurat$seurat_clusters))
write.csv(cluster_counts, file = "results/cell_counts_per_cluster.csv", row.names = FALSE)

# --- Save Primary UMAP Plot ---
umap_plot <- DimPlot(
  combined_seurat,
  reduction = "umap",
  label = TRUE
) + ggtitle("Merged 12-Sample Dataset — Primary UMAP")

ggsave(filename = "figures/04_primary_umap_clusters.png", plot = umap_plot, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Marker Gene Detection & Top Marker Visualizations
# ------------------------------------------------------------------------------
DefaultAssay(combined_seurat) <- "RNA"

all_markers <- FindAllMarkers(
  combined_seurat,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

signif_markers <- subset(all_markers, p_val_adj < 0.05)
write.csv(signif_markers, file = "results/cluster_marker_genes.csv", row.names = FALSE)

# Extract top 5 markers per cluster by log2FC
top5_markers <- signif_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# --- Save DotPlot of Top Markers ---
dot_plot <- DotPlot(combined_seurat, features = unique(top5_markers$gene)) +
  RotatedAxis() +
  ggtitle("Top 5 Marker Genes per Cluster")

ggsave(filename = "figures/05_top5_markers_dotplot.png", plot = dot_plot, width = 12, height = 8, dpi = 300)

# --- Save Heatmap of Top Markers ---
png("figures/06_top5_markers_heatmap.png", width = 3000, height = 2400, res = 300)
print(DoHeatmap(combined_seurat, features = top5_markers$gene) + NoLegend())
dev.off()

# ------------------------------------------------------------------------------
# 6. Subclustering Specific Cell Types using MNN (k = 5)
# ------------------------------------------------------------------------------
# Subset major cluster of interest (e.g., Cluster 0)
sub_celltype <- subset(combined_seurat, idents = "0")

# FastMNN integration via SeuratWrappers
sub_celltype <- RunFastMNN(
  object.list = SplitObject(sub_celltype, split.by = "orig.ident"),
  k = 5
)

# Run UMAP and Graph-based clustering on MNN dimensions
sub_celltype <- RunUMAP(sub_celltype, reduction = "mnn", dims = 1:10)
sub_celltype <- FindNeighbors(sub_celltype, reduction = "mnn", dims = 1:10)
sub_celltype <- FindClusters(sub_celltype, resolution = 0.3)

# Save Subcluster UMAP Plot
sub_umap_plot <- DimPlot(sub_celltype, reduction = "umap", label = TRUE) + ggtitle("Cluster 0 Subclusters (MNN)")
ggsave(filename = "figures/07_cluster0_subclusters_mnn.png", plot = sub_umap_plot, width = 8, height = 6, dpi = 300)

# Save final processed object with integrated reductions and clusters
saveRDS(combined_seurat, file = "results/postQC_integrated_clustered_12_samples.rds")