library(Seurat)
library(ggplot2)

# 1. Load Integrated Object
message("Loading integrated Seurat object...")
combined_seurat <- readRDS("results/02_integrated_seurat.rds")

# Marker panel (نفس اللي في Module 1)
cell_type_markers <- list(
  "T cell" = c("CD3D", "CD3E", "CD8A", "CD4"),
  "Epithelia" = c("KRT18", "EPCAM", "CD24", "KRT19"),
  "Monocytic" = c("CD14", "C1QA"),
  "Endothelia" = c("VCAM1", "PECAM1", "CLDN5"),
  "Cell cycle cell" = c("MKI67", "TOP2A"),
  "Fibroblast" = c("DCN", "OGN"),
  "B cell_plasma" = c("CD79A", "JCHAIN"),
  "SMC_myoFibroblast" = c("ACTA2", "MYH11", "TAGLN")
)

# 2. Fig. 1C
message("Generating Fig. 1C...")
fig1_genes <- unlist(cell_type_markers, use.names = FALSE)
fig1_genes <- intersect(unique(fig1_genes), rownames(combined_seurat))

p_1c <- DotPlot(combined_seurat, features = fig1_genes, group.by = "cell_type") +
  RotatedAxis() +
  labs(title = "Figure 1C: Canonical Marker Gene Expression",
       x = "Marker Genes", y = "Cell Lineage")

ggsave("figures/fig1c_markers_dotplot.png", plot = p_1c,
       width = 10, height = 5, dpi = 300, bg = "white")

# 3. Fig. 1E
message("Generating Fig. 1E...")
if ("stage" %in% colnames(combined_seurat@meta.data)) {
  p_1e <- DimPlot(combined_seurat, reduction = "umap",
                  split.by = "stage", group.by = "cell_type") +
    theme_classic() +
    ggtitle("Figure 1E: Cell Types Across Clinical Stages")
  ggsave("figures/fig1e_stage_umap.png", plot = p_1e,
         width = 14, height = 4, dpi = 300, bg = "white")
} else {
  warning("'stage' column not found in metadata. Skipping Fig. 1E export.")
}

# 4. Save
saveRDS(combined_seurat, "results/02_annotated_main_seurat.rds")
message("Module 2 complete.")