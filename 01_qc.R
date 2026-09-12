# ==============================================================================
# 01_qc.R
# Pipeline: Per-Sample QC -> Doublet Detection -> MT Removal -> Merge
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Packages & Set Environment
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}

library(SingleCellExperiment)
library(scDblFinder)

set.seed(42)

# Directory configurations and QC thresholds
data_dir <- "D:/Abdalfttah Academy Single Cell RNA-seq Project/GSE184880/"
out_dir  <- "results"
fig_dir  <- "figures"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

MIN_GENES_PER_CELL <- 200
MAX_MITO_PCT        <- 40

# Helper function to auto-fix GEO file names
fix_geo_filenames <- function(sample_path) {
  all_files <- list.files(sample_path, full.names = TRUE)
  for (f in all_files) {
    fname <- basename(f)
    if (grepl("matrix\\.mtx", fname) && fname != "matrix.mtx.gz" && fname != "matrix.mtx") {
      file.rename(f, file.path(sample_path, ifelse(grepl("\\.gz$", fname), "matrix.mtx.gz", "matrix.mtx")))
    }
    if (grepl("barcodes\\.tsv", fname) && fname != "barcodes.tsv.gz" && fname != "barcodes.tsv") {
      file.rename(f, file.path(sample_path, ifelse(grepl("\\.gz$", fname), "barcodes.tsv.gz", "barcodes.tsv")))
    }
    if ((grepl("features\\.tsv", fname) || grepl("genes\\.tsv", fname)) && 
        !fname %in% c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")) {
      file.rename(f, file.path(sample_path, ifelse(grepl("\\.gz$", fname), "features.tsv.gz", "features.tsv")))
    }
  }
}

# ------------------------------------------------------------------------------
# 2. Load Data
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1_path <- file.path(data_dir, "GSM5599220_Norm1")
fix_geo_filenames(GSM5599220_Norm1_path)
GSM5599220_Norm1 <- CreateSeuratObject(
  counts = Read10X(GSM5599220_Norm1_path),
  project = "GSM5599220_Norm1",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL
)

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1_path <- file.path(data_dir, "GSM5599225_Cancer1")
fix_geo_filenames(GSM5599225_Cancer1_path)
GSM5599225_Cancer1 <- CreateSeuratObject(
  counts = Read10X(GSM5599225_Cancer1_path),
  project = "GSM5599225_Cancer1",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL
)

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2_path <- file.path(data_dir, "GSM5599226_Cancer2")
fix_geo_filenames(GSM5599226_Cancer2_path)
GSM5599226_Cancer2 <- CreateSeuratObject(
  counts = Read10X(GSM5599226_Cancer2_path),
  project = "GSM5599226_Cancer2",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL
)

# ------------------------------------------------------------------------------
# 3. Quality Control (Mitochondrial Gene Percentage Calculation)
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1[["percent.mt"]] <- PercentageFeatureSet(GSM5599220_Norm1, pattern = "^MT-")

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1[["percent.mt"]] <- PercentageFeatureSet(GSM5599225_Cancer1, pattern = "^MT-")

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2[["percent.mt"]] <- PercentageFeatureSet(GSM5599226_Cancer2, pattern = "^MT-")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 3.1 QC Visualization Before Filtering
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1_vln_before <- VlnPlot(GSM5599220_Norm1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599220_Norm1_vln_before
ggsave(file.path(fig_dir, "GSM5599220_Norm1_vln_before.png"), GSM5599220_Norm1_vln_before, width = 10, height = 4, dpi = 300)

GSM5599220_Norm1_sc1 <- FeatureScatter(GSM5599220_Norm1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599220_Norm1_sc2 <- FeatureScatter(GSM5599220_Norm1, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599220_Norm1_scatter <- GSM5599220_Norm1_sc1 + GSM5599220_Norm1_sc2
GSM5599220_Norm1_scatter
ggsave(file.path(fig_dir, "GSM5599220_Norm1_scatter.png"), GSM5599220_Norm1_sc1 + GSM5599220_Norm1_sc2, width = 10, height = 4, dpi = 300)

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1_vln_before <- VlnPlot(GSM5599225_Cancer1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599225_Cancer1_vln_before
ggsave(file.path(fig_dir, "GSM5599225_Cancer1_vln_before.png"), GSM5599225_Cancer1_vln_before, width = 10, height = 4, dpi = 300)

GSM5599225_Cancer1_sc1 <- FeatureScatter(GSM5599225_Cancer1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599225_Cancer1_sc2 <- FeatureScatter(GSM5599225_Cancer1, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599225_Cancer1_scatter <- GSM5599225_Cancer1_sc1 + GSM5599225_Cancer1_sc2
GSM5599225_Cancer1_scatter
ggsave(file.path(fig_dir, "GSM5599225_Cancer1_scatter.png"), GSM5599225_Cancer1_sc1 + GSM5599225_Cancer1_sc2, width = 10, height = 4, dpi = 300)

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2_vln_before <- VlnPlot(GSM5599226_Cancer2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599226_Cancer2_vln_before
ggsave(file.path(fig_dir, "GSM5599226_Cancer2_vln_before.png"), GSM5599226_Cancer2_vln_before, width = 10, height = 4, dpi = 300)

GSM5599226_Cancer2_sc1 <- FeatureScatter(GSM5599226_Cancer2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599226_Cancer2_sc2 <- FeatureScatter(GSM5599226_Cancer2, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599226_Cancer2_scatter <- GSM5599226_Cancer2_sc1 + GSM5599226_Cancer2_sc2
GSM5599226_Cancer2_scatter
ggsave(file.path(fig_dir, "GSM5599226_Cancer2_scatter.png"), GSM5599226_Cancer2_sc1 + GSM5599226_Cancer2_sc2, width = 10, height = 4, dpi = 300)
# ------------------------------------------------------------------------------
# 3.2 Cell Filtering
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1 <- subset(GSM5599220_Norm1, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1 <- subset(GSM5599225_Cancer1, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2 <- subset(GSM5599226_Cancer2, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# ------------------------------------------------------------------------------
# 3.3 QC Visualization After Filtering
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1_vln_after <- VlnPlot(GSM5599220_Norm1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599220_Norm1_vln_after
ggsave(file.path(fig_dir, "GSM5599220_Norm1_vln_after.png"), GSM5599220_Norm1_vln_after, width = 10, height = 4, dpi = 300)

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1_vln_after <- VlnPlot(GSM5599225_Cancer1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599225_Cancer1_vln_after
ggsave(file.path(fig_dir, "GSM5599225_Cancer1_vln_after.png"), GSM5599225_Cancer1_vln_after, width = 10, height = 4, dpi = 300)

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2_vln_after <- VlnPlot(GSM5599226_Cancer2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599226_Cancer2_vln_after
ggsave(file.path(fig_dir, "GSM5599226_Cancer2_vln_after.png"), GSM5599226_Cancer2_vln_after, width = 10, height = 4, dpi = 300)
# ------------------------------------------------------------------------------
# 4. Doublet Detection Using scDblFinder
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
sce_GSM5599220_Norm1 <- suppressWarnings(as.SingleCellExperiment(GSM5599220_Norm1))
set.seed(100)
sce_GSM5599220_Norm1 <- scDblFinder(sce_GSM5599220_Norm1)
GSM5599220_Norm1$doublet_score <- colData(sce_GSM5599220_Norm1)$scDblFinder.score
GSM5599220_Norm1$doublet_class <- colData(sce_GSM5599220_Norm1)$scDblFinder.class

GSM5599220_Norm1_dbl1 <- ggplot(GSM5599220_Norm1@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599220_Norm1_dbl2 <- VlnPlot(GSM5599220_Norm1, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599220_Norm1_doublets <- GSM5599220_Norm1_dbl1 + GSM5599220_Norm1_dbl2
GSM5599220_Norm1_doublets
ggsave(file.path(fig_dir, "GSM5599220_Norm1_doublets.png"), GSM5599220_Norm1_doublets, width = 12, height = 5, dpi = 300)

# Sample 2: GSM5599225_Cancer1
sce_GSM5599225_Cancer1 <- suppressWarnings(as.SingleCellExperiment(GSM5599225_Cancer1))
set.seed(100)
sce_GSM5599225_Cancer1 <- scDblFinder(sce_GSM5599225_Cancer1)
GSM5599225_Cancer1$doublet_score <- colData(sce_GSM5599225_Cancer1)$scDblFinder.score
GSM5599225_Cancer1$doublet_class <- colData(sce_GSM5599225_Cancer1)$scDblFinder.class

GSM5599225_Cancer1_dbl1 <- ggplot(GSM5599225_Cancer1@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599225_Cancer1_dbl2 <- VlnPlot(GSM5599225_Cancer1, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599225_Cancer1_doublets <- GSM5599225_Cancer1_dbl1 + GSM5599225_Cancer1_dbl2
GSM5599225_Cancer1_doublets
ggsave(file.path(fig_dir, "GSM5599225_Cancer1_doublets.png"), GSM5599225_Cancer1_doublets, width = 12, height = 5, dpi = 300)

# Sample 3: GSM5599226_Cancer2
sce_GSM5599226_Cancer2 <- suppressWarnings(as.SingleCellExperiment(GSM5599226_Cancer2))
set.seed(100)
sce_GSM5599226_Cancer2 <- scDblFinder(sce_GSM5599226_Cancer2)
GSM5599226_Cancer2$doublet_score <- colData(sce_GSM5599226_Cancer2)$scDblFinder.score
GSM5599226_Cancer2$doublet_class <- colData(sce_GSM5599226_Cancer2)$scDblFinder.class

GSM5599226_Cancer2_dbl1 <- ggplot(GSM5599226_Cancer2@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599226_Cancer2_dbl2 <- VlnPlot(GSM5599226_Cancer2, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599226_Cancer2_doublets <- GSM5599226_Cancer2_dbl1 + GSM5599226_Cancer2_dbl2
GSM5599226_Cancer2_doublets
ggsave(file.path(fig_dir, "GSM5599226_Cancer2_doublets.png"), GSM5599226_Cancer2_doublets, width = 12, height = 5, dpi = 300)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 4.1 Record QC Summary Table & Cell Counts (Run BEFORE Section 4.1 subsetting)
# ------------------------------------------------------------------------------
library(dplyr)

qc_summary_table <- data.frame(
  Sample = c("GSM5599220_Norm1", "GSM5599225_Cancer1", "GSM5599226_Cancer2"),
  
  # Initial cell counts (before filtering doublets out)
  Initial_Cells = c(
    ncol(GSM5599220_Norm1), 
    ncol(GSM5599225_Cancer1), 
    ncol(GSM5599226_Cancer2)
  ),
  
  # Count of detected doublets
  Doublets_Found = c(
    sum(GSM5599220_Norm1$doublet_class == "doublet", na.rm = TRUE),
    sum(GSM5599225_Cancer1$doublet_class == "doublet", na.rm = TRUE),
    sum(GSM5599226_Cancer2$doublet_class == "doublet", na.rm = TRUE)
  ),
  
  # Singlet count (cells that will remain after subsetting)
  Final_Singlets = c(
    sum(GSM5599220_Norm1$doublet_class == "singlet", na.rm = TRUE),
    sum(GSM5599225_Cancer1$doublet_class == "singlet", na.rm = TRUE),
    sum(GSM5599226_Cancer2$doublet_class == "singlet", na.rm = TRUE)
  ),
  
  Median_Genes = c(median(GSM5599220_Norm1$nFeature_RNA), median(GSM5599225_Cancer1$nFeature_RNA), median(GSM5599226_Cancer2$nFeature_RNA)),
  Median_UMIs  = c(median(GSM5599220_Norm1$nCount_RNA), median(GSM5599225_Cancer1$nCount_RNA), median(GSM5599226_Cancer2$nCount_RNA)),
  Mean_Mito_Pct = c(round(mean(GSM5599220_Norm1$percent.mt), 2), round(mean(GSM5599225_Cancer1$percent.mt), 2), round(mean(GSM5599226_Cancer2$percent.mt), 2))
) %>%
  mutate(
    Cells_Removed = Initial_Cells - Final_Singlets,
    Retention_Rate_Pct = round((Final_Singlets / Initial_Cells) * 100, 2)
  )

print(qc_summary_table)
write.csv(qc_summary_table, file = file.path(out_dir, "qc_cell_counts_summary.csv"), row.names = FALSE)
# ------------------------------------------------------------------------------
# 4.2 Remove Doublets
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1 <- subset(GSM5599220_Norm1, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599220_Norm1)

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1 <- subset(GSM5599225_Cancer1, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599225_Cancer1)

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2 <- subset(GSM5599226_Cancer2, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599226_Cancer2)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 4.3 Record Post-Doublet Subset Metrics & Verification
# ------------------------------------------------------------------------------
library(dplyr)

# Extract confirmed post-subset cell numbers and verify zero doublets remain
post_subset_summary_table <- data.frame(
  Sample = c("GSM5599220_Norm1", "GSM5599225_Cancer1", "GSM5599226_Cancer2"),
  
  # Final verified singlet cell count retained in the active objects
  Post_Subset_Singlets = c(
    ncol(GSM5599220_Norm1), 
    ncol(GSM5599225_Cancer1), 
    ncol(GSM5599226_Cancer2)
  ),
  
  # Remaining Doublets Check (Must equal 0)
  Remaining_Doublets = c(
    sum(GSM5599220_Norm1$doublet_class == "doublet", na.rm = TRUE),
    sum(GSM5599225_Cancer1$doublet_class == "doublet", na.rm = TRUE),
    sum(GSM5599226_Cancer2$doublet_class == "doublet", na.rm = TRUE)
  ),
  
  # Summary Quality Metrics on clean singlets
  Median_Genes  = c(median(GSM5599220_Norm1$nFeature_RNA), median(GSM5599225_Cancer1$nFeature_RNA), median(GSM5599226_Cancer2$nFeature_RNA)),
  Median_UMIs   = c(median(GSM5599220_Norm1$nCount_RNA), median(GSM5599225_Cancer1$nCount_RNA), median(GSM5599226_Cancer2$nCount_RNA)),
  Mean_Mito_Pct = c(round(mean(GSM5599220_Norm1$percent.mt), 2), round(mean(GSM5599225_Cancer1$percent.mt), 2), round(mean(GSM5599226_Cancer2$percent.mt), 2))
)

# Print verification table to console
print(post_subset_summary_table)

# Ensure output directory exists and export to CSV
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(post_subset_summary_table, file = file.path(out_dir, "post_doublet_subset_summary.csv"), row.names = FALSE)
# ------------------------------------------------------------------------------
# 5. Remove Mitochondrial Genes
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1_mt <- grep("^MT-", rownames(GSM5599220_Norm1), value = TRUE)
GSM5599220_Norm1    <- subset(GSM5599220_Norm1, features = setdiff(rownames(GSM5599220_Norm1), GSM5599220_Norm1_mt))

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1_mt <- grep("^MT-", rownames(GSM5599225_Cancer1), value = TRUE)
GSM5599225_Cancer1    <- subset(GSM5599225_Cancer1, features = setdiff(rownames(GSM5599225_Cancer1), GSM5599225_Cancer1_mt))

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2_mt <- grep("^MT-", rownames(GSM5599226_Cancer2), value = TRUE)
GSM5599226_Cancer2    <- subset(GSM5599226_Cancer2, features = setdiff(rownames(GSM5599226_Cancer2), GSM5599226_Cancer2_mt))
# Check gene dimensions after MT removal
cat("Genes remaining per sample after MT removal:\n")
cat("GSM5599220_Norm1:  ", nrow(GSM5599220_Norm1), "genes\n")
cat("GSM5599225_Cancer1: ", nrow(GSM5599225_Cancer1), "genes\n")
cat("GSM5599226_Cancer2: ", nrow(GSM5599226_Cancer2), "genes\n")

# ------------------------------------------------------------------------------
# 6. Merge Samples & Save RDS Output
# ------------------------------------------------------------------------------
# Ensure output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Merge individual Seurat objects into a single combined object
merged_seurat <- merge(
  x = GSM5599220_Norm1,
  y = c(GSM5599225_Cancer1, GSM5599226_Cancer2),
  add.cell.ids = c("GSM5599220_Norm1", "GSM5599225_Cancer1", "GSM5599226_Cancer2")
)

# Export final merged post-QC Seurat object to RDS
rds_path <- file.path(out_dir, "postQC_GSM5599220-GSM5599225-GSM5599226.rds")
saveRDS(merged_seurat, file = rds_path)

# Console print status update
cat(sprintf(
  "\nSuccessfully completed Stage 1 QC pipeline!\nMerged object saved to: %s\nTotal dataset contains %d genes across %d cells.\n",
  rds_path,
  dim(merged_seurat)[1],
  dim(merged_seurat)[2]
))
# ------------------------------------------------------------------------------
# 6.1 Inspect Merged Seurat Object
# ------------------------------------------------------------------------------

# 1. Check overall object dimensions (Genes x Cells)
dim(merged_seurat)

# 2. Inspect active assay and Seurat v5 layer structure
# In Seurat v5, count matrices are stored as individual sample layers in RNA assay
merged_seurat[["RNA"]]
Layers(merged_seurat)

# 3. View metadata structure and available column names
head(merged_seurat@meta.data, 5)
colnames(merged_seurat@meta.data)

# 4. Check total cell counts per sample
table(merged_seurat$orig.ident)

# 5. Summarize key QC metrics across samples
aggregate(
  cbind(nCount_RNA, nFeature_RNA) ~ orig.ident, 
  data = merged_seurat@meta.data, 
  FUN = median
)
# Create the results directory if it does not exist
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# Save the merged Seurat object locally
saveRDS(merged_seurat, file = "results/postQC_GSM5599220-GSM5599225-GSM5599226.rds")
getwd()
# Opens your current working directory in your computer's file manager
shell.exec(getwd())  # For Windows
# system("open .")   # Use this instead if you are on Mac
file.edit(".gitignore")
write("*.rds", file = ".gitignore", append = TRUE)