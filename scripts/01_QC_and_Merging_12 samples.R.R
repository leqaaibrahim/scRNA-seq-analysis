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
# Sample 4: GSM5599224_Norm5
GSM5599224_Norm5_path <- file.path(data_dir, "GSM5599224_Norm5")
fix_geo_filenames(GSM5599224_Norm5_path)
GSM5599224_Norm5 <- CreateSeuratObject(
  counts = Read10X(GSM5599224_Norm5_path),
  project = "GSM5599224_Norm5",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL
)
# Sample 5: GSM5599231_Cancer7
GSM5599231_Cancer7_path <- file.path(data_dir, "GSM5599231_Cancer7")
fix_geo_filenames(GSM5599231_Cancer7_path)
GSM5599231_Cancer7 <- CreateSeuratObject(
  counts = Read10X(GSM5599231_Cancer7_path),
  project = "GSM5599231_Cancer7",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL
)
# Sample 6: GSM5599223_Norm4
GSM5599223_Norm4_path <- file.path(data_dir, "GSM5599223_Norm4")
fix_geo_filenames(GSM5599223_Norm4_path)
GSM5599223_Norm4 <- CreateSeuratObject(
  counts = Read10X(GSM5599223_Norm4_path),
  project = "GSM5599223_Norm4",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)

# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5_path <- file.path(data_dir, "GSM5599229_Cancer5")
fix_geo_filenames(GSM5599229_Cancer5_path)
GSM5599229_Cancer5 <- CreateSeuratObject(
  counts = Read10X(GSM5599229_Cancer5_path),
  project = "GSM5599229_Cancer5",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)

# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6_path <- file.path(data_dir, "GSM5599230_Cancer6")
fix_geo_filenames(GSM5599230_Cancer6_path)
GSM5599230_Cancer6 <- CreateSeuratObject(
  counts = Read10X(GSM5599230_Cancer6_path),
  project = "GSM5599230_Cancer6",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4_path <- file.path(data_dir, "GSM5599228_Cancer4")
fix_geo_filenames(GSM5599228_Cancer4_path)
GSM5599228_Cancer4 <- CreateSeuratObject(
  counts = Read10X(GSM5599228_Cancer4_path),
  project = "GSM5599228_Cancer4",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3_path <- file.path(data_dir, "GSM5599227_Cancer3")
fix_geo_filenames(GSM5599227_Cancer3_path)
GSM5599227_Cancer3 <- CreateSeuratObject(
  counts = Read10X(GSM5599227_Cancer3_path),
  project = "GSM5599227_Cancer3",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3_path <- file.path(data_dir, "GSM5599222_Norm3")
fix_geo_filenames(GSM5599222_Norm3_path)
GSM5599222_Norm3 <- CreateSeuratObject(
  counts = Read10X(GSM5599222_Norm3_path),
  project = "GSM5599222_Norm3",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2_path <- file.path(data_dir, "GSM5599221_Norm2")
fix_geo_filenames(GSM5599221_Norm2_path)
GSM5599221_Norm2 <- CreateSeuratObject(
  counts = Read10X(GSM5599221_Norm2_path),
  project = "GSM5599221_Norm2",
  min.cells = 3,
  min.features = MIN_GENES_PER_CELL)
# ------------------------------------------------------------------------------
# 3. Quality Control (Mitochondrial Gene Percentage Calculation)
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1[["percent.mt"]] <- PercentageFeatureSet(GSM5599220_Norm1, pattern = "^MT-")

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1[["percent.mt"]] <- PercentageFeatureSet(GSM5599225_Cancer1, pattern = "^MT-")

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2[["percent.mt"]] <- PercentageFeatureSet(GSM5599226_Cancer2, pattern = "^MT-")
# Sample 4: GSM5599224_Norm5
GSM5599224_Norm5[["percent.mt"]] <- PercentageFeatureSet(GSM5599224_Norm5, pattern = "^MT-")
# Sample 5: GSM5599231_Cancer7
GSM5599231_Cancer7[["percent.mt"]] <- PercentageFeatureSet(GSM5599231_Cancer7, pattern = "^MT-")
# Sample 6: GSM5599223_Norm4
GSM5599223_Norm4[["percent.mt"]] <- PercentageFeatureSet(GSM5599223_Norm4, pattern = "^MT-")
# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5[["percent.mt"]] <- PercentageFeatureSet(GSM5599229_Cancer5, pattern = "^MT-")
# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6[["percent.mt"]] <- PercentageFeatureSet(GSM5599230_Cancer6, pattern = "^MT-")
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4[["percent.mt"]] <- PercentageFeatureSet(GSM5599228_Cancer4, pattern = "^MT-")
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3[["percent.mt"]] <- PercentageFeatureSet(GSM5599227_Cancer3, pattern = "^MT-")
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3[["percent.mt"]] <- PercentageFeatureSet(GSM5599222_Norm3, pattern = "^MT-")
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2[["percent.mt"]] <- PercentageFeatureSet(GSM5599221_Norm2, pattern = "^MT-")
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

# Sample 4:  GSM5599224_Norm5
GSM5599224_Norm5_vln_before <- VlnPlot( GSM5599224_Norm5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599224_Norm5_vln_before
ggsave(file.path(fig_dir, " GSM5599224_Norm5_vln_before.png"),  GSM5599224_Norm5_vln_before, width = 10, height = 4, dpi = 300)

GSM5599224_Norm5_sc1 <- FeatureScatter( GSM5599224_Norm5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599224_Norm5_sc2 <- FeatureScatter( GSM5599224_Norm5, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599224_Norm5_scatter <-  GSM5599224_Norm5_sc1 +  GSM5599224_Norm5_sc2
GSM5599224_Norm5_scatter
ggsave(file.path(fig_dir, " GSM5599224_Norm5.png"),  GSM5599224_Norm5_sc1 +  GSM5599224_Norm5_sc2, width = 10, height = 4, dpi = 300)
# Sample 5:  GSM5599231_Cancer7
GSM5599231_Cancer7_vln_before <- VlnPlot( GSM5599231_Cancer7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599231_Cancer7_vln_before
ggsave(file.path(fig_dir, " GSM5599231_Cancer7_vln_before.png"),  GSM5599231_Cancer7_vln_before, width = 10, height = 4, dpi = 300)

GSM5599231_Cancer7_sc1 <- FeatureScatter( GSM5599231_Cancer7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599231_Cancer7_sc2 <- FeatureScatter( GSM5599231_Cancer7, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599231_Cancer7_scatter <-  GSM5599231_Cancer7_sc1 +  GSM5599231_Cancer7_sc2
GSM5599231_Cancer7_scatter
ggsave(file.path(fig_dir, " GSM5599231_Cancer7.png"),  GSM5599231_Cancer7_sc1 +  GSM5599231_Cancer7_sc2, width = 10, height = 4, dpi = 300)
# Sample 6:  GSM5599223_Norm4
GSM5599223_Norm4_vln_before <- VlnPlot( GSM5599223_Norm4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599223_Norm4_vln_before
ggsave(file.path(fig_dir, " GSM5599223_Norm4_vln_before.png"),  GSM5599223_Norm4_vln_before, width = 10, height = 4, dpi = 300)

GSM5599223_Norm4_sc1 <- FeatureScatter( GSM5599223_Norm4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599223_Norm4_sc2 <- FeatureScatter(GSM5599223_Norm4, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599223_Norm4_scatter <-  GSM5599223_Norm4_sc1 +  GSM5599223_Norm4_sc2
GSM5599223_Norm4_scatter
ggsave(file.path(fig_dir, " GSM5599223_Norm4.png"),  GSM5599223_Norm4_sc1 +  GSM5599223_Norm4_sc2, width = 10, height = 4, dpi = 300)
# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5_vln_before <- VlnPlot( GSM5599229_Cancer5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599229_Cancer5_vln_before
ggsave(file.path(fig_dir, " GSM5599229_Cancer5_vln_before.png"),  GSM5599229_Cancer5_vln_before, width = 10, height = 4, dpi = 300)

GSM5599229_Cancer5_sc1 <- FeatureScatter( GSM5599229_Cancer5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599229_Cancer5_sc2 <- FeatureScatter(GSM5599229_Cancer5, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599229_Cancer5_scatter <-  GSM5599229_Cancer5_sc1 +  GSM5599229_Cancer5_sc2
GSM5599229_Cancer5_scatter
ggsave(file.path(fig_dir, " GSM5599229_Cancer5.png"),  GSM5599229_Cancer5_sc1 +  GSM5599229_Cancer5_sc2, width = 10, height = 4, dpi = 300)
# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6_vln_before <- VlnPlot( GSM5599230_Cancer6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599230_Cancer6_vln_before
ggsave(file.path(fig_dir, " GSM5599230_Cancer6_vln_before.png"),  GSM5599230_Cancer6_vln_before, width = 10, height = 4, dpi = 300)

GSM5599230_Cancer6_sc1 <- FeatureScatter( GSM5599230_Cancer6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599230_Cancer6_sc2 <- FeatureScatter(GSM5599230_Cancer6, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599230_Cancer6_scatter <-  GSM5599230_Cancer6_sc1 +  GSM5599230_Cancer6_sc2
GSM5599230_Cancer6_scatter
ggsave(file.path(fig_dir, " GSM5599230_Cancer6.png"),  GSM5599230_Cancer6_sc1 +  GSM5599230_Cancer6_sc2, width = 10, height = 4, dpi = 300)
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4_vln_before <- VlnPlot( GSM5599228_Cancer4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599228_Cancer4_vln_before
ggsave(file.path(fig_dir, " GSM5599228_Cancer4_vln_before.png"),  GSM5599228_Cancer4_vln_before, width = 10, height = 4, dpi = 300)

GSM5599228_Cancer4_sc1 <- FeatureScatter( GSM5599228_Cancer4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599228_Cancer4_sc2 <- FeatureScatter(GSM5599228_Cancer4, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599228_Cancer4_scatter <-  GSM5599228_Cancer4_sc1 + GSM5599228_Cancer4_sc2
GSM5599228_Cancer4_scatter
ggsave(file.path(fig_dir, " GSM5599228_Cancer4.png"),  GSM5599228_Cancer4_sc1 +  GSM5599228_Cancer4_sc2, width = 10, height = 4, dpi = 300)
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3_vln_before <- VlnPlot( GSM5599227_Cancer3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599227_Cancer3_vln_before
ggsave(file.path(fig_dir, " GSM5599227_Cancer3_vln_before.png"),  GSM5599227_Cancer3_vln_before, width = 10, height = 4, dpi = 300)

GSM5599227_Cancer3_sc1 <- FeatureScatter( GSM5599227_Cancer3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599227_Cancer3_sc2 <- FeatureScatter(GSM5599227_Cancer3, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599227_Cancer3_scatter <-  GSM5599227_Cancer3_sc1 + GSM5599227_Cancer3_sc2
GSM5599227_Cancer3_scatter
ggsave(file.path(fig_dir, " GSM5599227_Cancer3.png"),  GSM5599227_Cancer3_sc1 +  GSM5599227_Cancer3_sc2, width = 10, height = 4, dpi = 300)
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3_vln_before <- VlnPlot( GSM5599222_Norm3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599222_Norm3_vln_before
ggsave(file.path(fig_dir, " GSM5599222_Norm3_vln_before.png"),  GSM5599222_Norm3_vln_before, width = 10, height = 4, dpi = 300)

GSM5599222_Norm3_sc1 <- FeatureScatter( GSM5599222_Norm3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599222_Norm3_sc2 <- FeatureScatter(GSM5599222_Norm3, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599222_Norm3_scatter <-  GSM5599222_Norm3_sc1 + GSM5599222_Norm3_sc2
GSM5599222_Norm3_scatter
ggsave(file.path(fig_dir, " GSM5599222_Norm3.png"),  GSM5599222_Norm3_sc1 +  GSM5599222_Norm3_sc2, width = 10, height = 4, dpi = 300)
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2_vln_before <- VlnPlot( GSM5599221_Norm2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599221_Norm2_vln_before
ggsave(file.path(fig_dir, " GSM5599221_Norm2_vln_before.png"),  GSM5599221_Norm2_vln_before, width = 10, height = 4, dpi = 300)

GSM5599221_Norm2_sc1 <- FeatureScatter( GSM5599221_Norm2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", slot = "counts")
GSM5599221_Norm2_sc2 <- FeatureScatter(GSM5599221_Norm2, feature1 = "nCount_RNA", feature2 = "percent.mt", slot = "counts")
GSM5599221_Norm2_scatter <-  GSM5599221_Norm2_sc1 + GSM5599221_Norm2_sc2
GSM5599221_Norm2_scatter
ggsave(file.path(fig_dir, " GSM5599221_Norm2.png"),  GSM5599221_Norm2_sc1 +  GSM5599221_Norm2_sc2, width = 10, height = 4, dpi = 300)
# ------------------------------------------------------------------------------
# 3.2 Cell Filtering
# ------------------------------------------------------------------------------
# Sample 1: GSM5599220_Norm1
GSM5599220_Norm1 <- subset(GSM5599220_Norm1, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)

# Sample 2: GSM5599225_Cancer1
GSM5599225_Cancer1 <- subset(GSM5599225_Cancer1, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)

# Sample 3: GSM5599226_Cancer2
GSM5599226_Cancer2 <- subset(GSM5599226_Cancer2, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 4: GSM5599224_Norm5
GSM5599224_Norm5 <- subset(GSM5599224_Norm5, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 5: GSM5599231_Cancer7
GSM5599231_Cancer7 <- subset(GSM5599231_Cancer7, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 6: GSM5599223_Norm4
GSM5599223_Norm4 <- subset(GSM5599223_Norm4, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5 <- subset(GSM5599229_Cancer5, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6 <- subset(GSM5599230_Cancer6, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4 <- subset(GSM5599228_Cancer4, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3 <- subset(GSM5599227_Cancer3, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3 <- subset(GSM5599222_Norm3, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2 <- subset(GSM5599221_Norm2, subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)
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
# Sample 4: GSM5599224_Norm5
GSM5599224_Norm5_vln_after <- VlnPlot(GSM5599224_Norm5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599224_Norm5_vln_after
ggsave(file.path(fig_dir, "GSM5599224_Norm5_vln_after.png"), GSM5599224_Norm5_vln_after, width = 10, height = 4, dpi = 300)
# Sample 5: GSM5599231_Cancer7
GSM5599231_Cancer7_vln_after <- VlnPlot(GSM5599231_Cancer7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599231_Cancer7_vln_after
ggsave(file.path(fig_dir, "GSM5599231_Cancer7_vln_after.png"), GSM5599226_Cancer2_vln_after, width = 10, height = 4, dpi = 300)
# Sample 6: GSM5599223_Norm4
GSM5599223_Norm4_vln_after <- VlnPlot(GSM5599223_Norm4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599223_Norm4_vln_after
ggsave(file.path(fig_dir, "GSM5599223_Norm4_vln_after.png"), GSM5599223_Norm4_vln_after, width = 10, height = 4, dpi = 300)
# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5_vln_after <- VlnPlot( GSM5599229_Cancer5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599229_Cancer5_vln_after
ggsave(file.path(fig_dir, " GSM5599229_Cancer5_vln_after.png"),  GSM5599229_Cancer5_vln_after, width = 10, height = 4, dpi = 300)
# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6_vln_after <- VlnPlot(GSM5599230_Cancer6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599230_Cancer6_vln_after
ggsave(file.path(fig_dir, "GSM5599230_Cancer6_vln_after.png"), GSM5599230_Cancer6_vln_after, width = 10, height = 4, dpi = 300)
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4_vln_after <- VlnPlot(GSM5599228_Cancer4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599228_Cancer4_vln_after
ggsave(file.path(fig_dir, "GSM5599228_Cancer4_vln_after.png"), GSM5599228_Cancer4_vln_after, width = 10, height = 4, dpi = 300)
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3_vln_after <- VlnPlot(GSM5599227_Cancer3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599227_Cancer3_vln_after
ggsave(file.path(fig_dir, "GSM5599227_Cancer3_vln_after.png"), GSM5599227_Cancer3_vln_after, width = 10, height = 4, dpi = 300)
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3_vln_after <- VlnPlot(GSM5599222_Norm3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599222_Norm3_vln_after
ggsave(file.path(fig_dir, "GSM5599222_Norm3_vln_after.png"), GSM5599222_Norm3_vln_after, width = 10, height = 4, dpi = 300)
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2_vln_after <- VlnPlot(GSM5599221_Norm2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts")
GSM5599221_Norm2_vln_after
ggsave(file.path(fig_dir, "GSM5599221_Norm2_vln_after.png"), GSM5599221_Norm2_vln_after, width = 10, height = 4, dpi = 300)
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
# Sample 4: GSM5599224_Norm5
sce_GSM5599224_Norm5 <- suppressWarnings(as.SingleCellExperiment(GSM5599224_Norm5))
set.seed(100)
sce_GSM5599224_Norm5 <- scDblFinder(sce_GSM5599224_Norm5)
GSM5599224_Norm5$doublet_score <- colData(sce_GSM5599224_Norm5)$scDblFinder.score
GSM5599224_Norm5$doublet_class <- colData(sce_GSM5599224_Norm5)$scDblFinder.class

GSM5599224_Norm5_dbl1 <- ggplot(GSM5599224_Norm5@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599224_Norm5_dbl2 <- VlnPlot(GSM5599224_Norm5, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599224_Norm5_doublets <- GSM5599224_Norm5_dbl1 + GSM5599224_Norm5_dbl2
GSM5599224_Norm5_doublets
ggsave(file.path(fig_dir, "GSM5599224_Norm5_doublets.png"), GSM5599224_Norm5_doublets, width = 12, height = 5, dpi = 300)
# Sample 5: GSM5599231_Cancer7
sce_GSM5599231_Cancer7 <- suppressWarnings(as.SingleCellExperiment(GSM5599231_Cancer7))
set.seed(100)
sce_GSM5599231_Cancer7 <- scDblFinder(sce_GSM5599231_Cancer7)
GSM5599231_Cancer7$doublet_score <- colData(sce_GSM5599231_Cancer7)$scDblFinder.score
GSM5599231_Cancer7$doublet_class <- colData(sce_GSM5599231_Cancer7)$scDblFinder.class

GSM5599231_Cancer7_dbl1 <- ggplot(GSM5599231_Cancer7@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599231_Cancer7_dbl2 <- VlnPlot(GSM5599231_Cancer7, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599231_Cancer7_doublets <- GSM5599231_Cancer7_dbl1 + GSM5599231_Cancer7_dbl2
GSM5599231_Cancer7_doublets
ggsave(file.path(fig_dir, "GSM5599231_Cancer7_doublets.png"), GSM5599231_Cancer7_doublets, width = 12, height = 5, dpi = 300)
# Sample 6: GSM5599223_Norm4
sce_GSM5599223_Norm4 <- suppressWarnings(as.SingleCellExperiment(GSM5599223_Norm4))
set.seed(100)
sce_GSM5599223_Norm4 <- scDblFinder(sce_GSM5599223_Norm4)
GSM5599223_Norm4$doublet_score <- colData(sce_GSM5599223_Norm4)$scDblFinder.score
GSM5599223_Norm4$doublet_class <- colData(sce_GSM5599223_Norm4)$scDblFinder.class

GSM5599223_Norm4_dbl1 <- ggplot(GSM5599223_Norm4@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599223_Norm4_dbl2 <- VlnPlot(GSM5599223_Norm4, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599223_Norm4_doublets <- GSM5599223_Norm4_dbl1 + GSM5599223_Norm4_dbl2
GSM5599223_Norm4_doublets
ggsave(file.path(fig_dir, "GSM5599223_Norm4_doublets.png"),GSM5599223_Norm4_doublets, width = 12, height = 5, dpi = 300)
# Sample 7: GSM5599229_Cancer5
sce_GSM5599229_Cancer5 <- suppressWarnings(as.SingleCellExperiment(GSM5599229_Cancer5))
set.seed(100)
sce_GSM5599229_Cancer5 <- scDblFinder(sce_GSM5599229_Cancer5)
GSM5599229_Cancer5$doublet_score <- colData(sce_GSM5599229_Cancer5)$scDblFinder.score
GSM5599229_Cancer5$doublet_class <- colData(sce_GSM5599229_Cancer5)$scDblFinder.class

GSM5599229_Cancer5_dbl1 <- ggplot(GSM5599229_Cancer5@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599229_Cancer5_dbl2 <- VlnPlot(GSM5599229_Cancer5, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599229_Cancer5_doublets <- GSM5599229_Cancer5_dbl1 + GSM5599229_Cancer5_dbl2
GSM5599229_Cancer5_doublets
ggsave(file.path(fig_dir, "GSM5599229_Cancer5_doublets.png"), GSM5599229_Cancer5_doublets, width = 12, height = 5, dpi = 300)
# Sample 8: GSM5599230_Cancer6
sce_GSM5599230_Cancer6 <- suppressWarnings(as.SingleCellExperiment(GSM5599230_Cancer6))
set.seed(100)
sce_GSM5599230_Cancer6 <- scDblFinder(sce_GSM5599230_Cancer6)
GSM5599230_Cancer6$doublet_score <- colData(sce_GSM5599230_Cancer6)$scDblFinder.score
GSM5599230_Cancer6$doublet_class <- colData(sce_GSM5599230_Cancer6)$scDblFinder.class

GSM5599230_Cancer6_dbl1 <- ggplot(GSM5599230_Cancer6@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599230_Cancer6_dbl2 <- VlnPlot(GSM5599230_Cancer6, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599230_Cancer6_doublets <- GSM5599230_Cancer6_dbl1 + GSM5599230_Cancer6_dbl2
GSM5599230_Cancer6_doublets
ggsave(file.path(fig_dir, "GSM5599230_Cancer6_doublets.png"), GSM5599230_Cancer6_doublets, width = 12, height = 5, dpi = 300)
# Sample 9: GSM5599228_Cancer4
sce_GSM5599228_Cancer4 <- suppressWarnings(as.SingleCellExperiment(GSM5599228_Cancer4))
set.seed(100)
sce_GSM5599228_Cancer4 <- scDblFinder(sce_GSM5599228_Cancer4)
GSM5599228_Cancer4$doublet_score <- colData(sce_GSM5599228_Cancer4)$scDblFinder.score
GSM5599228_Cancer4$doublet_class <- colData(sce_GSM5599228_Cancer4)$scDblFinder.class

GSM5599228_Cancer4_dbl1 <- ggplot(GSM5599228_Cancer4@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599228_Cancer4_dbl2 <- VlnPlot(GSM5599228_Cancer4, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599228_Cancer4_doublets <- GSM5599228_Cancer4_dbl1 + GSM5599228_Cancer4_dbl2
GSM5599228_Cancer4_doublets
ggsave(file.path(fig_dir, "GSM5599228_Cancer4_doublets.png"), GSM5599228_Cancer4_doublets, width = 12, height = 5, dpi = 300)
# Sample 10: GSM5599227_Cancer3
sce_GSM5599227_Cancer3 <- suppressWarnings(as.SingleCellExperiment(GSM5599227_Cancer3))
set.seed(100)
sce_GSM5599227_Cancer3 <- scDblFinder(sce_GSM5599227_Cancer3)
GSM5599227_Cancer3$doublet_score <- colData(sce_GSM5599227_Cancer3)$scDblFinder.score
GSM5599227_Cancer3$doublet_class <- colData(sce_GSM5599227_Cancer3)$scDblFinder.class

GSM5599227_Cancer3_dbl1 <- ggplot(GSM5599227_Cancer3@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599227_Cancer3_dbl2 <- VlnPlot(GSM5599227_Cancer3, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599227_Cancer3_doublets <- GSM5599227_Cancer3_dbl1 + GSM5599227_Cancer3_dbl2
GSM5599227_Cancer3_doublets
ggsave(file.path(fig_dir, "GSM5599227_Cancer3_doublets.png"), GSM5599227_Cancer3_doublets, width = 12, height = 5, dpi = 300)
# Sample 11: GSM5599222_Norm3
sce_GSM5599222_Norm3 <- suppressWarnings(as.SingleCellExperiment(GSM5599222_Norm3))
set.seed(100)
sce_GSM5599222_Norm3 <- scDblFinder(sce_GSM5599222_Norm3)
GSM5599222_Norm3$doublet_score <- colData(sce_GSM5599222_Norm3)$scDblFinder.score
GSM5599222_Norm3$doublet_class <- colData(sce_GSM5599222_Norm3)$scDblFinder.class

GSM5599222_Norm3_dbl1 <- ggplot(GSM5599222_Norm3@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599222_Norm3_dbl2 <- VlnPlot(GSM5599222_Norm3, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599222_Norm3_doublets <- GSM5599222_Norm3_dbl1 + GSM5599222_Norm3_dbl2
GSM5599222_Norm3_doublets
ggsave(file.path(fig_dir, "GSM5599222_Norm3_doublets.png"), GSM5599222_Norm3_doublets, width = 12, height = 5, dpi = 300)
# Sample 12: GSM5599221_Norm2
sce_GSM5599221_Norm2 <- suppressWarnings(as.SingleCellExperiment(GSM5599221_Norm2))
set.seed(100)
sce_GSM5599221_Norm2 <- scDblFinder(sce_GSM5599221_Norm2)
GSM5599221_Norm2$doublet_score <- colData(sce_GSM5599221_Norm2)$scDblFinder.score
GSM5599221_Norm2$doublet_class <- colData(sce_GSM5599221_Norm2)$scDblFinder.class

GSM5599221_Norm2_dbl1 <- ggplot(GSM5599221_Norm2@meta.data, aes(x = nCount_RNA, y = doublet_score, color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) + labs(x = "nCount_RNA", y = "Doublet Score") + theme_classic()
GSM5599221_Norm2_dbl2 <- VlnPlot(GSM5599221_Norm2, features = "doublet_score", group.by = "doublet_class", layer = "counts")
GSM5599221_Norm2_doublets <- GSM5599221_Norm2_dbl1 + GSM5599221_Norm2_dbl2
GSM5599221_Norm2_doublets
ggsave(file.path(fig_dir, "GSM5599221_Norm2_doublets.png"), GSM5599221_Norm2_doublets, width = 12, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 4.1 Record QC Summary Table & Cell Counts (Run BEFORE Section 4.1 subsetting)
# ------------------------------------------------------------------------------
library(dplyr)

sample_names <- c(
  "GSM5599220_Norm1", "GSM5599221_Norm2", "GSM5599222_Norm3", 
  "GSM5599223_Norm4", "GSM5599224_Norm5", "GSM5599225_Cancer1", 
  "GSM5599226_Cancer2", "GSM5599227_Cancer3", "GSM5599228_Cancer4", 
  "GSM5599229_Cancer5", "GSM5599230_Cancer6", "GSM5599231_Cancer7"
)

# Fetch all objects directly from memory matching sample_names
samples_list <- mget(sample_names, envir = .GlobalEnv)

qc_summary_table <- lapply(names(samples_list), function(s_name) {
  obj <- samples_list[[s_name]]
  
  data.frame(
    Sample          = s_name,
    Initial_Cells   = ncol(obj),
    Doublets_Found  = sum(obj$doublet_class == "doublet", na.rm = TRUE),
    Final_Singlets  = sum(obj$doublet_class == "singlet", na.rm = TRUE),
    Median_Genes    = median(obj$nFeature_RNA, na.rm = TRUE),
    Median_UMIs     = median(obj$nCount_RNA, na.rm = TRUE),
    Mean_Mito_Pct   = round(mean(obj$percent.mt, na.rm = TRUE), 2)
  )
}) %>% 
  bind_rows() %>%
  mutate(
    Cells_Removed      = Initial_Cells - Final_Singlets,
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
# Sample 4: GSM5599224_Norm5
GSM5599224_Norm5 <- subset(GSM5599224_Norm5, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599224_Norm5)
# Sample 5: GSM5599231_Cancer7
GSM5599231_Cancer7 <- subset(GSM5599231_Cancer7, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599231_Cancer7)
# Sample 6: GSM5599223_Norm4
GSM5599223_Norm4 <- subset(GSM5599223_Norm4, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599223_Norm4)
# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5 <- subset(GSM5599229_Cancer5, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599229_Cancer5)
# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6 <- subset(GSM5599230_Cancer6, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599230_Cancer6)
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4 <- subset(GSM5599228_Cancer4, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599228_Cancer4)
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3 <- subset(GSM5599227_Cancer3, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599227_Cancer3)
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3 <- subset(GSM5599222_Norm3, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599222_Norm3)
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2 <- subset(GSM5599221_Norm2, subset = doublet_class == "singlet")
# Check dimensions after doublet removal
dim(GSM5599221_Norm2)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 4.3 Record Post-Doublet Subset Metrics & Verification
# ------------------------------------------------------------------------------
library(dplyr)

# 1. Define sample names matching your environment objects
sample_names <- c(
  "GSM5599220_Norm1", "GSM5599221_Norm2", "GSM5599222_Norm3", 
  "GSM5599223_Norm4", "GSM5599224_Norm5", "GSM5599225_Cancer1", 
  "GSM5599226_Cancer2", "GSM5599227_Cancer3", "GSM5599228_Cancer4", 
  "GSM5599229_Cancer5", "GSM5599230_Cancer6", "GSM5599231_Cancer7"
)

# 2. Automatically fetch all 12 Seurat objects from memory
samples_list <- mget(sample_names, envir = .GlobalEnv)

# 3. Extract confirmed post-subset metrics across all samples
post_subset_summary_table <- lapply(names(samples_list), function(s_name) {
  obj <- samples_list[[s_name]]
  
  data.frame(
    Sample               = s_name,
    Post_Subset_Singlets = ncol(obj),
    Remaining_Doublets   = sum(obj$doublet_class == "doublet", na.rm = TRUE),
    Median_Genes         = median(obj$nFeature_RNA, na.rm = TRUE),
    Median_UMIs          = median(obj$nCount_RNA, na.rm = TRUE),
    Mean_Mito_Pct        = round(mean(obj$percent.mt, na.rm = TRUE), 2)
  )
}) %>% 
  bind_rows()

# 4. Print verification table to console
print(post_subset_summary_table)

# 5. Ensure output directory exists and export to CSV
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
# Sample 4: GSM5599224_Norm5
GSM5599224_Norm5_mt <- grep("^MT-", rownames(GSM5599224_Norm5), value = TRUE)
GSM5599224_Norm5    <- subset(GSM5599224_Norm5, features = setdiff(rownames(GSM5599224_Norm5), GSM5599224_Norm5_mt))
# Sample 5: GSM5599231_Cancer7
GSM5599231_Cancer7_mt <- grep("^MT-", rownames(GSM5599231_Cancer7), value = TRUE)
GSM5599231_Cancer7    <- subset(GSM5599231_Cancer7, features = setdiff(rownames(GSM5599231_Cancer7), GSM5599231_Cancer7_mt))
# Sample 6: GSM5599223_Norm4
GSM5599223_Norm4_mt <- grep("^MT-", rownames(GSM5599223_Norm4), value = TRUE)
GSM5599223_Norm4    <- subset(GSM5599223_Norm4, features = setdiff(rownames(GSM5599223_Norm4), GSM5599223_Norm4_mt))
# Sample 7: GSM5599229_Cancer5
GSM5599229_Cancer5_mt <- grep("^MT-", rownames(GSM5599229_Cancer5), value = TRUE)
GSM5599229_Cancer5    <- subset(GSM5599229_Cancer5, features = setdiff(rownames(GSM5599229_Cancer5), GSM5599229_Cancer5_mt))
# Sample 8: GSM5599230_Cancer6
GSM5599230_Cancer6_mt <- grep("^MT-", rownames(GSM5599230_Cancer6), value = TRUE)
GSM5599230_Cancer6    <- subset(GSM5599230_Cancer6, features = setdiff(rownames(GSM5599230_Cancer6), GSM5599230_Cancer6_mt))
# Sample 9: GSM5599228_Cancer4
GSM5599228_Cancer4_mt <- grep("^MT-", rownames(GSM5599228_Cancer4), value = TRUE)
GSM5599228_Cancer4    <- subset(GSM5599228_Cancer4, features = setdiff(rownames(GSM5599228_Cancer4), GSM5599228_Cancer4_mt))
# Sample 10: GSM5599227_Cancer3
GSM5599227_Cancer3_mt <- grep("^MT-", rownames(GSM5599227_Cancer3), value = TRUE)
GSM5599227_Cancer3    <- subset(GSM5599227_Cancer3, features = setdiff(rownames(GSM5599227_Cancer3), GSM5599227_Cancer3_mt))
# Sample 11: GSM5599222_Norm3
GSM5599222_Norm3_mt <- grep("^MT-", rownames(GSM5599222_Norm3), value = TRUE)
GSM5599222_Norm3    <- subset(GSM5599222_Norm3, features = setdiff(rownames(GSM5599222_Norm3), GSM5599222_Norm3_mt))
# Sample 12: GSM5599221_Norm2
GSM5599221_Norm2_mt <- grep("^MT-", rownames(GSM5599221_Norm2), value = TRUE)
GSM5599221_Norm2    <- subset(GSM5599221_Norm2, features = setdiff(rownames(GSM5599221_Norm2), GSM5599221_Norm2_mt))
# Check gene dimensions after MT removal
# 1. Define sample names matching your environment objects
sample_names <- c(
  "GSM5599220_Norm1", "GSM5599221_Norm2", "GSM5599222_Norm3", 
  "GSM5599223_Norm4", "GSM5599224_Norm5", "GSM5599225_Cancer1", 
  "GSM5599226_Cancer2", "GSM5599227_Cancer3", "GSM5599228_Cancer4", 
  "GSM5599229_Cancer5", "GSM5599230_Cancer6", "GSM5599231_Cancer7"
)

# 2. Fetch all 12 Seurat objects from memory
samples_list <- mget(sample_names, envir = .GlobalEnv)

# 3. Print gene dimensions after MT removal for all samples
cat("Genes remaining per sample after MT removal:\n")
gene_counts <- data.frame(
  Sample = names(samples_list),
  Genes_Remaining = sapply(samples_list, nrow)
)

invisible(lapply(1:nrow(gene_counts), function(i) {
  cat(sprintf("%-20s %d genes\n", paste0(gene_counts$Sample[i], ":"), gene_counts$Genes_Remaining[i]))
}))

# 4. Save to CSV file
write.csv(gene_counts, file = file.path(out_dir, "genes_remaining_after_mt_removal.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5.1 Export Cell Counts Before and After QC
# ------------------------------------------------------------------------------
library(dplyr)

# Create results directory if needed
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# 1. Define sample names matching your environment objects
sample_names <- c(
  "GSM5599220_Norm1", "GSM5599221_Norm2", "GSM5599222_Norm3", 
  "GSM5599223_Norm4", "GSM5599224_Norm5", "GSM5599225_Cancer1", 
  "GSM5599226_Cancer2", "GSM5599227_Cancer3", "GSM5599228_Cancer4", 
  "GSM5599229_Cancer5", "GSM5599230_Cancer6", "GSM5599231_Cancer7"
)

# 2. Fetch all 12 Seurat objects from memory
samples_list <- mget(sample_names, envir = .GlobalEnv)

# 3. If you have an initial counts lookup vector, define it here:
initial_counts <- c(
  "GSM5599220_Norm1" = 5831, "GSM5599221_Norm2" = 6000, "GSM5599222_Norm3" = 6000,
  "GSM5599223_Norm4" = 6000, "GSM5599224_Norm5" = 6000, "GSM5599225_Cancer1" = 6977,
  "GSM5599226_Cancer2" = 3473, "GSM5599227_Cancer3" = 6000, "GSM5599228_Cancer4" = 6000,
  "GSM5599229_Cancer5" = 6000, "GSM5599230_Cancer6" = 6000, "GSM5599231_Cancer7" = 6000
)

# 4. Generate post-QC cell counts summary across individual objects
qc_cell_counts <- data.frame(
  Sample = names(samples_list),
  Cells_After_QC = sapply(samples_list, ncol)
) %>%
  mutate(
    Cells_Before_QC    = initial_counts[Sample],
    Cells_Removed      = Cells_Before_QC - Cells_After_QC,
    Retention_Rate_Pct = round((Cells_After_QC / Cells_Before_QC) * 100, 2)
  ) %>%
  select(Sample, Cells_Before_QC, Cells_After_QC, Cells_Removed, Retention_Rate_Pct)

# View and export
print(qc_cell_counts)
write.csv(qc_cell_counts, file = "results/qc_cell_counts_per_sample.csv", row.names = FALSE)

# 6. Merge Samples & Save RDS Output
# ------------------------------------------------------------------------------
# Ensure output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Define sample names matching your environment objects
sample_names <- c(
  "GSM5599220_Norm1", "GSM5599221_Norm2", "GSM5599222_Norm3", 
  "GSM5599223_Norm4", "GSM5599224_Norm5", "GSM5599225_Cancer1", 
  "GSM5599226_Cancer2", "GSM5599227_Cancer3", "GSM5599228_Cancer4", 
  "GSM5599229_Cancer5", "GSM5599230_Cancer6", "GSM5599231_Cancer7"
)

# 2. Automatically retrieve all 12 objects from memory
samples_list <- mget(sample_names, envir = .GlobalEnv)

# 3. Merge all 12 Seurat objects into a single combined object
merged_seurat <- merge(
  x = samples_list[[1]],
  y = samples_list[-1],
  add.cell.ids = names(samples_list)
)

# 4. Export final merged post-QC Seurat object to RDS
rds_path <- file.path(out_dir, "postQC_merged_12_samples.rds")
saveRDS(merged_seurat, file = rds_path)

# 5. Console print status update
cat(sprintf(
  "\nSuccessfully completed Stage 1 QC pipeline!\nMerged object saved to: %s\nTotal dataset contains %d genes across %d cells.\n",
  rds_path,
  dim(merged_seurat)[1],
  dim(merged_seurat)[2]
))
# ------------------------------------------------------------------------------
# 6.1 Inspect Merged Seurat Object & Export Results
# ------------------------------------------------------------------------------

# 1. Check overall object dimensions (Genes x Cells)
dim(merged_seurat)

# 2. Inspect active assay and Seurat v5 layer structure
merged_seurat[["RNA"]]
Layers(merged_seurat)

# 3. View metadata structure and available column names
head(merged_seurat@meta.data, 5)
colnames(merged_seurat@meta.data)

# 4. Check total cell counts per sample across all 12 samples
table(merged_seurat$orig.ident)

# 5. Summarize key QC metrics across samples
aggregate(
  cbind(nCount_RNA, nFeature_RNA) ~ orig.ident, 
  data = merged_seurat@meta.data, 
  FUN = median
)

# 6. Create results directory if it does not exist
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# 7. Save the full 12-sample merged Seurat object locally
saveRDS(merged_seurat, file = "results/postQC_merged_12_samples.rds")

# 8. Save session info to a text file for reproducible research
writeLines(
  capture.output(sessionInfo()), 
  con = file.path("results", "session_info.txt")
)

# 9. Open working directory in file manager
if (.Platform$OS.type == "windows") {
  shell.exec(getwd())
} else if (Sys.info()["sysname"] == "Darwin") {
  system("open .")
} else {
  system("xdg-open .")
}

# 10. Update .gitignore to exclude large RDS files
if (!file.exists(".gitignore")) file.create(".gitignore")
write("*.rds", file = ".gitignore", append = TRUE)
write("results/*.rds", file = ".gitignore", append = TRUE)


# Define output directory and file path
env_dir <- "saved_environments"
if (!dir.exists(env_dir)) dir.create(env_dir, recursive = TRUE)

# Save the entire workspace environment
save.image(file = file.path(env_dir, "postQC_workspace.RData"))

# Update .gitignore to exclude .RData files (they are usually too large for GitHub)
if (!file.exists(".gitignore")) file.create(".gitignore")
write("*.RData", file = ".gitignore", append = TRUE)

cat("Environment successfully saved! You can safely close RStudio now.\n")

# Define directory
env_dir <- "saved_objects"
if (!dir.exists(env_dir)) dir.create(env_dir, recursive = TRUE)

# Save only key QC objects to a single file
save(
  merged_seurat, 
  qc_cell_counts, 
  post_subset_summary_table, 
  file = file.path(env_dir, "postQC_essential_objects.RData")
)

# Ensure .RData files are excluded from Git tracking
if (!file.exists(".gitignore")) file.create(".gitignore")
write("*.RData", file = ".gitignore", append = TRUE)

cat("Selected objects successfully saved!\n")