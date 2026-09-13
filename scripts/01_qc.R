# ==============================================================================
# 01_qc.R
# Paper thresholds (Xu et al. 2022, Materials and Methods):
#   - Keep cells with > 200 expressed genes/cell
#   - Keep cells with mitochondrial UMI rate < 40%
#   - Mitochondrial genes removed from the expression table 
#     after filtering (not used as HVGs/PCA input)
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# >>> EDIT THIS <<< Random seed — keep this the SAME value
# >>> across every script in the pipeline (01_qc.R, 02_..., etc.)
# ------------------------------------------------------------------------------
set.seed(42)

# ==============================================================================
# >>>>>>>>>>>>>>>>> USER CONFIGURATION — EDIT HERE <<<<<<<<<<<<<<<<<
# ==============================================================================

# Threshold parameters matching Xu et al. 2022 criteria
MIN_GENES_PER_CELL <- 500
MAX_MITO_PCT       <- 20

# >>> EDIT THIS <<<
# Main folder containing subfolders for each sample (each with matrix.mtx, barcodes.tsv, features.tsv)
data_dir <- "GSE184880"

# >>> EDIT THIS <<<
# Uncomment ONLY your assigned block from the QC Task distribution:

# --- HIBA (1 Normal + 2 Tumor) ---
# my_samples <- c("GSM5599220", "GSM5599225", "GSM5599226")

# --- TAGREED (1 Normal + 1 Tumor) ---
# my_samples <- c("GSM5599221", "GSM5599227")

# --- LEQAA (1 Normal + 1 Tumor) ---
# my_samples <- c("GSM5599222", "GSM5599228")

# --- MARIYAM (1 Normal + 2 Tumor) ---
# my_samples <- c("GSM5599223", "GSM5599229", "GSM5599230")

# --- NOUR (1 Normal + 1 Tumor) ---
 my_samples <- c("GSM5599224", "GSM5599231")

# Or set my_samples <- NULL to process all subfolders found inside data_dir
# my_samples <- NULL

# >>> EDIT THIS <<<
# Output directory definitions
out_dir <- "results"
fig_dir <- "figures"

# ==============================================================================
# >>>>>>>>>>>>>>>>> END OF USER CONFIGURATION <<<<<<<<<<<<<<<<<<<<<
# ==============================================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Determine sample folders to process
if (is.null(my_samples)) {
  sample_ids <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
} else {
  sample_ids <- my_samples
}

# ------------------------------------------------------------------------------
# 1. Load each sample and build a per-sample QC table
# ------------------------------------------------------------------------------
seurat_list <- list()
qc_summary <- data.frame(
  sample = character(),
  cells_before = numeric(),
  cells_after = numeric(),
  stringsAsFactors = FALSE
)

for (sid in sample_ids) {
  sample_path <- file.path(data_dir, sid)
  
  message(sprintf("Processing sample: %s", sid))
  
  # Load raw 10x matrix data
  counts <- Seurat::ReadMtx(
    mtx = "C:/Users/Lenovo/Downloads/GSE184880/GSM5599224/matrix.mtx",
    features = "C:/Users/Lenovo/Downloads/GSE184880/GSM5599224/features.tsv.tsv",
    cells = "C:/Users/Lenovo/Downloads/GSE184880/GSM5599224/barcodes.tsv.tsv",
    feature.column = 2
  )
  # Initialize Seurat Object with basic gene expression threshold
  so <- CreateSeuratObject(
    counts = counts,
    project = sid,
    min.cells = 3,
    min.features = MIN_GENES_PER_CELL
  )
  
  cells_before <- ncol(so)
  
  # Calculate percent mitochondrial reads per cell
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  
  # Apply paper's two cell-level filters (> 200 genes & < 40% mito)
  so <- subset(
    so,
    subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT
  )
  
  cells_after <- ncol(so)
  
  # Record QC log metrics
  qc_summary <- rbind(qc_summary, data.frame(    sample = sid,
    cells_before = cells_before,
    cells_after = cells_after
  ))
  
  seurat_list[[sid]] <- so
}

# Display and write out individual QC cell counts table
print(qc_summary)
write.csv(
  qc_summary,
  file.path(out_dir, paste0("qc_cell_counts_", paste(sample_ids, collapse = "-"), ".csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 2. Remove mitochondrial genes from the expression matrix
#    (paper: mito genes removed completely from matrix after QC filtering)
# ------------------------------------------------------------------------------
seurat_list <- lapply(seurat_list, function(so) {
  mt_genes   <- grep("^MT-", rownames(so), value = TRUE)
  keep_genes <- setdiff(rownames(so), mt_genes)
  subset(so, features = keep_genes)
})

# ------------------------------------------------------------------------------
# 3. Merge this person's samples into one object
#    (A teammate will merge everyone's outputs in Stage 3)
# ------------------------------------------------------------------------------
if (length(seurat_list) == 1) {
  merged <- seurat_list[[1]]
} else {
  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_ids
  )
}

# Export individual merged post-QC Seurat Object
saveRDS(merged, file.path(out_dir, paste0("postQC_", paste(sample_ids, collapse = "-"), ".rds")))

# ------------------------------------------------------------------------------
# 4. QC diagnostic plots
# ------------------------------------------------------------------------------
p1 <- VlnPlot(merged, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") +
  geom_hline(yintercept = MIN_GENES_PER_CELL, linetype = "dashed", color = "red") +
  ggtitle("Genes per cell (threshold = 200)")

p2 <- VlnPlot(merged, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  geom_hline(yintercept = MAX_MITO_PCT, linetype = "dashed", color = "red") +
  ggtitle("Mitochondrial % (threshold = 40%)")

p3 <- p1 + p2
ggsave(
  file.path(fig_dir, paste0("qc_violin_", paste(sample_ids, collapse = "-"), ".png")),
  p3, width = 12, height = 5, dpi = 300
)

# ------------------------------------------------------------------------------
# 5. Report totals against the paper's published numbers
# ------------------------------------------------------------------------------
cat(sprintf(
  "\nCells retained after QC for this subset (%s): %d\n",
  paste(sample_ids, collapse = ", "), ncol(merged)
))

# Output environment session info to guarantee reproducible tracking
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, paste0("sessionInfo_qc_", paste(sample_ids, collapse = "-"), ".txt"))
)