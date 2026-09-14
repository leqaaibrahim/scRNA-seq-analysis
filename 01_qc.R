# ==============================================================================
# 02_qc.R
# Pipeline: Per-Sample QC -> Doublet Detection -> MT Removal -> Merge
# ==============================================================================
## Quality control for GSE184880 (Xu et al. 2022, Clin Cancer Res)
## Filters replicate the ORIGINAL PAPER's thresholds, not our
## proposal's stricter draft thresholds. This is the deliberate
## "reproduce-first" step of the project (see README for rationale).
##
## Paper thresholds (Xu et al. 2022, Materials and Methods):
##   - Keep cells with > 200 expressed genes/cell
##   - Keep cells with mitochondrial UMI rate < 40%
##   - Mitochondrial genes removed from the expression table
##     after filtering (not used as HVGs/PCA input)
## ==============================================================
## >>>>>>>>>>>>>>>>>  1. Load Packages  <<<<<<<<<<<<<<<<<
## ==============================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
set.seed(42)

## ==============================================================
## >>>>>>>>>>>>>>>>>  USER CONFIGURATION  <<<<<<<<<<<<<<<<<
## ==============================================================

data_dir <- "data/GSE184880"
out_dir <- "results"
fig_dir <- "figures"

# Example for one team member's assigned samples:
## GSM5599221 (Normal_ovary_2), GSM5599227 (Cancer_HGSOC3)
my_samples <- c("GSM5599221", "GSM5599227")

# the paper — see README for why we default to the PAPER's values)
MIN_GENES_PER_CELL <- 200   # Xu et al. 2022: "> 200 expressed genes"
MAX_MITO_PCT       <- 40    # Xu et al. 2022: "mitochondrial UMI rate below 40%"

## ==============================================================
## >>>>>>>>>>>>>>>>>  END OF USER CONFIGURATION  <<<<<<<<<<<<<<<<<<<<<
## ==============================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Resolve which sample folders to load
## ------------------------------------------------------------
all_sample_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
all_sample_ids  <- basename(all_sample_dirs)

if (is.null(my_samples)) {
  sample_dirs <- all_sample_dirs
  sample_ids  <- all_sample_ids
} else {
  keep        <- all_sample_ids %in% my_samples
  sample_dirs <- all_sample_dirs[keep]
  sample_ids  <- all_sample_ids[keep]
  missing     <- setdiff(my_samples, all_sample_ids)
  if (length(missing) > 0) {
    warning("These sample folders were not found in data_dir: ", paste(missing, collapse = ", "))
  }
}

stopifnot("No sample folders matched — check data_dir and my_samples." = length(sample_dirs) > 0)
message("Processing ", length(sample_dirs), " sample(s): ", paste(sample_ids, collapse = ", "))

## ------------------------------------------------------------
## 2. Load each sample and build a per-sample QC table
## ------------------------------------------------------------
library(Matrix)

## Create a function to read GEO Sample
read_geo_sample <- function(sample, dir) {
  
  barcodes_file <- list.files(path = dir, pattern = "\\.barcodes.tsv.gz", full.names = TRUE)
  
  genes_file <- list.files(path = dir, pattern = "\\.genes.tsv.gz", full.names = TRUE)
  
  matrix_file <- list.files(path = dir, pattern = "\\.matrix.mtx.gz", full.names = TRUE)
  
  
  # Read matrix
  mat <- readMM(matrix_file)
  
  # Read genes
  genes <- read.delim(genes_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Read barcodes
  barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Assign row/column names
  rownames(mat) <- make.unique(genes[, 2])
  colnames(mat) <- barcodes[, 1]
  
  # Convert to sparse matrix
  mat <- as(mat, "CsparseMatrix")
  
  return(mat)
}
##-----------------------------------------------------------------------------##
seurat_list <- list()
qc_summary  <- data.frame()

for (i in seq_along(sample_dirs)) {

  sid <- sample_ids[i]
  message("Loading sample: ", sid)

  counts <- read_geo_sample(sample = sid, dir = sample_dirs[i]) # Read10X(data.dir = sample_dirs[i])
  so <- CreateSeuratObject(counts = counts, 
                           project = sid, 
                           min.cells = 3)    # standard sparse-gene filter, not a cell-level QC step
  

  # Percent mitochondrial reads per cell (before removing MT genes)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

  n_before <- ncol(so)

  # --- Apply the paper's two cell-level filters ---
  so <- subset(so, 
               subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT)

  n_after <- ncol(so)

  qc_summary <- rbind(qc_summary, data.frame(
    sample        = sid,
    cells_before  = n_before,
    cells_after   = n_after,
    pct_retained  = round(100 * n_after / n_before, 1)
  ))

  # >>> EDIT THIS <<<  (optional) tag stage/condition metadata once your
  # sample sheet is finalized, e.g.:
  # so$stage <- sample_metadata[sample_metadata$sample == sid, "stage"]

  seurat_list[[sid]] <- so
}

print(qc_summary)

#        sample cells_before cells_after pct_retained
# 1 GSM5599221         5174        5065         97.9
# 2 GSM5599227         5107        4744         92.9


write.csv(
  qc_summary,
  file.path(out_dir, paste0("qc_cell_counts_", paste(sample_ids, collapse = "-"), ".csv")),
  row.names = FALSE
)

## ------------------------------------------------------------
## 3. Remove mitochondrial genes from the expression matrix
##    (paper: mito genes removed from the table after QC filtering)
## ------------------------------------------------------------
seurat_list <- lapply(seurat_list, function(so) {
  mt_genes   <- grep("^MT-", rownames(so), value = TRUE)
  keep_genes <- setdiff(rownames(so), mt_genes)
  subset(so, features = keep_genes)
})

## ------------------------------------------------------------
## 4. Merge this person's samples into one object
##    (a teammate will merge everyone's outputs in Stage 3)
## ------------------------------------------------------------
if (length(seurat_list) == 1) {
  merged <- seurat_list[[1]]
} else {
  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_ids
  )
}

saveRDS(merged, file.path(out_dir, paste0("postQC_", paste(sample_ids, collapse = "-"), ".rds")))

## ------------------------------------------------------------
## 5. QC diagnostic plots
## ------------------------------------------------------------
p1 <- VlnPlot(merged, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") +
  geom_hline(yintercept = MIN_GENES_PER_CELL, linetype = "dashed", color = "red") +
  ggtitle("Genes per cell (threshold = 200)")

p2 <- VlnPlot(merged, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  geom_hline(yintercept = MAX_MITO_PCT, linetype = "dashed", color = "red") +
  ggtitle("Mitochondrial % (threshold = 40%)")

p3 <- p1 + p2

ggsave(
  file.path(fig_dir, paste0("qc_violin_", paste(sample_ids, collapse = "-"), ".png")),
  p3, width = 12, height = 5, dpi = 300)

## ------------------------------------------------------------
## 6. Report totals against the paper's published numbers
##    (paper: 59,324 cells total; 33,264 HGSOC / 26,060 nonmalignant,
##    across all 12 samples — your subset total will be smaller)
## ------------------------------------------------------------
cat(sprintf(
  "\nCells retained after QC for this subset (%s): %d\n",
  paste(sample_ids, collapse = ", "), ncol(merged)
))

### Cells retained after QC for this subset (GSM5599221, GSM5599227): 9809

## sessionInfo() saved for the record (guide: "record your versions")
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, paste0("sessionInfo_qc_", paste(sample_ids, collapse = "-"), ".txt"))
)
