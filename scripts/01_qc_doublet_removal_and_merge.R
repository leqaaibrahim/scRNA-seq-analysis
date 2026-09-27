# ==============================================================================
# 01. QUALITY CONTROL, DOUBLET REMOVAL, AND SAMPLE MERGING
#
# Input:  10x matrices for 12 GSE184880 samples
# Output: results/postQC_merged_12_samples.rds and QC audit tables/figures
# Design: apply identical QC and scDblFinder processing independently per sample
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# 0. Configuration -------------------------------------------------------------

RANDOM_SEED <- 42L
DOUBLET_SEED <- 100L
MIN_GENES_PER_CELL <- 200L
MAX_MITO_PCT <- 40

# Threshold rationale: remove very low-complexity/empty-droplet profiles while
# retaining stressed tumour-tissue cells, which can have elevated mitochondrial
# RNA. Before presentation, confirm this choice from the exported per-sample QC
# distributions and report both starting and retained cell counts.

# Keep data outside GitHub. Before running, either place the sample folders under
# data/GSE184880 or define SCRNA_DATA_DIR with the full local data directory.
DATA_DIR <- Sys.getenv("SCRNA_DATA_DIR", unset = file.path("data", "GSE184880"))
RESULT_DIR <- "results"
FIGURE_DIR <- "figures"

SAMPLE_IDS <- c(
  "GSM5599220_Norm1", "GSM5599221_Norm2", "GSM5599222_Norm3",
  "GSM5599223_Norm4", "GSM5599224_Norm5", "GSM5599225_Cancer1",
  "GSM5599226_Cancer2", "GSM5599227_Cancer3", "GSM5599228_Cancer4",
  "GSM5599229_Cancer5", "GSM5599230_Cancer6", "GSM5599231_Cancer7"
)

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(RANDOM_SEED)

required_packages <- c("Seurat", "SingleCellExperiment", "scDblFinder",
                       "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install required packages before running: ",
       paste(missing_packages, collapse = ", "))
}

# 1. Helpers -------------------------------------------------------------------

fix_geo_filenames <- function(sample_path) {
  if (!dir.exists(sample_path)) stop("Missing sample directory: ", sample_path)
  files <- list.files(sample_path, full.names = TRUE)

  rename_match <- function(pattern, allowed, target_stub) {
    hits <- files[grepl(pattern, basename(files)) & !basename(files) %in% allowed]
    for (old_file in hits) {
      target <- paste0(target_stub, ifelse(grepl("\\.gz$", old_file), ".gz", ""))
      file.rename(old_file, file.path(sample_path, target))
    }
  }

  rename_match("matrix\\.mtx", c("matrix.mtx", "matrix.mtx.gz"), "matrix.mtx")
  rename_match("barcodes\\.tsv", c("barcodes.tsv", "barcodes.tsv.gz"), "barcodes.tsv")
  rename_match("features\\.tsv|genes\\.tsv",
               c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"),
               "features.tsv")
}

save_plot <- function(plot, filename, width, height) {
  output <- file.path(FIGURE_DIR, filename)
  ggsave(output, plot = plot, width = width, height = height,
         dpi = 300, bg = "white")
  invisible(output)
}

qc_violin <- function(object) {
  VlnPlot(object,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3, layer = "counts")
}

qc_scatter <- function(object) {
  FeatureScatter(object, feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA", slot = "counts") +
    FeatureScatter(object, feature1 = "nCount_RNA",
                   feature2 = "percent.mt", slot = "counts")
}

detect_doublets <- function(object, sample_id) {
  set.seed(DOUBLET_SEED)
  sce <- suppressWarnings(as.SingleCellExperiment(object))
  sce <- scDblFinder(sce)
  object$doublet_score <- colData(sce)$scDblFinder.score
  object$doublet_class <- colData(sce)$scDblFinder.class

  score_plot <- ggplot(
    object[[]], aes(nCount_RNA, doublet_score, color = doublet_class)
  ) +
    geom_point(size = 1, alpha = 0.6) +
    labs(x = "nCount_RNA", y = "Doublet score", color = "Class") +
    theme_classic()

  violin_plot <- VlnPlot(
    object, features = "doublet_score", group.by = "doublet_class",
    layer = "counts"
  )
  save_plot(score_plot + violin_plot,
            paste0(sample_id, "_doublets.png"), 12, 5)
  object
}

# 2. Load and quality-filter each sample ---------------------------------------

samples <- setNames(vector("list", length(SAMPLE_IDS)), SAMPLE_IDS)
initial_counts <- setNames(integer(length(SAMPLE_IDS)), SAMPLE_IDS)

for (sample_id in SAMPLE_IDS) {
  message("Loading and processing: ", sample_id)
  sample_path <- file.path(DATA_DIR, sample_id)
  fix_geo_filenames(sample_path)

  counts <- Read10X(sample_path)
  # Record the downloaded matrix size before any Seurat filtering.
  initial_counts[sample_id] <- ncol(counts)
  object <- CreateSeuratObject(
    counts = counts, project = sample_id,
    min.cells = 3, min.features = MIN_GENES_PER_CELL
  )
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

  save_plot(qc_violin(object), paste0(sample_id, "_vln_before.png"), 10, 4)
  save_plot(qc_scatter(object), paste0(sample_id, "_scatter.png"), 10, 4)

  object <- subset(
    object,
    subset = nFeature_RNA > MIN_GENES_PER_CELL & percent.mt < MAX_MITO_PCT
  )
  save_plot(qc_violin(object), paste0(sample_id, "_vln_after.png"), 10, 4)
  samples[[sample_id]] <- detect_doublets(object, sample_id)
}

# 3. Audit and remove predicted doublets ---------------------------------------

qc_summary_table <- bind_rows(lapply(names(samples), function(sample_id) {
  object <- samples[[sample_id]]
  data.frame(
    Sample = sample_id,
    Initial_Cells = ncol(object),
    Doublets_Found = sum(object$doublet_class == "doublet", na.rm = TRUE),
    Final_Singlets = sum(object$doublet_class == "singlet", na.rm = TRUE),
    Median_Genes = median(object$nFeature_RNA, na.rm = TRUE),
    Median_UMIs = median(object$nCount_RNA, na.rm = TRUE),
    Mean_Mito_Pct = round(mean(object$percent.mt, na.rm = TRUE), 2)
  )
})) %>%
  mutate(
    Cells_Removed = Initial_Cells - Final_Singlets,
    Retention_Rate_Pct = round(100 * Final_Singlets / Initial_Cells, 2)
  )

write.csv(qc_summary_table,
          file.path(RESULT_DIR, "qc_cell_counts_summary.csv"), row.names = FALSE)

samples <- lapply(samples, function(object) {
  subset(object, subset = doublet_class == "singlet")
})

post_subset_summary_table <- bind_rows(lapply(names(samples), function(sample_id) {
  object <- samples[[sample_id]]
  data.frame(
    Sample = sample_id,
    Post_Subset_Singlets = ncol(object),
    Remaining_Doublets = sum(object$doublet_class == "doublet", na.rm = TRUE),
    Median_Genes = median(object$nFeature_RNA, na.rm = TRUE),
    Median_UMIs = median(object$nCount_RNA, na.rm = TRUE),
    Mean_Mito_Pct = round(mean(object$percent.mt, na.rm = TRUE), 2)
  )
}))
write.csv(post_subset_summary_table,
          file.path(RESULT_DIR, "post_doublet_subset_summary.csv"), row.names = FALSE)

# Remove mitochondrial genes only after their percentages have been recorded.
samples <- lapply(samples, function(object) {
  mt_genes <- grep("^MT-", rownames(object), value = TRUE)
  subset(object, features = setdiff(rownames(object), mt_genes))
})

gene_counts <- data.frame(
  Sample = names(samples),
  Genes_Remaining = vapply(samples, nrow, integer(1))
)
write.csv(gene_counts,
          file.path(RESULT_DIR, "genes_remaining_after_mt_removal.csv"),
          row.names = FALSE)

qc_cell_counts <- data.frame(
  Sample = names(samples),
  Cells_Before_QC = unname(initial_counts[names(samples)]),
  Cells_After_QC = vapply(samples, ncol, integer(1))
) %>%
  mutate(
    Cells_Removed = Cells_Before_QC - Cells_After_QC,
    Retention_Rate_Pct = round(100 * Cells_After_QC / Cells_Before_QC, 2)
  )
write.csv(qc_cell_counts,
          file.path(RESULT_DIR, "qc_cell_counts_per_sample.csv"), row.names = FALSE)

# 4. Merge, save, and record reproducibility ----------------------------------

merged_seurat <- merge(
  x = samples[[1]], y = samples[-1], add.cell.ids = names(samples)
)
saveRDS(merged_seurat, file.path(RESULT_DIR, "postQC_merged_12_samples.rds"))

sample_qc_medians <- aggregate(
  cbind(nCount_RNA, nFeature_RNA) ~ orig.ident,
  data = merged_seurat[[]], median
)
write.csv(sample_qc_medians,
          file.path(RESULT_DIR, "merged_sample_qc_medians.csv"), row.names = FALSE)

parameters <- data.frame(
  parameter = c("random_seed", "doublet_seed", "min_genes_per_cell",
                "max_mito_percent", "min_cells_per_gene"),
  value = c(RANDOM_SEED, DOUBLET_SEED, MIN_GENES_PER_CELL,
            MAX_MITO_PCT, 3)
)
write.csv(parameters, file.path(RESULT_DIR, "01_qc_parameters.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(RESULT_DIR, "01_session_info.txt"))

message("Stage 01 complete: ", nrow(merged_seurat), " genes and ",
        ncol(merged_seurat), " cells saved to results/postQC_merged_12_samples.rds")
