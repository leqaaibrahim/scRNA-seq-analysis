# ==============================================================================
# 00. PROJECT SETUP AND PACKAGE CHECK
# ==============================================================================
# Run from the repository root. This script creates the expected folders and
# checks dependencies; it does not install packages or alter analytical output.

set.seed(42)

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "ggrepel", "patchwork", "SingleCellExperiment",
  "scDblFinder", "batchelor", "BiocSingular", "pheatmap", "decoupleR",
  "dorothea", "limma", "CellChat"
)

optional_packages <- c("ragg", "ggtext")

for (directory in c("data", "figures", "results")) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

missing_required <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
missing_optional <- optional_packages[
  !vapply(optional_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required)) {
  stop(
    "Missing required packages: ", paste(missing_required, collapse = ", "),
    ". Install them before running Stage 01. See README.md."
  )
}

if (length(missing_optional)) {
  message(
    "Optional packages not installed: ",
    paste(missing_optional, collapse = ", "),
    ". The scripts retain documented fallback plotting behaviour."
  )
}

writeLines(capture.output(sessionInfo()), "results/00_session_info.txt")
message("Project folders and required R packages are ready.")
