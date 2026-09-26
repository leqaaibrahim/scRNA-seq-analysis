# ==============================================================================
# 04. SAMPLE-AWARE TF ACTIVITY AND STAGE-STRATIFIED CELLCHAT ANALYSIS
#     GSE184880 / Xu et al. HGSOC reanalysis
# ==============================================================================
# Biological question:
# Which transcription-factor activities and cell-communication programmes
# distinguish evidence-based myeloid/macrophage states across Non-malignant, Early,
# and Advanced HGSOC tissue groups?
#
# Design principles from the project guide and proposal:
#   * use internally validated macrophage annotations; no external atlas labels;
#   * preserve samples/patients as biological replicates;
#   * use non-integrated RNA counts/expression for inference;
#   * infer TF activity with curated DoRothEA A/B/C regulons and ULM;
#   * build consistently processed Non-malignant, Early, and Advanced CellChat
#     models;
#   * use epithelial and CD8 T-cell compartments as the two primary
#     macrophage communication axes;
#   * treat CLEC10A+ myeloid APC-like communication as supplementary;
#   * represent every other eligible macrophage-compartment axis with separate
#     exploratory bubble plots and retain all interactions in screening tables;
#   * report predictions as computational hypotheses, not mechanistic proof.
#
# Required upstream files:
#   results/03_clean_macrophages_Xu_annotated.rds
#     or results/03_clean_macrophages_subclustered.rds
#   results/02.1_integrated_seurat_refined.rds
#     or results/02_integrated_seurat.rds as a fallback
#
# Main output directories:
#   results/04_downstream/
#   figures/04_downstream/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(decoupleR)
  library(CellChat)
  library(limma)
})
install.packages("ggtext")
set.seed(42)

# Package installation, if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("decoupleR", "OmnipathR", "limma",
#                        "ComplexHeatmap"))
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
# remotes::install_github("jinworks/CellChat")

# ------------------------------------------------------------------------------
# 0. Parameters and output helpers
# ------------------------------------------------------------------------------

MACROPHAGE_RDS_CANDIDATES <- c(
  "results/03_clean_macrophages_Xu_annotated.rds",
  "results/03_clean_macrophages_subclustered.rds"
)

GLOBAL_RDS_CANDIDATES <- c(
  "results/02.1_integrated_seurat_refined.rds",
  "results/02_integrated_seurat.rds"
)

RESULT_DIR <- "results/04_downstream"
FIGURE_DIR <- "figures/04_downstream"
MIN_CELLS_PER_SAMPLE_STATE <- 20L
MIN_SAMPLES_PER_STAGE <- 2L
CELLCHAT_MIN_CELLS <- 10L
CELLCHAT_MAX_CELLS_PER_IDENTITY_STAGE <- 500L
BUBBLE_PVALUE_THRESHOLD <- 0.05
BUBBLE_MIN_HEIGHT <- 12
BUBBLE_MAX_HEIGHT <- 24
BUBBLE_WIDTH <- 22
# These pathways are emphasized visually because they are central to the
# macrophage-focused biological question. Emphasis does not change filtering.
PRIORITY_SIGNAL_PATTERNS <- c(
  "CXCL", "CCL", "TGFB", "TGF", "VEGF", "TNF", "IL1", "SPP1",
  "MIF", "GALECTIN"
)
# Display limit only; the complete TF statistics are still written to CSV.
TOP_TFS_HEATMAP_PER_STATE <- 5L
RANDOM_SEED <- 42L

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

first_existing_file <- function(paths, description) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(description, " was not found. Checked: ", paste(paths, collapse = ", "))
  }
  hit[1]
}

save_png <- function(plot, filename, width, height, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width, height = height, units = "in",
                  res = dpi, background = "white")
    print(plot)
    grDevices::dev.off()
  } else {
    ggsave(filename, plot, width = width, height = height, dpi = dpi,
           bg = "white", device = "png", limitsize = FALSE)
  }
}

save_base_png <- function(filename, width, height, code, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width, height = height, units = "in",
                  res = dpi, background = "white")
  } else {
    grDevices::png(filename, width = width, height = height, units = "in",
                   res = dpi, bg = "white")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

# Concise labels used only in figures. Full evidence-based annotations remain
# unchanged in Seurat metadata, CellChat objects and exported CSV files.
DISPLAY_STATE_LABELS <- c(
  "Mac-FOLR2-LYVE1-resident" = "Resident",
  "Mac-Cycling" = "Cycling",
  "Mac-OLFML3-CX3CR1-resident" = "OLFML3 resident",
  "Mac-TNF-CCL4-inflammatory" = "Inflammatory",
  "C7-APOBEC3A-like" = "APOBEC3A",
  "Mono-FCN1-S100A8-S100A9" = "FCN1 monocytes",
  "C6-SLC2A1-like" = "SLC2A1",
  "Mac-NUPR1-IGF2-stress" = "Stress-adapted",
  "Myeloid-CLEC10A-APC" = "CLEC10A APC-like",
  "Mono-CDKN1C" = "CDKN1C monocytes"
)

compact_cell_labels <- function(x) {
  output <- as.character(x)
  for (full_label in names(DISPLAY_STATE_LABELS)) {
    output <- gsub(
      full_label,
      DISPLAY_STATE_LABELS[[full_label]],
      output,
      fixed = TRUE
    )
  }
  output
}

wrap_plot_labels <- function(x, width = 24L) {
  vapply(
    as.character(x),
    function(label) paste(strwrap(label, width = width), collapse = "\n"),
    character(1)
  )
}

compact_wrap_labels <- function(x, width = 24L) {
  wrap_plot_labels(compact_cell_labels(x), width = width)
}

get_rna_layer <- function(object, layer = c("counts", "data")) {
  layer <- match.arg(layer)
  DefaultAssay(object) <- "RNA"
  available_layers <- Layers(object[["RNA"]])

  if (layer %in% available_layers) {
    return(GetAssayData(object, assay = "RNA", layer = layer))
  }

  matching_layers <- grep(paste0("^", layer, "\\."), available_layers,
                          value = TRUE)
  if (length(matching_layers) > 1) {
    object <- JoinLayers(object, assay = "RNA")
    return(GetAssayData(object, assay = "RNA", layer = layer))
  }

  stop("RNA layer `", layer, "` was not found. Available layers: ",
       paste(available_layers, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 1. Load validated macrophage and refined global objects
# ------------------------------------------------------------------------------

macrophage_rds <- first_existing_file(
  MACROPHAGE_RDS_CANDIDATES,
  "A cleaned macrophage object"
)
global_rds <- first_existing_file(
  GLOBAL_RDS_CANDIDATES,
  "A global integrated Seurat object"
)

message("Loading macrophage object: ", macrophage_rds)
macrophage_seurat <- readRDS(macrophage_rds)
message("Loading global object: ", global_rds)
global_seurat <- readRDS(global_rds)

if (!"orig.ident" %in% colnames(macrophage_seurat[[]]) ||
    !"orig.ident" %in% colnames(global_seurat[[]])) {
  stop("Both objects must contain `orig.ident` metadata.")
}

mac_label_candidates <- c(
  "final_macrophage_label",
  "final_annotation",
  "macrophage_subcluster"
)
mac_label_col <- mac_label_candidates[
  mac_label_candidates %in% colnames(macrophage_seurat[[]])
][1]
if (is.na(mac_label_col)) {
  stop("No validated macrophage-state column was found. Checked: ",
       paste(mac_label_candidates, collapse = ", "))
}

global_label_candidates <- c("cell_type_refined", "cell_type", "cell_type_Xu")
global_label_col <- global_label_candidates[
  global_label_candidates %in% colnames(global_seurat[[]])
][1]
if (is.na(global_label_col)) {
  stop("No global cell-type annotation column was found.")
}

message("Macrophage state column: ", mac_label_col)
message("Global cell-type column: ", global_label_col)

# ------------------------------------------------------------------------------
# 2. Harmonise clinical-stage metadata from Xu supplementary sample mapping
# ------------------------------------------------------------------------------

stage_map <- c(
  "norm1" = "Non-malignant", "norm2" = "Non-malignant",
  "norm3" = "Non-malignant", "norm4" = "Non-malignant",
  "norm5" = "Non-malignant",
  "cancer1" = "Advanced", "cancer2" = "Advanced",
  "cancer3" = "Early", "cancer4" = "Early",
  "cancer5" = "Advanced", "cancer6" = "Advanced",
  "cancer7" = "Early"
)

figo_map <- c(
  "norm1" = "Non-malignant", "norm2" = "Non-malignant",
  "norm3" = "Non-malignant", "norm4" = "Non-malignant",
  "norm5" = "Non-malignant",
  "cancer1" = "IIIB", "cancer2" = "IIB", "cancer3" = "IC2",
  "cancer4" = "IC2", "cancer5" = "IIB", "cancer6" = "IIIC",
  "cancer7" = "IC2"
)

extract_paper_sample_id <- function(x) {
  cleaned <- tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
  matched <- grepl("norm[1-5]|cancer[1-7]", cleaned)
  output <- rep(NA_character_, length(cleaned))
  output[matched] <- sub(".*(norm[1-5]|cancer[1-7]).*", "\\1",
                         cleaned[matched])
  output
}

add_stage_metadata <- function(object) {
  object$paper_sample_id <- extract_paper_sample_id(object$orig.ident)
  object$clinical_stage_group <- factor(
    unname(stage_map[object$paper_sample_id]),
    levels = c("Non-malignant", "Early", "Advanced")
  )
  object$figo_stage <- unname(figo_map[object$paper_sample_id])
  if (anyNA(object$clinical_stage_group)) {
    unresolved <- unique(as.character(object$orig.ident)[
      is.na(object$clinical_stage_group)
    ])
    stop("Stage mapping failed for orig.ident values: ",
         paste(unresolved, collapse = ", "))
  }
  object
}

macrophage_seurat <- add_stage_metadata(macrophage_seurat)
global_seurat <- add_stage_metadata(global_seurat)

stage_audit <- unique(data.frame(
  orig.ident = as.character(global_seurat$orig.ident),
  paper_sample_id = global_seurat$paper_sample_id,
  clinical_stage_group = as.character(global_seurat$clinical_stage_group),
  figo_stage = global_seurat$figo_stage
))
write.csv(stage_audit, file.path(RESULT_DIR, "04_sample_stage_audit.csv"),
          row.names = FALSE)

# Retain all Non-malignant, Early, and Advanced samples for comparisons.
macrophage_stage <- subset(
  macrophage_seurat,
  subset = clinical_stage_group %in% c("Non-malignant", "Early", "Advanced")
)
macrophage_stage$clinical_stage_group <- factor(
  macrophage_stage$clinical_stage_group,
  levels = c("Non-malignant", "Early", "Advanced")
)

# ------------------------------------------------------------------------------
# 3. Sample/state audit before inference
# ------------------------------------------------------------------------------

macrophage_stage$analysis_state <- as.character(
  macrophage_stage[[mac_label_col, drop = TRUE]]
)
if (anyNA(macrophage_stage$analysis_state) ||
    any(macrophage_stage$analysis_state == "")) {
  stop("Validated macrophage labels contain missing or empty values.")
}

state_sample_counts <- macrophage_stage[[]] %>%
  transmute(
    sample = as.character(orig.ident),
    stage = as.character(clinical_stage_group),
    macrophage_state = analysis_state
  ) %>%
  count(sample, stage, macrophage_state, name = "n_cells")

write.csv(state_sample_counts,
          file.path(RESULT_DIR, "04_macrophage_state_cells_by_sample.csv"),
          row.names = FALSE)

p_state_counts <- ggplot(
  state_sample_counts,
  aes(x = sample, y = macrophage_state, fill = n_cells)
) +
  geom_tile(color = "white") +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "C") +
  scale_y_discrete(labels = function(x) compact_wrap_labels(x, 22L)) +
  labs(
    title = "Macrophage-state representation by biological sample",
    x = NULL, y = NULL, fill = "Cells"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 25, 35, 10)
  ) +
  coord_cartesian(clip = "off")
save_png(p_state_counts,
         file.path(FIGURE_DIR, "04_macrophage_state_sample_audit.png"),
         14, 8)

# ------------------------------------------------------------------------------
# 4. Sample-aware pseudobulk expression for decoupleR
# ------------------------------------------------------------------------------

DefaultAssay(macrophage_stage) <- "RNA"
rna_counts <- get_rna_layer(macrophage_stage, "counts")

pb_metadata <- macrophage_stage[[]] %>%
  transmute(
    cell = rownames(macrophage_stage[[]]),
    sample = as.character(orig.ident),
    stage = as.character(clinical_stage_group),
    macrophage_state = analysis_state,
    pseudobulk_id = paste(sample, macrophage_state, sep = "___")
  )
pb_metadata <- pb_metadata[match(colnames(rna_counts), pb_metadata$cell), ]
if (anyNA(pb_metadata$cell)) stop("RNA count cells and metadata do not match.")

pb_factor <- factor(pb_metadata$pseudobulk_id)
design_cells <- Matrix::sparse.model.matrix(~ 0 + pb_factor)
colnames(design_cells) <- levels(pb_factor)
pseudobulk_counts <- rna_counts %*% design_cells

pb_info <- pb_metadata %>%
  count(pseudobulk_id, sample, stage, macrophage_state, name = "n_cells") %>%
  arrange(match(pseudobulk_id, colnames(pseudobulk_counts)))
if (!identical(pb_info$pseudobulk_id, colnames(pseudobulk_counts))) {
  stop("Pseudobulk column metadata order does not match the count matrix.")
}

eligible_pb <- pb_info$n_cells >= MIN_CELLS_PER_SAMPLE_STATE
pseudobulk_counts <- pseudobulk_counts[, eligible_pb, drop = FALSE]
pb_info <- pb_info[eligible_pb, , drop = FALSE]

write.csv(pb_info,
          file.path(RESULT_DIR, "04_eligible_pseudobulk_samples.csv"),
          row.names = FALSE)

if (ncol(pseudobulk_counts) < 4) {
  stop("Fewer than four eligible sample-state pseudobulks remain. Reduce ",
       "MIN_CELLS_PER_SAMPLE_STATE only after reviewing the sample audit.")
}

library_size <- Matrix::colSums(pseudobulk_counts)
if (any(library_size == 0)) stop("A pseudobulk library has zero counts.")
log_cpm <- log2(t(t(as.matrix(pseudobulk_counts)) / library_size) * 1e6 + 1)

# ------------------------------------------------------------------------------
# 5. DoRothEA/decoupleR ULM transcription-factor activity
# ------------------------------------------------------------------------------

# Load a versioned local copy of the human DoRothEA regulon. This avoids
# analysis failures caused by temporary OmniPath or Ensembl outages.
if (!requireNamespace("dorothea", quietly = TRUE)) {
  stop(
    "The `dorothea` package is required. Install it with: ",
    "BiocManager::install('dorothea')"
  )
}

data("dorothea_hs", package = "dorothea", envir = environment())

required_dorothea_columns <- c("tf", "target", "mor", "confidence")
missing_dorothea_columns <- setdiff(
  required_dorothea_columns,
  colnames(dorothea_hs)
)
if (length(missing_dorothea_columns) > 0) {
  stop(
    "Required DoRothEA columns are missing: ",
    paste(missing_dorothea_columns, collapse = ", ")
  )
}

dorothea_net <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  transmute(
    source = tf,
    target = target,
    mor = mor,
    confidence = confidence
  ) %>%
  distinct(source, target, .keep_all = TRUE)

if (nrow(dorothea_net) == 0) {
  stop("The offline DoRothEA A/B/C regulon is empty.")
}
message(
  "Loaded offline DoRothEA A/B/C regulon: ",
  nrow(dorothea_net), " TF-target interactions and ",
  n_distinct(dorothea_net$source), " transcription factors."
)

write.csv(dorothea_net,
          file.path(RESULT_DIR, "04_dorothea_ABC_network_used.csv"),
          row.names = FALSE)

set.seed(RANDOM_SEED)
ulm_long <- decoupleR::run_ulm(
  mat = log_cpm,
  net = dorothea_net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
) %>%
  rename(pseudobulk_id = condition, TF = source, activity = score) %>%
  left_join(pb_info, by = "pseudobulk_id")

write.csv(ulm_long,
          file.path(RESULT_DIR, "04_ULM_TF_activity_by_sample_state.csv"),
          row.names = FALSE)

# Sample-level pairwise TF comparisons are performed separately within each
# macrophage state. Cells are never treated as independent replicates.
tf_contrasts <- data.frame(
  numerator = c("Early", "Advanced", "Advanced"),
  denominator = c("Non-malignant", "Non-malignant", "Early"),
  contrast = c(
    "Early_vs_Non-malignant",
    "Advanced_vs_Non-malignant",
    "Advanced_vs_Early"
  ),
  stringsAsFactors = FALSE
)

test_tf_contrast <- function(state_name, numerator, denominator,
                             contrast_name, activity_table) {
  state_df <- activity_table %>%
    filter(
      macrophage_state == state_name,
      stage %in% c(numerator, denominator),
      is.finite(activity)
    )

  # Count independent sample-state pseudobulks, not repeated TF rows.
  sample_information <- state_df %>%
    distinct(pseudobulk_id, sample, stage)
  stage_n <- table(factor(
    sample_information$stage,
    levels = c(denominator, numerator)
  ))

  if (length(stage_n) != 2L || any(stage_n < MIN_SAMPLES_PER_STAGE)) {
    message(
      "Skipping ", state_name, " | ", contrast_name,
      ": insufficient independent pseudobulks (",
      denominator, "=", unname(stage_n[denominator]),
      ", ", numerator, "=", unname(stage_n[numerator]), ")."
    )
    return(NULL)
  }

  activity_wide <- state_df %>%
    select(TF, pseudobulk_id, activity) %>%
    distinct() %>%
    pivot_wider(names_from = pseudobulk_id, values_from = activity)
  activity_matrix <- as.matrix(activity_wide[, -1, drop = FALSE])
  rownames(activity_matrix) <- activity_wide$TF
  storage.mode(activity_matrix) <- "double"

  finite_rows <- apply(
    activity_matrix,
    1,
    function(x) all(is.finite(x))
  )
  activity_matrix <- activity_matrix[finite_rows, , drop = FALSE]
  if (nrow(activity_matrix) == 0) {
    message(
      "Skipping ", state_name, " | ", contrast_name,
      ": no TFs with complete finite activity values."
    )
    return(NULL)
  }

  column_information <- sample_information %>%
    slice(match(colnames(activity_matrix), pseudobulk_id))
  if (anyNA(column_information$pseudobulk_id)) {
    stop("TF activity columns could not be matched to pseudobulk metadata.")
  }

  stage_factor <- factor(
    column_information$stage,
    levels = c(denominator, numerator)
  )
  design <- model.matrix(~ stage_factor)
  residual_df <- nrow(design) - qr(design)$rank
  if (residual_df < 1L) {
    message(
      "Skipping ", state_name, " | ", contrast_name,
      ": no residual degrees of freedom."
    )
    return(NULL)
  }

  fit <- limma::lmFit(activity_matrix, design)
  finite_sigma <- is.finite(fit$sigma)
  if (!any(finite_sigma)) {
    message(
      "Skipping ", state_name, " | ", contrast_name,
      ": no finite residual standard deviations."
    )
    return(NULL)
  }
  fit <- fit[finite_sigma, ]

  fit <- tryCatch(
    limma::eBayes(fit),
    error = function(e) {
      warning(
        "limma eBayes skipped for ", state_name, " | ", contrast_name,
        ": ", conditionMessage(e)
      )
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  limma::topTable(
    fit, coef = 2, number = Inf, sort.by = "P"
  ) %>%
    rownames_to_column("TF") %>%
    transmute(
      macrophage_state = state_name,
      contrast = contrast_name,
      numerator = numerator,
      denominator = denominator,
      TF,
      activity_difference = logFC,
      t,
      p_value = P.Value,
      FDR = adj.P.Val,
      n_numerator = unname(stage_n[numerator]),
      n_denominator = unname(stage_n[denominator])
    )
}

tf_differential <- bind_rows(lapply(
  sort(unique(ulm_long$macrophage_state)),
  function(state_name) {
    bind_rows(lapply(seq_len(nrow(tf_contrasts)), function(i) {
      test_tf_contrast(
        state_name = state_name,
        numerator = tf_contrasts$numerator[i],
        denominator = tf_contrasts$denominator[i],
        contrast_name = tf_contrasts$contrast[i],
        activity_table = ulm_long
      )
    }))
  }
))

write.csv(tf_differential,
          file.path(RESULT_DIR,
                    "04_TF_activity_all_pairwise_comparisons_by_state.csv"),
          row.names = FALSE)

tf_testability <- pb_info %>%
  distinct(sample, stage, macrophage_state) %>%
  count(macrophage_state, stage, name = "n_samples") %>%
  complete(
    macrophage_state,
    stage = c("Non-malignant", "Early", "Advanced"),
    fill = list(n_samples = 0)
  ) %>%
  pivot_wider(names_from = stage, values_from = n_samples,
              names_prefix = "n_") %>%
  mutate(
    testable_Early_vs_Non_malignant =
      n_Early >= MIN_SAMPLES_PER_STAGE &
      `n_Non-malignant` >= MIN_SAMPLES_PER_STAGE,
    testable_Advanced_vs_Non_malignant =
      n_Advanced >= MIN_SAMPLES_PER_STAGE &
      `n_Non-malignant` >= MIN_SAMPLES_PER_STAGE,
    testable_Advanced_vs_Early =
      n_Advanced >= MIN_SAMPLES_PER_STAGE &
      n_Early >= MIN_SAMPLES_PER_STAGE
  )
write.csv(tf_testability,
          file.path(RESULT_DIR, "04_TF_state_testability_audit.csv"),
          row.names = FALSE)

if (nrow(tf_differential) > 0) {
  selected_tfs <- tf_differential %>%
    group_by(macrophage_state) %>%
    arrange(FDR, desc(abs(activity_difference)), .by_group = TRUE) %>%
    slice_head(n = TOP_TFS_HEATMAP_PER_STATE) %>%
    ungroup() %>%
    distinct(TF) %>%
    pull(TF)

  tf_heatmap_df <- ulm_long %>%
    filter(TF %in% selected_tfs) %>%
    group_by(TF, macrophage_state, stage) %>%
    summarise(mean_activity = mean(activity), .groups = "drop") %>%
    mutate(
      state_short = compact_cell_labels(macrophage_state),
      stage = factor(
        stage,
        levels = c("Non-malignant", "Early", "Advanced")
      )
    )

  state_order <- unique(c(
    unname(DISPLAY_STATE_LABELS),
    sort(setdiff(
      unique(tf_heatmap_df$state_short),
      unname(DISPLAY_STATE_LABELS)
    ))
  ))
  tf_heatmap_df$state_short <- factor(
    tf_heatmap_df$state_short,
    levels = state_order
  )

  p_tf_heatmap <- ggplot(
    tf_heatmap_df,
    aes(x = state_short, y = TF, fill = mean_activity)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(
      . ~ stage,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_x_discrete(labels = function(x) wrap_plot_labels(x, 18L)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0) +
    labs(
      title = "Sample-aware DoRothEA/ULM transcription-factor activity",
      subtitle = paste0(
        "Mean activity of the top ", TOP_TFS_HEATMAP_PER_STATE,
        " stage-associated TFs per macrophage state"
      ),
      x = NULL, y = NULL, fill = "ULM activity"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", size = 10),
      panel.spacing.x = grid::unit(0.8, "lines"),
      plot.margin = margin(12, 25, 55, 12),
      legend.position = "right"
    ) +
    coord_cartesian(clip = "off")
  save_png(p_tf_heatmap,
           file.path(FIGURE_DIR, "04_TF_activity_state_stage_heatmap.png"),
           20, 12)
}

# ------------------------------------------------------------------------------
# 6. Build CellChat identities from refined global types + macrophage states
# ------------------------------------------------------------------------------

DefaultAssay(global_seurat) <- "RNA"
global_labels <- as.character(global_seurat[[global_label_col, drop = TRUE]])
global_seurat$cellchat_identity <- global_labels

# Replace only validated macrophage/myeloid cells with their high-resolution
# biological state. Cell barcodes must overlap the global object.
shared_macrophage_cells <- intersect(
  colnames(macrophage_seurat),
  colnames(global_seurat)
)
if (length(shared_macrophage_cells) == 0) {
  stop("No cell barcodes overlap between macrophage and global objects.")
}

mac_state_by_cell <- setNames(
  as.character(macrophage_seurat[[mac_label_col, drop = TRUE]]),
  colnames(macrophage_seurat)
)
global_seurat$cellchat_identity[
  match(shared_macrophage_cells, colnames(global_seurat))
] <- mac_state_by_cell[shared_macrophage_cells]

# Refine broad T cells into CD8 T and Other T using positive multi-gene evidence.
# This is used only for CellChat targeting and does not overwrite cell_type.
t_cell_candidates <- grepl("(^|[ _-])T([ _-]|$)|T cell",
                           global_labels, ignore.case = TRUE)
rna_data_global <- get_rna_layer(global_seurat, "data")
cd8_genes <- intersect(c("CD8A", "CD8B", "CCL5", "NKG7"),
                       rownames(rna_data_global))
cd4_genes <- intersect(c("CD4", "IL7R", "LTB", "MAL"),
                       rownames(rna_data_global))

if (any(t_cell_candidates) && length(cd8_genes) >= 2 && length(cd4_genes) >= 2) {
  candidate_cells <- colnames(global_seurat)[t_cell_candidates]
  cd8_score <- Matrix::colMeans(
    rna_data_global[cd8_genes, candidate_cells, drop = FALSE] > 0
  )
  cd4_score <- Matrix::colMeans(
    rna_data_global[cd4_genes, candidate_cells, drop = FALSE] > 0
  )
  global_seurat$cellchat_identity[
    match(candidate_cells[cd8_score > cd4_score & cd8_score >= 0.25],
          colnames(global_seurat))
  ] <- "CD8 T"
  remaining_t <- setdiff(candidate_cells,
                         candidate_cells[cd8_score > cd4_score &
                                           cd8_score >= 0.25])
  global_seurat$cellchat_identity[
    match(remaining_t, colnames(global_seurat))
  ] <- "Other T"
}

cellchat_stage <- subset(
  global_seurat,
  subset = clinical_stage_group %in%
    c("Non-malignant", "Early", "Advanced")
)
cellchat_stage$clinical_stage_group <- factor(
  cellchat_stage$clinical_stage_group,
  levels = c("Non-malignant", "Early", "Advanced")
)

cellchat_count_audit <- cellchat_stage[[]] %>%
  transmute(
    sample = as.character(orig.ident),
    stage = as.character(clinical_stage_group),
    identity = as.character(cellchat_identity)
  ) %>%
  count(sample, stage, identity, name = "n_cells")
write.csv(cellchat_count_audit,
          file.path(RESULT_DIR, "04_CellChat_cells_by_sample_identity.csv"),
          row.names = FALSE)

# Retain identities present with at least CELLCHAT_MIN_CELLS in all groups.
identity_stage_totals <- cellchat_count_audit %>%
  group_by(stage, identity) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  pivot_wider(names_from = stage, values_from = n_cells, values_fill = 0) %>%
  mutate(
    retained_for_comparison = `Non-malignant` >= CELLCHAT_MIN_CELLS &
      Early >= CELLCHAT_MIN_CELLS & Advanced >= CELLCHAT_MIN_CELLS,
    balanced_n = pmin(`Non-malignant`, Early, Advanced,
                      CELLCHAT_MAX_CELLS_PER_IDENTITY_STAGE)
  )
write.csv(identity_stage_totals,
          file.path(RESULT_DIR, "04_CellChat_identity_stage_audit.csv"),
          row.names = FALSE)

shared_identities <- identity_stage_totals$identity[
  identity_stage_totals$retained_for_comparison
]
if (length(shared_identities) < 2) {
  stop("Fewer than two CellChat identities are represented in all groups.")
}

# Balance each retained identity across groups to reduce abundance-driven
# differences and make CellChat processing directly comparable.
set.seed(RANDOM_SEED)
balanced_cells <- unlist(lapply(shared_identities, function(identity_name) {
  target_n <- identity_stage_totals$balanced_n[
    identity_stage_totals$identity == identity_name
  ]
  unlist(lapply(c("Non-malignant", "Early", "Advanced"), function(stage_name) {
    eligible <- colnames(cellchat_stage)[
      cellchat_stage$cellchat_identity == identity_name &
        cellchat_stage$clinical_stage_group == stage_name
    ]
    sample(eligible, size = min(target_n, length(eligible)), replace = FALSE)
  }), use.names = FALSE)
}), use.names = FALSE)

cellchat_balanced <- subset(cellchat_stage, cells = balanced_cells)
DefaultAssay(cellchat_balanced) <- "RNA"

# ------------------------------------------------------------------------------
# 7. Consistently processed Non-malignant, Early, and Advanced CellChat models
# ------------------------------------------------------------------------------

create_cellchat_object <- function(seurat_object, stage_name) {
  stage_object <- subset(
    seurat_object,
    subset = clinical_stage_group == stage_name
  )
  stage_object$cellchat_identity <- droplevels(
    factor(stage_object$cellchat_identity, levels = shared_identities)
  )

  expression_data <- get_rna_layer(stage_object, "data")
  cell_metadata <- stage_object[[]] %>%
    transmute(
      labels = as.character(cellchat_identity),
      sample = as.character(orig.ident),
      stage = as.character(clinical_stage_group)
    )
  rownames(cell_metadata) <- rownames(stage_object[[]])
  cell_metadata <- cell_metadata[colnames(expression_data), , drop = FALSE]

  cellchat <- createCellChat(
    object = expression_data,
    meta = cell_metadata,
    group.by = "labels"
  )
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  set.seed(RANDOM_SEED)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = CELLCHAT_MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat
}

message("Building Non-malignant CellChat model...")
cellchat_non_malignant <- create_cellchat_object(
  cellchat_balanced,
  "Non-malignant"
)
message("Building Early CellChat model...")
cellchat_early <- create_cellchat_object(cellchat_balanced, "Early")
message("Building Advanced CellChat model...")
cellchat_advanced <- create_cellchat_object(cellchat_balanced, "Advanced")

saveRDS(cellchat_non_malignant,
        file.path(RESULT_DIR, "04_CellChat_Non_malignant_balanced.rds"))
saveRDS(cellchat_early,
        file.path(RESULT_DIR, "04_CellChat_Early_balanced.rds"))
saveRDS(cellchat_advanced,
        file.path(RESULT_DIR, "04_CellChat_Advanced_balanced.rds"))

non_malignant_lr <- subsetCommunication(cellchat_non_malignant)
early_lr <- subsetCommunication(cellchat_early)
advanced_lr <- subsetCommunication(cellchat_advanced)
write.csv(non_malignant_lr,
          file.path(RESULT_DIR,
                    "04_CellChat_Non_malignant_interactions.csv"),
          row.names = FALSE)
write.csv(early_lr,
          file.path(RESULT_DIR, "04_CellChat_Early_interactions.csv"),
          row.names = FALSE)
write.csv(advanced_lr,
          file.path(RESULT_DIR, "04_CellChat_Advanced_interactions.csv"),
          row.names = FALSE)

cellchat_list <- list(
  `Non-malignant` = cellchat_non_malignant,
  Early = cellchat_early,
  Advanced = cellchat_advanced
)
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
saveRDS(cellchat_merged,
        file.path(RESULT_DIR,
                  "04_CellChat_Non_malignant_Early_Advanced_merged.rds"))

# ------------------------------------------------------------------------------
# 8. Comparative and macrophage-focused CellChat outputs
# ------------------------------------------------------------------------------

p_number <- compareInteractions(
  cellchat_merged, show.legend = FALSE, group = c(1, 2, 3)
) + ggtitle("Number of inferred interactions")
p_weight <- compareInteractions(
  cellchat_merged, show.legend = FALSE, group = c(1, 2, 3),
  measure = "weight"
) + ggtitle("Aggregate interaction strength")
save_png(p_number + p_weight,
         file.path(FIGURE_DIR, "04_CellChat_global_interaction_comparison.png"),
         10, 5)

cellchat_pair_indices <- list(
  Early_vs_Non_malignant = c(1, 2),
  Advanced_vs_Non_malignant = c(1, 3),
  Advanced_vs_Early = c(2, 3)
)
for (contrast_name in names(cellchat_pair_indices)) {
  comparison_indices <- cellchat_pair_indices[[contrast_name]]
  save_base_png(
    file.path(
      FIGURE_DIR,
      paste0("04_CellChat_diff_network_", contrast_name, ".png")
    ),
    16, 12,
    netVisual_diffInteraction(
      cellchat_merged,
      comparison = comparison_indices,
      weight.scale = TRUE,
      measure = "weight"
    )
  )
}

p_role_non_malignant <- netAnalysis_signalingRole_scatter(
  cellchat_non_malignant, title = "Non-malignant: signalling roles"
)
p_role_early <- netAnalysis_signalingRole_scatter(
  cellchat_early, title = "Early: signalling roles"
)
p_role_advanced <- netAnalysis_signalingRole_scatter(
  cellchat_advanced, title = "Advanced: signalling roles"
)
save_png(p_role_non_malignant + p_role_early + p_role_advanced,
         file.path(FIGURE_DIR, "04_CellChat_signaling_roles.png"),
         21, 8)

macrophage_identities <- intersect(
  unique(as.character(macrophage_seurat[[mac_label_col, drop = TRUE]])),
  shared_identities
)
epithelial_targets <- shared_identities[
  grepl("epithel|tumou?r", shared_identities, ignore.case = TRUE)
]
cd8_targets <- shared_identities[
  grepl("CD8", shared_identities, ignore.case = TRUE)
]
apc_targets <- shared_identities[
  grepl("CLEC10A|APC|dendritic|(^|[ _-])DC([ _-]|$)",
        shared_identities, ignore.case = TRUE)
]

# Generate compartment-specific bubbles rather than one overcrowded plot.
# Epithelial and CD8 figures are primary presentation outputs. APC-like figures
# are supplementary because Xu et al. did not resolve a separate DC/APC group.

priority_lr_labels <- function(x) {
  labels <- wrap_plot_labels(x, 30L)
  is_priority <- grepl(
    paste(PRIORITY_SIGNAL_PATTERNS, collapse = "|"),
    x,
    ignore.case = TRUE
  )

  # ggtext allows selected pathway labels to be bold and dark red. If it is not
  # installed, all displayed significant LR labels remain bold in black.
  if (requireNamespace("ggtext", quietly = TRUE)) {
    labels[is_priority] <- paste0(
      "<b><span style='color:#9E2A2B'>", labels[is_priority], "</span></b>"
    )
  }
  labels
}

bubble_y_text <- function() {
  if (requireNamespace("ggtext", quietly = TRUE)) {
    ggtext::element_markdown(size = 9, lineheight = 1.05)
  } else {
    element_text(size = 9, face = "bold", lineheight = 1.05)
  }
}

estimate_bubble_height <- function(plot) {
  n_labels <- tryCatch({
    built <- ggplot_build(plot)
    panel <- built$layout$panel_params[[1]]
    labels <- panel$y$get_labels()
    length(unique(labels[!is.na(labels) & nzchar(labels)]))
  }, error = function(e) 20L)

  max(BUBBLE_MIN_HEIGHT, min(BUBBLE_MAX_HEIGHT, 5 + 0.34 * n_labels))
}

format_cellchat_bubble <- function(plot, title_text) {
  emphasis_note <- if (requireNamespace("ggtext", quietly = TRUE)) {
    "; bold dark-red labels identify priority signalling families"
  } else {
    "; ligand–receptor labels are enlarged and bold for readability"
  }
  plot +
    labs(
      title = title_text,
      subtitle = paste0(
        "Displayed interactions: CellChat p < ", BUBBLE_PVALUE_THRESHOLD,
        emphasis_note
      ),
      x = "Sender → receiver and disease-stage comparison",
      y = "Ligand–receptor interaction"
    ) +
    scale_x_discrete(
      labels = function(x) compact_wrap_labels(x, 22L),
      guide = guide_axis(n.dodge = 2)
    ) +
    scale_y_discrete(labels = priority_lr_labels) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"
      ),
      axis.text.y = bubble_y_text(),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 12)),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 12)),
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10, margin = margin(b = 10)),
      plot.margin = margin(15, 45, 95, 20),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25)
    ) +
    coord_cartesian(clip = "off")
}

save_bidirectional_bubbles <- function(targets, compartment_key,
                                       compartment_title,
                                       output_prefix = "04_CellChat_primary") {
  if (length(macrophage_identities) == 0 || length(targets) == 0) {
    warning("Skipping ", compartment_title,
            " bubbles: the required identities were not retained in all stages.")
    return(invisible(FALSE))
  }

  outgoing_saved <- tryCatch({
    p_outgoing <- netVisual_bubble(
      cellchat_merged,
      sources.use = macrophage_identities,
      targets.use = targets,
      comparison = c(1, 2, 3),
      thresh = BUBBLE_PVALUE_THRESHOLD,
      angle.x = 45,
      remove.isolate = TRUE
    )
    p_outgoing <- format_cellchat_bubble(
      p_outgoing,
      paste0(
        "Macrophage to ", compartment_title,
        ": outgoing ligand-receptor signalling"
      )
    )

    save_png(
      p_outgoing,
      file.path(
        FIGURE_DIR,
        paste0(output_prefix, "_macrophage_to_", compartment_key,
               "_outgoing_bubble.png")
      ),
      BUBBLE_WIDTH, estimate_bubble_height(p_outgoing)
    )
    TRUE
  }, error = function(e) {
    warning("Outgoing bubble skipped for ", compartment_title,
            ": ", conditionMessage(e))
    FALSE
  })

  incoming_saved <- tryCatch({
    p_incoming <- netVisual_bubble(
      cellchat_merged,
      sources.use = targets,
      targets.use = macrophage_identities,
      comparison = c(1, 2, 3),
      thresh = BUBBLE_PVALUE_THRESHOLD,
      angle.x = 45,
      remove.isolate = TRUE
    )
    p_incoming <- format_cellchat_bubble(
      p_incoming,
      paste0(
        compartment_title,
        " to macrophage: incoming ligand-receptor signalling"
      )
    )

    save_png(
      p_incoming,
      file.path(
        FIGURE_DIR,
        paste0(output_prefix, "_", compartment_key,
               "_to_macrophage_incoming_bubble.png")
      ),
      BUBBLE_WIDTH, estimate_bubble_height(p_incoming)
    )
    TRUE
  }, error = function(e) {
    warning("Incoming bubble skipped for ", compartment_title,
            ": ", conditionMessage(e))
    FALSE
  })

  invisible(outgoing_saved || incoming_saved)
}

# PRIMARY AXIS 1: tumour-directed macrophage communication.
save_bidirectional_bubbles(
  epithelial_targets,
  compartment_key = "epithelial",
  compartment_title = "epithelial cells",
  output_prefix = "04_CellChat_primary"
)

# PRIMARY AXIS 2: immune-directed macrophage communication.
save_bidirectional_bubbles(
  cd8_targets,
  compartment_key = "CD8_T",
  compartment_title = "CD8 T cells",
  output_prefix = "04_CellChat_primary"
)

# SUPPLEMENTARY AXIS: include only if the CLEC10A+/DC/APC-like identity passed
# the shared-stage and minimum-cell filters. Interpret as myeloid APC-like
# communication unless independent marker validation supports a definitive DC
# identity.
save_bidirectional_bubbles(
  apc_targets,
  compartment_key = "myeloid_APC_like",
  compartment_title = "CLEC10A+ myeloid APC-like cells",
  output_prefix = "04_CellChat_supplementary"
)

# EXPLORATORY AXES: create a separate pair of bubble plots for each remaining
# eligible compartment. Separate files preserve individual ligand-receptor
# detail without allowing these results to visually compete with the two
# explicitly labelled primary axes.
other_compartment_identities <- setdiff(
  shared_identities,
  unique(c(
    macrophage_identities,
    epithelial_targets,
    cd8_targets,
    apc_targets
  ))
)

safe_filename_key <- function(x) {
  key <- gsub("[^A-Za-z0-9]+", "_", x)
  key <- gsub("^_+|_+$", "", key)
  ifelse(nchar(key) == 0, "unnamed_compartment", key)
}

if (length(other_compartment_identities) > 0) {
  for (other_identity in other_compartment_identities) {
    save_bidirectional_bubbles(
      targets = other_identity,
      compartment_key = safe_filename_key(other_identity),
      compartment_title = other_identity,
      output_prefix = "04_CellChat_exploratory"
    )
  }
}

# Summarise communication involving macrophage states in machine-readable form.
tag_stage <- function(table, stage_name) {
  table$stage <- stage_name
  table
}
all_lr <- bind_rows(
  tag_stage(non_malignant_lr, "Non-malignant"),
  tag_stage(early_lr, "Early"),
  tag_stage(advanced_lr, "Advanced")
)
macrophage_lr <- all_lr %>%
  filter(source %in% macrophage_identities |
           target %in% macrophage_identities)
write.csv(macrophage_lr,
          file.path(RESULT_DIR,
                    "04_CellChat_all_macrophage_incoming_outgoing_LR.csv"),
          row.names = FALSE)

# Export the two primary axes, the supplementary APC-like axis and the broad
# all-other-compartment screen separately. Other compartments also receive
# separate exploratory bubble plots; no aggregated heatmap is produced.
extract_bidirectional_axis <- function(targets) {
  if (length(targets) == 0) return(macrophage_lr[0, , drop = FALSE])
  macrophage_lr %>%
    filter(
      (source %in% macrophage_identities & target %in% targets) |
        (source %in% targets & target %in% macrophage_identities)
    )
}

epithelial_macrophage_lr <- extract_bidirectional_axis(epithelial_targets)
cd8_macrophage_lr <- extract_bidirectional_axis(cd8_targets)
apc_macrophage_lr <- extract_bidirectional_axis(apc_targets)
main_and_supplementary_targets <- unique(c(
  macrophage_identities, epithelial_targets, cd8_targets, apc_targets
))
other_compartment_lr <- macrophage_lr %>%
  filter(
    (source %in% macrophage_identities &
       !(target %in% main_and_supplementary_targets)) |
      (target %in% macrophage_identities &
         !(source %in% main_and_supplementary_targets))
  )

write.csv(
  epithelial_macrophage_lr,
  file.path(RESULT_DIR,
            "04_CellChat_PRIMARY_macrophage_epithelial_LR.csv"),
  row.names = FALSE
)
write.csv(
  cd8_macrophage_lr,
  file.path(RESULT_DIR,
            "04_CellChat_PRIMARY_macrophage_CD8_T_LR.csv"),
  row.names = FALSE
)
write.csv(
  apc_macrophage_lr,
  file.path(RESULT_DIR,
            "04_CellChat_SUPPLEMENTARY_macrophage_myeloid_APC_like_LR.csv"),
  row.names = FALSE
)
write.csv(
  other_compartment_lr,
  file.path(RESULT_DIR,
            "04_CellChat_EXPLORATORY_macrophage_other_compartments_LR.csv"),
  row.names = FALSE
)

# Pathway-level summary for prioritising CXCL, CCL, TGF-beta and VEGF outputs.
pathway_summary <- macrophage_lr %>%
  group_by(stage, pathway_name, source, target) %>%
  summarise(
    n_LR_pairs = n(),
    summed_probability = sum(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(pathway_name, source, target, stage)
write.csv(pathway_summary,
          file.path(RESULT_DIR,
                    "04_CellChat_macrophage_pathway_summary.csv"),
          row.names = FALSE)

priority_pathways <- pathway_summary %>%
  filter(grepl("CXCL|CCL|TGF|VEGF", pathway_name, ignore.case = TRUE))
write.csv(priority_pathways,
          file.path(RESULT_DIR,
                    "04_CellChat_priority_pathways_CXCL_CCL_TGF_VEGF.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. Analysis manifest and reproducibility record
# ------------------------------------------------------------------------------

analysis_manifest <- c(
  paste0("macrophage_input=", macrophage_rds),
  paste0("global_input=", global_rds),
  paste0("macrophage_label_column=", mac_label_col),
  paste0("global_label_column=", global_label_col),
  "comparison=Non-malignant_vs_Early_vs_Advanced",
  "normal_tissue_label=Non-malignant",
  "non_malignant_samples_used_in_comparison=TRUE",
  paste0("min_cells_per_sample_state=", MIN_CELLS_PER_SAMPLE_STATE),
  paste0("min_samples_per_stage=", MIN_SAMPLES_PER_STAGE),
  "TF_method=DoRothEA_ABC_ULM_sample_state_pseudobulk",
  paste0("CellChat_min_cells=", CELLCHAT_MIN_CELLS),
  paste0("CellChat_max_cells_per_identity_stage=",
         CELLCHAT_MAX_CELLS_PER_IDENTITY_STAGE),
  "CellChat_identity_counts_balanced_between_stages=TRUE",
  "primary_CellChat_axes=macrophage_epithelial;macrophage_CD8_T",
  "supplementary_CellChat_axis=macrophage_CLEC10A_myeloid_APC_like",
  "other_compartments=separate_exploratory_bubble_plots_and_screening_table",
  "all_compartment_heatmap=FALSE",
  paste0("random_seed=", RANDOM_SEED),
  "interpretation=computational_hypotheses_not_mechanistic_proof"
)
writeLines(analysis_manifest,
           file.path(RESULT_DIR, "04_analysis_manifest.txt"))
writeLines(capture.output(sessionInfo()),
           file.path(RESULT_DIR, "04_session_info.txt"))

message("Downstream analysis completed.")
message("TF results: ", RESULT_DIR)
message("CellChat results: ", RESULT_DIR)
message("Figures: ", FIGURE_DIR)
