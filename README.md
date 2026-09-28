# HGSOC single-cell reanalysis (GSE184880): transcription-factor activity and cell communication in myeloid/macrophage states across stages

Reanalysis of the Xu et al. dataset GSE184880 to identify evidence-supported myeloid/macrophage states and determine how their transcription-factor activity and predicted communication programs differ across HGSOC stages.

## The question

Which transcription-factor activities and cell-communication programs distinguish evidence-supported myeloid/macrophage states across Non-malignant, Early, and Advanced HGSOC tissue groups?

## The data

- **Study:** Xu J, Fang Y, Chen K, et al. *Single-cell RNA sequencing reveals the tissue architecture in human high-grade serous ovarian cancer.* Clin Cancer Res. 2022;28(16):3590–3602. doi:[10.1158/1078-0432.CCR-22-0296](https://doi.org/10.1158/1078-0432.CCR-22-0296)
- **Public accession:** [GSE184880](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880) (NCBI GEO; 10x Genomics single-cell RNA-seq).
- **Input:** Processed 10x matrices (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) for twelve samples.
- **How to download:** Download the authors' processed matrices from the GSE184880 page, extract the 12 sample folders, and place them under `data/GSE184880/`. The full sample list and expected folder structure are in [`data/README.md`](data/README.md).
- **Data policy:** Expression matrices and `.rds` objects are excluded from GitHub; the repository contains code, small result tables, and final figures only.

### Sample design

| Group | Samples | n |
|---|---|---:|
| Non-malignant | Norm1–Norm5 | 5 |
| Early HGSOC | Cancer3, Cancer4, Cancer7 | 3 |
| Advanced HGSOC | Cancer1, Cancer2, Cancer5, Cancer6 | 4 |

## What we did

1. Applied sample-level quality control and `scDblFinder`, then merged the twelve samples.
2. Normalized the data, selected highly variable genes, scaled, ran PCA, assessed sample mixing, integrated by sample identity, and annotated the global atlas.
3. Reproduced the published global landscape as a pipeline-validation baseline and refined the broad myeloid annotation.
4. Reclustered the myeloid/macrophage compartment and assigned biological states using multiple markers and cluster differential expression.
5. Estimated sample-aware transcription-factor activity using DoRothEA/ULM.
6. Compared stage-stratified CellChat predictions, focusing on communication involving macrophages, epithelial cells, and CD8 T cells.

## Quality control summary

Doublets were identified with `scDblFinder` per sample. Of 59,661 cells before doublet removal, 4,558 (7.6%) were removed as doublets, leaving **55,103 singlets** (92.4% retained). No doublets remained after subsetting.

| Group | Samples | Cells before doublet removal | Doublets removed | Singlets retained |
|---|---:|---:|---:|---:|
| Non-malignant | 5 | 26,115 | 1,940 | 24,175 |
| Early | 3 | 11,975 | 797 | 11,178 |
| Advanced | 4 | 21,571 | 1,821 | 19,750 |
| **Total** | **12** | **59,661** | **4,558** | **55,103** |

Across samples, the median genes per cell in the retained singlets ranged from 1,313 to 3,313, the median UMIs per cell from 3,381 to 12,457, and the mean mitochondrial read percentage from 9.1% to 19.9% (highest in Cancer4). After removal of mitochondrial genes, between 17,549 and 20,041 genes remained per sample. Per-sample counts are in `results/qc_cell_counts_summary.csv` and `results/post_doublet_subset_summary.csv`.

## How to run it

Run R from the repository root so all relative paths resolve correctly.

| Order | Script | Main output |
|---:|---|---|
| 0 | `scripts/00_project_setup.R` | Folder/package checks and `results/00_session_info.txt` |
| 1 | `scripts/01_qc_doublet_removal_and_merge.R` | `results/postQC_merged_12_samples.rds` and QC audit tables |
| 2 | `scripts/02_preprocessing_integration_global_annotation.R` | `results/02_integrated_seurat.rds` |
| 2.1 | `scripts/02.1_refine_global_annotation.R` | `results/02.1_integrated_seurat_refined.rds` |
| 2.2 | `scripts/02.2_reproduce_fig1C_fig1E.R` | Reproduction figures and sample-stage mapping tables |
| 3 | `scripts/03_macrophage_subclustering_annotation.R` | `results/03_clean_macrophages_Xu_annotated.rds` |
| 4 | `scripts/04_downstream_TF_activity_CellChat.R` | `results/04_downstream/` and `figures/04_downstream/` |

```r
source("scripts/00_project_setup.R")
source("scripts/01_qc_doublet_removal_and_merge.R")
source("scripts/02_preprocessing_integration_global_annotation.R")
source("scripts/02.1_refine_global_annotation.R")
source("scripts/02.2_reproduce_fig1C_fig1E.R")
source("scripts/03_macrophage_subclustering_annotation.R")
source("scripts/04_downstream_TF_activity_CellChat.R")
```

Stage 04 is computationally expensive and can be run separately. Intermediate `.rds` files prevent unnecessary reruns.

### Repository layout

- `scripts/`: numbered analysis code, in the order it runs
- `results/`: output tables and session information
- `figures/`: high-resolution figures used in this README and the presentation
- `data/`: accession and download instructions only (no data files)

## Requirements

- **R version:** 4.6.1 (2026-06-24), run on Windows 11 (x86_64-w64-mingw32).
- **Core packages and versions:**

| Package | Version | Package | Version |
|---|---|---|---|
| Seurat | 5.5.1 | decoupleR | 2.17.0 |
| SeuratObject | 5.4.0 | dorothea | 1.23.0 |
| scDblFinder | 1.26.7 | CellChat | 2.2.0.9001 |
| SingleCellExperiment | 1.34.0 | limma | 3.68.5 |
| batchelor | 1.28.0 | BiocSingular | 1.28.0 |
| dplyr | 1.2.1 | tidyr | 1.3.2 |
| tibble | 3.3.1 | ggplot2 | 4.0.3 |
| patchwork | 1.3.2 | pheatmap | 1.0.13 |

- The setup script checks that the required R packages are installed. The full `sessionInfo()` output for each stage (01, 02, 02.1, 02.2, 03, and 04) is saved as a stage-specific session-information file in `results/`.

Key settings are more than 200 detected genes and less than 40% mitochondrial reads; `LogNormalize` with scale factor 10,000; 2,000 `vst` HVGs; 10 dimensions; CCA integration by `orig.ident`; UMAP with 30 neighbors, minimum distance 0.3, and cosine metric; and clustering resolution 0.08. Random seeds are fixed throughout (`set.seed(42)`).

## Results

### Reproduced global cell landscape

![Reproduced global UMAP of the major cell populations](figures/Figure1B_reproduced_umap.png)

The reanalysis recovered the eight principal populations reported in the baseline study: T cells, epithelial cells, fibroblasts, monocytic/myeloid cells, endothelial cells, cycling cells, B/plasma cells, and smooth-muscle/myofibroblast cells. This supports the validity of the preprocessing and annotation workflow, although exact UMAP geometry and cluster proportions were not expected to match because of software-version and parameter differences.

### Evidence-supported myeloid/macrophage states

![UMAP of evidence-supported myeloid and macrophage states](figures/03_macrophage_subcluster_umap_biological_labels.png)

Ten states were resolved within the cleaned myeloid/macrophage compartment:

- FOLR2/LYVE1 resident macrophages
- OLFML3/CX3CR1 resident/APC macrophages
- cycling macrophages
- TNF/CCL4 inflammatory macrophages
- APOBEC3A-like interferon-responsive macrophages
- FCN1/S100A8/S100A9 classical monocytes
- SLC2A1-like hypoxic/glycolytic macrophages
- NUPR1/IGF2 stress-adapted macrophages
- CLEC10A APC-like myeloid cells
- CDKN1C non-classical/intermediate monocytes

Labels were deliberately descriptive and evidence-based rather than forced into the paper's single-marker taxonomy.

### Stage-associated transcription-factor activity

![Heatmap of transcription-factor activity by macrophage state and clinical stage](figures/04_TF_activity_state_stage_heatmap.png)

The heatmap summarizes sample-aware DoRothEA/ULM activity patterns across macrophage states and Non-malignant, Early, and Advanced groups.

### Advanced versus Non-malignant communication

![CellChat differential network comparing Advanced and Non-malignant samples](figures/04_CellChat_diff_network_Advanced_vs_Non_malignant.png)

The differential CellChat network highlights predicted communication changes between Advanced HGSOC and Non-malignant tissue.

## Interpretation and limitations

The reanalysis reproduced the principal cellular architecture of GSE184880 and recovered recognizable macrophage biology, including resident, inflammatory/interferon-responsive, hypoxic/glycolytic, cycling, and monocyte-like programs. Rather than treating all globally monocytic cells as macrophages, an explicit lineage audit separated B-cell, cDC, and pDC signals before macrophage-specific inference. This improved interpretability but changed the analyzed cell set and therefore the UMAP geometry and cluster boundaries relative to the original study. Differences from Xu et al. are expected and scientifically explainable. Xu et al. used Seurat v3.1.4, MNN subclustering, and a 2021–2022 software environment, while this project used Seurat v5, current fastMNN, explicit random seeds, contemporary UMAP settings, and additional lineage-cleaning rules. Xu et al. labeled ten macrophage clusters with representative markers; this project used multi-gene programs and negative-lineage evidence, retaining Xu-like names only where support was specific. UMAP rotation, compactness, and intercluster distances are stochastic visual properties and were not interpreted as biological measurements.

TF-activity and cell-communication results are computational predictions, not evidence of direct TF binding, ligand-receptor engagement, or causality. The study does not establish a single linear M1-to-M2 or "malignant" macrophage transition; instead, it identifies coexisting resident, inflammatory, interferon-responsive, metabolic, stress, and recruited monocyte programs whose prevalence and predicted interactions vary across tissue groups.

**Key limitations:**
- **Small cohort** (5 non-malignant, 3 early, 4 advanced samples), limiting statistical power and increasing sensitivity to individual donors.
- **Stage and BRCA/HRR status are partly confounded**: all early tumors were BRCA/HRR wild-type, while 3 of 4 advanced tumors carried BRCA1/BRCA2/ATM/BRIP1 alterations, so observed stage effects may partly reflect genotype.
- **One advanced sample is a recurrence**, adding heterogeneity to that group.
- **Cross-sectional design**: cannot show individual macrophages transitioning from early to advanced states over time.
- **Analyst-defined lineage thresholds**; excluded cells remain available in a broader immune object for audit.
- **CellChat differences may be influenced by cell-type abundance**, not only by expression changes.

## Team

- **Hiba Hussein**: [GitHub profile](https://github.com/hibakhair24)
- **Leqaa Ibrahim**: [GitHub profile](https://github.com/leqaaibrahim)
- **Nour Ibrahim**: [GitHub profile](https://github.com/noour1762006)
