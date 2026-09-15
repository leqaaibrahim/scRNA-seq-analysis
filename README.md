# Data Directory & Access Instructions

This directory contains metadata, accession details, and download instructions for the single-cell RNA sequencing (scRNA-seq) datasets analyzed in this repository.

> **CRITICAL NOTICE:** In accordance with GitHub file size limits and project guidelines, raw count matrices, temporary cache files, and computed Seurat objects (`.rds` / `.RData`) are excluded from Git tracking via `.gitignore`[cite: 1]. **Do not commit raw or processed matrix files to this directory**[cite: 1].

---

## Dataset Overview

* **Dataset Accession:** [GSE184880](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880)[cite: 1]
* **Source:** NCBI Gene Expression Omnibus (GEO)[cite: 1]
* **Platform:** 10x Genomics Single-Cell RNA-Seq
* **Organism:** *Homo sapiens*

---

## Sample Breakdown (12 Samples)

The analysis encompasses 12 total samples across two biological conditions:

| Sample ID | Condition | Tissue / Group | GEO Accession |
| :--- | :--- | :--- | :--- |
| `GSM5599220_Norm1` | Normal | Healthy Control 1 | [GSM5599220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599220) |
| `GSM5599221_Norm2` | Normal | Healthy Control 2 | [GSM5599221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599221) |
| `GSM5599222_Norm3` | Normal | Healthy Control 3 | [GSM5599222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599222) |
| `GSM5599223_Norm4` | Normal | Healthy Control 4 | [GSM5599223](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599223) |
| `GSM5599224_Norm5` | Normal | Healthy Control 5 | [GSM5599224](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599224) |
| `GSM5599225_Cancer1` | Cancer | Tumor Sample 1 | [GSM5599225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599225) |
| `GSM5599226_Cancer2` | Cancer | Tumor Sample 2 | [GSM5599226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599226) |
| `GSM5599227_Cancer3` | Cancer | Tumor Sample 3 | [GSM5599227](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599227) |
| `GSM5599228_Cancer4` | Cancer | Tumor Sample 4 | [GSM5599228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599228) |
| `GSM5599229_Cancer5` | Cancer | Tumor Sample 5 | [GSM5599229](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599229) |
| `GSM5599230_Cancer6` | Cancer | Tumor Sample 6 | [GSM5599230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599230) |
| `GSM5599231_Cancer7` | Cancer | Tumor Sample 7 | [GSM5599231](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5599231) |

---

## How to Obtain Raw Data

### Option A: Direct Download from GEO
1. Navigate to the [GSE184880 GEO Accession Page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880)[cite: 1].
2. Scroll to the **Supplementary file** section at the bottom.
3. Download the supplementary count matrices or filtered 10x feature-barcode matrices for the 12 samples listed above.
4. Extract the downloaded files directly into this local `data/` folder on your computer.

### Option B: Shared Team HPC / Cloud Storage
For team members with access to our shared drive, pre-downloaded raw count matrices and intermediate post-QC Seurat `.rds` objects can be fetched from:
* **Shared Storage Path:** `[INSERT YOUR HPC DIRECTORY / CLOUD STORAGE LINK HERE]`

---

## Local Directory Layout (Untracked)

Once downloaded, your local `data/` directory should look like this (all data subfolders are ignored by Git)[cite: 1]:

```text
data/
├── README.md                          # Tracked in Git
├── raw_matrices/                      # Untracked by Git
│   ├── GSM5599220_Norm1/
│   ├── GSM5599221_Norm2/
│   ├── ...
│   └── GSM5599231_Cancer7/
└── processed_rds/                     # Untracked by Git
    └── postQC_merged_12_samples.rds