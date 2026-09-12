# Dataset Information

This directory contains metadata and references for the single-cell RNA sequencing (scRNA-seq) datasets analyzed in this repository.

## Overview
- **Dataset Accession Number:** GSE184880
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Download Link:** [GSE184880 GEO Accession Page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880)

## Sample Breakdown
The analysis focuses on three specific 10x Genomics scRNA-seq samples from this dataset:
1. **GSM5599220_Norm1** (Normal Tissue Sample 1)
2. **GSM5599225_Cancer1** (Cancer Tissue Sample 1)
3. **GSM5599226_Cancer2** (Cancer Tissue Sample 2)

## Data Handling & Storage Notice
- Large raw expression matrices and computed Seurat objects (`.rds` files) are excluded from Git tracking via `.gitignore` to maintain efficient repository size.
- Processed summary statistics and quality control metrics are exported to the `results/` folder, while generated figures are organized in `figures/`.
