# Dataset Information

This directory contains the data-access instructions, sample metadata, and
storage guidance for the single-cell RNA sequencing dataset analysed in this
repository.

## Overview

- **Dataset accession:** GSE184880
- **Repository:** NCBI Gene Expression Omnibus (GEO)
- **Study:** Xu et al., human high-grade serous ovarian cancer (HGSOC)
- **Technology:** 10x Genomics single-cell RNA sequencing
- **Download link:** [GSE184880 GEO accession page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880)
- **Samples included:** 12 samples—5 non-malignant tissue samples and 7 HGSOC
  tissue samples

## Complete Sample Breakdown

| GEO sample ID | Sample name | Biological group | Analysis group |
|---|---|---|---|
| GSM5599220 | Norm1 | Non-malignant tissue | Non-malignant |
| GSM5599221 | Norm2 | Non-malignant tissue | Non-malignant |
| GSM5599222 | Norm3 | Non-malignant tissue | Non-malignant |
| GSM5599223 | Norm4 | Non-malignant tissue | Non-malignant |
| GSM5599224 | Norm5 | Non-malignant tissue | Non-malignant |
| GSM5599225 | Cancer1 | HGSOC tissue | Advanced |
| GSM5599226 | Cancer2 | HGSOC tissue | Advanced |
| GSM5599227 | Cancer3 | HGSOC tissue | Early |
| GSM5599228 | Cancer4 | HGSOC tissue | Early |
| GSM5599229 | Cancer5 | HGSOC tissue | Advanced |
| GSM5599230 | Cancer6 | HGSOC tissue | Advanced |
| GSM5599231 | Cancer7 | HGSOC tissue | Early |


## Download and Directory Structure

1. Download the authors' processed 10x Genomics matrices from the GSE184880
   accession page.
2. Extract all 12 sample directories.
3. Place them under `data/GSE184880/` using the following structure:

```text
data/
└── GSE184880/
    ├── GSM5599220_Norm1/
    ├── GSM5599221_Norm2/
    ├── GSM5599222_Norm3/
    ├── GSM5599223_Norm4/
    ├── GSM5599224_Norm5/
    ├── GSM5599225_Cancer1/
    ├── GSM5599226_Cancer2/
    ├── GSM5599227_Cancer3/
    ├── GSM5599228_Cancer4/
    ├── GSM5599229_Cancer5/
    ├── GSM5599230_Cancer6/
    └── GSM5599231_Cancer7/
```

Each sample directory should contain the corresponding `matrix.mtx`,
`barcodes.tsv`, and `features.tsv` files, optionally gzip-compressed.

If the matrices are stored outside the repository, define their location before
running Stage 01:

```r
Sys.setenv(SCRNA_DATA_DIR = "path/to/GSE184880")
```

## Data Handling and Storage Notice

- Raw or processed expression matrices and computed Seurat objects (`.rds`)
  are excluded from Git tracking through `.gitignore`to maintain efficient repository size. 
- .rds, .h5, .h5ad, .mtx, and other large data files are not committed to the repository.
- Processed summary statistics, marker tables, quality-control reports, and
  other compact outputs are written to `results/`.
- High-resolution analysis figures are written to `figures/`.

