# Reproducibility Guide

## Required Software
- **R**: >= 4.1.0
- **R Packages**: Seurat, DESeq2, GSVA, clusterProfiler, ggplot2, tidyverse.
- **Python**: >= 3.8 (for formatting and supplementary scripts).

## Execution Order
1. Execute scripts/01_Bulk_RNAseq/step1_raw_counts.R through step9.
2. Execute scripts/02_Single_Cell_RNAseq/PhaseA_preprocessing.R through PhaseF.
3. Execute scripts/03_GSVA_Analysis/ scripts.
4. Execute scripts/04_Figure_Generation/ scripts to generate final publication panels.

## Known Limitations
- Raw FASTQ/Count matrices are assumed to be present in data/. Due to size limits, they must be downloaded from GEO independently.
- Manual intervention may be required if system-specific memory limits are reached during Seurat integration (Phase D).
