# Integrated Multi-Omics Analysis of Zika and Dengue Infection

## Scientific Overview
This repository contains a comprehensive, reproducible computational biology pipeline for the integrated analysis of Zika (ZIKV) and Dengue (DENV) virus infections. By synthesizing single-cell RNA-sequencing (scRNA-seq), bulk RNA-sequencing, and miRNA regulatory data, we investigate the conserved and virus-specific host responses to flavivirus infection.

## Research Objective
The primary objective is to identify shared host regulatory modules, key driving pathways, and potential therapeutic targets by triangulating evidence across multiple omics layers and independent experimental models (e.g., hNPCs, Huh7, and moDCs).

## Multi-Omics Strategy
1. **Bulk RNA-seq**: Establishes global transcriptomic shifts and high-confidence differentially expressed genes (DEGs).
2. **scRNA-seq**: Resolves cell-type-specific responses, infection trajectories, and bystander effects at single-cell resolution.
3. **miRNA Analysis**: Integrates post-transcriptional regulatory networks to identify upstream modulators of the conserved transcriptomic signature.
4. **Integration**: Cross-references these layers using Gene Set Variation Analysis (GSVA) and network topology to prioritize biological candidates.

## Repository Structure
`
.
|-- config/                  # Centralized configuration parameters
|-- data/                    # Raw and intermediate data (not tracked in git)
|-- docs/                    # Detailed scientific documentation and methods
|-- figures/                 # Generated publication-quality figures
|-- logs/                    # Execution logs
|-- results/                 # Statistical outputs and tables
|-- scripts/
|   |-- 01_Bulk_RNAseq/      # Steps 1-9 for bulk transcriptomics
|   |-- 02_Single_Cell_RNAseq/ # Phases A-F for scRNA-seq processing
|   |-- 03_GSVA_Analysis/    # Pathway enrichment and integration
|   |-- 04_Figure_Generation/ # 600 DPI publication figure scripts
|   |-- utils/               # Shared functions and visualization themes
|-- tables/                  # Supplementary and final result tables
`

## Reproducibility
This pipeline is designed for independent reproduction. See docs/REPRODUCIBILITY.md for detailed instructions on software environments, execution order, and required input data.
