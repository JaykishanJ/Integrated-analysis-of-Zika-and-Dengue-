# Computational Methods

## 1. Bulk RNA-seq Analysis
- **Data Source**: High-throughput sequencing datasets (e.g., GSE110512).
- **Preprocessing**: Raw counts were quantified and filtered for low-expression genes.
- **Statistical Model**: DESeq2 was used for differential expression, modeling condition while accounting for appropriate covariates (batch, tissue).
- **Multiple-Testing**: Benjamini-Hochberg (FDR < 0.05).

## 2. Single-Cell RNA-seq Analysis
- **Data Source**: ViscRNA-seq data (e.g., GSE110496).
- **QC**: Cells were filtered based on mitochondrial fraction, ERCC spike-ins, and detected features.
- **Normalization**: SCTransform / LogNormalize depending on the phase.
- **Integration**: Harmony was applied where batch integration was biologically required.
- **Differential Expression**: Wilcoxon Rank Sum test with Bonferroni correction for cell-state markers.

## 3. GSVA Pathway Analysis
- **Method**: Gene Set Variation Analysis (GSVA) using MSigDB hallmark, KEGG, and GO sets.
- **Biological Rationale**: Transforms gene-level data to pathway-level enrichment scores, enabling robust cross-dataset comparisons independent of exact gene-level noise.

## 4. Multi-Omics Integration and Candidate Prioritization
- **Method**: Hub genes were prioritized based on network centrality (PPI networks), consistent differential expression across datasets, and regulatory targeting by shared miRNAs.
- **Output**: A prioritized list of therapeutic targets (Drug Repurposing Phase F).
