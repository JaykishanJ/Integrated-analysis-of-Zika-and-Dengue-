<div align="center">
  <img src="https://raw.githubusercontent.com/JaykishanJ/Integrated-analysis-of-Zika-and-Dengue-/main/figures/workflow.png" alt="Integrated Multi-Omics Pipeline" width="800" onerror="this.style.display='none'"/>
  
  # Integrated Multi-Omics Analysis of Zika and Dengue Infection
  
  **A Comprehensive, Reproducible Computational Biology Pipeline**

  [![R Version](https://img.shields.io/badge/R-%3E%3D%204.1.0-blue?logo=R)](https://www.r-project.org/)
  [![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
  [![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
  [![Publication](https://img.shields.io/badge/Publication-Pending-orange)](#)
  [![Reproducible](https://img.shields.io/badge/Reproducible-Yes-success)](#reproducibility)
</div>

<br/>

> **Abstract:** By synthesizing single-cell RNA-sequencing (scRNA-seq), bulk RNA-sequencing, and miRNA regulatory data, we systematically investigate the conserved and virus-specific host responses to flavivirus infection (ZIKV & DENV). This repository contains the complete analytical workflow from raw count processing to final publication-ready 600 DPI figure generation.

<br/>

## 📑 Table of Contents
- [🎯 Research Objective](#-research-objective)
- [🧬 Multi-Omics Strategy](#-multi-omics-strategy)
- [📂 Repository Structure](#-repository-structure)
- [⚙️ Reproducibility & Execution](#️-reproducibility--execution)
- [📚 Documentation Directory](#-documentation-directory)
- [🤝 Contributing](#-contributing)

---

## 🎯 Research Objective
The primary objective of this project is to identify shared host regulatory modules, key driving pathways, and potential therapeutic targets by triangulating computational evidence across:
- **Multiple Omics Layers:** Transcriptomics, Single-Cell resolving, and Regulatory miRNAs.
- **Independent Experimental Models:** Human Neural Progenitor Cells (hNPCs), Hepatocytes (Huh7), and Monocyte-derived Dendritic Cells (moDCs).

---

## 🧬 Multi-Omics Strategy
Our unified approach breaks down into four sequential stages:

1. 📊 **Bulk RNA-seq Analysis**  
   Establishes global transcriptomic shifts and generates high-confidence differentially expressed genes (DEGs) while accounting for experimental covariates.
2. 🔬 **scRNA-seq Analysis**  
   Resolves cell-type-specific responses, charts infection trajectories, and separates true cellular responses from bystander effects at single-cell resolution.
3. 🕸️ **miRNA Regulatory Analysis**  
   Integrates post-transcriptional regulatory networks to identify upstream modulators driving the conserved transcriptomic signatures.
4. 🔗 **Integration & Prioritization**  
   Cross-references evidence using Gene Set Variation Analysis (GSVA) and Protein-Protein Interaction (PPI) topology to rank biological candidates for drug repurposing.

---

## 📂 Repository Structure

The workflow is completely modularized for ease of execution and review.

```text
📦 Integrated-analysis-of-Zika-and-Dengue
 ┣ 📂 config/                 # Centralized configuration parameters
 ┣ 📂 data/                   # Raw matrices & metadata (add locally)
 ┣ 📂 docs/                   # Scientific methodology and dictionaries
 ┣ 📂 figures/                # Output dir for SVGs/TIFFs
 ┣ 📂 logs/                   # Analytical execution logs
 ┣ 📂 results/                # Statistical outputs, RDS objects, tables
 ┣ 📂 scripts/
 ┃ ┣ 📂 01_Bulk_RNAseq/       # Step-by-step DESeq2 bulk pipeline
 ┃ ┣ 📂 02_Single_Cell_RNAseq/# Phase A-F Seurat scRNA-seq workflow
 ┃ ┣ 📂 03_GSVA_Analysis/     # Pathway scoring and integration logic
 ┃ ┣ 📂 04_Figure_Generation/ # Final 600 DPI publication figure code
 ┃ ┗ 📂 utils/                # Shared visualization themes (inject_themes.R)
 ┣ 📂 tables/                 # Final aggregated Excel/CSV tables
 ┣ 📜 CHANGELOG.md            # Version history and updates
 ┣ 📜 pipeline_manifest.yaml  # Complete snapshot of expected outputs
 ┗ 📜 README.md               # You are here
```

<details>
<summary><b>Click here to view detailed script order</b></summary>

- **Bulk**: `step1_raw_counts.R` ➔ `step9_diagnostics.R`
- **Single-Cell**: `PhaseA_preprocessing.R` ➔ `PhaseF_Drug_repurpusing.R`
- **GSVA**: `assemble_figure5_gsva.R` / `run_true_pathway_gsva.py`
- **Figures**: `Figure1_BulkRNAseq.R` ➔ `Figure6_Integration.R`
</details>

---

## ⚙️ Reproducibility & Execution

This pipeline is engineered for independent reproduction. We strongly enforce robust QC and strict dependency management.

### 🚀 Quick Start
1. Clone the repository:
   ```bash
   git clone https://github.com/JaykishanJ/Integrated-analysis-of-Zika-and-Dengue-.git
   cd Integrated-analysis-of-Zika-and-Dengue-
   ```
2. Download the required raw data matrices from GEO (e.g., `GSE110496`, `GSE110512`) into the `data/` directory.
3. Execute the pipeline sequentially, starting from `scripts/01_Bulk_RNAseq/`.

> [!NOTE]  
> Please refer to [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for detailed instructions on software environments, execution order, and required inputs.

---

## 📚 Documentation Directory

We have provided extensive documentation so any computational biologist can instantly understand our scientific rationale and thresholds.

| Document | Description |
|---|---|
| 📖 [**METHODS**](docs/METHODS.md) | Exhaustive scientific methodology detailing statistical models and QC cutoffs. |
| 🔄 [**REPRODUCIBILITY**](docs/REPRODUCIBILITY.md) | Software versions, environment setup, and execution order. |
| 📓 [**DATA DICTIONARY**](docs/DATA_DICTIONARY.md) | Glossary of sample traits, metadata fields, and dataset identifiers. |
| 🎨 [**FIGURE STYLE**](docs/FIGURE_STYLE.md) | Standardized typography, colors, and layout guidelines for visuals. |
| 📊 [**ANALYSIS SUMMARY**](docs/ANALYSIS_SUMMARY.md) | Plain-english translation of computational evidence into biological hypotheses. |
| ⚠️ [**LIMITATIONS**](docs/LIMITATIONS.md) | Scientific caveats, computational boundaries, and experimental caveats. |

<br/>

<div align="center">
  <sub>Built with ❤️ for Reproducible Bioinformatics.</sub>
</div>
