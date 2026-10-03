# Figure 6: Integrating Bulk and Single-Cell Data (GSVA Signature Overlap)
source("utils.R")
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR <- "."
BULK_DIR <- file.path(BASE_DIR, "ZIKA_Bulk_Som", "ZIKA_Bulk_Som")
GSVA_DIR <- file.path(BASE_DIR, "output_manuscript_figures")
OUT_DIR  <- file.path(BASE_DIR, "publication_figures_scripts", "output")
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 6A. GSVA Enrichment Boxplot — NES across bulk datasets (real data)
# ---------------------------------------------------------------------------
cat("[Fig6] Loading GSVA integration data...\n")

nes_path <- file.path(GSVA_DIR, "GSVA_sample_scores_all_datasets.csv")
if (!file.exists(nes_path)) {
  stop("GSVA sample scores not found at ", nes_path,
       ". Run the GSVA bulk-SC integration script first.")
}

nes_df <- read.csv(nes_path)

# Verify columns
if (!all(c("dataset", "sample", "condition", "term", "NES") %in% names(nes_df))) {
  stop("GSVA CSV missing expected columns. Got: ", paste(names(nes_df), collapse = ", "))
}

# Extract condition from sample name if not already present
if (!"condition" %in% names(nes_df) || all(is.na(nes_df$condition))) {
  nes_df$condition <- "Unknown"
  # Try to infer Infected vs Mock from dataset column (e.g. "Dengue" = infected)
}

# ---------------------------------------------------------------------------
# 6A. Boxplot: NES by condition, faceted by dataset
# ---------------------------------------------------------------------------
# Focus on the "Verified_23_All" core gene set
nes_23 <- filter(nes_df, term == "Verified_23_All")

# Reconstruct infection state from dataset + condition
# Known mappings: Dengue=Infected, ZIKVC=Infected, ZIKVM=Infected
nes_23$Infection <- factor(nes_23$condition,
                           levels = c("Mock", "Dengue", "ZIKVC", "ZIKVM",
                                      "Bystander", "Zika_bystander"),
                           labels = c("Mock", "Infected", "Infected", "Infected",
                                      "Infected", "Infected"))

nes_23$Infection <- factor(ifelse(nes_23$Infection == "Mock", "Mock", "Infected"),
                           levels = c("Mock", "Infected"))

# Clean dataset names for display
nes_23$Dataset <- gsub("GSE110512 \\(Huh7\\): DENV", "GSE110512 (DENV)", nes_23$dataset)
nes_23$Dataset <- gsub("GSE110512: ", "", nes_23$Dataset)

p_nes_box <- ggplot(nes_23, aes(x = Infection, y = NES, fill = Infection)) +
  geom_boxplot(alpha = 0.7, color = "black", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               fill = "white", color = "black", stroke = 1) +
  scale_fill_manual(values = c("Mock" = "#377eb8", "Infected" = "#e41a1c")) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 4) +
  labs(
    x = "Condition",
    y = "GSVA NES (Verified 23 Core Genes)",
    title = "Bulk–SC Consensus Signature Enrichment Across Datasets"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# 6B. NES Shift Barplot (Nested) — real contrast differential data
# ---------------------------------------------------------------------------
contrast_path <- file.path(GSVA_DIR, "GSVA_contrast_differential_analysis.csv")
if (file.exists(contrast_path)) {
  contrast_df <- read.csv(contrast_path)

  # Keep "Verified_23_All" across all datasets
  contrast_23 <- filter(contrast_df, term == "Verified_23_All") %>%
    mutate(
      Direction = ifelse(NES_shift > 0, "Up in Infected", "Down in Infected"),
      Dataset = gsub("GSE110512 \\(Huh7\\): DENV", "GSE110512 (DENV)", dataset),
      Dataset = gsub("GSE110512: ", "", Dataset)
    )

  p_shift <- ggplot(contrast_23, aes(x = Dataset, y = NES_shift, fill = Direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Up in Infected" = "#e41a1c",
                                  "Down in Infected" = "#377eb8")) +
    labs(
      x = "Dataset",
      y = "NES Shift (Infected - Mock)",
      title = "GSVA Signature Shift (Verified 23 Core Genes)"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
} else {
  warning("GSVA contrast differential CSV not found — skipping shift plot")
  p_shift <- ggplot() + labs(title = "GSVA Shift (data unavailable)") + theme_publication()
}

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
fig6_layout <- (p_nes_box / p_shift) + plot_annotation(tag_levels = 'A')
fig6_layout <- fig6_layout & theme(text = element_text(family = "sans", color = "black"))

ggsave(file.path(OUT_DIR, "Figure6_Integration.png"), plot = fig6_layout,
       width = FIG_WIDTH_2COL, height = 9, dpi = PUB_DPI)
ggsave(file.path(OUT_DIR, "Figure6_Integration.pdf"), plot = fig6_layout,
       width = FIG_WIDTH_2COL, height = 9, dpi = PUB_DPI)
ggsave(file.path(OUT_DIR, "Figure6_Integration.tiff"), plot = fig6_layout,
       width = FIG_WIDTH_2COL, height = 9, dpi = PUB_DPI, device = "tiff", compression = "lzw")
ggsave(file.path(OUT_DIR, "Figure6_Integration.svg"), plot = fig6_layout,
       width = FIG_WIDTH_2COL, height = 9, dpi = PUB_DPI, device = "svg")
cat("Figure 6 (Integration) generated and saved successfully.\n")
