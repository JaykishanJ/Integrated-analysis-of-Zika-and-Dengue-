# Figure 3: Cell-Type-Specific Transcriptional Responses (ZIKV vs Mock)
source("utils.R")
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR <- "."
RES_DIR  <- file.path(BASE_DIR, "results")
OUT_DIR  <- file.path(BASE_DIR, "publication_figures_scripts", "output")
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 3A. Gene Selection Funnel
# ---------------------------------------------------------------------------
cat("[Fig3] Creating Gene Selection Funnel...\n")
funnel_data <- data.frame(
  Stage = factor(c("1. Bulk DENV DEGs", "1. Bulk ZIKV DEGs", "2. Shared Core", 
                   "3. scRNA-seq DENV", "3. scRNA-seq ZIKV", "4. Final Hub Genes"),
                 levels = rev(c("1. Bulk DENV DEGs", "1. Bulk ZIKV DEGs", "2. Shared Core", 
                            "3. scRNA-seq DENV", "3. scRNA-seq ZIKV", "4. Final Hub Genes"))),
  Count = c(260, 207, 23, 138, 25, 8)
)

p_funnel <- ggplot(funnel_data, aes(x = Count, y = Stage)) +
  geom_bar(stat = "identity", fill = "#4daf4a", color = "black", width = 0.6) +
  geom_text(aes(label = Count), hjust = -0.2, fontface = "bold", size = 4) +
  labs(x = "Number of Genes", y = "Filtering Stage", title = "Gene Selection Pipeline") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_publication()

# ---------------------------------------------------------------------------
# 3B. Volcano — Differential expression in ZIKV-Infected cells vs Mock
# ---------------------------------------------------------------------------
# Use the Step11d all-timepoint DEG table with real per-cell-type DGEs
# from the cluster-specific DE table (Step11a within-plate Wilcoxon)
cat("[Fig3] Loading real DEG data...\n")

# Primary: ZIKV timepoint DEG table from phase D (Step11d)
# This has per-timepoint DEGs: zika_4h, zika_12h, zika_24h, zika_48h
# (plus dengue timepoints for comparison)
deg_table_path <- file.path(RES_DIR, "phaseD/tables/Step11d_All_Timepoint_DEGs.csv")
if (!file.exists(deg_table_path)) {
  stop("DEG table not found at ", deg_table_path)
}

volc_cluster <- read.csv(deg_table_path)
cat("[Fig3] DEG table loaded:", nrow(volc_cluster), "rows\n")
cat("[Fig3] Columns:", paste(names(volc_cluster), collapse = ", "), "\n")

# Identify the correct column names
logfc_col   <- intersect(c("avg_log2FC", "log2FoldChange", "avg_lfc"), names(volc_cluster))[1]
padj_col    <- intersect(c("p_val_adj", "p_adj", "padj", "fdr"), names(volc_cluster))[1]
gene_col    <- intersect(c("symbol", "gene", "Gene"), names(volc_cluster))[1]
cluster_col <- intersect(c("cluster", "Cluster", "cell_type"), names(volc_cluster))[1]

cat(sprintf("[Fig3] Column mapping: logfc=%s, padj=%s, gene=%s, cluster=%s\n",
            logfc_col, padj_col, gene_col, cluster_col))

if (is.na(logfc_col) || is.na(padj_col)) {
  stop("Required columns (log2FC, padj) not found in DEG table. Got: ",
       paste(names(volc_cluster), collapse = ", "))
}

# Handle NAs
volc_cluster[[padj_col]][is.na(volc_cluster[[padj_col]])] <- 1
volc_cluster[[logfc_col]][is.na(volc_cluster[[logfc_col]])] <- 0

# Filter for ZIKV conditions
if ("condition" %in% names(volc_cluster)) {
  cat("[Fig3] Unique conditions:\n")
  print(sort(unique(volc_cluster$condition)))

  # Select ZIKV conditions (24h peak response by default, fallback to first zika timepoint)
  zika_conds <- grep("^zika|ZIKV|zika_", unique(volc_cluster$condition), value = TRUE)
  if (length(zika_conds) == 0) {
    zika_conds <- grep("zika", unique(volc_cluster$condition), value = TRUE, ignore.case = TRUE)
  }

  if (length(zika_conds) == 0) {
    warning("No ZIKV conditions found, using condition with most significant genes")
    cond_counts <- volc_cluster %>%
      group_by(condition) %>%
      summarise(n_sig = sum(.data[[padj_col]] < 0.05 & abs(.data[[logfc_col]]) > 0.5),
                .groups = "drop") %>%
      arrange(desc(n_sig))
    target_condition <- cond_counts$condition[1]
  } else {
    # Prefer 24h timepoint (peak response shown in Figure 2)
    target_condition <- if (any(grepl("24", zika_conds))) {
      grep("24", zika_conds, value = TRUE)[1]
    } else {
      zika_conds[1]
    }
  }

  cat(sprintf("[Fig3] Using ZIKV condition: %s (from: %s)\n",
              target_condition, paste(zika_conds, collapse = ", ")))

  volc_df <- volc_cluster %>% filter(condition == target_condition)
} else {
  volc_df <- volc_cluster
  target_condition <- "ZIKV vs Mock"
}

volc_df$sig <- ifelse(volc_df[[padj_col]] < 0.05 & abs(volc_df[[logfc_col]]) > 0.5,
                      "Significant", "Not Significant")
top_genes_df <- volc_df %>% arrange(.data[[padj_col]]) %>% slice_head(n = 10)

# Plot volcano
p_volcano_cell <- ggplot(volc_df,
                         aes(x = .data[[logfc_col]],
                             y = -log10(.data[[padj_col]]),
                             color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "#e41a1c",
                                "Not Significant" = "grey80")) +
  ggrepel::geom_text_repel(
    data = top_genes_df,
    aes(label = .data[[gene_col]]),
    size = 3.5,
    max.overlaps = 15,
    color = "black"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  labs(
    title = sprintf("DEGs (%s)", target_condition),
    x = bquote(~Log[2]~"Fold Change"),
    y = bquote(~-Log[10]~"Adjusted P-value")
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# 3B. Heatmap — Top DEGs across clusters (real cluster signature scores)
# ---------------------------------------------------------------------------
# Use cluster signature scores from Phase_D_ZIKV_Sig_Matrix.rds
# This is a real matrix: clusters × genes, or genes × clusters
cat("[Fig3] Loading cluster signature matrix...\n")

sig_mat_path <- file.path(RES_DIR, "phaseD/objects/Phase_D_ZIKV_Sig_Matrix.rds")
if (!file.exists(sig_mat_path)) {
  stop("Cluster signature matrix not found at ", sig_mat_path)
}

sig_mat <- readRDS(sig_mat_path)

# Coerce to numeric matrix (cluster x gene or gene x cluster)
heat_mat <- as.matrix(sig_mat)
storage.mode(heat_mat) <- "numeric"

# Standardize across clusters so colors are comparable
heat_mat <- t(scale(t(heat_mat)))

# Limit to top 30 most-variable genes
gene_var <- apply(heat_mat, 1, var, na.rm = TRUE)
heat_mat <- heat_mat[order(-gene_var)[1:min(30, length(gene_var))], , drop = FALSE]

# Ensure column names are clusters
if (is.null(colnames(heat_mat))) {
  colnames(heat_mat) <- paste0("Cluster_", seq_len(ncol(heat_mat)))
}

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

ht <- Heatmap(
  heat_mat,
  name = "Z-Score",
  col = col_fun,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Cluster-Type-Specific ZIKV Infection Signatures",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title = paste0("Top DEGs (n=", nrow(heat_mat), ")"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title = "Z-Score",
    title_gp = gpar(fontface = "bold", fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
# Heatmap goes directly to PDF (ComplexHeatmap doesn't play well with ggsave)
pdf(file.path(OUT_DIR, "Figure3_Heatmap.pdf"), width = FIG_WIDTH_1.5COL, height = 8)
draw(ht)
dev.off()

tiff(file.path(OUT_DIR, "Figure3_Heatmap.tiff"), width = FIG_WIDTH_1.5COL, height = 8, units = "in", res = PUB_DPI, compression = "lzw")
draw(ht)
dev.off()
cat("[Fig3] Heatmap PDF and TIFF saved.\n")

# Funnel + Volcano as PNG/PDF
fig3_layout <- (p_funnel / p_volcano_cell) + plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = 'A')

ggsave(file.path(OUT_DIR, "Figure3_CellSpecific.png"), plot = fig3_layout,
       width = FIG_WIDTH_1COL + 2, height = 8, dpi = PUB_DPI)
ggsave(file.path(OUT_DIR, "Figure3_CellSpecific.pdf"), plot = fig3_layout,
       width = FIG_WIDTH_1COL + 2, height = 8, dpi = PUB_DPI)
ggsave(file.path(OUT_DIR, "Figure3_CellSpecific.tiff"), plot = fig3_layout,
       width = FIG_WIDTH_1COL + 2, height = 8, dpi = PUB_DPI, device = "tiff", compression = "lzw")
ggsave(file.path(OUT_DIR, "Figure3_CellSpecific.svg"), plot = fig3_layout,
       width = FIG_WIDTH_1COL + 2, height = 8, dpi = PUB_DPI, device = "svg")
cat("Figure 3 generated and saved successfully.\n")
