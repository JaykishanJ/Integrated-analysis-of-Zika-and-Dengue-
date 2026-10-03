# Figure 2: Single-Cell Atlas Overview
source("utils.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR <- "."
RES_DIR  <- file.path(BASE_DIR, "results")
OUT_DIR  <- file.path(BASE_DIR, "publication_figures_scripts", "output")
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load real data
# ---------------------------------------------------------------------------
cat("[Fig2] Loading single-cell objects...\n")
qc_cells  <- readRDS(file.path(RES_DIR, "phaseB/objects/Phase_B_cell_qc.rds"))
clust_obj <- readRDS(file.path(RES_DIR, "phaseD/objects/Phase_D_clustering.rds"))
log_mat   <- readRDS(file.path(RES_DIR, "phaseC/objects/Phase_C_log_matrix.rds"))

# Extract UMAP coordinates
umap_coords <- as.data.frame(clust_obj$embedding)
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Combine with metadata
umap_data <- cbind(umap_coords, qc_cells[rownames(umap_coords), ])

# Safely extract Cluster and Condition
if ("clusters" %in% names(clust_obj)) {
  umap_data$CellType <- as.character(clust_obj$clusters)
} else if ("seurat_clusters" %in% names(clust_obj)) {
  umap_data$CellType <- as.character(clust_obj$seurat_clusters)
} else {
  umap_data$CellType <- "Unknown"
}

if ("infection_class" %in% colnames(umap_data)) {
  umap_data$Condition <- as.character(umap_data$infection_class)
} else if ("virus" %in% colnames(umap_data)) {
  umap_data$Condition <- as.character(umap_data$virus)
} else {
  umap_data$Condition <- "Unknown"
}

# ---------------------------------------------------------------------------
# 2C. Marker Gene DotPlot — computed from real log-normalized matrix
# ---------------------------------------------------------------------------
cat("[Fig2] Computing marker gene expression from real data...\n")

cluster_markers_path <- file.path(RES_DIR, "phaseE/tables/Step17_TopClusterMarkers.csv")
if (file.exists(cluster_markers_path)) {
  # Use pre-computed cluster markers from phaseE
  marker_table <- read.csv(cluster_markers_path)

  # Get top markers per cluster (by avg_log2FC, significant only)
  top_markers <- marker_table %>%
    filter(significant == TRUE) %>%
    group_by(cluster1) %>%
    slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
    ungroup()

  # Extract gene symbols
  marker_genes <- unique(top_markers$symbol)
  cat(sprintf("[Fig2] Using %d pre-computed marker genes from %d clusters\n",
              length(marker_genes), length(unique(top_markers$cluster1))))

} else {
  # Fallback: compute markers from log matrix using Seurat-style logic
  warning("Step17_TopClusterMarkers.csv not found, computing markers from log matrix")

  gene_vars <- apply(log_mat, 1, var, na.rm = TRUE)
  hv_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]

  # Compute mean expression per cluster
  clusters_vec <- umap_data$CellType[match(colnames(log_mat), rownames(umap_data))]
  marker_expr <- sapply(split(seq_len(ncol(log_mat)), clusters_vec), function(idx) {
    rowMeans(log_mat[hv_genes, idx, drop = FALSE], na.rm = TRUE)
  })

  marker_genes <- rownames(marker_expr)
}

# Compute actual dotplot data from log matrix
# For each marker gene, compute avg expression and % expressing per cluster
marker_genes_valid <- intersect(marker_genes, rownames(log_mat))

if (length(marker_genes_valid) < 3) {
  warning("Too few valid marker genes in log matrix. Falling back to top variable genes.")
  gene_vars <- apply(log_mat, 1, var, na.rm = TRUE)
  marker_genes_valid <- names(sort(gene_vars, decreasing = TRUE))[1:15]
}

clusters_present <- sort(unique(umap_data$CellType))
cell_clusters_map <- umap_data$CellType[match(colnames(log_mat), rownames(umap_data))]

dotplot_data <- data.frame()
for (gene in marker_genes_valid) {
  gene_expr <- log_mat[gene, ]

  for (clust in clusters_present) {
    cell_idx <- which(cell_clusters_map == clust)
    if (length(cell_idx) == 0) next

    expr_vals <- as.numeric(gene_expr[cell_idx])
    avg_exp   <- mean(expr_vals, na.rm = TRUE)
    pct_exp   <- mean(expr_vals > 0, na.rm = TRUE) * 100

    # Scale avg expression relative to overall mean for visualization
    overall_mean <- mean(as.numeric(gene_expr), na.rm = TRUE)
    overall_sd   <- sd(as.numeric(gene_expr), na.rm = TRUE)
    avg_scaled   <- ifelse(overall_sd > 0, (avg_exp - overall_mean) / overall_sd, 0)

    dotplot_data <- rbind(dotplot_data, data.frame(
      Features.plot = gene,
      id = paste0("Cluster ", clust),
      avg.exp.scaled = avg_scaled,
      pct.exp = pct_exp
    ))
  }
}

dotplot_data$Features.plot <- factor(dotplot_data$Features.plot,
                                     levels = unique(marker_genes_valid))

# Limit to top 15 genes for readability
if (length(marker_genes_valid) > 15) {
  # Keep genes with highest variance in expression across clusters
  gene_var <- dotplot_data %>%
    group_by(Features.plot) %>%
    summarise(var_scaled = var(avg.exp.scaled), .groups = "drop") %>%
    arrange(desc(var_scaled))
  keep_genes <- head(gene_var$Features.plot, 15)
  dotplot_data <- filter(dotplot_data, Features.plot %in% keep_genes)
  dotplot_data$Features.plot <- factor(dotplot_data$Features.plot, levels = keep_genes)
}

# ---------------------------------------------------------------------------
# 2A. Master DEG Bar Chart
# ---------------------------------------------------------------------------
cat("[Fig2] Creating Master DEG Bar Chart...\n")
deg_summary <- data.frame(
  Contrast = c("hNPC DENV-2", "hNPC ZIKV-M", "hNPC ZIKV-C", "Huh7 DENV-2", "Huh7.5 ZIKV"),
  Upregulated = c(3251, 627, 99, 613, 160),
  Downregulated = c(3102, 842, 32, 203, 2)
)
deg_long <- pivot_longer(deg_summary, cols = c(Upregulated, Downregulated), 
                         names_to = "Direction", values_to = "Count")
# Make counts negative for Downregulated to create a bidirectional bar chart
deg_long$PlotCount <- ifelse(deg_long$Direction == "Downregulated", -deg_long$Count, deg_long$Count)
deg_long$Contrast <- factor(deg_long$Contrast, levels = deg_summary$Contrast)

p_deg_bar <- ggplot(deg_long, aes(x = Contrast, y = PlotCount, fill = Direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Upregulated" = "#e41a1c", "Downregulated" = "#377eb8")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(x = "Contrast (Infected vs Mock)", y = "Number of DEGs", title = "Master DEG Summary") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.title = element_blank()) +
  scale_y_continuous(labels = abs)

# ---------------------------------------------------------------------------
# 2B. Global UMAP colored by Cell Type
# ---------------------------------------------------------------------------
p_umap_ct <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Single-Cell Atlas", x = "UMAP 1", y = "UMAP 2") +
  theme_publication() +
  theme(legend.position = "right")

# ---------------------------------------------------------------------------
# 2C. UMAP Split by Condition
# ---------------------------------------------------------------------------
p_umap_cond <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Condition)) +
  geom_point(size = 0.5, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = col_palette_condition) +
  facet_wrap(~Condition) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "By Condition", x = "UMAP 1", y = "UMAP 2") +
  theme_publication() +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# 2C. Marker Gene DotPlot
# ---------------------------------------------------------------------------
p_dotplot <- ggplot(dotplot_data,
                    aes(x = Features.plot, y = id,
                        size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradientn(
    colors = c("#2166AC", "#F7F7F7", "#B2182B"),
    name = "Scaled\nExpression"
  ) +
  scale_size_continuous(range = c(1, 6), name = "% Expressing") +
  labs(x = "Marker Gene", y = "Cell Type",
       title = "Canonical Marker Expression") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
    axis.text.y = element_text(size = 9)
  )

# ---------------------------------------------------------------------------
# 2E. Compositional Analysis
# ---------------------------------------------------------------------------
comp_data <- umap_data %>%
  group_by(Condition, CellType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Condition) %>%
  mutate(Proportion = Count / sum(Count))

p_comp <- ggplot(comp_data, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  labs(y = "Relative Proportion", title = "Cellular Composition") +
  theme_publication()

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
fig2_top <- p_deg_bar + p_umap_ct + p_umap_cond + plot_layout(widths = c(1, 1, 1.2))
fig2_bottom <- p_dotplot + p_comp + plot_layout(widths = c(1.5, 1))
fig2_layout <- (fig2_top / fig2_bottom) + plot_annotation(tag_levels = 'A')

ggsave(file.path(OUT_DIR, "Figure2_SingleCell.png"), plot = fig2_layout,
       width = FIG_WIDTH_2COL + 2, height = 10, dpi = PUB_DPI)
ggsave(file.path(OUT_DIR, "Figure2_SingleCell.pdf"), plot = fig2_layout,
       width = FIG_WIDTH_2COL + 2, height = 10, dpi = PUB_DPI)
ggsave(file.path(OUT_DIR, "Figure2_SingleCell.tiff"), plot = fig2_layout,
       width = FIG_WIDTH_2COL + 2, height = 10, dpi = PUB_DPI, device = "tiff", compression = "lzw")
ggsave(file.path(OUT_DIR, "Figure2_SingleCell.svg"), plot = fig2_layout,
       width = FIG_WIDTH_2COL + 2, height = 10, dpi = PUB_DPI, device = "svg")
cat("Figure 2 (scRNA-seq Overview) generated and saved successfully.\n")
