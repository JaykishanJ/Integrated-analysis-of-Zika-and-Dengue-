# Figure 4: Cellular Trajectories and Dynamics (Monocle3)
source("utils.R")
library(ggplot2)
library(patchwork)
library(dplyr)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR <- "."
RES_DIR  <- file.path(BASE_DIR, "results")
OUT_DIR  <- file.path(BASE_DIR, "publication_figures_scripts", "output")
CACHE_DIR <- file.path(RES_DIR, "phaseD", "cache", "trajectory")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)

TRAJ_CACHE <- file.path(CACHE_DIR, "monocle3_trajectory.rds")
GENE_TREND_CACHE <- file.path(CACHE_DIR, "gene_trends.rds")

# ---------------------------------------------------------------------------
# Load required data
# ---------------------------------------------------------------------------
cat("[Fig4] Loading data objects...\n")

log_mat <- readRDS(file.path(RES_DIR, "phaseC/objects/Phase_C_log_matrix.rds"))
clust_obj  <- readRDS(file.path(RES_DIR, "phaseD/objects/Phase_D_clustering.rds"))
qc_cells   <- readRDS(file.path(RES_DIR, "phaseB/objects/Phase_B_cell_qc.rds"))

# Extract cluster labels
if ("clusters" %in% names(clust_obj)) {
  cell_clusters <- clust_obj$clusters
} else if ("seurat_clusters" %in% names(clust_obj)) {
  cell_clusters <- clust_obj$seurat_clusters
} else {
  stop("No cluster labels found in Phase_D_clustering.rds")
}
# Ensure cell_clusters has true cell IDs from the embedding
if (is.null(names(cell_clusters)) || !identical(names(cell_clusters), rownames(clust_obj$embedding))) {
  names(cell_clusters) <- rownames(clust_obj$embedding)
}

# Extract condition labels (and explicitly assign names)
if ("infection_class" %in% colnames(qc_cells)) {
  cell_condition <- as.character(qc_cells[names(cell_clusters), "infection_class"])
} else if ("virus" %in% colnames(qc_cells)) {
  cell_condition <- as.character(qc_cells[names(cell_clusters), "virus"])
} else {
  stop("No infection_class or virus column in Phase_B_cell_qc.rds")
}
names(cell_condition) <- names(cell_clusters)

# Identify root cluster (cluster 1 by default)
root_cluster_id <- 1
root_cells <- names(cell_clusters)[cell_clusters == root_cluster_id]
if (length(root_cells) == 0) {
  warning("Cluster ", root_cluster_id, " not found. Using lowest available cluster as root.")
  root_cluster_id <- sort(unique(cell_clusters))[1]
}

cat(sprintf("[Fig4] Using cluster %s as root (%d cells)\n",
             root_cluster_id, length(root_cells)))

# ---------------------------------------------------------------------------
# Monocle3 trajectory inference (with caching)
# ---------------------------------------------------------------------------
cat("[Fig4] Running trajectory inference...\n")

if (requireNamespace("monocle3", quietly = TRUE)) {
  library(monocle3)
} else {
  stop("monocle3 not installed. Run: BiocManager::install('monocle3')")
}

if (file.exists(TRAJ_CACHE)) {
  cat("[Fig4] Loading cached trajectory from", TRAJ_CACHE, "\n")
  cds <- readRDS(TRAJ_CACHE)
} else {
  # Build CellDataSet from log matrix
  # log_mat is genes x cells (rows = genes, cols = cells)
  feature_metadata <- data.frame(
    gene_short_name = rownames(log_mat),
    row.names = rownames(log_mat),
    stringsAsFactors = FALSE
  )
  cell_metadata <- data.frame(
    row.names = colnames(log_mat),
    cluster = factor(cell_clusters[colnames(log_mat)]),
    condition = factor(cell_condition[colnames(log_mat)]),
    stringsAsFactors = FALSE
  )
  expression_matrix <- as(as.matrix(log_mat), "sparseMatrix")

  cds <- new_cell_data_set(
    expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = feature_metadata
  )

  # Preprocess
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")

  # Cluster (reuse existing)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  
  # Inject Seurat clusters safely into Monocle3 UMAP clusters
  clusters_umap <- cell_metadata$cluster
  names(clusters_umap) <- rownames(cell_metadata)
  cds@clusters$UMAP$clusters <- clusters_umap
  
  # Learn trajectory, starting from root cluster
  cds <- learn_graph(cds, close_loop = FALSE)
  cds <- order_cells(cds, reduction_method = "UMAP",
                     root_cells = root_cells)
  
  # Save cache
  saveRDS(cds, TRAJ_CACHE)
  cat("[Fig4] Trajectory cached to", TRAJ_CACHE, "\n")
}

# Extract UMAP coordinates
umap_df <- as.data.frame(reducedDim(cds, "UMAP"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$CellID <- rownames(umap_df)
umap_df$Pseudotime <- pseudotime(cds, reduction_method = "UMAP")
umap_df$Cluster <- factor(cds@clusters$UMAP$clusters[rownames(umap_df)])
umap_df$Condition <- factor(cell_condition[umap_df$CellID])

# ---------------------------------------------------------------------------
# 4A. UMAP trajectory with pseudotime overlay
# ---------------------------------------------------------------------------
p_traj <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  labs(
    x = "UMAP 1", y = "UMAP 2",
    title = "Pseudotime Trajectory"
  ) +
  theme_publication(base_size = 11) +
  theme(legend.position = "right")

# Overlay trajectory graph on UMAP (simplified: cluster centroids + lineage)
centroids <- umap_df %>%
  group_by(Cluster) %>%
  summarise(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# 4B. Gene expression dynamics over pseudotime
# ---------------------------------------------------------------------------
cat("[Fig4] Computing gene expression trends...\n")

if (file.exists(GENE_TREND_CACHE)) {
  cat("[Fig4] Loading cached gene trends from", GENE_TREND_CACHE, "\n")
  trend_df <- readRDS(GENE_TREND_CACHE)
} else {
  # Select top DEGs by variance across pseudotime from the scaled matrix
  # Use all_timepoint DEGs table if available to pick biologically relevant genes
  deg_table_path <- file.path(RES_DIR,
    "phaseD/tables/Step11d_All_Timepoint_DEGs.csv")
  trend_genes <- c()

  if (file.exists(deg_table_path)) {
    top_deg <- tryCatch({
      read.csv(deg_table_path, nrows = 5000) %>%
        group_by(symbol) %>%
        summarise(max_lfc = max(abs(avg_log2FC), na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_lfc)) %>%
        slice_head(n = 8) %>%
        pull(symbol)
    }, error = function(e) character(0))

    trend_genes <- unique(c(top_deg, trend_genes))
  }

  # Fallback: top variable genes
  if (length(trend_genes) < 3) {
    gene_vars <- apply(log_mat, 1, var, na.rm = TRUE)
    trend_genes <- c(trend_genes, names(sort(gene_vars, decreasing = TRUE))[1:8])
  }
  trend_genes <- unique(head(trend_genes, 8))

  # Build long-format trend data frame
  valid_genes <- intersect(trend_genes, rownames(log_mat))
  if (length(valid_genes) == 0) {
    stop("No valid trend genes found in log matrix")
  }

  cells_with_pseudo <- names(pseudotime(cds, reduction_method = "UMAP"))
  trend_mat <- t(as.matrix(log_mat[valid_genes, cells_with_pseudo, drop = FALSE]))
  cell_pseudo <- pseudotime(cds, reduction_method = "UMAP")

  trend_df <- data.frame()
  for (g in valid_genes) {
    gene_vals <- as.numeric(trend_mat[, g])
    trend_df <- rbind(trend_df, data.frame(
      Gene = g,
      Pseudotime = as.numeric(cell_pseudo),
      Expression = gene_vals
    ))
  }
  trend_df$Gene <- factor(trend_df$Gene, levels = valid_genes)

  saveRDS(trend_df, GENE_TREND_CACHE)
  cat("[Fig4] Gene trends cached to", GENE_TREND_CACHE, "\n")
}

p_dyn <- ggplot(trend_df, aes(x = Pseudotime, y = Expression)) +
  geom_point(alpha = 0.25, color = "grey50", size = 0.8) +
  geom_smooth(
    aes(color = Gene),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    linewidth = 1,
    se = TRUE,
    alpha = 0.15
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Pseudotime",
    y = "Scaled Expression",
    title = "Expression Dynamics over Pseudotime"
  ) +
  theme_publication(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  )

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
fig4_layout <- (p_traj | p_dyn) + plot_annotation(tag_levels = 'A')
fig4_layout <- fig4_layout & theme(text = element_text(family = "sans", color = "black"))

ggsave(
  file.path(OUT_DIR, "Figure4_Trajectory.png"),
  plot = fig4_layout,
  width = FIG_WIDTH_2COL, height = 6, dpi = PUB_DPI
)
ggsave(
  file.path(OUT_DIR, "Figure4_Trajectory.pdf"),
  plot = fig4_layout,
  width = FIG_WIDTH_2COL, height = 6, dpi = PUB_DPI
)
ggsave(
  file.path(OUT_DIR, "Figure4_Trajectory.tiff"),
  plot = fig4_layout,
  width = FIG_WIDTH_2COL, height = 6, dpi = PUB_DPI, device = "tiff", compression = "lzw"
)
ggsave(
  file.path(OUT_DIR, "Figure4_Trajectory.svg"),
  plot = fig4_layout,
  width = FIG_WIDTH_2COL, height = 6, dpi = PUB_DPI, device = "svg"
)
cat("Figure 4 (Trajectory) generated and saved successfully.\n")
