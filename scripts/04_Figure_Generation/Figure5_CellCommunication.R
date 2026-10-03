# Figure 5: Intercellular Communication (CellChat)
source("utils.R")
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR <- "."
RES_DIR  <- file.path(BASE_DIR, "results")
OUT_DIR  <- file.path(BASE_DIR, "publication_figures_scripts", "output")
CACHE_DIR <- file.path(RES_DIR, "phaseD", "cache", "cellchat")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)

CELLCHAT_CACHE <- file.path(CACHE_DIR, "cellchat_object.rds")
SIG_LR_CACHE   <- file.path(CACHE_DIR, "significant_LR_pairs.rds")
NUM_GENES_FOR_CC <- 2000   # top N genes to use for CellChat

# ---------------------------------------------------------------------------
# Load required data
# ---------------------------------------------------------------------------
cat("[Fig5] Loading data objects...\n")

log_mat   <- readRDS(file.path(RES_DIR, "phaseC/objects/Phase_C_log_matrix.rds"))
clust_obj <- readRDS(file.path(RES_DIR, "phaseD/objects/Phase_D_clustering.rds"))
qc_cells  <- readRDS(file.path(RES_DIR, "phaseB/objects/Phase_B_cell_qc.rds"))

# Extract cluster labels — use clusters as cell types
if ("clusters" %in% names(clust_obj)) {
  cell_clusters <- clust_obj$clusters
} else if ("seurat_clusters" %in% names(clust_obj)) {
  cell_clusters <- clust_obj$seurat_clusters
} else {
  stop("No cluster labels found in Phase_D_clustering.rds")
}

# Extract condition labels
if ("infection_state" %in% colnames(qc_cells)) {
  cell_condition <- qc_cells[names(cell_clusters), "infection_state"]
} else if ("condition" %in% colnames(qc_cells)) {
  cell_condition <- qc_cells[names(cell_clusters), "condition"]
} else {
  stop("No infection_state or condition column in Phase_B_cell_qc.rds")
}

# Build metadata data.frame
meta_df <- data.frame(
  cell_id    = names(cell_clusters),
  cell_type  = factor(paste0("Cluster_", cell_clusters)),
  condition  = factor(cell_condition),
  stringsAsFactors = FALSE
)
rownames(meta_df) <- meta_df$cell_id

# Intersect log_mat cells with metadata
common_cells <- intersect(rownames(log_mat), meta_df$cell_id)
if (length(common_cells) < 100) {
  stop("Too few overlapping cells between log matrix and metadata (",
       length(common_cells), ")")
}
log_mat   <- log_mat[common_cells, ]
meta_df   <- meta_df[common_cells, ]

cat(sprintf("[Fig5] %d cells, %d genes loaded. %d cell types.\n",
            nrow(log_mat), ncol(log_mat), length(unique(meta_df$cell_type))))

# ---------------------------------------------------------------------------
# Subset to top-expressed genes (CellChat does not need all genes)
# ---------------------------------------------------------------------------
gene_means <- rowMeans(log_mat, na.rm = TRUE)
top_genes  <- names(sort(gene_means, decreasing = TRUE))[seq_len(min(NUM_GENES_FOR_CC, ncol(log_mat)))]
expr_input <- as.matrix(log_mat[, top_genes, drop = FALSE])

# ---------------------------------------------------------------------------
# Build CellChat object (with caching)
# ---------------------------------------------------------------------------
cat("[Fig5] Building CellChat object...\n")

if (requireNamespace("CellChat", quietly = TRUE)) {
  library(CellChat)
} else {
  stop("CellChat not installed. Run: BiocManager::install('CellChat')")
}

if (file.exists(CELLCHAT_CACHE)) {
  cat("[Fig5] Loading cached CellChat object from", CELLCHAT_CACHE, "\n")
  cellchat_obj <- readRDS(CELLCHAT_CACHE)
} else {
  cellchat_obj <- createCellChat(
    object = as.data.frame(t(expr_input)),
    meta = meta_df,
    group.by = "cell_type"
  )

  # Use human CellChat database
  cellchat_obj <- setDatabase(cellchat_obj, CellChatDB.human)

  # Subset to expressed genes
  cellchat_obj <- subsetData(cellchat_obj, features = rownames(expr_input))

  # Compute communication probability using "trppmi" (triplet permutation)
  cat("[Fig5] Computing communication probabilities...\n")
  cellchat_obj <- computeCommunProb(cellchat_obj,
                                    raw.use = TRUE,
                                    population.size = TRUE,
                                    type = "trppmi",
                                    min.nCells = 5)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 5)

  # Compute pathway activity
  cat("[Fig5] Computing pathway-level signaling...\n")
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)

  # Compute aggregated network
  cat("[Fig5] Computing aggregated network...\n")
  cellchat_obj <- aggregateNet(cellchat_obj)

  # Cache
  saveRDS(cellchat_obj, CELLCHAT_CACHE)
  cat("[Fig5] CellChat object cached to", CELLCHAT_CACHE, "\n")
}

# ---------------------------------------------------------------------------
# Extract significant ligand-receptor pairs
# ---------------------------------------------------------------------------
cat("[Fig5] Extracting significant LR interactions...\n")

if (file.exists(SIG_LR_CACHE)) {
  cat("[Fig5] Loading cached LR table from", SIG_LR_CACHE, "\n")
  sig_lr <- readRDS(SIG_LR_CACHE)
} else {
  sig_lr <- tryCatch({
    extractLR <- extractEnrichedLR(cellchat_obj, signific = TRUE)
    # Get the mean expression of LR pairs per cell type
    lr_mat <- extractEnrichedLR(cellchat_obj, signific = TRUE, slot = "net") %>%
      select(Ligand, Receptor, Source, Target, pval, prob) %>%
      distinct()

    # Add pathway information where available
    lr_mat
  }, error = function(e) {
    warning("Could not extract enriched LR: ", e$message)
    NULL
  })

  if (is.null(sig_lr) || nrow(sig_lr) == 0) {
    stop("[Fig5] No significant LR pairs found. Check min.cells threshold.")
  }
  saveRDS(sig_lr, SIG_LR_CACHE)
  cat("[Fig5] LR table cached to", SIG_LR_CACHE, "\n")
}

# ---------------------------------------------------------------------------
# 5A. Ligand-Receptor Dot Plot
# ---------------------------------------------------------------------------
# Build a cleaner LR pair label
sig_lr$LR_Pair <- paste0(sig_lr$Ligand, " - ", sig_lr$Receptor)
sig_lr$Comm    <- paste0(sig_lr$Source, " -> ", sig_lr$Target)

# Keep top N interactions by communication probability
sig_lr_top <- sig_lr %>%
  arrange(desc(prob)) %>%
  slice_head(n = 50) %>%
  mutate(
    LR_Pair = factor(LR_Pair, levels = rev(unique(LR_Pair))),
    Comm    = factor(Comm, levels = unique(Comm))
  )

p_lr <- ggplot(sig_lr_top, aes(x = Comm, y = LR_Pair)) +
  geom_point(aes(size = prob, color = -log10(pval + 1e-50))) +
  scale_color_gradientn(
    colors = c("#2166AC", "#F7F7F7", "#B2182B"),
    name = expression(-log[10]~"(p-value)")
  ) +
  scale_size_continuous(range = c(2, 8), name = "Communication\nProbability") +
  labs(
    x = "Source -> Target Cell Type",
    y = "Ligand - Receptor Pair",
    title = "Significant Intercellular Communication (Infected vs Mock)"
  ) +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(face = "bold.italic", size = 8),
    legend.position = "right"
  )

# ---------------------------------------------------------------------------
# 5B. Circle plot of incoming/outgoing signaling
# ---------------------------------------------------------------------------
cat("[Fig5] Building signaling circle plots...\n")

outgoing <- tryCatch({
  netVisual_aggregate(cellchat_obj, signaling = names(cellchat_obj@LR)$LRsig[1:5],
                      layout = "circle")
}, error = function(e) {
  warning("Could not build circle plot: ", e$message)
  NULL
})

if (!is.null(outgoing)) {
  # The circle plot is rendered directly by CellChat; save via pdf()
  sig_circle_pdf <- file.path(OUT_DIR, "Figure5_SignalingCircle.pdf")
  pdf(sig_circle_pdf, width = 8, height = 8)
  netVisual_aggregate(cellchat_obj,
                      signaling = names(cellchat_obj@LR)$LRsig[1:5],
                      layout = "circle",
                      edge.width.max = 10)
  dev.off()
  cat("[Fig5] Signaling circle saved to", sig_circle_pdf, "\n")
}

# ---------------------------------------------------------------------------
# 5C. Heatmap of pathway-level communication strength
# ---------------------------------------------------------------------------
cat("[Fig5] Building pathway communication heatmap...\n")

pathway_comm <- tryCatch({
  cc <- aggregateNet(cellchat_obj)
  netAnalysis_contribution(cellchat_obj, slot = "netP")
  cellchat_obj@netP$weight
}, error = function(e) {
  warning("Could not extract pathway communication: ", e$message)
  NULL
})

if (!is.null(pathway_comm) && nrow(pathway_comm) > 1) {
  p_path_hmap <- ComplexHeatmap::Heatmap(
    pathway_comm,
    name = "Strength",
    col = circlize::colorRamp2(c(0, max(pathway_comm, na.rm = TRUE) / 2,
                                  max(pathway_comm, na.rm = TRUE)),
                               c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "Incoming Signaling",
    row_title = "Outgoing Signaling",
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )

  pdf(file.path(OUT_DIR, "Figure5_PathwayHeatmap.pdf"),
      width = FIG_WIDTH_2COL, height = FIG_WIDTH_2COL)
  ComplexHeatmap::draw(p_path_hmap, heatmap_legend_param = list(title = "Strength"))
  dev.off()
  cat("[Fig5] Pathway heatmap saved.\n")
}

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
fig5_layout <- p_lr + plot_annotation(tag_levels = 'A')
fig5_layout <- fig5_layout & theme(text = element_text(family = "sans", color = "black"))

ggsave(
  file.path(OUT_DIR, "Figure5_Communication.png"),
  plot = fig5_layout,
  width = FIG_WIDTH_2COL, height = 7, dpi = PUB_DPI
)
ggsave(
  file.path(OUT_DIR, "Figure5_Communication.pdf"),
  plot = fig5_layout,
  width = FIG_WIDTH_2COL, height = 7, dpi = PUB_DPI
)
ggsave(
  file.path(OUT_DIR, "Figure5_Communication.tiff"),
  plot = fig5_layout,
  width = FIG_WIDTH_2COL, height = 7, dpi = PUB_DPI, device = "tiff", compression = "lzw"
)
ggsave(
  file.path(OUT_DIR, "Figure5_Communication.svg"),
  plot = fig5_layout,
  width = FIG_WIDTH_2COL, height = 7, dpi = PUB_DPI, device = "svg"
)
cat("Figure 5 (Cell Communication) generated and saved successfully.\n")
