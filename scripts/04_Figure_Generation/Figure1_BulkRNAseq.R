# ============================================================================
# FIGURE 1 — INTEGRATED BULK RNA-seq ANALYSIS
# ============================================================================
#
# Panels:
#   A — Study overview (dataset summary tile)
#   B — PCA across the three bulk datasets
#   C — Volcano plots for each dataset
#   D — GO Biological Process enrichment (dot plots)
#   E — Cross-dataset conserved upregulated genes
#
# Output: PNG (600 DPI), PDF (vector, cairo), SVG (vector)
# ============================================================================

source("utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)

  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)

  library(readr)
  library(scales)
  library(grid)
  library(gridExtra)
})

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

BASE_DIR <- "."
BULK_DIR <- file.path(BASE_DIR, "ZIKA_Bulk_Som", "ZIKA_Bulk_Som")
OUT_DIR  <- file.path(BASE_DIR, "publication_figures_scripts", "output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_WIDTH  <- 180   # mm, two-col
FIG_HEIGHT <- 220   # mm, tall enough for 5 panels
DPI        <- 600

# ============================================================================
# 2. COLOR PALETTE
# ============================================================================

COLORS <- list(
  ZIKV_AFRICAN = "#B2182B",
  DENV         = "#2166AC",
  ZIKV_ASIAN   = "#762A83",
  UP           = "#B2182B",
  DOWN         = "#2166AC",
  NS           = "#D9D9D9",
  CONTROL      = "#7F7F7F",
  TEXT         = "#222222",
  GRID         = "#E6E6E6"
)

# ============================================================================
# 3. DATASETS
# ============================================================================

DATASETS <- list(
  list(
    id       = "GSE80434_hNPC_ZIKVM",
    short    = "ZIKV-M",
    cell     = "hNPC",
    virus    = "African ZIKV",
    color    = COLORS$ZIKV_AFRICAN,
    dds_path = file.path(BULK_DIR, "results", "GSE80434", "step6_deseq", "dds_ZIKV_fitted.rds"),
    deg_csv  = file.path(BULK_DIR, "results", "GSE80434", "step8_degs", "ZIKVM_vs_Mock_DEGs_all.csv"),
    up_file  = file.path(BULK_DIR, "results", "GSE80434", "step8_degs", "ZIKVM_vs_Mock_gene_list_up_entrez.txt"),
    down_file= file.path(BULK_DIR, "results", "GSE80434", "step8_degs", "ZIKVM_vs_Mock_gene_list_down_entrez.txt"),
    padj_cut = 0.05, lfc_cut = 1
  ),
  list(
    id       = "GSE110512_Huh7_DENV",
    short    = "DENV-2",
    cell     = "Huh7",
    virus    = "DENV-2",
    color    = COLORS$DENV,
    dds_path = file.path(BULK_DIR, "results", "GSE110512", "step6_deseq", "dds_fitted.rds"),
    deg_csv  = file.path(BULK_DIR, "results", "GSE110512", "step8_degs", "DEGs_all.csv"),
    up_file  = file.path(BULK_DIR, "results", "GSE110512", "step8_degs", "gene_list_up_entrez.txt"),
    down_file= file.path(BULK_DIR, "results", "GSE110512", "step8_degs", "gene_list_down_entrez.txt"),
    padj_cut = 0.05, lfc_cut = 1
  ),
  list(
    id       = "GSE161783_Huh7.5_Zika",
    short    = "ZIKV-C",
    cell     = "Huh7.5 / moDC",
    virus    = "Asian ZIKV",
    color    = COLORS$ZIKV_ASIAN,
    dds_path = file.path(BULK_DIR, "results", "GSE161783", "step6_deseq", "dds_Huh7_5_fitted.rds"),
    deg_csv  = file.path(BULK_DIR, "results", "GSE161783", "step8_degs", "Huh7_5__Zika_vs_Mock_DEGs_all.csv"),
    up_file  = file.path(BULK_DIR, "results", "GSE161783", "step8_degs", "Huh7_5__Zika_vs_Mock_gene_list_up_entrez.txt"),
    down_file= file.path(BULK_DIR, "results", "GSE161783", "step8_degs", "Huh7_5__Zika_vs_Mock_gene_list_down_entrez.txt"),
    padj_cut = 0.05, lfc_cut = 1
  )
)

# ============================================================================
# 4. PUBLICATION THEME
# ============================================================================

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text          = element_text(family = "Arial", colour = COLORS$TEXT),
      axis.title    = element_text(size = 8, face = "plain"),
      axis.text     = element_text(size = 7, colour = COLORS$TEXT),
      axis.line     = element_line(linewidth = 0.35, colour = COLORS$TEXT),
      axis.ticks    = element_line(linewidth = 0.3,  colour = COLORS$TEXT),
      legend.title  = element_text(size = 7, face = "bold"),
      legend.text   = element_text(size = 6.5),
      legend.key.size = unit(0.35, "cm"),
      plot.title    = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7, colour = "#555555"),
      plot.margin   = margin(4, 5, 4, 5),
      strip.background = element_blank(),
      strip.text    = element_text(size = 7, face = "bold")
    )
}

# ============================================================================
# 5. HELPERS
# ============================================================================

panel_label <- function(x) {
  annotate("text", x = -Inf, y = Inf, label = x,
           hjust = -0.15, vjust = 1.35, fontface = "bold", size = 4)
}

find_group_col <- function(cd) {
  possible <- c("condition", "treatment", "group", "infection",
                "infected", "virus", "status", "sample_type", "arm")
  hit <- possible[tolower(possible) %in% tolower(colnames(cd))]
  if (length(hit) > 0) {
    return(colnames(cd)[match(tolower(hit[1]), tolower(colnames(cd)))])
  }
  candidates <- colnames(cd)[vapply(cd, function(x) {
    ux <- unique(x[!is.na(x)])
    length(ux) >= 2 && length(ux) <= 10
  }, logical(1))]
  if (length(candidates) == 0) stop("No suitable grouping variable found.")
  candidates[1]
}

# ---------------------------------------------------------------------------
# DEG loader — returns tibble with gene, log2FC, padj, ENTREZID
# ---------------------------------------------------------------------------
load_degs <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  lfc_col  <- intersect(c("log2FC_shrunk", "log2FoldChange", "log2FC_MLE"), colnames(df))[1]
  padj_col <- intersect(c("padj", "p_adj", "padj_BH"), colnames(df))[1]
  sym_col  <- intersect(c("symbol", "gene_symbol", "Gene"), colnames(df))[1]
  id_col   <- intersect(c("gene_id", "entrez", "ENTREZID"), colnames(df))[1]

  if (is.na(lfc_col) || is.na(padj_col))
    stop("Could not identify LFC/padj columns in: ", path)

  out <- tibble(
    gene      = if (!is.na(sym_col)) as.character(df[[sym_col]]) else rownames(df),
    log2FC    = suppressWarnings(as.numeric(df[[lfc_col]])),
    padj      = suppressWarnings(as.numeric(df[[padj_col]])),
    ENTREZID  = if (!is.na(id_col)) as.character(df[[id_col]]) else NA_character_
  )
  out <- out |> filter(is.finite(log2FC), !is.na(padj))
  out
}

# ---------------------------------------------------------------------------
# Safe enrichGO with progressive cutoff relaxation
# ---------------------------------------------------------------------------
run_go_safe <- function(genes, universe) {
  genes <- unique(as.character(genes)); genes <- genes[is.finite(as.numeric(genes))]
  universe <- unique(as.character(universe)); universe <- universe[is.finite(as.numeric(universe))]
  if (length(genes) < 3 || length(universe) < 10) return(NULL)

  result <- NULL
  for (cut in list(c(0.05, 0.01), c(0.10, 0.05), c(0.20, 0.10))) {
    qc <- cut[1]; pc <- cut[2]
    res <- tryCatch({
      enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID", ont = "BP",
               pAdjustMethod = "BH", qvalueCutoff = qc, pvalueCutoff = pc,
               readable = TRUE)
    }, error = function(e) NULL, warning = function(w) NULL)

    if (!is.null(res)) {
      df <- as.data.frame(res)
      if (nrow(df) > 0) {
        cat(sprintf("  [Fig1] enrichGO: %d terms (q<%.2f, p<%.2f)\n", nrow(df), qc, pc))
        result <- df
        break
      }
    }
  }
  result
}

# ============================================================================
# 6. PANEL A — STUDY OVERVIEW
# ============================================================================

make_overview <- function() {
  overview_df <- tibble(
    Dataset    = c("GSE80434", "GSE110512", "GSE161783"),
    Cell       = c("hNPC", "Huh7", "Huh7.5 / moDC"),
    Infection  = c("African ZIKV", "DENV-2", "Asian ZIKV"),
    Comparison = c("ZIKV-M vs Mock", "DENV-2 vs Mock", "ZIKV-C vs Mock")
  )

  ggplot(overview_df, aes(x = Dataset, y = Infection, fill = Infection)) +
    geom_tile(alpha = 0.9, height = 0.8) +
    geom_text(aes(label = paste0(Cell, "\n", Comparison)),
              colour = "white", size = 2.8, lineheight = 0.9) +
    scale_fill_manual(values = c(
      "African ZIKV" = COLORS$ZIKV_AFRICAN,
      "DENV-2"       = COLORS$DENV,
      "Asian ZIKV"   = COLORS$ZIKV_ASIAN
    )) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank(),
          legend.position = "none",
          panel.border = element_blank())
}

# ============================================================================
# 7. PANEL B — PCA
# ============================================================================

make_pca <- function(ds) {
  dds <- readRDS(ds$dds_path)
  vsd <- vst(dds, blind = FALSE)
  group_col <- find_group_col(colData(dds))
  pca <- plotPCA(vsd, intgroup = group_col, returnData = TRUE)

  df <- tibble(
    PC1       = pca$PC1,
    PC2       = pca$PC2,
    Condition = as.factor(pca$name)
  )

  ggplot(df, aes(PC1, PC2, colour = Condition)) +
    geom_point(size = 2.4, alpha = 0.85) +
    scale_color_manual(values = setNames(
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nlevels(df$Condition)),
      levels(df$Condition))) +
    guides(colour = guide_legend(title = NULL)) +
    labs(x = sprintf("PC1 (%.1f%%)", pca$percentVar[1]),
         y = sprintf("PC2 (%.1f%%)", pca$percentVar[2]),
         title = ds$short) +
    theme_pub() +
    theme(legend.position = "right", legend.key.size = unit(0.35, "cm"))
}

# ============================================================================
# 8. PANEL C — VOLCANO
# ============================================================================

make_volcano <- function(ds, label_n = 8) {
  res <- load_degs(ds$deg_csv)

  res <- res |>
    mutate(
      significance = case_when(
        padj < ds$padj_cut & log2FC >= ds$lfc_cut  ~ "Up",
        padj < ds$padj_cut & log2FC <= -ds$lfc_cut ~ "Down",
        TRUE                                        ~ "NS"
      ),
      neglog10padj = -log10(pmax(padj, 1e-300))
    )

  labels <- res |>
    filter(significance != "NS") |>
    arrange(desc(abs(log2FC) * neglog10padj)) |>
    slice_head(n = label_n)

  ggplot(res, aes(log2FC, neglog10padj)) +
    geom_point(data = subset(res, significance == "NS"),
               colour = COLORS$NS, size = 0.8, alpha = 0.55) +
    geom_point(data = subset(res, significance != "NS"),
               aes(colour = significance), size = 1.0, alpha = 0.75) +
    scale_color_manual(values = c(Up = COLORS$UP, Down = COLORS$DOWN)) +
    geom_hline(yintercept = -log10(ds$padj_cut),
               linetype = "dashed", linewidth = 0.3, colour = "#555555") +
    geom_vline(xintercept = c(-ds$lfc_cut, ds$lfc_cut),
               linetype = "dashed", linewidth = 0.3, colour = "#555555") +
    geom_text_repel(data = labels, aes(label = gene),
                    size = 2.5, box.padding = 0.35, point.padding = 0.2,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.linewidth = 0.25, seed = 42) +
    scale_x_continuous(limits = function(x) {
      rng <- range(x, na.rm = TRUE)
      c(min(rng[1], -ds$lfc_cut - 0.5), max(rng[2], ds$lfc_cut + 0.5))
    }) +
    labs(x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"adjusted"~italic(P)),
         title = ds$short) +
    guides(colour = guide_legend(title = NULL)) +
    theme_pub() +
    theme(legend.position = "none")
}

# ============================================================================
# 9. PANEL D — GO ENRICHMENT DOT PLOT
# ============================================================================

make_enrichment <- function(ds) {
  up   <- suppressWarnings(as.numeric(read_lines(ds$up_file)));    up   <- up[is.finite(up)]
  down <- suppressWarnings(as.numeric(read_lines(ds$down_file)));  down <- down[is.finite(down)]

  deg <- load_degs(ds$deg_csv)
  universe <- deg$ENTREZID[!is.na(deg$ENTREZID)]

  up_go   <- run_go_safe(up,   universe)
  down_go <- run_go_safe(down, universe)

  results <- list()
  if (!is.null(up_go))   { df <- as.data.frame(up_go);   df$Direction <- "Upregulated";   results[[1]] <- df }
  if (!is.null(down_go)) { df <- as.data.frame(down_go); df$Direction <- "Downregulated"; results[[length(results)+1]] <- df }

  if (length(results) == 0) {
    return(ggplot() + theme_pub() +
             annotate("text", 0, 0, label = sprintf("No significant GO:BP\n(q < 0.20)")) +
             labs(title = sprintf("GO:BP enrichment — %s", ds$short)))
  }

  df <- bind_rows(results) |>
    group_by(Direction) |>
    arrange(p.adjust, .by_group = TRUE) |>
    slice_head(n = 8) |>
    ungroup() |>
    mutate(Description = reorder(Description, -log10(p.adjust)))

  ggplot(df, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, colour = p.adjust), alpha = 0.85) +
    scale_color_continuous(low = "#B2182B", high = "#2166AC",
                           trans = "reverse", name = "Adj. P") +
    scale_size_continuous(range = c(1.5, 5), name = "Genes") +
    facet_wrap(~Direction, scales = "free_y") +
    labs(x = "Gene ratio", y = NULL,
         title = sprintf("GO:BP enrichment — %s", ds$short)) +
    theme_pub(base_size = 7) +
    theme(strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 6),
          legend.position = "right")
}

# ============================================================================
# 10. PANEL E — SHARED GENES
# ============================================================================

make_shared_genes <- function() {
  up_lists <- lapply(DATASETS, function(ds) {
    x <- suppressWarnings(as.numeric(read_lines(ds$up_file)))
    x[is.finite(x)]
  })
  names(up_lists) <- vapply(DATASETS, `[[`, character(1), "short")

  shared_df <- tibble(Gene = unique(unlist(up_lists)))
  shared_df$Dataset_Count <- vapply(shared_df$Gene, function(g) {
    sum(vapply(up_lists, function(x) g %in% x, logical(1)))
  }, numeric(1))
  shared_df <- shared_df |> filter(Dataset_Count >= 2)

  if (nrow(shared_df) == 0) {
    return(ggplot() + theme_pub() +
             annotate("text", 0, 0, label = "No shared upregulated genes") +
             labs(title = "Cross-dataset conservation"))
  }

  shared_df$Symbol <- mapIds(org.Hs.eg.db,
                             keys = as.character(shared_df$Gene),
                             column = "SYMBOL", keytype = "ENTREZID",
                             multiVals = "first")
  shared_df$Symbol[is.na(shared_df$Symbol)] <- as.character(shared_df$Gene[is.na(shared_df$Symbol)])

  shared_summary <- shared_df |> dplyr::count(Dataset_Count, name = "n")

  ggplot(shared_summary, aes(x = factor(Dataset_Count), y = n)) +
    geom_col(width = 0.65, fill = "#4D4D4D") +
    geom_text(aes(label = n), vjust = -0.4, size = 3) +
    labs(x = "Number of datasets", y = "Shared upregulated genes",
         title = "Cross-dataset conservation") +
    theme_pub()
}

# ============================================================================
# 11. GENERATE ALL PANELS
# ============================================================================

message("\n========== Generating Figure 1 ==========\n")

pca_list      <- list()
volcano_list  <- list()
enrich_list   <- list()

for (ds in DATASETS) {
  message("Processing: ", ds$id)
  pca_list[[ds$short]]     <- make_pca(ds)
  volcano_list[[ds$short]] <- make_volcano(ds)
  enrich_list[[ds$short]]  <- make_enrichment(ds)
}

# ============================================================================
# 12. ASSEMBLE FIGURE
# ============================================================================

p_overview <- make_overview()
p_shared   <- make_shared_genes()

pca_panel <- wrap_plots(pca_list, nrow = 1) +
  plot_annotation(title = "Sample-level transcriptomic structure",
                  theme = theme(plot.title = element_text(family = "Arial", size = 9, face = "bold")))

volcano_panel <- wrap_plots(volcano_list, nrow = 1) +
  plot_annotation(title = "Differential gene expression",
                  theme = theme(plot.title = element_text(family = "Arial", size = 9, face = "bold")))

enrichment_panel <- wrap_plots(enrich_list, nrow = 1) +
  plot_annotation(title = "Functional enrichment of differentially expressed genes",
                  theme = theme(plot.title = element_text(family = "Arial", size = 9, face = "bold")))

FIGURE_1 <-
  (p_overview  + panel_label("A")) /
  (pca_panel   + panel_label("B")) /
  (volcano_panel + panel_label("C")) /
  (enrichment_panel + panel_label("D")) /
  (p_shared    + panel_label("E")) +
  plot_layout(heights = c(0.6, 1.3, 1.3, 1.5, 0.8))

# ============================================================================
# 13. EXPORT
# ============================================================================

outfile <- file.path(OUT_DIR, "Figure_1_Bulk_RNAseq_Publication")

ggsave(paste0(outfile, ".png"), FIGURE_1,
       width = FIG_WIDTH, height = FIG_HEIGHT, units = "mm",
       dpi = DPI, bg = "white")

ggsave(paste0(outfile, ".pdf"), FIGURE_1,
       width = FIG_WIDTH, height = FIG_HEIGHT, units = "mm",
       device = cairo_pdf, dpi = DPI, bg = "white")

ggsave(paste0(outfile, ".tiff"), FIGURE_1,
       width = FIG_WIDTH, height = FIG_HEIGHT, units = "mm",
       device = "tiff", dpi = DPI, compression = "lzw", bg = "white")

tryCatch({
  ggsave(paste0(outfile, ".svg"), FIGURE_1,
         width = FIG_WIDTH, height = FIG_HEIGHT, units = "mm",
         device = "svg", bg = "white")
}, error = function(e) {
  warning("SVG export failed (install svglite): ", e$message)
})

# Save shared gene table
if (exists("shared_df") && nrow(shared_df) > 0) {
  write_csv(shared_df, file.path(OUT_DIR, "Figure_1_shared_upregulated_genes.csv"))
}

message("\n=================================================\n",
        "Figure 1 successfully generated.\n",
        "  PNG: ", outfile, ".png  (", DPI, " DPI)\n",
        "  PDF: ", outfile, ".pdf  (vector)\n",
        "=================================================\n")
