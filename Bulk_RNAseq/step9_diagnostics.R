#!/usr/bin/env Rscript
# STEP 9 -- Diagnostic visualisations. GSE80434.
#
#   Rscript scripts/GSE80434/step9_diagnostics.R
#
#   Volcano : -log10(padj) vs log2FC
#   MA      : log2FC vs baseMean
#   Heatmap : z-score scaled VST counts of top DEGs across samples
#
# For all three contrasts. Choices that decide whether these plots are honest:
#
#  * -log10(padj) is infinite wherever padj underflows to 0. Those are the strongest hits, so
#    they are capped at a marked ceiling rather than dropped, and the cap is stated on the plot.
#
#  * ZIKV heatmaps show ALL THREE arms (Mock, ZIKVC, ZIKVM), not just the pair being contrasted.
#    Given how nested the two strains turned out to be in step 8 (75% of ZIKVC's up-genes are
#    also ZIKVM's), seeing both on the same rows is the point.
#
#  * ZIKV heatmaps also annotate PLATFORM. It is a real batch effect in that sub-experiment
#    (step 5: median R2 0.449 on condition-residuals; NextSeq500 confined to Mock), so a reader
#    needs to be able to see whether a row's pattern tracks condition or instrument.
#
#  * Ranking by padj alone tends to return one direction, so a direction-balanced heatmap is
#    produced alongside.
#
#  * The heatmap uses the blind VST from step 5. Gene selection already depends on the contrast;
#    a design-aware transform on top would compound that circularity.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(ggplot2); library(ggrepel)
  library(pheatmap); library(RColorBrewer)
})

GSE      <- "GSE80434"
PADJ_CUT <- 0.05
LFC_CUT  <- 1
TOP_N    <- 50
DPI      <- 600
res_dir  <- file.path("results", GSE, "step7_results")
qc_dir   <- file.path("results", GSE, "step5_qc")
cd_in    <- file.path("results", GSE, "step2_metadata", "colData.rds")
out_dir  <- file.path("results", GSE, "step9_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("═══════════ STEP 9: DIAGNOSTIC VISUALISATIONS -- ", GSE, " ═══════════\n")
cd_all <- readRDS(cd_in)

PAL <- c(up = "#C0392B", down = "#2471A3", `sig, |LFC|<=1` = "#F4B942", ns = "grey80")
COND_COL <- c(Mock = "#5DADE2", DENV = "#8E44AD", ZIKVC = "#E74C3C", ZIKVM = "#E67E22")
PLAT_COL <- c(MiSeq = "#7F8C8D", HiSeq2000 = "#2C3E50", NextSeq500 = "#BDC3C7")

CONTRASTS <- list(
  list(arm = "DENV",  tag = "DENV_vs_Mock",  vst = "DENV_vst_matrix.csv"),
  list(arm = "ZIKVC", tag = "ZIKVC_vs_Mock", vst = "ZIKV_vst_matrix.csv"),
  list(arm = "ZIKVM", tag = "ZIKVM_vs_Mock", vst = "ZIKV_vst_matrix.csv"))

info <- list()

for (C in CONTRASTS) {
  f <- file.path(res_dir, paste0(C$tag, "_results_all_genes.csv"))
  if (!file.exists(f)) stop("missing input: ", f, " -- run step 7 first")
  res <- fread(f)
  cat("\n════════ ", C$arm, " vs Mock ════════\n", sep = "")

  res[, class := ifelse(is.na(padj) | padj >= PADJ_CUT, "ns",
               ifelse(log2FC_shrunk >  LFC_CUT, "up",
               ifelse(log2FC_shrunk < -LFC_CUT, "down", "sig, |LFC|<=1")))]
  res[, label := ifelse(is.na(symbol) | symbol == "", as.character(gene_id), symbol)]
  n_up <- sum(res$class == "up"); n_dn <- sum(res$class == "down")
  cat("  up ", n_up, ", down ", n_dn, " at padj < ", PADJ_CUT,
      " & |log2FC| > ", LFC_CUT, "\n", sep = "")

  # ══════════════════════════════════════════ 1. VOLCANO on padj
  v <- res[!is.na(padj)]
  n_zero <- sum(v$padj == 0)
  fin_max <- if (any(v$padj > 0)) max(-log10(v$padj[v$padj > 0])) else 1
  cap <- ceiling(fin_max / 25) * 25 + 25
  v[, negLog10padj := ifelse(padj == 0, cap, -log10(padj))]
  v[, capped := padj == 0]
  cat("  volcano: padj = 0 for ", n_zero, " gene(s) | largest finite -log10(padj) ",
      round(fin_max, 1), " | cap ", cap, "\n", sep = "")

  lb <- rbind(head(v[class == "up"][order(padj)], 12), head(v[class == "down"][order(padj)], 12))
  p <- ggplot(v, aes(log2FC_shrunk, negLog10padj, colour = class)) +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = 2, colour = "grey40") +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, colour = "grey40") +
    geom_point(aes(shape = capped), alpha = 0.6, size = 1.3) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_colour_manual(values = PAL, breaks = c("up", "down", "sig, |LFC|<=1", "ns")) +
    labs(title = paste0(GSE, " - hNPC: ", C$arm, " vs Mock"),
         subtitle = paste0("padj < ", PADJ_CUT, " and |log2FC| > ", LFC_CUT,
                           "  |  up ", n_up, ", down ", n_dn, "  |  apeglm-shrunken LFC"),
         x = expression(shrunken~log[2]~fold~change),
         y = expression(-log[10](adjusted~italic(p))), colour = NULL) +
    theme_bw(base_size = 12)
  if (nrow(lb)) p <- p + geom_text_repel(data = lb, aes(label = label), size = 2.8,
                                         colour = "black", max.overlaps = 40, show.legend = FALSE)
  if (n_zero > 0) p <- p + annotate("text", x = min(v$log2FC_shrunk), y = cap, hjust = 0,
        vjust = 1.5, size = 3, colour = "grey30",
        label = paste0("triangles: padj = 0, capped at ", cap))
  ggsave(file.path(out_dir, paste0(C$tag, "_volcano_padj.png")), p, width = 9, height = 7.5, dpi = DPI)

  # ══════════════════════════════════════════ 2. MA
  m <- res[baseMean > 0]
  lm_ <- rbind(head(m[class == "up"][order(-abs(log2FC_shrunk))], 10),
               head(m[class == "down"][order(-abs(log2FC_shrunk))], 10))
  p2 <- ggplot(m, aes(baseMean, log2FC_shrunk, colour = class)) +
    geom_hline(yintercept = 0, colour = "grey30") +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, colour = "grey55") +
    geom_point(alpha = 0.5, size = 1.05) +
    scale_x_log10(labels = scales::comma) +
    scale_colour_manual(values = PAL, breaks = c("up", "down", "sig, |LFC|<=1", "ns")) +
    labs(title = paste0(GSE, " - MA: hNPC, ", C$arm, " vs Mock"),
         subtitle = "shrunken log2FC; low-expression genes are pulled toward 0",
         x = "mean of normalised counts (baseMean, log10)",
         y = expression(shrunken~log[2]~fold~change), colour = NULL) +
    theme_bw(base_size = 12)
  if (nrow(lm_)) p2 <- p2 + geom_text_repel(data = lm_, aes(label = label), size = 2.7,
                                            colour = "black", max.overlaps = 30, show.legend = FALSE)
  ggsave(file.path(out_dir, paste0(C$tag, "_MA_plot.png")), p2, width = 9, height = 6.5, dpi = DPI)
  cat("  MA: baseMean ", round(min(m$baseMean), 2), " - ",
      format(round(max(m$baseMean)), big.mark = ","), "\n", sep = "")

  # ══════════════════════════════════════════ 3. HEATMAP
  vst <- fread(file.path(qc_dir, C$vst))
  vm <- as.matrix(vst[, -1, with = FALSE]); rownames(vm) <- as.character(vst[[1]])
  cd <- cd_all[colnames(vm), , drop = FALSE]
  stopifnot(identical(colnames(vm), rownames(cd)))
  ord <- order(factor(as.character(cd$condition), levels = c("Mock", "DENV", "ZIKVC", "ZIKVM")),
               cd$platform)

  deg <- res[class %in% c("up", "down")]
  fallback <- FALSE
  if (nrow(deg) == 0) {
    deg <- res[!is.na(padj) & padj < PADJ_CUT]; fallback <- TRUE
    cat("  heatmap: 0 genes clear the LFC cut -> falling back to ", nrow(deg),
        " padj-significant genes\n", sep = "")
  }
  if (nrow(deg) == 0) {
    cat("  heatmap: nothing to plot -- skipped\n")
    info[[C$tag]] <- data.table(arm = C$arm, n_up = n_up, n_down = n_dn,
                                padj_zero = n_zero, heatmap = "skipped")
    next
  }
  deg[, .k := -abs(log2FC_shrunk)]
  setorderv(deg, c("padj", ".k"), na.last = TRUE); deg[, .k := NULL]

  draw_hm <- function(sel, file, title) {
    sel <- sel[gene_id %in% rownames(vm)]
    if (!nrow(sel)) return(invisible(NULL))
    mat <- vm[as.character(sel$gene_id), , drop = FALSE]
    keepv <- apply(mat, 1, sd) > 0          # zero-variance rows would z-score to NaN
    mat <- mat[keepv, , drop = FALSE]; sel <- sel[keepv]
    if (!nrow(mat)) return(invisible(NULL))
    rownames(mat) <- make.unique(ifelse(is.na(sel$symbol) | sel$symbol == "",
                                        as.character(sel$gene_id), sel$symbol))
    z <- t(scale(t(mat)))
    ann_c <- data.frame(condition = as.character(cd$condition),
                        platform = as.character(cd$platform), row.names = rownames(cd))
    ann_r <- data.frame(direction = ifelse(sel$log2FC_shrunk > 0, "up", "down"),
                        row.names = rownames(z))
    png(file.path(out_dir, file), width = 9.5, height = max(5, 0.2 * nrow(z)),
        units = "in", res = DPI)
    pheatmap(z[, ord, drop = FALSE], annotation_col = ann_c, annotation_row = ann_r,
             annotation_colors = list(condition = COND_COL, platform = PLAT_COL,
                                      direction = c(up = "#C0392B", down = "#2471A3")),
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
             breaks = seq(-max(abs(z)), max(abs(z)), length.out = 256),
             cluster_rows = TRUE, cluster_cols = FALSE,
             show_rownames = nrow(z) <= 80, fontsize_row = 7, fontsize_col = 8,
             border_color = NA, main = title)
    dev.off()
    fwrite(data.table(gene = rownames(z), as.data.table(round(z, 4))),
           file.path(out_dir, sub("[.]png$", "_zscores.csv", file)))
    nrow(z)
  }

  arms_shown <- paste(levels(droplevels(factor(as.character(cd$condition)))), collapse = ", ")
  n1 <- draw_hm(head(deg, TOP_N), paste0(C$tag, "_heatmap_top50_zscore.png"),
                paste0(GSE, " - hNPC: ", C$arm, " vs Mock\n",
                       if (fallback) "top padj-significant genes (none clear |LFC|>1)"
                       else paste0("top ", min(TOP_N, nrow(deg)), " DEGs"),
                       ", z-scored VST | arms shown: ", arms_shown))
  half <- TOP_N %/% 2
  bal <- rbind(head(deg[log2FC_shrunk > 0][order(padj)], half),
               head(deg[log2FC_shrunk < 0][order(padj)], half))
  n2 <- draw_hm(bal, paste0(C$tag, "_heatmap_balanced_zscore.png"),
                paste0(GSE, " - hNPC: ", C$arm, " vs Mock\ntop ", half, " up + top ", half,
                       " down, z-scored VST"))
  n_up_top <- sum(head(deg, TOP_N)$log2FC_shrunk > 0)
  cat("  heatmap: top ", min(TOP_N, nrow(deg)), " contains ", n_up_top, " up / ",
      min(TOP_N, nrow(deg)) - n_up_top, " down; balanced version ", n2, " rows",
      " | arms shown: ", arms_shown, "\n", sep = "")

  fwrite(head(deg, TOP_N)[, .(gene_id, symbol, baseMean = round(baseMean, 1),
         log2FC_shrunk = round(log2FC_shrunk, 3), padj, description)],
         file.path(out_dir, paste0(C$tag, "_heatmap_top50_gene_table.csv")))

  info[[C$tag]] <- data.table(arm = C$arm, n_up = n_up, n_down = n_dn, padj_zero = n_zero,
                              heatmap = paste0(n1, " rows",
                                               if (fallback) " (padj fallback)" else ""),
                              arms_in_heatmap = arms_shown)
}

cat("\n════════════════════ SUMMARY ════════════════════\n")
I <- rbindlist(info, fill = TRUE); print(I, row.names = FALSE)
fwrite(I, file.path(out_dir, "diagnostics_summary.csv"))

cat("\n── Reading the ZIKV heatmaps ──\n")
cat("  They carry a platform annotation row because platform is a genuine batch effect in that\n")
cat("  sub-experiment and NextSeq500 sits only in Mock. If a row's colour tracks the platform\n")
cat("  band rather than the condition band, treat that gene with suspicion even though the\n")
cat("  model adjusted for platform.\n")

stamp <- file.path("scripts", "set_png_dpi.R")
if (file.exists(stamp)) {
  cat("\n── Figure resolution ──\n")
  cat(system2("Rscript", c(stamp, shQuote(out_dir), DPI), stdout = TRUE, stderr = TRUE), sep = "\n")
}

cat("\n── Output ──\n")
for (f in sort(list.files(out_dir))) cat("  ", f, "\n", sep = "")
cat("\n═══════════ RESULT ═══════════\n")
cat("  Volcano, MA and z-scored heatmaps written for all three contrasts at ", DPI, " dpi.\n", sep = "")
