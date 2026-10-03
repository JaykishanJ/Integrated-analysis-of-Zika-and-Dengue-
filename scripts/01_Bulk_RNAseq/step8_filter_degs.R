#!/usr/bin/env Rscript
# STEP 8 -- Filter DEGs and split up / down. GSE80434.
#
#   Rscript scripts/GSE80434/step8_filter_degs.R
#
#   padj < 0.05  AND  |log2FC| > 1
#   up   : log2FC >  1
#   down : log2FC < -1
#
# Applied to all three contrasts: DENV, ZIKVC and ZIKVM, each against its own mock set, all in
# hNPC. Thresholding uses the apeglm-SHRUNKEN log2FC -- the MLE inflates fold changes for
# low-count genes, so an |LFC| > 1 cut on MLE values admits counting artefacts.
#
# This is the dataset where a DENV-vs-ZIKV comparison is genuinely tissue-matched, so the
# cross-virus overlap is computed here -- restricted to genes both contrasts actually tested,
# because step 7 showed independent filtering retained different numbers per contrast
# (21,573 / 17,656 / 21,903).
#
# Two things to keep in view when reading the counts:
#  * DENV is ~3x more precisely measured than either ZIKV arm (median lfcSE 0.0895 vs 0.25-0.28),
#    for technical reasons established in step 5. Its larger DEG list is partly that.
#  * ZIKVC has n=2 against a control group whose internal spread exceeds the effect size, so its
#    list is the least reliable of the three.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(ggplot2); library(ggrepel)
})

GSE      <- "GSE80434"
PADJ_CUT <- 0.05
LFC_CUT  <- 1
DPI      <- 600
res_dir  <- file.path("results", GSE, "step7_results")
dds_dir  <- file.path("results", GSE, "step6_deseq")
out_dir  <- file.path("results", GSE, "step8_degs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-62s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 8: FILTER DEGs -- ", GSE, " ═══════════\n\n")
cat("── Cutoffs ──\n")
cat("  padj     < ", PADJ_CUT, "\n", sep = "")
cat("  |log2FC| > ", LFC_CUT, "   (on apeglm-shrunken LFC)\n\n", sep = "")

CONTRASTS <- list(
  list(arm = "DENV",  tag = "DENV_vs_Mock",  dds = "dds_DENV_fitted.rds"),
  list(arm = "ZIKVC", tag = "ZIKVC_vs_Mock", dds = "dds_ZIKV_fitted.rds"),
  list(arm = "ZIKVM", tag = "ZIKVM_vs_Mock", dds = "dds_ZIKV_fitted.rds"))

cascades <- list(); summ <- list(); mle_cmp <- list(); thr_cmp <- list(); sets <- list()

for (C in CONTRASTS) {
  f <- file.path(res_dir, paste0(C$tag, "_results_all_genes.csv"))
  if (!file.exists(f)) stop("missing input: ", f, " -- run step 7 first")
  res <- fread(f)

  cat("════════ ", C$arm, " vs Mock ════════\n", sep = "")
  n_tested <- sum(!is.na(res$padj))
  sig_padj <- res[!is.na(padj) & padj < PADJ_CUT]
  deg  <- sig_padj[abs(log2FC_shrunk) > LFC_CUT]
  up   <- deg[log2FC_shrunk >  LFC_CUT]
  down <- deg[log2FC_shrunk < -LFC_CUT]

  casc <- data.table(arm = C$arm,
    stage = c("genes in results table", "tested (padj not NA)",
              paste0("padj < ", PADJ_CUT),
              paste0("padj < ", PADJ_CUT, " AND |log2FC| > ", LFC_CUT),
              paste0("  up   (log2FC >  ", LFC_CUT, ")"),
              paste0("  down (log2FC < -", LFC_CUT, ")")),
    genes = c(nrow(res), n_tested, nrow(sig_padj), nrow(deg), nrow(up), nrow(down)))
  casc[, pct_of_tested := round(100 * genes / n_tested, 2)]
  print(casc[, -1], row.names = FALSE)
  cascades[[C$tag]] <- casc

  check(paste0(C$tag, ": up + down == total"), nrow(up) + nrow(down) == nrow(deg),
        sprintf("%d + %d = %d", nrow(up), nrow(down), nrow(deg)))
  check(paste0(C$tag, ": up and down disjoint"),
        length(intersect(up$gene_id, down$gene_id)) == 0)
  check(paste0(C$tag, ": every DEG meets padj cutoff"),
        nrow(deg) == 0 || all(deg$padj < PADJ_CUT))
  check(paste0(C$tag, ": every DEG meets LFC cutoff"),
        nrow(deg) == 0 || all(abs(deg$log2FC_shrunk) > LFC_CUT))
  check(paste0(C$tag, ": DEGs subset the padj-significant set"),
        all(deg$gene_id %in% sig_padj$gene_id))
  if (nrow(sig_padj))
    cat("  the |log2FC| > ", LFC_CUT, " cut removes ", nrow(sig_padj) - nrow(deg), " of ",
        nrow(sig_padj), " padj-significant genes (",
        round(100 * (1 - nrow(deg) / nrow(sig_padj)), 1), "%)\n", sep = "")

  # ---------------------------- shrunken vs MLE thresholding
  deg_mle <- sig_padj[!is.na(log2FC_MLE) & abs(log2FC_MLE) > LFC_CUT]
  only_mle <- setdiff(deg_mle$gene_id, deg$gene_id)
  mc <- data.table(arm = C$arm, thresholded_on = c("log2FC_shrunk (used)", "log2FC_MLE"),
    DEGs = c(nrow(deg), nrow(deg_mle)),
    median_baseMean = round(c(if (nrow(deg)) median(deg$baseMean) else NA_real_,
                              if (nrow(deg_mle)) median(deg_mle$baseMean) else NA_real_), 1))
  print(mc[, -1], row.names = FALSE)
  if (length(only_mle))
    cat("  genes admitted by MLE but not shrunken LFC: ", length(only_mle),
        "   their median baseMean ",
        round(median(sig_padj[gene_id %in% only_mle]$baseMean), 1), "\n", sep = "")
  mle_cmp[[C$tag]] <- mc

  # ---------------------------- formal lfcThreshold, for reference
  dds <- readRDS(file.path(dds_dir, C$dds))
  rt <- results(dds, name = paste0("condition_", C$arm, "_vs_Mock"),
                lfcThreshold = LFC_CUT, altHypothesis = "greaterAbs", alpha = PADJ_CUT)
  n_thr <- sum(rt$padj < PADJ_CUT, na.rm = TRUE)
  ov <- length(intersect(deg$gene_id, rownames(rt)[which(rt$padj < PADJ_CUT)]))
  cat("  reference: formal lfcThreshold=", LFC_CUT, " test gives ", n_thr,
      " genes (overlap with the ", nrow(deg), " filtered: ", ov, ")\n\n", sep = "")
  thr_cmp[[C$tag]] <- data.table(arm = C$arm, post_hoc_filter = nrow(deg),
                                 formal_test = n_thr, overlap = ov)
  fwrite(data.table(gene_id = rownames(rt), baseMean = rt$baseMean,
                    log2FC = rt$log2FoldChange, pvalue = rt$pvalue, padj = rt$padj
                    )[!is.na(padj) & padj < PADJ_CUT],
         file.path(out_dir, paste0(C$tag, "_reference_formal_lfcThreshold.csv")))

  # ---------------------------- export
  if (nrow(up))   setorder(up,   padj, -log2FC_shrunk)
  if (nrow(down)) setorder(down, padj,  log2FC_shrunk)
  if (nrow(deg))  setorderv(deg, "padj", na.last = TRUE)
  fwrite(deg,  file.path(out_dir, paste0(C$tag, "_DEGs_all.csv")))
  fwrite(up,   file.path(out_dir, paste0(C$tag, "_DEGs_upregulated.csv")))
  fwrite(down, file.path(out_dir, paste0(C$tag, "_DEGs_downregulated.csv")))
  writeLines(as.character(up$symbol[!is.na(up$symbol) & nzchar(up$symbol)]),
             file.path(out_dir, paste0(C$tag, "_gene_list_up_symbols.txt")))
  writeLines(as.character(down$symbol[!is.na(down$symbol) & nzchar(down$symbol)]),
             file.path(out_dir, paste0(C$tag, "_gene_list_down_symbols.txt")))
  writeLines(as.character(up$gene_id),   file.path(out_dir, paste0(C$tag, "_gene_list_up_entrez.txt")))
  writeLines(as.character(down$gene_id), file.path(out_dir, paste0(C$tag, "_gene_list_down_entrez.txt")))

  sets[[C$arm]] <- list(up = up$gene_id, down = down$gene_id,
                        universe = res$gene_id[!is.na(res$padj)])
  summ[[C$tag]] <- data.table(arm = C$arm, genes_tested = n_tested, padj_sig = nrow(sig_padj),
    DEGs = nrow(deg), up = nrow(up), down = nrow(down),
    up_down_ratio = if (nrow(down)) round(nrow(up) / nrow(down), 2) else NA_real_,
    median_absLFC = if (nrow(deg)) round(median(abs(deg$log2FC_shrunk)), 2) else NA_real_,
    max_absLFC = if (nrow(deg)) round(max(abs(deg$log2FC_shrunk)), 2) else NA_real_)

  # ---------------------------- volcano with cutoffs
  vp <- res[!is.na(pvalue)]
  vp[, class := ifelse(is.na(padj) | padj >= PADJ_CUT, "ns",
                ifelse(log2FC_shrunk >  LFC_CUT, "up",
                ifelse(log2FC_shrunk < -LFC_CUT, "down", "sig, |LFC|<=1")))]
  vp[, label := ifelse(is.na(symbol) | symbol == "", as.character(gene_id), symbol)]
  lb <- rbind(head(vp[class == "up"][order(padj)], 12), head(vp[class == "down"][order(padj)], 12))
  p <- ggplot(vp, aes(log2FC_shrunk, -log10(pvalue), colour = class)) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, colour = "grey45") +
    geom_point(alpha = 0.55, size = 1.2) +
    scale_colour_manual(values = c(up = "#C0392B", down = "#2471A3",
                                   `sig, |LFC|<=1` = "#F4B942", ns = "grey80"),
                        breaks = c("up", "down", "sig, |LFC|<=1", "ns")) +
    labs(title = paste0(GSE, " - hNPC: ", C$arm, " vs Mock"),
         subtitle = paste0("padj < ", PADJ_CUT, " and |log2FC| > ", LFC_CUT,
                           "   |   up ", nrow(up), ", down ", nrow(down),
                           "   (amber = significant but below the LFC cut)"),
         x = "shrunken log2 fold change", y = "-log10(p)", colour = NULL) +
    theme_bw(base_size = 12)
  if (nrow(lb)) p <- p + geom_text_repel(data = lb, aes(label = label), size = 2.8,
                                         colour = "black", max.overlaps = 30, show.legend = FALSE)
  ggsave(file.path(out_dir, paste0(C$tag, "_volcano_with_cutoffs.png")), p,
         width = 9, height = 7, dpi = DPI)
}

cat("════════════════════ SUMMARY ════════════════════\n")
S <- rbindlist(summ); print(S, row.names = FALSE)
fwrite(S, file.path(out_dir, "DEG_summary.csv"))
fwrite(rbindlist(cascades), file.path(out_dir, "filtering_cascade.csv"))
fwrite(rbindlist(mle_cmp),  file.path(out_dir, "shrunk_vs_mle_thresholding.csv"))
fwrite(rbindlist(thr_cmp),  file.path(out_dir, "formal_lfcThreshold_comparison.csv"))

ggsave(file.path(out_dir, "DEG_counts_all_contrasts.png"),
  ggplot(melt(S[, .(arm = factor(arm, levels = c("DENV", "ZIKVC", "ZIKVM")), up, down)],
              id.vars = "arm", variable.name = "direction", value.name = "n"),
         aes(arm, n, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.35, size = 3.4) +
    scale_fill_manual(values = c(up = "#C0392B", down = "#2471A3")) +
    labs(title = paste0(GSE, " - hNPC DEGs per contrast (padj < ", PADJ_CUT,
                        ", |log2FC| > ", LFC_CUT, ")"), x = NULL, y = "genes") +
    theme_bw(base_size = 12),
  width = 8, height = 5.5, dpi = DPI)

# ═══════════════ the cross-virus comparison this series is for
cat("\n── DENV vs ZIKV in the SAME tissue (shared tested genes only) ──\n")
ov <- list()
pairs <- list(c("DENV", "ZIKVC"), c("DENV", "ZIKVM"), c("ZIKVC", "ZIKVM"))
for (pr in pairs) {
  A <- sets[[pr[1]]]; B <- sets[[pr[2]]]
  shared_u <- intersect(A$universe, B$universe)
  for (dir in c("up", "down")) {
    a <- intersect(A[[dir]], shared_u); b <- intersect(B[[dir]], shared_u)
    ov[[length(ov) + 1]] <- data.table(virus_a = pr[1], virus_b = pr[2], direction = dir,
      shared_universe = length(shared_u), n_a = length(a), n_b = length(b),
      n_shared = length(intersect(a, b)),
      a_specific = length(setdiff(a, b)), b_specific = length(setdiff(b, a)),
      jaccard = if (length(union(a, b))) round(length(intersect(a, b)) / length(union(a, b)), 3) else NA_real_)
  }
}
OV <- rbindlist(ov); print(OV, row.names = FALSE)
fwrite(OV, file.path(out_dir, "cross_virus_overlap_shared_universe.csv"))

# Write the shared and virus-specific lists for the DENV-vs-ZIKVM pair, the best-powered one.
ann_map <- fread(file.path(res_dir, "DENV_vs_Mock_results_all_genes.csv"))[, .(gene_id, symbol)]
for (pr in list(c("DENV", "ZIKVM"), c("ZIKVC", "ZIKVM"))) {
  A <- sets[[pr[1]]]; B <- sets[[pr[2]]]; su <- intersect(A$universe, B$universe)
  for (dir in c("up", "down")) {
    a <- intersect(A[[dir]], su); b <- intersect(B[[dir]], su)
    fwrite(merge(data.table(gene_id = intersect(a, b)), ann_map, by = "gene_id", all.x = TRUE),
           file.path(out_dir, sprintf("shared_%s_%s_and_%s.csv", dir, pr[1], pr[2])))
    fwrite(merge(data.table(gene_id = setdiff(a, b)), ann_map, by = "gene_id", all.x = TRUE),
           file.path(out_dir, sprintf("%s_specific_%s_vs_%s.csv", pr[1], dir, pr[2])))
  }
}
check("cross-virus overlap computed on shared universes", nrow(OV) == 6)

cat("\n── Caveat that must travel with these counts ──\n")
cat("  DENV's list is larger partly because it is ~3x more precisely measured (median lfcSE\n")
cat("  0.0895 vs 0.25-0.28), being 3v3 on one platform against tight controls. ZIKVC is the\n")
cat("  weakest arm: n=2 against a control group more internally variable than the effect.\n")
cat("  Shared-versus-specific gene identity is the robust readout here; raw DEG counts are not.\n")

stamp <- file.path("scripts", "set_png_dpi.R")
if (file.exists(stamp)) {
  cat("\n── Figure resolution ──\n")
  cat(system2("Rscript", c(stamp, shQuote(out_dir), DPI), stdout = TRUE, stderr = TRUE), sep = "\n")
}

cat("\n═══════════ RESULT ═══════════\n")
if (length(fails)) {
  cat("  ", length(fails), " CHECK(S) FAILED:\n", sep = "")
  for (f in fails) cat("    - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("  All checks passed.\n")
for (i in seq_len(nrow(S)))
  cat("  ", S$arm[i], ": ", S$DEGs[i], " DEGs (", S$up[i], " up, ", S$down[i], " down)\n", sep = "")
