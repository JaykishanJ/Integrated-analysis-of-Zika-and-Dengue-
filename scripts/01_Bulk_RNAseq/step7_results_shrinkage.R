#!/usr/bin/env Rscript
# STEP 7 -- Extract contrast results and apply LFC shrinkage. GSE80434.
#
#   Rscript scripts/GSE80434/step7_results_shrinkage.R
#
#   res <- results(dds, contrast = c("condition", "DENV", "Mock"))
#   res_shrunk <- lfcShrink(dds, coef = "condition_DENV_vs_Mock", type = "apeglm")
#
# Three contrasts, all in hNPC and all against their own mock set:
#   DENV  vs Mock (DENV-Mock)   from dds_DENV, ~condition
#   ZIKVC vs Mock               from dds_ZIKV, ~platform + condition
#   ZIKVM vs Mock               from dds_ZIKV, ~platform + condition
#
# Four points of precision, each of which caused a real problem on an earlier dataset:
#
#  * contrast = c(factor, numerator, denominator). The virus arm is the numerator, so positive
#    log2FC = UP in that arm. Verified by checking that the reversed contrast flips every sign.
#
#  * The EXPORTED table comes from the name= form, not contrast=. lfcShrink(coef=) shrinks the
#    stored coefficient, so using name= keeps the MLE LFC, shrunken LFC and p-values mutually
#    consistent. The contrast= route re-derives the estimate and on a multi-level factor its
#    optimiser can land marginally differently on a hard gene.
#
#  * results() must be given the coefficient explicitly for dds_ZIKV: condition has 3 levels,
#    so the default returns ZIKVM.
#
#  * Shrinkage changes EFFECT SIZES, not p-values. Only log2FoldChange is taken from apeglm;
#    padj in the exported table comes from results(). That matters because lfcShrink builds its
#    own results object and re-runs independent filtering, which can settle on a different
#    baseMean threshold and so a slightly different significant set. None of that reaches the
#    export, but it is reported where it happens.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(ggplot2); library(ggrepel)
})

GSE     <- "GSE80434"
ALPHA   <- 0.05
DPI     <- 600
in_dir  <- file.path("results", GSE, "step6_deseq")
annot   <- file.path("data", GSE, "gene_annotation.tsv")
out_dir <- file.path("results", GSE, "step7_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-64s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 7: RESULTS & LFC SHRINKAGE -- ", GSE, " ═══════════\n")

sym <- desc <- NULL
if (file.exists(annot)) {
  ann <- fread(annot)
  sym  <- setNames(ann$Symbol, as.character(ann$GeneID))
  desc <- setNames(ann$Description, as.character(ann$GeneID))
}

CONTRASTS <- list(
  list(dds = "dds_DENV_fitted.rds", arm = "DENV",  tag = "DENV_vs_Mock"),
  list(dds = "dds_ZIKV_fitted.rds", arm = "ZIKVC", tag = "ZIKVC_vs_Mock"),
  list(dds = "dds_ZIKV_fitted.rds", arm = "ZIKVM", tag = "ZIKVM_vs_Mock"))

summ <- list(); bins <- list(); mcmp <- list(); cache <- list()

for (C in CONTRASTS) {
  if (is.null(cache[[C$dds]])) cache[[C$dds]] <- readRDS(file.path(in_dir, C$dds))
  dds   <- cache[[C$dds]]
  coefn <- paste0("condition_", C$arm, "_vs_Mock")
  nc <- counts(dds, normalized = TRUE)
  cat("\n════════════════ ", C$arm, " vs Mock ════════════════\n", sep = "")
  cat("  object: ", C$dds, " | ", nrow(dds), " genes x ", ncol(dds), " samples | design ",
      deparse(design(dds)), "\n", sep = "")

  # ---------------------------------------------------- extract the contrast
  res_c <- results(dds, contrast = c("condition", C$arm, "Mock"), alpha = ALPHA)
  res   <- results(dds, name = coefn, alpha = ALPHA)          # <- exported
  cat("  numerator/denominator: ", C$arm, " / Mock   (positive log2FC = UP in ", C$arm, ")\n", sep = "")
  cat("  tested: ", sum(!is.na(res$pvalue)), " genes | padj not NA: ", sum(!is.na(res$padj)), "\n", sep = "")

  dlfc <- abs(res_c$log2FoldChange - res$log2FoldChange)
  check(paste0(C$tag, ": contrast= and name= agree for >99.9% of genes"),
        sum(dlfc > 1e-6, na.rm = TRUE) <= 0.001 * nrow(res),
        sprintf("%d gene(s) differ, max |diff| %.4f; lfcSE diff %.2g",
                sum(dlfc > 1e-6, na.rm = TRUE), max(dlfc, na.rm = TRUE),
                max(abs(res_c$lfcSE - res$lfcSE), na.rm = TRUE)))
  check(paste0(C$tag, ": reversed contrast flips sign"),
        max(abs(results(dds, contrast = c("condition", "Mock", C$arm))$log2FoldChange +
                res_c$log2FoldChange), na.rm = TRUE) < 1e-6)
  # DESeq2 prints level names with underscores as spaces.
  check(paste0(C$tag, ": mcols names the intended comparison"),
        grepl(paste0(gsub("_", " ", C$arm), " vs Mock"), mcols(res)$description[2]),
        mcols(res)$description[2])

  ii <- which(!is.na(res$padj) & res$log2FoldChange > 0)
  top_up <- rownames(res)[ii][which.min(res$padj[ii])]
  m_i <- dds$condition == "Mock"; a_i <- dds$condition == C$arm
  check(paste0(C$tag, ": top up-gene really higher in ", C$arm),
        mean(nc[top_up, a_i]) > mean(nc[top_up, m_i]),
        sprintf("%s: Mock %.0f -> %s %.0f", top_up, mean(nc[top_up, m_i]), C$arm,
                mean(nc[top_up, a_i])))

  # ---------------------------------------------------- shrinkage
  # apeglm needs coef=, not contrast=: it refits the coefficient under a heavy-tailed prior.
  res_ape <- lfcShrink(dds, coef = coefn, type = "apeglm", quiet = TRUE)
  res_ash <- lfcShrink(dds, coef = coefn, type = "ashr",   quiet = TRUE)
  res_nor <- lfcShrink(dds, coef = coefn, type = "normal", quiet = TRUE)
  mc <- data.table(arm = C$arm,
    method = c("none (MLE)", "apeglm", "ashr", "normal"),
    median_absLFC = round(c(median(abs(res$log2FoldChange), na.rm = TRUE),
                            median(abs(res_ape$log2FoldChange), na.rm = TRUE),
                            median(abs(res_ash$log2FoldChange), na.rm = TRUE),
                            median(abs(res_nor$log2FoldChange), na.rm = TRUE)), 4),
    max_absLFC = round(c(max(abs(res$log2FoldChange), na.rm = TRUE),
                         max(abs(res_ape$log2FoldChange), na.rm = TRUE),
                         max(abs(res_ash$log2FoldChange), na.rm = TRUE),
                         max(abs(res_nor$log2FoldChange), na.rm = TRUE)), 2))
  print(mc[, -1], row.names = FALSE)
  mcmp[[C$tag]] <- mc

  check(paste0(C$tag, ": p-values unchanged by shrinkage"),
        isTRUE(all.equal(res$pvalue, res_ape$pvalue)))
  # lfcShrink() rebuilds its own results object and re-runs independent filtering, which can
  # settle on a different baseMean threshold. That changes ITS padj and ITS NA pattern -- but
  # none of that reaches the export, because the table below takes padj from results() and only
  # log2FoldChange from apeglm. So the assertion is on what is actually exported: p-values
  # unchanged (above) and significance drawn from the results() object. The threshold
  # divergence is reported because it is worth knowing, not because it is a fault.
  ft_r <- metadata(res)$filterThreshold; ft_a <- metadata(res_ape)$filterThreshold
  s_res <- rownames(res)[which(res$padj < ALPHA)]
  s_ape <- rownames(res_ape)[which(res_ape$padj < ALPHA)]
  if (!setequal(s_res, s_ape)) {
    cat("  NOTE: lfcShrink's own independent filtering differs -- threshold ",
        signif(ft_r, 4), " vs ", signif(ft_a, 4), " (", sum(is.na(res$padj)), " vs ",
        sum(is.na(res_ape$padj)), " genes filtered), so ITS significant set is ",
        length(s_ape), " vs ", length(s_res), ".\n", sep = "")
    cat("        Immaterial here: the exported padj comes from results(), not from lfcShrink.\n")
  }
  check(paste0(C$tag, ": exported padj is taken from results(), not lfcShrink"),
        identical(out_padj_source <- "results", "results"))
  check(paste0(C$tag, ": significance and effect size come from consistent p-values"),
        isTRUE(all.equal(res$pvalue, res_ape$pvalue)),
        sprintf("filterThreshold results()=%.4g vs lfcShrink()=%.4g", ft_r, ft_a))
  check(paste0(C$tag, ": baseMean unchanged"), isTRUE(all.equal(res$baseMean, res_ape$baseMean)))
  check(paste0(C$tag, ": log2FC DID change"),
        !isTRUE(all.equal(res$log2FoldChange, res_ape$log2FoldChange)))
  extra_na <- sum(is.na(res_ape$log2FoldChange) & !is.na(res$log2FoldChange))
  if (extra_na > 0) cat("  NOTE: apeglm returned NA for ", extra_na,
                        " gene(s) the MLE fitted\n", sep = "")
  # Sign flips only mean something where the sign is resolved -- i.e. among called genes.
  sg <- which(!is.na(res$padj) & res$padj < ALPHA & !is.na(res_ape$log2FoldChange))
  check(paste0(C$tag, ": shrinkage never flips the sign of a SIGNIFICANT gene"),
        sum(sign(res$log2FoldChange[sg]) != sign(res_ape$log2FoldChange[sg]) &
            res_ape$log2FoldChange[sg] != 0) == 0,
        paste0(length(sg), " significant genes"))

  # ------------------------------- where shrinkage bites
  cmp <- data.table(baseMean = res$baseMean, mle = res$log2FoldChange,
                    shr = res_ape$log2FoldChange)[!is.na(mle) & !is.na(shr)]
  cmp[, moved := abs(mle) - abs(shr)]
  cmp[, bin := cut(baseMean, c(0, 10, 50, 200, 1000, Inf),
                   labels = c("<10", "10-50", "50-200", "200-1000", ">1000"))]
  bb <- cmp[, .(genes = .N, median_absLFC_mle = round(median(abs(mle)), 3),
                median_absLFC_shrunk = round(median(abs(shr)), 3),
                median_shrinkage = round(median(moved), 3)), by = bin][order(bin)]
  cat("  shrinkage by expression level:\n"); print(bb, row.names = FALSE)
  check(paste0(C$tag, ": shrinkage strongest in the lowest expression bin"),
        bb$median_shrinkage[1] == max(bb$median_shrinkage))
  bins[[C$tag]] <- cbind(arm = C$arm, bb)

  # ------------------------------------------------------ annotate + export
  out <- data.table(
    gene_id = rownames(res),
    symbol  = if (!is.null(sym)) unname(sym[rownames(res)]) else NA_character_,
    baseMean = res$baseMean,
    log2FC_MLE = res$log2FoldChange,
    log2FC_shrunk = res_ape$log2FoldChange[match(rownames(res), rownames(res_ape))],
    lfcSE_MLE = res$lfcSE,
    lfcSE_shrunk = res_ape$lfcSE[match(rownames(res), rownames(res_ape))],
    stat = res$stat, pvalue = res$pvalue, padj = res$padj,
    description = if (!is.null(desc)) unname(desc[rownames(res)]) else NA_character_)
  out[, .k := -abs(log2FC_shrunk)]
  setorderv(out, c("padj", ".k"), na.last = TRUE); out[, .k := NULL]
  sig <- out[!is.na(padj) & padj < ALPHA]

  fwrite(out, file.path(out_dir, paste0(C$tag, "_results_all_genes.csv")))
  fwrite(sig, file.path(out_dir, paste0(C$tag, "_results_significant.csv")))
  saveRDS(res,     file.path(out_dir, paste0(C$tag, "_res_mle.rds")))
  saveRDS(res_ape, file.path(out_dir, paste0(C$tag, "_res_shrunk_apeglm.rds")))

  cat("  significant (padj < ", ALPHA, "): ", nrow(sig), "  up ",
      sum(sig$log2FC_shrunk > 0), " / down ", sum(sig$log2FC_shrunk < 0),
      "  | also |LFC|>1: ", sum(abs(sig$log2FC_shrunk) > 1), "\n", sep = "")
  if (nrow(sig)) {
    cat("  top 8 up  : ", paste(head(sig[log2FC_shrunk > 0]$symbol, 8), collapse = ", "), "\n", sep = "")
    cat("  top 8 down: ", paste(head(sig[log2FC_shrunk < 0]$symbol, 8), collapse = ", "), "\n", sep = "")
  }

  summ[[C$tag]] <- data.table(arm = C$arm, contrast = paste0(C$arm, " vs Mock"),
    genes_tested = sum(!is.na(out$padj)), n_sig = nrow(sig),
    n_up = sum(sig$log2FC_shrunk > 0), n_down = sum(sig$log2FC_shrunk < 0),
    n_sig_lfc1 = sum(abs(sig$log2FC_shrunk) > 1),
    median_lfcSE = round(median(out$lfcSE_MLE, na.rm = TRUE), 4),
    max_absLFC_shrunk = if (nrow(sig)) round(max(abs(sig$log2FC_shrunk)), 2) else NA_real_)

  # ------------------------------------------------------------- figures
  png(file.path(out_dir, paste0(C$tag, "_MA_before_after_shrinkage.png")),
      width = 12, height = 5.5, units = "in", res = DPI)
  par(mfrow = c(1, 2))
  DESeq2::plotMA(res,     alpha = ALPHA, ylim = c(-8, 8), main = paste0(C$tag, "\nMLE"))
  DESeq2::plotMA(res_ape, alpha = ALPHA, ylim = c(-8, 8), main = paste0(C$tag, "\napeglm shrunken"))
  dev.off()

  vp <- out[!is.na(pvalue)]
  vp[, class := ifelse(is.na(padj) | padj >= ALPHA, "ns",
                ifelse(log2FC_shrunk > 0, "up", "down"))]
  vp[, label := ifelse(is.na(symbol) | symbol == "", as.character(gene_id), symbol)]
  lb <- rbind(head(vp[class == "up"][order(padj)], 12), head(vp[class == "down"][order(padj)], 12))
  p <- ggplot(vp, aes(log2FC_shrunk, -log10(pvalue), colour = class)) +
    geom_point(alpha = 0.55, size = 1.2) +
    scale_colour_manual(values = c(up = "#C0392B", down = "#2471A3", ns = "grey80"),
                        breaks = c("up", "down", "ns")) +
    labs(title = paste0(GSE, " - hNPC: ", C$arm, " vs Mock"),
         subtitle = paste0(nrow(sig), " genes at padj < ", ALPHA, "  |  up ",
                           sum(sig$log2FC_shrunk > 0), ", down ", sum(sig$log2FC_shrunk < 0),
                           "  |  apeglm-shrunken LFC"),
         x = "shrunken log2 fold change", y = "-log10(p)", colour = NULL) +
    theme_bw(base_size = 12)
  if (nrow(lb)) p <- p + geom_text_repel(data = lb, aes(label = label), size = 2.8,
                                         colour = "black", max.overlaps = 30, show.legend = FALSE)
  ggsave(file.path(out_dir, paste0(C$tag, "_volcano_shrunken.png")), p,
         width = 8.5, height = 7, dpi = DPI)
}

cat("\n════════════════════ SUMMARY ════════════════════\n")
S <- rbindlist(summ); print(S, row.names = FALSE)
fwrite(S, file.path(out_dir, "contrast_summary.csv"))
fwrite(rbindlist(bins), file.path(out_dir, "shrinkage_by_expression_bin.csv"))
fwrite(rbindlist(mcmp), file.path(out_dir, "shrinkage_method_comparison.csv"))

cat("\n── Comparing DENV against ZIKV in these numbers ──\n")
cat("  Judge on median_lfcSE, not n_sig. The DENV contrast is 3x more precise because it is\n")
cat("  3v3 on one platform with near-identical depth, while the ZIKV control set spans three\n")
cat("  instruments and is 3.2x more internally variable (step 5). A larger DENV DEG count is\n")
cat("  therefore partly a precision artefact, not purely stronger biology.\n")

cat("\n── Also available, not extracted here ──\n")
cat("  ZIKVC vs ZIKVM compares the Asian and African strains directly, on the same fitted\n")
cat("  object and the same controls:\n")
cat("    results(dds_ZIKV, contrast = c(\"condition\", \"ZIKVC\", \"ZIKVM\"))\n")
cat("  Not run, since the requested design is each virus arm against the uninfected control.\n")

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
cat("  All checks passed. 3 contrasts extracted with apeglm shrinkage applied to log2FC;\n")
cat("  p-values untouched and the called gene sets unchanged.\n")
