#!/usr/bin/env Rscript
# STEP 6 -- Run the DESeq2 pipeline. GSE80434.
#
#   Rscript scripts/GSE80434/step6_deseq.R
#
#   dds <- DESeq(dds)
#
# DESeq() wraps three stages, and each object is fitted BOTH ways -- once via the wrapper and
# once stage by stage -- then checked for identity. That makes the "three automated steps"
# verifiable rather than asserted.
#
#   1. estimateSizeFactors  -- median-of-ratios library size normalisation
#   2. estimateDispersions  -- gene-wise MLE, parametric trend fit, empirical-Bayes shrinkage
#   3. nbinomWaldTest       -- negative binomial GLM fit + Wald test per coefficient
#
# Designs differ between the two objects, for reasons established in steps 2-5:
#
#   dds_DENV   ~ condition              all 6 samples on MiSeq, so no batch term exists.
#                                       Step 5: PC1 = 99.5% pure condition, separation 6.18x.
#   dds_ZIKV   ~ platform + condition   3 platforms; step 5 measured platform at median R2
#                                       0.449 on condition-residuals, and the Mock set is 3.2x
#                                       more internally variable than DENV-Mock. The batch term
#                                       is still checked empirically below.
#
# WATCH OUT: dds_ZIKV has TWO condition coefficients, so results() with no arguments returns
# the LAST one -- conditionZIKVM, the African strain. Both are named explicitly.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(ggplot2)
})

GSE     <- "GSE80434"
DPI     <- 600
in_dir  <- file.path("results", GSE, "step4_filtered")
out_dir <- file.path("results", GSE, "step6_deseq")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(warn = 1)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-62s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 6: DESeq2 PIPELINE -- ", GSE, " ═══════════\n")

objs <- c(DENV = "dds_DENV_filtered.rds", ZIKV = "dds_ZIKV_filtered.rds")
sf_tab <- list(); disp_tab <- list(); sens_tab <- list(); summ <- list()

for (ex in names(objs)) {
  path <- file.path(in_dir, objs[[ex]])
  if (!file.exists(path)) stop("missing input: ", path, " -- run step 4 first")
  dds <- readRDS(path)

  cat("\n════════════════════ ", ex, " ════════════════════\n", sep = "")
  mm <- model.matrix(design(dds), colData(dds))
  cat("  ", nrow(dds), " genes x ", ncol(dds), " samples\n", sep = "")
  cat("  design      : ", deparse(design(dds)),
      if (length(attr(terms(design(dds)), "term.labels")) > 1) "   (batch first, condition last)"
      else "   (no batch term -- platform is constant)", "\n", sep = "")
  cat("  reference   : ", levels(dds$condition)[1], "\n", sep = "")
  cat("  coefficients: ", ncol(mm), " -> residual df = ", ncol(dds) - ncol(mm), "\n", sep = "")
  check(paste0(ex, ": starts unfitted"), is.null(sizeFactors(dds)))
  check(paste0(ex, ": model matrix full rank"), qr(mm)$rank == ncol(mm))

  # ══════════════════════════════════════ the wrapper
  cat("\n  ── DESeq(dds) ──\n")
  t0 <- proc.time()[["elapsed"]]
  dds <- DESeq(dds, quiet = TRUE)
  cat("     completed in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep = "")

  # ══════════════════════════════════════ stage 1
  cat("\n  ── Stage 1: size factor normalisation ──\n")
  sf <- sizeFactors(dds)
  s1 <- data.table(experiment = ex, sample = names(sf),
                   condition = as.character(dds$condition),
                   platform = as.character(dds$platform),
                   size_factor = round(sf, 4),
                   library_M = round(colSums(counts(dds)) / 1e6, 2))
  print(s1[, -1], row.names = FALSE)
  cat("     method: median-of-ratios | range ", paste(round(range(sf), 3), collapse = "-"),
      " | geometric mean ", round(exp(mean(log(sf))), 4), "\n", sep = "")
  check(paste0(ex, ": size factors present and finite"), !anyNA(sf) && all(is.finite(sf)))
  check(paste0(ex, ": size factors positive"), all(sf > 0))
  check(paste0(ex, ": size factors centred near 1"), abs(exp(mean(log(sf))) - 1) < 0.05)
  check(paste0(ex, ": normalised counts == raw / size factor"),
        max(abs(counts(dds, normalized = TRUE) - sweep(counts(dds), 2, sf, "/"))) < 1e-8)
  sf_tab[[ex]] <- s1

  # ══════════════════════════════════════ stage 2
  cat("\n  ── Stage 2: dispersion estimation & shrinkage ──\n")
  md <- mcols(dds)
  cat("     gene-wise MLE : ", sum(!is.na(md$dispGeneEst)), " genes (median ",
      round(median(md$dispGeneEst, na.rm = TRUE), 5), ")\n", sep = "")
  cat("     fitted trend  : ", sum(!is.na(md$dispFit)), " genes (median ",
      round(median(md$dispFit, na.rm = TRUE), 5), ")\n", sep = "")
  cat("     final (MAP)   : ", sum(!is.na(md$dispersion)), " genes (median ",
      round(median(md$dispersion, na.rm = TRUE), 5), ")\n", sep = "")
  cat("     trend type    : ", attr(dispersionFunction(dds), "fitType"),
      " | dispersion outliers left unshrunk: ", sum(md$dispOutlier, na.rm = TRUE), "\n", sep = "")
  check(paste0(ex, ": gene-wise dispersions estimated"),
        sum(!is.na(md$dispGeneEst)) > 0.9 * nrow(dds))
  check(paste0(ex, ": final dispersions present for all genes"), !anyNA(md$dispersion))
  check(paste0(ex, ": dispersions positive"), all(md$dispersion > 0, na.rm = TRUE))
  sd_pre <- sd(log(md$dispGeneEst), na.rm = TRUE); sd_post <- sd(log(md$dispersion), na.rm = TRUE)
  check(paste0(ex, ": shrinkage reduced dispersion spread"), sd_post < sd_pre,
        sprintf("sd(log disp) %.3f -> %.3f", sd_pre, sd_post))
  disp_tab[[ex]] <- data.table(experiment = ex, gene_id = rownames(dds),
                               dispGeneEst = md$dispGeneEst, dispFit = md$dispFit,
                               dispersion = md$dispersion, dispOutlier = md$dispOutlier,
                               betaConv = md$betaConv)

  png(file.path(out_dir, paste0(ex, "_dispersion_estimates.png")),
      width = 7.5, height = 6, units = "in", res = DPI)
  plotDispEsts(dds, main = paste0(GSE, " - ", ex, " dispersion\n", deparse(design(dds))))
  dev.off()

  # ══════════════════════════════════════ stage 3
  cat("\n  ── Stage 3: negative binomial GLM + Wald test ──\n")
  rn <- resultsNames(dds)
  for (n in rn) cat("     - ", n, "\n", sep = "")
  cond_coefs <- grep("^condition", rn, value = TRUE)
  cat("     condition coefficients: ", paste(cond_coefs, collapse = ", "), "\n", sep = "")
  if (length(cond_coefs) > 1)
    cat("     results() default would return: ", tail(rn, 1),
        "  <- NOT necessarily the arm of interest\n", sep = "")
  conv <- sum(md$betaConv, na.rm = TRUE)
  cat("     GLM converged : ", conv, " of ", nrow(dds), " genes\n", sep = "")
  check(paste0(ex, ": GLM converged for >99% of genes"), conv > 0.99 * nrow(dds))
  check(paste0(ex, ": Wald statistics computed"), any(grepl("^WaldStatistic", names(md))))

  for (cf in cond_coefs) {
    r <- results(dds, name = cf, alpha = 0.05)
    cat("     ", cf, ": ", sum(r$padj < 0.05, na.rm = TRUE), " genes at padj<0.05  (up ",
        sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE), ", down ",
        sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE), ")\n", sep = "")
    summ[[paste(ex, cf)]] <- data.table(
      experiment = ex, coefficient = cf, design = paste(deparse(design(dds)), collapse = ""),
      residual_df = ncol(dds) - ncol(mm), genes_tested = sum(!is.na(r$padj)),
      n_sig = sum(r$padj < 0.05, na.rm = TRUE),
      n_up = sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
      n_down = sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE),
      median_lfcSE = round(median(r$lfcSE, na.rm = TRUE), 4))
  }

  # ══════════════════════════ wrapper == manual stages?
  cat("\n  ── Wrapper vs manual stages ──\n")
  m <- readRDS(path)
  m <- estimateSizeFactors(m); m <- estimateDispersions(m, quiet = TRUE); m <- nbinomWaldTest(m)
  check(paste0(ex, ": size factors identical"),
        isTRUE(all.equal(sizeFactors(m), sizeFactors(dds))))
  check(paste0(ex, ": dispersions identical"),
        isTRUE(all.equal(mcols(m)$dispersion, mcols(dds)$dispersion)))
  check(paste0(ex, ": log2FC identical"),
        isTRUE(all.equal(results(m, name = cond_coefs[1])$log2FoldChange,
                         results(dds, name = cond_coefs[1])$log2FoldChange)))
  check(paste0(ex, ": p-values identical"),
        isTRUE(all.equal(results(m, name = cond_coefs[1])$pvalue,
                         results(dds, name = cond_coefs[1])$pvalue)))

  # ══════════════════════════ sensitivity: does the batch term earn its df?
  has_batch <- length(attr(terms(design(dds)), "term.labels")) > 1
  if (!has_batch) {
    cat("\n  ── Sensitivity: not applicable ──\n")
    cat("     platform is constant across these 6 samples, so there is no alternative\n")
    cat("     model to compare against. The contrast is platform-matched by construction.\n")
  } else {
    cat("\n  ── Sensitivity: ~ condition (platform dropped) ──\n")
    u <- readRDS(path); design(u) <- ~ condition
    u <- DESeq(u, quiet = TRUE)
    mm_u <- model.matrix(design(u), colData(u))
    for (cf in cond_coefs) {
      rb <- results(dds, name = cf, alpha = 0.05); ru <- results(u, name = cf, alpha = 0.05)
      sens_tab[[paste(ex, cf)]] <- data.table(
        experiment = ex, coefficient = cf,
        blocked_resid_df = ncol(dds) - ncol(mm), unblocked_resid_df = ncol(u) - ncol(mm_u),
        blocked_median_lfcSE = round(median(rb$lfcSE, na.rm = TRUE), 4),
        unblocked_median_lfcSE = round(median(ru$lfcSE, na.rm = TRUE), 4),
        blocked_n_sig = sum(rb$padj < 0.05, na.rm = TRUE),
        unblocked_n_sig = sum(ru$padj < 0.05, na.rm = TRUE))
    }
    st <- rbindlist(sens_tab[grep(paste0("^", ex, " "), names(sens_tab))])
    print(st[, -1], row.names = FALSE)
    better <- all(st$blocked_median_lfcSE < st$unblocked_median_lfcSE)
    cat("     -> on median lfcSE alone, keeping platform ",
        if (better) "improved" else "did NOT improve", " every coefficient\n", sep = "")

    # lfcSE is the right criterion only when the batch factor is BALANCED across conditions.
    # If some batch level occurs in just one arm, its effect cannot be separated from that arm
    # unless it is in the model -- and then omitting it biases the estimate rather than merely
    # costing efficiency. Balance therefore decides first, and it overrides lfcSE.
    tb <- table(droplevels(colData(dds)$condition), droplevels(colData(dds)$platform))
    confined <- colnames(tb)[colSums(tb > 0) == 1]
    cat("\n     batch balance across conditions:\n")
    print(tb)
    if (length(confined)) {
      for (pl in confined) {
        arm <- rownames(tb)[tb[, pl] > 0]
        cat("     ", pl, " occurs ONLY in ", arm, " (n=", tb[arm, pl], " of ", sum(tb[arm, ]),
            " in that arm)\n", sep = "")
      }
      cat("     -> platform is UNBALANCED. Dropping it would fold the ", paste(confined, collapse = "/"),
          "\n        effect into the ", rownames(tb)[tb[, confined[1]] > 0],
          " baseline and bias every contrast against it.\n", sep = "")
      cat("        The batch term is kept on bias grounds, whatever median lfcSE says.\n")
    }
    check(paste0(ex, ": batch decision is defensible (balanced -> lfcSE, unbalanced -> keep)"),
          length(confined) > 0 || better,
          if (length(confined)) paste0("unbalanced: ", paste(confined, collapse = ", "),
                                       " confined to one arm -> keep platform")
          else "balanced, and lfcSE favours blocking")
    png(file.path(out_dir, paste0(ex, "_dispersion_design_comparison.png")),
        width = 12, height = 5.5, units = "in", res = DPI)
    par(mfrow = c(1, 2))
    plotDispEsts(dds, main = paste0(ex, ": ~ platform + condition  (primary)"))
    plotDispEsts(u,   main = paste0(ex, ": ~ condition  (sensitivity)"))
    dev.off()
    saveRDS(u, file.path(out_dir, paste0("dds_", ex, "_fitted_unblocked_sensitivity.rds")))
  }

  saveRDS(dds, file.path(out_dir, paste0("dds_", ex, "_fitted.rds")))
}

# ------------------------------------------------------------------------ summary
cat("\n════════════════════ SUMMARY ════════════════════\n")
S <- rbindlist(summ); print(S, row.names = FALSE)
fwrite(S, file.path(out_dir, "fit_summary.csv"))
fwrite(rbindlist(sf_tab), file.path(out_dir, "size_factors.csv"))
fwrite(rbindlist(disp_tab), file.path(out_dir, "dispersions.csv"))
if (length(sens_tab)) fwrite(rbindlist(sens_tab), file.path(out_dir, "design_sensitivity.csv"))

cat("\n── Reading these numbers ──\n")
cat("  The DENV and ZIKV contrasts are not equally powered, and the gap is partly technical.\n")
cat("  Step 5 found the ZIKV Mock set 3.2x more internally variable than DENV-Mock, because it\n")
cat("  spans three sequencing platforms while DENV-Mock is all MiSeq. Compare median lfcSE\n")
cat("  above rather than DEG counts alone when judging DENV against ZIKV here.\n")

cat("\n── Note on results() defaults ──\n")
cat("  dds_ZIKV has 3 condition levels, so resultsNames ends with conditionZIKVM.\n")
cat("  results(dds) with no arguments returns the ZIKVM (African strain) contrast.\n")
cat("  Step 7 must name each coefficient explicitly.\n")

stamp <- file.path("scripts", "set_png_dpi.R")
if (file.exists(stamp)) {
  cat("\n── Figure resolution ──\n")
  cat(system2("Rscript", c(stamp, shQuote(out_dir), DPI), stdout = TRUE, stderr = TRUE), sep = "\n")
}

cat("\n── Output ──\n")
for (f in sort(list.files(out_dir))) cat("  ", f, "\n", sep = "")

cat("\n═══════════ RESULT ═══════════\n")
if (length(fails)) {
  cat("  ", length(fails), " CHECK(S) FAILED:\n", sep = "")
  for (f in fails) cat("    - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("  All checks passed. Three stages completed and verified against manual execution\n")
cat("  for both sub-experiments. Baseline Mock throughout.\n")
