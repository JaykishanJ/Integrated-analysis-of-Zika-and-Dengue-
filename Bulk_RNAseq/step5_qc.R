#!/usr/bin/env Rscript
# STEP 5 -- Exploratory QC: VST / rlog, PCA, sample-to-sample distances. GSE80434.
#
#   Rscript scripts/GSE80434/step5_qc.R
#
# Transforms here are for VISUALISATION ONLY. DESeq2 tests raw counts and normalises
# internally; nothing produced in this step feeds the differential expression model. Both use
# blind = TRUE, so the transform ignores the design and any structure seen is genuinely in the
# data rather than something the transform was told to expect.
#
# The decisive question for THIS series is whether sequencing platform is a real batch effect
# in the ZIKV sub-experiment, which spans three instruments. Step 3 carried platform as a batch
# term because it is estimable; this step measures whether it deserves to be there.
#
# That measurement is made on condition-RESIDUALS, not on raw sample distances. When one axis
# holds most of the variance it swamps everything else and a raw-distance check returns a false
# negative -- exactly what happened on GSE110512 before the check was corrected.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
})

GSE     <- "GSE80434"
DPI     <- 600
in_dir  <- file.path("results", GSE, "step4_filtered")
out_dir <- file.path("results", GSE, "step5_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

notes <- character(0)
note <- function(...) { m <- paste0(...); cat(m, "\n"); notes <<- c(notes, m) }
gg <- function(p, f, w = 9, h = 6.5) ggsave(file.path(out_dir, f), p, width = w, height = h, dpi = DPI)

COND_COL <- c(Mock = "#2471A3", DENV = "#8E44AD", ZIKVC = "#C0392B", ZIKVM = "#E67E22")
PLAT_SHP <- c(MiSeq = 16, HiSeq2000 = 17, NextSeq500 = 15)

cat("═══════════ STEP 5: EXPLORATORY QC -- ", GSE, " ═══════════\n")

# ═══════════════════════════════════════════════ A. POOLED (all 16 samples)
cat("\n════════ A. All 16 samples (pooled) ════════\n")
dds <- estimateSizeFactors(readRDS(file.path(in_dir, "dds_all_filtered.rds")))
cat("  ", nrow(dds), " genes x ", ncol(dds), " samples | size factors ",
    paste(round(range(sizeFactors(dds)), 3), collapse = "-"), "\n", sep = "")

vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)
ntd <- normTransform(dds)

sd_trend <- function(m, lab) {
  mu <- rowMeans(m); s <- apply(m, 1, sd)
  q <- cut(mu, quantile(mu, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  data.table(transform = lab, decile = seq_along(tapply(s, q, median)),
             median_sd = round(as.numeric(tapply(s, q, median, na.rm = TRUE)), 3))
}
trend <- rbind(sd_trend(assay(ntd), "log2(norm+1)"), sd_trend(assay(vsd), "VST"),
               sd_trend(assay(rld), "rlog"))
cat("\n── Variance stabilisation check (median SD by expression decile) ──\n")
print(dcast(trend, decile ~ transform, value.var = "median_sd"), row.names = FALSE)
note(sprintf("  low/high SD ratio: log2+1 = %.2f, VST = %.2f (closer to 1 is better)",
             with(trend[transform == "log2(norm+1)"], median_sd[1] / median_sd[10]),
             with(trend[transform == "VST"], median_sd[1] / median_sd[10])))
fwrite(trend, file.path(out_dir, "variance_stabilisation_check.csv"))

pca_of <- function(x, ntop = 500) {
  rv <- rowVars(assay(x)); sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  p <- prcomp(t(assay(x)[sel, ]), center = TRUE, scale. = FALSE)
  list(sc = as.data.table(p$x[, 1:min(5, ncol(p$x))], keep.rownames = "sample"),
       pct = 100 * p$sdev^2 / sum(p$sdev^2))
}
P <- pca_of(vsd); d <- copy(P$sc)
idx <- match(d$sample, colnames(dds))
FACS <- c("condition", "experiment", "platform", "mock_set")
for (cc in FACS) set(d, j = cc, value = as.character(colData(dds)[[cc]])[idx])
set(d, j = "depth_M", value = round(colSums(counts(dds))[idx] / 1e6, 1))
cat("\n── Pooled PCA variance explained ──\n")
cat("  ", paste0("PC", 1:5, "=", round(P$pct[1:5], 1), "%", collapse = "  "), "\n", sep = "")

cat("\n── What each pooled PC tracks (R2) ──\n")
r2 <- sapply(paste0("PC", 1:4), function(p)
  sapply(FACS, function(v) {
    y <- d[[p]]; g <- d[[v]]; ok <- !is.na(g)
    if (length(unique(g[ok])) < 2) return(NA_real_)
    summary(lm(y[ok] ~ factor(g[ok])))$r.squared }))
print(round(r2, 3))
note(sprintf("  pooled PC1 tracks %s most strongly (R2 %.3f); condition R2 on PC1 = %.3f",
             rownames(r2)[which.max(r2[, "PC1"])], max(r2[, "PC1"], na.rm = TRUE),
             r2["condition", "PC1"]))
fwrite(as.data.table(round(r2, 4), keep.rownames = "factor"),
       file.path(out_dir, "pooled_pc_factor_R2.csv"))
fwrite(d, file.path(out_dir, "pooled_pca_scores.csv"))

gg(ggplot(d, aes(PC1, PC2, colour = condition, shape = platform)) +
     geom_point(size = 4.5, alpha = 0.9) +
     geom_text_repel(aes(label = sample), size = 2.3, colour = "grey35", max.overlaps = 25) +
     scale_colour_manual(values = COND_COL) + scale_shape_manual(values = PLAT_SHP) +
     labs(title = paste0(GSE, " - pooled PCA (VST, top 500 variable genes)"),
          subtitle = "all hNPC; two sub-experiments with separate mock sets",
          x = sprintf("PC1 (%.1f%%)", P$pct[1]), y = sprintf("PC2 (%.1f%%)", P$pct[2])) +
     theme_bw(base_size = 12), "pooled_PCA_by_condition.png", 9.5, 7)

gg(ggplot(d, aes(PC1, PC2, colour = experiment, shape = platform)) +
     geom_point(size = 4.5, alpha = 0.9) +
     geom_text_repel(aes(label = sample), size = 2.3, colour = "grey35", max.overlaps = 25) +
     scale_shape_manual(values = PLAT_SHP) +
     labs(title = paste0(GSE, " - pooled PCA coloured by sub-experiment"),
          subtitle = "if the two experiments separate here, pooling their mocks would be wrong",
          x = sprintf("PC1 (%.1f%%)", P$pct[1]), y = sprintf("PC2 (%.1f%%)", P$pct[2])) +
     theme_bw(base_size = 12), "pooled_PCA_by_experiment.png", 9.5, 7)

gg(ggplot(data.table(PC = factor(paste0("PC", 1:6), levels = paste0("PC", 1:6)),
                     pct = P$pct[1:6]), aes(PC, pct)) +
     geom_col(fill = "#34495E") + geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3.4) +
     labs(title = paste0(GSE, " - pooled variance explained"), y = "% variance", x = NULL) +
     ylim(0, max(P$pct) * 1.15) + theme_bw(base_size = 12), "pooled_PCA_scree.png", 7, 5)

# Do the two mock sets actually differ? That is the premise of splitting them.
cat("\n── Do the two mock sets differ from each other? ──\n")
mk <- colnames(dds)[colData(dds)$condition == "Mock"]
dm_all <- as.matrix(dist(t(assay(vsd))))
set_a <- mk[colData(dds)[mk, "mock_set"] == "DENV-Mock"]
set_b <- mk[colData(dds)[mk, "mock_set"] == "Mock"]
wa <- dm_all[set_a, set_a][lower.tri(dm_all[set_a, set_a])]
wb <- dm_all[set_b, set_b][lower.tri(dm_all[set_b, set_b])]
bw <- as.vector(dm_all[set_a, set_b])
note(sprintf("  within DENV-Mock %.1f | within Mock %.1f | between sets %.1f",
             mean(wa), mean(wb), mean(bw)))
# A between/within ratio is the obvious summary but it is misleading here: the two mock sets
# have very different internal spread, so pooling them into one "within" denominator hides
# that. Report the asymmetry, which is the actionable fact.
note(sprintf("  the Mock set is %.1fx more internally variable than DENV-Mock (%.1f vs %.1f)",
             mean(wb) / mean(wa), mean(wb), mean(wa)))
plat_mock <- table(as.character(colData(dds)[set_b, "platform"]))
note(sprintf("  and it spans %d platforms (%s) whereas DENV-Mock is all %s",
             length(plat_mock), paste(names(plat_mock), plat_mock, sep = "x", collapse = ", "),
             unique(as.character(colData(dds)[set_a, "platform"]))))
note(if (mean(wb) / mean(wa) > 2)
  paste0("  -> splitting them was NECESSARY: pooling would have imported the Mock set's\n",
         "     platform-driven spread into the otherwise tight DENV contrast")
  else "  -> the two sets are comparably tight; the split is conservative but not harmful")

lab <- paste0(colData(dds)$condition, " | ", colnames(dds))
dml <- dm_all; rownames(dml) <- colnames(dml) <- lab
ann <- data.frame(condition = colData(dds)$condition, experiment = colData(dds)$experiment,
                  platform = colData(dds)$platform, row.names = lab)
png(file.path(out_dir, "pooled_sample_distance_heatmap.png"), width = 11, height = 9,
    units = "in", res = DPI)
pheatmap(dml, clustering_distance_rows = as.dist(dm_all), clustering_distance_cols = as.dist(dm_all),
         annotation_row = ann, annotation_colors = list(condition = COND_COL),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = paste0(GSE, " - Euclidean distance, all samples (VST)"))
dev.off()
fwrite(as.data.table(dml, keep.rownames = "sample"), file.path(out_dir, "pooled_sample_distances.csv"))
fwrite(data.table(gene_id = rownames(vsd), as.data.table(assay(vsd))),
       file.path(out_dir, "pooled_vst_matrix.csv"))
fwrite(data.table(gene_id = rownames(rld), as.data.table(assay(rld))),
       file.path(out_dir, "pooled_rlog_matrix.csv"))

# ═══════════════════════════════════════════════ B. PER SUB-EXPERIMENT
cat("\n════════ B. Within each sub-experiment ════════\n")
outl <- list()
for (ex in c("DENV", "ZIKV")) {
  cat("\n──────── ", ex, " ────────\n", sep = "")
  dsub <- estimateSizeFactors(readRDS(file.path(in_dir, paste0("dds_", ex, "_filtered.rds"))))
  v <- vst(dsub, blind = TRUE)
  cat("  ", nrow(dsub), " genes x ", ncol(dsub), " samples | design ",
      deparse(design(dsub)), " | size factors ",
      paste(round(range(sizeFactors(dsub)), 3), collapse = "-"), "\n", sep = "")

  Pt <- pca_of(v); dt2 <- copy(Pt$sc)
  i2 <- match(dt2$sample, colnames(dsub))
  for (cc in c("condition", "platform")) set(dt2, j = cc, value = as.character(colData(dsub)[[cc]])[i2])
  cat("  variance explained: ", paste0("PC", 1:min(4, length(Pt$pct)), "=",
      round(Pt$pct[1:min(4, length(Pt$pct))], 1), "%", collapse = "  "), "\n", sep = "")
  r2t <- sapply(paste0("PC", 1:min(3, length(Pt$pct))), function(p)
    sapply(c("condition", "platform"), function(vv) {
      g <- dt2[[vv]]
      if (length(unique(g)) < 2) return(NA_real_)
      summary(lm(dt2[[p]] ~ factor(g)))$r.squared }))
  cat("  what each PC tracks (R2):\n"); print(round(r2t, 3))

  gg(ggplot(dt2, aes(PC1, PC2, colour = condition, shape = platform)) +
       geom_point(size = 5, alpha = 0.9) +
       geom_text_repel(aes(label = sample), size = 2.5, colour = "grey35", max.overlaps = 25) +
       scale_colour_manual(values = COND_COL) + scale_shape_manual(values = PLAT_SHP) +
       labs(title = paste0(GSE, " - ", ex, " sub-experiment PCA (VST)"),
            subtitle = "coloured by condition, shaped by sequencing platform",
            x = sprintf("PC1 (%.1f%%)", Pt$pct[1]), y = sprintf("PC2 (%.1f%%)", Pt$pct[2])) +
       theme_bw(base_size = 12), paste0(ex, "_PCA.png"), 9, 6.5)

  # ---- the platform question, measured on condition-residuals
  cat("\n  ── platform / batch check (on condition-residuals) ──\n")
  if (nlevels(droplevels(colData(dsub)$platform)) < 2) {
    note(sprintf("  [%s] single platform (%s) -- no batch effect possible, contrast is matched",
                 ex, levels(droplevels(colData(dsub)$platform))))
  } else {
    vm <- assay(v); rv <- rowVars(vm)
    top <- vm[order(rv, decreasing = TRUE)[seq_len(min(500, nrow(vm)))], ]
    cf <- colData(dsub)$condition; pf <- droplevels(colData(dsub)$platform)
    resid <- t(apply(top, 1, function(g) residuals(lm(g ~ cf))))
    r2p <- apply(resid, 1, function(g) summary(lm(g ~ pf))$r.squared)
    note(sprintf("  [%s] median R2 of platform on condition-residuals: %.3f (%d/500 genes > 0.5)",
                 ex, median(r2p), sum(r2p > 0.5)))
    note(sprintf("  [%s] -> platform batch term is %s", ex,
                 if (median(r2p) > 0.4) "well supported by the data"
                 else "weakly supported; step 6 must test whether it earns its df"))
  }

  # ---- distances + arm separation
  sd_t <- dist(t(assay(v))); dmt <- as.matrix(sd_t)
  labt <- paste0(colData(dsub)$condition, " | ", colData(dsub)$platform)
  dl <- dmt; rownames(dl) <- colnames(dl) <- make.unique(labt)
  png(file.path(out_dir, paste0(ex, "_sample_distance_heatmap.png")),
      width = 8.5, height = 7, units = "in", res = DPI)
  pheatmap(dl, clustering_distance_rows = sd_t, clustering_distance_cols = sd_t,
           annotation_row = data.frame(condition = colData(dsub)$condition,
                                       row.names = make.unique(labt)),
           annotation_colors = list(condition = COND_COL),
           col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
           display_numbers = TRUE, number_format = "%.0f", fontsize_number = 7,
           main = paste0(GSE, " - ", ex, ": Euclidean distance (VST)"))
  dev.off()

  cat("\n  ── separation between arms ──\n")
  cf <- colData(dsub)$condition
  for (lv in setdiff(levels(droplevels(cf)), "Mock")) {
    i1 <- which(cf == "Mock"); i2b <- which(cf == lv)
    wi <- c(dmt[i1, i1][lower.tri(dmt[i1, i1])], dmt[i2b, i2b][lower.tri(dmt[i2b, i2b])])
    bwn <- as.vector(dmt[i1, i2b])
    note(sprintf("  [%s] %-6s vs Mock: within %.1f, between %.1f, ratio %.2fx",
                 ex, lv, mean(wi), mean(bwn), mean(bwn) / mean(wi)))
  }

  med <- apply(dmt, 1, function(x) median(x[x > 0]))
  o <- data.table(experiment = ex, sample = colnames(dsub),
                  condition = as.character(cf),
                  platform = as.character(colData(dsub)$platform),
                  median_dist = round(med, 1),
                  depth_M = round(colSums(counts(dsub)) / 1e6, 2))
  o[, z := round((median_dist - mean(median_dist)) / sd(median_dist), 2)]
  print(o[, -1], row.names = FALSE)
  note(sprintf("  [%s] max |z| of median distance: %.2f -- %s", ex, max(abs(o$z)),
               if (max(abs(o$z)) < 2) "no outlier flagged" else "inspect flagged sample"))
  outl[[ex]] <- o

  fwrite(data.table(gene_id = rownames(v), as.data.table(assay(v))),
         file.path(out_dir, paste0(ex, "_vst_matrix.csv")))
  fwrite(dt2, file.path(out_dir, paste0(ex, "_pca_scores.csv")))
}

fwrite(rbindlist(outl), file.path(out_dir, "outlier_screen.csv"))
writeLines(notes, file.path(out_dir, "qc_summary.txt"))

stamp <- file.path("scripts", "set_png_dpi.R")
if (file.exists(stamp)) {
  cat("\n── Figure resolution ──\n")
  cat(system2("Rscript", c(stamp, shQuote(out_dir), DPI), stdout = TRUE, stderr = TRUE), sep = "\n")
}

cat("\n── Output ──\n")
for (f in sort(list.files(out_dir))) cat("  ", f, "\n", sep = "")
cat("\n═══════════ RESULT ═══════════\n")
cat("  VST and rlog computed (blind = TRUE); pooled and per-sub-experiment PCA plus distance\n")
cat("  heatmaps written. Visualisation ONLY -- DESeq2 uses the raw counts.\n")
