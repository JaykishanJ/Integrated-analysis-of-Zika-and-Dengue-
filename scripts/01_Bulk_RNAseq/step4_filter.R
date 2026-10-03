#!/usr/bin/env Rscript
# STEP 4 -- Low-count filtering. GSE80434.
#
#   Rscript scripts/GSE80434/step4_filter.R
#
#   keep <- rowSums(counts(dds) >= 10) >= min_samples
#   dds  <- dds[keep, ]
#
# min_samples is set PER OBJECT to that object's smallest condition-group size:
#
#   dds_DENV   Mock=3, DENV=3                 -> min_samples = 3
#   dds_ZIKV   Mock=4, ZIKVC=2, ZIKVM=4       -> min_samples = 2  (ZIKVC is the constraint)
#   dds_all    QC object                      -> min_samples = 2
#
# The ZIKV object is worth a word. It carries two contrasts with different group sizes, and the
# threshold has to follow the SMALLER one (ZIKVC, n=2) or a gene expressed only in ZIKVC would
# be deleted before it could ever be tested. That makes the filter more permissive than a
# 4-sample rule would be for the ZIKVM contrast -- the right trade, because DESeq2's own
# independent filtering inside results() then trims low-power genes per contrast anyway.
#
# The filter reads raw counts and never looks at condition labels, so it stays independent of
# the test statistic and does not bias the p-value distribution.
#
# This series exists to compare DENV against ZIKV in one tissue, so the two tested objects
# ending up with different gene universes matters for that comparison. It is quantified below.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table)
})

GSE       <- "GSE80434"
MIN_COUNT <- 10
in_dir    <- file.path("results", GSE, "step3_dds")
out_dir   <- file.path("results", GSE, "step4_filtered")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-60s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 4: LOW-COUNT FILTERING -- ", GSE, " ═══════════\n")

objs <- list(
  dds_all  = list(file = "dds_all_unfitted.rds",  min_samples = 2L,  role = "QC only"),
  dds_DENV = list(file = "dds_DENV_unfitted.rds", min_samples = NA,  role = "tested"),
  dds_ZIKV = list(file = "dds_ZIKV_unfitted.rds", min_samples = NA,  role = "tested"))

summ <- list(); kept <- list()

for (nm in names(objs)) {
  spec <- objs[[nm]]
  path <- file.path(in_dir, spec$file)
  if (!file.exists(path)) stop("missing input: ", path, " -- run step 3 first")
  dds <- readRDS(path)
  n_before <- nrow(dds)
  grp <- table(dds$condition)
  ms  <- if (is.na(spec$min_samples)) min(grp) else spec$min_samples

  cat("\n════════ ", nm, "  (", spec$role, ") ════════\n", sep = "")
  cat("  samples     : ", ncol(dds), "   ",
      paste(names(grp), grp, sep = "=", collapse = ", "), "\n", sep = "")
  cat("  min_count   : ", MIN_COUNT, "\n", sep = "")
  cat("  min_samples : ", ms,
      if (is.na(spec$min_samples)) paste0("  (= smallest condition group, ",
        names(grp)[which.min(grp)], ")") else "  (QC object)", "\n", sep = "")

  keep_g <- rowSums(counts(dds) >= MIN_COUNT) >= ms
  dds_f  <- dds[keep_g, ]
  removed <- counts(dds)[!keep_g, , drop = FALSE]

  cat("  genes       : ", format(n_before, big.mark = ","), " -> ",
      format(sum(keep_g), big.mark = ","),
      sprintf("  (%.1f%% kept)", 100 * mean(keep_g)), "\n", sep = "")
  cat("  removed     : ", format(sum(!keep_g), big.mark = ","),
      "   all-zero: ", format(sum(rowSums(removed) == 0), big.mark = ","),
      " | max count among removed: ", max(removed),
      sprintf(" | reads discarded: %.3f%%", 100 * sum(removed) / sum(counts(dds))), "\n", sep = "")

  # Single-arm genes are the whole reason min_samples tracks the group size.
  arm_only <- rep(FALSE, n_before)
  for (lv in levels(dds$condition)) {
    inn  <- counts(dds)[, dds$condition == lv, drop = FALSE]
    outn <- counts(dds)[, dds$condition != lv, drop = FALSE]
    this <- rowSums(inn >= MIN_COUNT) >= ms & rowSums(outn >= MIN_COUNT) == 0
    if (sum(this)) cat("  expressed only in ", lv, ": ", sum(this), "\n", sep = "")
    arm_only <- arm_only | this
  }
  check(paste0(nm, ": all single-arm genes retained"), all(keep_g[arm_only]),
        paste0(sum(arm_only), " single-arm genes"))

  check(paste0(nm, ": gene count reduced"), nrow(dds_f) < n_before)
  check(paste0(nm, ": sample count unchanged"), ncol(dds_f) == ncol(dds))
  check(paste0(nm, ": counts integer"), is.integer(counts(dds_f)))
  check(paste0(nm, ": counts unchanged for retained genes"),
        identical(counts(dds_f), counts(dds)[keep_g, , drop = FALSE]))
  check(paste0(nm, ": every retained gene meets the rule"),
        all(rowSums(counts(dds_f) >= MIN_COUNT) >= ms))
  check(paste0(nm, ": no all-zero genes remain"), all(rowSums(counts(dds_f)) > 0))
  check(paste0(nm, ": design preserved"),
        identical(deparse(design(dds_f)), deparse(design(dds))))
  check(paste0(nm, ": reference level preserved"),
        levels(dds_f$condition)[1] == levels(dds$condition)[1],
        paste0("baseline = ", levels(dds_f$condition)[1]))
  check(paste0(nm, ": colData still aligned"),
        identical(colnames(counts(dds_f)), rownames(colData(dds_f))))
  check(paste0(nm, ": still unfitted"), is.null(sizeFactors(dds_f)))

  saveRDS(dds_f, file.path(out_dir, paste0(nm, "_filtered.rds")))
  fwrite(data.table(gene_id = rownames(dds_f), as.data.table(counts(dds_f))),
         file.path(out_dir, paste0(nm, "_filtered_counts.tsv")), sep = "\t")
  fwrite(data.table(gene_id = rownames(dds), kept = keep_g,
                    n_samples_ge_10 = rowSums(counts(dds) >= MIN_COUNT),
                    total_count = rowSums(counts(dds))),
         file.path(out_dir, paste0(nm, "_filter_decision_per_gene.csv")))

  kept[[nm]] <- rownames(dds_f)
  summ[[nm]] <- data.table(object = nm, role = spec$role, samples = ncol(dds_f),
                           min_samples = ms, genes_before = n_before,
                           genes_kept = nrow(dds_f),
                           pct_kept = round(100 * nrow(dds_f) / n_before, 1),
                           single_arm_genes = sum(arm_only))
}

cat("\n── Summary ──\n")
S <- rbindlist(summ); print(S, row.names = FALSE)
fwrite(S, file.path(out_dir, "filtering_summary.csv"))

# ------------------------------- the DENV vs ZIKV comparison this series is for
cat("\n── Gene universes: DENV object vs ZIKV object ──\n")
a <- kept$dds_DENV; b <- kept$dds_ZIKV
sh <- intersect(a, b)
cat("  DENV object tested genes : ", format(length(a), big.mark = ","), "\n", sep = "")
cat("  ZIKV object tested genes : ", format(length(b), big.mark = ","), "\n", sep = "")
cat("  shared                   : ", format(length(sh), big.mark = ","), "\n", sep = "")
cat("  DENV only                : ", format(length(setdiff(a, b)), big.mark = ","), "\n", sep = "")
cat("  ZIKV only                : ", format(length(setdiff(b, a)), big.mark = ","), "\n", sep = "")
cat("  Jaccard                  : ", round(length(sh) / length(union(a, b)), 3), "\n", sep = "")
cat("  The point of this series is DENV vs ZIKV in one tissue, so that comparison must be\n")
cat("  made on the ", format(length(sh), big.mark = ","),
    " shared genes -- otherwise a gene missing from one universe reads as a\n", sep = "")
cat("  virus-specific negative when it was simply never tested there.\n")
fwrite(data.table(gene_id = union(a, b), in_DENV = union(a, b) %in% a,
                  in_ZIKV = union(a, b) %in% b),
       file.path(out_dir, "gene_universe_by_subexperiment.csv"))
check("shared universe covers most of both objects",
      length(sh) / min(length(a), length(b)) > 0.85,
      sprintf("%.1f%% of the smaller set", 100 * length(sh) / min(length(a), length(b))))

# The ZIKV object's threshold is driven by ZIKVC (n=2); show what a stricter rule would cost.
cat("\n── Why dds_ZIKV uses min_samples = 2, not 4 ──\n")
z <- readRDS(file.path(in_dir, "dds_ZIKV_unfitted.rds"))
alt <- rbindlist(lapply(c(2L, 3L, 4L), function(k) {
  kk <- rowSums(counts(z) >= MIN_COUNT) >= k
  zc <- counts(z)[, z$condition == "ZIKVC", drop = FALSE]
  other <- counts(z)[, z$condition != "ZIKVC", drop = FALSE]
  zc_only <- rowSums(zc >= MIN_COUNT) >= 2 & rowSums(other >= MIN_COUNT) == 0
  data.table(min_samples = k, genes_kept = sum(kk),
             ZIKVC_only_genes_lost = sum(zc_only & !kk),
             note = if (k == 2L) "<- used" else "")
}))
print(alt, row.names = FALSE)
check("min_samples = 2 loses no ZIKVC-specific gene",
      alt$ZIKVC_only_genes_lost[alt$min_samples == 2L] == 0)

cat("\n── Note ──\n")
cat("  Pre-filter only. DESeq2 applies its own independent filtering inside results(),\n")
cat("  per contrast, on top of this. Both are expected and complementary.\n")

cat("\n── Output ──\n")
for (f in sort(list.files(out_dir))) cat("  ", f, "\n", sep = "")

cat("\n═══════════ RESULT ═══════════\n")
if (length(fails)) {
  cat("  ", length(fails), " CHECK(S) FAILED -- do not proceed:\n", sep = "")
  for (f in fails) cat("    - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("  All checks passed.\n")
for (i in seq_len(nrow(S)))
  cat("  ", S$object[i], ": ", format(S$genes_before[i], big.mark = ","), " -> ",
      format(S$genes_kept[i], big.mark = ","), " genes (>=", MIN_COUNT, " in >=",
      S$min_samples[i], " samples)\n", sep = "")
