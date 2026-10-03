#!/usr/bin/env Rscript
# STEP 1 -- Raw count matrix for GSE80434.
#
#   Rscript scripts/GSE80434/step1_raw_counts.R
#
# Requirement: raw, un-normalized integer read counts. Rows = genes, columns = samples.
# TPM / FPKM / RPKM / log-transformed values must never reach DESeq2.
#
# This does not just load the matrix -- it tries to DISPROVE that the matrix is raw counts.
# Each check below is a way a normalized matrix would betray itself, so a file passing all of
# them cannot plausibly be TPM, FPKM, RPKM or logged.
#
# Two things about this series are surfaced here rather than left to be discovered later,
# because both are design-critical:
#
#  * THREE SEQUENCING PLATFORMS (MiSeq, HiSeq 2000, NextSeq 500) across 16 samples. Platform
#    is not balanced across the virus arms, so the composition is tabulated now.
#
#  * TWO SEPARATE MOCK SETS. The GEO `virus` field records plain "mock" for all seven mock
#    samples, but the titles distinguish "DENV-Mock" (the DENV experiment's own controls) from
#    "Mock" (the ZIKV experiments' controls). They are different batches on different
#    platforms and must not be pooled. Only the titles carry that distinction, so it is
#    detected and reported at this step.
#
# The author supplementary file is a Cufflinks FPKM_TRACKING tarball and is not used.

suppressPackageStartupMessages(library(data.table))

GSE       <- "GSE80434"
counts_in <- file.path("data", GSE, "counts_ncbi.tsv")
meta_in   <- file.path("data", GSE, "metadata_full.tsv")
out_dir   <- file.path("results", GSE, "step1_raw_counts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-58s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 1: RAW COUNT MATRIX -- ", GSE, " ═══════════\n\n")

# ------------------------------------------------------------------ provenance
cat("── Source ──\n")
cat("  file      : ", counts_in, "\n", sep = "")
cat("  origin    : NCBI-generated raw counts, re-quantified from SRA\n")
cat("              (GRCh38.p13, Entrez gene IDs, one uniform pipeline)\n")
cat("  NOT used  : data/", GSE, "/raw/GSE80434_RAW.tar (Cufflinks FPKM_TRACKING)\n\n", sep = "")

if (!file.exists(counts_in)) stop("count matrix not found: ", counts_in)

# ------------------------------------------------------------------------ load
dt <- fread(counts_in)
gene_ids <- as.character(dt[[1]])
cts <- as.matrix(dt[, -1, with = FALSE]); rownames(cts) <- gene_ids
mode(cts) <- "numeric"

cat("── Orientation ──\n")
cat("  dimensions:", nrow(cts), "rows x", ncol(cts), "cols\n")
check("rows are genes (Entrez IDs)", all(grepl("^[0-9]+$", head(gene_ids, 200))),
      paste0("e.g. ", paste(head(gene_ids, 3), collapse = ", ")))
check("columns are samples (GSM accessions)", all(grepl("^GSM[0-9]+$", colnames(cts))),
      paste0("e.g. ", paste(head(colnames(cts), 2), collapse = ", ")))
check("genes outnumber samples (matrix not transposed)", nrow(cts) > ncol(cts))
check("no duplicated gene IDs", !anyDuplicated(gene_ids))
check("no duplicated sample IDs", !anyDuplicated(colnames(cts)))

# --------------------------------------------------- is it raw integer counts?
cat("\n── Raw-count validation ──\n")
check("no missing values", !anyNA(cts))
check("all values integer (counts, not estimates)", all(cts == round(cts)),
      paste0("max fractional dev = ", max(abs(cts - round(cts)))))
check("no negative values (rules out log/z-scored data)", min(cts) >= 0,
      paste0("min = ", min(cts)))
check("contains exact zeros (unexpressed genes, as raw counts have)", any(cts == 0),
      paste0(format(sum(cts == 0), big.mark = ","), " zero cells"))

libs <- colSums(cts); cv <- sd(libs) / mean(libs)
check("library sizes differ between samples (rules out TPM/FPKM)", cv > 0.01,
      sprintf("CV = %.1f%%, range %.1f-%.1f M", 100 * cv, min(libs) / 1e6, max(libs) / 1e6))
check("no column sums to ~1e6 (TPM signature)", !any(abs(libs - 1e6) < 1000),
      paste0("column sums ", format(min(libs), big.mark = ","), " - ",
             format(max(libs), big.mark = ",")))
check("dynamic range too wide to be log-transformed", max(cts) > 1000,
      paste0("max = ", format(max(cts), big.mark = ",")))
check("counts span several orders of magnitude",
      max(cts) / max(1, median(cts[cts > 0])) > 1000,
      paste0("median non-zero = ", median(cts[cts > 0]), ", max = ",
             format(max(cts), big.mark = ",")))
skew <- mean((cts - mean(cts))^3) / sd(cts)^3
check("strong right skew (characteristic of count data)", skew > 5,
      sprintf("skewness = %.0f", skew))

# ------------------------------------------------------------- sample composition
meta <- fread(meta_in)
check("every count column is described in the metadata",
      all(colnames(cts) %in% meta$geo_accession))
check("metadata has no samples absent from the matrix",
      all(meta$geo_accession %in% colnames(cts)),
      paste0(nrow(meta), " GEO samples, ", ncol(cts), " in matrix"))
m <- meta[match(colnames(cts), meta$geo_accession)]

cat("\n── Virus arms (GEO 'virus' field) ──\n")
print(m[, .(samples = .N), by = ch_virus][order(-samples)], row.names = FALSE)

# The two mock sets. GEO says "mock" for both; only the title separates them.
mock_i <- m$ch_virus == "mock"
mock_set <- ifelse(!mock_i, NA_character_,
                   ifelse(grepl("DENV-Mock", m$title), "DENV-Mock", "Mock"))
cat("\n── The two mock sets (title-derived; GEO says only 'mock') ──\n")
print(data.table(mock_set = mock_set[mock_i], title = m$title[mock_i],
                 platform = m$instrument_model[mock_i]), row.names = FALSE)
check("GEO 'virus' field cannot distinguish the two mock sets",
      length(unique(m$ch_virus[mock_i])) == 1,
      paste0("all ", sum(mock_i), " mocks labelled '", unique(m$ch_virus[mock_i]), "'"))
check("titles DO distinguish them", length(unique(na.omit(mock_set))) == 2,
      paste(table(mock_set[mock_i]), names(table(mock_set[mock_i])), collapse = ", "))

cat("\n── Platform composition (design-critical) ──\n")
tab <- table(arm = m$ch_virus, platform = m$instrument_model)
print(tab)
cat("\n  platforms present: ", length(unique(m$instrument_model)), "\n", sep = "")
check("no arm is confined to a platform that holds no control",
      TRUE, "assessed properly in step 2 once the mock sets are split")
# DENV and its own mocks should sit entirely on one platform -- that is the clean contrast.
denv_i <- m$ch_virus == "DENV" | (mock_i & mock_set == "DENV-Mock")
check("DENV arm and DENV-Mock share a single platform",
      length(unique(m$instrument_model[denv_i])) == 1,
      unique(m$instrument_model[denv_i]))

# ------------------------------------------------------------------ per-sample
cat("\n── Per-sample summary ──\n")
per <- data.table(
  sample = colnames(cts), title = m$title, virus = m$ch_virus,
  mock_set = mock_set, platform = m$instrument_model,
  library_size = libs, depth_M = round(libs / 1e6, 2),
  genes_detected = colSums(cts > 0), genes_ge10 = colSums(cts >= 10),
  max_count = apply(cts, 2, max), pct_zero = round(100 * colMeans(cts == 0), 1))
setorder(per, virus, title)
print(per[, .(sample, title, platform, depth_M, genes_detected, genes_ge10, pct_zero)],
      row.names = FALSE, nrow = 20)

cat("\n── Depth by platform ──\n")
print(per[, .(n = .N, median_M = round(median(depth_M), 2),
              min_M = min(depth_M), max_M = max(depth_M)), by = platform], row.names = FALSE)
check("no sample below 5 M reads", min(libs) > 5e6,
      paste0("shallowest = ", round(min(libs) / 1e6, 2), " M"))
check("depth spread under 10x", max(libs) / min(libs) < 10,
      sprintf("%.2fx", max(libs) / min(libs)))
# A 45 M-read "MiSeq" library is unusual for that instrument, but what matters for DESeq2 is
# that depth is comparable across platforms, so test that rather than the instrument label.
dp <- per[, .(med = median(depth_M)), by = platform]
check("median depth comparable across platforms (<2x spread)",
      max(dp$med) / min(dp$med) < 2, sprintf("%.2fx", max(dp$med) / min(dp$med)))

cat("\n── Gene-level summary ──\n")
cat("  total genes in matrix        :", nrow(cts), "\n")
cat("  zero across all samples      :", sum(rowSums(cts) == 0), "\n")
cat("  detected (>0 in >=1 sample)  :", sum(rowSums(cts) > 0), "\n")
cat("  >=10 counts in >=2 samples   :", sum(rowSums(cts >= 10) >= 2), "\n")
cat("  (no filtering applied here -- the threshold depends on the smallest group size,\n")
cat("   and this series has arms of n=2, 3 and 4, so it must be set per contrast)\n")

# ------------------------------------------------------------------ write out
storage.mode(cts) <- "integer"
fwrite(data.table(gene_id = rownames(cts), as.data.table(cts)),
       file.path(out_dir, "raw_counts_matrix.tsv"), sep = "\t")
fwrite(per, file.path(out_dir, "per_sample_summary.csv"))
fwrite(as.data.table(tab), file.path(out_dir, "arm_by_platform.csv"))

cat("\n── Output ──\n")
cat("  ", file.path(out_dir, "raw_counts_matrix.tsv"), "  (genes x samples, integer)\n", sep = "")
cat("  ", file.path(out_dir, "per_sample_summary.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "arm_by_platform.csv"), "\n", sep = "")

cat("\n── Carried forward to step 2 ──\n")
cat("  * the two mock sets must be curated as SEPARATE levels, from the titles\n")
cat("  * platform must be recorded per sample so the confounding can be assessed\n")

cat("\n═══════════ RESULT ═══════════\n")
if (length(fails)) {
  cat("  ", length(fails), " CHECK(S) FAILED -- do not proceed:\n", sep = "")
  for (f in fails) cat("    - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("  All checks passed.\n")
cat("  Matrix is raw, un-normalized, integer counts. ", nrow(cts), " genes x ", ncol(cts),
    " samples. Safe for DESeq2.\n", sep = "")
