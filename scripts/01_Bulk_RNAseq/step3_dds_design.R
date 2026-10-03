#!/usr/bin/env Rscript
# STEP 3 -- Build the DESeqDataSet(s) and fix the design formulae. GSE80434.
#
#   Rscript scripts/GSE80434/step3_dds_design.R
#
# Convention: confounding / batch terms first, biological condition LAST.
#
# This series splits into two objects because it contains two sub-experiments with SEPARATE
# control sets (step 2). Fitting them together would mean pooling the two mock groups, which
# sit on different sequencing platforms -- that would import batch variation into the DENV
# contrast, which is otherwise perfectly platform-matched.
#
#   dds_all    ~ experiment + condition   assembled for QC / normalisation overview ONLY
#   dds_DENV   ~ condition                 3 DENV vs 3 DENV-Mock, all MiSeq -> platform is
#                                          constant, so there is no batch term to add and the
#                                          contrast is matched by construction
#   dds_ZIKV   ~ platform + condition      ZIKVC (2) + ZIKVM (4) vs Mock (4) across 3
#                                          platforms; step 2 showed platform is estimable
#                                          here, not aliased with condition
#
# dds_all is deliberately NOT used for testing: it would force one shared condition effect
# across two experiments run at different times on different instruments.

suppressPackageStartupMessages({
  library(DESeq2); library(data.table)
})

GSE       <- "GSE80434"
CONTROL   <- "Mock"
MIN_RESID <- 1L
counts_in <- file.path("results", GSE, "step1_raw_counts", "raw_counts_matrix.tsv")
cd_in     <- file.path("results", GSE, "step2_metadata", "colData.rds")
out_dir   <- file.path("results", GSE, "step3_dds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-60s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 3: DESeqDataSet & DESIGN -- ", GSE, " ═══════════\n\n")

for (p in c(counts_in, cd_in)) if (!file.exists(p)) stop("missing input: ", p)
dt  <- fread(counts_in)
cts <- as.matrix(dt[, -1, with = FALSE]); rownames(cts) <- as.character(dt[[1]])
storage.mode(cts) <- "integer"
# Read the .rds, not the .csv: the factor level order carries the reference level, and here a
# CSV round-trip would reset it to alphabetical -- making DENV the baseline and flipping signs.
cd <- readRDS(cd_in)

cat("── Inputs ──\n")
cat("  countData : ", nrow(cts), " genes x ", ncol(cts), " samples\n", sep = "")
cat("  colData   : ", nrow(cd), " samples x ", ncol(cd), " factors\n", sep = "")
cat("  reference : ", paste(levels(cd$condition), collapse = ", "),
    " -> baseline \"", levels(cd$condition)[1], "\"\n\n", sep = "")

cat("── Re-verify alignment before assembly ──\n")
check("all(colnames(counts) == rownames(colData))", all(colnames(cts) == rownames(cd)))
check("identical(colnames(counts), rownames(colData))", identical(colnames(cts), rownames(cd)))
check("reference level survived the round-trip", levels(cd$condition)[1] == CONTROL,
      paste0("baseline = ", levels(cd$condition)[1]))
check("counts still integer", is.integer(cts))

# ------------------------------------- why the two sub-experiments must be split
cat("\n── Why testing is split by sub-experiment ──\n")
cat("  the 7 mock samples do not form one group:\n")
print(as.data.table(cd)[condition == "Mock",
      .(n = .N, platforms = paste(sort(unique(as.character(platform))), collapse = ", ")),
      by = .(mock_set, experiment)], row.names = FALSE)
check("the two mock sets sit on different platform mixes",
      !identical(sort(unique(as.character(cd$platform[cd$mock_set == "DENV-Mock" &
                                                     !is.na(cd$mock_set)]))),
                 sort(unique(as.character(cd$platform[cd$mock_set == "Mock" &
                                                     !is.na(cd$mock_set)])))))
cat("  -> pooling them would put NextSeq/HiSeq variation into the all-MiSeq DENV contrast.\n")

# --------------------------------------------------------------- helper: build + report
constant_cols <- function(d) names(d)[vapply(d, function(x) length(unique(x)) == 1L, logical(1))]

build <- function(label, keep_rows, formula) {
  sub_cd <- droplevels(cd[keep_rows, , drop = FALSE])
  # droplevels can silently move the reference if the baseline level were dropped; it is not
  # here, but assert rather than trust.
  if ("condition" %in% names(sub_cd) && CONTROL %in% levels(sub_cd$condition))
    sub_cd$condition <- relevel(sub_cd$condition, ref = CONTROL)
  sub_ct <- cts[, rownames(sub_cd), drop = FALSE]
  stopifnot(identical(colnames(sub_ct), rownames(sub_cd)))
  mm <- model.matrix(formula, sub_cd)
  rdf <- nrow(sub_cd) - ncol(mm)

  cat("\n── ", label, " ──\n", sep = "")
  cat("  samples      : ", nrow(sub_cd), "   ",
      paste(names(table(sub_cd$condition)), table(sub_cd$condition), sep = "=",
            collapse = ", "), "\n", sep = "")
  cat("  design       : ", deparse(formula), "\n", sep = "")
  cat("  coefficients : ", ncol(mm), " -> residual df = ", rdf, "\n", sep = "")
  cat("  constant     : ", paste(constant_cols(sub_cd), collapse = ", "), "\n", sep = "")

  if (qr(mm)$rank < ncol(mm)) {
    check(paste0(label, ": model matrix full rank"), FALSE,
          sprintf("rank %d of %d", qr(mm)$rank, ncol(mm)))
    return(NULL)
  }
  if (rdf < MIN_RESID) {
    cat("  SKIP         : ", rdf, " residual df -- dispersion not estimable\n", sep = "")
    return(NULL)
  }

  dds <- DESeqDataSetFromMatrix(countData = sub_ct, colData = sub_cd, design = formula)
  check(paste0(label, ": model matrix full rank"), TRUE, sprintf("rank %d of %d", qr(mm)$rank, ncol(mm)))
  check(paste0(label, ": counts unchanged by assembly"), identical(counts(dds), sub_ct))
  check(paste0(label, ": colData in count-column order"),
        identical(colnames(counts(dds)), rownames(colData(dds))))
  check(paste0(label, ": reference is the control"), levels(dds$condition)[1] == CONTROL)
  check(paste0(label, ": unfitted (no size factors)"), is.null(sizeFactors(dds)))
  check(paste0(label, ": condition term is LAST"),
        tail(attr(terms(formula), "term.labels"), 1) == "condition")
  cat("  model matrix : ", paste(colnames(mm), collapse = ", "), "\n", sep = "")
  cond_coefs <- grep("^condition", colnames(mm), value = TRUE)
  cat("  condition coefficients to test: ", paste(cond_coefs, collapse = ", "), "\n", sep = "")
  if (length(cond_coefs) > 1)
    cat("  NOTE: >1 condition coefficient -- results() with no arguments would return \"",
        tail(colnames(mm), 1), "\". Name each coefficient explicitly in step 7.\n", sep = "")
  dds
}

# ═══════════════════════════════════════════ A. global object, QC use only
cat("\n════════ A. Global object (QC / normalisation only, NOT for testing) ════════")
dds_all <- build("dds_all", rep(TRUE, nrow(cd)), ~ experiment + condition)
cat("  NOTE: this model forces one shared condition effect across two experiments run on\n")
cat("  different instruments. Used only to give step 5 a single VST across all 16 samples.\n")

# ═══════════════════════════════════════════ B. per-sub-experiment objects
cat("\n════════ B. Per-sub-experiment objects (these get tested) ════════")
dds_list <- list()
for (ex in levels(cd$experiment)) {
  rows <- cd$experiment == ex
  n_plat <- length(unique(cd$platform[rows]))
  # platform enters as the batch term only where it actually varies within the sub-experiment
  f <- if (n_plat > 1) ~ platform + condition else ~ condition
  d <- build(paste0("dds_", ex), rows, f)
  if (!is.null(d)) dds_list[[ex]] <- d
}
check("DENV object built", "DENV" %in% names(dds_list))
check("ZIKV object built", "ZIKV" %in% names(dds_list))
check("DENV design has no batch term (platform constant)",
      identical(deparse(design(dds_list$DENV)), "~condition"))
check("ZIKV design carries platform as the batch term",
      identical(deparse(design(dds_list$ZIKV)), "~platform + condition"))

# ------------------------------------------------------------------------ summary
cat("\n── Summary of objects ──\n")
rows <- rbindlist(lapply(names(dds_list), function(n) {
  d <- dds_list[[n]]; mm <- model.matrix(design(d), colData(d))
  data.table(object = paste0("dds_", n), samples = ncol(d),
             arms = paste(levels(d$condition), collapse = "/"),
             design = paste(deparse(design(d)), collapse = ""),
             coefficients = ncol(mm), residual_df = ncol(d) - ncol(mm),
             genes = nrow(d), status = "will be tested")
}))
mm_all <- model.matrix(design(dds_all), colData(dds_all))
rows <- rbind(data.table(object = "dds_all", samples = ncol(dds_all),
                         arms = paste(levels(dds_all$condition), collapse = "/"),
                         design = paste(deparse(design(dds_all)), collapse = ""),
                         coefficients = ncol(mm_all),
                         residual_df = ncol(dds_all) - ncol(mm_all),
                         genes = nrow(dds_all), status = "QC only"), rows)
print(rows, row.names = FALSE)

cat("\n── Contrasts these objects support ──\n")
cat("  dds_DENV : DENV  vs Mock (DENV-Mock)      3 v 3, all MiSeq, platform-matched\n")
cat("  dds_ZIKV : ZIKVC vs Mock                  2 v 4, Asian strain\n")
cat("  dds_ZIKV : ZIKVM vs Mock                  4 v 4, African strain\n")
cat("  All three hold tissue constant (hNPC), which is what makes a later DENV-vs-ZIKV\n")
cat("  comparison in this series meaningful.\n")

# ------------------------------------------------------------------------ output
saveRDS(dds_all, file.path(out_dir, "dds_all_unfitted.rds"))
for (n in names(dds_list))
  saveRDS(dds_list[[n]], file.path(out_dir, paste0("dds_", n, "_unfitted.rds")))
fwrite(rows, file.path(out_dir, "design_summary.csv"))

cat("\n── Output ──\n")
for (f in sort(list.files(out_dir))) cat("  ", f, "\n", sep = "")

cat("\n═══════════ RESULT ═══════════\n")
if (length(fails)) {
  cat("  ", length(fails), " CHECK(S) FAILED -- do not proceed:\n", sep = "")
  for (f in fails) cat("    - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("  All checks passed.\n")
cat("  ", length(dds_list), " testable objects: dds_DENV (~condition) and dds_ZIKV",
    " (~platform + condition), baseline Mock.\n", sep = "")
cat("  Plus dds_all for QC.\n")
