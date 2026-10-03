#!/usr/bin/env Rscript
# STEP 2 -- Metadata (colData) and sample matching for GSE80434.
#
#   Rscript scripts/GSE80434/step2_metadata.R
#
# Builds colData from the GEO metadata, proves the sample names line up with the count matrix
# columns exactly, and sets the control level explicitly as the reference.
#
# The hard part of this series is that it contains TWO SUB-EXPERIMENTS with SEPARATE CONTROLS:
#
#   experiment = DENV : DENV-1/2/3        vs DENV-Mock-1/2/3   -- all on MiSeq
#   experiment = ZIKV : ZIKVC-1/2,
#                       ZIKVM-1/2/3/4     vs Mock-1/2/3/4      -- spread over 3 platforms
#
# GEO records plain "mock" for all seven control samples, so the split is only recoverable
# from the sample titles. Pooling the two mock sets would push MiSeq/NextSeq/HiSeq batch
# variation straight into the otherwise clean DENV contrast, so they are kept apart via an
# `experiment` column. That column then routes each mock to its own sub-experiment, which lets
# `condition` stay a single clean factor (Mock / DENV / ZIKVC / ZIKVM) with one reference.

suppressPackageStartupMessages(library(data.table))

GSE       <- "GSE80434"
CONTROL   <- "Mock"
counts_in <- file.path("results", GSE, "step1_raw_counts", "raw_counts_matrix.tsv")
meta_in   <- file.path("data", GSE, "metadata_full.tsv")
out_dir   <- file.path("results", GSE, "step2_metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fails <- character(0)
check <- function(label, ok, detail = "") {
  cat(sprintf("  %-58s %s%s\n", label, if (ok) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) fails <<- c(fails, label)
}

cat("═══════════ STEP 2: METADATA & SAMPLE MATCHING -- ", GSE, " ═══════════\n\n")

if (!file.exists(counts_in)) stop("step 1 output missing -- run step1_raw_counts.R first")
dt   <- fread(counts_in)
cts  <- as.matrix(dt[, -1, with = FALSE]); rownames(cts) <- as.character(dt[[1]])
meta <- fread(meta_in)
m <- meta[match(colnames(cts), meta$geo_accession)]

# ------------------------------------------------------------------ build colData
cat("── Factor derivation ──\n")
cat("  condition  <- ch_virus            'mock' -> Mock; DENV / ZIKVC / ZIKVM as given\n")
cat("  tissue     <- source_name_ch1     'hNPC' (constant)\n")
cat("  platform   <- instrument_model    MiSeq / HiSeq 2000 / NextSeq 500\n")
cat("  experiment <- TITLE               which sub-experiment, i.e. which mock set applies\n")
cat("                                    (GEO says only 'mock' -- titles are the sole source)\n\n")

cond <- ifelse(m$ch_virus == "mock", "Mock", m$ch_virus)
# The routing rule: a DENV-Mock title belongs to the DENV sub-experiment, everything DENV
# does too; all ZIKV arms and the plain Mock titles belong to the ZIKV sub-experiment.
experiment <- ifelse(grepl("DENV", m$title), "DENV", "ZIKV")
mock_set <- ifelse(cond != "Mock", NA_character_,
                   ifelse(grepl("DENV-Mock", m$title), "DENV-Mock", "Mock"))

# DESeq2 warns that design factor levels should stick to letters, numbers, '_' and '.',
# because level names become coefficient names. "Illumina MiSeq" would give the coefficient
# `platformIllumina MiSeq`, which is awkward to pass to results(name=) and easy to mistype.
# Sanitise the levels and keep the verbatim instrument string in its own column.
plat_map <- c("Illumina MiSeq" = "MiSeq", "Illumina HiSeq 2000" = "HiSeq2000",
              "NextSeq 500" = "NextSeq500")
plat <- unname(plat_map[m$instrument_model])

cd <- data.table(
  sample_id  = as.character(m$geo_accession),
  condition  = cond,
  experiment = experiment,
  platform   = plat,
  instrument = m$instrument_model,
  tissue     = "hNPC",
  mock_set   = mock_set,
  replicate  = sub("^.*-([0-9]+)$", "\\1", m$title),
  title      = m$title)

check("all GEO virus values mapped", all(cond %in% c("Mock", "DENV", "ZIKVC", "ZIKVM")),
      paste(unique(cond), collapse = ", "))
check("all instrument names mapped to safe platform labels", !anyNA(plat),
      paste(sort(unique(plat)), collapse = ", "))
check("platform labels are R-safe (no spaces or punctuation)",
      all(grepl("^[A-Za-z0-9_.]+$", unique(plat))))
check("every sample assigned to exactly one sub-experiment",
      all(experiment %in% c("DENV", "ZIKV")) && !anyNA(experiment))
check("DENV sub-experiment holds only DENV and its own mocks",
      setequal(unique(cd$condition[cd$experiment == "DENV"]), c("DENV", "Mock")))
check("ZIKV sub-experiment holds ZIKVC, ZIKVM and its own mocks",
      setequal(unique(cd$condition[cd$experiment == "ZIKV"]), c("ZIKVC", "ZIKVM", "Mock")))
check("the 7 mocks split 3 / 4 between sub-experiments",
      identical(as.integer(table(cd$experiment[cd$condition == "Mock"])), c(3L, 4L)),
      paste(names(table(cd$experiment[cd$condition == "Mock"])),
            table(cd$experiment[cd$condition == "Mock"]), sep = "=", collapse = ", "))

# ------------------------------------------------- sample matching (critical part)
cat("\n── Sample matching ──\n")
check("same number of samples", ncol(cts) == nrow(cd),
      paste0(ncol(cts), " columns vs ", nrow(cd), " rows"))
check("same set of sample names (order aside)", setequal(colnames(cts), cd$sample_id))
check("no count column absent from colData", length(setdiff(colnames(cts), cd$sample_id)) == 0)
check("no colData row absent from count matrix", length(setdiff(cd$sample_id, colnames(cts))) == 0)

# Reorder colData to the counts, never the reverse: the counts define the truth, and aligning
# by name is what makes the identity test below meaningful rather than circular.
setDF(cd); rownames(cd) <- cd$sample_id
cd <- cd[colnames(cts), , drop = FALSE]

check("all(colnames(counts) == rownames(metadata))", all(colnames(cts) == rownames(cd)))
check("identical(colnames(counts), rownames(metadata))", identical(colnames(cts), rownames(cd)))
cat("\n  verification (all ", nrow(cd), " rows):\n", sep = "")
print(data.frame(count_column = colnames(cts), metadata_rowname = rownames(cd),
                 match = colnames(cts) == rownames(cd), row.names = NULL))

# ---------------------------------------------------------- factors + reference level
cat("\n── Factor levels and reference ──\n")
cd$condition <- factor(cd$condition)
cat("  before relevel : ", paste(levels(cd$condition), collapse = ", "),
    "  -> reference \"", levels(cd$condition)[1], "\"\n", sep = "")
# Alphabetically DENV precedes Mock, so the default baseline here would be DENV and every
# fold change would come out sign-flipped. relevel() is doing real work in this dataset.
cd$condition <- relevel(cd$condition, ref = CONTROL)
cat("  after relevel  : ", paste(levels(cd$condition), collapse = ", "),
    "  -> reference \"", levels(cd$condition)[1], "\"\n", sep = "")
cat("  (alphabetical order would have made DENV the baseline -- relevel matters here)\n")
check(paste0("reference level is \"", CONTROL, "\""), levels(cd$condition)[1] == CONTROL)
check("4 condition levels", nlevels(cd$condition) == 4,
      paste(levels(cd$condition), collapse = ", "))

cd$experiment <- factor(cd$experiment, levels = c("DENV", "ZIKV"))
cd$platform   <- factor(cd$platform)
cd$tissue     <- factor(cd$tissue)
cd$replicate  <- factor(cd$replicate)

cat("\n── colData ──\n")
print(cd[, c("condition", "experiment", "platform", "tissue", "mock_set")])

cat("\n── Design tables ──\n")
cat("  condition x experiment:\n"); print(table(condition = cd$condition, experiment = cd$experiment))
cat("\n  condition x platform:\n");  print(table(condition = cd$condition, platform = cd$platform))

# --------------------------------------------------------------- design feasibility
cat("\n── Replication per contrast ──\n")
feas <- as.data.table(cd)[, .(n = .N), by = .(experiment, condition)]
print(dcast(feas, experiment ~ condition, value.var = "n", fill = 0), row.names = FALSE)
for (ex in levels(cd$experiment)) {
  sub <- cd[cd$experiment == ex, , drop = FALSE]
  n_mock <- sum(sub$condition == "Mock")
  for (arm in setdiff(unique(as.character(sub$condition)), "Mock")) {
    n_arm <- sum(sub$condition == arm)
    check(paste0(ex, ": ", arm, " vs Mock has >=2 per arm"), n_arm >= 2 && n_mock >= 2,
          paste0(n_arm, " vs ", n_mock))
  }
}

# ------------------------------------------------- can platform be modelled per experiment?
cat("\n── Platform: confounded or estimable, per sub-experiment? ──\n")
for (ex in levels(cd$experiment)) {
  sub <- droplevels(cd[cd$experiment == ex, , drop = FALSE])
  npl <- nlevels(sub$platform)
  if (npl == 1) {
    cat("  ", ex, ": single platform (", levels(sub$platform),
        ") -> constant, not a model term. Contrast is platform-matched.\n", sep = "")
    check(paste0(ex, ": platform constant, no confounding possible"), TRUE)
    next
  }
  mm_c  <- model.matrix(~ condition, sub)
  mm_pc <- model.matrix(~ platform + condition, sub)
  full  <- qr(mm_pc)$rank == ncol(mm_pc)
  cat("  ", ex, ": ", npl, " platforms | ~condition df=", nrow(sub) - ncol(mm_c),
      " | ~platform+condition df=", nrow(sub) - ncol(mm_pc),
      " | full rank: ", full, "\n", sep = "")
  print(table(condition = sub$condition, platform = sub$platform))
  check(paste0(ex, ": platform is estimable alongside condition (not aliased)"), full,
        sprintf("rank %d of %d", qr(mm_pc)$rank, ncol(mm_pc)))
  cat("  -> platform can be carried as a batch term; step 6 should test whether it earns\n")
  cat("     its degrees of freedom rather than assuming so.\n")
}

cat("\n── Factors that cannot be model terms ──\n")
for (v in c("tissue", "condition", "experiment", "platform", "replicate")) {
  nl <- nlevels(cd[[v]])
  cat(sprintf("  %-11s %d level(s)%s\n", v, nl, if (nl == 1) "  <- constant" else ""))
}
cat("  replicate numbering restarts within each arm, so it encodes no shared batch and is\n")
cat("  kept for reference only.\n")

# ------------------------------------------------------------------------ write out
fwrite(data.table(sample_id = rownames(cd),
                  cd[, c("condition", "experiment", "platform", "instrument", "tissue",
                         "mock_set", "replicate", "title")]),
       file.path(out_dir, "colData.csv"))
saveRDS(cd, file.path(out_dir, "colData.rds"))
fwrite(dcast(feas, experiment ~ condition, value.var = "n", fill = 0),
       file.path(out_dir, "design_feasibility.csv"))

cat("\n── Output ──\n")
cat("  ", file.path(out_dir, "colData.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "colData.rds"), "   (factor levels and reference preserved)\n", sep = "")
cat("  ", file.path(out_dir, "design_feasibility.csv"), "\n", sep = "")

cat("\n═══════════ RESULT ═══════════\n")
if (length(fails)) {
  cat("  ", length(fails), " CHECK(S) FAILED -- do not proceed:\n", sep = "")
  for (f in fails) cat("    - ", f, "\n", sep = "")
  quit(status = 1L)
}
cat("  All checks passed.\n")
cat("  colData aligned to count columns by name; reference level = \"", CONTROL, "\".\n", sep = "")
cat("  Positive log2FC will mean UP in the virus arm relative to its OWN mock set.\n")
