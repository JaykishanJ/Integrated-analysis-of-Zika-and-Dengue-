# ==============================================================================
# PHASE A - INPUT & PREPROCESSING 
# ==============================================================================
# GSE110496 viscRNA-seq: Count matrix assembly, feature routing, Ensembl mapping.
#
# Outputs:
# results/phaseA/objects/Phase_A_counts_host_ercc.rds -> Phase B
# results/phaseA/objects/Phase_A_cell_qc.rds -> Phase B
# results/phaseA/objects/Phase_A_gene_map.rds -> Phase B/C/D
# results/phaseA/cache/Phase_A_counts_raw.rds -> re-run cache
#
# NOTE: all figures are written as vector SVG only (results/phaseA/figures/*.svg)
# ==============================================================================

library(vroom)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
  source("E:/Zika/scripts/inject_themes.R")
library(httr)
library(tools)

# svglite gives a much more faithful/lightweight SVG device than the base
# grDevices::svg() backend, so we use it for every plot in this script.
if (!requireNamespace("svglite", quietly = TRUE))
  stop("svglite is not installed.\n",
       " install.packages('svglite')")
library(svglite)

# Guard against S4 generics masking dplyr verbs
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  stop("org.Hs.eg.db is not installed.\n",
       " if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
       " BiocManager::install('org.Hs.eg.db')")

.verbs <- c("select", "rename", "filter", "mutate", "slice", "desc", "first")
.owner <- vapply(.verbs, function(f) find(f)[1], character(1))
.masked <- .owner[!.owner %in% c("package:dplyr", "package:tidyr")]
if (length(.masked))
  stop("dplyr verbs are masked by another package: ",
       paste(sprintf("%s() <- %s", names(.masked), .masked), collapse = ", "),
       "\nRestart R and source this script in a clean session.")
rm(.verbs, .owner, .masked)

# ==============================================================================
# Configuration
# ==============================================================================
data_dir <- "E:/Zika"
raw_dir <- file.path(data_dir, "GSE110496_RAW")
matrix_file <- file.path(data_dir, "GSE110496_series_matrix.txt")

save_count_matrix <- TRUE
run_ensembl_check <- TRUE
use_cache <- TRUE
migrate_legacy <- TRUE
viral_bg_quantile <- 1.00 # 1.00 = highest MOI-0 count (zero FPs); 0.99 for higher sensitivity

setwd(data_dir)

# ==============================================================================
# Output layout
# ==============================================================================
results_root <- file.path(data_dir, "results", "phaseA")
dir_fig <- file.path(results_root, "figures")
dir_tab <- file.path(results_root, "tables")
dir_obj <- file.path(results_root, "objects")
dir_cch <- file.path(results_root, "cache")
for (d in c(dir_fig, dir_tab, dir_obj, dir_cch))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

fig <- function(name) file.path(dir_fig, name)
tab <- function(name) file.path(dir_tab, name)
obj <- function(name) file.path(dir_obj, name)

raw_cache <- file.path(dir_cch, "Phase_A_counts_raw.rds")

# One-time tidy-up for legacy flat-layout outputs
if (migrate_legacy) {
  legacy <- list(
    figures = c("Verification_Files_Found", "Feature_Class_Composition",
                "Verification_Viral_Scatter", "Viral_Load_By_Design",
                "Infected_Fraction", "Verification_Bucket_Counts",
                "Verification_Library_Complexity", "Verification_Plate_Position",
                "Verification_Highest_Expressed", "Verification_Alignment_Outcome",
                "Verification_Mapping_Outcome", "Verification_Unmapped_Expr"),
    tables = c("Sample_Metadata.csv", "Batch_File_Check.csv",
               "MD5_Subsample.csv", "Cell_QC_Table.csv",
               "Viral_qPCR_Concordance.csv", "Viral_Load_Summary.csv",
               "Bucket_Counts.csv", "Ensembl_Spotcheck.csv",
               "Gene_Mapping_Table.csv"),
    objects = c("Phase_A_counts_host_ercc.rds", "Phase_A_cell_qc.rds",
                "Phase_A_gene_map.rds"),
    cache = c("Phase_A_counts_raw.rds")
  )
  moved <- 0L
  move_one <- function(from, to) {
    if (file.exists(from) && !file.exists(to)) {
      if (file.rename(from, to)) return(1L)
    }
    0L
  }
  for (stem in legacy$figures)
    for (ext in c(".pdf", ".png", ".svg"))
      moved <- moved + move_one(file.path(data_dir, paste0(stem, ext)),
                                file.path(dir_fig, paste0(stem, ext)))
  for (f in legacy$tables)
    moved <- moved + move_one(file.path(data_dir, f), file.path(dir_tab, f))
  for (f in legacy$objects)
    moved <- moved + move_one(file.path(data_dir, f), file.path(dir_obj, f))
  for (f in legacy$cache)
    moved <- moved + move_one(file.path(data_dir, f), file.path(dir_cch, f))
  if (moved > 0L)
    cat(sprintf(" -> Migrated %d legacy output file(s) into results/phaseA/.\n", moved))
}

# --------------------------------------------------------------------------
# show_and_save(): Apply custom theme settings, print, and save as SVG.
# The global theme is set by inject_themes.R -> theme_pub().
# This function applies any plot-specific overrides before saving.
# --------------------------------------------------------------------------
show_and_save <- function(p, file, width = 7, height = 5, base_size = 10) {
  print(p)
  ggsave(fig(paste0(file, ".svg")), p, width = width, height = height,
         bg = "transparent", device = svglite::svglite)
  ggsave(fig(paste0(file, ".png")), p, width = width, height = height,
         dpi = 600, bg = "white")
  cat(sprintf(" [plot] %s.svg and .png\n", file))
  invisible(p)
}

min_pos <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x)) min(x) else 1
}

cat("Starting Phase A Pre-processing...\n")
cat(sprintf("Output directory: %s\n\n", getwd()))

# ==============================================================================
# Ingest & Merge
# ==============================================================================
cat("Ingest & Merge\n")

## 0a. Parse metadata from GEO series matrix
cat(" -> Parsing metadata from GEO series matrix...\n")
lines <- readLines(matrix_file, warn = FALSE)

split_row <- function(line) {
  v <- str_split(line, "\t")[[1]][-1]
  str_remove_all(v, '^"|"$')
}

acc_line <- lines[str_detect(lines, "^!Sample_geo_accession")]
if (length(acc_line) != 1L) stop("Could not find a unique !Sample_geo_accession row.")
gsm_ids <- split_row(acc_line[1])

title_line <- lines[str_detect(lines, "^!Sample_title")]
sample_titles <- if (length(title_line)) split_row(title_line[1]) else rep(NA_character_, length(gsm_ids))

metadata <- data.frame(
  Sample = gsm_ids,
  title = sample_titles,
  stringsAsFactors = FALSE
)

# Parse characteristics rows (key can differ per sample within one row)
char_lines <- lines[str_detect(lines, "^!Sample_characteristics_ch1")]
cat(sprintf(" %d characteristics rows, %d samples\n", length(char_lines), length(gsm_ids)))

char_long <- lapply(seq_along(char_lines), function(i) {
  vals <- split_row(char_lines[i])
  parts <- str_split_fixed(vals, ": ", 2)
  data.frame(
    Sample = gsm_ids,
    row_i = i,
    key = str_trim(parts[, 1]),
    value = str_trim(parts[, 2]),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  filter(key != "", !is.na(key))

keys_per_row <- char_long %>%
  group_by(row_i) %>%
  summarise(keys = paste(sort(unique(key)), collapse = " | "),
            n = n(), .groups = "drop")
cat(" Keys found per characteristics row:\n")
for (i in seq_len(nrow(keys_per_row)))
  cat(sprintf(" row %d (n=%4d): %s\n", keys_per_row$row_i[i],
              keys_per_row$n[i], keys_per_row$keys[i]))

char_wide <- char_long %>%
  mutate(col = str_replace_all(str_replace_all(tolower(key), "[^a-z0-9]+", "_"),
                               "^_|_$", "")) %>%
  select(Sample, col, value) %>%
  distinct(Sample, col, .keep_all = TRUE) %>%
  pivot_wider(names_from = col, values_from = value)

metadata <- left_join(metadata, char_wide, by = "Sample")

# Unify viral-load keys
num_or_na <- function(x) suppressWarnings(as.numeric(x))
get_col <- function(df, nm)
  if (nm %in% names(df)) df[[nm]] else rep(NA_character_, nrow(df))

metadata$n_virus_molecules <- num_or_na(dplyr::coalesce(
  get_col(metadata, "n_dengue_molecules"),
  get_col(metadata, "n_zika_molecules")
))
metadata$time_h <- num_or_na(get_col(metadata, "time_h"))
metadata$moi <- num_or_na(get_col(metadata, "moi"))
if ("virus_qpcr" %in% names(metadata))
  metadata$virus_qpcr <- num_or_na(metadata$virus_qpcr)

cat(sprintf(" Metadata: %d samples x %d fields\n", nrow(metadata), ncol(metadata)))
cat(sprintf(" virus: %s\n",
            paste(sprintf("%s=%d", names(table(metadata$virus)),
                          as.integer(table(metadata$virus))), collapse = ", ")))
cat(sprintf(" n_virus_molecules: min=%.0f median=%.0f max=%.0f (NA=%d)\n",
            min(metadata$n_virus_molecules, na.rm = TRUE),
            median(metadata$n_virus_molecules, na.rm = TRUE),
            max(metadata$n_virus_molecules, na.rm = TRUE),
            sum(is.na(metadata$n_virus_molecules))))

write.csv(metadata, tab("Sample_Metadata.csv"), row.names = FALSE)

## 0b. Locate GSM files
gsm_files <- list.files(raw_dir, pattern = "_counts\\.tsv\\.gz$", full.names = TRUE)
cat(sprintf(" -> Found %d GSM count files.\n", length(gsm_files)))

file_gsm <- str_match(basename(gsm_files), "^(GSM\\d+)_")[, 2]
file_batch <- str_match(basename(gsm_files), "^GSM\\d+_(\\d+)_")[, 2]
file_well <- str_match(basename(gsm_files), "^GSM\\d+_\\d+_([A-Z]\\d+)_")[, 2]
file_row <- str_match(file_well, "^([A-Z])")[, 2]
file_col <- as.integer(str_match(file_well, "^[A-Z](\\d+)$")[, 2])
stopifnot(!any(is.na(file_gsm)), !any(is.na(file_batch)), !any(is.na(file_well)))
cat(sprintf(" Plate layout: rows %s-%s, columns %d-%d (%d distinct wells used).\n",
            min(file_row), max(file_row), min(file_col), max(file_col),
            length(unique(file_well))))

# Plate = replicate unit. Check â‰¥3 plates per condition.
plate_counts <- metadata %>%
  filter(!is.na(virus), !is.na(moi)) %>%
  distinct(batch = file_batch[match(title, metadata$title)], virus, moi) %>%
  filter(!is.na(batch)) %>%
  count(virus, moi, name = "n_plates")
cat(" Plates per condition (virus x MOI):\n")
print(as.data.frame(plate_counts))
insufficient <- plate_counts %>% filter(n_plates < 3)
if (nrow(insufficient) > 0) {
  cat(sprintf(" WARNING: %d condition(s) have < 3 plates. Pseudobulk () will be skipped for these.\n",
              nrow(insufficient)))
  write.csv(insufficient, tab("Insufficient_Plates.csv"), row.names = FALSE)
}

meta_batch <- str_match(metadata$title, "^(\\d+)_")[, 2]

batch_check <- full_join(
  as.data.frame(table(batch = file_batch), stringsAsFactors = FALSE) %>%
    rename(Found = Freq),
  as.data.frame(table(batch = meta_batch), stringsAsFactors = FALSE) %>%
    rename(Expected = Freq),
  by = "batch"
) %>%
  mutate(across(c(Found, Expected), ~ replace_na(as.integer(.), 0L)),
         Missing = Expected - Found)

cat(sprintf(" Batches: %d | total expected %d | total found %d | mismatched: %d\n",
            nrow(batch_check), sum(batch_check$Expected), sum(batch_check$Found),
            sum(batch_check$Missing != 0)))
write.csv(batch_check, tab("Batch_File_Check.csv"), row.names = FALSE)

p1 <- batch_check %>%
  pivot_longer(c(Found, Expected), names_to = "Source", values_to = "n") %>%
  ggplot(aes(x = batch, y = n, fill = Source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black") +
  scale_fill_manual(values = c(Expected = "#888888", Found = "#0072B2")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Files found vs expected per batch",
       subtitle = "Expected = samples per plate in the series matrix; Found = .tsv.gz files on disk",
       x = "Batch (plate) ID", y = "Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
show_and_save(p1, "Verification_Files_Found", width = 9, height = 5)

## MD5 checksums on a 5% subsample
cat(" -> Computing MD5 checksums for a 5% subsample...\n")
set.seed(42)
sub_idx <- sample(seq_along(gsm_files), max(1, floor(length(gsm_files) * 0.05)))
md5_sums <- md5sum(gsm_files[sub_idx])
md5_df <- data.frame(file = basename(names(md5_sums)), md5 = unname(md5_sums),
                     stringsAsFactors = FALSE)
write.csv(md5_df, tab("MD5_Subsample.csv"), row.names = FALSE)
cat(sprintf(" %d checksums computed | %d unique | %d NA (unreadable)\n",
            nrow(md5_df), length(unique(md5_df$md5)), sum(is.na(md5_df$md5))))

## 0c. Load and merge
read_counts <- function(f) {
  vroom(f, delim = "\t", col_names = c("feature", "count"), skip = 1,
        col_types = cols(feature = col_character(), count = col_integer()),
        altrep = FALSE, progress = FALSE)
}

cache_ok <- FALSE
if (use_cache && file.exists(raw_cache)) {
  cat(sprintf(" -> Reusing cached count matrix from %s ...\n", raw_cache))
  count_mat <- readRDS(raw_cache)
  idx <- match(colnames(count_mat), file_gsm)
  if (!anyNA(idx) && nrow(count_mat) > 0) {
    file_gsm <- file_gsm[idx]
    file_batch <- file_batch[idx]
    file_well <- file_well[idx]
    file_row <- file_row[idx]
    file_col <- file_col[idx]
    genes <- rownames(count_mat)
    n_gene <- length(genes)
    cache_ok <- TRUE
    cat(sprintf(" Cache accepted: %d features x %d cells.\n",
                nrow(count_mat), ncol(count_mat)))
  } else {
    cat(" Cache does not match the files on disk; re-reading.\n")
    rm(count_mat)
  }
}

if (!cache_ok) {
  cat(" -> Reading files (header row is skipped; counts read as integer)...\n")
  first <- read_counts(gsm_files[1])
  genes <- first$feature
  n_gene <- length(genes)
  cat(sprintf(" %d features per file; %d unique.\n", n_gene, length(unique(genes))))
  stopifnot(!anyDuplicated(genes))
  
  count_mat <- matrix(NA_integer_, nrow = n_gene, ncol = length(gsm_files),
                      dimnames = list(genes, file_gsm))
  count_mat[, 1] <- first$count
  
  ok <- rep(TRUE, length(gsm_files))
  for (j in seq_along(gsm_files)) {
    if (j == 1) next
    d <- tryCatch(read_counts(gsm_files[j]), error = function(e) NULL)
    if (is.null(d) || nrow(d) != n_gene || !identical(d$feature, genes)) {
      ok[j] <- FALSE
      next
    }
    count_mat[, j] <- d$count
    if (j %% 250 == 0) cat(sprintf(" ... %d / %d\n", j, length(gsm_files)))
  }
  rm(list = intersect(c("first", "d"), ls())); invisible(gc())
  
  if (any(!ok)) {
    cat(sprintf(" WARNING: dropping %d file(s) that were unreadable or mis-ordered:\n",
                sum(!ok)))
    cat(paste0(" ", head(basename(gsm_files[!ok]), 10), collapse = "\n"), "\n")
    count_mat <- count_mat[, ok, drop = FALSE]
    file_batch <- file_batch[ok]
    file_gsm <- file_gsm[ok]
    file_well <- file_well[ok]
    file_row <- file_row[ok]
    file_col <- file_col[ok]
  }
  
  if (use_cache) {
    cat(sprintf(" Caching the assembled matrix to %s ...\n", raw_cache))
    saveRDS(count_mat, raw_cache)
  }
}
cat(sprintf(" Count matrix: %d features x %d cells.\n", nrow(count_mat), ncol(count_mat)))

## Row-name token composition
gene_tokens <- data.frame(gene = genes, stringsAsFactors = FALSE) %>%
  mutate(type = case_when(
    str_detect(gene, "^ENSG\\d+") ~ "Host (ENSG)",
    str_detect(gene, "^ERCC-\\d+") ~ "ERCC spike-in",
    str_detect(gene, "^(DENV|ZIKV)") ~ "Viral",
    str_detect(gene, "^__") ~ "HTSeq QC row",
    gene == "NIST_ConsensusVector" ~ "NIST vector",
    TRUE ~ "Other"
  ))
tok_tab <- as.data.frame(table(type = gene_tokens$type), stringsAsFactors = FALSE)
cat(" Feature-name tokens:\n")
for (i in seq_len(nrow(tok_tab)))
  cat(sprintf(" %-16s %6d\n", tok_tab$type[i], tok_tab$Freq[i]))

p2 <- ggplot(tok_tab, aes(x = reorder(type, -Freq), y = Freq, fill = type)) +
  geom_col(colour = "black", linewidth = 0.3) +
  geom_text(aes(label = Freq), vjust = -0.4, size = 3.2) +
  scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = feature_colors, name = NULL) +
  labs(title = "Row-name token distribution",
       subtitle = sprintf("%d features per file, identical order in all %d files",
                          n_gene, ncol(count_mat)),
       x = "Feature type", y = "Number of rows (log10)")
show_and_save(p2, "Feature_Class_Composition", width = 7.5, height = 5)

# ==============================================================================
# Route Features + Infection Calling + Viral load
# ==============================================================================
cat("\nRoute Features, Viral Load, and Infection Calling\n")

## 1a. Regex-split
idx_host <- str_detect(genes, "^ENSG\\d+")
idx_ercc <- str_detect(genes, "^ERCC-\\d+")
idx_viral <- str_detect(genes, "^(DENV|ZIKV)")
idx_qc <- str_detect(genes, "^__")
idx_other <- !(idx_host | idx_ercc | idx_viral | idx_qc)

cat(sprintf(" host=%d ERCC=%d viral=%d HTSeq-QC=%d other=%d\n",
            sum(idx_host), sum(idx_ercc), sum(idx_viral), sum(idx_qc), sum(idx_other)))

## 1b. Per-cell bucket sums
chunked_colsum <- function(idx, fun = identity) {
  out <- numeric(ncol(count_mat))
  if (!any(idx)) return(out)
  rows <- which(idx)
  for (chunk in split(rows, ceiling(seq_along(rows) / 5000)))
    out <- out + colSums(fun(count_mat[chunk, , drop = FALSE]))
  out
}
bucket_sum <- function(idx) chunked_colsum(idx)
host_sum <- bucket_sum(idx_host)
ercc_sum <- bucket_sum(idx_ercc)
viral_sum <- bucket_sum(idx_viral)
qc_sum <- bucket_sum(idx_qc)
other_sum <- bucket_sum(idx_other)

htseq_mat <- count_mat[idx_qc, , drop = FALSE]
htseq_by_cat <- as.data.frame(t(htseq_mat))
colnames(htseq_by_cat) <- str_remove(rownames(htseq_mat), "^__")
htseq_by_cat$Sample <- colnames(count_mat)

cell_qc <- data.frame(
  Sample = colnames(count_mat), batch = file_batch,
  well = file_well, plate_row = file_row, plate_col = file_col,
  host_counts = host_sum, ercc_counts = ercc_sum, viral_counts = viral_sum,
  htseq_qc_counts = qc_sum, other_counts = other_sum,
  stringsAsFactors = FALSE
) %>%
  mutate(n_genes_detected = chunked_colsum(idx_host, function(m) m > 0),
         frac_assigned = (host_counts + ercc_counts + other_counts) /
           (host_counts + ercc_counts + other_counts + htseq_qc_counts)) %>%
  left_join(metadata, by = "Sample")
write.csv(cell_qc, tab("Cell_QC_Table.csv"), row.names = FALSE)

## 1c. Viral load validated against orthogonal qPCR
if (sum(idx_viral) > 0 && sum(viral_sum) > 0)
  stop("Viral feature rows appeared in the count files;  needs revisiting.")

vdf <- cell_qc %>% filter(!is.na(n_virus_molecules), !is.na(virus_qpcr))
cat(sprintf(" -> %d cells carry both measurements (%d qPCR>0, %d molecules>0).\n",
            nrow(vdf), sum(vdf$virus_qpcr > 0), sum(vdf$n_virus_molecules > 0)))

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}
cor_by_virus <- vdf %>%
  filter(virus_qpcr > 0, n_virus_molecules > 0) %>%
  group_by(virus) %>%
  summarise(
    n = dplyr::n(),
    pearson = safe_cor(log10(virus_qpcr), log10(n_virus_molecules)),
    spearman = safe_cor(virus_qpcr, n_virus_molecules, method = "spearman"),
    .groups = "drop"
  )
for (i in seq_len(nrow(cor_by_virus)))
  cat(sprintf(" %-7s n=%4d log-Pearson r=%.3f Spearman rho=%.3f\n",
              cor_by_virus$virus[i], cor_by_virus$n[i],
              cor_by_virus$pearson[i], cor_by_virus$spearman[i]))
write.csv(cor_by_virus, tab("Viral_qPCR_Concordance.csv"), row.names = FALSE)

cor_lab <- cor_by_virus %>%
  mutate(lab = sprintf("n = %d\nlog-Pearson r = %.3f\nSpearman = %.3f",
                       n, pearson, spearman))

p3 <- ggplot(vdf, aes(x = virus_qpcr + min_pos(vdf$virus_qpcr),
                      y = n_virus_molecules + 1)) +
  geom_point(aes(colour = factor(moi)), alpha = 0.55, size = 1.4) +
  geom_text(data = cor_lab, aes(x = Inf, y = 1, label = lab),
            hjust = 1.05, vjust = -0.2, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ virus, scales = "free_x") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_brewer(palette = "Set1", name = "MOI") +
  labs(title = "Viral molecules vs the independent per-cell qPCR assay",
       subtitle = sprintf(
         "Two separate assays on the same cells; %d of %d cells carry both. Axes offset by the smallest positive value to keep zeros visible.",
         nrow(vdf), nrow(cell_qc)),
       x = "virus-qpcr (offset, log10)", y = "viral molecules + 1 (log10)")
show_and_save(p3, "Verification_Viral_Scatter", width = 9, height = 5.5)

## 1c-bis. Infection calling from empirical background
cat(" -> Calling infection status from MOI 0 controls...\n")
bg <- cell_qc %>%
  filter(moi == 0, !is.na(n_virus_molecules)) %>%
  group_by(virus) %>%
  summarise(n_control = n(),
            control_pos = sum(n_virus_molecules > 0),
            bg_threshold = quantile(n_virus_molecules, viral_bg_quantile, names = FALSE),
            .groups = "drop")

if (nrow(bg) == 0) {
  cat(" WARNING: no MOI 0 control cells found; falling back to threshold = 0.\n")
  bg <- cell_qc %>% distinct(virus) %>%
    mutate(n_control = 0L, control_pos = 0L, bg_threshold = 0)
}

cat(" Empirical background from MOI 0 controls:\n")
for (i in seq_len(nrow(bg)))
  cat(sprintf(" %-7s %4d controls, %3d with any viral molecule, threshold > %g\n",
              bg$virus[i], bg$n_control[i], bg$control_pos[i], bg$bg_threshold[i]))

cell_qc <- cell_qc %>%
  left_join(bg %>% select(virus, bg_threshold), by = "virus") %>%
  mutate(
    infected = case_when(
      is.na(n_virus_molecules) ~ NA,
      n_virus_molecules > bg_threshold ~ TRUE,
      TRUE ~ FALSE
    ),
    infection_class = case_when(
      is.na(n_virus_molecules) ~ NA_character_,
      moi == 0 ~ "control (MOI 0)",
      infected ~ "infected",
      TRUE ~ "exposed, uninfected"
    ),
    # Viral share of the cell's transcriptome (approximate: numerator from
    # Stampy viral pipeline, denominator from HTSeq host counting)
    viral_frac = n_virus_molecules / (n_virus_molecules + host_counts)
  )

viral_summary <- cell_qc %>%
  filter(!is.na(infection_class)) %>%
  group_by(virus, moi, time_h) %>%
  summarise(n_cells = n(),
            n_infected = sum(infected, na.rm = TRUE),
            pct_infected = 100 * mean(infected, na.rm = TRUE),
            median_virus = median(n_virus_molecules, na.rm = TRUE),
            max_virus = max(n_virus_molecules, na.rm = TRUE),
            .groups = "drop")
write.csv(viral_summary, tab("Viral_Load_Summary.csv"), row.names = FALSE)

cat(sprintf(" Infection calls: %d infected, %d exposed-uninfected, %d control.\n",
            sum(cell_qc$infection_class == "infected", na.rm = TRUE),
            sum(cell_qc$infection_class == "exposed, uninfected", na.rm = TRUE),
            sum(cell_qc$infection_class == "control (MOI 0)", na.rm = TRUE)))

p3b <- cell_qc %>%
  filter(!is.na(n_virus_molecules)) %>%
  ggplot(aes(x = factor(moi), y = n_virus_molecules + 1)) +
  geom_hline(data = bg, aes(yintercept = bg_threshold + 1),
             linetype = "dashed", colour = "red") +
  geom_jitter(aes(colour = infection_class), width = 0.22, alpha = 0.45, size = 1) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black", width = 0.5) +
  facet_grid(virus ~ time_h, labeller = label_both) +
  scale_y_log10() +
  scale_colour_manual(values = c("control (MOI 0)" = "#808080",
                                 "exposed, uninfected" = "#E69F00",
                                 "infected" = "#D55E00"),
                      name = NULL) +
  theme(legend.position = "bottom") +
  labs(title = "Viral load per cell across the experimental design",
       subtitle = "Dashed red line = empirical background threshold from the MOI 0 controls of that virus",
       x = "MOI", y = "viral molecules + 1 (log10)")
show_and_save(p3b, "Viral_Load_By_Design", width = 10, height = 6)

p3c <- viral_summary %>%
  filter(moi > 0) %>%
  ggplot(aes(x = factor(time_h), y = pct_infected, fill = factor(moi))) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "black") +
  geom_text(aes(label = sprintf("%d/%d", n_infected, n_cells)),
            position = position_dodge(width = 0.85), vjust = -0.35, size = 2.7) +
  facet_wrap(~ virus) +
  scale_fill_brewer(palette = "Blues", name = "MOI", direction = -1) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Infected fraction by virus, time and MOI",
       subtitle = "Labels show infected / total cells; MOI 0 controls excluded",
       x = "Hours post-infection", y = "% cells called infected")
show_and_save(p3c, "Infected_Fraction", width = 9, height = 5)

write.csv(cell_qc, tab("Cell_QC_Table.csv"), row.names = FALSE)

## 1d. Global counts per bucket + per-cell composition
bucket_df <- data.frame(
  Bucket = c("Host (ENSG)", "ERCC spike-in", "Viral", "HTSeq QC rows", "Other (NIST)"),
  Total_Counts = c(sum(host_sum), sum(ercc_sum), sum(viral_sum), sum(qc_sum), sum(other_sum)),
  stringsAsFactors = FALSE
) %>%
  mutate(Bucket = factor(Bucket, levels = rev(Bucket)),
         pct = 100 * Total_Counts / sum(Total_Counts))
cat(" Global counts per bucket:\n")
for (i in seq_len(nrow(bucket_df)))
  cat(sprintf(" %-15s %15.0f (%5.2f%%)\n", as.character(bucket_df$Bucket)[i],
              bucket_df$Total_Counts[i], bucket_df$pct[i]))
write.csv(bucket_df, tab("Bucket_Counts.csv"), row.names = FALSE)

bucket_cell <- cell_qc %>%
  transmute(Sample,
            `Host (ENSG)` = host_counts,
            `ERCC spike-in` = ercc_counts,
            `HTSeq QC rows` = htseq_qc_counts,
            `Other (NIST)` = other_counts) %>%
  pivot_longer(-Sample, names_to = "Bucket", values_to = "counts") %>%
  group_by(Sample) %>%
  mutate(pct = 100 * counts / sum(counts)) %>%
  ungroup()

p4 <- ggplot(bucket_cell, aes(x = reorder(Bucket, -pct, FUN = median), y = pct,
                              fill = Bucket)) +
  geom_violin(scale = "width", colour = "#4D4D4D", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white", linewidth = 0.3) +
  theme(legend.position = "none") +
  labs(title = "Per-cell composition across feature buckets",
       subtitle = sprintf("%d cells. Pooled totals: HTSeq QC rows carry %.1f%% of all counts, host %.1f%%.",
                          ncol(count_mat),
                          bucket_df$pct[bucket_df$Bucket == "HTSeq QC rows"],
                          bucket_df$pct[bucket_df$Bucket == "Host (ENSG)"]),
       x = NULL, y = "% of counts in that cell")
show_and_save(p4, "Verification_Bucket_Counts", width = 8, height = 5.5)

## 1d-bis. NEW: Library complexity, plate position, highest expressed, alignment outcome
p4b <- cell_qc %>%
  ggplot(aes(x = host_counts + 1, y = n_genes_detected, colour = frac_assigned)) +
  geom_point(alpha = 0.6, size = 1.1) +
  scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_colour_viridis_c(name = "fraction\nassigned", option = "viridis") +
  facet_wrap(~ virus) +
  labs(title = "Library complexity: host counts vs genes detected",
       subtitle = "One point per cell. Colour is the fraction of counts HTSeq could assign to a feature.",
       x = "Host counts + 1 (log10)", y = "Genes detected")
show_and_save(p4b, "Verification_Library_Complexity", width = 9, height = 5)

plate_df <- cell_qc %>%
  mutate(plate_row = factor(plate_row, levels = rev(sort(unique(plate_row)))))

p4c <- ggplot(plate_df, aes(x = factor(plate_col), y = plate_row,
                            fill = log10(host_counts + 1))) +
  geom_tile(colour = "#F0F0F0", linewidth = 0.15) +
  facet_wrap(~ batch, ncol = 4) +
  scale_fill_viridis_c(name = "log10\nhost counts", option = "magma") +
  theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 5, angle = 90),
        strip.text.y = element_text(size = 6)) +
  labs(title = "Library size by plate position",
       subtitle = "384-well layout per plate. Row or column gradients indicate handling artefacts, not biology.",
       x = "Column", y = "Row")
show_and_save(p4c, "Verification_Plate_Position", width = 13, height = 9)

plate_stats <- cell_qc %>%
  group_by(batch) %>%
  summarise(
    row_kw_p = tryCatch(kruskal.test(host_counts, factor(plate_row))$p.value,
                        error = function(e) NA_real_),
    col_kw_p = tryCatch(kruskal.test(host_counts, factor(plate_col))$p.value,
                        error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(row_fdr = p.adjust(row_kw_p, "BH"), col_fdr = p.adjust(col_kw_p, "BH"))
write.csv(plate_stats, tab("Plate_Position_Tests.csv"), row.names = FALSE)
cat(sprintf(" Plate-position effect on library size: %d/%d plates row-associated,\n",
            sum(plate_stats$row_fdr < 0.05, na.rm = TRUE), nrow(plate_stats)))
cat(sprintf(" %d/%d column-associated (Kruskal-Wallis, BH-adjusted p < 0.05).\n",
            sum(plate_stats$col_fdr < 0.05, na.rm = TRUE), nrow(plate_stats)))

## 1c-ter. NEW: Viral-load spatial analysis (cross-well contamination check)
cat(" -> Testing viral-load spatial clustering by plate position...\n")
spatial_viral <- cell_qc %>%
  filter(moi > 0, !is.na(n_virus_molecules)) %>%
  mutate(plate_row = factor(plate_row, levels = LETTERS[1:16]),
         plate_col = factor(plate_col)) %>%
  group_by(virus, batch) %>%
  summarise(
    row_kw_p = tryCatch(kruskal.test(n_virus_molecules, plate_row)$p.value,
                        error = function(e) NA_real_),
    col_kw_p = tryCatch(kruskal.test(n_virus_molecules, plate_col)$p.value,
                        error = function(e) NA_real_),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(row_fdr = p.adjust(row_kw_p, "BH"), col_fdr = p.adjust(col_kw_p, "BH"))
write.csv(spatial_viral, tab("Viral_Load_Spatial_Tests.csv"), row.names = FALSE)

for (v in unique(spatial_viral$virus))
  cat(sprintf(" %-7s: row-assoc %d/%d plates, col-assoc %d/%d plates (FDR<0.05)\n",
              v,
              sum(spatial_viral$virus == v & spatial_viral$row_fdr < 0.05, na.rm = TRUE),
              sum(spatial_viral$virus == v),
              sum(spatial_viral$virus == v & spatial_viral$col_fdr < 0.05, na.rm = TRUE),
              sum(spatial_viral$virus == v)))

p4d <- cell_qc %>%
  filter(moi > 0, !is.na(n_virus_molecules)) %>%
  ggplot(aes(x = plate_col, y = reorder(plate_row, as.numeric(plate_row)),
             fill = log10(n_virus_molecules + 1))) +
  geom_tile(colour = "grey90", linewidth = 0.15) +
  facet_grid(virus ~ batch, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(name = "log10\nviral mol+1", option = "inferno") +
  theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 5, angle = 90),
        strip.text.y = element_text(size = 6)) +
  labs(title = "Viral load by plate position",
       subtitle = "384-well layout. Spatial clustering (significant row/col gradient) indicates cross-well contamination; scattered pattern indicates true biological variation.",
       x = "Column", y = "Row")
show_and_save(p4d, "Verification_Viral_Spatial", width = 20, height = 12)

## 1c-quater. NEW: Dengue-only spatial auto-correlation (Moran's I on significant hits)

top_n_feat <- 30
row_tot <- numeric(nrow(count_mat))
for (ch in split(seq_len(nrow(count_mat)), ceiling(seq_len(nrow(count_mat)) / 5000)))
  row_tot[ch] <- rowSums(count_mat[ch, , drop = FALSE])
names(row_tot) <- genes
grand_total <- sum(row_tot)

top_ix <- order(row_tot, decreasing = TRUE)[seq_len(top_n_feat)]
top_feat <- data.frame(
  feature = genes[top_ix],
  total = row_tot[top_ix],
  pct = 100 * row_tot[top_ix] / grand_total,
  stringsAsFactors = FALSE
) %>%
  mutate(class = case_when(str_detect(feature, "^__") ~ "HTSeq QC row",
                           str_detect(feature, "^ERCC-") ~ "ERCC spike-in",
                           feature == "NIST_ConsensusVector" ~ "NIST vector",
                           TRUE ~ "Host (ENSG)"))
write.csv(top_feat, tab("Highest_Expressed_Features.csv"), row.names = FALSE)
cat(sprintf(" Top %d features account for %.1f%% of all counts.\n",
            top_n_feat, sum(top_feat$pct)))

p4d <- ggplot(top_feat, aes(x = pct, y = reorder(feature, pct), fill = class)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("HTSeq QC row" = "#D55E00", "ERCC spike-in" = "#E69F00",
                               "NIST vector" = "#8B3A8B", "Host (ENSG)" = "#A6A6A6"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = sprintf("Top %d features by share of total counts", top_n_feat),
       subtitle = sprintf("Across all %d cells.", ncol(count_mat)),
       x = "% of all counts", y = NULL)
show_and_save(p4d, "Verification_Highest_Expressed", width = 9, height = 7)

htseq_long <- htseq_by_cat %>%
  pivot_longer(-Sample, names_to = "category", values_to = "counts") %>%
  left_join(cell_qc %>% select(Sample, virus, host_counts, ercc_counts, other_counts,
                               htseq_qc_counts), by = "Sample") %>%
  mutate(total = host_counts + ercc_counts + other_counts + htseq_qc_counts,
         pct = 100 * counts / total)

p4e <- ggplot(htseq_long, aes(x = reorder(category, -pct, FUN = median), y = pct,
                              fill = category)) +
  geom_violin(scale = "width", colour = "#4D4D4D", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white", linewidth = 0.3) +
  facet_wrap(~ virus) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Per-cell mapping outcome across the five HTSeq categories",
       subtitle = "Percent of all counted reads in each cell. Explains where the unassigned majority goes.",
       x = NULL, y = "% of counted reads in that cell")
show_and_save(p4e, "Verification_Alignment_Outcome", width = 10, height = 5.5)

## 1e. Ensembl REST spot-check
if (run_ensembl_check) {
  cat(" -> Spot-checking 20 ENSG IDs via Ensembl REST API...\n")
  set.seed(42)
  ensg_20 <- sample(genes[idx_host], 20)
  chk <- vapply(ensg_20, function(g) {
    res <- tryCatch(GET(paste0("https://rest.ensembl.org/lookup/id/", g),
                        accept("application/json"), timeout(10)),
                    error = function(e) NULL)
    Sys.sleep(0.15)
    if (is.null(res)) "Unreachable (network/timeout)"
    else if (status_code(res) == 200) "Found"
    else if (status_code(res) == 400 || status_code(res) == 404) "Not found in Ensembl"
    else sprintf("HTTP %d", status_code(res))
  }, character(1))
  n_found <- sum(chk == "Found")
  n_absent <- sum(chk == "Not found in Ensembl")
  n_unreach <- sum(chk != "Found" & chk != "Not found in Ensembl")
  cat(sprintf(" %d found | %d genuinely absent | %d unreachable\n",
              n_found, n_absent, n_unreach))
  if (n_unreach > 0)
    cat(sprintf(" NOTE: %d calls did not reach the server; inconclusive for those IDs.\n",
                n_unreach))
  write.csv(data.frame(ensg = names(chk), status = unname(chk)),
            tab("Ensembl_Spotcheck.csv"), row.names = FALSE)
} else {
  cat(" -> Ensembl REST spot-check skipped.\n")
}

# ==============================================================================
# Map Ensembl IDs to Gene Symbols
# ==============================================================================
cat("\nMap Ensembl IDs\n")

ensg_all <- genes[idx_host]
ensg_clean <- str_remove(ensg_all, "\\..*$")
cat(sprintf(" %d host IDs (%d carried a version suffix).\n",
            length(ensg_all), sum(ensg_all != ensg_clean)))

cat(" -> Querying org.Hs.eg.db for SYMBOL...\n")
keys_u <- unique(ensg_clean)
sym_list <- suppressMessages(
  AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                        keys = keys_u, column = "SYMBOL",
                        keytype = "ENSEMBL", multiVals = "list")
)
n_sym <- vapply(sym_list, function(x) sum(!is.na(x)), integer(1))
first_ <- vapply(sym_list, function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA_character_ },
                 character(1))

m <- match(ensg_clean, keys_u)
map_df <- data.frame(
  ENSG = ensg_all,
  ENSG_Clean = ensg_clean,
  Symbol = unname(first_[m]),
  N_Symbols = unname(n_sym[m]),
  stringsAsFactors = FALSE
) %>%
  mutate(Status = if_else(is.na(Symbol), "Unmapped", "Mapped (1:1:many handled)"))

write.csv(map_df, tab("Gene_Mapping_Table.csv"), row.names = FALSE)
cat(sprintf(" Mapped %d | Unmapped %d | of mapped, %d are 1:many (first symbol kept)\n",
            sum(map_df$Status != "Unmapped"), sum(map_df$Status == "Unmapped"),
            sum(map_df$N_Symbols > 1)))

map_tab <- as.data.frame(table(Status = map_df$Status), stringsAsFactors = FALSE)
p5 <- ggplot(map_tab, aes(x = Status, y = Freq, fill = Status)) +
  geom_col(colour = "black", width = 0.6) +
  geom_text(aes(label = Freq), vjust = -0.4, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c("Mapped (1:1:many handled)" = "#009E73",
                               "Unmapped" = "#D55E00")) +
  theme(legend.position = "none") +
  labs(title = "Ensembl ID mapping outcome",
       subtitle = sprintf("org.Hs.eg.db ENSEMBL -> SYMBOL over %d host IDs", nrow(map_df)),
       x = NULL, y = "Number of genes")
show_and_save(p5, "Verification_Mapping_Outcome", width = 7, height = 5)

## Expression of unmapped IDs
gene_totals <- rowSums(count_mat)
names(gene_totals) <- genes
unmapped_genes <- map_df$ENSG[map_df$Status == "Unmapped"]

if (length(unmapped_genes)) {
  expr_df <- data.frame(
    gene = map_df$ENSG,
    Status = if_else(map_df$Status == "Unmapped", "Unmapped", "Mapped"),
    Expr = unname(gene_totals[map_df$ENSG]),
    stringsAsFactors = FALSE
  )
  un <- expr_df$Expr[expr_df$Status == "Unmapped"]
  cat(sprintf(" Unmapped IDs: %d total, %d with zero counts across all cells (%.1f%%)\n",
              length(un), sum(un == 0), 100 * mean(un == 0)))
  
  p6 <- ggplot(expr_df, aes(x = log10(Expr + 1), fill = Status)) +
    geom_histogram(bins = 60, colour = "black",
                   position = "identity", alpha = 0.65) +
    scale_fill_manual(values = c(Unmapped = "#D55E00", Mapped = "#A6A6A6"), name = NULL) +
    labs(title = "Total expression of unmapped vs mapped ENSG IDs",
         subtitle = "Summed across all cells; unmapped IDs are overwhelmingly low-expression",
         x = "log10(total counts + 1)", y = "Number of genes")
  show_and_save(p6, "Verification_Unmapped_Expr", width = 7.5, height = 5)
}

# ==============================================================================
# Persist for Phase B
# ==============================================================================
if (save_count_matrix) {
  cat("\n -> Saving host + ERCC count matrix to .rds...\n")
  count_mat_save <- count_mat[idx_host | idx_ercc, , drop = FALSE]
  invisible(gc())
  saveRDS(count_mat_save, obj("Phase_A_counts_host_ercc.rds"))
  saveRDS(cell_qc, obj("Phase_A_cell_qc.rds"))
  saveRDS(map_df, obj("Phase_A_gene_map.rds"))
  cat(sprintf(" Saved %d features x %d cells.\n", nrow(count_mat_save), ncol(count_mat_save)))
}

cat("\nPhase A Pre-processing completed.\n")
cat(sprintf("All figures (.svg) and tables (.csv) are in: %s\n", getwd()))