# ==============================================================================
# PHASE B - QUALITY CONTROL 
# ==============================================================================
# Per-cell QC metrics, per-stratum filtering, feature space cleaning.
#
# Inputs (from Phase A):
# results/phaseA/objects/Phase_A_counts_host_ercc.rds
# results/phaseA/objects/Phase_A_cell_qc.rds
# results/phaseA/objects/Phase_A_gene_map.rds
#
# Outputs (to Phase C):
# results/phaseB/objects/Phase_B_counts_active.rds
# results/phaseB/objects/Phase_B_cell_qc.rds
# results/phaseB/objects/Phase_B_gene_sets.rds
#
# NOTE: all figures are written as vector SVG only (results/phaseB/figures/*.svg)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  source("E:/Zika/scripts/inject_themes.R")
  library(scales)
  library(httr)
  library(jsonlite)
  library(patchwork)
  library(svglite)
})

if (!requireNamespace("svglite", quietly = TRUE))
  stop("svglite is not installed.\n",
       " install.packages('svglite')")

# Guard against S4 generics masking dplyr verbs
.verbs <- c("select", "rename", "filter", "mutate", "slice", "desc", "first")
.owner <- vapply(.verbs, function(f) find(f)[1], character(1))
.masked <- .owner[!.owner %in% c("package:dplyr", "package:tidyr")]
if (length(.masked)) {
  stop("dplyr verbs are masked: ",
       paste(sprintf("%s() <- %s", names(.masked), .masked), collapse = ", "),
       "\nIf the owner is .GlobalEnv, remove that object (or restart R) and re-source.")
}
rm(.verbs, .owner, .masked)

# ==============================================================================
# Configuration and output layout
# ==============================================================================
data_dir <- "E:/Zika"
setwd(data_dir)

nmads <- 3
min_cells_detected <- 3
hvg_n <- 2000
run_hgnc_check <- TRUE

FLOOR_HOST_COUNTS <- 20000
FLOOR_N_GENES <- 2000

results_root <- file.path(data_dir, "results", "phaseB")
dir_fig <- file.path(results_root, "figures")
dir_tab <- file.path(results_root, "tables")
dir_obj <- file.path(results_root, "objects")
for (d in c(dir_fig, dir_tab, dir_obj))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

fig <- function(n) file.path(dir_fig, n)
tab <- function(n) file.path(dir_tab, n)
obj <- function(n) file.path(dir_obj, n)

# --------------------------------------------------------------------------
# show_and_save(): Apply publication theme, print, and save as SVG.
# --------------------------------------------------------------------------
show_and_save <- function(p, file, width = 8, height = 5, base_size = 10) {
  p <- p + theme_pub(base_size = base_size)
  print(p)
  ggsave(fig(paste0(file, ".svg")), p, width = width, height = height,
         bg = "transparent", device = svglite::svglite)
  ggsave(fig(paste0(file, ".png")), p, width = width, height = height,
         dpi = 600, bg = "white")
  cat(sprintf(" [plot] %s.svg and .png\n", file))
  invisible(p)
}
write_tab <- function(x, name) {
  write.csv(x, tab(name), row.names = FALSE)
  cat(sprintf(" [table] %s\n", name))
}

min_pos <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) min(x) else 1 }

mad_threshold <- function(x, nmads = 3, side = c("lower", "higher"),
                          log = FALSE, pseudo = NULL) {
  side <- match.arg(side)
  x <- x[is.finite(x)]
  if (is.null(pseudo)) pseudo <- min_pos(x) / 2
  v <- if (log) log10(x + pseudo) else x
  med <- median(v); dev <- mad(v, center = med)
  thr <- if (side == "lower") med - nmads * dev else med + nmads * dev
  if (log) 10^thr - pseudo else thr
}

find_file <- function(bn, subdirs = c("results/phaseA/objects", "results/phaseA/tables",
                                      "results/phaseA/cache", ".")) {
  for (s in subdirs) {
    p <- file.path(data_dir, s, bn)
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NA_character_
}

cat("======================================================================\n")
cat("PHASE B - QUALITY CONTROL \n")
cat("======================================================================\n")

# ------------------------------------------------------------------------------
# Load and validate Phase A outputs
# ------------------------------------------------------------------------------
p_counts <- find_file("Phase_A_counts_host_ercc.rds")
p_qc <- find_file("Phase_A_cell_qc.rds")
p_map <- find_file("Phase_A_gene_map.rds")
if (is.na(p_counts) || is.na(p_qc))
  stop("Missing Phase A outputs. Need Phase_A_counts_host_ercc.rds and ",
       "Phase_A_cell_qc.rds under results/phaseA/objects/ or the project root.")

cat(sprintf("Counts : %s\n", p_counts))
cat(sprintf("Cell QC: %s\n", p_qc))
counts <- readRDS(p_counts)
cell_qc <- readRDS(p_qc)

gene_map <- if (!is.na(p_map)) readRDS(p_map) else NULL
if (!is.null(gene_map)) {
  ensg_col <- intersect(c("ENSG", "ENSG_Original", "ENSG_Clean"), names(gene_map))[1]
  if (is.na(ensg_col) || !"Symbol" %in% names(gene_map)) gene_map <- NULL else
    cat(sprintf("Gene map: %s (ID column '%s')\n", basename(p_map), ensg_col))
}
sym_of <- function(ids) {
  if (is.null(gene_map)) return(ids)
  s <- gene_map$Symbol[match(ids, gene_map[[ensg_col]])]
  ifelse(is.na(s), ids, s)
}

stopifnot(identical(colnames(counts), cell_qc$Sample))
cat(sprintf("Loaded %d features x %d cells\n", nrow(counts), ncol(counts)))

need <- c("Sample", "host_counts", "ercc_counts", "virus", "time_h", "moi",
          "n_virus_molecules", "batch", "plate_row", "plate_col")
lack <- setdiff(need, colnames(cell_qc))
if (length(lack))
  stop("Phase_A_cell_qc.rds lacks: ", paste(lack, collapse = ", "),
       "\nRe-run Phase A with the updated script.")

# ==============================================================================
# COMPUTE QC DISTRIBUTIONS
# ==============================================================================
cat("\ncomputing QC distributions\n")

feat <- rownames(counts)
is_ercc <- str_detect(feat, "^ERCC-\\d+")
is_host <- str_detect(feat, "^ENSG\\d+")

mt_chrM <- c("ENSG00000210049", "ENSG00000211459", "ENSG00000210077", "ENSG00000210082",
             "ENSG00000209082", "ENSG00000198888", "ENSG00000210100", "ENSG00000210107",
             "ENSG00000210112", "ENSG00000198763", "ENSG00000210117", "ENSG00000210127",
             "ENSG00000210135", "ENSG00000210140", "ENSG00000210144", "ENSG00000198804",
             "ENSG00000210151", "ENSG00000210154", "ENSG00000198712", "ENSG00000210156",
             "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000210164",
             "ENSG00000198840", "ENSG00000210174", "ENSG00000212907", "ENSG00000198886",
             "ENSG00000210176", "ENSG00000210184", "ENSG00000210191", "ENSG00000198786",
             "ENSG00000198695", "ENSG00000210194", "ENSG00000198727", "ENSG00000210195",
             "ENSG00000210196")
mt_genes <- intersect(mt_chrM, feat)
cat(sprintf(" MT genes: %d present\n", length(mt_genes)))
if (length(mt_genes) == 0) stop("No MT genes found; %MT would be identically zero.")
is_mt <- feat %in% mt_genes

csum <- function(idx, fun = identity) {
  out <- numeric(ncol(counts))
  if (!any(idx)) return(out)
  for (ch in split(which(idx), ceiling(seq_len(sum(idx)) / 5000)))
    out <- out + colSums(fun(counts[ch, , drop = FALSE]))
  out
}

## 3a. Per-cell percentages
## Ungroup first: Phase A RDS may carry a grouped_df structure that breaks
## downstream summarise / across calls.
qc <- cell_qc %>% ungroup() %>%
  mutate(
    mt_counts = csum(is_mt),
    n_genes_detected = csum(is_host, function(m) m > 0),
    n_virus_molecules = replace_na(n_virus_molecules, 0),
    pct_mt = 100 * mt_counts / host_counts,
    pct_ercc = 100 * ercc_counts / (host_counts + ercc_counts),
    pct_virus = 100 * n_virus_molecules / (n_virus_molecules + host_counts)
  )

## 3c. Metadata completeness
miss <- qc %>% mutate(m_hpi = is.na(time_h), m_moi = is.na(moi),
                      m_vir = is.na(virus) | virus == "")
n_miss <- sum(miss$m_hpi | miss$m_moi | miss$m_vir)
cat(sprintf(" Metadata completeness: %d missing hpi, %d missing MOI, %d missing virus\n",
            sum(miss$m_hpi), sum(miss$m_moi), sum(miss$m_vir)))
if (n_miss > 0) {
  cat(" WARNING: incomplete cells written to Missing_Metadata.csv\n")
  write_tab(miss %>% filter(m_hpi | m_moi | m_vir) %>%
              select(Sample, virus, time_h, moi), "Missing_Metadata.csv")
} else {
  cat(" PASS: every cell carries hpi, MOI and virus.\n")
}

qc_summary <- qc %>%
  summarise(across(c(host_counts, n_genes_detected, pct_mt, pct_ercc, pct_virus),
                   list(min = ~min(.x), p25 = ~quantile(.x, .25, names = FALSE),
                        med = ~median(.x), p75 = ~quantile(.x, .75, names = FALSE),
                        p99 = ~quantile(.x, .99, names = FALSE), max = ~max(.x)))) %>%
  pivot_longer(everything(), names_to = c("metric", "stat"),
               names_pattern = "(.*)_(min|p25|med|p75|p99|max)$") %>%
  pivot_wider(names_from = stat, values_from = value)
print(as.data.frame(qc_summary), digits = 4)
write_tab(qc_summary, "Distribution_Summary.csv")
write_tab(qc %>% select(Sample, virus, time_h, moi, batch,
                        host_counts, ercc_counts, mt_counts, n_genes_detected,
                        n_virus_molecules, pct_mt, pct_ercc, pct_virus),
          "Cell_Distributions.csv")

## 3d. Violins of %MT / %virus / %ERCC by hpi and MOI
qc_long <- qc %>%
  filter(!is.na(time_h), !is.na(moi)) %>%
  select(Sample, virus, time_h, moi, pct_mt, pct_ercc, pct_virus) %>%
  pivot_longer(c(pct_mt, pct_ercc, pct_virus), names_to = "metric", values_to = "pct") %>%
  mutate(metric = recode(metric, pct_mt = "% mitochondrial",
                         pct_ercc = "% ERCC spike-in", pct_virus = "% viral"))

off <- min_pos(qc_long$pct)
p3a <- ggplot(qc_long, aes(x = factor(time_h), y = pct + off, fill = factor(moi))) +
  geom_violin(scale = "width", position = position_dodge(width = 0.85),
              colour = "#4D4D4D", linewidth = 0.3) +
  geom_boxplot(aes(group = interaction(time_h, moi)),
               position = position_dodge(width = 0.85), width = 0.15, outlier.size = 0.4,
               fill = "white", linewidth = 0.3) +
  facet_grid(metric ~ virus, scales = "free_y") +
  scale_y_log10() + scale_fill_brewer(palette = "Set2", name = "MOI") +
  labs(title = "QC percentages by hours post-infection and MOI",
       subtitle = sprintf("All %d cells, unfiltered. Zeros offset by %.2g so they remain visible on the log axis.", nrow(qc), off),
       x = "Hours post-infection", y = "Percent (log10)")
show_and_save(p3a, "QC_Violins", width = 10, height = 8)

## 3e. %virus vs %MT, coloured by hpi
off_virus <- min_pos(qc$pct_virus)
p3b <- qc %>% filter(!is.na(time_h)) %>%
  ggplot(aes(x = pct_virus + off_virus, y = pct_mt, colour = factor(time_h))) +
  geom_point(alpha = 0.55, size = 1.2) +
  facet_wrap(~ virus) + scale_x_log10() +
  scale_colour_viridis_d(name = "hpi", option = "plasma", end = 0.9) +
  labs(title = "viral load vs mitochondrial fraction per cell",
       subtitle = sprintf("Unfiltered. x offset by %g so uninfected cells remain visible.", off_virus),
       x = "% viral (offset, log10)", y = "% mitochondrial")
show_and_save(p3b, "Virus_vs_MT", width = 9, height = 5)

## 3f. Spatial plate layout for Dengue (cross-contamination check)
dengue_qc <- qc %>% filter(virus == "dengue", !is.na(batch)) %>%
  mutate(
    plate_row_fac = factor(plate_row, levels = rev(LETTERS[1:16])),
    plate_col_num = as.integer(plate_col)
  )

if (nrow(dengue_qc) > 0) {
  deng_qc <- dengue_qc %>%
    mutate(concordance = case_when(
      n_virus_molecules > 0 ~ "Virus detected (seq +ve)",
      TRUE ~ "No virus molecules (seq -ve)"
    ))
  
  p3c_a <- ggplot(deng_qc, aes(x = plate_col_num, y = plate_row_fac,
                               fill = log10(n_virus_molecules + 1))) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "log10(vRNA + 1)", option = "plasma") +
    scale_x_continuous(breaks = seq(2, 24, 2)) +
    facet_wrap(~ batch) +
    theme(panel.grid = element_blank()) +
    labs(title = "Dengue spatial distribution (sequencing viral load)",
         subtitle = "Tile colour = sequencing reads.",
         x = "Plate Column", y = "Plate Row")
  
  discord_only <- deng_qc %>% filter(concordance == "No virus molecules (seq -ve)")
  p3c_b <- if (nrow(discord_only) > 0) {
    ggplot(discord_only, aes(x = plate_col_num, y = plate_row_fac)) +
      geom_tile(fill = "#D55E00", color = "white", alpha = 0.85) +
      scale_x_continuous(breaks = seq(2, 24, 2)) +
      scale_y_discrete(limits = rev(levels(deng_qc$plate_row_fac))) +
      facet_wrap(~ batch) +
      theme(panel.grid = element_blank()) +
      labs(title = paste0("Cells with no virus molecules (n = ", nrow(discord_only), ")"),
           subtitle = "Red tiles = cells with zero sequencing reads from virus.",
           x = "Plate Column", y = "Plate Row")
  } else {
    ggplot() + labs(title = "No cells without virus molecules")
  }
  
  p3c <- (p3c_a | p3c_b) +
    plot_annotation(
      title = "Dengue Sequencing Viral Load by Well Position",
      subtitle = paste0("Left: sequencing viral load. Right: ", nrow(discord_only),
                        " cells with zero virus molecules detected."),
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
  show_and_save(p3c, "Dengue_Spatial_Plate_Layout", width = 14, height = 8)
  
  concord_tbl <- deng_qc %>% count(concordance)
  write_tab(concord_tbl, "Dengue_Sequencing_Concordance.csv")
  cat(sprintf(" Dengue concordance: %s\n",
              paste(concord_tbl$concordance, concord_tbl$n, collapse = "; ")))
}

# ==============================================================================
# EMPIRICAL CELL FILTERING (per virus x hpi stratum)
# ==============================================================================
cat("\nEmpirical cell filtering\n")

# Per-stratum cutoffs with a floor that prevents relaxation below a sane minimum.
# The stratum may only RELAX the pooled cutoff, never tighten it.
pooled_lib <- mad_threshold(qc$host_counts, nmads, "lower", log = TRUE)
pooled_genes <- mad_threshold(qc$n_genes_detected, nmads, "lower", log = TRUE)
## Upper-tail thresholds (ERCC %, %MT): use raw scale. Log transform compresses
## the high-value tail and produces near-zero cutoffs that flag every cell.
thr_ercc <- mad_threshold(qc$pct_ercc, nmads, "higher", log = FALSE)

strata_thr <- qc %>% filter(!is.na(time_h), !is.na(virus)) %>%
  group_by(virus, time_h) %>%
  summarise(n = dplyr::n(),
            median_mt = median(pct_mt),
            thr_mt = mad_threshold(pct_mt, nmads, "higher", log = FALSE),
            median_lib = median(host_counts),
            mad_lib = mad_threshold(host_counts, nmads, "lower", log = TRUE),
            mad_genes = mad_threshold(n_genes_detected, nmads, "lower", log = TRUE),
            .groups = "drop") %>%
  mutate(pooled_lib_ref = .env$pooled_lib,
         pooled_genes_ref = .env$pooled_genes,
         thr_lib = pmin(pooled_lib_ref, pmax(mad_lib, FLOOR_HOST_COUNTS)),
         thr_genes = pmin(pooled_genes_ref, pmax(mad_genes, FLOOR_N_GENES)),
         lib_bound_by = case_when(
           abs(thr_lib - FLOOR_HOST_COUNTS) < 1e-6 ~ "floor",
           abs(thr_lib - pooled_lib_ref) < 1e-6 ~ "pooled",
           TRUE ~ "stratum MAD"),
         gene_bound_by = case_when(
           abs(thr_genes - FLOOR_N_GENES) < 1e-6 ~ "floor",
           abs(thr_genes - pooled_genes_ref) < 1e-6 ~ "pooled",
           TRUE ~ "stratum MAD"))

cat(sprintf(" Pooled cutoffs (REPORTED, NOT APPLIED): library < %.0f | genes < %.0f | %%ERCC < %.3f\n",
            pooled_lib, pooled_genes, thr_ercc))
cat(sprintf(" Applied floors: library >= %d | genes >= %d\n", FLOOR_HOST_COUNTS, FLOOR_N_GENES))
cat(" Per-stratum cutoffs applied:\n")
cat(sprintf(" %-7s %3s %5s %10s %10s %-12s %8s %8s %-12s %8s\n",
            "virus", "hpi", "n", "med lib", "lib cutoff", "set by", "med %MT", "MT cutoff", "gene set by", "gene cut"))
for (i in seq_len(nrow(strata_thr)))
  cat(sprintf(" %-7s %3d %5d %10.0f %10.0f %-12s %8.2f %8.2f %-12s %8.0f\n",
              strata_thr$virus[i], strata_thr$time_h[i], strata_thr$n[i],
              strata_thr$median_lib[i], strata_thr$thr_lib[i], strata_thr$lib_bound_by[i],
              strata_thr$median_mt[i], strata_thr$thr_mt[i],
              strata_thr$gene_bound_by[i], strata_thr$thr_genes[i]))
write_tab(strata_thr, "Thresholds_Per_Stratum.csv")
write_tab(strata_thr %>% select(virus, time_h, n, median_mt, thr_mt),
          "MT_Thresholds_Per_Stratum.csv")

qc <- qc %>%
  left_join(strata_thr %>% select(virus, time_h, thr_mt, thr_lib, thr_genes),
            by = c("virus", "time_h")) %>%
  mutate(fail_lib = !is.na(thr_lib) & host_counts < thr_lib,
         fail_genes = !is.na(thr_genes) & n_genes_detected < thr_genes,
         fail_ercc = pct_ercc > thr_ercc,
         fail_mt = !is.na(thr_mt) & pct_mt > thr_mt,
         keep = !(fail_lib | fail_genes | fail_ercc | fail_mt))

## Diagnostic: show how many cells fail each criterion and the threshold values
cat(sprintf("\n--- QC threshold diagnostics ---\n"))
cat(sprintf(" fail_lib : %5d cells (thr = %.0f)\n", sum(qc$fail_lib), unique(qc$thr_lib[!is.na(qc$thr_lib)][1])))
cat(sprintf(" fail_genes: %5d cells (thr = %.0f)\n", sum(qc$fail_genes), unique(qc$thr_genes[!is.na(qc$thr_genes)][1])))
cat(sprintf(" fail_ercc : %5d cells (thr = %.3f%%)\n", sum(qc$fail_ercc), thr_ercc))
cat(sprintf(" fail_mt : %5d cells (median thr = %.2f%%)\n", sum(qc$fail_mt), median(qc$thr_mt, na.rm = TRUE)))
cat(sprintf(" KEPT : %5d of %d cells (%.1f%%)\n",
            sum(qc$keep), nrow(qc), 100 * mean(qc$keep)))
cat(sprintf("--------------------------------\n\n"))

cat(sprintf("\n Applied: counts & genes per stratum (floored) | %%ERCC global | %%MT per stratum\n"))
cat(sprintf(" Removed by criterion: lib %d, genes %d, ERCC %d, MT %d\n",
            sum(qc$fail_lib), sum(qc$fail_genes), sum(qc$fail_ercc), sum(qc$fail_mt)))
cat(sprintf(" RETAINED %d of %d cells (%.1f%%)\n",
            sum(qc$keep), nrow(qc), 100 * mean(qc$keep)))
write_tab(qc %>% select(Sample, virus, time_h, moi, host_counts, n_genes_detected,
                        pct_mt, pct_ercc, pct_virus, thr_mt, thr_lib, thr_genes,
                        fail_lib, fail_genes, fail_ercc, fail_mt, keep),
          "Cell_Filter_Calls.csv")

## 4d. Before/after histograms per criterion
crit <- bind_rows(
  qc %>% transmute(metric = "Library size (host counts)", value = host_counts,
                   cutoff = thr_lib, side = "lower"),
  qc %>% transmute(metric = "Genes detected", value = n_genes_detected,
                   cutoff = thr_genes, side = "lower"),
  qc %>% transmute(metric = "% ERCC", value = pct_ercc,
                   cutoff = thr_ercc, side = "higher")) %>%
  mutate(status = if_else((side == "lower" & value < cutoff) |
                            (side == "higher" & value > cutoff), "Removed", "Retained"))

applied_lines <- crit %>% distinct(metric, cutoff) %>% filter(!is.na(cutoff))
ref_lines <- bind_rows(
  data.frame(metric = "Library size (host counts)", x = FLOOR_HOST_COUNTS, kind = "floor"),
  data.frame(metric = "Genes detected", x = FLOOR_N_GENES, kind = "floor"),
  data.frame(metric = "Library size (host counts)", x = pooled_lib, kind = "pooled (not applied)"),
  data.frame(metric = "Genes detected", x = pooled_genes, kind = "pooled (not applied)"))

p4a <- ggplot(crit, aes(x = value, fill = status)) +
  geom_histogram(bins = 60, colour = "black", linewidth = 0.2) +
  geom_vline(data = applied_lines, aes(xintercept = cutoff),
             linetype = "dashed", colour = "red", linewidth = 0.3, alpha = 0.75) +
  geom_vline(data = ref_lines, aes(xintercept = x, colour = kind), linewidth = 0.7) +
  scale_colour_manual(values = c("floor" = "black", "pooled (not applied)" = "#737373"), name = NULL) +
  facet_wrap(~ metric, scales = "free") + scale_x_log10() +
  scale_fill_manual(values = c(Retained = "#0072B2", Removed = "#D55E00"), name = NULL) +
  theme(legend.position = "bottom") +
  labs(title = "distribution and applied cutoffs per criterion",
       subtitle = sprintf("Thin dashed red = per-stratum cutoffs (median - %d MAD). Black = floor. Grey = pooled (not applied).", nmads),
       x = "Value (log10)", y = "Cells")
show_and_save(p4a, "Filtering_Cutoffs", width = 11, height = 5)

p4b <- qc %>% filter(!is.na(time_h)) %>%
  ggplot(aes(x = pct_mt, fill = fail_mt)) +
  geom_histogram(bins = 50, colour = "black", linewidth = 0.2) +
  geom_vline(data = strata_thr, aes(xintercept = thr_mt), linetype = "dashed", colour = "red") +
  facet_grid(virus ~ time_h, labeller = label_both, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"),
                    labels = c("Retained", "Removed"), name = NULL) +
  theme(legend.position = "bottom") +
  labs(title = "%MT thresholded within each virus x hpi stratum",
       subtitle = "Each panel gets its own cutoff, so late-hpi cells with genuinely elevated %MT are preserved.",
       x = "% mitochondrial", y = "Cells")
show_and_save(p4b, "MT_Per_Stratum", width = 11, height = 5.5)

## 4e. Cells retained per timepoint
retained <- qc %>% filter(!is.na(time_h)) %>%
  group_by(virus, time_h, moi) %>%
  summarise(before = dplyr::n(), after = sum(keep),
            pct_kept = 100 * mean(keep), .groups = "drop")
write_tab(retained, "Cells_Retained.csv")

p4c <- retained %>%
  pivot_longer(c(before, after), names_to = "stage", values_to = "n") %>%
  mutate(stage = factor(recode(stage, before = "Before QC", after = "After QC"),
                        levels = c("Before QC", "After QC"))) %>%
  ggplot(aes(x = factor(time_h), y = n, fill = stage)) +
  geom_col(position = position_dodge(0.85), width = 0.78, colour = "black") +
  facet_grid(virus ~ moi, labeller = label_both) +
  scale_fill_manual(values = c("Before QC" = "#CCCCCC", "After QC" = "#0072B2"), name = NULL) +
  theme(legend.position = "bottom") +
  labs(title = "cells retained per timepoint",
       subtitle = sprintf("%d of %d retained overall (%.1f%%)",
                          sum(qc$keep), nrow(qc), 100 * mean(qc$keep)),
       x = "Hours post-infection", y = "Cells")
show_and_save(p4c, "Cells_Retained", width = 10, height = 6)

## 4f. Sensitivity against generic thresholds (reference only)
sens <- data.frame(
  Threshold = c("Per-stratum + floor, APPLIED",
                sprintf("Pooled MAD on depth (not applied): counts < %.0f or genes < %.0f",
                        pooled_lib, pooled_genes),
                "Generic: host counts < 20000",
                "Generic: %MT > 15",
                "Generic: %MT > 5",
                "Generic: nFeature > 10000 ceiling"),
  Removed = c(sum(!qc$keep),
              sum(qc$host_counts < pooled_lib | qc$n_genes_detected < pooled_genes),
              sum(qc$host_counts < 20000),
              sum(qc$pct_mt > 15),
              sum(qc$pct_mt > 5),
              sum(qc$n_genes_detected > 10000)),
  stringsAsFactors = FALSE) %>%
  mutate(Retained = nrow(qc) - Removed, Pct_removed = 100 * Removed / nrow(qc))
cat(" Sensitivity vs generic thresholds (reference only, none applied):\n")
for (i in seq_len(nrow(sens)))
  cat(sprintf(" %-50s %4d (%5.1f%%)\n", sens$Threshold[i],
              sens$Removed[i], sens$Pct_removed[i]))
cat(sprintf(" Max genes detected in any cell: %d. The 10000 ceiling is inert here.\n",
            max(qc$n_genes_detected)))
write_tab(sens, "Sensitivity_Generic_Thresholds.csv")

p4d <- ggplot(sens, aes(x = reorder(Threshold, Removed), y = Removed,
                        fill = Threshold == "Per-stratum + floor, APPLIED")) +
  geom_col(colour = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", Removed, Pct_removed)),
            hjust = -0.08, size = 3) +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#B0B0B0"), guide = "none") +
  labs(title = "sensitivity check, empirical vs generic thresholds",
       subtitle = sprintf("Out of %d cells. Generic values are shown for comparison and are NOT applied.", nrow(qc)),
       x = NULL, y = "Cells removed")
show_and_save(p4d, "Sensitivity", width = 9, height = 4.5)

# ==============================================================================
# CLEAN FEATURE SPACE
# ==============================================================================
cat("\nCleaning the feature space\n")

keep_cells <- qc$Sample[qc$keep]
counts_f <- counts[, keep_cells, drop = FALSE]
rm(counts); invisible(gc())

detected <- rowSums(counts_f > 0) > min_cells_detected
mat_before <- counts_f[detected, , drop = FALSE]
cat(sprintf(" %d cells retained; %d features detected in > %d cells\n",
            ncol(mat_before), nrow(mat_before), min_cells_detected))

## 5a. Ribosomal gene set, verified against HGNC
ribo_regex <- if (is.null(gene_map)) character(0) else
  intersect(gene_map[[ensg_col]][!is.na(gene_map$Symbol) &
                                   str_detect(gene_map$Symbol, "^RP[LS]\\d")],
            rownames(mat_before))
ribo_genes <- ribo_regex
ribo_source <- "symbol regex (UNVERIFIED)"

if (run_hgnc_check && requireNamespace("jsonlite", quietly = TRUE)) {
  cat(" -> Fetching HGNC gene groups 728 / 729...\n")
  hg <- lapply(c(728, 729), function(gid) {
    r <- tryCatch(GET(paste0("https://rest.genenames.org/fetch/gene_group_id/", gid),
                      add_headers(Accept = "application/json"), timeout(30)),
                  error = function(e) NULL)
    Sys.sleep(1)
    if (is.null(r) || status_code(r) != 200) return(NULL)
    d <- jsonlite::fromJSON(content(r, as = "text", encoding = "UTF-8"))$response$docs
    if (is.null(d$ensembl_gene_id)) NULL else
      data.frame(symbol = d$symbol, ensembl = d$ensembl_gene_id, group = gid,
                 stringsAsFactors = FALSE)
  })
  if (!any(vapply(hg, is.null, logical(1)))) {
    ribo_ref <- bind_rows(hg) %>% filter(!is.na(ensembl), ensembl != "")
    write_tab(ribo_ref, "HGNC_Ribosomal_Reference.csv")
    hgnc_in_mat <- intersect(ribo_ref$ensembl, rownames(mat_before))
    only_hgnc <- setdiff(hgnc_in_mat, ribo_regex)
    only_regex <- setdiff(ribo_regex, ribo_ref$ensembl)
    cat(sprintf(" HGNC %d curated | %d survive detection filter | regex found %d\n",
                nrow(ribo_ref), length(hgnc_in_mat), length(ribo_regex)))
    cat(sprintf(" regex misses %d HGNC genes (%s); regex adds %d not in HGNC\n",
                length(only_hgnc), paste(head(sym_of(only_hgnc), 8), collapse = ", "),
                length(only_regex)))
    dropped <- setdiff(ribo_ref$ensembl, rownames(mat_before))
    if (length(dropped))
      cat(sprintf(" %d curated gene(s) below detection filter: %s\n",
                  length(dropped), paste(sym_of(dropped), collapse = ", ")))
    write_tab(data.frame(
      ensembl = c(only_hgnc, only_regex, dropped),
      symbol = sym_of(c(only_hgnc, only_regex, dropped)),
      issue = c(rep("in HGNC, missed by regex", length(only_hgnc)),
                rep("matched by regex, not in HGNC", length(only_regex)),
                rep("in HGNC, below detection filter", length(dropped)))),
      "Ribosomal_Crosscheck.csv")
    ribo_genes <- hgnc_in_mat
    ribo_source <- "HGNC groups 728/729"
  } else {
    cat(" HGNC unreachable; falling back to the regex set (UNVERIFIED).\n")
  }
}
cat(sprintf(" Ribosomal set in use: %d genes from %s. RETAINED in the active matrix.\n",
            length(ribo_genes), ribo_source))
if (length(ribo_genes) == 0)
  cat(" WARNING: no ribosomal genes identified; 5d will be skipped.\n")

## 5b. Remove ERCC / viral / MT from the active matrix
rn <- rownames(mat_before)
d_ercc <- str_detect(rn, "^ERCC-\\d+")
d_mt <- rn %in% mt_genes
d_viral <- str_detect(rn, "^(DENV|ZIKV)")
mat_after <- mat_before[!(d_ercc | d_mt | d_viral), , drop = FALSE]
cat(sprintf(" Removed %d ERCC + %d MT + %d viral -> active matrix %d x %d\n",
            sum(d_ercc), sum(d_mt), sum(d_viral), nrow(mat_after), ncol(mat_after)))
stopifnot(!any(rownames(mat_after) %in% mt_genes))

## 5c. PCA loadings before vs after removal
gene_stats <- function(m, lib, chunk = 5000) {
  ng <- nrow(m); mu <- numeric(ng); vr <- numeric(ng)
  for (ch in split(seq_len(ng), ceiling(seq_len(ng) / chunk))) {
    x <- log2(t(t(m[ch, , drop = FALSE]) / lib) * 1e6 + 1)
    mu[ch] <- rowMeans(x)
    vr[ch] <- rowSums((x - mu[ch])^2) / (ncol(x) - 1)
  }
  list(mean = mu, var = vr)
}
run_pca <- function(m, label) {
  lib <- colSums(m); lib[lib == 0] <- 1
  gs <- gene_stats(m, lib)
  sel <- order(gs$var, decreasing = TRUE)[seq_len(min(hvg_n, nrow(m)))]
  x <- log2(t(t(m[sel, , drop = FALSE]) / lib) * 1e6 + 1)
  pc <- prcomp(t(x), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  cat(sprintf(" PCA %-6s: %d HVGs | PC1 %.1f%% PC2 %.1f%%\n",
              label, length(sel), ve[1], ve[2]))
  list(pc = pc, ve = ve, label = label)
}
cat(" -> PCA before removal (ERCC and MT still present)...\n")
pca_before <- run_pca(mat_before, "before")
cat(" -> PCA after removal...\n")
pca_after <- run_pca(mat_after, "after")

top_load <- function(res, npc = 2, ntop = 15) {
  do.call(rbind, lapply(seq_len(npc), function(k) {
    ld <- res$pc$rotation[, k]
    ix <- order(abs(ld), decreasing = TRUE)[seq_len(ntop)]
    g <- names(ld)[ix]
    data.frame(stage = res$label, PC = paste0("PC", k), gene = g,
               symbol = sym_of(g), loading = unname(ld[ix]),
               class = case_when(g %in% mt_genes ~ "MT",
                                 str_detect(g, "^ERCC-") ~ "ERCC",
                                 g %in% ribo_genes ~ "Ribosomal",
                                 TRUE ~ "Other"),
               stringsAsFactors = FALSE)
  }))
}
load_df <- bind_rows(top_load(pca_before), top_load(pca_after)) %>%
  mutate(stage = factor(stage, levels = c("before", "after")))
write_tab(load_df, "PCA_Top_Loadings.csv")
cat(" Feature classes among the top loadings:\n")
print(as.data.frame(load_df %>% count(stage, PC, class) %>%
                      pivot_wider(names_from = class, values_from = n, values_fill = 0)))

p5a <- load_df %>%
  mutate(key = paste(stage, PC, symbol, sep = "|")) %>%
  ggplot(aes(x = loading, y = reorder(key, loading), fill = class)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_y_discrete(labels = function(z) sub("^[^|]*\\|[^|]*\\|", "", z)) +
  facet_wrap(~ stage + PC, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(MT = "#D55E00", ERCC = "#E69F00",
                               Ribosomal = "#009E73", Other = "#A6A6A6")) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "top PCA loadings before vs after removing ERCC and MT",
       subtitle = sprintf("Top 15 loadings on PC1/PC2. Same %d cells and same %d-HVG budget.", ncol(mat_after), hvg_n),
       x = "Loading", y = NULL)
show_and_save(p5a, "PCA_Loadings", width = 11, height = 9)

var_cmp <- data.frame(PC = rep(1:10, 2),
                      pct_variance = c(pca_before$ve[1:10], pca_after$ve[1:10]),
                      stage = rep(c("before", "after"), each = 10))
write_tab(var_cmp, "PCA_Variance_Explained.csv")

## 5d. Ribosomal expression vs hpi
if (length(ribo_genes)) {
  lib_a <- colSums(mat_after); lib_a[lib_a == 0] <- 1
  ribo_df <- data.frame(Sample = colnames(mat_after),
                        pct_ribo = 100 * colSums(mat_after[ribo_genes, , drop = FALSE]) / lib_a,
                        stringsAsFactors = FALSE) %>%
    left_join(qc %>% select(Sample, virus, time_h, moi), by = "Sample")
  write_tab(ribo_df, "Ribosomal_Expression.csv")
  cat(" Ribosomal fraction by stratum:\n")
  print(as.data.frame(ribo_df %>% filter(!is.na(time_h)) %>%
                        group_by(virus, time_h) %>%
                        summarise(n = dplyr::n(), median_pct = median(pct_ribo),
                                  .groups = "drop")), digits = 4)
  
  p5b <- ribo_df %>% filter(!is.na(time_h)) %>%
    ggplot(aes(x = factor(time_h), y = pct_ribo, fill = factor(moi))) +
    geom_violin(scale = "width", position = position_dodge(width = 0.85),
                colour = "#4D4D4D", linewidth = 0.3) +
    geom_boxplot(aes(group = interaction(time_h, moi)),
                 position = position_dodge(width = 0.85), width = 0.13,
                 outlier.size = 0.3, fill = "white", linewidth = 0.3) +
    facet_wrap(~ virus) +
    scale_fill_brewer(palette = "Set2", name = "MOI") +
    theme(legend.position = "bottom") +
    labs(title = "ribosomal protein gene expression vs hours post-infection",
         subtitle = sprintf("%d genes from %s, retained in the active matrix.",
                            length(ribo_genes), ribo_source),
         x = "Hours post-infection", y = "% of active-matrix counts")
  show_and_save(p5b, "Ribosomal_vs_hpi", width = 9, height = 5.5)
}

# ==============================================================================
# Persist for Phase C
# ==============================================================================
saveRDS(mat_after, obj("Phase_B_counts_active.rds"))
saveRDS(qc, obj("Phase_B_cell_qc.rds"))
saveRDS(list(mt = mt_genes, ribosomal = ribo_genes,
             ercc = rownames(mat_before)[d_ercc],
             ribo_source = ribo_source), obj("Phase_B_gene_sets.rds"))

cat("\n======================================================================\n")
cat("PHASE B COMPLETE\n")
cat(sprintf(" Active matrix : %d features x %d cells\n", nrow(mat_after), ncol(mat_after)))
cat(sprintf(" Objects : %s\n", dir_obj))
cat(sprintf(" Figures/tables: %s\n", results_root))
cat("======================================================================\n")