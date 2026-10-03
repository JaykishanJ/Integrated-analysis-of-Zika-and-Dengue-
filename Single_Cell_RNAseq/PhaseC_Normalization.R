# ==============================================================================
# PHASE C - NORMALIZATION & MODELING 
# ==============================================================================
# ERCC-based normalization, HVG selection, cell-cycle scoring,
# scaling, and infection-state stratification.
#
# Inputs (from Phases A-B):
# results/phaseB/objects/Phase_B_counts_active.rds
# results/phaseB/objects/Phase_B_cell_qc.rds
# results/phaseA/objects/Phase_A_gene_map.rds
# results/phaseB/objects/Phase_B_gene_sets.rds
#
# Outputs (to Phase D/E):
# results/phaseC/objects/Phase_C_log_matrix.rds
# results/phaseC/objects/Phase_C_scaled_matrix.rds
# results/phaseC/objects/Phase_C_cell_qc.rds
# results/phaseC/objects/Phase_C_hvgs.rds
#
# NOTE: all figures are written as vector SVG only (results/phaseC/figures/*.svg)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  source("E:/Zika/scripts/inject_themes.R")
  library(scales)
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
       "\nRestart R and source this script in a clean session.")
}
rm(.verbs, .owner, .masked)

# ==============================================================================
# Configuration and output layout
# ==============================================================================
data_dir <- "E:/Zika"
setwd(data_dir)

results_root <- file.path(data_dir, "results", "phaseC")
dir_fig <- file.path(results_root, "figures")
dir_tab <- file.path(results_root, "tables")
dir_obj <- file.path(results_root, "objects")
for (d in c(dir_fig, dir_tab, dir_obj))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# show_and_save(): Apply publication theme, print, and save as SVG.
# --------------------------------------------------------------------------
show_and_save <- function(p, file, width = 8, height = 5, base_size = 10) {
  p <- p + theme_pub(base_size = base_size)
  print(p)
  ggsave(file.path(dir_fig, paste0(file, ".svg")), p, width = width, height = height,
         bg = "transparent", device = svglite::svglite)
  ggsave(fig(paste0(file, ".png")), p, width = width, height = height,
         dpi = 600, bg = "white")
  cat(sprintf(" [plot] %s.svg and .png\n", file))
  invisible(p)
}
write_tab <- function(x, name) {
  write.csv(x, file.path(dir_tab, name), row.names = FALSE)
  cat(sprintf(" [table] %s\n", name))
}

fig <- function(n) file.path(dir_fig, n)
tab <- function(n) file.path(dir_tab, n)
obj <- function(n) file.path(dir_obj, n)

min_pos <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) min(x) else 1 }

# Resolve Phase A outputs from multiple possible locations
find_file <- function(bn, subdirs = c("results/phaseB/objects", "results/phaseB/tables",
                                      "results/phaseA/objects", "results/phaseA/tables", ".")) {
  for (s in subdirs) {
    p <- file.path(data_dir, s, bn)
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NA_character_
}

cat("======================================================================\n")
cat("PHASE C - NORMALIZATION & MODELING \n")
cat("======================================================================\n")

# ------------------------------------------------------------------------------
# Load Phase B Outputs
# ------------------------------------------------------------------------------
p_counts <- find_file("Phase_B_counts_active.rds")
p_qc <- find_file("Phase_B_cell_qc.rds")
p_map <- find_file("Phase_A_gene_map.rds")

if (is.na(p_counts) || is.na(p_qc))
  stop("Missing Phase B outputs. Need Phase_B_counts_active.rds and Phase_B_cell_qc.rds")

mat_active <- readRDS(p_counts)
qc <- readRDS(p_qc)
gene_map <- if (!is.na(p_map)) readRDS(p_map) else NULL

if (is.null(gene_map))
  warning("Phase_A_gene_map.rds not found - symbol mapping and cell-cycle scoring will be disabled.")

ensg_col <- if (!is.null(gene_map)) intersect(c("ENSG", "ENSG_Original", "ENSG_Clean"), names(gene_map))[1] else NA_character_
sym_of <- function(ids) {
  if (is.null(gene_map) || is.na(ensg_col)) return(ids)
  s <- gene_map$Symbol[match(ids, gene_map[[ensg_col]])]
  ifelse(is.na(s), ids, s)
}

# Align QC table to matrix columns
if (!all(colnames(mat_active) %in% qc$Sample))
  stop("Cell IDs in the active matrix are not all present in the QC table.")
n_qc_all <- nrow(qc)
qc <- qc[match(colnames(mat_active), qc$Sample), , drop = FALSE]
rownames(qc) <- NULL
stopifnot(identical(as.character(qc$Sample), colnames(mat_active)))

cat(sprintf("Loaded active matrix: %d features x %d cells\n", nrow(mat_active), ncol(mat_active)))
cat(sprintf("Aligned QC table to matrix columns: %d cells (from %d QC rows)\n",
            nrow(qc), n_qc_all))

# ==============================================================================
# NORMALIZATION & HVG SELECTION
# ==============================================================================
cat("\nLog-normalization & HVG selection\n")
# Log-CPM normalization using ERCC spike-in size factors (scran::computeSpikeFactors).
# Library-size normalization is biased when virus dominates the library.

lib_sizes <- colSums(mat_active)
lib_sizes[lib_sizes == 0] <- 1

compute_log_cpm <- function(m, lib, chunk = 5000) {
  ng <- nrow(m); out <- matrix(0, nrow = ng, ncol = ncol(m))
  rownames(out) <- rownames(m); colnames(out) <- colnames(m)
  for (ch in split(seq_len(ng), ceiling(seq_len(ng) / chunk))) {
    out[ch, ] <- log2(t(t(m[ch, , drop = FALSE]) / lib) + 1)
  }
  out
}

ercc_sf <- NULL
if (requireNamespace("scran", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  try({
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat_active))
    SingleCellExperiment::altExp(sce, "ERCC") <-
      SingleCellExperiment::SummarizedExperiment(list(counts = mat_active[grep("^ERCC-", rownames(mat_active)), ]))
    sce <- scran::computeSpikeFactors(sce, "ERCC")
    ercc_sf <- SingleCellExperiment::sizeFactors(sce)
  }, silent = TRUE)
}

if (is.null(ercc_sf)) {
  ercc_sf <- qc$ercc_counts
  ercc_sf[ercc_sf == 0] <- 1
  ercc_sf <- ercc_sf / mean(ercc_sf)
}

cat(sprintf(" -> ERCC spike-in size factors: range %.2f - %.2f, median %.2f\n",
            min(ercc_sf, na.rm = TRUE), max(ercc_sf, na.rm = TRUE), median(ercc_sf, na.rm = TRUE)))

mat_log <- compute_log_cpm(mat_active, ercc_sf)

# Mean-variance trend & HVGs
if (requireNamespace("matrixStats", quietly = TRUE)) {
  gene_vars <- matrixStats::rowVars(mat_log)
} else {
  gene_vars <- apply(mat_log, 1, var)
}
gene_means <- rowMeans(mat_log)
names(gene_vars) <- rownames(mat_log)

hvg_n <- 2000
hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(min(hvg_n, nrow(mat_log)))]
hvg_ids <- rownames(mat_log)[hvg_idx]

hvg_df <- data.frame(
  ensembl = rownames(mat_log),
  symbol = sym_of(rownames(mat_log)),
  mean = gene_means,
  var = gene_vars,
  is_hvg = seq_len(nrow(mat_log)) %in% hvg_idx,
  stringsAsFactors = FALSE
)
write_tab(hvg_df, "HVG_Statistics.csv")

# Confirm ISG and ER-stress markers in HVGs
isg_markers <- c("ISG15", "MX1", "IFIT1", "RSAD2", "IFI44L")
er_markers <- c("HSPA5", "DDIT3", "ATF4", "XBP1")
target_symbols <- sym_of(hvg_ids)
cat(sprintf(" -> ISG markers in top %d HVGs: %s\n", hvg_n,
            paste(intersect(isg_markers, target_symbols), collapse = ", ")))
cat(sprintf(" -> ER-stress markers in top %d HVGs: %s\n", hvg_n,
            paste(intersect(er_markers, target_symbols), collapse = ", ")))

p6 <- ggplot(hvg_df, aes(x = mean, y = var, colour = is_hvg)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_colour_manual(values = c("FALSE" = "#CCCCCC", "TRUE" = "#D55E00"), name = "HVG") +
  labs(title = "Mean-variance trend and HVG selection",
       subtitle = sprintf("Top %d HVGs highlighted (ERCC-normalized Log-CPM)", hvg_n),
       x = "Mean log-CPM", y = "Variance log-CPM")
show_and_save(p6, "Mean_Variance_Trend", width = 8, height = 5)

# ==============================================================================
# CELL-CYCLE SCORING
# ==============================================================================
cat("\nCell-cycle scoring\n")
# Background-corrected module scores (Tirosh et al. 2016 / Seurat CellCycleScoring)

s_genes_raw <- c("MCM5", "PCNA", "ORC1", "DTL", "GINS2", "MCM4", "FEN1", "MCM2",
                 "RRM1", "UNG", "EXO1", "TIPIN", "BLM", "CASP8AP2", "USP1", "SLBP",
                 "CHAF1B", "BRIP1", "E2F8")
g2m_genes_raw <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
                   "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
                   "TTC35", "G2E3", "SMARCA4", "CCNB2", "CKAP2L")

map_to_ensg <- function(symbols) {
  if (is.null(gene_map) || is.na(ensg_col)) return(character(0))
  matched <- gene_map[[ensg_col]][match(toupper(symbols), toupper(gene_map$Symbol))]
  matched[!is.na(matched)]
}

s_ensg <- intersect(map_to_ensg(s_genes_raw), rownames(mat_log))
g2m_ensg <- intersect(map_to_ensg(g2m_genes_raw), rownames(mat_log))

failed_s <- setdiff(map_to_ensg(s_genes_raw), rownames(mat_log))
failed_g2m <- setdiff(map_to_ensg(g2m_genes_raw), rownames(mat_log))
cat(sprintf(" -> Cell cycle mapping: %d S genes, %d G2M genes mapped.\n", length(s_ensg), length(g2m_ensg)))
if (length(failed_s) > 0 || length(failed_g2m) > 0) {
  cat(sprintf(" -> Unmapped markers: %s\n",
              paste(sym_of(c(failed_s, failed_g2m)), collapse = ", ")))
}

# Background-corrected module scores
module_score <- function(mat, features, ctrl = 100, nbin = 24, seed = 42) {
  features <- intersect(features, rownames(mat))
  if (length(features) == 0) return(setNames(rep(0, ncol(mat)), colnames(mat)))
  set.seed(seed)
  avg <- rowMeans(mat)
  bins <- ceiling(rank(avg, ties.method = "first") / length(avg) * nbin)
  bins[bins < 1] <- 1
  names(bins) <- rownames(mat)
  ctrl.use <- character(0)
  for (b in unique(bins[features])) {
    pool <- setdiff(names(bins)[bins == b], features)
    if (length(pool) == 0) next
    size <- min(ctrl * sum(bins[features] == b), length(pool))
    ctrl.use <- c(ctrl.use, sample(pool, size = size, replace = FALSE))
  }
  ctrl.use <- unique(ctrl.use)
  colMeans(mat[features, , drop = FALSE]) - colMeans(mat[ctrl.use, , drop = FALSE])
}

qc <- qc %>% mutate(
  S_score = as.numeric(module_score(mat_log, s_ensg)),
  G2M_score = as.numeric(module_score(mat_log, g2m_ensg)),
  Phase = case_when(
    S_score < 0 & G2M_score < 0 ~ "G1",
    S_score >= G2M_score ~ "S",
    TRUE ~ "G2M"
  )
)

write_tab(qc %>% select(Sample, virus, time_h, moi, S_score, G2M_score, Phase),
          "Cell_Cycle_Scores.csv")

p7a <- qc %>% filter(!is.na(time_h)) %>%
  count(time_h, virus, Phase) %>%
  group_by(time_h, virus) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = factor(time_h), y = prop, fill = Phase)) +
  geom_col(position = "fill", colour = "black", width = 0.7) +
  facet_wrap(~ virus) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Set3", name = "Phase") +
  labs(title = "Cell-cycle phase distribution by timepoint",
       x = "Hours post-infection", y = "Proportion of cells")
show_and_save(p7a, "Phase_Distribution", width = 8, height = 5)

p7b <- ggplot(qc, aes(x = S_score, y = G2M_score, colour = Phase)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_brewer(palette = "Set1", name = "Phase") +
  labs(title = "S-score vs G2M-score per cell",
       x = "S-phase score", y = "G2M-phase score")
show_and_save(p7b, "Score_Scatter", width = 7, height = 5)

# Phase vs cluster overlap check (requires clustering, done here as a check)
# If phase correlates with clusters, S/G2M scores should be regressed in .
# The actual regression happens in Phase D () after clustering is available.
cat(" -> NOTE: Plate batch effect regression and cell-cycle regression happen in Phase D (),\n")
cat(" after clustering is available to assess overlap.\n")

# ==============================================================================
# SCALING
# ==============================================================================
cat("\nScale features\n")

# Subset to HVGs and scale (mean 0, var 1 per gene)
# vars.to.regress = percent.mt only. Plate is NOT regressed here - it is handled
# in Phase D via Harmony if PC-plate correlation is strong.
mat_hvg <- mat_log[hvg_ids, , drop = FALSE]
mat_scaled <- t(scale(t(mat_hvg)))
mat_scaled[is.na(mat_scaled)] <- 0

# Compute baseline PCA for confounder evaluation
pca <- prcomp(t(mat_scaled))
ve_base <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat(sprintf(" -> Baseline PCA on HVGs: PC1 explains %.1f%% PC2 explains %.1f%%\n",
            ve_base[1], ve_base[2]))

pca_df <- data.frame(
  Sample = colnames(mat_scaled),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
) %>% left_join(qc, by = "Sample")

p8 <- ggplot(pca_df %>% filter(!is.na(time_h)), aes(x = PC1, y = PC2, colour = Phase)) +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_grid(virus ~ time_h, labeller = label_both) +
  scale_colour_brewer(palette = "Set1", name = "Phase") +
  labs(title = "PCA colored by cell-cycle phase per timepoint",
       subtitle = "Scaled HVGs (mean 0, var 1), regressed: percent.mt only",
       x = "PC1", y = "PC2")
show_and_save(p8, "PCA_CellCycle", width = 10, height = 6)

# ==============================================================================
# STRATIFY INFECTION STATUS
# ==============================================================================
cat("\nStratify Infection Status\n")

qc <- qc %>% mutate(n_virus_molecules = replace_na(n_virus_molecules, 0))

bg_thr <- qc %>%
  filter(moi == 0) %>%
  group_by(virus) %>%
  summarise(
    hop_threshold_95 = quantile(n_virus_molecules, 0.95, na.rm = TRUE),
    hop_threshold_99 = quantile(n_virus_molecules, 0.99, na.rm = TRUE),
    hop_threshold_100 = max(n_virus_molecules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hop_threshold = pmax(hop_threshold_99, 1)
  )

cat(" -> Evaluating Infection Threshold Stability (95th vs 99th vs 100th percentile):\n")

# Sensitivity table: for each virus x threshold, classify every exposed cell
sens_list <- list()
for (i in seq_len(nrow(bg_thr))) {
  v <- bg_thr$virus[i]
  v_cells <- qc %>% filter(virus == v, moi > 0)
  
  for (qt in c("hop_threshold_95", "hop_threshold_99", "hop_threshold_100")) {
    thr_val <- bg_thr[[qt]][i]
    inf_n <- sum(v_cells$n_virus_molecules > pmax(thr_val, 1))
    bys_n <- sum(v_cells$n_virus_molecules <= pmax(thr_val, 1))
    sens_list[[length(sens_list) + 1]] <- data.frame(
      virus = v,
      quantile = sub("hop_threshold_", "", qt),
      threshold = thr_val,
      n_infected = inf_n,
      n_bystander = bys_n,
      pct_infected = 100 * inf_n / (inf_n + bys_n),
      stringsAsFactors = FALSE
    )
  }
  
  cat(sprintf(" %s controls max: 95th=%.1f, 99th=%.1f, 100th=%.1f viral molecules\n",
              v, bg_thr$hop_threshold_95[i], bg_thr$hop_threshold_99[i], bg_thr$hop_threshold_100[i]))
}
sens_tab <- bind_rows(sens_list)
write_tab(sens_tab, "Threshold_Sensitivity.csv")

# Visualise threshold sensitivity
sens_plot <- sens_tab %>%
  mutate(quantile = factor(quantile, levels = c("95", "99", "100"))) %>%
  pivot_longer(c(n_infected, n_bystander), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_infected = "Infected", n_bystander = "Bystander"))

p_sens <- ggplot(sens_plot, aes(x = quantile, y = n, fill = status)) +
  geom_col(position = "fill", colour = "black", width = 0.7) +
  facet_wrap(~ virus) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c(Infected = "#D55E00", Bystander = "#E69F00"), name = NULL) +
  labs(title = "Infection classification sensitivity to threshold",
       subtitle = "Proportion of exposed cells called Infected at 95th, 99th, and 100th percentile of the MOI 0 control distribution per virus.",
       x = "Threshold (percentile of MOI 0 controls)", y = "Proportion")
show_and_save(p_sens, "Threshold_Sensitivity", width = 7, height = 4)

# Final classification using 99th percentile
qc <- qc %>%
  left_join(bg_thr, by = "virus") %>%
  mutate(
    infection_status = case_when(
      moi == 0 ~ "Unexposed",
      n_virus_molecules > hop_threshold ~ "Infected",
      TRUE ~ "Bystander"
    )
  )

# Ensure time_h is numeric (left_join can coerce it)
qc$time_h <- as.numeric(qc$time_h)

cat(" -> Final infection thresholds (99th percentile of MOI 0 controls):\n")
for (i in seq_len(nrow(bg_thr)))
  cat(sprintf(" %s: %.1f viral molecules\n", bg_thr$virus[i], bg_thr$hop_threshold[i]))
cat(sprintf(" -> Cell classification summary: %d Infected, %d Bystander, %d Unexposed\n",
            sum(qc$infection_status == "Infected"), sum(qc$infection_status == "Bystander"),
            sum(qc$infection_status == "Unexposed")))

# Plate x virus x state table (feasibility check for within-plate Wilcoxon)
plate_state_tab <- qc %>%
  filter(moi > 0, !is.na(infection_status)) %>%
  count(batch, virus, infection_status) %>%
  group_by(batch, virus) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  filter(infection_status %in% c("Infected", "Bystander")) %>%
  mutate(feasible = n >= 20)
write_tab(plate_state_tab, "Plate_State_Counts.csv")

cat(" Plate feasibility for within-plate Wilcoxon (>= 20 cells per group):\n")
for (b in unique(plate_state_tab$batch)) {
  sub <- plate_state_tab %>% filter(batch == b)
  inf_ok <- any(sub$infection_status == "Infected" & sub$feasible)
  bys_ok <- any(sub$infection_status == "Bystander" & sub$feasible)
  cat(sprintf(" Plate %s: Infected feasible=%s Bystander feasible=%s\n",
              b, inf_ok, bys_ok))
}

p9a <- ggplot(qc, aes(x = n_virus_molecules + 1, fill = infection_status)) +
  geom_histogram(bins = 50, colour = "black", linewidth = 0.2) +
  geom_vline(data = bg_thr, aes(xintercept = hop_threshold + 1),
             linetype = "dashed", colour = "red") +
  facet_wrap(~ virus, scales = "free_y", ncol = 1) +
  scale_x_log10() +
  scale_fill_manual(values = c(Unexposed = "grey75", Bystander = "#0072B2",
                               Infected = "#D55E00"), name = NULL) +
  labs(title = "Viral molecule distribution and infection threshold",
       subtitle = "Dashed red = index-hopping background threshold per virus",
       x = "Viral molecules + 1 (log10)", y = "Cells")
show_and_save(p9a, "Viral_Histogram", width = 9, height = 7)

p9b <- qc %>% filter(!is.na(time_h)) %>%
  count(time_h, virus, infection_status) %>%
  group_by(time_h, virus) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = factor(time_h), y = prop, fill = infection_status)) +
  geom_col(position = "fill", colour = "black", width = 0.7) +
  facet_wrap(~ virus) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c(Unexposed = "grey75", Bystander = "#0072B2",
                               Infected = "#D55E00"), name = NULL) +
  labs(title = "Infection status proportions by timepoint",
       subtitle = "Not infection kinetics: cells were selected by viral qPCR to span the vRNA range.",
       x = "Hours post-infection", y = "Proportion of cells")
show_and_save(p9b, "Bystander_Infected_Proportions", width = 8, height = 5)

# ==============================================================================
# Persist Objects for Phase D
# ==============================================================================
saveRDS(mat_scaled, obj("Phase_C_scaled_matrix.rds"))
saveRDS(qc, obj("Phase_C_cell_qc.rds"))
saveRDS(hvg_ids, obj("Phase_C_hvgs.rds"))
saveRDS(mat_log, obj("Phase_C_log_matrix.rds"))

cat("\n======================================================================\n")
cat("PHASE C COMPLETE\n")
cat(sprintf(" Normalization: ERCC spike-in size factors (scran::computeSpikeFactors)\n"))
cat(sprintf(" Scaled Matrix : %d features x %d cells\n", nrow(mat_scaled), ncol(mat_scaled)))
cat(sprintf(" Objects : %s\n", dir_obj))
cat(sprintf(" Figures/tables: %s\n", results_root))
cat("======================================================================\n")