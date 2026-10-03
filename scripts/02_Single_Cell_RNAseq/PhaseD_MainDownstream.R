# ==============================================================================
# PHASE D - MAIN DOWNSTREAM ANALYSIS 
# ==============================================================================
# PCA, clustering, within-plate Wilcoxon DE, pseudobulk (feasibility-gated),
# and continuous vRNA association.
#
# Inputs (from Phases A-C):
# results/phaseC/objects/Phase_C_log_matrix.rds
# results/phaseC/objects/Phase_C_scaled_matrix.rds
# results/phaseC/objects/Phase_C_cell_qc.rds
# results/phaseC/objects/Phase_C_hvgs.rds
# results/phaseB/objects/Phase_B_counts_active.rds
# results/phaseA/objects/Phase_A_gene_map.rds
# results/phaseB/objects/Phase_B_gene_sets.rds
#
# Outputs:
# results/phaseD/objects/Phase_D_clustering.rds -> Phase E
# results/phaseD/objects/Phase_D_ZIKV_Sig_Matrix.rds -> heatmaps
# results/phaseD/objects/Phase_D_DENV_Sig_Matrix.rds -> heatmaps
# results/phaseD/objects/Phase_D_Common_Sig_Matrix.rds -> heatmaps
# results/phaseD/tables/*.csv
# results/phaseD/tables/DE_*.csv -> Phase E
# results/phaseD/tables/Pseudobulk_*.csv -> Phase E
# results/phaseD/tables/vRNA_*.csv -> Phase E
# results/phaseD/tables/*.csv -> Phase E (timepoint DE + DEG lists)
# results/phaseD/figures/Volcano_*.svg
# results/phaseD/figures/MA_Plot_*.svg
# results/phaseD/supplementary/ -> Phase D supplementary
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  source("./scripts/inject_themes.R")
  library(scales)
  library(Matrix)
  library(igraph)
  library(FNN)
  library(Rtsne)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
  library(svglite)
})

set.seed(42)

# Fail loudly on masked dplyr
.verbs <- c("select", "filter", "rename", "count", "slice", "desc", "first")
.owner <- vapply(.verbs, function(f) find(f)[1], character(1))
.masked <- .owner[!.owner %in% c("package:dplyr", "package:tidyr")]
if (length(.masked))
  stop("dplyr verbs are masked: ",
       paste(sprintf("%s() <- %s", names(.masked), .masked), collapse = ", "),
       "\nRestart R and source in a clean session.")
rm(.verbs, .owner, .masked)

# uwot supplies UMAP; fall back to t-SNE if it fails
.uwot_err <- NA_character_
HAS_UWOT <- tryCatch({ suppressMessages(library(uwot)); TRUE },
                     error = function(e) { .uwot_err <<- conditionMessage(e); FALSE })

# Preflight: clusterProfiler needed by Phase E enrichment
.cp_ok <- suppressWarnings(tryCatch({
  suppressMessages(loadNamespace("clusterProfiler"))
  TRUE
}, error = function(e) conditionMessage(e)))
if (!isTRUE(.cp_ok)) {
  stop(sprintf(" needs clusterProfiler: %s\n",
               "Install with: if (!requireNamespace('BiocManager')) install.packages('BiocManager')\n",
               "BiocManager::install('clusterProfiler')\n",
               "Then restart R and re-source."), call. = FALSE)
}

# Preflight: svglite needed for SVG outputs
.svg_ok <- suppressWarnings(tryCatch({
  suppressMessages(loadNamespace("svglite"))
  TRUE
}, error = function(e) conditionMessage(e)))
if (!isTRUE(.svg_ok)) {
  stop(sprintf("SVG output needs svglite: %s\n",
               "Install with: install.packages('svglite')\n",
               "Then restart R and re-source."), call. = FALSE)
}

# ==============================================================================
# Configuration and output layout
# ==============================================================================
data_dir <- "."
setwd(data_dir)

results_root <- file.path(data_dir, "results", "phaseD")
dir_fig <- file.path(results_root, "figures")
dir_tab <- file.path(results_root, "tables")
dir_obj <- file.path(results_root, "objects")
dir_supp <- file.path(results_root, "supplementary")
dir_hm_cache <- file.path(results_root, "cache", "harmony")
for (d in c(dir_fig, dir_tab, dir_obj, dir_supp, dir_hm_cache))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Publication-ready plotting helpers
# - No grey panel/grid background
# - Larger, journal-friendly typography
# - Transparent SVG background
# - Consistent line weights and spacing
PUB_BASE_SIZE <- 13

theme_pub <- function(base_size = PUB_BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = base_size + 3, face = "bold",
                                colour = "#111111", hjust = 0, margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size - 1, colour = "#444444",
                                   hjust = 0, margin = margin(b = 8)),
      axis.title = element_text(size = base_size, face = "bold", colour = "#222222"),
      axis.text = element_text(size = base_size - 1, colour = "#333333"),
      axis.text.x = element_text(margin = margin(t = 4)),
      axis.text.y = element_text(margin = margin(r = 4)),
      axis.line = element_line(linewidth = 0.55, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.45, colour = "#333333"),
      axis.ticks.length = grid::unit(2.5, "pt"),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      panel.grid = element_blank(),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.position = "bottom",
      legend.margin = margin(t = 2),
      plot.margin = margin(10, 12, 10, 12)
    )
}

show_and_save <- function(p, file, width = 8, height = 5) {
  p <- p + theme_pub()
  print(p)
  ggsave(file.path(dir_fig, paste0(file, ".svg")), p,
         width = width, height = height, bg = "transparent")
  ggsave(file.path(dir_fig, paste0(file, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
  cat(sprintf(" [plot] %s.svg and .png\n", file))
  invisible(p)
}

supp_and_save <- function(p, file, width = 8, height = 5) {
  p <- p + theme_pub()
  print(p)
  ggsave(file.path(dir_supp, paste0(file, ".svg")), p,
         width = width, height = height, bg = "transparent")
  cat(sprintf(" [supp] %s.svg\n", file))
  invisible(p)
}
save_pheatmap <- function(ph, file, width = 8, height = 6) {
  svg(file.path(dir_fig, paste0(file, ".svg")), width = width, height = height,
      bg = "transparent")
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
  png(file.path(dir_fig, paste0(file, ".png")), width = width, height = height, units = "in", res = 600)
  grid::grid.draw(ph$gtable)
  dev.off()
  cat(sprintf(" [plot] %s.svg and .png\n", file))
}
write_tab <- function(x, name) {
  write.csv(x, file.path(dir_tab, name), row.names = FALSE)
  cat(sprintf(" [table] %s\n", name))
}

fig <- function(n) file.path(dir_fig, n)
tab <- function(n) file.path(dir_tab, n)
obj <- function(n) file.path(dir_obj, n)

# Harmony caching
hm_cache_path <- file.path(dir_hm_cache, "Phase_D_harmony_adjusted.rds")

find_file <- function(bn, subdirs = c("results/phaseC/objects","results/phaseB/objects",
                                      "results/phaseA/objects","results/phaseB/tables",
                                      "results/phaseD/objects","results/phaseD/tables",
                                      "results/phaseD_1/objects","results/phaseD_1/tables",
                                      ".")) {
  for (s in subdirs) {
    p <- file.path(data_dir, s, bn)
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NA_character_
}

cat("======================================================================\n")
cat("PHASE D - MAIN DOWNSTREAM ANALYSIS \n")
cat(sprintf("Embedding engine: %s\n",
            if (HAS_UWOT) "uwot (UMAP)"
            else sprintf("Rtsne (t-SNE) - uwot did not load: %s", .uwot_err)))
cat("======================================================================\n")

# ==============================================================================
# Load inputs from Phase C/B/A
# ==============================================================================
mat_scaled <- readRDS(find_file("Phase_C_scaled_matrix.rds"))
qc <- readRDS(find_file("Phase_C_cell_qc.rds"))
hvg_ids <- readRDS(find_file("Phase_C_hvgs.rds"))
mat_counts <- readRDS(find_file("Phase_B_counts_active.rds"))
gene_map <- readRDS(find_file("Phase_A_gene_map.rds"))
gene_sets <- tryCatch(readRDS(find_file("Phase_B_gene_sets.rds")), error = function(e) NULL)
mat_log <- tryCatch(readRDS(find_file("Phase_C_log_matrix.rds")), error = function(e) NULL)
if (is.null(mat_log)) {
  cat(" WARNING: Phase_C_log_matrix.rds not found; module scores will use mat_scaled.\n")
  mat_log <- mat_scaled
}

ensg_col <- intersect(c("ENSG", "ENSG_Original", "ENSG_Clean"), names(gene_map))[1]
if (is.na(ensg_col) || !"Symbol" %in% names(gene_map))
  stop("Phase_A_gene_map.rds has no recognisable ID/Symbol columns.")
sym_of <- function(ids) {
  s <- gene_map$Symbol[match(ids, gene_map[[ensg_col]])]
  ifelse(is.na(s), ids, s)
}

stopifnot(identical(colnames(mat_scaled), as.character(qc$Sample)))
stopifnot(identical(colnames(mat_counts), as.character(qc$Sample)))
qc$hpi <- qc$time_h
qc$log_viral <- log10(qc$n_virus_molecules + 1)

# Guard: MOI 0 controls must never be labelled Bystander
n_moi0_bys <- sum(qc$moi == 0 & qc$infection_status == "Bystander", na.rm = TRUE)
if (n_moi0_bys > 0)
  stop(sprintf("%d MOI 0 cells are labelled 'Bystander'. Unexposed controls are not bystanders.", n_moi0_bys))
cat(sprintf("PASS: %d unexposed cells held out of Bystander group.\n", sum(qc$moi == 0)))

# ==============================================================================
# PCA + PLATE CHECK + (OPTIONAL) HARMONY INTEGRATION
# ==============================================================================
cat("\nPCA + plate check + integration\n")

n_pc_max <- 50
pca <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE, rank. = n_pc_max)
ve <- 100 * pca$sdev[seq_len(n_pc_max)]^2 / sum(pca$sdev^2)

## PCA elbow
knee_point <- function(y) {
  n <- length(y); x <- seq_len(n)
  num <- abs((y[n] - y[1]) * x - (x[n] - x[1]) * y + x[n] * y[1] - y[n] * x[1])
  which.max(num / sqrt((y[n] - y[1])^2 + (x[n] - x[1])^2))
}
knee_ref <- knee_point(ve)
rel_drop <- c(-diff(ve)) / head(ve, -1)
FLAT_TOL <- 0.04; FLAT_RUN <- 3
flat_onset <- NA_integer_
for (i in seq_len(length(rel_drop) - FLAT_RUN + 1))
  if (all(rel_drop[i:(i + FLAT_RUN - 1)] < FLAT_TOL)) { flat_onset <- i; break }
n_pc_auto <- if (is.na(flat_onset)) knee_ref else flat_onset

cat(sprintf(" -> #PCs: chord knee = PC%d; flat tail = PC%s; auto-select = %d PCs (%.1f%% variance)\n",
            knee_ref, ifelse(is.na(flat_onset), "NA", flat_onset), n_pc_auto, sum(ve[seq_len(n_pc_auto)])))

## Print variance table so you can pick
cat("\n PC var% cum%\n")
cat(" ---- ----- ------\n")
for (k in seq_len(n_pc_max)) {
  cum <- sum(ve[seq_len(k)])
  marker <- if (k == n_pc_auto) " <-- auto" else ""
  cat(sprintf(" PC%-3d %5.2f %6.1f%s\n", k, ve[k], cum, marker))
  if (k %% 10 == 0 && k < n_pc_max)
    cat(" ---- ----- ------\n")
}

## Interactive override via environment variable (set before sourcing) or auto
## For pasting: set N_PC_OVERRIDE in your R session first, e.g.:
## N_PC_OVERRIDE <- 15
## source("scripts/phaseD/PhaseD_MainDownstream.R")
user_choice <- Sys.getenv("N_PC_OVERRIDE", unset = "")
if (nzchar(user_choice)) {
  n_pc <- as.integer(user_choice)
  if (is.na(n_pc) || n_pc < 1 || n_pc > n_pc_max)
    stop(sprintf("N_PC_OVERRIDE='%s' invalid. Must be between 1 and %d.", user_choice, n_pc_max))
  cat(sprintf(" -> N_PC_OVERRIDE set: using %d PCs (%.1f%% cumulative variance)\n", n_pc, sum(ve[seq_len(n_pc)])))
} else if (exists("N_PC_OVERRIDE", envir = .GlobalEnv)) {
  n_pc <- as.integer(get("N_PC_OVERRIDE", envir = .GlobalEnv))
  if (is.na(n_pc) || n_pc < 1 || n_pc > n_pc_max)
    stop(sprintf("N_PC_OVERRIDE=%s invalid. Must be between 1 and %d.", n_pc, n_pc_max))
  cat(sprintf(" -> N_PC_OVERRIDE set: using %d PCs (%.1f%% cumulative variance)\n", n_pc, sum(ve[seq_len(n_pc)])))
} else {
  n_pc <- n_pc_auto
  cat(sprintf(" -> Using %d PCs (%.1f%% cumulative variance)\n", n_pc, sum(ve[seq_len(n_pc)])))
  cat(sprintf(" (To override: N_PC_OVERRIDE <- %d; source(...))\n", n_pc))
}

## Build elbow data frame (must come before plot)
elbow_df <- data.frame(PC = seq_len(n_pc_max), var = ve)

## Threshold: variance < 0.25% = noise (visualise where PCs stop explaining anything)
NOISE_THR <- 0.25
noise_ix <- which(ve < NOISE_THR)[1]
if (is.na(noise_ix)) noise_ix <- n_pc_max + 1

p_elbow <- ggplot(elbow_df, aes(x = PC, y = var)) +
  # Noise region shading (where variance < threshold â†’ effectively zero)
  annotate("rect", xmin = noise_ix - 0.5, xmax = n_pc_max + 0.5, ymin = -Inf, ymax = Inf,
           alpha = 0.06, fill = "#D55E00") +
  # Retained region shading (up to cutoff)
  annotate("rect", xmin = 0.5, xmax = n_pc + 0.5, ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "#0072B2") +
  # Variance line + retained/dropped colouring
  geom_line(colour = "#808080", linewidth = 0.6) +
  geom_point(data = elbow_df, aes(colour = PC <= n_pc), size = 2.2, alpha = 0.9) +
  # Zero-variance threshold line
  geom_hline(yintercept = NOISE_THR, linetype = "dashed", colour = "#D55E00", linewidth = 0.5) +
  annotate("text", x = n_pc_max - 1, y = NOISE_THR + 0.15, hjust = 1, colour = "#D55E00", size = 2.8,
           label = sprintf("noise threshold = %.2f%%", NOISE_THR)) +
  # Cutoff marker
  geom_vline(xintercept = n_pc, linetype = "dashed", colour = "#0072B2", linewidth = 0.8) +
  annotate("text", x = n_pc + 1.5, y = max(ve) * 0.72, hjust = 0, colour = "#0072B2", size = 3.5,
           label = sprintf("PC1-%d retained\n(%.1f%% cum.)", n_pc, sum(ve[seq_len(n_pc)]))) +
  # Label the last retained PC
  annotate("text", x = n_pc, y = ve[n_pc] + 0.3, colour = "#0072B2", size = 3,
           label = sprintf("PC%d\n%.2f%%", n_pc, ve[n_pc]), vjust = -0.5) +
  # Label where noise begins
  if (!is.na(noise_ix) && noise_ix > n_pc + 2)
    annotate("text", x = noise_ix, y = NOISE_THR * 0.7, colour = "#D55E00", size = 2.8,
             label = sprintf("noise from\nPC%d (%.2f%%)", noise_ix, ve[noise_ix]), hjust = 0) +
  # Chord knee reference
  geom_vline(xintercept = knee_ref, linetype = "dotted", colour = "grey55", linewidth = 0.6) +
  annotate("text", x = knee_ref + 0.4, y = max(ve) * 0.90, hjust = 0, size = 3, colour = "#737373",
           label = sprintf("chord knee = PC%d", knee_ref)) +
  scale_colour_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#B0B0B0"),
                      labels = c("Retained", "Noise"), name = NULL) +
  scale_x_continuous(breaks = seq(0, n_pc_max, 5), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_pub() + theme(legend.position = "bottom") +
  labs(title = sprintf("PCA elbow â€” retaining PC1â€“%d (%.1f%% variance)", n_pc, sum(ve[seq_len(n_pc)])),
       subtitle = sprintf("Auto = PC%d | noise threshold = %.2f%% | chord knee = PC%d",
                          n_pc_auto, NOISE_THR, knee_ref),
       x = "Principal component", y = "% variance explained")
show_and_save(p_elbow, "PCA_Elbow", width = 11, height = 5)
write_tab(elbow_df, "PCA_Variance_Explained.csv")

pcs <- pca$x[, seq_len(n_pc), drop = FALSE]

## NEW: PC-plate correlation check + Harmony integration
# This is the critical addition that was missing in the previous workflow.
# If plate (batch) drives PC1/PC2, clusters will track plates rather than
# biology. Harmony removes this effect before clustering.
# Plate is confounded with condition (each plate = one virus/MOI/hpi combo),
# so we must check BEFORE assuming clusters are biological.

plate_num <- as.numeric(factor(qc$batch))

pc_cov <- data.frame(
  PC = paste0("PC", seq_len(min(50, n_pc_max))),
  var_pct = 100 * pca$sdev[seq_len(min(50, n_pc_max))]^2 / sum(pca$sdev^2)
)

corr_mat <- sapply(seq_len(min(50, n_pc_max)), function(k) {
  pc_k <- pca$x[, k]
  c(
    plate = suppressWarnings(cor(pc_k, plate_num, method = "spearman")),
    nCount = suppressWarnings(cor(pc_k, qc$host_counts, method = "spearman")),
    log_viral = suppressWarnings(cor(pc_k, qc$log_viral, method = "spearman")),
    S_score = suppressWarnings(cor(pc_k, qc$S_score, method = "spearman"))
  )
})
corr_mat <- t(corr_mat)
pc_cov <- cbind(pc_cov, as.data.frame(corr_mat))
pc_cov$plate_abs <- abs(pc_cov$plate)
pc_cov$nCount_abs <- abs(pc_cov$nCount)
pc_cov$log_viral_abs <- abs(pc_cov$log_viral)
pc_cov$S_score_abs <- abs(pc_cov$S_score)

pc_cov$plate_driven <- pc_cov$plate_abs > 0.3
n_plate_driven <- sum(pc_cov$plate_driven)
cat(sprintf(" -> PCs with |r(plate)| > 0.3: %d / %d\n", n_plate_driven, nrow(pc_cov)))
if (n_plate_driven > 0)
  cat(sprintf(" Plate-driven PCs: %s\n",
              paste(pc_cov$PC[pc_cov$plate_driven], collapse = ", ")))

write_tab(pc_cov, "PC_Covariate_Correlations.csv")

# Heatmap of correlations
corr_long <- pc_cov %>%
  select(PC, var_pct, plate, nCount, log_viral, S_score) %>%
  pivot_longer(-c(PC, var_pct), names_to = "covariate", values_to = "rho")

p_pc_corr <- ggplot(corr_long,
                    aes(x = covariate, y = reorder(PC, as.numeric(str_remove(PC, "PC"))), fill = rho)) +
  geom_tile(colour = "grey80", linewidth = 0.15) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "Spearman rho") +
  scale_y_discrete(limits = rev) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 10)) +
  labs(title = "PC-covariate correlations (Spearman rho)",
       subtitle = sprintf("%d of %d PCs have |r(plate)| > 0.3. %s",
                          n_plate_driven, nrow(pc_cov),
                          if (n_plate_driven > 0) "Plate-driven PCs will trigger Harmony integration."
                          else "No plate-driven PCs; standard clustering is adequate."),
       x = "Covariate", y = "Principal component")
show_and_save(p_pc_corr, "PC_Correlation_Heatmap", width = 8.5, height = 11)

# Scree with plate overlay
p_scree <- ggplot(pc_cov, aes(x = as.numeric(str_remove(PC, "PC")), y = var_pct)) +
  geom_line(colour = "#808080") +
  geom_point(aes(size = plate_abs, colour = plate_driven), alpha = 0.85) +
  scale_colour_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"), name = "Plate-driven") +
  scale_size_continuous(range = c(1, 4), name = expression(paste("|", rho, "(plate)|"))) +
  theme_pub() +
  labs(title = "PCA scree with plate-correlation overlay",
       subtitle = "Firebrick = PCs where |r(plate)| > 0.3",
       x = "Principal component", y = "% variance explained")
show_and_save(p_scree, "PCA_Scree_Plate", width = 9, height = 5)

# PC1 vs PC2 coloured by plate
p_pc12 <- ggplot(data.frame(PC1 = pcs[, 1], PC2 = pcs[, 2], batch = qc$batch),
                 aes(PC1, PC2, colour = batch)) +
  geom_point(alpha = 0.72, size = 1.5) +
  scale_colour_viridis_d(name = "Plate") +
  theme_pub() +
  labs(title = "PC1 vs PC2 coloured by plate",
       subtitle = sprintf("r(PC1, plate) = %.3f | r(PC2, plate) = %.3f",
                          pc_cov$plate[1], pc_cov$plate[2]),
       x = sprintf("PC1 (%.1f%%)", pc_cov$var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", pc_cov$var_pct[2]))
show_and_save(p_pc12, "PC1_PC2_Plate", width = 8, height = 6)

## Harmony integration (if plate-driven)
USE_HARMONY <- n_plate_driven > 0 && requireNamespace("harmony", quietly = TRUE)

if (!USE_HARMONY && n_plate_driven > 0) {
  cat(" Harmony not installed; falling back to standard PCA.\n",
      " Plate-driven PCs: ", paste(pc_cov$PC[pc_cov$plate_driven], collapse = ", "),
      " will be retained. Install: BiocManager::install('Harmony')", sep = "")
}

if (USE_HARMONY) {
  cat(" -> Plate is driving PCs. Running Harmony integration (group.by = 'batch')...\n")
  
  # Check cache
  if (file.exists(hm_cache_path)) {
    cat(sprintf(" -> Loading cached Harmony results from %s\n", hm_cache_path))
    hm_pcs <- readRDS(hm_cache_path)
    if (nrow(hm_pcs) == nrow(pcs) && identical(rownames(hm_pcs), rownames(pcs))) {
      cat(" -> Cache valid.\n")
    } else {
      cat(" -> Cache stale; re-running Harmony.\n")
      hm_pcs <- NULL
    }
  } else {
    hm_pcs <- NULL
  }
  
  if (is.null(hm_pcs)) {
    harmony_obj <- harmony::RunHarmony(
      data_mat = pcs,
      meta_data = qc %>% select(batch),
      vars_use = "batch",
      theta = 2,
      lambda = 1,
      verbose = FALSE
    )
    colnames(harmony_obj) <- paste0("harmonyPC_", seq_len(ncol(harmony_obj)))
    hm_pcs <- as.data.frame(harmony_obj, row.names = rownames(pcs))
    rownames(hm_pcs) <- rownames(pcs)
    saveRDS(hm_pcs, hm_cache_path)
    cat(sprintf(" -> Harmony complete: %d PCs adjusted. Saved to cache.\n", ncol(hm_pcs)))
  }
  
  # Evaluate post-Harmony plate correlation
  plate_corr_post <- suppressWarnings(cor(hm_pcs[, 1], plate_num, method = "spearman"))
  cat(sprintf(" -> Post-Harmony: r(PC1, plate) = %.3f (was %.3f)\n",
              plate_corr_post, pc_cov$plate[1]))
  
  if (abs(plate_corr_post) > 0.15) {
    cat(" WARNING: Harmony did not fully remove the plate effect on PC1.\n")
  }
  
  # Post-Harmony scree
  ve_hm <- 100 * apply(as.matrix(hm_pcs), 2, var) / sum(apply(as.matrix(hm_pcs), 2, var))
  hm_pc_keep <- sum(ve_hm[seq_len(n_pc)] > 0.01)
  hm_pc_keep <- max(hm_pc_keep, min(10, n_pc))
  hm_pc_keep <- min(hm_pc_keep, ncol(hm_pcs))
  cat(sprintf(" -> Using first %d Harmony PCs (%.1f%% cumulative variance)\n",
              hm_pc_keep, sum(ve_hm[seq_len(hm_pc_keep)])))
  
  pcs_final <- as.matrix(hm_pcs[, seq_len(hm_pc_keep), drop = FALSE])
  pcs_source <- "Harmony-integrated"
} else {
  cat(" -> No plate correction applied. Using standard PCA.\n")
  hm_pc_keep <- n_pc
  pcs_final <- pcs
  pcs_source <- "Standard PCA"
}


cat(sprintf(" -> PCs for downstream: %d %s PCs (%.1f%% cumulative)\n",
            hm_pc_keep, pcs_source, sum(ve[seq_len(hm_pc_keep)])))
if (!USE_HARMONY && n_plate_driven > 0)
  cat(sprintf(" WARNING: %d plate-correlated PCs included. Consider installing Harmony.\n",
              n_plate_driven))

# ==============================================================================
# CLUSTERING
# ==============================================================================
cat(sprintf("\nSNN clustering (%s, k=20)\n", pcs_source))

k_nn <- 20
knn <- FNN::get.knn(pcs_final, k = k_nn)$nn.index
N <- nrow(pcs_final)
M <- sparseMatrix(i = rep(seq_len(N), each = k_nn),
                  j = as.vector(t(knn)), x = 1, dims = c(N, N))
diag(M) <- 1
inter <- as(M %*% t(M), "TsparseMatrix")
deg <- (k_nn + 1)
jac <- inter@x / (2 * deg - inter@x)
keep <- jac >= (1/15) & inter@i != inter@j
edges <- data.frame(from = inter@i[keep] + 1, to = inter@j[keep] + 1, weight = jac[keep])
edges <- edges[edges$from < edges$to, ]
g_snn <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                       vertices = data.frame(name = seq_len(N)))

# Resolution sweep + clustree-style decision
res_grid <- seq(0.1, 1.2, by = 0.1)
res_tab <- lapply(res_grid, function(r) {
  cl <- igraph::membership(igraph::cluster_louvain(g_snn,
                                                   weights = E(g_snn)$weight, resolution = r))
  data.frame(resolution = r, n_clusters = length(unique(cl)),
             modularity = igraph::modularity(g_snn, cl, weights = E(g_snn)$weight))
}) %>% bind_rows()
write_tab(res_tab, "Resolution_Sensitivity.csv")

# Prefer the resolution with highest modularity, but ensure minimum useful granularity
p_res <- ggplot(res_tab, aes(resolution, n_clusters)) +
  geom_line(colour = "#808080") + geom_point(aes(size = modularity), colour = "#0072B2") +
  geom_text(aes(label = n_clusters), vjust = -1.1, size = 3) +
  scale_x_continuous(breaks = res_grid) +
  theme_pub() +
  labs(title = "Louvain resolution sensitivity",
       subtitle = "Cluster count vs resolution (point size = graph modularity)",
       x = "Louvain resolution", y = "Number of clusters")
show_and_save(p_res, "Resolution_Sensitivity", width = 8, height = 5)

res_pick <- res_tab$resolution[which.max(res_tab$modularity)]
clusters <- factor(igraph::membership(
  igraph::cluster_louvain(g_snn, weights = E(g_snn)$weight, resolution = res_pick)))
qc$cluster <- clusters
cat(sprintf(" -> Resolution %.1f (modularity = %.3f) -> %d clusters\n",
            res_pick, max(res_tab$modularity), nlevels(clusters)))

# Cluster stability (50 subsamples, Jaccard > 0.6 threshold)
STAB_B <- 50; STAB_FRAC <- 0.80; STAB_TH <- 0.6
build_snn <- function(P, k = k_nn) {
  nn <- FNN::get.knn(P, k = k)$nn.index; nr <- nrow(P)
  Mb <- sparseMatrix(i = rep(seq_len(nr), each = k), j = as.vector(t(nn)), x = 1, dims = c(nr, nr))
  diag(Mb) <- 1
  it <- as(Mb %*% t(Mb), "TsparseMatrix")
  jj <- it@x / (2 * (k + 1) - it@x)
  kp <- jj >= (1/15) & it@i != it@j
  e <- data.frame(from = it@i[kp] + 1, to = it@j[kp] + 1, weight = jj[kp])
  e <- e[e$from < e$to, ]
  igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = seq_len(nr)))
}
calc_jaccard <- function(a, b) {
  u <- length(base::union(a, b))
  if (u == 0) NA_real_ else length(base::intersect(a, b)) / u
}

set.seed(42)
stab <- do.call(rbind, lapply(seq_len(STAB_B), function(b) {
  ix <- sort(sample(seq_len(N), floor(STAB_FRAC * N)))
  gb <- build_snn(pcs_final[ix, , drop = FALSE])
  nc <- igraph::membership(igraph::cluster_louvain(gb, weights = E(gb)$weight,
                                                   resolution = res_pick))
  vapply(levels(clusters), function(cl) {
    ref <- which(clusters[ix] == cl)
    if (!length(ref)) return(NA_real_)
    max(vapply(unique(nc), function(k) calc_jaccard(ref, which(nc == k)), numeric(1)))
  }, numeric(1))
}))
stab_mean <- colMeans(stab, na.rm = TRUE)
stab_df <- data.frame(
  cluster = levels(clusters),
  n_cells = as.integer(table(clusters)),
  mean_jaccard = round(as.numeric(stab_mean), 3),
  stable = ifelse(stab_mean >= STAB_TH, "stable", "unstable")
)
write_tab(stab_df, "Cluster_Stability.csv")
n_unstable <- sum(stab_df$stable == "unstable", na.rm = TRUE)
cat(sprintf(" -> Cluster stability (%d subsamples at %.0f%%): %d of %d clusters below Jaccard %.2f\n",
            STAB_B, 100 * STAB_FRAC, n_unstable, nlevels(clusters), STAB_TH))

# ==============================================================================
# UMAP / t-SNE EMBEDDING
# ==============================================================================
cat(sprintf("\nNon-linear embedding (%s)\n", pcs_source))

emb_name <- if (HAS_UWOT) "UMAP" else "t-SNE"
if (HAS_UWOT) {
  emb <- uwot::umap(pcs_final, n_neighbors = k_nn, min_dist = 0.3, verbose = FALSE)
} else {
  cat(sprintf(" uwot unavailable; using t-SNE (%s)\n", .uwot_err))
  emb <- Rtsne::Rtsne(pcs_final, dims = 2, pca = FALSE, check_duplicates = FALSE,
                      perplexity = 30, verbose = FALSE)$Y
}
qc$emb1 <- emb[, 1]; qc$emb2 <- emb[, 2]

emb_panel <- function(colour_expr, title, legend_name) {
  ggplot(qc, aes(emb1, emb2, colour = {{ colour_expr }})) +
    geom_point(size = 1.0, alpha = 0.82) + theme_pub() +
    labs(title = title, x = paste0(emb_name, " 1"), y = paste0(emb_name, " 2"),
         colour = legend_name) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
}

p_umap4 <- (emb_panel(cluster, "Cluster", "cluster") |
              emb_panel(infection_status, "Infection status", "status")) /
  (emb_panel(factor(hpi), "Hours post-infection", "hpi") |
     emb_panel(factor(moi), "MOI", "MOI")) +
  patchwork::plot_annotation(
    title = sprintf("%s embedding (n = %d cells, %d PCs, resolution %.1f, %s)",
                    emb_name, N, hm_pc_keep, res_pick, pcs_source),
    subtitle = sprintf("%d Louvain clusters; %d fall below Jaccard %.2f over %d subsamples.",
                       nlevels(clusters), n_unstable, STAB_TH, STAB_B))
show_and_save(p_umap4, "Embedding_4panel", width = 12, height = 10)

# Embedding agreement: t-SNE vs primary embedding
if (HAS_UWOT) {
  ts <- Rtsne::Rtsne(pcs_final, dims = 2, pca = FALSE, check_duplicates = FALSE,
                     perplexity = 30, verbose = FALSE)$Y
  knn_overlap <- function(A, B, k = 30) {
    ia <- FNN::get.knn(A, k = k)$nn.index
    ib <- FNN::get.knn(B, k = k)$nn.index
    mean(vapply(seq_len(nrow(A)), function(i) length(intersect(ia[i, ], ib[i, ])) / k, numeric(1)))
  }
  kov <- knn_overlap(pcs_final, emb, 30)
  p_agree <- (ggplot(qc, aes(ts[, 1], ts[, 2], colour = cluster)) +
                geom_point(size = .9, alpha = .82) + theme_pub() + guides(colour = "none") +
                labs(title = "t-SNE", x = "t-SNE 1", y = "t-SNE 2")) |
    (ggplot(data.frame(PC1 = pcs_final[, 1], PC2 = pcs_final[, 2], cluster = clusters),
            aes(PC1, PC2, colour = cluster)) + geom_point(size = .9, alpha = .82) +
       theme_pub() + guides(colour = "none") + labs(title = "PC1/PC2"))
  p_agree <- p_agree + patchwork::plot_annotation(
    title = "Embedding agreement",
    subtitle = sprintf("Mean 30-NN preservation, %s vs PCA = %.2f", emb_name, kov))
  show_and_save(p_agree, "Embedding_Agreement", width = 10, height = 5)
}

# ==============================================================================
# WITHIN-PLATE WILCOXON (Infected vs Bystander)
# ==============================================================================
cat("\nWithin-plate Wilcoxon rank-sum (Infected vs Bystander)\n")

MIN_PER_GROUP <- 20

# Map symbols to active matrix features
feat_map <- setNames(rownames(mat_log), sym_of(rownames(mat_log)))

within_plate_wilcox <- function(plate_cells, mat_sub) {
  sc1 <- plate_cells$Sample[plate_cells$infection_status == "Infected"]
  sc2 <- plate_cells$Sample[plate_cells$infection_status == "Bystander"]
  if (length(sc1) < MIN_PER_GROUP || length(sc2) < MIN_PER_GROUP) return(NULL)
  
  idx1 <- match(sc1, colnames(mat_sub))
  idx2 <- match(sc2, colnames(mat_sub))
  if (anyNA(idx1) || anyNA(idx2)) return(NULL)
  
  R <- matrixStats::rowRanks(mat_sub[, c(idx1, idx2), drop = FALSE],
                             ties.method = "average")
  n1 <- length(sc1); n2 <- length(sc2); n <- n1 + n2
  grp1 <- c(rep(TRUE, n1), rep(FALSE, n2))
  U1 <- rowSums(R[, grp1, drop = FALSE]) - n1 * (n1 + 1) / 2
  mu <- n1 * n2 / 2
  varU <- (n1 * n2 / (n * (n - 1))) *
    (rowSums(R^2) - n * ((n + 1) / 2)^2)
  z <- (U1 - mu - sign(U1 - mu) * 0.5) / sqrt(varU)
  p <- 2 * pnorm(-abs(z))
  fc <- rowMeans(mat_sub[, grp1, drop = FALSE]) - rowMeans(mat_sub[, !grp1, drop = FALSE])
  pct1 <- rowMeans(mat_sub[, grp1, drop = FALSE] > 0)
  pct2 <- rowMeans(mat_sub[, !grp1, drop = FALSE] > 0)
  
  data.frame(
    ensembl = rownames(mat_sub),
    symbol = sym_of(rownames(mat_sub)),
    avg_log2FC = fc,
    pct_in_infected = pct1,
    pct_in_bystander = pct2,
    p_value = p,
    p_adj = p.adjust(p),
    n_cells = n1 + n2,
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(p_value), pct_in_infected > 0 | pct_in_bystander > 0) %>%
    mutate(stratum = paste0(plate_cells$virus[1], "_MOI", plate_cells$moi[1],
                            "_", plate_cells$time_h[1], "h_batch", plate_cells$batch[1],
                            "_", plate_cells$infection_status[1], "vs", plate_cells$infection_status[2])) %>%
    arrange(p_adj)
}

# Run within-plate DE for each plate that has both groups >= 20 cells
contrast_plates <- qc %>%
  filter(infection_status %in% c("Infected", "Bystander")) %>%
  count(batch, virus, moi, time_h, infection_status) %>%
  pivot_wider(names_from = infection_status, values_from = n, values_fill = 0) %>%
  filter(Infected >= MIN_PER_GROUP, Bystander >= MIN_PER_GROUP) %>%
  mutate(stratum = sprintf("%s_MOI%g_%dh", virus, moi, time_h))

cat(sprintf(" %d plates support the Infected vs Bystander contrast (>= %d per group).\n",
            nrow(contrast_plates), MIN_PER_GROUP))
if (nrow(contrast_plates) == 0) {
  cat(" No plates meet the feasibility threshold.  results will be empty.\n")
}

de_list <- list()
for (i in seq_len(nrow(contrast_plates))) {
  bt <- contrast_plates$batch[i]
  v <- contrast_plates$virus[i]
  cells <- qc$Sample[qc$batch == bt & qc$infection_status %in% c("Infected", "Bystander")]
  plate_qc <- qc %>% filter(Sample %in% cells)
  
  res <- within_plate_wilcox(plate_qc, mat_log)
  if (!is.null(res) && nrow(res) > 0) {
    de_list[[i]] <- res
    cat(sprintf(" Plate %s (%s, %s): %d significant genes (FDR < 0.05)\n",
                bt, v, contrast_plates$stratum[i], sum(res$p_adj < 0.05)))
  } else {
    cat(sprintf(" Plate %s: no significant genes or insufficient cells.\n", bt))
  }
}
de_all <- bind_rows(de_list)
write_tab(de_all, "DE_WithinPlate_Wilcoxon.csv")

# Summary: how many significant genes per stratum
if (nrow(de_all) > 0) {
  sig_sum <- de_all %>%
    group_by(stratum) %>%
    summarise(n_sig = sum(p_adj < 0.05 & abs(avg_log2FC) >= 0.25),
              n_tested = n(),
              pct_sig = 100 * n_sig / n_tested, .groups = "drop")
  write_tab(sig_sum, "DE_Summary.csv")
  cat(" DE summary per plate:\n")
  print(as.data.frame(sig_sum), row.names = FALSE)
}

# ==============================================================================
# ARTIFACT FILTER on within-plate DE
# ==============================================================================
cat("\nArtifact filtering on within-plate DE results\n")

artifact_ensg <- character(0)

## MT genes
if (!is.null(gene_sets) && !is.null(gene_sets$mt)) {
  artifact_ensg <- union(artifact_ensg, gene_sets$mt)
} else {
  mt_chrM <- c("ENSG00000210049","ENSG00000211459","ENSG00000210077","ENSG00000210082",
               "ENSG00000209082","ENSG00000198888","ENSG00000210100","ENSG00000210107",
               "ENSG00000210112","ENSG00000198763","ENSG00000210117","ENSG00000210127",
               "ENSG00000210135","ENSG00000210140","ENSG00000210144","ENSG00000198804",
               "ENSG00000210151","ENSG00000210154","ENSG00000198712","ENSG00000210156",
               "ENSG00000228253","ENSG00000198899","ENSG00000198938","ENSG00000210164",
               "ENSG00000198840","ENSG00000210174","ENSG00000212907","ENSG00000198886",
               "ENSG00000210176","ENSG00000210184","ENSG00000210191","ENSG00000198786",
               "ENSG00000198695","ENSG00000210194","ENSG00000198727","ENSG00000210195",
               "ENSG00000210196")
  artifact_ensg <- union(artifact_ensg, intersect(mt_chrM, rownames(mat_counts)))
}

## Ribosomal genes
if (!is.null(gene_sets) && !is.null(gene_sets$ribosomal)) {
  artifact_ensg <- union(artifact_ensg, gene_sets$ribosomal)
} else {
  ensg_col_local <- intersect(c("ENSG", "ENSG_Original", "ENSG_Clean"), names(gene_map))[1]
  if (!is.na(ensg_col_local) && !is.null(gene_map$Symbol)) {
    rp_ensg <- gene_map[[ensg_col_local]][!is.na(gene_map$Symbol) &
                                            str_detect(gene_map$Symbol, "^RP[LS]\\d")]
    artifact_ensg <- union(artifact_ensg, intersect(rp_ensg, rownames(mat_counts)))
  }
}

## ERCC features
ercc_ensg <- rownames(mat_counts)[str_detect(rownames(mat_counts), "^ERCC-")]
artifact_ensg <- union(artifact_ensg, ercc_ensg)

## ISG (interferon-stimulated gene) antiviral response genes
## These are expected to be upregulated in any infected cell â€” they are
## canonical antiviral responders, not condition-specific DEGs.
## Removed from all DE tables to focus on non-ISG biology.
isg_symbols_raw <- c(
  ## RNA sensors / PRRs
  "DDX58", "IFIH1", "OAS1", "OAS2", "OAS3", "OASL", "ZCCHC3", "ZNFX1",
  ## RNA degradation / translation inhibition
  "EIF2AK2", "RNASEL", "ISG20", "ZC3HAV1", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
  "ADAR1", "PARP12", "PARP10", "PARP11", "SHFL", "MOV10",
  ## GTPases / helicases
  "MX1", "MX2", "GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7", "IFI16",
  ## Membrane / entry / egress blockade
  "BST2", "IFITM1", "IFITM2", "IFITM3", "PLSCR1", "C6orf150", "HPSE",
  ## Protein modification / ubiquitination
  "ISG15", "TRIM25", "TRIM22", "TRIM5A", "TRIM32", "TRIM38", "TRIM62",
  "FAM110A", "FCAR", "FGL2",
  ## Metabolic / other effectors
  "RSAD2", "IFI44L", "IFI44", "CH25H", "IDO1", "LY6E", "CD74", "NCOA7", "SAMHD1",
  "APOBEC3F", "APOBEC3G", "APOBEC3B", "SLFN12", "C11orf83", "CMPK2",
  "DAXX", "MUC1", "MUC4", "MUC5AC", "MUC5B", "MUC6", "MUC12",
  ## Signaling / transcription factors
  "IRF1", "IRF7", "IRF9", "STAT1", "STAT2", "IFNB1",
  "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8",
  "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21",
  "IFNL1", "IFNL2", "IFNL3", "IFNL4",
  ## Additional ISG Atlas confirmed genes
  "CFB", "ADAMDEC1", "RTCB", "GPR141", "POLI", "LGALS9B", "ATXN1",
  "EHD4", "NOS2A", "TLK2", "TNFRSF10A", "PSMB8", "PSMB9", "IFI6",
  "RTP4", "CCDC75", "IFI27", "NAMPT", "CASP2", "PML", "SP110",
  "SERPINA5", "SERPINA10", "GPR146", "MITD1", "TAP1", "RHOH", "PVRL4", "MGAT1",
  "HERC5", "USP18", "UBE2L6", "PARP9", "MB21D1", "CGAS", "DTX3L", "LGP2",
  "CCL2", "CCL5", "CXCL10", "CXCL11")
isg_ensg <- gene_map[[ensg_col]][match(toupper(isg_symbols_raw), toupper(gene_map$Symbol))]
isg_ensg <- isg_ensg[!is.na(isg_ensg)]
isg_ensg <- intersect(isg_ensg, rownames(mat_counts))
artifact_ensg <- union(artifact_ensg, isg_ensg)
cat(sprintf(" -> ISG genes added to artifact filter: %d of %d found in matrix\n",
            length(isg_ensg), length(isg_symbols_raw)))

# Artifact filter deferred to after mock_sig is computed in -B


# ==============================================================================
# BYSTANDER VS MOCK (artifact check, )
# ==============================================================================
cat("
Bystander vs Mock (artifact check)\n")

mock_plates <- qc %>%
  filter(infection_status %in% c("Bystander", "Unexposed")) %>%
  count(batch, virus, moi, time_h, infection_status) %>%
  pivot_wider(names_from = infection_status, values_from = n, values_fill = 0) %>%
  filter(Bystander >= MIN_PER_GROUP, Unexposed >= MIN_PER_GROUP)

cat(sprintf(" %d plates support Bystander vs Unexposed (>= %d per group).\n",
            nrow(mock_plates), MIN_PER_GROUP))

mock_de_list <- list()
for (i in seq_len(nrow(mock_plates))) {
  bt <- mock_plates$batch[i]
  cells <- qc$Sample[qc$batch == bt & qc$infection_status %in% c("Bystander", "Unexposed")]
  plate_qc <- qc %>% filter(Sample %in% cells)
  res <- within_plate_wilcox(plate_qc, mat_log)
  if (!is.null(res) && nrow(res) > 0) {
    res$stratum <- paste0(plate_qc$virus[1], "_MOI", plate_qc$moi[1],
                          "_", plate_qc$time_h[1], "h_batch", bt, "_BystandervsUnexposed")
    mock_de_list[[i]] <- res
    cat(sprintf(" Plate %s: %d sig genes (FDR<0.05, |log2FC|>=0.25)\n",
                bt, sum(res$p_adj < 0.05 & abs(res$avg_log2FC) >= 0.25)))
  }
}
mock_de_all <- bind_rows(mock_de_list)
if (nrow(mock_de_all) > 0) {
  write_tab(mock_de_all, "Bystander_vs_Unexposed.csv")
  mock_sig <- mock_de_all %>%
    filter(p_adj < 0.05 & abs(avg_log2FC) >= 0.25) %>%
    pull(ensembl) %>% unique()
  cat(sprintf(" -> Bystander vs Unexposed: %d significant genes (artifact candidates).\n",
              length(mock_sig)))
} else {
  cat(" -> No plates meet Bystander vs Unexposed feasibility. Skipping.\n")
  mock_sig <- character(0)
}

## Augment with mock-responders and apply full filter
artifact_ensg <- union(artifact_ensg, mock_sig)
if (length(mock_sig) > 0)
  cat(sprintf(" Artifact list now includes %d mock-responder genes.\n",
              length(mock_sig)))

if (nrow(de_all) > 0 && length(artifact_ensg) > 0) {
  before <- nrow(de_all)
  de_all <- de_all %>% filter(!(ensembl %in% artifact_ensg))
  removed <- before - nrow(de_all)
  cat(sprintf("  filter: removed %d rows (%d -> %d)\n",
              removed, before, nrow(de_all)))
  
  ## Re-write cleaned DE table and summary
  write_tab(de_all, "DE_WithinPlate_Wilcoxon.csv")
  
  if (nrow(de_all) > 0) {
    sig_sum <- de_all %>%
      group_by(stratum) %>%
      summarise(n_sig = sum(p_adj < 0.05 & abs(avg_log2FC) >= 0.25),
                n_tested = n(),
                pct_sig = 100 * n_sig / n_tested, .groups = "drop")
    write_tab(sig_sum, "DE_Summary.csv")
  }
}

# ==============================================================================
# PSEUDOBULK (only if >= 3 plates per condition)
# ==============================================================================
cat("\nPseudobulk differential expression\n")
# Pseudobulk requires >=3 plates for adequate residual degrees of freedom.
# Only run for conditions with >= 3 plates carrying both Infected and Bystander.
# Viral rows are EXCLUDED from the pseudobulk matrix.

plate_feas <- qc %>%
  filter(infection_status %in% c("Infected", "Bystander")) %>%
  count(batch, virus, moi, time_h, infection_status) %>%
  pivot_wider(names_from = infection_status, values_from = n, values_fill = 0)

pseudobulk_eligible <- plate_feas %>%
  filter(Infected >= MIN_PER_GROUP, Bystander >= MIN_PER_GROUP) %>%
  count(virus, moi, time_h, name = "n_plates") %>%
  filter(n_plates >= 3)

cat(sprintf(" Conditions with >= 3 plates and >= %d per group: %d\n",
            MIN_PER_GROUP, nrow(pseudobulk_eligible)))
if (nrow(pseudobulk_eligible) > 0)
  print(as.data.frame(pseudobulk_eligible), row.names = FALSE)

if (nrow(pseudobulk_eligible) > 0) {
  # Load edgeR and limma if available
  if (requireNamespace("edgeR", quietly = TRUE) && requireNamespace("limma", quietly = TRUE)) {
    for (i in seq_len(nrow(pseudobulk_eligible))) {
      v <- pseudobulk_eligible$virus[i]
      mo <- pseudobulk_eligible$moi[i]
      th <- pseudobulk_eligible$time_h[i]
      cat(sprintf(" -> Pseudobulk: %s MOI%g %dh\n", v, mo, th))
      
      eligible_plates <- plate_feas %>%
        filter(virus == v, moi == mo, time_h == th,
               Infected >= MIN_PER_GROUP, Bystander >= MIN_PER_GROUP)
      
      # Get sample indices for each plate x status
      inf_samples <- unlist(lapply(eligible_plates$batch[eligible_plates$Infected >= MIN_PER_GROUP],
                                   function(bt) qc$Sample[qc$batch == bt & qc$infection_status == "Infected"]))
      bys_samples <- unlist(lapply(eligible_plates$batch[eligible_plates$Bystander >= MIN_PER_GROUP],
                                   function(bt) qc$Sample[qc$batch == bt & qc$infection_status == "Bystander"]))
      
      if (length(inf_samples) < MIN_PER_GROUP || length(bys_samples) < MIN_PER_GROUP) {
        cat(" Insufficient cells after filtering; skipping.\n")
        next
      }
      
      # Aggregate by plate (do NOT include viral rows)
      inf_batches <- eligible_plates$batch[eligible_plates$Infected >= MIN_PER_GROUP]
      bys_batches <- eligible_plates$batch[eligible_plates$Bystander >= MIN_PER_GROUP]
      
      # Remove viral genes
      d_viral <- str_detect(rownames(mat_counts), "^(DENV|ZIKV)")
      pb_mat <- mat_counts[!d_viral, , drop = FALSE]
      
      agg_fun <- function(samples) {
        if (length(samples) == 0) return(NULL)
        rowMeans(pb_mat[, samples, drop = FALSE])
      }
      
      inf_agg <- sapply(inf_batches, function(bt)
        agg_fun(qc$Sample[qc$batch == bt & qc$infection_status == "Infected"]))
      bys_agg <- sapply(bys_batches, function(bt)
        agg_fun(qc$Sample[qc$batch == bt & qc$infection_status == "Bystander"]))
      
      if (is.null(dim(inf_agg))) inf_agg <- matrix(inf_agg, ncol = 1)
      if (is.null(dim(bys_agg))) bys_agg <- matrix(bys_agg, ncol = 1)
      
      pb_counts <- cbind(inf_agg, bys_agg)
      colnames(pb_counts) <- c(
        paste0(inf_batches, "_Infected"),
        paste0(bys_batches, "_Bystander")
      )
      pb_group <- rep(c("Infected", "Bystander"), c(ncol(inf_agg), ncol(bys_agg)))
      pb_plate <- rep(c(inf_batches, bys_batches), c(ncol(inf_agg), ncol(bys_agg)))
      
      dge <- edgeR::DGEList(counts = pb_counts, group = pb_group)
      dge <- edgeR::filterByExpr(dge)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      
      design <- model.matrix(~ pb_group + pb_plate)
      vfit <- edgeR::voom(dge, design, plot = FALSE)
      fit <- limma::lmFit(vfit, design)
      fit <- limma::eBayes(fit)
      
      ct <- limma::makeContrasts(Infected_vs_Bystander = pb_groupInfected, levels = design)
      cfit <- limma::contrasts.fit(fit, ct)
      cfit <- limma::eBayes(cfit)
      
      res <- limma::topTable(cfit, n = Inf, adjust.method = "BH") %>%
        as.data.frame() %>%
        mutate(ensembl = rownames(.),
               symbol = sym_of(rownames(.))) %>%
        select(ensembl, symbol, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
        rename(avg_log2FC = logFC, p_adj = adj.P.Val, p_value = P.Value) %>%
        filter(is.finite(p_value), !is.na(symbol))
      
      write_tab(res, sprintf("Pseudobulk_%s_MOI%g_%dh.csv", v, mo, th))
      cat(sprintf(" Pseudobulk: %d significant genes (FDR < 0.05, |log2FC| >= %.2f)\n",
                  sum(res$p_adj < 0.05 & abs(res$avg_log2FC) >= 0.25), 0.25))
    }
  } else {
    cat(" edgeR/limma not available; pseudobulk skipped.\n")
    cat(" Install with: install.packages(c('edgeR', 'limma', 'statmod'))\n")
    write_tab(data.frame(
      message = "edgeR/limma not installed. Pseudobulk requires these packages.",
      stringsAsFactors = FALSE
    ), "Pseudobulk_Skipped.csv")
  }
} else {
  cat(" No conditions have >= 3 plates with both groups >= 20 cells.\n")
  cat(" Pseudobulk skipped (feasibility gate).\n")
  write_tab(data.frame(
    condition = "none",
    message = "No condition has >= 3 plates with >= 20 cells per group.",
    stringsAsFactors = FALSE
  ), "Pseudobulk_Skipped.csv")
}

# ==============================================================================
# CONTINUOUS vRNA ASSOCIATION (PRIMARY RESULT)
# ==============================================================================
cat("\nContinuous vRNA association (primary result)\n")

# Exposed cells only (Infected + Bystander, excluding MOI 0 controls)
exposed <- qc %>% filter(moi > 0, infection_status %in% c("Infected", "Bystander"))
cat(sprintf(" Exposed cells: %d (Infected + Bystander)\n", nrow(exposed)))

# Spearman rho: gene expression vs log_viral, separately per virus
vrna_list <- list()
for (v in unique(exposed$virus)) {
  v_cells <- exposed %>% filter(virus == v)
  if (nrow(v_cells) < 10) {
    cat(sprintf(" %s: only %d exposed cells, skipping.\n", v, nrow(v_cells)))
    next
  }
  
  v_ix <- match(v_cells$Sample, colnames(mat_log))
  lv <- v_cells$log_viral
  rho_v <- vapply(seq_len(nrow(mat_log)), function(g) {
    x <- mat_log[g, v_ix]
    ok <- is.finite(x) & is.finite(lv)
    if (sum(ok) < 10) return(NA_real_)
    if (sd(x[ok]) == 0 || sd(lv[ok]) == 0) return(NA_real_)
    cor(x[ok], lv[ok], method = "spearman")
  }, numeric(1))
  
  p_v <- vapply(seq_len(nrow(mat_log)), function(g) {
    x <- mat_log[g, v_ix]
    ok <- is.finite(x) & is.finite(lv)
    if (sum(ok) < 10) return(NA_real_)
    cor.test(x[ok], lv[ok], method = "spearman")$p.value
  }, numeric(1))
  
  res <- data.frame(
    ensembl = rownames(mat_log),
    symbol = sym_of(rownames(mat_log)),
    virus = v,
    rho = rho_v,
    p_value = p_v,
    p_adj = p.adjust(p_v),
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(rho), !is.na(p_adj)) %>%
    arrange(p_adj)
  
  vrna_list[[v]] <- res
  cat(sprintf(" %s: %d genes tested, %d significant (BH FDR < 0.05, |rho| > 0.15)\n",
              v, nrow(res), sum(res$p_adj < 0.05 & abs(res$rho) > 0.15)))
}

vrna <- bind_rows(vrna_list)
write_tab(vrna, "vRNA_Continuous_Association.csv")

# Volcano plot
if (nrow(vrna) > 0) {
  for (v in unique(vrna$virus)) {
    sub <- vrna %>% filter(virus == v)
    sub$neg_log10_p <- -log10(sub$p_value + 1e-300)
    top_pos <- sub %>% filter(rho > 0) %>% arrange(desc(rho)) %>% slice_head(n = 10)
    top_neg <- sub %>% filter(rho < 0) %>% arrange(rho) %>% slice_head(n = 10)
    
    p_vl <- ggplot(sub, aes(x = rho, y = neg_log10_p)) +
      geom_point(aes(colour = rho > 0), alpha = 0.4, size = 0.8) +
      geom_text_repel(data = bind_rows(top_pos, top_neg), aes(label = symbol),
                      size = 2.5, max.overlaps = 15, min.segment.length = 0, colour = "black") +
      scale_colour_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2")) +
      theme_pub() +
      labs(title = sprintf("Continuous vRNA association (%s)", v),
           subtitle = "Exposed cells only. Spearman rho vs log10(viral molecules + 1). BH-corrected.",
           x = "Spearman rho", y = "-log10(p-value)")
    show_and_save(p_vl, sprintf("vRNA_Volcano_%s", v), width = 9, height = 7)
  }
}

# ==============================================================================
# TIMEPOINT-SPECIFIC DE, VOLCANO & MA PLOTS, & DEG MATRICES
# ==============================================================================
cat("\nTimepoint-specific DE (Infected vs Bystander), Plots, and DEG Matrices\n")

# 1. Fast Wilcoxon for pooled cells across all plates of a given virus+timepoint
run_pooled_wilcox <- function(meta_sub, mat_sub) {
  inf_cells <- meta_sub$Sample[meta_sub$infection_status == "Infected"]
  bys_cells <- meta_sub$Sample[meta_sub$infection_status == "Bystander"]
  if (length(inf_cells) < 10 || length(bys_cells) < 10) return(NULL)
  
  idx_inf <- match(inf_cells, colnames(mat_sub))
  idx_bys <- match(bys_cells, colnames(mat_sub))
  all_idx <- c(idx_inf, idx_bys)
  
  ave_expr <- rowMeans(mat_sub[, all_idx, drop = FALSE])
  fc <- rowMeans(mat_sub[, idx_inf, drop = FALSE]) - rowMeans(mat_sub[, idx_bys, drop = FALSE])
  
  R <- matrixStats::rowRanks(mat_sub[, all_idx, drop = FALSE], ties.method = "average")
  n1 <- length(idx_inf); n2 <- length(idx_bys); n <- n1 + n2
  grp1 <- c(rep(TRUE, n1), rep(FALSE, n2))
  U1 <- rowSums(R[, grp1, drop = FALSE]) - n1 * (n1 + 1) / 2
  mu <- n1 * n2 / 2
  varU <- (n1 * n2 / (n * (n - 1))) * (rowSums(R^2) - n * ((n + 1) / 2)^2)
  z <- ifelse(varU > 0, (U1 - mu - sign(U1 - mu) * 0.5) / sqrt(varU), 0)
  p <- 2 * pnorm(-abs(z))
  
  data.frame(
    ensembl = rownames(mat_sub),
    symbol = sym_of(rownames(mat_sub)),
    AveExpr = ave_expr,
    avg_log2FC = fc,
    p_value = p,
    p_adj = p.adjust(p, method = "BH"),
    stringsAsFactors = FALSE
  )
}

# 2. Iterate through each virus and timepoint
viruses <- sort(unique(qc$virus))
timepoints <- sort(unique(qc$time_h[qc$time_h > 0]))
all_timepoint_de <- list()

for (v in viruses) {
  for (t in timepoints) {
    meta_sub <- qc %>% filter(virus == v, time_h == t,
                              infection_status %in% c("Infected", "Bystander"))
    
    res <- run_pooled_wilcox(meta_sub, mat_log)
    if (is.null(res)) {
      cat(sprintf(" Skipping %s %dh: insufficient cells.\n", v, t))
      next
    }
    
    res$virus <- v
    res$time_h <- t
    res$condition <- paste0(v, "_", t, "h")
    res$is_sig <- !is.na(res$p_adj) & res$p_adj < 0.05 & abs(res$avg_log2FC) >= 0.25
    all_timepoint_de[[res$condition[1]]] <- res
    
    top_genes <- res %>% filter(is_sig) %>% arrange(p_adj) %>% slice_head(n = 10)
    
    # Volcano plot
    p_volc <- ggplot(res, aes(x = avg_log2FC, y = -log10(p_adj + 1e-300), colour = is_sig)) +
      geom_point(alpha = 0.62, size = 1.15) +
      scale_colour_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#CCCCCC"), name = NULL) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "black", linewidth = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black", linewidth = 0.5) +
      theme_pub() +
      labs(title = sprintf("Volcano Plot (%s %dh)", v, t),
           subtitle = sprintf("%d Significant DEGs (Infected vs Bystander)", sum(res$is_sig, na.rm = TRUE)),
           x = "Log2 Fold Change", y = "-log10(FDR)") +
      theme(legend.position = "none")
    
    if (nrow(top_genes) > 0)
      p_volc <- p_volc + geom_text_repel(data = top_genes, aes(label = symbol),
                                         colour = "black", size = 3, max.overlaps = 15)
    show_and_save(p_volc, sprintf("Volcano_%s_%dh", v, t), width = 8.5, height = 7)
    
    # MA plot
    p_ma <- ggplot(res, aes(x = AveExpr, y = avg_log2FC, colour = is_sig)) +
      geom_point(alpha = 0.62, size = 1.15) +
      scale_colour_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#CCCCCC"), name = NULL) +
      geom_hline(yintercept = 0, colour = "black", linewidth = 0.8) +
      theme_pub() +
      labs(title = sprintf("MA Plot (%s %dh)", v, t),
           subtitle = "Average Expression vs Log2 Fold Change",
           x = "Average Log-Expression", y = "Log2 Fold Change") +
      theme(legend.position = "none")
    
    if (nrow(top_genes) > 0)
      p_ma <- p_ma + geom_text_repel(data = top_genes, aes(label = symbol),
                                     colour = "black", size = 3, max.overlaps = 15)
    show_and_save(p_ma, sprintf("MA_Plot_%s_%dh", v, t), width = 8.5, height = 7)
    
    cat(sprintf(" %s %dh: %d sig DEGs\n", v, t, sum(res$is_sig)))
  }
}

# 3. Save aggregated DEG lists and expression matrices
if (length(all_timepoint_de) > 0) {
  full_de_df <- bind_rows(all_timepoint_de)
  write_tab(full_de_df, "All_Timepoint_DEGs.csv")
  
  zika_sig <- full_de_df %>% filter(virus == "zika", is_sig) %>% pull(ensembl) %>% unique()
  dengue_sig <- full_de_df %>% filter(virus == "dengue", is_sig) %>% pull(ensembl) %>% unique()
  common_sig <- intersect(zika_sig, dengue_sig)
  
  # --- Artifact filter: remove MT-, RP-, ERCC, ISG, and Bystander-vs-Unexposed artifact genes ---
  artifact_genes <- rownames(mat_log)[
    str_detect(rownames(mat_log), "^MT-|^MT\\.|^RP[LS]\\d|^RP[LS]\\d+[A-Z]|ERCC|^\\.[0-9]+$")]
  artifact_genes <- union(artifact_genes, isg_ensg) # ISGs added in 
  n_isg_hit <- length(intersect(zika_sig, isg_ensg))
  clean_sig <- function(genes) setdiff(setdiff(genes, artifact_genes), mock_sig)
  zika_sig_clean <- clean_sig(zika_sig)
  dengue_sig_clean <- clean_sig(dengue_sig)
  common_sig_clean <- clean_sig(common_sig)
  n_rm <- length(zika_sig) - length(zika_sig_clean)
  n_rm_mock <- length(intersect(zika_sig, mock_sig))
  cat(sprintf(" -> Artifact filter removed %d ZIKV DEGs (%d MT-/RP-/ERCC, %d ISG, %d Mock-responding)\n",
              n_rm, length(intersect(zika_sig, artifact_genes)) - n_isg_hit,
              n_isg_hit, n_rm_mock))
  
  cat(sprintf("\n -> ZIKV: %d | DENV: %d | Common: %d DEGs across all timepoints.\n",
              length(zika_sig), length(dengue_sig), length(common_sig)))
  
  
  # --- ZIKV DEG table with stats ---
  zika_deg <- full_de_df %>%
    filter(virus == "zika", is_sig) %>%
    group_by(ensembl) %>%
    summarise(symbol = first(symbol), avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
              p_adj = min(p_adj, na.rm = TRUE), .groups = "drop") %>%
    mutate(Regulation = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DOWN", "NS"))) %>%
    select(ensembl, symbol, avg_log2FC, p_adj, Regulation) %>%
    arrange(desc(abs(avg_log2FC)))
  write_tab(zika_deg, "ZIKV_Sig_Genes.csv")
  cat(sprintf(" -> ZIKV DEGs: %d genes\n", nrow(zika_deg)))
  
  # --- DENV DEG table with stats ---
  denv_deg <- full_de_df %>%
    filter(virus == "dengue", is_sig) %>%
    group_by(ensembl) %>%
    summarise(symbol = first(symbol), avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
              p_adj = min(p_adj, na.rm = TRUE), .groups = "drop") %>%
    mutate(Regulation = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DOWN", "NS"))) %>%
    select(ensembl, symbol, avg_log2FC, p_adj, Regulation) %>%
    arrange(desc(abs(avg_log2FC)))
  write_tab(denv_deg, "DENV_Sig_Genes.csv")
  cat(sprintf(" -> DENV DEGs: %d genes\n", nrow(denv_deg)))
  
  # --- Common DEG table with stats ---
  common_deg <- full_de_df %>%
    filter(ensembl %in% common_sig, is_sig) %>%
    group_by(ensembl, virus) %>%
    summarise(symbol = first(symbol), avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
              p_adj = min(p_adj, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = virus, values_from = c(avg_log2FC, p_adj),
                names_glue = "{.value}_{virus}") %>%
    mutate(Regulation_Dengue = ifelse(avg_log2FC_dengue > 0, "UP", ifelse(avg_log2FC_dengue < 0, "DOWN", "NS")),
           Regulation_Zika = ifelse(avg_log2FC_zika > 0, "UP", ifelse(avg_log2FC_zika < 0, "DOWN", "NS"))) %>%
    rename(avg_log2FC_Dengue = avg_log2FC_dengue, p_adj_Dengue = p_adj_dengue,
           avg_log2FC_Zika = avg_log2FC_zika, p_adj_Zika = p_adj_zika) %>%
    select(ensembl, symbol, avg_log2FC_Dengue, p_adj_Dengue, Regulation_Dengue,
           avg_log2FC_Zika, p_adj_Zika, Regulation_Zika) %>%
    arrange(ensembl)
  write_tab(common_deg, "Common_Sig_Genes.csv")
  cat(sprintf(" -> Common DEGs: %d genes\n", nrow(common_deg)))
  
  # --- Merged wide-format DEG comparison table ---
  virus_summary <- full_de_df %>%
    group_by(ensembl, virus) %>%
    summarise(
      avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
      p_adj = min(p_adj, na.rm = TRUE),
      .groups = "drop") %>%
    group_by(ensembl) %>%
    mutate(Regulation = ifelse(avg_log2FC > 0, "UP", ifelse(avg_log2FC < 0, "DOWN", "NS"))) %>%
    ungroup()
  
  virus_wide <- virus_summary %>%
    select(ensembl, virus, avg_log2FC, p_adj, Regulation) %>%
    pivot_wider(names_from = virus,
                values_from = c(avg_log2FC, p_adj, Regulation),
                names_glue = "{.value}_{virus}") %>%
    rename(avg_log2FC_Dengue = avg_log2FC_dengue,
           p_adj_Dengue = p_adj_dengue,
           Regulation_Dengue = Regulation_dengue,
           avg_log2FC_Zika = avg_log2FC_zika,
           p_adj_Zika = p_adj_zika,
           Regulation_Zika = Regulation_zika) %>%
    mutate(symbol = sym_of(ensembl)) %>%
    select(ensembl, symbol, avg_log2FC_Dengue, p_adj_Dengue, Regulation_Dengue,
           avg_log2FC_Zika, p_adj_Zika, Regulation_Zika) %>%
    arrange(ensembl)
  
  # Save full table (all genes)
  write_tab(virus_wide, "DEG_Virus_Comparison.csv")
  
  # Save filtered table: significant in at least one virus (p_adj < 0.05, |log2FC| >= 0.25)
  sig_wide <- virus_wide %>%
    filter((p_adj_Dengue < 0.05 & abs(avg_log2FC_Dengue) >= 0.25) |
             (p_adj_Zika < 0.05 & abs(avg_log2FC_Zika) >= 0.25))
  write_tab(sig_wide, "DEG_Virus_Comparison_Sig.csv")
  cat(sprintf(" -> Saved merged DEG comparison: %d total, %d significant genes.\n", nrow(virus_wide), nrow(sig_wide)))
  
  mat_for_save <- tryCatch(mat_log, error = function(e) mat_scaled)
  if (length(zika_sig) > 0) {
    zika_mat <- mat_for_save[rownames(mat_for_save) %in% zika_sig, , drop = FALSE]
    saveRDS(zika_mat, obj("Phase_D_ZIKV_Sig_Matrix.rds"))
    cat(" -> Saved ZIKV sig DEG matrix (%d genes).\n", nrow(zika_mat))
  }
  if (length(dengue_sig) > 0) {
    dengue_mat <- mat_for_save[rownames(mat_for_save) %in% dengue_sig, , drop = FALSE]
    saveRDS(dengue_mat, obj("Phase_D_DENV_Sig_Matrix.rds"))
    cat(" -> Saved DENV sig DEG matrix (%d genes).\n", nrow(dengue_mat))
  }
  if (length(common_sig) > 0) {
    common_mat <- mat_for_save[rownames(mat_for_save) %in% common_sig, , drop = FALSE]
    saveRDS(common_mat, obj("Phase_D_Common_Sig_Matrix.rds"))
    cat(" -> Saved Common sig DEG matrix (%d genes).\n", nrow(common_mat))
  }
}

# ==============================================================================
# Persist clustering for Phase E
# ==============================================================================
saveRDS(list(pca = pca, n_pc = n_pc, clusters = clusters, res_pick = res_pick,
             stability = stab_df, embedding = emb, emb_name = emb_name,
             use_harmony = USE_HARMONY, hm_pc_keep = hm_pc_keep,
             pcs_source = pcs_source, plate_driven_pcs = n_plate_driven),
        obj("Phase_D_clustering.rds"))

cat("\n======================================================================\n")
cat("PHASE D MAIN COMPLETE\n")
cat(sprintf(" Clustering: %d clusters (%s)\n", nlevels(clusters), pcs_source))
cat(sprintf(" Within-plate DE: %d strata, %d genes tested\n",
            nlevels(factor(de_all$stratum)), nrow(de_all)))
cat(sprintf(" Timepoint DE: %d conditions tested\n", length(all_timepoint_de)))
cat(sprintf(" vRNA continuous: %d genes across %d viruses\n",
            nrow(vrna), length(unique(vrna$virus))))
cat(sprintf(" DEG matrices: ZIKV + DENV + Common saved to %s\n", dir_obj))
cat(" Objects: ", dir_obj, "\n")
cat(" Figures/tables: ", dir_fig, "\n")
cat("======================================================================\n")