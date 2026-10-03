# ==============================================================================
# PHASE D - SUPPLEMENTARY ANALYSIS 
# ==============================================================================
# Cluster characterisation, cross-virus concordance (H2 test), cluster markers.
#
# Inputs:
# results/phaseD/objects/Phase_D_clustering.rds
# results/phaseD/tables/DE_WithinPlate_Wilcoxon.csv
# results/phaseD/tables/vRNA_Continuous_Association.csv
# results/phaseC/objects/Phase_C_log_matrix.rds
# results/phaseC/objects/Phase_C_scaled_matrix.rds
# results/phaseC/objects/Phase_C_cell_qc.rds
# results/phaseB/objects/Phase_B_counts_active.rds
# results/phaseA/objects/Phase_A_gene_map.rds
# results/phaseB/objects/Phase_B_gene_sets.rds
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
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
})

set.seed(42)

# dplyr mask guard
.verbs <- c("select", "filter", "rename", "count", "slice", "desc", "first")
.owner <- vapply(.verbs, function(f) find(f)[1], character(1))
.masked <- .owner[!.owner %in% c("package:dplyr", "package:tidyr")]
if (length(.masked))
  stop("dplyr verbs are masked: ",
       paste(sprintf("%s() <- %s", names(.masked), .masked), collapse = ", "),
       "\nRestart R and source in a clean session.")
rm(.verbs, .owner, .masked)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
data_dir <- "."
setwd(data_dir)

results_root <- file.path(data_dir, "results", "phaseD")
dir_fig <- file.path(results_root, "figures")
dir_tab <- file.path(results_root, "tables")
dir_supp <- file.path(results_root, "supplementary")
dir.create(dir_supp, recursive = TRUE, showWarnings = FALSE)

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

show_and_save <- function(p, file, width = 8, height = 5) {
  p <- p + theme_pub()
  print(p)
  ggsave(file.path(dir_fig, paste0(file, ".pdf")), p, width = width, height = height)
  ggsave(file.path(dir_fig, paste0(file, ".png")), p, width = width, height = height, dpi = 600)
  cat(sprintf(" [plot] %s.pdf / .png\n", file))
  invisible(p)
}
supp_and_save <- function(p, file, width = 8, height = 5) {
  p <- p + theme_pub()
  print(p)
  ggsave(file.path(dir_supp, paste0(file, ".pdf")), p, width = width, height = height)
  ggsave(file.path(dir_supp, paste0(file, ".png")), p, width = width, height = height, dpi = 600)
  cat(sprintf(" [supp] %s.pdf / .png\n", file))
  invisible(p)
}
write_tab <- function(x, name) {
  write.csv(x, file.path(dir_tab, name), row.names = FALSE)
  cat(sprintf(" [table] %s\n", name))
}

cat("======================================================================\n")
cat("PHASE D SUPPLEMENTARY \n")
cat("======================================================================\n")

# ==============================================================================
# Load inputs
# ==============================================================================
mat_log <- readRDS(find_file("Phase_C_log_matrix.rds"))
mat_scaled <- readRDS(find_file("Phase_C_scaled_matrix.rds"))
qc <- readRDS(find_file("Phase_C_cell_qc.rds"))
gene_map <- readRDS(find_file("Phase_A_gene_map.rds"))
gene_sets <- tryCatch(readRDS(find_file("Phase_B_gene_sets.rds")), error = function(e) NULL)
mat_counts <- readRDS(find_file("Phase_B_counts_active.rds"))
clustering <- tryCatch(readRDS(find_file("Phase_D_clustering.rds")), error = function(e) NULL)

ensg_col <- intersect(c("ENSG", "ENSG_Original", "ENSG_Clean"), names(gene_map))[1]
sym_of <- function(ids) {
  if (is.null(gene_map) || is.na(ensg_col)) return(ids)
  s <- gene_map$Symbol[match(ids, gene_map[[ensg_col]])]
  ifelse(is.na(s), ids, s)
}

cat(sprintf("Loaded: mat_log %d x %d | mat_scaled %d x %d | QC %d rows\n",
            nrow(mat_log), ncol(mat_log), nrow(mat_scaled), ncol(mat_scaled), nrow(qc)))
if (!is.null(clustering))
  cat(sprintf("Clusters: %d levels\n", nlevels(clustering$clusters)))

# Load main Phase D results
de_all <- tryCatch(
  read.csv(find_file("DE_WithinPlate_Wilcoxon.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
vrna <- tryCatch(
  read.csv(find_file("vRNA_Continuous_Association.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)

# Load  DEG tables (primary significant gene lists)
sig_zika <- tryCatch(
  read.csv(find_file("ZIKV_Sig_Genes.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
sig_denv <- tryCatch(
  read.csv(find_file("DENV_Sig_Genes.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
sig_common <- tryCatch(
  read.csv(find_file("Common_Sig_Genes.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
sig_all <- tryCatch(
  read.csv(find_file("All_Timepoint_DEGs.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)

if (!is.null(de_all))
  cat(sprintf("DE table: %d rows, %d strata\n", nrow(de_all), nlevels(factor(de_all$stratum))))
if (!is.null(vrna))
  cat(sprintf("vRNA association: %d rows, %d strata\n", nrow(vrna), nlevels(factor(vrna$stratum))))
if (!is.null(sig_zika))
  cat(sprintf(" DEGs: ZIKV %d | DENV %d | Common %d\n",
              nrow(sig_zika), nrow(sig_denv), nrow(sig_common)))

# ==============================================================================
# SUPP_1 : Augment QC table
# ==============================================================================
cat("\n[Supp_1] Augmenting QC table with log_viral and pct_ribo\n")

qc$log_viral <- log10(qc$n_virus_molecules + 1)

ribo_genes <- character(0)
if (!is.null(gene_sets) && !is.null(gene_sets$ribosomal)) {
  ribo_genes <- gene_sets$ribosomal
} else {
  ribo_genes <- unique(gene_map[[ensg_col]][!is.na(gene_map$Symbol) &
                                              str_detect(gene_map$Symbol, "^RP[LS]\\d")])
}
ribo_in_mat <- intersect(ribo_genes, rownames(mat_counts))
cat(sprintf(" Ribosomal genes: %d available, %d in active matrix\n",
            length(ribo_genes), length(ribo_in_mat)))

if (length(ribo_in_mat) > 0) {
  ribo_counts <- colSums(mat_counts[ribo_in_mat, , drop = FALSE])
  lib <- colSums(mat_counts); lib[lib == 0] <- 1
  qc$pct_ribo <- 100 * ribo_counts / lib
} else {
  cat(" WARNING: no ribosomal genes in matrix; pct_ribo set to NA\n")
  qc$pct_ribo <- NA_real_
}

if (!is.null(clustering)) {
  qc$cluster <- clustering$clusters
  qc$cluster_name <- clustering$cluster_names %||% as.character(clustering$clusters)
}

cat(sprintf(" log_viral : range [%.3f, %.3f], mean = %.3f\n",
            min(qc$log_viral, na.rm = TRUE), max(qc$log_viral, na.rm = TRUE),
            mean(qc$log_viral, na.rm = TRUE)))
cat(sprintf(" pct_ribo : range [%.3f, %.3f], mean = %.3f\n",
            min(qc$pct_ribo, na.rm = TRUE), max(qc$pct_ribo, na.rm = TRUE),
            mean(qc$pct_ribo, na.rm = TRUE)))

saveRDS(qc, file.path(results_root, "objects", "Phase_D_qc_augmented.rds"))
write_tab(qc %>% select(Sample, virus, time_h, moi, batch, host_counts, ercc_counts,
                        n_genes_detected, pct_mt, pct_ercc, pct_virus, pct_ribo,
                        log_viral, infection_status, cluster, cluster_name),
          "Supp1_QC_Augmented.csv")

# ==============================================================================
# SUPP_2 : PC-plate correlation check ( verification)
# ==============================================================================
cat("\n[Supp_2] PC-plate correlation check (verification)\n")

if (!is.null(clustering) && !is.null(clustering$pca)) {
  pca <- clustering$pca
} else {
  pca <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE)
}
n_pc_used <- if (!is.null(clustering) && !is.null(clustering$n_pc)) clustering$n_pc else min(50, ncol(pca$x))

plate_num <- as.numeric(factor(qc$batch))

pc_cov <- data.frame(
  PC = paste0("PC", seq_len(min(50, n_pc_used))),
  var_pct = 100 * pca$sdev[seq_len(min(50, n_pc_used))]^2 / sum(pca$sdev^2)
)

corr_mat <- sapply(seq_len(min(50, n_pc_used)), function(k) {
  pc_k <- pca$x[, k]
  c(
    plate = suppressWarnings(cor(pc_k, plate_num, method = "spearman")),
    nCount = suppressWarnings(cor(pc_k, qc$host_counts, method = "spearman")),
    log_viral = suppressWarnings(cor(pc_k, qc$log_viral, method = "spearman")),
    S_score = if (is.null(qc$S_score)) NA_real_ else suppressWarnings(cor(pc_k, qc$S_score, method = "spearman"))
  )
})
corr_mat <- t(corr_mat)
pc_cov <- cbind(pc_cov, as.data.frame(corr_mat))
pc_cov$plate_abs <- abs(pc_cov$plate)
pc_cov$nCount_abs <- abs(pc_cov$nCount)
pc_cov$log_viral_abs <- abs(pc_cov$log_viral)
pc_cov$S_score_abs <- abs(pc_cov$S_score)

pc_cov$plate_driven <- pc_cov$plate_abs > 0.3
cat(sprintf(" PCs with |r(plate)| > 0.3: %d / %d\n",
            sum(pc_cov$plate_driven), nrow(pc_cov)))
if (any(pc_cov$plate_driven)) {
  cat(sprintf(" Plate-driven PCs: %s\n",
              paste(pc_cov$PC[pc_cov$plate_driven], collapse = ", ")))
}

write_tab(pc_cov, "Supp2_PC_Plate_Correlation.csv")

corr_long <- pc_cov %>%
  select(PC, var_pct, plate, nCount, log_viral, S_score) %>%
  pivot_longer(-c(PC, var_pct), names_to = "covariate", values_to = "rho")

p_pc_corr <- ggplot(corr_long,
                    aes(x = covariate, y = reorder(PC, as.numeric(str_remove(PC, "PC"))), fill = rho)) +
  geom_tile(colour = "grey80", linewidth = 0.15) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "Spearman rho") +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  labs(title = "PC-covariate correlations (Spearman rho)",
       subtitle = paste0(
         sum(pc_cov$plate_driven), " of ", nrow(pc_cov),
         " PCs have |r(plate)| > 0.3. ",
         if (any(pc_cov$plate_driven)) "Plate-driven PCs trigger integration." else
           "No plate-driven PCs; standard clustering is adequate."),
       x = "Covariate", y = "Principal component")
show_and_save(p_pc_corr, "Supp2_PC_Correlation_Heatmap", width = 7, height = 10)

# ==============================================================================
# SUPP_3 : Cluster characterisation ()
# ==============================================================================
cat("\n[Supp_3] Cluster characterisation\n")

if (is.null(clustering)) {
  cat(" WARNING: clustering object not found; skipping Supp_3.\n")
} else {
  clusters <- clustering$clusters
  qc$cluster <- clusters
  n_clust <- nlevels(clusters)
  cat(sprintf(" %d clusters, %d cells\n", n_clust, length(clusters)))
  
  # 12a. Viral load per cluster
  p_vcl <- ggplot(qc, aes(x = cluster, y = log_viral, fill = cluster)) +
    geom_violin(scale = "width", colour = "#4D4D4D", linewidth = 0.3) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", linewidth = 0.3) +
    scale_fill_viridis_d(guide = "none") +
    labs(title = "Viral load per cluster",
         x = "Cluster", y = expression(log[10]~"(viral molecules + 1)"))
  supp_and_save(p_vcl, "Supp3a_ViralLoad_ByCluster", width = 9, height = 5)
  
  # 12b. Infection state composition per cluster
  cl_state <- qc %>%
    filter(!is.na(infection_status)) %>%
    count(cluster, infection_status) %>%
    group_by(cluster) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  p_csb <- ggplot(cl_state, aes(x = cluster, y = prop, fill = infection_status)) +
    geom_col(colour = "black", linewidth = 0.2, width = 0.7) +
    scale_fill_manual(values = c(Unexposed = "#A6A6A6", Bystander = "#0072B2",
                                 Infected = "#D55E00")) +
    scale_y_continuous(labels = percent_format()) +
    labs(title = "Infection state composition per cluster",
         x = "Cluster", y = "Proportion", fill = "Infection status")
  supp_and_save(p_csb, "Supp3b_StateComposition_ByCluster", width = 9, height = 5)
  
  # 12c. Module scores: UPR, apoptosis, cell cycle (ISG removed â€” canonical antiviral responders)
  upr_genes_raw <- c("HSPA5", "DDIT3", "ATF4", "XBP1", "ATF6", "EIF2AK3", "ERN1",
                     "DNAJB9", "HYOU1", "CALR")
  apo_genes_raw <- c("CASP8", "CASP9", "CASP3", "FAS", "FASLG", "BAX", "BAK1",
                     "BCL2", "MCL1", "BID", "BIM", "PUMA", "NOXA", "XIAP",
                     "BIRC2", "BIRC3", "CFLAR", "TRAIL", "TNFRSF10A", "TNFRSF10B")
  s_genes_raw <- c("MCM5", "PCNA", "ORC1", "DTL", "GINS2", "MCM4", "FEN1", "MCM2",
                   "RRM1", "UNG", "EXO1", "TIPIN", "BLM", "USP1", "SLBP",
                   "CHAF1B", "BRIP1", "E2F8")
  g2m_genes_raw <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
                     "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
                     "CCNB2", "CKAP2L")
  
  map_to_ensg <- function(symbols) {
    if (is.null(gene_map) || is.na(ensg_col)) return(character(0))
    matched <- gene_map[[ensg_col]][match(toupper(symbols), toupper(gene_map$Symbol))]
    matched[!is.na(matched)]
  }
  
  gene_sets_list <- list(
    UPR = intersect(map_to_ensg(upr_genes_raw), rownames(mat_log)),
    Apoptosis = intersect(map_to_ensg(apo_genes_raw), rownames(mat_log)),
    CellCycle_S = intersect(map_to_ensg(s_genes_raw), rownames(mat_log)),
    CellCycle_G2M = intersect(map_to_ensg(g2m_genes_raw), rownames(mat_log))
  )
  
  cat(" Module gene counts:\n")
  for (nm in names(gene_sets_list))
    cat(sprintf(" %-15s %3d genes\n", nm, length(gene_sets_list[[nm]])))
  
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
  
  score_df <- data.frame(Sample = colnames(mat_log), stringsAsFactors = FALSE)
  for (nm in names(gene_sets_list)) {
    if (length(gene_sets_list[[nm]]) > 0) {
      score_df[[nm]] <- as.numeric(module_score(mat_log, gene_sets_list[[nm]]))
    }
  }
  
  qc <- left_join(qc, score_df, by = "Sample")
  
  # Dot plot: module scores by cluster
  score_long <- score_df %>%
    pivot_longer(-Sample, names_to = "module", values_to = "score") %>%
    left_join(qc %>% select(Sample, cluster), by = "Sample")
  
  dot_stats <- score_long %>%
    group_by(cluster, module) %>%
    summarise(avg_score = mean(score, na.rm = TRUE),
              pct_pos = mean(score > 0, na.rm = TRUE),
              .groups = "drop")
  
  p_dot <- ggplot(dot_stats, aes(x = module, y = factor(cluster), fill = avg_score)) +
    geom_tile(colour = "grey80", linewidth = 0.15) +
    geom_text(aes(label = sprintf("%.2f", avg_score)), size = 2.5) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, name = "Mean\nmodule score") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = "Module scores per cluster",
         x = "Module", y = "Cluster")
  supp_and_save(p_dot, "Supp3c_ModuleScores_DotPlot", width = 9, height = 6)
  
  # Violin per module
  for (nm in names(gene_sets_list)) {
    if (!(nm %in% names(score_df))) next
    p_mod <- ggplot(score_long %>% filter(module == nm),
                    aes(x = cluster, y = score, fill = cluster)) +
      geom_violin(scale = "width", colour = "#4D4D4D", linewidth = 0.3) +
      scale_fill_viridis_d(guide = "none") +
      labs(title = sprintf("Supp_3c: %s module score by cluster", nm),
           x = "Cluster", y = "Module score")
    supp_and_save(p_mod, sprintf("Supp3c_ModuleScore_%s", nm), width = 9, height = 5)
  }
  
  # 12d. Descriptive cluster names
  cl_desc <- dot_stats %>%
    pivot_wider(names_from = module, values_from = c(avg_score, pct_pos)) %>%
    mutate(
      upr_hi = avg_score_UPR > quantile(avg_score_UPR, 0.75, na.rm = TRUE),
      apo_hi = avg_score_Apoptosis > quantile(avg_score_Apoptosis, 0.75, na.rm = TRUE),
      cc_hi = (avg_score_CellCycle_S + avg_score_CellCycle_G2M) >
        quantile(avg_score_CellCycle_S + avg_score_CellCycle_G2M, 0.75, na.rm = TRUE),
      # Secondary prominence: apoptosis is visible in cluster 10 (0.461) even though
      # it falls below the top-quartile threshold â€” use median split as a backstop.
      apo_med = avg_score_Apoptosis > median(avg_score_Apoptosis, na.rm = TRUE)
    )
  
  vl_by_cl <- qc %>%
    group_by(cluster) %>%
    summarise(vl_med = median(log_viral, na.rm = TRUE),
              vl_q3 = quantile(log_viral, 0.75, na.rm = TRUE),
              .groups = "drop")
  cl_desc <- left_join(cl_desc, vl_by_cl, by = "cluster")
  cl_desc$virus_hi <- cl_desc$vl_q3 > quantile(cl_desc$vl_q3, 0.75, na.rm = TRUE)
  
  cl_desc <- cl_desc %>%
    mutate(
      name = case_when(
        upr_hi & apo_med ~ "UPR+Apoptosis",
        upr_hi ~ "UPR-high",
        apo_hi ~ "Apoptosis-high",
        cc_hi ~ "CellCycle-high",
        TRUE ~ "Basal"
      ),
      name = ifelse(virus_hi & vl_med > 1, paste0(name, "+Viral"), name)
    )
  
  write_tab(cl_desc %>%
              mutate(n_cells = NA_integer_) %>%
              select(cluster, n_cells, starts_with("avg_score_"), name),
            "Supp3d_Cluster_Descriptions.csv")
  
  cat(" Cluster descriptions:\n")
  print(as.data.frame(cl_desc %>% select(cluster, name, vl_med,
  avg_score_UPR, avg_score_Apoptosis,
  avg_score_CellCycle_S, avg_score_CellCycle_G2M)),
  digits = 2, row.names = FALSE)
  
  cl_names <- setNames(cl_desc$name, cl_desc$cluster)
  if (!is.null(clustering$clusters)) {
    qc$cluster_name <- cl_names[as.character(clusters)]
  }
  saveRDS(qc, file.path(results_root, "objects", "Phase_D_qc_augmented.rds"))
  
  # 12e. clustree: resolution sweep trajectory
  res_sens <- tryCatch(
    read.csv(find_file("Resolution_Sensitivity.csv"), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(res_sens) && nrow(res_sens) > 1) {
    cat(sprintf(" -> Building clustree from %d resolution points...\n", nrow(res_sens)))
    
    # Reconstruct cluster assignments at each resolution by re-running Louvain
    # on the SNN graph (rebuilt from Harmony PCs stored in clustering$embedding)
    pcs_for_clustree <- if (clustering$use_harmony) {
      # We don't have raw Harmony PCs saved, use PCA as approximation
      # Rebuild SNN from PCA (acceptable for clustree visualization)
      suppressWarnings(tryCatch({
        pca_full <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE, rank. = clustering$n_pc)
        pca_full$x[, seq_len(clustering$hm_pc_keep), drop = FALSE]
      }, error = function(e) NULL))
    } else {
      suppressWarnings(tryCatch({
        pca_full <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE, rank. = clustering$n_pc)
        pca_full$x[, seq_len(clustering$n_pc), drop = FALSE]
      }, error = function(e) NULL))
    }
    
    if (!is.null(pcs_for_clustree)) {
      knn_cl <- FNN::get.knn(pcs_for_clustree, k = 20)$nn.index
      nr <- nrow(pcs_for_clustree)
      Mb <- sparseMatrix(i = rep(seq_len(nr), each = 20), j = as.vector(t(knn_cl)), x = 1, dims = c(nr, nr))
      diag(Mb) <- 1
      it <- as(Mb %*% t(Mb), "TsparseMatrix")
      jj <- it@x / (2 * 21 - it@x)
      kp <- jj >= (1/15) & it@i != it@j
      e <- data.frame(from = it@i[kp] + 1, to = it@j[kp] + 1, weight = jj[kp])
      e <- e[e$from < e$to, ]
      g_cl <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = seq_len(nr)))
      
      # Assign clusters at each resolution
      res_vals <- res_sens$resolution
      clust_assign <- data.frame(cell = seq_len(nr), stringsAsFactors = FALSE)
      for (r in res_vals) {
        cl_r <- igraph::membership(igraph::cluster_louvain(g_cl, weights = E(g_cl)$weight, resolution = r))
        clust_assign[[sprintf("res_%.1f", r)]] <- as.character(cl_r)
      }
      
      # Clustree-style heatmap: for each resolution, show cluster sizes
      clustree_long <- clust_assign %>%
        pivot_longer(-cell, names_to = "resolution", values_to = "cluster") %>%
        mutate(resolution = as.numeric(str_remove(resolution, "res_"))) %>%
        count(resolution, cluster) %>%
        group_by(resolution) %>%
        mutate(prop = n / sum(n)) %>%
        ungroup()
      
      p_clustree <- ggplot(clustree_long, aes(x = resolution, y = cluster, fill = prop)) +
        geom_tile(colour = "grey80", linewidth = 0.15) +
        scale_fill_viridis_c(option = "D", name = "Proportion", labels = percent_format()) +
        scale_x_continuous(breaks = res_vals) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Clustree â€” cluster composition across resolutions",
             x = "Louvain resolution", y = "Cluster")
      supp_and_save(p_clustree, "Supp3e_Clustree", width = 10, height = 6)
      
      # Cluster count vs resolution (concise summary)
      p_res_sum <- ggplot(res_sens, aes(x = resolution, y = n_clusters)) +
        geom_line(colour = "#0072B2", linewidth = 1) +
        geom_point(aes(size = modularity), colour = "#0072B2") +
        geom_text(aes(label = n_clusters), vjust = -1.2, size = 3) +
        scale_x_continuous(breaks = res_vals) +
        labs(title = "Cluster count vs resolution",
             subtitle = sprintf("Selected resolution = %.1f (modularity = %.3f)",
                                res_sens$resolution[which.max(res_sens$modularity)],
                                max(res_sens$modularity)),
             x = "Louvain resolution", y = "Number of clusters",
             size = "Modularity")
      supp_and_save(p_res_sum, "Supp3e_Resolution_Sensitivity", width = 8, height = 5)
      
    } else {
      cat(" Could not rebuild SNN graph for clustree; skipping.\n")
    }
  }
}

cat("\n[Supp_3-F]  DEG summary\n")

if (!is.null(sig_zika) && !is.null(sig_denv)) {
  cat(sprintf(" ZIKV: %d | DENV: %d | Common: %d\n",
              nrow(sig_zika), nrow(sig_denv),
              if (!is.null(sig_common)) nrow(sig_common) else 0))
  
  venn_summary <- data.frame(
    comparison = c("ZIKV only", "DENV only", "Shared (both viruses)"),
    n_genes = c(
      nrow(sig_zika) - if (!is.null(sig_common)) nrow(sig_common) else 0,
      nrow(sig_denv) - if (!is.null(sig_common)) nrow(sig_common) else 0,
      if (!is.null(sig_common)) nrow(sig_common) else 0
    )
  ) %>%
    mutate(pct = 100 * n_genes / sum(n_genes))
  write_tab(venn_summary, "DEG_Venn_Summary.csv")
  
  p_venn <- ggplot(venn_summary, aes(x = comparison, y = n_genes, fill = comparison)) +
    geom_col(colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d (%.0f%%)", n_genes, pct)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("ZIKV only" = "#0072B2",
                                 "DENV only" = "#D55E00",
                                 "Shared (both viruses)" = "forestgreen")) + 
    theme(legend.position = "none") +
    labs(title = "DEG summary (Infected vs Bystander)",
         subtitle = sprintf("ZIKV %d | DENV %d | Shared %d",
                            nrow(sig_zika), nrow(sig_denv),
                            if (!is.null(sig_common)) nrow(sig_common) else 0),
         x = NULL, y = "Number of significant DEGs")
  supp_and_save(p_venn, "DEG_Summary_Barplot", width = 8, height = 6)
  
  if (!is.null(sig_all)) {
    cond_breakdown <- sig_all %>%
      filter(is_sig) %>%
      group_by(virus, time_h, condition) %>%
      summarise(n_sig = n(), .groups = "drop") %>%
      arrange(desc(n_sig))
    write_tab(cond_breakdown, "DEG_ByCondition.csv")
    cat(sprintf(" Per-condition breakdown saved (%d conditions with sig DEGs)\n", nrow(cond_breakdown)))
    
    p_cond <- ggplot(cond_breakdown, aes(x = factor(condition, levels = rev(condition)), y = n_sig, fill = virus)) +
      geom_col(colour = "black", linewidth = 0.2) +
      scale_fill_manual(values = c(zika = "#0072B2", dengue = "#D55E00")) +
      coord_flip() +
      labs(title = "Significant DEGs per condition",
           x = "Condition", y = "Number of DEGs", fill = "Virus")
    supp_and_save(p_cond, "DEG_ByCondition_Barplot", width = 9, height = 6)
  }
}

# ==============================================================================
# SUPP_4 : Cross-virus concordance ( - H2 test)
# ==============================================================================
cat("\n[Supp_4] Cross-virus concordance (ZIKV rho vs DENV rho)\n")

if (is.null(vrna) || nrow(vrna) == 0) {
  cat(" WARNING: no vRNA association data found; skipping Supp_4.\n")
} else {
  # Determine virus from stratum or virus column
  if ("virus" %in% names(vrna)) {
    vr_wide <- vrna %>%
      select(ensembl, symbol, virus, rho, p_adj) %>%
      pivot_wider(names_from = virus, values_from = c(rho, p_adj),
                  names_sep = "_")
  } else {
    # Parse from stratum column
    vr_wide <- vrna %>%
      mutate(virus = ifelse(grepl("zika", stratum, ignore.case = TRUE), "zika",
                            ifelse(grepl("dengue|denv", stratum, ignore.case = TRUE), "dengue", "unknown"))) %>%
      select(ensembl, symbol, virus, rho, p_adj) %>%
      pivot_wider(names_from = virus, values_from = c(rho, p_adj), names_sep = "_")
  }
  
  if (("rho_zika" %in% names(vr_wide) || "rho_dengue" %in% names(vr_wide)) &&
      any(grepl("^rho_", names(vr_wide)))) {
    vr_wide <- vr_wide %>%
      select(ensembl, symbol, starts_with("rho_")) %>%
      select(-ends_with("_NA")) %>%
      distinct()
    
    # Ensure we have both viruses
    zika_col <- grep("^rhosig_zika$", names(vr_wide), value = TRUE)
    denv_col <- grep("^rho_dengue$|^rhosig_denv$", names(vr_wide), value = TRUE)
    
    if (length(zika_col) > 0 && length(denv_col) > 0) {
      vr_wide <- vr_wide %>%
        rename(rho_zika = !!sym(zika_col[1]), rho_dengue = !!sym(denv_col[1])) %>%
        filter(is.finite(rhosig_zika), is.finite(rho_dengue))
      
      cat(sprintf(" Genes with both ZIKV and DENV rho: %d\n", nrow(vr_wide)))
      
      p_conc <- ggplot(vr_wide, aes(x = rhosig_zika, y = rho_dengue)) +
        geom_point(alpha = 0.3, size = 0.8, colour = "#0072B2") +
        geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
        geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
        geom_density_2d(colour = "#D55E00", linewidth = 0.4) +
        labs(title = "Cross-virus concordance (H2 test)",
             subtitle = sprintf(
               "Spearman rho of the two rho vectors = %.3f (p = %.3g). Each point = one gene.",
               cor(vr_wide$rhosig_zika, vr_wide$rho_dengue, method = "spearman"),
               suppressWarnings(cor.test(vr_wide$rhosig_zika, vr_wide$rho_dengue,
                                         method = "spearman"))$p.value),
             x = expression(paste("ZIKV: ", rho, " (expr vs viral load)")),
             y = expression(paste("DENV: ", rho, " (expr vs viral load)")))
      supp_and_save(p_conc, "Supp4_CrossVirus_Concordance", width = 8, height = 7)
      
      # Label top discordant genes
      vr_wide$rho_diff <- vr_wide$rho_zika - vr_wide$rho_dengue
      vr_wide$rho_abs_sum <- abs(vr_wide$rhosig_zika) + abs(vr_wide$rho_dengue)
      top_disc <- vr_wide %>%
        filter(rho_abs_sum > 0.5) %>%
        arrange(desc(abs(rho_diff))) %>%
        slice_head(n = 15)
      
      p_conc_lab <- p_conc +
        geom_text_repel(data = top_disc, aes(label = symbol),
                        size = 2.5, max.overlaps = Inf, min.segment.length = 0,
                        colour = "#D55E00")
      supp_and_save(p_conc_lab, "Supp4_CrossVirus_Concordance_Labelled", width = 9, height = 8)
      
      # Concordance stats
      conc_stats <- data.frame(
        metric = c("n_genes", "spearman_rho", "spearman_p",
                   "mean_abs_rho_zika", "mean_abs_rho_dengue",
                   "pct_concordant", "pct_strong_concordant", "pct_discordant",
                   "mean_rho_difference"),
        value = c(
          nrow(vr_wide),
          cor(vr_wide$rhosig_zika, vr_wide$rho_dengue, method = "spearman"),
          suppressWarnings(cor.test(vr_wide$rhosig_zika, vr_wide$rho_dengue,
                                    method = "spearman"))$p.value,
          mean(abs(vr_wide$rhosig_zika), na.rm = TRUE),
          mean(abs(vr_wide$rho_dengue), na.rm = TRUE),
          100 * mean(sign(vr_wide$rhosig_zika) == sign(vr_wide$rho_dengue), na.rm = TRUE),
          100 * mean(sign(vr_wide$rhosig_zika) == sign(vr_wide$rho_dengue) &
                       abs(vr_wide$rhosig_zika) > 0.3 & abs(vr_wide$rho_dengue) > 0.3, na.rm = TRUE),
          100 * mean(sign(vr_wide$rhosig_zika) != sign(vr_wide$rho_dengue) &
                       abs(vr_wide$rhosig_zika) > 0.2 & abs(vr_wide$rho_dengue) > 0.2, na.rm = TRUE),
          mean(abs(vr_wide$rho_zika - vr_wide$rho_dengue), na.rm = TRUE)
        ),
        stringsAsFactors = FALSE
      )
      write_tab(conc_stats, "Supp4_Concordance_Stats.csv")
      
      # Quadrant breakdown
      vr_wide$quadrant <- case_when(
        vr_wide$rho_zika > 0 & vr_wide$rho_dengue > 0 ~ "Both +",
        vr_wide$rho_zika < 0 & vr_wide$rho_dengue < 0 ~ "Both -",
        vr_wide$rho_zika > 0 & vr_wide$rho_dengue < 0 ~ "ZIKV+ / DENV-",
        TRUE ~ "ZIKV- / DENV+"
      )
      quad_tab <- vr_wide %>% count(quadrant) %>%
        mutate(pct = 100 * n / sum(n))
      write_tab(quad_tab, "Supp4_Quadrant_Breakdown.csv")
      
      cat(sprintf(" Spearman rho (ZIKV rho vs DENV rho) = %.3f (p = %.3g)\n",
                  conc_stats$value[2], conc_stats$value[3]))
      cat(sprintf(" Concordant: %.1f%% | Discordant: %.1f%%\n",
                  conc_stats$value[6], conc_stats$value[8]))
      
      # Permutation test: is concordance greater than expected by chance?
      n_perm <- 500
      cat(" -> Running permutation test (%d iterations, matrix-vectorized)...\n", n_perm)
      exposed_idx <- which(qc$moi > 0 & qc$infection_status %in% c("Infected", "Bystander"))
      virus_labels <- qc$virus[exposed_idx]
      null_rhos <- numeric(n_perm)
      null_ok <- logical(n_perm)
      
      # Pre-extract exposed-cell expression and rank-transform ONCE (Spearman = Pearson on ranks)
      mat_exp <- mat_log[, exposed_idx, drop = FALSE]
      lv_exp <- qc$log_viral[exposed_idx]
      n_cells_exp <- length(exposed_idx)
      
      # Rank-transform expression matrix: genes x cells (done once, not per iteration)
      expr_ranked <- apply(mat_exp, 1, rank, ties.method = "average")
      lv_ranked <- rank(lv_exp, ties.method = "average")
      
      for (p in seq_len(n_perm)) {
        perm_virus <- sample(virus_labels)
        msk_z <- perm_virus == "zika"
        msk_d <- perm_virus == "denv"
        if (sum(msk_z) < 10 || sum(msk_d) < 10) { null_rhos[p] <- NA; null_ok[p] <- FALSE; next }
        
        # Vectorized Spearman: cor(t(ranked_genes), ranked_lv) â€” single C call per virus
        rho_z <- suppressWarnings(cor(t(expr_ranked[, msk_z, drop = FALSE]),
                                      lv_ranked[msk_z], method = "pearson"))
        rho_d <- suppressWarnings(cor(t(expr_ranked[, msk_d, drop = FALSE]),
                                      lv_ranked[msk_d], method = "pearson"))
        
        ok <- is.finite(rho_z) & is.finite(rho_d)
        null_ok[p] <- (sum(ok) > 10)
        if (sum(ok) > 10) {
          null_rhos[p] <- cor(rho_z[ok], rho_d[ok], method = "spearman")
        } else {
          null_rhos[p] <- NA_real_
        }
        if (p %% 100 == 0) cat(sprintf(" permutation %d/%d...\n", p, n_perm))
      }
      
      obs_rho <- conc_stats$value[2]
      valid_rhos <- null_rhos[null_ok & !is.na(null_rhos)]
      
      if (length(valid_rhos) >= 10) {
        emp_p <- mean(valid_rhos >= obs_rho)
        perm_stats <- data.frame(
          metric = c("observed_rho", "permutation_p", "null_mean", "null_sd",
                     "n_perm_valid", "n_perm_na"),
          value = c(
            obs_rho, emp_p, mean(valid_rhos), sd(valid_rhos),
            length(valid_rhos), sum(!null_ok, na.rm = TRUE)
          ),
          stringsAsFactors = FALSE
        )
        write_tab(perm_stats, "Supp4_Permutation_Test.csv")
        
        cat(sprintf(" Permutation test: observed rho = %.3f, empirical p = %.3g (null mean = %.3f, sd = %.3f, n_valid = %d)\n",
                    obs_rho, emp_p, mean(valid_rhos), sd(valid_rhos), length(valid_rhos)))
        
        # Null distribution plot
        p_perm <- ggplot(data.frame(null = valid_rhos), aes(x = null)) +
          geom_histogram(bins = 40, fill = "#0072B2", colour = "black", linewidth = 0.2) +
          geom_vline(xintercept = obs_rho, colour = "#D55E00", linewidth = 1.2) +
          annotate("text", x = obs_rho, y = Inf, label = sprintf("Observed = %.3f", obs_rho),
                   vjust = -1.5, colour = "#D55E00", size = 3.5) + 
          labs(title = "Permutation test for concordance",
               subtitle = sprintf("Null distribution (n = %d), observed rho = %.3f, emp. p = %.3g",
                                  length(valid_rhos), obs_rho, emp_p),
               x = "Spearman rho (ZIKV rho vs DENV rho) under null", y = "Frequency")
        supp_and_save(p_perm, "Supp4_Permutation_NullDist", width = 8, height = 5)
      } else {
        cat(sprintf(" WARNING: only %d valid permutations; skipping permutation test.\n", length(valid_rhos)))
        perm_stats <- data.frame(
          metric = c("observed_rho", "permutation_p", "null_mean", "null_sd",
                     "n_perm_valid", "n_perm_na"),
          value = c(obs_rho, NA, NA, NA, length(valid_rhos), sum(!null_ok, na.rm = TRUE)),
          stringsAsFactors = FALSE
        )
        write_tab(perm_stats, "Supp4_Permutation_Test.csv")
      }
    } else {
      cat(" Only one virus present; cross-virus concordance not computable.\n")
    }
  }
}

# ==============================================================================
# SUPP_4-B : Meta-analysis across plates (Stouffer's method, )
# ==============================================================================
cat("\n[Supp_4-B] Meta-analysis across plates (Stouffer's method)\n")

if (is.null(de_all) || nrow(de_all) == 0) {
  cat(" WARNING: no per-plate DE results; skipping meta-analysis.\n")
} else {
  # Extract per-stratum z-scores from within-plate Wilcoxon
  # z = sign(FC) * Phi^{-1}(1 - p/2) (one-tailed conversion)
  de_all <- de_all %>%
    mutate(z_stratum = sign(avg_log2FC) * qnorm(pmin(p_value / 2, 0.4999), lower.tail = FALSE),
           w_stratum = 1 / pmax(as.numeric(n_cells), 1))
  
  # Stouffer's method: Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
  meta_stouffer <- de_all %>%
    group_by(ensembl, symbol) %>%
    summarise(
      Z_stouffer = sum(w_stratum * z_stratum, na.rm = TRUE) /
        sqrt(sum(w_stratum^2, na.rm = TRUE)),
      p_meta = 2 * pnorm(-abs(Z_stouffer)),
      n_strata = n(),
      mean_log2FC = mean(avg_log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_meta, method = "BH")) %>%
    arrange(p_adj)
  
  write_tab(meta_stouffer, "Supp4b_Stouffer_MetaAnalysis.csv")
  
  n_sig_meta <- sum(meta_stouffer$p_adj < 0.05 & abs(meta_stouffer$mean_log2FC) >= 0.25, na.rm = TRUE)
  cat(sprintf(" Meta-analysis: %d genes across %d strata, %d significant (FDR<0.05, |log2FC|>=0.25)\n",
              nrow(meta_stouffer), length(unique(de_all$stratum)), n_sig_meta))
  
  # Volcano plot for meta-analysis
  p_meta_volc <- ggplot(meta_stouffer, aes(x = mean_log2FC, y = -log10(p_adj + 1e-300))) +
    geom_point(alpha = 0.4, size = 0.8, colour = "#0072B2") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black", linewidth = 0.5) +
    geom_text_repel(data = meta_stouffer %>% filter(p_adj < 0.001 & abs(mean_log2FC) >= 0.5) %>% slice_head(n = 15),
                    aes(label = symbol), size = 2.5, max.overlaps = 15) + 
    labs(title = "Stouffer meta-analysis (across plates)",
         subtitle = sprintf("%d genes, %d strata, %d significant",
                            nrow(meta_stouffer), length(unique(de_all$stratum)), n_sig_meta),
         x = "Mean log2 fold change", y = "-log10(FDR)")
  supp_and_save(p_meta_volc, "Supp4b_Stouffer_Volcano", width = 8, height = 6)
  
  # Top up- and down-regulated genes
  top_up <- meta_stouffer %>% filter(p_adj < 0.05, mean_log2FC >= 0.25) %>% arrange(desc(mean_log2FC)) %>% slice_head(n = 20)
  top_dn <- meta_stouffer %>% filter(p_adj < 0.05, mean_log2FC <= -0.25) %>% arrange(mean_log2FC) %>% slice_head(n = 20)
  write_tab(bind_rows(
    top_up %>% mutate(direction = "up"),
    top_dn %>% mutate(direction = "down")
  ), "Supp4b_Stouffer_TopGenes.csv")
  
  cat(sprintf(" Top upregulated: %s\n",
              paste(top_up$symbol[1:5], collapse = ", ")))
  cat(sprintf(" Top downregulated: %s\n",
              paste(top_dn$symbol[1:5], collapse = ", ")))
}

# ==============================================================================
# SUPP_5 : Cluster pairwise markers ()
# ==============================================================================
cat("\n[Supp_5] Cluster pairwise markers\n")

if (is.null(clustering) || !requireNamespace("Matrix", quietly = TRUE)) {
  cat(" WARNING: clustering not available; skipping Supp_5.\n")
} else {
  clusters <- clustering$clusters
  n_clust <- nlevels(clusters)
  
  # Background-corrected Wilcoxon for cluster markers
  find_markers <- function(cl1, cl2, label1, label2) {
    idx1 <- which(clusters == cl1)
    idx2 <- which(clusters == cl2)
    if (length(idx1) < 5 || length(idx2) < 5) return(NULL)
    
    sub <- mat_log[, c(idx1, idx2), drop = FALSE]
    keep_g <- rowSums(sub > 0) >= 3
    sub <- sub[keep_g, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    
    is_grp1 <- c(rep(TRUE, length(idx1)), rep(FALSE, length(idx2)))
    R <- matrixStats::rowRanks(sub, ties.method = "average")
    n1 <- length(idx1); n2 <- length(idx2); n <- n1 + n2
    U1 <- rowSums(R[, is_grp1, drop = FALSE]) - n1 * (n1 + 1) / 2
    mu <- n1 * n2 / 2
    varU <- (n1 * n2 / (n * (n - 1))) *
      (rowSums(R^2) - n * ((n + 1) / 2)^2)
    z <- (U1 - mu - sign(U1 - mu) * 0.5) / sqrt(varU)
    p <- 2 * pnorm(-abs(z))
    fc <- rowMeans(sub[, is_grp1, drop = FALSE]) - rowMeans(sub[, !is_grp1, drop = FALSE])
    pct1 <- rowMeans(sub[, is_grp1, drop = FALSE] > 0)
    pct2 <- rowMeans(sub[, !is_grp1, drop = FALSE] > 0)
    
    res <- data.frame(
      cluster1 = label1, cluster2 = label2,
      ensembl = rownames(sub), symbol = sym_of(rownames(sub)),
      avg_log2FC = fc, pct_in_cl1 = pct1, pct_in_cl2 = pct2,
      p_value = p, p_adj = p.adjust(p),
      stringsAsFactors = FALSE
    )
    res <- res[is.finite(res$p_value), ]
    res <- res[order(res$p_adj), ]
    res$comparison <- sprintf("%s vs %s", label1, label2)
    res
  }
  
  cl_order <- levels(clusters)
  marker_list <- list()
  for (i in seq_along(cl_order)) {
    for (j in seq(i + 1, length(cl_order))) {
      mr <- find_markers(cl_order[i], cl_order[j], cl_order[i], cl_order[j])
      if (!is.null(mr) && nrow(mr) > 0) {
        marker_list[[length(marker_list) + 1]] <- mr
      }
    }
  }
  
  if (length(marker_list)) {
    all_markers <- bind_rows(marker_list)
    all_markers$direction <- ifelse(all_markers$avg_log2FC > 0,
                                    sprintf("Up in %s", all_markers$cluster1),
                                    sprintf("Up in %s", all_markers$cluster2))
    all_markers$significant <- all_markers$p_adj < 0.05 & abs(all_markers$avg_log2FC) >= 0.25
    
    write_tab(all_markers, "Supp5_ClusterPairwise_Markers.csv")
    
    sig_summary <- all_markers %>%
      group_by(comparison, direction) %>%
      summarise(n_sig = sum(significant), .groups = "drop")
    write_tab(sig_summary, "Supp5_Marker_Summary.csv")
    cat(sprintf(" Pairwise comparisons: %d\n", length(marker_list)))
    cat(sprintf(" Total significant markers (FDR<0.05, |log2FC|>=0.25): %d\n",
                sum(all_markers$significant)))
    print(as.data.frame(sig_summary), row.names = FALSE)
  } else {
    cat(" No pairwise comparisons produced results.\n")
  }
}

# ==============================================================================
# Update QC with all augmentations
# ==============================================================================
saveRDS(qc, file.path(results_root, "objects", "Phase_D_qc_augmented.rds"))

cat("\n======================================================================\n")
cat("PHASE D SUPPLEMENTARY COMPLETE\n")
cat(sprintf(" Augmented QC : %d cells x %d columns\n", nrow(qc), ncol(qc)))
cat(sprintf(" Key additions: log_viral, pct_ribo, cluster, cluster_name\n"))
cat(sprintf(" Figures: %s\n", dir_supp))
cat(sprintf(" Tables: %s\n", dir_tab))
cat("======================================================================\n")