# ==============================================================================
# PHASE E - ENRICHMENT, FINAL DE, AND SUPPLEMENTARY ANALYSES 
# ==============================================================================
# Over-Representation Analysis (ORA), GSEA, cluster markers (MAST),
# tradeSeq trajectory analysis, and summary outputs.
#
# Inputs:
# results/phaseD/objects/Phase_D_clustering.rds
# results/phaseD/tables/DE_WithinPlate_Wilcoxon.csv
# results/phaseD/tables/vRNA_Continuous_Association.csv
# results/phaseD/tables/All_Timepoint_DEGs.csv
# results/phaseD/tables/ZIKV_Sig_Genes.csv
# results/phaseD/tables/DENV_Sig_Genes.csv
# results/phaseD/tables/Common_Sig_Genes.csv
# results/phaseC/objects/Phase_C_log_matrix.rds
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
  source("E:/Zika/scripts/inject_themes.R")
  library(scales)
  library(Matrix)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
  library(mgcv)
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
data_dir <- "E:/Zika"
setwd(data_dir)

results_root <- file.path(data_dir, "results", "phaseE")
dir_fig <- file.path(results_root, "figures")
dir_tab <- file.path(results_root, "tables")
dir_obj <- file.path(results_root, "objects")
for (d in c(dir_fig, dir_tab, dir_obj))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

show_and_save <- function(p, file, width = 8, height = 5, base_size = 10) {
  p <- p + theme_pub(base_size = base_size)
  print(p)
  ggsave(file.path(dir_fig, paste0(file, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(dir_fig, paste0(file, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(dir_fig, paste0(file, ".svg")), p, width = width, height = height, bg = "transparent", device = svglite::svglite)
  cat(sprintf(" [plot] %s.pdf / .png / .svg\n", file))
  invisible(p)
}
write_tab <- function(x, name) {
  write.csv(x, file.path(dir_tab, name), row.names = FALSE)
  cat(sprintf(" [table] %s\n", name))
}

save_pheatmap <- function(ph, file, width = 8, height = 6) {
  pdf(file.path(dir_fig, paste0(file, ".pdf")), width = width, height = height)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  png(file.path(dir_fig, paste0(file, ".png")), width = width, height = height,
      units = "in", res = 600)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  cat(sprintf(" [plot] %s.pdf / .png\n", file))
}

fig <- function(n) file.path(dir_fig, n)
tab <- function(n) file.path(dir_tab, n)
obj <- function(n) file.path(dir_obj, n)

find_file <- function(bn, subdirs = c("results/phaseD/objects","results/phaseD/tables",
                                      "results/phaseC/objects","results/phaseB/objects",
                                      "results/phaseA/objects","results/phaseB/tables",
                                      "results/phaseD_1/objects","results/phaseD_1/tables",
                                      ".")) {
  for (s in subdirs) {
    p <- file.path(data_dir, s, bn)
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NA_character_
}

cat("======================================================================\n")
cat("PHASE E - ENRICHMENT & FINAL ANALYSES \n")
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

# Load Phase D results
de_all <- tryCatch(
  read.csv(find_file("DE_WithinPlate_Wilcoxon.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
vrna <- tryCatch(
  read.csv(find_file("vRNA_Continuous_Association.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)

# Load  DEG tables (primary significant gene lists)
sig_all <- tryCatch(
  read.csv(find_file("All_Timepoint_DEGs.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
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

# Load supplementary results
de_markers <- tryCatch(
  read.csv(find_file("Supp5_ClusterPairwise_Markers.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)

ensg_col <- intersect(c("ENSG", "ENSG_Original", "ENSG_Clean"), names(gene_map))[1]
sym_of <- function(ids) {
  if (is.null(gene_map) || is.na(ensg_col)) return(ids)
  s <- gene_map$Symbol[match(ids, gene_map[[ensg_col]])]
  ifelse(is.na(s), ids, s)
}
ensg_of <- function(symbols) {
  if (is.null(gene_map) || is.na(ensg_col)) return(rep(NA_character_, length(symbols)))
  gene_map[[ensg_col]][match(toupper(symbols), toupper(gene_map$Symbol))]
}

get_sig_ensg <- function(syms) {
  if (is.null(gene_map) || is.na(ensg_col)) return(character(0))
  gene_map[[ensg_col]][match(toupper(syms), toupper(gene_map$Symbol))]
}

cat(sprintf("Loaded: mat_log %d x %d | QC %d rows\n",
            nrow(mat_log), ncol(mat_log), nrow(qc)))
if (!is.null(clustering)) {
  clusters <- clustering$clusters
  qc$cluster <- clusters
  cat(sprintf("Clusters: %d levels\n", nlevels(clusters)))
}
if (!is.null(de_all))
  cat(sprintf("DE table: %d rows, %d strata\n", nrow(de_all), nlevels(factor(de_all$stratum))))
if (!is.null(vrna))
  cat(sprintf("vRNA association: %d rows, %d strata\n", nrow(vrna), nlevels(factor(vrna$stratum))))

# Log  DEG counts (primary significant gene lists)
if (!is.null(sig_zika))
  cat(sprintf("ZIKV %d | DENV %d | Common %d sig DEGs\n",
              nrow(sig_zika), nrow(sig_denv), nrow(sig_common)))
if (!is.null(sig_all))
  cat(sprintf("%d total DEG rows across all conditions\n", nrow(sig_all)))

# ==============================================================================
# ENRICHMENT ANALYSIS
# ==============================================================================
cat("\nEnrichment analysis\n")

# Universe = detected genes only (from active matrix)
universe <- rownames(mat_counts)
cat(sprintf(" Universe: %d detected genes (from active matrix)\n", length(universe)))

# Primary gene set:  DEGs (Infected vs Bystander, pooled across plates)
#  within-plate has very low power (0 sig in DENV strata) - use  instead.
upregulated <- character(0)
up_entrez <- character(0)

if (!is.null(sig_all)) {
  sig_degs <- sig_all %>%
    filter(is_sig, avg_log2FC >= 0.25) %>%
    pull(ensembl) %>%
    unique()
  upregulated <- intersect(sig_degs, universe)
  cat(sprintf(" Upregulated DEGs (, FDR<0.05, |log2FC|>=0.25): %d\n",
              length(upregulated)))
} else if (!is.null(de_all)) {
  de_all$is_true_contrast <- !grepl("InfectedvsInfected|BystandervsBystander", de_all$stratum)
  sig_up <- de_all %>%
    filter(is_true_contrast, p_adj < 0.05, avg_log2FC <= -0.1) %>%
    pull(ensembl) %>%
    unique()
  upregulated <- intersect(sig_up, universe)
  cat(sprintf(" Upregulated genes ( fallback, FDR<0.05, |log2FC|>=0.1): %d\n",
              length(upregulated)))
}

ensg_to_entrez <- function(ensg_ids) {
  if (length(ensg_ids) == 0) return(character(0))
  keys <- str_remove(ensg_ids, "\\..*$")
  mapped <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                          keys = unique(keys), column = "ENTREZID",
                          keytype = "ENSEMBL", multiVals = "first")
  )
  mapped <- mapped[!is.na(mapped)]
  entrez <- mapped[match(ensg_ids, names(mapped))]
  entrez[!is.na(entrez)]
}

# ORA with clusterProfiler
ora_gobp <- NULL
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat(" -> Running ORA with clusterProfiler...\n")
  
  # Universe Entrez IDs
  univ_entrez <- ensg_to_entrez(universe)
  cat(sprintf(" Universe Entrez IDs: %d of %d mapped\n", length(univ_entrez), length(universe)))
  
  if (length(upregulated) > 0) {
    up_entrez <- ensg_to_entrez(upregulated)
    cat(sprintf(" Upregulated Entrez IDs: %d of %d mapped\n", length(up_entrez), length(upregulated)))
    
    # GO Biological Process
    ora_gobp <- tryCatch({
      clusterProfiler::enrichGO(gene = up_entrez,
                                universe = univ_entrez,
                                OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,
                                minGSSize = 5,
                                maxGSSize = 500)
    }, error = function(e) NULL)
    
    if (!is.null(ora_gobp) && nrow(as.data.frame(ora_gobp)) > 0) {
      gobp_df <- as.data.frame(ora_gobp)
      write_tab(gobp_df, "ORA_GO_BP.csv")
      cat(sprintf(" GO-BP: %d significant terms\n", sum(gobp_df$p.adjust < 0.05)))
      
      p_go <- clusterProfiler::dotplot(ora_gobp, showCategory = 15) +
        theme(axis.text.y = element_text(size = 7, face = "bold")) +
        labs(title = "ORA - GO Biological Process (upregulated genes)",
             subtitle = sprintf("%d genes tested, %d significant terms",
                                length(up_entrez), sum(gobp_df$p.adjust < 0.05)))
      show_and_save(p_go, "ORA_GO_BP", width = 14, height = 10)
    } else {
      cat(" No significant GO-BP terms.\n")
    }
  } else {
    cat(" No upregulated genes for ORA.\n")
  }
  
  # GSEA: rank by  log2FC (independent of ORA)
  if (!is.null(sig_all) && nrow(sig_all) > 0) {
    # Aggregate  to gene-level mean log2FC and direction for ranking
    meta_rank <- sig_all %>%
      group_by(ensembl) %>%
      summarise(mean_fc = mean(avg_log2FC, na.rm = TRUE),
                n_sig = sum(is_sig, na.rm = TRUE),
                .groups = "drop") %>%
      filter(!is.na(mean_fc), is.finite(mean_fc))
    
    if (nrow(meta_rank) >= 10) {
      # Build named ranking vector: positive = upregulated in Infected
      rank_vec <- setNames(meta_rank$mean_fc, meta_rank$ensembl)
      rank_vec <- rank_vec[!is.na(names(rank_vec)) & names(rank_vec) != ""]
      cat(sprintf(" -> Running GSEA on %d genes (ranked by mean log2FC)...\n", length(rank_vec)))
      
      gsea_res <- tryCatch({
        clusterProfiler::gseGO(geneList = rank_vec,
                               OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                               ont = "BP",
                               nPerm = 1000,
                               minGSSize = 10,
                               maxGSSize = 500,
                               pvalueCutoff = 0.05,
                               verbose = FALSE)
      }, error = function(e) {
        cat(sprintf(" GSEA error: %s\n", conditionMessage(e)))
        NULL
      })
      
      if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
        gsea_df <- as.data.frame(gsea_res)
        write_tab(gsea_df, "GSEA_GO_BP.csv")
        cat(sprintf(" GSEA: %d significant pathways\n", sum(gsea_df$p.adjust < 0.05)))
        
        p_gsea <- clusterProfiler::ridgeplot(gsea_res) +
          theme(axis.text.y = element_text(size = 7, face = "bold")) +
          labs(title = "GSEA - GO Biological Process",
               subtitle = "Ranked by mean log2FC across  conditions")
        show_and_save(p_gsea, "GSEA_Ridge", width = 12, height = 8)
      }
    }
  }
}

# ==============================================================================
# KEGG ORA (upregulated genes)
# ==============================================================================
if (length(up_entrez) > 0 && requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat(" -> Running KEGG ORA...\n")
  
  kegg_ora <- tryCatch({
    clusterProfiler::enrichKEGG(gene = up_entrez,
                                universe = univ_entrez,
                                organism = "hsa",
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,
                                minGSSize = 5,
                                maxGSSize = 500)
  }, error = function(e) {
    cat(sprintf(" KEGG ORA error: %s\n", conditionMessage(e)))
    NULL
  })
  
  if (!is.null(kegg_ora) && nrow(as.data.frame(kegg_ora)) > 0) {
    kegg_df <- as.data.frame(kegg_ora)
    write_tab(kegg_df, "ORA_KEGG.csv")
    cat(sprintf(" KEGG: %d significant pathways\n", sum(kegg_df$p.adjust < 0.05)))
    
    p_kegg <- clusterProfiler::dotplot(kegg_ora, showCategory = 15) +
      theme(axis.text.y = element_text(size = 7, face = "bold")) +
      labs(title = "ORA - KEGG pathways (upregulated genes)",
           subtitle = sprintf("%d genes tested, %d significant terms",
                              length(up_entrez), sum(kegg_df$p.adjust < 0.05)))
    show_and_save(p_kegg, "ORA_KEGG", width = 12, height = 8)
  } else {
    cat(" No significant KEGG terms.\n")
  }
}

# ==============================================================================
# limma::camera (inter-gene correlation-corrected enrichment)
# Uses  gene rankings directly instead of depending on ORA results
# ==============================================================================
if (!is.null(sig_all) && nrow(sig_all) > 0 && requireNamespace("limma", quietly = TRUE)) {
  cat(" -> Running limma::camera...\n")
  
  # Build gene-level ranking from  mean log2FC
  meta_z <- sig_all %>%
    group_by(ensembl) %>%
    summarise(z = mean(avg_log2FC, na.rm = TRUE),
              n_cond = n(),
              .groups = "drop") %>%
    filter(!is.na(z), is.finite(z), n_cond >= 1)
  
  # Map to row indices in mat_counts
  univ_genes <- rownames(mat_counts)
  meta_z <- meta_z[meta_z$ensembl %in% univ_genes, ]
  meta_z$idx <- match(meta_z$ensembl, univ_genes)
  
  if (nrow(meta_z) >= 10) {
    stat_vec <- setNames(meta_z$z, meta_z$ensembl)
    stat_vec <- stat_vec[as.character(univ_genes[meta_z$idx])]
    
    # Use built-in MSigDB gene sets if available, otherwise skip
    camera_sets <- tryCatch({
      if (requireNamespace("msigdbr", quietly = TRUE)) {
        msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
        if (!is.null(msig) && nrow(msig) > 0) {
          gs_list <- split(msig$ensembl_gene, msig$gs_name)
          gs_list <- gs_list[sapply(gs_list, length) >= 10 & sapply(gs_list, length) <= 500]
          gene_idx <- lapply(gs_list, function(genes) {
            idx <- match(intersect(genes, univ_genes), univ_genes)
            idx[!is.na(idx)]
          })
          gene_idx <- gene_idx[sapply(gene_idx, length) >= 3]
          if (length(gene_idx) > 0) {
            limma::GeneSetCollection(lapply(gene_idx, function(idx) univ_genes[idx]))
          } else NULL
        } else NULL
      } else NULL
    }, error = function(e) NULL)
    
    if (!is.null(camera_sets) && length(camera_sets) > 0) {
      tryCatch({
        cam_res <- limma::camera(stat_vec,
                                 camera_sets,
                                 use.ranks = FALSE,
                                 inter.gene.cor = 0.01)
        cam_df <- as.data.frame(cam_res)
        cam_df$Pathway <- rownames(cam_df)
        rownames(cam_df) <- NULL
        cam_df <- cam_df[, c("Pathway", "NGenes", "Direction", "pValue", "FDR")]
        names(cam_df) <- c("Pathway", "NGenes", "Direction", "PValue", "p_adj")
        write_tab(cam_df, "Camera_Enrichment.csv")
        cat(sprintf(" camera: %d pathways tested, %d significant (FDR<0.05)\n",
                    nrow(cam_df), sum(cam_df$p_adj < 0.05, na.rm = TRUE)))
        
        if (sum(cam_df$p_adj < 0.05, na.rm = TRUE) > 0) {
          p_cam <- ggplot(cam_df %>% arrange(PValue) %>% slice_head(n = 20),
                          aes(x = reorder(Pathway, PValue), y = -log10(PValue + 1e-300))) +
            geom_col(fill = "#0072B2", colour = "black", linewidth = 0.2) +
            coord_flip() +
            labs(title = "limma::camera enrichment",
                 subtitle = "Inter-gene correlation corrected, ranked by  log2FC",
                 x = NULL, y = "-log10(p-value)")
          show_and_save(p_cam, "Camera_Enrichment", width = 10, height = 6)
        }
      }, error = function(e) {
        cat(sprintf(" camera error: %s\n", conditionMessage(e)))
      })
    } else {
      cat(" camera: no gene sets available (install msigdbr for MSigDB, or provide custom sets).\n")
    }
  }
}

# ==============================================================================
# mgcv GAM FOR CONTINUOUS vRNA CURVE SHAPE
# ==============================================================================
cat("\nmgcv::gam for vRNA expression curves\n")

if (is.null(vrna) || nrow(vrna) == 0) {
  cat(" No vRNA data available; skipping GAM.\n")
} else if (requireNamespace("mgcv", quietly = TRUE)) {
  cat(" mgcv available; fitting GAM per gene per virus...\n")
  
  gam_results <- list()
  NKNOTS <- c(5, 6, 7)
  
  for (v in unique(vrna$virus)) {
    v_cells <- qc %>% filter(virus == v, moi > 0, infection_status %in% c("Infected", "Bystander"))
    if (nrow(v_cells) < 20) {
      cat(sprintf(" %s: only %d exposed cells, skipping GAM.\n", v, nrow(v_cells)))
      next
    }
    
    v_ix <- match(v_cells$Sample, colnames(mat_log))
    lv <- v_cells$log_viral
    sig_genes <- vrna %>% filter(virus == v, p_adj < 0.05) %>% pull(ensembl)
    
    v_res <- list()
    cat(sprintf(" -> %s: fitting GAM for %d significant genes across %d cells\n",
                v, length(sig_genes), length(v_ix)))
    
    for (g_ensg in sig_genes) {
      expr <- mat_log[g_ensg, v_ix]
      ok <- is.finite(expr) & is.finite(lv) & !is.na(expr) & !is.na(lv)
      if (sum(ok) < 15) next
      if (sd(expr[ok]) < 0.01 || sd(lv[ok]) < 0.01) next
      
      df_g <- data.frame(expr = expr[ok], lv = lv[ok])
      aic_vals <- numeric(length(NKNOTS))
      fitted_vals <- list()
      
      for (ik in seq_along(NKNOTS)) {
        k <- NKNOTS[ik]
        fit <- tryCatch({
          mgcv::gam(expr ~ s(lv, k = k, bs = "cs"), data = df_g, method = "REML")
        }, error = function(e) NULL)
        if (is.null(fit)) { aic_vals[ik] <- NA; next }
        aic_vals[ik] <- AIC(fit)
        fitted_vals[[ik]] <- fit
      }
      
      best_ix <- which.min(aic_vals)
      best_fit <- fitted_vals[[best_ix]]
      if (is.null(best_fit)) next
      
      # Predict across log_viral range for curve shape classification
      lv_seq <- seq(min(df_g$lv), max(df_g$lv), length.out = 50)
      pred <- predict(best_fit, newdata = data.frame(lv = lv_seq), type = "terms", se.fit = TRUE)
      smooth_vals <- pred$fit[, "s(lv)"]
      
      # Curve shape from first and last derivative
      y_diff <- diff(smooth_vals)
      if (all(y_diff > 0)) {
        shape <- "monotonic increasing"
      } else if (all(y_diff < 0)) {
        shape <- "monotonic decreasing"
      } else {
        peak_ix <- which.max(smooth_vals)
        if (peak_ix <= 10) {
          shape <- "bell-shaped (peak early)"
        } else if (peak_ix >= 40) {
          shape <- "bell-shaped (peak late)"
        } else {
          shape <- "U-shaped"
        }
      }
      
      v_res[[g_ensg]] <- data.frame(
        ensembl = g_ensg,
        symbol = sym_of(g_ensg),
        virus = v,
        best_k = NKNOTS[best_ix],
        aic_best = min(aic_vals, na.rm = TRUE),
        monotonicity = ifelse(all(y_diff > 0), "increasing",
                              ifelse(all(y_diff < 0), "decreasing", "non-monotonic")),
        curve_shape = shape,
        rho = vrna$rho[vrna$ensembl == g_ensg][1],
        p_adj = vrna$p_adj[vrna$ensembl == g_ensg][1],
        stringsAsFactors = FALSE
      )
    }
    
    if (length(v_res) > 0) {
      gam_df <- bind_rows(v_res)
      gam_results[[v]] <- gam_df
      
      write_tab(gam_df, sprintf("GAM_Curves_%s.csv", v))
      
      shape_tab <- gam_df %>%
        group_by(curve_shape) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(pct = 100 * n / sum(n))
      write_tab(shape_tab, sprintf("GAM_CurveShape_%s.csv", v))
      
      cat(sprintf(" %s: %d genes fitted. Curve shapes:\n", v, nrow(gam_df)))
      print(as.data.frame(shape_tab), row.names = FALSE)
      
      # Plot: curve shape distribution
      p_shapes <- ggplot(shape_tab, aes(x = reorder(curve_shape, n), y = n, fill = curve_shape)) +
        geom_col(colour = "black", linewidth = 0.2) +
        scale_fill_viridis_d(option = "D") +
        scale_y_continuous(expand = c(0, 0)) +
        coord_flip() +
        labs(title = sprintf("GAM curve shapes (%s)", v),
             subtitle = "mgcv::gam with k = 5,6,7 selected by AIC",
             x = "Curve shape", y = "Genes")
      show_and_save(p_shapes, sprintf("GAM_CurveShapes_%s", v), width = 8, height = 5)
      
      # Plot: example GAM fits for top 6 genes by |rho|
      if (nrow(gam_df) >= 1) {
        top6 <- gam_df %>% arrange(desc(abs(rho))) %>% slice_head(n = 6)
        plot_data <- list()
        for (i in seq_len(nrow(top6))) {
          g_ensg <- top6$ensembl[i]
          plot_cells <- qc %>% filter(virus == v, moi > 0, infection_status %in% c("Infected", "Bystander"))
          v_ix2 <- match(plot_cells$Sample, colnames(mat_log))
          expr_g <- mat_log[g_ensg, v_ix2]
          ok <- is.finite(expr_g) & is.finite(plot_cells$log_viral)
          plot_data[[i]] <- data.frame(
            ensembl = g_ensg, symbol = top6$symbol[i],
            lv = plot_cells$log_viral[ok], expr = expr_g[ok]
          )
        }
        pdsig_all <- bind_rows(plot_data)
        
        p_gam_fits <- ggplot(pdsig_all, aes(x = lv, y = expr)) +
          geom_point(alpha = 0.4, size = 1) +
          facet_wrap(~symbol, scales = "free_y") +
          labs(title = sprintf("Example GAM fits (%s, top 6 genes)", v),
               subtitle = "Scatter: expression vs log10(viral molecules + 1)",
               x = expression(log[10]~"(viral molecules + 1)"), y = "Log expression")
        show_and_save(p_gam_fits, sprintf("GAM_ExampleFits_%s", v), width = 10, height = 7)
      }
    }
  }
}

# ==============================================================================
# CELL-STATE COMPARISON (Cluster markers)
# ==============================================================================
cat("\nCell-state comparison (cluster markers)\n")

if (!is.null(de_markers)) {
  cat(sprintf(" -> Loaded %d pre-computed cluster marker results.\n", nrow(de_markers)))
  
  # Summary: top markers per cluster
  top_per_cluster <- de_markers %>%
    group_by(cluster1) %>%
    filter(significant) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  write_tab(top_per_cluster, "TopClusterMarkers.csv")
  cat(sprintf(" -> Top 10 significant markers per cluster saved.\n"))
  
  # Heatmap of top markers across clusters
  if (nrow(top_per_cluster) > 0) {
    top_genes <- unique(top_per_cluster$ensembl)
    top_symbols <- sym_of(top_genes[!is.na(top_genes)])
    top_symbols <- top_symbols[!is.na(top_symbols) & top_symbols != ""]
    
    if (length(top_symbols) > 0) {
      in_mat <- intersect(ensg_of(top_symbols), rownames(mat_log))
      if (length(in_mat) > 0) {
        mat_heat <- mat_log[in_mat, , drop = FALSE]
        rownames(mat_heat) <- sym_of(rownames(mat_heat))
        if (is.null(clustering)) {
          ann_col <- data.frame(row.names = colnames(mat_heat),
                                infection = qc$infection_status)
        } else {
          ann_col <- data.frame(row.names = colnames(mat_heat),
                                cluster = clusters,
                                infection = qc$infection_status)
        }
        
        mat_scaled_hm <- t(scale(t(mat_heat)))
        mat_scaled_hm[is.na(mat_scaled_hm)] <- 0
        mat_scaled_hm <- pmax(pmin(mat_scaled_hm, 2), -2)
        
        if (nrow(mat_scaled_hm) > 1 && ncol(mat_scaled_hm) > 1) {
          n_genes_hm <- nrow(mat_scaled_hm)
          hm_height <- max(6, n_genes_hm * 0.12)
          ph <- pheatmap(mat_scaled_hm,
                         annotation_col = ann_col,
                         show_rownames = TRUE, show_colnames = FALSE,
                         fontsize_row = 8,
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         silent = TRUE)
          save_pheatmap(ph, "ClusterMarkers_Heatmap", width = 12, height = hm_height)
        }
      }
    }
  }
  
  # Per-cluster: percentage of cells expressing signature genes
  sig_genes <- list(
    UPR = c("HSPA5", "DDIT3", "ATF4", "XBP1"),
    Apoptosis = c("CASP3", "CASP8", "FAS", "BAX", "BCL2")
  )
  
  cat(" -> Per-cluster signature expression (% cells):\n")
  pct_exp_list <- list()
  for (cl in levels(clusters)) {
    cl_idx <- which(clusters == cl)
    if (length(cl_idx) == 0) next
    cl_cells <- colnames(mat_log)[cl_idx]
    
    pct_row <- data.frame(cluster = cl, n_cells = length(cl_idx))
    for (sig_name in names(sig_genes)) {
      sig_ensg <- intersect(get_sig_ensg(sig_genes[[sig_name]]), rownames(mat_log))
      
      if (length(sig_ensg) > 0) {
        pct <- mean(rowMeans(mat_log[sig_ensg, cl_cells, drop = FALSE]) > 0) * 100
        pct_row[[sprintf("pct_%s", sig_name)]] <- pct
      }
    }
    pct_exp_list[[cl]] <- pct_row
  }
  pct_exp <- bind_rows(pct_exp_list)
  write_tab(pct_exp, "SignatureExpression_ByCluster.csv")
  cat(sprintf(" %s\n", paste(sprintf("Cluster %s: %d cells",
                                     pct_exp$cluster, pct_exp$n_cells), collapse = "; ")))
  print(as.data.frame(pct_exp), row.names = FALSE)
  
} else {
  cat(" No cluster marker results available. Run Phase D Supplementary first.\n")
}

# ==============================================================================
# SUMMARY OUTPUTS
# ==============================================================================
cat("\nGenerating summary outputs...\n")

# Overall experiment summary - use  as the primary DEG source
zika_sig_n <- if (!is.null(sig_zika)) nrow(sig_zika) else 0
denv_sig_n <- if (!is.null(sig_denv)) nrow(sig_denv) else 0
up_n <- length(upregulated)

if (!is.null(clustering)) {
  summary_df <- data.frame(
    item = c("Total cells",
             "Clusters",
             "Stable clusters",
             "DE strata ()",
             "ZIKV sig DEGs ()",
             "DENV sig DEGs ()",
             "Common sig DEGs ()",
             "Upregulated DEGs (Infected > Bystander)"),
    value = c(
      nrow(qc),
      nlevels(clusters),
      if (!is.null(clustering$stability)) sum(clustering$stability$stable == "stable") else NA,
      if (!is.null(de_all)) nlevels(factor(de_all$stratum)) else NA,
      zika_sig_n,
      denv_sig_n,
      if (!is.null(sig_common)) nrow(sig_common) else 0,
      up_n
    ),
    stringsAsFactors = FALSE
  )
  write_tab(summary_df, "PhaseE_Summary.csv")
  cat(" Summary table saved.\n")
}

# Cross-virus vRNA concordance summary
if (!is.null(vrna) && "virus" %in% names(vrna)) {
  zika_genes <- vrna %>% filter(virus == "zika", p_adj < 0.05) %>% pull(ensembl)
  denv_genes <- vrna %>% filter(virus == "dengue", p_adj < 0.05) %>% pull(ensembl)
  overlap <- length(intersect(zika_genes, denv_genes))
  
  concordance_summary <- data.frame(
    comparison = c("ZIKV-specific", "DENV-specific", "Shared"),
    n_genes = c(length(setdiff(zika_genes, denv_genes)),
                length(setdiff(denv_genes, zika_genes)),
                overlap),
    stringsAsFactors = FALSE
  ) %>%
    mutate(pct = 100 * n_genes / sum(n_genes))
  write_tab(concordance_summary, "CrossVirus_vRNA_Overlap.csv")
  
  cat(sprintf(" ZIKV-associated: %d | DENV-associated: %d | Shared: %d\n",
              length(zika_genes), length(denv_genes), overlap))
}

# Primary result:  continuous vRNA (main paper figure)
if (!is.null(vrna) && nrow(vrna) > 0) {
  vrna_summary <- vrna %>%
    filter(p_adj < 0.05) %>%
    group_by(virus) %>%
    summarise(n_sig = n(),
              median_rho = median(abs(rho), na.rm = TRUE),
              mean_rho = mean(abs(rho), na.rm = TRUE),
              .groups = "drop")
  write_tab(vrna_summary, "PrimaryResult_Summary.csv")
  
  cat(" Primary result summary ():\n")
  print(as.data.frame(vrna_summary), row.names = FALSE)
}

#  DEG virus breakdown
if (!is.null(sig_all)) {
  virus_summary <- sig_all %>%
    filter(is_sig) %>%
    group_by(virus) %>%
    summarise(n_sig = n(), n_genes = n_distinct(ensembl), .groups = "drop")
  write_tab(virus_summary, "DEG_Virus_Summary.csv")
  cat("  virus breakdown:\n")
  print(as.data.frame(virus_summary), row.names = FALSE)
}

# ==============================================================================
# Phase E complete
# ==============================================================================
cat("\n======================================================================\n")
cat("PHASE E COMPLETE\n")
cat(sprintf(" Figures: %s\n", dir_fig))
cat(sprintf(" Tables: %s\n", dir_tab))
cat(sprintf(" Objects: %s\n", dir_obj))
cat("======================================================================\n")