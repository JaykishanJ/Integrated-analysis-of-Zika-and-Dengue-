#!/usr/bin/env Rscript
# CROSS-VIRUS COMPARISON -- Zika vs Dengue, tissue held constant.
#
#   Rscript compare/compare_zika_dengue.R
#
# Takes the per-study DESeq2 results already produced (steps 1-9 of each series) and asks how
# far the two viruses agree, using four independent statistical readouts on the genes that move
# the SAME WAY in both:
#
#   1. Hypergeometric test          -- is the DEG overlap bigger than chance?
#   2. Sign test + RRHO             -- do shared genes agree in direction? (thresholded and
#                                     threshold-free)
#   3. Effect-size correlation      -- do the log2 fold changes track each other?
#   4. Pathway analysis             -- what biology do the concordant genes represent?
#
# EVERY comparison is restricted to genes BOTH contrasts actually tested. DESeq2's independent
# filtering retains different gene sets per contrast, so an unrestricted intersection would
# treat "never tested here" as "not changed here" and understate agreement.
#
# Comparison A is the scientifically clean one: DENV and ZIKV in the SAME cell type, from the
# SAME study, processed together. Comparison B pairs two different studies and is reported with
# its confounds stated -- Huh7.5 carries a RIG-I T55I defect that Huh7 does not, so any
# interferon-related difference there is uninterpretable.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
  library(pheatmap); library(RColorBrewer)
})

DPI      <- 600
PADJ_CUT <- 0.05
LFC_CUT  <- 1
ROOT     <- "compare"
dir.create(ROOT, showWarnings = FALSE)

# ------------------------------------------------------------------ input locations
# GSE110512 has a single contrast so its file carries no prefix; the others are prefixed.
res_path <- function(gse, tag = NULL) {
  d <- file.path("results", gse, "step7_results")
  if (is.null(tag)) file.path(d, "results_all_genes.csv")
  else file.path(d, paste0(tag, "_results_all_genes.csv"))
}

COMPARISONS <- list(
  list(id = "A_hNPC_GSE80434_DENV_vs_ZIKVM",
       tissue = "hNPC", quality = "PRIMARY -- same tissue, same study, processed together",
       dengue = list(label = "DENV (GSE80434, hNPC)",  file = res_path("GSE80434", "DENV_vs_Mock")),
       zika   = list(label = "ZIKVM (GSE80434, hNPC)", file = res_path("GSE80434", "ZIKVM_vs_Mock")),
       caveat = "None material. Both arms are hNPC from one study; each is tested against its own mock set."),
  list(id = "B_hNPC_GSE80434_DENV_vs_ZIKVC",
       tissue = "hNPC", quality = "PRIMARY -- same tissue, same study; ZIKV Asian strain",
       dengue = list(label = "DENV (GSE80434, hNPC)",  file = res_path("GSE80434", "DENV_vs_Mock")),
       zika   = list(label = "ZIKVC (GSE80434, hNPC)", file = res_path("GSE80434", "ZIKVC_vs_Mock")),
       caveat = "ZIKVC has n=2 against a control group whose internal spread exceeds the effect; it is the weakest arm in the project."),
  list(id = "C_Huh7lineage_crossstudy_DENV_vs_ZIKV",
       tissue = "Huh7 / Huh7.5", quality = "SECONDARY -- cross-study, confounded",
       dengue = list(label = "DENV (GSE110512, Huh7)",     file = res_path("GSE110512")),
       zika   = list(label = "ZIKV (GSE161783, Huh7.5)",   file = res_path("GSE161783", "Huh7_5__Zika_vs_Mock")),
       caveat = "Huh7.5 carries a RIG-I T55I defect Huh7 does not, so interferon differences are confounded with genotype. Strain, timepoint (18h vs 24h), MOI and study also differ.")
)

log_lines <- character(0)
say <- function(...) { m <- paste0(...); cat(m, "\n"); log_lines <<- c(log_lines, m) }

# ═══════════════════════════════════════════════════════════ RRHO implementation
# Rank-rank hypergeometric overlap. RRHO/RRHO2 are not installed, and the algorithm is short
# enough to implement transparently: rank both lists by SIGNED significance, then at every
# point on a grid of rank cut-offs ask whether the two top-lists overlap more than chance.
# Threshold-free, so it does not depend on the padj/LFC cut-offs used elsewhere.
rrho <- function(g, s1, s2, steps = 50) {
  N <- length(g)
  o1 <- order(s1, decreasing = TRUE); o2 <- order(s2, decreasing = TRUE)
  r1 <- integer(N); r1[o1] <- seq_len(N)
  r2 <- integer(N); r2[o2] <- seq_len(N)
  cuts <- unique(round(seq(max(10, floor(N / steps)), N, length.out = steps)))
  M <- matrix(NA_real_, length(cuts), length(cuts),
              dimnames = list(cuts, cuts))
  for (a in seq_along(cuts)) for (b in seq_along(cuts)) {
    k1 <- cuts[a]; k2 <- cuts[b]
    ov <- sum(r1 <= k1 & r2 <= k2)
    # -log10 hypergeometric p for enrichment of the overlap
    p <- phyper(ov - 1L, k1, N - k1, k2, lower.tail = FALSE)
    M[a, b] <- -log10(max(p, .Machine$double.xmin))
  }
  # Signed: positive where the two agree in direction, negative where they oppose.
  list(matrix = M, cuts = cuts, N = N,
       max_up = M[1, 1],                       # both most-up
       max_overall = max(M, na.rm = TRUE),
       rank1 = r1, rank2 = r2)
}

# ═══════════════════════════════════════════════════════════ per-comparison work
summary_rows <- list()

for (CMP in COMPARISONS) {
  out <- file.path(ROOT, CMP$id)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  say("\n", strrep("═", 78))
  say("  ", CMP$id)
  say("  tissue : ", CMP$tissue)
  say("  status : ", CMP$quality)
  say("  Dengue : ", CMP$dengue$label)
  say("  Zika   : ", CMP$zika$label)
  say("  caveat : ", CMP$caveat)
  say(strrep("═", 78))

  for (f in c(CMP$dengue$file, CMP$zika$file))
    if (!file.exists(f)) stop("missing results file: ", f)
  D <- fread(CMP$dengue$file); Z <- fread(CMP$zika$file)

  # ---------------------------------------------------- shared tested universe
  uD <- D$gene_id[!is.na(D$padj)]; uZ <- Z$gene_id[!is.na(Z$padj)]
  U  <- intersect(uD, uZ)
  say("\n── Shared tested universe ──")
  say("  Dengue tested : ", format(length(uD), big.mark = ","))
  say("  Zika tested   : ", format(length(uZ), big.mark = ","))
  say("  shared (used) : ", format(length(U), big.mark = ","),
      sprintf("  (%.1f%% of the smaller set)", 100 * length(U) / min(length(uD), length(uZ))))

  d <- D[gene_id %in% U]; z <- Z[gene_id %in% U]
  M <- merge(d[, .(gene_id, symbol, lfc_D = log2FC_shrunk, padj_D = padj,
                   p_D = pvalue, bm_D = baseMean)],
             z[, .(gene_id, lfc_Z = log2FC_shrunk, padj_Z = padj, p_Z = pvalue,
                   bm_Z = baseMean)], by = "gene_id")
  M[, `:=`(sig_D = !is.na(padj_D) & padj_D < PADJ_CUT & abs(lfc_D) > LFC_CUT,
           sig_Z = !is.na(padj_Z) & padj_Z < PADJ_CUT & abs(lfc_Z) > LFC_CUT)]
  M[, dir_D := ifelse(lfc_D > 0, "up", "down")]
  M[, dir_Z := ifelse(lfc_Z > 0, "up", "down")]
  M[, concordant := sig_D & sig_Z & dir_D == dir_Z]
  M[, discordant := sig_D & sig_Z & dir_D != dir_Z]

  nD <- sum(M$sig_D); nZ <- sum(M$sig_Z)
  say("\n── DEGs within the shared universe (padj < ", PADJ_CUT, ", |log2FC| > ", LFC_CUT, ") ──")
  say("  Dengue DEGs : ", nD, "  (up ", sum(M$sig_D & M$dir_D == "up"),
      ", down ", sum(M$sig_D & M$dir_D == "down"), ")")
  say("  Zika DEGs   : ", nZ, "  (up ", sum(M$sig_Z & M$dir_Z == "up"),
      ", down ", sum(M$sig_Z & M$dir_Z == "down"), ")")
  say("  significant in BOTH : ", sum(M$sig_D & M$sig_Z))
  say("  CONCORDANT (same direction) : ", sum(M$concordant),
      "   DISCORDANT : ", sum(M$discordant))

  concordant <- M[concordant == TRUE]
  fwrite(M, file.path(out, "merged_gene_table.csv"))
  fwrite(concordant[order(-abs(lfc_D + lfc_Z) / 2)], file.path(out, "concordant_genes.csv"))
  fwrite(concordant[dir_D == "up"],   file.path(out, "concordant_up.csv"))
  fwrite(concordant[dir_D == "down"], file.path(out, "concordant_down.csv"))
  fwrite(M[discordant == TRUE],       file.path(out, "discordant_genes.csv"))

  # ══════════════════════════════════ 1. HYPERGEOMETRIC TEST
  say("\n═══ 1. HYPERGEOMETRIC TEST ═══")
  say("  Is the DEG overlap larger than expected if the two viruses were independent?")
  hyper <- function(lbl, a, b) {
    N <- nrow(M); k1 <- sum(a); k2 <- sum(b); ov <- sum(a & b)
    expct <- k1 * k2 / N
    p <- phyper(ov - 1L, k1, N - k1, k2, lower.tail = FALSE)
    data.table(comparison = lbl, universe = N, dengue_n = k1, zika_n = k2,
               observed_overlap = ov, expected_overlap = round(expct, 2),
               fold_enrichment = round(ov / max(expct, 1e-12), 2),
               p_value = p, jaccard = round(ov / length(union(which(a), which(b))), 4))
  }
  HY <- rbindlist(list(
    hyper("all DEGs (any direction)", M$sig_D, M$sig_Z),
    hyper("upregulated in both",   M$sig_D & M$dir_D == "up",   M$sig_Z & M$dir_Z == "up"),
    hyper("downregulated in both", M$sig_D & M$dir_D == "down", M$sig_Z & M$dir_Z == "down")))
  HY[, p_adj := p.adjust(p_value, "BH")]
  print(HY, row.names = FALSE)
  log_lines <<- c(log_lines, capture.output(print(HY, row.names = FALSE)))
  fwrite(HY, file.path(out, "1_hypergeometric_test.csv"))
  say("  -> the one-sided test asks for ENRICHMENT; fold_enrichment is observed/expected.")

  # ══════════════════════════════════ 2. SIGN TEST + RRHO
  say("\n═══ 2. SIGN TEST and RRHO ═══")
  both <- M[sig_D & sig_Z]
  if (nrow(both) >= 5) {
    agree <- sum(both$dir_D == both$dir_Z)
    bt <- binom.test(agree, nrow(both), p = 0.5, alternative = "two.sided")
    say("  Sign test on the ", nrow(both), " genes significant in BOTH:")
    say("    agreeing in direction : ", agree, " / ", nrow(both),
        sprintf("  (%.1f%%)", 100 * agree / nrow(both)))
    say("    binomial p vs 50%%    : ", format.pval(bt$p.value, digits = 3))
    say("    95%% CI on agreement  : ", sprintf("%.3f - %.3f", bt$conf.int[1], bt$conf.int[2]))
    ST <- data.table(n_both = nrow(both), n_agree = agree,
                     pct_agree = round(100 * agree / nrow(both), 2),
                     binom_p = bt$p.value,
                     ci_low = bt$conf.int[1], ci_high = bt$conf.int[2])
  } else {
    say("  Sign test skipped: only ", nrow(both), " genes significant in both.")
    ST <- data.table(n_both = nrow(both), n_agree = NA_integer_, pct_agree = NA_real_,
                     binom_p = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  }
  fwrite(ST, file.path(out, "2a_sign_test.csv"))

  # RRHO -- threshold-free, so it does not inherit the padj/LFC cut-offs
  say("\n  RRHO (threshold-free rank-rank hypergeometric overlap):")
  M[, sig_score_D := -log10(pmax(p_D, .Machine$double.xmin)) * sign(lfc_D)]
  M[, sig_score_Z := -log10(pmax(p_Z, .Machine$double.xmin)) * sign(lfc_Z)]
  RR <- rrho(M$gene_id, M$sig_score_D, M$sig_score_Z, steps = 50)
  say("    genes ranked        : ", RR$N)
  say("    max -log10(p) overall: ", round(RR$max_overall, 1))
  say("    concordant-up corner : ", round(RR$max_up, 1), " (-log10 p)")
  # Also the down-down corner: reverse both rankings.
  RRd <- rrho(M$gene_id, -M$sig_score_D, -M$sig_score_Z, steps = 50)
  say("    concordant-down corner: ", round(RRd$max_up, 1), " (-log10 p)")
  # And the discordant corners, which should be weak if the viruses agree.
  RRx <- rrho(M$gene_id, M$sig_score_D, -M$sig_score_Z, steps = 50)
  say("    DISCORDANT corner    : ", round(RRx$max_up, 1),
      " (-log10 p)  <- should be far lower if the viruses agree")
  RRtab <- data.table(corner = c("up-up (concordant)", "down-down (concordant)",
                                 "up-down (discordant)"),
                      neglog10_p = round(c(RR$max_up, RRd$max_up, RRx$max_up), 3),
                      max_in_map = round(c(RR$max_overall, RRd$max_overall, RRx$max_overall), 3))
  print(RRtab, row.names = FALSE)
  log_lines <<- c(log_lines, capture.output(print(RRtab, row.names = FALSE)))
  fwrite(RRtab, file.path(out, "2b_rrho_corners.csv"))
  fwrite(as.data.table(RR$matrix, keep.rownames = "dengue_rank_cut"),
         file.path(out, "2b_rrho_matrix_up.csv"))

  pdf(file.path(out, "2b_RRHO_heatmap.pdf"), width = 8, height = 8)
  par(family = "Arial")
  pheatmap(RR$matrix, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("white", "#FDD49E", "#FC8D59", "#B30000"))(255),
           main = "Rank Rank Hypergeometric Overlap",
           show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
           fontsize = 8, fontsize_row = 8, fontsize_col = 8, fontfamily = "Arial")
  dev.off()

  # ══════════════════════════════════ 3. EFFECT-SIZE CORRELATION
  say("\n═══ 3. EFFECT-SIZE CORRELATION (log2FC) ═══")
  cor_on <- function(lbl, idx) {
    if (sum(idx) < 10) return(data.table(set = lbl, n = sum(idx), pearson_r = NA_real_,
                                         pearson_p = NA_real_, spearman_rho = NA_real_,
                                         spearman_p = NA_real_, slope = NA_real_))
    x <- M$lfc_D[idx]; y <- M$lfc_Z[idx]
    pe <- cor.test(x, y, method = "pearson"); sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
    data.table(set = lbl, n = sum(idx), pearson_r = round(pe$estimate, 4), pearson_p = pe$p.value,
               spearman_rho = round(sp$estimate, 4), spearman_p = sp$p.value,
               slope = round(coef(lm(y ~ x))[2], 4))
  }
  CO <- rbindlist(list(
    cor_on("all shared tested genes", rep(TRUE, nrow(M))),
    cor_on("DEG in either virus",     M$sig_D | M$sig_Z),
    cor_on("DEG in both viruses",     M$sig_D & M$sig_Z),
    cor_on("concordant genes only",   M$concordant)))
  print(CO, row.names = FALSE)
  log_lines <<- c(log_lines, capture.output(print(CO, row.names = FALSE)))
  fwrite(CO, file.path(out, "3_effect_size_correlation.csv"))
  say("  Spearman is the safer headline: log2FC distributions are heavy-tailed, so Pearson is")
  say("  sensitive to a handful of extreme genes.")

  M[, panel := ifelse(concordant, "concordant", ifelse(discordant, "discordant",
                ifelse(sig_D | sig_Z, "DEG in one only", "not DEG")))]
  lab_pts <- M[concordant == TRUE][order(-abs(lfc_D) - abs(lfc_Z))][1:min(15, sum(M$concordant))]
  ggsave(file.path(out, "3_log2FC_scatter.pdf"),
    ggplot(M, aes(lfc_D, lfc_Z, colour = panel)) +
      geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(alpha = 0.45, size = 1.1) +
      scale_colour_manual(values = c(concordant = "#1B7837", discordant = "#762A83",
                                     `DEG in one only` = "#F4B942", `not DEG` = "grey82"),
                          breaks = c("concordant", "discordant", "DEG in one only", "not DEG")) +
      geom_text_repel(data = lab_pts, aes(label = symbol), size = 2.7, colour = "black",
                      max.overlaps = 30, show.legend = FALSE) +
      labs(title = paste0("Effect-size agreement: ", CMP$tissue),
           subtitle = paste0(CMP$dengue$label, "  vs  ", CMP$zika$label,
                             "\nSpearman rho = ", CO[set == "all shared tested genes"]$spearman_rho,
                             " (all shared), ", CO[set == "DEG in either virus"]$spearman_rho,
                             " (DEG in either)"),
           x = paste0("Dengue shrunken log2FC"), y = paste0("Zika shrunken log2FC"),
           colour = NULL) +
      theme_bw(base_size = 12),
    width = 8.5, height = 7, dpi = DPI)

  summary_rows[[CMP$id]] <- data.table(
    comparison = CMP$id, tissue = CMP$tissue, status = sub(" --.*", "", CMP$quality),
    shared_universe = nrow(M), dengue_DEGs = nD, zika_DEGs = nZ,
    both = sum(M$sig_D & M$sig_Z), concordant = sum(M$concordant),
    discordant = sum(M$discordant),
    hyper_p_all = HY[comparison == "all DEGs (any direction)"]$p_value,
    hyper_fold_all = HY[comparison == "all DEGs (any direction)"]$fold_enrichment,
    sign_pct_agree = ST$pct_agree, sign_p = ST$binom_p,
    rrho_up = RRtab[corner == "up-up (concordant)"]$neglog10_p,
    rrho_discordant = RRtab[corner == "up-down (discordant)"]$neglog10_p,
    spearman_all = CO[set == "all shared tested genes"]$spearman_rho,
    spearman_DEG_either = CO[set == "DEG in either virus"]$spearman_rho)

  # ══════════════════════════════════ 4. PATHWAY ANALYSIS
  say("\n═══ 4. PATHWAY ANALYSIS on the concordant genes ═══")
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    say("  clusterProfiler / org.Hs.eg.db unavailable -- skipped")
  } else if (nrow(concordant) < 10) {
    say("  only ", nrow(concordant), " concordant genes -- too few for enrichment, skipped")
  } else {
    suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
    # The universe MUST be the shared tested set, not all human genes: using the whole genome
    # would report enrichment of anything expressed in this cell type at all.
    univ <- as.character(M$gene_id)
    say("  universe = ", length(univ), " shared tested genes (not the whole genome)")
    for (dir in c("up", "down")) {
      gl <- as.character(concordant[dir_D == dir]$gene_id)
      say("\n  -- concordant ", toupper(dir), " (", length(gl), " genes) --")
      if (length(gl) < 10) { say("     too few, skipped"); next }
      go <- tryCatch(enrichGO(gene = gl, universe = univ, OrgDb = org.Hs.eg.db,
                              keyType = "ENTREZID", ont = "BP",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              readable = TRUE), error = function(e) NULL)
      if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
        gdf <- as.data.frame(go)
        say("     GO:BP terms at padj<0.05: ", nrow(gdf))
        say("     top 8:")
        for (i in seq_len(min(8, nrow(gdf))))
          say("       ", sprintf("%-58s p.adj=%.2e  %s", substr(gdf$Description[i], 1, 58),
                                 gdf$p.adjust[i], gdf$GeneRatio[i]))
        fwrite(gdf, file.path(out, paste0("4_GO_BP_concordant_", dir, ".csv")))
        ggsave(file.path(out, paste0("4_GO_BP_concordant_", dir, "_dotplot.png")),
               enrichplot::dotplot(go, showCategory = 15) +
                 ggtitle(paste0("GO:BP -- concordant ", dir, " genes (", CMP$tissue, ")")),
               width = 9.5, height = 8, dpi = DPI)
      } else say("     no GO:BP term reached padj < 0.05")

      kg <- tryCatch(enrichKEGG(gene = gl, universe = univ, organism = "hsa",
                                pAdjustMethod = "BH", pvalueCutoff = 0.05),
                     error = function(e) NULL)
      if (!is.null(kg) && nrow(as.data.frame(kg)) > 0) {
        kdf <- as.data.frame(kg)
        say("     KEGG pathways at padj<0.05: ", nrow(kdf))
        for (i in seq_len(min(6, nrow(kdf))))
          say("       ", sprintf("%-50s p.adj=%.2e  %s", substr(kdf$Description[i], 1, 50),
                                 kdf$p.adjust[i], kdf$GeneRatio[i]))
        fwrite(kdf, file.path(out, paste0("4_KEGG_concordant_", dir, ".csv")))
      } else say("     no KEGG pathway reached padj < 0.05 (or KEGG lookup unavailable offline)")
    }
  }
}

# ═══════════════════════════════════════════════════════════ overall summary
say("\n\n", strrep("═", 78))
say("  SUMMARY ACROSS COMPARISONS")
say(strrep("═", 78))
S <- rbindlist(summary_rows)
print(S, row.names = FALSE)
log_lines <- c(log_lines, capture.output(print(S, row.names = FALSE)))
fwrite(S, file.path(ROOT, "SUMMARY_all_comparisons.csv"))

say("\n── How to read this ──")
say("  The four tests answer different questions and can disagree:")
say("   1. hypergeometric  -- do the DEG LISTS overlap more than chance? (set-level)")
say("   2. sign test/RRHO  -- do shared genes AGREE IN DIRECTION? (sign-level; RRHO is")
say("                         threshold-free, so it does not inherit the padj/LFC cut-offs)")
say("   3. correlation     -- do the EFFECT SIZES track each other? (magnitude-level)")
say("   4. pathways        -- what BIOLOGY the concordant genes represent")
say("  A large overlap with poor correlation would mean the same genes move but by different")
say("  amounts; strong correlation with small overlap would mean agreement that the DEG")
say("  cut-offs are hiding.")

stamp <- file.path("scripts", "set_png_dpi.R")
if (file.exists(stamp)) {
  say("\n── Figure resolution ──")
  say(paste(system2("Rscript", c(stamp, shQuote(ROOT), DPI), stdout = TRUE, stderr = TRUE),
            collapse = "\n"))
}

writeLines(log_lines, file.path(ROOT, "comparison_log.txt"))
say("\nDone -> ", normalizePath(ROOT))
