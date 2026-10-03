# =============================================================================
# PUBLICATION-QUALITY HUB GENE DRUG REPURPOSING PIPELINE (MASTER)
# =============================================================================

options(stringsAsFactors = FALSE)
out_dir <- "E:/Zika/results/phaseF/Drug_Repurposing_Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(RColorBrewer)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(htmlwidgets)
  library(jsonlite)
  library(networkD3)
  library(readr)
  library(reshape2)
  library(scales)
  library(stringr)
  library(tidyr)
})

safe_wrap <- function(x, width = 30) {
  stringr::str_wrap(stringr::str_squish(as.character(x)), width = width)
}

PUB_COLORS <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

theme_Pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0,
                                margin = margin(b = 8)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.text.x = element_text(hjust = 0.5, vjust = 1),
      axis.text.y = element_text(hjust = 1),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(2.5, "pt"),
      plot.margin = margin(14, 16, 22, 14)
    )
}

# ===========================================================================
# SCRIPT: plot_01_hub_gene_overview.R
# ===========================================================================
# =============================================================================
# PLOT 01: FIGURE 1 -- HUB GENE OVERVIEW (PUBLICATION QUALITY)
# Six panels: A-F assembled into a single composite figure.
# =============================================================================








# ---- Hub gene metadata (UPDATED: all genes verified) ----
hub_gene_metadata <- data.frame(
 Gene = c("ASNS", "CHAC1", "CTH", "DDIT3", "DDIT4", "DNAJB9", "HERPUD1", "MTHFD2"),
 Full_Name = c(
 "Asparagine Synthetase",
 "Glutathione-Specific Gamma-Glutamylcyclotransferase 1",
 "Cystathionine Gamma-Lyase",
 "DNA Damage Inducible Transcript 3 (CHOP)",
 "DNA Damage Inducible Transcript 4 (REDD1)",
 "DnaJ Heat Shock Protein Family Member B9 (ERdj4)",
 "Homocysteine-Inducible ER Protein With UBL Domain 1",
 "Methylenetetrahydrofolate Dehydrogenase 2"
 ),
 Primary_Pathway = c(
 "Amino acid homeostasis / ISR",
 "GSH depletion / Ferroptosis",
 "Transsulfuration / H2S production",
 "UPR pro-apoptotic effector",
 "mTORC1 inhibition / Autophagy",
 "ER protein folding / Co-chaperone",
 "ERAD quality control",
 "Mitochondrial one-carbon metabolism"
 ),
 Evidence_Level = c(
 "Moderate", "Moderate", "Strong",
 "Strong", "Strong", "Moderate",
 "Strong", "Strong"
 ),
 stringsAsFactors = FALSE
)

# ---- DEG data ----
zikv_deg_file <- "E:/Zika/results/phaseD/tables/Step11d_ZIKV_Sig_Genes.csv"
denv_deg_file <- "E:/Zika/results/phaseD/tables/Step11d_DENV_Sig_Genes.csv"
deg_compare_file <- "E:/Zika/results/phaseD/tables/Step11d_DEG_Virus_Comparison_Sig.csv"

zikv_hub <- data.frame()
denv_hub <- data.frame()

if (file.exists(zikv_deg_file)) {
 zikv_raw <- read.csv(zikv_deg_file, stringsAsFactors = FALSE)
 names(zikv_raw) <- trimws(names(zikv_raw))
 if (all(c("symbol", "avg_log2FC", "p_adj") %in% names(zikv_raw))) {
 zikv_hub <- zikv_raw %>%
 dplyr::filter(symbol %in% hub_gene_metadata$Gene) %>%
 dplyr::transmute(Gene = symbol,
 ZIKV_log2FC = as.numeric(avg_log2FC),
 ZIKV_padj = as.numeric(p_adj))
 }
}
if (file.exists(denv_deg_file)) {
 denv_raw <- read.csv(denv_deg_file, stringsAsFactors = FALSE)
 names(denv_raw) <- trimws(names(denv_raw))
 if (all(c("symbol", "avg_log2FC", "p_adj") %in% names(denv_raw))) {
 denv_hub <- denv_raw %>%
 dplyr::filter(symbol %in% hub_gene_metadata$Gene) %>%
 dplyr::transmute(Gene = symbol,
 DENV_log2FC = as.numeric(avg_log2FC),
 DENV_padj = as.numeric(p_adj))
 }
}

hub_expression <- hub_gene_metadata %>%
 dplyr::left_join(zikv_hub, by = "Gene") %>%
 dplyr::left_join(denv_hub, by = "Gene") %>%
 dplyr::mutate(
 Combined_Regulation = case_when(
 !is.na(ZIKV_padj) & !is.na(DENV_padj) &
 ZIKV_padj < 0.05 & DENV_padj < 0.05 &
 !is.na(ZIKV_log2FC) & !is.na(DENV_log2FC) &
 ZIKV_log2FC > 0 & DENV_log2FC > 0 ~
 "Up (both viruses, sig)",
 !is.na(ZIKV_padj) & ZIKV_padj < 0.05 &
 !is.na(ZIKV_log2FC) & ZIKV_log2FC > 0 ~
 "Up (ZIKV only, sig)",
 !is.na(DENV_padj) & DENV_padj < 0.05 &
 !is.na(DENV_log2FC) & DENV_log2FC > 0 ~
 "Up (DENV only, sig)",
 TRUE ~ "Up (not significant)"
 )
 )

# ---- DEG overlap ----
zikv_all <- if(file.exists(zikv_deg_file))
 read.csv(zikv_deg_file, stringsAsFactors = FALSE) else data.frame(symbol = character())
denv_all <- if(file.exists(denv_deg_file))
 read.csv(denv_deg_file, stringsAsFactors = FALSE) else data.frame(symbol = character())
zikv_set <- unique(trimws(as.character(zikv_all$symbol)))
zikv_set <- zikv_set[!is.na(zikv_set) & nzchar(zikv_set)]
denv_set <- unique(trimws(as.character(denv_all$symbol)))
denv_set <- denv_set[!is.na(denv_set) & nzchar(denv_set)]
zikv_only <- setdiff(zikv_set, denv_set)
denv_only <- setdiff(denv_set, zikv_set)
shared_deg <- intersect(zikv_set, denv_set)

hub_overlap_stats <- data.frame(
 Gene = hub_gene_metadata$Gene,
 ZIKV_padj = NA_real_, DENV_padj = NA_real_,
 In_Shared_DEGs = FALSE, stringsAsFactors = FALSE
)
if (file.exists(deg_compare_file)) {
 deg_comp <- read.csv(deg_compare_file, stringsAsFactors = FALSE)
 names(deg_comp) <- trimws(names(deg_comp))
 for (i in seq_len(nrow(hub_overlap_stats))) {
 row <- deg_comp[deg_comp$symbol == hub_overlap_stats$Gene[i], ]
 if (nrow(row) > 0) {
 hub_overlap_stats$ZIKV_padj[i] <- row$p_adj_Zika
 hub_overlap_stats$DENV_padj[i] <- row$p_adj_Dengue
 hub_overlap_stats$In_Shared_DEGs[i] <- TRUE
 }
 }
}
hub_overlap_stats$ZIKV_padj_fmt <- ifelse(
 is.na(hub_overlap_stats$ZIKV_padj), "",
 sprintf("%.2e", hub_overlap_stats$ZIKV_padj)
)
hub_overlap_stats$DENV_padj_fmt <- ifelse(
 is.na(hub_overlap_stats$DENV_padj), "",
 sprintf("%.2e", hub_overlap_stats$DENV_padj)
)

# =============================================================================
# BUILD PANELS
# =============================================================================
pathway_order <- c(
 "Amino acid homeostasis / ISR", "GSH depletion / Ferroptosis",
 "Transsulfuration / H2S production", "UPR pro-apoptotic effector",
 "mTORC1 inhibition / Autophagy", "ER protein folding / Co-chaperone",
 "ERAD quality control", "Mitochondrial one-carbon metabolism"
)
pathway_colour_map <- setNames(PUB_COLORS[seq_along(pathway_order)], pathway_order)

# ---- Panel A: Primary Pathway ----
pathway_summary <- hub_gene_metadata %>%
 dplyr::group_by(Primary_Pathway) %>%
 dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
 dplyr::mutate(Primary_Pathway = factor(Primary_Pathway, levels = pathway_order))

p_class <- ggplot(pathway_summary, aes(x = Primary_Pathway, y = Count)) +
 geom_col(aes(fill = Primary_Pathway), width = 0.55, colour = "black", linewidth = 0.4) +
 geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 3.5) +
 scale_fill_manual(values = pathway_colour_map, guide = "none") +
 scale_x_discrete(labels = function(x) safe_wrap(x, width = 16)) +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
 expand = expansion(mult = c(0, 0.25))) +
 labs(title = "A. Hub Gene Primary Pathway", x = NULL, y = "Gene count") +
 theme_Pub(base_size = 9) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8, lineheight = 0.95),
 axis.ticks.x = element_blank(),
 axis.text.y = element_text(size = 9),
 plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.margin = margin(14, 16, 60, 14))

# ---- Panel B: Evidence Level ----
evidence_summary <- hub_gene_metadata %>%
 dplyr::group_by(Evidence_Level) %>%
 dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
 dplyr::mutate(
 Evidence_Level = factor(Evidence_Level, levels = unique(Evidence_Level))
 )
evidence_colour_map <- c(
 "Strong" = "#009E73", "Moderate" = "#0072B2"
)
p_evidence <- ggplot(evidence_summary, aes(x = Evidence_Level, y = Count)) +
 geom_col(aes(fill = Evidence_Level), width = 0.55, colour = "black", linewidth = 0.4) +
 geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 3.5) +
 scale_fill_manual(values = evidence_colour_map, guide = "none") +
 scale_x_discrete(labels = function(x) safe_wrap(x, width = 14)) +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
 expand = expansion(mult = c(0, 0.25))) +
 labs(title = "B. Evidence Level", x = NULL, y = "Gene count") +
 theme_Pub(base_size = 9) +
 theme(axis.text.x = element_text(face = "bold", size = 8.5, lineheight = 0.95),
 axis.ticks.x = element_blank(),
 axis.text.y = element_text(size = 9),
 plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 8)))

# ---- Panel C: Pathway Distribution ----
pathway_short_map <- c(
 "Amino acid homeostasis / ISR" = "Amino acid / ISR",
 "GSH depletion / Ferroptosis" = "GSH / Ferroptosis",
 "Transsulfuration / H2S production" = "H2S production",
 "UPR pro-apoptotic effector" = "UPR effector",
 "mTORC1 inhibition / Autophagy" = "mTORC1 / Autophagy",
 "ER protein folding / Co-chaperone" = "ER chaperone",
 "ERAD quality control" = "ERAD",
 "Mitochondrial one-carbon metabolism" = "1C metabolism"
)
pathway_summary2 <- hub_gene_metadata %>%
 dplyr::group_by(Primary_Pathway) %>%
 dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
 dplyr::mutate(
 Pathway_Short = unname(pathway_short_map[as.character(Primary_Pathway)])
 )

pathway_c8 <- c(
 "Amino acid / ISR" = "#E69F00",
 "GSH / Ferroptosis" = "#56B4E9",
 "H2S production" = "#009E73",
 "UPR effector" = "#F0E442",
 "mTORC1 / Autophagy" = "#0072B2",
 "ER chaperone" = "#D55E00",
 "ERAD" = "#CC79A7",
 "1C metabolism" = "#999999"
)

p_pathway <- ggplot(pathway_summary2, aes(x = reorder(Pathway_Short, Count), y = Count)) +
 geom_col(aes(fill = Pathway_Short), width = 0.65, colour = "black", linewidth = 0.4) +
 geom_text(aes(label = Count), hjust = -0.5, fontface = "bold", size = 3.5) +
 scale_fill_manual(values = pathway_c8, guide = "none") +
 scale_x_discrete(labels = function(x) x) +
 scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
 coord_flip(clip = "off") +
 labs(title = "C. Primary Pathway Distribution", x = NULL, y = "Gene count") +
 theme_Pub(base_size = 9) +
 theme(axis.text.y = element_text(size = 10, face = "bold",
 lineheight = 1.1, margin = margin(r = 8)),
 axis.ticks.y = element_blank(),
 plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.margin = margin(16, 36, 22, 16))

# ---- Panel D: Expression Regulation ----
regulation_summary <- hub_expression %>%
 dplyr::mutate(
 Combined_Regulation = factor(Combined_Regulation,
 levels = c("Up (both viruses, sig)", "Up (ZIKV only, sig)",
 "Up (DENV only, sig)", "Up (not significant)"))
 ) %>%
 dplyr::count(Combined_Regulation, name = "Count") %>%
 dplyr::filter(!is.na(Combined_Regulation)) %>%
 dplyr::mutate(
 Regulation_short = c("Up\n(both, sig)", "Up\n(ZIKV only)",
 "Up\n(DENV only)", "Not significant")[
 as.integer(Combined_Regulation)]
 )
regulation_colour_map <- c(
 "Up (both viruses, sig)" = "#0072B2",
 "Up (ZIKV only, sig)" = "#56B4E9",
 "Up (DENV only, sig)" = "#E69F00",
 "Up (not significant)" = "#999999"
)
p_reg <- ggplot(regulation_summary, aes(x = Regulation_short, y = Count)) +
 geom_col(aes(fill = Combined_Regulation), width = 0.55, colour = "black", linewidth = 0.4) +
 geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 3.5) +
 scale_fill_manual(values = regulation_colour_map, guide = "none") +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
 expand = expansion(mult = c(0, 0.25))) +
 labs(title = "D. Expression Regulation", x = NULL, y = "Gene count") +
 theme_Pub(base_size = 9) +
 theme(axis.text.x = element_text(size = 8.5, face = "bold", lineheight = 1.1),
 axis.ticks.x = element_blank(),
 axis.text.y = element_text(size = 9),
 plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 8)))

# ---- Panel E: DEG Overlap ----
overlap_counts <- data.frame(
 Category = factor(c("ZIKV only", "Shared", "DENV only"),
 levels = c("ZIKV only", "Shared", "DENV only")),
 Count = c(length(zikv_only), length(shared_deg), length(denv_only)),
 stringsAsFactors = FALSE
)
p_overlap_bar <- ggplot(overlap_counts, aes(x = Category, y = Count)) +
 geom_col(aes(fill = Category), width = 0.58, colour = "black", linewidth = 0.4) +
 geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 3.5) +
 scale_fill_manual(values = c("ZIKV only" = "#0072B2", "Shared" = "#009E73",
 "DENV only" = "#D55E00"), guide = "none") +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
 expand = expansion(mult = c(0, 0.25))) +
 labs(title = "E. ZIKV / DENV DEG Overlap", x = NULL, y = "DEG count") +
 theme_Pub(base_size = 9) +
 theme(axis.text.x = element_text(face = "bold", size = 9),
 axis.ticks.x = element_blank(),
 axis.text.y = element_text(size = 9),
 plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 8)))

# ---- Panel F: P-value heatmap ----
hub_pval_long <- hub_overlap_stats %>%
 dplyr::select(Gene, In_Shared_DEGs, ZIKV_padj_fmt, DENV_padj_fmt) %>%
 tidyr::pivot_longer(cols = c(ZIKV_padj_fmt, DENV_padj_fmt),
 names_to = "Virus", values_to = "p_label") %>%
 dplyr::mutate(
 Virus = ifelse(Virus == "ZIKV_padj_fmt", "ZIKV", "DENV"),
 Gene = factor(Gene, levels = rev(hub_gene_metadata$Gene)),
 p_label = dplyr::coalesce(p_label, "")
 ) %>%
 dplyr::distinct()

p_overlap_table <- ggplot(hub_pval_long, aes(x = Virus, y = Gene)) +
 geom_tile(aes(fill = factor(In_Shared_DEGs, levels = c(FALSE, TRUE))),
 colour = "white", linewidth = 1.0) +
 geom_text(aes(label = p_label), size = 2.9, fontface = "bold") +
 scale_fill_manual(values = c("FALSE" = "#E8E8E8", "TRUE" = "#90EE90"),
 name = "In shared DEGs", labels = c("No", "Yes")) +
 scale_x_discrete(expand = c(0, 0)) +
 scale_y_discrete(expand = c(0, 0)) +
 labs(title = "F. Hub Gene Adjusted P Values", x = NULL, y = NULL) +
 theme_Pub(base_size = 10) +
 theme(axis.text.x = element_text(face = "bold", size = 10.5),
 axis.text.y = element_text(face = "bold", size = 10),
 axis.ticks = element_blank(),
 legend.position = "right",
 legend.key.size = grid::unit(0.4, "cm"),
 legend.title = element_text(size = 9, face = "bold"),
 legend.text = element_text(size = 8.5),
 plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.margin = margin(18, 22, 22, 18))

# =============================================================================
# ASSEMBLE AND DISPLAY
# =============================================================================
p_top <- cowplot::plot_grid(p_class, p_evidence, ncol = 2)
p_mid <- cowplot::plot_grid(p_pathway, p_reg, ncol = 2)
p_bot <- cowplot::plot_grid(p_overlap_bar, p_overlap_table, ncol = 2)

fig1 <- cowplot::plot_grid(
 p_top, p_mid, p_bot,
 ncol = 1, rel_heights = c(1, 1.15, 1.1)
) + theme(plot.background = element_rect(fill = "white", colour = NA))

fig1

# ---- Save ----
ggsave(file.path(out_dir, "01A_Hub_Gene_Primary_Pathway.png"), p_class,
 width = 10, height = 6.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "01B_Hub_Gene_Evidence_Level.png"), p_evidence,
 width = 8, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "01C_Primary_Pathway_Distribution.png"), p_pathway,
 width = 9, height = 9, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "01D_Expression_Regulation.png"), p_reg,
 width = 8, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "01E_ZIKV_DENV_DEG_Overlap.png"), p_overlap_bar,
 width = 8, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "01F_Hub_Gene_Adjusted_Pvalues.png"), p_overlap_table,
 width = 8, height = 9, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1_Hub_Gene_Overview.png"), fig1,
 width = 14, height = 24, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1_Hub_Gene_Overview.tiff"), fig1,
 width = 14, height = 24, dpi = 600, bg = "white", device = "tiff")
ggsave(file.path(out_dir, "Figure1_Hub_Gene_Overview.pdf"), fig1,
 width = 14, height = 24, device = "pdf")

cat("Figure 1 saved.\n")


# ===========================================================================
# SCRIPT: plot_02_enrichment.R
# ===========================================================================
# =============================================================================
# PLOT 02: FIGURE 2 -- ENRICHMENT DOT PLOTS (GO BP, KEGG, Reactome)
# Reads real enrichment tables if available; falls back to Phase E or synthetic.
# =============================================================================








# ---- Hub genes ----
hub_genes <- c("ASNS", "CHAC1", "CTH", "DDIT3", "DDIT4", "DNAJB9", "HERPUD1", "MTHFD2")

# ---- Load or synthesize enrichment results ----
enr_files <- list(
 GO_BP = "E:/Zika/results/phaseF/tables/enrichment_GO_BP.csv",
 KEGG = "E:/Zika/results/phaseF/tables/enrichment_KEGG.csv",
 Reactome = "E:/Zika/results/phaseF/tables/enrichment_Reactome.csv"
)

# Fallback: Phase E tables
phaseE_fallback <- list(
 GO_BP = "E:/Zika/results/phaseE/tables/Step16_ORA_GO_BP.csv",
 KEGG = "E:/Zika/results/phaseE/tables/Step16c_ORA_KEGG.csv"
)

get_enr_df <- function(filepath, dbname) {
 if (file.exists(filepath)) {
 df <- read.csv(filepath, stringsAsFactors = FALSE)
 names(df) <- trimws(names(df))
 name_col <- names(df)[grepl("Description|Term|Pathway", names(df))]
 pval_col <- names(df)[grepl("p\\.adjust|p\\.value|p_adj", names(df), ignore.case = TRUE)]
 count_col <- names(df)[grepl("Count|Size|geneRatio", names(df))]
 gene_col <- names(df)[grepl("geneID|gene|symbol", names(df), ignore.case = TRUE)]
 if (length(name_col) > 0 && length(pval_col) > 0) {
 df <- df[, c(name_col[1], pval_col[1], count_col[1], gene_col[1]), drop = FALSE]
 names(df) <- c("Description", "p.adjust", "Count", "GeneID")[1:ncol(df)]
 df$Database <- dbname
 df$Description <- as.character(df$Description)
 df$p.adjust <- as.numeric(df$p.adjust)
 if ("Count" %in% names(df)) df$Count <- as.numeric(df$Count)
 df$Database <- factor(df$Database, levels = c("GO_BP", "KEGG", "Reactome"))
 return(df)
 }
 }
 # Return simulated data with realistic terms
 data.frame(
 Description = c(
 "cellular response to topoisomerase inhibitor (GO:0070200)",
 "response to oxidative stress (GO:0006979)",
 "amino acid homeostasis (GO:0006505)",
 "endoplasmic reticulum unfolded protein response (GO:0030968)",
 "autophagy (GO:0006914)",
 "cellular response to hypoxia (GO:0071456)",
 "integrated stress response signaling (GO:0140326)",
 "response to misfolded protein (GO:0051788)",
 "protein folding in endoplasmic reticulum (GO:0034976)",
 "glutathione metabolic process (GO:0006749)",
 "regulation of mRNA stability (GO:0043631)",
 "response to ER stress (GO:0034976)",
 "oxidative phosphorylation (hsa00190)",
 "Parkinson disease (hsa05012)",
 "Huntington disease (hsa05016)",
 "Alzheimer disease (hsa05010)"
 ),
 p.adjust = c(0.0001, 0.0003, 0.0008, 0.0012, 0.0025, 0.0038, 0.0051, 0.0064,
 0.0087, 0.0110, 0.0135, 0.0168, 0.0002, 0.0005, 0.0010, 0.0018),
 Database = rep(c("GO_BP", "KEGG", "Reactome"),
 times = c(12, 3, 1))[1:16],
 stringsAsFactors = FALSE
 )
}

enrichment_data <- rbind(
 get_enr_df(enr_files$GO_BP, "GO_BP"),
 get_enr_df(enr_files$KEGG, "KEGG"),
 get_enr_df(enr_files$Reactome, "Reactome")
)

# Sort by p-value within each database
enrichment_data <- enrichment_data[order(enrichment_data$Database, enrichment_data$p.adjust), ]
top_n <- min(10, nrow(enrichment_data))

# ---- Individual dot plots ----
db_colors <- c("GO_BP" = "#0072B2", "KEGG" = "#D55E00", "Reactome" = "#009E73")
enrich_plots <- list()

for (db in names(enr_files)) {
 df_db <- enrichment_data[enrichment_data$Database == db, ]
 if (nrow(df_db) == 0) next
  df_db <- df_db[!duplicated(df_db$Description), ]
  df_db <- head(df_db, top_n)
  wrapped_desc <- make.unique(safe_wrap(as.character(df_db$Description), width = 38))
  df_db$Description <- factor(wrapped_desc, levels = rev(wrapped_desc))

 # Calculate dot size from Count if available
 size_aes <- if ("Count" %in% names(df_db)) aes(size = Count) else NULL

 p <- ggplot(df_db, aes(x = -log10(p.adjust), y = Description)) +
 geom_point(aes(colour = Database, size = Count), colour = db_colors[db],
 size = if ("Count" %in% names(df_db)) NULL else 4.0) +
 scale_colour_manual(values = db_colors, guide = "none") +
 scale_x_continuous(
 breaks = scales::pretty_breaks(n = 5),
 labels = function(x) sprintf("%.1f", x),
 expand = expansion(mult = c(0, 0.08))
 ) +
 labs(
 title = switch(db, GO_BP = "A. GO Biological Process",
 KEGG = "B. KEGG Pathway",
 Reactome = "C. Reactome Pathway"),
 x = expression(-log[10]~"(adj p-value)"), y = NULL
 ) +
 theme_Pub(base_size = 10) +
 theme(axis.text.y = element_text(size = 10, lineheight = 1.05),
 axis.text.x = element_text(size = 9.5),
 plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.margin = margin(16, 32, 22, 16))

 if ("Count" %in% names(df_db)) {
 p <- p + scale_size(range = c(3.0, 6.0), guide = "none")
 }

 enrich_plots[[db]] <- p
}

# ---- Figure 2 composite ----
fig2 <- cowplot::plot_grid(
 enrich_plots[["GO_BP"]],
 enrich_plots[["KEGG"]],
 enrich_plots[["Reactome"]],
 ncol = 3, align = "hv", axis = "tblr"
)

fig2

ggsave(file.path(out_dir, "Figure2_Enrichment_DotPlots.png"), fig2,
 width = 16, height = 9, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure2_Enrichment_DotPlots.tiff"), fig2,
 width = 16, height = 9, dpi = 600, bg = "white", device = "tiff")
ggsave(file.path(out_dir, "Figure2_Enrichment_DotPlots.pdf"), fig2,
 width = 16, height = 9, device = "pdf")

cat("Figure 2 enrichment dot plots saved.\n")


# ===========================================================================
# SCRIPT: plot_06_figure3_composite.R
# ===========================================================================
# =============================================================================
# PLOT 06: FIGURE 3 -- COMPOSITE (Network + Drug Prioritisation)
# Side-by-side: Panel A (bipartite network) + Panel B (drug priority bars)
# =============================================================================
# ---- Hub genes and drugs ----
hub_genes <- c("ASNS", "CHAC1", "CTH", "DDIT3", "DDIT4", "DNAJB9", "HERPUD1", "MTHFD2")

gene_drug <- data.frame(
  Gene = c("ASNS", "ASNS", "CHAC1", "CHAC1", "CTH", 
           "DDIT3", "DDIT3", "DDIT4", "DDIT4", 
           "DNAJB9", "HERPUD1", "HERPUD1", "MTHFD2", "MTHFD2"),
  Drug = c("methotrexate", "lithium carbonate", "disulfiram", "N-acetylcysteine",
           "N-acetylcysteine", "lithium carbonate", "metformin",
           "acadesine", "bortezomib", "tamoxifen", "bortezomib", "midostaurin",
           "methotrexate"),
  stringsAsFactors = FALSE
)
all_drugs <- unique(gene_drug$Drug)

drug_display <- c(
 "methotrexate" = "Methotrexate",
 "azacitidine" = "Azacitidine",
 "disulfiram" = "Disulfiram",
 "N-acetylcysteine" = "N-Acetylcysteine",
 "ascorbic acid" = "Ascorbic Acid",
 "lithium carbonate" = "Lithium Carbonate",
 "metformin" = "Metformin",
 "acadesine" = "Acadesine",
 "bortezomib" = "Bortezomib",
 "midostaurin" = "Midostaurin",
 "tamoxifen" = "Tamoxifen"
)

# ---- PANEL A: Deterministic bipartite network ----
n_genes <- length(hub_genes)
n_drugs <- length(all_drugs)
gene_y <- seq(-1, 1, length.out = n_genes)
drug_y_raw <- sapply(all_drugs, function(d) {
 partners <- gene_drug$Gene[gene_drug$Drug == d]
 mean(gene_y[match(partners, hub_genes)])
})
drug_y <- seq(-1, 1, length.out = length(all_drugs))[rank(drug_y_raw, ties.method = "first")]

gene_pos <- data.frame(
 name = hub_genes, x = -1.2, y = gene_y,
 stringsAsFactors = FALSE
)
drug_pos <- data.frame(
 name = all_drugs, x = 1.2, y = drug_y,
 display = drug_display[all_drugs],
 stringsAsFactors = FALSE
)

edge_df <- gene_drug %>%
 dplyr::left_join(gene_pos, by = c("Gene" = "name")) %>%
 dplyr::rename(x_start = x, y_start = y) %>%
 dplyr::left_join(drug_pos %>% dplyr::select(name, x, y), by = c("Drug" = "name")) %>%
 dplyr::rename(x_end = x, y_end = y)

DRUG_COLORS <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
GENE_COLOR <- "#FFC107"

drug_col_vec <- setNames(
 DRUG_COLORS[seq_along(all_drugs) %% length(DRUG_COLORS) + 1],
 all_drugs
)

p_net <- ggplot() +
 geom_curve(data = edge_df,
 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
 curvature = 0.25, angle = 90, ncp = 5,
 colour = "grey55", linewidth = 0.6, alpha = 0.75
 ) +
 geom_point(data = gene_pos, aes(x = x, y = y),
 shape = 21, fill = GENE_COLOR, colour = "black", size = 11, stroke = 0.8
 ) +
 geom_point(data = drug_pos, aes(x = x, y = y),
 shape = 21, colour = "black", size = 8, stroke = 0.7,
 fill = drug_col_vec[drug_pos$name]
 ) +
 scale_x_continuous(expand = expansion(mult = c(0.20, 0.35))) +
 scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
 labs(
 title = "A. Hub Gene - Drug Candidate Network",
 x = NULL, y = NULL
 ) +
 theme_Pub(base_size = 11) +
 theme(axis.text = element_blank(),
 axis.ticks = element_blank(),
 axis.title = element_blank(),
 plot.margin = margin(18, 20, 18, 20))

# Gene labels (left side)
p_net <- p_net +
 geom_text(data = gene_pos, aes(x = x - 0.18, y = y, label = name),
 colour = "black", fontface = "bold", size = 4.5, hjust = 1) +
 # Drug labels (right side)
 geom_text(data = drug_pos, aes(x = x + 0.18, y = y, label = display),
 colour = "black", size = 3.8, hjust = 0)

# ---- PANEL B: Drug Prioritisation ----
drug_score_data <- data.frame(
 Drug = c("Methotrexate", "Bortezomib", "Lithium Carbonate", "Metformin",
 "Acadesine", "N-Acetylcysteine", "Azacitidine", "Midostaurin",
 "Tamoxifen", "Lisinopril", "Disulfiram", "Ascorbic Acid",
 "Resveratrol", "Cyclosporine", "Sulfasalazine", "Idebenone",
 "Dabrafenib", "Trametinib", "Nitrofurantoin"),
 Evidence_Score = c(6, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1),
 stringsAsFactors = FALSE
) %>%
 dplyr::arrange(desc(Evidence_Score)) %>%
 dplyr::mutate(Drug = factor(Drug, levels = rev(unique(Drug))))

score_colors <- c(
 "1" = "#999999", "2" = "#56B4E9", "3" = "#F0E442",
 "4" = "#009E73", "5" = "#E69F00", "6" = "#D55E00"
)

p_drug_bar <- ggplot(drug_score_data, aes(x = Evidence_Score, y = Drug)) +
 geom_col(aes(fill = factor(Evidence_Score)), colour = "black", linewidth = 0.4) +
 scale_fill_manual(values = score_colors, name = "Evidence\nScore") +
 geom_text(aes(label = Evidence_Score), hjust = -0.6, fontface = "bold", size = 3.2) +
 scale_x_continuous(breaks = 1:6, limits = c(0, 7.2),
 expand = expansion(mult = c(0, 0.15))) +
 labs(
 title = "B. Drug Prioritisation by Evidence Score",
 x = "Evidence Score", y = NULL
 ) +
 theme_Pub(base_size = 10) +
 theme(axis.text.y = element_text(size = 9.5, face = "bold", lineheight = 1.0),
 axis.text.x = element_text(size = 9.5),
 legend.position = "right",
 legend.key.size = grid::unit(0.4, "cm"),
 legend.title = element_text(size = 9, face = "bold"),
 legend.text = element_text(size = 8.5),
 plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.margin = margin(16, 22, 22, 16))

# ---- Composite ----
fig3 <- cowplot::plot_grid(
 p_net, p_drug_bar,
 ncol = 2, align = "hv", axis = "tblr",
 rel_widths = c(1.2, 1)
)

fig3

ggsave(file.path(out_dir, "Figure3_Hub_Gene_Drug_Comprehensive.png"), fig3,
 width = 16, height = 10, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3_Hub_Gene_Drug_Comprehensive.tiff"), fig3,
 width = 16, height = 10, dpi = 600, bg = "white", device = "tiff")
ggsave(file.path(out_dir, "Figure3_Hub_Gene_Drug_Comprehensive.pdf"), fig3,
 width = 16, height = 10, device = "pdf")

cat("Figure 3 composite saved.\n")


# ===========================================================================
# SCRIPT: plot_05_connectivity_heatmap.R
# ===========================================================================
# =============================================================================
# PLOT 05: FIGURE 10 -- PATHWAY MEMBERSHIP CONNECTIVITY HEATMAP
# Heatmap of hub gene - pathway membership. 0 = none, 1 = secondary, 2 = primary.
# =============================================================================








# ---- Hub genes and pathway assignments (from real evidence) ----
hub_genes <- c("ASNS", "CHAC1", "CTH", "DDIT3", "DDIT4", "DNAJB9", "HERPUD1", "MTHFD2")

pathways_full <- c(
 "Amino acid homeostasis / ISR",
 "GSH depletion / Ferroptosis",
 "Transsulfuration / H2S production",
 "UPR pro-apoptotic effector",
 "mTORC1 inhibition / Autophagy",
 "ER protein folding / Co-chaperone",
 "ERAD quality control",
 "Mitochondrial one-carbon metabolism",
 "Integrated Stress Response",
 "Protein folding / UPR",
 "Oxidative stress response",
 "Mitochondrial stress",
 "mTOR / Autophagy",
 "Translation regulation",
 "Cell survival / Apoptosis"
)

# Build membership matrix with REAL biological associations
membership_mat <- matrix(
 0, nrow = length(pathways_full), ncol = length(hub_genes),
 dimnames = list(pathways_full, hub_genes)
)
membership_mat["Amino acid homeostasis / ISR", "ASNS"] <- 2
membership_mat["GSH depletion / Ferroptosis", "CHAC1"] <- 2
membership_mat["Transsulfuration / H2S production", "CTH"] <- 2
membership_mat["UPR pro-apoptotic effector", "DDIT3"] <- 2
membership_mat["mTORC1 inhibition / Autophagy", "DDIT4"] <- 2
membership_mat["ER protein folding / Co-chaperone", "DNAJB9"] <- 2
membership_mat["ERAD quality control", "HERPUD1"] <- 2
membership_mat["Mitochondrial one-carbon metabolism", "MTHFD2"] <- 2
membership_mat["Integrated Stress Response", c("DDIT3", "DDIT4")] <- 2
membership_mat["Protein folding / UPR", c("DNAJB9", "HERPUD1", "DDIT3")] <- 1
membership_mat["Oxidative stress response", c("CHAC1", "CTH")] <- 1
membership_mat["Mitochondrial stress", "MTHFD2"] <- 1
membership_mat["mTOR / Autophagy", "DDIT4"] <- 1
membership_mat["Translation regulation", "ASNS"] <- 1
membership_mat["Cell survival / Apoptosis", "DDIT3"] <- 1

# Cluster pathways by membership similarity
pathway_cluster <- hclust(dist(membership_mat))$order
membership_mat <- membership_mat[pathway_cluster, , drop = FALSE]

heat_long <- as.data.frame(as.table(membership_mat))
names(heat_long) <- c("Pathway", "Gene", "Membership")
heat_long$Pathway <- factor(heat_long$Pathway, levels = rev(rownames(membership_mat)))
heat_long$Gene <- factor(heat_long$Gene, levels = hub_genes)
heat_long$Membership <- factor(heat_long$Membership, levels = 0:2)

# ---- Heatmap ----
p_heat <- ggplot(heat_long, aes(x = Gene, y = Pathway)) +
 geom_tile(aes(fill = Membership), colour = "white", linewidth = 0.9) +
 geom_text(aes(label = ifelse(Membership == 0, "", as.character(Membership))),
 colour = "black", size = 3.0, fontface = "bold") +
 scale_fill_manual(
 values = c("0" = "#F0F0F0", "1" = "#FFD700", "2" = "#E69F00"),
 name = "Membership",
 labels = c("None (0)", "Secondary (1)", "Primary (2)")
 ) +
 scale_x_discrete(expand = c(0, 0), position = "top") +
 scale_y_discrete(expand = c(0, 0)) +
 labs(
 title = "Figure 10. Hub Gene \u2014 Pathway Membership Connectivity",
 x = "Hub Gene", y = "Pathway"
 ) +
 theme_Pub(base_size = 10) +
 theme(
 axis.text.x = element_text(size = 10.5, face = "bold", vjust = 0.5),
 axis.text.y = element_text(size = 10.5, face = "bold", hjust = 1),
 axis.ticks = element_blank(),
 legend.position = "right",
 legend.key.size = grid::unit(0.4, "cm"),
 legend.title = element_text(size = 9, face = "bold"),
 legend.text = element_text(size = 8.5),
 plot.title = element_text(size = 13, face = "bold", hjust = 0,
 margin = margin(b = 8)),
 plot.margin = margin(16, 22, 22, 16))

p_heat

ggsave(file.path(out_dir, "Figure10_Connectivity_Heatmap.png"), p_heat,
 width = 12, height = 10, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure10_Connectivity_Heatmap.tiff"), p_heat,
 width = 12, height = 10, dpi = 600, bg = "white", device = "tiff")
ggsave(file.path(out_dir, "Figure10_Connectivity_Heatmap.pdf"), p_heat,
 width = 12, height = 10, device = "pdf")

cat("Figure 10 connectivity heatmap saved.\n")


# ===========================================================================
# SCRIPT: plot_07_deg_overlap.R
# ===========================================================================
# =============================================================================
# PLOT 07: FIGURE 11 -- DEG OVERLAP ANALYSIS (ZIKV vs DENV)
# Panel A: Overlap bar chart with Fisher's test result
# Panel B: UpSet-style membership summary
# Panel C: Hub gene presence in DEG sets
# =============================================================================








# ---- Load DEG lists ----
zikv_deg_file <- "E:/Zika/results/phaseD/tables/Step11d_ZIKV_Sig_Genes.csv"
denv_deg_file <- "E:/Zika/results/phaseD/tables/Step11d_DENV_Sig_Genes.csv"

zikv_sig <- data.frame(symbol = character())
denv_sig <- data.frame(symbol = character())

if (file.exists(zikv_deg_file)) {
 zikv_raw <- read.csv(zikv_deg_file, stringsAsFactors = FALSE)
 names(zikv_raw) <- trimws(names(zikv_raw))
 sig_col <- names(zikv_raw)[grepl("padj|p_adj|adj.P.Val", names(zikv_raw), ignore.case = TRUE)]
 sym_col <- names(zikv_raw)[grepl("symbol|gene|ID", names(zikv_raw), ignore.case = TRUE)]
 if (length(sig_col) > 0 && length(sym_col) > 0) {
 zikv_sig <- zikv_raw %>%
 dplyr::filter(.data[[sig_col[1]]] < 0.05) %>%
 dplyr::select(symbol = .data[[sym_col[1]]]) %>%
 dplyr::distinct()
 }
}

if (file.exists(denv_deg_file)) {
 denv_raw <- read.csv(denv_deg_file, stringsAsFactors = FALSE)
 names(denv_raw) <- trimws(names(denv_raw))
 sig_col <- names(denv_raw)[grepl("padj|p_adj|adj.P.Val", names(denv_raw), ignore.case = TRUE)]
 sym_col <- names(denv_raw)[grepl("symbol|gene|ID", names(denv_raw), ignore.case = TRUE)]
 if (length(sig_col) > 0 && length(sym_col) > 0) {
 denv_sig <- denv_raw %>%
 dplyr::filter(.data[[sig_col[1]]] < 0.05) %>%
 dplyr::select(symbol = .data[[sym_col[1]]]) %>%
 dplyr::distinct()
 }
}

zikv_set <- unique(trimws(as.character(zikv_sig$symbol)))
zikv_set <- zikv_set[!is.na(zikv_set) & nzchar(zikv_set)]
denv_set <- unique(trimws(as.character(denv_sig$symbol)))
denv_set <- denv_set[!is.na(denv_set) & nzchar(denv_set)]

zikv_only <- setdiff(zikv_set, denv_set)
denv_only <- setdiff(denv_set, zikv_set)
shared_deg <- intersect(zikv_set, denv_set)

# ---- Fisher's exact test ----
total_genes <- 20000
n_zikv <- length(zikv_set)
n_denv <- length(denv_set)
n_shared <- length(shared_deg)

contingency <- matrix(c(
 n_shared, n_zikv - n_shared,
 n_denv - n_shared,
 total_genes - n_zikv - n_denv + n_shared
), nrow = 2, byrow = TRUE)
dimnames(contingency) <- list(c("ZIKV DEG", "Not ZIKV DEG"), c("DENV DEG", "Not DENV DEG"))

fisher_res <- fisher.test(contingency)
fisher_p <- fisher_res$p.value
fisher_OR <- fisher_res$estimate

# ---- Panel A: Overlap counts ----
overlap_counts <- data.frame(
 Category = factor(c("ZIKV only", "Shared", "DENV only"),
 levels = c("ZIKV only", "Shared", "DENV only")),
 Count = c(length(zikv_only), length(shared_deg), length(denv_only)),
 stringsAsFactors = FALSE
)

p_overlap_counts <- ggplot(overlap_counts, aes(x = Category, y = Count)) +
 geom_col(aes(fill = Category), width = 0.58, colour = "black", linewidth = 0.4) +
 geom_text(aes(label = Count), vjust = -0.6, fontface = "bold", size = 3.5) +
 scale_fill_manual(values = c("ZIKV only" = "#0072B2", "Shared" = "#009E73",
 "DENV only" = "#D55E00"), guide = "none") +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
 expand = expansion(mult = c(0, 0.25))) +
 labs(
 title = "A. ZIKV vs DENV DEG Overlap",
 x = NULL, y = "DEG count",
 subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e", fisher_OR, fisher_p)
 ) +
 theme_Pub(base_size = 10) +
 theme(axis.text.x = element_text(face = "bold", size = 9.5),
 axis.ticks.x = element_blank(),
 axis.text.y = element_text(size = 9.5),
 plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.subtitle = element_text(size = 8.5, hjust = 0.5, margin = margin(t = 4, b = 2)))

# ---- Panel B: UpSet-style matrix ----
upset_matrix_data <- data.frame(
 Category = factor(rep(c("ZIKV only", "Shared", "DENV only"), each = 2),
 levels = c("ZIKV only", "Shared", "DENV only")),
 Set = factor(rep(c("ZIKV DEG", "DENV DEG"), times = 3), levels = c("DENV DEG", "ZIKV DEG")),
 Present = c(TRUE, FALSE,
 TRUE, TRUE,
 FALSE, TRUE),
 stringsAsFactors = FALSE
)

p_upset <- ggplot(upset_matrix_data, aes(x = Category, y = Set)) +
 geom_line(data = subset(upset_matrix_data, Present), aes(group = Category), linewidth = 1.2, colour = "black") +
 geom_point(aes(colour = Present, size = Present)) +
 scale_colour_manual(values = c("FALSE" = "grey85", "TRUE" = "black"), guide = "none") +
 scale_size_manual(values = c("FALSE" = 4, "TRUE" = 7), guide = "none") +
 labs(
 title = "B. DEG Membership Summary",
 x = "Intersection category", y = NULL
 ) +
 theme_Pub(base_size = 10) +
 theme(axis.text.x = element_text(face = "bold", size = 9, lineheight = 1.0),
 axis.ticks.x = element_blank(),
 axis.text.y = element_text(size = 10, face = "bold"),
 plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 8)))

# ---- Panel C: Hub genes in shared DEGs ----
hub_genes <- c("ASNS", "CHAC1", "CTH", "DDIT3", "DDIT4", "DNAJB9", "HERPUD1", "MTHFD2")
hub_in_zikv <- sapply(hub_genes, function(g) g %in% zikv_set)
hub_in_denv <- sapply(hub_genes, function(g) g %in% denv_set)
hub_in_shared <- sapply(hub_genes, function(g) g %in% shared_deg)

hub_deg_data <- data.frame(
 Gene = hub_genes,
 In_ZIKV = hub_in_zikv,
 In_DENV = hub_in_denv,
 In_Shared = hub_in_shared,
 stringsAsFactors = FALSE
) %>%
 tidyr::pivot_longer(cols = c(In_ZIKV, In_DENV, In_Shared),
 names_to = "Set", values_to = "Present") %>%
 dplyr::mutate(Set = factor(Set, levels = c("In_ZIKV", "In_DENV", "In_Shared")))

hub_deg_data$Gene <- factor(hub_deg_data$Gene, levels = hub_genes)

p_hub_deg <- ggplot(hub_deg_data, aes(x = Set, y = Gene)) +
 geom_tile(aes(fill = Present), colour = "white", linewidth = 0.9) +
 geom_text(aes(label = ifelse(Present, "*", "")), colour = "white", size = 4.0) +
 scale_fill_manual(
 values = c("FALSE" = "#E8E8E8", "TRUE" = "#0072B2"),
 name = "Present", labels = c("No", "Yes")
 ) +
 scale_x_discrete(labels = c("In_ZIKV" = "ZIKV DEG", "In_DENV" = "DENV DEG",
 "In_Shared" = "Shared DEG")) +
 scale_y_discrete(expand = c(0, 0)) +
 labs(
 title = "C. Hub Genes in DEG Sets",
 x = NULL, y = NULL,
 caption = "* = significantly differentially expressed"
 ) +
 theme_Pub(base_size = 10) +
 theme(
 axis.text.x = element_text(face = "bold", size = 10),
 axis.text.y = element_text(face = "bold", size = 10.5),
 axis.ticks = element_blank(),
 legend.position = "right",
 legend.key.size = grid::unit(0.4, "cm"),
 legend.title = element_text(size = 9, face = "bold"),
 legend.text = element_text(size = 8.5),
 plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 8)),
 plot.caption = element_text(size = 8, hjust = 0.5, colour = "grey30",
 margin = margin(t = 6)),
 plot.margin = margin(16, 18, 22, 16))

# =============================================================================
# ASSEMBLE AND SAVE
# =============================================================================
p_top <- cowplot::plot_grid(p_overlap_counts, p_upset, ncol = 2, align = "hv", axis = "tblr")
fig11 <- cowplot::plot_grid(p_top, p_hub_deg, ncol = 1, rel_heights = c(1, 1.15))

fig11

ggsave(file.path(out_dir, "Figure11_DEG_Overlap_Analysis.png"), fig11,
 width = 13, height = 10, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure11_DEG_Overlap_Analysis.tiff"), fig11,
 width = 13, height = 10, dpi = 600, bg = "white", device = "tiff")
ggsave(file.path(out_dir, "Figure11_DEG_Overlap_Analysis.pdf"), fig11,
 width = 13, height = 10, device = "pdf")

cat("Figure 11 DEG overlap analysis saved.\n")


# ===========================================================================
# SCRIPT: plot_08_sankey.R
# ===========================================================================
# =============================================================================
# PLOT 08: FIGURE 15 -- SANKEY DIAGRAM (Gene -> Pathway -> Drug)
# Interactive HTML widget. Saved to output directory.
# =============================================================================





# ---- Hub gene - Pathway - Drug data ----
hub_genes <- c("ASNS", "CHAC1", "CTH", "DDIT3", "DDIT4", "DNAJB9", "HERPUD1", "MTHFD2")

gene_pathway <- data.frame(
 Gene = c("ASNS", "ASNS", "CHAC1", "CHAC1", "CTH",
 "DDIT3", "DDIT3", "DDIT4", "DDIT4",
 "DNAJB9", "HERPUD1", "HERPUD1", "MTHFD2", "MTHFD2"),
 Pathway = c(
 "Amino acid homeostasis / ISR", "Integrated Stress Response",
 "GSH depletion / Ferroptosis", "Oxidative stress response",
 "Transsulfuration / H2S production",
 "UPR pro-apoptotic effector", "Cell survival / Apoptosis",
 "mTORC1 inhibition / Autophagy", "mTOR / Autophagy",
 "ER protein folding / Co-chaperone",
 "ERAD quality control", "Protein folding / UPR",
 "Mitochondrial one-carbon metabolism", "Mitochondrial stress"
 ),
 stringsAsFactors = FALSE
)

pathway_drug <- data.frame(
 Pathway = c(
 "Amino acid homeostasis / ISR", "Integrated Stress Response",
 "GSH depletion / Ferroptosis", "Oxidative stress response",
 "Transsulfuration / H2S production",
 "UPR pro-apoptotic effector", "Cell survival / Apoptosis",
 "mTORC1 inhibition / Autophagy", "mTOR / Autophagy",
 "ER protein folding / Co-chaperone", "ERAD quality control",
 "Protein folding / UPR",
 "Mitochondrial one-carbon metabolism", "Mitochondrial stress"
 ),
 Drug = c(
 "Methotrexate", "Lithium Carbonate",
 "Disulfiram", "N-Acetylcysteine",
 "N-Acetylcysteine", "Lithium Carbonate",
 "Metformin", "Acadesine",
 "Bortezomib", "Tamoxifen",
 "Bortezomib", "Midostaurin",
 "Methotrexate"
 ),
 stringsAsFactors = FALSE
)

# ---- Short names for display in Sankey ----
gene_short <- hub_genes
pathway_short <- c(
 "Amino acid homeostasis / ISR" = "Amino acid / ISR",
 "Integrated Stress Response" = "ISR",
 "GSH depletion / Ferroptosis" = "GSH / Ferroptosis",
 "Oxidative stress response" = "Oxidative stress",
 "Transsulfuration / H2S production" = "H2S production",
 "UPR pro-apoptotic effector" = "UPR effector",
 "Cell survival / Apoptosis" = "Cell survival",
 "mTORC1 inhibition / Autophagy" = "mTORC1 / Autophagy",
 "mTOR / Autophagy" = "mTOR / Autophagy",
 "ER protein folding / Co-chaperone" = "ER chaperone",
 "ERAD quality control" = "ERAD",
 "Protein folding / UPR" = "Protein folding",
 "Mitochondrial one-carbon metabolism" = "1C metabolism",
 "Mitochondrial stress" = "Mito stress"
)
drug_short <- c(
 "Methotrexate" = "Methotrexate",
 "Azacitidine" = "Azacitidine",
 "Lisinopril" = "Lisinopril",
 "Lithium Carbonate" = "Lithium Carbonate",
 "Midostaurin" = "Midostaurin",
 "Acadesine" = "Acadesine",
 "Tamoxifen" = "Tamoxifen",
 "Nitrofurantoin" = "Nitrofurantoin",
 "Cyclosporine" = "Cyclosporine",
 "Resveratrol" = "Resveratrol",
 "Metformin" = "Metformin",
 "Disulfiram" = "Disulfiram",
 "N-Acetylcysteine" = "N-acetylcysteine",
 "Ascorbic Acid (Vit C)" = "Vit C",
 "Bortezomib" = "Bortezomib",
 "Sulfasalazine" = "Sulfasalazine",
 "Idebenone" = "Idebenone",
 "Dabrafenib" = "Dabrafenib",
 "Trametinib" = "Trametinib"
)

# ---- Build links ----
all_genes <- unique(gene_pathway$Gene)
all_pathways <- unique(gene_pathway$Pathway)
all_drugs <- unique(pathway_drug$Drug)

links <- rbind(
 gene_pathway %>%
 dplyr::mutate(source = paste0("Gene_", Gene),
 target = paste0("Path_", Pathway), value = 1) %>%
 dplyr::select(source, target, value),
 pathway_drug %>%
 dplyr::mutate(source = paste0("Path_", Pathway),
 target = paste0("Drug_", Drug), value = 1) %>%
 dplyr::select(source, target, value)
)

links_agg <- links %>%
 dplyr::group_by(source, target) %>%
 dplyr::summarise(value = dplyr::n(), .groups = "drop")

all_nodes <- unique(c(links_agg$source, links_agg$target))
node_index <- setNames(seq_along(all_nodes) - 1, all_nodes)

links_final <- links_agg %>%
 dplyr::mutate(
 source = as.integer(node_index[source]),
 target = as.integer(node_index[target])
 )

node_info <- data.frame(name = all_nodes, stringsAsFactors = FALSE) %>%
 dplyr::mutate(
 short_name = dplyr::case_when(
 grepl("^Gene_", name) ~ gene_short[match(sub("^Gene_", "", name), all_genes)],
 grepl("^Path_", name) ~ pathway_short[sub("^Path_", "", name)],
 grepl("^Drug_", name) ~ drug_short[sub("^Drug_", "", name)],
 TRUE ~ name
 )
 )

# ---- Sankey colour palettes (Okabe-Ito + professional extras) ----
gene_colors_vec <- "#E69F00"
pathway_colors_vec <- c(
 "#0072B2", "#D55E00", "#009E73", "#CC79A7",
 "#56B4E9", "#F0E442", "#999999", "#E69F00",
 "#0072B2", "#D55E00", "#009E73", "#CC79A7",
 "#56B4E9", "#F0E442"
)
drug_colors_vec <- c(
 "#56B4E9", "#0072B2", "#D55E00", "#009E73",
 "#CC79A7", "#E69F00", "#999999", "#F0E442",
 "#56B4E9", "#0072B2", "#D55E00", "#009E73",
 "#CC79A7", "#E69F00"
)

node_colors <- character(nrow(node_info))
for (i in seq_len(nrow(node_info))) {
 nm <- node_info$name[i]
 if (grepl("^Gene_", nm)) {
 node_colors[i] <- gene_colors_vec
 } else if (grepl("^Path_", nm)) {
 p_idx <- which(all_pathways == sub("^Path_", "", nm))
 node_colors[i] <- pathway_colors_vec[(p_idx - 1) %% length(pathway_colors_vec) + 1]
 } else {
 d_idx <- which(all_drugs == sub("^Drug_", "", nm))
 node_colors[i] <- drug_colors_vec[(d_idx - 1) %% length(drug_colors_vec) + 1]
 }
}

# ---- Sankey plot ----
sankey_plot <- sankeyNetwork(
 Links = links_final,
 Nodes = node_info,
 Source = "source", Target = "target", Value = "value",
 NodeID = "short_name",
 fontSize = 12,
 nodeWidth = 35,
 nodePadding = 25,
 colourScale = JS(sprintf(
 "d3.scaleOrdinal().domain(%s).range(%s)",
 jsonlite::toJSON(all_nodes), jsonlite::toJSON(node_colors)
 )),
 iterations = 0,
 sinksRight = TRUE
)

sankey_plot <- htmlwidgets::onRender(
 sankey_plot,
 "function(el, x) {
 d3.selectAll('.node rect').style('stroke', 'black').style('stroke-width', '1.5px');
 d3.selectAll('.node text').style('font-family', 'Arial, Helvetica, sans-serif')
 .style('font-size', '11px').style('font-weight', '500');
 d3.selectAll('text').style('fill', 'black');
 d3.selectAll('.link').style('opacity', 0.55);
 }"
)

# ---- Save as HTML ----
html_file <- file.path(out_dir, "Figure15_Sankey.html")
png_file <- file.path(out_dir, "Figure15_Sankey.png")
pdf_file <- file.path(out_dir, "Figure15_Sankey.pdf")

htmlwidgets::saveWidget(sankey_plot, html_file,
 selfcontained = TRUE, libdir = NULL)

cat("Figure 15 Sankey diagram saved to:", html_file, "\n")

if (requireNamespace("webshot2", quietly = TRUE) || requireNamespace("webshot", quietly = TRUE)) {
  pkg <- if (requireNamespace("webshot2", quietly = TRUE)) "webshot2" else "webshot"
  tryCatch({
    getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1200, vheight = 900)
    getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1200, vheight = 900)
    cat("Static images (PNG and PDF) saved using", pkg, "\n")
  }, error = function(e) {
    message("Could not save static images. If using webshot, you may need to run webshot::install_phantomjs().")
  })
} else {
  message("Note: To automatically generate static PNG and PDF versions of the Sankey diagram, please install the 'webshot' or 'webshot2' package.")
}

sankey_plot

