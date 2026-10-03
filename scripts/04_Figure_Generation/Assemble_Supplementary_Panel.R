# =============================================================================
# Script to Assemble a Supplementary Multi-Panel Figure
# =============================================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

# Define paths to supplementary plots (loading the native 600 DPI PNGs)
file_A <- "results/phaseD/supplementary/Supp3a_ViralLoad_ByCluster.png"
file_B <- "results/phaseD/supplementary/Supp3b_StateComposition_ByCluster.png"
file_C <- "results/phaseD/supplementary/Supp3c_ModuleScores_DotPlot.png"
file_D <- "results/phaseD/supplementary/Supp3e_Clustree.png"
file_E <- "results/phaseF/Drug_Repurposing_Results/Figure3_Hub_Gene_Drug_Comprehensive.png"
file_F <- "results/phaseD/supplementary/DEG_Summary_Barplot.png"

cat("Reading 6 supplementary PNG images at native 600 DPI resolution...\n")
pA <- ggdraw() + draw_image(file_A)
pB <- ggdraw() + draw_image(file_B)
pC <- ggdraw() + draw_image(file_C)
pD <- ggdraw() + draw_image(file_D)
pE <- ggdraw() + draw_image(file_E)
pF <- ggdraw() + draw_image(file_F)

cat("Arranging 2x3 supplementary panel...\n")
# Combine into a 2-row, 3-column grid, maintaining journal standards
row1 <- plot_grid(pA, pB, pC, labels = c("A", "B", "C"), label_size = 24, label_fontfamily = "sans", label_fontface = "bold", ncol = 3)
row2 <- plot_grid(pD, pE, pF, labels = c("D", "E", "F"), label_size = 24, label_fontfamily = "sans", label_fontface = "bold", ncol = 3)

# Stack the two rows
mega_panel <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))

# Ensure the plots directory exists for the final output
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

out_file_png <- "plots/Supplementary_Workflow_Panel.png"
out_file_pdf <- "plots/Supplementary_Workflow_Panel.pdf"
out_file_tiff <- "plots/Supplementary_Workflow_Panel.tiff"

cat("Saving massive supplementary panel to PDF, PNG, and TIFF...\n")

# Save as vector format
ggsave(out_file_pdf, mega_panel, width = 20, height = 14, bg = "white")

# Save as raster formats (Nature/Cell requested 600 DPI)
ggsave(out_file_png, mega_panel, width = 20, height = 14, dpi = 600, bg = "white")
ggsave(out_file_tiff, mega_panel, width = 20, height = 14, dpi = 600, bg = "white", device = "tiff")

cat("Done! Rendering in RStudio...\n")
print(mega_panel)

