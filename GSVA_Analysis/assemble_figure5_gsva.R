# Assemble Figure 5 from GSVA plots
library(magick)
library(cowplot)
library(ggplot2)

in_dir <- "E:/Zika/output_manuscript_figures"
out_dir <- "E:/Zika/figures_final"

cat("Reading PDFs...\n")
p1 <- image_read_pdf(file.path(in_dir, "GSVA_enrichment_shifts_heatmap.pdf"), density = 300)
p2 <- image_read_pdf(file.path(in_dir, "GSVA_enrichment_boxplots.pdf"), density = 300)

cat("Converting to ggplot objects...\n")
g1 <- ggdraw() + draw_image(p1)
g2 <- ggdraw() + draw_image(p2)

cat("Assembling final plot...\n")
final_plot <- plot_grid(g1, g2, labels = c("A", "B"), ncol = 1, label_size = 20, rel_heights=c(1, 1.2))

cat("Saving...\n")
ggsave(file.path(out_dir, "Figure5_GSVA_Pathway_Analysis.png"), final_plot, width = 16, height = 18, dpi = 300, bg="white")
ggsave(file.path(out_dir, "Figure5_GSVA_Pathway_Analysis.pdf"), final_plot, width = 16, height = 18, bg="white")

cat("Done!\n")
