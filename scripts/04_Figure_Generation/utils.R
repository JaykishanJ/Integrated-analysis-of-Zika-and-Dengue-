# Utility functions for publication-ready figures
library(ggplot2)

# Define standard DPI and dimensions for publication
PUB_DPI <- 600
FIG_WIDTH_1COL <- 3.54 # inches (90 mm)
FIG_WIDTH_1.5COL <- 5.5  # inches (140 mm)
FIG_WIDTH_2COL <- 7.48 # inches (190 mm)
FIG_HEIGHT_STD <- 5.0 # inches (adjust as needed)

# Define a publication-ready ggplot theme
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(color = "black", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 2),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size, face = "bold"),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# Define custom color palettes
col_palette_condition <- c("Mock" = "#377eb8", "Infected" = "#e41a1c")
# Add more cell type palettes as needed
