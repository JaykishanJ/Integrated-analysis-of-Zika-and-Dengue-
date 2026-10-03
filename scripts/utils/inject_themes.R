# Bridge script: inject_themes.R
# This sources 92_figure_theme.R and provides backward-compatible aliases
# so existing scripts that source this file continue to work.
#
# Usage in scripts: source("./scripts/inject_themes.R")

source("./scripts/92_figure_theme.R")

# ---- Backward-compatible aliases ----
my_pub_theme <- function(base_size = 10, base_family = "Helvetica") {
  theme_pub(base_size = base_size, base_family = base_family)
}

apply_pub_theme <- apply_pub_theme

save_pub_plot <- save_pub_plot

set_pub_geom_defaults <- set_pub_geom_defaults

scale_x_log10_pub <- scale_x_log10_pub
scale_y_log10_pub <- scale_y_log10_pub
