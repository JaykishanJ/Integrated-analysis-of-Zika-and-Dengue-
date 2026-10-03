import os

filepath = 'E:/Zika/scripts/Assemble_Workflow_Panel.R'
with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# 1. Update density to 600
new_content = content.replace('density = 300', 'density = 600')

# 2. Update plot_grid for high impact journal standards
new_content = new_content.replace(
    'label_size = 18, ncol = 3)',
    'label_size = 24, label_fontfamily = "sans", label_fontface = "bold", ncol = 3)'
)

# 3. Add PNG and TIFF exports at 600 DPI
add_exports = '''
out_file_png <- "E:/Zika/plots/Comprehensive_A_to_F_Workflow_Panel.png"
out_file_tiff <- "E:/Zika/plots/Comprehensive_A_to_F_Workflow_Panel.tiff"

# Save as raster formats
ggsave(out_file_png, mega_panel, width = 20, height = 14, dpi = 600, bg = "white")
ggsave(out_file_tiff, mega_panel, width = 20, height = 14, dpi = 600, bg = "white", device = "tiff")
'''

if 'out_file_png' not in new_content:
    new_content = new_content.replace(
        'ggsave(out_file_pdf, mega_panel, width = 20, height = 14, bg = "white")',
        'ggsave(out_file_pdf, mega_panel, width = 20, height = 14, bg = "white")\n' + add_exports
    )

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(new_content)
