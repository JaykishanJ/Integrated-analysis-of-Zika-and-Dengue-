import os
import glob
import re

directory = './scripts'
files = ['PhaseD_MainDownstream.R', 'PhaseD_Supplementary.R']

for fn in files:
    filepath = os.path.join(directory, fn)
    if not os.path.exists(filepath): continue
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    new_content = content
    
    # 1. Clean up weird left-overs like '## 10a.', '## 15a.'
    new_content = re.sub(r'## \d+[a-z]\. ', '## ', new_content)
    
    # 2. Clean up '-B:' -> ''
    new_content = re.sub(r'# -B: ', '# ', new_content)
    new_content = re.sub(r'\\n-B: ', '\\n', new_content)
    
    # 3. Add 600 DPI PNG saving to show_and_save
    if 'bg = "transparent")\n  cat(sprintf(" [plot] %s.svg' in new_content and 'dpi = 600' not in new_content:
        new_content = new_content.replace(
            'bg = "transparent")\n  cat(sprintf(" [plot] %s.svg\\n", file))',
            'bg = "transparent")\n  ggsave(file.path(dir_fig, paste0(file, ".png")), p, width = width, height = height, dpi = 600, bg = "white")\n  cat(sprintf(" [plot] %s.svg and .png\\n", file))'
        )
        
    # 4. Add 600 DPI PNG saving to supp_and_save
    if 'bg = "transparent")\n  cat(sprintf(" [supp] %s.svg' in new_content and 'dpi = 600' not in new_content:
        new_content = new_content.replace(
            'bg = "transparent")\n  cat(sprintf(" [supp] %s.svg\\n", file))',
            'bg = "transparent")\n  ggsave(file.path(dir_supp, paste0(file, ".png")), p, width = width, height = height, dpi = 600, bg = "white")\n  cat(sprintf(" [supp] %s.svg and .png\\n", file))'
        )

    # 5. Add 600 DPI PNG saving to save_pheatmap (which uses grid.draw)
    if 'grid::grid.draw(ph$gtable)\n  dev.off()\n  cat(sprintf(" [plot] %s.svg' in new_content and 'png(' not in new_content:
        new_content = new_content.replace(
            'grid::grid.draw(ph$gtable)\n  dev.off()\n  cat(sprintf(" [plot] %s.svg\\n", file))',
            'grid::grid.draw(ph$gtable)\n  dev.off()\n  png(file.path(dir_fig, paste0(file, ".png")), width = width, height = height, units = "in", res = 600)\n  grid::grid.draw(ph$gtable)\n  dev.off()\n  cat(sprintf(" [plot] %s.svg and .png\\n", file))'
        )
        
    if new_content != content:
        print(f'Updated {fn}')
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
