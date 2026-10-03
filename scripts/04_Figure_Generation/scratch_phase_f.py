import os, re

files = ['PhaseF_Drug_repurpusing.R', 'PhaseF_1_Drug_repurpusing.R']
directory = './scripts'

for fn in files:
    filepath = os.path.join(directory, fn)
    if not os.path.exists(filepath): continue
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    new_content = content
    # Remove Step X:
    new_content = re.sub(r'(?i)step\s*\d+[a-z]*:\s*', '', new_content)
    new_content = re.sub(r'(?i)step\s*\d+[a-z]*\b', '', new_content)
    new_content = re.sub(r'(?i)\(?steps?\s*\d+\s*-\s*\d+\)?', '', new_content)
    
    # 600 DPI fix
    new_content = new_content.replace('units = "in", res = 150)', 'units = "in", res = 600)')
    if 'bg = "transparent")\n  cat(sprintf(" [plot] %s.svg' in new_content and 'dpi = 600' not in new_content:
        new_content = new_content.replace(
            'bg = "transparent")\n  cat(sprintf(" [plot] %s.svg\\n", file))',
            'bg = "transparent")\n  ggsave(file.path(dir_fig, paste0(file, ".png")), p, width = width, height = height, dpi = 600, bg = "white")\n  cat(sprintf(" [plot] %s.svg and .png\\n", file))'
        )

    if new_content != content:
        print(f'Updated {fn}')
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
