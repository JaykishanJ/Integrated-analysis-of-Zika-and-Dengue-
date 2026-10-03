import os
import re

filepath = './scripts/PhaseE_Enrichment.R'
with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

new_content = content

# 1. Fix save_pheatmap to 600 DPI
new_content = new_content.replace('units = "in", res = 150)', 'units = "in", res = 600)')

# 2. Fix '-B:' leftover
new_content = new_content.replace('# -B: ', '# ')

# 3. Fix missing pluses before theme()
new_content = re.sub(r'(\n\s*)theme\(', r' + \1theme(', new_content)
# If it resulted in double pluses like `+  + \n`, clean it up
new_content = re.sub(r'\+\s*\+', '+', new_content)
# Wait, let's just do precise replacements to be safe
new_content = content # reset

new_content = new_content.replace('units = "in", res = 150)', 'units = "in", res = 600)')
new_content = new_content.replace('# -B: ', '# ')

lines = new_content.split('\n')
for i in range(len(lines)):
    if 'theme(axis.text.y =' in lines[i] and 'clusterProfiler::' in lines[i-1]:
        if lines[i-1].rstrip().endswith(')'):
            lines[i-1] = lines[i-1].rstrip() + ' +'

new_content = '\n'.join(lines)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(new_content)
