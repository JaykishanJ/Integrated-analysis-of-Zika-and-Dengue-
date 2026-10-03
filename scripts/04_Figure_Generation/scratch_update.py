import os
import glob
import re

directory = './scripts'
r_files = glob.glob(os.path.join(directory, '*.R'))

for filepath in r_files:
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    new_content = content
    
    # 1. Remove "Step \d+[a-zA-Z]*: " from comments and cat() strings
    new_content = re.sub(r'(cat\([^"\'\n]*["\'])Step \d+[a-zA-Z]*:\s*', r'\1', new_content)
    new_content = re.sub(r'(#+)\s*Step \d+[a-zA-Z]*:\s*', r'\1 ', new_content)
    
    # Also catch "(Steps X-Y)"
    new_content = re.sub(r'\(Steps? \d+-\d+\)', '', new_content)
    
    # 2. Add 600 DPI PNG save to show_and_save()
    # Find the ggsave for .svg and add another one for .png
    if 'device = svglite::svglite)' in new_content and 'dpi = 600' not in new_content:
        png_save = '''device = svglite::svglite)
  ggsave(fig(paste0(file, ".png")), p, width = width, height = height,
         dpi = 600, bg = "white")'''
        new_content = new_content.replace('device = svglite::svglite)', png_save)
        
        # also update the cat() message to reflect .png
        new_content = new_content.replace('[plot] %s.svg', '[plot] %s.svg and .png')
        
    # 3. Update SingleCell_FigurePanel.R dpi=300 to dpi=600
    if 'dpi = 300' in new_content:
        new_content = new_content.replace('dpi = 300', 'dpi = 600')

    if new_content != content:
        print(f'Updating {os.path.basename(filepath)}')
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
