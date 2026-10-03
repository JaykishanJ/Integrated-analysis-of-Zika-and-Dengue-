import os
import glob

directory = 'E:/Zika/scripts'
r_files = glob.glob(os.path.join(directory, 'Phase*.R'))
r_files.append(os.path.join(directory, 'SingleCell_FigurePanel.R'))

for filepath in r_files:
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    new_content = content
    
    if 'device = svglite::svglite)' in new_content and 'dpi = 600' not in new_content:
        png_save = '''device = svglite::svglite)
  ggsave(fig(paste0(file, ".png")), p, width = width, height = height,
         dpi = 600, bg = "white")'''
        new_content = new_content.replace('device = svglite::svglite)', png_save)
        
        # also update the cat() message to reflect .png
        new_content = new_content.replace('[plot] %s.svg', '[plot] %s.svg and .png')
        
    if new_content != content:
        print(f'Updating DPI in {os.path.basename(filepath)}')
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
