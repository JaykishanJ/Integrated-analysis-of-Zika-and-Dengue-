import os
import glob
import re

directory = 'E:/Zika/scripts'
r_files = glob.glob(os.path.join(directory, '*.R'))

for filepath in r_files:
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Replace 'title = "Step XX: ' with 'title = "'
    new_content = re.sub(r'(title\s*=\s*[\"\'])Step[^:\"]*:\s*', r'\1', content)
    
    if new_content != content:
        print(f'Updating {os.path.basename(filepath)}')
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
