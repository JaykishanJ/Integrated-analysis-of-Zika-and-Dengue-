import os
import re

filepath = 'E:/Zika/scripts/PhaseD_Supplementary.R'
with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

new_content = content

# 1. Clean up StepSupp
new_content = re.sub(r'StepSupp', 'Supp', new_content)
new_content = re.sub(r'Step Supp', 'Supp', new_content)

# 2. Fix missing pluses
new_content = new_content.replace('name = "Mean\\nmodule score") \n    theme(axis.text.x', 'name = "Mean\\nmodule score") + \n    theme(axis.text.x')
new_content = new_content.replace('forestgreen")) \n    theme(legend.position', 'forestgreen")) + \n    theme(legend.position')
new_content = new_content.replace('size = 3.5) \n          labs(title = "Permutation', 'size = 3.5) + \n          labs(title = "Permutation')
new_content = new_content.replace('max.overlaps = 15) \n    labs(title = "Stouffer', 'max.overlaps = 15) + \n    labs(title = "Stouffer')

# 3. Add p <- p + theme_pub() to save functions if not there
if 'p <- p + theme_pub()' not in new_content:
    new_content = new_content.replace(
        'show_and_save <- function(p, file, width = 8, height = 5) {\n  print(p)',
        'show_and_save <- function(p, file, width = 8, height = 5) {\n  p <- p + theme_pub()\n  print(p)'
    )
    new_content = new_content.replace(
        'supp_and_save <- function(p, file, width = 8, height = 5) {\n  print(p)',
        'supp_and_save <- function(p, file, width = 8, height = 5) {\n  p <- p + theme_pub()\n  print(p)'
    )

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(new_content)
