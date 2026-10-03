import re

def fix_python_script(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except:
        return
    
    # Remove thioridazine from lists
    content = content.replace(", 'thioridazine'", "")
    content = content.replace("'thioridazine',", "")
    content = content.replace(", 'Thioridazine'", "")
    content = content.replace("'Thioridazine',", "")
    
    # Remove from dictionaries
    content = re.sub(r"'thioridazine':\s*[^,}]+,?", "", content)
    
    # DDIT3 mapping (remove thioridazine from DDIT3)
    content = content.replace("'DDIT3': ['thioridazine', 'lithium carbonate']", "'DDIT3': ['lithium carbonate']")
    
    # Figure 9 -> Supplementary Figure S7 in text and titles
    content = content.replace('Figure 9:', 'Supplementary Figure S7:')
    content = content.replace('Figure 9A', 'Supplementary Figure S7A')
    content = content.replace('Figure 9B', 'Supplementary Figure S7B')
    content = content.replace('Hub gene', 'module gene')
    content = content.replace('Hub Gene', 'module gene')
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

fix_python_script(r'E:\Zika\all_Python_scripts\build_figure9_drug_network.py')
fix_python_script(r'E:\Zika\output_manuscript_figures\build_figure9_drug_network.py')

def fix_r_script(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except:
        return
        
    content = content.replace('"thioridazine", ', '')
    content = content.replace('"thioridazine"', '')
    content = content.replace('"Thioridazine", ', '')
    content = content.replace(', "Thioridazine"', '')
    content = content.replace('"Thioridazine"', '')
    
    # Remove the row from gene_drug dataframe
    content = re.sub(r'\"DDIT3\"(?:,\s*\"DDIT3\")?', '\"DDIT3\"', content)
    content = content.replace('"thioridazine", "lithium carbonate"', '"lithium carbonate"')
    
    # Also adjust score list for Thioridazine (it had score 1 at the end)
    content = re.sub(r'Evidence_Score = c\(6, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1\)', 
                     'Evidence_Score = c(6, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1)', content)
    
    content = content.replace('Hub gene', 'module gene')
    content = content.replace('Hub Gene', 'module gene')
    
    # xlim fix for truncating labels
    content = content.replace('expand = expansion(mult = c(0.40, 0.45))', 'expand = expansion(mult = c(0.50, 0.55))')
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

for r_file in [r'E:\Zika\scripts\PhaseF_1_Drug_repurpusing.R', 
               r'E:\Zika\scripts\PhaseF_Drug_repurpusing.R',
               r'E:\Zika\output_manuscript_figures\R_scripts_consolidated\PhaseF_1_Drug_repurpusing.R',
               r'E:\Zika\output_manuscript_figures\R_scripts_consolidated\PhaseF_Drug_repurpusing.R',
               r'E:\Zika\all_R_scripts\PhaseF_1_Drug_repurpusing.R',
               r'E:\Zika\all_R_scripts\PhaseF_Drug_repurpusing.R']:
    fix_r_script(r_file)

print("done")
