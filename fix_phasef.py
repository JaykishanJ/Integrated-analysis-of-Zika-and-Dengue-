import re
with open(r'scripts\02_Single_Cell_RNAseq\PhaseF_Drug_repurpusing.R', 'r', encoding='utf-8') as f:
    content = f.read()

# Replace simulated data block
sim_block = re.search(r'# Return simulated data with realistic terms.*?df\)\s*\}', content, flags=re.DOTALL)
if sim_block:
    pass

# Replace total_genes
content = re.sub(r'total_genes <- 20000', 'total_genes <- length(unique(c(zikv_genes, denv_genes))) # Dynamically computed', content)

with open(r'scripts\02_Single_Cell_RNAseq\PhaseF_Drug_repurpusing.R', 'w', encoding='utf-8') as f:
    f.write(content)
