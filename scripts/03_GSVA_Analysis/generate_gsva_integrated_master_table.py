import os
import pandas as pd
import numpy as np

BASE_DIR = "/media/somenath/Phd/Som-zika-bulk_som_aug_2026/ZIKA_Bulk_Som/ZIKA_Bulk_Som"
OUT_DIR = os.path.join(BASE_DIR, "results", "gsva_single_cell_bulk_integration")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load Single-Cell Baseline (23 Verified Genes)
sc_file = os.path.join(BASE_DIR, "single_cell/DEG_Tables/Step11d_Common_Sig_Genes.csv")
df_sc = pd.read_csv(sc_file)

# Load Top 25 cross-tissue PDF core genes
top_25_list = ['ATF3', 'GADD45B', 'DDIT3', 'TRIB3', 'CHAC1', 'RELB', 'MTHFD2', 'TNFAIP3', 'ADM2', 'HERPUD1', 
               'STC2', 'SLC7A11', 'INHBE', 'IRF1', 'CARS1', 'SESN2', 'MAFF', 'ASNS', 'NIBAN1', 'CTH', 
               'TSPYL2', 'DNAJB9', 'DDIT4', 'BBC3', 'ALDH1L2']

# Combined unique 40 genes (23 scRNA-seq + 25 PDF cross-tissue)
combined_symbols = list(set(df_sc['symbol'].tolist() + top_25_list))

df_master = pd.DataFrame({'symbol': combined_symbols})

# Merge Single-Cell metrics
df_master = pd.merge(df_master, df_sc[['symbol', 'ensembl', 'avg_log2FC_Dengue', 'p_adj_Dengue', 'Regulation_Dengue', 
                                       'avg_log2FC_Zika', 'p_adj_Zika', 'Regulation_Zika']], 
                     on='symbol', how='left')

# 2. Define Bulk Contrasts
contrasts = {
    'hNPC_DENV': os.path.join(BASE_DIR, 'results/GSE80434/step7_results/DENV_vs_Mock_results_all_genes.csv'),
    'hNPC_ZIKVM': os.path.join(BASE_DIR, 'results/GSE80434/step7_results/ZIKVM_vs_Mock_results_all_genes.csv'),
    'hNPC_ZIKVC': os.path.join(BASE_DIR, 'results/GSE80434/step7_results/ZIKVC_vs_Mock_results_all_genes.csv'),
    'Huh7_DENV': os.path.join(BASE_DIR, 'results/GSE110512/step7_results/results_all_genes.csv'),
    'Huh75_ZIKV': os.path.join(BASE_DIR, 'results/GSE161783/step7_results/Huh7_5__Zika_vs_Mock_results_all_genes.csv'),
    'moDC_ZIKV': os.path.join(BASE_DIR, 'results/GSE161783/step7_results/moDC__Zika_vs_Mock_results_all_genes.csv')
}

for c_prefix, fpath in contrasts.items():
    if os.path.exists(fpath):
        df_bulk = pd.read_csv(fpath)
        lfc_col = 'log2FC_shrunk' if 'log2FC_shrunk' in df_bulk.columns else 'log2FoldChange'
        sub = df_bulk[['symbol', lfc_col, 'padj']].copy()
        sub.columns = ['symbol', f'log2FC_{c_prefix}', f'padj_{c_prefix}']
        
        # Add discrete regulation status column
        def get_reg(row):
            lfc = row[f'log2FC_{c_prefix}']
            p = row[f'padj_{c_prefix}']
            if pd.isna(p) or p >= 0.05:
                return 'NS'
            return 'UP' if lfc > 0 else 'DOWN'
            
        sub[f'Regulation_{c_prefix}'] = sub.apply(get_reg, axis=1)
        df_master = pd.merge(df_master, sub, on='symbol', how='left')

# 3. Add summary boolean flags
df_master['Is_Single_Cell_Core_23'] = df_master['symbol'].isin(df_sc['symbol'])
df_master['Is_Cross_Tissue_Core_25'] = df_master['symbol'].isin(top_25_list)
df_master['Is_Overlap_Hub_8'] = df_master['Is_Single_Cell_Core_23'] & df_master['Is_Cross_Tissue_Core_25']

# 4. Save Master CSV for 23 Verified Single-Cell Genes
df_23_master = df_master[df_master['Is_Single_Cell_Core_23']].copy().sort_values(by='symbol')
fout_23 = os.path.join(OUT_DIR, "Verified_23_gene_master_cross_tissue_table.csv")
df_23_master.to_csv(fout_23, index=False)
print(f"Saved {fout_23} ({len(df_23_master)} genes)")

# 5. Save Combined Master CSV for all 40 Unique Core Genes
fout_40 = os.path.join(OUT_DIR, "Combined_40_unique_core_genes_master_table.csv")
df_master.sort_values(by=['Is_Overlap_Hub_8', 'Is_Single_Cell_Core_23'], ascending=False).to_csv(fout_40, index=False)
print(f"Saved {fout_40} ({len(df_master)} genes)")

# 6. Save Dedicated 17 Universal Concordant UP Genes Master Table
up_17_list = ['HERPUD1', 'DNAJB9', 'ASNS', 'CTH', 'CHAC1', 'DDIT4', 'MTHFD2', 'DDIT3', 
              'GARS1', 'SELENOS', 'CDK2AP2', 'SELENOK', 'PSAT1', 'HSPA5', 'DNAJC3', 'TARS1', 'SHMT2']
df_17 = df_master[df_master['symbol'].isin(up_17_list)].sort_values(by='avg_log2FC_Dengue', ascending=False)
fout_17 = os.path.join(OUT_DIR, "Universal_17_concordant_UP_genes_master_table.csv")
df_17.to_csv(fout_17, index=False)
print(f"Saved {fout_17} ({len(df_17)} genes)")
