import os
import glob
import pandas as pd
import numpy as np
import gseapy as gp
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = "E:/Zika/ZIKA_Bulk_Som/ZIKA_Bulk_Som"
OUT_DIR = "E:/Zika/figures_final/gsva_true_pathways"
os.makedirs(OUT_DIR, exist_ok=True)

# We need the VST matrices for the 3 datasets
vst_files = {
    "GSE110512": os.path.join(BASE_DIR, "results/GSE110512/step5_qc/vst_matrix.csv"),
    "GSE80434": os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
    "GSE161783_Huh7.5": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
    "GSE161783_moDC": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv")
}

# The target pathways
target_pathways = {
    "Response to ER stress": ["DDIT3", "HERPUD1", "DNAJB9", "ATF3", "ERN1", "EIF2AK3", "ATF6", "HSPA5", "XBP1", "DDIT4", "CHOP"],
    "Response to unfolded protein": ["DDIT3", "HERPUD1", "DNAJB9", "HSPA5", "XBP1", "ATF6", "ERN1", "EIF2AK3", "SEC24D"],
    "Integrated stress response": ["DDIT3", "DDIT4", "ATF3", "ATF4", "EIF2S1", "EIF2AK3", "EIF2AK4", "EIF2AK1", "EIF2AK2", "PPP1R15A"],
    "One-carbon metabolism": ["MTHFD2", "MTHFD1", "MTHFD1L", "SHMT1", "SHMT2", "MTHFR", "TYMS", "DHFR"],
    "Aminoacyl-tRNA biosynthesis": ["ASNS", "CARS1", "AARS1", "DARS1", "EARS1", "FARSB", "GARS1", "HARS1", "IARS1", "KARS1", "LARS1", "MARS1"],
    "Glutathione metabolism": ["CHAC1", "CTH", "GCLM", "GCLC", "GSS", "GSR", "GPX1", "GPX4", "GSTP1"],
    "Apoptosis": ["DDIT3", "DDIT4", "CHAC1", "GADD45B", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS", "APAF1"],
    "mTORC1 signaling": ["DDIT4", "SLC7A11", "MTOR", "RPTOR", "AKT1", "TSC1", "TSC2", "RHEB", "EIF4EBP1", "RPS6KB1"]
}

# 1. Read metadata mapping from all experiments
meta_data = []

# GSE110512
meta1 = pd.read_csv(os.path.join(BASE_DIR, "results/GSE110512/step2_metadata/colData.csv"))
meta1['dataset'] = 'GSE110512'
meta1['cell_type'] = 'Huh7'
meta_data.append(meta1)

# GSE80434
meta2 = pd.read_csv(os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"))
meta2['dataset'] = 'GSE80434'
meta2['cell_type'] = 'hNPC'
meta_data.append(meta2)

# GSE161783
meta3 = pd.read_csv(os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/design_feasibility.csv"))
# Wait, GSE161783 colData might not be combined easily, let's just infer from PCA scores
pca_huh = pd.read_csv(os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_pca_scores.csv"))
pca_huh['dataset'] = 'GSE161783'
pca_huh['cell_type'] = 'Huh7.5'
pca_huh['sample_id'] = pca_huh['sample']

pca_modc = pd.read_csv(os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_pca_scores.csv"))
pca_modc['dataset'] = 'GSE161783'
pca_modc['cell_type'] = 'moDC'
pca_modc['sample_id'] = pca_modc['sample']

# 2. Process each dataset
all_ssgsea = []

def run_ssgsea(vst_path, dataset_name, df_inf_path):
    if not os.path.exists(vst_path):
        return
    print(f"Running ssGSEA for {dataset_name}...")
    vst = pd.read_csv(vst_path, index_col=0)
    
    # Map Entrez to Symbol using the DEG results file
    df_inf = pd.read_csv(df_inf_path)
    mapping = dict(zip(df_inf['gene_id'], df_inf['symbol']))
    
    # Convert index to symbols
    vst.index = vst.index.map(mapping)
    vst = vst.loc[~vst.index.isna()]
    vst = vst[~vst.index.duplicated(keep='first')]
    
    # gseapy can sometimes misread index, make an explicit NAME column
    vst = vst.reset_index()
    vst = vst.rename(columns={'gene_id': 'NAME'})
    
    # Run ssGSEA
    ss = gp.ssgsea(data=vst, gene_sets=target_pathways, outdir=None, min_size=1, max_size=2000, n_jobs=4)
    res = ss.res2d
    
    # gseapy ss.res2d is already in long format! Columns: ['Name', 'Term', 'ES', 'NES']
    res_melt = res.rename(columns={'Name': 'sample_id', 'NES': 'NES_score'})
    res_melt = res_melt[['Term', 'sample_id', 'NES_score']]
    res_melt['dataset_group'] = dataset_name
    all_ssgsea.append(res_melt)

run_ssgsea(vst_files["GSE110512"], "GSE110512_Huh7", os.path.join(BASE_DIR, "results/GSE110512/step7_results/results_all_genes.csv"))
run_ssgsea(vst_files["GSE80434"], "GSE80434_hNPC", os.path.join(BASE_DIR, "results/GSE80434/step7_results/DENV_vs_Mock_results_all_genes.csv"))
run_ssgsea(vst_files["GSE161783_Huh7.5"], "GSE161783_Huh7.5", os.path.join(BASE_DIR, "results/GSE161783/step7_results/Huh7_5__Zika_vs_Mock_results_all_genes.csv"))
run_ssgsea(vst_files["GSE161783_moDC"], "GSE161783_moDC", os.path.join(BASE_DIR, "results/GSE161783/step7_results/moDC__Zika_vs_Mock_results_all_genes.csv"))

final_df = pd.concat(all_ssgsea, ignore_index=True)

# Merge conditions
def get_condition(row):
    samp = row['sample_id']
    if samp in meta1['sample_id'].tolist(): return meta1[meta1['sample_id']==samp]['condition'].values[0]
    if samp in meta2['sample_id'].tolist(): return meta2[meta2['sample_id']==samp]['condition'].values[0]
    if samp in pca_huh['sample_id'].tolist(): return pca_huh[pca_huh['sample_id']==samp]['condition'].values[0]
    if samp in pca_modc['sample_id'].tolist(): return pca_modc[pca_modc['sample_id']==samp]['condition'].values[0]
    return "Unknown"

print(final_df['sample_id'].head().tolist())
final_df['condition'] = final_df.apply(get_condition, axis=1)

print("Before replace:")
print(final_df['condition'].value_counts())

# Clean up conditions
final_df['condition'] = final_df['condition'].replace({'ZIKVC': 'Zika', 'ZIKVM': 'Zika', 'DENV': 'Dengue'})
print("After replace:")
print(final_df['condition'].value_counts())

final_df = final_df[final_df['condition'].isin(['Mock', 'Dengue', 'Zika'])]

# Map unified Infected vs Mock
final_df['infection_status'] = final_df['condition'].apply(lambda x: 'Mock' if x == 'Mock' else 'Infected')

final_df.to_csv(os.path.join(OUT_DIR, "True_GSVA_Pathway_Scores.csv"), index=False)
print("Finished ssGSEA. Saved to True_GSVA_Pathway_Scores.csv")
