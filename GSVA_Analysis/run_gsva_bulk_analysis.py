import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy.stats import ttest_ind

# Set style — font config MUST come before sns.set_theme to override seaborn's default
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['figure.dpi'] = 300
sns.set_theme(style="whitegrid", font="Arial")

# Base directories
BASE_DIR = "/sessions/trusting-loving-gauss/mnt/Zika/bulk_results/ZIKA_Bulk_Som"
OUT_DIR = os.path.join(BASE_DIR, "gsva_single_cell_bulk_integration")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load the 23 verified genes from single-cell table
sc_file = "/sessions/trusting-loving-gauss/mnt/Zika/results/phaseD/tables/Step11d_Common_Sig_Genes.csv"
sc_df = pd.read_csv(sc_file)
print(f"Loaded {len(sc_df)} verified common genes from {sc_file}")

# Gene lists
up_18_genes = sc_df[(sc_df['Regulation_Dengue'] == 'UP') & (sc_df['Regulation_Zika'] == 'UP')]['symbol'].tolist()
down_4_genes = sc_df[(sc_df['Regulation_Dengue'] == 'DOWN') & (sc_df['Regulation_Zika'] == 'DOWN')]['symbol'].tolist()
discordant_genes = sc_df[sc_df['Regulation_Dengue'] != sc_df['Regulation_Zika']]['symbol'].tolist()
all_23_genes = sc_df['symbol'].tolist()

print(f"Verified 23 All: {len(all_23_genes)}")
print(f"Verified 18 UP: {len(up_18_genes)}")
print(f"Verified 4 DOWN: {len(down_4_genes)}")
print(f"Discordant 1: {len(discordant_genes)}")

# Gene sets dict for GSVA / ssGSEA
gene_sets = {
    "Verified_23_All": all_23_genes,
    "Verified_18_UP": up_18_genes,
    "Verified_4_DOWN": down_4_genes,
    "Discordant_1": discordant_genes
}

# 2. Define Contrast Files mapping across datasets
contrast_files = {
    "GSE110512__DENV_vs_Mock": {
        "dataset": "GSE110512",
        "file": os.path.join(BASE_DIR, "results/GSE110512/step7_results/results_all_genes.csv"), # Wait, let's check exact path for GSE110512
        "vst": os.path.join(BASE_DIR, "results/GSE110512/step5_qc/vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE110512/step2_metadata/colData.csv"),
        "contrast_name": "GSE110512 (Huh7): DENV vs Mock"
    },
    "GSE80434__DENV_vs_Mock": {
        "dataset": "GSE80434",
        "file": os.path.join(BASE_DIR, "results/GSE80434/step7_results/DENV_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
        "contrast_name": "GSE80434 (hNPC): DENV vs Mock"
    },
    "GSE80434__ZIKVM_vs_Mock": {
        "dataset": "GSE80434",
        "file": os.path.join(BASE_DIR, "results/GSE80434/step7_results/ZIKVM_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
        "contrast_name": "GSE80434 (hNPC): ZIKV-M vs Mock"
    },
    "GSE80434__ZIKVC_vs_Mock": {
        "dataset": "GSE80434",
        "file": os.path.join(BASE_DIR, "results/GSE80434/step7_results/ZIKVC_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
        "contrast_name": "GSE80434 (hNPC): ZIKV-C vs Mock"
    },
    "GSE161783__moDC_Zika_vs_Mock": {
        "dataset": "GSE161783",
        "file": os.path.join(BASE_DIR, "results/GSE161783/step7_results/moDC__Zika_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (moDC): ZIKV Inf vs Mock"
    },
    "GSE161783__moDC_Zika_bystander_vs_Mock": {
        "dataset": "GSE161783",
        "file": os.path.join(BASE_DIR, "results/GSE161783/step7_results/moDC__Zika_bystander_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (moDC): ZIKV Bystander vs Mock"
    },
    "GSE161783__Huh7_5_Zika_vs_Mock": {
        "dataset": "GSE161783",
        "file": os.path.join(BASE_DIR, "results/GSE161783/step7_results/Huh7_5__Zika_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (Huh7.5): ZIKV Inf vs Mock"
    },
    "GSE161783__Huh7_5_Zika_bystander_vs_Mock": {
        "dataset": "GSE161783",
        "file": os.path.join(BASE_DIR, "results/GSE161783/step7_results/Huh7_5__Zika_bystander_vs_Mock_results_all_genes.csv"),
        "vst": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (Huh7.5): ZIKV Bystander vs Mock"
    }
}

# 3. Process Gene-by-Gene Table across Bulk Contrasts
gene_contrast_rows = []
tally_rows = []

for c_key, c_info in contrast_files.items():
    fpath = c_info['file']
    if not os.path.exists(fpath):
        print(f"Warning: file {fpath} does not exist!")
        continue
    
    df_res = pd.read_csv(fpath)
    # Check symbol column name
    sym_col = None
    for candidate in ['symbol', 'gene_symbol', 'Gene', 'gene_name', 'GeneSymbol']:
        if candidate in df_res.columns:
            sym_col = candidate
            break
    if sym_col is None:
        # Check index or Ensembl
        if 'ensembl' in df_res.columns:
            sym_col = 'ensembl' # fallback
        else:
            sym_col = df_res.columns[0]
            
    # Filter for the 23 verified genes
    # Some tables use symbol column or index
    df_sub = df_res[df_res[sym_col].isin(all_23_genes)].copy()
    
    # Ensure log2FC and padj are extracted
    lfc_col = [c for c in df_res.columns if 'log2FoldChange' in c or 'log2FC' in c or 'logFC' in c][0]
    padj_col = [c for c in df_res.columns if 'padj' in c or 'p_adj' in c or 'FDR' in c][0]
    pval_col = [c for c in df_res.columns if 'pvalue' in c or 'p_val' in c or 'PValue' in c][0] if any('pval' in c or 'p_val' in c for c in df_res.columns) else padj_col

    n_up_pval = 0
    n_down_pval = 0
    n_ns_pval = 0
    
    n_up_strict = 0
    n_down_strict = 0
    n_ns_strict = 0

    for gene in all_23_genes:
        g_match = df_sub[df_sub[sym_col] == gene]
        if len(g_match) > 0:
            lfc = g_match[lfc_col].values[0]
            padj = g_match[padj_col].values[0]
            pval = g_match[pval_col].values[0]
            
            # Direction (p_adj < 0.05)
            if pd.notnull(padj) and padj < 0.05:
                reg = 'UP' if lfc > 0 else 'DOWN'
            else:
                reg = 'NS'
                
            # Strict direction (p_adj < 0.05 and |log2FC| > 1)
            if pd.notnull(padj) and padj < 0.05 and abs(lfc) > 1:
                reg_strict = 'UP' if lfc > 0 else 'DOWN'
            else:
                reg_strict = 'NS'
        else:
            lfc = np.nan
            padj = np.nan
            pval = np.nan
            reg = 'NOT_TESTED'
            reg_strict = 'NOT_TESTED'

        # Tally counts
        if reg == 'UP': n_up_pval += 1
        elif reg == 'DOWN': n_down_pval += 1
        elif reg == 'NS': n_ns_pval += 1
        
        if reg_strict == 'UP': n_up_strict += 1
        elif reg_strict == 'DOWN': n_down_strict += 1
        elif reg_strict == 'NS': n_ns_strict += 1

        gene_contrast_rows.append({
            "contrast_key": c_key,
            "contrast_name": c_info['contrast_name'],
            "dataset": c_info['dataset'],
            "symbol": gene,
            "log2FoldChange": lfc,
            "pvalue": pval,
            "padj": padj,
            "Regulation_pval_0.05": reg,
            "Regulation_strict_lfc1": reg_strict
        })
        
    tally_rows.append({
        "contrast_key": c_key,
        "contrast_name": c_info['contrast_name'],
        "dataset": c_info['dataset'],
        "total_verified_genes": len(all_23_genes),
        "UP_padj_0.05": n_up_pval,
        "DOWN_padj_0.05": n_down_pval,
        "NS_padj_0.05": n_ns_pval,
        "UP_strict_lfc1": n_up_strict,
        "DOWN_strict_lfc1": n_down_strict,
        "NS_strict_lfc1": n_ns_strict
    })

# Save Gene-by-Gene Master Table
master_gene_df = pd.DataFrame(gene_contrast_rows)
master_gene_df.to_csv(os.path.join(OUT_DIR, "Verified_23_gene_bulk_quantification.csv"), index=False)

# Save Tally Table
tally_df = pd.DataFrame(tally_rows)
tally_df.to_csv(os.path.join(OUT_DIR, "Verified_23_gene_direction_tally_per_dataset.csv"), index=False)
print("Saved gene-level quantification and tally table successfully.")

# 4. GSVA / ssGSEA Enrichment Analysis on VST matrices
all_gsva_res = []
gsva_contrast_stats = []

datasets_vst = {
    "GSE110512": {
        "vst": os.path.join(BASE_DIR, "results/GSE110512/step5_qc/vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE110512/step2_metadata/colData.csv")
    },
    "GSE80434": {
        "vst": os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv")
    },
    "GSE161783_moDC": {
        "vst": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv")
    },
    "GSE161783_Huh7_5": {
        "vst": os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv")
    }
}

for ds_name, ds_paths in datasets_vst.items():
    if not os.path.exists(ds_paths['vst']) or not os.path.exists(ds_paths['coldata']):
        print(f"Skipping {ds_name}, VST or coldata missing.")
        continue
        
    vst_df = pd.read_csv(ds_paths['vst'], index_col=0)
    col_df = pd.read_csv(ds_paths['coldata'], index_col=0)
    
    # Gene annotation / symbol check
    ann_file = os.path.join(BASE_DIR, f"data/{ds_name.split('_')[0]}/gene_annotation.tsv")
    if os.path.exists(ann_file):
        ann_df = pd.read_csv(ann_file, sep='\t', low_memory=False)
        id_col = [c for c in ann_df.columns if 'id' in c.lower() or 'ensembl' in c.lower() or 'geneid' in c.lower()][0]
        sym_col = [c for c in ann_df.columns if 'symbol' in c.lower() or 'name' in c.lower()][0]
        id_to_sym = dict(zip(ann_df[id_col].astype(str), ann_df[sym_col].astype(str)))
        # Map VST index safely
        vst_df.index = [id_to_sym.get(str(idx), str(idx)) for idx in vst_df.index]
    
    # Remove duplicate gene symbols by keeping max expression
    vst_df = vst_df.groupby(vst_df.index).max()
    
    # Run ssGSEA/GSVA
    try:
        ss_res = gp.ssgsea(data=vst_df, gene_sets=gene_sets, min_size=1, outdir=None)
        es_matrix = ss_res.res2d.pivot(index='Name', columns='Term', values='NES')
        
        # Merge with coldata metadata
        es_matrix['dataset'] = ds_name
        merged_es = es_matrix.join(col_df, how='inner')
        merged_es.to_csv(os.path.join(OUT_DIR, f"GSVA_scores_{ds_name}.csv"))
        
        for term in gene_sets.keys():
            for sample in merged_es.index:
                condition = merged_es.loc[sample, 'condition'] if 'condition' in merged_es.columns else merged_es.loc[sample, 'group'] if 'group' in merged_es.columns else 'Unknown'
                all_gsva_res.append({
                    "dataset": ds_name,
                    "sample": sample,
                    "condition": condition,
                    "term": term,
                    "NES": merged_es.loc[sample, term]
                })
        print(f"Successfully computed GSVA for {ds_name}")
    except Exception as e:
        print(f"Error computing GSVA for {ds_name}: {e}")

gsva_df = pd.DataFrame(all_gsva_res)
gsva_df.to_csv(os.path.join(OUT_DIR, "GSVA_sample_scores_all_datasets.csv"), index=False)

# Compute differential statistics per dataset & contrast
diff_stats_rows = []
for ds_name in gsva_df['dataset'].unique():
    sub_df = gsva_df[gsva_df['dataset'] == ds_name]
    mock_conds = [c for c in sub_df['condition'].unique() if 'mock' in c.lower() or 'control' in c.lower()]
    if not mock_conds:
        continue
    mock_cond = mock_conds[0]
    
    inf_conds = [c for c in sub_df['condition'].unique() if c != mock_cond]
    for inf_cond in inf_conds:
        for term in gene_sets.keys():
            mock_vals = sub_df[(sub_df['condition'] == mock_cond) & (sub_df['term'] == term)]['NES'].values
            inf_vals = sub_df[(sub_df['condition'] == inf_cond) & (sub_df['term'] == term)]['NES'].values
            
            if len(mock_vals) > 0 and len(inf_vals) > 0:
                mean_mock = np.mean(mock_vals)
                mean_inf = np.mean(inf_vals)
                nes_shift = mean_inf - mean_mock
                
                # t-test if samples >= 2
                if len(mock_vals) >= 2 and len(inf_vals) >= 2:
                    t_stat, p_val = ttest_ind(inf_vals, mock_vals, equal_var=False)
                else:
                    t_stat, p_val = np.nan, np.nan
                    
                diff_stats_rows.append({
                    "dataset": ds_name,
                    "infected_condition": inf_cond,
                    "mock_condition": mock_cond,
                    "term": term,
                    "mean_mock_NES": mean_mock,
                    "mean_infected_NES": mean_inf,
                    "NES_shift": nes_shift,
                    "t_statistic": t_stat,
                    "pvalue": p_val
                })

diff_stats_df = pd.DataFrame(diff_stats_rows)
diff_stats_df.to_csv(os.path.join(OUT_DIR, "GSVA_contrast_differential_analysis.csv"), index=False)
print("Saved GSVA contrast differential statistics table successfully.")

# 5. Generate Heatmap of Log2FC across 23 Verified Genes
heatmap_pivot = master_gene_df.pivot(index='symbol', columns='contrast_name', values='log2FoldChange')
hub_8 = ['DDIT3', 'DDIT4', 'CTH', 'CHAC1', 'ASNS', 'DNAJB9', 'HERPUD1', 'MTHFD2']
remaining_up = [g for g in up_18_genes if g in heatmap_pivot.index and g not in hub_8]
ordered_genes = [g for g in hub_8 if g in heatmap_pivot.index] + remaining_up + [g for g in down_4_genes if g in heatmap_pivot.index] + [g for g in discordant_genes if g in heatmap_pivot.index]
heatmap_pivot = heatmap_pivot.reindex(ordered_genes)

plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(12, 10))
ax = sns.heatmap(heatmap_pivot, cmap="coolwarm", center=0, annot=True, fmt=".2f",
                 cbar_kws={'label': 'Expression Log2FC'}, linewidths=0.5, linecolor="white",
                 annot_kws={"size": 8})

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Expression Log2FC', fontproperties={'family': 'Arial', 'size': 9})

for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color("white")

plt.title("Expression Profile of 23 Verified Genes Across Bulk RNA-seq Contrasts", fontsize=10, pad=15, color='black')
plt.xlabel("Bulk RNA-seq Contrast", fontsize=9, labelpad=10, color='black')
plt.ylabel("Verified Gene Symbol", fontsize=9, labelpad=10, color='black')
plt.xticks(rotation=30, ha='right', fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "Verified_23_gene_heatmap_singlecell_and_bulk.pdf"), format='pdf', dpi=600)
plt.close()

# 6. Generate Dedicated GSVA Pathway Enrichment NES Shift Heatmap
if not diff_stats_df.empty:
    diff_stats_df['contrast'] = diff_stats_df['dataset'].str.replace('GSE161783_', '') + ': ' + diff_stats_df['infected_condition']
    gsva_shift_pivot = diff_stats_df.pivot(index='term', columns='contrast', values='NES_shift')
    plt.figure(figsize=(9.5, 4.5), dpi=300)
    ax_gsva = sns.heatmap(gsva_shift_pivot, cmap="coolwarm", center=0, annot=True, fmt=".3f",
                          cbar_kws={'label': 'GSVA NES Shift (Infected - Mock)'}, linewidths=0.75, linecolor="white",
                          annot_kws={"size": 10, "weight": "bold"})
    plt.title("GSVA Pathway Module Enrichment Shifts Across Bulk Infection Contrasts", fontsize=12, pad=15, fontweight='bold', color='#1A365D')
    plt.xlabel("Bulk RNA-Seq Infection Contrast", fontsize=10, labelpad=10, fontweight='bold', color='#2B6CB0')
    plt.ylabel("Core Gene Set Module", fontsize=10, labelpad=10, fontweight='bold', color='#2B6CB0')
    plt.xticks(rotation=25, ha='right', fontsize=9.5, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "GSVA_enrichment_shifts_heatmap.pdf"), format='pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved GSVA enrichment shifts heatmap to: {os.path.join(OUT_DIR, 'GSVA_enrichment_shifts_heatmap.png')}")

# 7. Generate Boxplots of GSVA Enrichment Scores
if not gsva_df.empty:
    plt.figure(figsize=(14, 8))
    g = sns.catplot(
        data=gsva_df[gsva_df['term'].isin(['Verified_18_UP', 'Verified_4_DOWN', 'Verified_23_All'])],
        x='condition', y='NES', hue='term', col='dataset', kind='box',
        palette='Set2', sharey=False, height=5, aspect=1.2
    )
    g.fig.subplots_adjust(top=0.85)
    g.fig.suptitle("GSVA Normalized Enrichment Scores (NES) Across Bulk RNA-seq Conditions", fontsize=16, fontweight='bold')
    g.set_xticklabels(rotation=45, ha='right')
    g.set_axis_labels("Experimental Condition", "GSVA Enrichment Score (NES)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "GSVA_enrichment_boxplots.pdf"), format='pdf')
    plt.close()

print("All GSVA bulk analysis steps completed successfully!")
