import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sns.set_theme(style="whitegrid")
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.edgecolor'] = '#CBD5E0'
plt.rcParams['axes.linewidth'] = 0.8

BASE_DIR = "/media/somenath/Phd/Som-zika-bulk_som_aug_2026/ZIKA_Bulk_Som/ZIKA_Bulk_Som"
CSV_PATH = os.path.join(BASE_DIR, "results/gsva_single_cell_bulk_integration/GSVA_sample_scores_all_datasets.csv")
OUT_DIR = os.path.join(BASE_DIR, "results/gsva_single_cell_bulk_integration")

df = pd.read_csv(CSV_PATH)

# Filter for primary modules
df_sub = df[df['term'].isin(['Verified_18_UP', 'Verified_4_DOWN', 'Verified_23_All'])].copy()

# Dataset configurations (display title and output filename)
dataset_configs = {
    "GSE110512": {
        "title": "GSVA Module Enrichment: Human Neural Progenitor Cells (hNPC)\nDengue Virus Infection (GSE110512)",
        "outfile": "GSVA_boxplot_GSE110512_hNPC_DENV.png",
        "palette": ["#4299E1", "#E53E3E", "#DD6B20"]
    },
    "GSE80434": {
        "title": "GSVA Module Enrichment: Human Neural Progenitor Cells (hNPC)\nZika & Dengue Multi-Strain Infection (GSE80434)",
        "outfile": "GSVA_boxplot_GSE80434_hNPC_ZIKV.png",
        "palette": ["#4299E1", "#E53E3E", "#DD6B20", "#319795"]
    },
    "GSE161783_moDC": {
        "title": "GSVA Module Enrichment: Monocyte-Derived Dendritic Cells (moDC)\nZika Virus Infection & Bystander Cells (GSE161783)",
        "outfile": "GSVA_boxplot_GSE161783_moDC_ZIKV.png",
        "palette": ["#4299E1", "#E53E3E", "#DD6B20"]
    },
    "GSE161783_Huh7_5": {
        "title": "GSVA Module Enrichment: Huh7.5 Hepatocytes\nZika Virus Infection & Bystander Cells (GSE161783)",
        "outfile": "GSVA_boxplot_GSE161783_Huh7_5_ZIKV.png",
        "palette": ["#4299E1", "#E53E3E", "#DD6B20"]
    }
}

print("Generating 4 separate GSVA enrichment boxplot figures...")

for ds_name, config in dataset_configs.items():
    sub_ds = df_sub[df_sub['dataset'] == ds_name]
    if sub_ds.empty:
        continue
        
    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=300)
    
    # Create boxplot with clean formatting
    sns.boxplot(
        data=sub_ds,
        x='condition',
        y='NES',
        hue='term',
        palette=["#319795", "#DD6B20", "#718096"],
        ax=ax,
        fliersize=0,
        width=0.6,
        linewidth=1.2,
        boxprops=dict(alpha=0.85)
    )
    
    # Overlay individual data points
    sns.stripplot(
        data=sub_ds,
        x='condition',
        y='NES',
        hue='term',
        dodge=True,
        jitter=0.12,
        marker='o',
        color='#1A202C',
        alpha=0.75,
        size=6,
        linewidth=0.5,
        edgecolor='white',
        ax=ax,
        legend=False
    )
    
    # Expand y-axis limits to prevent touching top/bottom borders
    y_min, y_max = sub_ds['NES'].min(), sub_ds['NES'].max()
    y_range = y_max - y_min if y_max != y_min else 0.5
    ax.set_ylim(y_min - 0.15 * y_range, y_max + 0.20 * y_range)
    
    # Customize Aesthetics
    ax.set_title(config['title'], fontsize=12, fontweight='bold', pad=18, color='#1A365D')
    ax.set_xlabel("Experimental Condition", fontsize=10.5, fontweight='bold', labelpad=10, color='#2B6CB0')
    ax.set_ylabel("GSVA Enrichment Score (NES)", fontsize=10.5, fontweight='bold', labelpad=10, color='#2B6CB0')
    
    plt.xticks(fontsize=10, fontweight='bold', rotation=0, ha='center')
    plt.yticks(fontsize=9.5, fontweight='bold')
    
    # Place Legend OUTSIDE the plot area to guarantee zero collision/overlap
    handles, labels = ax.get_legend_handles_labels()
    # Map term names to publication-friendly labels
    clean_labels = {
        'Verified_23_All': 'All 23 Core Genes',
        'Verified_18_UP': '18 UP Core Genes',
        'Verified_4_DOWN': '4 DOWN Core Genes'
    }
    legend_labels = [clean_labels.get(lbl, lbl) for lbl in labels[:3]]
    
    ax.legend(
        handles[:3], 
        legend_labels, 
        title="Core Gene Module", 
        title_fontsize='9.5', 
        fontsize='9', 
        bbox_to_anchor=(1.02, 1.0), 
        loc='upper left', 
        frameon=True,
        facecolor='#F7FAFC',
        edgecolor='#CBD5E0'
    )
    
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, config['outfile'])
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved cleanly formatted plot to: {out_path}")

print("All 4 GSVA boxplot figures created successfully!")
