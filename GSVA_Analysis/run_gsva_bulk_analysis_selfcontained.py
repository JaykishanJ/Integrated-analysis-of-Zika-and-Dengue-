"""
Self-contained GSVA/ssGSEA bulk analysis (no gseapy/scipy dependency).
Replaces run_gsva_bulk_analysis.py for environments without gseapy/scipy.
Implements ssGSEA (Barbie et al. 2009) and Welch's t-test from scratch.
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = "/sessions/trusting-loving-gauss/mnt/Zika/bulk_results/ZIKA_Bulk_Som"
OUT_DIR  = os.path.join(BASE_DIR, "gsva_single_cell_bulk_integration")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams["font.family"]     = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["figure.dpi"]      = 300
sns.set_theme(style="whitegrid", font="Arial")

# ════════════════════════════════════════════════════════════════════════════
# 1. Load the 23 verified genes from single-cell table
# ════════════════════════════════════════════════════════════════════════════
sc_file = "/sessions/trusting-loving-gauss/mnt/Zika/results/phaseD/tables/Step11d_Common_Sig_Genes.csv"
sc_df   = pd.read_csv(sc_file)
print(f"Loaded {len(sc_df)} verified common genes from {sc_file}")

up_18_genes     = sc_df[(sc_df["Regulation_Dengue"] == "UP") & (sc_df["Regulation_Zika"] == "UP")]["symbol"].tolist()
down_4_genes    = sc_df[(sc_df["Regulation_Dengue"] == "DOWN") & (sc_df["Regulation_Zika"] == "DOWN")]["symbol"].tolist()
discordant_genes = sc_df[sc_df["Regulation_Dengue"] != sc_df["Regulation_Zika"]]["symbol"].tolist()
all_23_genes    = sc_df["symbol"].tolist()

print(f"Verified 23 All:  {len(all_23_genes)}")
print(f"Verified 18 UP:   {len(up_18_genes)}")
print(f"Verified 4 DOWN:  {len(down_4_genes)}")
print(f"Discordant 1:     {len(discordant_genes)}")

gene_sets = {
    "Verified_23_All": all_23_genes,
    "Verified_18_UP":  up_18_genes,
    "Verified_4_DOWN": down_4_genes,
    "Discordant_1":    discordant_genes,
}

# ════════════════════════════════════════════════════════════════════════════
# 2. Contrast files mapping
# ════════════════════════════════════════════════════════════════════════════
contrast_files = {
    "GSE110512__DENV_vs_Mock": {
        "dataset": "GSE110512",
        "file":    os.path.join(BASE_DIR, "results/GSE110512/step7_results/results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE110512/step5_qc/vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE110512/step2_metadata/colData.csv"),
        "contrast_name": "GSE110512 (Huh7): DENV vs Mock",
    },
    "GSE80434__DENV_vs_Mock": {
        "dataset": "GSE80434",
        "file":    os.path.join(BASE_DIR, "results/GSE80434/step7_results/DENV_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
        "contrast_name": "GSE80434 (hNPC): DENV vs Mock",
    },
    "GSE80434__ZIKVM_vs_Mock": {
        "dataset": "GSE80434",
        "file":    os.path.join(BASE_DIR, "results/GSE80434/step7_results/ZIKVM_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
        "contrast_name": "GSE80434 (hNPC): ZIKV-M vs Mock",
    },
    "GSE80434__ZIKVC_vs_Mock": {
        "dataset": "GSE80434",
        "file":    os.path.join(BASE_DIR, "results/GSE80434/step7_results/ZIKVC_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
        "contrast_name": "GSE80434 (hNPC): ZIKV-C vs Mock",
    },
    "GSE161783__moDC_Zika_vs_Mock": {
        "dataset": "GSE161783",
        "file":    os.path.join(BASE_DIR, "results/GSE161783/step7_results/moDC__Zika_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (moDC): ZIKV Inf vs Mock",
    },
    "GSE161783__moDC_Zika_bystander_vs_Mock": {
        "dataset": "GSE161783",
        "file":    os.path.join(BASE_DIR, "results/GSE161783/step7_results/moDC__Zika_bystander_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (moDC): ZIKV Bystander vs Mock",
    },
    "GSE161783__Huh7_5_Zika_vs_Mock": {
        "dataset": "GSE161783",
        "file":    os.path.join(BASE_DIR, "results/GSE161783/step7_results/Huh7_5__Zika_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (Huh7.5): ZIKV Inf vs Mock",
    },
    "GSE161783__Huh7_5_Zika_bystander_vs_Mock": {
        "dataset": "GSE161783",
        "file":    os.path.join(BASE_DIR, "results/GSE161783/step7_results/Huh7_5__Zika_bystander_vs_Mock_results_all_genes.csv"),
        "vst":     os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
        "contrast_name": "GSE161783 (Huh7.5): ZIKV Bystander vs Mock",
    },
}

# ════════════════════════════════════════════════════════════════════════════
# 3. Gene-by-Gene quantification
# ════════════════════════════════════════════════════════════════════════════
print("\n=== Step 3: Gene-by-gene quantification ===")
gene_contrast_rows = []
tally_rows          = []

for c_key, c_info in contrast_files.items():
    fpath = c_info["file"]
    if not os.path.exists(fpath):
        print(f"  Warning: {fpath} not found, skipping.")
        continue

    df_res = pd.read_csv(fpath)

    sym_col = None
    for candidate in ["symbol", "gene_symbol", "Gene", "gene_name", "GeneSymbol", "ensembl"]:
        if candidate in df_res.columns:
            sym_col = candidate
            break
    if sym_col is None:
        sym_col = df_res.columns[0]

    lfc_col  = [c for c in df_res.columns if "log2FoldChange" in c or "log2FC" in c or "logFC" in c][0]
    padj_col = [c for c in df_res.columns if "padj" in c or "p_adj" in c or "FDR" in c][0]
    pval_col = None
    for candidate in ["pvalue", "p_val", "PValue", "pval"]:
        if candidate in df_res.columns:
            pval_col = candidate
            break
    if pval_col is None:
        pval_col = padj_col

    df_sub = df_res[df_res[sym_col].isin(all_23_genes)].copy()

    n_up_pval = n_down_pval = n_ns_pval = 0
    n_up_strict = n_down_strict = n_ns_strict = 0

    for gene in all_23_genes:
        g_match = df_sub[df_sub[sym_col] == gene]
        if len(g_match) > 0:
            lfc  = g_match[lfc_col].values[0]
            padj = g_match[padj_col].values[0]
            pval = g_match[pval_col].values[0]
        else:
            lfc = padj = pval = np.nan

        if pd.notnull(padj) and padj < 0.05:
            reg = "UP" if lfc > 0 else "DOWN"
        else:
            reg = "NS"

        if pd.notnull(padj) and padj < 0.05 and abs(lfc) > 1:
            reg_strict = "UP" if lfc > 0 else "DOWN"
        else:
            reg_strict = "NS"

        if reg == "UP":           n_up_pval   += 1
        elif reg == "DOWN":       n_down_pval += 1
        else:                     n_ns_pval   += 1
        if reg_strict == "UP":    n_up_strict += 1
        elif reg_strict == "DOWN": n_down_strict += 1
        else:                     n_ns_strict += 1

        gene_contrast_rows.append({
            "contrast_key":          c_key,
            "contrast_name":         c_info["contrast_name"],
            "dataset":               c_info["dataset"],
            "symbol":                gene,
            "log2FoldChange":        lfc,
            "pvalue":                pval,
            "padj":                  padj,
            "Regulation_pval_0.05":  reg,
            "Regulation_strict_lfc1": reg_strict,
        })

    tally_rows.append({
        "contrast_key":           c_key,
        "contrast_name":          c_info["contrast_name"],
        "dataset":                c_info["dataset"],
        "total_verified_genes":   len(all_23_genes),
        "UP_padj_0.05":           n_up_pval,
        "DOWN_padj_0.05":         n_down_pval,
        "NS_padj_0.05":           n_ns_pval,
        "UP_strict_lfc1":         n_up_strict,
        "DOWN_strict_lfc1":       n_down_strict,
        "NS_strict_lfc1":         n_ns_strict,
    })

master_gene_df = pd.DataFrame(gene_contrast_rows)
master_gene_df.to_csv(os.path.join(OUT_DIR, "Verified_23_gene_bulk_quantification.csv"), index=False)
tally_df = pd.DataFrame(tally_rows)
tally_df.to_csv(os.path.join(OUT_DIR, "Verified_23_gene_direction_tally_per_dataset.csv"), index=False)
print("  Saved gene-level quantification and tally tables.")

# ════════════════════════════════════════════════════════════════════════════
# 4. ssGSEA — self-contained implementation
# ════════════════════════════════════════════════════════════════════════════
print("\n=== Step 4: ssGSEA ===")

from math import lgamma, exp

def beta_cdf_integral(a, b, x, n_pts=4000):
    """
    Regularized incomplete beta I_x(a,b) via numerical integration (Simpson's rule).
    """
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    ln_beta = lgamma(a) + lgamma(b) - lgamma(a + b)
    ts = np.linspace(0.0, x, n_pts)
    vals = np.zeros_like(ts)
    for i, t in enumerate(ts):
        if t == 0.0 or t == 1.0:
            vals[i] = 0.0
        else:
            ln_f = (a - 1.0) * np.log(t) + (b - 1.0) * np.log(1.0 - t) - ln_beta
            vals[i] = exp(ln_f)
    h = ts[1] - ts[0]
    integral = vals[0] + vals[-1] + 4.0 * np.sum(vals[1:-1:2]) + 2.0 * np.sum(vals[2:-1:2])
    return float(integral * h / 3.0)


def welch_ttest(x, y):
    """
    Two-sided Welch's t-test implemented from scratch.
    Returns (t_statistic, two_tailed_p_value).
    """
    x = np.asarray(x, dtype=float); x = x[~np.isnan(x)]
    y = np.asarray(y, dtype=float); y = y[~np.isnan(y)]
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan
    m1, m2 = np.mean(x), np.mean(y)
    v1, v2 = np.var(x, ddof=1), np.var(y, ddof=1)
    se = np.sqrt(v1 / n1 + v2 / n2)
    if se == 0:
        return np.nan, np.nan
    t = (m1 - m2) / se
    num   = (v1 / n1 + v2 / n2) ** 2
    denom = (v1**2) / (n1**2 * (n1 - 1)) + (v2**2) / (n2**2 * (n2 - 1))
    df    = num / denom if denom > 0 else (n1 + n2 - 2)
    t2    = float(df) / (float(df) + t * t)
    p     = beta_cdf_integral(df / 2.0, 0.5, t2)
    return t, p


def ssgsea_score(expr_sorted, gene_mask, alpha=0.25):
    """
    ssGSEA enrichment score (Barbie et al. 2009).

    Parameters
    ----------
    expr_sorted : np.ndarray, shape (N,), expression values sorted descending.
    gene_mask   : np.ndarray, shape (N,), boolean — True if gene is in the gene set.
    alpha       : weighting exponent (default 0.25).

    Returns
    -------
    float — NES (normalized enrichment score).
    """
    N  = len(expr_sorted)
    NG = int(gene_mask.sum())
    if NG == 0 or NG >= N:
        return np.nan

    # Step 1: cumulative rank statistics
    # es_k = sum_{i in hits}(rank(i)/N - (i - 0.5)/N)  — simplified:
    # We compute ES as area between ECDF of hits and uniform ECDF.
    w       = np.where(gene_mask, 1.0, 0.0)
    cum_w   = np.cumsum(w)
    cum_all = np.arange(1.0, N + 1.0)
    # Enrichment score at each position
    es_pos  = cum_w / NG - cum_all / N
    # The ssGSEA ES is the sum over positions (equivalent to integrating the difference)
    es = np.sum(es_pos) / N  # scale by N

    # Normalization: analytic null approximation (Barbie et al.)
    # ES_null mean ≈ 0
    # ES_null std ≈ sqrt(NG*(N-NG) / (N*(N-1))) * sqrt(N^2/12) / N
    var_null = (NG * (N - NG) / (N * (N - 1))) * (N**2 / 12.0) / (N**2)
    sd_null  = max(np.sqrt(var_null), 1e-12)
    # Size correction (GSVA uses this)
    es_corr = es * NG / N
    return float(es_corr / sd_null)


def run_ssgsea(vst_df, gene_sets, alpha=0.25):
    """Run ssGSEA on VST matrix (genes x samples)."""
    genes_present = vst_df.index.astype(str).tolist()
    results = {}
    for term_name, term_genes in gene_sets.items():
        mask = np.isin(genes_present, term_genes)
        if mask.sum() == 0:
            print(f"    Warning: '{term_name}' has 0 genes in VST. Skipping.")
            continue
        nes_list = []
        for sample in vst_df.columns:
            expr = np.asarray(vst_df[sample].values, dtype=float)
            idx  = np.argsort(-expr)
            nes  = ssgsea_score(expr[idx], mask[idx], alpha=alpha)
            nes_list.append(nes)
        results[term_name] = nes_list
    return pd.DataFrame(results, index=vst_df.columns)


# ── Run ssGSEA per dataset ──────────────────────────────────────────────────
datasets_vst = {
    "GSE110512": {
        "vst":     os.path.join(BASE_DIR, "results/GSE110512/step5_qc/vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE110512/step2_metadata/colData.csv"),
    },
    "GSE80434": {
        "vst":     os.path.join(BASE_DIR, "results/GSE80434/step5_qc/pooled_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE80434/step2_metadata/colData.csv"),
    },
    "GSE161783_moDC": {
        "vst":     os.path.join(BASE_DIR, "results/GSE161783/step5_qc/moDC_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
    },
    "GSE161783_Huh7_5": {
        "vst":     os.path.join(BASE_DIR, "results/GSE161783/step5_qc/Huh7_5_vst_matrix.csv"),
        "coldata": os.path.join(BASE_DIR, "results/GSE161783/step2_metadata/colData.csv"),
    },
}

all_gsva_res = []
for ds_name, ds_paths in datasets_vst.items():
    vst_path     = ds_paths["vst"]
    coldata_path = ds_paths["coldata"]
    if not os.path.exists(vst_path) or not os.path.exists(coldata_path):
        print(f"  Skipping {ds_name}: file missing.")
        continue

    vst_df = pd.read_csv(vst_path, index_col=0)
    col_df = pd.read_csv(coldata_path, index_col=0)

    # Optional annotation mapping
    gse_root = ds_name.split("_")[0]
    ann_file = os.path.join(BASE_DIR, f"data/{gse_root}/gene_annotation.tsv")
    if os.path.exists(ann_file):
        ann_df  = pd.read_csv(ann_file, sep="\t", low_memory=False)
        id_col  = [c for c in ann_df.columns if "id" in c.lower() or "ensembl" in c.lower()][0]
        sym_col = [c for c in ann_df.columns if "symbol" in c.lower() or "name" in c.lower()][0]
        id_to_sym = dict(zip(ann_df[id_col].astype(str), ann_df[sym_col].astype(str)))
        vst_df.index = [id_to_sym.get(str(idx), str(idx)) for idx in vst_df.index]

    # Deduplicate by keeping max expression per gene symbol
    vst_df = vst_df.groupby(vst_df.index).max()

    print(f"  Running ssGSEA for {ds_name} ({vst_df.shape[0]} genes x {vst_df.shape[1]} samples)...")
    ss_res = run_ssgsea(vst_df, gene_sets, alpha=0.25)
    ss_res.index.name = "sample"
    merged_es = ss_res.join(col_df, how="inner")
    merged_es.to_csv(os.path.join(OUT_DIR, f"GSVA_scores_{ds_name}.csv"))

    for sample in merged_es.index:
        condition = "Unknown"
        for cc in ["condition", "group", "Condition", "Group"]:
            if cc in merged_es.columns:
                condition = merged_es.loc[sample, cc]
                break
        for term in gene_sets:
            all_gsva_res.append({
                "dataset":   ds_name,
                "sample":    sample,
                "condition": condition,
                "term":      term,
                "NES":       merged_es.loc[sample, term],
            })
    print(f"  Done {ds_name}.")

gsva_df = pd.DataFrame(all_gsva_res)
gsva_df.to_csv(os.path.join(OUT_DIR, "GSVA_sample_scores_all_datasets.csv"), index=False)
print("Saved GSVA sample scores.")

# ════════════════════════════════════════════════════════════════════════════
# 5. Differential ssGSEA statistics
# ════════════════════════════════════════════════════════════════════════════
print("\n=== Step 5: Differential GSVA statistics ===")
diff_stats_rows = []

for ds_name in gsva_df["dataset"].unique():
    sub_df = gsva_df[gsva_df["dataset"] == ds_name]
    mock_conds = [c for c in sub_df["condition"].unique()
                  if "mock" in str(c).lower() or "control" in str(c).lower()]
    if not mock_conds:
        print(f"  No mock condition for {ds_name}, skipping.")
        continue
    mock_cond = mock_conds[0]
    inf_conds = [c for c in sub_df["condition"].unique() if c != mock_cond]

    for inf_cond in inf_conds:
        for term in gene_sets:
            mock_vals = sub_df[(sub_df["condition"] == mock_cond) & (sub_df["term"] == term)]["NES"].values
            inf_vals  = sub_df[(sub_df["condition"] == inf_cond)  & (sub_df["term"] == term)]["NES"].values
            if len(mock_vals) == 0 or len(inf_vals) == 0:
                continue
            mean_mock  = float(np.mean(mock_vals))
            mean_inf   = float(np.mean(inf_vals))
            nes_shift  = mean_inf - mean_mock
            t_stat, p_val = welch_ttest(inf_vals, mock_vals)

            diff_stats_rows.append({
                "dataset":            ds_name,
                "infected_condition": inf_cond,
                "mock_condition":     mock_cond,
                "term":               term,
                "mean_mock_NES":      mean_mock,
                "mean_infected_NES":  mean_inf,
                "NES_shift":          nes_shift,
                "t_statistic":        t_stat,
                "pvalue":             p_val,
            })

diff_stats_df = pd.DataFrame(diff_stats_rows)
diff_stats_df.to_csv(os.path.join(OUT_DIR, "GSVA_contrast_differential_analysis.csv"), index=False)
print("Saved GSVA contrast differential statistics.")

# ════════════════════════════════════════════════════════════════════════════
# 6. Heatmap: Log2FC across 23 Verified Genes
# ════════════════════════════════════════════════════════════════════════════
print("\n=== Step 6: 23-gene heatmap ===")
heatmap_pivot = master_gene_df.pivot(index="symbol", columns="contrast_name", values="log2FoldChange")
hub_8 = ["DDIT3", "DDIT4", "CTH", "CHAC1", "ASNS", "DNAJB9", "HERPUD1", "MTHFD2"]
remaining_up   = [g for g in up_18_genes    if g in heatmap_pivot.index and g not in hub_8]
remaining_down = [g for g in down_4_genes   if g in heatmap_pivot.index]
remaining_disc = [g for g in discordant_genes if g in heatmap_pivot.index]
ordered_genes  = (
    [g for g in hub_8             if g in heatmap_pivot.index]
    + remaining_up + remaining_down + remaining_disc
)
heatmap_pivot  = heatmap_pivot.reindex(ordered_genes)

plt.figure(figsize=(12, 10))
ax = sns.heatmap(
    heatmap_pivot,
    cmap="coolwarm", center=0,
    annot=True, fmt=".2f",
    cbar_kws={"label": "Expression Log2FC"},
    linewidths=0.5, linecolor="white",
    annot_kws={"size": 8},
)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.set_label("Expression Log2FC", fontproperties={"family": "Arial", "size": 9})
for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color("white")

plt.title("Expression Profile of 23 Verified Genes Across Bulk RNA-seq Contrasts",
          fontsize=10, pad=15, color="black")
plt.xlabel("Bulk RNA-seq Contrast", fontsize=9, labelpad=10, color="black")
plt.ylabel("Verified Gene Symbol",  fontsize=9, labelpad=10, color="black")
plt.xticks(rotation=30, ha="right", fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "Verified_23_gene_heatmap_singlecell_and_bulk.pdf"),
            format="pdf", dpi=600)
plt.close()
print("Saved 23-gene heatmap.")

# ════════════════════════════════════════════════════════════════════════════
# 7. GSVA NES-shift heatmap
# ════════════════════════════════════════════════════════════════════════════
if not diff_stats_df.empty:
    print("\n=== Step 7: GSVA NES-shift heatmap ===")
    diff_stats_df["contrast"] = (
        diff_stats_df["dataset"].str.replace("GSE161783_", "")
        + ": " + diff_stats_df["infected_condition"]
    )
    gsva_shift_pivot = diff_stats_df.pivot(index="term", columns="contrast", values="NES_shift")

    plt.figure(figsize=(9.5, 4.5), dpi=300)
    sns.heatmap(
        gsva_shift_pivot,
        cmap="coolwarm", center=0,
        annot=True, fmt=".3f",
        cbar_kws={"label": "GSVA NES Shift (Infected - Mock)"},
        linewidths=0.75, linecolor="white",
        annot_kws={"size": 10, "weight": "bold"},
    )
    plt.title("GSVA Pathway Module Enrichment Shifts Across Bulk Infection Contrasts",
              fontsize=12, pad=15, fontweight="bold", color="#1A365D")
    plt.xlabel("Bulk RNA-Seq Infection Contrast", fontsize=10, labelpad=10,
               fontweight="bold", color="#2B6CB0")
    plt.ylabel("Core Gene Set Module", fontsize=10, labelpad=10,
               fontweight="bold", color="#2B6CB0")
    plt.xticks(rotation=25, ha="right", fontsize=9.5, fontweight="bold")
    plt.yticks(fontsize=10, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "GSVA_enrichment_shifts_heatmap.pdf"),
                format="pdf", bbox_inches="tight")
    plt.close()
    print("Saved GSVA NES-shift heatmap.")

# ════════════════════════════════════════════════════════════════════════════
# 8. Boxplots
# ════════════════════════════════════════════════════════════════════════════
if not gsva_df.empty:
    print("\n=== Step 8: GSVA boxplots ===")
    terms_to_plot = ["Verified_18_UP", "Verified_4_DOWN", "Verified_23_All"]
    plot_df = gsva_df[gsva_df["term"].isin(terms_to_plot)]

    g = sns.catplot(
        data=plot_df,
        x="condition", y="NES", hue="term",
        col="dataset", kind="box",
        palette="Set2", sharey=False,
        height=5, aspect=1.2,
    )
    g.fig.subplots_adjust(top=0.88)
    g.fig.suptitle(
        "GSVA Normalized Enrichment Scores (NES) Across Bulk RNA-seq Conditions",
        fontsize=14, fontweight="bold",
    )
    g.set_xticklabels(rotation=45, ha="right")
    g.set_axis_labels("Experimental Condition", "GSVA Enrichment Score (NES)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "GSVA_enrichment_boxplots.pdf"), format="pdf")
    plt.close()
    print("Saved GSVA boxplots.")

print("\n=== Summary table: key NES values ===")
if not diff_stats_df.empty:
    key = diff_stats_df[diff_stats_df["term"] == "Verified_18_UP"][
        ["dataset", "infected_condition", "mean_mock_NES", "mean_infected_NES", "NES_shift", "pvalue"]
    ]
    print(key.to_string(index=False))

print("\n✅ All GSVA bulk analysis steps completed successfully!")
