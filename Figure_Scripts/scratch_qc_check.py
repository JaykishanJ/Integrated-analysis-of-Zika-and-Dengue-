import pandas as pd
import numpy as np

# Load the Step 3 distribution (all cells, has n_virus_molecules)
df_dist = pd.read_csv("E:/Zika/results/phaseB/tables/Cell_Distributions.csv")

# Load the Step 4 filter calls (has the 'keep' flag)
df_filter = pd.read_csv("E:/Zika/results/phaseB/tables/Cell_Filter_Calls.csv")

# Join on Sample
df = pd.merge(df_dist, df_filter[['Sample', 'keep']], on='Sample', how='left')

df['n_virus_molecules'] = df['n_virus_molecules'].fillna(0)

# Calculate 99th percentile of MOI 0 controls per virus
bg_thr = df[df['moi'] == 0].groupby('virus')['n_virus_molecules'].quantile(0.99).reset_index()
bg_thr.rename(columns={'n_virus_molecules': 'hop_threshold'}, inplace=True)
bg_thr['hop_threshold'] = np.maximum(bg_thr['hop_threshold'], 1)

# Join threshold back
df = pd.merge(df, bg_thr, on='virus', how='left')

# Classify infection status
def classify(row):
    if pd.isna(row['virus']) or row['virus'] == "":
        return "Unknown"
    if row['moi'] == 0:
        return "Unexposed"
    if row['n_virus_molecules'] > row['hop_threshold']:
        return "Infected"
    return "Bystander"

df['infection_status'] = df.apply(classify, axis=1)

# Calculate pass rates
summary = df.groupby('infection_status').agg(
    Total_Cells=('Sample', 'count'),
    Passed_QC=('keep', 'sum')
).reset_index()

summary['Failed_QC'] = summary['Total_Cells'] - summary['Passed_QC']
summary['Pass_Rate_Pct'] = (100 * summary['Passed_QC'] / summary['Total_Cells']).round(2)

print(summary.to_string(index=False))

# Breakdown of the failed cells by infection status
failed_breakdown = summary[['infection_status', 'Failed_QC']].copy()
failed_breakdown['Pct_of_Failed'] = (100 * failed_breakdown['Failed_QC'] / failed_breakdown['Failed_QC'].sum()).round(2)
print("\nBreakdown of failed cells:")
print(failed_breakdown.to_string(index=False))
