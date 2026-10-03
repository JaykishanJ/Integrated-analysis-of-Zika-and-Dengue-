import os
import shutil

source_dir = "e:/Zika"
dest_dir = "e:/Zika/Figure_Scripts_Package"

scripts = [
    "publication_figures_scripts/Figure1_BulkRNAseq.R",
    "publication_figures_scripts/Figure2_SingleCellOverview.R",
    "publication_figures_scripts/Figure3_CellSpecificResponses.R",
    "publication_figures_scripts/Figure4_Trajectories.R",
    "publication_figures_scripts/Figure5_CellCommunication.R",
    "publication_figures_scripts/Figure6_Integration.R",
    "publication_figures_scripts/utils.R"
]

data_files = [
    "results/phaseB/objects/Phase_B_cell_qc.rds",
    "results/phaseC/objects/Phase_C_log_matrix.rds",
    "results/phaseD/objects/Phase_D_clustering.rds",
    "results/phaseD/objects/Phase_D_ZIKV_Sig_Matrix.rds",
    "results/phaseD/tables/Step11d_All_Timepoint_DEGs.csv",
    "results/phaseE/tables/Step17_TopClusterMarkers.csv",
    "results/phaseD/cache/trajectory/monocle3_trajectory.rds",
    "results/phaseD/cache/trajectory/gene_trends.rds",
    "results/phaseD/cache/cellchat/cellchat_object.rds",
    "results/phaseD/cache/cellchat/significant_LR_pairs.rds",
    "ZIKA_Bulk_Som/ZIKA_Bulk_Som/results/GSE80434/step6_deseq/dds_ZIKV_fitted.rds",
    "ZIKA_Bulk_Som/ZIKA_Bulk_Som/results/GSE80434/step8_degs/ZIKVM_vs_Mock_DEGs_all.csv",
    "ZIKA_Bulk_Som/ZIKA_Bulk_Som/results/GSE110512/step6_deseq/dds_fitted.rds",
    "ZIKA_Bulk_Som/ZIKA_Bulk_Som/results/GSE110512/step8_degs/DEGs_all.csv",
    "ZIKA_Bulk_Som/ZIKA_Bulk_Som/results/GSE161783/step6_deseq/dds_Huh7_5_fitted.rds",
    "ZIKA_Bulk_Som/ZIKA_Bulk_Som/results/GSE161783/step8_degs/Huh7_5__Zika_vs_Mock_DEGs_all.csv",
    "output_manuscript_figures/GSVA_sample_scores_all_datasets.csv",
    "output_manuscript_figures/GSVA_contrast_differential_analysis.csv"
]

os.makedirs(dest_dir, exist_ok=True)

for sf in scripts:
    s_path = os.path.join(source_dir, sf)
    d_path = os.path.join(dest_dir, os.path.basename(sf))
    if os.path.exists(s_path):
        with open(s_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        if "utils.R" not in sf:
            # Replace BASE_DIR
            content = content.replace('BASE_DIR <- "e:/Zika"', 'BASE_DIR <- "."')
            # Replace absolute util source
            content = content.replace('source("e:/Zika/publication_figures_scripts/utils.R")', 'source("utils.R")')
            
        with open(d_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Copied script {sf}")
    else:
        print(f"Missing script {sf}")

for df in data_files:
    s_path = os.path.join(source_dir, df)
    d_path = os.path.join(dest_dir, df)
    os.makedirs(os.path.dirname(d_path), exist_ok=True)
    if os.path.exists(s_path):
        shutil.copy2(s_path, d_path)
        print(f"Copied data {df}")
    else:
        print(f"Missing data {df}")

print("Done grouping files.")
