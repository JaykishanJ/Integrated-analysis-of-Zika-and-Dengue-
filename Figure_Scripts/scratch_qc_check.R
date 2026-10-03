# R script to calculate QC pass rate by infection status
library(dplyr)
library(tidyr)

# Load Phase A QC data (pre-filter)
qc_a <- readRDS("E:/Zika/results/phaseA/objects/Phase_A_cell_qc.rds")

# Load Phase B QC filter calls
qc_b <- read.csv("E:/Zika/results/phaseB/tables/Step4_Cell_Filter_Calls.csv")

# Join to get 'keep' status for all Phase A cells
qc_all <- qc_a %>% left_join(qc_b %>% select(Sample, keep), by="Sample")
qc_all$keep[is.na(qc_all$keep)] <- FALSE # Just in case

# Make sure n_virus_molecules is numeric and not NA
qc_all$n_virus_molecules[is.na(qc_all$n_virus_molecules)] <- 0

# Calculate background thresholds (from Phase C logic)
bg_thr <- qc_all %>%
  filter(moi == 0) %>%
  group_by(virus) %>%
  summarise(
    hop_threshold = pmax(quantile(n_virus_molecules, 0.99, na.rm = TRUE), 1),
    .groups = "drop"
  )

# Assign infection status to ALL cells
qc_all <- qc_all %>%
  left_join(bg_thr, by = "virus") %>%
  mutate(
    infection_status = case_when(
      is.na(virus) | virus == "" ~ "Unknown",
      moi == 0 ~ "Unexposed",
      n_virus_molecules > hop_threshold ~ "Infected",
      TRUE ~ "Bystander"
    )
  )

# Calculate Pass Rate by Infection Status
summary_stats <- qc_all %>%
  group_by(infection_status) %>%
  summarise(
    Total_Cells = n(),
    Passed_QC = sum(keep),
    Failed_QC = Total_Cells - Passed_QC,
    Pass_Rate_Pct = round(100 * Passed_QC / Total_Cells, 2)
  )

print(summary_stats)
