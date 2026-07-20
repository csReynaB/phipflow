#!/usr/bin/env bash
set -euo pipefail

module purge
module load Nextflow

# nextflow run main.nf \
#   -profile lisc,apptainer \
#   -c conf/params.example.config \
#   -resume

nextflow run main.nf -profile lisc,apptainer   --base_dir /lisc/data/scratch/ccr/GI_projects/analysis/phiper   --project_name PCa_Innsbruck   --metadata_file PC
a_Innsbruck_metadata.csv   --group_cols "Prostate_cancer"   --report_group_cols "Prostate_cancer"    --rank_cols "genus,species,protein_seq_id"   --aggregate_stat maxmean   --delta_min_m_eff 20  --results_name result
s_pept7-93prev_maxmean_test   --force
