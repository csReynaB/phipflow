#!/usr/bin/env bash
set -euo pipefail

module purge
module load Nextflow
# Nextflow-native test run using the built-in test profile.
# The profile sets:
#   project_name  = IBD-Berlin
#   metadata_file = IBD-Berlin_metadata.csv
#   group_cols    = group_test

nextflow run main.nf -profile lisc,apptainer,test -resume
