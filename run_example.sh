#!/usr/bin/env bash
set -euo pipefail

module purge
module load Nextflow
module load Apptainer

nextflow run main.nf \
  -profile lisc,apptainer \
  -c conf/params.example.config \
  -resume
