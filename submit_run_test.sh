#!/usr/bin/env bash
#SBATCH --job-name=phipflow
#SBATCH --output=logs/phipflow_%j.out
#SBATCH --error=logs/phipflow_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB

set -euo pipefail

cd /lisc/data/scratch/ccr/CR_projects/phiper/phipflow

module load Nextflow

nextflow run main.nf \
  -profile lisc,apptainer,test \
  -resume
