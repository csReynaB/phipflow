#!/usr/bin/env bash
#SBATCH --job-name=phipflow
#SBATCH --output=logs/phipflow_%j.out
#SBATCH --error=logs/phipflow_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

config=$1
run_name=$2

phipflow_dir="/lisc/data/scratch/ccr/CR_projects/phiper/phipflow"
run_dir="${phipflow_dir}/runs/${run_name}"

mkdir -p "${run_dir}/logs"
cd "${run_dir}"

module load Nextflow

nextflow run "${phipflow_dir}/main.nf" \
  -profile lisc,apptainer \
  -c "${config}" \
  -resume
