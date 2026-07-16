# phipflow 

Version 0.9.0

`phipflow` is a Nextflow DSL2 wrapper for running the PHIPER analysis workflow on a server.

It automates the full PHIPER workflow:

1. Create a PHIPER parquet object from `exist.csv`, `fold.csv`, and metadata.
2. Run PHIPER analyses for one or more group columns.
3. Render automated Quarto HTML reports.

The workflow is designed to keep reusable workflow code separate from project data and results.

---

## Recommended use

For production runs, use the **LiSC + Apptainer** profile:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

This runs the workflow on SLURM and uses the reproducible Apptainer/SIF image containing R, Quarto, Pandoc, `phiper`, `phiperio`, and the required R packages.

The fallback `lisc` profile can also be used, but then the `base_dir` must already contain a working `renv` environment with all required R packages installed. This is currently true for:

```text
/lisc/data/scratch/ccr/CR_projects/phiper
```

---

## Repository layout

A typical setup on LiSC is:

```text
/lisc/data/scratch/ccr/CR_projects/
├── phiper/
│   ├── renv/
│   ├── IBD-Berlin/
│   │   ├── Data/
│   │   │   ├── exist.csv
│   │   │   ├── fold.csv
│   │   │   └── IBD-Berlin.parquet
│   │   ├── Metadata/
│   │   │   └── IBD-Berlin_metadata.csv
│   │   ├── R/
│   │   │   └── group_config.R
│   │   └── results/
│   └── other-projects/
│
└── phipflow/
    ├── conf/
    ├── containers/
    │   └── phipflow-r4.6.0-quarto1.9.37.sif
    ├── docker/
    │   └── install_r_packages.R
    ├── modules/
    │   └── local/
    ├── src/
    │   ├── R/
    │   │   ├── 01-create_phiper_object.R
    │   │   ├── 02-run_phiper_analysis.R
    │   │   ├── 03-render_phiper_reports.R
    │   │   └── helper_functions.R
    │   └── template/
    │       └── phiper_summary_report.qmd
    ├── main.nf
    ├── nextflow.config
    ├── Dockerfile
    ├── README.md
    ├── run_example.sh
    ├── run_test.sh
    └── submit_run_test.sh
```

The important distinction is:

```text
base_dir   = /lisc/data/scratch/ccr/CR_projects/phiper
projectDir = the phipflow repository directory
```

`base_dir` contains project data, metadata, project-specific `group_config.R`, results, and optionally the LiSC `renv` environment used by the non-container profile.

`projectDir` contains the reusable Nextflow workflow, R scripts, Quarto template, Dockerfile, and Apptainer/SIF image location.

---

## Required project structure

Each PHIPER project should have this layout:

```text
<base_dir>/<project_name>/
├── Data/
│   ├── exist.csv
│   └── fold.csv
├── Metadata/
│   └── <metadata_file>.csv
├── R/
│   └── group_config.R
└── results/
```

For example:

```text
/lisc/data/scratch/ccr/CR_projects/phiper/IBD-Berlin/
├── Data/
│   ├── exist.csv
│   └── fold.csv
├── Metadata/
│   └── IBD-Berlin_metadata.csv
├── R/
│   └── group_config.R
└── results/
```

The file `<project_name>/R/group_config.R` defines the PHIPER group configurations used by the analysis.

---

## Input data format

`phipflow` expects PHIPER input matrices and metadata in a consistent format.

### `exist.csv`

`exist.csv` should be a **wide peptide-by-sample matrix**:

```text
peptide_name,sample_1,sample_2,sample_3,...
pep_0001,0,1,0,...
pep_0002,1,1,0,...
pep_0003,0,0,1,...
```

Rows are peptides/features. Columns are samples.

By default, the peptide identifier column is expected to be:

```text
peptide_name
```

This can be changed with:

```bash
--peptide_col your_peptide_column
```

The sample columns are detected from the first column name starting with the configured sample prefix. By default:

```bash
--sample_prefix R
```

### `fold.csv`

`fold.csv` should have the same wide structure as `exist.csv`:

```text
peptide_name,sample_1,sample_2,sample_3,...
pep_0001,0,2.3,0,...
pep_0002,1.8,4.1,0,...
pep_0003,0,0,3.5,...
```

It must contain the same peptide identifier column and matching sample columns.

`exist.csv` is used as the binary enrichment matrix, where:

```text
1 = enriched / present
0 = not enriched / absent
```

`fold.csv` provides the corresponding fold-change values. Infinite values can be replaced by the maximum finite fold-change value using:

```bash
--replace_inf TRUE
```

This is the default behavior.

### Metadata file

The metadata file should be in **long/sample-level format**:

```text
SampleName,Sex,Age,group_test,group_IBD,subject_id,Smoker,Timepoint,Smoker_Timepoint
R001,F,34,Control,Control,S01,Smoker,BL,Smoker_BL
R002,F,34,Disease,IBD,S01,Smoker,FU,Smoker_FU
R003,M,41,Control,Control,S02,NonSmoker,BL,NonSmoker_BL
```

The first column is assumed to contain the sample name / sample ID. During object creation, the first metadata column is renamed internally to:

```text
sample_id
```

The sample IDs in the metadata must match the sample column names in `exist.csv` and `fold.csv`.

### Longitudinal data

For longitudinal or paired analyses, the metadata should include a subject-level identifier column named exactly:

```text
subject_id
```

This column is used to match repeated samples from the same individual.

For now, longitudinal group columns should be created by concatenating the biological group column and the timepoint column:

```text
<group_col>_<timepoint_col>
```

For example, if the metadata has:

```text
Smoker
Timepoint
```

where `Timepoint` contains values such as:

```text
BL, M3, FU
```

then the longitudinal analysis column should be named:

```text
Smoker_Timepoint
```

and should contain combined labels such as:

```text
Smoker_BL
Smoker_M3
Smoker_FU
NonSmoker_BL
NonSmoker_FU
```

This combined column is the value that should be used in `group_config.R` and passed to `--group_cols` when running a longitudinal timepoint analysis.

Example:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols Smoker_Timepoint \ # this can be named different than the group_col "Smoker_Timepoint" in metadata but for consistency it's better to name it the same
  -resume
```

In `group_config.R`, the same name should be present in `group_definitions`, for example:

```r
group_definitions <- list(
  Smoker_Timepoint = list(
    group_col = "Smoker_Timepoint",
    groups = c("Smoker_BL", "Smoker_FU", "NonSmoker_BL", "NonSmoker_FU"),
    comparisons = list(
      c("Smoker_BL", "Smoker_FU"),
      c("NonSmoker_BL", "NonSmoker_FU")
    ),
    longitudinal = c(TRUE, TRUE)
  )
)
```

---

## Runtime profiles

`phipflow` supports two main runtime modes.

### 1. Recommended: LiSC + Apptainer

Use this for reproducible production runs:

```bash
nextflow run main.nf -profile lisc,apptainer,test -resume
```

or for a real project:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

This mode uses:

```text
containers/phipflow-r4.6.0-quarto1.9.37.sif
```

The image contains:

```text
R 4.6.0
Quarto 1.9.37
Pandoc
phiper
phiperio
dplyr
duckdb
ggplot2
ggpubr
ggtext
plotly
vegan
openxlsx
DT
and other workflow dependencies
```

In this mode, the workflow does **not** load LiSC R/Quarto/Pandoc modules inside the SLURM jobs. The R environment comes from the Apptainer image.

### 2. Fallback: LiSC modules + base_dir renv

Use this only if you want to run with LiSC modules and the `renv` environment under `base_dir`:

```bash
nextflow run main.nf -profile lisc,test -resume
```

or:

```bash
nextflow run main.nf \
  -profile lisc \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

This mode assumes that:

```text
<base_dir>/renv/
```

exists and already contains all required R packages. For the current LiSC setup, this is true for:

```text
/lisc/data/scratch/ccr/CR_projects/phiper
```

If you use another `--base_dir` with the `lisc` profile only, make sure that directory has a working `renv` environment and that R/Quarto/Pandoc modules are available.

---

## First-time setup: build the Apptainer image on LiSC

The large `.sif` file is not directly given in this GitHub repo. Instead, it can be built on LiSC from the Docker image hosted on GitHub Container Registry.

The Docker image is:

```text
ghcr.io/csreynab/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

Build the Apptainer/SIF image:

```bash
cd /lisc/data/scratch/ccr/CR_projects/phipflow/containers


apptainer build phipflow-r4.6.0-quarto1.9.37.sif \
  docker://ghcr.io/csreynab/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

If the file already exists and you want to rebuild it:

```bash
apptainer build --force phipflow-r4.6.0-quarto1.9.37.sif \
  docker://ghcr.io/csreynab/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

Test the image:

```bash
apptainer exec phipflow-r4.6.0-quarto1.9.37.sif R --version

apptainer exec phipflow-r4.6.0-quarto1.9.37.sif quarto --version

apptainer exec phipflow-r4.6.0-quarto1.9.37.sif \
  Rscript -e "library(phiper); library(phiperio); library(duckdb); library(DT); library(quarto); cat('OK\\n')"
```

The `.sif` file is expected by default at:

```text
<phipflow repo>/containers/phipflow-r4.6.0-quarto1.9.37.sif
```

---

## Quick start on LiSC

Go to the workflow directory:

```bash
cd /lisc/data/scratch/ccr/CR_projects/phipflow
```

Load Nextflow:

```bash
module load Nextflow
```

Run the default test case:

```bash
nextflow run main.nf -profile lisc,apptainer,test -resume
```

The `test` profile uses:

```groovy
params.project_name  = 'IBD-Berlin'
params.metadata_file = 'IBD-Berlin_metadata.csv'
params.group_cols    = 'group_test'
```

If this test succeeds, the workflow, SLURM profile, Apptainer image, and project structure are working.

---

## Running your own project

Example using the recommended Apptainer profile:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

Example using LiSC modules and the `base_dir` renv environment:

```bash
nextflow run main.nf \
  -profile lisc \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

When using `-profile lisc` without `apptainer`, make sure `base_dir` has a valid `renv` environment with all required packages.

---

## Running several group columns

Pass comma-separated group columns:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test,group_combined,group_3 \
  -resume
```

Each group column is passed to `02-run_phiper_analysis.R` as:

```text
ACTIVE_GROUP=<group_col>
```

For example:

```text
ACTIVE_GROUP=group_test
ACTIVE_GROUP=group_combined
ACTIVE_GROUP=group_3
```

The analyses are submitted as independent SLURM jobs and can run in parallel.

In a config/profile, use one comma-separated string:

```groovy
params.group_cols = 'group_test,group_combined,group_3'
```

not multiple quoted strings.

---

## Report rendering

By default, reports are rendered for the same groups listed in `--group_cols`.

For example:

```bash
--group_cols group_test,group_combined,group_3
```

will render reports for:

```text
group_test
group_combined
group_3
```

If you want to analyze several groups but render reports only for a subset, use `--report_group_cols`:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test,group_combined,group_3 \
  --report_group_cols group_test,group_combined \
  -resume
```

Reports are written to:

```text
<base_dir>/<project_name>/results/<group_col>/summary_report_<group_col>.html
```

For example:

```text
/lisc/data/scratch/ccr/CR_projects/phiper/IBD-Berlin/results/group_test/summary_report_group_test.html
```

---

## Automatic group discovery

If `--group_cols` is not provided, `phipflow` tries to discover all available groups from:

```text
<base_dir>/<project_name>/R/group_config.R
```

specifically from the object:

```r
group_definitions <- list(...)
```

Example:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  -resume
```

This will run all group definitions found in `group_config.R`.

---

## Main parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--base_dir` | no | `/lisc/data/scratch/ccr/CR_projects/phiper` | Root directory containing PHIPER projects. |
| `--project_name` | yes | `null` | Project folder name inside `base_dir`. |
| `--metadata_file` | yes | `null` | Metadata file inside `<project_name>/Metadata/`. |
| `--group_cols` | no | `null` | Comma-separated group columns to analyze. If null, groups are discovered from `group_config.R`. |
| `--report_group_cols` | no | `null` | Comma-separated group columns to render. If null, uses `group_cols`. |
| `--exist_file` | no | `exist.csv` | Enrichment matrix inside `<project_name>/Data/`. |
| `--fold_file` | no | `fold.csv` | Fold-change matrix inside `<project_name>/Data/`. |
| `--out_parquet` | no | `<project_name>.parquet` | Output parquet file created in `<project_name>/Data/`. |
| `--all` | no | `FALSE` | Passed to `02-run_phiper_analysis.R`. |
| `--default_longitudinal` | no | `FALSE` | Passed to `02-run_phiper_analysis.R`. |
| `--manual_comparison_file` | no | `NULL` | Optional manual comparison file inside `<project_name>/R/`. |
| `--container` | no | set by `apptainer` profile | Path to the local `.sif` image. |

---

## Workflow steps

### 1. Create PHIPER object

Process:

```text
CREATE_PHIPER_OBJECT
```

Runs:

```text
src/R/01-create_phiper_object.R
```

Creates:

```text
<base_dir>/<project_name>/Data/<out_parquet>
```

---

### 2. Run PHIPER analysis

Process:

```text
RUN_PHIPER_ANALYSIS
```

Runs:

```text
src/R/02-run_phiper_analysis.R
```

One job is launched per group column.

Outputs are written to:

```text
<base_dir>/<project_name>/results/<group_col>/
```

---

### 3. Render reports

Process:

```text
RENDER_PHIPER_REPORTS
```

Runs:

```text
src/R/03-render_phiper_reports.R
```

Uses template:

```text
src/template/phiper_summary_report.qmd
```

Outputs:

```text
<base_dir>/<project_name>/results/<group_col>/summary_report_<group_col>.html
```

---

## SLURM resources

Default resources are defined in `nextflow.config`:

```groovy
withLabel: create_object {
    cpus = 1
    memory = '25 GB'
    time = '30 min'
}

withLabel: run_analysis {
    cpus = 25
    memory = '50 GB'
    time = '8 h'
}

withLabel: render_reports {
    cpus = 1
    memory = '10 GB'
    time = '30 min'
}
```

The number of CPUs is passed to the R analysis as:

```text
N_CORES=<task.cpus>
```

Memory is converted to:

```text
MAX_GB
```

and passed to the R script.

---

## Running workflows

For short tests, it is fine to run Nextflow interactively from a `tmux` session or `screen` in your local computer.  

For longer production runs, it is recommended to submit the Nextflow controller itself as a small SLURM job.

Nextflow will then submit the heavy workflow steps as separate SLURM jobs.

### Option 1: run inside `tmux`

For test/short runs, launch Nextflow inside `tmux`:

```bash
tmux new -s phipflow

cd /lisc/data/scratch/ccr/CR_projects/phipflow

module load Nextflow

nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

Detach:

```text
Ctrl+b, then d
```

Reconnect:

```bash
tmux attach -t phipflow
```

Nextflow submits the workflow steps as SLURM jobs. The Nextflow command itself acts as the controller and should remain running until the workflow finishes.

### Option 2: submit the Nextflow controller as a SLURM job

```bash
#!/usr/bin/env bash
#SBATCH --job-name=phipflow
#SBATCH --output=logs/phipflow_%j.out
#SBATCH --error=logs/phipflow_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

cd /lisc/data/scratch/ccr/CR_projects/phipflow 

mkdir -p logs

module load Nextflow

nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test,group_combined \
  -resume
```

Submit it with:

```text
sbatch submit_run.sh
```

---

## Checking SLURM jobs

```bash
squeue -u $USER
```

---

## Resuming runs

Use `-resume` to continue from completed steps:

```bash
nextflow run main.nf \
  -profile lisc,apptainer \
  --project_name IBD-Berlin \
  --metadata_file IBD-Berlin_metadata.csv \
  --group_cols group_test \
  -resume
```

If the parquet object already exists and the corresponding process hash has not changed, Nextflow can reuse the cached result.

---

## Rebuilding and publishing the Docker image

The Apptainer image is built from this GHCR Docker image:

```text
ghcr.io/csreynab/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

To rebuild locally, a Dockerfile and docker/install_r_packages.R are provides so you can run:

```bash
docker build -t phipflow-r4.6.0-quarto1.9.37 . # or whathever label you want
```

Tag for GitHub Container Registry:

```bash
docker tag phipflow-r4.6.0-quarto1.9.37 \
  ghcr.io/user/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

Push:

```bash
docker push ghcr.io/user/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

Then rebuild the `.sif` on LiSC:

```bash
cd /lisc/data/scratch/ccr/CR_projects/phipflow/containers

apptainer build --force phipflow-r4.6.0-quarto1.9.37.sif \
  docker://ghcr.io/csreynab/phipflow-r4.6.0-quarto1.9.37:0.1.0
```

Note: overwriting the same `0.1.0` tag is convenient, but less reproducible than using versioned tags such as `0.1.1`, `0.1.2`, etc.

---

## Notes

- `phipflow` expects project-specific settings in `<project_name>/R/group_config.R`.
- Reusable R scripts and the Quarto template live in `phipflow/src/`.
- Project data and results live outside `phipflow`, under `base_dir`.
- Use `-profile lisc,apptainer` for reproducible production runs.
- Use `-profile lisc` only when the `base_dir` `renv` environment is correctly set up.
- Use `-profile lisc,apptainer,test` for the default container-based test case.
