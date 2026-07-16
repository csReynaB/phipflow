process RENDER_PHIPER_REPORTS {
    tag "${project_name}:${group_cols}"

    label 'render_reports'

    input:
    path analysis_markers
    val base_dir
    val project_name
    val group_cols
    val output_group_mode
    val workflow_src_dir
    val workflow_template
    val use_modules

    output:
    path "render_phiper_reports.${project_name}.done", emit: report_marker

    script:
    """
    set -euo pipefail

    workdir=\$PWD

    if [[ "${use_modules}" == "true" ]]; then
      module purge
      module load R
      module load Quarto
      module load Pandoc/3.9 || true
    fi

    echo "Runtime check"
    echo "R path    : \$(which R || true)"
    echo "Rscript   : \$(which Rscript || true)"
    echo "Quarto    : \$(which quarto || true)"
    R --version | head -n 1 || true
    quarto --version || true

    cd "${base_dir}"

    echo "Rendering PHIPER reports"
    echo "Base dir          : ${base_dir}/${project_name}/results"
    echo "Active groups     : ${group_cols}"
    echo "Output group mode : ${output_group_mode}"
    echo "Workflow src      : ${workflow_src_dir}"
    echo "Template          : ${workflow_template}"
    echo "Use modules       : ${use_modules}"

    Rscript --vanilla "${workflow_src_dir}/03-render_phiper_reports.R" \\
      BASE_DIR="${base_dir}/${project_name}/results" \\
      GROUP_COLS="${group_cols}" \\
      OUTPUT_GROUP_MODE="${output_group_mode}" \\
      PHIPFLOW_SRC="${workflow_src_dir}" \\
      TEMPLATE="${workflow_template}"

    echo "PHIPER reports finished for project: ${project_name}"
    echo "DONE: ${project_name} reports" > "\${workdir}/render_phiper_reports.${project_name}.done"
    """
}