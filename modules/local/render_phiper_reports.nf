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
    val results_name
    val use_modules
    val delta_min_m_eff
    val rank_cols

    output:
    path "render_phiper_reports.${project_name}.done", emit: report_marker

    script:
    """
    set -euo pipefail

    workdir=\$PWD
    results_dir="${base_dir}/${project_name}/${results_name}"

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
    echo "Results dir   : \${results_dir}"
    echo "Active groups     : ${group_cols}"
    echo "Output group mode : ${output_group_mode}"
    echo "Workflow src      : ${workflow_src_dir}"
    echo "Template          : ${workflow_template}"
    echo "Use modules       : ${use_modules}"
    echo "Delta minimum effect size : ${delta_min_m_eff}"
    echo "Rank columns used : ${rank_cols}"

    Rscript --vanilla "${workflow_src_dir}/03-render_phiper_reports.R" \\
      BASE_DIR="\${results_dir}" \\
      GROUP_COLS="${group_cols}" \\
      OUTPUT_GROUP_MODE="${output_group_mode}" \\
      PHIPFLOW_SRC="${workflow_src_dir}" \\
      TEMPLATE="${workflow_template}" \\
      DELTA_MIN_M_EFF="${delta_min_m_eff}" \\
      RANK_COLS="${rank_cols}"

    echo "PHIPER reports finished for project: ${project_name}"
    echo "DONE: \${results_dir} reports" > "\${workdir}/render_phiper_reports.${project_name}.done"
    """
}