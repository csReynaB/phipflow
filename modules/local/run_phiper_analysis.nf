process RUN_PHIPER_ANALYSIS {
    tag "${project_name}:${active_group}"

    label 'run_analysis'

    input:
    tuple path(parquet_marker), val(active_group)
    val base_dir
    val project_name
    val parquet_name
    val all_flag
    val default_longitudinal
    val manual_comparison_file
    val force
    val output_group_mode
    val workflow_src_dir
    val peptide_library
    val rank_cols
    val aggregate_stat
    val delta_min_m_eff
    val results_name
    val use_modules
  
    output:
    path "run_phiper_analysis.${project_name}.${active_group}.done", emit: analysis_marker

    script:
    """
    set -euo pipefail

    workdir=\$PWD
    results_dir="${base_dir}/${project_name}/${results_name}"

    if [[ "${use_modules}" == "true" ]]; then
      module purge
      module load R
      module load Pandoc/3.9 || true
    fi

    echo "Runtime check"
    echo "R path    : \$(which R || true)"
    echo "Rscript   : \$(which Rscript || true)"
    R --version | head -n 1 || true

    cd "${base_dir}"

    if [[ -n "\${SLURM_MEM_PER_NODE:-}" ]]; then
      max_gb=\$(( SLURM_MEM_PER_NODE / 1024 - 1 ))
      if [[ "\${max_gb}" -lt 1 ]]; then max_gb=1; fi
    else
      max_gb=${task.memory.toGiga()}
    fi

    echo "Running PHIPER analysis"
    echo "Project                : ${project_name}"
    echo "Active group           : ${active_group}"
    echo "Output group mode      : ${output_group_mode}"
    echo "Parquet                : ${parquet_name}"
    echo "Peptide library        : ${peptide_library}"
    echo "ALL                    : ${all_flag}"
    echo "DEFAULT_LONGITUDINAL   : ${default_longitudinal}"
    echo "MANUAL_COMPARISON_FILE : ${manual_comparison_file}"
    echo "Workflow src           : ${workflow_src_dir}"
    echo "Use modules            : ${use_modules}"
    echo "N_CORES                : ${task.cpus}"
    echo "MAX_GB                 : \${max_gb}"

    Rscript --vanilla "${workflow_src_dir}/02-run_phiper_analysis.R" \\
      N_CORES="${task.cpus}" \\
      MAX_GB="\${max_gb}" \\
      LOG=true \\
      FORCE="${force}" \\
      ACTIVE_GROUP="${active_group}" \\
      OUTPUT_GROUP_MODE="${output_group_mode}" \\
      ALL="${all_flag}" \\
      DEFAULT_LONGITUDINAL="${default_longitudinal}" \\
      PROJECT_DIR="${project_name}" \\
      PARQUET_NAME="${parquet_name}" \\
      MANUAL_COMPARISON_FILE="${manual_comparison_file}" \\
      PHIPFLOW_SRC="${workflow_src_dir}" \\
      PEPTIDE_LIBRARY="${peptide_library}" \\
      RANK_COLS="${rank_cols}" \\
      AGGREGATE_STAT="${aggregate_stat}" \\
      DELTA_MIN_M_EFF="${delta_min_m_eff}" \\
      RESULTS_NAME="${results_name}"

    echo "Checking result directory: \${results_dir}"
    echo "PHIPER analysis finished for active group: ${active_group}"
    echo "DONE: \${results_dir} ${active_group}" > "\${workdir}/run_phiper_analysis.${project_name}.${active_group}.done"
    """
}