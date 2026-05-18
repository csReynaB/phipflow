process CREATE_PHIPER_OBJECT {
    tag "${project_name}"

    label 'create_object'

    input:
    val base_dir
    val project_name
    val exist_file
    val fold_file
    val metadata_file
    val out_parquet
    val sample_prefix
    val peptide_col
    val replace_inf
    val workflow_src_dir
    val use_modules

    output:
    path "create_phiper_object.${project_name}.done", emit: parquet_marker

    script:
    """
    set -euo pipefail

    workdir=\$PWD

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

    echo "Running PHIPER object creation"
    echo "Base dir       : ${base_dir}"
    echo "Project        : ${project_name}"
    echo "Exist file     : ${exist_file}"
    echo "Fold file      : ${fold_file}"
    echo "Metadata file  : ${metadata_file}"
    echo "Output parquet : ${out_parquet}"
    echo "Workflow src   : ${workflow_src_dir}"
    echo "Use modules    : ${use_modules}"

    Rscript --vanilla "${workflow_src_dir}/01-create_phiper_object.R" \\
      --base_dir "${base_dir}" \\
      --project_name "${project_name}" \\
      --exist_file "${exist_file}" \\
      --fold_file "${fold_file}" \\
      --metadata_file "${metadata_file}" \\
      --out_parquet "${out_parquet}" \\
      --sample_prefix "${sample_prefix}" \\
      --peptide_col "${peptide_col}" \\
      --replace_inf "${replace_inf}"

    test -s "${project_name}/Data/${out_parquet}"

    echo "${base_dir}/${project_name}/Data/${out_parquet}" > "\${workdir}/create_phiper_object.${project_name}.done"
    """
}