nextflow.enable.dsl = 2

include { CREATE_PHIPER_OBJECT } from './modules/local/create_phiper_object'
include { RUN_PHIPER_ANALYSIS } from './modules/local/run_phiper_analysis'
include { RENDER_PHIPER_REPORTS } from './modules/local/render_phiper_reports'

/*
 * PHIPER end-to-end workflow
 *
 * projectDir:
 *   The phipflow pipeline root.
 *   Example:
 *     /lisc/data/scratch/ccr/CR_projects/phipflow
 *
 * params.base_dir:
 *   The PHIPER project/data root.
 *   Example:
 *     /lisc/data/scratch/ccr/CR_projects/phiper
 *
 * Meaning of group parameters:
 *   params.group_cols =
 *     group columns to run through 02-run_phiper_analysis.R as ACTIVE_GROUP.
 *     Example:
 *       group_test,group_combined,group_subtype
 *
 *     If null/empty, all entries from:
 *       <base_dir>/<project_name>/R/group_config.R::group_definitions
 *     are used.
 *
 *   params.report_group_cols =
 *     optional override for 03-render_phiper_reports.R as GROUP_COLS.
 *     If not provided, defaults to the final group_cols list.
 */

def parseCsvParam(value) {
    if (value == null) return []

    def s = value.toString().trim()

    if (!s || s.equalsIgnoreCase('null')) return []

    return s
        .split(',')
        .collect { it.trim() }
        .findAll { it }
}

def discoverGroupsFromConfig(baseDir, projectName) {
    def cfg = file("${baseDir}/${projectName}/R/group_config.R")

    if (!cfg.exists()) {
        error "--group_cols was not provided and group_config.R was not found: ${cfg}"
    }

    def text = cfg.text
    def start = text.indexOf('group_definitions')

    if (start < 0) {
        error "--group_cols was not provided, but no group_definitions object was found in ${cfg}"
    }

    def block = text.substring(start)
    def groups = []

    def matcher = block =~ /(?m)^\s*([A-Za-z][A-Za-z0-9_.]*)\s*=\s*list\s*\(/

    matcher.each { m ->
        def name = m[1].toString()

        // Avoid accidentally capturing nested keys inside each group definition.
        if (!(name in ['groups', 'comparisons', 'longitudinal'])) {
            groups << name
        }
    }

    groups = groups.unique()

    if (groups.isEmpty()) {
        error "--group_cols was not provided and no group names could be parsed from ${cfg}"
    }

    return groups
}

def requireExistingPaths(Map paths) {
    paths.each { label, pathObj ->
        if (!pathObj.exists()) {
            error "${label} does not exist: ${pathObj}"
        }
    }
}

workflow {

    /*
     * Required user/project parameters
     */

    if (!params.base_dir) {
        error """
Missing required parameter: --base_dir

Example:
  --base_dir /lisc/data/scratch/ccr/CR_projects/phiper
"""
    }

    if (!params.project_name) {
        error """
Missing required parameter: --project_name

Example:
  nextflow run main.nf \\
    -profile lisc \\
    --project_name IBD-Berlin \\
    --metadata_file IBD-Berlin_metadata.csv
"""
    }

    if (!params.metadata_file) {
        error """
Missing required parameter: --metadata_file

metadata_file should be the file name inside:
  <project_name>/Metadata/

Example:
  --metadata_file IBD-Berlin_metadata.csv
"""
    }

    /*
     * Pipeline-level reusable paths.
     * These are inside phipflow/, not inside params.base_dir.
     */

    def workflow_src_dir  = "${projectDir}/src/R"
    def workflow_template = "${projectDir}/src/template/phiper_summary_report.qmd"
    def peptide_library = params.peptide_library ?: "${projectDir}/peplib/peptide_library.rds"
    
    /*
     * Project-specific paths.
     * These remain inside params.base_dir/project_name.
     */

    def base_dir_path      = file(params.base_dir)
    def project_dir_path   = file("${params.base_dir}/${params.project_name}")
    def data_dir_path      = file("${params.base_dir}/${params.project_name}/Data")
    def metadata_dir_path  = file("${params.base_dir}/${params.project_name}/Metadata")
    def project_r_dir_path = file("${params.base_dir}/${params.project_name}/R")
    
    def metadata_path      = file("${params.base_dir}/${params.project_name}/Metadata/${params.metadata_file}")
    def exist_path         = file("${params.base_dir}/${params.project_name}/Data/${params.exist_file}")
    def fold_path          = file("${params.base_dir}/${params.project_name}/Data/${params.fold_file}")
    def group_config_path  = file("${params.base_dir}/${params.project_name}/R/group_config.R")

    /*
     * Optional peptide library path, defaulting to phipflow/peplib/peptide_library.rds
     */
    def peptide_library_path = file(peptide_library)
    
    /*
     * Reusable R scripts and QMD template live in phipflow/src.
     */

    def r_create_path   = file("${workflow_src_dir}/01-create_phiper_object.R")
    def r_analysis_path = file("${workflow_src_dir}/02-run_phiper_analysis.R")
    def r_render_path   = file("${workflow_src_dir}/03-render_phiper_reports.R")
    def r_helper_path   = file("${workflow_src_dir}/helper_functions.R")
    def template_path   = file(workflow_template)

    /*
     * Fail early before submitting SLURM jobs.
     */

    requireExistingPaths([
        'base_dir'                   : base_dir_path,
        'Project directory'          : project_dir_path,
        'Project Data directory'     : data_dir_path,
        'Project Metadata directory' : metadata_dir_path,
        'Project R directory'        : project_r_dir_path,
        'Metadata file'              : metadata_path,
        'Exist matrix'               : exist_path,
        'Fold matrix'                : fold_path,
        'group_config.R'             : group_config_path,
        'Peptide library RDS'        : peptide_library_path,
        '01-create_phiper_object.R'  : r_create_path,
        '02-run_phiper_analysis.R'   : r_analysis_path,
        '03-render_phiper_reports.R' : r_render_path,
        'helper_functions.R'         : r_helper_path,
        'Quarto template'            : template_path
    ])

   

    /*
     * Resolve groups.
     */

    def group_cols = parseCsvParam(params.group_cols)

    if (group_cols.isEmpty()) {
        group_cols = discoverGroupsFromConfig(params.base_dir, params.project_name)
        log.info "No --group_cols supplied; using all group_definitions from ${params.project_name}/R/group_config.R"
    }

    def report_group_cols = parseCsvParam(params.report_group_cols)

    if (report_group_cols.isEmpty()) {
        report_group_cols = group_cols
    }

    /*
     * Resolve parquet output name.
     */

    def parquet_name = params.out_parquet ?: "${params.project_name}.parquet"

    /*
     * Logging.
     */

    log.info "PHIPER workflow"
    log.info "  phipflow projectDir : ${projectDir}"
    log.info "  workflow_src_dir    : ${workflow_src_dir}"
    log.info "  workflow_template   : ${workflow_template}"
    log.info "  base_dir            : ${params.base_dir}"
    log.info "  project_name        : ${params.project_name}"
    log.info "  metadata_file       : ${params.metadata_file}"
    log.info "  parquet_name        : ${parquet_name}"
    log.info "  group_cols          : ${group_cols.join(', ')}"
    log.info "  report_group_cols   : ${report_group_cols.join(', ')}"
    log.info "  peptide_library     : ${peptide_library}"
    log.info "  report BASE_DIR     : ${params.base_dir}/${params.project_name}/results"
    log.info "  use_modules         : ${params.use_modules}"
    log.info "  container           : ${params.container ?: 'none'}"
    /*
     * 01-create_phiper_object.R
     */

    CREATE_PHIPER_OBJECT(
        params.base_dir,
        params.project_name,
        params.exist_file,
        params.fold_file,
        params.metadata_file,
        parquet_name,
        params.sample_prefix,
        params.peptide_col,
        params.replace_inf,
        workflow_src_dir,
        params.use_modules
    )

    /*
     * 02-run_phiper_analysis.R
     *
     * The object creation marker is combined with the list of group columns,
     * so analysis starts only after the parquet creation step has completed.
     */

    def analysis_input_ch = CREATE_PHIPER_OBJECT.out.parquet_marker.combine(Channel.fromList(group_cols))
    RUN_PHIPER_ANALYSIS(
        analysis_input_ch,
        params.base_dir,
        params.project_name,
        parquet_name,
        params.all,
        params.default_longitudinal,
        params.manual_comparison_file,
        params.force,
        workflow_src_dir,
        peptide_library,
        params.use_modules
    )

    /*
     * 03-render_phiper_reports.R
     *
     * collect() makes rendering wait until all requested group-column analyses finish.
     */

    RENDER_PHIPER_REPORTS(
        RUN_PHIPER_ANALYSIS.out.analysis_marker.collect(),
        params.base_dir,
        params.project_name,
        report_group_cols.join(','),
        workflow_src_dir,
        workflow_template,
        params.use_modules
    )
}