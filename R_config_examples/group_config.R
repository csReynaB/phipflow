# ------------------------------------------------------------------------------
# Group definitions and ordering
# ------------------------------------------------------------------------------

# for now this filed must be called group_config.R per project

# `group_definitions` is a named list. Each top-level entry defines one analysis
# group to be run by the workflow.
#
# 1. Each entry must have a name, e.g. `Treatment_Timepoint`.
#    This name is used internally by the workflow and for organizing outputs.
#    It can be different from `group_col`, but using the same name is recommended
#    to make results easier to track.
#
# 2. `group_col` must match a column name in the metadata.
#
# 3. `groups` defines the levels from `group_col` to include in the analysis.
#    These do not need to include all values present in the metadata column.
#
# 4. `comparisons` is optional.
#    - If omitted, all possible pairwise comparisons among `groups` are generated.
#    - If provided, it must be a list of pairs.
#    - Each label used in `comparisons` must also be present in `groups`.
#
# 5. `longitudinal` is optional.
#    - If omitted, it defaults to FALSE.
#    - This default can be changed globally with the config parameter:
#        default_longitudinal = TRUE
#    - Set to TRUE for paired/longitudinal comparisons, where samples are matched
#      by `subject_id`.

group_definitions <- list(
  Treatment_Timepoint = list(
    group_col = "Treatment_Timepoint",
    groups = c("Syst_BL", "Syst_FU"),
    comparisons = list(    
      c("Syst_BL", "Syst_FU") 
    ),
    longitudinal = c(TRUE)
  ),
  Treatment_type_Timepoint = list(
    group_col = "Treatment_type_Timepoint",
    groups = c("CHT_BL", "IMM_BL", "CHT_FU", "IMM_FU"),
    comparisons = list(
      c("CHT_BL", "CHT_FU"), 
      c("IMM_BL", "IMM_FU") 
    ),
    longitudinal = c(TRUE, TRUE)
  ),
  Cohort_Timepoint = list(
    group_col = "Cohort_Timepoint",
    groups = c("Cis_BL", "Carbo_BL", "ICI_BL", "EVP_BL", "Cis_FU", "Carbo_FU", "ICI_FU", "EVP_FU" ),
    comparisons = list(
      c("Cis_BL", "Cis_FU"),
      #c("Carbo_BL", "Carbo_FU"),
      c("ICI_BL", "ICI_FU"),
      c("EVP_BL", "EVP_FU")
    ),
    longitudinal = c(TRUE, TRUE, TRUE)
  ),
  PFS_12mo_Timepoint = list(
    group_col = "PFS_12mo_Timepoint",
    groups = c("S_BL", "P_BL", "S_FU", "P_FU"),
    comparisons = list(
      c("S_BL", "S_FU"),
      c("P_BL", "P_FU")
    ),
    longitudinal = c(TRUE, TRUE)
  ),
  Setting_stratified_Timepoint = list(
    group_col = "Setting_stratified_Timepoint",
    groups = c("Adjuvant_BL", "Neoadjuvant_BL", "Metastatic_BL", "Adjuvant_FU", "Neoadjuvant_FU", "Metastatic_FU"),
    comparisons = list(
      #c("Adjuvant_BL", "Adjuvant_FU"),
      #c("Neoadjuvant_BL", "Neoadjuvant_FU"),
      c("Metastatic_BL", "Metastatic_FU")
    ),
    longitudinal = c(TRUE)
  ),
  VitalStatus_stratified_Timepoint = list(
    group_col = "VitalStatus_stratified_Timepoint",
    groups = c("Alive_BL", "Alive_FU"),
    comparisons = list(
      c("Alive_BL", "Alive_FU")
    ),
    longitudinal = c(TRUE)
  ),
  OS_group_12mo_Timepoint = list(
    group_col = "OS_group_12mo_Timepoint",
    groups = c(">12 mo_BL", "≤12 mo_BL", ">12 mo_FU", "≤12 mo_FU"),
    comparisons = list(
      #c(">12 mo_BL", ">12 mo_FU"),
      c("≤12 mo_BL", "≤12 mo_FU")
    ),
    longitudinal = c(TRUE)
  ),
  ResponseToChemo_stratified_Timepoint = list(
    group_col = "ResponseToChemo_stratified_Timepoint",
    groups = c("Progression under chemo_BL", "No progression under chemo_BL", "Progression under chemo_FU", "No progression under chemo_FU"),
    comparisons = list(
      c("Progression under chemo_BL", "Progression under chemo_FU"),
      c("No progression under chemo_BL", "No progression under chemo_FU")
    ),
    longitudinal = c(TRUE, TRUE)
  ),
  DiseaseStage_stratified_Timepoint = list(
    group_col = "DiseaseStage_stratified_Timepoint",
    groups = c("Nodal_BL", "Local_BL", "Metastatic_BL", "Nodal_FU", "Local_FU", "Metastatic_FU"),
    comparisons = list(
      c("Nodal_BL", "Nodal_FU"),
      c("Local_BL", "Local_FU"),
      c("Metastatic_BL", "Metastatic_FU")
    ),
    longitudinal = c(TRUE, TRUE, TRUE)
  ),
  group_test_Timepoint = list(
    group_col = "group_test_Timepoint",
    groups = c("BC_BL", "BC_FU"),
    comparisons = list(
      c("BC_BL","BC_FU")
    ),
   longitudinal = c(TRUE)
  )
)
