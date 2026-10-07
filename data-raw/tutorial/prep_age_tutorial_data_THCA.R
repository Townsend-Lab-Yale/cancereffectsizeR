# Prepare packaged data for the age-dependent selection vignette
#
# This is a maintainer/development script, not part of the user-facing API.
#
# Expected output:
#   inst/continuous_selection_tutorial/TCGA_THCA_age.maf.gz
#   inst/continuous_selection_tutorial/TCGA_THCA_age_clinical.txt

library(cancereffectsizeR)
library(data.table)
library(dplyr)
library(readr)

source_maf <- "../../data/TCGA/TCGA-THCA.maf.gz"
source_clinical <- "../../data/TCGA/TCGA-THCA_clinical.tsv"

dir.create(dirname(source_maf), recursive = TRUE, showWarnings = FALSE)

dir.create(
  "inst/continuous_selection_tutorial",
  recursive = TRUE,
  showWarnings = FALSE
)

#-------------------------------------------------------------------------------
# TCGA MAF
#-------------------------------------------------------------------------------

# Download TCGA-THCA somatic variant data from the current GDC release.
# get_TCGA_project_MAF() excludes non-primary TCGA tumors by default.
if (!file.exists(source_maf)) {
  get_TCGA_project_MAF(
    project = "THCA",
    filename = source_maf
  )
}

TCGA_maf <- fread(source_maf, header = TRUE, sep = "\t")

length(unique(TCGA_maf$Tumor_Sample_Barcode))
length(unique(TCGA_maf$Unique_Patient_Identifier))

# Tumor_Sample_Barcode is a sample-level barcode, while
# Unique_Patient_Identifier is the patient-level ID.
# For TCGA, the first 12 characters of the sample barcode identify the patient.
patient_id_from_barcode <- substr(TCGA_maf$Tumor_Sample_Barcode, 1, 12)

if (!all(patient_id_from_barcode == TCGA_maf$Unique_Patient_Identifier)) {
  stop(
    "Unique_Patient_Identifier does not match the patient ID encoded in ",
    "Tumor_Sample_Barcode."
  )
}

#-------------------------------------------------------------------------------
# TCGA clinical data
#-------------------------------------------------------------------------------

# Preferred option:
# Download clinical data for TCGA-THCA from the GDC Data Portal and save the
# TSV as ../../data/TCGA/TCGA-THCA_clinical.tsv.
#
# Alternatively, if that file is not present and TCGAbiolinks is installed,
# retrieve the same indexed GDC clinical information through the GDC API.

if (file.exists(source_clinical)) {
  TCGA_sample <- read_tsv(
    source_clinical,
    show_col_types = FALSE
  )
} else {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop(
      "Clinical data were not found at ", source_clinical,
      ". Download TCGA-THCA clinical data from the GDC Data Portal, ",
      "or install TCGAbiolinks and rerun this script."
    )
  }

  TCGA_sample <- TCGAbiolinks::GDCquery_clinic(
    project = "TCGA-THCA",
    type = "clinical"
  )
}

TCGA_sample <- as.data.frame(TCGA_sample)

# GDC clinical downloads may name the patient ID column either
# case_submitter_id or submitter_id.
if ("submitter_id" %in% names(TCGA_sample) &&
    !"case_submitter_id" %in% names(TCGA_sample)) {
  names(TCGA_sample)[names(TCGA_sample) == "submitter_id"] <-
    "case_submitter_id"
}

if (!"case_submitter_id" %in% names(TCGA_sample)) {
  stop("Clinical data are missing the patient identifier column.")
}

if (!"age_at_index" %in% names(TCGA_sample)) {
  stop("Clinical data are missing age_at_index.")
}

# The GDC clinical export may contain multiple rows for the same case.
TCGA_sample <- TCGA_sample |>
  dplyr::distinct(case_submitter_id, .keep_all = TRUE)

# Match clinical data to the MAF by patient ID.
TCGA_sample <- TCGA_sample[
  TCGA_sample$case_submitter_id %in%
    TCGA_maf$Unique_Patient_Identifier,
]

# Clean age
TCGA_sample$age_at_index <- as.character(
  TCGA_sample$age_at_index
)

TCGA_sample$age_at_index <- trimws(
  TCGA_sample$age_at_index
)

TCGA_sample$age_at_index[
  TCGA_sample$age_at_index %in%
    c("'--", "--", "", "NA", "N/A", "Not Reported", "not reported")
] <- NA

TCGA_sample$age_at_index <- gsub(
  "[^0-9.]",
  "",
  TCGA_sample$age_at_index
)

TCGA_sample$age_at_index[
  TCGA_sample$age_at_index == ""
] <- NA

TCGA_sample$age_at_index <- as.numeric(
  TCGA_sample$age_at_index
)

# Exclude patients without age.
TCGA_sample <- TCGA_sample[
  !is.na(TCGA_sample$age_at_index),
]

TCGA_sample_id <- unique(
  TCGA_sample$case_submitter_id
)

TCGA_maf <- TCGA_maf[
  TCGA_maf$Unique_Patient_Identifier %in%
    TCGA_sample_id,
]

#-------------------------------------------------------------------------------
# Build the clinical table used by cancereffectsizeR
#-------------------------------------------------------------------------------

TCGA_sample <- TCGA_sample |>
  dplyr::select(
    case_submitter_id,
    age_at_index
  )

colnames(TCGA_sample) <- c(
  "Unique_Patient_Identifier",
  "AGE_AT_SEQ_REPORT"
)

# Final check: MAF and clinical table should contain the same patients.
stopifnot(
  setequal(
    unique(TCGA_maf$Unique_Patient_Identifier),
    unique(TCGA_sample$Unique_Patient_Identifier)
  ),
  nrow(TCGA_sample) ==
    length(unique(TCGA_maf$Unique_Patient_Identifier))
)

#-------------------------------------------------------------------------------
# Write tutorial data
#-------------------------------------------------------------------------------

fwrite(
  TCGA_maf,
  "inst/continuous_selection_tutorial/TCGA_THCA_age.maf.gz",
  sep = "\t",
  quote = FALSE,
  na = "",
  compress = "gzip"
)

fwrite(
  TCGA_sample,
  "inst/continuous_selection_tutorial/TCGA_THCA_age_clinical.txt",
  sep = "\t",
  quote = FALSE,
  na = ""
)

message(
  "Wrote TCGA-THCA tutorial data: ",
  length(unique(TCGA_maf$Unique_Patient_Identifier)),
  " patients."
)
