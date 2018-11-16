# load in packages required for cancereffectsizeR ----
library(cancereffectsizeR)


# Load in mutation data ----

NCI_data <- read.delim(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf",stringsAsFactors = F,skip = 5,header = T)

YG_data <- read.delim(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/mutationsTN_62_Lung_Adenocarcinoma_(Yale___MD_Anderson).maf",header = T,stringsAsFactors = F)

# Preprocess for analysis ----
# These functions will be automatically loaded into the package space,
# but let's source them into this working file.


# TODO: check all these preprocessing steps to make sure they are necessary and efficient

# convert NCI data to hg19 ----
source("R/hg_converter.R")

NCI_data <- hg_converter(chain = "~/Box Sync/ML_cancereffectsizeR_data/hg38ToHg19.over.chain",maf_to_convert =  NCI_data)

# merge the data frames along common headers ----
source("R/merging_NCI_and_local_MAF_files.R")

MAF_input <- merging_TCGA_and_local_MAFdata_function(NCI_data = NCI_data,Local_data = YG_data,check_for_same_tumor = T)

# add unique tumor names ----
# is this necessary anymore? Want to remove superfluous functions and code.
source("R/unique_tumor_addition.R")

MAF_input <- unique_tumor_addition_function(MAF.file = MAF_input,non.TCGA.characters.to.keep = "all")

# tumor allele adder ----
source("R/tumor_allele_adder.R")

MAF_input <- tumor_allele_adder(MAF = MAF_input)

# remove DNP, TNP ----
source("R/DNP_TNP_remover.R")

MAF_input <- DNP_TNP_remover(MAF = MAF_input,delete_recur = T)

# save for easy recall -----
save(MAF_input, file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/LUAD_MAF_for_analysis.RData")


# End of preprocessing ----

# Start of analysis ----
