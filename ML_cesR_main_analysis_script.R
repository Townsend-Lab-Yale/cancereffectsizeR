# load in packages required for cancereffectsizeR ----
library(cancereffectsizeR)


# Load in mutation data ----

NCI_data <- read.delim(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf",stringsAsFactors = F,skip = 5,header = T)

YG_data <- read.delim(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/mutationsTN_62_Lung_Adenocarcinoma_(Yale___MD_Anderson).maf",header = T,stringsAsFactors = F)

# Preprocess for analysis ----
# These functions will be automatically loaded into the package space,
# but let's source them into this working file.


# TODO: check all these preprocessing steps to make sure they are necessary and efficient

## convert NCI data to hg19 ----
source("R/hg_converter.R")

NCI_data <- hg_converter(chain = "~/Box Sync/ML_cancereffectsizeR_data/hg38ToHg19.over.chain",maf_to_convert =  NCI_data)

## merge the data frames along common headers ----
source("R/merging_NCI_and_local_MAF_files.R")

MAF_input <- merging_TCGA_and_local_MAFdata_function(NCI_data = NCI_data,Local_data = YG_data,check_for_same_tumor = T)

## add unique tumor names ----
# is this necessary anymore? Want to remove superfluous functions and code.
source("R/unique_tumor_addition.R")

MAF_input <- unique_tumor_addition_function(MAF.file = MAF_input,non.TCGA.characters.to.keep = "all")

## tumor allele adder ----
source("R/tumor_allele_adder.R")

MAF_input <- tumor_allele_adder(MAF = MAF_input)

## remove DNP, TNP ----
source("R/DNP_TNP_remover.R")

MAF_input <- DNP_TNP_remover(MAF = MAF_input,delete_recur = T)

## save for easy recall -----
save(MAF_input, file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/LUAD_MAF_for_analysis.RData")


# End of preprocessing ----

# load("~/Box Sync/ML_cancereffectsizeR_data/LUAD/LUAD_MAF_for_analysis.RData")

# Start of analysis ----

# will use functions and sub-functions to complete each step
# have clean and consistent inputs/outputs at each step
# will eventually wrap this all in a single function, that takes preprocessed MAF and outputs selection results.

# 1. Find trinucleotide signature weights ----
## Use deconstructSigs
### Keep formatting consistent (i.e. A[A>T]G) with deconstructSigs output (not original cancereffectsizeR
###

source("R/deconstructSigs_input_preprocess.R")

MAF_input_deconstructSigs_preprocessed <- deconstructSigs_input_preprocess(MAF = MAF_input)

# This provides us with the number of each trinucleotide in each tumor
# Even the ones with less than 50 variants
# Useful for the nearest neighbor calculation
trinuc_breakdown_per_tumor <- deconstructSigs::mut.to.sigs.input(mut.ref =
                                                         MAF_input_deconstructSigs_preprocessed,
                                                       sample.id = "Unique_patient_identifier",
                                                       chr = "Chromosome",
                                                       pos = "Start_Position",
                                                       ref = "Reference_Allele",
                                                       alt = "Tumor_allele")



# 2. Find nearest neighbor to tumors with < 50 mutations, assign identical weights as neighbor ----

# 3. Assign genes to MAF_input ----


# 4. Calculate mutation rates at gene-level (dndscv) ----

# 5. For each substitution, calculate the gene- and tumor- and mutation-specific mutation rate----


# 6. Maximum-likelihood framework

















# test <- cancereffectsizeR::trinuc_profile_function_with_weights(input.MAF = MAF_input)






