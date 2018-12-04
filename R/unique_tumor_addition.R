#' Unique tumor addition function
#'
#' Adds a column called "Unique_patient_identifier" to your MAF file,
#' useful if TCGA codes have more data than sample and you want summary
#' statistics of the number of mutations per tumor.
#'
#' @param MAF_file MAF file with "Tumor_Sample_Barcode" column
#' @param non_TCGA_characters_to_keep numeric number of characters to keep
#' in "Tumor_Sample_Barcode" column for sample names that
#' do not start with "TCGA". If non-numeric object then defaults to
#' all characters
#' @param sum_stats Boolean, T prints summary() of number of mutations per tumor
#' @param figures Boolean, T saves histogram of sum.stats and saves as "Figures/mutations_per_tumor.pdf"
#' @export



# Adds column to the MAF that contains the unique tumor ID from the $Tumor_Sample_Barcode

unique_tumor_addition_function <- function(MAF_file,non_TCGA_characters_to_keep="all",sum_stats=T,figures=F){
  this_maf <- MAF_file
  TCGA_tumors <- grep(pattern = "TCGA",x = this_maf$Tumor_Sample_Barcode)
  other_tumors <- setdiff(1:nrow(this_maf),TCGA_tumors)

  if(class(non_TCGA_characters_to_keep)=="numeric"){
    this_maf$Unique_patient_identifier <- NA
    this_maf$Unique_patient_identifier[TCGA_tumors] <- substr(this_maf$Tumor_Sample_Barcode[TCGA_tumors],1,12)
    this_maf$Unique_patient_identifier[other_tumors] <- substr(this_maf$Tumor_Sample_Barcode[other_tumors],1,non_TCGA_characters_to_keep)
  }else{
    this_maf$Unique_patient_identifier <- NA
    this_maf$Unique_patient_identifier[TCGA_tumors] <- substr(this_maf$Tumor_Sample_Barcode[TCGA_tumors],1,12)
    this_maf$Unique_patient_identifier[other_tumors] <- this_maf$Tumor_Sample_Barcode[other_tumors]
  }

  if(sum_stats){
    unique_patients <- unique(this_maf$Unique_patient_identifier)
    tumor_mutation_number <- NULL
    for(i in 1:length(unique_patients)){
      tumor_mutation_number[i] <- length(which(this_maf$Unique_patient_identifier==unique_patients[i]))
    }

    if(figures){
      pdf(file = "Figures/mutations_per_tumor.pdf")
      hist(tumor_mutation_number,breaks=100,xlab="Number of mutations per tumor",
           main="Histogram of the number of mutations in each tumor")
      dev.off()
    }

    message("Summary statistics of the number of mutations per unique tumor:")
    print(summary(tumor_mutation_number))
  }

  return(this_maf)
}
