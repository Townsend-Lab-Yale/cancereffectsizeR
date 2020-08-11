#' Reads TCGA Tumor Sample Barcodes and consolidates barcodes into patient-specific
#' identifiers, such patients with multiple distinct tumors in the analysis share
#' share just one identifer. Sample barcodes that don't start with "TCGA" are unaltered.
#'
#' @param tumor_sample_barcodes an MAF-style Tumor_Sample_Barcode vector, of which some or all are TCGA samples
#' @param non_TCGA_characters_to_keep numeric number of characters to keep
#' in "Unique_Patient_Identifier" column for sample names that
#' do not start with "TCGA". If non-numeric object then defaults to
#' all characters
#' @param sum_stats Boolean, T prints summary() of number of mutations per tumor
#' @param figures Boolean, T saves histogram of sum.stats and saves as "Figures/mutations_per_tumor.pdf"
#' @keywords internal


# Adds column to the MAF that contains the unique tumor ID from the $Unique_Patient_Identifier
consolidate_tcga_tumors_by_patient <- function(tumor_sample_barcodes, non_TCGA_characters_to_keep="all", sum_stats=T,figures=F){
  TCGA_tumors <- grep(pattern = "TCGA",x = tumor_sample_barcodes)
  other_tumors <- setdiff(1:length(tumor_sample_barcodes),TCGA_tumors)

  if (length(TCGA_tumors) > 1) {
    msg = paste0("TCGA samples appear to be present (tumor IDs start with \"TCGA\"). ",
             "Distinct tumors that by TCGA naming conventions appear to derive from ",
             "the same patient will be consolidated under one ID per patient.")
    message(msg)
  } else {
    return(tumor_sample_barcodes)
  }


  if(class(non_TCGA_characters_to_keep)=="numeric"){
    tumor_sample_barcodes[TCGA_tumors] <- substr(tumor_sample_barcodes[TCGA_tumors],1,12)
    tumor_sample_barcodes[other_tumors] <- substr(tumor_sample_barcodes[other_tumors],1,non_TCGA_characters_to_keep)
  }else{
    tumor_sample_barcodes[TCGA_tumors] <- substr(tumor_sample_barcodes[TCGA_tumors],1,12)
    tumor_sample_barcodes[other_tumors] <- tumor_sample_barcodes[other_tumors]
  }

  if(sum_stats){
    tumor_mutation_number = as.numeric(table(tumor_sample_barcodes))

    if(figures){
      pdf(file = "Figures/mutations_per_tumor.pdf")
      hist(tumor_mutation_number,breaks=100,xlab="Number of mutations per tumor",
           main="Histogram of the number of mutations in each tumor")
      dev.off()
    }

    message("Summary statistics of the number of mutations per unique tumor:")
    print(summary(tumor_mutation_number))
  }

  return(tumor_sample_barcodes)
}
