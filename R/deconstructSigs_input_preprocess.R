#' Preprocess for deconstructSigs
#'
#' This function removes recurrent variants from the MAF file and returns a MAF of non-recurrent variants.
#'
#' @param MAF MAF file to have recurrent variants removed
#' @param sample_ID_column column in MAF with sample ID data
#' @param ref_column column in MAF with reference allele data
#' @param alt_column column in MAF with alternative allele data
#' @param pos_column column in MAF with chromosome nucleotide location data
#' @param chr_column column in MAF with chromosome data
#'
#' @return a MAF of non-recurrent variants
#' @export
#'
#' @examples
deconstructSigs_input_preprocess <- function(MAF,
                                             sample_ID_column="Unique_patient_identifier",
                                             chr_column = "Chromosome",
                                             pos_column = "Start_Position",
                                             ref_column = "Reference_Allele",
                                             alt_column = "Tumor_allele") {

# function to prepare for deconstructSigs

# need to:
## 1. make sure all input is SNV
## 2. remove recurrent varaints
## 3. and "chr" to chromosome

  # keeping just single nucleotide variants

  MAF <- MAF[which(MAF[, ref_column] %in% c('A', 'T', 'C', 'G') & MAF[, alt_column] %in% c('A', 'T', 'C', 'G')),]

  # removing recurrent variants

  duplicated_vec_first <- duplicated(MAF[,c(pos_column,chr_column,alt_column)])
  duplicated_vec_last <- duplicated(MAF[,c(pos_column,chr_column,alt_column)],fromLast=T)
  duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)

  MAF <- MAF[-duplicated_vec_pos,]

  MAF[,chr_column] <- paste("chr",trimws(MAF[,chr_column]),sep="")

  return(MAF)
}
