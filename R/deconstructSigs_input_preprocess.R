#' deconstructSigs_input_preprocess
#'
#' Gathers all SNV and removes recurrent variants
#'
#' @param MAF
#' @param sample_ID_column
#' @param chr_column
#' @param pos_column
#' @param ref_column
#' @param alt_column
#'
#' @return
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
