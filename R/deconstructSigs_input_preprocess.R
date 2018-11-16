#' deconstructSigs_input_preprocess
#'
#' Gathers all SNV and removes recurrent variants
#'
#' @param MAF
#' @param sample.ID.column
#' @param chr.column
#' @param pos.column
#' @param ref.column
#' @param alt.column
#'
#' @return
#' @export
#'
#' @examples
deconstructSigs_input_preprocess <- function(MAF,
                                             sample.ID.column="Unique_patient_identifier",
                                             chr.column = "Chromosome",
                                             pos.column = "Start_Position",
                                             ref.column = "Reference_Allele",
                                             alt.column = "Tumor_allele"){

# function to prepare for deconstructSigs

# need to:
## 1. make sure all input is SNV
## 2. remove recurrent varaints
## 3. and "chr" to chromosome

  # keeping just single nucleotide variants

  MAF <- MAF[which(MAF[, ref.column] %in% c('A', 'T', 'C', 'G') & MAF[, alt.column] %in% c('A', 'T', 'C', 'G')),]

  # removing recurrent variants


  duplicated.vec.first <- duplicated(MAF[,c(pos.column,chr.column,alt.column)])
  duplicated.vec.last <- duplicated(MAF[,c(pos.column,chr.column,alt.column)],fromLast=T)

  duplicated.vec.pos <- which(duplicated.vec.first | duplicated.vec.last)


  MAF <- MAF[-duplicated.vec.pos,]

  MAF[,chr.column] <- paste("chr",trimws(MAF[,chr.column]),sep="")

  return(MAF)

}
