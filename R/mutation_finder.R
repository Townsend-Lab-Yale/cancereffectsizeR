#' mutation_finder
#'
#' Finds original trinucleotides of variants near splice sites.
#'
#' @param RefCDS_instance Slice of the RefCDS list for this gene
#' @param MAF_input_row Row from the MAF being analyzed with the mutation in question.
#'
#' @return
#' @export
#'
#' @examples
mutation_finder <- function(RefCDS_instance,MAF_input_row){
  AA_translations$Nucs <- as.character(AA_translations$Nucs)

  AA_translations$AA_short <- as.character(AA_translations$AA_short)

  # data structure for the codon

  codon_matrix <- matrix(nrow=3,ncol=3)
  rownames(codon_matrix) <- c("pos1","pos2","pos3")
  colnames(codon_matrix) <- c("upstream","main","downstream")
  shift <- MAF_input_row$codon_pos - 1
  RefCDS_instance$seq_cds <- as.character(RefCDS_instance$seq_cds)
  RefCDS_instance$seq_cds1up <- as.character(RefCDS_instance$seq_cds1up)
  RefCDS_instance$seq_cds1down <- as.character(RefCDS_instance$seq_cds1down)

  codon_matrix[c("pos1", "pos2", "pos3"),"main"] <- unlist(strsplit(RefCDS_instance$seq_cds, split=""))[(MAF_input_row$nuc_position-shift):(MAF_input_row$nuc_position-shift+2)]
  codon_matrix[c("pos1", "pos2", "pos3"),"upstream"] <- unlist(strsplit(RefCDS_instance$seq_cds1up, split=""))[(MAF_input_row$nuc_position-shift):(MAF_input_row$nuc_position-shift+2)]
  codon_matrix[c("pos1", "pos2", "pos3"),"downstream"] <- unlist(strsplit(RefCDS_instance$seq_cds1down, split=""))[(MAF_input_row$nuc_position-shift):(MAF_input_row$nuc_position-shift+2)]

  nt <- c("A","T","G","C")
  trinucs_for_rate <- NULL
  original_codon <- codon_matrix[,2]

  for(position in 1:3){
    these_nucs <- setdiff(nt,codon_matrix[position,2])
   for(nucleotide in 1:3){

     new_codon <- original_codon
     new_codon[position] <- these_nucs[nucleotide]

     if(AA_translations[which(AA_translations$Nucs==paste(new_codon,collapse = "")),"AA_short"] == MAF_input_row$coding_variant_AA_mut){

      trinucs_for_rate <- c(trinucs_for_rate,trinuc_translator[paste(codon_matrix[position,"upstream"],codon_matrix[position,"main"],codon_matrix[position,"downstream"],":",these_nucs[nucleotide],sep="",collapse = ""),"deconstructSigs_format"])
     }

   }

  }

return(trinucs_for_rate)
}

# function to find trinucleotides for when a variant is up against an exon edge.
