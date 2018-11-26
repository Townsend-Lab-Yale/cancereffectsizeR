# function to find trinucleotides for when a variant is up against an exon edge.


AA_translations$Nucs <- as.character(AA_translations$Nucs)

AA_translations$AA_short <- as.character(AA_translations$AA_short)



mutation_finder <- function(RefCDS_instance,MAF_input_row){

  # data structure for the codon

  codon_matrix <- matrix(nrow=3,ncol=3)
  rownames(codon_matrix) <- c("pos1","pos2","pos3") ; colnames(codon_matrix) <- c("upstream","main","downstream")


  if(MAF_input_row$codon_pos == 1){
    codon_matrix["pos1","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
    codon_matrix["pos2","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position+1,stop = MAF_input_row$nuc_position+1))
    codon_matrix["pos3","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position+2,stop = MAF_input_row$nuc_position+2))

    codon_matrix["pos1","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
    codon_matrix["pos2","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position+1,stop = MAF_input_row$nuc_position+1))
    codon_matrix["pos3","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position+2,stop = MAF_input_row$nuc_position+2))

    codon_matrix["pos1","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
    codon_matrix["pos2","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position+1,stop = MAF_input_row$nuc_position+1))
    codon_matrix["pos3","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position+2,stop = MAF_input_row$nuc_position+2))
  }

  if(MAF_input_row$codon_pos == 2){
    codon_matrix["pos1","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position-1,stop = MAF_input_row$nuc_position-1))
    codon_matrix["pos2","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
    codon_matrix["pos3","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position+1,stop = MAF_input_row$nuc_position+1))

    codon_matrix["pos1","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position-1,stop = MAF_input_row$nuc_position-1))
    codon_matrix["pos2","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
    codon_matrix["pos3","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position+1,stop = MAF_input_row$nuc_position+1))

    codon_matrix["pos1","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position-1,stop = MAF_input_row$nuc_position-1))
    codon_matrix["pos2","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
    codon_matrix["pos3","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position+1,stop = MAF_input_row$nuc_position+1))
  }

  if(MAF_input_row$codon_pos == 3){
    codon_matrix["pos1","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position-2,stop = MAF_input_row$nuc_position-2))
    codon_matrix["pos2","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position-1,stop = MAF_input_row$nuc_position-1))
    codon_matrix["pos3","main"] <- as.character(substr(RefCDS_instance$seq_cds,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))

    codon_matrix["pos1","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position-2,stop = MAF_input_row$nuc_position-2))
    codon_matrix["pos2","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position-1,stop = MAF_input_row$nuc_position-1))
    codon_matrix["pos3","upstream"] <- as.character(substr(RefCDS_instance$seq_cds1up,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))

    codon_matrix["pos1","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position-2,stop = MAF_input_row$nuc_position-2))
    codon_matrix["pos2","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position-1,stop = MAF_input_row$nuc_position-1))
    codon_matrix["pos3","downstream"] <- as.character(substr(RefCDS_instance$seq_cds1down,start = MAF_input_row$nuc_position,stop = MAF_input_row$nuc_position))
  }

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

