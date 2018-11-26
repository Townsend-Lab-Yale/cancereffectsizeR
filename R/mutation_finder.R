# function to find trinucleotides for when a variant is up against an exon edge.

# OR2L13 M1T

# > MAF_input[i,]
# Unique_patient_identifier Chromosome Start_Position
# 128852              TCGA-MP-A4TA          1      248262679
# Reference_Allele Tumor_allele Gene_name strand triseq
# 128852                T            C    OR2L13      +    ATG
# unique_variant_ID is_coding Tumor_allele_correct_strand trinuc_dcS
# 128852     1 248262679 C      TRUE                           C    A[T>C]G
# nuc_variant coding_variant nuc_position codon_pos
# 128852         T2C            M1T            2         2
# amino_acid_context unique_variant_ID_AA coding_variant_AA_mut
# 128852            TCATGGA                  M1T                   Thr

# MAF_input_row <- MAF_input[118273,]
# RefCDS_instance <- RefCDS[["AL161645.2"]]

AA_translations$Nucs <- as.character(AA_translations$Nucs)

AA_translations$AA_short <- as.character(AA_translations$AA_short)


# any(abs(RefCDS[["OR2L13"]]$intervals_cds - 248262679) <= 3)

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

mutation_finder(RefCDS_instance = RefCDS[["OR2L13"]],MAF_input_row <- MAF_input[118273,])


#
# substr(as.character(RefCDS[["NFE2L2"]]$seq_cds),start = 593,stop = 595)
#
# RefCDS[["NFE2L2"]]$intervals_cds
#
# chr2cds(pos = 178098733,cds_int = RefCDS[["NFE2L2"]]$intervals_cds,strand = -1)
#
# chr2cds(pos = 178097120,cds_int = RefCDS[["NFE2L2"]]$intervals_cds,strand = -1)
#
# substr(as.character(RefCDS[["NFE2L2"]]$seq_cds),start = 593,stop = 595)
# substr(as.character(RefCDS[["NFE2L2"]]$seq_cds1down),start = 592,stop = 594)
#
#
#
