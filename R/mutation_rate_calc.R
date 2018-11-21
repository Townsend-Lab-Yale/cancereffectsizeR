# function to calculate mutation rate
# given inputMAF, gene, gene mutation rate, and trinucs of all tumors.
# only thing specified is the gene

# requires "data/gene_trinuc_comp.RData"
# requires "data/"data/AA_mutation_list.RData"

# MAF <- MAF_input
# gene <- "BRINP3"
# gene_mut_rate<- mutrates
# tumor_trinucs <- trinuc_proportion_matrix

mutation_rate_calc <- function(MAF, gene, gene_mut_rate, tumor_trinucs){

  mutation_rate_nucs <- matrix(nrow=nrow(trinuc_proportion_matrix),ncol=ncol(trinuc_proportion_matrix),data = NA)
  rownames(mutation_rate_nucs) <- rownames(trinuc_proportion_matrix); colnames(mutation_rate_nucs) <- colnames(trinuc_proportion_matrix)

  if(0 %in% gene_trinuc_comp[[gene]]$gene_trinuc$count){
    gene_trinuc_comp[[gene]]$gene_trinuc$count <- gene_trinuc_comp[[gene]]$gene_trinuc$count + 1
  }

  for(i in 1:nrow(mutation_rate_nucs)){
    mutation_rate_nucs[i,] <- ((gene_trinuc_comp[[gene]]$gene_trinuc$count * trinuc_proportion_matrix[i,] / mean(gene_trinuc_comp[[gene]]$gene_trinuc$count * trinuc_proportion_matrix[i,]))) * gene_mut_rate[gene]

  }

  # mutation_rate_nucs is now the rate of each trinucleotide in each tumor for this gene

  # need to find unique variants and then rates


  this_MAF <- subset(MAF, Gene_name==gene & Reference_Allele %in% c("A","T","G","C") & Tumor_allele %in% c("A","T","G","C")) # subset the MAF into just this gene
  # this_MAF <- this_MAF[!duplicated(this_MAF[,c("Start_Position","Tumor_allele")]),] # find all unique variants
  # this_MAF$unique_mutation <- paste(this_MAF$Gene_name, this_MAF$Chromosome, this_MAF$Start_Position, this_MAF$Tumor_allele)





  this_MAF <- this_MAF[!duplicated(this_MAF[,c("unique_variant_ID_AA")]),]

  # Need to account for different nucleotide changes giving the same amino acid
  # Assign amino acids here
  # Need to give this information back to the main function to count total variants in population


  # as.numeric(gsub("\\D", "", dndscvout$annotmuts$ntchange))

  mutation_rate_matrix <- matrix(nrow=nrow(trinuc_proportion_matrix), ncol=nrow(this_MAF))
  rownames(mutation_rate_matrix) <- rownames(trinuc_proportion_matrix); colnames(mutation_rate_matrix) <- this_MAF$unique_variant_ID_AA


  for(i in 1:nrow(mutation_rate_matrix)){
    for(j in 1:ncol(mutation_rate_matrix)){
      mutation_rate_matrix[i,j] <- ifelse(this_MAF$is_coding[j],sum(mutation_rate_nucs[i,as.character(unlist(AA_mutation_list[[this_MAF$amino_acid_context[j]]][this_MAF$coding_variant_AA_mut[j]]))]),mutation_rate_nucs[i,this_MAF$trinuc_dcS[j]])
    }
  }


 return(mutation_rate_matrix)

}

#
# dndscvout$annotmuts$nuc_pos <- as.numeric(gsub("\\D", "", dndscvout$annotmuts$ntchange))
# dndscvout$annotmuts$codon_pos <- (dndscvout$annotmuts$nuc_pos %% 3)
#
# head(dndscvout$annotmuts)
#
#


