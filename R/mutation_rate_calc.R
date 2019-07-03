#' Mutation rate calculation
#'
#' This function calculates the mutation rates of all the unique varaints within a gene given gene-level mutation rate and tumor-specific trinucleotide contexts.
#'
#' @param this_MAF The subset MAF (just the gene to be analyzed) file to extract mutational data from
#' @param gene Gene name in question
#' @param trinuc_proportion_matrix matrix constructed from deconstructSigs output, containing proportion of each trinucleotide mutated in each tumor
#' @param gene_trinuc_comp list containing matrices of counts of each trinuc in every gene
#' @param RefCDS RefCDS loaded in from data("RefCDS_TP53splice",package = "cancereffectsizeR")
#' @param progression CESProgressions
#' @param gene_mut_rate mutation rate at the gene-level
#'
#' @return
#' @export
#'
#' @examples
mutation_rate_calc <- function(this_MAF,
                               gene, gene_mut_rate,
                               trinuc_proportion_matrix,
                               gene_trinuc_comp,
                               RefCDS,
                               relative_substitution_rates,
                               tumor_specific_rate=F,
                               all_tumors,
                               progressions
                               ){

  mutation_rate_nucs <- matrix(nrow=nrow(trinuc_proportion_matrix),ncol=ncol(trinuc_proportion_matrix),data = NA)
  rownames(mutation_rate_nucs) <- rownames(trinuc_proportion_matrix)
  colnames(mutation_rate_nucs) <- colnames(trinuc_proportion_matrix)

  # if there are no substitutions within a tumor after our preprocessing,
  # do not calculate the mutation rate within that tumor

  if(length(which(!rownames(mutation_rate_nucs) %in% all_tumors)) > 0){
    mutation_rate_nucs <- as.matrix(mutation_rate_nucs[-which(!rownames(mutation_rate_nucs) %in% all_tumors),])
    colnames(mutation_rate_nucs) <- colnames(trinuc_proportion_matrix)
  }


  if(0 %in% gene_trinuc_comp[[gene]]$gene_trinuc$count){
    gene_trinuc_comp[[gene]]$gene_trinuc$count <- gene_trinuc_comp[[gene]]$gene_trinuc$count + 1
  }



  norm_constant = sum(gene_trinuc_comp[[gene]]$gene_trinuc[,"count"])
  calc_normalizers = function(x) {
    return(sum(gene_trinuc_comp[[gene]]$gene_trinuc$count * x))
  }
  normalizers = apply(trinuc_proportion_matrix, 1, calc_normalizers) / norm_constant

  for(i in 1:nrow(mutation_rate_nucs)){
    mutation_rate_nucs[i,] <- (trinuc_proportion_matrix[i,] / normalizers[i]) * gene_mut_rate[[get_progression_number(progressions, rownames(mutation_rate_nucs)[i])]][gene]
  }

  # mutation_rate_nucs is now the rate of each trinucleotide in each tumor for this gene
  # need to find unique variants and then rates


  # frequency of all variants in the different subsets.
  variant_freq_list <- vector(mode = "list", length = length(progressions@order))
  names(variant_freq_list) <- names(progressions@order)
  for(this_level in 1:length(progressions@order)){
    variant_freq_list[[this_level]] <- c(
      table(this_MAF[this_MAF$Unique_Patient_Identifier %in% get_progression_tumors(progressions, this_level),
                     "unique_variant_ID_AA"]),
      setNames(rep(0,length(base::setdiff(unique(this_MAF[,
                                                          "unique_variant_ID_AA"]),
                                          unique(this_MAF[this_MAF$Unique_Patient_Identifier %in% get_progression_tumors(progressions, this_level),
                                                          "unique_variant_ID_AA"])))),
               base::setdiff(unique(this_MAF[,
                                             "unique_variant_ID_AA"]),
                             unique(this_MAF[this_MAF$Unique_Patient_Identifier %in% get_progression_tumors(progressions, this_level),
                                                          "unique_variant_ID_AA"])))
      )
  }

  this_MAF <- this_MAF[!duplicated(this_MAF[,c("unique_variant_ID_AA")]),]

  # Need to account for different nucleotide changes giving the same amino acid
  # Assign amino acids here
  # Need to give this information back to the main function to count total variants in population


  # as.numeric(gsub("\\D", "", dndscvout$annotmuts$ntchange))

  mutation_rate_matrix <- matrix(nrow=nrow(trinuc_proportion_matrix), ncol=nrow(this_MAF))
  rownames(mutation_rate_matrix) <- rownames(trinuc_proportion_matrix)
  colnames(mutation_rate_matrix) <- this_MAF$unique_variant_ID_AA

  if(length(which(!rownames(mutation_rate_matrix) %in% all_tumors)) > 0){
    mutation_rate_matrix <- as.matrix(mutation_rate_matrix[-which(!rownames(mutation_rate_matrix) %in% all_tumors),])
    colnames(mutation_rate_matrix) <- this_MAF$unique_variant_ID_AA
  }



  for(i in 1:nrow(mutation_rate_matrix)){
    for(j in 1:ncol(mutation_rate_matrix)){

      if(this_MAF$next_to_splice[j]){

        if(this_MAF$is_coding[j]){
          mutation_rate_matrix[i,j] <- sum(mutation_rate_nucs[i,cancereffectsizeR::mutation_finder(RefCDS_instance = RefCDS[[gene]],MAF_input_row = this_MAF[j,])])
        } else {
          mutation_rate_matrix[i,j] <- mutation_rate_nucs[i,this_MAF$trinuc_dcS[j]]
        }

      } else {

        if(this_MAF$is_coding[j]) {
          mutation_rate_matrix[i,j] <- sum(mutation_rate_nucs[i,as.character(unlist(AA_mutation_list[[this_MAF$amino_acid_context[j]]][this_MAF$coding_variant_AA_mut[j]]))])
        } else {
          mutation_rate_matrix[i,j] <- mutation_rate_nucs[i,this_MAF$trinuc_dcS[j]]
        }

      }
    }
  }

  # adjusting tumor-specific rate to reflect the tumor-specific non-recurrent substitution load

  # mutation_rate_matrix_rel <- (mutation_rate_matrix) * t(relative_substitution_rates[rownames(mutation_rate_matrix)])
  if(tumor_specific_rate){
    for(this_row in 1:nrow(mutation_rate_matrix)){
      mutation_rate_matrix[this_row,] <-  mutation_rate_matrix[this_row,] * relative_substitution_rates[rownames(mutation_rate_matrix)[this_row]]
    }
  }
  unsure_genes_vec <- this_MAF$unsure_gene_name



  return(list(mutation_rate_matrix=mutation_rate_matrix,
              unsure_genes_vec=unsure_genes_vec,
              variant_freq = variant_freq_list,
              unique_variant_ID_vec = this_MAF$unique_variant_ID))

}
