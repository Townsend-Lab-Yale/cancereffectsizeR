#' Mutation rate calculation
#'
#' This function calculates the mutation rates of all the unique varaints within a gene given gene-level mutation rate and tumor-specific trinucleotide contexts.
#'
#' @param this_MAF The subset MAF (just the gene to be analyzed) file to extract mutational data from
#' @param gene current gene
#' @param trinuc_proportion_matrix matrix constructed from deconstructSigs output, containing proportion of each trinucleotide mutated in each tumor
#' @param progression CESProgressions
#' @param gene_mut_rate mutation rate at the gene-level
#'
#' @return
#'
#' @examples
mutation_rate_calc <- function(this_MAF,
                               gene,
                               gene_mut_rate,
                               trinuc_proportion_matrix,
                               all_tumors,
                               progressions
                               ){

  # trinuc_proportion_matrix: rows = samples, columns = relative frequency of trinucleotide-context-specific mutation (sums to 1)
  # mutation_rate_nucs: rows = samples, columns = gene-specific relative frequency of trinuclotide-context-specific
  ##                    mutations, calculated by normalizing trinuc_proportion_matrix by gene content
  mutation_rate_nucs <- matrix(nrow=nrow(trinuc_proportion_matrix),ncol=ncol(trinuc_proportion_matrix),data = NA)
  rownames(mutation_rate_nucs) <- rownames(trinuc_proportion_matrix)
  colnames(mutation_rate_nucs) <- colnames(trinuc_proportion_matrix)

  # if there are no substitutions within a tumor after our preprocessing,
  # do not calculate the mutation rate within that tumor

  if(length(which(!rownames(mutation_rate_nucs) %in% all_tumors)) > 0){
    mutation_rate_nucs <- as.matrix(mutation_rate_nucs[-which(!rownames(mutation_rate_nucs) %in% all_tumors),])
    colnames(mutation_rate_nucs) <- colnames(trinuc_proportion_matrix)
  }

  trinuc_counts = gene_trinuc_comp[[gene]]
  if (0 %in% trinuc_counts) {
    trinuc_counts = trinuc_counts + 1
  }
  norm_constant = sum(trinuc_counts)

  calc_normalizers = function(x) {
    return(sum(trinuc_counts * x))
  }
  normalizers = apply(trinuc_proportion_matrix, 1, calc_normalizers) / norm_constant



  for(i in 1:nrow(mutation_rate_nucs)){
    mutation_rate_nucs[i,] <- (trinuc_proportion_matrix[i,] / normalizers[i]) * gene_mut_rate[[get_progression_number(progressions, rownames(mutation_rate_nucs)[i])]][gene]
  }

  # mutation_rate_nucs is now the rate of each trinucleotide in each tumor for this gene
  # need to find unique variants and then rates


  this_MAF <- this_MAF[!duplicated(this_MAF[,c("unique_variant_ID_AA")]),]

  # Need to account for different nucleotide changes giving the same amino acid
  # Assign amino acids here
  # Need to give this information back to the main function to count total variants in population



  # mutation rate matrix: rows = tumors, columns = expected relative rates of amino acid changes in gene for tumor
  ## this is determined by summing over the rates for different nucleotides mutations that create same amino acid changes
  mutation_rate_matrix <- matrix(nrow=nrow(trinuc_proportion_matrix), ncol=nrow(this_MAF))
  rownames(mutation_rate_matrix) <- rownames(trinuc_proportion_matrix)
  colnames(mutation_rate_matrix) <- this_MAF$unique_variant_ID_AA

  if(length(which(!rownames(mutation_rate_matrix) %in% all_tumors)) > 0){
    mutation_rate_matrix <- as.matrix(mutation_rate_matrix[-which(!rownames(mutation_rate_matrix) %in% all_tumors),])
    colnames(mutation_rate_matrix) <- this_MAF$unique_variant_ID_AA
  }


  for(i in 1:nrow(mutation_rate_matrix)){
    for(j in 1:ncol(mutation_rate_matrix)){
      mutation_rate_matrix[i, j] = sum(mutation_rate_nucs[i, unlist(this_MAF$equivalent_aa_muts[j])])
    }
  }

  return(mutation_rate_matrix)
}
