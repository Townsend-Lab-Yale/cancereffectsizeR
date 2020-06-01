#' Mutation rate calculation
#'
#' This function calculates the mutation rates of all the unique varaints within a gene given gene-level mutation rate and tumor-specific trinucleotide contexts.
#'
#' @param this_MAF The subset MAF (just the gene to be analyzed) file to extract mutational data from
#' @param gene current gene
#' @param trinuc_proportion_matrix matrix constructed from deconstructSigs output, containing proportion of each trinucleotide mutated in each tumor
#' @param samples data table with sample info
#' @param gene_mut_rate mutation rate at the gene-level
#' @param gene_trinuc_comp genome/transcriptome-specific gene trinuc composition
#'
#' @return
#' @export
#' @keywords internal
mutation_rate_calc <- function(this_MAF,
                               gene,
                               gene_mut_rate,
                               trinuc_proportion_matrix,
                               gene_trinuc_comp,
                               samples
                               ){

  # trinuc_proportion_matrix: rows = samples, columns = relative frequency of trinucleotide-context-specific mutation (sums to 1)
  # mutation_rate_nucs: rows = samples, columns = gene-specific relative frequency of trinuclotide-context-specific
  ##                    mutations, calculated by normalizing trinuc_proportion_matrix by gene content
  mutation_rate_nucs <- matrix(nrow=nrow(trinuc_proportion_matrix),ncol=ncol(trinuc_proportion_matrix),data = NA)
  rownames(mutation_rate_nucs) <- rownames(trinuc_proportion_matrix)
  colnames(mutation_rate_nucs) <- colnames(trinuc_proportion_matrix)

  trinuc_comp = gene_trinuc_comp[[gene]]
  calc_normalizers = function(x) {
    return(sum(trinuc_comp * x))
  }
  normalizers = apply(trinuc_proportion_matrix, 1, calc_normalizers)

  current_gene = gene
  gene_rates_by_state = as.numeric(gene_mut_rate[gene == current_gene, -"gene"])

  gene_rates_by_sample = gene_rates_by_state[samples[rownames(mutation_rate_nucs), progression_index]]
  mutation_rate_nucs = (trinuc_proportion_matrix / normalizers) * gene_rates_by_sample
  

  # mutation_rate_nucs is now the rate of each trinucleotide in each tumor for this gene
  # need to find unique variants and then rates
  this_MAF <- this_MAF[!duplicated(this_MAF[,c("nt_mut_id")]),]


  # mutation rate matrix: rows = tumors, columns = expected relative rates of amino acid changes in gene for tumor
  ## this is determined by summing over the rates for different nucleotides mutations that create same amino acid changes
  mutation_rate_matrix <- matrix(nrow=nrow(trinuc_proportion_matrix), ncol=nrow(this_MAF))
  rownames(mutation_rate_matrix) <- rownames(trinuc_proportion_matrix)
  colnames(mutation_rate_matrix) <- this_MAF$nt_mut_id


  for(i in 1:nrow(mutation_rate_matrix)){
    for(j in 1:ncol(mutation_rate_matrix)){
      mutation_rate_matrix[i, j] = sum(mutation_rate_nucs[i, unlist(this_MAF$equivalent_aa_muts[j])])
    }
  }

  return(mutation_rate_matrix)
}
