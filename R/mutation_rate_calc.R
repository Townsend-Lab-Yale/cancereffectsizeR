#' Mutation rate calculation
#'
#' This function calculates the mutation rates of all the unique varaints within a gene given gene-level mutation rate and tumor-specific trinucleotide contexts.
#'
#' @param this_MAF The subset MAF (just the gene to be analyzed) file to extract mutational data from
#' @param gene current gene
#' @param trinuc_proportion_matrix matrix constructed from deconstructSigs output, containing proportion of each trinucleotide mutated in each tumor
#' @param gene_trinuc_comp list containing matrices of counts of each trinuc in every gene
#' @param gene_refcds information from RefCDS data set for the current gene
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
                               gene_refcds,
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

      if(this_MAF$next_to_splice[j]){

        if(this_MAF$is_coding[j]){
          mutation_rate_matrix[i,j] <- sum(mutation_rate_nucs[i,cancereffectsizeR::mutation_finder(RefCDS_instance = gene_refcds, MAF_input_row = this_MAF[j,])])
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

  unsure_genes_vec <- this_MAF$unsure_gene_name



  return(list(mutation_rate_matrix=mutation_rate_matrix,
              unsure_genes_vec=unsure_genes_vec,
              unique_variant_ID_vec = this_MAF$unique_variant_ID))

}
