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


#' baseline_mutation_rates
#' 
#' @export
baseline_mutation_rates = function(cesa, aac_ids = NULL, snv_ids = NULL) {
  
  samples = cesa@samples
  mutations = cesa@mutations
  relevant_genes = union(mutations$amino_acid_change[aac_id %in% aac_ids, gene], mutations$snv[snv_id %in% snv_ids, unlist(genes)])
  
  
  sample_gene_rates = as.data.table(expand.grid(gene = relevant_genes, Unique_Patient_Identifier = samples$Unique_Patient_Identifier, 
                                                stringsAsFactors = F),key = "Unique_Patient_Identifier")
  
  sample_gene_rates = sample_gene_rates[samples[, .(Unique_Patient_Identifier, progression_name)]]
  
  melted_mutrates = melt.data.table(cesa@mutrates[gene %in% relevant_genes], id.vars = c("gene"))
  setnames(melted_mutrates, c("variable", "value"), c("progression_name", "raw_rate"))
  
  sample_gene_rates = melted_mutrates[sample_gene_rates, , on = c("gene", "progression_name")]
  
  gene_trinuc_comp = get_genome_data(cesa, "gene_trinuc_comp")
  sample_gene_rates[, trinuc_comp := lapply(gene, function(x) gene_trinuc_comp[[x]])]
  trinuc_mat = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  sample_gene_rates[, aggregate_rate := raw_rate / sum(unlist(trinuc_comp) * trinuc_mat[Unique_Patient_Identifier, ]), by = c("gene", "Unique_Patient_Identifier")]
  
  
  baseline_rates = samples[, .(Unique_Patient_Identifier)] # pre-keyed
  trinuc_mut_by_aac = mutations$amino_acid_change[aac_id %in% aac_ids, .(trinuc_mut = list(mutations$snv[unlist(all_snv_ids), trinuc_mut])), by = "aac_id"]
  setkey(sample_gene_rates, "gene")
  
  gene_by_aac = mutations$amino_acid_change[aac_ids, gene]
  aac_genes = unique(gene_by_aac)
  tmp = sample_gene_rates[aac_genes, .(list(aggregate_rate)), by = c("gene"), allow.cartesian = T, nomatch = NULL]
  agg_rates_by_gene = tmp$V1
  names(agg_rates_by_gene) = tmp$gene
  agg_rates_by_gene = agg_rates_by_gene[gene_by_aac]
  agg_rates_by_gene = lapply(agg_rates_by_gene, unlist)
  
  get_baseline = function(aac, agg_rates) {
    trinuc_mut = trinuc_mut_by_aac[aac, unlist(trinuc_mut)]
    rs = rowSums(trinuc_mat[, trinuc_mut, drop = F])
    return(agg_rates * rs)
  }
  tmp = mapply(get_baseline, aac_ids, agg_rates_by_gene, SIMPLIFY = F)
  # Watch out for changes to setalloccol in future versions of data.table
  # Need to raise number of allocated columns high enough to add all columns at once
  if(length(tmp) + length(baseline_rates) > truelength(baseline_rates)) {
    invisible(suppressWarnings(setalloccol(baseline_rates, length(tmp))))
  }
  baseline_rates[, (aac_ids) := tmp]
  
  for (snv in snv_ids) {
    # occasionally more than one gene annotated; in this case take the average
    curr_genes = mutations$snv[snv, unlist(genes)]
    trinuc_mut = mutations$snv[snv, trinuc_mut]
    tmp = trinuc_mat[, trinuc_mut]
    #tmp = dt[, sum(as.numeric(.SD)), .SDcols = trinuc_mut, by = "Unique_Patient_Identifier"]
    tmp =  sample_gene_rates[gene %in% curr_genes][ , V1 := tmp, by = "gene"]
    tmp[, per_gene_rate := tmp$aggregate_rate * tmp$V1]
    averaged_over_genes = tmp[, mean(per_gene_rate), by = "Unique_Patient_Identifier"]
    baseline_rates[, (snv) := averaged_over_genes$V1]
  }
  return(baseline_rates)
}
