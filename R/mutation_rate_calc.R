#' Baseline mutation rate calculation
#' 
#' Caculates neutral mutation rates at specific sites based on gene mutation rates and the relative
#' trinucleotide-context-specific SNV mutation rates of each sample
#' 
#' @param cesa CESAnalysis with gene mutation rates and tumor-specific trinucleotide-context-specific mutation rates already calculated
#' @param aac_ids vector of IDs for amino acid change variants
#' @param snv_ids vector of IDs for SNVs
#' @param samples vector of sample IDs (Unique_Patient_Identifier) to include in mutation rate table (defaults to all samples)
#' @return a data table of mutation rates with one column per variant, and a Unique_Patient_Identifier column identifying each row
#' @export
baseline_mutation_rates = function(cesa, aac_ids = NULL, snv_ids = NULL, samples = NULL) {
  
  # Let user specify a subset of samples to calculate rates (or, by default, use all samples)
  if(! is.null(samples)) {
    if (! is.character(samples)) {
      stop("Samples should be character vector", call. = F)
    }
    samples = unique(samples)
    subsetted_samples = cesa@samples[samples]
    if(subsetted_samples[, .N] != length(samples)) {
      stop("One or more of the requested samples isn't in the CESAnalysis", call. = F)
    }
    samples = subsetted_samples
  } else {
    samples = cesa@samples
  }
  mutations = cesa@mutations
  
  
  # produce a table with all pairwise combinations of Unique_Patient_Identifier and relevant genes
  # relevant genes are those associated with one of the AACs/SNVs of interest
  
  relevant_genes = union(mutations$amino_acid_change[aac_id %in% aac_ids, gene], mutations$snv[snv_id %in% snv_ids, unlist(genes)])
  sample_gene_rates = as.data.table(expand.grid(gene = relevant_genes, Unique_Patient_Identifier = samples$Unique_Patient_Identifier, 
                                                stringsAsFactors = F),key = "Unique_Patient_Identifier")
  
  # add gene mutation rates to the table by using @mutrates and the progressions of each samples
  sample_gene_rates = sample_gene_rates[samples[, .(Unique_Patient_Identifier, progression_name)]]
  melted_mutrates = melt.data.table(cesa@mutrates[gene %in% relevant_genes], id.vars = c("gene"))
  setnames(melted_mutrates, c("variable", "value"), c("progression_name", "raw_rate"))
  sample_gene_rates = melted_mutrates[sample_gene_rates, , on = c("gene", "progression_name")]
  
  # add a list column with the trinuc composition of each gene (composition is a 96-length numeric, deconstructSigs order)
  gene_trinuc_comp = get_genome_data(cesa, "gene_trinuc_comp")
  sample_gene_rates[, trinuc_comp := lapply(gene, function(x) gene_trinuc_comp[[x]])]
  trinuc_mat = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  
  # dot product of trinuc comp and patient's expected relative trinuc rates yields the denominator
  # for site-specific mutation rate calculation; numerator is raw gene rate multipled by the patient's relative rate for the site's trinuc context
  # (this last value gets multipled in by get_baseline functions below)
  sample_gene_rates[, aggregate_rate := raw_rate / sum(unlist(trinuc_comp) * trinuc_mat[Unique_Patient_Identifier, ]), by = c("gene", "Unique_Patient_Identifier")]
  setkey(sample_gene_rates, "gene")
  

  

  
  if(length(aac_ids > 0)) {
    trinuc_mut_by_aac = mutations$amino_acid_change[aac_id %in% aac_ids, .(trinuc_mut = list(mutations$snv[unlist(all_snv_ids), trinuc_mut])), by = "aac_id"]
    gene_by_aac = mutations$amino_acid_change[aac_ids, gene]
    aac_genes = unique(gene_by_aac)
    
    # build vector agg_rates (see above comment) corresponding to each aac
    tmp = sample_gene_rates[aac_genes, .(list(aggregate_rate)), by = c("gene")]
    agg_rates_by_gene = tmp$V1
    names(agg_rates_by_gene) = tmp$gene
    agg_rates_by_gene = agg_rates_by_gene[gene_by_aac]
    agg_rates_by_gene = lapply(agg_rates_by_gene, unlist)
    
    # for each sample, multiply agg_rate by relative trinuc rate of the context for all AAC sites
    get_baseline_aac = function(aac, agg_rates) {
      trinuc_mut = trinuc_mut_by_aac[aac, unlist(trinuc_mut)]
      sample_rates = rowSums(trinuc_mat[, trinuc_mut, drop = F])
      return(agg_rates * sample_rates)
    }
    aac_rate_list = mapply(get_baseline_aac, aac_ids, agg_rates_by_gene, SIMPLIFY = F)
  } else {
    aac_rate_list = NULL
  }
  

  # repeat with SNVs (slightly different handling since SNVs can have more than one gene associated)
  trinuc_mut_by_snv = mutations$snv[snv_ids, trinuc_mut, keyby="snv_id"]
  gene_by_snv = mutations$snv[snv_ids, genes, by = "snv_id"]
  snv_genes = unique(unlist(gene_by_snv$genes))
  tmp = sample_gene_rates[snv_genes, .(list(aggregate_rate)), by = c("gene")]
  agg_rates_by_gene = tmp$V1
  names(agg_rates_by_gene) = tmp$gene
  
  get_baseline_snv = function(snv, genes) {
    trinuc_mut = trinuc_mut_by_snv[snv, trinuc_mut]
    sample_rates = trinuc_mat[, trinuc_mut, drop = F]
    if(length(genes) > 1) {
      agg_rates = rowMeans(simplify2array(agg_rates_by_gene[genes]))
    } else {
      agg_rates = agg_rates_by_gene[[genes]]
    }
    return(agg_rates * sample_rates)
  }
  snv_rate_list = mapply(get_baseline_snv, snv_ids, gene_by_snv$genes, SIMPLIFY = F)
  
  
  # Combine all rates and create table
  # Need to raise number of allocated columns high enough to add all columns at once (possibly thousands)
  # Watch out for changes to setalloccol in future versions of data.table
  all_rates_list = c(aac_rate_list, snv_rate_list)
  baseline_rates = samples[, .(Unique_Patient_Identifier)] # pre-keyed
  columns_needed = length(all_rates_list) + 1
  if(truelength(baseline_rates) < columns_needed) {
    invisible(suppressWarnings(setalloccol(baseline_rates, columns_needed)))
  }
  baseline_rates[, (names(all_rates_list)) := all_rates_list]
  return(baseline_rates)
}
