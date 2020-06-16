#' Calculate gene-level selection intensity under an assumption of pairwise gene epistasis
#' @param cesa CESAnalysis object
#' @param genes which genes to calculate effect sizes within
#' @param cores number of cores to use
#' @param optimx_args list of optimization parameters to feed into optimx::opm (see optimx docs); if you use this argument,
#'                    none of the CES default values will be used (opm defaults will be used for parameters you don't provide)
#' @param  return_all_opm_output whether to include full parameter optimization report with output (default false)
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export


ces_gene_epistasis = function(cesa = NULL, genes = character(), cores = 1, optimx_args = ces_gene_epistasis_opm_args(), return_all_opm_output = FALSE)
{
  setkey(cesa@samples, "Unique_Patient_Identifier") # in case dt has forgotten its key
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  # Some optimx::opm optimization methods crash for unclear reasons: Rnmin, nmkb, newuoa, Rtnmin
  # Others should be skipped automatically because they aren't appropriate, but it's
  # necessary to skip them manually: "subplex", "snewtonm", "snewton", "CG", "BFGS","Nelder-Mead", "nlm", "lbfgs"
  working_methods = c("L-BFGS-B", "nlminb", "Rcgmin", "Rvmmin","spg", "bobyqa", "hjkb", "hjn")
  if (is(optimx_args, "list")) {
      if(is.null(optimx_args[["method"]])) {
        warning(paste0("You specified custom optimization parameters with optimx_args, but you didn't include a \"method\" argument ",
                        "(see optimx docs), which means all optimx methods will be used. However, several of these tend to crash this function, ",
                        "and running many methods will substantially slow runtime. The following methods, while not all useful, are believed ",
                        "to not crash: ", paste(working_methods, collapse = ", "), "."))
      } 
  } else {
      stop("optimx_args should of type list (typically, leave out this argument to use CES default parameters)")
  }
  # check for custom optimx arguments
  # can't combine epistasis analysis with multi-stage yet
	if (length(cesa@progressions) > 1) {
	  stop("Epistasis analysis is not compatible yet with multi-stage analyses. You'll have to create a new single-stage analysis.")
	}

  if (length(genes) < 2) {
    stop("Supply at least two genes to analyze.")
  }
  
  if(return_all_opm_output) {
    message(silver("FYI, you can access full parameter optimization output in [CESAnalysis]@advanced$opm_output."))
  }
  
  maf = cesa@maf[Variant_Type == "SNV"]
	genes = unique(genes)
	genes_in_dataset = unique(unlist(maf$genes))
	genes_to_analyze = genes[genes %in% genes_in_dataset]

	recurrent_aac_id = maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut))][, .N, by = aac_id][N > 1, aac_id]
	recurrent_snv_id = maf[! is.na(snv_id), .(snv_id)][, .N, by = "snv_id"][N > 1, snv_id]

	genes_with_recurrent_variants = unique(c(cesa@mutations$snv[recurrent_snv_id, unlist(genes)], cesa@mutations$amino_acid_change[recurrent_aac_id, gene]))
	
	passing_genes = genes_with_recurrent_variants[genes_with_recurrent_variants %in% genes_to_analyze]
	num_genes = length(passing_genes)
  if(num_genes == 0) {
    stop("None of your requested genes have any recurrent variants (see docs for why this is required).", call. = F)
  }
	if(num_genes != length(genes_to_analyze)) {
	  num_missing = length(genes_to_analyze) - num_genes
	  message(paste0(num_genes, " of your requested genes had no recurrent variants, so will not be included (see docs for why this is required)."))
	}
  genes_to_analyze <- passing_genes

    selection_epistasis_results <- t(utils::combn(genes_to_analyze,2))
    selection_epistasis_results <- data.frame(t(selection_epistasis_results),stringsAsFactors=F)
    rownames(selection_epistasis_results) <- c("Variant_1","Variant_2")
    selection_epistasis_results_list <- as.list(selection_epistasis_results)
    selection_results = pbapply::pblapply(X = selection_epistasis_results_list,
                                           FUN = epistasis_gene_level,
                                           cesa=cesa,
                                           optimx_args = optimx_args,
                                           cl = cores)
    cesa@gene_epistasis_results = data.table::rbindlist(lapply(selection_results, function(x) x[[1]]))
    if(return_all_opm_output) {
      opm_output = lapply(selection_results, function(x) x[[3]])
      names(opm_output) = lapply(selection_results, function(x) x[[2]])
      # need to take the rownames (which are the methods used) and put them as a column for conversion to data.table
      opm_output = lapply(opm_output, function(x) {x$method = rownames(x); return(x)})
      opm_dt = data.table::rbindlist(opm_output, idcol = "genes")
      cesa@advanced[["opm_output"]] = opm_dt
    }
    cesa@status[["gene epistasis"]] = "view pairwise gene-level epistatic effect sizes with gene_epistasis_results()"
    return(cesa)
}

epistasis_gene_level = function(cesa, genes_to_analyze, optimx_args) {
  genes_to_analyze = sort(genes_to_analyze) # varies by locale, but hopefully makes them alphabetical
  variant1 = genes_to_analyze[1]
  variant2 = genes_to_analyze[2]

  
  # when including target gene sequencing data, need to throw out any samples that don't have coverage at ALL variant sites in these genes
  # this could be a problem if some of the "exome" samples are actually whole-genome if they haven't been trimmed strictly enough
  
  # Collect mutations present in MAF in the two genes (note that some AACs may be same site, different gene, but okay for coverage check)
  maf = cesa@maf
  maf[, `:=`(in_v1 = variant1 %in% genes, in_v2 = variant2 %in% genes), by = "snv_id"]
  if(maf[in_v1 == T & in_v2 == T, .N] > 0) {
    ## To-do: test behavior
    warning(sprintf("Genes %s and %s having overlapping positions in the MAF data, so the gene pair will be skipped.",
                    variant1, variant2))
    return(NULL)
  }
  
  
  
  # get all amino acid changes in either gene in the data set, then subset to get just the IDs of recurrent ones
  maf = maf[in_v1 | in_v2]
  aac_table = maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut), in_v1, in_v2)]
  aac_table = unique(aac_table[, .(in_v1, in_v2, .N), by = aac_id][N > 1])
  aac_v1 = aac_table[in_v1 == T, aac_id]
  aac_v2 = aac_table[in_v2 == T, aac_id]
  aac = c(aac_v1, aac_v2)
  
  # repeat with noncoding SNVs
  noncoding_table =unique(maf[is.na(assoc_aa_mut), .(in_v1, in_v2, .N), by = "snv_id"][N > 1])
  # remove intergenic records
  noncoding_table[cesa@mutations$snv[, .(intergenic, snv_id)], on = "snv_id", nomatch = NULL][intergenic == F]
  noncoding_v1 = noncoding_table[in_v1 == T, snv_id]
  noncoding_v2 = noncoding_table[in_v2 ==T, snv_id]
  noncoding_snv = c(noncoding_v1, noncoding_v2)
  
  ## restric analysis to samples that have all variants of both genes covered (that is, all remaining recurrent variants)
  eligible_samples = unique(maf$Unique_Patient_Identifier) # these are samples in the v1/v2-only maf
  all_region_sets_in_data = cesa@samples[eligible_samples, unique(covered_regions)]
  eligible_regions = all_region_sets_in_data
  num_aac_sites = length(aac)
  num_noncoding_sites = length(noncoding_snv)
  for (region_set in all_region_sets_in_data) {
    if(cesa@mutations$amino_acid_change[aac, region_set %in% unlist(covered_in), by = "aac_id"][, sum(V1)] < num_aac_sites) {
      eligible_regions = setdiff(all_region_sets_in_data, region_set)
    } else if(cesa@mutations$snv[noncoding_snv, region_set %in% unlist(covered_in), by = "snv_id"][, sum(V1)] < num_noncoding_sites) {
      eligible_regions = setdiff(all_region_sets_in_data, region_set)
    }
  }
  if(length(eligible_regions) == 0) {
    # To-do: test behavior
    warning(sprintf("There are no tumors that have coverage at all recurrent variant sites in genes %s and %s, so the pair must be skipped.",
                    variant1, variant2))
    return(NULL)
  }
  eligible_tumors = cesa@samples[covered_regions %in% eligible_regions, Unique_Patient_Identifier]
  
  # calculate rates across all sites and sum
  baseline_rates_v1 = baseline_mutation_rates(cesa, aac_ids = aac_v1, snv_ids = noncoding_v1, samples = eligible_tumors)
  baseline_rates_v1 = rowSums(as.matrix(baseline_rates_v1[, -"Unique_Patient_Identifier"], rownames = baseline_rates_v1$Unique_Patient_Identifier))
  baseline_rates_v2 = baseline_mutation_rates(cesa, aac_ids = aac_v2, snv_ids = noncoding_v2, samples = eligible_tumors)
  baseline_rates_v2 = rowSums(as.matrix(baseline_rates_v2[, -"Unique_Patient_Identifier"], rownames = baseline_rates_v2$Unique_Patient_Identifier))

  # calculate tumors with and without recurrent mutations in the two genes
  # only the tumors containing a recurrent variant factor into the selection analysis
  ## MAF_input1 and MAF_input2 are already restricted to just eligible tumors
  maf = maf[Unique_Patient_Identifier %in% eligible_tumors]
  all_v1_snv = c(cesa@mutations$amino_acid_change[aac_v1, unlist(all_snv_ids)], noncoding_v1)
  tumors_with_variant1_mutated = maf[snv_id %in% all_v1_snv, unique(Unique_Patient_Identifier)]
  
  all_v2_snv = c(cesa@mutations$amino_acid_change[aac_v2, unlist(all_snv_ids)], noncoding_v2)
  tumors_with_variant2_mutated = maf[snv_id %in% all_v2_snv, unique(Unique_Patient_Identifier)]
  
  tumors_with_both_mutated = intersect(tumors_with_variant1_mutated,tumors_with_variant2_mutated)
  tumors_with_ONLY_variant1_mutated <- tumors_with_variant1_mutated[which(!tumors_with_variant1_mutated %in% tumors_with_variant2_mutated)]
  tumors_with_ONLY_variant2_mutated <- tumors_with_variant2_mutated[which(!tumors_with_variant2_mutated %in% tumors_with_variant1_mutated)]
  tumors_with_neither_mutated <- setdiff(eligible_tumors, c(tumors_with_both_mutated,tumors_with_ONLY_variant1_mutated,tumors_with_ONLY_variant2_mutated))
  
  get_summed_mut_rates = function(tumors) {
    if(length(tumors) == 0) {
      return(NULL)
    }
    rates1 = baseline_rates_v1[tumors]
    rates2 = baseline_rates_v2[tumors]
    return(list(rates1, rates2))
  }
  
  with_just_1 = get_summed_mut_rates(tumors_with_ONLY_variant1_mutated)
  with_just_2 = get_summed_mut_rates(tumors_with_ONLY_variant2_mutated)
  with_both = get_summed_mut_rates(tumors_with_both_mutated)
  with_neither = get_summed_mut_rates(tumors_with_neither_mutated)
  
  par = 1000:1003 # initialized values (rumor has it some methods come up with trivial optimizations when all parameters start the same?)
  args = c(list(par, fn = ml_objective_epistasis_full_gene,with_just_1 = with_just_1,
                with_just_2 = with_just_2, with_both = with_both, with_neither = with_neither), optimx_args)
  
  # suppress warnings because opm complains too much about everything
  opm_output = suppressWarnings(do.call(optimx::opm, args = args))
  
  #TODO: explore the best possible optimization algorithm, fnscale, etc.
  opm_output$value <- -opm_output$value # we did a minimization of the negative, so need to take the negative again to find the maximum
  
  # each row corresponds to the output from one optimization method; sort from best to worst result
  optimization <- opm_output[order(opm_output$value,decreasing = T),]
  params = as.numeric(optimization[1,1:4]) # take best set of parameters

  gene_pair = paste(variant1, variant2, sep=",")
  ces_results = list(gene_1 = variant1, gene_2 = variant2, ces_g1 = params[1], ces_g2 = params[2],
              ces_g1_after_g2 = params[3], ces_g2_after_g1 = params[4], 
              tumors_with_recurrent_muts_only_g1 = length(tumors_with_ONLY_variant1_mutated),
              tumors_with_recurrent_muts_only_g2 = length(tumors_with_ONLY_variant2_mutated),
              tumors_with_recurrent_muts_in_both = length(tumors_with_both_mutated),
              tumors_with_recurrent_muts_in_neither = length(tumors_with_neither_mutated))
  
  return(list(ces_results, gene_pair, optimization))
}


#' Get list of arguments to feed into optimx::opm for parameter optimization during gene-level epistasis
#' @export
#' @keywords internal
ces_gene_epistasis_opm_args = function() {
  # L-BFGS-B is one of the fastest methods and had very good performance in testing
  return(list(gr = "grfwd", lower=1e-3, upper=1e20, method= "L-BFGS-B"))
}
