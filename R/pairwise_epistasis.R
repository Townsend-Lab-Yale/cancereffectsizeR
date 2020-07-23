


#' Gene-level epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise gene-level epistasis. 
#' 
#' @param cesa CESAnalysis object
#' @param genes vector of gene names; SIs will be calculated for all gene pairs 
#' @param cores number of cores for parallel processing of gene pairs
#' @param optimx_args list of optimization parameters to feed into optimx::opm (see optimx docs); if you use this argument,
#'                    none of the CES default values will be used (opm defaults will be used for parameters you don't provide)
#' @param  return_all_opm_output whether to include full parameter optimization report with output (default false)
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export
ces_gene_epistasis = function(cesa = NULL, genes = character(), cores = 1, optimx_args = ces_epistasis_opm_args(), return_all_opm_output = FALSE)
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
	  message(paste0(num_missing, " of your requested genes had no recurrent variants, so will not be included (see docs for why this is required)."))
	}
  genes_to_analyze <- passing_genes
  if (length(genes_to_analyze) == 1) {
    stop(sprintf("Only 1 requested gene (%s) has recurrent variants, so epistasis can't be run.", genes_to_analyze), call. = F)
  }
  gene_pairs <- utils::combn(sort(genes_to_analyze),2, simplify = F)
  selection_results = pbapply::pblapply(X = gene_pairs,
                                         FUN = pairwise_gene_epistasis,
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
  return(cesa)
}

#' Variant-level pairwise epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise epistasis between pairs of variants
#' 
#' @param cesa CESAnalysis
#' @param variants list where each element is a character-type pair of variant IDs (either amino-acid-change (coding) variant IDs or SNV IDs)
#' @param cores number of cores for parallel processing of gene pairs
#' @return a data.table with pairwise-epistatic selection intensities and variant frequency info
#' @export
ces_epistasis = function(cesa = NULL, variants = NULL, cores = 1) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should specify a CESAnalysis object", call. = F)
  }
  
  # can't combine epistasis analysis with multi-stage yet
  if (length(cesa@progressions) > 1) {
    stop("Epistasis analysis is not compatible yet with multi-stage analyses. You'll have to create a new single-stage analysis.", call. = F)
  }
  
  if (is.null(cesa@mutations$snv)) {
    stop("The CESAnalysis must have mutation annotations.", call. = F)
  } 
  
  if(! is(variants, "list")) {
    stop("variants should be of list type")
  }
  if (length(variants) < 1) {
    stop("Variants list is empty.", call. = F)
  }
  
  setkey(cesa@samples, "Unique_Patient_Identifier") # in case dt has forgotten its key
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  # validate that all variant IDs appear in mutation annotations
  for (pair in variants) {
    if (! is(pair, "character") || length(pair) != 2) {
      stop("All elements of variants list must be character vectors of length 2 (for two variant IDs)", call. = F)
    }
    for (variant in pair) {
      if (cesa@mutations$snv[variant, .N, nomatch = NULL] + cesa@mutations$amino_acid_change[variant, .N, nomatch = NULL] == 0) {
        stop(sprintf("Requested variant %s is not in the mutation annotation tables.", variant), call. = F)
      }
    }
  }
  selection_results = rbindlist(pbapply::pblapply(X = variants, FUN = pairwise_variant_epistasis, cesa=cesa, cl = cores))
  return(selection_results)
  
}

#' Calculate SIs at variant level under pairwise epistasis model
#' @keywords internal
pairwise_variant_epistasis = function(cesa, variant_pair, variant_types, optimx_args = ces_epistasis_opm_args()) {
  v1 = variant_pair[1]
  v2 = variant_pair[2]
  
  v1_coverage = c(cesa$mutations$amino_acid_change[v1, unlist(covered_in), nomatch = NULL], 
                  cesa$mutations$snv[v1, unlist(covered_in), nomatch = NULL])
  v2_coverage = c(cesa$mutations$amino_acid_change[v2, unlist(covered_in), nomatch = NULL], 
                  cesa$mutations$snv[v2, unlist(covered_in), nomatch = NULL])
  
  joint_coverage = intersect(v1_coverage, v2_coverage)
  
  eligible_tumors = cesa@samples[covered_regions %in% joint_coverage, Unique_Patient_Identifier]
  
  if (length(eligible_tumors) == 0) {
    warning(sprintf("No samples have coverage at both %s and %s, so this variant pair had to be skipped.", v1, v2), immediate. = T, call. = F)
  }
  
  all_rates = baseline_mutation_rates(cesa = cesa, variant_ids = c(v1, v2), samples = eligible_tumors)
  setcolorder(all_rates, c("Unique_Patient_Identifier", v1, v2))
  setkey(all_rates, "Unique_Patient_Identifier")
  
  tumors_with_v1 = cesa@maf[, v1 %in% c(snv_id, assoc_aa_mut), by = "Unique_Patient_Identifier"][V1 == T, Unique_Patient_Identifier]
  tumors_with_v2 = cesa@maf[, v2 %in% c(snv_id, assoc_aa_mut), by = "Unique_Patient_Identifier"][V1 == T, Unique_Patient_Identifier]
  tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
  tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
  tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
  tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
  
  # 2-item lists: first item is baseline rates for v1; second for v2
  with_just_1 = as.list(all_rates[tumors_just_v1])[2:3]
  with_just_2 = as.list(all_rates[tumors_just_v1])[2:3]
  with_both = as.list(all_rates[tumors_with_both])[2:3]
  with_neither = as.list(all_rates[tumors_with_neither])[2:3]
  
  par = 1000:1003 # initialized values
  args = c(list(par, fn = ml_objective_pairwise_epistasis,with_just_1 = with_just_1,
                with_just_2 = with_just_2, with_both = with_both, with_neither = with_neither), optimx_args)
  
  # suppress warnings because opm complains too much about everything
  opm_output = suppressWarnings(do.call(optimx::opm, args = args))
  opm_output$value <- -opm_output$value # we did a minimization of the negative, so need to take the negative again to find the maximum
  
  # each row corresponds to the output from one optimization method; sort from best to worst result
  optimization <- opm_output[order(opm_output$value,decreasing = T),]
  params = as.numeric(optimization[1,1:4]) # take best set of parameters
  
  ces_results = list(variant1 = v1, variant2 = v2, ces_v1 = params[1], ces_v2 = params[2],
                     ces_v1_after_v2 = params[3], ces_v2_after_v1 = params[4], 
                     tumors_just_v1 = length(tumors_just_v1),
                     tumors_just_v2 = length(tumors_just_v2),
                     tumors_with_both = length(tumors_with_both),
                     tumors_with_neither = length(tumors_with_neither))
  return(ces_results)
}


#' Calculate SIs at gene level under pairwise epistasis model
#' @keywords internal
pairwise_gene_epistasis = function(cesa, genes, optimx_args) {
  gene1 = genes[1]
  gene2 = genes[2]

  # when including target gene sequencing data, we need to throw out any samples
  # that don't have coverage at ALL recurrent variant sites in these genesthis
  # could be a problem if some of the "exome" samples are actually whole-genome
  # if they haven't been trimmed strictly enough
  
  # Collect mutations present in MAF in the two genes (note that some AACs may
  # be same site, different gene, but okay for coverage check)
  maf = cesa@maf
  maf[, `:=`(in_g1 = gene1 %in% genes, in_g2 = gene2 %in% genes), by = "snv_id"]
  if(maf[in_g1 == T & in_g2 == T, .N] > 0) {
    ## To-do: test behavior
    warning(sprintf("Genes %s and %s having overlapping positions in the MAF data, so the gene pair will be skipped.",
                    gene1, gene2))
    return(NULL)
  }
  
  # get all amino acid changes in either gene in the data set, then subset to get just the IDs of recurrent ones
  maf = maf[in_g1 | in_g2]
  aac_table = maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut), in_g1, in_g2)]
  aac_table = unique(aac_table[, .(in_g1, in_g2, .N), by = aac_id][N > 1])
  aac_v1 = aac_table[in_g1 == T, aac_id]
  aac_v2 = aac_table[in_g2 == T, aac_id]
  aac = c(aac_v1, aac_v2)
  
  # repeat with noncoding SNVs
  noncoding_table = unique(maf[is.na(assoc_aa_mut), .(in_g1, in_g2, .N), by = "snv_id"][N > 1])
  # remove intergenic records
  noncoding_table[cesa@mutations$snv[, .(intergenic, snv_id)], on = "snv_id", nomatch = NULL][intergenic == F]
  noncoding_v1 = noncoding_table[in_g1 == T, snv_id]
  noncoding_v2 = noncoding_table[in_g2 ==T, snv_id]
  noncoding_snv = c(noncoding_v1, noncoding_v2)
  
  ## restrict analysis to samples that have all variants of both genes covered (that is, all remaining recurrent variants)
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
                    gene1, gene2))
    return(NULL)
  }
  eligible_tumors = cesa@samples[covered_regions %in% eligible_regions, Unique_Patient_Identifier]
  
  # calculate rates across all sites and sum
  baseline_rates_g1 = baseline_mutation_rates(cesa, aac_ids = aac_v1, snv_ids = noncoding_v1, samples = eligible_tumors)
  baseline_rates_g1 = rowSums(as.matrix(baseline_rates_g1[, -"Unique_Patient_Identifier"], rownames = baseline_rates_g1$Unique_Patient_Identifier))
  baseline_rates_g2 = baseline_mutation_rates(cesa, aac_ids = aac_v2, snv_ids = noncoding_v2, samples = eligible_tumors)
  baseline_rates_g2 = rowSums(as.matrix(baseline_rates_g2[, -"Unique_Patient_Identifier"], rownames = baseline_rates_g2$Unique_Patient_Identifier))

  # calculate tumors with and without recurrent mutations in the two genes
  # only the tumors containing a recurrent variant factor into the selection analysis
  ## MAF_input1 and MAF_input2 are already restricted to just eligible tumors
  maf = maf[Unique_Patient_Identifier %in% eligible_tumors]
  all_g1_snv = c(cesa@mutations$amino_acid_change[aac_v1, unlist(all_snv_ids)], noncoding_v1)
  tumors_with_gene1_mutated = maf[snv_id %in% all_g1_snv, unique(Unique_Patient_Identifier)]
  
  all_g2_snv = c(cesa@mutations$amino_acid_change[aac_v2, unlist(all_snv_ids)], noncoding_v2)
  tumors_with_gene2_mutated = maf[snv_id %in% all_g2_snv, unique(Unique_Patient_Identifier)]
  
  tumors_with_both_mutated = intersect(tumors_with_gene1_mutated,tumors_with_gene2_mutated)
  tumors_with_ONLY_gene1_mutated <- tumors_with_gene1_mutated[which(!tumors_with_gene1_mutated %in% tumors_with_gene2_mutated)]
  tumors_with_ONLY_gene2_mutated <- tumors_with_gene2_mutated[which(!tumors_with_gene2_mutated %in% tumors_with_gene1_mutated)]
  tumors_with_neither_mutated <- setdiff(eligible_tumors, c(tumors_with_both_mutated,tumors_with_ONLY_gene1_mutated,tumors_with_ONLY_gene2_mutated))
  
  get_summed_mut_rates = function(tumors) {
    if(length(tumors) == 0) {
      return(NULL)
    }
    rates1 = unname(baseline_rates_g1[tumors])
    rates2 = unname(baseline_rates_g2[tumors])
    return(list(rates1, rates2))
  }
  
  with_just_1 = get_summed_mut_rates(tumors_with_ONLY_gene1_mutated)
  with_just_2 = get_summed_mut_rates(tumors_with_ONLY_gene2_mutated)
  with_both = get_summed_mut_rates(tumors_with_both_mutated)
  with_neither = get_summed_mut_rates(tumors_with_neither_mutated)
  
  par = 1000:1003 # initialized values (rumor has it some methods come up with trivial optimizations when all parameters start the same?)
  args = c(list(par, fn = ml_objective_pairwise_epistasis,with_just_1 = with_just_1,
                with_just_2 = with_just_2, with_both = with_both, with_neither = with_neither), optimx_args)
  
  # suppress warnings because opm complains too much about everything
  opm_output = suppressWarnings(do.call(optimx::opm, args = args))
  opm_output$value <- -opm_output$value # we did a minimization of the negative, so need to take the negative again to find the maximum
  
  # each row corresponds to the output from one optimization method; sort from best to worst result
  optimization <- opm_output[order(opm_output$value,decreasing = T),]
  params = as.numeric(optimization[1,1:4]) # take best set of parameters

  gene_pair = paste(gene1, gene2, sep=",")
  ces_results = list(gene_1 = gene1, gene_2 = gene2, ces_g1 = params[1], ces_g2 = params[2],
              ces_g1_after_g2 = params[3], ces_g2_after_g1 = params[4], 
              tumors_with_recurrent_muts_only_g1 = length(tumors_with_ONLY_gene1_mutated),
              tumors_with_recurrent_muts_only_g2 = length(tumors_with_ONLY_gene2_mutated),
              tumors_with_recurrent_muts_in_both = length(tumors_with_both_mutated),
              tumors_with_recurrent_muts_in_neither = length(tumors_with_neither_mutated))
  
  return(list(ces_results, gene_pair, optimization))
}


#' Get recommended optimx optimization parameters for epistasis functions
#' @export
#' @keywords internal
ces_epistasis_opm_args = function() {
  # L-BFGS-B is one of the fastest methods and had very good performance in testing
  return(list(gr = "grfwd", lower=1e-3, upper=1e20, method= "L-BFGS-B"))
}
