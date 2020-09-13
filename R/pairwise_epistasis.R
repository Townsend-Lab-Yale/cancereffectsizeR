#' Gene-level epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise gene-level epistasis.
#' More specifically, selection at the gene level is assumed to act through all recurrent
#' variants (variants that appear in just one sample are ignored), with selection
#' intensity of each gene assumed to vary based on mutational status of the other gene.
#' Returns a table of results, \strong{NOT} a CESAnalysis, so be careful not to overwrite
#' your CESAnalysis.
#' 
#' Only samples that have coverage at all recurrent sites can be included in analysis
#' since samples lacking full coverage may or may not have mutations at the uncovered
#' sites
#' 
#' 
#' @param cesa CESAnalysis object
#' @param genes vector of gene names; SIs will be calculated for all gene pairs 
#' @param cores number of cores for parallel processing of gene pairs
#' @param optimx_args list of optimization parameters to feed into optimx::opm (see optimx docs); if you use this argument,
#'                    none of the CES default values will be used (opm defaults will be used for parameters you don't provide)
#' @param include_opm_report dev option to include full parameter optimization report with output (default false)
#' @return a table of gene-level epistasis results
#' @export
ces_gene_epistasis = function(cesa = NULL, genes = character(), cores = 1, optimx_args = ces_epistasis_opm_args(), include_opm_report = FALSE)
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
	if (length(cesa@groups) > 1) {
	  stop("Epistasis analysis is not compatible yet with multi-stage analyses. You'll have to create a new single-stage analysis.")
	}

  if (length(genes) < 2) {
    stop("Supply at least two genes to analyze.")
  }
  
  maf = cesa@maf[variant_type == "snv"]
	
	
	# subset input genes to those with data
  genes = unique(genes)
	genes_to_analyze = intersect(unique(unlist(maf$genes)), genes)
	recurrent_variants = select_variants(cesa = cesa, genes = genes_to_analyze, min_freq = 2, collapse_lists = F)
	recurrent_snv_id = recurrent_variants[variant_type == "snv", variant_id]
	recurrent_aac_id = recurrent_variants[variant_type == "aac", variant_id]
	
	# use "all_genes" column because SNVs might have multiple genes, and only one appears in the "gene" column
	genes_with_recurrent_variants = unique(unlist(recurrent_variants$all_genes)) 

	genes_to_analyze = intersect(genes_with_recurrent_variants, genes_to_analyze)
	num_passing_genes = length(genes_to_analyze)
  if(num_passing_genes == 0) {
    stop("None of your requested genes have any recurrent variants (see docs for why this is required).", call. = F)
  }
	if(num_passing_genes != length(genes)) {
	  num_missing = length(genes) - num_passing_genes
	  pretty_message(paste0(num_missing, " of your requested genes had no recurrent variants, so will not be included (see docs for why this is required)."))
	}

	if (length(genes_to_analyze) == 1) {
	  stop(sprintf("Only 1 requested gene (%s) has recurrent variants, so epistasis can't be run.", genes_to_analyze), call. = F)
	}
	
	gr_table = as.data.table(get_ref_data(cesa@ref_data_dir, "gr_genes"))
	setkey(gr_table, "names")
	passing_gene_ranges = gr_table[genes_to_analyze]
	setkey(passing_gene_ranges, "seqnames", "start", "end")
	
  gene_pairs <- utils::combn(sort(genes_to_analyze), 2, simplify = F)
  bad_pairs = numeric()
  for (i in 1:length(gene_pairs)) {
    if (foverlaps(passing_gene_ranges[names == gene_pairs[[i]][1]], 
              passing_gene_ranges[names == gene_pairs[[i]][2]], nomatch = NULL)[, .N] > 0) {
      bad_pairs = c(bad_pairs, i)
    }
  }
  if (length(bad_pairs) > 0) {
    msg = paste0("The following pairs of genes were skipped because their ranges overlap in reference data:\n",
                paste(sapply(gene_pairs[bad_pairs], function(x) paste(x, collapse = "/")), collapse = ", "))
    pretty_message(msg, black = F)
    gene_pairs = gene_pairs[- bad_pairs]
    if (length(gene_pairs) == 0) {
      stop("There are no remaining gene pairs to test for epistasis.")
    }
  }
  
  selection_results = pbapply::pblapply(X = gene_pairs,
                                         FUN = pairwise_gene_epistasis,
                                         cesa=cesa,
                                         optimx_args = optimx_args,
                                         cl = cores)
  results = data.table::rbindlist(lapply(selection_results, function(x) x[[1]]))
  if(include_opm_report) {
    opm_output = lapply(selection_results, function(x) x[[3]])
    names(opm_output) = lapply(selection_results, function(x) x[[2]])
    # need to take the rownames (which are the methods used) and put them as a column for conversion to data.table
    opm_output = lapply(opm_output, function(x) {x$method = rownames(x); return(x)})
    opm_dt = data.table::rbindlist(opm_output, idcol = "genes")
    results = list(results = results, opm_report = opm_dt)
  }
  return(results)
}

#' Variant-level pairwise epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise epistasis between pairs of variants.
#' Return a table, \strong{NOT} a CESAnalysis, so be careful not to overwrite your analysis.
#' 
#' @param cesa CESAnalysis
#' @param variants list where each element is a character-type pair of variant IDs (either
#'   amino-acid-change (coding) variant IDs or SNV IDs)
#' @param cores number of cores for parallel processing of gene pairs
#' @return a data table with pairwise-epistatic selection intensities and variant
#'   frequencies for in tumors that have coverage at both variants in each pair
#' @export
ces_epistasis = function(cesa = NULL, variants = NULL, cores = 1) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should specify a CESAnalysis object", call. = F)
  }
  
  # can't combine epistasis analysis with multi-stage yet
  if (length(cesa@groups) > 1) {
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
  validated_variants = character()
  for (i in 1:length(variants)) {
    pair = variants[[i]]
    if (! is(pair, "character") || length(pair) != 2) {
      stop("All elements of variants list must be character vectors of length 2 (for two variant IDs)")
    }
    # this is to allow short IDs in the variants list
    variants[[i]] = unlist(suppressMessages(select_variants(cesa, variant_ids = pair, ids_only = T)), use.names = F)
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
  
  # Samples have to have v1 and v2 coverage (and samples with covered_regions == "genome" always have coverage)
  joint_coverage = c("genome", intersect(v1_coverage, v2_coverage))
  
  eligible_tumors = cesa@samples[covered_regions %in% joint_coverage, Unique_Patient_Identifier]
  
  if (length(eligible_tumors) == 0) {
    warning(sprintf("No samples have coverage at both %s and %s, so this variant pair had to be skipped.", v1, v2), immediate. = T, call. = F)
  }
  
  all_rates = baseline_mutation_rates(cesa = cesa, variant_ids = c(v1, v2), samples = eligible_tumors)
  setcolorder(all_rates, c("Unique_Patient_Identifier", v1, v2))
  setkey(all_rates, "Unique_Patient_Identifier")
  
  tumors_with_v1 = cesa@maf[, v1 %in% c(variant_id, assoc_aac), by = "Unique_Patient_Identifier"][V1 == T, Unique_Patient_Identifier]
  tumors_with_v2 = cesa@maf[, v2 %in% c(variant_id, assoc_aac), by = "Unique_Patient_Identifier"][V1 == T, Unique_Patient_Identifier]
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
                     covered_tumors_just_v1 = length(tumors_just_v1),
                     covered_tumors_just_v2 = length(tumors_just_v2),
                     covered_tumors_with_both = length(tumors_with_both),
                     covered_tumors_with_neither = length(tumors_with_neither))
  return(ces_results)
}


#' Calculate SIs at gene level under pairwise epistasis model
#' 
#' The genes are assumed to not overlap in any ranges (the calling
#' function checks for this).
#' @keywords internal
pairwise_gene_epistasis = function(cesa, genes, optimx_args) {
  gene1 = genes[1]
  gene2 = genes[2]

  # When including targeted gene sequencing data, we need to throw out any samples
  # that don't have coverage at ALL recurrent variant sites in these genes. This
  # could be a problem if some of the "exome" samples are actually whole-genome
  # if they haven't been trimmed strictly enough, especially if some of those
  # sites are low-complexity and filled with calling errors

  # Get the IDs of recurrent variants in the genes
  # Note that intergenic SNVs won't be included, even if they're close to a gene
  # Doing each gene individually to facilitate finding t
  rec_variants_g1 = select_variants(cesa, genes = gene1, min_freq = 2, collapse_lists = F)
  aac_g1 = rec_variants_g1[variant_type == "aac", variant_id]
  noncoding_g1 = rec_variants_g1[variant_type == "snv", variant_id]
  all_g1_snv = unique(c(noncoding_g1, na.omit(unlist(rec_variants_g1$constituent_snvs))))
  
  rec_variants_g2 = select_variants(cesa, genes = gene2, min_freq = 2, collapse_lists = F)
  aac_g2 = rec_variants_g2[variant_type == "aac", variant_id]
  noncoding_g2 = rec_variants_g2[variant_type == "snv", variant_id]
  all_g2_snv = unique(c(noncoding_g1, na.omit(unlist(rec_variants_g2$constituent_snvs))))
  
  rec_variants = rbind(rec_variants_g1, rec_variants_g2)
  aac = c(aac_g1, aac_g2)
  noncoding_snv = c(noncoding_g1, noncoding_g2)
  
  # Can only use samples that have coverage at all sites
  num_covered_by_panel = table(unlist(rec_variants$covered_in))
  eligible_regions = c(names(which(num_covered_by_panel == rec_variants[, .N])), "genome")
  eligible_tumors = cesa@samples[covered_regions %in% eligible_regions, Unique_Patient_Identifier]
  
  if(length(eligible_tumors) == 0) {
    no_tumors_with_coverage_warning = sprintf("There are no tumors that have coverage at all recurrent variant sites in genes %s and %s, so the pair must be skipped.",
                                              gene1, gene2)
    warning(no_tumors_with_coverage_warning)
    return(NULL)
  }

  # calculate rates across all sites and sum
  baseline_rates_g1 = baseline_mutation_rates(cesa, aac_ids = aac_g1, snv_ids = noncoding_g1, samples = eligible_tumors)
  baseline_rates_g1 = rowSums(as.matrix(baseline_rates_g1[, -"Unique_Patient_Identifier"], rownames = baseline_rates_g1$Unique_Patient_Identifier))
  baseline_rates_g2 = baseline_mutation_rates(cesa, aac_ids = aac_g2, snv_ids = noncoding_g2, samples = eligible_tumors)
  baseline_rates_g2 = rowSums(as.matrix(baseline_rates_g2[, -"Unique_Patient_Identifier"], rownames = baseline_rates_g2$Unique_Patient_Identifier))

  # Find recurrent mutation status of all tumors in the two genes
  maf = cesa@maf[Unique_Patient_Identifier %in% eligible_tumors]
  setkey(maf, "variant_id")
  tumors_with_gene1_mutated = maf[all_g1_snv, unique(Unique_Patient_Identifier), nomatch = NULL]
  tumors_with_gene2_mutated = maf[all_g2_snv, unique(Unique_Patient_Identifier), nomatch = NULL]
  
  tumors_with_both_mutated = intersect(tumors_with_gene1_mutated,tumors_with_gene2_mutated)
  tumors_with_ONLY_gene1_mutated <- setdiff(tumors_with_gene1_mutated, tumors_with_both_mutated)
  tumors_with_ONLY_gene2_mutated <- setdiff(tumors_with_gene2_mutated, tumors_with_both_mutated)
  tumors_with_neither_mutated <- setdiff(eligible_tumors, c(tumors_with_gene1_mutated, tumors_with_gene2_mutated))
  
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
