#' Gene-level epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise gene-level epistasis.
#' Selection at the gene level is assumed to act through all 
#' specified variants (see options), with the selection for each gene's variants allowed
#' to vary based on the mutational status of the other gene's variants.
#' 
#' Only samples that have coverage at all included sites in both genes can be included in the
#' inference since samples lacking full coverage may or may not have mutations at the
#' uncovered sites.
#' 
#' @param cesa CESAnalysis object
#' @param genes Vector of gene names; SIs will be calculated for all gene pairs.
#' @param variants Which variants to include in inference for each gene. Either
#'   "recurrent" for all variants present in two or more samples (across all MAF data),
#'   "nonsilent" for all coding variants except non-splice-site synonymous variants, or a
#'   data.table containing all variants to include (as returned by
#'   \code{select_variants()} or by subsetting\code{[CESAnalysis\]$variants}). For
#'   noncoding variants with multiple gene annotations, the one listed in the "gene"
#'   column is used. In the recurrent method, nearby noncoding variants may be included.
#' @param samples Which samples to include in inference. Defaults to all samples. Can be a
#'   vector of Unique_Patient_Identifiers, or a data.table containing rows from the
#'   CESAnalysis sample table.
#' @param run_name Optionally, a name to identify the current run.
#' @param conf Confidence interval size from 0 to 1 (.95 -> 95\%). NULL skips calculation,
#'   which may be helpful to reduce runtime when analyzing many gene pairs.
#' @param cores number of cores for parallel processing of gene pairs
#' @return CESAnalysis with epistasis analysis results added
#' @export
ces_gene_epistasis = function(cesa = NULL, genes = NULL, variants = NULL,
                              samples = character(), run_name = "auto", cores = 1, conf = .95)
{
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  
  samples = select_samples(cesa, samples)
  num_excluded = cesa@samples[, .N] - samples[, .N]
  if(num_excluded > 0) {
    pretty_message(paste0("Note that ", num_excluded, " samples are being excluded from selection inference."))
  }
  
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if(run_name %in% names(cesa@epistasis)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@epistasis) + 1
    run_name = paste0('gene_epistasis_', run_number)
    while(run_name %in% names(cesa@epistasis)) {
      run_number = run_number + 1
      run_name = paste0('gene_epistasis_', run_number)
    }
  }
  
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }

  # genes can be list of character pairs or a character vector, from which all possible pairs will be taken
  gene_pairs = NULL
  if(is(genes, "list") && all(sapply(genes, is.character)) && all(sapply(genes, length) ==2)) {
    gene_pairs = unique(genes)
    genes = unique(unlist(genes))
  } else if (is(genes, "character")) {
    if (length(genes) < 2) {
      stop("Supply at least two genes to analyze.")
    }
  } else{
    stop("genes should be character (at least two gene names) or a list of gene pairs.")
  }
  
  invalid_genes = setdiff(genes, .ces_ref_data[[cesa@ref_key]]$gene_names)
  num_invalid = length(invalid_genes)
  if(num_invalid > 0) {
    if(num_invalid > 19) {
      invalid_genes = c(invalid_genes[1:15], paste0("and ", num_invalid - 15, " more"))
    }
    stop("Input genes are not present in reference data: ", paste(invalid_genes, collapse = ", "))
  } 
  
  maf = cesa@maf[variant_type == "snv"]
  
  if(is.null(variants)) {
    variants = 'recurrent'
    pretty_message(paste0("As of cancereffectsizeR 2.4.0, there is a \"variants\" argument that lets you control ",
      "which variants from each gene are used in epistasis inference. Defaulting to \"recurrent\" strategy to match ",
      "previous behavior. See docs for other options."))
  }
  if(is.character(variants) && length(variants) == 1) {
    if(variants == 'recurrent') {
      variants_to_use = select_variants(cesa = cesa, genes = genes, min_freq = 2)
    } else if (variants == 'nonsilent') {
      variants_to_use = select_variants(cesa  = cesa, genes = genes, min_freq = 1)
      snv_to_use = variants_to_use[variant_type == 'snv']
      variants_to_use = variants_to_use[variant_type == 'aac'][aa_ref != aa_alt | essential_splice == TRUE]
      variants_to_use = rbind(variants_to_use, snv_to_use)
    } else {
      stop("variants should be \"recurrent\", \"nonsilent\", or a data.table.")
    }
  } else if (is.data.table(variants)) {
    if (! 'variant_id' %in% names(variants)) {
      stop("Input variants table lacks variant_id column")
    }
    message("Verifying input variant table....", appendLF = F)
    variants_to_use = select_variants(cesa = cesa, variant_ids = variants$variant_id)
    if(attr(variants_to_use, 'nonoverlapping') == FALSE) {
      msg = paste0("Input variants table contains overlapping variant records (e.g., amino acid substitutions that ",
           "contain the same constituent SNVs). To avoid confusion, ",
           "all variants in the inference should be nonoverlapping.")
      stop(pretty_message(msg, emit = F))
    }
    message(" Done.")
  } else if (is.null(variants)) {
    stop("Argument variants should be set to \"recurrent\", \"nonsilent\", or a data.table of variants to include for each gene.")
  } else {
    stop("variants should be \"recurrent\", \"nonsilent\", or a data.table.")
  }
	
	
	# use "all_genes" column because SNVs might have multiple genes, and only one appears in the "gene" column
	genes_present = unique(variants_to_use$gene)
	genes_to_analyze = intersect(genes_present, genes)
	num_passing_genes = length(genes_to_analyze)
  if(num_passing_genes == 0) {
    stop("None of your requested genes have any eligible variants to include in selection inference. ",
         "Changing the \"variants\" argument might help.")
  }
	if(num_passing_genes != length(genes)) {
	  num_missing = length(genes) - num_passing_genes
	  pretty_message(paste0(num_missing, " of your input genes had no eligible variants in input, so they will not be included in pairwise epistasis inference.",
	                        "You may want to verify that your \"variants\" input is what you intended."))
	}

	if (num_passing_genes == 1) {
	  msg = paste0("Only 1 requested gene (", genes_to_analyze, " has eligible variants  in input, so epistasis can't be tested. ",
	               "You may want to verify that your \"variants\" input is what you intended.")
	  stop(pretty_message(msg, emit = F))
	}
	
	if(is.null(gene_pairs)) {
	  gene_pairs = utils::combn(sort(genes_to_analyze), 2, simplify = F)
	} else {
	  gene_pairs = gene_pairs[sapply(gene_pairs, function(x) all(x %in% genes_to_analyze))]
	}
	
	gr_table = as.data.table(.ces_ref_data[[cesa@ref_key]]$gr_genes)
	# Handle pid-based gr_genes
	if('gene' %in% names(gr_table)) {
	  gr_table[, names := gene]
	  gr_table[, gene := NULL]
	}
	setkey(gr_table, "names")
	passing_gene_ranges = gr_table[genes_to_analyze]
	setkey(passing_gene_ranges, "seqnames", "start", "end")
  
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
    gene_pairs = gene_pairs[-bad_pairs]
    if (length(gene_pairs) == 0) {
      stop("There are no remaining gene pairs to test for epistasis.")
    }
  }
  
  # use separate CompoundVariantSets for each pair to avoid possible gene-overlap issues
  setkey(variants_to_use, 'gene')
  results = rbindlist(pbapply::pblapply(X = gene_pairs, FUN = 
    function(x) {
      variant_ids = setNames(list(variants_to_use[x[1], variant_id], variants_to_use[x[2], variant_id]), x)
      comp = CompoundVariantSet(cesa, variant_ids)
      results = pairwise_variant_epistasis(cesa = cesa, variant_pair = c(1, 2), samples = samples, compound_variants = comp,
                                           conf = conf)
    }, cl = cores
  ))
  
  # pairwise epistasis function uses v1/v2 in parameter names (as in, variants); sub in g1/g2 for gene
  colnames(results) = gsub("_v1", "_g1", colnames(results))
  colnames(results) = gsub("_v2", "_g2", colnames(results))
  
  if (results[(joint_cov_samples_just_g1 == 0 | joint_cov_samples_just_g2 == 0) & 
              joint_cov_samples_with_both == 0, .N] > 0) {
    pretty_message(paste0("Some gene pairs had no eligible variants in one or both genes of jointly-covering samples. ",
                   "Epistatic selection intensities are all NAs for these pairs."))
  }
  # Make column names match historical gene-based output.
  setnames(results, c("variant1", "variant2", "joint_cov_samples_just_g1", "joint_cov_samples_just_g2", 
                      "joint_cov_samples_with_both", "joint_cov_samples_with_neither"),
           c("gene_1", "gene_2", "joint_cov_samples_just_g1_mut", "joint_cov_samples_just_g2_mut", 
             "joint_cov_samples_both_mut", "joint_cov_samples_no_mut"))
           
  
  results = list(results)
  names(results) = run_name
  cesa@epistasis = c(cesa@epistasis, results)
  return(cesa)
}

#' Variant-level pairwise epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise epistasis between pairs of variants.
#' CompoundVariantSets are supported.
#' 
#' @param cesa CESAnalysis
#' @param variants To test pairs of variants, supply a list where each element is a
#'   2-length vector of CES-style variant IDs. Alternatively (and often more usefully),
#'   supply a CompoundVariantSet (see \code{define_compound_variants}) to test all pairs
#'   of compound variants in the set.
#' @param samples Which samples to include in inference. Defaults to all samples.
#'   Can be a vector of Unique_Patient_Identifiers, or a data.table containing rows from
#'   the CESAnalysis sample table.
#' @param run_name Optionally, a name to identify the current run.
#' @param conf confidence interval size from 0 to 1 (.95 -> 95\%); NULL skips calculation,
#'   reducing runtime.
#' @param cores number of cores for parallel processing of variant pairs
#' @return a data table with pairwise-epistatic selection intensities and variant
#'   frequencies for in tumors that have coverage at both variants in each pair
#' @export
ces_epistasis = function(cesa = NULL, variants = NULL, samples = character(), run_name = "auto", cores = 1, conf = .95) {
  
  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }
  
  if (! is.null(variants) && length(variants) < 1) {
    stop("variants is 0-length.", call. = F)
  }
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should specify a CESAnalysis object", call. = F)
  }
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  
  samples = select_samples(cesa, samples)
  num_excluded_samples = cesa@samples[, .N] - samples[, .N]
  if(num_excluded_samples > 0) {
    message("Note that ", num_excluded_samples, " are being excluded from selection inference.")
  }
  
  if(run_name %in% names(cesa@epistasis)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@epistasis) + 1
    run_name = paste0('variant_epistasis_', run_number)
    while(run_name %in% names(cesa@epistasis)) {
      run_number = run_number + 1
      run_name = paste0('variant_epistasis_', run_number)
    }
  }
  
  if(is(variants, "CompoundVariantSet")) {
    index_pairs = utils::combn(1:length(variants), 2, simplify = F)
    results = rbindlist(pbapply::pblapply(X = index_pairs, FUN = pairwise_variant_epistasis, samples = samples, compound_variants = variants, cesa=cesa, conf = conf, cl = cores))
  } else if (is(variants, "list")) {
    setkey(cesa@mutations$amino_acid_change, "aac_id")
    setkey(cesa@mutations$snv, "snv_id")
    
    # validate that all variant IDs appear in mutation annotations
    aac_names = cesa@mutations$amino_acid_change[, .(aac_id, variant_name = paste(gene, aachange, sep = "_"))]
    notify_name_conversion = F
    for (i in 1:length(variants)) {
      pair = variants[[i]]
      if (! is(pair, "character") || length(pair) != 2) {
        stop("All elements of variants list must be character vectors of length 2 (for two variant IDs)")
      }
      
      for (j in 1:2) {
        variant = pair[j]
        if (cesa@mutations$snv[variant, .N, nomatch = NULL] == 0 && aac_names[variant == aac_id, .N, nomatch = NULL] == 0) {
          variant = gsub(' ', '_', variant)
          record_by_name = aac_names[variant, on = "variant_name", nomatch = NULL]
          num_records = record_by_name[, .N]
          if (num_records == 0) {
            stop("Variant ", variant, " could not be found. Either it's not present in the CESAnalysis, or it's a bad ID.")
          } else if (num_records > 1) {
            stop("More than one variant in the CESAnalysis has the name \"", variant, "\". Use a full variant_id to specify which to test.")
          } else {
            notify_name_conversion = T
            variants[[i]][j] = record_by_name[, aac_id]
          }
        }
      }
    }
    if (notify_name_conversion) {
      pretty_message("Shorthand variant names were recognized and uniquely paired with cancereffectsizeR's aac_ids.")
    }
    results = rbindlist(pbapply::pblapply(X = variants, FUN = pairwise_variant_epistasis, cesa=cesa, samples = samples, conf = conf, cl = cores))
  } else {
    stop("variants should be of type list or CompoundVariantSet")
  }
  
  if(results[(joint_cov_samples_just_v1 == 0 | joint_cov_samples_just_v2 == 0) & 
             joint_cov_samples_with_both == 0, .N] > 0) {
    pretty_message(paste0("Some variant pairs had prevalence of 0 in one or both variants across jointly-covering input samples. ",
                          "Epistatic selection intensities are all NAs for these pairs."))
  }
  results = list(results)
  names(results) = run_name
  cesa@epistasis = c(cesa@epistasis, results)
  return(cesa)
  
}

#' Calculate SIs at variant level under pairwise epistasis model
#' @param cesa CESAnalysis
#' @param samples Validated samples data.table (as from select_samples())
#' @param variant_pair 2-length character of variant IDs, or 2-length numeric giving
#'   indices of CompoundVariantSet for the current two compound variants
#' @param compound_variants If testing a pair of compound variants, the CompoundVariantSet defining them
#' @keywords internal
pairwise_variant_epistasis = function(cesa, variant_pair, samples, conf, compound_variants = NULL) {
  running_compound = FALSE
  if (is(compound_variants, "CompoundVariantSet")) {
    running_compound = TRUE
    compound_variants = compound_variants[variant_pair]
    joint_coverage = c("genome", intersect(compound_variants@compounds$shared_cov[[1]], compound_variants@compounds$shared_cov[[2]]))
    v1 = compound_variants@compounds$compound_name[[1]]
    v2 = compound_variants@compounds$compound_name[[2]]
    v1_ids = compound_variants@snvs[compound_name == v1, snv_id]
    v2_ids = compound_variants@snvs[compound_name == v2, snv_id]
    variant_ids = c(v1_ids, v2_ids)
  } else {
    v1 = variant_pair[1]
    v2 = variant_pair[2]
    
    v1_coverage = c(cesa@mutations$amino_acid_change[v1, unlist(covered_in), nomatch = NULL], 
                    cesa@mutations$snv[v1, unlist(covered_in), nomatch = NULL])
    v2_coverage = c(cesa@mutations$amino_acid_change[v2, unlist(covered_in), nomatch = NULL], 
                    cesa@mutations$snv[v2, unlist(covered_in), nomatch = NULL])
    
    # Samples have to have v1 and v2 coverage (and samples with covered_regions == "genome" always have coverage)
    joint_coverage = c("genome", intersect(v1_coverage, v2_coverage))
    variant_ids = c(v1, v2)
  }

  eligible_tumors = samples[covered_regions %in% joint_coverage, Unique_Patient_Identifier]
  if (length(eligible_tumors) == 0) {
    warning(sprintf("No samples have coverage at both %s and %s, so this variant pair had to be skipped.", v1, v2), immediate. = T, call. = F)
  }
  all_rates = baseline_mutation_rates(cesa = cesa, variant_ids = variant_ids, samples = eligible_tumors)
  
  if (running_compound) {
    # sample_calls could include tumors that don't have coverage across both in pair
    tumors_with_v1 = intersect(compound_variants@sample_calls[[v1]], eligible_tumors)
    tumors_with_v2 = intersect(compound_variants@sample_calls[[v2]], eligible_tumors)
    tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
    tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
    tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
    tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
    
    v1_rates = all_rates[, ..v1_ids]
    v1_rates = rowSums(v1_rates)
    names(v1_rates) = all_rates$Unique_Patient_Identifier
    v2_rates = all_rates[, ..v2_ids]
    v2_rates = rowSums(v2_rates)
    names(v2_rates) = all_rates$Unique_Patient_Identifier
    
    # 2-item lists: first item is baseline rates for v1; second for v2
    with_just_1 = NULL; with_just_2 = NULL; with_both = NULL; with_neither = NULL
    if(length(tumors_just_v1) > 0){
      with_just_1 = list(v1_rates[tumors_just_v1], v2_rates[tumors_just_v1])
    }
    if(length(tumors_just_v2) > 0) {
      with_just_2 = list(v1_rates[tumors_just_v2], v2_rates[tumors_just_v2])
    }
    if (length(tumors_with_both) > 0) {
      with_both = list(v1_rates[tumors_with_both], v2_rates[tumors_with_both])
    }
    if(length(tumors_with_neither) > 0) {
      with_neither = list(v1_rates[tumors_with_neither], v2_rates[tumors_with_neither])
    }
  } else {
    setcolorder(all_rates, c("Unique_Patient_Identifier", v1, v2))
    setkey(all_rates, "Unique_Patient_Identifier")
    covered_maf = cesa@maf[eligible_tumors, on = "Unique_Patient_Identifier", nomatch = NULL]
    tumors_with_v1 = intersect(samples_with(cesa, v1), eligible_tumors)
    tumors_with_v2 = intersect(samples_with(cesa, v2), eligible_tumors)
    tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
    tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
    tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
    tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
    
    # 2-item lists: first item is baseline rates for v1; second for v2
    with_just_1 = as.list(all_rates[tumors_just_v1])[2:3]
    with_just_2 = as.list(all_rates[tumors_just_v2])[2:3]
    with_both = as.list(all_rates[tumors_with_both])[2:3]
    with_neither = as.list(all_rates[tumors_with_neither])[2:3]
  }

  # call factory function to get variant-specific likelihood function
  fn = pairwise_epistasis_lik(with_just_1, with_just_2, with_both, with_neither)
  par_init = formals(fn)[[1]]
  names(par_init) = bbmle::parnames(fn)
  
  # No point of testing epistasis if either variant doesn't appear
  if (length(tumors_with_v1) == 0 || length(tumors_with_v2) == 0) {
    early_output = list(variant1 = v1, variant2 = v2, ces_v1 = NA_real_, ces_v2 = NA_real_,
                        ces_v1_after_v2 = NA_real_, ces_v2_after_v1 = NA_real_, 
                        joint_cov_samples_just_v1 = length(tumors_just_v1),
                        joint_cov_samples_just_v2 = length(tumors_just_v2),
                        joint_cov_samples_with_both = length(tumors_with_both),
                        joint_cov_samples_with_neither = length(tumors_with_neither))
    if (! is.null(conf)) {
     si_names = names(par_init)
     ci_high_colname = paste0("ci_high_", conf * 100)
     ci_low_colname = paste0("ci_low_", conf * 100)
     low_colnames = paste(ci_low_colname, si_names, sep = "_")
     high_colnames = paste(ci_high_colname, si_names, sep = "_")
     ci_colnames = unlist(S4Vectors::zipup(low_colnames, high_colnames))
     ci_output = as.list(rep.int(NA_real_, 8))
     names(ci_output) = ci_colnames
     early_output = c(early_output, ci_output) # low/high for 4 parameters
    }
    return(early_output)
  }
  
  # find optimized selection intensities
  # the selection intensity when some group has 0 variants will be on the lower boundary; will muffle the associated warning
  withCallingHandlers(
    {
      fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = T, lower=1e-3, upper=1e9)
    },
    warning = function(w) {
      if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  params = bbmle::coef(fit)
  
  variant_ep_results = list(variant1 = v1, variant2 = v2, ces_v1 = params[1], ces_v2 = params[2],
                     ces_v1_after_v2 = params[3], ces_v2_after_v1 = params[4], 
                     joint_cov_samples_just_v1 = length(tumors_just_v1),
                     joint_cov_samples_just_v2 = length(tumors_just_v2),
                     joint_cov_samples_with_both = length(tumors_with_both),
                     joint_cov_samples_with_neither = length(tumors_with_neither))
  if(! is.null(conf)) {
    variant_ep_results = c(variant_ep_results, univariate_si_conf_ints(fit, fn, .001, 1e9, conf))
  }
  return(variant_ep_results)
}

