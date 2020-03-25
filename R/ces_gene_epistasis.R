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
  # Some optimx::opm optimization methods crash for unclear reasons: Rnmin, nmkb, newuoa, Rtnmin
  # Others should be skipped automatically because they aren't appropriate, but it's
  # necessary to skip them manually: "subplex", "snewtonm", "snewton", "CG", "BFGS","Nelder-Mead", "nlm", "lbfgs"
  working_methods = c("L-BFGS-B", "nlminb", "lbfgsb3", "Rcgmin", "Rvmmin","spg", "bobyqa", "hjkb", "hjn")
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
	genes_in_dataset = unique(maf$Gene_name)
	genes_to_analyze = genes[genes %in% genes_in_dataset]


	cesa_subset <- maf[Gene_name %in% genes_to_analyze,]
	cesa_subset$identifier <- cesa_subset$unique_variant_ID_AA
	cesa_subset$identifier[
	  which(sapply(strsplit(cesa_subset$identifier,split = " "),
	               function(x) length(x))==1)
	  ] <- paste(cesa_subset$Gene_name[
	    which(sapply(strsplit(cesa_subset$identifier,split = " "),
	                 function(x) length(x))==1)],
	    cesa_subset$identifier[
	      which(sapply(strsplit(cesa_subset$identifier,split = " "),
	                   function(x) length(x))==1)
	      ],sep=" ")

	cesa_subset_table <- table(cesa_subset$identifier)
	cesa_subset$recurrent_val <- cesa_subset_table[cesa_subset$identifier]

	recurrently_subed_genes <- unique(cesa_subset$Gene_name[cesa_subset$recurrent_val > 1])

    if(length(recurrently_subed_genes) < 2){
     stop("Less than 2 of the 'genes' you specified have recurrent variants.
          Choose 2 genes or more with recurrent variants to measure epistasis
          among those recurrently substituted variants. ")
    }

    genes_to_analyze <- recurrently_subed_genes

    selection_epistasis_results <- t(utils::combn(genes_to_analyze,2))
    selection_epistasis_results <- data.frame(t(selection_epistasis_results),stringsAsFactors=F)
    rownames(selection_epistasis_results) <- c("Variant_1","Variant_2")
    selection_epistasis_results_list <- as.list(selection_epistasis_results)
    gene_trinuc_comp = get_genome_data(cesa, "gene_trinuc_comp")
    selection_results = pbapply::pblapply(X = selection_epistasis_results_list,
                                           FUN = epistasis_gene_level,
                                           cesa=cesa,
                                           gene_trinuc_comp = gene_trinuc_comp,
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

epistasis_gene_level = function(genes_to_analyze, 
                                gene_trinuc_comp,
                                cesa,
                                optimx_args) {
  mutrates_list = cesa@mutrates_list
  MAF = cesa@maf[Variant_Type == "SNV"]
  trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  all_tumors = unique(MAF$Unique_Patient_Identifier)

  variant1 = genes_to_analyze[1]
  variant2 = genes_to_analyze[2]
  MAF_input1= MAF[Gene_name == variant1]
  MAF_input2= MAF[Gene_name == variant2]
  
  # when including target gene sequencing data, need to throw out any samples that don't have coverage at ALL variant sites in these genes
  # this could be a problem if some of the "exome" samples are actually whole-genome if they haven't been trimmed strictly enough
  
  # this is clunky but okay for now and will change with next refactor
  eligible_tumors = all_tumors # samples without coverage at all sites will get intersected out
  for (maf in list(MAF_input1, MAF_input2)) {
    for (i in 1:nrow(maf)) {
      site_coverage = unlist(maf[i, covered_in]) # this returns a character vector naming the covered regions with coverage
      tumors_covering_locus = cesa@samples[covered_regions %in% site_coverage, Unique_Patient_Identifier]
      eligible_tumors = intersect(eligible_tumors, tumors_covering_locus)
    }
  }
  
  # restrict MAFs to eligible tumors
  # temporarily make data.frame for compatibility
  MAF_input1 = data.frame(MAF_input1[Unique_Patient_Identifier %in% eligible_tumors,])
  MAF_input2 = data.frame(MAF_input2[Unique_Patient_Identifier %in% eligible_tumors,])
  
  variant_freq_1 <- table(MAF_input1$unique_variant_ID_AA)
  variant_freq_2 <- table(MAF_input2$unique_variant_ID_AA)

  # only run the selection algorithm if there are 2 or more tumors with
  # recurrent variants of each gene present.
  these_mutation_rates1 <-
    mutation_rate_calc(
      this_MAF = MAF_input1,
      gene = variant1,
      gene_mut_rate = mutrates_list,
      trinuc_proportion_matrix = trinuc_proportion_matrix,
      gene_trinuc_comp = gene_trinuc_comp,
      all_tumors = eligible_tumors,
      samples = cesa@samples)

  these_mutation_rates2 <-
    mutation_rate_calc(
      this_MAF = MAF_input2,
      gene = variant2,
      gene_mut_rate = mutrates_list,
      trinuc_proportion_matrix = trinuc_proportion_matrix,
      gene_trinuc_comp = gene_trinuc_comp,
      all_tumors = eligible_tumors,
      samples = cesa@samples)


  # since we are looking at selection at the gene level,
  # we only consider selection at sites that are recurrently
  # substituted.

  these_mutation_rates1 <- as.matrix(these_mutation_rates1[,colnames(these_mutation_rates1) %in% names(variant_freq_1[variant_freq_1>1])])
  these_mutation_rates2 <- as.matrix(these_mutation_rates2[,colnames(these_mutation_rates2) %in% names(variant_freq_2[variant_freq_2>1])])

  # calculate tumors with and without recurrent mutations in the two genes
  # only the tumors containing a recurrent variant factor into the selection analysis
  ## MAF_input1 and MAF_input2 are already restricted to just eligible tumors
  tumors_with_variant1_mutated <- unique(MAF_input1[which(MAF_input1$unique_variant_ID_AA %in% names(which(variant_freq_1>1))),"Unique_Patient_Identifier"])
  tumors_with_variant2_mutated <- unique(MAF_input2[which(MAF_input2$unique_variant_ID_AA %in% names(which(variant_freq_2>1))),"Unique_Patient_Identifier"])
  tumors_with_both_mutated <- base::intersect(tumors_with_variant1_mutated,tumors_with_variant2_mutated)
  tumors_with_ONLY_variant1_mutated <- tumors_with_variant1_mutated[which(!tumors_with_variant1_mutated %in% tumors_with_variant2_mutated)]
  tumors_with_ONLY_variant2_mutated <- tumors_with_variant2_mutated[which(!tumors_with_variant2_mutated %in% tumors_with_variant1_mutated)]
  tumors_with_neither_mutated <- setdiff(eligible_tumors, c(tumors_with_both_mutated,tumors_with_ONLY_variant1_mutated,tumors_with_ONLY_variant2_mutated))
  
  specific_mut_rates1=these_mutation_rates1
  specific_mut_rates2=these_mutation_rates2
  get_summed_mut_rates = function(tumors) {
    if(length(tumors) == 0) {
      return(NULL)
    }
    rates1 = specific_mut_rates1[tumors, , drop=F]
    rates2 = specific_mut_rates2[tumors, , drop=F]
    return(list(rowSums(rates1), rowSums(rates2)))
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
