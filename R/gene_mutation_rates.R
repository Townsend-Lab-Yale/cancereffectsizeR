#' Use dNdScv with tissue-specific covariates to calculate gene-level mutation rates
#'
#' This function calculates gene-level neutral mutation rates based on counts
#' of nonsynonymous and synonymous mutations per gene under the dNdScv package's model, 
#' as described in Martincorena et al. (https://doi.org/10.1016/j.cell.2017.09.042).
#'
#' @param cesa CESAnalysis object
#' @param covariates Name of gene mutation rate covariates to use (run
#'   list_ces_covariates() to see choices). For hg19 data only, you can also use "default"
#'   for dNdScv's non-tissue-specific covariates. If no covariates are available, set NULL
#'   to run without.
#' @param sample_group Which sample groups to include in the gene rate calculation;
#'   defaults to all groups. (To calculate different rates for different groups, you'll 
#'   run this function multiple times, changing this argument each time.)
#' @param save_all_dndscv_output default false; when true, saves all dndscv output, not
#'   just what's needed by cancereffectsizeR (will make object size very large)
#' @return CESAnalysis object with gene-level mutation rates calculated
#' @export
# don't change this function at all without being sure you're not messing up tests
gene_mutation_rates <- function(cesa, covariates = NULL, sample_group = NULL, save_all_dndscv_output = FALSE){
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa expected to be a CESAnalysis object", call. = F)
  }
  
  if(is.null(sample_group)) {
    sample_group = cesa@groups
  }
  
  if(! is(sample_group, "character")) {
    stop("sample_group should be character")
  }
  sample_group = unique(sample_group)
  if (! all(sample_group %in% cesa@groups)) {
    possible_groups = setdiff(cesa@groups, "stageless") # historical
    if (length(possible_groups) == 0) {
      stop("The CESAnalysis has no user-defined sample groups, so you can't use the sample_group parameter.")
    }
    stop("Unrecognized sample groups. Those defined in the CESAnalysis are ", paste(possible_groups, sep = ", "), ".")
  }
  
  cesa = update_cesa_history(cesa, match.call())

  RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS
  gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes
  
  
  # select MAF records to use for gene mutation rate calculation (via dNdScv)
  # for now, records from TGS samples are kept out; in the future, we could check out dNdScv's panel sequencing features
  dndscv_samples = cesa@samples[(coverage == "exome" | coverage == "genome") & group %in% sample_group, Unique_Patient_Identifier]
  
  if(length(dndscv_samples) == 0) {
    if (length(sample_group) != length(cesa@groups)) {
      stop("Cannot run dNdScv on sample group(s) because there are no whole-exome or whole-genome samples.")
    } else {
      stop("Cannot run dNdScv because the CESAnalysis has no whole-exome or whole-genome samples.")
    }
  }
  
  
  # Don't allow rates to be overwritten (too confusing for now to risk overlapping sample group specification in sequential calls)
  if (! is.null(cesa@samples$gene_rate_grp)) {
    if (! all(is.na(cesa@samples[group %in% sample_group, gene_rate_grp]))) {
      stop("Gene mutation rates have already been calculated for some or all samples specified.")
    }
  }
  
  
  if(is.null(covariates)){
    cv = NULL
    genes_in_pca = NULL
    warning("Calculating gene mutation rates with no covariate data; stop and re-run with covariates if available.\n",
            "(Check with list_ces_covariates(); for hg19 only, you can also specify \"default\" for dNdScv default\n",
            "covariates.)", call. = F, immediate. = T)
  } else if (is(covariates, "character") && covariates[1] == "default") {
    if(cesa@ref_key == "ces_hg19_v1") {
      pretty_message("Loading dNdScv default covariates for hg19 (stop and re-run with tissue-specific covariates if available)...")
      data("covariates_hg19", package = "dndscv", envir = environment())
      genes_in_pca <- rownames(covs)
      cv = "hg19"
    } else {
      stop("There is no default covariates data for this genome build, so you'll need to supply your own\n",
           "or run without by setting covariates = NULL.", call. = F)
    }
  } else if (! is(covariates, "character") || length(covariates) != 1) {
    stop("covariates expected to be 1-length character. Check available covariates with list_ces_covariates()")
  } else {
    covariates = paste0("covariates/", covariates)
    if(! check_for_ref_data(cesa, covariates)) {
      stop("Covariates could not be found. Check available covariates with list_ces_covariates().")
    }
    this_cov_pca <- get_ref_data(cesa, covariates) 
    cv = this_cov_pca$rotation
    genes_in_pca = rownames(cv)
  }
  
  mutations = cesa@maf[Unique_Patient_Identifier %in% dndscv_samples & variant_type == "snv", 
                       .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
  if(mutations[, .N] == 0) {
    stop("Can't run dNdScv because there are no usable SNV mutations in the input samples.")
  }
  
  # Run in separate function for quick unit tests of gene_mutation_rates
  dndscv_output = run_dndscv(mutations = mutations, gene_list = genes_in_pca, cv = cv, refdb = RefCDS, gr_genes = gr_genes)
  
  
  # Get RefCDS data on number of synonymous mutations possible at each site
  # Per dNdScv docs, L matrices list "number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context"
  num_syn = sapply(RefCDS, function(x) colSums(x$L)[1])
  names(num_syn) = sapply(RefCDS, function(x) x$gene_name)
  
  dndscv_genes = dndscv_output$genemuts$gene_name # dndscv uses same set of genes for each stage
  num_syn = num_syn[names(num_syn) %in% dndscv_genes]
  
  
  message("Using dNdScv output to calculate gene-level mutation rates...")
  number_of_tumors_in_this_subset <- length(unique(dndscv_output$annotmuts$sampleID))
  if(dndscv_output$nbreg$theta>1){
    # see page e4 of dNdScv paper (Martincorena 2017, Cell)
    mutrates_vec <- ((dndscv_output$genemuts$n_syn +
                        dndscv_output$nbreg$theta -
                        1) /
                       (1 +
                          (dndscv_output$nbreg$theta /
                             dndscv_output$genemuts$exp_syn_cv)
                       )
    ) /
      num_syn /
      number_of_tumors_in_this_subset
    
  } else{
    mutrates_vec <- rep(NA,length(dndscv_output$genemuts$exp_syn_cv))
    syn_sites <- num_syn
    for(i in 1:length(mutrates_vec)){
      mutrates_vec[i] <-  max(dndscv_output$genemuts$exp_syn_cv[i],  ((dndscv_output$genemuts$n_syn[i] +
                                                                         dndscv_output$nbreg$theta -
                                                                         1) /
                                                                        (1 +
                                                                           (dndscv_output$nbreg$theta /
                                                                              dndscv_output$genemuts$exp_syn_cv[i])
                                                                        ))) /
        num_syn[i] /
        number_of_tumors_in_this_subset
    }
  }
  
  if (ncol(cesa@mutrates) == 0) {
    curr_rate_group = 'rate_grp_1'
  } else {
    curr_rate_group = paste0('rate_grp_', ncol(cesa@mutrates))
  }
  mutrates_dt = data.table(gene = dndscv_output$genemuts$gene_name, rate = mutrates_vec)
  setnames(mutrates_dt, "rate", curr_rate_group)
  
  
  # Genes in gr_genes that are not present in the covariates data (a few are missing, usually),
  # won't get rates calcualted by dNdScv. Here, we assign them the rate of the nearest gene,
  # as measured by center-to-center distance.
  missing_genes = setdiff(GenomicRanges::mcols(gr_genes)["names"][,1], mutrates_dt$gene)
  
  tmp = as.data.table(gr_genes)
  gene_intervals = tmp[, .(chr = as.character(seqnames[1]), start = min(start), end = max(end)), by = "names"]
  setnames(gene_intervals, "names", "gene")
  gene_intervals[, center := start + trunc((end - start)/2)]
  gene_intervals = gene_intervals[order(chr, center)]
  present_genes = gene_intervals[! gene %in% missing_genes]
  not_present = gene_intervals[gene %in% missing_genes]
  
  # find nearest genes by chr/center to the missing genes from the present genes,
  # and eliminate the rare tie by taking the first match for each gene
  nearest_genes = present_genes[not_present, on = c("chr", "center"), roll = "nearest"]
  nearest_genes = nearest_genes[! duplicated(i.gene)]
  
  setkey(mutrates_dt, "gene")
  nearest_rates = mutrates_dt[nearest_genes$gene, -"gene"]
  missing_rates = cbind(gene = nearest_genes$i.gene, nearest_rates)
  mutrates_dt = rbind(mutrates_dt, missing_rates)
  
  # keep just the main gene-level selection output from dNdScv, unless user wanted everything
  if(! save_all_dndscv_output) {
    dndscv_output = dndscv_output$sel_cv
  }
  cesa@samples[group %in% sample_group, gene_rate_grp := curr_rate_group]
  
  if (ncol(cesa@mutrates) > 0) {
    cesa@mutrates = cesa@mutrates[mutrates_dt, on = "gene"]
  } else {
    cesa@mutrates = mutrates_dt
  }
  
  dndscv_output = list(dndscv_output)
  names(dndscv_output) = curr_rate_group
  cesa@dndscv_out_list = c(cesa@dndscv_out_list, dndscv_output)
  cesa@advanced$locked = T
  
  if (! any(is.na(cesa@samples$gene_rate_grp))) {
    cesa@advanced$gene_rates_done = T
  }
  return(cesa)
}

#' Internal function to run dNdScv
#' 
#' @keywords internal
run_dndscv = function(mutations, gene_list, cv, refdb, gr_genes) {
  # hacky way of forcing an object of name gr_genes into the dndscv::dndscv function environment,
  # since the object is required by dndscv but there's no argument to supply your own copy of it
  our_env = new.env(parent = environment(dndscv::dndscv))
  our_env$gr_genes = gr_genes
  our_dndscv = dndscv::dndscv
  environment(our_dndscv) = our_env
  
  message("Running dNdScv...")
  withCallingHandlers(
    {
      dndscv_raw_output = our_dndscv(mutations = mutations, gene_list = gene_list, cv = cv, refdb = refdb)
    }, error = function(e) {
      if (startsWith(conditionMessage(e), "bad 'file' argument"))  {
        stop("You need to update dNdScv. Try running \"remotes::update_packages(packages = \"dndscv\")\".")
      }
    }, warning = function(w) {
      dndscv_msg = conditionMessage(w)
      # Recurrent mutations are presumed real in CES
      # As of 09/03/20, dNdScv's MNV check is flawed (and CES has already handled them in load_maf())
      if (startsWith(dndscv_msg, "Same mutations observed in different sampleIDs") || 
          startsWith(dndscv_msg, "Mutations observed in contiguous sites")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  return(dndscv_raw_output)
}



#' Fill missing genes
#' 
#' For each gene in gr_genes missing from gene rate table, assign the 
#' rate of the nearest non-missing gene. Gene rate table is assumed
#' two-column (gene and rate).
#' 
#' @param gr_genes reference data set style gr_genes
#' @param gene_rates two-column data.table (gene, rate)
#' @keywords internal
fill_missing_genes = function(gr_genes, gene_rates) {
  missing_genes = setdiff(GenomicRanges::mcols(gr_genes)["names"][,1], gene_rates$gene)
  tmp = as.data.table(gr_genes)
  gene_intervals = tmp[, .(chr = as.character(seqnames[1]), start = min(start), end = max(end)), by = "names"]
  setnames(gene_intervals, "names", "gene")
  gene_intervals[, center := start + trunc((end - start)/2)]
  gene_intervals = gene_intervals[order(chr, center)]
  present_genes = gene_intervals[! gene %in% missing_genes]
  not_present = gene_intervals[gene %in% missing_genes]
  
  # find nearest genes by chr/center to the missing genes from the present genes,
  # and eliminate the rare tie by taking the first match for each gene
  nearest_genes = present_genes[not_present, on = c("chr", "center"), roll = "nearest"]
  nearest_genes = nearest_genes[! duplicated(i.gene)]
  
  setkey(gene_rates, "gene")
  nearest_rates = gene_rates[nearest_genes$gene, -"gene"]
  missing_rates = cbind(gene = nearest_genes$i.gene, nearest_rates)
  gene_rates = rbind(gene_rates, missing_rates)
  return(gene_rates)
}




