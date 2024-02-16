#' Use dNdScv with tissue-specific covariates to calculate gene-level mutation rates
#'
#' This function calculates gene-level neutral mutation rates based on counts
#' of nonsynonymous and synonymous mutations per gene under the dNdScv package's model, 
#' as described in \href{https://doi.org/10.1016/j.cell.2017.09.042}{Martincorena et al.}
#'
#' @param cesa CESAnalysis object
#' @param covariates Tissue-specific mutation rate covariates. Typically, supply the
#'   covariates object from your refset (e.g., ces.refset.hg19$covariates$bladder), or the
#'   object name ("bladder"). Run list_ces_covariates() to see choices. For hg19 data
#'   only, set to "hg19" to use dNdScv's non-tissue-specific covariates. If no appropriate
#'   covariates data are available, set to NULL to run without. Finally, you can also
#'   supply custom covariates data in the form of a matrix or prcomp object (see website
#'   for details).
#' @param samples Which samples to include in the current run. Defaults to all samples. Can be a
#'   vector of patient_ids, or a data.table containing rows from the CESAnalysis
#'   sample table.
#' @param dndscv_args Custom arguments to pass to dNdScv. (The arguments \code{mutations},
#'   \code{gene_list}, \code{cv}, and \code{refdb} are supplied by cancereffectsizeR and can't be
#'   substituted.)
#' @param save_all_dndscv_output Default FALSE; when TRUE, saves all dndscv output, not
#'   just what's needed by cancereffectsizeR. (Full output can be very large, in the gigabytes.)
#' @return CESAnalysis object with gene-level mutation rates calculated
#' @export
# don't change this function at all without being sure you're not messing up tests
gene_mutation_rates <- function(cesa, covariates = NULL, samples = character(), dndscv_args = list(), 
                                save_all_dndscv_output = FALSE) {
  
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa expected to be a CESAnalysis object", call. = F)
  }
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  
  sample_info = select_samples(cesa, samples)
  
  # select MAF records to use for gene mutation rate calculation (via dNdScv)
  # for now, records from TGS samples are kept out; in the future, we could check out dNdScv's panel sequencing features
  curr_run_grp_samples = sample_info$patient_id
  sample_info = sample_info[coverage %in% c("exome", "genome")]
  if(sample_info[, .N] == 0) {
    stop("Cannot run dNdSCv because there are no whole-exome or whole-genome samples among the input samples.")
  } 
  
  # Don't allow rates to be overwritten (too confusing for now to risk overlapping sample group specification in sequential calls)
  if (! all(is.na(sample_info$gene_rate_grp))) {
    msg = paste0("Gene mutation rates have already been calculated for some or all samples specified. If you ",
                 "want to re-run, run clear_gene_rates() first.")
    stop(pretty_message(msg, emit = F))
  }
  dndscv_samples = sample_info$patient_id
  
  if (! is(dndscv_args, "list") || uniqueN(names(dndscv_args)) != length(dndscv_args)) {
    stop("dndscv_args should a named list of arguments to pass.")
  }
  
  reserved_args = c("mutations", "gene_list", "cv", "refdb")
  if(any(reserved_args %in% names(dndscv_args))) {
    stop("The following arguments are passed to dNdScv automatically and cannot be specified via ",
         "dndscv_args:\n", paste(reserved_args, collapse = ", "), '.')
  }

  RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS.dndscv
  gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes.dndscv
  if(is.null(RefCDS)) {
    RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS
    gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes
    full_gr_genes = gr_genes
  } else {
    full_gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes
  }
  using_cds_rates = "gene" %in% names(GenomicRanges::mcols(gr_genes))
  
  skip_covariate_validation = FALSE
  if(is.null(covariates)){
    skip_covariate_validation = TRUE
    cv = NULL
    genes_in_pca = NULL
    warning("Calculating gene mutation rates with no covariate data; stop and re-run with covariates if available.\n",
            "(Check with list_ces_covariates(); for hg19 only, you can also specify \"default\" for dNdScv default\n",
            "covariates.)", call. = F, immediate. = T)
  } else if (is(covariates, "character") && (covariates[1] == "default" || covariates[1] == "hg19")) {
    skip_covariate_validation = ! using_cds_rates # will need to handle gene/CDS issue if using gene-based dNdScv covariates with CDS-based refset
    if(cesa@advanced$genome_info$build_name == 'hg19') {
      pretty_message("Loading dNdScv default covariates for hg19 (stop and re-run with tissue-specific covariates if available)...")
      data("covariates_hg19", package = "dndscv", envir = environment())
      cv = covs # dNdScv object from data() call
      genes_in_pca <- rownames(cv)
    } else {
      stop("There is no default covariates data for this genome build, so you'll need to supply your own\n",
           "or run without by setting covariates = NULL.", call. = F)
    }
  } else if (is(covariates, "prcomp") || is(covariates, 'matrix')) {
      # To-do: validate custom covariates
      cv = covariates
      if(is(cv, "prcomp")) {
        cv = cv$rotation
      }
      genes_in_pca = rownames(cv) 
      
  } else if (! is(covariates, "character") || length(covariates) != 1) {
    stop("covariates expected to be 1-length character. Check available covariates with list_ces_covariates()")
  } else {
    covariates = paste0("covariates/", covariates)
    if(! check_for_ref_data(cesa, covariates)) {
      stop("Covariates could not be found. Check available covariates with list_ces_covariates().")
    }
    this_cov_pca = get_ref_data(cesa, covariates) 
    cv = this_cov_pca$rotation
    genes_in_pca = rownames(cv)
  }
  
  # When gr_genes is protein-based, need to track which pids go with each gene identifier
  if(using_cds_rates) {
    # If rates are by CDS, synch-up gene-based covariates by assuming same covariates across each gene.
    # In the future, CDS-based (or, maybe better, coordinate-based) covariates may be supported.
    pid_to_gene = unique(as.data.table(GenomicRanges::mcols(gr_genes)))
    setnames(pid_to_gene, 'names', 'pid')
  }
  
  # Validate covariates
  if (! skip_covariate_validation) {
    if (using_cds_rates) {
      # For compatibility with cancereffectsizeR covariates that split isoforms for CDKN2A 
      # (The isoforms use the same covariates, so we can just copy one.)
      current_cv_genes = rownames(cv)
      if (pid_to_gene["CDKN2A", .N, on = "gene"] > 0 && ! "CDKN2A" %in% current_cv_genes & "CDKN2A.p14arf" %in% current_cv_genes) {
        to_copy = cv["CDKN2A.p14arf", , drop = F]
        rownames(to_copy) = "CDKN2A"
        cv = rbind(cv, to_copy)
      }
      present_genes = intersect(pid_to_gene$gene, rownames(cv)) # some may be missing, which gets dealt with later
      present_pid_to_gene = pid_to_gene[present_genes, on = "gene"]
      cv = cv[present_pid_to_gene$gene, ]
      rownames(cv) = genes_in_pca = present_pid_to_gene$pid
    }
    
    if(is.null(genes_in_pca)) {
      stop("Covariates matrix must have gene names as rownames.")
    }
    
    if (length(genes_in_pca) != uniqueN(genes_in_pca)) {
      stop("Some gene names are repeated in covariates matrix rownames.")
    }
    bad_genes = setdiff(genes_in_pca, get_ref_data(cesa, "gene_names"))
    
    # Deal with ces.refset.hg19's special handling of CDKN2A
    if ("CDKN2A" %in% bad_genes && cesa@ref_key == 'ces.refset.hg19' && 
        ! any(c("CDKN2A.p14arf", "CDKN2A.p16INK4a") %in% genes_in_pca)) {
      to_insert = cv[c("CDKN2A", "CDKN2A"), ]
      rownames(to_insert) = c("CDKN2A.p14arf", "CDKN2A.p16INK4a")
      cv = cv[! rownames(cv) == "CDKN2A",]
      cv = rbind(cv, to_insert)
      bad_genes = setdiff(bad_genes, "CDKN2A")
      genes_in_pca = rownames(cv)
    }
    
    # On CDS refsets, pid and gene names are not really properly validated yet
    prop_bad = length(bad_genes)/length(genes_in_pca)
    if(prop_bad > 0.2 && ! using_cds_rates) {
      if(prop_bad == 1) {
        stop('None of the gene names in the reference data are present in the covariates data.')
      }
      if(length(bad_genes) > 30) {
        bad_genes = c(bad_genes[1:30], '...')
      }
      warning("Rownames of covariates matrix include genes not present in reference data:\n",
           paste(bad_genes, collapse = ', '), '.')
    }
    if(length(bad_genes) > 0 && ! using_cds_rates) {
      good_genes = intersect(setdiff(genes_in_pca, bad_genes), names(RefCDS))
      genes_in_pca = good_genes
      cv = cv[genes_in_pca, ]
    }
  }
 
  
  
  mutations = cesa@maf[patient_id %in% dndscv_samples & variant_type == "sbs", 
                       .(patient_id, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
  if(mutations[, .N] == 0) {
    stop("Can't run dNdScv because there are no usable sbs mutations in the input samples.")
  }
  
  # Run in separate function for quick unit tests of gene_mutation_rates
  dndscv_args = c(list(mutations = mutations, gene_list = genes_in_pca, cv = cv, refdb = RefCDS, 
                       gr_genes = gr_genes),
                  dndscv_args)
  dndscv_output = do.call(run_dndscv, dndscv_args)
  
  # Will want to get CIs on gene mutation rates. (Encapsulating use of predict for testing purposes.)
  fit = get_dndscv_model_fit(dndscv_output)
  
  theta = dndscv_output$nbreg$theta # parameter controlling variability in rates across genes
  lower = exp(fit$fit - 1.96 * fit$se.fit)
  upper = exp(fit$fit + 1.96 * fit$se.fit)
  
  # Note exp(fit$fit) will be equivalent to dndscv_output$genemuts$exp_syn_cv
  mle_rate = dndscv_output$genemuts$exp_syn_cv
  
  # Get RefCDS data on number of synonymous mutations possible at each site
  # Per dNdScv docs, L matrices list "number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context"
  message("Extracting gene mutation rates from dNdScv output...")
  dndscv_sample_num = uniqueN(dndscv_output$annotmuts$sampleID)
  actual_syn_counts = dndscv_output$genemuts$n_syn
  dndscv_gene_names = dndscv_output$genemuts$gene_name
  nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])
  if(dndscv_output$nbreg$theta>1){
    # see page e4 of dNdScv paper (Martincorena 2017, Cell)
    mutrates_vec <- ((actual_syn_counts + theta - 1) /
                    (1 + (theta / mle_rate))) /
                    nsyn_sites / dndscv_sample_num
    lower_rates = ((actual_syn_counts + theta - 1) /
                      (1 + (theta / lower))) /
                    nsyn_sites / dndscv_sample_num
    upper_rates = ((actual_syn_counts + theta - 1) /
                     (1 + (theta / upper))) /
                    nsyn_sites / dndscv_sample_num

  } else{
    # theta should definitely be >1, so none of this should really happen
    mutrates_vec = pmax(mle_rate,
                        (actual_syn_counts + theta - 1) / (1 +(theta / mle_rate))) /
      nsyn_sites / dndscv_sample_num
    
    # No CIs when theta<1
    upper_rates = rep(NA, length(mutrates_vec))
    lower_rates = upper_rates
  }
  
  curr_rate_group = as.integer(max(c(0, na.omit(cesa@samples$gene_rate_grp))) + 1)
  rate_grp_colname = paste0("rate_grp_", curr_rate_group)
  mutrates_dt = data.table(gene = dndscv_output$genemuts$gene_name, rate = mutrates_vec,
                           low = lower_rates, high = upper_rates)
  
  ci_low_name = paste0(rate_grp_colname, '_95_low')
  ci_high_name = paste0(rate_grp_colname, '_95_high')
  setnames(mutrates_dt, c("rate", "low", "high"), c(rate_grp_colname, ci_low_name, ci_high_name))

  # Genes in gr_genes that are not present in the covariates data (a few are missing,
  # usually), won't get rates calculated by dNdScv. Here, we assign them the rate of the
  # nearest gene, as measured by center-to-center distance. (And, in the case of CDS
  # rates, this same logic works.)
  
  if(using_cds_rates || is.null(full_gr_genes$gene)) {
    missing_genes = setdiff(GenomicRanges::mcols(full_gr_genes)["names"][,1], mutrates_dt$gene)
  } else {
    missing_genes = setdiff(GenomicRanges::mcols(full_gr_genes)["gene"][,1], mutrates_dt$gene)
  }
  
  
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
  
  # convert dNdScv's main output to data.table
  setDT(dndscv_output$sel_cv)

  cesa@samples[curr_run_grp_samples, gene_rate_grp := curr_rate_group, on = "patient_id"]
  
  if(using_cds_rates) {
    setnames(mutrates_dt, 'gene', 'pid')
    mutrates_dt[pid_to_gene, gene := gene, on = 'pid']
    setcolorder(mutrates_dt, c('pid', 'gene'))
    setnames(dndscv_output$sel_cv, 'gene_name', 'pid')
    dndscv_output$sel_cv[mutrates_dt, gene := gene, on = 'pid']
  }
  
  # keep just the main gene-level selection output from dNdScv, unless user wanted everything
  if(! save_all_dndscv_output) {
    dndscv_output = dndscv_output$sel_cv
  }
  
  
  if (curr_rate_group > 1) {
    if (using_cds_rates) {
      cesa@mutrates = cesa@mutrates[mutrates_dt[, -"gene"], on = "pid"]
    } else {
      cesa@mutrates = cesa@mutrates[mutrates_dt, on = "gene"]
    }
  } else {
    cesa@mutrates = mutrates_dt
  }
  
  dndscv_output = list(dndscv_output)
  names(dndscv_output) = rate_grp_colname
  cesa@dndscv_out_list = c(cesa@dndscv_out_list, dndscv_output)
  return(cesa)
}

#' This little function called by gene_mutation_rates() is separated for testing purposes.
#'
#' @keywords internal
get_dndscv_model_fit = function(dndscv_output) {
  predict(dndscv_output$nbreg, se.fit = T)
}

#' Internal function to run dNdScv
#' 
#' @keywords internal
run_dndscv = function(mutations, gene_list, cv, refdb, gr_genes, ...) {
  # hacky way of forcing an object of name gr_genes into the dndscv::dndscv function environment,
  # since the object is required by dndscv but there's no argument to supply your own copy of it
  our_env = new.env(parent = environment(dndscv::dndscv))
  our_env$gr_genes = gr_genes
  our_dndscv = dndscv::dndscv
  environment(our_dndscv) = our_env
  
  message("Running dNdScv...")
  withCallingHandlers(
    {
      dndscv_raw_output = our_dndscv(mutations = mutations, gene_list = gene_list, cv = cv, refdb = refdb, ...)
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


#' Clear regional mutation rates
#' 
#' Remove all gene/coding region mutation rates, usually in order to re-run with different
#' parameters without having to create a new CESAnalysis.
#' 
#' @param cesa CESAnalysis
#' @return The CESAnalysis with rates cleared
#' @export
clear_gene_rates = function(cesa = NULL) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  
  cesa = copy_cesa(cesa)
  cesa@dndscv_out_list = list()
  cesa@mutrates = data.table()
  cesa@samples[, gene_rate_grp := NA_integer_]
  
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}

