#' Create a cancereffectsizeR analysis
#' 
#' Creates a CESAnalysis object, the central data structure of cancereffectsizeR.
#' @param ref_set name of reference data set to use; run \code{list_ces_ref_sets()} for
#'   available ref sets. Alternatively, the path to a custom reference data directory.
#' @param sample_groups Optionally, supply labels identifying different groups of samples.
#'   To be able to perform analyses that require ordered groups of samples, give the
#'   labels in proper order: e.g., c("Primary", "Metastatic"). If you will not be doing
#'   such analyses, ordering doesn't matter.
#' @return CESAnalysis object
#' @export
CESAnalysis = function(ref_set = "ces_hg19_v1", sample_groups = NULL) {
  
  # Check for and load reference data for the chosen genome/transcriptome data
  ref_set_name = ref_set
  if(is.null(ref_set_name) || ! is(ref_set_name, "character") || length(ref_set_name) != 1) {
    stop("Expected reference data source to be given as character. Run list_ces_ref_sets() to see available reference data.", call. = F)
  }
  
  ref_set_dirs = get_ref_set_dirs()
  if(ref_set_name %in% names(ref_set_dirs)) {
    data_dir = ref_set_dirs[ref_set_name]
    # avoid weird edge case
    if (dir.exists(ref_set_name)) {
      stop("There's a folder in your working directory with the same name as your chosen reference data set.\n",
           "Change your working directory or rename it, please.")
    }
  } else {
    if (! dir.exists(ref_set_name)) {
      if (grepl('/', ref_set_name)) {
        stop("Could not find reference data at ", ref_set_name)
      } else {
        stop("Invalid reference set name. Check spelling, or view available data sets with list_ces_ref_sets().")
      }
    }
    data_dir = ref_set_name
    ref_set_name = basename(ref_set_name)
    
    # To avoid confusion, you can't create a custom ref set with the same name as a built-in one
    builtin_sets = list.dirs(system.file("ref_sets/", package = "cancereffectsizeR"), full.names = F, recursive = F)
    if (ref_set_name %in% builtin_sets) {
      stop("The name of your reference data set (", ref_set_name, ") exactly matches a built-in CES reference set. Please rename it.")
    }
  }
  
  # preload some reference data, which will get stored in the .ces_ref_data env under ref_set_name
  preload_ref_data(data_dir)
  bsg = .ces_ref_data[[ref_set_name]]$genome
  
  # Validate sample_groups
  if (is.null(sample_groups)) {
    sample_groups = c("stageless")
  }
  
  if (is.numeric(sample_groups)) {
    sample_groups = as.character(sample_groups)
  }

  if (! is(sample_groups, "character")) {
    stop("sample_groups should be a character vector of chronological tumor states (e.g., Primary, Metastatic)")
  }
  
  # advanced is a grab bag of additional stuff to keep track of
  ## annotated: whether loaded MAF records are annotated
  ## using_exome_plus: whether previously loaded and any future generic exome data uses the "exome+" coverage option 
  ##  (either all generic data must, or none of it, based on choice of enforce_generic_exome_coverage on first load_maf call)
  ## recording: whether "run_history" is currently being recorded (gets set to false during some internal steps for clarity)
  ## locked: whether load_maf can still be used (can't load more data after trinuc_mutation_rates or gene_mutation_rates)
  advanced = list("version" = packageVersion("cancereffectsizeR"), annotated = F, using_exome_plus = F, recording = T, locked = F)
  cesa = new("CESAnalysis", run_history = character(),  ref_key = ref_set_name, maf = data.table(), excluded = data.table(),
             groups = sample_groups, mutrates = data.table(),
             selection_results = data.table(), ref_data_dir = data_dir,
             advanced = advanced, samples = data.table(), mutations = list())
  cesa = update_cesa_history(cesa, match.call())
  
  msg = paste0("This CESAnalysis will use ", ref_set_name, " reference data and the ", tolower(BSgenome::commonName(bsg)),
               " genome, assembly ", BSgenome::providerVersion(bsg), '.')
  pretty_message(msg)
  return(cesa)
}



#' Load a previously saved CESAnalysis object
#' 
#' @param file filename/path of CESAnalysis that has been saved with saveRDS, expected to end in .rds
#' @export
load_cesa = function(file) {
  if (! endsWith(file, '.rds')) {
    stop("Expected filename to end in .rds (because saveRDS() is the recommended way to save a CESAnalysis).", call. = F)
  }
  
  cesa = readRDS(file)
  
  if (! .hasSlot(cesa, "ref_data_dir")) {
    ## TEMPORARY
    cesa@ref_data_dir = c(ces_hg19_v1 = "/Users/Jeff/cancereffectsizeR/inst/ref_sets/ces_hg19_v1")
    cesa@ref_key = "ces_hg19_v1"
  }
  if(.hasSlot(cesa, "progressions")) {
    cesa@groups = cesa@progressions
    if(cesa@samples[, .N] > 0) {
      setnames(cesa@samples, c("progression_name", "progression_index"), c("group", "group_index"))
    }
  }
  
  available_ref_sets = get_ref_set_dirs()
  ref_key = cesa@ref_key
  if (! ref_key %in% names(available_ref_sets)) {
    if (! dir.exists(cesa@ref_data_dir)) {
      cesa@ref_data_dir = NA_character_
      msg = paste0("Reference data associated with the CESAnalysis (", ref_key, ") not found. You can view the data in this CESAnalysis, ",
                   "but many functions will not work as expected. If this is a custom reference data set, ",
                   "you can fix the issue by using set_ces_ref_set_dir() to associate the path to your data with the analyis.")
      warning(paste0(strwrap(msg), collapse = "\n"))
    } 
  } else {
    cesa@ref_data_dir = available_ref_sets[ref_key]
  }
  
  if (! is.na(cesa@ref_data_dir)) {
    preload_ref_data(cesa@ref_data_dir)
  }
  
  # Allow back-compatibility with column name changes
  if(! is.null(cesa@mutations$amino_acid_change)) {
    setnames(cesa@mutations$amino_acid_change, 'all_snv_ids', 'constituent_snvs', skip_absent = T)
  }
  if(! is.null(cesa@mutations$snv)) {
    setnames(cesa@mutations$snv, 'assoc_aa_mut', 'assoc_aac', skip_absent = T)
    setnames(cesa@maf, 'assoc_aa_mut', 'assoc_aac', skip_absent = T)
  }
  setnames(cesa@maf, c("Variant_Type", "snv_id"), c("variant_type", "variant_id"), skip_absent = T)
  cesa@maf[variant_type == "SNV", variant_type := "snv"]
  
  if (! .hasSlot(cesa, "run_history")) {
    cesa@run_history = character()
  }
  
  cesa = update_cesa_history(cesa, match.call())
  
  # temporary
  if (is.null(cesa@advanced$using_exome_plus)) {
    cesa@advanced$using_exome_plus = F
  }
  if (is.null(cesa@advanced$locked)) {
    # not generally desired behavior, but temporary
    cesa@advanced$locked = cesa@advanced$annotated
  }
  return(cesa)
}


#' View data loaded into CESAnalysis
#' 
#' returns a data.table containing MAF records used in the given CESAnalysis
#' @param cesa CESAnalysis object
#' @export
maf_records = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: maf_records(cesa), where cesa is a CESAnalysis")
  }
  if(cesa@maf[,.N] == 0) {
    stop("No MAF data has been loaded")
  }
  return(cesa@maf)
}

#' View excluded MAF data
#' 
#' returns a data.table containing MAF records that were excluded from the given CESAnalysis
#' @param cesa CESAnalysis object
#' @export
excluded_maf_records = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: excluded_maf_records(cesa), where cesa is a CESAnalysis")
  }
  if(cesa@maf[,.N] == 0) {
    stop("No MAF data has been loaded yet, so naturally no records have been excluded.")
  }
  if(cesa@excluded[,.N] == 0) {
    message("Returned an empty data table since no records have been excluded.")
  }
  return(cesa@excluded)
}

#' View sample metadata
#' 
#' returns a data.table with info on all samples in the CESAnalysis
#' @param cesa CESAnalysis object
#' @export
get_sample_info = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_sample_info(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@samples[,.N] == 0) {
    stop("No MAF data has been loaded yet, so naturally there is no sample data.")
  }
  
  # user doesn't need group columns for single-group analyses
  if(length(cesa@groups) == 1) {
    return(cesa@samples[, -c("group_index", "group")])
  } else {
    return(cesa@samples)
  }
}

#' Get expected relative trinucleotide-specific SNV mutation rates
#' 
#' @param cesa CESAnalysis object
#' @export
get_trinuc_rates = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_trinuc_rates(cesa), where cesa is a CESAnalysis")
  }
  return(as.data.table(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix, keep.rownames = "Unique_Patient_Identifier"))
}

#' Get table of signature weights by tumor
#' 
#' @param cesa CESAnalysis object
#' @param artifacts_zeroed If TRUE, return weights table in which artifact signatures have
#'   been set to zero and remaining weights renormalized to sum to 1, reflecting relative
#'   contributions of biological processes to mutations. If FALSE, return a table with
#'   artifact weights unaltered.
#' @export
get_signature_weights = function(cesa = NULL, artifacts_zeroed = T) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_signature_weights(cesa), where cesa is a CESAnalysis")
  }
  if (artifacts_zeroed) {
    return(cesa@trinucleotide_mutation_weights$signature_weight_table)
  } else {
    return(cesa@trinucleotide_mutation_weights$signature_weight_table_with_artifacts)
  }
}

#' Get table of neutral gene mutation rates
#' 
#' @param cesa CESAnalysis object
#' @export
get_gene_rates = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_gene_rates(cesa), where cesa is a CESAnalysis")
  }
  gene_rates = cesa@mutrates
  if (ncol(gene_rates) == 2) {
    colnames(gene_rates) = c("gene", "rate")
  }
  return(gene_rates)
}


#' View results from ces_snv
#' 
#' returns a data table of SNV effect sizes generated with ces_snv
#' @param cesa CESAnalysis object
#' @export
snv_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: snv_results(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@selection_results[, .N] == 0) {
    stop("No results yet from ces_snv in this CESAnalysis")
  }
  annotations = suppressMessages(select_variants(cesa, variant_ids = cesa@selection_results$variant_id, min_freq = 0))
  results = cesa@selection_results[annotations, on = c("variant_id", "variant_type")]
  results_cols = colnames(results)
  # try to flip variant_name and variant_id columns
  if(results_cols[1] == "variant_id") {
    name_col = which(results_cols == "variant_name")[1]
    setcolorder(results, c("variant_name", 
                           results_cols[2:(name_col - 1)], "variant_id", 
                           results_cols[(name_col + 1):length(results_cols)]))
  }
  return(results)
}





#' clean_granges_for_bsg
#' 
#' Tries to format an input GRanges object to be compatible with a CESAnalysis's reference
#' genome. Optionally also applies padding to start and end positions of ranges, stopping
#' at chromosome ends. Either stops with an error or returns a clean granges object.
#' 
#' @param bsg BSgenome object (CES-formatted, typically)
#' @param gr GRanges object
#' @param padding How many bases to expand start and end of each position
#' @keywords internal
clean_granges_for_bsg = function(bsg = NULL, gr = NULL, padding = 0) {
  stopifnot(is(padding, "numeric"),
            length(padding) == 1,
            padding >= 0,
            padding - as.integer(padding) == 0)

  # try to make gr style/seqlevels match bsg (if this fails, possibly the genome build does not match)
  GenomeInfoDb::seqlevelsStyle(gr) = "NCBI"
  
  tryCatch({
    msg = paste0("An input granges (or converted BED file) does't seem compatible with the current reference genome.\n",
                 "Make sure it uses the same genome assembly. It may also help to subset to just the\n",
                 "primary chromosomes, if any obscure contigs are present in your regions.\n",
                 "Original warning/error:")
    GenomeInfoDb::seqlevels(gr) = GenomeInfoDb::seqlevels(bsg)
    GenomeInfoDb::seqinfo(gr) = GenomeInfoDb::seqinfo(bsg)
  }, error = function(e) {
    message(msg)
    stop(conditionMessage(e))
  }, warning = function(w) {
    message(msg)
    stop(conditionMessage(w))
  })
  
  # drop any metadata
  GenomicRanges::mcols(gr) = NULL
  
  # sort, reduce, unstrand
  gr = GenomicRanges::reduce(GenomicRanges::sort(gr), drop.empty.ranges = T)
  GenomicRanges::strand(gr) = "*"
  
  # require genome name to match the reference genome (too many potential errors if we allow anonymous or mismatched genome)
  expected_genome = GenomeInfoDb::genome(bsg)[1]
  gr_genome = GenomeInfoDb::genome(gr)[1]
  if (expected_genome != gr_genome) {
    stop(paste0("The genome name of an input granges object (", gr_genome, ") does not match the current reference genome (",
                expected_genome, ")."))
  }
  
  
  if (padding > 0) {
    # Suppress the out-of-range warning since we'll trim afterwards
    withCallingHandlers(
      {
        GenomicRanges::start(gr) = GenomicRanges::start(gr) - padding
        GenomicRanges::end(gr) = GenomicRanges::end(gr) + padding
      }, warning = function(w) 
      {
        if (grepl("out-of-bound range", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
    gr = GenomicRanges::reduce(GenomicRanges::trim(gr))
  }
  return(gr)
}

update_cesa_history = function(cesa, comm) {
  if (identical(cesa@advanced$recording, TRUE)) {
    cesa@run_history =  c(cesa@run_history, deparse(comm, width.cutoff = 500))
  }
  return(cesa)
}
