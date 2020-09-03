#' Create a cancereffectsizeR analysis
#' 
#' Creates a CESAnalysis object, the central data structure of cancereffectsizeR.
#' @param ref_set reference data to use (e.g. "ces_hg19_v1")
#' @param progression_order evolutionary order of tumor stage progression (e.g. c("Primary", "Metastatic"))
#' @return CESAnalysis object
#' @export
CESAnalysis = function(ref_set = "ces_hg19_v1", progression_order = NULL) {
  
  # Check for and load reference data for the chosen genome/transcriptome data
  ref_key = ref_set
  preload_ref_data(ref_key)
  ref_data_dir = get_ref_set_dirs()[ref_key]
  bsg = .ces_ref_data[[ref_key]]$genome
  message(crayon::black(paste0("Okay, this CES analysis will use the ", 
                               tolower(BSgenome::commonName(bsg)), " genome (", BSgenome::releaseName(bsg), ").")))
  
  
  # Validate progression_order
  if (is.null(progression_order)) {
    progression_order = c("stageless")
  }
  
  if (is.numeric(progression_order)) {
    progression_order = as.character(progression_order)
  }

  if (! is(progression_order, "character")) {
    stop("progression_order should be a character vector of chronological tumor states (e.g., Primary, Metastatic)")
  }
  
  
  # advanced is a grab bag of additional stuff to keep track of
  ## annotated: whether loaded MAF records are annotated
  ## using_exome_plus: whether previously loaded and any future generic exome data uses the "exome+" coverage option 
  ##  (either all generic data must, or none of it, based on choice of enforce_generic_exome_coverage on first load_maf call)
  ## recording: whether "run_history" is currently being recorded (gets set to false during some internal steps for clarity)
  advanced = list("version" = packageVersion("cancereffectsizeR"), annotated = F, using_exome_plus = F, recording = T)
  cesa = new("CESAnalysis", run_history = character(),  ref_key = ref_key, maf = data.table(), excluded = data.table(),
             progressions = progression_order, mutrates = data.table(),
             gene_epistasis_results = data.table(), selection_results = data.table(), ref_data_dir = ref_data_dir,
             advanced = advanced, samples = data.table(), mutations = list())
  cesa = update_cesa_history(cesa, match.call())
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
  }
  
  
  ref_key = names(cesa@ref_data_dir)[1]
  
  available_ref_sets = get_ref_set_dirs()
  if (! ref_key %in% names(available_ref_sets)) {
    warning("Reference data for ", ref_key, " not found. You can view data in this CESAnalysis,\n",
            "but most functions will not work as expected.")
    return(cesa)
  }
  
  cesa@ref_data_dir = available_ref_sets[ref_key]
  preload_ref_data(ref_key)
  cesa@ref_key = ref_key
  
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
  
  # user doesn't need progression columns for single-progression-state analyses
  if(length(cesa@progressions) == 1) {
    return(cesa@samples[, -c("progression_index", "progression_name")])
  } else {
    return(cesa@samples)
  }
  return(cesa@samples)
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
#' @export
get_signature_weights = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_signature_weights(cesa), where cesa is a CESAnalysis")
  }
  return(cesa@trinucleotide_mutation_weights$signature_weight_table)
}

#' Get table of neutral gene mutation rates by progression state
#' 
#' @param cesa CESAnalysis object
#' @export
get_gene_rates = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_gene_rates(cesa), where cesa is a CESAnalysis")
  }
  return(cesa@mutrates)
}

#' Get lists of mutations and annotations
#' 
#' @param cesa CESAnalysis object
#' @export
get_mutations = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_mutations(cesa), where cesa is a CESAnalysis")
  }
  return(cesa@mutations)
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
  return(cesa@selection_results)
}

#' View results from gene-level epistasis analysis
#' 
#' returns a data table of pairwise gene epistasis effect sizes generated with ces_gene_epistasis
#' @param cesa CESAnalysis object
#' @export
gene_epistasis_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: gene_epistasis_results(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@gene_epistasis_results[, .N] == 0) {
    stop("No results yet from ces_gene_epistasis in this CESAnalysis")
  }
  return(cesa@gene_epistasis_results)
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
