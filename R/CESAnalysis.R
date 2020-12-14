#' Create a cancereffectsizeR analysis
#' 
#' Creates a CESAnalysis object, the central data structure of cancereffectsizeR.
#' 
#' @param refset name of reference data set (refset) to use; run \code{list_ces_refsets()} for
#'   available refsets. Alternatively, the path to a custom reference data directory.
#' @param sample_groups Optionally, supply labels identifying different groups of samples.
#'   Each designated group of samples can be run independently through functions like
#'   \code{trinuc_mutation_rates()} and \code{gene_mutation_rates()}, and some selection
#'   models (such as sswm_sequential) require multiple sample groups.
#' @return CESAnalysis object
#' @export
CESAnalysis = function(refset = "ces.refset.hg19", sample_groups = NULL) {
  refset_name = refset
  if (! is(refset_name, "character")) {
    stop("refset should be type character (either the name of a supported reference data set package or a path to a custom reference data directory.")
  }
  
  # Check for and load reference data for the chosen genome/transcriptome data
  refset_name = refset
  if (refset_name %in% names(.official_refsets)) {
    if(file.exists(refset_name)) {
      stop("You've given the name of a CES reference data set package, but a file/folder with the same name is in your working directory. Stopping to avoid confusion.")
    }
    if(! require(refset_name, character.only = T, quietly = T)) {
      if(refset_name == "ces.refset.hg19") {
        message("Install ces.refset.hg19 like this:\n",
                "remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg19@*release\")")
      }
      stop("CES reference data set ", refset_name, " not installed.")
    }
    req_version = .official_refsets[[refset_name]]
    actual_version = packageVersion(refset_name)
    if (actual_version < req_version) {
      stop("CES reference data set ", refset_name, " is version ", actual_version, ", but your version of cancereffectsizeR requires at least ",
           "version ", req_version, ".\nRun this to update:\n",
           "remotes::install_github(\"Townsend-Lab-Yale/ces-reference-data/", refset_name, "\")")
    }
    ref_data_version = actual_version
    data_dir = system.file("refset", package = refset_name)
  } else {
    if (! dir.exists(refset_name)) {
      if (grepl('/', refset_name)) {
        stop("Could not find reference data at ", refset_name)
      } else {
        stop("Invalid reference set name. Check spelling, or view available data sets with list_ces_refsets().")
      }
    }
    
    data_dir = refset_name
    refset_name = basename(refset_name)
    if(refset_name %in% names(.official_refsets)) {
      stop("Your custom reference data set has the same name (", refset_name, ") as a CES reference data package. Please rename it.")
    }
  }
  
  # load reference data
  if (! refset_name %in% ls(.ces_ref_data)) {
    message("Loading reference data set ", refset_name, "...")
    .ces_ref_data[[refset_name]] = preload_ref_data(data_dir)
    .ces_ref_data[[refset_name]][["data_dir"]] = data_dir
  }
  
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
  
  # Enforce format for sample group names
  # Among other things, disallowing commas/whitespace will allow comma-delimiting groups within data table columns
  # Unliked covered_regions_name, can start with letter or number (not just letter)
  legal_name = '^[a-z0-9][0-9a-z\\_\\-\\.]*$'
  if (! all(grepl(legal_name, tolower(sample_groups), perl = T))) {
    stop("Invalid sample group names. Start with a letter/number and use only letters, numbers, '-', '_', '.'.")
  }
    
  # advanced is a grab bag of additional stuff to keep track of
  ## annotated: whether loaded MAF records are annotated
  ## using_exome_plus: whether previously loaded and any future generic exome data uses the "exome+" coverage option 
  ##  (either all generic data must, or none of it, based on choice of enforce_generic_exome_coverage on first load_maf call)
  ## recording: whether "run_history" is currently being recorded (gets set to false during some internal steps for clarity)
  ## locked: whether load_maf can still be used (can't load more data after trinuc_mutation_rates or gene_mutation_rates)
  ## trinuc_done: have all trinuc mutation rates been calculated?
  ## gene_rates_done: have all samples been through gene_mutation_rates?
  ## uid: a unique-enough identifier for the CESAnalysis (just uses epoch time)
  ## genome_info: environment with stuff like genome build name, species, name of associated BSgenome
  genome_info = get_ref_data(data_dir, "genome_build_info")
  ces_version = packageVersion("cancereffectsizeR")
  advanced = list("version" = ces_version, annotated = F, using_exome_plus = F, 
                  recording = T, locked = F, trinuc_done = F, gene_rates_done = F,
                  uid = unclass(Sys.time()), genome_info = genome_info)
  cesa = new("CESAnalysis", run_history = character(),  ref_key = refset_name, maf = data.table(), excluded = data.table(),
             groups = sample_groups, mutrates = data.table(),
             selection_results = list(), ref_data_dir = data_dir, epistasis = list(),
             advanced = advanced, samples = data.table(), mutations = list())
  cesa@run_history = c(paste0("[Version: cancereffectsizeR ", ces_version, "]" ))
  cesa = update_cesa_history(cesa, match.call())

  msg = paste0("This CESAnalysis will use ", refset_name, " reference data and the ", genome_info$species,
               " genome, assembly ", genome_info$build_name, '.')
  pretty_message(msg)
  return(cesa)
}

#' Save a CESAnalysis in progress
#' 
#' Saves a CESAnalysis to a file by calling using base R's saveRDS function. Also updates
#' run history for reproducibility. Files saved should be reloaded with \code{load_cesa()}.
#' 
#' Note that the genome reference data associated with a CESAnalysis (refset) is not
#' actually part of the CESAnalysis, so it is not saved here. (Saving this data with the
#' analysis would make file sizes too large.) When you reload the CESAnalysis, you can
#' re-associate the correct reference data.
#' @param cesa CESAnalysis to save
#' @param file filename to save to (must end in .rds)
#' @export 
save_cesa = function(cesa, file) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  if(! is(file, "character") || length(file) != 1) {
    stop("file should be a valid file path (1-length character)")
  }
  if(! grepl('\\.rds$', file)) {
    stop("filename should end in .rds (indicates R data serialized format)")
  }
  cesa_to_save = update_cesa_history(cesa, match.call())
  saveRDS(cesa_to_save, file)
}

#' Load a previously saved CESAnalysis
#' 
#' @param file filename/path of CESAnalysis that has been saved with saveRDS, expected to end in .rds
#' @export
load_cesa = function(file) {
  if (! endsWith(tolower(file), '.rds')) {
    stop("Expected filename to end in .rds (because saveRDS() is the recommended way to save a CESAnalysis).", call. = F)
  }
  cesa = readRDS(file)
  
  # Data tables must be reset to get them working properly
  cesa@samples = setDT(cesa@samples)
  cesa@maf = setDT(cesa@maf)
  if (! is.null(cesa@mutations$amino_acid_change)) {
    cesa@mutations$amino_acid_change = setDT(cesa@mutations$amino_acid_change, key = "aac_id")
  }
  if (! is.null(cesa@mutations$snv)) {
    cesa@mutations$snv = setDT(cesa@mutations$snv, key = "snv_id")
  }
  if (! is.null(cesa@trinucleotide_mutation_weights[["signature_weight_table"]])) {
    cesa@trinucleotide_mutation_weights[["signature_weight_table"]] = setDT(cesa@trinucleotide_mutation_weights[["signature_weight_table"]])
  }
  if (! is.null(cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]])) {
    cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]] = setDT(cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]])
  }
  cesa@mutrates = setDT(cesa@mutrates)
  cesa@selection_results = lapply(cesa@selection_results, setDT)
  cesa@epistasis = lapply(cesa@epistasis, setDT)
  
  refset_name = cesa@ref_key
  if (refset_name %in% names(.official_refsets)) {
    if(! require(refset_name, character.only = T, quietly = T)) {
      stop("CES reference data set ", refset_name, " not installed. Run this to install:\n", 
           "remotes::install_github(\"Townsend-Lab-Yale/ces-reference-data/", refset_name, "\")")
    }
    req_version = .official_refsets[[refset_name]]
    actual_version = packageVersion(refset_name)
    if (actual_version < req_version) {
      stop("CES reference data set ", refset_name, " is version ", actual_version, ", but your version of cancereffectsizeR requires at least ",
           "version ", req_version, ".\nRun this to update:\n",
           "remotes::install_github(\"Townsend-Lab-Yale/ces-reference-data/", refset_name, "\")")
    }
    cesa@ref_data_dir = system.file("refset", package = refset_name)
  } else {
    if (! dir.exists(cesa@ref_data_dir)) {
      cesa@ref_data_dir = NA_character_
      msg = paste0("Reference data not found at ", cesa@ref_data_dir, ". You can view the data in this CESAnalysis, ",
                   "but many functions will not work as expected. If this is a custom reference data set, ",
                   "you can fix the issue by using set_refset_dir() to associate the path to your data with the analyis.")
      warning(pretty_message(msg, emit = F))
    }
  }
  
  if (! is.na(cesa@ref_data_dir) & ! cesa@ref_key %in% ls(.ces_ref_data)) {
    .ces_ref_data[[cesa@ref_key]] = preload_ref_data(cesa@ref_data_dir)
  }
  current_version = packageVersion("cancereffectsizeR")
  previous_version = cesa@advanced$version
  if (as.character(current_version) != as.character(previous_version)) {
    warning("Version change: CESAnalysis previously created in CES ", previous_version, ".\n",
            "Currently running version ", current_version, '.', call. = F)
    cesa@advanced$version = paste(previous_version, current_version, sep = '/' )
    cesa@run_history = c(cesa@run_history, paste0("\n[Now running CES ", current_version, ']'))
  }
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}


#' Set reference data directory
#'
#' When working with custom reference data or loading a previously saved CESAnalysis in a
#' new environment, use this function to reassociate the location of reference data with
#' the analysis. (If \code{load_cesa()} didn't give you a warning when loading your
#' analysis, you probably don't need to use this function.)
#' 
#' @param cesa CESAnalysis
#' @param dir path to data directory
#' @export
set_refset_dir = function(cesa, dir) {
  if(cesa@ref_key %in% names(.official_refsets)) {
    stop("You can't set the reference directory on a built-in CES reference data set.")
  }
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  if(! is(dir, "character") || length(dir) != 1) {
    stop("dir should be a path to a directory")
  }
  if(! dir.exists(dir)) {
    stop("Directory not found at ", dir)
  }
  dir = normalizePath(dir)
  
  cesa@ref_data_dir = dir
  .ces_ref_data[[cesa@ref_key]] = preload_ref_data(dir)
  .ces_ref_data[[cesa@ref_key]][["data_dir"]] = dir
  cesa = update_cesa_history(cesa, match.call())
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
    return(cesa@samples[, -c("group")])
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


#' View results from ces_variant
#' 
#' returns a list of data tables with results from ces_variant(), plus annotations
#' @param cesa CESAnalysis object
#' @export
snv_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: snv_results(cesa), where cesa is a CESAnalysis")
  }
  if (length(cesa@selection_results) == 0) {
    stop("No results yet from ces_variant in this CESAnalysis")
  }
  
  output = list()
  for (i in 1:length(cesa@selection_results)) {
    curr_selection = cesa@selection_results[[i]]
    
    # compound variant runs have already have all available annotations
    if (identical(attr(curr_selection, "is_compound"), TRUE)) {
      output = c(output, list(curr_selection))
      next
    }
    annotations = suppressMessages(select_variants(cesa, variant_passlist = curr_selection$variant_id))
    results = curr_selection[annotations, on = c("variant_id", "variant_type")]
    results_cols = colnames(results)
    # try to flip variant_name and variant_id columns
    if(results_cols[1] == "variant_id") {
      name_col = which(results_cols == "variant_name")[1]
      setcolorder(results, c("variant_name", 
                             results_cols[2:(name_col - 1)], "variant_id", 
                             results_cols[(name_col + 1):length(results_cols)]))
    }
    setattr(results, "si_cols", attr(curr_selection, "si_cols", T))
    
    output = c(output, list(results))
  }
  names(output) = names(cesa@selection_results)
  
  return(output)
}

#' View output from epistasis functions
#' 
#' returns a list of data tables with results from epistasis functions
#' @param cesa CESAnalysis object
#' @export
epistasis_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: epistasis_results(cesa), where cesa is a CESAnalysis")
  }
  if (length(cesa@epistasis) == 0) {
    stop("No results yet from epistasis functions in this CESAnalysis")
  }
  return(cesa@epistasis)
}



#' clean_granges_for_cesa
#' 
#' Tries to format an input GRanges object to be compatible with a CESAnalysis's reference
#' genome. Optionally also applies padding to start and end positions of ranges, stopping
#' at chromosome ends. Either stops with an error or returns a clean granges object.
#' 
#' @param cesa CESAnalysis
#' @param gr GRanges object
#' @param padding How many bases to expand start and end of each position
#' @keywords internal
clean_granges_for_cesa = function(cesa = NULL, gr = NULL, padding = 0) {
  bsg = get_cesa_bsg(cesa)
  
  stopifnot(is(padding, "numeric"),
            length(padding) == 1,
            padding >= 0,
            padding - as.integer(padding) == 0)

  # try to make gr style/seqlevels match bsg (if this fails, possibly the genome build does not match)
  # For now, suppressing "more than one best sequence renaming map"; tends to appear on single-chr inputs
  withCallingHandlers(
    { 
      GenomeInfoDb::seqlevelsStyle(gr) = "NCBI" 
    },
    warning = function(w) {
      if (grepl("more than one best sequence renaming map", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }        
  )

  
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
  
  # subset to just supported contigs
  supported_chr = cesa@advanced$genome_info$supported_chr
  gr = gr[as.character(GenomeInfoDb::seqnames(gr)) %in% supported_chr]
  gr = GenomeInfoDb::keepSeqlevels(gr, supported_chr)
  
  
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
