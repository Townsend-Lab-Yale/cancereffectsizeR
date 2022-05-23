#' Create a cancereffectsizeR analysis
#' 
#' Creates a CESAnalysis, the central data structure of cancereffectsizeR.
#' 
#' @param refset Name of reference data set (refset) to use; run \code{list_ces_refsets()} for
#'   available refsets. Alternatively, the path to a custom reference data directory.
#' @param sample_groups (Deprecated; no longer necessary.) Optionally, supply labels
#'   identifying different groups of samples. Each designated group of samples can be run
#'   independently through functions like \code{trinuc_mutation_rates()} and
#'   \code{gene_mutation_rates()}, and some selection models (such as sswm_sequential)
#'   require multiple sample groups.
#' @return CESAnalysis object
#' @export
CESAnalysis = function(refset = NULL, sample_groups = NULL) {
  if(is.null(refset)) {
    msg = paste0("Required argument refset: Supply a reference data package (e.g., ces.refset.hg38 or ces.refset.hg19), or ",
                 "custom reference data (see docs).")
    stop(paste0(strwrap(msg, exdent = 2), collapse = "\n"))
  }
  if(is(refset, "environment")) {
    refset_name = as.character(substitute(refset))
  } else {
    refset_name = refset
  }
  # Check for and load reference data for the chosen genome/transcriptome data
  if (! is(refset_name, "character")) {
    stop("refset should be a refset object, the name of an installed refset package, or a path to a custom refset directory.")
  }
  using_custom_refset = TRUE
  if (refset_name %in% names(.official_refsets)) {
    using_custom_refset = FALSE
    if(file.exists(refset_name)) {
      stop("You've given the name of a CES reference data set package, but a file/folder with the same name is in your working directory. Stopping to avoid confusion.")
    }
    if(! require(refset_name, character.only = T)) {
      if(refset_name == "ces.refset.hg19") {
        message("Install ces.refset.hg19 like this:\n",
                "options(timeout = 600)\n",
                "remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg19@*release\")")
      } else if(refset_name == "ces.refset.hg38") {
        message("Install ces.refset.hg38 like this:\n",
                "options(timeout = 600)\n",
                "remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg38@*release\")")
      }
      stop("CES reference data set ", refset_name, " not installed.")
    }
    req_version = .official_refsets[[refset_name]]
    actual_version = packageVersion(refset_name)
    if (actual_version < req_version) {
      stop("CES reference data set ", refset_name, " is version ", actual_version, ", but your version of cancereffectsizeR requires at least ",
           "version ", req_version, ".\nRun this to update:\n",
           "remotes::install_github(\"Townsend-Lab-Yale/", refset_name, "\")")
    }
    refset_version = actual_version
    data_dir = system.file("refset", package = refset_name)
  } else {
    refset_version = NA_character_
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
    if (using_custom_refset) {
      .ces_ref_data[[refset_name]] = preload_ref_data(data_dir)
    } else {
      .ces_ref_data[[refset_name]] = get(refset_name, envir = as.environment(paste0('package:', refset_name)))
    }
    .ces_ref_data[[refset_name]][["data_dir"]] = data_dir
  }
  
  # Validate sample_groups
  if (is.null(sample_groups)) {
    sample_groups = c("stageless")
  } else {
    warning("sample_groups is deprecated and will be removed in a future update. No downstream functions require declaration of sample_groups.")
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
  ## using_exome_plus: whether previously loaded and any future generic exome data uses the "exome+" coverage option 
  ##  (either all generic data must, or none of it, based on choice of enforce_default_exome_coverage on first load_maf call)
  ## recording: whether "run_history" is currently being recorded (gets set to false during some internal steps for clarity)
  ## uid: a unique-enough identifier for the CESAnalysis (just uses epoch time)
  ## genome_info: environment with stuff like genome build name, species, name of associated BSgenome
  ## snv_signatures: List CES signature sets used in the analysis
  ## cached_variants (not populated here): output of select_variants() run with default arguments
  ##      (automatically updated as needed by load_cesa/update_covered_in)
  genome_info = get_ref_data(data_dir, "genome_build_info")
  ces_version = packageVersion("cancereffectsizeR")
  advanced = list("version" = ces_version, using_exome_plus = F, 
                  recording = T, uid = unclass(Sys.time()), genome_info = genome_info,
                  snv_signatures = list(), refset_version = refset_version)
  
  # Mutation table specifications (see template tables declared in imports.R)
  mutation_tables = list(amino_acid_change = copy(aac_annotation_template), 
                         snv = copy(snv_annotation_template), aac_snv_key = copy(aac_snv_key_template))
  
  cesa = new("CESAnalysis", run_history = character(),  ref_key = refset_name, maf = data.table(), excluded = data.table(),
             groups = sample_groups, mutrates = data.table(),
             selection_results = list(), ref_data_dir = data_dir, epistasis = list(),
             advanced = advanced, samples = copy(sample_table_template), mutations = mutation_tables,
             coverage = list())

  cesa@run_history = c(paste0("[Version: cancereffectsizeR ", ces_version, "]" ))
  if (! is.na(refset_version)) {
    cesa@run_history= c(cesa@run_history, 
                        paste0("[Refset: ", refset_name, " ", refset_version, "]"))
  }
  cesa = update_cesa_history(cesa, match.call())

  msg = paste0("This CESAnalysis will use ", refset_name, " reference data and the ", genome_info$species,
               " genome, assembly ", genome_info$build_name, '.')
  pretty_message(msg)
  return(cesa)
}

#' Create an independent copy of a CESAnalysis
#' 
#' Used internally to "copy" CESAnalysis objects while keeping memory use to a minimum.
#'
#' The trick is to use data.table's copy function on all data.tables (and lists of
#' data.tables) within CESAnalysis slots. (If you just call copy on the whole object, the
#' data tables won't be handled in a memory-efficient way. And if you call copy on
#' non-data.tables, it's actually less efficient since it forces an immediate full copy
#' instead of the usual copy-on-modify.)
#' 
#' @param cesa CESAnalysis
#' @keywords internal
copy_cesa = function(cesa) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  trinuc_orig = cesa@trinucleotide_mutation_weights
  trinuc_copy = list()
  if(length(trinuc_orig) > 0) {
    for (i in c("trinuc_snv_counts", "trinuc_proportion_matrix", "group_average_dS_output")) {
      trinuc_copy[[i]] = trinuc_orig[[i]]
    }
    trinuc_copy[c("signature_weight_table", "signature_weight_table_with_artifacts",
                  "raw_signature_weights")] = copy(trinuc_orig[c("signature_weight_table", "signature_weight_table_with_artifacts",
                                                                 "raw_signature_weights")])
  }

  
  cesa = new("CESAnalysis", run_history = cesa@run_history,  groups = cesa@groups,
             ref_key = cesa@ref_key, ref_data_dir = cesa@ref_data_dir, dndscv_out_list = copy(cesa@dndscv_out_list),
             maf = copy(cesa@maf), excluded = copy(cesa@excluded),
             mutrates = copy(cesa@mutrates), selection_results = copy(cesa@selection_results),
             epistasis = copy(cesa@epistasis), samples = copy(cesa@samples), mutations = copy(cesa@mutations),
             advanced = copy(cesa@advanced), coverage = copy(cesa@coverage),
             trinucleotide_mutation_weights = trinuc_copy)
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
  # Note updating of history before copying, so save_cesa call is recorded both in the
  # user's active analysis and in the saved file.
  cesa_to_save = update_cesa_history(cesa, match.call())
  cesa_to_save = copy_cesa(cesa)
  cesa_to_save@advanced$cached_variants = NULL # reduce file size (data recalculated on reload)
  saveRDS(cesa_to_save, file)
}

#' Load a previously saved CESAnalysis
#' 
#' @param file filename/path of CESAnalysis that has been saved with saveRDS, expected to end in .rds
#' @export
load_cesa = function(file) {
  if (! endsWith(tolower(file), '.rds')) {
    stop("Expected filename to end in .rds.", call. = F)
  }
  cesa = readRDS(file)
  
  # Data tables must be reset to get them working properly
  cesa@samples = setDT(cesa@samples)
  
  # Convert old versions of sample table
  if(cesa@samples[, .N] == 0) {
    cesa@samples = copy(sample_table_template)
  }
  if(is.null(cesa@samples$sig_analysis_grp)) {
    cesa@samples[, sig_analysis_grp := NA_integer_]
    if(! is.null(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix)) {
      cesa@samples[rownames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix), sig_analysis_grp := 0]
    }
  }
  if(is.null(cesa@samples$gene_rate_grp)) {
    cesa@samples[, gene_rate_grp := NA_integer_]
  } else if(is.character(cesa@samples$gene_rate_grp)) {
    # previously there was rate_grp_1, etc. here (now just 1, 2, ...)
    cesa@samples[, gene_rate_grp := as.integer(sub('.*_', '', gene_rate_grp))]
  }
  cesa@maf = setDT(cesa@maf)
  
  # Older versions lack annotation table templates
  if (is.null(cesa@mutations$snv)) {
    cesa@mutations = list(amino_acid_change = aac_annotation_template, snv = snv_annotation_template)
  } else {
    cesa@mutations$amino_acid_change = setDT(cesa@mutations$amino_acid_change, key = "aac_id")
    cesa@mutations$snv = setDT(cesa@mutations$snv, key = "snv_id")
    if (! is.null(cesa@mutations$aac_snv_key)) {
      # if it is NULL, gets handled later
      cesa@mutations$aac_snv_key = setDT(cesa@mutations$aac_snv_key, key = "aac_id")
    }
  }
  

  if (! is.null(cesa@trinucleotide_mutation_weights[["signature_weight_table"]])) {
    cesa@trinucleotide_mutation_weights[["signature_weight_table"]] = setDT(cesa@trinucleotide_mutation_weights[["signature_weight_table"]])
  }
  if (! is.null(cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]])) {
    cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]] = setDT(cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]])
  }
  
  if (! is.null(cesa@trinucleotide_mutation_weights[["raw_signature_weights"]])) {
    cesa@trinucleotide_mutation_weights[["raw_signature_weights"]] = setDT(cesa@trinucleotide_mutation_weights[["raw_signature_weights"]])
  }
  
  cesa@mutrates = setDT(cesa@mutrates)
  cesa@selection_results = lapply(cesa@selection_results, setDT)
  cesa@epistasis = lapply(cesa@epistasis, setDT)
  
  # Now a list of signature sets. Formerly just 1 set, so put in enclosing list if necessary.
  used_sig_sets = cesa@advanced$snv_signatures
  if (! is.null(used_sig_sets) && length(used_sig_sets) > 0) {
    if (! is.list(used_sig_sets[[1]])) {
      cesa@advanced$snv_signatures = list(cesa@advanced$snv_signatures)
      names(cesa@advanced$snv_signatures) = cesa@advanced$snv_signatures[[1]]$name
    }
    
    # Get each signature set's meta data.table and call setDT
    lapply(lapply(cesa@advanced$snv_signatures, '[[', 'meta'), setDT)
  } else {
    cesa@advanced$snv_signatures = list()
  }
  
  refset_name = cesa@ref_key
  if (refset_name %in% names(.official_refsets)) {
    if(! require(refset_name, character.only = T, quietly = T)) {
      stop("CES reference data set ", refset_name, " not installed. Run this to install:\n", 
           "remotes::install_github(\"Townsend-Lab-Yale/", refset_name, "\")")
    }
    req_version = .official_refsets[[refset_name]]
    actual_version = packageVersion(refset_name)
    if (actual_version < req_version) {
      stop("CES reference data set ", refset_name, " is version ", actual_version, ", but your version of cancereffectsizeR requires at least ",
           "version ", req_version, ".\nRun this to update:\n",
           "remotes::install_github(\"Townsend-Lab-Yale/", refset_name, "\")")
    }
    previous_version = cesa@advanced$refset_version
    if (is.null(previous_version) && refset_name == 'ces.refset.hg38') {
      msg = paste0("This CESAnalysis was likely created with an older version of ces.refset.hg38 with known issues. ",
              "You should create a new CESAnalysis and start over if you want to continue analysis.")
      warning(pretty_message(msg, emit = F))
    } else if(! is.null(previous_version)) {
      previous_version = as.package_version(previous_version)
      if (previous_version$major != actual_version$major | previous_version$minor != actual_version$minor) {
        msg = paste0("This CESAnalysis was annotated with data from ", refset_name, ' ', previous_version, " and you currently have ",
                     "version ", actual_version, ". You should create a new CESAnalysis and start over if you want to continue analysis.")
        warning(pretty_message(msg, emit = F))
      }
    }
    cesa@ref_data_dir = system.file("refset", package = refset_name)
  } else {
    if (! dir.exists(cesa@ref_data_dir)) {
      cesa@ref_data_dir = NA_character_
      msg = paste0("Reference data not found at ", cesa@ref_data_dir, ". You can view the data in this CESAnalysis, ",
                   "but many functions will not work as expected. If this is a custom reference data set, ",
                   "you can fix the issue by using set_refset_dir() to associate the path to your data with the analysis.")
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
  
  # Before v2.2, nearest_pid not annotated; add it if needed
  # This would crash on custom refsets with missing data directories, but unlikely
  # situation will ever arise with a pre-2.2 refset
  if (! "nearest_pid" %in% names(cesa@mutations$snv)) {
    snv_and_gene = cesa@mutations$snv[, .(gene = unlist(genes)), by = "snv_id"]
    to_lookup = snv_and_gene[, .(gene = unique(gene))]
    RefCDS = get_ref_data(cesa, "RefCDS")
    to_lookup[, pid := sapply(RefCDS[gene], '[[', 'protein_id')]
    snv_and_gene[to_lookup, pid := pid, on = 'gene']
    snv_and_pid = snv_and_gene[, .(nearest_pid = list(pid)), by = "snv_id"]
    cesa@mutations$snv[snv_and_pid, nearest_pid := nearest_pid, on = "snv_id"]
  }
  
  if (is.null(cesa@mutations$aac_snv_key)) {
    if (cesa@mutations$amino_acid_change[, .N] == 0) {
      cesa@mutations$aac_snv_key = copy(aac_snv_key_template)
    } else {
      aac_snv_key = cesa@mutations$snv[, .(aac_id = unlist(assoc_aac)), by = 'snv_id'][! is.na(aac_id)]
      setcolorder(aac_snv_key, c('aac_id', 'snv_id'))
      aac_snv_key[, multi_anno_site := uniqueN(aac_id) > 1, by = 'snv_id']
      setkey(aac_snv_key, 'aac_id')
      cesa@mutations$aac_snv_key = aac_snv_key
      cesa@mutations$snv[, assoc_aac := NULL]
      cesa@maf[, assoc_aac := NULL]
    }
  }

  # cache variant table for easy user access
  if(length(cesa@mutations$snv[, .N]) > 0) {
    cesa@advanced$cached_variants = select_variants(cesa)
    
    # Annotate top coding implications of each variant, based on select_variants() tiebreakers
    consequences = cesa@advanced$cached_variants[variant_type != 'snv', .(snv_id = unlist(constituent_snvs), variant_name, gene), by = 'variant_id']
    cesa@maf[consequences, c("top_gene", "top_consequence") := list(gene, variant_name), on = c(variant_id = 'snv_id')]
  }
  
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
  return(copy(cesa@maf))
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
  return(copy(cesa@excluded))
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
  
  to_include = names(cesa@samples)
  for(col in c("sig_analysis_grp", "gene_rate_grp")) {
    if(all(is.na(cesa@samples[[col]]))) {
      to_include = setdiff(to_include, col)
    }
  }
  if(length(cesa@groups) == 1) {
    to_include = setdiff(to_include, "group")
  }
  return(copy(cesa@samples[, .SD, .SDcols = to_include]))
  
}

#' Get estimated relative rates of trinucleotide-specific SNV mutation
#' 
#' @param cesa CESAnalysis object
#' @export
get_trinuc_rates = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_trinuc_rates(cesa), where cesa is a CESAnalysis")
  }
  return(as.data.table(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix, keep.rownames = "Unique_Patient_Identifier"))
}

#' Get table of signature attributions
#' 
#' View SNV signature attributions associated with CESAnalysis samples.
#' 
#' 
#' Use raw = TRUE to get signature attributions as produced by the signature extraction
#' tool (or as provided by the user with set_signature_weights()), without any of the
#' adjustments that are made by cancereffectsizeR's trinuc_mutation_rates().
#' 
#' @param cesa CESAnalysis object
#' @param raw Default FALSE. When TRUE, return raw signature attributions as found by the
#'   signature extraction tool. Format may vary by tool.
#' @param artifacts_zeroed Deprecated.
#' @return A data.table of signature attributions for each sample. By default, these are
#'   estimated relative weights of biologically-associated signatures (i.e., non-artifact
#'   signatures) that sum to 1.
#' @export
get_signature_weights = function(cesa = NULL, raw = F, artifacts_zeroed = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_signature_weights(cesa), where cesa is a CESAnalysis")
  }
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  
  if(! is.logical(raw) || length(raw) != 1) {
    stop("raw should be T/F.")
  }
  
  if(! is.null(artifacts_zeroed)) {
    stop("artifacts_zeroed is deprecated. Consider using raw = TRUE.")
  }

  if (raw == TRUE) {
    return(copy(cesa@trinucleotide_mutation_weights$raw_signature_weights))
  } else {
    return(copy(cesa@trinucleotide_mutation_weights$signature_weight_table))
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
  gene_rates = copy(cesa@mutrates)
  if (cesa@samples[, identical(unique(gene_rate_grp), 1L)]) {
    setnames(gene_rates, 'rate_grp_1', 'rate')
  }
  return(gene_rates)
}


#' View results from ces_variant
#' 
#' returns a list of ces_variant() results tables, with variant annotations added
#' 
#' @param cesa CESAnalysis object
#' @export
snv_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: snv_results(cesa), where cesa is a CESAnalysis")
  }

  return(copy(cesa@selection_results))
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
  return(copy(cesa@epistasis))
}



#' clean_granges_for_cesa
#' 
#' Tries to format an input GRanges object to be compatible with a CESAnalysis reference
#' genome. Optionally, also applies padding to start and end positions of ranges, stopping
#' at chromosome ends. Either stops with an error or returns a clean granges object.
#' 
#' @param cesa CESAnalysis
#' @param gr GRanges object
#' @param padding How many bases to expand start and end of each position
#' @keywords internal
clean_granges_for_cesa = function(cesa = NULL, gr = NULL, padding = 0, refset_env = NULL) {
  if(is.null(cesa)) {
    bsg = refset_env$genome
    supported_chr = refset_env$supported_chr
  } else {
    bsg = get_cesa_bsg(cesa)
    supported_chr = cesa@advanced$genome_info$supported_chr
  }
  
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
      } else if(grepl("cannot switch some of.*to NCBI style", conditionMessage(w))) {
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
    cesa@run_history = c(cesa@run_history, deparse(comm, width.cutoff = 500))
  }
  return(cesa)
}
