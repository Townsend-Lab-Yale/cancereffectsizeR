#' get_ref_set_dirs
#' 
#' returns a character vector mapping ref set names to their data directories
#' @keywords internal
get_ref_set_dirs = function() {
  ref_set_super_dir = system.file("ref_sets", package = "cancereffectsizeR")
  ref_sets = list.dirs(ref_set_super_dir, recursive = F)
  names(ref_sets) = basename(ref_sets)
  
  if (exists(".ces_ref_data")) {
    custom_set_names = setdiff(ls(.ces_ref_data), names(ref_sets))
    if (length(custom_set_names) > 0) {
      custom_set_dirs = sapply(custom_set_names, function(x) .ces_ref_data[[x]][["data_dir"]])
      names(custom_set_dirs) = custom_set_names
      ref_sets = c(ref_sets, custom_set_dirs)
    }
  }
  
  return(ref_sets)
}

#' get_ref_data
#' 
#' reads in the requested reference data for the ref set associated with the CESAnalysis
#' @keywords internal
get_ref_data = function(data_dir_or_cesa, datatype) {
  data_dir = data_dir_or_cesa
  if (is(data_dir_or_cesa, "CESAnalysis")) {
    data_dir = data_dir_or_cesa@ref_data_dir
  }
  path = paste0(data_dir, "/", datatype, ".rds")
  if (check_for_ref_data(data_dir_or_cesa, datatype)) {
    return(readRDS(path))
  } else {
    stop(paste0("Something's wrong with the reference data installation:\n",
                "Expected to find data at ", path, "."))
  }
}

#' check_for_ref_data
#' 
#' checks if the requested reference data exists and returns T/F
#' @param data_dir_or_cesa CESAnalysis, or file path for reference data directory
#' @keywords internal
check_for_ref_data = function(data_dir_or_cesa, datatype) {
  data_dir = data_dir_or_cesa
  if (is(data_dir_or_cesa, "CESAnalysis")) {
    data_dir = data_dir_or_cesa@ref_data_dir
  }
  if (! dir.exists(data_dir)) {
    error_message = paste0("Error: Could not locate genome data directory. This can happen when a CESAnalysis is saved and\n",
                           "reloaded in a different computing environment that does not have the data available.")
    stop(error_message)
  }
  path = paste0(data_dir, "/", datatype, ".rds")
  return(file.exists(path))
}




#' preload_ref_data
#' 
#' Used when loading or creating a CESAnalysis to load reference into an environment for quick access
#' @keywords internal
preload_ref_data = function(data_dir) {
  
  # by design, name of ref set is always the directory name
  ref_set_name = basename(data_dir) 
  
  if (! ref_set_name %in% ls(.ces_ref_data)) {
    message("Loading reference data set ", ref_set_name, "..." )
    .ces_ref_data[[ref_set_name]] = new.env()
    .ces_ref_data[[ref_set_name]][["RefCDS"]] = get_ref_data(data_dir, "RefCDS")
    .ces_ref_data[[ref_set_name]][["gr_genes"]] = get_ref_data(data_dir, "gr_genes")
    .ces_ref_data[[ref_set_name]][["gene_names"]] = get_ref_data(data_dir, "gene_names")
    
    
    genome_info = get_ref_data(data_dir, "genome_build_info")
    bsg = BSgenome::getBSgenome(genome_info$BSgenome)
    GenomeInfoDb::seqlevelsStyle(bsg) = "NCBI"
    .ces_ref_data[[ref_set_name]][["genome"]] = bsg

    if(check_for_ref_data(data_dir, "generic_exome_gr")) {
      .ces_ref_data[[ref_set_name]][["default_exome"]] = get_ref_data(data_dir, "generic_exome_gr")
    }
    # Load trinuc composition of each gene (composition is a 96-length numeric, deconstructSigs order)
    .ces_ref_data[[ref_set_name]][["gene_trinuc_comp"]] = get_ref_data(data_dir, "gene_trinuc_comp")
    .ces_ref_data[[ref_set_name]][["data_dir"]] = data_dir
  }
}


#' list_ces_ref_sets
#' 
#' Prints names of built-in reference data sets
#' @export
list_ces_ref_sets = function() {
  genome_names = names(get_ref_set_dirs())
  if(length(genome_names) == 0) {
    cat("No genome data found. See documentation for how to get genome data before running cancereffectsizeR.\n")
  } else {
    cat(paste0("Available reference data sets: ", paste(genome_names, collapse = ", "), ".\n"))
  }
}

#' list_ces_covariates
#'
#' Prints names of available built-in covariate data sets for all loaded CES genomes
#' @export
list_ces_covariates = function() {
  ref_set_dirs = get_ref_set_dirs()
  ref_set_names = names(ref_set_dirs)
  if (length(ref_set_dirs) == 0) {
    pretty_message("No convariates data available.\n")
    return(invisible())
  }
  
  longest_length = max(sapply(ref_set_names, nchar))
  for (i in 1:length(ref_set_dirs)) {
    ref_dir = ref_set_dirs[i]
    ref_set = ref_set_names[i]
    cov_files = list.files(paste0(ref_dir, "/covariates/"))
    cov_files = gsub("\\.rds$", '', cov_files)
    
    print_pad = paste0(rep(" ", longest_length - nchar(ref_set)), collapse = "")
    initial = paste0(print_pad, ref_set, ": ")
    exdent = nchar(initial)
    covs = ifelse(length(cov_files) > 0, paste(cov_files, collapse = ", "), "(no covariate data available)")
    msg = paste0(strwrap(initial = initial, x = covs, exdent = exdent), collapse = "\n")
    message(crayon::black(msg))
  }
}

#' list_ces_signature_sets
#'
#' Prints names of available mutational signature sets. Just to be clear, we're calling
#' them ces_signature_sets because they're ready to use with cancereffectsizeR. We didn't
#' derive any of these signature sets.
#' @export
list_ces_signature_sets = function() {
  ref_set_dirs = get_ref_set_dirs()
  ref_set_names = names(ref_set_dirs)
  if (length(ref_set_dirs) == 0) {
    pretty_message("No signature data available.\n")
    return(invisible())
  }
  for (i in 1:length(ref_set_dirs)) {
    genome = ref_set_names[i]
    ref_dir = ref_set_dirs[i]
    sig_files = list.files(paste0(ref_dir, "/signatures/"))
    sig_sets = gsub("_signatures\\.rds$", '', sig_files)
    if (length(sig_sets) == 0) {
      pretty_message(paste0(genome, ": (no signature sets available)"))
    } else {
      initial = paste0(genome, ": ")
      msg = strwrap(initial = initial, x = paste(sig_sets, collapse = ", "), exdent = nchar(initial))
      pretty_message(msg)
    }
  }
  pretty_message("[Plug signature set names into the signature_set option of trinuc_mutation_rates().]\n")
}


#' get_ces_signature_set
#'
#' For a given CES reference data collection and signature set name, returns
#' cancereffectsizeR's internal data for the signature set in a three-item list: 
#' the signature set name, a data table of signature metadata, and a signature 
#' definition data frame
#' @export
get_ces_signature_set = function(ref_set, name) {
  ref_set_dirs = get_ref_set_dirs()
  if (! is(ref_set, "character") | length(ref_set) != 1) {
    stop("ref_set should be 1-length character")
  }
  if (! ref_set %in% names(ref_set_dirs)) {
    stop("Could not find input reference data set (see list_ces_ref_sets())")
  }
  ref_set_dir = ref_set_dirs[ref_set]

  sig_file = paste0(ref_set_dir, "/signatures/", name, "_signatures.rds")
  if (! file.exists(sig_file)) {
    stop("Couldn't find signature data; expected to find it at ", sig_file)
  }
  return(readRDS(sig_file))
}

#' get_cesa_bsg 
#' 
#' Loads the right BSgenome for a CESAnalysis, formatted NCBI-style
#' @param cesa CESAnalysis
#' @keywords internal
get_cesa_bsg = function(cesa) {
  ref_key = cesa@ref_key
  if (ref_key %in% ls(.ces_ref_data)) {
    return(.ces_ref_data[[ref_key]]$genome)
  }
  bsg = BSgenome::getBSgenome(cesa@advanced$genome_info$BSgenome)
  GenomeInfoDb::seqlevelsStyle(bsg) = "NCBI"
  return(bsg)
}

#' validate_signature_set
#' 
#' Checks if a custom CES signature is properly formatted; stops with an error if not
#' @param signature_set signature set list (see docs for format)
#' @export
validate_signature_set = function(signature_set) {
  signature_set_data = signature_set
  
  if (! is(signature_set_data, "list")) {
    stop("Invalid signature set: type should be list")
  }
  signature_set_name = signature_set_data$name
  signatures = signature_set_data$signatures
  signature_metadata = signature_set_data$meta
  if (! is(signature_set_name, "character") || length(signature_set_name) != 1 || ! is(signatures, "data.frame") ||
      ! is(signature_metadata, "data.table")) {
    stop("Improperly formatted custom signature set: name should be 1-length character, signatures should be data.frame, meta should be data.table", call. = F)
  }
  if(is(signatures, "data.table") || is(signatures, "tbl")) {
    stop("For compatibility with deconstructSigs, signature definitions must be given as a pure data.frame (see docs).", call. = F)
  }
  if(! identical(sort(deconstructSigs_trinuc_string), sort(colnames(signatures)))) {
    tmp = paste0("\n", 'c("', paste(deconstructSigs_trinuc_string, collapse = '", "'), '")')
    message("Expected signature definition column names:")
    message(strwrap(tmp, indent = 4, exdent = 4))
    stop("Your signature definition data frame has improper column names.", call. = F)
  }
  # Validate signature metadata if it's not empty
  if (signature_metadata[, .N] > 0) {
    if (is.null(signature_metadata$Signature)) {
      stop("Signature metadata incorrectly formatted (see docs).")
    }
    if (any(! rownames(signatures) %in% signature_metadata$Signature)) {
      stop("Improperly formatted signature set: Some signatures in your signature definitions are missing from the metadata table.")
    }
    if(length(signature_metadata$Signature) != length(unique(signature_metadata$Signature))) {
      stop("Improperly formatted signature set: Some signatures are repeated in your signature metadata table")
    }
    
    meta_cols = colnames(signature_metadata)
    if ("Likely_Artifact" %in% meta_cols && ! is(signature_metadata$Likely_Artifact, "logical")) {
      stop("Improperly formatted signature set metadata: column Likely_Artifact should be logical.")
    }
    
    num_hm_cols = sum(c("Exome_Min", "Genome_Min") %in% meta_cols)
    if (num_hm_cols == 2) {
      if (! signature_metadata[, all(sapply(.SD, is.numeric)), .SDcols = c("Exome_Min", "Genome_Min")]) {
        stop("Invalid signature set: metadata columns Exome_Min and Genome_Min should be numeric.")
      }
      if (signature_metadata[, any(xor(is.na(Exome_Min), is.na(Genome_Min)))]) {
        stop("Invalid signature set metadata: For any signature, Exome_Min/Genome_Min must be both defined or both NA.")
      }
      if (signature_metadata[! is.na(Exome_Min), any(Exome_Min >= Genome_Min)]) {
        stop("Invalid signature set metadata: Exome_Min must be <= Genome_Min in all signatures where they're defined.")
      }
    } else if (num_hm_cols != 0) {
      stop("Invalid signature set metadata: there should be both Exome_Min and Genome_Min numeric columns, or neither.")
    }
  }
}
