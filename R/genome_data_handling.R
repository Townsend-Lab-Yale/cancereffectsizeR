#' get_refset_dirs
#' 
#' returns a character vector mapping ref set names to their data directories
#' @keywords internal
get_refset_dirs = function() {
  refsets = sapply(names(.official_refsets), function(x) system.file("refset", package = x))
  refsets = refsets[refsets != ""]
  
  if (exists(".ces_ref_data")) {
    custom_set_names = setdiff(ls(.ces_ref_data), names(refsets))
    if (length(custom_set_names) > 0) {
      custom_set_dirs = sapply(custom_set_names, function(x) .ces_ref_data[[x]][["data_dir"]])
      names(custom_set_dirs) = custom_set_names
      refsets = c(refsets, custom_set_dirs)
    }
  }
  
  return(refsets)
}

#' get_ref_data
#' 
#' reads in the requested reference data for the ref set associated with the CESAnalysis
#' @keywords internal
get_ref_data = function(data_dir_or_cesa, datatype) {
  data_dir = data_dir_or_cesa
  if (is(data_dir_or_cesa, "CESAnalysis")) {
    # Usually, the data should already be loaded into the package environment
    if(data_dir_or_cesa@ref_key %in% ls(.ces_ref_data)) {
      cached_data = .ces_ref_data[[data_dir_or_cesa@ref_key]][[datatype]]
      if(! is.null(cached_data)) {
        return(cached_data)
      }
    }
    data_dir = data_dir_or_cesa@ref_data_dir
  }
  path = paste0(data_dir, "/", datatype, ".rds")
  
  if (check_for_ref_data(data_dir_or_cesa, datatype)) {
    return(readRDS(path))
  } else {
    stop(paste0("Something's wrong with the reference data set:\n",
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
  ref = new.env()
  ref[["RefCDS"]] = get_ref_data(data_dir, "RefCDS")
  
  if(check_for_ref_data(data_dir, "RefCDS.dndscv")) {
    ref[["RefCDS.dndscv"]] = get_ref_data(data_dir, "RefCDS.dndscv")
    ref[["gr_genes.dndscv"]] = get_ref_data(data_dir, "gr_genes.dndscv")
  }
  
  if(check_for_ref_data(data_dir, "transcript_info")) {
    ref[["transcripts"]] = get_ref_data(data_dir, "transcript_info")
  }
  
  ref[["gr_genes"]] = get_ref_data(data_dir, "gr_genes")
  ref[["gene_names"]] = get_ref_data(data_dir, "gene_names")
  
  genome_info = get_ref_data(data_dir, "genome_build_info")
  bsg = BSgenome::getBSgenome(genome_info$BSgenome)
  tryCatch(
    # new versions complain about not all seqlevels being switchable
    suppressWarnings({GenomeInfoDb::seqlevelsStyle(bsg) = seqlevelsStyle(ref[["gr_genes"]])}),
    error = function(e) {
      if(conditionMessage(e) %like% "cannot open URL.*chromInfo" && 
         check_for_ref_data(data_dir, 'cached_chromInfo')) {
        cached_chromInfo = get_ref_data(data_dir, "cached_chromInfo")
        ucsc_info = cached_chromInfo$UCSC
        ncbi_info = cached_chromInfo$NCBI
        assign(ucsc_info$name, ucsc_info$value, envir = get(".UCSC_cached_chrom_info", envir = asNamespace('GenomeInfoDb')))
        assign(ncbi_info$name, ncbi_info$value, envir = get(".NCBI_cached_chrom_info", envir = asNamespace('GenomeInfoDb')))
        suppressWarnings({GenomeInfoDb::seqlevelsStyle(bsg) = seqlevelsStyle(ref[["gr_genes"]])})
      } else {
        stop(e)
      }
    }
  )
   
  ref[["genome"]] = bsg
  ref[["supported_chr"]] = genome_info$supported_chr

  if(check_for_ref_data(data_dir, "default_exome_gr")) {
    ref[["default_exome"]] = get_ref_data(data_dir, "default_exome_gr")
  }
  

  
  # Load trinuc composition of each gene (composition is a 96-length numeric, deconstructSigs order)
  ref[["gene_trinuc_comp"]] = get_ref_data(data_dir, "gene_trinuc_comp")
  
  cov_files = list.files(paste0(data_dir, "/covariates"), full.names = T)
  if (length(cov_files) > 0) {
    cov_list = list()
    for (cov_file in cov_files) {
      cov_name = gsub('.rds$', '', base::basename(cov_file))
      cov_list[[cov_name]] = readRDS(cov_file)
    }
    ref[["covariates"]] = cov_list
  }
  
  signature_sets = list.files(paste0(data_dir, "/signatures/"), pattern =  "_signatures.rds$", full.names = T)
  if (length(signature_sets) > 0) {
    sig_set_list = list()
    for (sig_file in signature_sets) {
      sig_set_name = gsub('_signatures.rds$', '', base::basename(sig_file))
      sig_set_list[[sig_set_name]] = readRDS(sig_file)
    }
    ref[["signatures"]] = sig_set_list
  }
  ref[['data_dir']] = data_dir
  
  preload_anno_files = list.files(paste0(data_dir, "/maf_preload_anno"), pattern = '\\.rds$', full.names = T)
  if(length(preload_anno_files) > 0) {
    ref[['preload_anno']] = lapply(preload_anno_files, readRDS)
  }
  lockEnvironment(ref, bindings = TRUE)
  return(ref)
}


#' list_ces_refsets
#' 
#' Prints names of built-in reference data sets
#' @export
list_ces_refsets = function() {
  genome_names = names(get_refset_dirs())
  if(length(genome_names) == 0) {
    cat("No reference data found. See documentation for how to get reference data before running cancereffectsizeR.\n")
  } else {
    cat(paste0("Available reference data sets: ", paste(genome_names, collapse = ", "), ".\n"))
  }
}

#' list_ces_covariates
#'
#' Prints names of available built-in covariate data sets for all loaded CES genomes
#' @export
list_ces_covariates = function() {
  refset_dirs = get_refset_dirs()
  refset_names = names(refset_dirs)
  if (length(refset_dirs) == 0) {
    pretty_message("No convariates data available.\n")
    return(invisible())
  }
  
  longest_length = max(sapply(refset_names, nchar))
  for (i in 1:length(refset_dirs)) {
    ref_dir = refset_dirs[i]
    refset = refset_names[i]
    cov_files = list.files(paste0(ref_dir, "/covariates/"))
    cov_files = gsub("\\.rds$", '', cov_files)
    
    print_pad = paste0(rep(" ", longest_length - nchar(refset)), collapse = "")
    initial = paste0(print_pad, refset, ": ")
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
  refset_dirs = get_refset_dirs()
  refset_names = names(refset_dirs)
  if (length(refset_dirs) == 0) {
    pretty_message("No signature data available.\n")
    return(invisible())
  }
  for (i in 1:length(refset_dirs)) {
    genome = refset_names[i]
    ref_dir = refset_dirs[i]
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
#' @param refset name of refset (if using a custom refset, it must be loaded into a CESAnalysis already)
#' @param name name of signature set
#' @export
get_ces_signature_set = function(refset, name) {
  refset_dirs = get_refset_dirs()
  if (! is(refset, "character") | length(refset) != 1) {
    stop("refset should be 1-length character")
  }
  if (! refset %in% names(refset_dirs)) {
    stop("Could not find input reference data set (see list_ces_refsets())")
  }
  refset_dir = refset_dirs[refset]

  sig_file = paste0(refset_dir, "/signatures/", name, "_signatures.rds")
  if (! file.exists(sig_file)) {
    no_spaces = paste0(refset_dir, "/signatures/", gsub(' ', '_', name), "_signatures.rds")
    if (! file.exists(no_spaces)) {
      stop("Couldn't find signature data; expected to find it at ", sig_file)
    }
    sig_file = no_spaces
  }
  return(readRDS(sig_file))
}

#' get_cesa_bsg 
#' 
#' Loads the right BSgenome for a CESAnalysis
#' @param cesa CESAnalysis
#' @keywords internal
get_cesa_bsg = function(cesa) {
  ref_key = cesa@ref_key
  if (ref_key %in% ls(.ces_ref_data)) {
    return(.ces_ref_data[[ref_key]]$genome)
  }
  bsg = BSgenome::getBSgenome(cesa@advanced$genome_info$BSgenome)
  
  # Chromosome naming style now present in genome_info. On older versions, it's not,
  # but we can pull the style from gr_genes.
  chromosome_style = cesa@advanced$genome_info$chromosome_style
  if (is.null(chromosome_style)) {
    chromosome_style = seqlevelsStyle(get_ref_data(cesa, 'gr_genes'))[1]
  }
  suppressWarnings({GenomeInfoDb::seqlevelsStyle(bsg) = chromosome_style}) # avoid unswitchable seqlevels complaints
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
    stop("For compatibility with signature extractors, signature definitions must be given as a pure data.frame (see docs).", call. = F)
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
