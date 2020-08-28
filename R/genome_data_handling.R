#' get_genome_dirs
#' 
#' returns a character vector mapping genome names to their data directories
#' @param full_paths default TRUE; returns full paths to genome data instead of just genome names 
#' @keywords internal
get_genome_dirs = function() {
  genome_superdir = system.file("genomes", package = "cancereffectsizeR")
  genomes = list.dirs(genome_superdir, recursive = F)
  names(genomes) = basename(genomes)
  return(genomes)
}

#' get_genome_data
#' 
#' reads in the requested reference data for the genome build associated with the CESAnalysis
#' @keywords internal
get_genome_data = function(data_dir_or_cesa, datatype) {
  data_dir = data_dir_or_cesa
  if (is(data_dir_or_cesa, "CESAnalysis")) {
    data_dir = data_dir_or_cesa@genome_data_dir
  }
  path = paste0(data_dir, "/", datatype, ".rds")
  if (check_for_genome_data(data_dir_or_cesa, datatype)) {
    return(readRDS(path))
  } else {
    stop(paste0("Something's wrong with the genome data installation:\n",
                "Expected to find data at ", path, "."))
  }
}

#' check_genome_data
#' 
#' checks if the requested reference data exists and returns T/F
#' @keywords internal
check_for_genome_data = function(data_dir_or_cesa, datatype) {
  data_dir = data_dir_or_cesa
  if (is(data_dir_or_cesa, "CESAnalysis")) {
    data_dir = data_dir_or_cesa@genome_data_dir
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
preload_ref_data = function(ref_key) {
  if(is.null(ref_key) || ! is(ref_key, "character") || length(ref_key) != 1) {
    stop("Expected genomic data source to be given as character. Run list_ces_genomes() to see available genomes.", call. = F)
  }
  genome_dirs = get_genome_dirs()
  if(ref_key %in% names(genome_dirs)) {
    data_dir = genome_dirs[ref_key]
  } else {
    stop("Unrecognized genomic data source. Run list_ces_genomes() to see available genomes.", call. = F)
  }
  
  genome_path = paste0(data_dir, "/genome_package_name.rds")
  if(! file.exists(genome_path)) {
    stop(paste0("Something is wrong with the genome data installation.\n",
                "Expected to find a reference genome at ", genome_path, "."))
  }
  
  if (! ref_key %in% ls(.ces_ref_data)) {
    genome_package = readRDS(genome_path)
    bsg = BSgenome::getBSgenome(genome_package)
    GenomeInfoDb::seqlevelsStyle(bsg) = "NCBI"
    
    message("Loading reference data for ", ref_key, "..." )
    .ces_ref_data[[ref_key]] = new.env()
    .ces_ref_data[[ref_key]][["RefCDS"]] = get_genome_data(data_dir, "RefCDS")
    .ces_ref_data[[ref_key]][["gr_genes"]] = get_genome_data(data_dir, "gr_genes")
    .ces_ref_data[[ref_key]][["genome"]] = bsg
    if(check_for_genome_data(data_dir, "generic_exome_gr")) {
      .ces_ref_data[[ref_key]][["default_exome"]] = get_genome_data(data_dir, "generic_exome_gr")
    }
  }
}


#' list_ces_genomes
#' 
#' Prints names of reference data collections that are ready for use with cancereffectsizeR
#' @export
list_ces_genomes = function() {
  genome_names = names(get_genome_dirs())
  if(length(genome_names) == 0) {
    cat("No genome data found. See documentation for how to get genome data before running cancereffectsizeR.\n")
  } else {
    cat(paste0("Available genomes: ", paste(genome_names, collapse = ", "), ".\n"))
  }
}

#' list_ces_covariates
#'
#' Prints names of available built-in covariate data sets for all loaded CES genomes
#' @export
list_ces_covariates = function() {
  loaded_genomes = ls(.ces_ref_data)
  if (length(loaded_genomes) == 0) {
    cat("No covariates data available. Create a CESAnalysis first.\n")
    return(invisible())
  }
  for (genome in loaded_genomes) {
    # To-do: will fail for custom genomes
    cov_files = list.files(system.file(paste0("genomes/", genome, "/covariates/"), package = "cancereffectsizeR"))
    cov_files = gsub("\\.rds$", '', cov_files)
    
    if (length(cov_files) == 0) {
      cat(genome, ": (no covariate data available)\n", sep = "")
    } else {
      initial = paste0(genome, ": ")
      msg = strwrap(initial = initial, x = paste(cov_files, collapse = ", "), exdent = nchar(initial))
      writeLines(msg)
    }
  }
}

#' list_ces_signature_sets
#'
#' Prints names of available mutational signature sets. Just to be clear, we're calling
#' them ces_signature_sets because they're ready to use with cancereffectsizeR. We didn't
#' derive any of these signature sets.
#' @export
list_ces_signature_sets = function() {
  loaded_genomes = ls(.ces_ref_data)
  if (length(loaded_genomes) == 0) {
    cat("No signature data available. Create a CESAnalysis first.\n")
    return(invisible())
  }
  for (genome in loaded_genomes) {
    sig_files = list.files(system.file(paste0("genomes/", genome, "/signatures/"), package = "cancereffectsizeR"))
    sig_sets = gsub("_signatures\\.rds$", '', sig_files)
    if (length(sig_sets) == 0) {
      cat(genome, ": (no signature sets available)\n", sep = "")
    } else {
      initial = paste0(genome, ": ")
      msg = strwrap(initial = initial, x = paste(sig_sets, collapse = ", "), exdent = nchar(initial))
      writeLines(msg)
    }
  }
  cat("[Plug signature set names into the signature_set option of trinuc_mutation_rates().]\n")
}


#' get_ces_signature_set
#'
#' For a given CES reference data collection and signature set name, returns
#' cancereffectsizeR's internal data for the signature set in a three-item list: 
#' the signature set name, a data table of signature metadata, and a signature 
#' definition data frame
#' @export
get_ces_signature_set = function(genome, name) {
  genome_dirs = get_genome_dirs()
  if (! is(genome, "character") | length(genome) != 1) {
    stop("genome should be 1-length character")
  }
  if (! genome %in% names(genome_dirs)) {
    stop("Could not find genome (see list_ces_genomes())")
  }
  genome_dir = genome_dirs[genome]

  sig_file = paste0(genome_dir, "/signatures/", name, "_signatures.rds")
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
  bsg_pkg_name = get_genome_data(cesa, "genome_package_name")
  bsg = BSgenome::getBSgenome(bsg_pkg_name)
  GenomeInfoDb::seqlevelsStyle(bsg) = "NCBI"
  return(bsg)
}




