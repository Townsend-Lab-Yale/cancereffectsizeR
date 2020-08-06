#' get_genome_dirs
#' 
#' returns a character vector mapping genome names to their data directories
#' @param full_paths default TRUE; returns full paths to genome data instead of just genome names 
get_genome_dirs = function() {
  genome_superdir = system.file("genomes", package = "cancereffectsizeR")
  genomes = list.dirs(genome_superdir, recursive = F)
  names(genomes) = basename(genomes)
  return(genomes)
}

#' get_genome_data
#' 
#' reads in the requested reference data for the genome build associated with the CESAnalysis
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
check_for_genome_data = function(data_dir_or_cesa, datatype) {
  data_dir = data_dir_or_cesa
  if (is(data_dir_or_cesa, "CESAnalysis")) {
    data_dir = cesa@genome_data_dir
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
preload_ref_data = function(ref_key) {
  if(is.null(ref_key) || ! is(ref_key, "character") || length(ref_key) != 1) {
    stop("Expected genomic data source to be given as character. Run list_genomes() to see available genomes.", call. = F)
  }
  genome_dirs = get_genome_dirs()
  if(ref_key %in% names(genome_dirs)) {
    data_dir = genome_dirs[ref_key]
  } else {
    stop("Unrecognized genomic data source. Run list_genomes() to see available genomes.", call. = F)
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


#' list_genomes
#' 
#' Prints names of all reference data collections that are ready for use with cancereffectsizeR
#' @export
list_genomes = function() {
  genome_names = names(get_genome_dirs())
  if(length(genome_names) == 0) {
    cat("No genome data found. See documentation for how to get genome data before running cancereffectsizeR.\n")
  } else {
    cat(paste0("Available genomes: ", paste(genome_names, collapse = ", "), ".\n"))
  }
}
#' 
#' #' list_gene_mutation_covariates (coming soon)
#' #' 
#' #' for a given CESAnalysis (that is, the genome build it uses), prints available covariates for use with gene mutation rate calculation
#' #' @export
#' list_gene_mutation_covariates = function() {
#'   
#' }


