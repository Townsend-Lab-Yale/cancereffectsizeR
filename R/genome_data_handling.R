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
get_genome_data = function(cesa, datatype) {
  path = paste0(cesa@genome_data_dir, "/", datatype, ".rds")
  if (check_for_genome_data(cesa, datatype)) {
    return(readRDS(path))
  } else {
    stop(paste0("Something's wrong with the genome data installation:\n",
                "Expected to find data at ", path, "."))
  }
}

#' check_genome_data
#' 
#' checks if the requested reference data exists and returns T/F
check_for_genome_data = function(cesa, datatype) {
  data_dir = cesa@genome_data_dir
  if (! dir.exists(data_dir)) {
    error_message = paste0("Error: Could not locate genome data directory. This can happen when a CESAnalysis is saved and\n",
                           "reloaded in a different computing environment. Assuming you have the necessary data installed,\n",
                           "you can try running set_genome_data_directory() to fix the problem.") 
    stop(error_message)
  }
  path = paste0(data_dir, "/", datatype, ".rds")
  return(file.exists(path))
}

#' get_genome_data_directory
#' 
#' takes in a genome name (e.g., "hg38"), searches for an associated data directory, and returns the path
get_genome_data_directory = function(genome) {
  if(is.null(genome) || ! is(genome, "character") || length(genome) != 1) {
    stop("You must specify a genome build (e.g., \"hg38\"). Run list_genomes() to see available genomes.")
  }
  genome_dirs = get_genome_dirs()
  if(genome %in% names(genome_dirs)) {
    genome_dir = genome_dirs[genome]
  } else {
    stop("Unrecognized genome. Run list_genomes() to see available genomes.")
  }
  return(genome_dir)
}

#' set_genome_data_directory
#' 
#' assigns a genome data directory to a CESAnalysis; can be used when loading a saved CESAnalysis in a new envrionment
#' @export
set_genome_data_directory = function(cesa, genome) {
  genome_dir = get_genome_data_directory(genome)
  cesa@genome_data_dir = genome_dir
  return(cesa)
}


#' list_genomes
#' 
#' prints names of all genome builds that are ready for use with cancereffectsizeR
#' @export
list_genomes = function() {
  genome_names = names(get_genome_dirs())
  if(length(genome_names) == 0) {
    message(silver("No genome data found. See documentation for how to get genome data before running cancereffectsizeR."))
  } else {
    message(silver(paste0("Available genomes: ", paste(genome_names, collapse = ", "), ".")))
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


