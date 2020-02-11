#' get_genome_dirs
#' returns a character vector mapping genome names to their data directories
#' @param full_paths default TRUE; returns full paths to genome data instead of just genome names 
get_genome_dirs = function() {
  genome_superdir = system.file("genomes", package = "cancereffectsizeR")
  genomes = list.dirs(genome_superdir, recursive = F)
  names(genomes) = basename(genomes)
  return(genomes)
}

#' get_genome_data
#' reads in the requested reference data for the genome build associated with the CESAnalysis
get_genome_data = function(cesa, datatype) {
  data_dir = cesa@genome_data_dir
  path = paste0(data_dir, "/", datatype, ".rds")
  if (! file.exists(path)) {
    stop(paste0("Something's wrong with the genome data installation\n",
                 "Expected to find data at ", path, "."))
  }
  return(readRDS(path))
}

#' prints names of  all genome builds that are ready for use with cancereffectsizeR
#' list_genomes
#' @export
list_genomes = function() {
  genome_names = names(get_genome_dirs())
  if(length(genome_names) == 0) {
    message(silver("No genome data found. See documentation for how to get genome data before running cancereffectsizeR."))
  } else {
    message(silver(paste0("Available genomes: ", paste(genome_names, collapse = ", "), ".")))
  }
}


