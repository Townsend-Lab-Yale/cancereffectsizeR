#' get_genome_dirs
#' returns a character vector with the paths of genome data available to use with CES (or just genome names)
#' does not verify that all required files are present in the genome directories
#' @param full_paths default TRUE; returns full paths to genome data instead of just genome names 
get_genome_dirs = function(full_paths = TRUE) {
  genome_superdir = system.file("genomes", package = "cancereffectsizeR")
  genomes = list.dirs(genome_superdir,recursive = F, full.names = full_paths)
  return(genomes)
}


#' prints names of  all genome builds that are ready for use with cancereffectsizeR
#' list_genomes
#' @export
list_genomes = function() {
  genome_names = get_genome_dirs(full_paths = FALSE)
  if(length(genome_names) == 0) {
    message(silver("No genome data found. See documentation for how to get genome data before running cancereffectsizeR."))
  } else {
    message(silver(paste0("Available genomes: ", paste(genome_names, collapse = ", "), ".")))
  }
}

