#' CESAnalysis
#' @description Creates a CESAnalysis object (the data structure used by cancereffectsizeR), and optionally loads in MAF data
#' @param progression_order evolutionary order of tumor stage progression (e.g. c("Primary", "Metastatic"))
#' @param genome_build Genome build of MAF data (currently, just hg19 supported)
#' @param maf optional path of tab-delimited text file in MAF format, or an MAF in data.frame format
#' @param ... when loading MAF data, supply additional arguments for load_maf
#' @return CESAnalysis object
#' @export

CESAnalysis = function(progression_order = NULL, maf = NULL, genome_build = "hg19", ...) {
  supported_genomes = c("hg19")
  
  if (length(genome_build) != 1) {
    stop("Genome build should be a 1-length vector.")
  }
  if (! genome_build %in% supported_genomes) {
    stop("Sorry, your chosen genome build is not supported. Let us know you want it!")
  }
  
  if (is.null(progression_order)) {
    progression_order = c("1")
  }
  
  cesa = new("CESAnalysis", genome_build = genome_build, progressions = CESProgressions(order = progression_order))
  
  if(! is.null(maf)) {
    cesa = load_maf(cesa = cesa, maf = maf, ...)
  }
  return(cesa)
  
}

