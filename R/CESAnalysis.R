#' CESAnalysis
#' @description Creates a CESAnalysis object, the central data structure of cancereffectsizeR
#' @param progression_order evolutionary order of tumor stage progression (e.g. c("Primary", "Metastatic"))
#' @param genome Genome build of MAF data (currently, just hg19 supported)
#' @return CESAnalysis object
#' @export
CESAnalysis = function(genome = NULL, progression_order = NULL) {
  genome_dir = get_genome_data_directory(genome)
  genome_path = paste0(genome_dir, "/genome_package_name.rds")
  if(! file.exists(genome_path)) {
    stop(paste0("Something is wrong with the genome data installation.\n",
                "Expected to find a reference genome at ", genome_path, "."))
  }
  genome_package = readRDS(genome_path)
  genome = BSgenome::getBSgenome(genome_package)
  GenomeInfoDb::seqlevelsStyle(genome) = "NCBI"
  message(crayon::black(paste0("Okay, this CES analysis will use the ", 
                               tolower(BSgenome::commonName(genome)), " genome (", BSgenome::releaseName(genome), ").")))
  message(crayon::black("Note: We'll be using NCBI-style chromosome names (i.e., no \"chr\" prefixes)."))
  
  if (is.null(progression_order)) {
    progression_order = c("1")
  }

  status = list("genome" = GenomeInfoDb::providerVersion(genome),
                "progressions" = paste0(as.character(progression_order), collapse = ", "),
                "MAF data" = "none so far (run load_maf)",
                "trinucleotide mutation rates" = "uncalculated (run trinucleotide_mutation_weights)",
                "gene mutation rates" = "uncalculated (need trinucleotide rates first)"
                )
  if(length(progression_order) == 1) {
    status[["progressions"]] = NULL
  }
  advanced = list("version" = packageVersion("cancereffectsizeR"))
  cesa = new("CESAnalysis", status = status, genome = genome, maf = data.table(), excluded = data.table(),
             progressions = CESProgressions(order = progression_order), 
             gene_epistasis_results = data.table(), selection_results = data.table(), genome_data_dir = genome_dir,
             advanced = advanced)
  return(cesa)
}

#' maf
#' 
#' returns a data.table containing MAF records used in the given CESAnalysis
#' @param cesa CESAnalysis object
#' @export
maf = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: maf(cesa), where cesa is a CESAnalysis")
  }
  if(cesa@maf[,.N] == 0) {
    stop("No MAF data has been loaded")
  }
  return(cesa@maf)
}

#' excluded_maf_records
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

