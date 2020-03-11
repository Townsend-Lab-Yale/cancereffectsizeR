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
  
  
  # Validate progression_order
  if (is.null(progression_order)) {
    progression_order = c("stageless")
  }
  
  if (is.numeric(progression_order)) {
    progression_order = as.character(progression_order)
  }

  if (! is(progression_order, "character")) {
    stop("progression_order should be a character vector of chronological tumor states (e.g., Primary, Metastatic)")
  }

  
  # simple analysis status tracking; used to guide user in show(CESAnalysis)
  status = list("genome" = GenomeInfoDb::providerVersion(genome),
                "progressions" = paste0(progression_order, collapse = ", "),
                "MAF data" = "none so far (run load_maf)",
                "trinucleotide mutation rates" = "uncalculated (run trinucleotide_mutation_weights)",
                "gene mutation rates" = "uncalculated (need trinucleotide rates first)"
                )
  if(length(progression_order) == 1) {
    status[["progressions"]] = NULL
  }
  advanced = list("version" = packageVersion("cancereffectsizeR"))
  cesa = new("CESAnalysis", status = status, genome = genome, maf = data.table(), excluded = data.table(),
             progressions = progression_order,
             gene_epistasis_results = data.table(), selection_results = data.table(), genome_data_dir = genome_dir,
             advanced = advanced, samples = data.table())
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

#' samples
#' 
#' returns a data.table with info on all samples in the CESAnalysis (at least, all samples with any valid mutations)
#' @param cesa CESAnalysis object
#' @export
samples = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: samples(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@samples[,.N] == 0) {
    stop("No MAF data has been loaded yet, so naturally there is no sample data.")
  }
  
  # user doesn't need progression columns for single-progression-state analyses
  if(length(cesa@progressions) == 1) {
    return(cesa@samples[, -c("progression_index", "progression_name")])
  } else {
    return(cesa@samples)
  }
  return(cesa@samples)
}

#' snv_results
#' 
#' returns a data table of SNV effect sizes generated with ces_snv
#' @param cesa CESAnalysis object
#' @export
snv_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: snv_results(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@selection_results[, .N] == 0) {
    stop("No results yet from ces_snv in this CESAnalysis")
  }
  return(cesa@selection_results)
}

#' gene_epistasis_results
#' 
#' returns a data table of pairwise gene epistasis effect sizes generated with ces_gene_epistasis
#' @param cesa CESAnalysis object
#' @export
gene_epistasis_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: gene_epistasis_results(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@gene_epistasis_results[, .N] == 0) {
    stop("No results yet from ces_gene_epistasis in this CESAnalysis")
  }
  return(cesa@gene_epistasis_results)
}

