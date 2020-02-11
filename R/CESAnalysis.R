#' CESAnalysis
#' @description Creates a CESAnalysis object (the data structure used by cancereffectsizeR)
#' @param progression_order evolutionary order of tumor stage progression (e.g. c("Primary", "Metastatic"))
#' @param genome Genome build of MAF data (currently, just hg19 supported)
#' @return CESAnalysis object
#' @export

CESAnalysis = function(genome = NULL, progression_order = NULL) {
  if(is.null(genome) || ! is(genome, "character") || length(genome) != 1) {
    stop("You must specify a genome build (e.g., \"hg38\"). Run list_genomes() to see available genomes.")
  }
  genome_dirs = get_genome_dirs()
  if(genome %in% names(genome_dirs)) {
    genome_dir = genome_dirs[genome]
  } else {
    stop("Unrecognized genome. Run list_genomes() to see available genomes.")
  }
  
  genome_path = paste0(genome_dir, "/genome.rds")
  if(! file.exists(genome_path)) {
    stop(paste0("Something is wrong with the genome data installation.\n",
                "Expected to find a reference genome at ", genome_path, "."))
  }
  genome = readRDS(genome_path)
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
  cesa = new("CESAnalysis", status = status, genome = genome, progressions = CESProgressions(order = progression_order), 
             gene_epistasis_results = data.table(), selection_results = data.table(), genome_data_dir = genome_dir)
  return(cesa)
}

