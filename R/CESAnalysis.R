#' CESAnalysis
#' @description Creates a CESAnalysis object (the data structure used by cancereffectsizeR)
#' @param progression_order evolutionary order of tumor stage progression (e.g. c("Primary", "Metastatic"))
#' @param genome Genome build of MAF data (currently, just hg19 supported)
#' @return CESAnalysis object
#' @export

CESAnalysis = function(genome = NULL, progression_order = NULL) {
  if(is.null(genome)) {
    stop("You must specify a genome build: either \"hg19\", \"hg38\", or any BSgenome object.")
  }
  if("BSgenome" %in% class(genome)) {
    message(crayon::black(paste0("Okay, this CES analysis will use the ", 
                                 tolower(BSgenome::commonName(genome)), "genome (", BSgenome::releaseName(genome), ").")))
  } else if ("character" %in% class(genome) && length(genome) == 1) {
    genome = tolower(genome)
    if(genome %in% c("hg19", "grch37")) {
      if(! requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = T)) {
        stop("You need to install BSgenome.Hsapiens.UCSC.hg19 from Bioconductor.")
      }
      message(crayon::black("This CES analysis will use the hg19/GRCh37 build of the human genome."))
      genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    } else if (genome %in% c("hg38", "grch38")) {
      if(! requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = T)) {
        stop("You need to install BSgenome.Hsapiens.UCSC.hg38 from Bioconductor.")
      }
      message(crayon::black("This CES analysis will use the hg38/GRCh38 build of the human genome."))
      genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    } else {
      stop("Unrecognized genome. Choose \"hg19\" or \"hg38\", or supply any BSgenome object.")
    }
  }
  GenomeInfoDb::seqlevelsStyle(genome) = "NCBI"
  message(crayon::black("Note: We'll be using NCBI-style chromosome names (i.e., no \"chr\" prefixes)."))
  
  if (is.null(progression_order)) {
    progression_order = c("1")
  }
  
  cesa = new("CESAnalysis", genome = genome, progressions = CESProgressions(order = progression_order), 
             gene_epistasis_results = data.table(), selection_results = data.table())
  return(cesa)
}

