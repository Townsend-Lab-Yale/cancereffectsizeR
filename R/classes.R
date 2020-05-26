setClass("CESAnalysis", representation(maf = "data.table", trinucleotide_mutation_weights = "list",
          progressions = "character", mutrates_list = "list", dndscv_out_list = "list",
          excluded = "data.table", selection_results = "data.table", gene_epistasis_results = "data.table", coverage = "list",
          genome = "BSgenome", advanced = "list", genome_data_dir = "character", status = "list", samples = "data.table", 
          mutations = "list"))

setMethod("$", "CESAnalysis",
  function(x, name)
  {
    if(name == "maf") {
      if(x@maf[, .N] > 0) {
        return(maf_records(x))
      } 
    } else if (name == "status") {
      return(x) # calls show method
    } else if (name == "samples") {
      return(samples(x))
    } else if (name == "excluded") {
      return(excluded_maf_records(x))
    } else if (name == "trinuc_rates") {
      return(x@trinucleotide_mutation_weights$trinuc_proportion_matrix)
    } else if (name == "mutational_signatures") {
      return(x@trinucleotide_mutation_weights$signature_weight_table)
    }
  }
)

#' @export
.DollarNames.CESAnalysis <- function(x, pattern = "") {
  features = c("status")
  if(x@maf[, .N] > 0) {
    features = c(features, "maf")
  }
  if(x@samples[, .N] > 0) {
    features = c(features, "samples")
  }
  if(x@excluded[, .N] > 0) {
    features = c(features, "excluded")
  }
  if(length(x@trinucleotide_mutation_weights) > 0) {
    features = c(features, "trinuc_rates")
    if ("signature_weight_table" %in% names(x@trinucleotide_mutation_weights)) {
      features = c(features, "mutational_signatures")
    }
  }
  grep(pattern, features, value=TRUE)
}


setMethod("show", "CESAnalysis", 
  function(object) {
    steps = names(object@status)
    for (step in steps) {
      cat(paste0(step,": ", object@status[[step]], "\n"))
    }
    cat(paste0("[Created in cancereffectsizeR, version ", object@advanced$version, ".]"))
  }
)



setValidity("CESAnalysis",
    function(object) {
     if(object@maf[, .N] > 0) {
       if(! identical(colnames(object@maf)[1:5], c("Unique_Patient_Identifier", "Chromosome", "Start_Position",
                                                   "Reference_Allele", "Tumor_Allele") )) {
         return("First MAF columns not Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele")
       }
     } 
     TRUE
    }
)
