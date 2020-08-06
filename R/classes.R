setClass("CESAnalysis", representation(maf = "data.table", trinucleotide_mutation_weights = "list",
          progressions = "character", mutrates = "data.table", dndscv_out_list = "list",
          excluded = "data.table", selection_results = "data.table", gene_epistasis_results = "data.table", coverage = "list",
          ref_key = "character", advanced = "list", genome_data_dir = "character", status = "list", samples = "data.table", 
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
      return(get_sample_info(x))
    } else if (name == "excluded") {
      return(excluded_maf_records(x))
    } else if (name == "trinuc_rates") {
      return(get_trinuc_rates(x))
    } else if (name == "mutational_signatures") {
      return(get_signature_weights(x))
    } else if (name == "gene_rates") {
      return(get_gene_rates(x))
    } else if (name == "mutations") {
      return(get_mutations(x))
    } else if (name == "selection") {
      return(snv_results(x))
    } else if (name == "epistasis") {
      return(gene_epistasis_results(x))
    } else if (name == "reference_data") {
      return(list(snv_signatures = x@advanced$snv_signatures))
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
  if(length(x@mutations) > 0) {
    features = c(features, "mutations")
  }
  if(x@selection_results[, .N] > 0) {
    features = c(features, "selection")
  }
  if(x@gene_epistasis_results[, .N] > 0) {
    features = c(features, "epistasis")
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
  if(x@mutrates[, .N] > 0) {
    features = c(features, "gene_rates")
  }
  if("snv_signatures" %in% names(x@advanced)) {
    features = c(features, "reference_data")
  }

  grep(pattern, features, value=TRUE)
}


setMethod("show", "CESAnalysis", 
  function(object) {
    steps = names(object@status)
    if(object@maf[, .N] > 0) {
      cat("Samples:\n")
      print(object@samples[, .(num_samples = .N), by = c("progression_name", "progression_index", "coverage")][order(progression_index)], row.names = F)
    }
    if(length(object@mutations) > 0) {
      cat("\nAnnotated mutations:\n")
      print(object@mutations)
    }
    if(! is.null(object@trinucleotide_mutation_weights$signature_weight_table)) {
      cat("\nSNV signatures:\n")
      print(object@trinucleotide_mutation_weights$signature_weight_table, topn = 5)
    }
    if(object@mutrates[, .N] > 0) {
      cat("\nGene mutation rates:\n")
      print(object@mutrates)
    }
    if(object@selection_results[, .N] > 0) {
      cat("\nSelection intensities of single variants:\n")
      print(object@selection_results, topn = 5)
    }
    if(object@gene_epistasis_results[, .N] > 0) {
      cat("\nGene-level recurrent variant epistasis:\n")
      print(object@gene_epistasis_results)
    }
    cat("\nRun summary:\n")
    for (step in steps) {
      cat(paste0(step,": ", object@status[[step]], "\n"))
    }
    cat(paste0("\n[Created in cancereffectsizeR, version ", object@advanced$version, ".]"))
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
