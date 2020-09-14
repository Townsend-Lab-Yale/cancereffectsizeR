setClass("CESAnalysis", representation(maf = "data.table", trinucleotide_mutation_weights = "list",
          groups = "character", mutrates = "data.table", dndscv_out_list = "list",
          excluded = "data.table", selection_results = "data.table", coverage = "list",
          ref_key = "character", advanced = "list", ref_data_dir = "character", run_history = "character", samples = "data.table", 
          mutations = "list"))


setMethod("$", "CESAnalysis",
  function(x, name)
  {
    if(name == "maf") {
      if(x@maf[, .N] > 0) {
        return(maf_records(x))
      } 
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
    } else if (name == "variants") {
      return(select_variants(cesa = x, min_freq = 0))
    } else if (name == "selection") {
      return(snv_results(x))
    } else if (name == "reference_data") {
      if (! x@ref_key %in% ls(.ces_ref_data)) {
        preload_ref_data(x@ref_key)
      }
      ref_data = list(RefCDS = .ces_ref_data[[x@ref_key]]$RefCDS, gene_ranges = .ces_ref_data[[x@ref_key]]$gr_genes,
                      genome = .ces_ref_data[[x@ref_key]]$genome)
      snv_signatures = x@advanced$snv_signatures
      if (! is.null(snv_signatures)) {
        ref_data = c(ref_data, list(snv_signatures = snv_signatures))
      }
      return(invisible(ref_data))
    } else if (name == "coverage_ranges") {
        return(x@coverage)
    } else if (name == "run_history") {
      ces_version = paste0("[Version: cancereffectsizeR ", as.character(x@advanced$version), ']')
      run_history = c(x@run_history, "", ces_version)
      CES_Run_History(run_history)
    }
  }
)

#' @export
.DollarNames.CESAnalysis <- function(x, pattern = "") {
  features = character()
  if(x@maf[, .N] > 0) {
    features = c(features, "maf")
  }
  if(x@samples[, .N] > 0) {
    features = c(features, "samples")
  }
  if(length(x@mutations) > 0) {
    features = c(features, "variants")
  }
  if(x@selection_results[, .N] > 0) {
    features = c(features, "selection")
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
  features = c(features, c("reference_data", "coverage_ranges", "run_history"))
  grep(pattern, features, value=TRUE)
}


setMethod("show", "CESAnalysis", 
  function(object) {
    if(object@maf[, .N] > 0) {
      cat("Samples:\n")
      print(object@samples[, .(num_samples = .N), by = c("group", "group_index", "coverage")][order(group_index)], row.names = F)
      num_snvs = object@maf[variant_type == "snv", .N]
      cat("\nMAF data: ", num_snvs, " SNVs loaded", sep = "")
      if (identical(object@advanced$annotated, TRUE)) {
        cat(" and annotated.\n")
      } else {
        cat(" (but not annotated).\n")
      }
    }
    signature_set = object@advanced$snv_signatures
    if (! is.null(signature_set)) {
      signature_set_name = signature_set$name
      cat("Mutational processes: Sample-level extraction of ", signature_set_name, " SNV signatures.\n", sep = "")
    }
    cat("CES reference data set: ", object@ref_key, "\n", sep = "")
    cat("Run history: See [CESAnalysis]$run_history.\n")
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

setClass("CES_Run_History", representation(history = "character"))
CES_Run_History = function(history) {
  new("CES_Run_History", history = history)
  
}

setMethod("show", "CES_Run_History",
  function(object) {
    run_history = strwrap(object@history, exdent = 4)
    writeLines(run_history)
  }
)

as.character.CES_Run_History = function(object) {
  history = object@history[object@history != ""]
}

