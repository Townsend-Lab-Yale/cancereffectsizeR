setClass("CESAnalysis", representation(maf = "data.table", trinucleotide_mutation_weights = "list",
          mutrates = "data.table", dndscv_out_list = "list",
          excluded = "data.table", selection_results = "list", coverage = "list", epistasis = "list",
          ref_key = "character", advanced = "list", ref_data_dir = "character", run_history = "character", samples = "data.table", 
          mutations = "list"))

#' @export
.DollarNames.CESAnalysis <- function(x, pattern = "") {
  features = character()
  if(x@maf[, .N] > 0) {
    features = c(features, "maf")
  }
  if(x@samples[, .N] > 0) {
    features = c(features, "samples")
  }
  if(sum(sapply(x@mutations, function(y) y[, .N])) > 0) {
    features = c(features, "variants")
  }
  if(length(x@selection_results) > 0) {
    features = c(features, "selection")
  }
  if(length(x@epistasis) > 0) {
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
  if(length(x@dndscv_out_list) > 0) {
    features = c(features, "dNdScv_results")
  }
  features = c(features, "reference_data")
  if(length(x@coverage) > 0) {
    features = c(features, "coverage_ranges")
  }
  features = c(features, "run_history")
  grep(pattern, features, value=TRUE)
}

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
              return(list(snv_counts = x@trinucleotide_mutation_weights$trinuc_snv_counts,
                          raw_attributions = get_signature_weights(x, raw = T),
                          biological_weights = get_signature_weights(x, raw = F),
                          help = function() cat("See Value in ?trinuc_mutation_rates for methods and tips.\n")))
            } else if (name == "gene_rates") {
              return(get_gene_rates(x))
            } else if (name == "variants") {
              cached_variants = copy(x@advanced$cached_variants)
              if (is.null(cached_variants)) {
                x@advanced$cached_variants = suppressMessages(select_variants(cesa = x))
                cached_variants = copy(x@advanced$cached_variants)
              }
              return(cached_variants)
            } else if (name == "selection") {
              return(snv_results(x))
            } else if (name == "epistasis") {
              return(epistasis_results(x))
            } else if (name == "reference_data") {
              ref_data = list(genome = get_cesa_bsg(x))
              snv_signatures = copy(x@advanced$snv_signatures)
              if (length(snv_signatures) > 0) {
                ref_data = c(ref_data, list(snv_signatures = snv_signatures))
              }
              return(invisible(ref_data))
            } else if (name == "coverage_ranges") {
              return(x@coverage) 
            } else if(name == "dNdScv_results") {
              return(lapply(x@dndscv_out_list, 
                            function(y) {
                              if(is.data.table(y)) {
                                return(copy(y))
                              } else {
                                # for the rare user who requests all dNdScv output
                                return(y)
                              }
                            }))
            }
              else if (name == "run_history") {
              CES_Run_History(x@run_history)
            }
          }
)


setMethod("show", "CESAnalysis", 
  function(object) {
    genome_name = object@advanced$genome_info$build_name
    cat("CESAnalysis of ", genome_name, " data\n", sep = "")
    refset_version = object@advanced$refset_version
    refset_version_msg = ''
    if(! is.null(refset_version) && ! is.na(refset_version)) {
      refset_version_msg = paste0(' (v', refset_version, ')')
    }
    cat("Reference data set: ", object@ref_key, refset_version_msg, "\n", sep = "")
    if(object@maf[, .N] > 0) {
      cat("Samples:\n")
      print(object@samples[, .(num_samples = .N), by = "coverage"], row.names = F)
      num_snvs = object@maf[variant_type == "snv", .N]
      cat("\nMAF data: ", num_snvs, " SNVs loaded.\n", sep = "")
    }
    signature_sets = object@advanced$snv_signatures
    if (length(signature_sets) > 0) {
      signature_set_names = paste(names(signature_sets), collapse = ', ')
      cat("Mutational processes: Sample-level extraction of ", signature_set_names, " SNV signatures.\n", sep = "")
    }
    cat("Run history: See [CESAnalysis]$run_history.\n")
    cat(paste0("\n[Created in cancereffectsizeR, version ", object@advanced$version, ".]"))
  }
)

setClass("CompoundVariantSet", representation(snvs = "data.table", compounds = "data.table", sample_calls = "list",
                                              cesa_uid = "numeric", cesa_num_samples = "integer"))
setMethod("show", "CompoundVariantSet", function(object) {
  num_compound = object@compounds[, .N]
  num_snv = object@snvs[, .N]
  plural = ifelse(num_compound == 1, '', 's')
  plural2 = ifelse(num_snv == 1, '', 's')
  msg = paste0("CompoundVariantSet with ", num_compound, " compound variant", plural, " consisting of ", 
                num_snv, " SNV", plural2, ":\n")
  cat(msg)
  print(object@compounds, topn = 3)
})

setMethod("length", "CompoundVariantSet", function(x) {
  return(x@compounds[, .N])
})

# allow subsetting CompoundVariantSet by index or name (that is, row number/compound_name entry in compound table)
# returns a new CompoundVariantSet; fine interactively but maybe not efficient enough for high-throughput use
setMethod("[", "CompoundVariantSet", function(x, i , j, ..., drop) {
  compounds = `[`(x@compounds, i, j, nomatch = NULL, ...)
  sample_calls = x@sample_calls[compounds$compound_name]
  snvs = x@snvs[compounds$compound_name, on = 'compound_name']
  return(new("CompoundVariantSet", compounds = compounds, snvs = snvs, sample_calls = sample_calls,
             cesa_uid = x@cesa_uid, cesa_num_samples = x@cesa_num_samples))
})

# thanks to https://stackoverflow.com/a/26080137
as.list.CompoundVariantSet <-function(x) {
  lapply(seq_along(x), function(i) x[i])
}
setMethod("as.list", "CompoundVariantSet", as.list.CompoundVariantSet)


#' @export
.DollarNames.CompoundVariantSet <- function(x, pattern = "") {
  features = c("compounds", "snv_info", "samples_with", "definitions")
  grep(pattern, features, value=TRUE)
}

setMethod("$", "CompoundVariantSet",
  function(x, name)
  {
    if(name == "snv_info") {
      return(x@snvs)
    } else if (name == "compounds") {
      return(x@compounds)
    } else if (name == "samples_with") {
      return(x@sample_calls)
    }else if (name == "definitions") {
      tmp = x@snvs[, .(snvs = list(snv_id)), by = "compound_name"]
      return(stats::setNames(tmp$snvs, tmp$compound_name))
    }
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

