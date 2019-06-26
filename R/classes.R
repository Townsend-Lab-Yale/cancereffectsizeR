setClass("CESProgressions", representation(order = "numeric", progression_by_tumor = "environment", 
                                            tumors_by_progression = "environment"))
setValidity("CESProgressions",
  function(object) {
    problems = character()
    if (length(object@order) == 0) {
      problems = c(problems, "Need at least one tumor progression stage for effect size analysis")
    }
    if (length(ls(object@progression_by_tumor)) == 0) {
      problems = c(problems, "No tumors in tumor progression dictionary")
    } else {
        invalid_progressions = character()
        for (tumor in ls(envir = object@progression_by_tumor)) {
            progression = get(tumor, envir = object@progression_by_tumor)
            if (! is.numeric(progression)) {
                problems = c(problems, "tumor progression dictionary contains non-numeric values")
                break
            }
            if(length(progression) != 1) {
              problems = c(problems, "tumor progression dictionary values should be numeric vectors of length 1")
              break
            }
            if (! progression %in% object@order) {
              invalid_progressions = c(invalid_progressions, progression)
            }
        }
        if (length(invalid_progressions) > 0) {
          if (length(invalid_progressions) > 10) {
            num_additional = length(invalid_progressions) - 10
            invalid_progressions = paste0(paste(invalid_progressions[1:10], collapse=", "), " (and ", num_additional, " more)") 
          } else {
            invalid_progressions = paste(invalid_progressions, collapse=", ")
          }
          problems = c(problems, paste0("Tumor progression dictionary includes progression stages not found in ",
                      "tumor progression order:\n", invalid_progressions))
        }
    }
    if (length(problems) > 0) {
      problems
    } else{
      TRUE
    }
  }
)

CESProgressions = function(samples, sample_progressions, ordered_progression_labels) {
  if(length(samples) != length(sample_progressions)) {
    stop("Samples and progressions vectors are not equal-length.")
  }
  progression_order = 1:length(ordered_progression_labels)
  names(progression_order) = ordered_progression_labels
  progression_by_tumor = new.env()
  tumors_by_progression = new.env()
  for (i in 1:length(samples)) {
    tumor = as.character(samples[i])
    progression_label = sample_progressions[i]
    progression_number = progression_order[[progression_label]]
    if (exists(tumor, progression_by_tumor) && progression_by_tumor[[tumor]] != progression_number) {
      stop(paste0("Tumor ", tumor, " has multiple progression states in the input MAF."))
    } else {
      progression_by_tumor[[tumor]] = progression_number
    }
    if (exists(as.character(progression_number), tumors_by_progression) && ! tumor %in% tumors_by_progression[[as.character(progression_number)]]) {
      tumors_by_progression[[as.character(progression_number)]] = c(tumor, tumors_by_progression[[as.character(progression_number)]])
    } else if (! exists(as.character(progression_number), tumors_by_progression)) {
      tumors_by_progression[[as.character(progression_number)]] = tumor
    }
  }

  progressions = new("CESProgressions", order = progression_order, progression_by_tumor = progression_by_tumor,
                          tumors_by_progression = tumors_by_progression) 
  return(progressions)
}

setGeneric("get_progression_number", function(a, b) standardGeneric("get_progression_number")) 
setMethod("get_progression_number", signature("CESProgressions", "character"), function(a, b) {
  progressions = a; tumor = b
  return(sapply(tumor, function(x) { progressions@progression_by_tumor[[x]]}))
})

setGeneric("get_progression_name", function(a, b) standardGeneric("get_progression_name"))
setMethod("get_progression_name", signature("CESProgressions", "character"), function(a, b) {
  progressions = a; tumor = b
  return(sapply(tumor, function(x) { names(progressions@order)[[progressions@progression_by_tumor[[x]]]]}))
})

setGeneric("get_progression_tumors", function(a, b) standardGeneric("get_progression_tumors"))
setMethod("get_progression_tumors", signature("CESProgressions", "numeric"), 
  function(a, b) {
    progressions = a; progression_number = b
    return(progressions@tumors_by_progression[[as.character(b)]])
  }
)
setMethod("get_progression_tumors", signature("CESProgressions", "character"),
  function(a, b) {
    progressions = a; progression_name = b
    return(progressions@tumors_by_progression[[as.character(progressions@order[[progression_name]])]])
  }
)



#' A class that extends data.frame to enforce consistent handling of MAF data
setClass("MAFdf", contains = "data.frame")
setValidity("MAFdf",
  function(object) {
    # add column validation later (e.g., valid chromosomes, alleles)
    maf_colnames = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
    if(! identical(colnames(object)[1:5], maf_colnames)) {
      "Illegal MAF column names"
    } else {
      TRUE
    }
  }
)


MAFdf = function(df, ...) {
  new("MAFdf", df, ...)
}

setClass("CESAnalysis", representation(main.maf = "MAFdf", annotated.snv.maf = "MAFdf", trinucleotide_mutation_weights = "list", 
          progressions = "CESProgressions", mutrates_list = "list", dndscv_out_list = "list",
          relative_substitution_rates = "table", refcds_data = "array", excluded = "environment", selection_results = "list"))

## To-do: Add generic methods for editing slots that enfore these validity checks
## Currently, validity is only checked on object creation, so downstream functions can mess up the object
setValidity("CESAnalysis", 
    function(object) {
      ## Can add more validation later
      problems = character()
      if (nrow(object@main.maf) == 0) {
        problems = c(problems, "CESAnalysis must be provided MAF data")
      }

      if (length(problems) > 0) {
        problems
      } else{
        TRUE
      }
    }
)
