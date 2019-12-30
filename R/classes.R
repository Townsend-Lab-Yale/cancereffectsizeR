setClass("CESProgressions", representation(order = "numeric", progression_by_tumor = "environment",
                                            tumors_by_progression = "environment", num_tumors = "numeric"))
setValidity("CESProgressions",
  function(object) {
    problems = character()
    if (length(object@order) == 0) {
      problems = c(problems, "Need at least one tumor progression stage for effect size analysis")
    }
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
      problems = c(problems, paste0("Input data lists tumor progression stages not found in ",
                  "CESAnalysis object's tumor progression order:\n", invalid_progressions))
    }
    if (length(problems) > 0) {
      problems
    } else{
      TRUE
    }
  }
)



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




setClass("CESAnalysis", representation(maf = "data.frame", annotated.snv.maf = "data.frame", trinucleotide_mutation_weights = "list",
          progressions = "CESProgressions", mutrates_list = "list", dndscv_out_list = "list",
          excluded = "data.frame", selection_results = "list", coverage = "list",
          genome = "BSgenome"))

setValidity("CESAnalysis",
    function(object) {
     # add validation later
     TRUE
    }
)
