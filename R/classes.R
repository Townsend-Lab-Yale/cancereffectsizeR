setClass("CESProgressions", representation(order = "character", by.tumor = "environment"))
setValidity("CESProgressions",
  function(object) {
    problems = character()
    if (length(object@order) == 0) {
      problems = c(problems, "Need at least one tumor progression stage for effect size analysis")
    }
    if (length(ls(object@by.tumor)) == 0) {
      problems = c(problems, "No tumors in tumor progression dictionary")
    } else {
        invalid_progressions = character()
        for (tumor in ls(envir = object@by.tumor)) {
            progression = get(tumor, envir = object@by.tumor)
            if (! is.character(progression)) {
                problems = c(problems, "tumor progression dictionary contains non-character values")
                break
            }
            if(length(progression) != 1) {
              problems = c(problems, "tumor progression dictionary values should be character vectors of length 1")
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

setClass("CESAnalysis", representation(mutations.maf = "MAFdf",  reference.mismatch.maf = "MAFdf", snv.maf = "MAFdf",
         pred.mnv.maf = "MAFdf", trinucleotide_mutation_weights = "list", tumor.progressions = "CESProgressions",
         mutrates_list = "list", dndscv_out_list = "list", relative_substitution_rates = "table", refcds_data = "array"))


## To-do: Add generic methods for editing slots that enfore these validity checks
## Currently, validity is only checked on object creation, so downstream functions can mess up the object
setValidity("CESAnalysis", 
    function(object) {
      ## Can add more validation later
      problems = character()
      if (nrow(object@mutations.maf) == 0) {
        problems = c(problems, "CESAnalysis must be provided MAF data")
      }

      if (length(problems) > 0) {
        problems
      } else{
        TRUE
      }
    }
)
