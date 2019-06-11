setClass("CESList", representation(maf = "data.frame",  reference.mismatch.maf = "data.frame", non.snv.maf = "data.frame",
         pred.mnv.maf = "data.frame", trinucleotide_mutation_weights = "list"), )

## To-do: Add generic methods for editing slots that enfore these validity checks
## Currently, validity is only checked on object creation, so downstream functions can mess up the object
setValidity("CESList", 
    function(object) {
      maf.colnames = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
      
      problems = character()
      if (nrow(object@maf) == 0) {
        problems = c(problems, "CESList must be provided MAF data")
      } else if (! identical(colnames(object@maf), maf.colnames)) {
        problems = c(problems, "Invalid column names in object@maf")
      }
      
      if (nrow(object@reference.mismatch.maf) != 0 && ! identical(colnames(object@reference.mismatch.maf), maf.colnames)) {
        problems = c(problems, "Invalid column names in object@reference.mismatch.maf")
      }
      if (nrow(object@non.snv.maf) != 0 && ! identical(colnames(object@non.snv.maf), maf.colnames)) {
        problems = c(problems, "Invalid column names in object@non.snv.maf")
      }
      if (nrow(object@pred.mnv.maf) != 0 && ! identical(colnames(object@pred.mnv.maf), maf.colnames)) {
        problems = c(problems, "Invalid column names in object@pred.mnv.maf")
      }
      
      if (length(problems) > 0) {
        problems
      } else{
        TRUE
      }
    }
)
