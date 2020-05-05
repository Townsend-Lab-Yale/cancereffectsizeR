#' suggest_cosmic_v3_signatures_to_remove
#'
#' Get a list of suggested signatures_to_remove for trinuc_mutation_rates.
#' For details, run vignette("cosmic_cancer_type_note").
#' 
#' @param cancer_type Run vignette("cosmic_cancer_type_note") for possible cancer type labels
#' @param treatment_naive give TRUE if samples were taken pre-treatment; FALSE or leave NULL otherwise
#' @param quiet (default false) for non-interactive use, suppress explanations and advice
#' @return a string of signatures to feed to signatures_to_remove
#' @export
suggest_cosmic_v3_signatures_to_remove = function(cancer_type = NULL, treatment_naive = NULL, quiet = FALSE) {
  data_source = paste0(system.file("extdata", package = "cancereffectsizeR"), '/pcawg_tcga_cancer_types.txt')
  dt = data.table::fread(data_source)
  to_remove = character()
  if(is.null(cancer_type)) {
    if(! quiet) {
      message("SBS25 (dubious, specific to Hodgkin's lymphoma cell lines) can usually be excluded.")
      message("Specify cancer_type for suggestions on which signatures do not appear in specific cancers.")
      message("(Run vignette(\"cosmic_cancer_type_note\") to see if there is this function has data for your cancer type.)")
    }
    to_remove = "SBS25"
  } 
  if(is.null(treatment_naive)) {
    if(! quiet) {
    message(paste0("If your samples were taken before treatment (e.g., TCGA data), re-run with treatment_naive = TRUE \nto get a list of ",
                   "signatures that are chemotherapy-associated, which you may want to exclude."))
    }
    treatment_naive = FALSE
  }
  
  if(! is.null(cancer_type)) {
    if(length(cancer_type) != 1 || ! is.character(cancer_type)) {
      stop("cancer_type should be a 1-length character vector")
    }
    # handle the only case where two different TCGA studies are the same PCAWG group
    if(cancer_type == "COAD" || cancer_type == "READ") {
      cancer_type = "ColoRect-AdenoCA"
    }
    index = which(dt$PCAWG == cancer_type)
    if (length(index) != 1) {
      index = which(dt$Applicable_TCGA == cancer_type)
    }
    if (length(index) != 1) {
      message(paste0("Input cancer_type not recognized.\n",
                    "Run vignette(\"cosmic_cancer_type_note\") and find your cancer type in the table. If ",
                    "it's not there, then there is no cancer-type-specific data available."))
      stop()
    }
    to_remove = names(which(unlist(dt[index,]) == 0))
    if(! quiet) {
      message(crayon::black(paste0("The following signatures were absent in ",
                                   dt[index, Number_of_tumors], " tumors in ", cancer_type, " in Alexandrov 2020:")))
      cat(to_remove, sep = ", ")
      cat("\n")
    }
  }
  
  treatment_sigs = c("SBS11", "SBS31", "SBS32", "SBS35")
  if(treatment_naive == TRUE) {
    if(! quiet) {
      message(crayon::black(paste0("\nThe following signatures were are associated with various cancer drugs:\n",
                                   "SBS11 (alkylating agents, like temozolomide), SBS31 (platinum drugs), SBS32 ",
                                   "(azathioprine, immunosuppressant), SBS35 (platinum drugs)")))     
    }
    to_remove = c(to_remove, treatment_sigs)
  }
  if(! quiet) {
    sig_string = paste0("signatures_to_remove = c(\"", paste(to_remove, collapse = "\", \""), "\")")
    message(crayon::black("\nIf you want to make all suggested exclusions: "))
    message(crayon::black(sig_string))   
  }
  return(invisible(to_remove))
}
