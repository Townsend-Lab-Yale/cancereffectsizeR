#' Tissue-specific mutational signature exclusions
#'
#' Get suggestions on signatures_to_remove for trinuc_mutation_rates for COSMIC signatures v3 and later.
#' For details, see [this article](https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cosmic_cancer_type_note.html)
#' on our website.
#' 
#' @param cancer_type See [here](https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cosmic_cancer_type_note.html) 
#'   for supported cancer type labels.
#' @param treatment_naive give TRUE if samples were taken pre-treatment; FALSE or leave
#'   NULL otherwise.
#' @param quiet (default false) for non-interactive use, suppress explanations and advice.
#' @return a vector of signatures to feed to the \code{trinuc_mutation_rates()}
#'   \code{signature_exclusions} argument.
#' @md
#' @export
suggest_cosmic_signature_exclusions = function(cancer_type = NULL, treatment_naive = NULL, quiet = FALSE) {
  data_source = paste0(system.file("extdata", package = "cancereffectsizeR"), '/COSMIC_v3.2_signatures_by_cancer_type.txt')
  dt = data.table::fread(data_source)
  to_remove = character()
  if(is.null(cancer_type)) {
    if(! quiet) {
      pretty_message(paste0("SBS25, SBS89, and SBS91 can usually be excluded."), black = F)
      message("Specify cancer_type (see documentation) for suggestions on which signatures do not appear in specific cancers.")
    }
    to_remove = c("SBS25", "SBS89", "SBS91")
  } 
  if(is.null(treatment_naive)) {
    if(! quiet) {
    message(paste0("If your samples were taken before treatment (e.g., TCGA data), re-run with treatment_naive = TRUE \nto get a list of ",
                   "signatures that are chemotherapy-associated, which you may want to exclude."))
    }
    treatment_naive = FALSE
  }
  
  if (! identical(quiet, TRUE) && ! identical(quiet, FALSE) ) {
    stop("argument quiet should be T/F")
  }
  
  sig_metadata = cosmic_signature_info()
  suppressWarnings({sig_metadata$short_name = NULL})
  original_sig_order = copy(sig_metadata$name)
  setkey(sig_metadata, "name")
  
  if(! is.null(cancer_type)) {
    if(length(cancer_type) != 1 || ! is.character(cancer_type)) {
      stop("cancer_type should be a 1-length character vector")
    }
    # handle the only case where two different TCGA studies are the same PCAWG group
    if(cancer_type == "COAD" || cancer_type == "READ") {
      cancer_type = "ColoRect-AdenoCA"
    }
    index = which(dt$PCAWG %ilike% cancer_type)
    if (length(index) != 1) {
      index = which(dt$Applicable_TCGA %ilike% cancer_type)
    }
    if (length(index) != 1) {
      message(paste0("Input cancer_type not recognized.\n",
                    "See \"Cancer type considerations for COSMIC signatures\" on the cancereffectsizeR website and find your cancer type in the table. If ",
                    "it's not there, then there is no cancer-type-specific data available."))
      stop()
    }
    to_remove = c(names(which(unlist(dt[index,]) == 0)))
    to_remove = intersect(to_remove, sig_metadata$name) # compatibility with v3.1, when using
    if(! quiet) {
      pretty_message(paste0("The following signatures are suggested absent in ", cancer_type, ", either by Alexandrov 2020 or the COSMIC signature website:"))
      print(sig_metadata[to_remove])
      cat("\n")
    }
  }
  
  treatment_sigs = c("SBS11", "SBS31", "SBS25", "SBS32", "SBS35", "SBS86", "SBS87", "SBS90", "SBS99")
  if(treatment_naive == TRUE) {
    if(! quiet) {
      cat("The following signatures are associated with various treatments:\n")
      print(sig_metadata[treatment_sigs])
    }
    to_remove = c(to_remove, treatment_sigs)
  }
  
  # COSMIC v3.4 split SBS22 into 22a,b and SBS40 into 40a,b,c. We'll apply the same tissue exclusions.
  if('SBS22' %in% to_remove) {
    to_remove = c(to_remove, c("SBS22a", 'SBS22b'))
  }
  if('SBS40' %in% to_remove) {
    to_remove = c(to_remove, c("SBS40a", 'SBS40b', 'SBS40c'))
  }
  
  # make unique and put signatures in numeric order
  to_remove = unique(to_remove)
  
  to_remove = original_sig_order[original_sig_order %in% to_remove]
  if(! quiet) {
    sig_string = paste0("signature_exclusions = c(\"", paste(to_remove, collapse = "\", \""), "\")")
    message(crayon::black("\nSilently returning the following suggested exclusions: "))
    message(crayon::black(sig_string))
  }
  return(invisible(to_remove))
}
