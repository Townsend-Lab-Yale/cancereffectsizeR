#' Tissue-specific mutational signature exclusions
#'
#' Get suggestions on signatures_to_remove for trinuc_mutation_rates for COSMIC signatures v3 and later.
#' For details, see \code{vignette("cosmic_cancer_type_note")}.
#' 
#' @param cancer_type See chart on website for supported cancer type labels.
#' @param treatment_naive give TRUE if samples were taken pre-treatment; FALSE or leave
#'   NULL otherwise.
#' @param quiet (default false) for non-interactive use, suppress explanations and advice.
#' @return a vector of signatures to feed to the \code{trinuc_mutation_rates()}
#'   \code{signature_exclusions} argument.
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
  
  if(require("ces.refset.hg38", character.only = T)) {
    signature_set_path = system.file("refset/signatures/COSMIC_v3.2_signatures.rds", package = "ces.refset.hg38")
  } else if(require("ces.refset.hg19", character.only = T)) {
    # Use v3.2 if refset package is up-to-date; otherwise, use v3.1.
    signature_set_path = system.file("refset/signatures/COSMIC_v3.2_signatures.rds", package = "ces.refset.hg19")
    if (signature_set_path == '') {
      signature_set_path = system.file("refset/signatures/COSMIC_v3.1_signatures.rds", package = "ces.refset.hg19")
    }
  } else {
    message("To use this function, install ces.refset.h38 or ces.refset.hg19:")
    message("remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg38@*release\")")
    message("remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg19@*release\")")
    stop("Required reference data not available; install a CES refset data package.")
  }
  
  sig_metadata = readRDS(signature_set_path)$meta
  original_sig_order = copy(sig_metadata$Signature)
  setkey(sig_metadata, "Signature")
  
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
                    "See \"Cancer type considerations for COSMIC signatures\" on the cancereffectsizeR website and find your cancer type in the table. If ",
                    "it's not there, then there is no cancer-type-specific data available."))
      stop()
    }
    to_remove = c(names(which(unlist(dt[index,]) == 0)))
    to_remove = intersect(to_remove, sig_metadata$Signature) # compatibility with v3.1, when using
    if(! quiet) {
      pretty_message(paste0("The following signatures are suggested absent in ", cancer_type, ", either by Alexandrov 2020 or the COSMIC signature website:"))
      print(sig_metadata[to_remove, .(Signature, Etiology)])
      cat("\n")
    }
  }
  
  treatment_sigs = c("SBS11", "SBS31", "SBS32", "SBS35", "SBS86", "SBS87", "SBS90")
  if(treatment_naive == TRUE) {
    if(! quiet) {
      cat("The following signatures are associated with various treatments:\n")
      print(sig_metadata[treatment_sigs, .(Signature, Etiology)])
    }
    to_remove = c(to_remove, treatment_sigs)
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

#' Previous name of suggest_cosmic_signature_exclusions
#' 
#' This function was renamed suggest_cosmic_signature_exclusions(). For now, calling it
#' with the old name will simply pass the arguments along to the new version. It will
#' eventually be removed.
#' 
#' @export
#' @keywords internal
suggest_cosmic_signatures_to_remove = function(...) {
  message("Calling suggest_cosmic_signature_exclusions().")
  message("The alias suggest_cosmic_signatures_to_remove() will itself be removed in a future update.")
  suggest_cosmic_signature_exclusions(...)
}