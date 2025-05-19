#' Parse user-supplied refset name
#' 
#' @param refset_name The name of a reference data package or directory.
#' @keywords internal
parse_refset_name = function(refset_name) {
  # Check for and load reference data for the chosen genome/transcriptome data
  if (! rlang::is_scalar_character(refset_name)) {
    stop("refset should be the name of an installed refset package or a path to a custom refset directory.")
  }
  using_custom_refset = TRUE
  if (refset_name %in% names(.official_refsets)) {
    using_custom_refset = FALSE
    if(file.exists(refset_name)) {
      stop("You've given the name of a CES reference data set package, but a file/folder with the same name is in your working directory. Stopping to avoid confusion.")
    }
    if(! require(refset_name, character.only = T)) {
      if(refset_name == "ces.refset.hg19") {
        message("Install ces.refset.hg19 like this:\n",
                "options(timeout = 600)\n",
                "remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg19@*release\")")
      } else if(refset_name == "ces.refset.hg38") {
        message("Install ces.refset.hg38 like this:\n",
                "options(timeout = 600)\n",
                "remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg38@*release\")")
      }
      stop("CES reference data set ", refset_name, " not installed.")
    }
    req_version = .official_refsets[[refset_name]]
    actual_version = packageVersion(refset_name)
    if (actual_version < req_version) {
      stop("CES reference data set ", refset_name, " is version ", actual_version, ", but your version of cancereffectsizeR requires at least ",
           "version ", req_version, ".\nRun this to update:\n",
           "remotes::install_github(\"Townsend-Lab-Yale/", refset_name, "\")")
    }
    refset_version = actual_version
    data_dir = system.file("refset", package = refset_name)
  } else {
    refset_version = NA_character_
    if (! dir.exists(refset_name)) {
      if (grepl('/', refset_name)) {
        stop("Could not find reference data at ", refset_name)
      } else {
        stop("Invalid reference set name. Check spelling, or view available data sets with list_ces_refsets().")
      }
    }
    
    data_dir = refset_name
    refset_name = basename(refset_name)
    if(refset_name %in% names(.official_refsets)) {
      stop("Your custom reference data set has the same name (", refset_name, ") as a CES reference data package. Please rename it.")
    }
  }
  return(list(refset_name = refset_name, data_dir = data_dir, refset_version = refset_version))
}
