#' CESProgressions
#'
#' @param samples a vector of unique patient identifiers
#' @param sample_progressions a vector of tumor progression stages corresonding to samples
#' @param ordered_progression_labels vector of progression stages in chronological order (e.g. c("1", "2", "3", ))
#'
#' @return a CESProgressions object for keeping track of tumor stages in a data set
#'
CESProgressions = function(samples = NULL, sample_progressions = NULL, ordered_progression_labels) {
  progression_order = 1:length(ordered_progression_labels)
  names(progression_order) = ordered_progression_labels
  progression_by_tumor = new.env()
  tumors_by_progression = new.env()
  num_tumors = rep(0, length(progression_order))
  
  progressions = new("CESProgressions", order = progression_order, progression_by_tumor = progression_by_tumor,
                     tumors_by_progression = tumors_by_progression, num_tumors = num_tumors)
  if (! is.null(samples)) {
    progressions = add_samples_to_CESProgressions(progressions = progressions, samples = samples, sample_progressions = sample_progressions)
  }
  
  return(progressions)
}

#' add_samples_to_CESProgressions
#' @param progressions a CESProgressions object
#' @param samples a vector of unique patient identifiers
#' @param sample_progressions a vector of tumor progression stages corresonding to samples
#' @return CESProgressions object with samples/progressions added

add_samples_to_CESProgressions = function(progressions, samples, sample_progressions) {
  tumors_by_progression = progressions@tumors_by_progression
  progression_by_tumor = progressions@progression_by_tumor
  progression_order = progressions@order
  
  if(length(samples) != length(sample_progressions)) {
    stop("Samples and progressions vectors are not equal-length.")
  }
  
  already_samples = sum(progressions@num_tumors) > 0 # for clearer error messages, check if progressions object already has samples
  for (i in 1:length(samples)) {
    tumor = as.character(samples[i])
    progression_label = sample_progressions[i]
    if (! progression_label %in% progression_order) {
      stop(paste0("Error: Unexpected progressions stage (\"", progression_label, "\") in MAF data"))
    }
    progression_number = progression_order[[progression_label]]
    if (exists(tumor, progression_by_tumor) && progression_by_tumor[[tumor]] != progression_number) {
      if (already_samples) {
        stop(paste0("Sample ", tumor, " has multiple progression stages in the data. It's possible that it's ",
            "listed as being different stages in different data sources."))
      } else {
        stop(paste0("Sample ", tumor, " has multiple progression stages in input data."))
      }
    } else {
      progression_by_tumor[[tumor]] = progression_number
    }
    if (exists(as.character(progression_number), tumors_by_progression) && ! tumor %in% tumors_by_progression[[as.character(progression_number)]]) {
      tumors_by_progression[[as.character(progression_number)]] = c(tumor, tumors_by_progression[[as.character(progression_number)]])
    } else if (! exists(as.character(progression_number), tumors_by_progression)) {
      tumors_by_progression[[as.character(progression_number)]] = tumor
    }
  }
  
  progressions@progression_by_tumor = progression_by_tumor
  progressions@tumors_by_progression = tumors_by_progression
  progressions@num_tumors = sapply(progressions@order, function(x) length(get_progression_tumors(progressions, x)))
  return(progressions)
}




