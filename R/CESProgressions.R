#' CESProgressions
#'
#' @param samples
#' @param sample_progressions
#' @param ordered_progression_labels
#'
#' @return
#' @export
#'
#' @examples
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
