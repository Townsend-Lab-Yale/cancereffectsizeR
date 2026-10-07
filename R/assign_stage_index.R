#' Assign samples to ordered progression stages
#'
#' Builds the \code{sample_index} table required by \code{step_selection_lik()} (and, in turn,
#' \code{ces_variant_step()}), associating each sample with an integer index and a display name
#' for its position in an ordered sequence of tumor progression stages (for example, normal
#' tissue, primary tumor, metastasis).
#'
#' @param cesa CESAnalysis object
#' @param stage_col Name of a sample-level data column (as added via \code{load_maf()} or
#'   \code{load_sample_data()}) that records each sample's stage.
#' @param stage_order Describes the values of \code{stage_col} in earliest-to-latest stage order.
#'   Either a character vector with one value of \code{stage_col} per stage (stage names default
#'   to these same values), or a named list where each element gives one-or-more values of
#'   \code{stage_col} that should be grouped into the same stage, with the element's name used as
#'   the stage's display name.
#' @param samples Which samples to include. Defaults to all samples in the CESAnalysis. Can be a
#'   vector of Unique_Patient_Identifiers, or a data.table containing rows from the CESAnalysis
#'   sample table.
#' @return A data.table with columns Unique_Patient_Identifier, group_index (1 = earliest stage),
#'   and group_name.
#' @export
assign_stage_index = function(cesa = NULL, stage_col = NULL, stage_order = NULL, samples = character()) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (! is.character(stage_col) || length(stage_col) != 1) {
    stop("stage_col should be 1-length character naming a sample data column.")
  }
  sample_info = select_samples(cesa, samples)
  if (! stage_col %in% names(sample_info)) {
    stop("stage_col \"", stage_col, "\" is not a column in the CESAnalysis samples table.")
  }
  if (is.null(stage_order) || length(stage_order) < 2) {
    stop("stage_order should give the values of stage_col for at least 2 stages, in ",
         "earliest-to-latest order (see documentation).")
  }

  # Normalize stage_order into a named list, one element per stage
  if (is.null(names(stage_order))) {
    if (length(unlist(stage_order)) == length(stage_order)) {
      # plain character vector: stage names default to the values themselves
      names(stage_order) = as.character(unlist(stage_order))
    } else {
      stop("If stage_order groups multiple values per stage, it must be a named list ",
           "(so each stage has a display name).")
    }
  } else if (any(names(stage_order) == "")) {
    stop("If stage_order is a named list, every element must be named.")
  }

  index_by_value = list()
  name_by_value = list()
  for (i in seq_along(stage_order)) {
    values_i = as.character(unlist(stage_order[[i]]))
    for (v in values_i) {
      if (v %in% names(index_by_value)) {
        stop("Value \"", v, "\" of stage_col appears in more than one stage of stage_order.")
      }
      index_by_value[[v]] = i
      name_by_value[[v]] = names(stage_order)[i]
    }
  }

  stage_values = as.character(sample_info[[stage_col]])
  unrecognized = setdiff(unique(stage_values), names(index_by_value))
  if (length(unrecognized) > 0) {
    stop("Some values of ", stage_col, " are not covered by stage_order: ",
         paste(unrecognized, collapse = ", "), ".")
  }

  sample_index = data.table::data.table(
    Unique_Patient_Identifier = sample_info$Unique_Patient_Identifier,
    group_index = as.integer(unlist(index_by_value[stage_values])),
    group_name = unlist(name_by_value[stage_values])
  )
  data.table::setkey(sample_index, "Unique_Patient_Identifier")
  return(sample_index)
}
