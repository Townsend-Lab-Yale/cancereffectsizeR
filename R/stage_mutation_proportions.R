#' Calculate per-gene, per-stage mutation rate proportions
#'
#' For use with \code{step_selection_lik()}/\code{ces_variant_step()}. Takes cumulative per-gene
#' mutation rates estimated separately for each ordered progression stage (each from a
#' \code{gene_mutation_rates()} call restricted to samples that have reached at least that stage;
#' see the step-specific selection vignette) and converts them into the proportion of a gene's
#' total mutation rate that is estimated to accumulate during each stage.
#'
#' Cumulative rates should be non-decreasing across stages, since later stages by definition
#' cover at least as much mutational "time" as earlier ones. When dN/dS-based rate estimates are
#' noisy, a later stage's estimated cumulative rate can come out below an earlier one for a given
#' gene, which would produce an invalid (negative) proportion. When this happens, affected genes
#' are flagged with a warning, and handled according to \code{on_invalid}.
#'
#' @param cesa CESAnalysis object with gene rates already calculated once per stage, via repeated
#'   calls to \code{gene_mutation_rates()} (or \code{set_gene_rates()}), each restricted to the
#'   samples that have reached at least that stage.
#' @param rate_cols Character vector of column names in \code{get_gene_rates(cesa)}, one per
#'   stage, in earliest-to-latest order (for example, \code{c("rate_grp_1", "rate_grp_2")}).
#'   Defaults to all \code{rate_grp_*} columns present, in ascending numeric order.
#' @param stage_names Optional character vector of display names for the stages, in the same
#'   order as \code{rate_cols}. Defaults to \code{rate_cols}.
#' @param on_invalid What to do when a gene's cumulative rates are not non-decreasing across
#'   stages: \code{"floor"} (default) floors negative/undefined stage contributions to
#'   \code{floor_prop} (of the gene's total rate) and rescales that gene's proportions to sum to
#'   1, with a warning; \code{"error"} stops; \code{"NA"} sets all of that gene's proportions to
#'   \code{NA_real_}.
#' @param floor_prop Minimum stage proportion used when \code{on_invalid = "floor"} (default
#'   1e-6).
#' @return A data.table with a gene (or pid) identifier column, then one proportion column per
#'   stage (named \code{p_<stage_name>}). Each row sums to 1, except rows set to NA under
#'   \code{on_invalid = "NA"}.
#' @export
stage_mutation_proportions = function(cesa = NULL, rate_cols = NULL, stage_names = NULL,
                                       on_invalid = "floor", floor_prop = 1e-6) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (! is.character(on_invalid) || length(on_invalid) != 1 || ! on_invalid %in% c("floor", "error", "NA")) {
    stop('on_invalid should be "floor", "error", or "NA".')
  }
  if (! is.numeric(floor_prop) || length(floor_prop) != 1 || is.na(floor_prop) || floor_prop < 0 || floor_prop >= 1) {
    stop("floor_prop should be a single number in [0, 1).")
  }

  gene_rates = get_gene_rates(cesa)
  id_col = if ("gene" %in% names(gene_rates)) "gene" else "pid"

  if (is.null(rate_cols)) {
    rate_cols = grep("^rate_grp_[0-9]+$", names(gene_rates), value = TRUE)
    if (length(rate_cols) > 0) {
      rate_cols = rate_cols[order(as.integer(sub("^rate_grp_", "", rate_cols)))]
    }
  } else if (! is.character(rate_cols)) {
    stop("rate_cols should be a character vector of column names.")
  }
  if (length(rate_cols) < 2) {
    stop("Need at least two stage rate columns to calculate stage proportions (one ",
         "gene_mutation_rates()/set_gene_rates() call per stage; see documentation).")
  }
  missing_cols = setdiff(rate_cols, names(gene_rates))
  if (length(missing_cols) > 0) {
    stop("rate_cols not present in get_gene_rates(cesa): ", paste(missing_cols, collapse = ", "), ".")
  }
  if (is.null(stage_names)) {
    stage_names = rate_cols
  } else if (! is.character(stage_names) || length(stage_names) != length(rate_cols)) {
    stop("stage_names must be a character vector the same length as rate_cols.")
  }

  num_stages = length(rate_cols)
  cum_rates = as.matrix(gene_rates[, rate_cols, with = FALSE])
  storage.mode(cum_rates) = "double"

  # successive differences give each stage's own (non-cumulative) contribution to the total rate;
  # the first stage's contribution is just its cumulative rate
  stage_rate = cum_rates
  if (num_stages > 1) {
    stage_rate[, 2:num_stages] = cum_rates[, 2:num_stages, drop = FALSE] -
      cum_rates[, 1:(num_stages - 1), drop = FALSE]
  }

  final_total = cum_rates[, num_stages]
  row_invalid = function(i) {
    row = stage_rate[i, ]
    any(is.na(row)) || any(row < 0) || is.na(final_total[i]) || final_total[i] <= 0
  }
  invalid = vapply(seq_len(nrow(cum_rates)), row_invalid, logical(1))

  if (any(invalid)) {
    bad_genes = gene_rates[[id_col]][invalid]
    bad_genes_shown = bad_genes
    if (length(bad_genes_shown) > 30) {
      bad_genes_shown = c(bad_genes_shown[1:30], "...")
    }
    if (on_invalid == "error") {
      stop("Cumulative mutation rates are not non-decreasing across stages for ", length(bad_genes),
           " gene(s): ", paste(bad_genes_shown, collapse = ", "),
           ". Pass on_invalid = \"floor\" or \"NA\" instead of \"error\" to handle rather than stop.")
    }
    msg = paste0(length(bad_genes), " gene(s) had cumulative mutation rates that were not ",
                 "non-decreasing across stages, which is required (later stages should reflect ",
                 "at least as much mutational time as earlier ones); this usually reflects ",
                 "estimation noise: ", paste(bad_genes_shown, collapse = ", "), ". ")
    if (on_invalid == "floor") {
      msg = paste0(msg, "Negative/undefined stage contributions were floored to ", floor_prop,
                   " of the gene's total rate, and that gene's proportions were rescaled to sum to 1.")
    } else {
      msg = paste0(msg, "Proportions for these genes were set to NA.")
    }
    warning(msg, call. = FALSE)
  }

  props = matrix(NA_real_, nrow = nrow(cum_rates), ncol = num_stages)
  for (i in seq_len(nrow(cum_rates))) {
    if (invalid[i]) {
      if (on_invalid != "floor") {
        next # NA (either on_invalid == "NA", or "error" already stopped above)
      }
      total = final_total[i]
      if (is.na(total) || total <= 0) {
        next # nothing sensible to floor to; leave as NA
      }
      row = stage_rate[i, ]
      row[is.na(row) | row < 0] = floor_prop * total
      props[i, ] = row / sum(row)
    } else {
      props[i, ] = stage_rate[i, ] / final_total[i]
    }
  }

  colnames(props) = paste0("p_", stage_names)
  result = data.table::data.table(gene_rates[[id_col]])
  data.table::setnames(result, "V1", id_col)
  result = cbind(result, data.table::as.data.table(props))
  return(result)
}
