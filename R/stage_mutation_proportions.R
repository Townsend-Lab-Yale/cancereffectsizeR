#' Calculate per-gene, per-stage mutation rate proportions
#'
#' For use with \code{step_selection_lik()}/\code{ces_variant_step()}. Takes cumulative per-gene
#' mutation rates estimated separately for each ordered progression stage (each from a
#' \code{gene_mutation_rates()} call restricted to the samples at that stage;
#' see the step-specific selection vignette) and converts them into the proportion of a gene's
#' total mutation rate that is estimated to accumulate during each stage.
#'
#' By default (\code{rate_source = "dndscv_expected"}), each stage's cumulative rate is dNdScv's
#' covariate-based expected synonymous mutation count for the gene (\code{exp_syn_cv}), divided by
#' the gene's number of synonymous sites and by the number of samples in the dNdScv run. These are
#' uncorrected rates: unlike the rates that \code{gene_mutation_rates()} assigns, they are not
#' adjusted toward each gene's observed synonymous mutation count. Using the same kind of estimate
#' at every stage keeps stage proportions internally consistent (the adjustment is skipped
#' altogether when dNdScv's overdispersion parameter theta is below 1, which is common for sparse
#' early-stage samples). This requires full dNdScv output, so run each stage's
#' \code{gene_mutation_rates()} call with \code{save_all_dndscv_output = TRUE}. Alternatively,
#' \code{rate_source = "gene_rates"} uses the rates in \code{get_gene_rates(cesa)} as-is (for
#' example, rates supplied with \code{set_gene_rates()}).
#'
#' \code{ces_variant_step()} expects every sample to carry the same, final-stage cumulative gene
#' rate, which the step-specific model then divides among stages using these proportions. After
#' calculating proportions, clear the per-stage rates with \code{clear_gene_rates()} and assign the
#' final-stage cumulative rate (the last \code{rate_<stage name>} column of the output) to all
#' samples with \code{set_gene_rates()}.
#'
#' Cumulative rates should be non-decreasing across stages, since later stages by definition
#' cover at least as much mutational "time" as earlier ones. When dN/dS-based rate estimates are
#' noisy, a later stage's estimated cumulative rate can come out below an earlier one for a given
#' gene, which would produce an invalid (negative) proportion. When this happens, affected genes
#' are flagged with a warning, and handled according to \code{on_invalid}.
#'
#' @param cesa CESAnalysis object with gene rates already calculated once per stage, via repeated
#'   calls to \code{gene_mutation_rates()} (or \code{set_gene_rates()}), each restricted to the
#'   samples at that stage.
#' @param rate_cols Character vector of gene rate groups, one per stage, in earliest-to-latest
#'   order (for example, \code{c("rate_grp_1", "rate_grp_2")}). These name dNdScv runs in
#'   \code{cesa$dNdScv_results} (when \code{rate_source = "dndscv_expected"}) or columns of
#'   \code{get_gene_rates(cesa)} (when \code{rate_source = "gene_rates"}). Defaults to all
#'   \code{rate_grp_*} groups present, in ascending numeric order.
#' @param rate_source \code{"dndscv_expected"} (default) for uncorrected rates calculated from
#'   full dNdScv output, or \code{"gene_rates"} to use \code{get_gene_rates(cesa)} as-is. See
#'   details.
#' @param stage_names Optional character vector of display names for the stages, in the same
#'   order as \code{rate_cols}. Defaults to \code{rate_cols}.
#' @param on_invalid What to do when a gene's cumulative rates are not non-decreasing across
#'   stages: \code{"floor"} (default) floors negative/undefined stage contributions to
#'   \code{floor_prop} (of the gene's total rate) and rescales that gene's proportions to sum to
#'   1, with a warning; \code{"error"} stops; \code{"NA"} sets all of that gene's proportions to
#'   \code{NA_real_}.
#' @param floor_prop Minimum stage proportion used when \code{on_invalid = "floor"} (default
#'   1e-6).
#' @return A data.table with a gene (or pid) identifier column, one proportion column per stage
#'   (named \code{p_<stage_name>}), and one cumulative rate column per stage (named
#'   \code{rate_<stage_name>}). Proportions in each row sum to 1, except rows set to NA under
#'   \code{on_invalid = "NA"}.
#' @export
stage_mutation_proportions = function(cesa = NULL, rate_cols = NULL, stage_names = NULL,
                                       rate_source = "dndscv_expected",
                                       on_invalid = "floor", floor_prop = 1e-6) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (! is.character(rate_source) || length(rate_source) != 1 || ! rate_source %in% c("dndscv_expected", "gene_rates")) {
    stop('rate_source should be "dndscv_expected" or "gene_rates".')
  }
  if (! is.character(on_invalid) || length(on_invalid) != 1 || ! on_invalid %in% c("floor", "error", "NA")) {
    stop('on_invalid should be "floor", "error", or "NA".')
  }
  if (! is.numeric(floor_prop) || length(floor_prop) != 1 || is.na(floor_prop) || floor_prop < 0 || floor_prop >= 1) {
    stop("floor_prop should be a single number in [0, 1).")
  }

  if (rate_source == "gene_rates") {
    gene_rates = get_gene_rates(cesa)
    available_groups = names(gene_rates)
    id_col = if ("gene" %in% names(gene_rates)) "gene" else "pid"
  } else {
    dndscv_out = cesa@dndscv_out_list
    available_groups = names(dndscv_out)
    id_col = if ("pid" %in% names(get_gene_rates(cesa))) "pid" else "gene"
  }

  if (is.null(rate_cols)) {
    rate_cols = grep("^rate_grp_[0-9]+$", available_groups, value = TRUE)
    if (length(rate_cols) > 0) {
      rate_cols = rate_cols[order(as.integer(sub("^rate_grp_", "", rate_cols)))]
    }
  } else if (! is.character(rate_cols)) {
    stop("rate_cols should be a character vector of gene rate group names.")
  }
  if (length(rate_cols) < 2) {
    stop("Need at least two stage rate groups to calculate stage proportions (one ",
         "gene_mutation_rates()/set_gene_rates() call per stage; see documentation).")
  }
  missing_cols = setdiff(rate_cols, available_groups)
  if (length(missing_cols) > 0) {
    where = if (rate_source == "gene_rates") "get_gene_rates(cesa)" else "cesa$dNdScv_results"
    stop("rate_cols not present in ", where, ": ", paste(missing_cols, collapse = ", "), ".")
  }

  if (rate_source == "dndscv_expected") {
    # Same RefCDS that gene_mutation_rates() passes to dNdScv
    RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS.dndscv
    if (is.null(RefCDS)) {
      RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS
    }
    if (is.null(RefCDS)) {
      stop("Reference data for ", cesa@ref_key, " is not loaded (load the CESAnalysis with load_cesa()).")
    }
    rates_by_stage = lapply(rate_cols, function(grp) {
      output = dndscv_out[[grp]]
      if (is.data.table(output) || is.null(output$genemuts) || is.null(output$annotmuts)) {
        stop("Full dNdScv output is not available for ", grp, ". Run gene_mutation_rates() with ",
             "save_all_dndscv_output = TRUE for each stage, or use rate_source = \"gene_rates\".")
      }
      nsyn_sites = sapply(RefCDS[output$genemuts$gene_name], function(x) colSums(x[["L"]])[1])
      num_samples = uniqueN(output$annotmuts$sampleID)
      out = data.table(id = output$genemuts$gene_name, rate = output$genemuts$exp_syn_cv / nsyn_sites / num_samples)
      setnames(out, "rate", grp)
      return(out)
    })
    gene_rates = Reduce(function(x, y) merge(x, y, by = "id"), rates_by_stage)
    setnames(gene_rates, "id", id_col)
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
  colnames(cum_rates) = paste0("rate_", stage_names)
  result = data.table::data.table(gene_rates[[id_col]])
  data.table::setnames(result, "V1", id_col)
  result = cbind(result, data.table::as.data.table(props), data.table::as.data.table(cum_rates))
  return(result)
}
