#' Plot step-specific cancer effects
#'
#' Visualize cancer effects estimated separately across ordered progression stages, as produced
#' by \code{ces_variant_step()}. Works for any number of stages, not just two: one panel is drawn
#' per variant (or gene, or other grouping column), with stage on the x-axis.
#'
#' @param effects Cancer effects table produced by \code{ces_variant_step()}. You can combine
#'   multiple such tables via rbind() to plot multiple runs together.
#' @param group_by Facet the plot by this column: "variant_name" (default), "gene", or another
#'   column present in \code{effects}.
#' @param stage_order Optional character vector giving stage display order (matching the
#'   \code{si_<stage>} column suffixes in \code{effects}). Defaults to the order those columns
#'   appear in the table, which matches the stage order used in the original
#'   \code{ces_variant_step()} call.
#' @param stage_labels Optional named character vector to relabel stages for display, e.g.
#'   \code{c(Pre = "Normal tissue", Pri = "Primary tumor")}. Names should match the stage names
#'   (i.e., \code{si_<stage>} suffixes, or \code{stage_order} if supplied).
#' @param show_ci TRUE/FALSE to depict confidence intervals (error bars) in the plot (default
#'   TRUE; ignored with a message if \code{effects} has no ci_low/ci_high columns).
#' @param lrt Optional output of \code{step_selection_LRT()}. When supplied, each facet is
#'   annotated with a significance code (\code{*}/\code{**}/\code{***}; \code{ns} if not
#'   significant) based on p_value. Requires \code{group_by} to identify a single row of
#'   \code{lrt} per facet (true when \code{group_by} is "variant_name", or "gene" when each gene
#'   has exactly one variant/compound in \code{effects}).
#' @param ncol Number of facet columns (default: ggplot chooses).
#' @param color_by Set to "stage" (default) to color points by stage (viridis discrete scale), a
#'   single R color to use throughout, or NULL to disable coloring.
#' @param viridis_option Viridis color map option, used when \code{color_by = "stage"}.
#' @param title Main plot title (default none).
#' @param x_title,y_title Axis titles.
#' @param legend.position Passed to ggplot's legend.position (none, left, right, top, bottom).
#' @return A ggplot
#' @export
plot_effects_step = function(effects, group_by = 'variant_name', stage_order = NULL, stage_labels = NULL,
                             show_ci = TRUE, lrt = NULL, ncol = NULL, color_by = 'stage',
                             viridis_option = 'cividis', title = '', x_title = NULL, y_title = NULL,
                             legend.position = 'right') {
  if (! require("ggplot2")) {
    stop("Package ggplot2 is needed for plotting. Run install.packages('ggplot2').")
  }
  if (! is(effects, 'data.table')) {
    stop('effects should be a data.table of cancer effects, as produced by ces_variant_step().')
  }
  effects = data.table::copy(effects)
  if (! is.character(group_by) || length(group_by) != 1 || ! group_by %in% names(effects)) {
    stop('group_by should be a 1-length character naming a column in effects.')
  }
  if (! identical(show_ci, TRUE) && ! identical(show_ci, FALSE)) {
    stop('show_ci should be TRUE/FALSE.')
  }

  si_cols = grep('^si_', names(effects), value = TRUE)
  if (length(si_cols) < 2) {
    stop('effects does not look like ces_variant_step() output (no si_<stage> columns found).')
  }
  stage_names = sub('^si_', '', si_cols)
  if (is.null(stage_order)) {
    stage_order = stage_names
  } else if (! setequal(stage_order, stage_names)) {
    stop('stage_order must contain exactly the stage names found in effects: ',
         paste(stage_names, collapse = ', '), '.')
  }

  ci_low_cols = grep('^ci_low_[0-9.]+_si_', names(effects), value = TRUE)
  has_ci = show_ci && length(ci_low_cols) > 0
  if (show_ci && length(ci_low_cols) == 0) {
    message('effects has no ci_low_*/ci_high_* columns (was conf set in ces_variant_step()?); omitting error bars.')
  }

  # reshape from one row per variant (wide, one column per stage) to one row per variant-stage (long)
  non_stage_cols = grep('^si_|^ci_low_|^ci_high_', names(effects), value = TRUE)
  id_cols = setdiff(names(effects), non_stage_cols)
  long_list = lapply(stage_names, function(st) {
    out = effects[, ..id_cols]
    out[, stage := st]
    out[, selection_intensity := effects[[paste0('si_', st)]]]
    if (has_ci) {
      low_col = grep(paste0('^ci_low_[0-9.]+_si_', st, '$'), names(effects), value = TRUE)
      high_col = grep(paste0('^ci_high_[0-9.]+_si_', st, '$'), names(effects), value = TRUE)
      out[, ci_low := if (length(low_col) == 1) effects[[low_col]] else NA_real_]
      out[, ci_high := if (length(high_col) == 1) effects[[high_col]] else NA_real_]
    }
    return(out)
  })
  long = data.table::rbindlist(long_list)

  if (! is.null(stage_labels)) {
    if (is.null(names(stage_labels))) {
      stop('stage_labels should be a named character vector (names = stage names, values = display labels).')
    }
    long[, stage := unname(ifelse(stage %in% names(stage_labels), stage_labels[stage], stage))]
    stage_order = unname(ifelse(stage_order %in% names(stage_labels), stage_labels[stage_order], stage_order))
  }
  long[, stage := factor(stage, levels = stage_order)]

  if (! is.null(lrt)) {
    if (! is(lrt, 'data.table') || ! all(c('variant_name', 'p_value') %in% names(lrt))) {
      stop('lrt should be output from step_selection_LRT() (a data.table with variant_name and p_value).')
    }
    long[lrt, p_value := i.p_value, on = 'variant_name']
    long[, significance := cut(p_value, breaks = c(-Inf, .001, .01, .05, Inf),
                               labels = c('***', '**', '*', 'ns'))]
  }

  p = ggplot2::ggplot(long, ggplot2::aes(x = stage, y = selection_intensity))
  if (identical(color_by, 'stage')) {
    p = p + ggplot2::geom_point(ggplot2::aes(color = stage), size = 2.5) +
      ggplot2::scale_color_viridis_d(option = viridis_option)
  } else if (! is.null(color_by)) {
    p = p + ggplot2::geom_point(color = color_by, size = 2.5)
  } else {
    p = p + ggplot2::geom_point(size = 2.5)
  }
  if (has_ci) {
    p = p + ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.25)
  }

  facet_formula = stats::as.formula(paste0('~', group_by))
  p = p + ggplot2::facet_wrap(facet_formula, ncol = ncol, scales = 'free_y')

  if (! is.null(lrt)) {
    label_data = long[, .(y = max(c(ci_high, selection_intensity), na.rm = TRUE) * 1.05,
                          label = as.character(significance[1])), by = group_by]
    label_data[, stage := stage_order[1]]
    p = p + ggplot2::geom_text(data = label_data, ggplot2::aes(x = stage, y = y, label = label),
                               inherit.aes = FALSE, size = 5)
  }

  x_axis_title = if (is.null(x_title)) 'Progression stage' else x_title
  y_axis_title = if (is.null(y_title)) 'Cancer effect (scaled selection coefficient)' else y_title

  p = p + ggplot2::labs(x = x_axis_title, y = y_axis_title, title = title, color = 'Stage') +
    ggplot2::theme_bw() +
    ggplot2::expand_limits(y = 0) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = 'lightgrey', linetype = 'dotted') +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = legend.position)

  return(p)
}
