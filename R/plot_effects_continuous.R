#' Plot continuous-covariate cancer effects
#'
#' Visualize cancer-effect estimates from the linear, logistic, and generalized
#' sigmoid continuous-selection models. The function can display all three
#' continuous models together, select the best-fitting model by AIC, or display
#' any user-specified subset of the three models. It can return the continuous
#' plot alone or place a discrete-cohort plot beside the continuous plot.
#'
#' Point-estimate runs are read directly from \code{cesa@selection_results}
#' using the supplied run names.
#'
#' The three continuous models are
#'
#' \deqn{\gamma(x) = \beta_0 + \beta_1 x}
#'
#' for the linear model,
#'
#' \deqn{\gamma(x) = \frac{\exp(L)}{1 + \exp[-k(x-m)]}}
#'
#' for the logistic model, and
#'
#' \deqn{
#' \gamma(x) =
#' \exp(C) + [\exp(L)-\exp(C)]
#' \frac{x^s}{x^s+m^s}
#' }
#'
#' for the generalized sigmoid model.
#'
#' When \code{continuous_models = "best"}, AIC is calculated separately for
#' each variant as \eqn{2k - 2\ell}, using 2, 3, and 4 fitted parameters for
#' the linear, logistic, and sigmoid models, respectively.
#'
#' When \code{output = "cohort_continuous"}, the cohort and continuous panels
#' use the same y-axis limits. The shared range is determined from the cohort
#' point estimates and confidence intervals and from the continuous curve(s)
#' that are actually displayed.
#'
#' @param cesa A \code{CESAnalysis} object containing the requested selection
#'   runs.
#' @param linear_run_name Run name for point estimates from
#'   \code{ces_variant_linear()}.
#' @param logistic_run_name Run name for point estimates from
#'   \code{ces_variant_logistic()}.
#' @param sigmoid_run_name Run name for point estimates from
#'   \code{ces_variant_sigmoid()}.
#' @param cohort_run_names Optional character vector of run names for discrete
#'   cohort-specific effect estimates. Required when
#'   \code{output = "cohort_continuous"}. If the vector is named, its names are
#'   used as the cohort labels and their order is preserved; otherwise the run
#'   names themselves are used as labels. Cohort runs must contain
#'   \code{ci_low_95} and \code{ci_high_95} for the plotted 95\% confidence
#'   intervals.
#' @param variants Optional character vector of variant IDs or variant names to
#'   plot. If \code{NULL}, all variants present in all three continuous-model
#'   runs are plotted.
#' @param covariate_col Name of the numeric continuous-covariate column in
#'   \code{cesa@samples}.
#' @param covariate_range Optional numeric vector of length two giving the
#'   minimum and maximum covariate values over which to draw the fitted curves.
#'   If \code{NULL}, the observed range of \code{covariate_col} is used.
#' @param x_range Numeric vector of length two giving the displayed x-axis
#'   limits for the continuous panel. Defaults to \code{c(0, 100)}.
#' @param n_covariate_points Number of points used to draw each continuous curve.
#' @param output Either \code{"continuous"} for the continuous plot alone or
#'   \code{"cohort_continuous"} for discrete-cohort and continuous panels side
#'   by side.
#' @param continuous_models Which continuous models to display. Use
#'   \code{"best"} for the lowest-AIC model, \code{"all"} for all three models,
#'   or a character vector containing one or more of \code{"linear"},
#'   \code{"logistic"}, and \code{"sigmoid"}.
#' @param title Optional plot title. With multiple variants, a supplied title is
#'   used as a prefix; otherwise each variant name is used.
#' @param x_title X-axis title for the continuous plot. If \code{NULL},
#'   \code{covariate_col} is used.
#' @param y_title Y-axis title used for both continuous and cohort plots.
#' @param cohort_x_title X-axis title for the discrete-cohort plot.
#' @param legend.position Position of the continuous-model line legend. The
#'   default places the legend inside the continuous panel.
#' @param cohort_colors Colors for the discrete cohorts. The default provides
#'   three colors; supply a vector with at least one color per cohort when more
#'   cohorts are plotted. A named vector may be used to assign colors by cohort
#'   label.
#' @param model_colors Named character vector giving colors for
#'   \code{"Linear"}, \code{"Logistic"}, and \code{"Sigmoid"}.
#'
#' @return For one variant, a ggplot object when \code{output = "continuous"}
#'   or a patchwork object when \code{output = "cohort_continuous"}. For
#'   multiple variants, a named list of such objects. The best-fitting model
#'   and AIC values are stored as attributes on each returned plot.
#'
#' @seealso \code{\link{plot_effects}}, \code{\link{ces_variant_linear}},
#'   \code{\link{ces_variant_logistic}}, \code{\link{ces_variant_sigmoid}}
#'
#' @family continuous selection models
#'
#' @export
plot_effects_continuous <- function(
    cesa,
    linear_run_name,
    logistic_run_name,
    sigmoid_run_name,
    cohort_run_names = NULL,
    variants = NULL,
    covariate_col,
    covariate_range = NULL,
    x_range = c(0, 100),
    n_covariate_points = 100,
    output = c("continuous", "cohort_continuous"),
    continuous_models = "best",
    title = NULL,
    x_title = NULL,
    y_title = "Strength of selection",
    cohort_x_title = "Cohort",
    legend.position = c(0.82, 0.58),
    cohort_colors = c("#E71F19", "#F4A016", "#2B6A99"),
    model_colors = c(
      "Linear" = "#C2A5CF",
      "Logistic" = "#2B6A99",
      "Sigmoid" = "#F4A016"
    )
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  output <- match.arg(output)

  if (output == "cohort_continuous" &&
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required for side-by-side plots. ",
      "Please install it with install.packages('patchwork')."
    )
  }

  if (!is(cesa, "CESAnalysis")) stop("cesa should be a CESAnalysis.")
  if (!is.character(covariate_col) || length(covariate_col) != 1L) {
    stop("covariate_col should be 1-length character.")
  }
  if (!is.numeric(n_covariate_points) || length(n_covariate_points) != 1L ||
      n_covariate_points < 2 || n_covariate_points != as.integer(n_covariate_points)) {
    stop("n_covariate_points should be an integer >= 2.")
  }

  valid_models <- c("linear", "logistic", "sigmoid")
  if (!is.character(continuous_models) || length(continuous_models) < 1L) {
    stop("continuous_models should be 'best', 'all', or one or more model names.")
  }
  continuous_models <- tolower(continuous_models)
  if (length(continuous_models) == 1L && continuous_models == "best") {
    model_choice <- "best"
  } else if (length(continuous_models) == 1L && continuous_models == "all") {
    model_choice <- "all"
  } else {
    if (any(!continuous_models %in% valid_models)) {
      stop(
        "continuous_models should be 'best', 'all', or a character vector containing ",
        "'linear', 'logistic', and/or 'sigmoid'."
      )
    }
    model_choice <- unique(continuous_models)
  }

  if (is.null(x_title)) x_title <- covariate_col

  get_run <- function(run_name, arg_name) {
    if (!is.character(run_name) || length(run_name) != 1L) {
      stop(arg_name, " should be 1-length character.")
    }
    if (!run_name %in% names(cesa@selection_results)) {
      stop(
        "Selection run '", run_name, "' supplied through ", arg_name,
        " was not found in cesa@selection_results."
      )
    }
    data.table::as.data.table(
      data.table::copy(cesa@selection_results[[run_name]])
    )
  }

  linear_results <- get_run(linear_run_name, "linear_run_name")
  logistic_results <- get_run(logistic_run_name, "logistic_run_name")
  sigmoid_results <- get_run(sigmoid_run_name, "sigmoid_run_name")

  required_point_cols <- list(
    linear = c("beta0", "beta1", "loglikelihood"),
    logistic = c("L", "k", "m", "loglikelihood"),
    sigmoid = c("C", "L", "s", "m", "loglikelihood")
  )
  point_tables <- list(
    linear = linear_results,
    logistic = logistic_results,
    sigmoid = sigmoid_results
  )
  for (model_name in names(point_tables)) {
    missing_cols <- setdiff(
      required_point_cols[[model_name]],
      names(point_tables[[model_name]])
    )
    if (length(missing_cols) > 0L) {
      stop(
        "Run for ", model_name, " model is missing required columns: ",
        paste(missing_cols, collapse = ", "), "."
      )
    }
  }

  if (is.null(covariate_range)) {
    if (!covariate_col %in% names(cesa@samples)) {
      stop(
        "Column '", covariate_col,
        "' was not found in cesa@samples. Supply covariate_col or an explicit covariate_range."
      )
    }
    observed_covariate <- suppressWarnings(
      as.numeric(as.character(cesa@samples[[covariate_col]]))
    )
    observed_covariate <- observed_covariate[is.finite(observed_covariate)]
    if (length(observed_covariate) == 0L) {
      stop("No finite numeric values were found in cesa@samples$", covariate_col, ".")
    }
    covariate_range <- range(observed_covariate)
  } else if (!is.numeric(covariate_range) || length(covariate_range) != 2L ||
             any(!is.finite(covariate_range)) || covariate_range[1] >= covariate_range[2]) {
    stop("covariate_range should contain two finite increasing numeric values.")
  }

  if (!is.numeric(x_range) || length(x_range) != 2L ||
      any(!is.finite(x_range)) || x_range[1] >= x_range[2]) {
    stop("x_range should contain two finite increasing numeric values.")
  }

  covariate_values <- seq(
    covariate_range[1],
    covariate_range[2],
    length.out = n_covariate_points
  )

  get_variant_row <- function(dt, target_variant) {
    out <- NULL
    if ("variant_id" %in% names(dt)) {
      out <- dt[dt[["variant_id"]] %in% target_variant]
    }
    if ((is.null(out) || nrow(out) == 0L) &&
        "variant_name" %in% names(dt)) {
      out <- dt[dt[["variant_name"]] %in% target_variant]
    }
    if (is.null(out) || nrow(out) == 0L) return(NULL)
    out[1]
  }

  resolve_variant <- function(x) {
    row <- get_variant_row(linear_results, x)
    if (is.null(row)) return(NULL)
    id <- if ("variant_id" %in% names(row) && !is.na(row$variant_id[1])) {
      as.character(row$variant_id[1])
    } else {
      as.character(row$variant_name[1])
    }
    label <- if ("variant_name" %in% names(row) && !is.na(row$variant_name[1])) {
      as.character(row$variant_name[1])
    } else {
      id
    }
    list(id = id, label = label)
  }

  if (is.null(variants)) {
    ids_for <- function(dt) {
      if ("variant_id" %in% names(dt)) {
        as.character(dt$variant_id)
      } else {
        as.character(dt$variant_name)
      }
    }
    variant_values <- Reduce(
      intersect,
      list(
        ids_for(linear_results),
        ids_for(logistic_results),
        ids_for(sigmoid_results)
      )
    )
    variants_resolved <- lapply(variant_values, resolve_variant)
  } else {
    if (!is.character(variants)) {
      stop("variants should be NULL or a character vector of variant IDs/names.")
    }
    variants_resolved <- lapply(variants, resolve_variant)
    missing_variants <- variants[
      vapply(variants_resolved, is.null, logical(1))
    ]
    if (length(missing_variants) > 0L) {
      stop(
        "The following variants were not found in the linear run: ",
        paste(missing_variants, collapse = ", "), "."
      )
    }
  }
  variants_resolved <- Filter(Negate(is.null), variants_resolved)
  if (length(variants_resolved) == 0L) {
    stop("No variants were available to plot.")
  }

  curve_values <- function(model, row, x) {
    if (model == "linear") {
      y <- row$beta0[1] + row$beta1[1] * x
    } else if (model == "logistic") {
      y <- exp(row$L[1]) / (1 + exp(-row$k[1] * (x - row$m[1])))
    } else {
      lower <- exp(row$C[1])
      upper <- exp(row$L[1])
      y <- lower + (upper - lower) * x^row$s[1] /
        (x^row$s[1] + row$m[1]^row$s[1])
    }
    y[!is.finite(y)] <- NA_real_
    y
  }

  model_aic <- function(linear_row, logistic_row, sigmoid_row) {
    loglikelihoods <- c(
      linear = as.numeric(linear_row$loglikelihood[1]),
      logistic = as.numeric(logistic_row$loglikelihood[1]),
      sigmoid = as.numeric(sigmoid_row$loglikelihood[1])
    )
    n_parameters <- c(
      linear = 2,
      logistic = 3,
      sigmoid = 4
    )
    aic <- 2 * n_parameters - 2 * loglikelihoods
    aic[!is.finite(aic)] <- Inf
    if (all(is.infinite(aic))) {
      stop("All three AIC values are non-finite for a plotted variant.")
    }
    list(values = aic, best = names(which.min(aic)))
  }

  get_cohort_data <- function(variant_id, variant_label) {
    if (is.null(cohort_run_names)) return(NULL)
    cohort_labels <- names(cohort_run_names)
    if (is.null(cohort_labels) || any(cohort_labels == "")) {
      cohort_labels <- cohort_run_names
    }

    out <- vector("list", length(cohort_run_names))
    for (ii in seq_along(cohort_run_names)) {
      dt <- get_run(cohort_run_names[ii], "cohort_run_names")
      row <- get_variant_row(dt, variant_id)
      if (is.null(row)) row <- get_variant_row(dt, variant_label)
      if (is.null(row)) next

      required_cohort_cols <- c(
        "selection_intensity",
        "ci_low_95",
        "ci_high_95"
      )
      missing_cohort_cols <- setdiff(
        required_cohort_cols,
        names(row)
      )
      if (length(missing_cohort_cols) > 0L) {
        stop(
          "Cohort run '", cohort_run_names[ii],
          "' is missing required columns: ",
          paste(missing_cohort_cols, collapse = ", "),
          "."
        )
      }

      prevalence <- NA_real_
      if (all(c("included_with_variant", "included_total", "held_out") %in% names(row))) {
        eligible <- as.numeric(row$included_total[1]) +
          as.numeric(row$held_out[1])
        if (is.finite(eligible) && eligible > 0) {
          prevalence <- as.numeric(row$included_with_variant[1]) / eligible
        }
      }

      out[[ii]] <- data.table::data.table(
        cohort = cohort_labels[ii],
        selection_intensity = as.numeric(row$selection_intensity[1]),
        ci_low = as.numeric(row$ci_low_95[1]),
        ci_high = as.numeric(row$ci_high_95[1]),
        prevalence = prevalence
      )
    }

    out <- data.table::rbindlist(out, fill = TRUE)
    if (!nrow(out)) return(NULL)
    out[, cohort := factor(cohort, levels = cohort_labels)]
    out
  }

  # Manuscript-style scientific notation, e.g. 2 x 10^6.
  pretty_sci <- function(x) {
    labels <- vapply(x, function(value) {
      if (!is.finite(value)) return("NA")
      if (value == 0) return("0")
      scientific <- format(
        value,
        scientific = TRUE,
        trim = TRUE,
        digits = 2
      )
      parts <- strsplit(scientific, "e", fixed = TRUE)[[1]]
      coefficient <- sub("\\.?0+$", "", parts[1])
      exponent <- as.integer(parts[2])
      paste0(coefficient, "%*%10^", exponent)
    }, character(1))
    parse(text = labels)
  }

  make_one_plot <- function(v) {
    variant_id <- v$id
    variant_label <- v$label

    linear_row <- get_variant_row(linear_results, variant_id)
    if (is.null(linear_row)) {
      linear_row <- get_variant_row(linear_results, variant_label)
    }
    logistic_row <- get_variant_row(logistic_results, variant_id)
    if (is.null(logistic_row)) {
      logistic_row <- get_variant_row(logistic_results, variant_label)
    }
    sigmoid_row <- get_variant_row(sigmoid_results, variant_id)
    if (is.null(sigmoid_row)) {
      sigmoid_row <- get_variant_row(sigmoid_results, variant_label)
    }

    if (any(vapply(
      list(linear_row, logistic_row, sigmoid_row),
      is.null,
      logical(1)
    ))) {
      stop(
        "Variant ", variant_label,
        " is not present in all three continuous-model runs."
      )
    }

    aic_info <- model_aic(
      linear_row,
      logistic_row,
      sigmoid_row
    )
    best_model <- aic_info$best

    if (identical(model_choice, "best")) {
      models_to_plot <- best_model
    } else if (identical(model_choice, "all")) {
      models_to_plot <- valid_models
    } else {
      models_to_plot <- model_choice
    }

    model_labels <- c(
      linear = "Linear",
      logistic = "Logistic",
      sigmoid = "Sigmoid"
    )
    point_rows <- list(
      linear = linear_row,
      logistic = logistic_row,
      sigmoid = sigmoid_row
    )

    curve_data <- data.table::rbindlist(
      lapply(models_to_plot, function(model_name) {
        data.table::data.table(
          covariate_value = covariate_values,
          effect = curve_values(
            model_name,
            point_rows[[model_name]],
            covariate_values
          ),
          model = model_labels[[model_name]]
        )
      })
    )

    cohort_data <- if (output == "cohort_continuous") {
      get_cohort_data(variant_id, variant_label)
    } else {
      NULL
    }
    if (output == "cohort_continuous" && is.null(cohort_data)) {
      stop("No cohort data were available for variant ", variant_label, ".")
    }

    y_values <- curve_data$effect
    if (!is.null(cohort_data)) {
      y_values <- c(
        y_values,
        cohort_data$selection_intensity,
        cohort_data$ci_low,
        cohort_data$ci_high
      )
    }
    y_values <- y_values[is.finite(y_values)]
    y_max <- max(y_values, na.rm = TRUE)
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    y_upper <- y_max
    y_breaks <- pretty(c(0, y_upper), n = 4)
    y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_upper]
    if (!0 %in% y_breaks) y_breaks <- c(0, y_breaks)

    plot_title <- if (is.null(title)) {
      variant_label
    } else if (length(variants_resolved) == 1L) {
      title
    } else {
      paste0(title, ": ", variant_label)
    }

    x_breaks <- pretty(x_range, n = 5)
    x_breaks <- x_breaks[
      x_breaks >= x_range[1] &
        x_breaks <= x_range[2]
    ]

    gg_cont <- ggplot2::ggplot()

    line_alpha <- c(
      "Linear" = 0.8,
      "Logistic" = 0.8,
      "Sigmoid" = 0.7
    )

    for (model_name in model_labels[models_to_plot]) {
      dd <- curve_data[curve_data$model == model_name, ]
      gg_cont <- gg_cont +
        ggplot2::geom_line(
          data = dd,
          ggplot2::aes(
            x = covariate_value,
            y = effect,
            color = model
          ),
          linewidth = 0.5,
          alpha = line_alpha[[model_name]],
          na.rm = TRUE,
          show.legend = TRUE
        )
    }

    gg_cont <- gg_cont +
      ggplot2::scale_color_manual(
        values = model_colors,
        breaks = model_labels[models_to_plot],
        name = NULL
      ) +
      ggplot2::scale_x_continuous(
        limits = x_range,
        breaks = x_breaks,
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, y_upper),
        breaks = y_breaks,
        labels = pretty_sci(y_breaks),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::labs(
        title = if (output == "continuous") plot_title else NULL,
        x = x_title,
        y = y_title
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_text(size = 15),
        plot.title = ggplot2::element_text(size = 15),
        text = ggplot2::element_text(size = 15),
        legend.position = legend.position,
        legend.justification = c(0, 0.5),
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 15)
      )

    if (output == "continuous") {
      attr(gg_cont, "best_model") <- best_model
      attr(gg_cont, "AIC") <- aic_info$values
      return(gg_cont)
    }

    cohort_levels <- levels(cohort_data$cohort)

    if (is.null(names(cohort_colors)) || any(names(cohort_colors) == "")) {
      if (length(cohort_colors) < length(cohort_levels)) {
        stop(
          "cohort_colors must provide at least one color per cohort."
        )
      }
      cohort_colors_use <- stats::setNames(
        cohort_colors[seq_along(cohort_levels)],
        cohort_levels
      )
    } else {
      missing_colors <- setdiff(cohort_levels, names(cohort_colors))
      if (length(missing_colors) > 0L) {
        stop(
          "cohort_colors is missing colors for: ",
          paste(missing_colors, collapse = ", "),
          "."
        )
      }
      cohort_colors_use <- cohort_colors[cohort_levels]
    }

    gg_cohort <- ggplot2::ggplot(
      cohort_data,
      ggplot2::aes(x = cohort, y = selection_intensity)
    )

    have_prevalence <- any(is.finite(cohort_data$prevalence))

    if (have_prevalence) {
      cohort_data[, prevalence_label :=
                    ifelse(
                      is.finite(prevalence),
                      paste0(format(prevalence * 100, digits = 2), "%"),
                      ""
                    )
      ]

      gg_cohort <- gg_cohort +
        ggplot2::geom_point(
          data = cohort_data,
          ggplot2::aes(size = prevalence, fill = cohort),
          shape = 21,
          color = "gray20",
          na.rm = TRUE,
          show.legend = FALSE
        ) +
        ggplot2::geom_text(
          data = cohort_data,
          ggplot2::aes(label = prevalence_label),
          vjust = -1.2,
          hjust = -0.4,
          size = 4,
          na.rm = TRUE,
          show.legend = FALSE
        ) +
        ggplot2::scale_size_continuous(
          range = c(5, 10),
          guide = "none"
        )
    } else {
      gg_cohort <- gg_cohort +
        ggplot2::geom_point(
          ggplot2::aes(fill = cohort),
          shape = 21,
          color = "gray20",
          size = 5,
          na.rm = TRUE,
          show.legend = FALSE
        )
    }

    # Draw cohort confidence intervals after the points so the CI lines
    # appear on top of the point estimates.
    gg_cohort <- gg_cohort +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = ci_low,
          ymax = ci_high
        ),
        width = 0.2,
        na.rm = TRUE
      )

    gg_cohort <- gg_cohort +
      ggplot2::scale_fill_manual(
        values = cohort_colors_use,
        guide = "none"
      ) +
      ggplot2::scale_x_discrete(labels = cohort_levels) +
      ggplot2::scale_y_continuous(
        limits = c(0, y_upper),
        breaks = y_breaks,
        labels = pretty_sci(y_breaks),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::labs(
        title = plot_title,
        x = cohort_x_title,
        y = y_title
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_text(size = 15),
        plot.title = ggplot2::element_text(size = 15),
        text = ggplot2::element_text(size = 15),
        legend.position = "none"
      )

    combined <- gg_cohort + gg_cont +
      patchwork::plot_layout(
        widths = c(1.15, 2.3),
        guides = "keep"
      )
    attr(combined, "best_model") <- best_model
    attr(combined, "AIC") <- aic_info$values
    combined
  }

  plots <- lapply(variants_resolved, make_one_plot)
  names(plots) <- vapply(
    variants_resolved,
    function(x) x$label,
    character(1)
  )
  if (length(plots) == 1L) return(plots[[1]])
  plots
}
