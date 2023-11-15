#' Plot cancer effects
#'
#' Visualize and compare cancer effects for variants of interest.
#' 
#' @param effects Cancer effects table, as produced by \code{ces_variant()}. You can combine multiple tables via rbind()
#' to plot multiple effects per variant, such as to compare effects across subgroups.
#' @param topn Include up to this many variants. The highest-effect variants are plotted. (Or, if
#'   \code{group_by} is gene, include up to this many groups. Groups are ranked by their
#'   highest-effect variants.)
#' @param group_by If 'variant' (the default), one variant per row in the plot. If "gene" or some other column name, 
#'   variants will be plotted together accordingly. When "gene", some 
#' @param y_label Y-axis labels for each group of variants. By default ("auto"), will be variant names
#' when \code{group_by = "variant"}, the values in group_by otherwise.
#' @param color_by A single color to use for geom_point fill (default "darkseagreen4"). Or, the name of
#'   a column that specifies color groupings. Can be used to distinguish points when multiple effects
#'   are plotted per variant (for example, when comparing effects between subgroups), or to
#'   highlight related groups of variants. A viridis color scale will be applied, unless ever single value
#'   in the color column is interpretable as an R color, in which case the given colors will be used.
#' @param prevalence_method Show each variant's prevalence as a raw mutation count ("count", the default), or as
#'   a percentage of samples with sequencing coverage at the site ("percent"). If the effects table
#'   has the same number of samples covering every inference, you can choose "both".
#' @param color_label If color_by is supplying color names for scale_color_identity(), optionally include color_label
#' so that colors can be labeled in the plot legend. 
#' @param legend.position Passed to ggplot's legend.position (none, left, right, top, bottom, or
#'   coordinates). Use "none" to eliminate the legend. Defaults to "right".
#' @param legend_size_name The title for the point size scale (larger points = more prevalent variants).
#' @param legend_size_breaks Vector of specific mutation counts (or percentages) to depict in the point size legend.
#'   Specify numeric values if you don't like what gets produced by the default ("auto"). Set to
#'   FALSE or to a single desired point size to turn of size scaling.
#' @param legend_color_name The title for the point fill color scale.
#' @param viridis_option If using \code{color_by}, this argument
#' specifies which viridis color map to use. Ignored if you specify your own colors.
#' @param label_individual_variants When TRUE (default), individual variants within groups will be
#'   labeled when group_by is not 'variant'. Set FALSE to not label variants, or specify a column
#'   name that supplies a label for each row in the effects table. By default, variant names will be
#'   used for labels. If group_by is exactly "gene", labels will be shortened to just the amino acid
#'   changes.
#' @param order_by_effect When TRUE (default), variants are plotted in order of effect. When FALSE,
#'   variants are plotted top-down in the order they are supplied.
#' @param show_ci TRUE/FALSE to depict confidence intervals in plot (default TRUE).
#' @param title Main title for the plot (by default, no title)
#' @param x_title Text for the X-axis label.
#' @param y_title Text for the Y-axis label.
#' @return A ggplot
#' @export
plot_effects = function(effects, topn = 30, group_by = 'variant',
                        title = '',
                        x_title = NULL, y_title = NULL,
                        y_label = 'auto',
                        color_by = 'darkseagreen4', color_label = NULL,
                        legend.position = 'right',
                        legend_size_name = 'auto',
                        legend_color_name = NULL,
                        viridis_option = 'cividis',
                        legend_size_breaks = 'auto', 
                        label_individual_variants = TRUE,
                        order_by_effect = TRUE,
                        prevalence_method = 'auto',
                        show_ci = TRUE) {
  # Verify ggplot2/ggrepel are installed, since they are not package dependencies
  if (! require("ggplot2") || ! require("ggrepel")) {
    stop("Packages needed for plotting are not installed. Run install.packages(c('ggplot2', 'ggrepel')).")
  }
  
  # Validate viridis_option
  if(! is.character(viridis_option) || length(viridis_option) != 1) {
    stop('viridis_option should be 1-length character.')
  }
  
  # Determine axis titles
  if(is.null(x_title)) {
    x_axis_title = 'Cancer effect (scaled selection coefficient)'
  } else if (is.character(x_title) && length(x_title) == 1) {
    x_axis_title = x_title
  } else {
    stop('x_title should be 1-length character (or NULL for default title).')
  }
  if(is.null(y_title)) {
    y_axis_title = 'Somatic variant'
  } else if (is.character(y_title) && length(y_title) == 1) {
    y_axis_title = y_title
  } else {
    stop('y_title should be 1-length character (or NULL for default title).')
  }
  if(! is.character(title) || length(title) != 1)  {
    stop('title should be 1-length character')
  }
  
  # Validate that effects is table with required columns
  if (! is(effects, 'data.table')) {
    stop('effects should be a data.table of cancer effects.')
  }
  effects = copy(effects)
  
  # For compound variants, there is no variant_id, just variant_name
  if(! 'variant_id' %in% names(effects)) {
    effects[, variant_id := variant_name]
  }
  
  required_cols = c('variant_name', 'variant_type', 'selection_intensity', 'included_with_variant', 'held_out')
  
  if(identical(show_ci, TRUE)) {
    required_cols = c(required_cols, c('ci_low_95', 'ci_high_95'))
  } else if (! identical(show_ci, FALSE)) {
    stop('Argument show_ci should be TRUE/FALSE.')
  }
  
  missing_cols = setdiff(required_cols, names(effects))
  if(length(missing_cols) > 0) {
    if(identical(missing_cols, 'held_out')) {
      msg('Column held_out is missing from effects table. This column was added to effects output in cancereffectsizeR v2.8.0. ',
          'If you are trying to use effects loaded from an older analysis, the solution is to re-run the ces_variant() call.')
      warning(pretty_message(msg, emit = F))
    }
    stop("Missing required columns in effects table: ", paste(missing_cols, collapse = ', '), '.')
  }
  
  if(effects[, .N] == 0) {
    stop('effects table has zero rows.')
  }

  # group_by can be variant (default), gene (also gets special behavior), or any other character/factor column name.
  if(! is.character(group_by) || length(group_by) != 1) {
    stop('group_by should be 1-length character.')
  }
  
  if(group_by == 'variant') {
    # When grouping by variant, there is only 1 variant per variant group
    effects[, variant_group := variant_id]
    effects[, top_by_group := max(selection_intensity, na.rm = T), by = 'variant_group']
  } else {
    if(! group_by %in% names(effects)) {
      stop('Specified group_by column ', group_by, ' is not present in effects table.')
    }
    setnames(effects, group_by, 'variant_group')
    if(! is.character(effects$variant_group) && ! is.factor(effects$variant_group)) {
      effects$variant_group = as.factor(effects$variant_group)
    }
    
    has_na = effects[, anyNA(variant_group)]
    if (has_na) {
      effects = effects[! is.na(variant_group)]
      msg = paste0('Some variants in effects table have NA values in group_by column ', group_by, '. ',
                   'These variants have been filtered out. If you want to include these in the plot, ',
                   'assign them non-NA.')
      warning(pretty_message(msg, emit = F))
    }
    effects[, top_by_group := max(selection_intensity, na.rm = T), by = 'variant_group']
  }
  
  # Remove variants (or variant groups) outside of topn.
  if(! is.null(topn)) {
    if(! is.numeric(topn) || length(topn) != 1 || as.integer(topn) != topn) {
      stop('topn should be a positive integer.')
    }
    if(uniqueN(effects$top_by_group) > topn) {
      lowest_passing = sort(unique(effects$top_by_group), decreasing = T)[topn]
      effects = effects[top_by_group >= lowest_passing]
    }
  }
  
  # Deal with NA selection/CI
  lowest_label = NULL
  lowest_real = effects[included_with_variant > 0, min(ci_low_95, na.rm = TRUE)]
  even_lower = 10^floor(log10(lowest_real)) # rounding down to next factor of 10 below any lower CI
  values_to_check = unlist(effects[, .(selection_intensity, ci_low_95)])
  if(anyNA(values_to_check) || any(values_to_check < lowest_real)) {
    lowest_label = paste0('   <', format(even_lower, scientific = F, big.mark = ',')) # whitespace for aesthetics
    effects[is.na(ci_low_95) | ci_low_95 < lowest_real, ci_low_95 := even_lower]
    effects[is.na(selection_intensity) | selection_intensity < lowest_real, selection_intensity := even_lower]
  }
  
  # Sort into desired plot order (top group will be top of plot, with variants ordered by selection within groups).
  if (identical(order_by_effect, TRUE)) {
    effects = effects[order(top_by_group, selection_intensity)]
  } else if (identical(order_by_effect, FALSE)) {
    effects = effects[.N:1] # Reverse given order so that they are plotted in top-down order.
  } else {
    stop('Argument order_by_effect should be TRUE/FALSE.')
  }
  
  # Use the chosen prevalence method to scale variant point sizes
  effects[, num_samples := included_total + held_out]
  if(identical(prevalence_method, 'auto')) {
    # As stated in docs, we'll use count if sample numbers are similar enough (20%).
    prevalence_method = ifelse(effects[, max(num_samples)/min(num_samples)] > 1.2, 'percent', 'count')
    if(prevalence_method == 'percent') {
      pretty_message("Depicting variant prevalence as percent of eligible samples that have mutation. If you prefer counts, set prevalence_method = \"count\".",
                     black = F)
    }
  }
  if(identical(prevalence_method, 'count') || identical(prevalence_method, 'both')) {
    effects[, prevalence := included_with_variant]
    if(prevalence_method == 'both' && uniqueN(effects$num_samples) != 1) {
      msg = paste0('Not all variants have sequencing coverage in the same number of samples, so ',
                   'prevalance_method \"both\" can\'t be used.')
      stop(pretty_message(msg, emit = F))
    }
  } else if (identical(prevalence_method, 'percent')) {
    # Rounding to match the rounding that will be applied to labels, to 
    # avoid possibly getting the same rounded label for two different breaks.
    effects[, prevalence := round(included_with_variant / num_samples, 3)]
  } else{
    stop('prevalence_method should be "count", "percent", "both", or "auto".')
  }
  
  # legend_size_breaks controls point size and what point sizes are displayed in legend
  if(identical(legend_size_breaks, 'auto')) {
    if(effects[, .N] < 6) {
      size_breaks = sort(unique(effects$prevalence))
    } else {
      ordered_prev = sort(unique(effects$prevalence))
      first_break = ordered_prev[1]
      last_break = max(ordered_prev)
      
      num_middle_breaks_left = min(3, length(setdiff(ordered_prev, c(first_break, last_break))))
      middle_breaks = numeric()
      while(num_middle_breaks_left > 0) {
        biggest_left = ordered_prev[length(ordered_prev)]
        next_biggest = ordered_prev[length(ordered_prev) - 1]
        ideal_spacing = (biggest_left - ordered_prev[1])/(num_middle_breaks_left + 1)
        ideal_next_break = biggest_left - ideal_spacing
        
        if(ideal_next_break > next_biggest) {
          next_break = next_biggest
        } else {
          next_index = which.min(abs(ordered_prev - ideal_next_break))
          next_break = ordered_prev[next_index]
          i = 1
          while(next_break/biggest_left > (1 + num_middle_breaks_left)/5 && 
                next_index > 1) {
            next_break = ordered_prev[next_index - i]
            i = i + 1
          }
        }
        middle_breaks = c(middle_breaks, next_break)
        ordered_prev = c(ordered_prev[ordered_prev < next_break], next_break)
        num_middle_breaks_left = num_middle_breaks_left - 1
        if(next_break == first_break) {
          middle_breaks_left = 0
        }
      }
      size_breaks = unique(sort(c(first_break, middle_breaks, last_break)))
    }
  } else if (identical(legend_size_breaks, FALSE)) {
    size_breaks = 1.5
    effects[, prevalence := as.numeric(prevalence)] # convert from integer to avoid warning
    effects$prevalence = 1.5 # medium-small; we're going to do scale_size_identity()
  } else if(is.numeric(legend_size_breaks)) {
    size_breaks = legend_size_breaks
    if(length(size_breaks) == 1) {
      effects$prevalence = size_breaks # for scale_size_identity()
    }
  } else {
    msg = paste0('legend_size_breaks should be "auto", numeric vector of prevalences to depict, or FALSE to make all points small, ',
                 'or a single numeric point size.')
    stop(pretty_message(msg, emit = F))
  }
  
  if(identical(legend_size_name, 'auto')) {
    if(prevalence_method == 'count'){
      legend_size_name = 'Variant prevalence'
    } else if (prevalence_method == 'percent') {
      legend_size_name = 'Variant frequency\n(within covering samples, \nper effect inference)'
    } else if (prevalence_method == 'both') {
      legend_size_name = 'Variant prevalence\n(percent of samples)'
    }
  } else if(! is.character(legend_size_name) || length(legend_size_name) != 1) {
    stop('legend_size_name should be 1-length character.')
  }
  
  # Handle variant (or variant group) labels
  if(identical(y_label, 'auto')) {
    if(group_by == 'variant') {
      if(uniqueN(effects$variant_name) == uniqueN(effects$variant_id)) {
        effects[, variant_group_label := gsub('_', ' ', variant_name)]
      } else {
        effects[, variant_group_label := variant_id] # unusual situation
      }
    } else {
      effects[, variant_group_label := variant_group]
    }
  } else {
    if(! y_label %in% names(effects)) {
      stop('Column ', y_label, ' does not exist in effects table.')
    }
    setnames(effects, y_label, 'variant_group_label')
  }
  
  # Validate color specification
  # Other nice choices: "darkseagreen3", "lightskyblue4"
  if(! is.character(color_by)) {
    stop('color_by should be type character')
  }
  if(length(color_by) == 1) {
    if(color_by %in% names(effects)) {
      if(color_by %in% colors()) {
        msg = paste0("This is kind of silly: You chose a value for ",
                     "color_by that is both a column name and an R color.")
        stop(pretty_message(msg, emit = F))
      }
      setnames(effects, color_by, 'point_fill')
      
      use_fill_identity = FALSE
      if(is.character(effects$point_fill)) {
        effects[, is_color := point_fill %in% colors()]
        
        # Help out user when only missing colors are NAs.
        if(effects[is_color == F, .N > 0 && all(is.na(point_fill))]) {
          stop("It seems like your color_by column is giving color names, but some rows have NA values.")
        }
        use_fill_identity = all(effects$is_color)
      }

    } else {
      effects$point_fill = color_by
      use_fill_identity = TRUE
    }
  } else {
    stop("color_by should be an R color name (\"purple4\") or the name of a column in effects")
  }

  variant_groups = uniqueN(effects$variant_group)
  group_labels = uniqueN(effects$variant_group_label)
  if(length(variant_groups) != length(group_labels)) {
    if(group_by == 'variant') {
      stop('There is not exactly one unique label per variant. Check your y_label.')
    } else {
      stop("There is not exactly one unique label per group of variants. Check y_label.")
    }
  }
  
  # When there are multiple variants with same y-position (that is, multiple variants in a variant
  # group), nudge y-position so that points/CIs don't overlap.
  effects[, y_nudge := scale((1:.N)/.N, center = T, scale = F), by = 'variant_group']
  effects[y_nudge != 0, y_nudge := (y_nudge * .15) / max(y_nudge), by = 'variant_group']
  
  x_limits = c(min(effects$ci_low_95, na.rm = T), max(effects$ci_high_95, na.rm = T))
  
  # Alternating rows of output will have darker and lighter dashed lines to connect to group names.
  # The dashed lines go until highest lower CI in each group.
  effects[, line_color := ifelse(.GRP %% 2 == 0, 'gray60', 'gray90'), by = 'variant_group']
  
  # When just one variant per row, dashed lines go to lower CI. With multiple variants, we'll 
  # do the dashed line all the way across
  if(identical(as.integer(effects[, .N, by = 'variant_group'][, unique(N)]), 1L)) {
    effects[, dash_end := ci_low_95]
  } else {
    effects[, dash_end := Inf] # to end of visible plot
  }
  
  # For aesthetics, we'll eliminate the CI crossbars for groups that have lots of variants.
  effects[, ci_width := ifelse(.N > 4, 0, .2), by = 'variant_group']
  
  
  x_labeler = function(x) {
    first_visible_label = which(! is.na(x))[1]
    x = format(x, scientific = F, big.mark = ',')
    if(! is.null(lowest_label)) {
      x[first_visible_label] = lowest_label 
    }
    return(x)
  }
  
  gg = ggplot(effects, aes(x = selection_intensity, y = variant_group)) +
    geom_segment(aes(x = x_limits[1], xend = dash_end, y = variant_group, yend = variant_group), 
                 color = effects$line_color, linetype = 'dotted', na.rm = T) +
    geom_errorbar(aes(xmin = ci_low_95, xmax = ci_high_95), color = "azure4", na.rm = T, 
                  position = position_nudge(x = 0, y = effects$y_nudge), width = effects$ci_width, linewidth = .25) +
    geom_point(shape = 21, color = 'gray20', aes(size = prevalence, fill = point_fill), na.rm = T,
               position = position_nudge(x = 0, y = effects$y_nudge)) +
    scale_x_log10(expand = expansion(mult = c(.01, .05)), labels = x_labeler) + 
    scale_y_discrete(limits = unique(effects$variant_group), labels = unique(effects$variant_group_label),
                     expand = expansion(add = 1))  +
    labs(title = title, x = x_axis_title, y = y_axis_title) +
    theme(axis.title.x = element_text(margin = margin(6, 0, 6, 0)),
          axis.title.y = element_text(margin = margin(0, 6, 0, 6)),
          axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 8),
          axis.text.x = element_text(size = 8),
          axis.ticks.x = element_line(color = 'gray50'),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = legend.position,
          legend.direction = 'vertical',
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          plot.margin = margin(l = 6, r = 15, b = 6, unit = 'pt'),
          plot.title = element_text(margin = margin(t = 6, b = 6)))
  
  # Put a border around the legend if it's within the plot space
  if(is.numeric(legend.position) && all(legend.position > 0) && all(legend.position < 1)) {
    gg = gg + theme(legend.background = element_rect(fill = alpha(c("white"), 0.9), linewidth = .2, color = 'gray20'))
  }
  
  if(length(size_breaks) > 1) {
    # If there is just one point fill color in the plot, make legend's size glyphs that color.
    size_override = list()
    if(use_fill_identity && uniqueN(effects$point_fill) == 1) {
      only_color = effects$point_fill[1]
      size_override = list(fill = only_color, alpha = 1)
    }
    size_labels = size_breaks
    if(prevalence_method == 'percent') {
      size_labels = scales::label_percent(accuracy = .1)(size_labels)
    }
    if(prevalence_method == 'both') {
      # Have already verified that all effects have same number of samples
      size_labels = paste0(size_breaks, ' (', scales::label_percent(accuracy = .1)(size_breaks/effects$num_samples[1]), ')')
    }
    size_range = c(1, 6)
    if (length(size_breaks) < 3) {
      size_range = c(1, length(size_breaks))
    }
    gg = gg + scale_size(breaks = size_breaks, labels = size_labels, 
                         limits = c(min(effects$prevalence), max(effects$prevalence)),
                         guide = guide_legend(title.position = 'top', override.aes = size_override),
                         name = legend_size_name, range = size_range)
  } else {
    gg = gg + scale_size_identity()
  }
  
  
  # Validate label_individual_variants and decide whether individual labels are happening.
  # Only worth labeling variants if there is more than one variant per variant_group/color grouping
  if(identical(label_individual_variants, TRUE) && 
     effects[, .N, by = c('variant_group', 'point_fill')][, all(N == 1)]) {
    label_individual_variants = FALSE
  }
  
  if(identical(label_individual_variants, TRUE)) {
    # Use variant_name (unless variant_id is necessary due to ambiguity) if nothing supplied,
    # unless grouping by gene, in which get aachange from gene name.
    if(group_by == 'gene') {
      effects[variant_type == 'aac', individual_label := gsub('.*_', '', variant_name)]
      effects[variant_type != 'aac', individual_label := gsub('_', ' ', variant_name)]
      
    } else if (uniqueN(effects$variant_id) != uniqueN(effects$variant_name)) {
      effects[, individual_label := variant_id]
    } else {
      effects[, individual_label := gsub('_', ' ', variant_name)]
    }
    
  } else if(is.character(label_individual_variants) && length(label_individual_variants) == 1 && 
            label_individual_variants %in% names(effects)) {
    setnames(effects, label_individual_variants, 'individual_label')
    if(! is.character(effects$individual_label)) {
      msg = paste0('Column specified for label_individual_variants (', label_individual_variants, ') is not type character.')
      stop(pretty_message(msg, emit = F))
    }
    label_individual_variants = TRUE
  } else if (! identical(label_individual_variants, FALSE)){
      stop('label_individual_variants should be TRUE/FALSE or the name of a column in the effects table.')
  }
    
  if(label_individual_variants) {
    gg = gg + geom_label_repel(aes(label = individual_label), size = 2.5, box.padding = .3, label.r = .2,
                               fill = alpha(c("white"), 0.9), label.size = .1, label.padding = .15,
                               segment.color = 'grey20', segment.size = .4,
                               position = position_nudge(x = 0, y = effects$y_nudge))
  }
  
  # Change axis label using group_by when not "variant" (unless user already explicitly specified
  # via y_title).
  if(group_by != 'variant' && is.null(y_title)) {
    if (group_by == 'gene') {
      gg = gg + labs(y = 'Gene') # auto-capitalize
    } else {
      gg = gg + labs(y = group_by)
    }
  }
  
  # Validate legend_color_name
  if(is.null(legend_color_name)) {
    legend_color_name = color_by
  } else if(! is.character(legend_color_name) || length(legend_color_name) != 1) {
    stop('legend_color_name should be NULL or 1-length character.')
  }
  
  # Handle legend labels for colors.
  if(! is.null(color_label)) {
    if(! is.character(color_label)) {
      stop('color_label should be type character.')
    }
    if(! color_label %in% names(effects)) {
      stop("color_label column ", color_label, ' not found in effects table.')
    }
    setnames(effects, color_label, 'fill_label')
    if(! is.character(effects$fill_label)) {
      stop("color_label column ", color_label, ' is not type character.')
    }
    if(uniqueN(effects[, .(point_fill, fill_label)]) != uniqueN(effects$fill_label)) {
      stop("There is not a one-to-one correspondence between color names in color_by and labels in color_label.")
    }
  }
  
  if(use_fill_identity) {
    if(is.null(color_label)) {
      gg = gg + scale_fill_identity()
    } else {
      gg = gg + scale_fill_identity(breaks = unique(effects$point_fill), labels = unique(effects$fill_label), 
                                    guide = guide_legend(), name = legend_color_name)
    }
  } else {
    if(is.null(color_label)) {
      if(is.numeric(effects$point_fill)) {
        gg  = gg + scale_fill_viridis_c(name = legend_color_name, option = viridis_option)
      } else {
        gg  = gg + scale_fill_viridis_d(name = legend_color_name, begin = .2, end = .9, option = viridis_option)
      }
    } else {
      unique_colors = unique(effects[, .(point_fill, fill_label)])
      labels_by_color = setNames(unique_colors$fill_label, unique_colors$point_fill)
      if(is.numeric(effects$point_fill)) {
        gg = gg + scale_fill_viridis_c(name = legend_color_name, breaks = labels_by_color, labels = labels_by_color)
      } else {
        gg = gg + scale_fill_viridis_d(name = legend_color_name, labels = labels_by_color)
      }
    }
  }
  gg
}
