#' Plot pairwise epistasis
#'
#' Visualize pairwise epistatic scaled selection coefficients, as calculated by
#' \code{ces_epistasis()} or \code{ces_gene_epistasis()}. For each variant, the isolated site
#' effect--the overall scaled selection at the site without regard for the mutation status of the
#' other site--is also depicted.
#' 
#' Variant pairs for which the epistatic selection model is not significantly better than a
#' no-epistasis model (in which \code{ces_A0 = ces_A_on_B} and \code{ces_B0 = ces_B_on_A}) are depicted
#' with faded arrows; the threshold is \code{p_epistasis < .05} (or, instead of p_epistasis,
#' whatever corrected significance column is specified with \code{significance_cols}).
#' 
#' @param epistatic_effects Epistatic effects inference table, as produced by \code{ces_epistasis()}
#'   or \code{ces_gene_epistasis()}.
#' @param pairs_per_row How many epistatic pairs to show in each plot row. The provided value is
#'   incremented if need to prevent the legend schematic from being isolated on its own row.
#' @param x_title X-axis label. Set NULL for no label.
#' @param variant_label_size Text size for the variant labels.
#' @param dodge_labels TRUE/FALSE on using n.dodge (height staggering) on variant labels. Defaults
#' to "auto"; you can try setting manually if labels are not looking good.
#' @param alternating_colors Colors, provided as character vector, to use on epistatic effect
#'   arrows. It's recommended to supply one or two colors, but more will work.
#' @param include_schematic TRUE/FALSE on whether to include the schematic that shows how to interpret the plot.
#' @param schematic_label_size Text size of labels in the schematic.
#' @param significance_levels A vector of 1-3 distinct numeric values on (0, 1) in descending order
#'   to use for significance annotations.
#' @param significance_cols A named list of column names that give significance values for nonzero
#'   change in selection for each pair of sites A and B, and for the performance of the epistatic
#'   selection model over a model that ignores epistatic interactions. List elements must be named
#'   A_change (default "p_A_change") B_change ("p_B_change"), and model ("p_epistasis"). The purpose
#'   of this option is to support the use of different columns when multiple testing correction is
#'   performed.
#' @param inference_floor Numeric value of the optimization floor used in epistatic effect
#'   inference. Typically, should be left at the default value, which matches cancereffectsizeR
#'   epistatic inference defaults. For plot legibility, there will be a dashed horizontal line in
#'   the output plot, higher than \code{inference_floor} and lower than any non-minimized parameter
#'   inference, which indicates that all arrows pointing below the line have estimates at the
#'   optimization floor.
#' @export
plot_epistasis = function(epistatic_effects, pairs_per_row = 8,
                          x_title = 'Site pairs',
                          variant_label_size = 6.5, dodge_labels = 'auto', alternating_colors = c("#7cb4de", "#7DD3AF"),
                          include_schematic = TRUE,
                          schematic_label_size = 2.5, 
                          significance_levels = c(.05, .01, .001), 
                          significance_cols = list(A_change = 'p_A_change', 
                                                   B_change = 'p_B_change', 
                                                   model = 'p_epistasis'), 
                          inference_floor = .001) {
  # old colors: c("#8dbee2", "#a3d8c2")
  if(! require('ggplot2')) {
    stop('ggplot2 must be installed.')
  }
  if(packageVersion('ggplot2') < as.package_version('3.5.0')) {
    stop('Update ggplot2: it must be version 3.5.0 or later.')
  }
  if(! require('ggrepel')) {
    stop('ggrepel must be installed.')
  }
  
  if(! identical(dodge_labels, 'auto')) {
    if(! rlang::is_bool(dodge_labels)) {
      stop('dodge_labels should be TRUE/FALSE, or \"auto\".')
    }
  }
  
  if(! rlang::is_bool(include_schematic)) {
    stop('include_schematic should be TRUE/FALSE')
  }
  if(! rlang::is_scalar_integerish(pairs_per_row) || pairs_per_row < 1) {
    stop('pairs_per_row should be a positive integer.')
  }
  
  if(! is.numeric(inference_floor) || length(inference_floor) != 1 || inference_floor <=0) {
    stop('inference_floor should be numeric and greater than 0.')
  }
  
  remove_x_title = is.null(x_title)
  if(is.null(x_title) || (is.na(x_title) && length(x_title) == 1)) {
    x_title = ''
  }
  if(! rlang::is_scalar_character(x_title)) {
    stop('x_title should be character (or NULL for no axis label).')
  }
  
  if(! is.numeric(variant_label_size) || length(variant_label_size) != 1 || variant_label_size <= 0) {
    stop('variant_label_size should be a positive numeric value.')
  }
  
  if(! is.numeric(schematic_label_size) || length(schematic_label_size) != 1 || schematic_label_size <= 0) {
    stop('schematic_label_size should be a positive numeric value.')
  }
  
  if(! is.character(alternating_colors) || length(alternating_colors) < 1) {
    stop('alternating_colors should be a character vector that specifies colors.')
  }
  
  if(! is.list(significance_cols) || length(significance_cols) != 3 || 
     ! all(c('A_change', 'B_change', 'model') %in% names(significance_cols))) {
    stop('significance_cols must be named list of length three (see docs).')
  }
  
  
  required_cols = c("ces_A0", "ces_B0", "ces_A_on_B", "ces_B_on_A", "ci_low_95_ces_A0", 
                    "ci_high_95_ces_A0", "ci_low_95_ces_B0", "ci_high_95_ces_B0", 
                    "ci_low_95_ces_A_on_B", "ci_high_95_ces_A_on_B", "ci_low_95_ces_B_on_A", 
                    "ci_high_95_ces_B_on_A", "ces_A_null", "ces_B_null",
                    "variant_A", "variant_B")
  if(! all(sapply(significance_cols, rlang::is_scalar_character))) {
    stop('All elements of the significance_cols list must be 1-length character vectors.')
  }
  if(is.null(significance_levels)) {
    required_cols = c(required_cols, significance_cols[['model']])
  } else {
    required_cols = c(required_cols, unlist(significance_cols))
  }
  
  effects = copy(epistatic_effects)
  if(! is.data.table(effects)) {
    stop('epistatic_effects should be a data.table (from ces_epistasis()/ces_gene_epistasis()).')
  }
  
  missing_cols = setdiff(required_cols, names(effects))
  if(length(missing_cols) > 0) {
    stop("Required columns missing from epistatic_effects: ", "\n", paste0(missing_cols, collapse = ', '), '.')
  }
  if(length(names(effects)[names(effects) %in% required_cols]) != length(required_cols)) {
    stop('epistatic_effects table has duplicate column names.')
  }
  
  effects$model_significance = effects[[significance_cols[['model']]]]
  effects[, is_signif := model_significance < .05]
  
  effects[, change_in_A := ces_A_on_B - ces_A0] 
  effects[change_in_A < 0, A_arrow := 'triangle down filled']
  effects[change_in_A >= 0, A_arrow := 'triangle filled']
  
  effects[, change_in_B := ces_B_on_A - ces_B0]
  effects[change_in_B < 0, B_arrow := 'triangle down filled']
  effects[change_in_B >= 0, B_arrow := 'triangle filled']
  
  # If less than 6 per row, keep adding more per row until legend is not overhanging
  legend_width = floor(pairs_per_row / 4) + 1
  n_grp = ceiling(effects[, .N]/pairs_per_row) 
  
  while(TRUE && include_schematic == TRUE) {
    if(pairs_per_row > effects[, .N]) {
      spaces_open = pairs_per_row - effects[, .N]
    } else {
      spaces_open = pairs_per_row - (effects[, .N] - (n_grp - 1) * pairs_per_row)
    }
    # Legend must have enough space and also not be on its own row
    if(spaces_open >= legend_width && spaces_open != pairs_per_row) break
    pairs_per_row = pairs_per_row + 1
    n_grp = ceiling(effects[, .N]/pairs_per_row) 
    legend_width = floor(pairs_per_row / 4) + 1
  }
  
  # Avoid empty space if someone sets pairs_per_row = 10 with only 4 pairs, for example
  while(pairs_per_row > effects[, .N] + legend_width) {
    pairs_per_row = pairs_per_row - 1
    legend_width = max(1, legend_width - 1)
  } 
  
  effects[, grp := ceiling(.I / pairs_per_row)]
  
  effects[, x := 1:.N]
  effects[, v1_x := x - .20]
  effects[, v2_x := x + .20]
  
  low_values = na.omit(effects[, unlist(.(ces_A0, ces_B0, ces_A_on_B, ces_B_on_A,
                                          ci_low_95_ces_A0, ci_low_95_ces_B0, ci_low_95_ces_A_on_B, 
                                          ci_low_95_ces_B_on_A))])
  data_min = min(low_values[low_values > inference_floor])
  
  threshold_line = 10^floor(log10(data_min))
  plot_ymin = threshold_line * .6
  #threshold_line = sqrt(plot_ymin * data_min) # log-space mean of plot_ymin and data_min
  
  # Enforce threshold to display, so all of the selection-to-zero arrows are cleaner
  cols_to_set = c('ci_low_95_ces_A0', 'ci_low_95_ces_A_on_B', 'ci_low_95_ces_B0', 'ci_low_95_ces_B_on_A',
                  'ces_A0', 'ces_B0', 'ces_A_on_B', 'ces_B_on_A')
  for(x in cols_to_set) {
    effects[[x]][is.na(effects[[x]])] = threshold_line
    effects[[x]][effects[[x]] < data_min] = threshold_line
  }
  setDT(effects) # avoid a harmless data.table warning caused by the above loop
  
  
  effects[, var1_signif := ""]
  effects[, var2_signif := ""]
  effects[, change_A_signif_y := ces_A_on_B]
  effects[, change_B_signif_y := ces_B_on_A]
  
  if(! is.null(significance_levels)) {
    if(! is.numeric(significance_levels) || length(significance_levels) < 1 || length(significance_levels) > 3 ||
       any(significance_levels <= 0 | significance_levels > 1) ||
       uniqueN(significance_levels) != length(significance_levels) || 
       ! identical(sort(significance_levels, decreasing = T), significance_levels)) {
      stop('significance_levels should be a vector of 1-3 distinct values in descending order (or NULL).')
    }
    effects[p_A_change < significance_levels[1], var1_signif := "*"]
    effects[p_B_change < significance_levels[1], var2_signif := "*"]
    
    if(length(significance_levels) > 1) {
      effects[p_A_change < significance_levels[2], var1_signif := "**"]
      effects[p_B_change < significance_levels[2], var2_signif := "**"]
    }
    if(length(significance_levels) > 2) {
      effects[p_A_change < significance_levels[3], var1_signif := "***"]
      effects[p_B_change < significance_levels[3], var2_signif := "***"]
    }
  }
  
  
  
  x_label_pos = unlist(S4Vectors::zipup(effects$v1_x, effects$v2_x))
  x_labels = unlist(S4Vectors::zipup(effects$variant_A, effects$variant_B))
  
  effects[, ymin := plot_ymin]
  effects[, ymax := max(ci_high_95_ces_A0, ci_high_95_ces_B0, ci_high_95_ces_A_on_B, ci_high_95_ces_B_on_A,
                        change_A_signif_y, change_B_signif_y, na.rm = T) * 2] # factor of 2 for breathing room on top
  effects[, bg_rect_xmin := v1_x - .27]
  effects[, bg_rect_xmax := v2_x + .27]
  
  if(effects[grp == 1, .N %% 2] == 0) {
    effects[grp %% 2 == 1, bg_color := rep_len(c('gray96', 'gray90'), .N)]
    effects[grp %% 2 == 0, bg_color := rep_len(c('gray90', 'gray96'), .N)]
    effects[grp %% 2 == 1, alt_color := rep_len(alternating_colors, .N)]
    effects[grp %% 2 == 0, alt_color := rep_len(rev(alternating_colors), .N)]
    
  } else {
    effects[, bg_color := rep_len(c('gray96', 'gray90'), .N)]
    effects[, alt_color := rep_len(alternating_colors, .N)]
  }
  
  main_line_width = 2
  row_boundaries = effects[, .(c(min(bg_rect_xmin), max(bg_rect_xmax))), by = 'grp'][, unlist(V1)]
  names(row_boundaries) = rep(' ', length(row_boundaries))
  names(x_label_pos) = scales::label_wrap(19)(x_labels)
  final_pos = sort(c(row_boundaries, x_label_pos))
  
  grouping_table = effects[, .(x = c(min(bg_rect_xmin), max(bg_rect_xmax)), y = threshold_line), by = 'grp']
  
  row_width = grouping_table$x[2] - grouping_table$x[1]
  grouping_table[.N, x := grouping_table[.N - 1, x] + row_width]
  
  if(identical(dodge_labels, 'auto')) {
    n_dodge = ifelse(pairs_per_row > 9, 2, 1)
  } else {
    n_dodge = ifelse(dodge_labels, 2, 1)
  }
  
  
  draw_signif_key = function(data, params, size) {
    grid::grobTree(grid::linesGrob(x = c(.1, .9), y = c(.35, .35), arrow = arrow(type = 'closed', length = unit(.15, 'npc')),
                                   gp = grid::gpar(col = alpha(alternating_colors[1], data$alpha), lwd = 2)),
                   grid::linesGrob(x = c(.1, .9), y = c(.65, .65), arrow = arrow(type = 'closed', length = unit(.15, 'npc')),
                                   gp = grid::gpar(col = alpha(alternating_colors[2], data$alpha), lwd = 2)))
  }
  
  y_labeler = function(x) {
    x = format(x, big.mark = ",", scientific = FALSE)
    which_threshold = suppressWarnings(which(as.numeric(x) == threshold_line))
    x[which_threshold] = format(paste0('<', as.numeric(x[which_threshold])), big.mark = ',', scientific = FALSE)
    return(x)
  }
  gg = ggplot(data = effects) + 
    
    # Alternating grey background boxes
    geom_rect(aes(xmin = bg_rect_xmin, xmax = bg_rect_xmax, ymin = plot_ymin, ymax = ymax), 
              fill = effects$bg_color) + 
    
    # Dashed line indicating estimates/CIs that go to floor
    geom_segment(aes(x = bg_rect_xmin, xend = bg_rect_xmax, y = threshold_line, yend = threshold_line), color = 'darkgrey', linetype = 'dashed') +
    
    # Plot v1 arrows and errors bars
    geom_segment(aes(x = v1_x - .1, xend = v1_x - .1, y = ci_low_95_ces_A0, yend = ci_high_95_ces_A0), color = "grey62") +
    geom_segment(aes(x = v1_x + .1, xend = v1_x + .1, y = ci_low_95_ces_A_on_B, yend = ci_high_95_ces_A_on_B), color = "grey62") +
    
    # Need to put background-colored circle/triangle behind each circle/triangle to
    # avoid alpha values stacking between those objects and the geom_segment on insignificant pairs.
    geom_segment(aes(x = v1_x, xend = v1_x, y = ces_A0, yend = ces_A_on_B, color = alt_color, alpha = is_signif), 
                 key_glyph = draw_signif_key, linewidth = main_line_width, lineend = "butt") +
    geom_point(aes(x = v1_x, y = ces_A0, color = effects$bg_color), size = 2.5) +
    geom_point(aes(x = v1_x, y = ces_A0, color = alt_color, alpha = is_signif), size = 2.5, show.legend = FALSE) +
    geom_point(aes(x = v1_x, y = ces_A_on_B, shape = A_arrow, fill = effects$bg_color, color = effects$bg_color), size = 2.5) +
    geom_point(aes(x = v1_x, y = ces_A_on_B, shape = A_arrow, fill = alt_color, color = alt_color, alpha = is_signif), 
               size = 2.5, show.legend = F) +  
    
    geom_point(aes(x = v1_x, y = ces_A_null), color = 'darkslateblue', shape = 'square') +
    geom_text(data = effects[var1_signif != ''], aes(x = v1_x - .06, y = change_A_signif_y, label = var1_signif, hjust = 'right', vjust = 'inward'), size = 5) + 
    
    # Plot v2 arrows and error bars
    geom_segment(aes(x = v2_x - .1, xend = v2_x - .1, y = ci_low_95_ces_B0, yend = ci_high_95_ces_B0), color = "grey62") + 
    geom_segment(aes(x = v2_x + .1, xend = v2_x + .1, y = ci_low_95_ces_B_on_A, yend = ci_high_95_ces_B_on_A), color = "grey62") +
    
    geom_segment(aes(x = v2_x, xend = v2_x, y = ces_B0, yend = ces_B_on_A, color = alt_color, alpha = is_signif), 
                 key_glyph = draw_signif_key, linewidth = main_line_width, lineend = "butt") +
    geom_point(aes(x = v2_x, y = ces_B0, color = effects$bg_color), size = 2.5) +
    geom_point(aes(x = v2_x, y = ces_B0, color = alt_color, alpha = is_signif), size = 2.5, show.legend = FALSE) + 
    geom_point(aes(x = v2_x, y = ces_B_on_A, shape = B_arrow, color = effects$bg_color, fill = effects$bg_color), size = 2.5) +
    geom_point(aes(x = v2_x, y = ces_B_on_A, shape = B_arrow, fill = alt_color, color = alt_color, alpha = is_signif), 
               size = 2.5, show.legend = FALSE) +              
    
    geom_point(aes(x = v2_x, y = ces_B_null), color = 'darkslateblue', shape = 'square') +
    geom_text(data = effects[var2_signif != ''], aes(x = v2_x - .06, y = change_B_signif_y, label = var2_signif, hjust = 'right', vjust = 'inward'), size = 5) + 
    geom_point(data = grouping_table, aes(x = x, y = y), alpha = 0) +
    scale_fill_identity(guide = guide_none(), aesthetics = c('color', 'fill')) +
    scale_shape_identity() +
    scale_alpha_manual(breaks = c(T, F), values = c(1, .25), labels = c('Significant', 'Not significant'),
                       guide = guide_legend(title = 'Pairwise epistasis', 
                                            override.aes = list(color = alternating_colors[1]))) + 
    scale_x_continuous(breaks = final_pos, labels = names(final_pos), expand = expansion(0),
                       guide = guide_axis(cap = T, n.dodge = n_dodge)) + 
    scale_y_log10(labels = y_labeler, limits = c(plot_ymin, NA), 
                  expand = expansion(0)) +
    xlab(x_title) + ylab('Scaled selection coefficients') + 
    theme_classic() + 
    theme(axis.ticks.length.x = unit(0, 'cm'), 
          axis.text.x = element_text(size= variant_label_size, margin = margin(t = 4, unit = 'pt')),
          plot.margin = margin(t = 12, r = 12, b= 24, l = 12)) + 
    
    facet_wrap(~grp, nrow = n_grp, scales = 'free') + theme(strip.text = element_blank(), panel.grid = element_blank())
  
  ymax = effects$ymax[1]
  
  if(include_schematic) {
    if(legend_width == 1) {
      legend_x = n_grp * pairs_per_row - .25
      legend_box_xmin = n_grp * pairs_per_row - .47
    } else if(legend_width == 2) {
      if(n_grp < 4) {
        legend_x = n_grp * pairs_per_row - 1
      } else {
        legend_x = n_grp * pairs_per_row - 1.25
      }
      legend_box_xmin = n_grp * pairs_per_row - 1.47 
    } else if(legend_width == 3) {
      if(n_grp < 4) {
        legend_x = n_grp * pairs_per_row - 1.75
      } else {
        legend_x = n_grp * pairs_per_row - 2
      }
      legend_box_xmin = n_grp * pairs_per_row - 2.47 
    } else {
      legend_x = n_grp * pairs_per_row - (legend_width + 1) * .5
      legend_box_xmin = n_grp * pairs_per_row - (legend_width - 1) - .47
    }
    legend_box_xmax = n_grp * pairs_per_row + .47
    legend_title_x = mean(c(legend_box_xmin, legend_box_xmax))
    
    log_range = log10(ymax) - log10(plot_ymin)
    legend_y_low = 10^(log10(plot_ymin) + .12 * log_range) # low end 25% from bottom of box
    legend_y_high = 10^(log10(plot_ymin) + .7 * log_range)
    legend_null_y = 10^(log10(plot_ymin) + .3 * log_range)
    legend_title_y = plot_ymin + 10^(log10(plot_ymin) + .9 * log_range)
    legend_data = data.table(grp = n_grp)
    legend_text_spacing = ifelse(legend_width == 1 || n_grp < 4, '\n', ' ')
    gg = gg + 
      geom_text(data = legend_data, aes(x = legend_title_x, y = legend_title_y, vjust = 'inward', hjust = 'center', 
                                        label = 'Types of effects'), size = schematic_label_size + 1, lineheight = .75) +
      geom_point(data = legend_data, aes(x = legend_x, y = legend_y_low), size = 2.5, color = 'goldenrod2') +
      geom_segment(data = legend_data, aes(x = legend_x, xend = legend_x, y = legend_y_low, yend = legend_y_high), color = 'goldenrod2', linewidth = 2, lineend = "butt") +
      geom_point(data = legend_data, aes(x = legend_x, y = legend_y_high), shape = 'triangle filled', size = 2.5, fill = 'goldenrod2', color = 'goldenrod2') +
      geom_text_repel(data = legend_data, aes(x = legend_x, y = legend_y_low, label = paste0('Epistatic effect:', legend_text_spacing, 'paired site wildtype')),
                      size = schematic_label_size, lineheight = .9, force = 0, min.segment.length = 0, point.padding = unit(.7 , 'lines'),
                      hjust = 'left', vjust = 'inward', nudge_x = legend_width * .15) +
      geom_text_repel(data = legend_data, aes(x = legend_x, y = legend_y_high, label = paste0('Epistatic effect:', legend_text_spacing, 'paired site mutated')),
                      size = schematic_label_size, lineheight = .9, force = 0, min.segment.length = 0, point.padding = unit(.7 , 'lines'),
                      hjust = 'left', vjust = 'inward', nudge_x = legend_width * .15) +
      geom_point(data= legend_data, aes(x = legend_x, y = legend_null_y), color = 'darkslateblue', shape = 'square' ) +
      geom_text_repel(data = legend_data, aes(x = legend_x, y = legend_null_y, label = paste0('\n\nIsolated site effect:', legend_text_spacing, 'paired site ignored\n')),
                      size = schematic_label_size, force = 0, lineheight = .9, direction = 'y', min.segment.length = 0, point.padding = unit(.7 , 'lines'),
                      hjust = 'left', vjust = 'inward', nudge_x = legend_width * .15) + 
      geom_rect(data = legend_data, aes(xmin = legend_box_xmin,
                                        xmax = legend_box_xmax,
                                        ymin = plot_ymin, ymax = ymax), fill = NA, color = 'black')
  }
  
  if(remove_x_title) {
    gg = gg + theme(axis.title.x = element_blank())
  }
  
  if(all(effects$is_signif)) {
    gg = gg + theme(legend.position = 'none')
  } else {
    gg = gg + theme(legend.direction = 'horizontal', 
                    legend.title = element_text(size = 10, lineheight = .5, hjust = 0, vjust = .6, margin = margin(b = 0, r = 12)),
                    legend.title.position = 'left',
                    legend.position = 'inside', legend.position.inside = c(0, 0), 
                    legend.justification = 'left', legend.margin = margin(t = 55 + 12 * (n_dodge > 1), r = 6),
                    legend.background = element_rect(fill = rgb(0, 0, 0, 0)))
  }
  return(list(gg, effects))
}

