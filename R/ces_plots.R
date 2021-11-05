#' Plot selection intensities
#' 
#' Lollipop plots can be used to compare highly selected variants
#' within or among groups of samples.
#' 
#' To keep the lollipops readable, no more than \code{max_sites} variants are shown for
#' any group. At the same time, all SIs that are within the limits of the output plot are
#' depicted (which means that specification of \code{max_sites} creates an implicit SI
#' floor). If the labels don't look good, you may have to reduce max_sites or break 
#' up the plot into pieces. (If you find a better way to get lots of variant labels to
#' display nicely, we would love to hear about it!)
#' 
#' 
#' @param si_list Selection output table or named list of selection output tables. Each
#'   table should have columns called "variant_name", "variant_type", and "gene", plus a
#'   selection intensity column (or columns, and then each column will have its values
#'   plotted on a separate lollipop.)
#' @param group_names Names to use for labeling each group of selection intensities.
#' @param si_col Vector giving names of selection intensity column in each SI table,
#'               either 1-length or same length as si_list.
#' @param max_sites Maximum number of variant sites to include per lollipop; if you try to
#'   include too many, you may have a challenge getting it to look good.
#' @param title Plot title to pass to ggplot.
#' @param ylab Y-axis label to pass to ggplot.
#' @param label_size Text size for labels, either 1-length or same length as si_list.
#' @param merge_dist How close points must be to be eligible to have their labels combined
#'   (.04 = 4 percent of plot space). Either 1-length or same length as si_list. Try tweaking if
#'   labels are not looking good. Set to zero to prevent any labels being combined.
#' @return ggplot object with lollipops
#' @export
lollipops = function(si_list, title = "My SIs", ylab = "selection intensity", 
                     max_sites = 50, label_size = 3.0, merge_dist = .04, 
                     si_col = "auto", group_names = "auto") {
  if (! require("ggplot2") || ! require("ggrepel")) {
    stop("To use lollipops(), install ggplot2 and ggrepel: install.packages(c('ggplot2', 'ggrepel')).")
  }
  # help out user if they called the function on cesa instead of [cesa]$seleciton
  if (is(si_list, "CESAnalysis")) {
    if(length(si_list@selection_results) > 0) {
      si_list = snv_results(si_list)
    } else {
      stop("This CESAnalysis has nothing in [cesa]$selection.")
    }
  }
  if (! is(si_list, "list") || ! all(sapply(si_list, is.data.table))) {
    if (is(si_list, "data.table")) {
      si_list = list(' ' = si_list) # may not need a display name when plotting one group
    } else {
      stop("si_list should be a list of selection output tables (or a single data table).")
    }
  }

  if(! is(title, "character") | length(title) > 1) {
    stop("title should be 1-length character")
  }
  
  if (! is(group_names, "character")) {
    stop("group_names should be character")
  }
  
  if(! is(max_sites, "numeric")) {
    stop("max_sites should be numeric")
  }
  
  if(! is(si_col, "character")) {
    stop("si_col should be character")
  }
  
  final_si_list = list()
  final_si_col = character()
  poss_auto_names = character() # build up names to use if auto-naming groups
  input_si_names = names(si_list)
  
  if (length(si_col) == 1 && si_col == "auto") {
    for (i in 1:length(si_list)) {
      curr_si = si_list[[i]]
      curr_name = input_si_names[i]
      curr_si_col = attr(curr_si, "si_cols", T)
      if(is.null(curr_si_col)) {
        stop("si_col could not be automatically determined for element ", i, " of si_list.")
      }
      num_in_group = length(curr_si_col)
      for (col in curr_si_col) {
        final_si_col = c(final_si_col, col)
        final_si_list = c(final_si_list, list(curr_si))
        # IF EGFR has two SIs, give labels like egfr_pre, egfr_met (otherwise, just call it egfr)
        if (num_in_group > 1) {
          poss_auto_names = c(poss_auto_names, paste(curr_name, col))
        } else {
          poss_auto_names = c(poss_auto_names, curr_name)
        }
      }
    }
  } else {
    if (length(si_col) != length(si_list)) {
      if (length(si_col) == 1) {
        si_col = as.list(rep(si_col, length(si_list)))
      } else {
        stop("si_col should match length of si_length or be 1-length")
      }
    }
    for (i in 1:length(si_list)) {
      curr_si = si_list[[i]]
      curr_si_col = si_col[[i]]
      if(! is(curr_si_col, "character")) {
        stop("All elements of list si_col must be type character")
      }
      if (any(duplicated(curr_si_col)) || ! all(curr_si_col %in% names(curr_si))) {
        stop("si_col column names not found in corresponding tables (or duplicates)")
      }
      final_si_col = c(final_si_col, curr_si_col)
      final_si_list = c(final_si_list, list(curr_si))
    }
  }
  si_list = final_si_list
  si_col = final_si_col
  
  num_groups = length(si_list)
  if (length(group_names) == 1 & group_names[1] == "auto") {
    if (is.null(input_si_names)) {
      stop("Can't determine group_names automatically because input si_list unnamed.")
    }
    group_names = poss_auto_names
  } else if (length(group_names) != num_groups) {
    stop("group_names length doesn't match number of SI groups to be plotted.")
  }

  
  if (! is(label_size, "numeric") || ! length(label_size) %in% c(1, num_groups) || any(label_size < 0)) {
    stop("Illegal label_size value.")
  }
  if (length(label_size) < num_groups) {
    label_size = rep(label_size, num_groups)
  }
  
  if (! is(merge_dist, "numeric") || ! length(merge_dist) %in% c(1, num_groups) || any(merge_dist < 0) || any(merge_dist > 1)) {
    stop("Illegal merge_dist value.")
  }
  if (length(merge_dist) < num_groups) {
    merge_dist = rep(merge_dist, num_groups)
  }
  
  
  # for each group, get the highest N SIs (set by max_sites), then pick the
  # greatest of all of these as the SI floor so that no more than max_sites SIs print per lollipop
  min_si = 0
  genes_available = TRUE
  all_genes = character()
  for (i in 1:num_groups) {
    curr = copy(si_list[[i]])
    curr_cols = colnames(curr)
    curr_si_col = si_col[i]
    if (! curr_si_col %in% curr_cols) {
      stop("Selection intensity column not found in input table ", i, " (set with si_col).")
    }
    if (! 'variant_id' %in% curr_cols) {
      stop("Missing required column variant_id in input table ", i, " (see help).")
    }
    curr[variant_type == 'aac', variant_name := gsub('_[^_]*$', '', variant_id)]
    curr[variant_type == 'snv', variant_name := variant_id]
    curr[variant_type == 'aac', gene := gsub('_.*', '', variant_name)]

    
    setnames(curr, curr_si_col, "selection_intensity")
    if ("gene" %in% curr_cols) {
      curr = curr[, .(variant_name, variant_type, gene, selection_intensity)]
    } else {
      curr = curr[, .(variant_name, variant_type, selection_intensity)]
      genes_available = FALSE
    }
    
  
    if (curr[, .N] > max_sites) {
      lowest_passing = curr[order(selection_intensity, decreasing = T), selection_intensity[max_sites]]
      if (lowest_passing > min_si) {
        min_si = lowest_passing
      }
    }
    si_list[[i]] = curr
  }
  
  si_list = lapply(si_list, function(x) x[selection_intensity >= min_si])
  
  # variants from genes with multiple variants will get colors in output
  if (genes_available) {
    all_genes = na.omit(unlist(lapply(si_list, function(x) x$gene)))
    multi_variant_genes = names(which(table(all_genes) > 1))
  } else {
    all_variants = na.omit(unlist(lapply(si_list, function(x) x$variant_name)))
    repeated_variants =  names(which(table(all_variants) > 1))
  }



  placeholder = data.table(group = numeric(), selection_intensity = numeric(), variant_name_print = character(),
                           display_color = factor())
  xlim_low = .8
  xlim_high = num_groups + 1
  lolli = ggplot(placeholder, mapping = aes(group, selection_intensity, label = variant_name_print, color = display_color)) + 
          scale_y_log10() + ylab(ylab) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            legend.position = "none",
            plot.title = element_text(size = 20, hjust = .5)
          ) + ggtitle(label = title) +  scale_colour_discrete() + 
          scale_x_continuous(breaks = 1:num_groups, labels = group_names, limits = c(xlim_low, xlim_high))

  for (i in 1:num_groups) {
    top_hits = si_list[[i]][order(selection_intensity, decreasing = T)]
    if (top_hits[, .N] == 0) {
      warning("Group ", i, " has nothing to show in the plot space.", call. = F)
      next
    }
    top_hits[, log_si := log10(selection_intensity)]
    top_hits[, variant_name_print := gsub('_', ' ', variant_name)]
    curr_label_size = label_size[i]
    curr_merge_dist = merge_dist[i]
    
    # leave out gene names if all AACs have same gene
    # otherwise, keep gene names and add them to SNVs as well
    # since intergenic variants have no gene, specify display color
    top_hits[, display_color := NA_character_]
    if (genes_available) {
      if (length(all_genes) == 1) {
        regex = paste0('^', top_hits[variant_type == "aac", gene], '')
        old_names = top_hits[variant_type == "aac", variant_name_print]
        new_names = mapply(function(r, n) gsub(r, '', n), regex, old_names)
        top_hits[variant_type == "aac", variant_name_print := new_names]
      } 
      top_hits[variant_type == "snv" & ! is.na(gene), variant_name_print := paste0('(', gene, ') ', variant_name_print)] 
      top_hits[gene %in% multi_variant_genes, display_color := factor(gene)]
    } else {
      top_hits[, variant_name_print := variant_name]
      top_hits[variant_name %in% repeated_variants, display_color := factor(variant_name)]
    }
  
    top_hits[, group := i]
    
    # Greedily merge labels of tightly clustered points for readability
    log_si_max = max(top_hits$log_si)
    log_si_range = log_si_max - min(top_hits$log_si)
    label_data = top_hits
    num_labels = label_data[, .N]
    wrap_width = (3 / curr_label_size) * 42
    min_wrap_width = wrap_width * .6 
    if (num_labels > 3) {
      label_data[, merge_eligible := T]
      label_data[, screen_y := (log_si_max - log_si) / log_si_range]
      while(TRUE) {
        label_data[, num_char := nchar(variant_name_print)]
        label_data[, num_near := sum(abs(label_data$screen_y - screen_y) < curr_merge_dist & merge_eligible == T), by = "variant_name"]
        l_center = label_data[merge_eligible == T & num_near == max(num_near), which = T][1]
        curr_screen_y = label_data[l_center, screen_y]
        l_to_add = setdiff(label_data[abs(screen_y - curr_screen_y) < curr_merge_dist & merge_eligible == T, which = T], l_center)
        if (length(l_to_add) == 0) {
          break
        }
        
        labels_to_merge = sort(c(l_center, l_to_add))
        labels_text = label_data$variant_name_print[labels_to_merge]
        num_to_add = length(labels_text)
        
        combined_label = paste(labels_text, collapse = ", ")
        combined_label = strwrap(combined_label, width = wrap_width)
        num_lines = length(combined_label)

        if (num_lines > 1 & nchar(combined_label[num_lines]) < min_wrap_width) {
          while(length(combined_label) == num_lines) {
            num_to_add = num_to_add - 1
            labels_to_merge = labels_to_merge[1:num_to_add]
            labels_text = label_data$variant_name_print[labels_to_merge]
            combined_label = paste(labels_text, collapse = ", ")
            combined_label = strwrap(combined_label, width = wrap_width)
          }
        }
        combined_label = paste0(combined_label, collapse = "\n")
        label_data[labels_to_merge[1], variant_name_print := combined_label]
        label_data[labels_to_merge[1], merge_eligible := F]
        
        # When merging labels fails due to unworkable nchar, you can end up with num_to_add == 1
        if(num_to_add > 1) {
          to_drop = labels_to_merge[2:num_to_add]
          label_data = label_data[! to_drop]
        }
      }
    }
    
    # each group i gets SIs plotted on the x = i line, with labels allowed roughly in the (i, i + 1) space
    lolli = lolli + 
            geom_vline(xintercept = i, color = "grey") +
            geom_label_repel(data = label_data, nudge_x = .1, direction = 'both', 
                             segment.size = .3, size = curr_label_size, color = "black",  xlim = c(i, i + .98))
    # silly: scale_color_discrete fails on all NAs despite NA color option
    if (all(is.na(top_hits$display_color))) {
      lolli = lolli + geom_point(data = top_hits, alpha = .7, color = "gray50")
    } else {
      lolli = lolli + geom_point(data = top_hits, alpha = .7)
    }
  }
  return(lolli)
}


