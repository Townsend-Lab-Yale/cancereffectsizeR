#' Plot mutational source and effect attributions
#'
#' Compare the extent to which mutational signatures contribute mutations (mutational source share)
#' to the degree to which they contribute high-effect mutations (cancer effect share).
#' 
#' @param mutational_effects Output from mutational_signature_effects(). To compare
#' groups of samples, supply a named list with each element corresponding to output
#' from a separate run of mutational_signature_effects().
#' @param signature_groupings A data.table of signature names and descriptions; signatures with
#'   identical descriptions are grouped together. Only signatures present in the data get displayed.
#'   Setting to "auto" (the default) uses the table returned by cosmic_sbs_etiologies() which only
#'   makes sense when using COSMIC signatures. A custom table should have columns "name",
#'   "short_name", and "description". Optionally, include a "color" column to manually specify colors
#'   for each group. Alternatively, setting to "cannataro"
#'   applies the same signature grouping and color palette as
#'   \href{https://academic.oup.com/mbe/article/39/5/msac084/6570859}{Cannataro et al. 2022}.
#' @param viridis_option A viridis color mapping, specified with a single letter ('A' to 'H'). By
#'   default, map 'G' (mako) is used.
#' @param num_sig_groups How many groups of signatures to display. Remaining signatures (from groups
#'   with lower effect shares, when averaged across sample groups) get lumped into an "Other
#'   signatures" group.

plot_signature_effects = function(mutational_effects = NULL,
                                               signature_groupings = 'auto',
                                               viridis_option = NULL,
                                               num_sig_groups = 7) {
  if(! require("ggplot2")) {
    stop('Package ggplot2 must be installed for plotting.')
  }
  
  if(! is.numeric(num_sig_groups) || length(num_sig_groups) != 1 || 
     num_sig_groups - as.integer(num_sig_groups) != 0 || num_sig_groups < 1) {
    stop('num_sig_groups must be a positive integer.')
  }
  
  running_cannataro = identical(tolower(signature_groupings), 'cannataro')
  using_cannataro_colors = is.null(viridis_option) && running_cannataro
  
  if(is.null(viridis_option)) {
    viridis_option = 'G'
  } else if(! is.character(viridis_option) || length(viridis_option) != 1 || ! nchar(viridis_option) == 1) {
    stop('Specify viridis color map with a single letter ("A"-"H")')
    # Will leave the specific viridis character unvalidated, in case more maps are created in the future.
  }
  viridis_option = toupper(viridis_option)
  
  # If num_sig_groups is not default, warn that num_sig_groups is ignored under cannataro
  if(num_sig_groups != 7 && running_cannataro) {
    warning('Ignoring num_sig_groups because signature_groupings = "cannataro".')
  }
  
  # assign signature_groupings to table according to signature_groupings parameter
  if (identical(signature_groupings, 'auto')) {
    signature_groupings = cosmic_sbs_etiologies()
  } 
  
  # Check that mutational_effects is mutational_signature_effects() output, or a list of such outputs.
  if (is(mutational_effects, 'list') && length(mutational_effects) == 2 &&
      identical(names(mutational_effects), c("mutational_sources", "effect_shares"))) {
    
    # Convert a single mutational effects output into 1-length list.
    mutational_effects = list(mutational_effects)
    sample_groupings = NA
    
  } else if (is(mutational_effects, 'list') && ! is.null(names(mutational_effects)) &&
             all(sapply(mutational_effects, function(x) length(x) == 2 &&
                        identical(names(x), c("mutational_sources", "effect_shares"))))) {
    # Check that mutational_effects is a named list containing expected outputs.
    # If so, extract names and put into sample_groupings.
    sample_groupings = names(mutational_effects)
    if(uniqueN(setdiff(sample_groupings, '')) != length(sample_groupings)) {
      stop('mutational_effects list entries should have unique names')
    }
  } else {
    stop('mutational_effects must be mutational_signature_effects() output or a named list of such outputs.')
  }
  
  # Create empty list to load data into
  final_df_list = list()
  
  df_weights = rbindlist(lapply(mutational_effects,
                                function(x) {
                                  weights = x$mutational_sources$average_source_shares
                                  data.table(type = 'SW', prop = weights, name = names(weights))
                                }), idcol = 'sample_group')
  df_effects = rbindlist(lapply(mutational_effects,
                                function(x) {
                                  effects = x$effect_shares$average_effect_shares
                                  data.table(type = 'CEW', prop = effects, name = names(effects))
                                }), idcol = 'sample_group')
  df_final = rbind(df_weights, df_effects)
  
  if (running_cannataro) {
    # group signatures into suggested categories as outlined in Cannataro et al.
    signature_groupings = list(
      "Deamination with age, clock-like (1)" = "SBS1",
      "Unknown, clock-like (5)" = "SBS5",
      "APOBEC (2, 13)" = c("SBS2", "SBS13"),
      "Defective homologous recombination (3)" = "SBS3",
      "Tobacco (4, 29)" = c("SBS4", "SBS29"),
      "UV light (7a–d, 38)" = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
      "Prior treatment (11, 31, 32, 35)" = c("SBS11", "SBS31", "SBS32", "SBS35"),
      "Mutagenic chemical exposure (22, 24, 42, 88)" = c("SBS22", "SBS24", "SBS42", "SBS88"),
      "Alcohol-associated (16)" = "SBS16"
    )
    other_label = "Non-actionable and unknown signatures"
    
  } else if(! is.data.table(signature_groupings)) {
    stop('signature_groupings should be "auto", "cannataro", or a data.table describing the signatures (see docs).')
  } else {
    
    if(length(signature_groupings) != uniqueN(names(signature_groupings))) {
      stop('Input signature_groupings table has repeated column names.')
    }
    required_cols = c('name', 'short_name', 'description')
    missing_cols = setdiff(required_cols, names(signature_groupings))
    if(length(missing_cols) > 0) {
      stop('Missing columns in signature_groupings table: ', paste0(missing_cols, collapse = ', '), '.')
    }
    
    if('other signatures' %in% signature_groupings$description) {
      stop('signature_groupings has signatures with description \"other signatures\", which is reserved.')
    }
    other_label = 'other signatures'
    
    all_effect_shares = rbindlist(lapply(mutational_effects, function(x) as.list(x$effect_shares$average_effect_shares)), idcol = 'cohort')
    shares_melted = melt.data.table(all_effect_shares, id.vars = 'cohort', variable.factor = F, variable.name = 'name')
    
    # Subset signature groupings to only include signatures that actually appear
    signature_groupings = signature_groupings[name %in% shares_melted$name]
    
    shares_melted[signature_groupings, c('description', 'short_name') := .(description, short_name), on = 'name']
    
    # Don't consider signatures in other groups (they'll be represented in the "other" group)
    shares_melted = shares_melted[! is.na(short_name)]
    shares_melted[, sig_grp_id := .GRP, by = 'description']
    summed_shares = shares_melted[, .(in_grp_sum = sum(value)), by = c('cohort', 'sig_grp_id')]
    
    top_grp_ids = summed_shares[order(in_grp_sum, decreasing = T)][! duplicated(sig_grp_id)][1:min(.N, num_sig_groups), sig_grp_id]
    
    final_groupings = list()
    for(id in top_grp_ids) {
      grp_info = shares_melted[sig_grp_id == id]
      short_names = paste(unique(grp_info$short_name), collapse = ', ')
      curr_descrip = paste0(grp_info$description[1], ' (', short_names, ')')
      curr_sigs = unique(grp_info$name)
      final_groupings[[curr_descrip]] = curr_sigs
    }
    
    if('color' %in% names(signature_groupings)) {
      if(! is.character(signature_groupings$color)) {
        stop('The optional color column in signature_groupings is expected to be type character.')
      }
      if(uniqueN(signature_groupings[, .(description, color)]) != uniqueN(signature_groupings$description)) {
        stop('Unusable color column in signature_groupings table: Exactly one color must be associated with each signature group.')
      }
      sig_to_color = unique(signature_groupings[, .(name, color)])
      df_final[sig_to_color, color := color, on = 'name']
      df_final[! name %in% unlist(final_groupings), color := NA]
    }
    signature_groupings = final_groupings
  }
  
  # Signatures present in data but not represented in a signature grouping join "other" signatures
  other_signatures = setdiff(df_effects$name, unlist(signature_groupings))
  
  if(all(unlist(signature_groupings) %like% '^SBS') && ! all(other_signatures %like% '^SBS')) {
    msg = paste0('Data contains some non-SBS signatures, but the signature_groupings table only describes SBS signatures. ',
                 'Non-SBS signatures will all be grouped with "', other_label, '". Consider updating ',
                 'the signature_groupings to include these signatures.')
    warning(pretty_message(msg, emit = F))
  }
  signature_groupings[[other_label]] = other_signatures
  
  
  # Unwind list to get a table matching signatures to their labels
  signature_labels = rbindlist(lapply(1:length(signature_groupings), 
                                      function(x) data.table(name = signature_groupings[[x]], 
                                                             label = names(signature_groupings)[x])))
  # create group column and assign signature_groupings
  # nested loop syntax is for cases when signature_groupings contains several signatures
  # Ex: "UV light" = SBS7a-d and SBS38
  df_final[signature_labels, group := label, on = 'name']
  
  # first, we reorder and rename some values to make the graph pretty.
  # reorder facets because the order is weird
  df_final$facet = factor(df_final$sample_group,
                          levels = sample_groupings)
  # rename values
  df_final$weight_type = ifelse(df_final$type == "SW",
                                "Source\nshare",
                                "Effect\nshare")
  # reorder weight bars manually to show Source Share, then Effect Share
  df_final$weight_type = factor(
    df_final$weight_type,
    levels = c("Effect\nshare", "Source\nshare")
  )
  # Order signature groups by mean effect share (or use the cannataro order)
  summed_shares = df_final[type == 'CEW' & group != other_label, .(summed_share = sum(prop)), by = c('group', 'sample_group')]
  mean_share_order = summed_shares[, .(mean_share = mean(summed_share)), by = 'group'][order(mean_share), group]
  final_order = c(other_label, mean_share_order)
  df_final$group = factor(df_final$group, levels = final_order)
  
  gg = ggplot(data = df_final) +
    geom_bar(
      mapping = aes(
        y = weight_type,
        fill = group,
        weight = prop
      ),
      position = 'fill',
      width = .9,
      color = 'black'
    ) +
    xlab('Proportion') + ylab('') +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      axis.line.y = element_blank(),
      legend.position = 'right') +
    scale_x_continuous(n.breaks = 10, limits = c(-0.001, 1.001), expand = expansion(add = 0))
  if(! any(is.na(df_final$facet))) {
    gg = gg + facet_wrap(~ facet, ncol = 1, strip.position = 'left') + 
      theme(strip.background = element_blank(), strip.placement = 'outside', 
            strip.text.y.left = element_text(size = 12, angle = 0))
  }
  sig_legend_name = 'Signatures'
  if(all(df_final$name %like% '^SBS')) {
    sig_legend_name = 'Signatures (COSMIC SBS)'
  }
  if(using_cannataro_colors) {
    # Use custom palette outlined in Cannataro et al.
    cannataro_colors = c(
      "Deamination with age, clock-like (1)" = "gray40",
      "Unknown, clock-like (5)" = "gray60",
      "APOBEC (2, 13)" = "#7570b3",
      "Defective homologous recombination (3)" = "#e7298a",
      "Tobacco (4, 29)" = "#a6761d",
      "UV light (7a–d, 38)" = "#e6ab02",
      "Prior treatment (11, 31, 32, 35)" = "#1b9e77",
      "Mutagenic chemical exposure (22, 24, 42, 88)" = "#66a61e",
      "Alcohol-associated (16)" = "#d95f02",
      "Non-actionable and unknown signatures" = "black"
    )
    gg = gg + scale_fill_manual(name = sig_legend_name, values = cannataro_colors)
  } else if(is.null(df_final$color)) {
    gg = gg + scale_fill_viridis_d(name = sig_legend_name, option = viridis_option, begin = .2, direction = -1)
  } else {
    df_final[group == other_label, color := 'white']
    group_to_color = setNames(unique(df_final$color), unique(df_final$group))
    gg = gg + scale_fill_manual(name = sig_legend_name, values = group_to_color)
  }
  
  gg = gg + guides(fill = guide_legend(reverse=TRUE))
  return(gg)
}

#' COSMIC signatures with known etiologies
#'
#' @export
cosmic_sbs_etiologies = function() {
  return(cosmic_sbs_signature_etiology)
}
