#' Plot heatmap
#'
#' Plots an hierarchial heatmap with the package 'ComplexHeatmap' for a given a phyloseq object
#' with 'otu_table', 'tax_table', 'sam_data' (required). If provided a 'tax_rank' the heatmap
#' will be subsetted to the taxonomic rank specified. Samples can be annotated by specifying one
#' or more categorical/factor variables (character) among sample data (in 'sample_variables(physeq)').
#' @param physeq phyloseq-class object with 'otu_table', 'tax_table', 'sam_data' slots (required).
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @param rm_na include ('TRUE') or not ('FALSE') NAs (logical), i.e., taxa without classification at
#' the taxonomic level specified at 'tax_rank'. It only works if 'tax_rank' was provided.
#' Otherwise it is ignored.
#' @param annot_samples one or more categorical/factor variables (character) among sample data
#' (in 'sample_variables(physeq)'). Default is 'NULL'.
#' @param norm_mth normalize/transform the feature table. Default is 'NULL'. One character to pass to
#' 'vegan::decostand()' - one of c("total", "max", "frequency", "normalize", "range", "rank",
#' "standardize", "pa", "chi", "hellinger", "log").
#' @param scale_by scale (Z-score) data by "samples" or "taxa" (character). Default is 'NULL'.
#' @param tr_heat transpose heatmap. Default is 'FALSE' (logical).
#' @param set_seed set seed to allow reproducibility. Default is '1024' (numeric). Set to 'NULL' to
#' turn it off.
#' @param ... parameters to pass to the function 'ComplexHeatmap::Heatmap()' with the exception of
#' 'matrix'. Also 'top_annotation' (samples in columns) or 'left_annotation' (samples in rows)
#' cannot be specified if 'annot_samples' was specified.
#' @return Heatmap (from 'ComplexHeatmap').
#' @export

plot_taxa_heatmap <- function(physeq, tax_rank = NULL,
                              rm_na = FALSE,
                              annot_samples = NULL,
                              norm_mth = NULL,
                              scale_by = NULL, # "samples" or "taxa"
                              tr_heat = FALSE,
                              set_seed = 1024,
                              ...) {

  # Written: 16/09/2020
  # Last update: 16/09/2020

  # packages
  require("phyloseq")
  require("dplyr")
  require("ComplexHeatmap")

  # check input
  if ( !is.null(scale_by) ) stopifnot( scale_by %in% c("samples", "taxa") & length(scale_by)==1 )
  if ( !is.null(norm_mth) ) stopifnot( norm_mth %in% c("total", "max", "frequency",
                                                       "normalize", "range", "rank",
                                                       "standardize", "pa", "chi",
                                                       "hellinger", "log") &
                                         length(norm_mth)==1 )
  if ( !is.null(annot_samples) ) stopifnot( annot_samples %in% sample_variables(physeq) )

  samples <- "rows"
  # Get data
  if ( !is.null(tax_rank) ) { # add features taxa to col names
    feature_w_tax = TRUE
  } else {
    feature_w_tax = FALSE
  }
  physeq_list <- get_physeq_tbls(physeq = physeq, tax_rank = tax_rank,
                                 rm_na = rm_na, feature_w_tax = feature_w_tax)


  # normalize
  if ( !is.null(norm_mth) ) physeq_list[["feature"]] <- vegan::decostand(x = physeq_list[["feature"]],
                                                                         method = norm_mth)

  # z-score
  if ( !is.null(scale_by) ) {
    if ( scale_by == "taxa" ) {
      physeq_list[["feature"]] <- scale(physeq_list[["feature"]])
      } else {
        physeq_list[["feature"]] <- t(scale(t(physeq_list[["feature"]])))
      }
    }

  if ( isTRUE(tr_heat) ) {
    physeq_list[["feature"]] <- t(physeq_list[["feature"]])
    samples <- "cols"
  }

  # annot heatmap samples
  if ( !is.null(annot_samples) ) {
    annot_df <- data.frame(physeq_list[["metadata"]][,annot_samples])
    colnames(annot_df) <- annot_samples; rownames(annot_df) <- rownames(physeq_list[["metadata"]]);
    if ( samples == "rows") {
      annot_samp <- ComplexHeatmap::HeatmapAnnotation(df = annot_df, which = "row")
      # set seed and plot heatmap
      set.seed(set_seed)
      heat_plot <- do.call(ComplexHeatmap::Heatmap,
                           list(matrix = physeq_list[["feature"]],
                                left_annotation = annot_samp,

                                ...))
    } else {
      annot_samp <- ComplexHeatmap::HeatmapAnnotation(df = annot_df)
      # set seed and plot heatmap
      set.seed(set_seed)
      heat_plot <- do.call(ComplexHeatmap::Heatmap,
                           list(matrix = physeq_list[["feature"]],
                                top_annotation = annot_samp,
                                ...))

    }
  } else {
    # set seed and plot heatmap
    set.seed(set_seed)
    heat_plot <- do.call(ComplexHeatmap::Heatmap,
                         list(matrix = physeq_list[["feature"]],
                              ...))
  }

  return(heat_plot)
}

#------------------------------------------------------------------------------------------------------------------------

#' Profile taxa by samples
#'
#' Given a phyloseq object (with absolute abundance counts) and a taxonomic
#' rank (character - one among 'rank_names(physeq)') it gives the taxa (by 'tax_rank')
#' profiled by samples.
#' @param physeq phyloseq-class object.
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @param count_type type of counts to perform 'rank_taxa()'. It can be absolute, i.e.,
#' 'abs', or percentage, i.e., 'perc' (character). Default is 'abs' - absolute counts.
#' @param col_bar column bar to annotate samples (character). One character/factor
#' variable among 'sample_variables(physeq)' to be used to annotate the samples.
#' The maximum number of levels that the character/factor variable can have is 10.
#' @param top_taxa top most abundant taxa to select by 'tax_rank' (numeric - coerced
#' to integer) in the overall data set. Default 'NULL'. If 'count_type = perc', the
#' percentage ranging from 0-100\% will be determined for the 'top_taxa' only.
#' @param show_top top n taxa to show (numeric - coerced to integer). Default is 'NULL'.
#' The difference to 'top_taxa' is that if 'count_type = "perc"', the percentage presented
#' still refers to the percentage of overall data before discarding the non-top n taxa.
#' @param group factorial variable to group, one among 'sample_variables(physeq)'.
#' Default is 'NULL'.
#' @param taxa_perc_cutoff percentage cutoff to remove less abundant taxa than
#' 'taxa_perc_cutoff'. Default is 'NULL'.
#' @param ord_by logical (or character). Default is 'FALSE'. Order the samples (x-axis)
#' by the levels from the character/factorial variable indicated at 'col_bar'. Only
#' applied if 'col_bar' parameter is used. If provided a character of the 'col_bar'
#' factor levels, the samples will be ordered by the 'col_bar' factor level provided.
#' @param rm_na include (TRUE) or not (FALSE) NAs, i.e., taxa without classification at
#' the taxonomic level specified at 'tax_rank'.
#' @param fill_other
#' @param ... parameters to be passed to the function 'filter_feature_table()'.
#' @return It returns a list with a tibble data frame ('$data') and ggplot ('$plot')
#' profiling the taxa picked at the taxonomic rank select at 'tax_rank' by samples.
#' @export

profile_taxa_by_samples <- function(physeq, tax_rank, count_type = "abs",
                                    col_bar = NULL, top_taxa = NULL, show_top = NULL,
                                    group = NULL, taxa_perc_cutoff = NULL,
                                    ord_by = FALSE, rm_na = FALSE, fill_other = FALSE,
                                    ...) {

  # Written: 27/09/2020
  # Last update: 27/09/2020

  # packages
  require("phyloseq")
  require("dplyr")
  require("ggplot2")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(deparse(substitute(physeq)), " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table") %in% slotNames(physeq) ) ) { # check the 2 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 2 slots: 'otu_table',
    'tax_table'!\nThe 2 slots are necessary to use the 'rank_taxa' function.\nAborting..."))
  }
  stopifnot( all(c(is.character(tax_rank), is.character(count_type))) )
  stopifnot( all(c(length(tax_rank) == 1, tax_rank %in% rank_names(physeq))) )
  stopifnot( count_type %in% c("abs", "perc") )
  stopifnot( !all(c(!is.null(show_top), !is.null(taxa_perc_cutoff))) )
  stopifnot( c(is.logical(rm_na), is.logical(fill_other)) )
  if ( !is.null(group) ) stopifnot( all(c(length(group) == 1, group %in% sample_variables(physeq))) )
  if ( !is.null(col_bar) ) stopifnot( all(c(length(col_bar) == 1, group %in% sample_variables(physeq))) )
  if ( !is.null(show_top) ) stopifnot( all(c(length(show_top)==1, is.numeric(show_top))) )
  if ( !is.null(taxa_perc_cutoff) ) stopifnot( all(c(length(taxa_perc_cutoff)==1, is.numeric(taxa_perc_cutoff)), taxa_perc_cutoff < 100) )
  checkInteger <- (physeq@otu_table@.Data%%1==0) # to deal with double like, e.g., 1 (by default - double)
  #and not as integer (coerce 1L)
  if ( !all(checkInteger == TRUE) ) stop(paste0("otu_table(", deparse(substitute(physeq)), ") is not integer!"))

  #-----------------------------------------------------------------------------------------------------
  # Future implementation should allow
  #different variables in the x-axis
  x_var = "Sample"
  #-----------------------------------------------------------------------------------------------------

  # tax glom by 'tax_rank'
  physeq_rank <- tax_glom(physeq = physeq, taxrank = tax_rank, NArm = rm_na)
  physeq_rank <- do.call(filter_feature_table, list(physeq = physeq_rank, ...))

  if ( !is.null(top_taxa) ) { # get 'top_taxa'
    stopifnot( is.numeric(top_taxa) )
    top_taxa <- as.integer(top_taxa)
    top_taxa <- names(sort(taxa_sums(x = physeq_rank), TRUE)[1:top_taxa])
    physeq_rank   <- prune_taxa(taxa = top_taxa, x = physeq_rank)
  }

  # melt data
  physeq_rank_melted <- melt_physeq(physeq = physeq_rank)
  # rank by abundance
  if ( is.null(group) ) { # plot rank abundance for the overall data set
    # data
    physeq_rank_melted_2_plot <-
      physeq_rank_melted %>%
      select(.data[[x_var]], .data[[tax_rank]], Abundance) %>%
      group_by(.data[[x_var]], .data[[tax_rank]]) %>%
      summarise(Abundance = sum(Abundance)) %>%
      filter(Abundance != 0) %>% # rm phyla without abundance
      arrange(desc(Abundance)) %>%
      ungroup(.) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(!!tax_rank := factor(.data[[tax_rank]],
                                  levels = unique(.data[[tax_rank]])))

    if ( count_type == "abs" ) {
      y_lab <- "Absolute abundance"
      abundance <- "Abundance"
      y_max <- max(physeq_rank_melted_2_plot[, "Abundance", drop = TRUE]) * 1.1
    } else {
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        group_by(.data[[x_var]]) %>%
        mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) )
      y_lab <- "Percentage (%)"
      abundance <- "Percentage"
      y_max <- max(physeq_rank_melted_2_plot[, "Percentage", drop = TRUE]) + 5
    }

    # show top taxa
    if ( !is.null(show_top) ) {
      show_top <- as.numeric(show_top)
      list_taxa <- physeq_rank_melted_2_plot %>%
        pull(.data[[tax_rank]]) %>%
        levels(.)
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        filter(.data[[tax_rank]] %in% list_taxa[1:show_top])
      if ( fill_other ) { # fill the remaining empth space with "Other"
        physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
          ungroup(.) %>%
          mutate_if(is.factor, as.character)
        physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
          group_by(.data[[x_var]]) %>%
          summarize("Percentage" = 100 - sum(Percentage)) %>%
          mutate(!!tax_rank := "Other") %>%
          bind_rows(., physeq_rank_melted_2_plot)
        }
      }
    # filter taxa based on percentage
    if ( !is.null(taxa_perc_cutoff) ) {
      if ( count_type == "abs" ) {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          group_by(.data[[x_var]]) %>%
          mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) ) %>%
          filter(Percentage >= taxa_perc_cutoff) %>%
          select(-Percentage)
      } else {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          filter(Percentage >= taxa_perc_cutoff)
      }
    }
    if ( is.null(col_bar) ) {
    # plot
    physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
      mutate(!!tax_rank := factor(.data[[tax_rank]],
                                  levels =
                                    physeq_rank_melted_2_plot %>% ungroup() %>%
                                    group_by(.data[[tax_rank]]) %>%
                                    summarise("Sum" = sum(.data[[abundance]])) %>%
                                    arrange(Sum) %>%
                                    pull(.data[[tax_rank]]) %>%
                                    as.character()
        ))
    taxa_barplot_ranked <-
      ggplot(data = physeq_rank_melted_2_plot, aes_string(x = x_var,
                                                          y = abundance,
                                                          fill = tax_rank)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab(y_lab) #+
      #geom_text(aes_string(label = abundance),
      #          angle = 45, vjust=0, hjust = 0, size=3) +
      #ylim(c(0, y_max))
    }
  } else { # plot rank abundance by group
    # samples by group
    n_samp_by_group <- physeq_rank_melted %>%
      group_by(.data[[group]]) %>%
      summarise("n" = n_distinct(Sample)) %>%
      arrange(desc(n)) %>%
      filter(!is.na(.data[[group]]))
    facet_label <- paste0(n_samp_by_group[,group, drop = TRUE],
                          " (n=", n_samp_by_group$n, ")")
    names(facet_label) <- n_samp_by_group[,group, drop = TRUE]
    #data
    physeq_rank_melted_2_plot <-
      physeq_rank_melted %>%
      group_by(.data[[group]], .data[[x_var]], .data[[tax_rank]]) %>%
      summarise(Abundance = sum(Abundance)) %>%
      filter(Abundance != 0) %>% # rm phyla without abundance
      arrange(.data[[group]], desc(Abundance)) %>%
      ungroup(.) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(!!tax_rank := factor(.data[[tax_rank]],
                                  levels = unique(.data[[tax_rank]]))) %>%
      mutate(!!group := factor(.data[[group]],
                               levels = unique(n_samp_by_group[,group, drop = TRUE])))

    if ( count_type == "abs" ) {
      y_lab <- "Absolute abundance"
      abundance <- "Abundance"
      #y_max <- max(physeq_rank_melted_2_plot[, "Abundance", drop = TRUE]) * 1.1
    } else {
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        group_by(.data[[group]], .data[[x_var]]) %>%
        mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) )
      y_lab <- "Percentage (%)"
      abundance <- "Percentage"
      #y_max <- max(physeq_rank_melted_2_plot[, "Percentage", drop = TRUE]) + 5
    }

    # show top taxa
    if ( !is.null(show_top) ) {
      show_top <- as.numeric(show_top)
      list_taxa <- physeq_rank_melted_2_plot %>%
        pull(.data[[tax_rank]]) %>%
        levels(.)
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        filter(.data[[tax_rank]] %in% list_taxa[1:show_top])
      if ( fill_other ) { # fill the remaining empth space with "Other"
        physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
          ungroup(.) %>%
          mutate_if(is.factor, as.character)
        physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
          group_by(.data[[group]], .data[[x_var]]) %>%
          summarize("Percentage" = 100 - sum(Percentage)) %>%
          mutate(!!tax_rank := "Other") %>%
          bind_rows(., physeq_rank_melted_2_plot) %>%
          ungroup(.) %>%
          mutate(!!group := factor(.data[[group]],
                                   levels = unique(n_samp_by_group[,group, drop = TRUE]))) %>%
          group_by(.data[[group]], .data[[x_var]])
      }
    }
    # filter taxa based on percentage
    if ( !is.null(taxa_perc_cutoff) ) {
      if ( count_type == "abs" ) {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          group_by(.data[[group]], .data[[x_var]]) %>%
          mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) ) %>%
          filter(Percentage >= taxa_perc_cutoff) %>%
          select(-Percentage)
      } else {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          filter(Percentage >= taxa_perc_cutoff)
      }
    }
    if ( is.null(col_bar) ) {
      # plot
      if ( fill_other & !is.null(show_top) ) {
        tax_level_order <- physeq_rank_melted_2_plot %>% ungroup() %>%
          group_by(.data[[tax_rank]]) %>%
          summarise("Sum" = sum(.data[[abundance]])) %>%
          arrange(Sum) %>%
          pull(.data[[tax_rank]]) %>%
          as.character()
        tax_level_order <- c(tax_level_order[tax_level_order!="Other"], "Other")
        physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
          mutate(!!tax_rank := factor(.data[[tax_rank]],
                                      levels = tax_level_order ))
      } else {
        physeq_rank_melted_2_plot <- physeq_rank_melted_2_plot %>%
          mutate(!!tax_rank := factor(.data[[tax_rank]],
                                      levels =
                                        physeq_rank_melted_2_plot %>% ungroup() %>%
                                        group_by(.data[[tax_rank]]) %>%
                                        summarise("Sum" = sum(.data[[abundance]])) %>%
                                        arrange(Sum) %>%
                                        pull(.data[[tax_rank]]) %>%
                                        as.character()
          ))
    }
      taxa_barplot_ranked <-
        ggplot(data = physeq_rank_melted_2_plot, aes_string(x = x_var,
                                                            y = abundance,
                                                            fill = tax_rank)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ylab(y_lab) +
        facet_wrap(as.formula(paste("~", group)), scales = "free",
                   labeller = labeller(.cols = facet_label),
                   ncol = 3)
    }
  }

  ## Save output or plot bar annotation (and save output)
  if ( is.null(col_bar) ) {
    # save output
    listOut <- list(data = physeq_rank_melted_2_plot,
                    plot = taxa_barplot_ranked)
  } else {  ## Add segment annotation bar
    if ( is.null(group) ) facet_label <- NULL
    listOut <- .profile_taxa_bar_plot(data2plot = physeq_rank_melted_2_plot,
                                      melted_physeq = physeq_rank_melted,
                                      x_var = x_var, col_bar = col_bar,
                                      y_var = tax_rank, group = group,
                                      count_type = count_type, ord_by = ord_by,
                                      facet_label = facet_label,
                                      fill_other = fill_other)
  }

  return(listOut)
}

# helper (hidden) function
.profile_taxa_bar_plot <- function(data2plot, melted_physeq, x_var, col_bar, y_var,
                                   count_type = "abs", group = NULL, ord_by = FALSE,
                                   facet_label = NULL, fill_other = FALSE) {

  ## packages
  require("ggplot2")
  require("dplyr")
  #require("gsci")

  # which annotation vars to select
  var_cols <- col_bar
  cols_names <- "col_bar"
  if ( !is.null(group) ) {
    var_cols <- c(group, col_bar)
    cols_names <- c("group", "col_bar")
  }

  # perc vs. abs y-axis column selection
  if ( count_type == "abs" ) {
    abundance = "Abundance"
    y_lab <- "Absolute abundance"
  } else {
    abundance = "Percentage"
    y_lab <- "Percentage (%)"
  }

  ## descending order "y_var"
  if ( fill_other ) {
    tax_level_order <- data2plot %>% ungroup() %>%
      group_by(.data[[y_var]]) %>%
      summarise("Sum" = sum(.data[[abundance]])) %>%
      #summarise("Sum" = sum(Abundance)) %>% # aggs: order fct var by abs. abundance
      arrange(Sum) %>%
      pull(.data[[y_var]]) %>%
      as.character()
    tax_level_order <- c(tax_level_order[tax_level_order!="Other"], "Other")
    data2plot <- data2plot %>%
      mutate(!!y_var := factor(.data[[y_var]], levels = tax_level_order ))
    } else {
      data2plot <- data2plot %>%
        mutate(!!y_var := factor(.data[[y_var]],
                                 levels =
                                   data2plot %>% ungroup() %>%
                                   group_by(.data[[y_var]]) %>%
                                   summarise("Sum" = sum(.data[[abundance]])) %>%
                                   arrange(Sum) %>%
                                   pull(.data[[y_var]]) %>%
                                   as.character()
        ))
    }

  ## Annotation df
  annot_df <- data.frame(sort(unique(data2plot[,x_var, drop=TRUE]) ))
  colnames(annot_df) <- x_var
  annot_df <- left_join(annot_df, distinct(melted_physeq[,c(x_var, var_cols)]),
                        by = x_var)
  annot_df <- annot_df %>%
    mutate_if(is.character,as.factor)
  colnames(annot_df)[ (ncol(annot_df)+1-length(var_cols)) : ncol(annot_df) ] <- cols_names;

  if ( !isFALSE(ord_by) ) { # order 'x_var' by 'col_bar'
    if( !is.null(group) ) {
      annot_df <- annot_df %>%
        group_by(group) %>%
        `if`( is.character(ord_by), mutate(., col_bar = factor(col_bar, levels = ord_by)), .) %>%
        arrange(col_bar)
      annot_df <- annot_df %>%
        group_by(group) %>%
        mutate("x_axis" = order(col_bar))
      x_var_levels <- as.character(annot_df[,x_var, drop = TRUE]) # order 'x_var' fct var
    } else {
      annot_df <- annot_df %>%
        arrange(col_bar) %>%
        mutate("x_axis" = order(col_bar))
      x_var_levels <- as.character(annot_df[,x_var, drop = TRUE]) # order 'x_var' fct var
    }
  } else if ( isFALSE(ord_by) & !is.null(col_bar) & !is.null(group) ) {
    annot_df <- annot_df %>%
      group_by(group) %>%
      mutate("x_axis" = order(col_bar))
    x_var_levels <- annot_df %>% pull(.data[[x_var]]) %>%
      levels(.) # order 'x_var' fct var
  } else {
    annot_df <- annot_df %>%
      mutate("x_axis" = order(.data[[x_var]]))
    x_var_levels <- annot_df %>% pull(.data[[x_var]]) %>%
      levels(.) # order 'x_var' fct var
  }
  data2plot <- left_join(data2plot, annot_df, by = x_var)
  data2plot <- data2plot %>% ungroup() %>%
    mutate(!!x_var := factor(.data[[x_var]], levels = x_var_levels))
  if ( is.null(group) ) {
    data2plot$y_axis <- data2plot %>% group_by(.data[[x_var]]) %>%
      summarise("Sum" = sum(.data[[abundance]])) %>% pull(Sum) %>% max(.)
    taxa_barplot_ranked  <- ggplot(data = data2plot,
                                   aes_string(x = x_var,
                                              y = abundance,
                                              fill = y_var)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab(y_lab) +
      geom_segment(aes(x = x_axis - 0.5, xend = x_axis + 0.5,
                       y = y_axis * 1.05, yend = y_axis * 1.05,
                       color = col_bar), size = 3) +
      scale_color_manual(name = col_bar, values = ggsci::pal_npg()(10))
  } else {
     y_axis_data <- data2plot %>% group_by(group, .data[[x_var]]) %>%
      summarise("Sum" = sum(.data[[abundance]])) %>%
      summarise("y_axis" = max(Sum))
     width_facets <- annot_df %>% group_by(group) %>%
       distinct(.data[[x_var]]) %>%  tally() %>%
       mutate("facet_width" = n / max(n)) %>%
       select(-n)
     y_axis_data <- left_join(y_axis_data, width_facets, by = "group")
     data2plot <- left_join(data2plot, y_axis_data, by = "group")
     taxa_barplot_ranked  <- ggplot(data = data2plot,
                                   aes(x = .data[[x_var]],
                                       y = .data[[abundance]],
                                       fill = .data[[y_var]],
                                       width = .9 * facet_width)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab(y_lab) +
      geom_segment(aes(x = x_axis - 0.5, xend = x_axis + 0.5,
                       y = y_axis * 1.05, yend = y_axis * 1.05,
                       color = col_bar), size = 3) +
      scale_color_manual(name = col_bar, values = ggsci::pal_npg()(10)) +
      facet_wrap(as.formula(paste("~", group)), scales = "free",
                 labeller = labeller(.cols = facet_label),
                 ncol = 3)
     data2plot <- data2plot %>% select(-facet_width)
  }

  ## Format and save result as a list
  data2plot <- data2plot %>%
    select(-c(x_axis, y_axis))
  colnames(data2plot)[ (ncol(data2plot)+1-length(var_cols)) : ncol(data2plot) ] <- var_cols
  data2plot <- data2plot[,!duplicated(colnames(data2plot))]
  # save output
  listOut <- list(data = data2plot,
                  plot = taxa_barplot_ranked)

  return(listOut)
}

# helper function to add extra layers of row annotation to the profile_taxa_bar_plot()
# .add_col_bar <- function(plot, physeq, extra_col_bar) {
#   data_from_plot <- plot$data
#   metadata <- data.frame("Sample" = rownames(sample_data(physeq)),
#                          "Extra_bar" = get_variable(physeq = physeq,
#                                                     varName = extra_col_bar))
#   data_to_plot <- left_join(data_from_plot, metadata, by = "Sample")
#   plot$data <- data_to_plot
#   Plot <- plot + geom_segment( aes(x = x_axis - 0.5, xend = x_axis + 0.5,
#                                    y = y_axis * 1.15, yend = y_axis * 1.15,
#                                    colour = Extra_bar), size = 3) +
#     scale_color_brewer(palette = "Spectral")
#   return(Plot)
# }

#----------------------------------------------------------------------------------------------------------------

#' Rank taxa
#'
#' Given a phyloseq object (with absolute abundance counts) and a taxonomic
#' rank (character - one among 'rank_names(physeq)') it gives the taxa (by 'tax_rank')
#' summarized for the whole data set ranked by abundance.
#' @param physeq phyloseq-class object.
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @param count_type type of counts to perform 'rank_taxa()'. It can be absolute, i.e.,
#' 'abs', or percentage, i.e., 'perc' (character). Default is 'abs' - absolute counts.
#' @param top_taxa top most abundant taxa to select by 'tax_rank' (numeric - coerced
#' to integer) in the overall data set. Default 'NULL'. If 'count_type = perc', the
#' percentage ranging from 0-100\% will be determined for the 'top_taxa' only.
#' @param show_top top n taxa to show (numeric - coerced to integer). Default is 'NULL'.
#' The difference to 'top_taxa' is that if 'count_type = "perc"', the percentage presented
#' still refers to the percentage of overall data before discarding the non-top n taxa.
#' @param group factorial variable to group, one among 'sample_variables(physeq)'.
#' Default is 'NULL'.
#' @param taxa_perc_cutoff percentage cutoff to remove less abundant taxa than
#' 'taxa_perc_cutoff'. Default is 'NULL'.
#' @param rm_na include (TRUE) or not (FALSE) NAs, i.e., taxa without classification at
#' the taxonomic level specified at 'tax_rank'.
#' @param ... parameters to be passed to the function 'filter_feature_table()'.
#' @return It returns a list with a tibble data frame ('$data') and ggplot ('$plot')
#' ranking the taxa picked at the taxonomic rank select at 'tax_rank'.
#' @export

rank_taxa <- function(physeq, tax_rank, count_type = "abs",
                      top_taxa = NULL, show_top = NULL,
                      group = NULL, taxa_perc_cutoff = NULL,
                      rm_na = FALSE, ...) {

  # packages
  require("phyloseq")
  require("dplyr")
  require("ggplot2")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(deparse(substitute(physeq)), " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table") %in% slotNames(physeq) ) ) { # check the 2 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 2 slots: 'otu_table',
    'tax_table'!\nThe 2 slots are necessary to use the 'rank_taxa' function.\nAborting..."))
  }
  stopifnot( all(c(is.character(tax_rank), is.character(count_type))) )
  stopifnot( all(c(length(tax_rank) == 1, tax_rank %in% rank_names(physeq))) )
  stopifnot( count_type %in% c("abs", "perc") )
  stopifnot( !all(c(!is.null(show_top), !is.null(taxa_perc_cutoff))) )
  stopifnot( is.logical(rm_na) )
  if ( !is.null(group) ) stopifnot( all(c(length(group) == 1, group %in% sample_variables(physeq))) )
  if ( !is.null(show_top) ) stopifnot( all(c(length(show_top)==1, is.numeric(show_top))) )
  if ( !is.null(taxa_perc_cutoff) ) stopifnot( all(c(length(taxa_perc_cutoff)==1, is.numeric(taxa_perc_cutoff)), taxa_perc_cutoff < 100) )
  checkInteger <- (physeq@otu_table@.Data%%1==0) # to deal with double like, e.g., 1 (by default - double)
  #and not as integer (coerce 1L)
  if ( !all(checkInteger == TRUE) ) stop(paste0("otu_table(", deparse(substitute(physeq)), ") is not integer!"))

  # tax glom by 'tax_rank'
  physeq_rank <- tax_glom(physeq = physeq, taxrank = tax_rank, NArm = rm_na)
  physeq_rank <- do.call(filter_feature_table, list(physeq = physeq_rank, ...))

  if ( !is.null(top_taxa) ) { # get 'top_taxa'
    stopifnot( is.numeric(top_taxa) )
    top_taxa <- as.integer(top_taxa)
    top_taxa <- names(sort(taxa_sums(x = physeq_rank), TRUE)[1:top_taxa])
    physeq_rank   <- prune_taxa(taxa = top_taxa, x = physeq_rank)
  }

  # melt data
  physeq_rank_melted <- melt_physeq(physeq = physeq_rank)
  # rank by abundance
  if ( is.null(group) ) { # plot rank abundance for the overall data set
    # data
    physeq_rank_melted_2_plot <-
      physeq_rank_melted %>%
      select(.data[[tax_rank]], Abundance) %>%
      group_by(.data[[tax_rank]]) %>%
      summarise(Abundance = sum(Abundance)) %>%
      filter(Abundance != 0) %>% # rm phyla without abundance
      arrange(desc(Abundance)) %>%
      ungroup(.) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(!!tax_rank := factor(.data[[tax_rank]],
                                  levels = unique(.data[[tax_rank]])))

    if ( count_type == "abs" ) {
      y_lab <- "Absolute abundance"
      abundance <- "Abundance"
      y_max <- max(physeq_rank_melted_2_plot[, "Abundance", drop = TRUE]) * 1.1
    } else {
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) )
      y_lab <- "Percentage (%)"
      abundance <- "Percentage"
      y_max <- max(physeq_rank_melted_2_plot[, "Percentage", drop = TRUE]) + 5
    }

    # show top taxa
    if ( !is.null(show_top) ) {
      show_top <- as.numeric(show_top)
      list_taxa <- physeq_rank_melted_2_plot %>%
        pull(.data[[tax_rank]]) %>%
        levels(.)
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        filter(.data[[tax_rank]] %in% list_taxa[1:show_top])
    }
    # filter taxa based on percentage
    if ( !is.null(taxa_perc_cutoff) ) {
      if ( count_type == "abs" ) {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) ) %>%
          filter(Percentage >= taxa_perc_cutoff) %>%
          select(-Percentage)
        } else {
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        filter(Percentage >= taxa_perc_cutoff)
        }
      }

    # plot
    taxa_barplot_ranked <-
      ggplot(data = physeq_rank_melted_2_plot, aes_string(x = tax_rank,
                                                          y = abundance)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab(y_lab) +
      geom_text(aes_string(label = abundance),
                angle = 45, vjust=0, hjust = 0, size=3) +
      ylim(c(0, y_max))
  } else { # plot rank abundance by group
    # samples by group
    n_samp_by_group <- physeq_rank_melted %>%
      group_by(.data[[group]]) %>%
      summarise("n" = n_distinct(Sample)) %>%
      arrange(desc(n)) %>%
      filter(!is.na(.data[[group]]))
    facet_label <- paste0(n_samp_by_group[,group, drop = TRUE],
                          " (n=", n_samp_by_group$n, ")")
    names(facet_label) <- n_samp_by_group[,group, drop = TRUE]
    #data
    physeq_rank_melted_2_plot <-
      physeq_rank_melted %>%
      group_by(.data[[group]], .data[[tax_rank]]) %>%
      summarise(Abundance = sum(Abundance)) %>%
      filter(Abundance != 0) %>% # rm phyla without abundance
      arrange(.data[[group]], desc(Abundance)) %>%
      ungroup(.) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(!!tax_rank := factor(.data[[tax_rank]],
                                  levels = unique(.data[[tax_rank]]))) %>%
      mutate(!!group := factor(.data[[group]],
                               levels = unique(n_samp_by_group[,group, drop = TRUE])))

    if ( count_type == "abs" ) {
      y_lab <- "Absolute abundance"
      abundance <- "Abundance"
      #y_max <- max(physeq_rank_melted_2_plot[, "Abundance", drop = TRUE]) * 1.1
    } else {
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        group_by(.data[[group]]) %>%
        mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) )
      y_lab <- "Percentage (%)"
      abundance <- "Percentage"
      #y_max <- max(physeq_rank_melted_2_plot[, "Percentage", drop = TRUE]) + 5
    }

    # show top taxa
    if ( !is.null(show_top) ) {
      list_taxa <- physeq_rank_melted_2_plot %>%
        pull(.data[[tax_rank]]) %>%
        levels(.)
      physeq_rank_melted_2_plot <-
        physeq_rank_melted_2_plot %>%
        filter(.data[[tax_rank]] %in% list_taxa[1:show_top])
    }
    # filter taxa based on percentage
    if ( !is.null(taxa_perc_cutoff) ) {
      if ( count_type == "abs" ) {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          group_by(.data[[group]]) %>%
          mutate("Percentage" = round( Abundance / sum( Abundance ) * 100, 2) ) %>%
          filter(Percentage >= taxa_perc_cutoff) %>%
          select(-Percentage)
        } else {
        physeq_rank_melted_2_plot <-
          physeq_rank_melted_2_plot %>%
          filter(Percentage >= taxa_perc_cutoff)
        }
      }

    # plot
    taxa_barplot_ranked <-
      ggplot(data = physeq_rank_melted_2_plot, aes_string(x = tax_rank,
                                                          y = abundance)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab(y_lab) +
      #geom_text(aes_string(label = abundance),
      #          angle = 45, vjust=0, hjust = 0, size=3) +
      facet_wrap(as.formula(paste("~", group)), scales = "free_y",
                 labeller = labeller(.cols = facet_label),
                 ncol = 3)
  }

  # save output
  listOut <- list(data = physeq_rank_melted_2_plot,
                  plot = taxa_barplot_ranked)

  return(listOut)
}

#---------------------------------------------------------------------------------------------------------------

#' Filter phyloseq object based on samples and taxa.
#'
#' Filter phyloseq object based on samples and taxa names or minimum number of reads.
#' @param physeq phyloseq-class object.
#' @param taxa2filter taxonomic names to exclude (character). It can include the
#' taxa names, i.e., OTU/ASV ids - 'taxa_names(physeq)'. Also it can be instead a logical
#' of the same length as OTU/ASV ids (i.e., 'length(taxa_names(physeq))'), where 'TRUE'
#' represents the taxa to keep and 'FALSE' to remove. Default 'NULL'.
#' @param samples2filter name of the samples to discard (character). Default 'NULL'.
#' @param min_sample_depth minimum number of reads per sample to keep a sample (numeric). Default '0'.
#' @param min_taxa_counts minimum number of reads per feature i.e., OTU/ASV, to keep
#' a feature (features sum across samples) (numeric). Default '0'.
#' @return A phyloseq object filtered.
#' @export

filter_feature_table <- function(physeq, taxa2filter = NULL,
                                 samples2filter = NULL, min_sample_depth = 0,
                                 min_taxa_counts = 0) {
  # packages
  require("phyloseq")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(deparse(substitute(physeq)), " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table") %in% slotNames(physeq) ) ) { # check the 2 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 2 slots: 'otu_table',
    'tax_table'!\nThe 2 slots are necessary to use the 'filter_feature_table' function.\nAborting..."))
  }
  stopifnot( is.numeric(c(min_sample_depth, min_taxa_counts)) ) # check if input is numeric

  # check if otu table is integer
  checkInteger <- (physeq@otu_table@.Data%%1==0) # to deal with double like, e.g., 1 (by default - double)
  #and not as integer (coerce 1L)
  if ( !all(checkInteger == TRUE) ) stop(paste0("otu_table(", deparse(substitute(physeq)), ") is not integer!"))

  ## Filter samples
  #
  if ( !is.null(samples2filter) ) {
    stopifnot( all(samples2filter %in% sample_names(physeq)) )
    physeq <- prune_samples(samples = !(sample_names(physeq) %in% samples2filter), x = physeq)
  }
  samples2keep <- sample_sums(physeq) > min_sample_depth
  physeq <- prune_samples(samples = samples2keep, x = physeq)

  ## Filter taxa
  #
  `%notin%` <- Negate(`%in%`)
  if ( !is.null(taxa2filter) ) {
    stopifnot( (class(taxa2filter) %in% c("character", "logical")) )
    if( is.logical(taxa2filter) & (length(taxa2filter) == length(taxa_names(physeq))) ) {
      physeq <- prune_taxa(taxa = taxa2keep, x = physeq)
    } else if ( is.character(taxa2filter) & all( taxa2filter %in% taxa_names(physeq)) ) {
      physeq <- prune_taxa(taxa = taxa2keep, x = physeq)
    } else if ( is.character(taxa2filter) & any(taxa2filter %notin% taxa_names(physeq)) ) {
      full_tax_tbl <- cbind(physeq@tax_table@.Data, rownames(physeq@tax_table@.Data))
      taxa2keep <- apply(full_tax_tbl, 1, function(rank)
        ifelse( length(grep(paste(taxa2filter, collapse = "|"), rank))>0, FALSE, TRUE) )
      physeq <- prune_taxa(taxa = taxa2keep, x = physeq)
    }
  }
  taxa2keep <- taxa_sums(physeq) > min_taxa_counts
  physeq <- prune_taxa(taxa = taxa2keep, x = physeq)

  return(physeq)
}
