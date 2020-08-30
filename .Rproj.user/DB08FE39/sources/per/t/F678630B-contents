#setClass("micrUBI", representation(name = "character", data = "data.frame", plot = c("gg", "ggplot")))

#' Within beta-diversity boxplot
#'
#' Given a phyloseq object, a beta-diversity metric and a group/factorial variable
#' to group the samples, it calculates the within group beta-diversity.
#' @param physeq phyloseq-class object.
#' @param method beta-diversity metric. Default "bray", i.e., Bray-Curtis dissimilarity.
#' @param group factorial variable to group.
#' @return A list with the data frame and with the box plot of within group beta-diversity.
#' @export

beta_boxplot <- function(physeq, method = "bray", group) {

  ## Packages
  require("phyloseq") # v.1.30.0
  require("ggplot2") # v.3.3.2

  ## Identify the correspondence: group and samples
  group2samp <- list() # list to save the correspondence between group <--> samples
  group_list <- get_variable(sample_data(physeq), group) # list of group elements
  for (groups in levels(group_list)) { # loop over the no. of group levels
    target_group <- which(group_list == groups) # vct pos of the curr group variable
    group2samp[[ groups ]] <- sample_names(physeq)[target_group] # matching samples: based on vct pos
  }

  ## Calculate beta-diversity
  beta_div_dist <- distance(physeq = physeq, method = method)
  beta_div_dist <- as(beta_div_dist, "matrix")

  ## Coerce distance mtx into a tidy data frame by group
  dist_df <- data.frame() # save results in df
  counter <- 1
  for (groups in names(group2samp)) { # loop over group fct levels
    sub_dist <- beta_div_dist[ group2samp[[groups]], group2samp[[groups]] ] # subset dist mtx to curr samples
    #print(sub_dist)
    no_samp_col <- ncol(sub_dist) # n cols: curr sub dist
    no_samp_row <- nrow(sub_dist) # n rows: curr sub dist
    for ( cols in seq(no_samp_col) ) { # loop over cols: curr sub_dist
      if ( cols > 1 ) {
        for ( rows in seq((cols-1)) ) { # loop over rows: curr sub_dist
          ## Save results
          dist_df[ counter, "sample_pair" ] <- paste0( colnames(sub_dist)[cols], "-",
                                                       rownames(sub_dist)[rows] ) # sample pair
          dist_df[ counter, "group" ] <- groups # group
          dist_df[ counter, "beta_div_method" ] <- method # method
          dist_df[ counter, "beta_div_value" ] <- sub_dist[rows, cols] # beta-diversity for the sample pair
          counter = counter + 1
        }
      }
    }
  }

  ## Create a ggplot2 boxplot
  plot_boxplot <- ggplot(data = dist_df, aes(x = group, y = beta_div_value, color = group)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter() +
    theme_bw() +
    xlab(group) + ylab(method) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  ## Save df and boxplot into a list
  list_Out <- list("data" = dist_df, "plot" = plot_boxplot)
  #micrUBI <- new(Class = "micrUBI", data = dist_df, plot = plot_boxplot)

  return(list_Out)
}
