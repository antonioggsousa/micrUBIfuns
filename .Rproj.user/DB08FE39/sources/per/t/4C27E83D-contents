#' Prevalence of features (OTUs/ASVs)
#'
#' Given a phyloseq object and a taxonomic rank (character) it gives the
#' prevalence of features (OTUs/ASVs), i.e., the number of samples where
#' each feature (OTUs/ASVs) appear (observed).
#' @param physeq phyloseq-class object.
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @return A data frame with Feature_ID (lists the unique features (OTUs/ASVs)),
#' "Feature_[tax_rank]" (lists the taxonomic affiliation at the taxonomic rank defined
#' by 'tax_rank', Feature_Abundance (absolute abundance of each feature),
#' Feature_Relative (relative abundance of each feature), Prevalence (prevalence of
#' each feature - the no. of samples where each feature is observed),
#' Prevalence_Relative (the same as Prevalence, but in relative abundance).
#' @export

prevalence <- function(physeq, tax_rank) {

  # packages
  require("phyloseq")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(deparse(substitute(physeq)), " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table") %in% slotNames(physeq) ) ) { # check the 2 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 2 slots: 'otu_table',
    'tax_table'!\nThe 2 slots are necessary to use the 'prevalence' function.\nAborting..."))
  }

  # retrieve otu tbl
  if ( taxa_are_rows(physeq) & all( rownames(otu_table(physeq)) == rownames(tax_table(physeq)) ) ) { # check if
    #taxa are rows & ASV/OTU names match between otu_table() and tax_table()
    otu_tbl <- physeq@otu_table@.Data
  } else { # retrieve mtx and if not 'taxa_are_rows' tranpose
    stopifnot( all( colnames(otu_table(physeq)) == rownames(tax_table(physeq)) ) )
    otu_tbl <- t(physeq@otu_table@.Data)
  }

  # check if otu table is integer
  #checkInteger <- apply(otu_tbl, 2, function(x) is.integer(x)) # check if mtx 'x' are intergers
  checkInteger <- (otu_tbl%%1==0) # to deal with double like, e.g., 1 (by default - double)
  #and not as integer (coerce 1L)
  if ( !all(checkInteger == TRUE) ) stop(paste0("otu_table(", deparse(substitute(physeq)), ") is not integer!"))

  # retrieve data for each row: OTU/ASV/feature
  features <- rownames(otu_tbl) # names/ids
  feature_abund <- rowSums(otu_tbl) # abundance (sum across all samples)
  total_reads <- sum(feature_abund) # total abundance
  feature_rel <- feature_abund / total_reads # relative abundance
  prevalence <- rowSums(otu_tbl>=1) # prevalence
  prevalence_rel <- prevalence / ncol(otu_tbl) # relative prevalence
  tax_names <- physeq@tax_table@.Data[,tax_rank] # taxonomic classification
  feature_tax <- paste0("Feature_", get("tax_rank"))

  # save result as df
  df_out <- data.frame("Feature_ID" = features,
                       "Feature_Taxa" = tax_names,
                       "Feature_Abundance" = feature_abund,
                       "Feature_Relative" = feature_rel,
                       "Prevalence" = prevalence,
                       "Prevalence_Relative" = prevalence_rel)
  colnames(df_out)[2] <- feature_tax
  return(df_out)
}

#-------------------------------------------------------------------------------------------------------------------

#' Good's coverage index
#'
#' The 'good_coverage' function estimates the Good's coverage index
#' based on a matrix of OTUs/ASVs/Species counts.
#' @param x numeric (of integers) matrix where rows are samples/observations and
#' columns are OTUs/ASVs/Species. If given a data frame or phyloseq-class object,
#' it will be coerced to a matrix.
#' @return A named (with samples) vector of Good's coverage estimates per sample.
#' @export

good_coverage <- function(x) {

  # implementation based on: (http://scikit-bio.org/docs/0.1.1/generated/skbio.math.diversity.alpha.gini_index.html#skbio.math.diversity.alpha.gini_index)
  # packages
  require("phyloseq")

  # check input
  if ( !(class(x) %in% c("matrix", "data.frame", "phyloseq")) ) { # class
    stop(paste("You provide an object of ", class(x),
               " class instead of a matrix, data.frame or phyloseq!"))
  }
  if ( is.data.frame(x) ) { # convert df into mtx
    x <-  as(x, "matrix")
  }
  if ( class(x) == "phyloseq" ) { # convert physeq into mtx
    if ( taxa_are_rows(x) ) { # if 'taxa_are_rows' tranpose
      x <- t(as(otu_table(x), "matrix"))
    } else {
      x <- as(otu_table(x), "matrix")
    }
  }
  #checkInteger <- apply(x, 2, function(x) is.integer(x)) # check if mtx 'x' are intergers
  checkInteger <- (x%%1==0) # to deal with double like, e.g., 1 and not as integer (coerce 1L)
  if ( !all(checkInteger == TRUE) ) stop(paste0(deparse(substitute(x)), " ", class(x), " is not integer!"))

  # calculate index
  singletons <- rowSums(x==1) # the no. of singletons OTUs/ASVs/Species x sample
  N <- rowSums(x) # the total no. of reads/individuals x sample
  GC_index <- 1 - (singletons / N)
  return(GC_index)
}

#---------------------------------------------------------------------------------------------------------------

#' Uquiness of features (OTUs/ASVs)
#'
#' Given a phyloseq object and a taxonomic rank (character) it gives the
#' uniquess of features (OTUs/ASVs), i.e., the number
#' of features (OTUs/ASVs) observed in one sample only.
#' @param physeq phyloseq-class object.
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @return A data frame with Taxa (lists the unique taxa that belongs to the
#' 'tax_rank' provided), Unique_Taxa (no. of features - OTUs/ASVs - that appear
#' only once, in one sample), Total_Taxa (no. of distinct features),
#' Percentage_unique _Taxa (percentage of Unique_Taxa across the Total_Taxa).
#' @export

uni_taxa_rank <- function(physeq, tax_rank) {

  # packages
  require("dplyr")
  require("phyloseq")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(deparse(substitute(physeq)), " is not a 'phyloseq' object!"))
  if( !is.character(tax_rank) | !(tax_rank %in% colnames(tax_table(physeq))) ) {
    stop(paste0(tax_rank, " is not a 'character' or/and not a taxonomic rank!"))
  }

  # metl physeq obj
  df <- melt_physeq(physeq) # tidy data.frame with physeq data

  # distinct taxa inside the input df
  taxa <- unique(df[,tax_rank]) %>%
    as.character(.)
  df_out <- data.frame("Taxa" = taxa,
                       "Unique_Taxa" = rep(0, length(taxa)),
                       "Total_Taxa" = rep(0, length(taxa)),
                       "Percentage_unique_Taxa" = rep(0, length(taxa)))

  for (i in 1:length(taxa)) {
    # retrieve a df by distinct taxa
    tax_level_df <- df %>%
      filter( (get(tax_rank) == as.vector(df_out$Taxa)[i]) & (Abundance != 0) )

    # count the no. of times that a taxa appear
    count_tax_level_df <- tax_level_df %>%
      group_by(OTU) %>%
      count(OTU)

    # df with the unique ASVs
    count_unique_taxa <- count_tax_level_df %>%
      filter(n == 1)

    ## add the output to the df_out
    df_out[i,"Unique_Taxa"] <- count_unique_taxa %>%
      nrow(.) # add the no. of unique taxa per rank to the df_out
    df_out[i,"Total_Taxa"] <- count_tax_level_df %>%
      nrow(.) # add the no. of taxa per rank to the df_out
  }
  # add the no. of % taxa per rank to the df_out
  df_out[,"Percentage_unique_Taxa"] <- df_out$Unique_Taxa / df_out$Total_Taxa * 100
  return(df_out)
}

#--------------------------------------------------------------------------------------------------------------

#' Uquiness of features (OTUs/ASVs) per sample
#'
#' Given a phyloseq object and a taxonomic rank (character), it gives the
#' uniquess of features (OTUs/ASVs), i.e., the number of features
#' (OTUs/ASVs) observed in one sample only, per sample per taxa.
#' @param physeq phyloseq-class object.
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @return A data frame with Taxa (lists the taxa that belongs to the
#' 'tax_rank' provided), Sample (lists the samples), Count (no. of unique
#' features belonging to that taxa - 'Taxa' - and found in only that sample - 'Sample').
#' @export

uni_taxa_per_sample <- function(physeq, tax_rank) {

  # packages
  require("dplyr")
  require("phyloseq")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(deparse(substitute(physeq)), " is not a 'phyloseq' object!"))
  if( !is.character(tax_rank) | !(tax_rank %in% colnames(tax_table(physeq))) ) {
    stop(paste0(tax_rank, " is not a 'character' or/and not a taxonomic rank!"))
  }

  # melt physeq obj
  df <- melt_physeq(physeq) # tidy data.frame with physeq data

  # retrieve uniq. samples and taxa
  sample <- unique(df$Sample) %>%
    as.character(.) # unique samples
  taxa <- unique(df[, tax_rank, drop = TRUE]) %>%
    as.character(.) # unique taxa

  # df to save out
  nr <- rep(length(taxa) * length(sample)) # size of nrows of df_out
  df_out <- data.frame("Taxa" = rep(0, nr),
                       "Sample" = rep(0, nr),
                       "Count" = rep(0, nr))
  counter <- seq(0, nr - length(sample), by = length(sample))
  for (i in 1:length(taxa)) { # loop over taxa
    taxa_df <- df %>%
      filter( get(tax_rank) == taxa[i] & Abundance!=0) # retrieve a df by distinct phylum
    count_taxa_df <- taxa_df %>%
      group_by(get(tax_rank)) %>%
      count(OTU) # count the no. of times that a feature appear across samples
    count_unique_taxa_df <- count_taxa_df %>%
      filter(n == 1) # df with the unique ASVs
    unique_taxa_names <- count_unique_taxa_df$OTU
    counter_x <- counter[i]
    for (j in 1:length(sample))  { # loop over samples
      x <- counter_x + j # move to the next taxa jumping from 0 to n_sample and so forth
      df_out[x,"Taxa"] <- taxa[i]
      df_out[x,"Sample"] <- sample[j]
      df_out[x,"Count"] <- taxa_df[taxa_df$Sample == sample[j] &
                                                 taxa_df$OTU %in% unique_taxa_names, ] %>%
        nrow() # add the count of unique taxa per taxa rank per sample to the df_out
    }
  }
  return(df_out)
}

#----------------------------------------------------------------------------------------------------------------

