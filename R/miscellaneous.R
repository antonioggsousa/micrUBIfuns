#' Get phyloseq tables
#'
#' Given a phyloseq object with 'otu_table', 'tax_table', 'sam_data' (required)
#' it gives a list with those tables.
#' @param physeq phyloseq-class object with 'otu_table', 'tax_table', 'sam_data' slots (required).
#' @param tax_rank taxonomic rank (character). One of the taxonomic ranks among
#' the column names of 'tax_table()' of the 'physeq' object given.
#' @param rm_na include (TRUE) or not (FALSE) NAs (logical), i.e., taxa without classification at
#' the taxonomic level specified at 'tax_rank'. It only works if 'tax_rank' was provided.
#' Otherwise it is ignored.
#' @param feature_w_tax substitute the feature labels (in columns) by the taxa names specified at
#' 'tax_rank' level (logical). Default is 'FALSE'. It only works if 'tax_rank' was provided.
#' Otherwise it is ignored.
#' @return A list with the 'feature' ('otu_table'), 'taxonomy' ('tax_table'), and 'metadata'
#' ('sample_data') tables. If the 'tax_rank' is specified, it will give 'feature' and 'taxonomy'
#' summarized at that 'tax_rank' level. If 'feature_w_tax' is equal 'TRUE' (only works together
#' with 'tax_rank'), the column names in the 'feature' table will be the taxa names at 'tax_rank'
#' level specied.
#' @export

get_physeq_tbls <- function(physeq, tax_rank = NULL,
                            rm_na = FALSE,
                            feature_w_tax = FALSE) {

  # Written: 15/09/2020
  # Last update: 16/09/2020

  # packages
  require("phyloseq")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(physeq, " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table", "sam_data") %in% slotNames(physeq) ) ) { # check the 3 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 3 slots: 'otu_table',
    'tax_table', 'sam_data'!\nThe 3 slots are necessary to use the 'melt_physeq' function.\nAborting..."))
  }
  if( !is.null(tax_rank) ) stopifnot( tax_rank %in% rank_names(physeq) )
  stopifnot( is.logical(feature_w_tax) )

  # get physeq by 'tax_rank'
  if ( !is.null(tax_rank) ) {
    physeq <- tax_glom(physeq = physeq,
                       taxrank = tax_rank,
                       NArm = rm_na)
  }

  # get data: otu, tax, sample tbls
  physeq_list <- list() # to save results
  if ( taxa_are_rows(physeq) ) {
    otu_tbl <- t(physeq@otu_table@.Data)
  } else {
    otu_tbl <- physeq@otu_table@.Data
  }
  physeq_list[["feature"]] <- otu_tbl
  physeq_list[["taxonomy"]] <- physeq@tax_table@.Data
  physeq_list[["metadata"]] <- data.frame(physeq@sam_data@.Data)
  colnames(physeq_list[["metadata"]]) <- physeq@sam_data@names
  rownames(physeq_list[["metadata"]]) <- physeq@sam_data@row.names

  # get 'feature_w_tax'
  if ( isTRUE(feature_w_tax) & !is.null(tax_rank) ) {
    if ( all(colnames(physeq_list[["feature"]]) == rownames(physeq_list[["taxonomy"]])) ) {
      colnames(physeq_list[["feature"]]) <- physeq_list[["taxonomy"]][,tax_rank]
    }
  }

  return(physeq_list)
}

#---------------------------------------------------------------------------------------------------------------------

#' Melt phyloseq object
#'
#' Given a phyloseq object with 'otu_table', 'tax_table', 'sam_data' (required)
#' it gives the data.frame melted with all the OTU/ASV (tidy), taxonomic and sample data.
#' Same functionality of the function 'psmelt()' from the 'phyloseq' R package, but
#' fastest as it uses 'data.table' R package for the task.
#' @param physeq phyloseq-class object with 'otu_table', 'tax_table', 'sam_data' slots (required).
#' @return A data frame with with all the OTU/ASV (tidy), taxonomic and sample data. The same
#' output to as 'psmelt()' function (from 'phyloseq').
#' @export

melt_physeq <- function(physeq) {

  require("phyloseq")
  library("data.table")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(physeq, " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table", "sam_data") %in% slotNames(physeq) ) ) { # check the 3 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 3 slots: 'otu_table',
    'tax_table', 'sam_data'!\nThe 3 slots are necessary to use the 'melt_physeq' function.\nAborting..."))
  }
  # combine: otu and tax tables
  if ( taxa_are_rows(physeq) ) { # check if
    #taxa are rows & ASV/OTU names match between otu_table() and tax_table()
    stopifnot( all( rownames(otu_table(physeq)) == rownames(tax_table(physeq)) ) )
    physeq_melt <- cbind.data.frame(physeq@otu_table@.Data, physeq@tax_table@.Data)
  } else { # transpose before cbind
    stopifnot( all( colnames(otu_table(physeq)) == rownames(tax_table(physeq)) ) )
    physeq_melt <- cbind.data.frame(t(physeq@otu_table@.Data), physeq@tax_table@.Data)
  }

  # melt and tidy: otu and tax tbls
  samples_names <- sample_names(physeq) # sample names
  tax_rank_names <- c(rank_names(physeq), "OTU") # tax rank names, inclusive, OTU/ASV ids
  physeq_melt$OTU <- rownames(physeq_melt) # save OTU (rownames) before convert to DT
  physeq_melt <- setDT(x = physeq_melt) # convert to DT
  physeq_melt <- melt(data = physeq_melt, id.vars = tax_rank_names,
                      measure.vars = samples_names, variable.name = "Sample",
                      value.name = "Abundance") # tidy read count only

  # add sample data to the previous 'physeq_melt'
  sample_dt <- setDT(x = physeq@sam_data@.Data)
  colnames(sample_dt) <- physeq@sam_data@names
  sample_dt[, Sample := physeq@sam_data@row.names]
  physeq_melt <- merge.data.table(x = physeq_melt, y = sample_dt, by = "Sample")

  # convert DT to DF
  physeq_melt <- setDF(x = physeq_melt)
  return(physeq_melt)
}
