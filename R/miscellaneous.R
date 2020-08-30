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
  require("data.table")

  # check input
  if(class(physeq) != "phyloseq") stop(paste0(physeq, " is not a 'phyloseq' object!"))
  if( !all( c("otu_table", "tax_table", "sam_data") %in% slotNames(physeq) ) ) { # check the 3 physeq slots
    stop(paste0(deparse(substitute(physeq)), " does not contain at least one of the 3 slots: 'otu_table',
    'tax_table', 'sam_data'!\nThe 3 slots are necessary to use the 'melt_physeq' function.\nAborting..."))
  }
  # combine: otu and tax tables
  if ( taxa_are_rows(physeq) &
       all( rownames(otu_table(physeq)) == rownames(tax_table(physeq)) ) ) { # check if
    #taxa are rows & ASV/OTU names match between otu_table() and tax_table()
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
