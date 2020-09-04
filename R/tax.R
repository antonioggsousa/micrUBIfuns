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
    'tax_table'!\nThe 2 slots are necessary to use the 'prevalence' function.\nAborting..."))
  }
  stopifnot( is.numeric(c(min_sample_depth, min_taxa_counts)) ) # check if input is numeric

  # check if otu table is integer
  checkInteger <- (physeq@otu_table@.Data%%1==0) # to deal with double like, e.g., 1 (by default - double)
  #and not as integer (coerce 1L)
  if ( !all(checkInteger == TRUE) ) stop(paste0("otu_table(", deparse(substitute(physeq)), ") is not integer!"))

  ## Filter samples
  #
  if ( !is.null(samples2filter) ) physeq <- prune_samples(samples = samples2filter, x = physeq)
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
