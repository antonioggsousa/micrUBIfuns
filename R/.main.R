## Define class 'micrUBI'
#
# TODO: make it functional
#
micrUBI <- setClass(Class = "micrUBI",
                    slots = list(physeq = "phyloseq",
                                 tax = "list",
                                 alpha_div = "list",
                                 beta_div = "list")
                    )

## Create micrUBI
#
# physeq
setGeneric("physeq", function(x) standardGeneric("physeq"))
setGeneric("physeq<-", function(x, value) standardGeneric("physeq<-"))
setMethod("physeq", "micrUBI", function(x) x@physeq)
setMethod("physeq<-", "micrUBI", function(x, value) {
  x@physeq <- value
  x
})

# tax
setGeneric("tax", function(x) standardGeneric("tax"))
setGeneric("tax<-", function(x, value) standardGeneric("tax<-"))
setMethod("tax", "micrUBI", function(x) x@tax)
setMethod("tax<-", "micrUBI", function(x, value) {
  x@tax <- value
  x
})

# alpha
setGeneric("alpha_div", function(x) standardGeneric("alpha_div"))
setGeneric("alpha_div<-", function(x, value) standardGeneric("alpha_div<-"))
setMethod("alpha_div", "micrUBI", function(x) x@alpha_div)
setMethod("alpha_div<-", "micrUBI", function(x, value) {
  x@alpha_div <- value
  x
})

# beta
setGeneric("beta_div", function(x) standardGeneric("beta_div"))
setGeneric("beta_div<-", function(x, value) standardGeneric("beta_div<-"))
setMethod("beta_div", "micrUBI", function(x) x@beta_div)
setMethod("beta_div<-", "micrUBI", function(x, value) {
  x@beta_div <- value
  x
})

#micrubi <- new(Class = "micrUBI", physeq = GlobalPatterns, tax = list(1, 2, 3))
