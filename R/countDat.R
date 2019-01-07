require(edgeR)

#' An S4 class that stores count data.
#' @slot data Contains data.
#' @slot RPKM Contains RPKM data.
#' @slot annotation Contains annotation column data.
#' @slot cellObservables Contains cellObservables data.
#' @slot densityFunction Contains densityFunction data.
#' @slot estProps Contains estProps data.
#' @slot groups Contains groups data.
#' @slot nullPosts Contains nullPosts data.
#' @slot orderings Contains orderings data.
#' @slot posteriors Contains posteriors data.
#' @slot priorModels Contains priorModels data.
#' @slot priorType Contains priorType data.
#' @slot priors Contains priors data.
#' @slot replicates Contains replicates data.
#' @slot rowObservables Contains rowObservables data.
#' @slot sampleObservables Contains sampleObservables data.
#' @export

setClass("countDat", representation(data = "array", RPKM = "array", replicates = "factor", groups = "list", rowObservables = "list", sampleObservables = "list", cellObservables = "list", annotation = "data.frame", priorModels = "list", priorType = "character", densityFunction = "list", priors = "list", posteriors = "matrix", nullPosts = "matrix", estProps = "numeric", orderings = "data.frame" ))

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods
#' @param x Value
#' @param value Value

setGeneric("libsizes<-", function(x, value) standardGeneric("libsizes<-"))

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods

setMethod("libsizes<-", signature = "countDat", function(x, value) {
  x@sampleObservables$libsizes <- value
  x
})

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods

setGeneric("libsizes", function(x) standardGeneric("libsizes"))

#' libsizes method for testClass
#'
#' @docType methods
#' @rdname libsizes-methods

setMethod("libsizes", signature = "countDat", function(x) {
  x@sampleObservables$libsizes
})

###############################

#' getLibsizes method for testClass
#'
#' @docType methods
#' @rdname getLibsizes-methods
#' @param cD Value
#' @param subset Value
#' @param estimationType Value
#' @param quantile Value
#' @param ... Value

'getLibsizes2' <- function(cD, subset = NULL, estimationType = c("quantile", "total", "edgeR"), quantile = 0.75, ...)
{

  data<-cD@data 				# NEED THIS
  replicates <- cD@replicates		# NEED THIS

  if(missing(subset)) subset <- NULL
  if(is.null(subset)) subset <- 1:nrow(data)
  estimationType = match.arg(estimationType)
  if(is.na(estimationType)) stop("'estimationType' not known")

  estLibs <- function(data, replicates)
  {
    libsizes <- switch(estimationType,
                       total = colSums(data[subset,,drop = FALSE], na.rm = TRUE),
                       quantile = apply(data[subset,, drop = FALSE], 2, function(z) {
                         x <- z[z > 0]
                         sum(x[x <= quantile(x, quantile, na.rm = TRUE)], na.rm = TRUE)
                       }),

                       edgeR = {
                         if(!("edgeR" %in% loadedNamespaces()))
                           requireNamespace("edgeR", quietly = TRUE)
                         d <- edgeR::DGEList(counts = data[subset,, drop = FALSE], group = replicates, lib.size = colSums(data, na.rm = TRUE))
                         d <- edgeR::calcNormFactors(d, ...)
                         d$samples$norm.factors * d$samples$lib.size
                       })
    names(libsizes) <- colnames(data)
    libsizes
  }

  if(length(dim(data)) == 2) estLibsizes <- estLibs(data, replicates)
  if(length(dim(data)) == 3) {
    combData <- do.call("cbind", lapply(1:dim(data)[3], function(kk) data[,,kk]))
    combReps <- paste(as.character(rep(replicates, dim(data)[3])), rep(c("a", "b"), each = ncol(data)), sep = "")
    estLibsizes <- estLibs(combData, combReps)
    estLibsizes <- do.call("cbind",
                           split(estLibsizes, cut(1:length(estLibsizes), breaks = dim(data)[3], labels = FALSE)))
  }

  if(!missing(cD))
    if(inherits(cD, what = "pairedData")) return(list(estLibsizes[1:ncol(cD)], estLibsizes[1:ncol(cD) + ncol(cD)]))

  if(length(dim(data)) > 2) estLibsizes <- array(estLibsizes, dim = dim(cD@data)[-1])

  return(estLibsizes)
}

###############################

#' summary method for testClass
#'
#' @docType methods
#' @rdname extract-methods
#' @param x Value
#' @param i Value
#' @param j Value
#' @param ... Value
#' @param drop Value

setMethod("[", "countDat", function(x, i, j, ..., drop = FALSE) {
  if(missing(j)) {
    j <- 1:ncol(x@data)
  } else {
    if(is.logical(j)) j <- which(j)

    if(!all(1:ncol(x@data) %in% j))
    {
      replicates(x) <- as.character(x@replicates[j])

      if(length(x@groups) > 0)
      {
        newgroups <- list()
        newgroups <- lapply(x@groups, function(x) {
          x[j]
          rep(1:length(unique(x[j])), sapply(unique(x[j]), function(z) sum(x[j] == z)))[unlist(sapply(unique(x[j]), function(z) which(x[j] == z)))]
        })
        x@groups <- newgroups[!duplicated(newgroups) | duplicated(x@groups)]
      }

      if(length(x@posteriors) > 0)
      {
        warning("Selection of samples (columns) will invalidate the values calculated in slot 'posteriors', and so these will be discarded.")
        x@posteriors <- matrix(nrow = 0, ncol = 0)
      }
      if(length(x@orderings) > 0)
      {
        warning("Selection of samples (columns) will invalidate the values calculated in slot 'orderings', and so these will be discarded.")
        x@orderings <- data.frame()
      }

    }
  }

  if(missing(i))
    i <- 1:nrow(x@data)
  if(is.logical(i)) i <- which(i)

  if(nrow(x@data) > 0)
    x@data <- .sliceArray2(list(i, j), x@data)
  x@RPKM <- .sliceArray2(list(i, j), x@RPKM)

  x@annotation <- x@annotation[i,, drop = FALSE]
  if(nrow(x@posteriors) > 0)
    x@posteriors <- x@posteriors[i,, drop = FALSE]
  if(nrow(x@orderings) > 0)
    x@orderings <- x@orderings[i,, drop = FALSE]
  if(length(x@nullPosts) > 0)
    x@nullPosts <- x@nullPosts[i,,drop = FALSE]

  x@rowObservables <- lapply(x@rowObservables, function(z) .sliceArray2(list(i),z, drop = FALSE))
  x@sampleObservables <- lapply(x@sampleObservables, function(z) .sliceArray2(list(j), z, drop = FALSE))
  x@cellObservables <- lapply(x@cellObservables, function(z) .sliceArray2(list(i,j), z, drop = FALSE))

  x
})

###############################

.sliceArray2 <- function(slices, array, drop = FALSE) {
  if((is.vector(array) & sum(!sapply(slices, is.null)) > 1) || (is.array(array) & length(slices) > length(dim(array)))) warning("dimensions of slice exceed dimensions of array")
  sarray <- abind::asub(array, slices, dims = 1:length(slices), drop = drop)
  sarray
}

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods
#' @param x Value
#' @param value Value

setGeneric("replicates<-", function(x, value) standardGeneric("replicates<-"))

#' replicates method for testClass
#'
#' @docType methods
#' @rdname replicates-methods

setMethod("replicates<-", signature = "countDat", function(x, value) {
  x@replicates <- as.factor(value)
  x
})

