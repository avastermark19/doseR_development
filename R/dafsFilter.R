#' @title dafsFilter	Function to filter expression data within a countData object.
#' @description This function filters the expression of the supplied countData object, by invoking the dafsFilter function. dafsFilter is a filtering function used to remove rows (genes) of various expression data.
#' @usage dafsFilter(cD, PLOT=TRUE)
#' @param cD A countData object.
#' @param PLOT Boolean, toggles plotting.
#' @details This function filters the expression of the supplied countData object using a Data Adaptive Flag filter. The internal function (dafs) uses a vector to store Kolmogorov Smirnov distance statistics, loops through cuts of the data to determine targeted K-S statistic, selects data greater than a quantile and runs Mclust on that data to determine theoretical distribution. The wrapper uses simpleFilter to determine first left-most local minima (using the Earth library).
#' @return Returns an invisible filtered countData object.
#' @author AJ Vaestermark, JR Walters.
#' @references BMC Bioinformatics, 2014, 15:92



dafsFilter <- function(cD, PLOT=TRUE) {

  #utils::globalVariables(dafsFilter)

  VEC1 <- rowMeans(cD@RPKM)

  cutoff <- dafs (VEC1, PLOT)

  invisible(   simpleFilter(cD, mean_cutoff=2^cutoff, counts=FALSE    )    )

}
