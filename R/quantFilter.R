#' @title quantFilter	Function to filter expression data of a countData object.
#' @description This function filters the expression of the supplied countData object; quantFilter is a filtering function used to remove rows (genes) of various expression data.
#' @usage quantFilter (cD, perc_cutoff=NULL, perc_cutoff2=0.05, MEDIAN=FALSE)
#' @param cD A countData object.
#' @param perc_cutoff	The lower cutoff, expressed as a percentage.
#' @param perc_cutoff2 The upper cutoff, expressed as a percentage.
#' @param MEDIAN Boolean, Calculate RowMeans or RowMedians.
#' @details This function filters the expression of the supplied countData object, based on a selected percentage cutoff.
#' @return Returns a filtered countData object.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

quantFilter <- function(cD, perc_cutoff=NULL, perc_cutoff2=0.05, MEDIAN=FALSE)
{
  cD@RPKM <- log2(cD@RPKM)
  ori_len<- nrow(cD@data)

  if(is.null(perc_cutoff)  ) {
    warning ('No filtering parameters specified... cancelling...')
    return (NULL)
  }

  if(length(cD@RPKM) == 0) {
    stop ('No RPKM data saved in count data object... cancelling...')
    return (NULL)
  }

  averages <- if(MEDIAN==FALSE) apply(cD@RPKM , 1, function(x) mean(x) ) else apply(cD@RPKM , 1, function(x) median(x) )
  averages <- averages[averages>0]

  mean_cutoff  <- as.numeric(quantile(averages,  perc_cutoff ))
  mean_cutoff2 <- as.numeric(quantile(averages,1-perc_cutoff2))

  cD <- if(MEDIAN==FALSE) cD[   apply(cD@RPKM , 1, function(x) mean(x) ) > mean_cutoff, ] else cD[   apply(cD@RPKM , 1, function(x) median(x) ) > mean_cutoff, ]
  cD <- if(MEDIAN==FALSE) cD[   apply(cD@RPKM , 1, function(x) mean(x) ) < mean_cutoff2,] else cD[   apply(cD@RPKM , 1, function(x) median(x) ) < mean_cutoff2,]

  if(!is.null(perc_cutoff) && length(cD@RPKM)>0 ) {

    percent <- function(x, digits = 2, format = "f", ...) { paste0(formatC(100 * x, format = format, digits = digits, ...), "%") }
    pct <- percent( 1 - ( nrow(cD@data) / ori_len )  )
    message( "Filtering removed ", ori_len-nrow(cD@data), " (", pct, ") of ", ori_len, " total loci." )

    cD@RPKM <- 2^(cD@RPKM)
    return(cD)

  }
}
