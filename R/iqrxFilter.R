#' @title iqrxFilter Function to filter expression data within a countData object.
#' @description This function filters the expression of the supplied countData object. iqrxFilter is a filtering function used to remove rows (genes) of various expression data.
#' @usage iqrxFilter(cD, perc_cutoff=0.25, perc_cutoff2=0.25, IQR_MULTIPLIER=1.5, MEDIAN=FALSE)
#' @param cD A countData object.
#' @param perc_cutoff	The lower cutoff, expressed as a percentage.
#' @param perc_cutoff2 The upper cutoff, expressed as a percentage.
#' @param IQR_MULTIPLIER Numeric multiplier; removes any outliers that are IQR_MULTIPLIER times the mid-50 percentile distance greater or less than the perc_cutoff and 1-perc_cutoff2 (25th and 75th percentiles, by default)
#' @param MEDIAN Boolean, Calculate RowMeans or RowMedians.
#' @details This function filters the expression of the supplied countData object, based on a selected percentage cutoff and selected interquartile range multiplier. The function iqrxFilter will: (1) log-base two transform all RPKM values (obligatory); (2) remove any outliers that were 1.5 times the mid-50 percentile distance greater or less than the 75th and 25th percentiles (by default), respectively; and (3) uses mean values and instead of median values (by default).
#' @return Returns a filtered countData object.
#' @author AJ Vaestermark, JR Walters.
#' @references Jue et al. BMC Genomics201314:150

iqrxFilter <- function(cD, perc_cutoff=0.25, perc_cutoff2=0.25, IQR_MULTIPLIER=1.5, MEDIAN=FALSE)
{
  cD@RPKM <- log2(cD@RPKM)
  ori_len<- nrow(cD@data)

  perc_cutoff <-perc_cutoff -((IQR_MULTIPLIER-1)/4)
  perc_cutoff2<-perc_cutoff2-((IQR_MULTIPLIER-1)/4)
  message( "new lb: ", perc_cutoff)
  message( "new ub: ", perc_cutoff2)

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

  if( length(cD@RPKM)>0 ) {

    percent <- function(x, digits = 2, format = "f", ...) { paste0(formatC(100 * x, format = format, digits = digits, ...), "%") }
    pct <- percent( 1 - ( nrow(cD@data) / ori_len )  )
    message( "Filtering removed ", ori_len-nrow(cD@data), " (", pct, ") of ", ori_len, " total loci." )

    cD@RPKM <- 2^(cD@RPKM)
    invisible(cD)

  }
} # iqrxFilter
