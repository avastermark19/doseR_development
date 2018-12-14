#' @title plotExpr Function to generate a boxplot (expression) for a countData object based on the replicate data.
#' @description This function boxplots the expression of the supplied countData object. plotExpr is a generic function used to produce boxplot summaries of various expression data. The function invokes particular plottings which depend on the selected grouping and treatment.
#' @usage plotExpr(cD, groupings, mode_mean, treatment, LOG2, clusterby_grouping, ...)
#' @param cD A countData object, an object for which a summary is desired.
#' @param groupings A grouping (annotation column), e.g. groupings="something".
#' @param mode_mean Boolean, Calculate RowMeans or RowMedians.
#' @param treatment list, indicating which treatments (columns) to be plotted.
#' @param LOG2 Boolean, Calculate LOG2.
#' @param clusterby_grouping	Boolean, Cluster boxplots by grouping (or treatment).
#' @param ...	Passthrough arguments to boxplot (additional arguments affecting the summary produced).
#' @details This function boxplots the expression of the supplied countData object, based on a selected annotation column and selected treatments.
#' @return Returns an invisible data frame containing values and labels.
#' @author AJ Vaestermark, JR Walters.
#' @references The "doseR" package, 2018 (in press).

## Factorized anntoation column input:
## cD@annotation$something <- factor(x = cD@annotation$something, levels = c("X", "e"))

plotExpr <- function (cD, groupings= NULL, mode_mean=TRUE, treatment=levels(cD@replicates), LOG2=TRUE, clusterby_grouping=TRUE, ...) {

  MyGroups<-cD@annotation[[groupings]]

  if(is.null(groupings)) {
    stop ('No groupings, e.g. groupings="something"...')
    return (NULL)
  }

  if(
    is.element( FALSE, treatment %in% levels(cD@replicates) )
  ) {
    stop ('Some treatment not in levels(cD@replicates), please check...')
    return (NULL)
  }

#  require(matrixStats)

  MyLabels<-NULL
  PLOT <- NULL
  NAMES <- NULL

  if ( is.factor(MyGroups) ) { MyGroups <- droplevels(MyGroups) }
  if ( is.factor(cD@replicates) )               { cD@replicates <- droplevels(cD@replicates) }

  Super_ch<- if(clusterby_grouping) unique(MyGroups) else treatment
  Super_dh<- if(clusterby_grouping) treatment else unique(MyGroups)

  for (ch in Super_ch) {
    for (dh in Super_dh) {

      actual_ch <- if(clusterby_grouping) dh else ch
      actual_dh <- if(clusterby_grouping) ch else dh

      if( is.element(actual_ch, treatment )) {

        NAMES <- c( NAMES, paste0(actual_ch, actual_dh))

        RM<- if(mode_mean) rowMeans(cD@RPKM[,cD@replicates==actual_ch]) else matrixStats::rowMedians(cD@RPKM[,cD@replicates==actual_ch])

        PLOT <-

          c( PLOT, RM[MyGroups == actual_dh ] )

        MyLabels<- c( MyLabels,                 rep(paste0(actual_ch, actual_dh),length(  RM[MyGroups == actual_dh ] ))        )

      }}}

  if(LOG2) { PLOT <-  log2(PLOT) }
  PLOT[is.infinite(PLOT)] <- NA
  MyLabels<- factor(MyLabels, levels=NAMES)

  boxplot(PLOT~MyLabels, ...)

  invisible( data.frame(values=PLOT,labels=MyLabels) )

}# plotExpr
