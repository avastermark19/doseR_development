#' @title plotExpr Function to generate a boxplot (expression) for a countData object based on the replicate data.
#' @description This function generates a boxplot of FPKM expression values from the supplied countData object. FPKM values are averaged across replicates and partitioned among groups of loci as specified in a selected column from the annotation slot of the provided countData object.
#' @usage plotExpr(cD, groupings, mode_mean, treatment, LOG2, clusterby_grouping, ...)
#' @param cD A countData object containing FPKM values and at least one annotation column.
#' @param groupings	Specifies which column in the dataframe of the annotation slot that will be used to group loci in the boxplot. Can provide either a character value matching the column name, or a single numerical value used as an index of dataframe columns.
#' @param mode_mean Logical. If TRUE then FPKM values are averaged by mean across replicates within treatment. If FALSE, values are averaged by median.
#' @param treatment A character vector indicating which treatments (i.e. levels in the replicates slot vector) will be plotted. Order matters, and controls the ordering of treatments represented in the boxplot.
#' @param LOG2 Logical. If TRUE then average FPKM values are Log2 transformed.
#' @param clusterby_grouping Logical. If TRUE then boxplots are arranged by locus annotation grouping. If FALSE they are arranged by treatment levels, as indicated in the treatment argument.
#' @param ... Additional named arguments and graphical parameters passed to the boxplot function.
#' @details This function generates boxplots to visualize the distribution of FPKM expression values provided in a countData object, arranged by selected treatments and locus annotations. FPKM values are averaged (mean or median) within selected treatments, to provide a single expression value per locus per treatment. Loci are partitioned into groupings based on a specified column in the dataframe of annotations slot of the countData object. Thus a box is drawn for each grouping of loci for each treatment indicated. Desired treatments and their ordering are specified by the treatment argument. Groupings are arranged by sort order of the annotation column indicated, and can thus be controlled by providing a factor with a pre-specified level order. By default (clusterby_grouping = TRUE), boxes are arranged by annotation group first, and then by treatment, but setting this option to FALSE arranges boxes by treatment and then annotation group. This function uses the base graphics “boxplot” function to generate the plot, so can accept all relevant graphical arguments for customizing the figure; see “boxplot” for details.
#' @return Returns an invisible data frame containing values and labels used to generate the figure.
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
