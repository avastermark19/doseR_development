#' @title dafsFilter
#' @description This function filters the expression.
#' @param VEC1 Vector 1.
#' @param PLOT Boolean, toggles plotting.
#' @details This function filters the expression.
#' @author AJ Vaestermark, JR Walters.
#' @references BMC Bioinformatics, 2014, 15:92

dafsFilter <- function(VEC1, PLOT) {

  #utils::globalVariables(c("mclust", "earth"))
  #utils::suppressForeignCheck("mclust")
  #requireNamespace(mclust)
 # require(earth)

  # remove 0 values
  #VEC1 <- VEC1[-which(VEC1==0)]
  VEC1 <- VEC1[VEC1!=0]

  #take log2 of data
  log2xx <- log2(VEC1)

  #vector to store Kolmogorov Smirnov distance statistics
  vv <- rep(0,0)
  vx <- rep(0,0)

  #select start point
  start <- length(log2xx[log2xx==min(log2xx)])/length(log2xx)

  #set sequence
  s <- seq(round(start,2),0.5,by=0.005)

  #loop through cuts of the data to determine targeted K-S statistic
  for(q in s) {

    #select data greater than a quantile and run Mclust on that data to determine theoretical distribution
    d <- log2xx[which(log2xx>quantile(log2xx,q,na.rm=T))]
    vx <- c(vx,quantile(log2xx,q,na.rm=T))
    out <- mclust::Mclust(d,G=1 , verbose=FALSE )
    ks <- suppressWarnings( ks.test(d,"pnorm",out$parameter$mean, sqrt(out$parameter$variance$sigmasq)) )
    vv <- c(vv,ks$statistic)

  }

  #determine first left-most local minima
  out <- earth::earth(s,vv,thresh=0.005)

  if(PLOT==TRUE) {
    plot(density(log2xx), ylim=c(0,0.5))
    lines(vx,vv,col="red")
  }

  cutv <- min(out$cuts[out$cuts>0])

  return( vx[which.min(vv)] )
} # dafsFilter
