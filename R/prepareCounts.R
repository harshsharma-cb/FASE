#' Prepare counts for ExonPointer and IntronPointer
#'
#' @description Internal function for \code{\link{EPrnaseq}}/\code{\link{iPrnaseq}}, not to be run separately.
#'
#' @import matrixStats
#'
#' @return
#' @references Henrik Bengtsson (2017). matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors). R package version 0.52.2. https://github.com/HenrikBengtsson/matrixStats
#'
#'
#'
.prepareCounts<-function(y,designM,Groups=NULL,threshold=NULL, ...)
{
  #require(matrixStats)
  y <- as.matrix(y)
  if(is.null(threshold)) threshold <- 2 else threshold<- threshold
  if(is.null(Groups)) {
    if(any(rowSums(designM) >1))
      Groups <- rowSums(designM %*% t(designM))  else Groups<-which(designM!=0,arr.ind=TRUE)[,2]
  } else Groups <- Groups
  #-------------------------------
  Gfactor<- as.factor(Groups)
  n<- as.numeric(levels(Gfactor))
  GIndex <- match(Gfactor, n)
  #Filter the probes
  rs<- lapply(1:length(n),function(x) rowMedians(y[,(Gfactor==n[x]),drop=FALSE]))
  rs<-do.call('cbind',rs)
  newy <- y[(rowMaxs(rs) >threshold),,drop=FALSE]
  return(newy)
}
