#' Cassette exon/intron retention event ranking by contrast
#'
#' @description Takes the annotated and fitted object of \code{\link{EPrnaseq}}/\code{\link{iPrnaseq}} and the name or number of contrast as given in contrast matrix as input and finds differentially alternatively spliced cassette exon/intron retention events for that contrast.
#'
#' @param fit output of addAnnotationRnaSeq function.
#' @param contrast contrast whose ranking is required, for example, 'NormalvsTumor' (as used in contrast matrix).
#'
#' @return Data frame that contains ranking of cassette exon/intron retention events in the given contrast or comparision with their annotation. The output can be saved as csv/xlsx file.
#' @export
#'

getPvaluesByContrast <- function(fit, contrast=NULL)
{
  if (attr(fit, "name") != "EPrnaSeq" & attr(fit, "name") != "iPrnaSeq") stop('Not an object of EPrnaSeq or iPrnaSeq')
  atNames<- attr(fit, "name")
  if (attr(fit, "fun") != 'addAnnotationRnaSeq') stop('Please run first addAnnotationRnaSeq')
  avContrast <- levels(fit[,7])

  if (is.null(contrast)) {
    cat('Available contrast :', length(avContrast),': Using only first one','\n',sep='')
    index<- match(as.vector(fit[,7]), avContrast[1])
    fit <- fit[!is.na(index),,drop=FALSE]
  } else {

    if(is.numeric(contrast)) {
      if(contrast > length(avContrast)) stop('Please select from available contrasts in contrast matrix')
      index<- match(as.vector(fit[,7]), avContrast[contrast])
      fit <- fit[!is.na(index),,drop=FALSE]
    }

    if(is.character(contrast)){
      if(!any(avContrast==contrast)) stop('Please select from available contrasts in contrast matrix')
      index<- match(as.vector(fit[,7]), contrast)
      fit <- fit[!is.na(index),,drop=FALSE]
    }

  }

  return(fit)

}


#' Sum of p-values
#'
#' @description Internal function of \code{\link{EPrnaseq}}/\code{\link{iPrnaseq}}, not to be called separately.
#'
#' @param x
#'
#' @import matrixStats
#' @import stats
#' 
#'
#' @return
#'

.sumPvalsMethod <-  function(x,n) {

  if (n < 10) {

    psumunif(x,n)
  } else {

    pnorm(x,n/2,sqrt(n/12),lower=TRUE)
  }}

psumunif <-  function(x,n) 1/factorial(n) * sum(sapply(0:n, function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)))

prepareCounts<-function(y,designM,Groups=NULL,threshold=NULL, ...)
{

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

#' Remove low-expressed reads
#'
#' @description Internal function of \code{\link{EPrnaseq}}/\code{\link{iPrnaseq}}, not to be called separately.
#'
#' @import matrixStats
#'
#'
#' @return
#'
#' @references Henrik Bengtsson (2017). matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors). R package version 0.52.2. https://github.com/HenrikBengtsson/matrixStats
#'
#'

.removeLECounts<-function(y,designM,Groups=NULL,threshold=6.32, ...)
{

  y <- as.matrix(y)

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
