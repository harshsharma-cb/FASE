
#' Junction Matrix Annotation
#'
#' @description This function is used for correcting the annotation of junction matrix in case junction.bed files are obtained from "BAM" files using "regtools" software. This is the case where "BAM" files are downloaded some repository instead of running tophat2 pipeline. This function then produces a correctly annotated junction matrix on the basis of chromosome nomenclature as used in standard "GTF" file.\cr
#'
#'
#' @param gtf gtf file of the organism.
#' @param JunctionMatrix matrix with read counts of junctions.
#'
#' @return Annotated junction matrix file: JunctionMatix.
#' @export
#'
GTFnomencJunctionM<- function(gtf,JunctionMatrix){
  #Due to change readcounts done by regtools for Junction extracted from Bamfiles
  gtf<- read.csv(gtf, sep='\t',header=FALSE,quote="")
  gtfnames<- levels(gtf[,1])
  lchanges<- needtochanges <- levels(JunctionMatrix[,1])
  indexTochange<- grep('chr',needtochanges)
  chrnames<- strsplit(needtochanges[indexTochange],'chr')
  chrnames<- unlist(lapply(chrnames, function(x) x[2]))
  needtochanges[indexTochange] <- chrnames

  indexTochange<- grep('_',needtochanges)
  chrnames<- strsplit(needtochanges[indexTochange],'_')
  chrnames<- unlist(lapply(chrnames, function(x) x[2]))
  needtochanges[indexTochange] <- chrnames

  indexTochange<- grep('v',needtochanges)
  chrnames<- strsplit(needtochanges[indexTochange],'v')
  chrnames<- unlist(lapply(chrnames,function(x) paste(x,collapse = '.')))
  needtochanges[indexTochange] <- chrnames

  finalchange <- as.vector(JunctionMatrix[,1])
  indexTochange<- match(finalchange,lchanges)
  finalchange <- needtochanges[indexTochange]
  JunctionMatrix[,1]<- as.factor(finalchange)
  save(JunctionMatrix,file="JunctionCounts.Rdata")

}
