###------------------------------------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------------------------------------###

#' correctJunctionCoordinate
#'
#' @description Internal function for \code{\link{getJunctionCountMatrix}}, not to be run separately.
#'
#'
#' @return
#'

.correctJunctionCoordinate<- function(bed){
  vec<- strsplit(as.vector(bed[,11]),',')
  vec<- do.call('rbind',vec)
  bed[,2] <- bed[,2] + as.numeric(vec[,1])
  bed[,3] <- bed[,3] - as.numeric(vec[,2])
  return(bed)
}


#' Generate junction count matrix
#'
#' @description This function combines tophat2 pipeline output junctions.bed files after mapping reads to genome/trancriptome. It can be called separately for combining junction.bed files for the FASE pipeline.
#'
#' @param files junction bed files.
#'
#' @import BioPhysConnectoR
#' @return Junction read counts matrix.
#'
#' @references \enumerate{
#' \item  F. Hoffgaard, P. Weil, K. Hamacher. BioPhysConnectoR: Connecting Sequence Information and Biophysical Models. BMC Bioinformatics volume 11, Article number: 199 (2010).
#' }
#' @export
#'

getJunctionCountMatrix<- function(files) {
  # library(BioPhysConnectoR)
  tstamp <- paste("[",system("date",intern=TRUE),"]",sep="",collapse="")
  cat(paste(tstamp," Reading bed files...",sep="",collapse=""))
  #require(BioPhysConnectoR)
  files_C <- paste(files, collapse = ",")
  if (nchar(files_C) == 0)
    stop("No read files provided!")
  #Reading Junction.BED files
  junction_f <- lapply(1:length(files),function(x) read.csv(files[x], sep='\t',header=FALSE,quote=""))
  junction_f <- lapply(1:length(files),function(x) junction_f[[x]][-1,,drop=FALSE])

  cat('done','\n',sep='')
  tstamp <- paste("[",system("date",intern=TRUE),"]",sep="",collapse="")
  cat(paste(tstamp," Correcting Junction coordinates...",sep="",collapse=""))
  ##Correcting the overhung position of junction to get the correct corrdinate of the Junctions
  junction_f <- lapply(1:length(files),function(x) .correctJunctionCoordinate(junction_f[[x]]))

  cat('done','\n',sep='')
  tstamp <- paste("[",system("date",intern=TRUE),"]",sep="",collapse="")
  cat(paste(tstamp," Running...",sep="",collapse=""))

  AJ<- do.call('rbind',junction_f)
  AJ<- mat.sort(AJ,c(1,2))

  AJpos<-  paste(AJ[,1],AJ[,2],AJ[,3],AJ[,6],sep=":")
  AJposU <- unique(AJpos)
  UJnames <- paste(rep('JUNC',length(AJposU)),seq(1,length(AJposU),1),sep='')
  AJpos<- strsplit(AJposU,':')
  AJpos<- do.call('rbind', AJpos)
  AJpos<- data.frame(as.vector(AJpos[,1]),as.numeric(AJpos[,2]),as.numeric(AJpos[,3]),as.vector(AJpos[,4]),UJnames, matrix(0,nrow=length(AJposU), ncol=length(files)))

  junction_pos <- lapply(1:length(files),function(x) paste(junction_f[[x]][,1], junction_f[[x]][,2], junction_f[[x]][,3],junction_f[[x]][,6],sep=":"))
  junction_index <- lapply(1:length(files),function(x) match(junction_pos[[x]],AJposU))
  junction_f<-  lapply(1:length(files),function(x) junction_f[[x]][,5])
  for (i in 1: length(files))  AJpos[(junction_index[[i]]),(5+i)] <- junction_f[[i]]
  #Make seperate vector
  AJposT <-  lapply(6:ncol(AJpos),function(x) as.numeric(AJpos[,x]))
  AJposT<- do.call('cbind',AJposT)
  AJpos<- data.frame(as.vector(AJpos[,1]),as.numeric(AJpos[,2]),as.numeric(AJpos[,3]),as.vector(AJpos[,4]),as.vector(AJpos[,5]),AJposT)
  colnames(AJpos)<- c('chr','start','end','strand','UniqueJun',paste(rep('sample',length(files)),1:length(files),sep=''))

  cat('done','\n',sep='')
  return(AJpos)
}
