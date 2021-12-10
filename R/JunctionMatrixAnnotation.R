#' Junction Matrix Annotation
#'
#' @description Annotation of junction matrix using gtf file.
#'
#' @param gtf gtf file of the organism.
#' @param JunctionMatrix matrix containing junction read counts.
#'
#' @return Annotated junction matrix file: JunctionMatixA.
#' @export
#'

JunctionMatrixAnnotation<- function(gtf,JunctionMatrix) {
  gtf<- read.csv(gtf, sep='\t',header=FALSE,quote="")
  #Filter genes only
  genes<- as.character(gtf[,3])
  index<- match(genes,"gene")
  gtf<- gtf[!is.na(index),]

  genes<- as.character(gtf[,9])
  gene_id <- grep("gene_id", (strsplit(genes[1],';')[[1]]))
  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],';')[[1]][gene_id]))
  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],'gene_id')[[1]][2]))
  sIndex<- sort(genes,index.return=TRUE)$ix
  gtf <- gtf[sIndex,,drop=FALSE]
  genes<- genes[sIndex]
  gtf[,9] <- as.vector(genes)
  indexJP<- lapply(1:nrow(JunctionMatrix),function(x)
    which((as.vector(JunctionMatrix[x,2]) >= as.vector(gtf[,4])) & (as.vector(JunctionMatrix[x,2]) < as.vector(gtf[,5])) & (as.vector(JunctionMatrix[x,1])== as.vector(gtf[,1]))==1))

  newJP<- lapply(indexJP, function(x) x[1][[1]])
  newJP<- unlist(newJP)
  genes<- as.vector(gtf[newJP,9])
  JunctionMatrixA <- cbind(JunctionMatrix,genes)
  save(JunctionMatrixA,file='JunctionMatrixAnnotation.Rdata')

}
