
#' Read Membership Matrix
#'
#' @description RMM describes association of each exon with meta-features (introns, skipping junctions and flanking junctions) of that gene. It can be generated using a gtf file and a combined junction matrix generated via \code{\link{getJunctionCountMatrix}}. RMM is a pre-requisite matrix for running \code{\link{EPrnaseq}}.
#'
#' @import BioPhysConnectoR
#'
#' @param gtf gtf file of the organism.
#' @param JunctionMatrix junction matrix contains read counts of each junction mapped by tophat2 alongwith their annotation.
#'
#' @return readMembershipMatrix creates a gene-wise list which is saved by default as RMM.Rdata. Each gene is represented by a matrix of meta-features times the number of exons in gene. \cr
#'         A number is assigned for each meta-feature association to exons in the gene as: \cr
#'         \itemize{
#'                    \item 0  : No association
#'                    \item 0.5: Skipping junction to the exon
#'                    \item 1  : Exon with itself
#'                    \item 2  : Flanking junction to the exon
#'                    \item 3  : Intron associated with the exon
#'
#'         }
#'
#' @references \enumerate{
#' \item F. Hoffgaard, P. Weil, K. Hamacher. BioPhysConnectoR: Connecting Sequence Information and Biophysical Models. BMC Bioinformatics volume 11, Article number: 199 (2010).
#' }
#' @export
#'

readMembershipMatrix <- function(gtf,JunctionMatrix=JunctionMatrix) {
  #library(BioPhysConnectoR)
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  cat(paste(tstamp," Reading GTF...",sep="",collapse=""))
  #   #gtf<- read.delim(gtf,sep='\t',header=FALSE)
  gtf<- read.csv(gtf, sep='\t',header=FALSE,quote="")
  # gtf<- read.csv('gene_corrected.gtf', sep='\t',header=FALSE,quote="")
  cat('done','\n',sep='')
  cat(paste(tstamp," Preprocessing...",sep="",collapse=""))

  ##Changing this to matrix creates problem for as.numeric(as.vector)
  #JunctionMatrix <- as.matrix(JunctionMatrix)

  #extract the genes names #extracted only exons and introns   # assigning unique exon and intron names
  indexExon<- match(as.vector(gtf[,3]),'exon')
  indexExon<- which(!is.na(indexExon))
  exGTF<- gtf[indexExon,]
  exGTF    <- exGTF[!duplicated(exGTF[,c(1,4,5)]),,drop=FALSE]
  exons <- paste('EX', seq(1,nrow(exGTF),1),sep='')
  exGTF<- cbind(exGTF,exons)
  indexIntron<- match(as.vector(gtf[,3]),'intron')
  indexIntron<- which(!is.na(indexIntron))
  inGTF<- gtf[indexIntron,]
  intron <- paste('INT', seq(1,length(indexIntron),1),sep='')
  inGTF <- cbind(inGTF,intron)
  colnames(exGTF) <- c(colnames(exGTF)[1:9],'EX_IN')
  colnames(inGTF) <- c(colnames(exGTF)[1:9],'EX_IN')
  gtf <- rbind(exGTF,inGTF)
  rm(exGTF,inGTF,indexIntron,indexExon); gc()

  #sort the gtf file
  genes<- as.character(gtf[,9])
  #Checking which is actually gene_id
  gene_id <- grep("gene_id", (strsplit(genes[1],';')[[1]]))
  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],';')[[1]][gene_id]))
  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],'"')[[1]][2]))

  #   #sort the gtf file
  #  genes<- as.character(gtf[,9])
  #  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],'"')[[1]][4]))
  #genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],'gene_id')[[1]][2]))
  sIndex<- sort(genes,index.return=TRUE)$ix
  gtf <- gtf[sIndex,,drop=FALSE]
  genes<- genes[sIndex]
  ugenes<- levels(factor(genes))
  lgene<- length(ugenes)
  cat('done','\n',sep='')
  rm(sIndex)

  #i<-1
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  cat(paste(tstamp," Parsing Gene Structure...",sep="",collapse=""))
  if(lgene < 5000) ntimes <-  seq(1, lgene, (lgene-1)/2) else ntimes <- seq(5000, lgene, 5000)
  RMM<- lapply(1:lgene, function(i) {
    #RMM<- lapply(1:100, function(i) {
    #i<-2
    if(any(ntimes==i)) cat(i,' ',sep='');
    # cat(i,'\n ',sep='');
    index<- match(genes,ugenes[i])
    sGTF<- gtf[!is.na(index),,drop=FALSE]
    if(nrow(sGTF) < 2) return(NA)
    sGTF<- mat.sort(sGTF,c(1,4))
    #Extract Junction reads
    fpmin <- min(as.vector(sGTF[,4]))
    fpmax <- max(as.vector(sGTF[,5])) #Changed to five since there could be a long exon defined but has some junction
    #indexJP<- (as.numeric(as.vector(JunctionMatrix[,2])) >= fpmin) & (as.numeric(as.vector(JunctionMatrix[,2])) < fpmax) & (as.vector(JunctionMatrix[,1])== as.vector(sGTF[1,1]))
    indexJP<- (JunctionMatrix[,2] >= fpmin) & (JunctionMatrix[,2] < fpmax) & (as.vector(JunctionMatrix[,1])== as.vector(sGTF[1,1])) #since i changed the code of JunctionMatrix
    #indexJP<- (as.numeric(JunctionMatrix[,2]) >= fpmin) & (as.numeric(JunctionMatrix[,2]) < fpmax) & (as.vector(JunctionMatrix[,1])== as.vector(sGTF[1,1]))
    sJM<- JunctionMatrix[indexJP,,drop=FALSE]
    #if(nrow(sJM) >0) cat('yes',ugenes[i],sep='')
    #image(t(.Gstructure(sGTF,sJM)))
    G<-  .Gstructure(sGTF,sJM)
    #flag0 return(G)
    Jun<- as.vector(sJM[,5])
    return(list(G=G,Jun=Jun))

  })
  cat('\n',sep='')
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  cat(paste(tstamp," Parsing Gene Structure...Done","\n",sep="",collapse=""))
  #names(RMM) <- ugenes #ncase need to remove the flag0 potion

  #flag0 .................................................................
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  cat(paste(tstamp," Post processing...","\n",sep="",collapse=""))
  Jnames <- lapply(RMM, function(x) x[2][[1]])
  RMM<- lapply(RMM, function(x) x[[1]])
  names(Jnames) <- ugenes
  names(RMM) <- ugenes
  jun<- JunctionMatrix[,c(1:5)]
  annotation <- gtf[,c(1,4,5,7,10)]
  rm(gtf, JunctionMatrix);gc();
  annotation<- cbind(annotation,genes)
  Jnames <- Jnames[!is.na(Jnames)]
  Jnames<- lapply(1: length(Jnames), function(x) cbind(Jnames[[x]],names(Jnames[x])))
  lenJnames<- lapply(Jnames,function(x) ncol(x))
  Jnames<- Jnames[lenJnames ==2]
  Jnames <- do.call('rbind',Jnames)
  index<- match(jun[,5],Jnames[,1])
  Jnames<- as.vector(Jnames[index,2])
  jun<- cbind(jun,Jnames)
  colnames(jun) <- colnames(annotation)
  annotation <- rbind(annotation,jun)
  save(annotation,file='Annotation.Rdata')
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  rm(annotation,jun);
  cat(paste(tstamp," Post processing...Done","\n",sep="",collapse=""))
  #flag0 ...................................................................
  return(RMM)
  save(RMM, file='RMM.Rdata')
}

#' Gene structure
#'
#' @description Internal function for \code{\link{readMembershipMatrix}}, not to be called separately.
#'
#'
#' @return
#'

.Gstructure <- function(sGTF,sJM) {

  indexExon<- match(as.vector(sGTF[,3]),'exon')
  exGTF<- sGTF[!is.na(indexExon),,drop=FALSE]
  indexIntron<- match(as.vector(sGTF[,3]),'intron')
  inGTF<- sGTF[!is.na(indexIntron),,drop=FALSE]
  nxlen <- nrow(exGTF)
  nilen<- nrow(inGTF)

  if(nxlen <2 ) return(NA)
  Gmatrix<- matrix(0, nrow=(nrow(sGTF)+nrow(sJM)),ncol=nrow(exGTF))
  colnames(Gmatrix) <- exGTF[,10]
  rownames(Gmatrix) <- c(as.vector(exGTF[,10]), as.vector(inGTF[,10]), as.vector(sJM[,5]))
  diag(Gmatrix[1:nrow(exGTF),]) <- 1

  #introns
  if(nilen> 0) {
    for (i in 1:nrow(inGTF)) {
      indexI<- ((exGTF[,5]+1)== inGTF[i,4] | (exGTF[,4]-1)== inGTF[i,5]	)
      Gmatrix[nxlen+i,indexI]<- 3
    }}

  if(nrow(sJM) <1) return(Gmatrix)
  #Junctions
  for (i in 1:nrow(sJM)) {
    #indexJ<-((exGTF[,5])== as.numeric(as.vector(sJM[i,2])) | (exGTF[,4]-1)== as.numeric(as.vector(sJM[i,3])))
    indexJ<-((exGTF[,5])== sJM[i,2] | (exGTF[,4]-1)== sJM[i,3]) #since i changed the code to gen JunctionMatrix
    #if(sum(indexJ)< 2) indexJ<- ((exGTF[,4]-5) <  as.numeric(as.vector(sJM[i,3])) & (exGTF[,5]+10) > as.numeric(as.vector(sJM[i,2])) )
    if(sum(indexJ)< 2) indexJ<- ((exGTF[,4]-5) < sJM[i,3] & (exGTF[,5]+10) > sJM[i,2] )
    Gmatrix[nxlen+nilen+i,indexJ]<- 2
  }

  #May be an extra check for the small exons or part of the exons that are not directly connected to Introns
  ## Adding Skipping Junction
  indSkipJ<- lapply(1:nrow(sJM), function(j) (sJM[j,2] < exGTF[,4]) &  (sJM[j,3] > exGTF[,5]))
  for (i in 1:nrow(sJM)) Gmatrix[(nxlen+nilen+i),indSkipJ[[i]]]<- 0.5

  return(Gmatrix)
  #End of function
}

##-------------------------------------------
## Update
#added support to generate skipping junction 21062016
##-------------------------------------------
