#' extractIntrons
#'
#' @description Internal function of \code{\link{intronGTFparser}}, not to be called separately.
#'
#'
#' @return
#'
.extractIntrons<- function(sGTF) {


  if(nrow(sGTF) < 2) return(NA);
  vec<- paste(as.vector(sGTF[,4]),as.vector(sGTF[,5]))
  index<- which(duplicated(vec))
  if(length(index)!=0)  sGTF<- sGTF[-index,,drop=FALSE];

  #Finding the overlapping Exons
  indOverlap<- lapply(1:nrow(sGTF), function(x) which((sGTF[x,4] >=  sGTF[,4]) &  (sGTF[x,5] <= sGTF[,5])))
  indOverlap<- lapply(1:length(indOverlap),function(x) cbind(indOverlap[[x]],x))
  indOverlap <- do.call('rbind',indOverlap)

  idmatrix<- sparseMatrix(i=indOverlap[,1], j=indOverlap[,2],x=1)
  rnames<- levels(as.factor(indOverlap[,1]))
  colnames(idmatrix) <- rnames
  rownames(idmatrix) <- rnames

  idmatrix <- as(idmatrix,"TsparseMatrix")
  gR <- ftM2graphNEL(cbind(1+idmatrix@i, 1+idmatrix@j))
  Salida <- connectedComp(gR)

  sGTF<- lapply(1: length(Salida), function(x) {
    indEX<- as.numeric(Salida[[x]])
    return(cbind(min(sGTF[indEX,4]), max(sGTF[indEX,5]), as.vector(sGTF[indEX[1],10]))) })
  sGTF <- do.call('rbind',sGTF)
  if(nrow(sGTF) < 2) return(NA);
  iGTF<- cbind((as.numeric(sGTF[1:(nrow(sGTF)-1),2])+1), (as.numeric(sGTF[2:(nrow(sGTF)),1])-1), as.vector(sGTF[1:(nrow(sGTF)-1),3]))
  return(iGTF)
}


#' intronGTFparser
#'
#' @description Parse intron location given in a gtf file and updated gtf will be written. Intron information can be used then for counting reads with Rsubread package (check wrapper functions: \code{\link{ppAuto}} and \code{\link{ppSumEIG}} for read count summarization). However, information associated with these introns (related to transcripts) can not be used as annotation since this transcript information is added in the corresponding field to avoid unnecessary errors.
#'
#' @param gtf gtf file of the organism.
#'
#' @import Matrix
#' @import RBGL
#' @import graph
#'
#' @return gtf file with intron information.
#'
#' @references \enumerate{
#'  \item http://Matrix.R-forge.R-project.org/
#'  \item Carey V, Long L, Gentleman R (2019). RBGL: An interface to the BOOST graph library. https://bioconductor.org/packages/RBGL/
#'  \item Gentleman R, Whalen E, Huber W, Falcon S (2019). graph: graph: A package to handle graph data structures. http://www.bioconductor.org/packages/release/bioc/html/graph.html
#'  }
#'
#' @export
intronGTFparser <- function(gtf){
  #require(Matrix); require(RBGL);require(graph)

  #Parse for the intron between exons in the genome
  #remove duplicate exons for two or more transcript but the same Gene
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  cat(paste(tstamp," Reading GTF...",sep="",collapse=""))
  #   #gtf<- read.delim(gtf,sep='\t',header=FALSE)
  #   gtf<- read.csv('genes_TAIR10.gtf', sep='\t',header=FALSE,quote="")
  gtfName<- gtf
  gtfName<- strsplit(gtfName,".gtf")[[1]]
  gtf<- read.csv(gtf, sep='\t',header=FALSE,quote="")
  if (ncol(gtf)==1) stop('Please manually remove headers(comments) from GTF and rerun')
  indexExon<- match(as.vector(gtf[,3]),'exon')
  indexExon<- which(!is.na(indexExon))
  exGTF<- gtf[indexExon,]
  cat('done','\n',sep='')

  ##Parsing the introns
  eNumber<- paste(1:nrow(exGTF),'e',sep='')
  exGTF<-cbind(exGTF,eNumber)
  genes<- as.character(exGTF[,9])
  #Checking which is actually gene_id
  gene_id <- grep("gene_id", (strsplit(genes[1],';')[[1]]))  #incase of gff file please use gene_id <- grep("ID=", (strsplit(genes[1],';')[[1]]))
  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],';')[[1]][gene_id]))
  # genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],'ID=')[[1]][2]))  #run incase use are using the GFF file also the previous line
  # genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],'"')[[1]][2]))
  #  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],';')[[1]][gene_id]))
  ugenes<- levels(factor(genes))
  lgene<- length(ugenes)
  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  #cat(paste(tstamp," Parsing possible introns...",sep="",collapse=""),"\n")
  cat(paste(tstamp," Parsing possible introns...",sep="",collapse=""));
  if(lgene < 5000) ntimes <-  seq(1, lgene, (lgene-1)/2) else ntimes <- seq(5000, lgene, 5000)


  iGTF<- lapply(1:lgene, function(i) {
    #iGTF<- lapply(1:10, function(i) {
    if(any(ntimes==i)) cat(i,' ',sep='');
    # cat(i,'\n ',sep='');
    index<- match(genes,ugenes[i])
    sGTF<- exGTF[!is.na(index),,drop=FALSE]
    #in case of single exon retun NA
    if(nrow(sGTF) < 2) return(NA)
    if(length(levels(factor(sGTF[,7])))==1) { #what about + and -, would not be needing ## there would be problem : i realized in case of both in one genes
      iGTF<- .extractIntrons(sGTF)
      return(iGTF)

    } else {
      psGTF <- sGTF[(as.character(sGTF[,7]) == '+'),,drop=FALSE]
      nsGTF <- sGTF[(as.character(sGTF[,7]) == '-'),,drop=FALSE]
      piGTF<- .extractIntrons(psGTF)
      niGTF<- .extractIntrons(nsGTF)
      iGTF<- rbind(piGTF, niGTF)
      return(iGTF)
    }
    #i can devide both of them and then by merging
  } )

  cat('Done','\n',sep='')
  iGTF<-iGTF[(!is.na(iGTF))]
  iGTF<-do.call('rbind',iGTF)
  gc() #make some free memory
  index<-match(iGTF[,3], exGTF[,10])
  iGTF<- cbind( exGTF[index,1], exGTF[index,2],rep('intron', nrow(iGTF)), as.vector(iGTF[,1]),as.vector(iGTF[,2]),exGTF[index,6:9])
  iGTF<- iGTF[(as.numeric(as.vector(iGTF[,5])) - as.numeric(as.vector(iGTF[,4])) != -1),,drop=FALSE]
  index<- which(as.numeric(as.vector(iGTF[,5])) -as.numeric(as.vector(iGTF[,4])) < 0)
  iGTF<- iGTF[-index,,drop=FALSE]
  igtfName <- paste(gtfName,'_temp_igtf.Rdata',sep='')
  #save(iGTF,file=igtfName)
  colnames(iGTF)<- colnames(gtf)
  finalGTF <- rbind(gtf,iGTF)

  tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
  #cat(paste(tstamp," Parsing possible introns...",sep="",collapse=""),"\n")
  cat(paste(tstamp," Writing Output files...",sep="",collapse=""))

  ##sort by Genes
  genes<- as.character(finalGTF[,9])
  genes<- unlist(lapply(1:length(genes),function(x) strsplit(genes[x],';')[[1]][2]))
  sIndex<- sort(genes,index.return=TRUE)$ix
  finalGTF <- finalGTF[sIndex,,drop=FALSE]

  finalGTF<-data.frame(as.vector(finalGTF[,1]),as.vector(finalGTF[,2]),as.vector(finalGTF[,3]),
                       as.vector(finalGTF[,4]),as.vector(finalGTF[,5]),as.vector(finalGTF[,6]),
                       as.vector(finalGTF[,7]),as.vector(finalGTF[,8]),as.vector(finalGTF[,9]))

  gtfName <- paste(gtfName,'_corrected.gtf',sep='')
  ##Write GTF
  write.table(finalGTF,file=gtfName,sep='\t',quote = FALSE, row.names = FALSE,
              col.names = FALSE)

  #cat('gene_corrected.gtf'," have been written successfully",sep='')
  cat('Done','\n',sep='')
  return(gtfName)
}

