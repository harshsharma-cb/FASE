#' Intron Membership Matrix
#'
#' @description iMM describes association of each intron with meta-features (exons, skipping junctions and flanking junctions) of that gene. It can be generated using a gtf file and a combined junction matrix generated via \code{\link{getJunctionCountMatrix}}. iMM is a pre-requisite matrix for running \code{\link{iPrnaseq}}. It should be run only after running \code{\link{readMembershipMatrix}}.
#'
#' @import BioPhysConnectoR
#'
#' @param verbose TRUE
#' @param annotation matrix; contains annotation of exons and introns, created using \code{\link{readMembershipMatrix}}.
#'
#' @return intronMembershipMatrix creates gene-wise list which is saved by default as iMM.Rdata. Each gene is represented by a matrix of meta-features times the number of introns in gene. \cr
#'         A number is assigned for each meta-feature association to introns in the gene as: \cr
#'         \itemize{
#'                    \item 0  : No association
#'                    \item 1  : Exon associated with the intron
#'                    \item 2  : Intron with itself
#'                    \item 3  : Junction associated with the intron
#'
#'         }
#'
#' @references \enumerate{
#' \item F. Hoffgaard, P. Weil, K. Hamacher. BioPhysConnectoR: Connecting Sequence Information and Biophysical Models. BMC Bioinformatics volume 11, Article number: 199 (2010).
#' }
#' @export
#'

intronMembershipMatrix <- function(verbose=TRUE, annotation) {

    # library(BioPhysConnectoR)

    tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
    cat(paste(tstamp," Preprocessing...",sep="",collapse=""))
    #load('Annotation.Rdata')
    genes<-as.vector(annotation[,6])

    #For the moment i remove the junction that are not mapped to any genes but later i would think why
    iRetain<- !is.na(genes)
    genes <- genes[iRetain]
    annotation <- annotation[iRetain,,drop=FALSE]
    sIndex<- sort(genes,index.return=TRUE)$ix
    annotation <- annotation[sIndex,,drop=FALSE]
    genes<- genes[sIndex]
    ugenes<- levels(factor(genes))
    lgene<- length(ugenes)
    cat('done','\n',sep='')
    rm(sIndex)



    tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
    cat(paste(tstamp," Parsing intron Structure...",sep="",collapse=""))
    if(lgene < 5000) ntimes <-  seq(1, lgene, (lgene-1)/2) else ntimes <- seq(5000, lgene, 5000)
    iMM<- lapply(1:lgene, function(i) {
        #RMM<- lapply(1:10, function(i) {
        #i<-2
        if(any(ntimes==i)) cat(i,' ',sep='');
        # cat(i,'\n ',sep='');
        index<- match(genes,ugenes[i])
        sGTF<- annotation[!is.na(index),,drop=FALSE]
        if(nrow(sGTF) < 2) return(NA)
        sGTF<- mat.sort(sGTF,c(1,4))
        G<-  .iStructure(sGTF)
        return(G)

    })
    cat('\n',sep='')
    tstamp <- paste("[",Sys.time(),"]",sep="",collapse="")
    cat(paste(tstamp," Parsing intron Structure...Done","\n",sep="",collapse=""))
    names(iMM) <- ugenes #incase need to remove the flag0 portion
    save(iMM,file='iMM.Rdata') #now user save the data

    #flag0 ...................................................................
    #save(iMM, file='iMM.Rdata')
    return(iMM)
}



#' intron Structure
#'
#' @description Internal function for \code{\link{intronMembershipMatrix}}, not to be called separately.
#'
#'
#' @return
#'

.iStructure<- function(sGTF) {
    eijnames<- as.vector(sGTF[,5])
    indexExon<-grep('EX',eijnames)
    exGTF<- sGTF[indexExon,,drop=FALSE]
    indexIntron<-grep('IN',eijnames)
    inGTF<- sGTF[indexIntron,,drop=FALSE]
    indexJUN<- grep('JUN',eijnames)
    junGTF<- sGTF[indexJUN,,drop=FALSE]
    nxlen <- nrow(exGTF)
    nilen<- nrow(inGTF)
    njlen<- nrow(junGTF)
    #introns
    if(nxlen <2 | nilen< 1 ) return(NA)
    Gmatrix<- matrix(0, nrow=nrow(sGTF),ncol=nrow(inGTF))
    colnames(Gmatrix) <- inGTF[,5]
    rownames(Gmatrix) <- c(as.vector(exGTF[,5]), as.vector(inGTF[,5]), as.vector(junGTF[,5]))
    if (nilen !=1 ) diag(Gmatrix[(nxlen+1):(nxlen+nilen),]) <- 2 else Gmatrix[(nxlen+1):(nxlen+nilen),] <- 2
    #Exons
    if(nxlen> 0) {
        for (i in 1:nrow(inGTF)) {
            indexE<- ((exGTF[,3]+1)== inGTF[i,2] | (exGTF[,2]-1)== inGTF[i,3])
            Gmatrix[(1:nxlen)[indexE],i]<- 1
        }}
    ## Adding Skipping Junction
    if(nrow(junGTF) <1) return(Gmatrix)
    indSkipJ<- lapply(1:nrow(junGTF), function(j) (junGTF[j,2] <= inGTF[,2]) &  (junGTF[j,3] >= inGTF[,3]))
    for (i in 1:nrow(junGTF)) Gmatrix[(nxlen+nilen+i),indSkipJ[[i]]]<- 3

    return(Gmatrix)
    #End of function
}

