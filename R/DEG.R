#' Differentially Expressed Genes
#'
#' @description A wrapper function of limma package to find differentially expressed genes, given summarized read counts of genes obtained from \code{\link[Rsubread]{featureCounts}} function or preprocessing intructions.
#'
#' @param geneCount summarized read counts of genes.
#' @param designM design matrix required by limma.
#' @param contrastM contrast matrix required by limma.
#' @param Groups list of sample groups. \cr
#'        Example: If there are two sample groups with three samples each, 'Groups' should be formed as:
#'        \enumerate{
#'        \item numeric: c(1, 1, 1, 2, 2, 2)
#'        }
#'
#' @import edgeR
#'
#' @return Saves raw gene counts and log2cpm expression for all genes. Meta-data generated through this function is saved in fit2.Rdata file. Further, annotation of this meta-data is performed by \code{\link{addAnnotationDEG}} function. Contrast-wise ranking of annotated differentially expressed genes can be obtained using \code{\link{cpmCountsDEG}} function.
#'
#' @references \enumerate{
#' \item Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139–140 (2009)
#' }
#' @export
#'
DEG<- function(geneCount, designM, contrastM, Groups)
{

        #load('counts_genes.Rdata')
        counts<- geneCount$counts #saving raw counts for all genes
        write.csv(counts,file='RawGeneCounts.csv')
        #load('DCmatrix.Rdata')

        #removing genes with low read counts less than 3
        counts <- prepareCounts(counts, designM, Groups, threshold=2)
        ncounts <- DGEList(counts, Groups)
        ncounts <- calcNormFactors(ncounts)
        fit<- voom(ncounts,designM)

        E <- ((counts !=0)*1)* fit$E
        write.csv(E,file='RawGeneCounts_log2cpm.csv')

        fit <- lmFit(E, designM)
        fit2 <- contrasts.fit(fit, contrastM)
        fit2 <- eBayes(fit2)
        return(fit2)

}

#' Annotation of differentially expressed genes
#'
#' @description Add gene annotation to ranked differentially expressed genes for a given contrast, using output of \code{\link{DEG}} function.
#'
#' @param geneCount summarized read counts of genes.
#' @param fit output of \code{\link{DEG}} function that contains ranking of differentially expressed genes.
#' @param contrast a contrast from contrast matrix, whose ranking is required.
#'
#' @return Annotated ranking of differentially expressed genes of given contrast. The output can be saved using write.csv or \code{\link[openxlsx]{write.xlsx}}.
#' @export

addAnnotationDEG <- function(geneCount,fit,contrast)
{
        test<- topTable(fit, coef=contrast, n = Inf, sort = "p")
        #load('counts_genes.Rdata')
        annotation <- geneCount$annotation
        if(ncol(test) > 6) # When fit$coefficients has duplicate geneIDs
                index<- match(as.vector(test[,1]),as.vector(annotation[,1]))
        else
                index<- match(as.vector(rownames(test)),as.vector(annotation[,1]))
        test<- cbind(test,annotation[index,,drop=FALSE])
        return(test)
}


#' cpmCountsDEG
#'
#' @description Generates read counts and log2CPM expression for differentially expressed genes for a given contrast, using output of \code{\link{addAnnotationDEG}} function.
#'
#' @param geneCount summarized read counts of genes.
#' @param filename filename in which output of addAnnotationDEG is saved.
#' @param designM design matrix required by limma
#' @param contrastM contrast matrix required by limma.
#' @param Groups list of sample groups. \cr
#'        Example: If there are two sample groups with three samples each, 'Groups' should be formed as:
#'        \enumerate{
#'        \item numeric: c(1, 1, 1, 2, 2, 2)
#'        }
#'
#' @return Read counts and log2cpm expression of given contrast for ranked differentially expressed genes.
#'
#' @references \enumerate{
#' \item Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139–140 (2009)
#' }
#' @export
#'

cpmCountsDEG<- function(geneCount, filename, designM, contrastM, Groups)
{

        #load('DCmatrix.Rdata')
        test <- read.csv(filename)
        filename<- strsplit(filename,'.',fixed=TRUE)[[1]][1]

        eventName<- as.vector(test$X)

        # load('counts_genes.Rdata')
        counts<- geneCount$counts #saving raw counts for all genes
        # load('DCmatrix.Rdata')
        index <- match(eventName,rownames(counts))
        counts <- counts[index,,drop=FALSE]

        countfilename <- paste(filename,'_count.csv',sep="",collapse="")
        write.csv(counts,file=countfilename)

        ncounts <- DGEList(counts=counts, group=Groups)
        ncounts <- calcNormFactors(ncounts)
        fit<- voom(ncounts,designM)
        E <- ((counts !=0)*1)* fit$E
        l_filename <- paste(filename,'_log2cpm.csv',sep="",collapse="")
        write.csv(E,file=l_filename)
}
