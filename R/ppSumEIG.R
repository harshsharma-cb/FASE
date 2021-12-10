#' RNA Seq Preprocessing Read Summarization
#'
#' @description ppSumEIG is a manual wrapper function that provides summarization of read counts for exons, introns and genes using \code{\link[Rsubread]{featureCounts}}. The gtf file passed to this function should first be passed to \code{\link{intronGTFparser}} to find the location of introns. Reads used for summarization by ppSumEIG should already be mapped, sorted and index using tophat2 and samtools or their wrapper function: \code{\link{ppRawData}}. The summarized counts produced by ppSumEIG can be further processed using \code{\link{ppFASE}}, which produces several matrices required by ExonPointer and IntronPointer algorithms for finding alternative splicing events. ppSumEIG need not be run if \code{\link{ppAuto}} has already been run.
#'
#' @param folderSRA path of directory containing aligned and indexed bam file folders. (default=current directory)
#' @param pairedend boolean, TRUE if reads are paired-end and FALSE if reads are single-end. (default=FALSE).
#' @param p number of threads to be utilized by Rsubread package. (default=1)
#' @param gtf intron parsed gtf file of the organism.
#' @param ... other parameters to be passed to \code{\link[Rsubread]{featureCounts}}.
#'
#' @import Rsubread
#'
#' @return Lists of gene counts, exon counts and intron counts saved in folderSRA directory as respective Rdata files. (Can be run separately using \code{\link[Rsubread]{featureCounts}} or automatic wrapper function for entire pre-processing: \code{\link{ppAuto}})
#'
#' @references \enumerate{
#' \item Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47, e47 (2019).
#' }
#' @export
#'

ppSumEIG <- function(folderSRA = FALSE, pairedend = FALSE, p=1, gtf = gtf, srlist = NULL, ...) {

  if(folderSRA == FALSE){
    folderSRA <- getwd();
  } else {folderSRA=folderSRA}

  start.time <- Sys.time()
  if(is.null(srlist)){
    load('srlist.Rdata')
  } else {
    srlist = srlist
    save(srlist, file = 'srlist.Rdata')
  }

  #load('srlist.Rdata')
  # if(files == 'sra'){
  # fastq_dump <- readline(prompt="Is fastq-dump installed in the system? y/n ")
  # if(fastq_dump == 'n'){
  # stop("Please install fastq-dump to continue.")
  # }
  # }
  #
  # tophat2 <- readline(prompt="Is tophat2 installed in the system? y/n ")
  # if(tophat2 == 'n'){
  # stop("Please install tophat2 to continue.")
  # }
  #
  # samtools <- readline(prompt="Is samtools installed in the system? y/n ")
  # if(samtools == 'n'){
  # stop("Please install samtools to continue.")
  # }

  GTF_corrected <- readline(prompt = "Has intronGTFparser() been run (to find introns in gtf)? y/n ")
  if(GTF_corrected == 'n'){
    stop('Please parse gtf using intronGTFparser() to continue and use _corrected.gtf for analysis.')
  }

  ### Summarization
  srlistbam <- paste(folderSRA,'/',srlist,'_tophat_out/accepted_hits_sorted.bam',sep='')


  ##Counting exon reads

  cat(paste('\n', "Task 1 out of 3:",'\n' ,"Summarization of reads for exons...", '\n' ))

  #require(Rsubread)
  tstamp <- Sys.time()
  cat(paste("[", tstamp,"]", " Counting exon reads... ", '\n', sep="", collapse=""))
  exonCount <- featureCounts(files= srlistbam, isPairedEnd = pairedend,requireBothEndsMapped=pairedend,GTF.featureType = "exon", GTF.attrType = "gene_id", useMetaFeatures = FALSE, isGTFAnnotationFile=TRUE, annot.ext=gtf,allowMultiOverlap=TRUE,nthreads=p)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Counting exon reads...done ", '\n', sep="", collapse=""))
  save(exonCount, file='counts_exons.Rdata'); rm(exonCount);gc();
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving exon reads count...done ", '\n', sep="", collapse=""))


  ##Counting intron reads

  cat(paste('\n', "Task 2 out of 3:",'\n' ,"Summarization of reads for introns...", '\n' ))

  tstamp <- Sys.time()
  cat(paste("[", tstamp,"]", " Counting intron reads... ", '\n', sep="", collapse=""))
  intronCount <-featureCounts(files=srlistbam, isPairedEnd=pairedend,requireBothEndsMapped=pairedend,GTF.featureType = "intron", GTF.attrType = "gene_id", useMetaFeatures = FALSE, isGTFAnnotationFile=TRUE, annot.ext=gtf,allowMultiOverlap=TRUE,nthreads=p)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Counting intron reads...done ", '\n', sep="", collapse=""))
  save(intronCount, file='counts_introns.Rdata');
  rm(intronCount);
  gc();
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving intron reads count...done ", '\n', sep="", collapse=""))



  ##Counting genes

  cat(paste('\n', "Task 3 out of 3:",'\n' ,"Summarization of reads for genes...", '\n' ))

  tstamp <- Sys.time()
  cat(paste("[", tstamp,"]", " Counting gene reads... ", '\n', sep="", collapse=""))
  geneCount <- featureCounts(files=srlistbam, isPairedEnd=pairedend, requireBothEndsMapped=pairedend,GTF.featureType = "gene", GTF.attrType = "gene_id", useMetaFeatures = FALSE, isGTFAnnotationFile=TRUE, annot.ext=gtf,allowMultiOverlap=TRUE,nthreads=p)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Counting gene reads...done ", '\n', sep="", collapse=""))
  save(geneCount, file='counts_genes.Rdata');
  rm(geneCount);
  gc();
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving gene reads count...done ", '\n', sep="", collapse=""))

}

