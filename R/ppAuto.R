#' RNA Seq and Alternative Splicing preprocessing function
#'
#' @description ppAuto is a wrapper function for several tools and functions that perform preprocessing of RNA Sequencing data. This function performs preprocessing that includes mapping of reads, sorting and indexing of bam files, to summarization of read counts for exons, introns, genes and junctions. ppAuto also creates several prerequisite matrices including junction matrix, ReadMembershipMatrix (RMM), IntronMembershipMatrix (iMM) and Gcount matrix in order to run ExonPointer and IntronPointer algorithms. \cr
#'
#' System requirements for ppAuto include: \cr
#'     \enumerate{
#'               \item fastq-dump (if files='SRA')
#'               \item tophat2
#'               \item samtools
#' }
#'
#' @param folderSRA path of directory containing fastq or SRA files. (default=current directory)
#' @param srlist list of unique sample names of fastq/SRA files created by default in the function. Please follow naming convention for the sample files: \cr
#'                   For SRA files  : "Sample-S1_1" "Sample-S1_2" (for paired-end reads) and "Sample-S1" (for single-end reads). \cr
#'                   For fastq files: "Sample-S1_1.fastq" "Sample-S1_2.fastq" (for paired-end reads) and "Sample-S1.fastq" (for single-end reads). \cr
#'
#' @param pairedend boolean, TRUE if reads are paired-end and FALSE if reads are single-end. All files should be either single-end or paired-end. (default=FALSE)
#' @param genomeBI path of genome build of the organism created using bowtie2-build command.
#' @param gtf intron parsed gtf file of the organism. Please check \code{\link{intronGTFparser}} to generate intron parsed gtf file (to generate intron read counts).
#' @param files type of raw read file: fastq or sra (downloaded from NCBI). All files should be in same format and have same read length. (default=fastq)
#' @param p number of threads to be utilized by samtools and Rsubread package. (default=1)
#' @param N accepted read mismatches. Reads with more than N mismatches are discarded. (default=6) [tophat2 parameter]
#' @param r expected inner distance between mate pair. (default=44) [tophat2 parameter]
#' @param mate_std_dev the standard deviation for the distribution on inner distances between mate pairs. (default=30) [tophat2 parameter]
#' @param read_edit_dist final read alignments having more than these many edit distance are discarded. (default=6) [tophat2 parameter]
#' @param max_intron_length when searching for junctions ab initio, TopHat2 will ignore donor/acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read. (default=10000) [tophat2 parameter]
#' @param min_intron_length topHat2 will ignore donor/acceptor pairs closer than this many bases apart. (default=50) [tophat2 parameter]
#' @param segment_length each read is divided into this length and mapped independently to find junctions. [tophat2 parameter]
#' @param ... other parameter to be passed to tophat2.
#'
#' @import Rsubread
#' @import parallel
#'
#' @return \enumerate{
#'         \item Mapped, sorted and indexed bam files. (Can be run separately using tophat2 and samtools or wrapper function: \code{\link{ppRawData}})
#'         \item Lists of gene counts, exon counts and intron counts saved in folderSRA directory as respective Rdata files. (Can be run separately using \code{\link[Rsubread]{featureCounts}} or wrapper function: \code{\link{ppSumEIG}})
#'         \item Junction Matrix: Matrix with annotated junction count reads. (Can be run separately using \code{\link{getJunctionCountMatrix}} or wrapper function: \code{\link{ppRawData}})
#'         \item RMM            : ReadMembershipMatrix. (Can be run separately using \code{\link{readMembershipMatrix}} or wrapper function: \code{\link{ppFASE}})
#'         \item iMM            : intronMembershipMatrix. (Can be run separately using \code{\link{intronMembershipMatrix}} or wrapper function: \code{\link{ppFASE}})
#'         \item Gcount         : A list of gene-wise read count summarization of meta-features times samples in the study. (Can be run separately using \code{\link{countMatrixGenes}} or wrapper function: \code{\link{ppFASE}})
#'         }
#'
#' @references \enumerate{
#' \item Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47, e47 (2019).
#' }
#'
#' @export
#'

ppAuto <- function(folderSRA = FALSE, srlist = NULL, pairedend = FALSE, genomeBI, gtf, files='fastq', p=1, N= 6, r=44, mate_std_dev=30, read_edit_dist= 6, max_intron_length=10000, min_intron_length=50, segment_length=NULL,...) {
  #devtools::load_all()

  if(folderSRA==FALSE){
    folderSRA <- getwd();
  } else {folderSRA=folderSRA}

  start.time <- Sys.time()
  
  paste('\n', "Assuming genome index has already been built using bowtie2-build")

  if(files == 'sra'){
    paste('\n', "Assuming fastq-dump is already installed and working.", '\n')
  }

  paste('\n', "Assuming tophat2 and samtools is already installed and working.", '\n')

  GTF_corrected <- readline(prompt = "Has intronGTFparser() been run (to find introns in gtf)? y/n ")
  if(GTF_corrected == 'n'){
    stop('Please parse gtf using intronGTFparser() to continue and use _corrected.gtf for analysis.')
  }

  if(files == 'sra'){
    cat(paste('\n', "Task 1 out of 11:",'\n' ,"Converting SRA files to fastq files..." ,'\n'))
  }

  if(is.null(srlist)){
    if (files == 'sra' && pairedend == TRUE) {
    #srlist <- dir(pattern = '.sra')
    srlist <- dir()
    if(length(srlist) == 0) stop('No SRA file found in SRA folder')

    ##Separate paired-end files using fastq-dump # I need to use the parallel for this
    tstamp <- Sys.time()
    cat(paste("[", tstamp,"]", " Converting paired-end sra files to fastq... ", '\n', sep="", collapse=""))
    new<- lapply(srlist, function(x) {
      system(paste('fastq-dump -I -v -v -v -v --split-files ', x, sep="", collapse=""), ignore.stdout = TRUE)
      return(TRUE)        })
    cat('\n')
    tstamp <- Sys.time()
    cat(paste("[",tstamp,"]", " Converting paired-end sra files to fastq...done ", '\n', sep="", collapse=""))
    #srlist<- unlist(strsplit(srlist, '.sra', fixed=TRUE))
    save(srlist, file='srlist.Rdata')

  } else if (files == 'sra' && pairedend == FALSE) {
    #srlist <- dir(pattern = '.sra')
    srlist <- dir()
    if(length(srlist) == 0) stop('No SRA file found in SRA folder')

    ##Separate paired-end files using fastq-dump # I need to use the parellel for this
    tstamp <- Sys.time()
    cat(paste("[", tstamp,"]", " Converting sra files to fastq... ", '\n', sep="", collapse=""))
    new<- lapply(srlist, function(x) {
      system(paste('fastq-dump -I -v ', x, sep="", collapse=""), ignore.stdout = TRUE)
      return(TRUE)        })
    cat('\n')
    tstamp <- Sys.time()
    cat(paste("[", tstamp,"]", " Converting sra files to fastq...done ", '\n', sep="", collapse=""))
    # srlist<- unlist(strsplit(srlist,'.sra',fixed=TRUE))
    #srlist <- unlist(strsplit(srlist, '.sra', fixed=TRUE))
    save(srlist, file='srlist.Rdata')
  } else if (files == 'fastq' && pairedend == TRUE) {
    srlist <- dir(pattern = '.fastq')
    if(length(srlist) == 0) stop('No fastq file found in the folder')
    srlist <- unlist(strsplit(srlist, '.fastq', fixed=TRUE))

    #i have to make something that should take into pairedend into account otherwise there are lots file now
    ##might not work for all purpose
    nsrlist <- strsplit(srlist, '_', TRUE)
    nsrlist <- do.call('rbind', nsrlist)
    srlist <- unique(nsrlist[,1])
    save(srlist, file='srlist.Rdata')  #only for test purpose, I should delete this later.
  } else if (files == 'fastq'&& pairedend == FALSE) {
    srlist <- dir(pattern = '.fastq')
    if(length(srlist) == 0) stop('No fastq file found in the folder')
    srlist <- unlist(strsplit(srlist, '.fastq', fixed=TRUE))
    save(srlist, file='srlist.Rdata')
  }
  }
  

  ##Running tophat

  #segment_length
  if(is.null(segment_length)){
    fsl <- dir(pattern = '*.fastq')
    fsl <- fsl[1]
    system(paste('a=`head -2 ', fsl, ' | tail -1 | wc -L` ; echo $a > .rl.csv' ))
    segment_length <- as.integer(as.vector(read.csv(file='.rl.csv', header = F))/2-5)
  } else {
    segment_length <- segment_length
  }


  if(files == 'sra'){
    cat(paste('\n', "Task 2 out of 11:",'\n' ,"Mapping reads to genome using tophat2...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 1 out of 10:",'\n' ,"Mapping reads to genome using tophat2...", '\n'))
  }

  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Running TopHat... ", '\n', sep="", collapse=""))

  if(pairedend == TRUE){
    new<- lapply(srlist, function(x) {
      system(paste("tophat2 -o ",paste(folderSRA,'/',x,'_tophat_out',sep="",collapse=""), paste(' -p ',p,' -N ',N,' -r ',r,' --segment-length ', segment_length, ' --max-intron-length ',max_intron_length,' --no-coverage-search --min-intron-length ',min_intron_length,' --mate-std-dev ', mate_std_dev,' --read-edit-dist ',read_edit_dist,sep="",collapse=""),
                   paste(' ',genomeBI,' ',sep='',collapse=""),paste(folderSRA,'/',x,'_1.fastq ',folderSRA,'/',x,'_2.fastq ',sep='',collapse=""),sep='',collapse=""),ignore.stdout = TRUE)
      return(TRUE)		  })
    cat('\n')
    tstamp <- Sys.time()
    cat(paste("[",tstamp,"]"," Running TopHat...done ",'\n',sep="",collapse=""))

  } else {
    new<- lapply(srlist, function(x) {
      system(paste("tophat2 -o ",paste(folderSRA,'/',x,'_tophat_out',sep="",collapse=""), paste(' -p ',p,' -N ',N,' -r ',r,' --segment-length ', segment_length,' --max-intron-length ',max_intron_length,' --no-coverage-search --min-intron-length ',min_intron_length,' --mate-std-dev ', mate_std_dev,' --read-edit-dist ',read_edit_dist,sep="",collapse=""),
                   paste(' ',genomeBI,' ',sep='',collapse=""),paste(folderSRA,'/',x,'.fastq ',sep='',collapse=""),sep='',collapse=""),ignore.stdout = TRUE)
      return(TRUE)		  })
  }
  cat('\n')
  tstamp <- Sys.time()
  cat(paste("[", tstamp,"]", " Running TopHat...done ", '\n', sep="", collapse=""))



  ##Sorting files generated by topHat using samtools # I need to use the parallel for this

  if(files == 'sra'){
    cat(paste('\n', "Task 3 out of 11:",'\n' ,"Sorting mapped bam files using samtools...", '\n'))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 2 out of 10:",'\n' ,"Sorting mapped bam files using samtools...", '\n'))
  }

    tstamp <- Sys.time()
    cat(paste("[",tstamp,"]"," Sorting files with samtools... ",'\n',sep="",collapse=""))
    new<- mclapply(srlist, function(x) {
      cat('Sorting ',x,' ',sep='')
      system(paste(' samtools sort ',folderSRA,'/',x,'_tophat_out/accepted_hits.bam',' ', '-o', ' ' ,folderSRA,'/',x,'_tophat_out/accepted_hits_sorted.bam',sep="",collapse=""),ignore.stdout = TRUE)
      return(TRUE)				   }, mc.cores = p)
    tstamp <- Sys.time()
    cat(paste("[",tstamp,"]"," Sorting files with samtools...done ",'\n',sep="",collapse=""))

  ##Make index of these files
  #samtools index accepted_hits_sorted.bam

  if(files == 'sra'){
    cat(paste('\n', "Task 4 out of 11:",'\n' ,"Indexing sorted bam files using samtools...", '\n', "Assuming samtools is already installed and working." , '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 3 out of 10:",'\n' ,"Indexing sorted bam files using samtools...", '\n', "Assuming samtools is already installed and working." , '\n' ))
  }

    tstamp <- Sys.time()
    cat(paste("[",tstamp,"]"," Indexing files with samtools... ",'\n',sep="",collapse=""))
    new <- mclapply(srlist, function(x) {
      cat('Indexing ', x,' ', sep='')
      system(paste('samtools index ',folderSRA,'/',x,'_tophat_out/accepted_hits_sorted.bam' ,sep="",collapse=""),ignore.stdout = TRUE)
      return(TRUE)				   }, mc.cores = p)
    tstamp <- Sys.time()
    cat(paste("[", tstamp, "]", " Indexing files with samtools...done ", '\n', sep="", collapse=""))



    ##Counting Junction reads
    #source('/home/harsh/Scripts_All/getJunctionCountMatrix.R')
    if(files == 'sra'){
      cat(paste('\n', "Task 5 out of 11:",'\n' ,"Preparing junction matrix...", '\n' ))
    } else if (files == 'fastq'){
      cat(paste('\n', "Task 4 out of 10:",'\n' ,"Preparing junction matrix...", '\n' ))
    }
    tstamp <- Sys.time()
    cat(paste("[", tstamp, "]", " Counting Junction reads... ", '\n', sep="", collapse=""))
    jfiles <- unlist(lapply(srlist, function(x) paste(folderSRA, '/' ,x ,'_tophat_out/junctions.bed' ,sep='', collapse="")))
    JunctionMatrix <- getJunctionCountMatrix(jfiles)
    tstamp <- Sys.time()
    cat(paste("[", tstamp, "]", " Counting Juntion reads...done ", '\n', sep="", collapse=""))
    ##Change name of sample here and will be same through out##
    colnames(JunctionMatrix) <- c(colnames(JunctionMatrix)[1:5],srlist)
    rownames(JunctionMatrix)<- JunctionMatrix[,5]
    save(JunctionMatrix, file='JunctionCounts.Rdata')
    tstamp <- Sys.time()
    cat(paste("[", tstamp, "]", " Saving Junction reads count...done ", '\n', sep="", collapse=""))


  ### Summarization
  srlistbam <- paste(folderSRA,'/',srlist,'_tophat_out/accepted_hits_sorted.bam',sep='')
  #srlistbam <- paste(folderSRA,'/',srlist,'_tophat_out/accepted_hits_sorted.bam',sep='')


  ##Counting exon reads
  if(files == 'sra'){
    cat(paste('\n', "Task 6 out of 11:",'\n' ,"Summarization of reads for exons...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 5 out of 10:",'\n' ,"Summarization of reads for exons...", '\n' ))
  }

  #require(Rsubread)
  tstamp <- Sys.time()
  cat(paste("[", tstamp,"]", " Counting exon reads... ", '\n', sep="", collapse=""))
  exonCount <- featureCounts(files= srlistbam, isPairedEnd= pairedend,requireBothEndsMapped=pairedend,GTF.featureType = "exon", GTF.attrType = "gene_id", useMetaFeatures = FALSE, isGTFAnnotationFile=TRUE, annot.ext=gtf,allowMultiOverlap=TRUE,nthreads=p)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Counting exon reads...done ", '\n', sep="", collapse=""))
  save(exonCount, file='counts_exons.Rdata'); rm(exonCount);gc();
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving exon reads count...done ", '\n', sep="", collapse=""))


  ##Counting intron reads
  if(files == 'sra'){
    cat(paste('\n', "Task 7 out of 11:",'\n' ,"Summarization of reads for introns...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 6 out of 10:",'\n' ,"Summarization of reads for introns...", '\n' ))
  }
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
  if(files == 'sra'){
    cat(paste('\n', "Task 8 out of 11:",'\n' ,"Summarization of reads for genes...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 7 out of 10:",'\n' ,"Summarization of reads for genes...", '\n' ))
  }

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


  ##Generating readMembershipMatrix
  # source('/home/harsh/Scripts_All/readMembershipMatrix.R')
  if(files == 'sra'){
    cat(paste('\n', "Task 9 out of 11:",'\n' ,"Creating Read Membership Matrix...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 8 out of 10:",'\n' ,"Creating Read Membership Matrix...", '\n' ))
  }

  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Getting readMembershipMatrix... ", '\n', sep="", collapse=""))
  RMM <- readMembershipMatrix(gtf, JunctionMatrix)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Getting readMembershipMatrix...done ", '\n', sep="", collapse=""))
  save(RMM, file='RMM.Rdata');
  rm(RMM);
  gc();
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving readMembershipMatrix...done ", '\n', sep="", collapse=""))


  #intronMembershipMatrix
  # source('/home/harsh/Scripts_All/intronMembershipMatrix.R')
  if(files == 'sra'){
    cat(paste('\n', "Task 10 out of 11:",'\n' ,"Creating Intron Membership Matrix...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 9 out of 10:",'\n' ,"Creating Intron Membership Matrix...", '\n' ))
  }

  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Getting intronMembershipMatrix... ", '\n', sep="", collapse=""))
  load('Annotation.Rdata')
  iMM <- intronMembershipMatrix(annotation = annotation)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Getting intronMembershipMatrix...done ", '\n', sep="", collapse=""))
  save(iMM, file='iMM.Rdata');
  rm(iMM);
  gc();
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving intronMembershipMatrix...done ", '\n', sep="", collapse=""))

  ##Generating genecountMatrix
  # source('/home/harsh/Scripts_All/countMatrixGenes.R')
  # source('/home/harsh/Scripts_All/countMatrixGenes.R')
  if(files == 'sra'){
    cat(paste('\n', "Task 11 out of 11:",'\n' ,"Creating gene count matrix...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 10 out of 10:",'\n' ,"Creating gene count matrix...", '\n' ))
  }
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Generating per Gene meta-features count... ", '\n', sep="", collapse=""))
  load('counts_exons.Rdata')
  load('counts_introns.Rdata')
  load("Annotation.Rdata")
  Gcount <- countMatrixGenes(JunctionMatrix, annotation=annotation, intronList=intronCount, exonList= exonCount)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Generating per Gene count...done ", '\n', sep="", collapse=""))
  save(Gcount, file='Gcount.Rdata')
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving per Gene count...done ", '\n', sep="", collapse=""))
  ttime<- tstamp -start.time
  cat('Total time of spent:', ttime,sep='')
  rm(list=ls())
  ;gc()
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Preprocessing has been done successfully! ", '\n', sep="", collapse=""))

  ##All finished
}

############################################################
##Junction matrix names of the chromosome are in different nomenclature than gtf so changing them
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
  save(JunctionMatrix,file="JunctionMatrix.Rdata")
}


