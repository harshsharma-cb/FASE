#' RNA Sequencing raw data preprocessing
#'
#' @description Manual function to map reads with the reference genome, given SRA/fastq files. It also sorts and indexes the mapped reads for further processing. Reads produced by ppRawData can be summarized for genes, exons and introns using \code{\link{ppSumEIG}}. \code{\link{ppAuto}} is not required if ppRawData has been called. \cr
#'
#'             System requirements for ppRawData include: \cr
#'     \enumerate{
#'               \item fastq-dump (if files='SRA')
#'               \item tophat2
#'               \item samtools
#'               }
#'
#' @param folderSRA path of directory containing fastq or SRA files. (default=current directory)
#'
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
#' @param ... other parameters to be passed to tophat2.
#'
#' @import parallel
#'
#' @return \enumerate{
#' \item Mapped, sorted and indexed bam files. (Can be run separately using tophat2 and samtools or automatic wrapper function: \code{\link{ppAuto}})
#' \item Junction Matrix: Matrix with junction count reads. (Can be run separately using \code{\link{getJunctionCountMatrix}} or wrapper function: \code{\link{ppAuto}})
#' }
#'
#' @references \enumerate{
#'  \item https://CRAN.R-project.org/view=HighPerformanceComputing
#'  \item Sequence Read Archive Submissions Staff. Using the SRA Toolkit to convert .sra files into other formats. In: SRA Knowledge Base [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2011-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK158900/.
#'  \item https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump
#'  \item Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL. TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biol. 25;14(4):R36 (2013 Apr). http://ccb.jhu.edu/software/tophat.
#'  \item Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9.
#'  }
#' @export
#'

ppRawData <- function(folderSRA=FALSE, srlist = NULL, pairedend = FALSE, genomeBI, files='fastq', p=1, N= 6, r=44, mate_std_dev=30, read_edit_dist= 6, max_intron_length=10000, min_intron_length=50, segment_length=NULL,...) {


  # if(files=='fastq'){
  # paste("In order to make srlist (if not provided by user), fastq file should be named as:", '\n', 'Sample-S1_1.fastq (for paired-end reads).', '\n', "Sample-S1.fastq (for single-end reads).")
  # } else {
  #   paste("In order to make srlist (if not provided by user), sra file should be named as:", '\n', 'Sample-S1_1.sra (for paired-end reads).', '\n', "Sample-S1.sra (for single-end reads).")
  # }
  if(folderSRA==FALSE){
    folderSRA <- getwd();
  } else {folderSRA=folderSRA}

  start.time <- Sys.time()

  if(files == 'sra'){
    paste('\n', "Assuming fastq-dump is already installed and working.", '\n')
  }

  paste('\n', "Assuming tophat2 and samtools is already installed and working.", '\n')


  if(files == 'sra'){
    cat(paste('\n', "Task 1 out of 5:",'\n' ,"Converting SRA files to fastq files..." ,'\n'))
  }

  #folderSRA <- getwd(); # for preparation of packages
  if(is.null(srlist)){
    if (files == 'sra' && pairedend == TRUE) {
    #srlist <- dir(pattern = '.sra')
    srlist <- dir()
    if(length(srlist) == 0) stop('No SRA file found in SRA folder')

    ##Separate paired-end files using fastq-dump # I need to use the parellel for this
    tstamp <- Sys.time()
    cat(paste("[", tstamp,"]", " Converting paired-end sra files to fastq... ", '\n', sep="", collapse=""))
    new<- lapply(srlist, function(x) {
      system(paste('fastq-dump -I -v -v -v -v --split-files ', x, sep="", collapse=""), ignore.stdout = TRUE)
      return(TRUE)        })
    cat('\n')
    tstamp <- Sys.time()
    cat(paste("[",tstamp,"]", " Converting pair-end sra files to fastq...done ", '\n', sep="", collapse=""))
    srlist<- unlist(strsplit(srlist, '.sra', fixed=TRUE))
    save(srlist, file='srlist.Rdata')

  } else if (files == 'sra' && pairedend == FALSE) {
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
    srlist <- unlist(strsplit(srlist, '.sra', fixed=TRUE))
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
    save(srlist, file='srlist.Rdata')
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
    cat(paste('\n', "Task 2 out of 5:",'\n' ,"Mapping reads to genome using tophat2...", '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 1 out of 4:",'\n' ,"Mapping reads to genome using tophat2...", '\n'))
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
    cat(paste('\n', "Task 3 out of 5:",'\n' ,"Sorting mapped bam files using samtools...", '\n'))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 2 out of 4:",'\n' ,"Sorting mapped bam files using samtools...", '\n'))
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
    cat(paste('\n', "Task 4 out of 5:",'\n' ,"Indexing sorted bam files using samtools...", '\n', "Assuming samtools is already installed and working." , '\n' ))
  } else if (files == 'fastq'){
    cat(paste('\n', "Task 3 out of 4:",'\n' ,"Indexing sorted bam files using samtools...", '\n', "Assuming samtools is already installed and working." , '\n' ))
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
      cat(paste('\n', "Task 5 out of 5:",'\n' ,"Preparing junction matrix...", '\n' ))
    } else if (files == 'fastq'){
      cat(paste('\n', "Task 4 out of 4:",'\n' ,"Preparing junction matrix...", '\n' ))
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
  }
