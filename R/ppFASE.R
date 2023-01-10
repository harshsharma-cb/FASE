#' Alternative Splicing Pre-processing
#'
#' @description Alternative Splicing preprocessing function. This function creates several prerequisite matrices for related to meta-features:junction matrix (by combining output of tophat2: junction.bed), ReadMembershipMatrix (RMM), IntronMembershipMatrix (IMM) and Gcount matrix in order to run ExonPointer and IntronPointer algorithms. \code{\link{ppFASE}} should be run only after tophat2 (or its wrapper function: \code{\link{ppRawData}}) has mapped all the raw read files and the reads have been summarized using \code{\link[Rsubread]{featureCounts}} (or its wrapper function: \code{\link{ppSumEIG}}).
#'
#' @param folderSRA directory containing fastq or SRA files.
#' @param gtf intron parsed gtf file of the organism.
#' @param exonCounts list of summarized exon counts. If ppSumEIG has been run, exonCounts are saved in counts_exons.Rdata.
#' @param intronCounts list of summarized intron counts. If ppSumEIG has been run, intronCounts are saved in counts_introns.Rdata.
#' @param JunctionMatrix junction count matrix. If ppRawData has been run, JunctionMatrix is saved in JunctionCounts.Rdata.
#'
#' @import Rsubread
#'
#' @return \enumerate{
#'         \item Junction Matrix: Matrix with Junction count reads and their annotation. (Can be run separately using \code{\link{getJunctionCountMatrix}})
#'         \item RMM            : ReadMembershipMatrix. (Can be run separately using \code{\link{readMembershipMatrix}} or wrapper function: \code{\link{ppAuto}})
#'         \item iMM            : intronMembershipMatrix. (Can be run separately using \code{\link{intronMembershipMatrix}} or wrapper function: \code{\link{ppAuto}})
#'         \item Gcount         : A list of gene-wise read count summarization of meta-features times samples in the study. (Can be run separately using \code{\link{countMatrixGenes}} or wrapper function: \code{\link{ppAuto}})
#' }
#'
#' @references Liao Y., Smyth G.K., Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47, e47 (2019).
#'
#' @export
#'

ppFASE <- function(folderSRA=FALSE, gtf, exonCount, intronCount, JunctionMatrix) {

  if(folderSRA==FALSE){
    folderSRA <- getwd();
  } else {folderSRA=folderSRA}

  start.time <- Sys.time()


  ##Generating readMembershipMatrix
  # source('/home/harsh/Scripts_All/readMembershipMatrix.R')

  cat(paste('\n', "Task 1 out of 3:",'\n' ,"Creating Read Membership Matrix...", '\n' ))

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

  cat(paste('\n', "Task 2 out of 3:",'\n' ,"Creating Intron Membership Matrix...", '\n' ))
  load('Annotation.Rdata')
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Getting intronMembershipMatrix... ", '\n', sep="", collapse=""))
  iMM <- intronMembershipMatrix(annotation=annotation)
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

  cat(paste('\n', "Task 3 out of 3:",'\n' ,"Creating gene count matrix...", '\n' ))

  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Generating per Gene meta-features count... ", '\n', sep="", collapse=""))
  #load('counts_exons.Rdata')
  #load('counts_introns.Rdata')
  #load("Annotation.Rdata")
  #load("JunctionCounts.Rdata")
  Gcount <- countMatrixGenes(JunctionMatrix, annotation=annotation, intronList=intronCount, exonList= exonCount)
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Generating per Gene count...done ", '\n', sep="", collapse=""))
  save(Gcount, file='Gcount.Rdata')
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Saving per Gene count...done ", '\n', sep="", collapse=""))
  ttime<- tstamp -start.time
  cat('Total time of spent:', ttime,sep='')
  rm(list=ls())
  gc()
  tstamp <- Sys.time()
  cat(paste("[", tstamp, "]", " Preprocessing has been done successfully! ", '\n', sep="", collapse=""))

  ##All finished
}


