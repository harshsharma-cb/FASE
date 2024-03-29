% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppFASE.R
\name{ppFASE}
\alias{ppFASE}
\title{Alternative Splicing Pre-processing}
\usage{
ppFASE(
  folderSRA = FALSE,
  gtf = gtf,
  exonCount = exonCount,
  intronCount = intronCount,
  JunctionMatrix = JunctionMatrix
)
}
\arguments{
\item{folderSRA}{directory containing fastq or SRA files.}

\item{gtf}{intron parsed gtf file of the organism.}

\item{JunctionMatrix}{junction count matrix. If ppRawData has been run, JunctionMatrix is saved in JunctionCounts.Rdata.}

\item{exonCounts}{list of summarized exon counts. If ppSumEIG has been run, exonCounts are saved in counts_exons.Rdata.}

\item{intronCounts}{list of summarized intron counts. If ppSumEIG has been run, intronCounts are saved in counts_introns.Rdata.}
}
\value{
\enumerate{
        \item Junction Matrix: Matrix with Junction count reads and their annotation. (Can be run separately using \code{\link{getJunctionCountMatrix}})
        \item RMM            : ReadMembershipMatrix. (Can be run separately using \code{\link{readMembershipMatrix}} or wrapper function: \code{\link{ppAuto}})
        \item iMM            : intronMembershipMatrix. (Can be run separately using \code{\link{intronMembershipMatrix}} or wrapper function: \code{\link{ppAuto}})
        \item Gcount         : A list of gene-wise read count summarization of meta-features times samples in the study. (Can be run separately using \code{\link{countMatrixGenes}} or wrapper function: \code{\link{ppAuto}})
}
}
\description{
Alternative Splicing preprocessing function. This function creates several prerequisite matrices for related to meta-features:junction matrix (by combining output of tophat2: junction.bed), ReadMembershipMatrix (RMM), IntronMembershipMatrix (IMM) and Gcount matrix in order to run ExonPointer and IntronPointer algorithms. \code{\link{ppFASE}} should be run only after tophat2 (or its wrapper function: \code{\link{ppRawData}}) has mapped all the raw read files and the reads have been summarized using \code{\link[Rsubread]{featureCounts}} (or its wrapper function: \code{\link{ppSumEIG}}).
}
\references{
Liao Y., Smyth G.K., Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47, e47 (2019).
}
