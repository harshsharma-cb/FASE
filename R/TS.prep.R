#' Title
#'
#' @import edgeR
#' @import matrixStats
#'
#' @param Gcount 
#' @param designM 
#' @param Groups 
#' @param threshold 
#' @param ... 
#'
#' @return

.removeLECountsTS <- function(Gcount, designM, Groups, threshold=6.32, ...){
  #removing LE count meta-features
  #require(edgeR)
  #require(matrixStats)
  counts <- Gcount
  ucounts <- counts
  counts <- DGEList(counts=counts, group=Groups)
  counts <- calcNormFactors(counts)
  try(fit <- voom(counts, designM), silent=TRUE)
  if(!exists("fit", inherits=FALSE)) return(NA) 
  ##removing zero counts false expression
  E <- ((ucounts !=0)*1)* fit$E
  
  ##Check low expressed reads #remove that have expression less than 6.32
  y <- as.matrix(E)
  if(is.null(Groups)){
    if(any(rowSums(designM) > 1)) 
    Groups <- rowSums(designM %*% t(designM))  else Groups <- which(designM != 0, arr.ind = TRUE)[ , 2]
    } else Groups <- Groups
  #-------------------------------
  Gfactor <- as.factor(Groups)
  n <- as.numeric(levels(Gfactor))
  GIndex <- match(Gfactor, n)
  #Filter the probes
  rs <- lapply(1:length(n), function(x) rowMedians(y[ , (Gfactor == n[x]), drop = FALSE]))
  rs <- do.call('cbind',rs)
  newy <- y[(rowMaxs(rs) > threshold),, drop=FALSE]
  return(newy)
}

#' Title
#'
#' @param geneID 
#' @param expression 
#' @param RMM 
#' @param iMM 
#'
#' @return

.cmm <- function(expression, RMM, iMM){
  #combining RMM and iMM to form combined membership matrix. Also to make RMM and iMM numbers for FJ and SJ constant.
  if(!is.null(RMM) && !is.null(dim(RMM))){
    rmm.index.rows <- match(rownames(expression), rownames(RMM))
    rmm.index.cols <- match(colnames(RMM), rownames(expression))
    RMM <- RMM[rmm.index.rows, !is.na(rmm.index.cols), drop=FALSE]  
  }
  #imm <- iMM[[geneID]]
  if(!is.null(iMM) && !is.null(dim(iMM))){
    imm.index.rows <- match(rownames(expression), rownames(iMM))
    imm.index.cols <- match(colnames(iMM), rownames(expression))
    iMM <- iMM[imm.index.rows, !is.na(imm.index.cols), drop=FALSE]
  }
  
  # making order of events consistent in rmm and imm and combining them into one matrix: CMM.
  if(is.null(iMM) || is.null(dim(iMM))){
    cmm <- RMM; return(cmm)
  }
  
  iMM <- ((iMM ==3)*0.5) + ((iMM == 2)*0) + ((iMM == 1)*3) #changing skipping junction to intron (3 to 0.5), intron to intron (2 to 0) and intron to exon (3 to 1). 

  if(is.null(RMM) || is.null(dim(RMM))){
    cmm <- iMM; return(cmm)
    }

  index <- match(rownames(iMM), rownames(RMM))
  cmm <- cbind(RMM[index,,drop=FALSE], iMM)
  
  return(cmm)
}

