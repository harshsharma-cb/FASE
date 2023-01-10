#' Transcript Structure: seed propagation function
#'
#' @param cmm 
#' @param ts.p 
#' @param cjmm 
#' @param se 
#' @param s_exon 
#' @param event 
#' @param eventtype 
#' @param connected.exon 
#' @param s.exon 
#' @param keep.intron 
#'
#' @return

.propagation <- function(cmm, ts.p, cjmm, se, s_exon, event, eventtype, connected.exon, s.exon = NULL, keep.intron = FALSE){
  #finding how many TS can be formed with one seed exon pair. It is possible for an event to not have both right and left flanking (seed) exon, so putting these conditions:
#browser()
  #adding seed junctions
  if(eventtype == 'ep.incl'){
    if(!is.na(s_exon[1])) lsj <- intersect(rownames(cjmm)[cjmm[ , s_exon[1]] == 2], rownames(cjmm)[cjmm[ , event] == 2]) else lsj <- NULL
    if(!is.na(s_exon[2])) rsj <- intersect(rownames(cjmm)[cjmm[ , s_exon[2]] == 2], rownames(cjmm)[cjmm[ , event] == 2]) else rsj <- NULL
    ts.p[, c(lsj, rsj)] <- 1
    } else if(eventtype == 'ip.excl' || eventtype == 'ep.excl'){
          #junc.mat <- junc.mat
          s_exon <- s_exon[!is.na(s_exon)]
          if(length(s_exon) == 1)
          index <- rownames(cjmm)[cjmm[ , s_exon] == 2]
          else
          index <- intersect(rownames(cjmm)[cjmm[ , s_exon[1]] == 2], rownames(cjmm)[cjmm[ , s_exon[2]] == 2])
          ts.p[, index] <- 1
        }

  ##adding intron(s) to TS in keep.intron=T for EP incl event:
  if(keep.intron == TRUE){
    if((length(lsj) == 0) && !is.na(s_exon[1])) int.incl.l <- intersect(rownames(cmm)[(cmm[ , event] == 3)], rownames(cmm)[(cmm[ , s_exon[1]] == 3)]) else int.incl.l <- NULL
    if((length(rsj) == 0) && !is.na(s_exon[2])) int.incl.r <- intersect(rownames(cmm)[(cmm[ , event] == 3)], rownames(cmm)[(cmm[ , s_exon[2]] == 3)]) else int.incl.r <- NULL
    int.incl <- c(int.incl.l, int.incl.r)
    if(length(int.incl) > 0) ts.p[ , int.incl] <- 1
  }

  #making combinations of possible left and right structures of each seed exon pair
  lse.struct <- connected.exon[[1]][[s_exon[1]]]
  rse.struct <- connected.exon[[2]][[s_exon[2]]]
  n_ts <- if(length(lse.struct) > 0 && length(rse.struct) > 0) n_ts <- length(lse.struct)*length(rse.struct) else if(length(lse.struct) > 0) n_ts <- length(lse.struct) else if(length(rse.struct) > 0) n_ts <- length(rse.struct) else return(ts.p)

  inc.exon.combinations <- NULL
  if(!is.null(lse.struct) && !is.null(rse.struct)){
    #inc.exon.combinations <- expand.grid(lse.struct, rse.struct)
    inc.exon.combinations <- expand.grid(lapply(1:length(lse.struct), function(x) lse.struct[[x]]), lapply(1:length(rse.struct), function(x) rse.struct[[x]]))
    inc.exon.combinations <- lapply(1:nrow(inc.exon.combinations), function(x) as.vector(unlist(inc.exon.combinations[x, ])))
    } else if(!is.null(lse.struct)){
      inc.exon.combinations <- lapply(1:length(lse.struct), function(x) lse.struct[[x]])
      } else if(!is.null(rse.struct)){
        inc.exon.combinations <- lapply(1:length(rse.struct), function(x) rse.struct[[x]])
      }


  #replicating ts.p row nrow(inc.exon.combinations)-1 times
  if(n_ts > 1)
  ts.p <- ts.p[rep(1, times = length(inc.exon.combinations)), ]
  
  #Propagation
  for(i in 1:nrow(ts.p)){
    ts.p[i, match(inc.exon.combinations[[i]], colnames(ts.p))] <- 1
  }
  

  ##adding junctions
  incl.exons.ts.p <- lapply(1:nrow(ts.p), function(x) colnames(ts.p)[ts.p[x, ] == 1]) #extracting included exons in each transcript structure
  incl.exons.ts.p <- lapply(1:length(incl.exons.ts.p ), function(x) incl.exons.ts.p[[x]][grep("EX", incl.exons.ts.p [[x]])]) # sometimes junctions are also included, removing them here
  
  if(eventtype == 'ep.incl' || eventtype == 'ep.excl'){
    l.exonlist <- colnames(ts.p)[grep("EX", colnames(ts.p))[1:((match(event, colnames(ts.p))) - 1)]]
    r.exonlist <- colnames(ts.p)[grep("EX", colnames(ts.p))[((match(event, colnames(ts.p))) + 1):ncol(ts.p)]]
    } else{
      if(length(s.exon$ls_exon.connected.exons) > 0) l.exonlist <- colnames(ts.p)[grep("EX", colnames(ts.p))[1:(match(s.exon$ls_exon.connected.exons, colnames(ts.p)))]] else l.exonlist <- NA
      if(length(s.exon$rs_exon.connected.exons) > 0) r.exonlist <- colnames(ts.p)[grep("EX", colnames(ts.p))[(match(s.exon$rs_exon.connected.exons, colnames(ts.p))):ncol(ts.p)]] else r.exonlist <- NA
    }

  #dividing included exons into left and right
  l.exons <- lapply(1:length(incl.exons.ts.p), function(x)
  {
    index <- match(incl.exons.ts.p[[x]], l.exonlist)
    l.exons <- incl.exons.ts.p[[x]][!is.na(index)]
    })
  
  r.exons <- lapply(1:length(incl.exons.ts.p), function(x)
  {
    index <- match(incl.exons.ts.p[[x]], r.exonlist)
    r.exons <- incl.exons.ts.p[[x]][!is.na(index)]
    })
  
  #adding junctions to left side  
  if(length(l.exons[[i]]) > 1){
    #ls_index <- which(colnames(cjmm) == s_exon[1])
    for(i in 1:length(l.exons)){
      for(j in 1:(length(l.exons[[i]])) - 1){
        fj1 <- rownames(cjmm)[which(cjmm[ , l.exons[[i]][j]] == 2)]
        fj2 <- rownames(cjmm)[which(cjmm[ , l.exons[[i]][j+1]] ==2)]
        fj <- intersect(fj1, fj2)
        ts.p[i, fj] <- 1
      }
    }
  }
  

  #adding junctions to right side
  if(length(r.exons[[i]]) > 1){
    for(i in 1:length(r.exons)){
      for(j in 1:(length(r.exons[[i]])) - 1){
        fj1 <- rownames(cjmm)[which(cjmm[,r.exons[[i]][j]] == 2)]
        fj2 <- rownames(cjmm)[which(cjmm[,r.exons[[i]][j+1]] == 2)]
        fj <- intersect(fj1, fj2)
        ts.p[i, fj] <- 1
      }
    }
  }


  return(ts.p)
}
