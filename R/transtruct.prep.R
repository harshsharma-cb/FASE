#' Transtruct preparation function
#'
#' @param cmm 
#' @param event 
#' @param eventlist 
#' @param se 
#' @param ls_exon 
#' @param rs_exon 
#' @param eventtype 
#' @param junc.mat 
#' @param p 
#' @param transtruct.ep.incl 
#' @param transtruct.ep.excl 
#' @param transtruct.ip.incl 
#' @param transtruct.ip.excl 
#' @param s.exon 
#' @param keep.intron 
#'
#' @return
#' 
#'
.transtruct.prep <- function(cmm, event, eventlist, se = s.exon$se, ls_exon = s.exon$ls_exon, rs_exon = s.exon$rs_exon, eventtype, junc.mat = NULL , p, transtruct.ep.incl = NULL, transtruct.ep.excl = NULL, transtruct.ip.incl = NULL, transtruct.ip.excl = NULL, s.exon = NULL, keep.intron){
  #browser()
  index <- grep('JUNC', rownames(cmm))
  if(length(index) == 1){
    cjmm <- t(cmm[index, ])
    rownames(cjmm) <- rownames(cmm)[index]
    } else
      cjmm <- cmm[index, ] #to reduce matrix size to junctions only

  if(!is.null(se)){ #if no seed exons, e.g. flanking introns not included, ts.p will not go through this
    se = s.exon$se; ls_exon = s.exon$ls_exon; rs_exon = s.exon$rs_exon
    
    #finding transcript structures starting from all seed exons
    connected.exon <- connected.exons.igraph(cmm = cmm, event = event, ls_exon = ls_exon, rs_exon = rs_exon, s.exon = s.exon, eventtype = eventtype)
    
    # #calculating total possible structures. If total_n_ts>10,000, ts would be discarded
    # total_n_ts <-  sum(unlist(lapply(1:length(connected.exon[[1]]), function(x) if(!is.null(connected.exon[[1]][[x]])) length(connected.exon[[1]][[x]]) else 1)))*sum(unlist(lapply(1:length(connected.exon[[2]]), function(x) if(!is.null(connected.exon[[2]][[x]])) length(connected.exon[[2]][[x]]) else 1)))
    # if(total_n_ts > 10000){
    #   cat(event, " has too many structures, possibly due to alternative 3'/5' splice site Alternative Splicing or gene fusion.")
    #   ts.result <- NULL
    #   return(ts.result)
    # }

    #Generating transcript structure
    for(z in 1:ncol(se)){
      ts.p <- matrix(0, ncol = nrow(cmm), nrow = 1)
      colnames(ts.p) <- c(colnames(cjmm), rownames(cjmm))
      
      #adding event to ts.p
      if(eventtype == 'ep.incl' || eventtype == 'ip.incl') ts.p[ , event] <- 1
      
      # #adding seed exons to ts.p
      s_exon <- as.vector(se[,z])
      index <- match(s_exon, colnames(ts.p))
      ts.p[, index] <- 1
      
      #Propagation
      ts.p <- .propagation(cmm = cmm, ts.p = ts.p, cjmm = cjmm, se = se, s_exon = s_exon, event = event, eventtype = eventtype, connected.exon = connected.exon, s.exon = s.exon, keep.intron = keep.intron)
      
      if(z == 1){
        ts <- ts.p
        } else if(z > 1){
          ts <- rbind(ts, ts.p)
        }

      #adding rownames to ts.p
      if(z == ncol(se)){
        ts.rownames <- paste0(event, "_", c(1:nrow(ts)))
        rownames(ts) <- ts.rownames
      }  
    }    
  } else{
      ts <- matrix(0, ncol = nrow(cmm), nrow = 1)
      colnames(ts) <- c(colnames(cjmm), rownames(cjmm))
      ts.rownames <- paste0(event, "_", c(1:nrow(ts)))
      rownames(ts) <- ts.rownames
    }


  #adding event to TS.
  if(eventtype == 'ep.incl' || eventtype == 'ip.incl') ts[ , event] <- 1
  
  #removing skipping junction to event (it maybe included in case seed exons are 3'/5' alt spliced)
  if(eventtype == 'ep.incl' || eventtype == 'ip.incl'){
    skipj <- rownames(cmm)[cmm[,event] == 0.5]
    ts[ , skipj] <- 0
  } else if(eventtype == 'ep.excl'){
    flankj <- rownames(cmm)[cmm[,event] == 2]
    ts[ , flankj] <- 0
  }
  
  
  #Combining TS for different events
  if(p > 1 && eventtype == 'ep.incl')
  ts.new <- transtruct.ep.incl$ts.new
  if(p > 1 && eventtype == 'ep.excl')
  ts.new <- transtruct.ep.excl$ts.new
  if(p > 1 && eventtype == 'ip.incl')
  ts.new <- transtruct.ip.incl$ts.new
  if(p > 1 && eventtype == 'ip.excl')
  ts.new <- transtruct.ip.excl$ts.new  
  
  if(length(eventlist) == 1){ 
    ts.new <- ts
    } else if(p == 1) {
      ts.new <- ts
      } else {
        ts.new <- rbind(ts.new, ts)
  } # ts.new is new ts
  
  #Each transcript should have at least 2 meta-features.
  if(p == length(eventlist)){
    #If there is only one TS for EP or IP, to convert it from vector to matrix and provide event name as its rownames
    if(nrow(ts.new) == 0){
      ts.new <- t(as.matrix(ts.new))
      rownames(ts.new)[1] <- event
    }
  }  
  
  ts.result <- list('ts.new' = ts.new)
  return(ts.result)
  
}
