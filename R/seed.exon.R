#' Seed exons
#'
#' @param cmm 
#' @param ep_event 
#' @param keep.intron 
#'
#' @return

.seed.exon.ep.incl <- function(cmm, ep_event, keep.intron){
  #making list of all exons
  exon.list <- colnames(cmm)[grep("EX", colnames(cmm))] 
  #checking flanking junctions to ep incl event
  s.junc <- rownames(cmm)[cmm[,ep_event] == 2]
  #if no FJ and keep.intron = FALSE, no need to continue
  if(length(s.junc) == 0 && keep.intron == FALSE){
    s.exon <- list('se' = NULL, 'n_ts' = NULL, 'ls_exon' = NULL, 'rs_exon' = NULL)
    warning('No seed exons found for', ep_event)
    } else{
      ls.exons <- exon.list[1:(match(ep_event, exon.list) - 1)]
      rs.exons <- exon.list[(match(ep_event, exon.list) + 1):length(exon.list)]
      if(length(s.junc) > 0){
        #exons associated with flanking junctions
        s.exon.junc <- unique(unlist(lapply(1:length(s.junc), function(x) colnames(cmm)[cmm[s.junc[x], exon.list] == 2])))
        #dividing seed exons into left and right seed exon
        ls_exon.junc <- ls.exons[!is.na(match(ls.exons, s.exon.junc))] #dividing seed exons into left seed
        rs_exon.junc <- rs.exons[!is.na(match(rs.exons, s.exon.junc))] #dividing seed exons into right seed
        } else{
          ls_exon.junc <- NULL
          rs_exon.junc <- NULL
        }
        #exons associated with flanking introns
        if(keep.intron == TRUE){
          s.int <- rownames(cmm)[cmm[,ep_event] == 3] #flanking intron(s)
          if(length(s.int) > 0) s.exon.int <- unique(unlist(lapply(1:length(s.int), function(x) colnames(cmm)[cmm[s.int[x], exon.list] == 3]))) else s.exon.int <- NULL #exon(s) associated with flanking intron(s)

          if(length(ls_exon.junc) == 0) ls_exon.int <- ls.exons[!is.na(match(ls.exons, s.exon.int))] else ls_exon.int <- NULL
          if(length(rs_exon.junc) == 0) rs_exon.int <- rs.exons[!is.na(match(rs.exons, s.exon.int))] else rs_exon.int <- NULL
        } else{
            ls_exon.int <- NULL
            rs_exon.int <- NULL
          } 

        #combining seed exons from flanking junction and flanking intron
        ls_exon <- unique(c(ls_exon.junc, ls_exon.int)) 
        rs_exon <- unique(c(rs_exon.junc, rs_exon.int)) 


        s.exon <- c(ls_exon, rs_exon)
        #s.exon <- s.exon.junc
        s.exon <- s.exon[lapply(s.exon,length) > 0] # removing empty elements
        # s.exon <- unique(unlist(s.exon))

        ## preparing matrix/list of left and right seed exon for each transcript
        if(length(ls_exon) == 0  && length(rs_exon) == 0){
          warning(ep_event, " does not have any flanking exons")
        } else if(length(ls_exon) == 0){
            n_ts <- length(rs_exon)
          } else if(length(rs_exon) == 0){
              n_ts <- length(ls_exon)
            } else{
                n_ts <- length(ls_exon)*length(rs_exon)
              }
      #initializing seed exon matrix. It will contain all possible combinations of left and right seed exons.
      se <- matrix(nrow=2, ncol=n_ts, dimnames=list(c('lse', 'rse')))
        
      # making combinations of left and right seed exons.
      if(length(ls_exon > 0) & length(rs_exon > 0)){
        se.p <- expand.grid('lse' = ls_exon, 'rse' = rs_exon)
        se[1, ] <- as.vector(se.p[, 1])
        se[2, ] <- as.vector(se.p[, 2])
      } else if(length(ls_exon) == 0){
          se[2, ] <- as.vector(rs_exon)
        } else if(length(rs_exon) == 0){
            se[1, ] <- as.vector(ls_exon)
          }
          s.exon <- list('se'=se, 'n_ts'=n_ts, 'ls_exon'=ls_exon, 'rs_exon'=rs_exon, 'keep.intron' = keep.intron)
      }
        return(s.exon)
}

#' Title
#'
#' @param ip_event 
#' @param cmm 
#' @param annotation 
#'
#' @return

.s.exon.ip.incl <- function(ip_event = ip_event, cmm = cmm, annotation = annotation){
  #annotation <- annotation[which(annotation$genes == geneID),]
  annotation <- annotation[order(annotation$V4), ]

  #finding seed exons to IP event
  s_exon <- colnames(cmm)[cmm[ ,ip_event] == 3]
  ls_exon <- s_exon[!is.na(match(s_exon, annotation[1:((match(ip_event, annotation[ , 'EX_IN'])) - 1), 'EX_IN']))]
  rs_exon <- s_exon[!is.na(match(s_exon, annotation[((match(ip_event, annotation[ , 'EX_IN'])) + 1):nrow(annotation), 'EX_IN']))]
  
  #required in connected exons for dividing metafeatures. It finds exonID of the two exons sharing boundary with ip event.
  index1 <- which(annotation$V5 == (annotation$V4[annotation$EX_IN == ip_event]-1))
  ls_exon.connected.exons <- annotation[index1, 'EX_IN']
  index2 <- which(annotation$V4 == (annotation$V5[annotation$EX_IN == ip_event]+1))
  rs_exon.connected.exons <- annotation[index2, 'EX_IN']
  if(length(ls_exon.connected.exons) > 1) ls_exon.connected.exons <- ls_exon.connected.exons[length(ls_exon.connected.exons)]
  if(length(rs_exon.connected.exons) > 1) rs_exon.connected.exons <- rs_exon.connected.exons[1]

  #finding max number of transcripts that can be formed
  if(length(ls_exon) == 0  && length(rs_exon) == 0){
    warning(ip_event, " does not have any flanking exons")
    } else if(length(ls_exon) == 0){
      n_ts <- length(rs_exon)
      } else if(length(rs_exon) == 0){
        n_ts <- length(ls_exon)
        } else{
          n_ts <- length(ls_exon)*length(rs_exon)
        }

 #initializing seed exon matrix. It will contain all possible combinations of left and right seed exons.
 se <- matrix(nrow=2, ncol=n_ts, dimnames=list(c('lse', 'rse')))

  # making combinations of left and right seed exons.
  if(length(ls_exon > 0) & length(rs_exon > 0)){
    se.p <- expand.grid('lse' = ls_exon, 'rse' = rs_exon)
    se[1, ] <- as.vector(se.p[, 1])
    se[2, ] <- as.vector(se.p[, 2])
    } else if(length(ls_exon) == 0){
      se[2, ] <- as.vector(rs_exon)
      } else if(length(rs_exon) == 0){
        se[1, ] <- as.vector(ls_exon)
        }
      s.exon <- list('se'=se, 'ls_exon'=ls_exon, 'rs_exon'=rs_exon, 'rs_exon.connected.exons' = rs_exon.connected.exons, 'ls_exon.connected.exons' = ls_exon.connected.exons)
      return(s.exon)
}

#' Title
#'
#' @param ep_event 
#' @param cmm 
#'
#' @return

.seed.exon.ep.excl <- function(ep_event = ep_event, cmm = cmm){
  # list of exons
  exon.list <- colnames(cmm)[grep("EX", colnames(cmm))]

  # finding seed skipping junction(s) to ep event
  s.mf <- rownames(cmm)[(cmm[,ep_event] == 0.5)]

  if(is.null(length(s.mf))){
    s.exon <- list('se'=NULL, 'n_ts'=NULL, 'ls_exon'=NULL, 'rs_exon'=NULL)
    warning('No skipping junction found for ', ep_event)
  } else{
  
    #seed exons for ep event
    s.exon <- lapply(1:length(s.mf), function(x) colnames(cmm)[cmm[s.mf[x], exon.list] == 2]) # exons associated with skipping junctions
    s.exon <- s.exon[lapply(s.exon,length) > 0] # removing empty elements
    s.exon <- unique(unlist(s.exon))
  
    # dividing seed exons into left and right seed exon
    ls.exons <- exon.list[1:(match(ep_event, exon.list) - 1)]
    rs.exons <- exon.list[(match(ep_event, exon.list) + 1):length(exon.list)]
    ls_exon <- ls.exons[!is.na(match(ls.exons, s.exon))] #dividing seed exons into left seed
    rs_exon <- rs.exons[!is.na(match(rs.exons, s.exon))] #dividing seed exons into right seed
  
  
    ## preparing matrix/list of left and right seed exon for each transcript
    if(length(ls_exon) == 0  && length(rs_exon) == 0){
      warning(ep_event, " does not have any flanking exons")
      } else if(length(ls_exon) == 0){
        n_ts <- length(rs_exon)
        } else if(length(rs_exon) == 0){
          n_ts <- length(ls_exon)
          } else{
            n_ts <- length(ls_exon)*length(rs_exon)
          }
  
    #initializing seed exon matrix. It will contain all possible combinations of left and right seed exons.
    se <- matrix(nrow=2, ncol=n_ts, dimnames=list(c('lse', 'rse')))
    
    # making combinations of left and right seed exons.
    if(length(ls_exon > 0) & length(rs_exon > 0)){
      se.p <- expand.grid('lse' = ls_exon, 'rse' = rs_exon)
      se[1, ] <- as.vector(se.p[, 1])
      se[2, ] <- as.vector(se.p[, 2])
      } else if(length(ls_exon) == 0){
        se[2, ] <- as.vector(rs_exon)
        } else if(length(rs_exon) == 0){
          se[1, ] <- as.vector(ls_exon)
        }
        s.exon <- list('se'=se, 'n_ts'=n_ts, 'ls_exon'=ls_exon, 'rs_exon'=rs_exon)
  }
      return(s.exon)

    }

#' Title
#'
#' @param ip_event 
#' @param cmm 
#' @param annotation 
#'
#' @return
.seed.exon.ip.excl <- function(ip_event = ip_event, cmm = cmm, annotation = annotation){
  exon.list <- colnames(cmm)[grep("EX", colnames(cmm))]

  # finding seed skipping junction(s) to ip event
  s.mf <- rownames(cmm)[(cmm[ , ip_event] == 0.5)]

  #seed exons for ip event
  s.exon <- lapply(1:length(s.mf), function(x) colnames(cmm)[cmm[s.mf[x], exon.list] == 2]) # exons associated with flanking junctions
  s.exon <- s.exon[lapply(s.exon,length) > 0] # removing empty elements
  s.exon <- unique(unlist(s.exon))

  # dividing seed exons into left and right seed exon
  annotation <- annotation[order(annotation$V4), ]
  ls_exon <- s.exon[!is.na(match(s.exon, annotation[1:((match(ip_event, annotation[ , 'EX_IN'])) - 1), 'EX_IN']))]
  rs_exon <- s.exon[!is.na(match(s.exon, annotation[((match(ip_event, annotation[ , 'EX_IN'])) + 1):nrow(annotation), 'EX_IN']))]
  
  #required in connected exons for dividing metafeatures. It finds exonID of the two exons sharing boundary with ip event.
  index1 <- which(annotation$V5 == (annotation$V4[annotation$EX_IN == ip_event]-1))
  ls_exon.connected.exons <- annotation[index1, 'EX_IN']
  index2 <- which(annotation$V4 == (annotation$V5[annotation$EX_IN == ip_event]+1))
  rs_exon.connected.exons <- annotation[index2, 'EX_IN']
  if(length(ls_exon.connected.exons) > 1) ls_exon.connected.exons <- ls_exon.connected.exons[length(ls_exon.connected.exons)]
  if(length(rs_exon.connected.exons) > 1) rs_exon.connected.exons <- rs_exon.connected.exons[1]

  #finding max number of transcripts that can be formed
  if(length(ls_exon) == 0  && length(rs_exon) == 0){
    warning(ip_event, " does not have any flanking exons")
    } else if(length(ls_exon) == 0){
      n_ts <- length(rs_exon)
      } else if(length(rs_exon) == 0){
        n_ts <- length(ls_exon)
        } else{
          n_ts <- length(ls_exon)*length(rs_exon)
        }

  #initializing seed exon matrix. It will contain all possible combinations of left and right seed exons.
  se <- matrix(nrow=2, ncol=n_ts, dimnames=list(c('lse', 'rse')))
  
  # making combinations of left and right seed exons.
  if(length(ls_exon > 0) & length(rs_exon > 0)){
    se.p <- expand.grid('lse' = ls_exon, 'rse' = rs_exon)
    se[1, ] <- as.vector(se.p[, 1])
    se[2, ] <- as.vector(se.p[, 2])
    } else if(length(ls_exon) == 0){
      se[2, ] <- as.vector(rs_exon)
      } else if(length(rs_exon) == 0){
        se[1, ] <- as.vector(ls_exon)
      }
      s.exon <- list('se' = se, 'ls_exon'=ls_exon, 'rs_exon'=rs_exon, 'rs_exon.connected.exons' = rs_exon.connected.exons, 'ls_exon.connected.exons' = ls_exon.connected.exons)
      return(s.exon)
    }
