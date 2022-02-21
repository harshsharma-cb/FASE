#' Connected Exons using igraph
#'
#' @import Matrix
#' @import igraph
#' @import parallel
#' @import methods
#' 
#' @param cmm 
#' @param event 
#' @param ls_exon 
#' @param rs_exon 
#' @param s.exon 
#' @param eventtype 
#' 
#' @return

connected.exons.igraph <- function(cmm = cmm, event = event, ls_exon = ls_exon, rs_exon = rs_exon, s.exon = NULL, eventtype = NULL){
  #making all FJs as 1 and the rest 0
  cmm <- (cmm == 2)*1
  #here i'll take both junctions and exons in adj mat
  adj.mat <- t(cmm) %*% cmm
  adj.mat <- as(adj.mat, "TsparseMatrix")
  #making adjacency matrix an upper right triangle matrix
  adj.mat[lower.tri(adj.mat, diag = TRUE)] <- 0
  igraph.result <- igraph::graph_from_adjacency_matrix(adj.mat, mode=c("directed"))
  
  #finding exons left and right to event
  if(eventtype == 'ep.incl' || eventtype == 'ep.excl'){
    if(length(ls_exon) > 0) ls_mf <- rownames(adj.mat)[1:((match(event, rownames(adj.mat))) - 1)] else ls_mf <- NULL
    if(length(rs_exon) > 0) rs_mf <- rownames(adj.mat)[((match(event, rownames(adj.mat))) + 1):nrow(adj.mat)] else rs_mf <- NULL
    } else{
      if(length(ls_exon) > 0) ls_mf <- rownames(adj.mat)[1:(match(s.exon$ls_exon.connected.exons, rownames(adj.mat)))] else ls_mf <- NULL
      if(length(rs_exon) > 0) rs_mf <- rownames(adj.mat)[(match(s.exon$rs_exon.connected.exons, rownames(adj.mat))):nrow(adj.mat)] else rs_mf <- NULL
    }
  
  #for left seed exons
  lse.struct <- NULL
  if(length(ls_exon) > 0){
    if((!is.null(ls_exon) && (length(ls_mf) > 1))){
      lse.struct <- lapply(1:length(ls_exon), function(x){
        if(!is.na(match(ls_exon[x], as_ids(V(igraph.result))))){
          lse.struct.n <- lapply((1:match(ls_exon[x], ls_mf)), function(z) all_simple_paths(graph = igraph.result, from = z, to = ls_exon[x]))
          lse.struct.n <- lse.struct.n[lengths(lse.struct.n) != 0]
          #removing subsets of long paths
          if(length(lse.struct.n) > 0){
            #extracting exonIDs from paths.
            for(l in 1:length(lse.struct.n)){
              if(!is.null(lse.struct.n[[l]])) 
                lse.struct.n[[l]] <- lapply(1:length(lse.struct.n[[l]]), function(x) as_ids(lse.struct.n[[l]][[x]]))
            }
            lse.struct.n <- unlist(lse.struct.n, recursive = F) #converting list of lists to list

            #removing sub sequences
            lse.struct.n.vec <- sapply(lse.struct.n, function(m) paste0(";", paste(m, collapse = ";"), ";"))
            # check each element, if this is found somewhere in any other element
            if(length(lse.struct.n.vec) > 1000){
              lse.struct.n <- NULL
              #paste0(event, " has too many structures, possibly due to 3'/5' Alternative Splicing or gene fusion.")
              return()
              } else{
                lse.struct.n <- lse.struct.n[sapply(seq_along(lse.struct.n.vec), function(m) !any(grepl(lse.struct.n.vec[m], lse.struct.n.vec[-m], fixed = TRUE)))] ##limiting step
                lse.struct.n <- lse.struct.n[lengths(lse.struct.n) != 0]
                }
          } else lse.struct.n <- NULL
        }
      })
    }
  }
  #naming sub-lists by seed exon
  if(!is.null(lse.struct)) names(lse.struct) <- ls_exon
  
  
  # for right seed exons
  rse.struct <- NULL
  if(length(rs_exon) > 0){
      rse.struct <- lapply(1:length(rs_exon), function(x){
        if(!is.null(rs_exon[x]) && (length(rs_mf) > 1) && !is.na(match(rs_exon[x], as_ids(V(igraph.result))))){
          rse.struct.n <- lapply((match(rs_exon[x], colnames(adj.mat))):ncol(adj.mat), function(z) all_simple_paths(graph = igraph.result, from = rs_exon[x], to = z))
          rse.struct.n <- rse.struct.n[lengths(rse.struct.n) != 0]
          #removing subsets of long paths
          if(length(rse.struct.n) > 0){
            for(l in 1:length(rse.struct.n)){
              if(!is.null(rse.struct.n[[l]])) 
                rse.struct.n[[l]] <- lapply(1:length(rse.struct.n[[l]]), function(x) as_ids(rse.struct.n[[l]][[x]]))
            }
            rse.struct.n <- unlist(rse.struct.n, recursive = F) #converting list of lists to list
            rse.struct.n.vec <- sapply(rse.struct.n, function(m) paste0(";", paste(m, collapse = ";"), ";"))
            # check each element, if this is found somewhere in any other element
            if(length(rse.struct.n.vec) > 1000){
              rse.struct.n <- NULL
              #paste0(event, " has too many structures, possibly due to 3'/5' Alternative Splicing or gene fusion.")
              return()
              } else{
            rse.struct.n <- rse.struct.n[sapply(seq_along(rse.struct.n.vec), function(m) !any(grepl(rse.struct.n.vec[m], rse.struct.n.vec[-m], fixed = TRUE)))] ##limiting step
            
            rse.struct.n <- rse.struct.n[lengths(rse.struct.n) != 0]
                }
            } else rse.struct.n <- NULL
          }
          })
  }
  if(!is.null(rse.struct)) names(rse.struct) <- rs_exon
  return(list( 'lse.struct' = lse.struct, 'rse.struct' = rse.struct))
}
