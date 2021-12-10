#' Heatscore using Splice Index.
#' 
#' @description heatscore uses expression of ranked EP and IP events 
#'
#' @param ep ExonPointer file containing cassette exon event ranking. 
#' @param ip IntronPointer file containing intron retention event ranking.
#' @param ep.exp ExonPointer cassete exon expression file.
#' @param ip.exp IntronPointer intron retention expression file.
#'
#' @return Heat Score/Splice Index required by HotNet software.
#' @export

heatscore <- function(ep, ip, ep.exp, ip.exp){
  #calculating event heatscore after combining EP and IP expression
  epip.exp <- .epip.exp(ep = ep, ep.exp = ep.exp, ip = ip, ip.exp = ip.exp)
  
  #calculating gene heatscore and normalizing it
  hs <- aggregate(evenths ~ genes, data = epip.exp, sum)
  hs$normHeatscore <- unlist(lapply(1:nrow(hs), function(x) hs[x,3] <- (hs[x , 2] - min(hs[ , 2]))/(max(hs[ , 2]) - min(hs[ , 2])) ))
  hs <- hs[ , c(1,3)]
  colnames(hs) <- c('gene', 'heatscore')
  return(hs)
}

#' Similarity Score for network edge weight
#'
#' @param ep ExonPointer file containing cassette exon event ranking. 
#' @param ip IntronPointer file containing intron retention event ranking.
#' @param ep.exp ExonPointer cassete exon expression file.
#' @param ip.exp IntronPointer intron retention expression file.
#' 
#' @return Similarity Score gives cosine similarity between the splice index of two genes. It can be used as edge weight for an unweighted network.
#' @export 
#'
simscore <- function(ep, ep.exp, ip, ip.exp){
  epip.exp <- .epip.exp(ep = ep, ep.exp = ep.exp, ip = ip, ip.exp = ip.exp)
  epip.exp$eventhssq <- epip.exp$evenths^2
  
  ##Aggregating heat scores by genes
  hss <- aggregate(evenths ~ genes, data = epip.exp, sum)
  hssq <- aggregate(eventhssq ~ genes, data = epip.exp, sum)
  
  ##Final hs
  heats <- cbind(hss, hssq[ , 2])
  heats[ , 3] <- (heats[ , 3]) ^ (1/2)
  rownames(heats) <- heats[ , 1]; heats <- heats[ , -1]
  colnames(heats) <- c('heat', 'rootsumsqheat')
  
  ##Generating similatity matrix
  sim.score <- cbind(expand.grid("gene1" = rownames(heats), "gene2" = rownames(heats)), expand.grid(heats[ , 1], heats[ , 1]), expand.grid(heats[ , 2], heats[ , 2]))
  sim.score <- cbind(sim.score[ , 1:2], (sim.score[ , 3] * sim.score[ , 4]) / (sim.score[ , 5] * sim.score[ , 6]))
  colnames(sim.score) <- c("gene1", "gene2", "sim.score")
  sim.score[,3] <- sim.score[,3] - min(sim.score[,3])
  sim.score[,3] <- ((sim.score[,3] - min(sim.score[,3])) / (max(sim.score[,3]) - min(sim.score[,3])))
  return(sim.score)
}

#' Prepare files for hs and simscore
#'
#' @param ep 
#' @param ip 
#' @param ep.exp 
#' @param ip 
#' @param exp 
#'
#' @return

.epip.exp <- function(ep, ip, ep.exp, ip.exp){
  ##Preparing ep and ip files
  rownames(ep.exp) <- ep.exp[ , 1]; ep.exp <- ep.exp[ , -1]
  rownames(ip.exp) <- ip.exp[ , 1]; ip.exp <- ip.exp[ , -1]
  
  ##Combining ep and ip expression and finding heat score
  epip.exp <- rbind(ep.exp, ip.exp)
  epip.exp$evenths <- rowMeans(epip.exp[1:(ncol(epip.exp)-6)]) - rowMeans(epip.exp[(ncol(epip.exp)-6): ncol(epip.exp)])
  epip.exp$genes <- c(as.vector(ep$gene), as.vector(ip$gene))
  return(epip.exp)
}
