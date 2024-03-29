#24.09.2020. modified:27.10.2020, 02.11.2020
#Transcript Concentration
#After running TS_6.2 and obtaining TS
#transtruct = output of function transtruct

#main function
#' Transcript Concentration
#' 
#' @description Transcript concentration gives the abundance of transcripts expressed in different samples and conditions. It uses transcript structure generated by transtruct method and superimposes expression values of corresponding meta-features.
#'
#' @import MASS
#'
#' @param transtruct transcript structure for EP/IP inclusion and exclusion for a gene.
#' @param designM design matrix for the two conditions in question. It should be a binary matrix where 1 represents samples of the condition and 0 represents the samples of other condition.
#'
#' @return transconc returns three matrices: \cr
#' \enumerate{
#'          \item transconc.samples: transcript concentration by samples.
#'          \item transconc.condition: transcript concentration by condition.
#'          \item transconc.TS: structures of transcripts and concentration of their meta-features, used for evaluating transcript abundance (transconc combines all structures by transtruct and retains only unique structures).
#' }
#' @export

transconc <- function(transtruct, designM){
	#require('MASS')
	transtruct <- .r.dup.ts(transtruct = transtruct)
	Y <- t(transtruct$expression) #samples*metafeatures
	G <- transtruct$ts #transcripts*metafeatures
	##colnames of Y and G should be in exactly same order
	T <- Y %*% ginv(G) #ginv(G) is samples*transcripts

	G1 <- G; T1 <- T
	 ts.rownames <- paste0('TS', 1:nrow(G1))
  rownames(G1) <- ts.rownames
  colnames(T1) <- ts.rownames
  T1 <- .r.a.bias(ts = transtruct$ts, tc = T1)
  #finding transcript concentration in different conditions:
  T2 <- t(T1)%*%designM
  T2 <- sweep(T2, 2, colSums(designM), FUN = '/')
  
	return(list('transconc.samples' = T1, 'transconc.conditions' = T2, 'TS' = transtruct$ts))
}

#Removing duplicate transcript(s), if any
#' Removing duplicate transcript
#'
#' @param transtruct 
#'
#' @return
.r.dup.ts <- function(transtruct){
	ts <- rbind(transtruct[[1]], transtruct[[2]])
	up<- lapply(1:nrow(ts),function(x) paste(ts[x,],collapse=""))
	uup<- unique(up)
	index_up<- unlist(lapply(1:length(uup),function(x) pmatch(uup[x],up)))
	ts<- ts[index_up,,drop=FALSE]
	#renaming rownames
	#ts.rownames <- sprintf('TS%03d', 1:nrow(ts))
	ts.rownames <- paste0('TS', 1:nrow(ts))
  	rownames(ts) <- ts.rownames
	return(list('ts' = ts,  'expression' = transtruct$expression))
}

#Removing algebraic bias in full length transcript(s)
#' Removing algebraic bias
#'
#' @param
#'
#' @return
.r.a.bias <- function(ts, tc){
	exons <- grep(colnames(ts), pattern = 'EX', fixed = T)
	nexons <- rowSums(ts[ , exons])
	flt <- names(which(nexons == max(nexons)))
	for(x in 1:length(flt)) {
		flt.sign <- sign(tc[,flt[x]])
		tc.temp <- abs(tc[ , flt[x]]) - rowMeans(tc)
		###if(any(tc[ , flt[x]]) < 0) temp <- (tc[ , flt[x]])*-1
		tc[ , flt[x]] <- flt.sign*tc.temp
		}
	return('T1' = tc)
}
