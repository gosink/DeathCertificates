########################################################
#
#  Methods for working with the garbage code death certificates from the IHME
#
#     -John Gosink
#      7/29/14
#

# ICD codes are almost, but not quite right:
# ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2015/ICD10CM_FY2015_code_descriptions.zip
# This is depredcated as I find that the CMS version of the ICD10 list is correct
makeICDList <- function(detailDat) {
	
	mapInfo  <- NULL
	temp     <- ddply(detailDat, .(gs_code55), summarize, description=unique(gs_text55))
	inx      <- !(temp$description %in% mapInfo)
	temp2    <- temp$description[inx];   names(temp2) <- temp$gs_code55[inx]
	mapInfo  <- c(mapInfo, temp2)
	
	temp     <- ddply(detailDat, .(gs_code46), summarize, description=unique(gs_text46))
	inx      <- !(temp$description %in% mapInfo)
	temp2    <- temp$description[inx];   names(temp2) <- temp$gs_code46[inx]
	mapInfo  <- c(mapInfo, temp2)
	
	temp     <- ddply(detailDat, .(gs_code34), summarize, description=unique(gs_text34))
	inx      <- !(temp$description %in% mapInfo)
	temp2    <- temp$description[inx];   names(temp2) <- temp$gs_code34[inx]
	mapInfo  <- c(mapInfo, temp2)
		
	return(mapInfo)
}


# -----------------------------------------
# ICD10 data as follows:
# 00003 A001    1 Cholera due to Vibrio cholerae 01, biovar eltor              Cholera due to Vibrio cholerae 01, biovar eltor
icdGuess <- function(x, icdVec) {
	if (missing(icdVec)) icdVec <- get(grep('icd.*Vec', ls(envir=globalenv()), value=TRUE, perl=TRUE, ignore.case=TRUE)[1])
	if (length(x) > 1000) {
		numSeq <- length(x)
		cat('  I expect it will take about', round(0.009*numSeq), 'seconds to process your', numSeq, 'entries...\n')
	}
	xMod <- cleanICDs(x, drop.empty=FALSE)  #toupper(gsub(' ','', x))
	ans  <- icdVec[xMod]
	ans[nchar(xMod)==0] <- ''
	inx  <- which(is.na(ans))
	for (i in inx) {
		nextCode <- which(xMod[i] < names(icdVec))
		ans[i]   <- paste('?', icdVec[nextCode-1][1], '|', icdVec[nextCode][1])
		xMod[i]  <- paste(names(icdVec)[nextCode-1][1], names(icdVec)[nextCode][1], sep='|')
	}
	
	names(ans)                    <- x
	attr(ans, 'closest.ICD.code') <- xMod
	
	return(ans)
}


# -----------------------------------------
# ICD10 data as follows:
readICDList <- function(aFile, fileFormat=10) {
	if (fileFormat == 10) {
		icdDat <- read.fwf(file=aFile, widths=c(5, -1, 7, -1, 1, -1, 60, -1, 206), 
							header=FALSE, n=-1, stringsAsFactors=FALSE, 
							strip.white=TRUE)
		icdVec        <- icdDat[,4]
		names(icdVec) <- icdDat[,2]
	} else {
		icdDat <- read.fwf(file=aFile, widths=c(5, 226), colClasses='character',
							header=FALSE, n=-1, stringsAsFactors=FALSE, 
							strip.white=TRUE)
		icdVec        <- icdDat[,2]
		names(icdVec) <- icdDat[,1]
	}
	return(icdVec)
}



# -----------------------------------------
# Simple helper function
cleanICDs <- function(x, drop.empty=TRUE) {
	x  <- gsub(' ', '', (unlist(x)))  
	if (drop.empty)  x <- x[nchar(x)>0] 
	x  <- toupper(x)
	return(x)
}


# -----------------------------------------
# Directed graph from the death certificate MCD list.  The assumption is that the
# input mcdCols are in order, the first term being the actual cause of death and 
# terms further to the right are increasingly distally antecedant.
# To do:  should I initialize a transition matrix with all possibilites?   Nah..
getMCDTransitionMatrix <- function(x, mcdCols=c("causadef", "codcau11", "codcau21", 
												"codcau31", "codcau41", "codcau51"), ...) {

	numInput      <- nrow(x)
	icdVec        <- sort(unique(cleanICDs(unlist(x[,mcdCols]), drop.empty=TRUE)))
	n             <- length(icdVec)
	codVec        <- rep(0, length(icdVec))
	names(codVec) <- icdVec
	
	# Q is the returned transition matrix for each step in the causal chains.
	#   - The structure of Q is:   Q[ from,  to ]
	# Qglobal is the "1-step" transition matrix from any causal event to 
	# the terminal (causadef) event.
	Q <- Qglobal <- matrix(rep(0, n^2), nrow=n, ncol=n, dimnames=list(icdVec, icdVec))
	for (i in 1:numInput) {
		rowVec    <- cleanICDs(x[i,mcdCols], drop.empty=TRUE)   # Delete empty cells
		numCauses <- length(rowVec)
		codVec[rowVec[1]] <- codVec[rowVec[1]] + 1
		if (numCauses==1) {
			Qglobal[rowVec[1], rowVec[1]] <- Qglobal[rowVec[1], rowVec[1]] + 1			
			Q[rowVec[1], rowVec[1]]       <- Q[rowVec[1], rowVec[1]] + 1
		} else {
			for (j in numCauses:2) {
				Qglobal[ rowVec[j], rowVec[1] ] <- Qglobal[ rowVec[j], rowVec[1] ] + 1
				Q[ rowVec[j], rowVec[j-1] ]     <- Q[ rowVec[j], rowVec[j-1] ] + 1
			}
		}
	}
	
	# Transition matrix rows must sum to one
	for (i in 1:n) {
		if (sum(Q[i,])==0) Q[i,i] <- 1
		Q[i,] <- Q[i,] / sum(Q[i,])
		# Likewise for Qglobal
		if (sum(Qglobal[i,])==0) Qglobal[i,i] <- 1
		Qglobal[i,] <- Qglobal[i,] / sum(Qglobal[i,])
	}
	
	# codVec probabilities must sum to one
	codVec   <- codVec / sum(codVec)
	probObj  <- list(Q=Q, Qglobal=Qglobal, codVec=codVec, dimQ=dim(Q), numInputExamples=numInput)
	
	return(probObj)
}




# -----------------------------------------
# probObj    - the probObj obj with  list(Q, Qglobal, codVec, n)
# evidence   - a vector of one or more MCD observations, ie  c("K709", "I608", "K550", "M259")
bayesCompute <- function(evidence, probObj, ...) {
#	print(class(x)); 
	evidence          <- cleanICDs(evidence, drop.empty=TRUE)
	goodEvInx         <- evidence %in% names(probObj$codVec)
	dropEv            <- evidence[!goodEvInx]
	if (length(dropEv) > 0) warning(paste('Couldnt find transition probabilities for', paste(dropEv, collapse=', ')))
	evidence          <- evidence[goodEvInx]
	obsVec            <- probObj$codVec * 0
	obsVec[evidence]  <- 1

	eStepVec    <- obsVec %*% (probObj$Q %^% length(evidence))
	eStepVec    <- eStepVec[1,] / sum(eStepVec[1,])
	eGlobalVec  <- obsVec %*% probObj$Qglobal 
	eGlobalVec  <- eGlobalVec[1,] / sum(eGlobalVec[1,])
	
	outcomes <- list(evidence=evidence, dropEv=dropEv, eStepVec=eStepVec, eGlobalVec=eGlobalVec)
	return(outcomes)
}


showBayesCompute <- function(bayesDatObj) {
	old.par <- par(mfrow=c(2,1))
	barplot(head(sort(bayesDatObj$eStepVec[bayesDatObj$eStepVec >0], decreasing=TRUE), 10),
			main='Top 10 predictions using the "Stepwise" approach', cex.names=0.8); grid()
	barplot(head(sort(bayesDatObj$eGlobalVec[bayesDatObj$eGlobalVec >0], decreasing=TRUE),10),
			main='Top 10 predictions using the "Global" approach', cex.names=0.8); grid()
	bringToTop()
	cat('\nThe MCD trail is as follows:\n----------------------------\n')
	for (i in 1:length(bayesDatObj$evidence)) {
		cat(i, bayesDatObj$evidence[i], icdGuess(bayesDatObj$evidence[i]), '\n')
	}
	if (length(bayesDatObj$dropEv) > 0) {
		cat('\n\nTerms in the input evidence but not found in the transition matrix include:\n')
		cat('---------------------------------------------------------------------------\n')
		for (i in 1:length(bayesDatObj$dropEv)) {
			cat(i, bayesDatObj$dropEv[i], icdGuess(bayesDatObj$dropEv[i]), '\n')
		}	
	}
	cat('\n')
	par(old.par)
}



# =====================  D E B R I S  ===================== #


# -----------------------------------------
# Wrap a 'try' around accessing an element in a transition matrix, Q.
# If you look for a cell that doesn't exist in Q, then it will return
# by default a very small number, the square of the smallest value in Q.
safeMatLookup <- function(Q, a, b, cell.default) {
	if (missing(cell.default)) cell.default = (min(Q[Q > 0], na.rm=TRUE))^2
	if (class(try(Q[a,b], silent=TRUE)) == 'try-error') return(cell.default)
	return(Q[a,b])
}



# -----------------------------------------
# Q is the transition matrix
# x is the line from the standard data file
evaluateSubject <- function(x, Q, codCol="causadef", 
							mcdCols=c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51"), ...) {
#	print(class(x)); 
#	mcdDat    <- gsub(' ', '', x[c(codCol, mcdCols)])
#	mcdDat    <- toupper(mcdDat[nchar(mcdDat)>0])
	mcdDat    <- cleanICDs(x[c(codCol, mcdCols)], drop.empty=TRUE)
	numCauses <- length(mcdDat)
	finalProb <- 1
	if (numCauses == 1) {
		finalProb <- finalProb * safeMatLookup(Q, mcdDat[1], mcdDat[1])
	} else {
		for (i in length(mcdDat):2) {
#		cat(i, finalProb, '\n')
			finalProb <- finalProb * safeMatLookup(Q, mcdDat[i], mcdDat[i-1])
		}
	}
	
	return(c(numCauses=numCauses, finalProb=finalProb))
}



# -----------------------------------------
# Directed graph from the death certificate MCD list
reportBigrams <- function(x, method=c('all', 'ordered', 'root1')) {
	bigramDat <- NULL
	rootList  <- list()
	x         <- gsub(' ', '', x)    # strip white
	x         <- x[ length(x) > 0 ]  # drop missing cells

	# Root ordering with causadef being the special root (should just sort & unique at the end)
	rootList[[ x[1] ]] <- sort(unique(c(rootList[[ x[1] ]], x[2:length(x)])))
	
	# Ordered return list with nothing special about causadef (except that it's always a edge recipient)
	edgeCount <- 0
	edgeVec   <- list()
	for (i in 1:(length(x)-1)) {
		edge          <- paste(x[i], x[i+1], sep='_')
#		print(str(edgeVec))
		edgeVec[[edge]] <- ifelse(is.null(edgeVec[edge]), 1, edgeVec[[edge]]+1)
		edgeCount     <- edgeCount+1
	}
	
	retList <- list(rootList=rootList, edgeVec=edgeVec, edgeCount=edgeCount)
	return(retList)
}







