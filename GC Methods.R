########################################################
#
#  Methods for working with the garbage code death certificates from the IHME
#
#     -John Gosink
#      7/29/14
#

require(edgeR)



# -----------------------------------------
# Maps y onto x in the sense that it returns a list of of y corresponding
# to sort(unique(x)) such that the y's are the first y's seen in sort(x)
map <- function(x, y) {
	temp   <- data.frame(x=x, y=y, stringsAsFactors=FALSE)
	result <- ddply(temp, .(x), summarize, y[1])[,2]
	return(result)
}

#str(map(x=c(1:6,1:6), y=LETTERS[c(1:12)]))



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
	if (drop.empty) {
		x <- x[!is.na(x)]
		x <- x[nchar(x)>0]
	}		
	x  <- toupper(x)
	return(x)
}


# -----------------------------------------
# Directed graph from the death certificate MCD list.  The assumption is that the
# input mcdCols are in order, the first term being the actual cause of death and 
# terms further to the right are increasingly distally antecedant.
# To do:  should I initialize a transition matrix with all possibilites?   Nah..
getMCDTransitionMatrix <- function(x, priorMethod='Good-Turing', finalCause="gs_text",
		mcdCols=c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51"),
		cofactors=c('gender', 'ageGroup', 'smokingGroup'), ...) {

	x[,finalCause] <- make.names(x[,finalCause])
	finalVec       <- sort(unique(x[, finalCause]))
	icdVec         <- sort(unique(cleanICDs(unlist(x[,mcdCols]), drop.empty=TRUE)))

	numInput       <- nrow(x)
	nRows          <- length(icdVec)
	nCols          <- length(finalVec)
	QnBayes        <- matrix(rep(0, nRows*nCols), nrow=nRows, ncol=nCols, dimnames=list(icdVec, finalVec))
	Q              <- QnBayes  # belay this for now
	
	for (i in 1:numInput) {
		rowVec    <- cleanICDs(x[i,mcdCols], drop.empty=TRUE)   # Delete empty cells
		if (length(rowVec)==0) next;
		for (j in 1:length(rowVec)) {
			cause  <- rowVec[j]
			result <- x[i,finalCause]
			QnBayes[cause, result] <- QnBayes[cause, result] + 1  
		}
	}

	# Now add on the transition matrices for the cofactors (eg. gender, ageGroup & etc.)
	# TO DO:  THINK ABOUT HOW TO DEAL WITH MISSING INFORMATION OR 'UNKNOWN' CATEGORIES.
	#	     THESE MAY REALLY BE SKEWED BY JUST A FEW RANDOM CASES IN THE TRAINING DATA.
	for (aCofactor in cofactors) {
#		Qcof    <- getCofactorTransitionMatrix(x, finalCause, cofactor=aCofactor, method=priorMethod)
		Qcof    <- table(x[,aCofactor], x[,finalCause])
		rownames(Qcof) <- cleanICDs(paste(aCofactor, rownames(Qcof), sep='.'))  # Do this for consistency later on
#		rownames(Qcof) <- paste(aCofactor, rownames(Qcof), sep='.')
		QnBayes <- rbind(QnBayes, Qcof)   # how to confirm the same columns and column order????
	}

	# codVec is the counts of each final cause of death by antecedant cause
	codVec <- rowSums(QnBayes)
	
	# Need to add a background or prior value on the transitions that are zero.
	# Also, transition matrix rows must sum to one.
	nRows <- nrow(QnBayes)
	for (i in 1:nRows) {
		if (sum(QnBayes[i,])==0) QnBayes[i,] <- 1/nCols
		QnBayes[i,] <- priorFunction(x=as.numeric(QnBayes[i,]), method=priorMethod, constant=.Machine$double.eps)		
#		QnBayes[i,] <- priorFunction(x=as.numeric(QnBayes[i,]), method=priorMethod, constant=1/(sum(QnBayes[i,])))
		QnBayes[i,] <- QnBayes[i,] / sum(QnBayes[i,])
	}

	probObj  <- list(Q=Q, QnBayes=QnBayes, icdVec=icdVec, codVec=codVec, dimQ=dim(Q), numInputExamples=numInput,
					priorMethod=priorMethod)
	return(probObj)
}




# -----------------------------------------
# A helper function for the transition matrix to account for events 
# that weren't quite observed.  Ie a zero in an observation doesn't
# mean that it can't happen, just that we didn't observe it yet...
# I've kinda implemented goodTuring from library(edgeR), but am *not quite*
# using it correctly...
# 10/18/14 - after a *lot* of trial and error I see that the Good-Turing method doesn't
#            always put a small number into the zero cells... urk.
priorFunction <- function(x, method='Good-Turing', constant) {
	zeroInx <- (x == 0)
	nZero   <- sum(zeroInx)
	
	if (method=='square.min') {  	
		fill       <- sqrt(min(x[x > 0], na.rm=TRUE)) / nZero
		x[zeroInx] <- fill
		if (!missing(constant)) x <- x + constant
		return(x)
		
	} else if (method=='Good-Turing') {
		fill       <- (goodTuring(x)$P0) / nZero   # not quite right...
		x[zeroInx] <- fill
		if (!missing(constant)) x <- x + constant
		return(x)
		
	} else {
		if (!missing(constant)) x <- x + constant
		return(x)
	}
}


# -----------------------------------------
# A helper function to visualize the transition matrix
visualizeTransitionMatrix <- function(Q) {
	topYLabs <- sample(rownames(Q), 30)
	xLabs    <- colnames(Q)
	yLabs    <- rownames(Q)
	yLabs[! yLabs %in% topYLabs] = ''
	
	heatmap(Q, scale='column', col=gray(32:0/32),  margins=c(5,5),
		labRow=yLabs, labCo=xLabs)
	bringToTop()
}





# -----------------------------------------
# probObj    - the probObj obj with  list(Q, QnBayes, codVec, n)
# evidence   - a vector of one or more MCD observations, ie  c("K709", "I608", "K550", "M259")
bayesCompute <- function(evidence, probObj, ...) {
	evidence          <- cleanICDs(evidence, drop.empty=TRUE)
	goodEvInx         <- evidence %in% names(probObj$codVec)
	dropEv            <- evidence[!goodEvInx]
	if (length(dropEv) > 0) warning(paste("Couldn't find transition probabilities for", paste(dropEv, collapse=', ')))
	evidence          <- evidence[goodEvInx]
	obsVec            <- probObj$codVec * 0
	for (i in 1:length(evidence)) {  obsVec[evidence[i]]  <-  obsVec[evidence[i]] + 1 }
	
	# 10/15/14 meeting with Abie indicates this is Naive Bayes and I should move forward with this.
	naiveBayesVec  <- obsVec %*% log(probObj$QnBayes) 
	naiveBayesVec  <- exp(naiveBayesVec[1,])
	naiveBayesVec  <- naiveBayesVec / sum(naiveBayesVec)
	
	outcomes <- list(evidence=evidence, dropEv=dropEv, naiveBayesVec=naiveBayesVec)
	return(outcomes)
}


showBayesCompute <- function(bayesDatObj) {
	old.par <- par(mfrow=c(2,1))
#	barplot(head(sort(bayesDatObj$eStepVec[bayesDatObj$eStepVec >0], decreasing=TRUE), 10),
#			main='Top 10 predictions using the "Stepwise" approach', cex.names=0.8); grid()
	barplot(head(sort(bayesDatObj$naiveBayesVec[bayesDatObj$naiveBayesVec >0], decreasing=TRUE),10),
			main='Top 10 predictions using the Naive Bayes (1-gram) approach', cex.names=0.8); grid()
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





# -----------------------------------------
# HAVEN'T STARTED THIS YET.   NEED TO BREAK IT OUT AS IT'S OWN FUNCTION.
# EVERYTHING HERE IS JUST CUT AND PASTED FROM MY SCRIPTING LINES..
#  TO DO:  
#    0. Correctly figure out how to group and panel arbitrary data
#    0.1 correctly incorporate cofactors into figures
#    1. Ordered causality
#    2. Visualize as a complete Naive Bayes  (eg. one step transition)
visualizeCausalChains <- function() {

	# A plot for each final cause of death  
	icdGuess(dat[sample(nrow(dat), 30), 'codcau21'])  # Here is a sample of ICD descriptions
	codBucket  <- with(temp2, unique(MCD[ grepl('_gc', icd_name) & causeOrder==1]))    # one way to search
	idBucket   <- with(temp2, unique(sid[ gs_text=='Homicide' ]))    # by gs_text
	descBucket <- with(temp2, unique(MCD[ grepl('infarction', MCDfull, ignore.case=TRUE) ]))  # more descriptive
	idBucket   <- with(temp2, unique(sid[ (MCD %in% descBucket) ]))
	temp3      <- subset(temp2, subset = ((sid %in% idBucket) & (as.numeric(MCD) > 1)))

	pdf(file='MCD Paths from the Mexico Data Set.pdf', width=11, height=8.5)
	for (aText in sort(unique(temp2$gs_text))) {
		idBucket   <- with(temp2, unique(sid[ gs_text==aText ]))    # everyone who has this gs_text
		temp3      <- subset(temp2, subset = ((sid %in% idBucket) & (as.numeric(MCD) > 1)))

		yLabels <- as.character(count(temp3$MCD)$x)
		yLabels <- sample(yLabels, size=min(5,length(yLabels)))  
		yAt     <- match(yLabels, levels(temp3$MCD)) 

		ageGroups   <- cut(temp3$age, breaks=c(-2,2,16,50,200), labels=c('Infant', 'Child', 'Adult', 'Senior'))  
		colGroups   <- c('red', 'blue')[as.numeric(as.factor(temp3$gender))]  # display.brewer.all()
		pointGroups <- c(49:52)[as.numeric(ageGroups)]                        # pch=48:57  are '0' -> '9'

		aPlot <- xyplot(as.numeric(MCD) ~ jitter(causeOrder, amount=0.1) | icd_name, groups=sid, data=temp3, type=c('b','g'),
						main=paste('All patients with "gs_text" (died from):', aText, '\n(panel by icd_name)'),
						sub=list('\nBlue=male, Red=Female\n1=Infant, 2=Child, 3=Adult, 4=Senior', cex=0.7),
						col=colGroups, pch=pointGroups, cex=1,
						scales=list(y=list(at=yAt, labels=yLabels, relation='same')),
						xlab='< < <== Direction of Causality', ylab='ICD-10 Code'); #bringToTop()
		print(aPlot)
	}					

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







