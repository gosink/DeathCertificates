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
# Grab all of the MCDs from a block of data and see if there
# is something obvious that distinguishes the groups.
# Useage:    tabulateTerms(dat)
tabulateTerms <- function(x, groups='gs_text', 
						mcdCols=c('causadef','codcau11','codcau21','codcau31','codcau41','codcau51')) {
						
	allGroups = sort(unique(x[,groups]))
	allTerms  = cleanICDs(unique(unlist(x[,mcdCols])))
	n = length(allTerms)
	p = length(allGroups)
	mat = matrix(rep(0, n*p), nrow=n, ncol=p, dimnames=list(allTerms, allGroups))
	for (aGroup in allGroups) {
		inx    = x[,groups] == aGroup
		values = table(cleanICDs(unlist(x[inx, mcdCols])))
		mat[names(values), aGroup] = values
	}
	mat = mat[order(mat[,1], decreasing=TRUE),]
	return(mat)
}





# -----------------------------------------
# Generate nGrams from incoming data
# Usage:   nGramFunction(x=c('Dog', 'Cat', 'Horse', NA, 'Flower'), n=2, sep='.')
nGramFunction <- function(x, n=1, sep='.', ...) {
	n    <- min(c(n, length(x)))
	if (n <= 1) return(x)
	n    <- n - 1
	temp <- NULL
	for (nn in 1:n) {
		for (i in 1:(length(x) - nn)) {
			temp <- c(temp, paste(x[i:(i+nn)], collapse=sep))
		}
	}
	return(c(x, temp))
}



# -----------------------------------------
# Directed graph from the death certificate MCD list.  The assumption is that the
# input mcdCols are in order, the first term being the actual cause of death and 
# terms further to the right are increasingly distally antecedant.
# To do:  should I initialize a transition matrix with all possibilites?   Nah..
getMCDTransitionMatrix <- function(x, priorMethod='Good-Turing', constant=.Machine$double.eps,
		finalCause="gs_text", nGrams=1, columnNormalize=FALSE,
		mcdCols=c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51"),
		cofactors=c('gender', 'ageGroup', 'smokingGroup'), ...) {

	x[,finalCause] <- make.names(x[,finalCause])
	finalVec       <- sort(unique(x[, finalCause]))
	icdVec         <- cleanICDs(unlist(t(x[,mcdCols])), drop.empty=TRUE)
	icdVec         <- nGramFunction(icdVec, n=nGrams, sep='.')
	icdVec         <- sort(unique(icdVec))

	numInput       <- nrow(x)
	nRows          <- length(icdVec)
	nCols          <- length(finalVec)
	QnBayes        <- matrix(rep(0, nRows*nCols), nrow=nRows, ncol=nCols, dimnames=list(icdVec, finalVec))
	
	for (i in 1:numInput) {
		rowVec    <- cleanICDs(x[i,mcdCols], drop.empty=TRUE)   # Delete empty cells
		rowVec    <- nGramFunction(rowVec, n=nGrams, sep='.')   # Add in tuples, if appropriate		
		if (length(rowVec)==0) next;
		QnBayes[rowVec, x[i,finalCause]] <- QnBayes[rowVec, x[i,finalCause]] + 1
	}
	
	# Because of the way I set up the nGram tuples, I may have a largish number of
	# rows in QnBayes at this moment that don't have any values.
	QnBayes <- QnBayes[rowSums(QnBayes) != 0,]
	
	# Now add on the transition matrices for the cofactors (eg. gender, ageGroup & etc.)
	# TO DO:  THINK ABOUT HOW TO DEAL WITH MISSING INFORMATION OR 'UNKNOWN' CATEGORIES.
	#	     THESE MAY REALLY BE SKEWED BY JUST A FEW RANDOM CASES IN THE TRAINING DATA.
	allCofRowNames <- NULL
	for (aCofactor in cofactors) {
#		Qcof     <- getCofactorTransitionMatrix(x, finalCause, cofactor=aCofactor, method=priorMethod)
		Qcof     <- table(x[,aCofactor], x[,finalCause])
		cofNames <- cleanICDs(paste(aCofactor, rownames(Qcof), sep='.'))  # Do this for consistency later on
		rownames(Qcof) <- cofNames
		allCofRowNames <- c(allCofRowNames, cofNames)
		QnBayes        <- rbind(QnBayes, Qcof)   
	}

	# codVec is the counts of each final cause of death by antecedant cause
	codVec  <- rowSums(QnBayes)
	rawMat  <- QnBayes
	
	# The training and test data sets may have vastly different sizes for the classification
	# populations.   Account for this by performing column-wise normalization.  Notice that 
	# I treat cofactor information separately from the other information as everyone will
	# have a value for every cofactor.   Not so with the mcdCols.
	if (columnNormalize) {
		cofInx            <- rownames(QnBayes) %in% allCofRowNames
		QnBayes[!cofInx,] <- apply(QnBayes[!cofInx,], 2, function(x) { x / sum(x) })
		QnBayes[cofInx,]  <- apply(QnBayes[cofInx,],  2, function(x) { x / sum(x) })
	}
	
	# Need to add a background or prior value on the transitions that are zero.
	# Also, transition matrix rows must sum to one.
	nRows <- nrow(QnBayes)
	for (i in 1:nRows) {
		if (sum(QnBayes[i,])==0) QnBayes[i,] <- 1/nCols
		QnBayes[i,] <- priorFunction(x=as.numeric(QnBayes[i,]), method=priorMethod, constant=constant)		
#		QnBayes[i,] <- priorFunction(x=as.numeric(QnBayes[i,]), method=priorMethod, constant=sum(QnBayes[i,])/ncol(QnBayes))
		QnBayes[i,] <- QnBayes[i,] / sum(QnBayes[i,])
	}

	probObj  <- list(QnBayes=QnBayes, rawMat=rawMat, icdVec=icdVec, codVec=codVec, 
					numInputExamples=numInput, priorMethod=priorMethod, constant=constant, 
					mcdCols=mcdCols, nGrams=nGrams, cofactors=cofactors, finalCause=finalCause, 
					columnNormalize=columnNormalize, createTime=Sys.time())
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
visualizeTransitionMatrix <- function(Q, ...) {
	sampleSize <- min(c(length(rownames(Q)), 30))
	topYLabs   <- sample(rownames(Q), sampleSize)
	xLabs      <- colnames(Q)
	yLabs      <- rownames(Q)
	yLabs[! yLabs %in% topYLabs] = ''
	
	heatmap(Q, scale='column', col=gray(32:0/32),  margins=c(5,5),
		labRow=yLabs, labCo=xLabs)
	bringToTop()
}




# -----------------------------------------
# probObj    - the probObj obj with  list(QnBayes, codVec, n, etc.)
# evidence   - a vector of one or more MCD observations, ie  c("K709", "I608", "K550", "M259")
bayesCompute <- function(evidence, probObj, nGrams=1, ...) {
	evidence          <- cleanICDs(evidence, drop.empty=TRUE)
	evidence          <- nGramFunction(evidence, n=nGrams, sep='.')

	goodEvInx         <- evidence %in% names(probObj$codVec)
	dropEv            <- evidence[!goodEvInx]
	if (length(dropEv) > 0) warning(paste("Couldn't find transition probabilities for", paste(dropEv, collapse=', ')))
	evidence          <- evidence[goodEvInx]
	if (length(evidence) == 0) 
		return(list(evidence=evidence, dropEv=dropEv, naiveBayesVec=colSums(probObj$QnBayes)/sum(probObj$QnBayes)))
	
	obsVec            <- probObj$codVec * 0
	for (i in 1:length(evidence)) {  obsVec[evidence[i]]  <-  obsVec[evidence[i]] + 1 }
	
	# 10/15/14 meeting with Abie indicates this is Naive Bayes and I should move forward with this.
	naiveBayesVec  <- obsVec %*% log(probObj$QnBayes) 
	naiveBayesVec  <- exp(naiveBayesVec[1,])
	naiveBayesVec  <- naiveBayesVec / sum(naiveBayesVec)
	QPartial       <- probObj$QnBayes[evidence, names(sort(naiveBayesVec, decreasing=TRUE))]   # interesting return obj
	outcomes <- list(evidence=evidence, dropEv=dropEv, naiveBayesVec=naiveBayesVec, QPartial=QPartial)
	return(outcomes)
}


showBayesCompute <- function(bayesDatObj) {
#	old.par <- par(mfrow=c(2,1))
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
	
	windows()
	visualizeTransitionMatrix(bayesDatObj$QPartial, Colv=bayesDatObj$naiveBayesVec) 
#	heatmap.2(bayesDatObj$QPartial)
	cat('\n')
#	par(old.par)
}





# -----------------------------------------
# HAVEN'T STARTED THIS YET.   NEED TO BREAK IT OUT AS IT'S OWN FUNCTION.
# EVERYTHING HERE IS JUST CUT AND PASTED FROM MY SCRIPTING LINES..
#  TO DO:  
#    0. Correctly figure out how to group and panel arbitrary data
#    0.1 correctly incorporate cofactors into figures
#    1. Ordered causality
#    2. Visualize as a complete Naive Bayes  (eg. one step transition)
visualizeCausalChains <- function(x, orderedCauses=c(11:6), colorGroups='gender', pchGroups='smokingGroup') {

	temp          <- reshape(x, idvar='sid', varying=orderedCauses, v.names='MCD', direction='long')
	temp$MCD      <- cleanICDs(temp$MCD, drop.empty=FALSE)
	temp$MCD      <- as.factor(temp$MCD)
	names(temp)   <- gsub('time', 'causeOrder', names(temp))
	temp[,colorGroups] <- factor(temp[,colorGroups])
	temp[,pchGroups]   <- factor(temp[,pchGroups])

	for (aText in sort(unique(temp$gs_text))) {
		idBucket   <- with(temp, unique(sid[ gs_text==aText ]))    # everyone who has this gs_text
		temp2      <- subset(temp, subset = ((sid %in% idBucket) & (as.numeric(MCD) > 1)))

		# Pick a sample of MCDs to decorate the y-axis
		yLabels <- as.character(count(temp2$MCD)$x)
		yLabels <- sample(yLabels, size=min(7, length(yLabels)))  
		yAt     <- match(yLabels, levels(temp2$MCD)) 

		colorChoices <- brewer.pal(n=8, 'Dark2')
		colGroups    <- colorChoices[as.numeric(temp2[,colorGroups])]  # display.brewer.all()
		pointGroups  <- c(49:62)[as.numeric(temp2[,pchGroups])]        # pch=48:57  are '0' -> '9'

		aPlot <- xyplot(as.numeric(MCD) ~ jitter(causeOrder, amount=0.01) | ageGroup, groups=sid, data=temp2, 
						type=c('b','g'),
						main=paste('All patients with "gs_text" (died from):', aText),
						col=colGroups, pch=pointGroups, cex=1,
						scales=list(y=list(at=yAt, labels=yLabels, relation='same'),
									x=list(at=(1:length(orderedCauses)), labels=colnames(x)[orderedCauses]), rot=45),
						xlab='Direction of Causality', ylab='ICD-10 Code',
						sub=list(paste('points=', pchGroups, '   colors=', colorGroups), cex=0.7))
		print(aPlot)
	}
}


# -----------------------------------------
# Tally up the counts of various pairings and make them into a matrix
#  counts     extra if I have a count of certain pairings.  (not implemented yet)
#  allVals    if provided, then force matrix to include thiese terms (and be square)
pairsToMatrix <- function(rowVals, colVals, counts, square=TRUE, allVals=NULL) {
	if (!is.null(allVals)) square=TRUE
	allVals <- sort(unique(c(rowVals, colVals, allVals)))
	if (length(rowVals) != length(colVals)) stop('rowVals and colVals need to be the same length')
	
	if (square) {
		n     <- length(allVals)
		mat   <- matrix(rep(0, n^2), nrow=n, ncol=n,
					dimnames=list(allVals, allVals))
	} else {
		nRows <- length(unique(rowVals))
		nCols <- length(unique(colVals))
		mat   <- matrix(rep(0, nRows*nCols), nrow=nRows, ncol=nCols,
					dimnames=list(sort(unique(rowVals)), sort(unique(colVals))))
	}
	
	for (i in 1:length(rowVals)) {
		mat[rowVals[i], colVals[i]] = mat[rowVals[i], colVals[i]] + 1
	}
	return(mat)
}



# -----------------------------------------
# Given a probObj and some data that has the true CoD, create a transition matrix
# that shows the distribution of CoD given an estimated CoD.  For example, a set
# of evidence may turn up 'Diabetes' as substantially the most likely CoD, but
# after looking at a couple hundred such estimations we find that 10% of such
# calls of 'Diabetes' are actually 'Lung.Cancer' and 4% are 'Falls'.  The result
# found here will account for that.  The transition matrix is typically added back 
# into the probObj for use in future calculations.
createDoubleDipMatrix <- function(probObj, evidence, nGrams=1, ...) {
	Q2 <- matrix(rep(0, nCols^2),   nrow=nCols, ncol=nCols, dimnames=list(finalVec, finalVec))
}



# -----------------------------------------
# A mini-call to for use in iteratively calling/guessing the CoD
# in a dataframe of MCD evidence.
bestCall.2 <- function(x, probObj, mcdCols, cofactors, nGrams=1, ...) {
	evidence = cleanICDs(unlist(x[ mcdCols ]))
	evidence = nGramFunction(x=evidence, n=nGrams, ...)
	if (! missing(cofactors)) {
		for (aCofactor in cofactors) {
			evidence <- c(evidence, cleanICDs(paste(aCofactor, unlist(x[aCofactor]), sep='.')))
		}
	}
	result  = sort(bayesCompute(evidence=evidence, probObj=probObj, ...)$naiveBayesVec, decreasing=TRUE)
	myCall  = names(result)[1]
	prob    = result[1]
	logOdds = log(prob/(sum(result[-1])))
	return(c(bestCall=myCall, logOdds=logOdds))
}

# OK, like above, but meant for use by the following call.   It will produce a matrix of allCauses x allSubjects.
#  apply(trainDat, 1, function(x) { bestCall.3(x=x, probObj=probObj, mcdCols=mcdCols, cofactors=cofactors, nGrams=1))
bestCall.3 <- function(x, probObj, mcdCols, cofactors, nGrams=1, allCauses, ...) {
	if (missing(allCauses)) allCauses = colnames(probObj$QnBayes)
	allCauses = sort(make.names(allCauses))
	evidence  = cleanICDs(x[ mcdCols ])
	evidence  = nGramFunction(x=evidence, n=nGrams, ...)
	for (aCofactor in cofactors) {
		evidence <- c(evidence, toupper(paste(aCofactor, unlist(x[aCofactor]), sep='.')))			
	}
	answer <- bayesCompute(evidence=evidence, probObj=probObj)
	return(answer$naiveBayesVec[allCauses])
}





# -----------------------------------------
# Cross validate with the various sample fractions
# x            the original data frame
# nSamples     number of cv reps you would like to run
# testFraction fraction that will be test data, the rest will be training data.
# Dirichlet sampling from gtools
cvRun <- function(x, nSamples=17, testFraction=0.25, nGrams=1, priorMethod='Good-Turing', constant=.Machine$double.eps, 
					mcdCols=c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51"),
					cofactorCols=c('gender', 'ageGroup'), useQ2=TRUE, columnNormalize=FALSE, 
					progress=TRUE, ...) {

	csmf         <- NULL    # return object
	csmf2        <- NULL    # return object
	pCCC         <- NULL    # return object
	bigTest      <- NULL    # return object
	bayesDat     <- NULL    # return object
	x$gs_text    <- make.names(x$gs_text)
	cutoff       <- round(nrow(x) * testFraction)    # size of the testDat
	allCauses    <- sort(unique(x$gs_text))
	N            <- length(allCauses)
	probList     <- NULL
	if (progress) cat('\n  Working on CV replicate: \n')
	for (aRep in 1:nSamples) {
		if (progress) iterationCount(i=aRep, reportFreq=1, last=nSamples)
		ttList              <- cvSample(x, testFraction)
		trainDat            <- ttList$trainDat
		testDat             <- ttList$testDat
		
		probObj  <- getMCDTransitionMatrix(x=trainDat, nGrams=nGrams,
						priorMethod=priorMethod, constant=constant,
						mcdCols=mcdCols, cofactors=cofactorCols,
						finalCause="gs_text", columnNormalize=columnNormalize, ...)
		
		if (progress) plot(table(testDat$gs_text), main=paste('Iteration', aRep))
		
		# Do the calculations for the training data
		trainResult <- apply(trainDat, 1, function(x) { 
						bestCall.3(x=x, probObj=probObj, mcdCols=mcdCols, cofactors=cofactorCols, 
									nGrams=nGrams, allCauses=allCauses, ...) })
		trainResult          <- data.frame(t(trainResult), stringsAsFactors=FALSE)
		trainResult$observed <- trainDat$gs_text
		if (!all(sort(unique(trainResult$observed)) == allCauses))
			stop(paste('\n Failure to have allCauses in the trainResult on rep', aRep, '\n'))
		
		# Do the calculations for the test data
		testResult <- apply(testDat, 1, function(x) { 
						bestCall.3(x=x, probObj=probObj, mcdCols=mcdCols, cofactors=cofactorCols, 
									nGrams=nGrams, allCauses=allCauses, ...) })									
		testResult <- t(testResult)

		# Correct for expected missassignment by reversing choice bias from Bayes estimate.
		# For example, deaths labelled 'Diabetes' are often from many differnt causes.
		if (useQ2) {
			Q2             <- ddply(trainResult, .(observed), colwise(sum))
			Q2             <- t(apply(Q2[,-1], 1, function(x) x/sum(x)))
			rownames(Q2)   <- colnames(Q2)
			probObj$Q2     <- t(Q2)     # Whoa, I think I have to do this so that columns are the TRUTH (gs_text)
			testResult     <- exp(testResult %*% log(probObj$Q2))  # 1/21/15 use the log(Q2) ???  
		}
		probList[[aRep]]  <- probObj
		
		testResult <- data.frame(testResult, stringsAsFactors=FALSE)
		prefix     <- data.frame(cvReplicate=aRep, i=1:nrow(testResult), sid=testDat$sid,
								inferred=NA, gs_text=testDat$gs_text, matchPosition=NA,
								stringsAsFactors=FALSE)
		testResult <- cbind(prefix, testResult)
		
		# Pick the most likely CoD for each patient.   Then calculate CSMFAccuracy in the standard way.
		# In a minute I'll go back and calculate CSMF & etc on the basis of partial calls (see csmf2).
		probCols            <- which(colnames(testResult) %in% allCauses)
		testResult$inferred <- apply(testResult, 1, function(z) {  
									names(z[probCols])[order(as.numeric(z[probCols]), decreasing=TRUE)][1] })
		csmfTemp    <- calculateCSMF(inferred=testResult$inferred, 
									observed=testResult$gs_text, 
									allCauses=allCauses)
									
		# Similar to above, but calculations based on partial CoD calls.
		csmf2Temp   <- csmfCount(testResult, allCauses=allCauses)
		
		csmfTemp$CSMFAccuracy <- calculateCSMFAccuracy(inferred=testResult$inferred, 
														observed=testResult$gs_text, 
														allCauses=allCauses)

		pCCCTemp              <- calculatePCCCk(x=testResult, K=1:4, allCauses=allCauses)
		pCCCTemp$cvReplicate  <- aRep
		csmfTemp$cvReplicate  <- aRep
		csmf2Temp$cvReplicate <- aRep
		csmf        		  <- rbind(csmf, csmfTemp)
		csmf2                 <- rbind(csmf2, csmf2Temp)
		pCCC                  <- rbind(pCCC, pCCCTemp)
		bayesDat              <- rbind(bayesDat, testResult)
	}
	
	# Then go back and recalculate CCC for the whole shebang on the basis of partial calls, 
	# ie not making a specific call for each patient.
	rownames(bayesDat) <- NULL
	retList <- list(csmf=csmf, csmf2=csmf2, pCCC=pCCC, bayesDat=bayesDat, probList=probList)
	attr(retList, 'description')   <- paste('Run time', Sys.time(), '  priorMethod=', priorMethod, 
											'  constant=', constant, '  useQ2=', useQ2,
											'  nGrams=', nGrams)  
	return(retList)
}




# -----------------------------------------
dirichletSample <- function(x, minValue=1) {
	returnVec  <- NULL
	n          <- length(x)
	uniTerms   <- unique(x)
	nTerms     <- length(unique(x))
	rowSamples <- round(c(rdirichlet(1, rep(1, nTerms))*n))  # won't always exactly add up to 384 rows
	rowSamples[rowSamples < minValue] <- minValue
#	print(rowSamples)
	for (i in 1:nTerms) {
		aTerm     <- uniTerms[i]
		returnVec <- c(returnVec, sample(which(x==aTerm), rowSamples[i], replace=TRUE))
	}
	return(returnVec)
}

#dirichletSample(datBayes$gs_text[1:30])
#dirichletSample(LETTERS[1:26])
#dirichletSample(LETTERS[rep(1:26, 3)])


# -----------------------------------------
# Sample the raw data frame into two sets.   One, the training
# data set is just a random sample (except that it will have at least
# one case of each term).   The other, the test data set, will have
# a Dirichlet sample from each type.
cvSample <- function(x, testFraction) {
	cutoff      <- round(nrow(x) * testFraction)    # size of the testDat
	allCauses   <- unique(x$gs_text)                # * gs_text currently hardcoded here!
	
	rowSample   <- sample(nrow(x), replace=FALSE)   # a radomization of the row numbers
	testDat     <- x[rowSample[1:cutoff], ]         # testDat gets the first 'testFraction' part...
	trainDat    <- x[-c(rowSample[1:cutoff]), ]     # ... trainDat gets the rest.
	trainCauses <- unique(trainDat$gs_text)         # some causes might be missing from train, so

	# Steal them from test and put them in train
	lostCauses  <- allCauses[!allCauses %in% trainCauses]  
	for (aCause in lostCauses) {
		inx      <- which(aCause == testDat$gs_text)[1]
		trainDat <- rbind(testDat[inx,], trainDat)
		testDat  <- testDat[-inx,]
	}

	# Dirichlet sample the testDat
	dSample  <- dirichletSample(testDat$gs_text)
	testDat  <- testDat[dSample,]

	cvList <- list(trainDat=trainDat, testDat=testDat)
	return(cvList)
}

#x = data.frame(gs_text=c(rep(LETTERS[1:5], each=5), 'Z'), score=1:26, stringsAsFactors=FALSE)
#a = cvSample(x, testFraction=0.25); a


# -----------------------------------------
# Cause Specific Mortality Fractions (CSMF)
#  See "Robust metrics..." p. 6
# See correction at:   http://www.pophealthmetrics.com/content/12/1/7
# inferred   vector of inferred causes
# observed   vector of gold standard causes
calculateCSMF <- function(inferred, observed, allCauses) {
	if (missing(allCauses))  allCauses <- sort(unique(c(inferred,observed)))
	csmf <- data.frame(cause=allCauses, obs=0, obsFrac=0,
						pred=0, predFrac=0, stringsAsFactors=FALSE)
	if (length(inferred) != length(observed)) stop('\n  inferred and observed need to be the same length!\n\n')
	N   <- length(allCauses)
	for (aCause in allCauses) {
		inx                   <- csmf$cause == aCause
		truePos               <- sum((inferred == observed) & (observed == aCause))
		falseNeg              <- sum((inferred != observed) & (observed == aCause)) 
		csmf[inx, 'obs']      <- sum(observed == aCause) 
		csmf[inx, 'obsFrac']  <- csmf[inx, 'obs'] / length(inferred)
		csmf[inx, 'pred']     <- sum(inferred == aCause)
		csmf[inx, 'predFrac'] <- csmf[inx, 'pred'] / length(inferred)
		ccc                   <- ( (truePos / (truePos+falseNeg)) - (1/N) )  / (1 - (1/N))
		# 2/6/15 Abie:  in the case of 0 predicted and 0 observed, the CCC_j should be 1.0
		ccc                   <- ifelse(is.finite(ccc), ccc, 1)
#		cat(aCause, truePos, falseNeg, ccc, '\n')
		csmf[inx, 'CCC']      <- ccc
	}
	
	return(csmf)
}


# -----------------------------------------
# First version deals with two text vectors:
# inferred   vector of inferred causes
# observed   vector of gold standard causes
#
# See "Robust metrics..." p. 7  http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3160921/
#
calculateCSMFAccuracy <- function(inferred, observed, allCauses) {
	if (missing(allCauses))  allCauses <- c(inferred, observed)
	allCauses   <- unique(allCauses)
	n           <- length(observed)
	numerator   <- 0
	if (length(inferred) != length(observed)) stop('\n  inferred and observed need to be the same length!\n\n')
	minCSMFtrue <- 1
	for (aCause in allCauses) {
		csmfTrue    <- sum(observed == aCause) / n
		csmfPred    <- sum(inferred == aCause) / n
		minCSMFtrue <- ifelse((minCSMFtrue > csmfTrue), csmfTrue, minCSMFtrue)
		numerator   <- numerator + abs(csmfTrue - csmfPred) 
	}
	csmfMaxError <- 2*(1-minCSMFtrue)
	csmfAcc      <- 1 - (numerator / csmfMaxError)
	return(csmfAcc)
}


# -----------------------------------------
#  The following function creates a data structure and calculates results
#  about partial CSMF estimates.  Ie statistics around assertions of 83% AIDS,
#  7% Pneumonia, 6% Falls,  etc.
#  See correction at:   http://www.pophealthmetrics.com/content/12/1/7
#
#   x     bayesDat obj.  Ie a column for each cause and couple of columns
#                        to indicated cvReplicate and gs_text, etc.
csmfCount <- function(x, allCauses) {
	csmf      = NULL
	if (missing(allCauses))  allCauses <- c(x$gs_text, x$inferred)
	N         = length(allCauses)
	n         = nrow(x)
	probCols  = which(colnames(x) %in% allCauses)

#	fn <- function(z, probCols) { 
#		z[probCols] = as.numeric(z[probCols])
#		z[probCols] = z[probCols] / sum(z[probCols], na.rm=TRUE)
#		return(z)
#	}
#	
#	x[,probCols] <- apply(as.matrix(x[,probCols]), 1, 
#							function(z) { z = z/sum(z, na.rm=TRUE); return(z) })
	
	for (aCause in allCauses) {
		obs      = sum(x$gs_text == aCause)
		obsFrac  = obs / n    
		pred     = sum(x[,aCause])
		predFrac = pred / n
#		cat(aCause, obs, obsFrac, pred, predFrac, '\n')
		truePos  = pred   
		falseNeg = obs - pred 
		CCC      = ( (truePos / (truePos+falseNeg)) - (1/N) )  / (1 - (1/N))
		# 2/6/15 Abie:  in the case of 0 predicted and 0 observed, the CCC_j should be 1.0
		CCC      = ifelse(is.finite(CCC), CCC, 1)  
		score    = pred / obs    # not needed??
		csmf     = rbind(csmf, data.frame(cause=aCause, obs=obs, obsFrac=obsFrac, 
									pred=pred, predFrac=predFrac,
									CCC=CCC, score=score, stringsAsFactors=FALSE))							
	}

	return(csmf)
}

#csmfCount(subset(retList$bayesDat, cvReplicate==3))


# -----------------------------------------
# Calculate Partial Cause Concordance (PCCCk).  That is the fraction of times
# that the correct CoD occurs in the 1st, 2nd, or 3rd (etc.) position. 
# See 'Robust metrics...' p. 6-7. 
# x     bayesDat object
#
 calculatePCCCk <- function(x, K=1:4, allCauses) {
	pccDat    = NULL   # Return object
	if (missing(allCauses))  allCauses <- c(x$gs_text, x$inferred)
	N         = length(allCauses)
	n         = nrow(x)
	probCols  = which(colnames(x) %in% allCauses)
	x$matchPosition = apply(x, 1, function(z) { match(z['gs_text'], 
				names(z[probCols])[order(as.numeric(z[probCols]), decreasing=TRUE)]) })

	for (k in K) {
		C    = sum(x$matchPosition <= k) / n
		PCk  = k / N
		pccDat = c(pccDat, (C - PCk) / (1 - PCk))
	}
	pccDat            <- data.frame(t(pccDat))
	colnames(pccDat)  <- paste('pCCC', K, sep='.')
	return(pccDat)
}
 
#calculatePCCCk(x=retList$bayesDat, K=1:4)
 
 
 
# -----------------------------------------
#  Idea:  could I linearly adjust the CSMFs to fall on the 45 deg line?
#  Answ:  Yes I can, but then I'll end up with many negative CSMFs....
#         Also, will have to do separate CV splits to get linear reg coeff.  
#    * The code below is not really correct!
#
linearAdjustCSMF <- function(csmfNew, csmfOld) {
	for (aCause in unique(csmfOld$cause)) {
#		print(aCause)
		inx1 <- csmfOld$cause == aCause
		fit  <- lm(predFrac ~ obsFrac, data=csmfOld[inx1,])
		if (is.na(coefficients(fit)[2])) next
		if (coefficients(fit)[2] == 0) next
		inx2 <- csmfNew$cause == aCause
		if (sum(inx2) == 0) next
#		print(str(fit))
#		csmfNew[inx2, 'predFrac'] <- predict(fit, newdata=csmfNew[inx2,])
		csmfNew[inx2, 'predFrac'] <- csmfNew[inx2, 'predFrac'] - coefficients(fit)[1]
		csmfNew[inx2, 'predFrac'] <- csmfNew[inx2, 'predFrac'] / coefficients(fit)[2]
	}
	csmfNew <- ddply(csmfNew, .(cvReplicate), mutate, predFrac = predFrac / sum(predFrac))
	return(csmfNew)
}

#junk <- linearAdjustCSMF(retList$csmf, retList$csmf)



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



