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
	codVec  <- rowSums(QnBayes)
	rawMat  <- QnBayes     # will I need this?

	
	# Need to add a background or prior value on the transitions that are zero.
	# Also, transition matrix rows must sum to one.
	nRows <- nrow(QnBayes)
	for (i in 1:nRows) {
		if (sum(QnBayes[i,])==0) QnBayes[i,] <- 1/nCols
		QnBayes[i,] <- priorFunction(x=as.numeric(QnBayes[i,]), method=priorMethod, constant=.Machine$double.eps)		
#		QnBayes[i,] <- priorFunction(x=as.numeric(QnBayes[i,]), method=priorMethod, constant=1/(sum(QnBayes[i,])))
		QnBayes[i,] <- QnBayes[i,] / sum(QnBayes[i,])
	}

	probObj  <- list(Q=Q, QnBayes=QnBayes, rawMat=rawMat, icdVec=icdVec, codVec=codVec, dimQ=dim(Q), numInputExamples=numInput,
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
	if (length(evidence) == 0) 
		return(list(evidence=evidence, dropEv=dropEv, naiveBayesVec=colSums(probObj$QnBayes)/sum(probObj$QnBayes)))
	
	
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



# -----------------------------------------
# A mini-call to for use in iteratively calling/guessing the CoD
# in a dataframe of MCD evidence.
bestCall.2 <- function(x, probObj, mcdCols, cofactors, ...) {
	evidence = cleanICDs(unlist(x[ mcdCols ]))
	if (! missing(cofactors)) {
		for (aCofactor in cofactors) {
			evidence <- c(evidence, cleanICDs(paste(aCofactor, as.character(x[,aCofactor]), sep='.')))
		}
	}
	result  = sort(bayesCompute(evidence=evidence, probObj=probObj, ...)$naiveBayesVec, decreasing=TRUE)
	myCall  = names(result)[1]
	prob    = result[1]
#	print(str(result))
	logOdds = log(prob/(sum(result[-1])))
	return(c(bestCall=myCall, logOdds=logOdds))
}



# -----------------------------------------
# Cause Specific Mortality Fractions (CSMF)
#  See "Robust metrics..." p. 6
# x   vector of inferred causes
# y   vector of gold standard causes
calculateCSMF <- function(x, y, allCauses) {
	if (missing(allCauses))  allCauses <- unique(c(x,y))
	csmf <- data.frame(cause=allCauses, obs=NA, obsFrac=NA,
						pred=NA, predFrac=NA, stringsAsFactors=FALSE)
	if (length(x) != length(y)) stop('\n  x and y need to be the same length!\n\n')
	N   <- length(allCauses)
#print(data.frame(x=x, y=y))
	for (aCause in allCauses) {
		inx                   <- csmf$cause == aCause
		truePos               <- sum(y == aCause)
		trueNeg               <- sum(y != aCause)
		csmf[inx, 'obs']      <- truePos
		csmf[inx, 'obsFrac']  <- csmf[inx, 'obs'] / length(x)
		csmf[inx, 'pred']     <- sum(x == aCause)
		csmf[inx, 'predFrac'] <- csmf[inx, 'pred'] / length(x)
		ccc   <- ( (truePos / (truePos+trueNeg)) - (1/N) )  / (1 - (1/N))
		csmf[inx, 'CCC']      <- ccc
	}
	
	return(csmf)
}


# -----------------------------------------
# See "Robust metrics..." p. 7
#   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3160921/
# x   vector of inferred causes
# y   vector of gold standard causes
calculateCSMFAccuracy <- function(x, y, allCauses) {
	if (missing(allCauses))  allCauses <- unique(c(x,y))
	denom        <- 0
	if (length(x) != length(y)) stop('\n  x and y need to be the same length!\n\n')
	minTruePos  <- length(allCauses)
	for (aCause in allCauses) {
		truePos  <- sum(y == aCause)
		truePred <- sum(x == aCause)
		minTruePos <- ifelse(minTruePos > truePos, truePos, minTruePos)
		denom      <- denom + abs(truePos - truePred)/length(x) 
	}
	minTruePos <- minTruePos / length(x) 
	csmfAcc    <- 1 - (denom / (2*(1-minTruePos)))
	
	return(csmfAcc)
}




# -----------------------------------------
# Don't know yet how to do Dirichlet sampling, instead for each category
# of interest, resample more or less items of that type and then fill in 
# the rest with a random sampling of items *not* in that category. 
# For each gs_text sample(gs_text, 1:2n)  + sample(everything - number already picked)
sampleFractions <- function(x, minNum, maxNum, categories=c('Stroke', 'AIDS', 'Diabetes')) {
	minNum <- ifelse(missing(minNum), 1, minNum)  # min number of a certain category
	maxNum <- ifelse(missing(maxNum), 3, maxNum)  # a multiplicative value
	choiceFrame <- NULL
	for (aCategory in categories) {
		xRow   <- which(x == aCategory)
		notX   <- which(x != aCategory)
		xCount <- length(xRow)
		# Could fail if a certain category starts off pretty big
		for (i in 1:25) {
			numKeepers  <- ceiling(runif(1, minNum, maxNum*xCount))
			keepers     <- sample(xRow, numKeepers, replace=TRUE)
			keepers     <- c(keepers, sample(notX, length(x) - length(keepers), replace=TRUE))
			choiceFrame <- rbind(choiceFrame, data.frame(category=aCategory, run=i, 
														numKeepers=numKeepers, keep=keepers, 
														stringsAsFactors=FALSE))
		}
	}
	return(choiceFrame)
}

#junk = sampleFractions(x=datBayes$gs_text); str(junk)

# -----------------------------------------
# Cross validate with the various sample fractions
# x            the original data frame
# numSamples   the number of fractional population sizes you would like to test
# maxMultiple  max multiple of the original obs population fraction you would like to test
# categories   which gs_text categories would you like to test
# testFraction fraction that will be test data, the rest will be training data.
# fractions    multiplicative values for the observed fractions of interesting categories.
#              That is, synthesize test data sets with these multiples of the observed data fraction.
# Dirichlet sampling at:  http://127.0.0.1:20951/library/MCMCpack/html/dirichlet.html
cvRun <- function(x, fractions=c(0,0.5,1,2), categories, testFraction=0.25, ...) {
	# foreach category and run create a new data frame via the selected rows.
	# partition it in 75% train,   25% test
	# create a probObj, then run it on the test data
	# calculateCSMF and store results 
	mcdCols      <- c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51")
	cofactorCols <- c('gender', 'ageGroup')
	csmf         <- NULL    # return object
	x$gs_text    <- make.names(x$gs_text)
	x$naiveBayes <- NA
	x$logOdds    <- NA
	for (aCategory in unique(categories)) {
#		for (i in 1:numSamples) {
		for (aFraction in fractions) {
			rowSample   <- sample(nrow(x), replace=FALSE)         # a radomization of the row numbers
			cutoff      <- round(nrow(x) * testFraction)          # size of the testDat
			obsCatFrac  <- sum(x$gs_text == aCategory) / nrow(x)  # observed CSMF for this category
			targetFrac  <- obsCatFrac * aFraction
#			targetFrac <- (obsCatFrac * maxMultiple) / (numSamples/i)
			targetNum  <- ceiling(targetFrac * cutoff)
#cat('\n  -------------------\n')
cat(aCategory, '  obsCatFrac=', obsCatFrac, '  targetNum=', targetNum, '  targetFrac=', targetFrac, '  cutoff=', cutoff, '\n')

			testDat   <- x[rowSample[1:cutoff], ]
			trainDat  <- x[-c(rowSample[1:cutoff]), ]
			catInx    <- which(testDat$gs_text == aCategory)
			if (length(catInx) == 0)  {
#cat('reworking\n')
				testDat <- rbind(x[sample(which(x$gs_text == aCategory), 1), ], testDat)
				cat('reworking\n')
				catInx <- 1
			}

#cat('\n  * Sum of catInx:', sum(catInx), '\n\n')
			x1        <- testDat[catInx,]    # all subjects in the test data set who ARE in the category
			x2        <- testDat[-catInx,]   # everything in testDat that is NOT the category of interest
			
			x1Rows   <- sample(nrow(x1), targetNum, replace=TRUE)
			testDat  <- x1[x1Rows,]                          # make the new testDat with a sample of the target category...
			if (targetNum == 0) testDat <- NULL
			x2Inx    <- min(c(nrow(x2), cutoff-targetNum))   # plus enough non-category to fill it out to expected length
			testDat  <- rbind(testDat, x2[1:x2Inx,])
			probObj  <- getMCDTransitionMatrix(x=trainDat, priorFunction='Good-Turing',
							mcdCols=mcdCols, cofactors=cofactorCols,
							finalCause="gs_text")
#print(dim(testDat))
			for (j in 1:nrow(testDat)) {
				testDat[j,c('naiveBayes', 'logOdds')] <- bestCall.2(x=testDat[j,], 
														probObj=probObj,
														mcdCols=mcdCols, 
														cofactors=cofactorCols)
			}
#print(tail(testDat[,1:10]))			
			csmfTemp            <- calculateCSMF(x=testDat$naiveBayes, y=testDat$gs_text)
			csmfTemp$targetNum  <- targetNum
			csmfTemp$targetFrac <- targetFrac
			csmfTemp$category   <- aCategory
			csmf                <- rbind(csmf, csmfTemp)
		}
	}
	return(csmf)
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







