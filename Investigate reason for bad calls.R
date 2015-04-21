##############################################
#
# Run big runs on the Amazon EC2 cluster and collect the data.  Material here
# is just for the (white)paper.
#
#   -John Gosink
#    3/24/15
#

rm(list=ls())
setwd('C:/Users/John/Desktop/DeathCertificates/')
library(lattice)       # graphics
library(latticeExtra)  # graphics
library(plyr)          # Map/Reduce type operations
library(gplots)        # more graphics
library(ggplot2)       # 

source('GC Methods.R')
source('C:/Users/John/Desktop/R_Scripts/gosink_methods.R')

# -----------------------------------------------
load(file='dat.RData')          # Load the (1587x67) data as processed in earlier
load(file='datBayes.RData')     # In development.   Most recent predictions
load(file='icd10Vec.RData')     # 91,737 ICD-10 definitions
load(file='colDesc.RData')      # Descriptions of the column headers (547 headers as a named vector)


load(file='Data/retList.RData')
load(file='Data/pCCCList.RData')
load(file='Data/testingGrid.RData')
load(file='Data/csmfList.RData')
load(file='Data/csmfList2.RData')

	
	
# ====================  C S M F  =================== #

#### SHOW THE STRUCTURE OF THE ORIGINAL DATA ####
source('GC Methods.R')
visualizeCausalChains(x=subset(dat, subset=(ageGroup=='Adult' & gs_text=='Falls')), maxLabels=41)
visualizeCausalChains(x=subset(dat, subset=(ageGroup=='Adult' & gs_text=='Pneumonia')))


goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))
dat1 <- subset(dat, gs_text=='Falls')[,goodCols]

#### SHOW RECOVERY OF _gc CODES AND CORRECT ASSIGNMENT ####
temp   = retList$bayesDat
gcTemp = NULL
temp   = merge(dat[,c('sid', 'icd_name')], temp, by='sid', sort=FALSE)
allCauses = colnames(temp)[8:34]

temp1 = subset(temp, gs_text=='Falls')

# Here are the people grouped by their difficulty to predict
table(temp1$sid, temp1$matchPosition)
a = table(temp$gs_text, temp$sid, temp$matchPosition)
a = dlply(temp, .(gs_text), summarise, table(sid, matchPosition))

###### ANSWER 1
#  - The CSMF Accuracy is not (necessarily) related to the size of the training data set.  
#    For example, Other Non-communicable Diseases and Pneumonia are amongst the most highly
#    represented CoD's but have very low CCCs.  Conversely Suicide and Prostate Cancer have very
#    small numbers but are well predicted.  This is broadly consistent with the 'Assessing quality' 
#    paper.

#####  ABOUT THE UNDERLYING DATA...
# ==============   MAKE A FRESH datBayes  AND probObj  ==============
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))
mcdCols  <- c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51")
 
probObj  <- getMCDTransitionMatrix(x=subset(dat, ageGroup=='Adult'), priorFunction='Good-Turing', 
				cofactors='gender', columnNormalize=FALSE, mcdCols=mcdCols, finalCause="gs_text")

visualizeTransitionMatrix(probObj$QnBayes)


# 1.	Number of training examples
x = sort(table(dat$gs_text), decreasing=TRUE)
xyplot(x ~ 1:length(x), type=c('h','g'),
	lwd=4, main='Number of examples of each class',
	ylab='Number of examples', xlab=NULL,
	scales=list(x=list(labels=names(x), rot=45, at=1:length(x))))
	
# 2.	Heterogeneity of the training samples
# calculate median percent novelty per new sequence
allCauses = sort(unique(subset(dat, ageGroup=='Adult')$gs_text))
wordList  = NULL
for (aCause in allCauses) {
	x     = dat[dat$gs_text==aCause, mcdCols]
	xList = apply(x, 1, function(z) nGramFunction(unlist(cleanICDs(z)), n=1))
	wordList[[aCause]] = as.character(unlist(xList))
}
	
	
temp = NULL
for (aCause in allCauses) {
	x        = dat[dat$gs_text==aCause, mcdCols]
	xList    = apply(x, 1, function(z) nGramFunction(unlist(cleanICDs(z)), n=1))
	wordSoup = unlist(wordList[-which(names(wordList)==aCause)])
	# I think the following loop is wrong as we will be considering the first couple
	# of patients many more times than the last couple of patients.  Fixed now?
	for (i in 2:length(xList)) {
		fracTotal = i / length(xList)
		soupSample   = sample(wordSoup, 0.75 * length(wordSoup))
		for (j in 1:3) {
			rowInx    = sample(length(xList))[1:i]
			icdsInRow = unlist(xList[rowInx[1]])
			icdsSeen  = unlist(xList[rowInx[-1]])
			fracSeen  = sum(icdsInRow %in% icdsSeen) / length(icdsInRow)
			# The following is probably wrong for these type of estimates.
			# table(match(LETTERS[c(1:20, 3:7)], LETTERS[1:5]))
			outsideHits         = sort(table(match(soupSample, icdsInRow)))
			leastOutsideHits    = outsideHits[1]
			medianOutsideCause  = median(outsideHits)
			temp = rbind(temp, data.frame(aCause=aCause, i=i, j=j, fracSeen=fracSeen, 
						fracTotal=fracTotal, leastOutsideHits=leastOutsideHits, 
						medianOutsideCause=medianOutsideCause,
						stringsAsFactors=FALSE))
		}
	}		
}


junk = ddply(temp, .(aCause), summarize, median(medianOutsideCause, na.rm=TRUE))
junk[order(junk[,2]),]

bwplot(medianOutsideCause ~ reorder(factor(aCause), medianOutsideCause, median, na.rm=TRUE), data=temp,
		main='Specificity of the ICD codes for each cause',
		ylab='Median number of times an ICD code\nis found in other causes',
		scales=list(x=list(rot=45)));  bringToTop()

bwplot(leastOutsideHits ~ reorder(factor(aCause), leastOutsideHits, median, na.rm=TRUE), data=temp,
		main='Specificity of the ICD codes for each cause',
		ylab='Minimum number of times an ICD code\nis found in other causes',
		scales=list(x=list(rot=45)));  bringToTop()

# NB, THERE IS NO LINEAR CORRELATION BETWEEN MEDIAN SPECIFICITY AND MEDIAN CCC....
#  HMM, NEED TO ABANDON THIS AS A PREDICTOR TOO...

		
		
xyplot(jitter(fracSeen) ~ i | aCause, data=temp, type=c('a','smooth','g'),
	as.table=TRUE,
	main='Fraction of ICD codes seen with increasing training size',
	xlab='Number of distinct cases in the training set',
	ylab='Fraction of ICD codes found in the\ntraining set for a new case')

xyplot(jitter(fracSeen) ~ i | aCause, data=temp, type=c('a','g'),
	as.table=TRUE, lwd=2,
	main='Fraction of ICD codes seen with increasing training size',
	xlab='Number of distinct cases in the training set',
	ylab='Fraction of ICD codes found already\n in the training set for a new case',
	panel=function(x,y,...) {
		panel.xyplot(x, y, ...)
		aFit = lm(y ~ x + I(x^0.5))
		newX = seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=100)
		newY = predict(aFit, newdata=data.frame(x=newX), type='response')
		panel.xyplot(newX, newY, col='red', type='l', lwd=1)
		panel.rug(x=newX[75], col='red', start=0.1)
		panel.rug(y=newY[75], col='red', end=0.1, lwd=2)
		cat(newY[75], '\n')
	})
	
bringToTop()
save(temp, file='seenICDRate.RData')
load(file='seenICDRate.RData')



# 3.	Similarity between the classes.
dmat <- dist(t(round(probObj$QnBayes, 3))) 
		
x    <- probObj$rawMat
x    <- x[!grepl('GENDER', rownames(x)),]
x    <- apply(x, 1, function(z) z/sum(z))
dmat <- dist(x)

hc   <- hclust(dmat, method='ward')
plot(hc); bringToTop()

visualizeCausalChains(x=subset(dat, subset=(ageGroup=='Adult' & gs_text=='Cirrhosis')))
x = probObj$rawMat[,c('Cirrhosis','Maternal', 'Prostate.Cancer', 'Diarrhea.Dysentery')]
x[rowSums(x) > 1,]


# 2D heatmap or graph
library(MASS)
dmat <- as.matrix(dist(t(probObj$QnBayes)))
heatmap.2(dmat, dendrogram='both', col=gray(32:0/32), 
		margins=c(10,14), trace='none', srtCol=45);  bringToTop();
		
fit <- cmdscale(dmat, eig=TRUE, k=2) # k is the number of dim
fit <- isoMDS(dmat, k=2)             # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Metric MDS", type="n", bty='n')
text(x, y, labels = row.names(dmat), cex=.7)


# Manually encode the medianCCC along with the data size and the
# ICDSeen75 score (see above) and the median number of ICD codes found
# outside of the specific CoD.   Model this for explained variance.
x = read.csv('Data/Underlying data structure.csv', stringsAsFactors=FALSE); str(x)
fit = lm(medianCCC ~ Size + ICDSeen75 + SpecICDs, data=x)
summary(fit)
par(mfrow=c(2,3))
plot(fit, which=c(1:6))
# The fit seems pretty good.   Other.Non.communicable.Diseases  is the biggest outlier, but its
# not immediately apparent that there is any reason to reject it.
# My CCC are broadly consistent with (Assessing quality... Hernández et al. Population Health Metrics 2011, 9:38)


library(visreg)
par(mfrow=c(1,2), cex.lab=0.9, mar=c(5, 5, 4, 2) + 0.1)
y = visreg(fit, xvar='Size', main='', points=list(cex=1, pch=1),
	xlab='Number of cases available in the original data', 
	ylab='medianCCC\n(other factors kept constant)'); grid()
identify(x$Size, y$res$visregRes, labels=x$cause, cex=0.75)

y = visreg(fit, xvar='ICDSeen75', main='', points=list(cex=1, pch=1),
	xlab='Fraction of ICD codes seen\nin a 75% training data', 
	ylab='medianCCC\n(other factors kept constant)'); grid()
identify(x$ICDSeen75, y$res$visregRes, labels=x$cause, cex=0.75)		
	


