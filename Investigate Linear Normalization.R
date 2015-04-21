##############################################
#
# Look into back calculation for my Naive Bayes machine.
# After 4-5 days of working on this and useQ2, I don't see
# them really helping.   In fact they seem to casue more problems
# than they solve.
#
#   -John Gosink
#    4/16/15
#

rm(list=ls())
setwd('C:/Users/John/Desktop/DeathCertificates/')
library(lattice)       # graphics
library(latticeExtra)  # graphics
library(plyr)          # Map/Reduce type operations
library(gplots)        # more graphics
library(ggplot2)       # 
library(gtools)        # smartbind  (to bind dfs with different columns)
library(edgeR)         # Good-Turing estimator background for transition matrix
library(MCMCpack)      # Dirichlet sampling

source('GC Methods.R')
source('C:/Users/John/Desktop/R_Scripts/gosink_methods.R')

# -----------------------------------------------
# Start by reading in the Mexico death certificate data.   This supposedly has the data from 1589 Mexican
# subjects.   These are 'gold standard' reviewed documents and 224 of them are noted as being garbage codes.
# Prep the data in the other script file, 'DeathCertificates_Basic_Data_Shaping_and_Storage.R'
load(file='dat.RData')          # Load the (1587x67) data as processed in earlier
load(file='datBayes.RData')     # In development.   Most recent predictions
load(file='icd10Vec.RData')     # 91,737 ICD-10 definitions
load(file='colDesc.RData')      # Descriptions of the column headers (547 headers as a named vector)
dat$gs_text <- make.names(dat$gs_text)


# ==============   MAKE A FRESH datBayes  AND probObj  ==============
source('GC Methods.R')
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))

someIds      <- dat[dat$ageGroup == 'Adult', 'sid']
mcdCols      <- c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51")
cofactorCols <- c('gender', 'ageGroup')

dat1         <- subset(dat[,goodCols], sid %in% someIds)
allCauses    <- sort(unique(dat1$gs_text)); 

probObj  <- getMCDTransitionMatrix(x=datBayes, priorFunction='Good-Turing', columnNormalize=TRUE,
				mcdCols=mcdCols, finalCause="gs_text", cofactors=cofactorCols)


# ==============   FIGURES AND RESULTS   ==============

source('GC Methods.R')
retList  <- cvRun(x=dat1, nSamples=15, testFraction=0.25,  nGrams=2, useQ2=TRUE, 
					mcdCols=c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51"), 
					useLinNorm=FALSE, columnNormalize=FALSE,  equalTrainSizes=TRUE,
					priorFunction='Good-Turing', constant=.Machine$double.eps)

str(retList[1:5]); names(retList)
plot(retList$csmf$predFrac, retList$csmfLinNorm$predFrac); grid(); abline(a=0, b=1)
plot(retList$csmf$CSMFAccuracy, retList$csmfLinNorm$CSMFAccuracy); grid(); abline(a=0, b=1)
xyplot(pred ~ obs | cause, data=retList$csmf, type=c('p', 'r', 'g'))


csmfNew = retList$csmf; csmfOld=retList$csmfTrain
junk = linearAdjustCSMF(csmfNew, csmfOld)

xyplot(pred ~ obs | cause, data=retList$csmf, 
	type=c('p', 'g'), pch=19, alpha=0.3, 
	main=as.character(retList$probList[[1]]$createTime),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='same', rot=45), y=list(relation='same', rot=45)),
	as.table=TRUE,  xlab='True Cause Fraction', ylab='Estimated Cause Fraction',
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=1)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	})


# Figuring out the linear regression	
x = subset(retList$csmf, cvReplicate==8)
fit = lm(x$pred ~ x$obs)	
plot(x$obs, x$pred); grid(); abline(reg=fit); abline(a=0,b=1,col='red', lty=2, lwd=2)
newPred <- (x$pred - coefficients(fit)[1]) / coefficients(fit)[2]
points(x$obs, newPred, col='blue', pch=19)	
abline(reg=lm(newPred ~ x$obs), col='blue')


junk <- with(retList, cbind(csmf[,c('cvReplicate', 'cause', 'obs', 'pred')], csmfTrain[,c('obs', 'pred')]))
keepCols <- c('cvReplicate', 'cause', 'obs', 'pred')
junk <- with(retList, rbind(cbind(csmf[,keepCols],group='csmf'), cbind(csmfTrain[,keepCols], group='csmfTrain')))
xyplot(pred ~ obs | cause, groups=group, data=junk,
	type=c('p', 'g', 'r'), pch=19, alpha=0.8,  
	main=as.character(retList$probList[[1]]$createTime),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='free', rot=45), y=list(relation='free', rot=45)),
	as.table=TRUE,  xlab='True Obs Count', ylab='Estimated Obs Count',
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=1)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	})

# Get the data from a 500x CV run on the cloud
load('Data/retListPart.RData')
temp = rbind(retListPart[[1]]$csmfTrain, retListPart[[1]]$csmf, retListPart[[1]]$csmfLinNorm)
temp$group = rep(c('train', 'csmf', 'linNorm'), each=nrow(retListPart[[1]]$csmf))

# Do linAdjust by hand
source('GC Methods.R')
junk = linearAdjustCSMF(csmfNew=retListPart[[1]]$csmf, csmfOld=retListPart[[1]]$csmfTrain)
#inx=sample(13500,15); junk[inx,], junk1[inx,]
temp = rbind(junk, retListPart[[1]]$csmf)
temp$group = rep(c('handLinNorm', 'csmf'), each=nrow(retListPart[[1]]$csmf))

xyplot(pred ~ obs | cause, groups=group, data=temp,
#	subset=cause %in% c('AIDS', 'Breast.Cancer', 'COPD', 'Diabetes'),
	type=c('p', 'g', 'r'), pch=1, alpha=0.5,  
	main=as.character(retList$probList[[1]]$createTime),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='same', rot=45), y=list(relation='free', rot=45)),
	as.table=TRUE,  xlab='True Obs Count', ylab='Estimated Obs Count',
	auto.key=list(space='bottom', columns=2, lines=TRUE),
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	})



