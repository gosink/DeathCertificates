##############################################
#
# Consider the MCD fields to be likely nodes in a network.  People often die
# from  disease A --> B --> C(ardiac arrest).   So there should be an emission pattern
# of  A -> B   and  B -> C  and a common grouping of A,B, & C.
# 
# Also the emission likelihood is influenced by age, gender, & etc.
# Also you can have hidden states: ie we may not have B on the certificate, but it can be inferred.
# Might also be value in stemming the codes or otherwise considering adjacent codes
# Might also be value in parsing the long text of codes to find associations
#
# 9/16 Email from Abie:  Ah, I think I see.  My goal is to predict the gs_text column, 
#  and for the non _gc columns causadef should be very useful.  For the _gc columns, 
#  predicting the gs_text value is the hard part.
#
#   -John Gosink
#    8/12/14
#


# Details of the ICD codes:
# http://en.wikipedia.org/wiki/ICD-10
# http://en.wikipedia.org/wiki/ICD-10_Chapter_XX:_External_causes_of_morbidity_and_mortality
# CMS32_DESC_LONG_DX.txt  ICD9 codes
# icd10cm_order_2015.txt  ICD10 codes


rm(list=ls())
setwd('C:/Users/John/Desktop/DeathCertificates/')
library(lattice)
library(latticeExtra)
library(plyr)
library(gplots)
library(ggplot2)
library(gtools)       # smartbind  (to bind dfs with different columns)
library(edgeR)        # Good-Turing estimator background for transition matrix
library(MCMCpack)     # Dirichlet sampling

source('GC Methods.R')

# -----------------------------------------------
# Start by reading in the Mexico death certificate data.   This supposedly has the data from 1589 Mexican
# subjects.   These are 'gold standard' reviewed documents and 224 of them are noted as being garbage codes.
# Prep the data in the other script file, 'DeathCertificates_Basic_Data_Shaping_and_Storage.R'
load(file='dat.RData')          # Load the (1587x67) data as processed in earlier
load(file='datBayes.RData')     # In development.   Most recent predictions
load(file='icd10Vec.RData')     # 91,737 ICD-10 definitions
load(file='colDesc.RData')      # Descriptions of the column headers (547 headers as a named vector)





		
# ==============   TESTING ON SMALLER SUBSETS OF DATA   ==============
source('GC Methods.R')
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))
someIds <- sample(dat$sid[dat$icd_name != '_gc'])
dat[dat$sid %in% someIds[10:15], goodCols]

#   First, build the transition matrix
x=dat[dat$sid %in% someIds,]
probObj  <- getMCDTransitionMatrix(x=x, priorFunction='Good-Turing',
				mcdCols=c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51"),
				finalCause="gs_text", cofactors=c('gender', 'ageGroup', 'smokingGroup'))

junk = probObj$QnBayes[sample(nrow(probObj$QnBayes),30), sample(ncol(probObj$QnBayes),7)]
str(probObj);  signif(junk,2);
probObj$rawMat[probObj$codVec > 30, 1:5]



visualizeTransitionMatrix(probObj$QnBayes)

#    Second, use the transition matrix to compute the probabe outcome of test data
source('GC Methods.R')
results  <- bayesCompute(evidence=c('J189', 'J969', 'AGEGROUP.ADULT', 'gender.Female'), probObj);  str(results)
head(sort(results$naiveBayesVec, decreasing=TRUE))
showBayesCompute(results)



# ==============   AUTO-RUN THE PROCEDURE ON THE WHOLE DATA SET   ==============
source('GC Methods.R')
inx      <- dat$icd_name != '_gc'
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))
													
mcdCols=c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51")
													
probObj  <- getMCDTransitionMatrix(x=dat[inx, ], priorFunction='Good-Turing',
				mcdCols=mcdCols, finalCause="gs_text", cofactors=c('gender', 'ageGroup', 'smokingGroup'))


# A minifunction for ddply that takes in 1 line of dat
bestCall.2 <- function(x, mcdCols, cofactors, ...) {
	evidence = cleanICDs(unlist(x[ mcdCols ]))
	if (! missing(cofactors)) {
		for (aCofactor in cofactors) {
			evidence <- c(evidence, cleanICDs(paste(aCofactor, as.character(x[,aCofactor]), sep='.')))
		}
	}
	result  = sort(bayesCompute(evidence=evidence, ...)$naiveBayesVec, decreasing=TRUE)
	myCall  = names(result)[1]
	prob    = result[1]
#	print(str(result))
	logOdds = log(prob/(sum(result[-1])))
	return(c(bestCall=myCall, logOdds=logOdds))
}
bestCall.2(x=dat[16,], mcdCols=mcdCols, cofactors=c('ageGroup', 'gender'), probObj=probObj)

#someRows <- sample(nrow(dat),50)         # sample just 50 rows
#someRows <- which(dat$icd_name=='_gc')   # try just the _gc data

mcdCols=c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51")
datBayes <- dat[,goodCols]   # process the whole shebang
datBayes <- cbind(datBayes[,1:3], 'naiveBayes'=NA, 'logOdds'=NA, datBayes[,4:ncol(datBayes)]) 
datBayes$naiveBayes <- NA
datBayes$logOdds    <- NA
for (i in 1:nrow(datBayes)) {
	datBayes[i,c('naiveBayes', 'logOdds')] <- bestCall.2(x=datBayes[i,], 
											mcdCols=mcdCols, 
											cofactors=c('gender', 'ageGroup'), 
											probObj=probObj)
}
datBayes$gs_text <- make.names(datBayes$gs_text)
datBayes$logOdds <- as.numeric(datBayes$logOdds)
head(datBayes)

glort = as.numeric(datBayes$gs_text == datBayes$naiveBayes)
xyplot(jitter(glort, amount=0.02) ~ logOdds | gs_text, data=datBayes,
	type=c('p', 'smooth'), span=1, 
	as.table=TRUE, 	par.strip.text=list(lines=0.8, cex=0.7),
	scales=list(x=list(relation='free')))



# Visualize Cause Specific Mortality Fraction (CSMF)
library(ggplot2)
a = data.frame(table(datBayes$naiveBayes)); a$source='Bayes'
b = data.frame(table(datBayes$gs_text));    b$source='gs_text'
ab = rbind(a,b)
ab$Var1 <- reorder(ab$Var1, ab$Freq, median)
ggplot(ab, aes(x=Var1, y=Freq/1587, fill=source)) +
	geom_bar(position='dodge', stat='identity') +
	theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
	labs(x='Cause of Death', y='CSMF', title='Cause Specific Mortality Counts (estimated and observed)')


	
	
#write.csv(datBayes, file='datBayes.csv', quote=FALSE)  # save it for perusal in Excel & etc
#save(datBayes, file='datBayes.RData')                  # save it for further R work

# BASIC SENSITIVITY AND SPECIFICITY
sum(datBayes$gs_text == datBayes$naiveBayes); nrow(datBayes)
accuCount = ddply(datBayes, .(gs_text), summarize, count=length(gs_text), trainingCount=sum(icd_name!='_gc'),
											correct=sum(gs_text==naiveBayes), sens=correct/count)
plot(accuCount$count, accuCount$sens, main='Sensitivty as a function of the number of cases')
abline(reg=lm(accuCount$sens ~ accuCount$count), col='red', lty=2)
write.csv(accuCount, file='accuCount.csv', quote=FALSE)  # save it for perusal in Excel & etc

# Graphically depict which choices were made in the transitions
temp1 = temp2 = datBayes
inx = sample(nrow(temp1))
labelCodes =  sort(unique(datBayes$gs_text, datBayes$naiveBayes))
temp1$x   = 1
temp1$y   = factor(datBayes$gs_text, labels=labelCodes, levels=labelCodes)
temp2$x   = 2
temp2$y   = factor(datBayes$naiveBayes, labels=labelCodes, levels=labelCodes)
temp      = rbind(temp1[inx,], temp2[inx,])

sugarDat           <- data.frame(id=sort(unique(temp$sid)), pch=NA, col=NA)
rownames(sugarDat) <- sugarDat$id
sugarDat$pch       <- c(70, 77)[ map(temp$sid, temp$gender=='Male')+1 ]   # M and F unicode symbols!
sugarDat$col       <- brewer.pal(n=3, name="Dark2")[ map(temp$sid, temp$icd_name=='_gc')+1 ]  

#x = 60:85;  plot(x, x, type='n'); grid(); points(x, x, pch=x)
 
pdf(file='Info_and_Results/Prediction_calls_gender_only_cofactor.pdf', height=11, width=17)
# First, panel by the gs_text listing
xyplot(jitter(as.numeric(y)) ~ jitter(x, amount=0.05) | gs_text, groups=sid, data=temp, type=c('p', 'l', 'g'), as.table=TRUE,
	scales=list(y=list(labels=as.character(temp$y), at=as.numeric(temp$y), alternating=3, cex=0.4),
	            x=list(draw=FALSE)),
	par.strip.text=list(lines=0.8, cex=0.7),
	alpha=0.3, pch=sugarDat$pch, col=sugarDat$col, 
	ylab='', xlab='', main='gs_text  --->  naiveBayes')
	
# Then, panel by the Bayes prediction
xyplot(jitter(as.numeric(y)) ~ x | naiveBayes, groups=sid, data=temp, type=c('p', 'l', 'g'), as.table=TRUE,
	scales=list(y=list(labels=as.character(temp$y), at=as.numeric(temp$y), alternating=3, cex=0.4),
	            x=list(draw=FALSE)),
	par.strip.text=list(lines=0.8, cex=0.7),
	alpha=0.3, pch=sugarDat$pch, col=sugarDat$col, 
	ylab='', xlab='', main='gs_text  --->  naiveBayes')	
dev.off()


# 10/26/14 -  Adding a prior constant of 2 to every cell *makes things worse* as opposed to 1/(Q row sum)
#          -  Adding a prior constant of 1/(2 Q rowsum) --> *a bit better & a bit worse* than 1/(Q row sum) 
#          -  Adding a prior constant of .Machine.double.eps -->  *a bit better & a bit worse* than 1/(Q row sum) 
#          - check issue with inability to find transition probablilities for various items
#          - check issue with transition probablilities for AGEGROUP.NA and other missing values
#          - set up automated sensitivity and specificity test for each category
	
	

	
# Take a closer look at pneumonia:
inx <- which(datBayes$gs_text == 'Pneumonia')      # gs_text says they're pneumonia
inx <- inx[ order(datBayes[inx,'naiveBayes']) ]    
datBayes[inx, ]
ddply(datBayes[inx,], .(naiveBayes), summarize, res=length(gs_text))   # here is what my code called them

# Following figure is NOT informative.   Columnwise, both the well predicted and poorly predicted
# data look similar to each other..
plot(sort(probObj$QnBayes[,'Pneumonia']), xlim=c(400, probObj$dimQ[1]+10)); grid()  # very poor sens
points(sort(probObj$QnBayes[,'Stillbirth']),                    col='green')  # very poor sens
points(sort(probObj$QnBayes[,'Other.Cardiovascular.Diseases']), col='red')    # very poor sens
points(sort(probObj$QnBayes[,'AIDS']), col='grey')       # very good sens
points(sort(probObj$QnBayes[,'Diabetes']), col='blue')   # very good sens
	
# See, sid= 2281  (amazing call),   515  bad call even if terms look like pnemonia

	
	
# ====================  C S M F  =================== #
# Several papers by Abie and Murray et. al. indicate that automated methods
# for determining mortality often fail.  I am now writing up functions to
# calculate the Cause Specific Mortality Fraction (CSMF) for a number of cross
# validation runs on the data.   As per the papers I create different subsets of
# data, each with bigger or smaller fractions of various causes of death (by
# resampling with replacement).   Then I split those data into training and
# test data sets and look at the relationship between expected and observed
# fractions for each mortality group.
source('GC Methods.R')
categories=c('Stroke', 'AIDS', 'Diabetes', 'Other.Non.communicable.Diseases', 'Breast.Cancer')
categories=c('Stroke', 'AIDS')
categories=make.names(c('Prostate Cancer'))
categories=make.names(unique(dat$gs_text))


csmf  <- cvRun(x=dat, fractions=seq(0, 2, length.out=17), categories=categories, testFraction=0.25)
str(csmf); csmf.bak = csmf
maxFrac <- max(c(csmf$predFrac, csmf$obsFrac), na.rm=TRUE)
xyplot(predFrac ~ obsFrac | cause, data=csmf, 
	type=c('p', 'g'), pch=19,
	subset=category==cause, #xlim=c(0,maxFrac), ylim=c(0,maxFrac),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='free', rot=45), y=list(relation='free', rot=45)),
	as.table=TRUE,  xlab='True Cause Fraction', ylab='Estimated Cause Fraction',
	panel = function(x,y, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=2)
		panel.lmline(c(x,y), c(x,y), col='red')
	})

	
subset(csmf, subset=(category=='Poisonings' & cause==category))
	
	
# TO DO:
# - Dirichlet sampling for cross validation
#  Next two from "Robust metrics for verbal autopsy"
# - Correcting for concordance p.5
# - CSMFAccuracy   p. 7
#
#
#



