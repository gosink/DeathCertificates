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

		
# ==============   FIGURES AND RESULTS   ==============


# Figure showing the transitions in the underlying CoD.  Panel by age group.
source('GC Methods.R')
pdf(file='CoD Chains.PDF', width=11, height=8.5)
visualizeCausalChains(x=datBayes, orderedCauses=c(11:6), colorGroups='gender', pchGroups='smokingGroup')
dev.off()


# Figure showing the cutoff for accurate prediction
xyplot(as.numeric(gs_text == naiveBayes) ~ logOdds | gs_text, data=datBayes,
#xyplot(as.numeric(gs_text == naiveBayes) ~ logOdds | gs_text, data=retList$bigTest,
#xyplot(as.numeric(gs_text == naiveBayes) ~ logOdds | naiveBayes, data=retList$bigTest,
	xlab='logOdds (logistic regression line in red),', ylab='False Negative or True Positive',
	main='Predicted CoD correctness vs. log odds of the prediction\n(paneled by gold standard classification)',
	scales=list(y=list(at=c(0,1), labels=c('FN', 'TP'), limits=c(-0.5, 1.5)), 
				x=list(relation='free', rot=45)), 
	as.table=TRUE, 	par.strip.text=list(lines=1, cex=0.8), alpha=0.8, #pch='.',
	panel=function(x,y,...) {
		panel.xyplot(x, jitter(y, amount=0.1), ...)
		if (length(x) > 3) {
			aFit = glm(y~x, family=binomial)
			newX = seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=100)
			newY = predict(aFit, newdata=data.frame(x=newX), type='response')
			panel.xyplot(newX, newY, col='red', type='l', lty=2)
		}
	})



# Visualize Cause Specific Mortality Fraction (CSMF) as a bargraph
library(ggplot2)
a = data.frame(table(datBayes$naiveBayes)); a$source='Bayes'
b = data.frame(table(datBayes$gs_text));    b$source='gs_text'
ab = rbind(a,b)
ab$Var1 <- reorder(ab$Var1, ab$Freq, median)
ggplot(ab, aes(x=Var1, y=Freq/1587, fill=source)) +
	geom_bar(position='dodge', stat='identity') +
	theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
	labs(x='Cause of Death', y='CSMF', title='Cause Specific Mortality Counts (estimated and observed)')


	
	
# ====================  C S M F  =================== #
# Several papers by Abie and Murray et. al. looking at CSMFAccuracy, CCC, pCCC(k) etc.
# from 500 resamplings of the data.   See the development of these metrics in:
# DeathCertificates via Bayes for paper v0.5.R

source('GC Methods.R')
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))


testingGrid <-  expand.grid(ageGroup=c('Adult', 'Child', 'Neonate'), gender=0:1,  nGrams=1:2,
							priorMethod=c('square.min', 'Good-Turing'), constant=c(.Machine$double.eps, 1),
							stringsAsFactors=FALSE)

exclusions <- c('Stillbirth')  #, 'Poisonings', 'Suicide', 'Diarrhea/Dysentery')
nSamples   <- 5

for (i in 1:nrow(testingGrid)) {
	cat('\n  ----------------------  \n\n'); print(testingGrid[i,])
	dat1 <- subset(dat[,goodCols], subset=(ageGroup==testingGrid[i,'ageGroup']) &
	                                      !(gs_text %in% exclusions))
	retList <- cvRun(x=dat1, nSamples=nSamples, testFraction=0.25, 
						nGrams=testingGrid[i,'nGrams'],
						priorMethod=testingGrid[i,'priorMethod'], 
						constant=testingGrid[i,'constant'])
						
print(xyplot(predFrac ~ obsFrac | cause, data=retList$csmf2, 
	type=c('p', 'g'), pch=19, alpha=0.5, #col=score, 
	main=as.character(retList$probObj$createTime),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='same', rot=45), y=list(relation='same', rot=45)),
	as.table=TRUE,  xlab='True Cause Fraction', ylab='Estimated Cause Fraction',
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=1)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	}))
bringToTop(); Sys.sleep(2)

}

dat1     <- subset(dat[,goodCols], subset=(!(gs_text %in% exclusions)) & (ageGroup == 'Adult'))
retList  <- cvRun(x=dat1, nSamples=5, testFraction=0.25,  nGrams=2, priorMethod='square.min', constant=.Machine$double.eps)


	
#score = brewer.pal(n=6, 'Spectral')[cut(retList$csmf$score, breaks=c(-1,1,2,4,8,16,2000))]
xyplot(predFrac ~ obsFrac | cause, data=retList$csmf2, 
	type=c('p', 'g'), pch=19, alpha=0.5, #col=score, 
	main=as.character(retList$probObj$createTime),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='same', rot=45), y=list(relation='same', rot=45)),
	as.table=TRUE,  xlab='True Cause Fraction', ylab='Estimated Cause Fraction',
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=1)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	})



# How do CSMF compare when using all of the info in the Naive Bayes estimate vs.
# picking just the top scoring CoD from the estimate? 
b=ddply(retList$csmf2, .(cause), summarize, median(CCC[is.finite(CCC)]))
a=ddply(retList$csmf, .(cause), summarize, median(CCC[is.finite(CCC)]))
temp = cbind(a[,1:2],b[,2])
colnames(temp) <- c('cause', 'csmf', 'csmf2')
temp$cause = factor(temp$cause, levels=order(temp$csmf2))
xyplot(csmf2 ~ csmf, groups=reorder(cause, -csmf2, median), data=temp, type=c('p','g'), 
	par.settings=simpleTheme(pch=1:25, cex=1.3), 
	main=paste('Comparison of CCC picking just the best CoD for each subject (x-axis)',
	           'or assigning partial fractions for each CoD for each subject (y-axis)', sep='\n'),
	xlab='Single best CoD assignment', ylab='Fractional CoD assignments',
	auto.key=list(space='right', points=TRUE, cex=0.8))



		
# In the final paper as per "Robust metrics..." p. 8+	
# 1. Average chance-corrected concordance of individual cause assignment (p. 8 Discussion)
ddply(retList$csmf, .(cause), summarize, median(CCC[is.finite(CCC)]))
with(retList$csmf, median(CCC[is.finite(CCC)]))   # overall or grand median



# 2. pCCC(k)  -  (p. 8 Discussion)
summary(retList$pCCC[,2:5])


# 3. Median CSMFAccuracy (middle left column p. 8, & Discussion)
median(retList$csmf$CSMFAccuracy)   # approx 0.71 for adults


	
# 3.5. Regression coefficients (upper right column p. 8)
fitParams <- function(x,y) {
	fit       <- lm(y ~ x)
	alpha     <- coefficients(summary(lm(y ~ x)))[1,1]
	beta1     <- coefficients(summary(lm(y ~ x)))[2,1]
	r.squared <- summary(fit)$r.squared
	sigma     <- summary(fit)$sigma
	return(c(alpha=alpha, beta=beta1, sigma=sigma, r.squared=r.squared))
}

a = ddply(retList$csmf, .(cause), function(x) fitParams(x=x$obsFrac, y=x$predFrac))
str(a); cbind(CoD=a[,1], signif(a[,2:5], 3))


# 4.  Full NxN matrix for the process
mat = pairsToMatrix(rowVals=retList$bayesDat$gs_text, colVals=retList$bayesDat$inferred, 
					square=TRUE) #, allVals=dat$gs_text) 
mat = t(apply(mat, 1, function(x) { x=x+1; log(x/sum(x))} ))
heatmap(mat, Rowv=NA, Colv=NA, main='Observed (rows) vs. Predicted (columns) CoD',
        col=gray(32:0/32),  margins=c(7,7))


