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


datBayes = subset(datBayes, subset=c(gs_text != 'Stillbirth'))



# ==============   MAKE A FRESH datBayes  AND probObj  ==============
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))
someIds  <- dat[dat$gs_text != 'Stillbirth', 'sid']
someIds  <- dat[dat$ageGroup == 'Child', 'sid']
mcdCols  <- c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51")

datBayes <- dat[dat$sid %in% someIds, goodCols]   
datBayes <- cbind(datBayes[,1:3], 'naiveBayes'=NA, 'logOdds'=NA, datBayes[,4:ncol(datBayes)]) 
probObj  <- getMCDTransitionMatrix(x=datBayes, priorFunction='Good-Turing', columnNormalize=TRUE,
				mcdCols=mcdCols, finalCause="gs_text", cofactors=c('gender', 'ageGroup'))

for (i in 1:nrow(datBayes)) {
	datBayes[i,c('naiveBayes', 'logOdds')] <- bestCall.2(x=datBayes[i,], 
											mcdCols=mcdCols, nGrams=2,
											cofactors=c('gender', 'ageGroup'), 
											probObj=probObj)
}
datBayes$gs_text <- make.names(datBayes$gs_text)
datBayes$logOdds <- as.numeric(datBayes$logOdds)
head(datBayes)
#write.csv(datBayes, file='datBayes.csv', quote=FALSE, row.names=FALSE)  # save it for perusal in Excel & etc
#save(datBayes, file='datBayes.RData')                                   # save it for further R work

# Diagnose failures of the system:
subject  <- datBayes[datBayes$sid==1503, ]
evidence <- unlist(c(subject[,mcdCols], as.character(interaction('AGEGROUP',subject$ageGroup)), 
					as.character(interaction('GENDER',subject$gender))))
results  <- bayesCompute(evidence=evidence, probObj)
head(sort(results$naiveBayesVec, decreasing=TRUE))
showBayesCompute(results)

source('GC Methods.R')
bestCall.3(x=datBayes[150,], mcdCols=mcdCols, cofactors=c('gender', 'ageGroup'), probObj=probObj)
probObj = retList$probList[[10]]
a=apply(subset(datBayes, gs_text=='Suicide'), 1, 
		function(x) bestCall.3(x=x, mcdCols=mcdCols, cofactors=c('gender', 'ageGroup'), probObj=probObj))
a; apply(a, 2, rank)
subset(retList$bayesDat, gs_text=='Suicide')[,1:8]
		
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




# 10/26/14 -  Adding a prior constant of 2 to every cell *makes things worse* as opposed to 1/(Q row sum)
#          -  Adding a prior constant of 1/(2 Q rowsum) --> *a bit better & a bit worse* than 1/(Q row sum) 
#          -  Adding a prior constant of .Machine.double.eps -->  *a bit better & a bit worse* than 1/(Q row sum) 
#          - check issue with inability to find transition probablilities for various items
#          - check issue with transition probablilities for AGEGROUP.NA and other missing values
#          - set up automated sensitivity and specificity test for each category
	

	
	
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
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))

exclusions  <- c('Stillbirth')  #, 'Poisonings', 'Suicide', 'Diarrhea/Dysentery')
dat1        <- subset(dat[,goodCols], subset=(!(gs_text %in% exclusions)) & (ageGroup == 'Adult'))

retList  <- cvRun(x=dat1, nSamples=5, testFraction=0.25,  nGrams=2, useQ2=TRUE, 
					columnNormalize=TRUE,
					priorFunction='Good-Turing', constant=.Machine$double.eps)
str(retList, max.level=2); retList.bak = retList



## See DeathCertificates_via_Bayes_for_paper_v_1.0.Rmd  - there I set up a loop to analyze many different
#  combinations of parameters.   I store the results in:
#  testingGrid - 1_13_15.RData    37 CVs of many different items  (adult, child, and neonate)
#  testingGrid.RData             101 CVs of gender, nGrams, prior, constant  (Adults only)
#
load('testingGrid.RData')
load('testingGrid - 1_21_15.RData')  # this is the better log(Q2) version
testingGrid[order(testingGrid$csmfAccuracy, decreasing=TRUE),]
#fit = lm(csmfAccuracy ~  gender + useQ2 + nGrams + useQ2, data=testingGrid); summary(fit)
fit = lm(csmfAccuracy ~  gender + columnNormalize + useCausadef, data=testingGrid); summary(fit)


#
# Results are:
#  - Adults are FAR BETTER predicted than neonates or children
#  - 2 vs. 1 nGrams give 5-6% improvement on CSMFAccuracy    
#  - including gender DROPS CSMFAccuracy by 3-4 percent                 (CHECK THIS!)
#  - using a Good-Turing prior is 4-5% CSMFAccuracy better than square.min
#  - constant of 1 (vs. e-16) drops CSMFAccuracy ~15%
#  - CSMFAccuracy is highly correlated to both medianCCC and medianpCCC1
#  - picking a single best assignment is a few percent better medianCCC than using fractional assignments
#  - the same patterns of 'from' and 'to' assignments occur for adult, child, and neonates.  
#    This is frustrating as I thought the multiplication by probObj$Q2 would eliminate this.  
#     useQ2  may DROP CSMFAccuracy by 3% 
#  - use of columnNormalization improves CSMFAccuracy by ~3-4% 
#  - use of causadef improves CSMFAccuracy by ~4-5%






# Two different ways of plotting the data
xyplot(pred ~ obs, groups=cause, data=retList$csmf, type=c('p', 'r', 'g'),
	par.settings=simpleTheme(col=brewer.pal(n=8, 'Dark2'), pch=19:23, lty=1:2),
	auto.key=list(space='right', lines=TRUE, cex=0.8))
	
#score = brewer.pal(n=6, 'Spectral')[cut(retList$csmf$score, breaks=c(-1,1,2,4,8,16,2000))]
xyplot(predFrac ~ obsFrac | cause, data=retList$csmf2, 
	type=c('p', 'g'), pch=19, alpha=0.3, #col=score, 
	main=as.character(retList$probList[[1]]$createTime),
	par.strip.text=list(lines=1, cex=0.7),
	scales=list(x=list(relation='same', rot=45), y=list(relation='same', rot=45)),
	as.table=TRUE,  xlab='True Cause Fraction', ylab='Estimated Cause Fraction',
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=1)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	})


# Look at various things
source('GC Methods.R')
# Work out process to calculate stats on partial deaths
bayesDat <- retList$bayesDat; #head(bayesDat)
csmf2    <- csmfCount(subset(bayesDat, cvReplicate < 5))
calculateCSMFAccuracy(inferred=datBayes$naiveBayes, observed=datBayes$gs_text, allCauses=datBayes$gs_text)


y_true = c('Stroke', 'Diabetes', 'Other')
y_pred = c('Stroke', 'Other', 'Diabetes')    # ccc = 0.0
calculateCSMF(inferred=y_pred, observed=y_true)

y_true = c('Stroke', 'Diabetes', 'Other')
y_pred = c('Diabetes', 'Other', 'Stroke')    # ccc = -0.5
calculateCSMF(inferred=y_pred, observed=y_true)

y_true = c(3, 4, 3, 4, 5, 5, 2, 5, 2, 5)
y_pred = c(2, 2, 4, 5, 3, 4, 3, 2, 4, 2)              # ccc = -0.333333333333
calculateCSMF(inferred=y_pred, observed=y_true)

y_true = c(15,  1, 16,  1,  6,  7, 18,  9,  3, 14)
y_pred = c(16,  3,  3,  7,  6, 11,  0,  2,  2, 16)    # ccc = 0.0181818181818
calculateCSMF(inferred=y_pred, observed=y_true)


# How many replicates do I need?
x = retList$csmf
plot(1,1, xlim=c(1, max(x$cvReplicate)), ylim=c(0.5,0.85), type='n'); grid()
cumMed = NULL
for (i in sort(unique(x$cvReplicate))) {
	inx = which(x$cvReplicate == i)[1]
	points(i, x[inx, 'CSMFAccuracy'])
	cumMed = c(cumMed, x[inx, 'CSMFAccuracy'])
	points(i+0.1, median(cumMed), pch=2)
}


plot(retList$csmf$CSMFAccuracy, retList$csmf$CCC)
subset(retList$csmf, subset=(category=='Poisonings' & cause==category))
temp = subset(retList$bigTest, cvReplicate==21)  # Poisonings diagnosed very poorly
temp1 = ddply(retList$bigTest, .(cvReplicate, naiveBayes), summarize, medLogOdds = c(median(logOdds)))
temp2 = ddply(temp1, .(naiveBayes), mutate, rankLoD = rank(medLogOdds))
plot(temp2$cvReplicate, temp2$rankLoD)
qqmath(~medLogOdds | naiveBayes, data=temp2)

# How do CSMF compare when using all of the info in the Naive Bayes estimate vs.
# picking just the top scoring CoD from the estimate? 
b = merge(retList$csmf, retList$csmf2, by=c('cause', 'cvReplicate'))
xyplot(predFrac.y ~ predFrac.x, groups=cause, data=b, type=c('p', 'g'),
		par.settings=simpleTheme(col=brewer.pal(n=8, 'Dark2'), pch=19:23))
b=ddply(retList.bak$csmf2, .(cause), summarize, median(CCC[is.finite(CCC)]))
a=ddply(retList.bak$csmf, .(cause), summarize, median(CCC[is.finite(CCC)]))
temp = cbind(a[,1:2],b[,2])
colnames(temp) <- c('cause', 'csmf', 'csmf2')
xyplot(csmf2 ~ csmf, groups=cause, data=temp, type=c('p','g'), 
	par.settings=simpleTheme(pch=1:32), auto.key=list(space='right', points=TRUE))

		
		
# Do the frequencies of misattrribution look the same irrespective of the training
# data set and the test data set?   If so I can implement the 'double-dipping'
# stragegy to get a better estimate of the death causes.   The answer is, yes,
# this looks like it could be a useful exercise.
a = pairsToMatrix(rowVals=datBayes$naiveBayes, colVals=datBayes$gs_text, allVals=unique(datBayes$gs_text))
temp = pairsToMatrix(rowVals=retList$bigTest$naiveBayes, colVals=retList$bigTest$gs_text, allVals=unique(datBayes$gs_text))
a = t(apply(a, 1, function(x) { x=x+1; log(x/sum(x))} ))
heatmap(a, Rowv=NA, Colv=NA)
for (i in unique(retList$bigTest$cvReplicate)) {
	x    = subset(retList$bigTest, cvReplicate==i)
	temp = pairsToMatrix(rowVals=x$naiveBayes, colVals=x$gs_text, allVals=unique(datBayes$gs_text))
	temp = t(apply(temp, 1, function(x) { x=x+1; log(x/sum(x))} ))
	heatmap(temp, Rowv=NA, Colv=NA, main=paste('cvReplicate', i))
	Sys.sleep(1)
}

# In terms of the 'double-dipping' method, can I just sum up the results from a bayesDat
# type object from the training data?   The answer is, 'yeah...', kinda looks the same.
plot(1,1, ylim=c(0,1), xlim=c(0,5), type='n'); grid(); bringToTop()
for (i in sample(nrow(retList$bayesDat), 70)) {
	points(sort(as.numeric(retList$bayesDat[i, colnames(retList$probObj$QnBayes)]), decreasing=TRUE), type='l')
}
x = retList$bayesDat[, colnames(retList$probObj$QnBayes)]
a = matrix(rep(0, ncol(x)^2), nrow=ncol(x), ncol=ncol(x), dimnames=list(colnames(x), colnames(x)))
for (i in 1:nrow(x)) {
	a[, retList$bayesDat[i,'gs_text']] = a[, retList$bayesDat[i,'gs_text']] + as.numeric(x[i,])
}
a = t(apply(a, 1, function(x) { x=x+1; log(x/sum(x))} ))
heatmap(a, Rowv=NA, Colv=NA)


# But is this double dipping (useQ2) approach correct?.. Maybe not??
subset(retList$bayesDat, inferred=='Diabetes',  c(1:8))  # The diagram suggests, YES this suggests NO.
subset(retList$bayesDat, inferred=='Pneumonia', c(1:8))  # The diagram suggests, YES this suggests NO.
sort(table(subset(retList$bayesDat, inferred=='Diabetes', 'gs_text')))
sort(table(subset(retList$bayesDat, inferred=='AIDS', 'gs_text')))
sort(table(subset(retList$bayesDat, inferred=='TB', 'gs_text')))
Q2 = retList$probList[[25]]$Q2
for (i in 1:25) { heatmap(log(retList$probList[[i]]$Q2), Rowv=NA, Colv=NA, col=gray(32:0/32), main=i); Sys.sleep(0.5); }

# What is the average Q2 matrix like?   Gah, this is hard to understand/interpret.   Basically it averages
# out such that almost everything is exactly on the diagonal.
superQ2 = retList$probList[[1]]$Q2 * 0
for (i in 1:25) superQ2 = superQ2 + retList$probList[[i]]$Q2


# Heatmap of Q2 and my understanding of Q2 tables are in alignment
# 

 
# Test out new method (bestCall.3) for assessing all probabilities from a set of data.  Takes ~4 seconds.
trainDat=dat1[sample(nrow(dat1)),]; nGrams=1;  prior='Good-Turing'; constant=.Machine$double.eps 
mcdCols=c("causadef", "codcau11", "codcau21", "codcau31", "codcau41", "codcau51"); cofactorCols=c('gender', 'ageGroup')
allCauses = sort(unique(dat1$gs_text)); probObj=retList$probObj
a=apply(trainDat, 1, function(x) { bestCall.3(x=x, probObj=probObj, mcdCols=mcdCols, cofactors=cofactors, nGrams=1, allCauses=allCauses) })
a=data.frame(t(a), stringsAsFactors=FALSE)
a$observed=make.names(trainDat$gs_text); str(a)



		
# In the final paper as per "Robust metrics..." p. 8+	
# 1. Average chance-corrected concordance of individual cause assignment (p. 8 Discussion)
ddply(retList$csmf, .(cause), summarize, mean(CCC[is.finite(CCC)]))
with(retList$csmf, mean(CCC[is.finite(CCC)]))   # overall or grand median



# 2. pCCC(k)  -  (p. 8 Discussion)
summary(retList$pCCC)


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
mat = pairsToMatrix(rowVals=retList$bayesDat$inferred, colVals=retList$bayesDat$gs_text, 
					square=TRUE, allVals=datBayes$gs_text) 
#mat = t(apply(mat, 1, function(x) { x=x+1; log(x/sum(x))} ))
mat = log(mat+1)
heatmap(mat, Rowv=NA, Colv=NA, col=gray(32:0/32), main='Predicted (rows) vs. Observed (columns) CoD')



