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


datBayes = subset(datBayes, subset=c(gs_text != 'Stillbirth'))


		
# ==============   TESTING ON SMALLER SUBSETS OF DATA   ==============


# ==============   AUTO-RUN THE PROCEDURE ON THE WHOLE DATA SET   ==============
head(datBayes)

# Figure showing the cutoff for accurate prediction
glort = as.numeric(datBayes$gs_text == datBayes$naiveBayes)
xyplot(jitter(glort, amount=0.02) ~ logOdds | gs_text, data=datBayes,
	type=c('p', 'smooth'), span=1, 
	as.table=TRUE, 	par.strip.text=list(lines=0.8, cex=0.7),
	scales=list(x=list(relation='free')))



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
		panel.lmline(x, y, col='red', lty=1)
		panel.lmline(c(x,y), c(x,y), col='green', lty=2)
	})

	

# Calculate the overall CSMFAccuracy for a run.   See "Robust metrics..." p. 7+	
# In the final paper I should calculate the the median CSMFAccuracy
# for the 500 CV samplings.
calculateCSMFAccuracy(datBayes$naiveBayes, datBayes$gs_text)  

	
subset(csmf, subset=(category=='Poisonings' & cause==category))
	
	



