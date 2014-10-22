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
#library(jsonlite)     # .JSON files
#library(RJSONIO)
#library(gdata)        # Excel files
library(lattice)
library(latticeExtra)
library(plyr)
library(gplots)
library(gtools)       # smartbind  (to bind dfs with different columns)
#library(RWeka)        # Various machine learning implementations
library(expm)         # matrix exponentiation
library(edgeR)        # Good-Turing estimator background for transition matrix

source('GC Methods.R')

# -----------------------------------------------
# Start by reading in the Mexico death certificate data.   This supposedly has the data from 1589 Mexican
# subjects.   These are 'gold standard' reviewed documents and 224 of them are noted as being garbage codes.
# Prep the data in the other script file, 'DeathCertificates_Basic_Data_Shaping_and_Storage.R'
load(file='dat.RData')          # Load the (1587x67) data as processed in earlier
load(file='icd10Vec.RData')     # 91,737 ICD-10 definitions
load(file='colDesc.RData')      # Descriptions of the column headers (547 headers as a named vector)


# So now take a look at the data before and after classifications.
# There are about 35 distinct gs_text values 
#                 52 distinct gs_text55 values
#                124 distinct icd_name values
table(dat$gs_text); table(dat$gs_text55); table(dat$icd_name)
dat[sample(nrow(dat), 20), c('sid', 'gs_text', 'gs_text55', 'icd_name')]


# What are typical causal chains for cause of death?
causeCols <- c('codcau11', 'codcau21', 'codcau31', 'codcau41', 'codcau51')
temp      <- ddply(detailInfo, ~gs_code55, summarize, description=unique(gs_text55))
map_icd2  <- c(temp$description); names(map_icd2) <- temp$gs_code55
apply(dat[529, causeCols], 2, function(x) map_icd2[x])



# --- ASSEMBLE AND STRUCTURE THE DATA TO A TALL/THIN DATAFRAME WITH KEY INFORMATION ---
# Note that I fold in causadef and codcau1..5 into the causal chain.   As we see in 
# a minute, that is fine.   That is, 'causadef' terms are NOT their own unique category
# of terms.  Note a little below that, however, that the _gc category of causadef terms
# is disjoint from the non _gc set of terms.   See ... communication with Abie.
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", "gender", "education", "cigsPerDay")))
temp    <- dat[, goodCols]
temp2   <- reshape(temp, idvar='sid', varying=list(4:9), v.names='MCD', direction='long')

temp2$MCDfull <- icdGuess(temp2$MCD)
names(temp2)  <- c('sid', 'gs_text', 'icd_name', 'age', 'gender', 'education', 'cigsPerDay', 'causeOrder', 'MCD', 'MCDfull')
temp2$MCD     <- factor(gsub(' ', '', temp2$MCD))
str(temp2); temp[temp$sid==12,]; temp2[temp2$sid==12,] 



# Compare usage of codes by position in the causal chain:
tail(sort(table(temp2$MCD[temp2$causeOrder==1])), 20)
tail(sort(table(temp2$MCD[temp2$causeOrder==2])), 20)
tail(sort(table(temp2$MCD[temp2$causeOrder==3])), 20)
tail(sort(table(temp2$MCD[temp2$causeOrder>1])),  20)

# Do the names in earlier causal states even show up in the final state?
#  Yes 23%, 40%, 49%, 35%, 25% of final cause of death terms appear in positions 2-6
#  Thus 75% of the final cause of death terms show up as one or more antecedant causes.
sum(unique(as.character(temp2$MCD[temp2$causeOrder==1])) %in% unique(as.character(temp2$MCD[temp2$causeOrder==2])))
sum(unique(as.character(temp2$MCD[temp2$causeOrder==1])) %in% unique(as.character(temp2$MCD[temp2$causeOrde>1])))

# AHA!!!! JUNK DIAGNOSES ARE COMPLETELY DIFFERENT FROM GOOD DIAGNOSES.  THIS IS THE ARBITRATOR
# OF WHETHER SOMETHING IS JUNK!!!  THIS COMPLETELY SEPARATES THE TWO CLASSES!!!!
junk <- unique(as.character(temp2$MCD[temp2$causeOrder==1 & temp2$icd_name=='_gc']))  # junk diagnoses
good <- unique(as.character(temp2$MCD[temp2$causeOrder==1 & temp2$icd_name!='_gc']))  # good diagnoses
junk %in% good   # 0 of 73 junk codes are found in the good codes

# Double check, do this another way:
inx <- dat$icd_name=='_gc'
junk <- unique(gsub(' ', '', dat$causadef[inx]))
good <- unique(gsub(' ', '', dat$causadef[!inx]))


# 75% of ICD codes that lead to junk diagnoses are found in good chains
junk2 <- unique(as.character(temp2$MCD[temp2$causeOrder>1 & temp2$icd_name=='_gc']))  # junk diagnoses
good2 <- unique(as.character(temp2$MCD[temp2$causeOrder>1 & temp2$icd_name!='_gc']))  # good diagnoses
junk2 %in% good2
sum(junk2 %in% good2)/length(junk2)  



# --- GRAPH THE DATA IN VARIOUS WAYS ---   
#  Too many causes of death to see them all listed as labels on the y-axis.  
#  Instead, pick a representative set.
yLabels = names(table(temp2$MCD)[ table(temp2$MCD) > 60 ])[-1]
yAt     = match(yLabels, levels(temp2$MCD)) 
yAxis   = list(at=yAt, labels=yLabels)

# Plot 1.  More or less the raw data
groupings <- cut(temp2$age, 4)  
groupings <- temp2$gender
groupings <- temp2$education
groupings <- factor(temp2$icd_name == '_gc', labels=c('Not_GC', '_GC'))  
colGroups <- brewer.pal(n=5, name="Dark2")[as.numeric(groupings)]
xyplot(as.numeric(MCD) ~ jitter(causeOrder, amount=0.01) | groupings, groups=sid, data=temp2, type=c('b'),
		subset=(as.numeric(MCD) > 1) & sid < 400,
		par.settings=simpleTheme(col=unique(colGroups), col.line=colGroups, 
								col.points=colGroups, border='white', alpha=0.3),
		scales=list(y=list(at=yAt, labels=yLabels, relation='same')),
		xlab='< < <== Direction of Causality', ylab='ICD-10 Code',
		main='Multiple Causes of Death as observed in the data',
		auto.key=list(text=levels(groupings), space='right', rectangles=TRUE))



# !!HERE IS A KEY FIGURE!!
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
dev.off()	
		
		
				
# Plot 2.  Condense all 'cause' lines to one bucket
temp3 <- NULL
for (aCol in 4:8) {
	temp3 <- rbind(temp3, data.frame(causeOrder=1, MCD=gsub(' ', '', temp[,'causadef']), 
								sid=temp[,'sid']+(aCol-3)/10,  icd_name=temp$icd_name,
								age=temp$age, gender=temp$g1_05, stringsAsFactors=TRUE))
	temp3 <- rbind(temp3, data.frame(causeOrder=2, MCD=gsub(' ', '', temp[, aCol]),      
								sid=temp[,'sid']+(aCol-3)/10, icd_name=temp$icd_name,
								age=temp$age, gender=temp$g1_05, stringsAsFactors=TRUE))
}

temp3$MCD    <- factor(gsub(' ', '', temp3$MCD))
temp3$gender <- factor(temp3$gender)
str(temp3); temp[temp$sid==12,]; temp3[round(temp3$sid)<12,] 

groupings <- as.factor(temp3$gender)
groupings <- cut(temp3$age, breaks=c(-2,2,16,50,200), labels=c('Infant', 'Child', 'Adult', 'Senior'))  
groupings <- factor(temp3$icd_name == '_gc', labels=c('Not_GC', '_GC'))
colGroups <- brewer.pal(n=5, name="Dark2")[as.numeric(groupings)]
temp3$groupings <- groupings
temp3$colGroups <- colGroups; rm(groupings, colGroups)
xyplot(as.numeric(MCD) ~ jitter(causeOrder, amount=0.01) | groupings, groups=sid, data=temp3, type=c('b'),
		subset=(as.numeric(MCD) > 1) & sid < 4000, 
		as.table=TRUE, col.line='black', alpha=0.2,
		scales=list(y=list(at=yAt, labels=yLabels, relation='same'), x=list(draw=FALSE)),
		xlab='< < <== Direction of Causality', ylab='ICD-10 Code',
		main='Multiple Causes of Death as with "codcau" data grouped together')



#--- CAUSAL CHAINS ---#
# See 'discrete hidden Markov model (DHMM)   http://www.robots.ox.ac.uk/~vgg/rg/slides/hmm.pdf
# See also:  http://en.wikipedia.org/wiki/Bayesian_inference#Multiple_observations
#   Q - each step in the causal chain is included in the matrix  [from, to]
#   QnBayes - only the terminal events are 

source('GC Methods.R')
Q = getMCDTransitionMatrix(x=dat[sample(nrow(dat),40),])$QnBayes
Q = getMCDTransitionMatrix(x=dat[dat$icd_name=='_gc',])$Q
heatmap.2(apply(Q, 1:2, asinh), trace='none', col=rev(gray(0:18/18)), ylab='Antecedant cause...', xlab='...Leads to')
heatmap.2(Q %^% 5, trace='none', ylab='Antecedant cause...', xlab='...Leads to')


# Explore around for a while and come to the conclusion that the 
# spaces in the causa- chains are irrelevant.   
dat[sample(nrow(dat),20), c(1, 10, 12:17, 28:29, 67:70)]  # Most useful information


Q = getMCDTransitionMatrix(x=dat[dat$icd_name!='_gc',])$Q  # Good causal chains only
Q[sample(nrow(Q),20), sample(nrow(Q),20)]                  # Take a look at a part of it  


# What is the distribution of probabilites of all observed non garbage paths?
# What is the probability of an observed causal path for a single patient?
#   - for these types of stats the distributions are highly skewed.  
#       - Use rank order & Wilcox & etc?
source('GC Methods.R')
dat[1504, c(1, 10, 12:17, 28:29, 67:70)]  # Lung cancer patient
dat[1571, c(1, 10, 12:17, 28:29, 67:70)]  # CVD patient
safeMatLookup(Q, 'J180', 'J960')
evaluateSubject(x=dat[1504,], Q=Q) 
goodProbs <- data.frame(t(apply(dat[dat$icd_name!='_gc',], 1, function(x) evaluateSubject(x=x, Q=Q))))
# As expected, the longer the chain, the less probable the result.  Also, the chains
# are, as expected, high skewed.
densityplot(~finalProb^0.1 | factor(numCauses), data=goodProbs, 
			type=c('g'), ref=FALSE, as.table=TRUE, layout=c(1,6),
			scales=list(y=list(relation='free')))



# 10/8/14 THE FOLLOWING SEEMS TO WORK!!!			
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", "gender", "education", "cigsPerDay")))
source('GC Methods.R')
someIds <- sample(dat$sid, 100)
dat[dat$sid %in% someIds[10:15], goodCols]
probObj <- getMCDTransitionMatrix(x=dat[dat$sid %in% someIds,])
str(probObj)
plot(sort(probObj$codVec))

evidence <- c('G936', 'G610', 'B24')  # sid 612 died of AIDS
evidence <- c('P285', 'P280', 'P220', 'P072')  # sid 1361 died of Preterm delivery
evidence <- c('I509')  # sid 643 died of IHD - Acute Myocardial Infarction
evidence <- c('J960', 'J81', 'N189', 'E119', 'I10')
results  <- bayesCompute(evidence, probObj);  str(results)
tail(sort(results$eStepVec))
tail(sort(results$naiveBayesVec))

# FULL SCALE RUN.
source('GC Methods.R')
goodCols <- c(1,4,10,12:17, which(names(dat) %in% c("age", 'ageGroup',"gender", "education", 
													"cigsPerDay", 'smokingGroup')))
someIds <- sample(dat$sid[dat$icd_name != '_gc'])
dat[dat$sid %in% someIds[10:15], goodCols]


# Note here I am also ingesting the 'causadef' field.   Above I saw that causadef for _gc and non_gc
# patients were completely disjoint sets.
dat[1:5, goodCols]
source('GC Methods.R')
x=dat[dat$sid %in% someIds,]
probObj  <- getMCDTransitionMatrix(x=x, priorFunction='Good-Turing',
				mcdCols=c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51"),
				finalCause="gs_text", cofactors=c('gender', 'ageGroup', 'smokingGroup'))

junk = probObj$QnBayes[sample(nrow(probObj$QnBayes),30), sample(ncol(probObj$QnBayes),7)]
str(probObj);  signif(junk,2);
visualizeTransitionMatrix(probObj$QnBayes)

source('GC Methods.R')
results  <- bayesCompute(evidence=c("C924","J960","R048","C924", 'gender.Female'), probObj);  str(results)
head(sort(results$naiveBayesVec, decreasing=TRUE))
showBayesCompute(results)

source('GC Methods.R')
mcdCols=c("codcau11", "codcau21", "codcau31", "codcau41", "codcau51")

# A minifunction for ddply that takes in 1 line of dat
bestCall.2 <- function(x, mcdCols, cofactors, ...) {
	evidence = unlist(x[ mcdCols ])
	for (aCofactor in cofactors) {
		evidence <- c(evidence, cleanICDs(paste(aCofactor, as.character(x[,aCofactor]), sep='.')))
	}
	result = bayesCompute(evidence=evidence, ...)$naiveBayesVec
	return(names(sort(result, decreasing=TRUE))[1])
}
bestCall.2(x=dat[16,], mcdCols=mcdCols, cofactors=c('ageGroup', 'gender'), probObj=probObj)



# Again, note here I'm doing predictions with the causadef field
someRows <- sample(nrow(dat),50)
someRows <- which(dat$icd_name=='_gc')

datBayes <- dat[,goodCols]
datBayes$naiveBayes <- NA
for (i in 1:nrow(datBayes)) {
	datBayes[i,'naiveBayes'] <- bestCall.2(x=datBayes[i,], mcdCols=mcdCols, cofactors=c('ageGroup', 'gender'), probObj=probObj)
}										
datBayes <- cbind(datBayes[,1:3], naiveBayes=datBayes[,ncol(datBayes)], datBayes[,4:(ncol(datBayes)-1)])
head(datBayes)

write.csv(datBayes, file='datBayes.csv', quote=FALSE)






