##############################################
#
# As with the work on the Mexico death certificates, now make
# causal chain graphs for the New York data set
#
# HERE IS SOMEONE WHO HAS MOST OF AN R PACKAGE FOR WORKING WITH ICD CODES
#    http://rpackages.ianhowson.com/cran/icd9/
#
#   -John Gosink
#    9/22/14
#


# Details of the ICD codes:
# http://en.wikipedia.org/wiki/ICD-10
# http://en.wikipedia.org/wiki/ICD-10_Chapter_XX:_External_causes_of_morbidity_and_mortality
# CMS32_DESC_LONG_DX.txt  ICD9 codes
# icd10cm_order_2015.txt  ICD10 codes


rm(list=ls())
setwd('C:/Users/John/Desktop/DeathCertificates/')
library(gdata)        # Excel files
library(lattice)
library(latticeExtra)
library(plyr)
library(gplots)
library(gtools)       # smartbind  (to bind dfs with different columns)

source('GC Methods.R')

# -----------------------------------------------
# Read in the 30 or so csv files with the New York death certificates
# for ~498,000 people.
nyDat <- NULL
for (aFile in dir('./Data/', pattern=glob2rx('36-*.csv'), full.names=TRUE)) {
	temp    <- read.csv(aFile, header=TRUE, stringsAsFactors=FALSE, colClasses='character')
	cat(aFile, dim(temp), '\n')
	nyDat   <- smartbind(nyDat, temp, verbose=FALSE)
}
nyDat <- cbind(id=1:nrow(nyDat), nyDat)
#ny.bak = nyDat
for (i in c(1:14,16,37,58,59)) {
	nyDat[,i] <- as.numeric(nyDat[,i])
}
str(nyDat)  # 498,218  x 58   !

#save(nyDat, file='nyDat.RData')                 # Saving the data at this step.
load('nyDat.RData')
smallDat <- nyDat[sample(nrow(nyDat), 2000),c(1:27, 58:59)]   # work with just a part of it
str(smallDat)

# ---------------------------------------
# Read in the ICD-9 info.  
icd9Vec <- readICDList(aFile='./Data/icd9_CMS32_DESC_LONG_DX.txt', fileFormat=9) 



# --- ASSEMBLE AND STRUCTURE THE DATA TO A TALL/THIN DATAFRAME WITH KEY INFORMATION ---
keepCols   <- c('id', 'data_year', 'age', 'sex', 'occurence_county_fips', 'marital', 'underlying_cause')
entityCols <- grep('cause_entity', names(smallDat), value=TRUE)
keepCols   <- c(keepCols, entityCols)
temp       <- smallDat[, which(names(smallDat) %in% c(keepCols, entityCols))]
tallDat    <- reshape(temp, idvar='id', varying=entityCols, v.names='MCD', direction='long')
tallDat    <- tallDat[nchar(tallDat$MCD)>0,]
tallDat$MCD <- factor(tallDat$MCD)

# To do, look up the ICD codes.   Note that underlying_cause needs to have a zero added.
aPatient = 139963
tallDat[sample(nrow(tallDat),20),]; smallDat[smallDat$id==aPatient,]; tallDat[tallDat$id==aPatient,]
icdGuess(tallDat$MCD[tallDat$id==aPatient])


# Make graphs with most common causes
sort(table(tallDat$underlying_cause), decreasing=TRUE)
pdf(file='MCD Paths from the New York Data Set.pdf', width=11, height=8.5)
i <- 1
plotLayout <- data.frame(aCol=rep(1:2, times=20), aRow=rep(c(1,1,2,2), times=5), newpage=rep(c(T,F,F,F), times=5))
for (aCause in names(sort(table(tallDat$underlying_cause), decreasing=TRUE))[1:20]) {
#	aCause = '4280'
	causeBucket <- names(sort(table(tallDat$underlying_cause), decreasing=TRUE)[1:20])
	idBucket   <- with(smallDat, unique(id[ underlying_cause %in% causeBucket ]))    # everyone who has this gs_text
	idBucket   <- with(smallDat, unique(id[ underlying_cause == aCause ]))    # everyone who has this gs_text
	temp3      <- subset(tallDat, subset = ((id %in% idBucket)))

	yLabels <- as.character(count(temp3$MCD)$x)
	yLabels <- sample(yLabels, size=min(5,length(yLabels)))  
	yAt     <- match(yLabels, levels(temp3$MCD)) 
	
	ageGroups   <- cut(temp3$age, breaks=c(-2,2,16,50,200), labels=c('Infant', 'Child', 'Adult', 'Senior'))  
	colGroups   <- c('blue', 'red')[temp3$sex] 
	pointGroups <- c(49:52)[as.numeric(ageGroups)]     # pch=48:57  are '0' -> '9'
	codName     <- paste(substr(unlist(strsplit(icdGuess(aCause), split=' \\| ')), 1, 50), collapse=' -or-\n')
	codName     <- gsub('\\? ', '', codName)
	codName     <- substr(unlist(icdGuess(temp3$underlying_cause)), 1, 15)

	xyplot(as.numeric(MCD) ~ jitter(time, amount=0.1) | codName, groups=id, data=temp3, 
				type=c('b','g'),
				sub=list('\nBlue=male, Red=Female\n1=Infant, 2=Child, 3=Adult, 4=Senior', cex=0.7),
				col=colGroups, pch=pointGroups, cex=0.8, alpha=0.3,
				scales=list(y=list(at=yAt, labels=yLabels, relation='same')),
				xlab='< < <== Direction of Causality', ylab='ICD-9 Code?'); # bringToTop()

	
	aPlot = xyplot(as.numeric(MCD) ~ jitter(time, amount=0.1), groups=id, data=temp3, type=c('b','g'),
				main=list(paste('Cases with underlying_cause:', codName, sep='\n'), cex=0.8),
				sub=list('\nBlue=male, Red=Female\n1=Infant, 2=Child, 3=Adult, 4=Senior', cex=0.7),
				col=colGroups, pch=pointGroups, cex=0.8, alpha=0.6,
				scales=list(y=list(at=yAt, labels=yLabels, relation='same')),
				xlab='< < <== Direction of Causality', ylab='ICD-9 Code?'); # bringToTop()
	print(aPlot, split=list(plotLayout[i,'aCol'], plotLayout[i,'aRow'], 2, 2),
			newpage=plotLayout[i,'newpage'])
	i = i + 1
}
dev.off()


# --- GRAPH THE DATA IN VARIOUS WAYS ---   

#pdf(file='MCD Paths from the Mexico Data Set.pdf', width=11, height=8.5)
for (aText in sort(unique(temp2$gs_text))) {
	idBucket   <- with(temp2, unique(sid[ gs_text==aText ]))    # everyone who has this gs_text
	temp3      <- subset(temp2, subset = ((sid %in% idBucket) & (as.numeric(MCD) > 1)))

	yLabels <- as.character(count(temp3$MCD)$x)
	yLabels <- sample(yLabels, size=min(5,length(yLabels)))  
	yAt     <- match(yLabels, levels(temp3$MCD)) 

	ageGroups   <- cut(temp3$age, breaks=c(-2,2,16,50,200), labels=c('Infant', 'Child', 'Adult', 'Senior'))  
	colGroups   <- c('blue', 'red')[as.numeric(as.factor(temp3$gender))]  # display.brewer.all()
	pointGroups <- c(49:52)[as.numeric(ageGroups)]                        # pch=48:57  are '0' -> '9'

	aPlot <- xyplot(as.numeric(MCD) ~ jitter(causeOrder, amount=0.1) | icd_name, groups=sid, data=temp3, type=c('b','g'),
					main=paste('All patients with "gs_text" (died from):', aText, '\n(panel by icd_name)'),
					sub=list('\nBlue=male, Red=Female\n1=Infant, 2=Child, 3=Adult, 4=Senior', cex=0.7),
					col=colGroups, pch=pointGroups, cex=1,
					scales=list(y=list(at=yAt, labels=yLabels, relation='same')),
					xlab='< < <== Direction of Causality', ylab='ICD-10 Code'); #bringToTop()
	print(aPlot)
}		
#dev.off()	
		
		
