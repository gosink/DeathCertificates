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
load(file='dat.RData')
load(file='Data/retList.RData')
load(file='Data/pCCCList.RData')
load(file='Data/testingGrid.RData')
load(file='Data/csmfList.RData')
load(file='Data/csmfList2.RData')

	
	
# ====================  C S M F  =================== #

#### SHOW THE STRUCTURE OF THE ORIGINAL DATA ####
source('GC Methods.R')
visualizeCausalChains(x=subset(dat, subset=(ageGroup=='Adult' & gs_text=='AIDS')))
visualizeCausalChains(x=subset(dat, subset=(ageGroup=='Adult' & gs_text=='Pneumonia')))


#### SHOW RECOVERY OF _gc CODES AND CORRECT ASSIGNMENT ####
temp   = retList$bayesDat
gcTemp = NULL
temp = merge(dat[,c('sid', 'icd_name')], temp, by='sid', sort=FALSE)
allCauses = colnames(temp)[8:34]
for (aCause in allCauses) {
	inx = temp$gs_text == aCause
	cat(aCause, sum(inx), sum(temp$matchPosition[inx]==1), '\n')
	gcTemp = rbind(gcTemp, data.frame(cause=aCause, count=sum(inx), 
									goodCallGC    = sum(temp$matchPosition == 1 & inx & temp$icd_name == '_gc'),
									goodCallNotGC = sum(temp$matchPosition == 1 & inx & temp$icd_name != '_gc'),
									badCallGC     = sum(temp$matchPosition != 1 & inx & temp$icd_name == '_gc'),
									badCallNotGC  = sum(temp$matchPosition != 1 & inx & temp$icd_name != '_gc'),
									stringsAsFactors=FALSE))
}
gcTemp$fracGCcorrected = with(gcTemp, signif(goodCallGC / (goodCallGC + badCallGC), 1))
gcTemp



#### SHOW HOW CSMF AND CSMF2 COMPARE ####
summary(lm(I(csmfAccuracy2 - csmfAccuracy) ~ ageGroup + useCausadef, data=testingGrid))
# RESULTS:   No significant difference.   Should carry on with CSMFAccuracy
	"
	csmf  <- csmfList[[3]]
	csmf2 <- csmfList2[[3]]
	temp <- data.frame(cvReplicate=csmf$cvReplicate, cause=csmf$cause, 
						CSMFAccuracy=csmf$CSMFAccuracy, CSMFAccuracy2=csmf2$CSMFAccuracy, stringsAsFactors=FALSE)
						
	axis.line = trellis.par.get()$axis.line
	axis.line$col = "transparent"
	trellis.par.set(axis.line=axis.line)
						
	xyplot(CSMFAccuracy2 ~ CSMFAccuracy, data=temp, subset=cause=='AIDS', 
		type=c('p', 'r', 'g'), pch=19, alpha=0.8, bty='n',
		main='CSMF Accuracy comparing 2 methods\nblue line = best fit,  red line = 45 equiprobability line',
		xlab='CSMF Accuracy using only the highest probability match', 
		ylab='CSMF Accuracy using fractional probability matches',
		panel = function(x,y, score, ...){
			panel.xyplot(x, y, ...)
			panel.lmline(x, y, col='blue', lty=1)
			panel.lmline(c(x,y), c(x,y), col='red', lty=2)
		})
	"
	

#### SEE IMPACT OF INCLUDING Causadef ####
summary(lm(csmfAccuracy ~ useCausadef, subset=ageGroup=='Adult', data=testingGrid))   # 2.1% boost for adults
summary(lm(csmfAccuracy ~ useCausadef, subset=ageGroup=='Child', data=testingGrid))   # no effect for child
summary(lm(csmfAccuracy ~ useCausadef, subset=ageGroup=='Neonate', data=testingGrid)) # 1.9% boost for neonates

	
str(csmfList[[3]])   # Basic CSMF, Adults with useCausadef=TRUE

# Plot the overall results	
temp <- csmfList2[[3]]
temp$cause <- gsub('\\.+', ' ', temp$cause, perl=TRUE); sort(unique(temp$cause))
temp$cause <- gsub('IHD Acute Myocardial', '       IHD Acute\nMyocardial', temp$cause, perl=TRUE); #sort(unique(temp$cause))
temp$cause <- gsub('Other Cardiovascular ', 'Other Cardiovascular\n       ', temp$cause, perl=TRUE); #sort(unique(temp$cause))
temp$cause <- gsub('Other Infectious Diseases', 'Other Infectious\n    Diseases', temp$cause, perl=TRUE); #sort(unique(temp$cause))
temp$cause <- gsub('Other Non communicable Diseases', 'Other Non-comm\n    Diseases', temp$cause, perl=TRUE); #sort(unique(temp$cause))

xyplot(predFrac ~ obsFrac | cause, data=temp, 
	type=c('p', 'g'), pch=1, alpha=0.3, #col=score, 
	main='Observed vs. Predicted CSMF in Adults',
	par.strip.text=list(lines=2.5, cex=0.6),
	scales=list(x=list(relation='same', rot=45, limits=c(-0.01,0.2)), 
				y=list(relation='same', rot=45, limits=c(-0.01,0.2))),
	as.table=TRUE,  xlab='True Cause Fraction', ylab='Estimated Cause Fraction',
	panel = function(x,y, score, ...){
		panel.xyplot(x, y, ...)
		panel.lmline(x, y, col='blue', lty=1)
		panel.lmline(c(x,y), c(x,y), col='red', lty=2)
	})



#### HOW MANY REPLICATES DO I NEED?  ####
# Basic CSMF, Adults with useCausadef=TRUE
x = subset(csmfList[[3]], cause=='AIDS')$CSMFAccuracy
cumMean = data.frame(i=1:500, mu=NA, ll=NA, ul=NA) 
for (i in 1:500) {
	x1 = sort(sample(x, i, replace=TRUE))
	cumMean[i, 'mu'] = median(x1)
	x2 = rep(NA, i)
	for (j in 1:i) {
		x2[j] = median(sample(x, i, replace=TRUE))
	}
	cumMean[i, c('ll', 'ul')] = quantile(x2, probs=c(0.05, 0.95))
}
par(bty='n')
plot(cumMean$i, cumMean$mu, ylim=c(0.75, 0.83), type='l', col='blue',
	main='CMSF Accuracy for CV sampling of adults\n(includes Causadef)',
	xlab='Number of samples', ylab='Median CSMF Accuracy'); grid()
lines(1:500, cumMean$ul, col='red', lty=2)
lines(1:500, cumMean$ll, col='red', lty=2)
abline(h=median(x), col='slategray')

		

# 1. AVERAGE CHANCE-CORRECTED CONCORDANCE OF INDIVIDUAL CAUSE ASSIGNMENT (P. 8 DISCUSSION)
# In the final paper as per "Robust metrics..." p. 8+	
# On p. 6, we should use the median.
temp = ddply(csmfList[[2]], .(cause), summarize, CCC=signif(median(CCC[is.finite(CCC)]), 2))
temp[order(temp$CCC, decreasing=TRUE),]
with(csmfList[[1]], signif(median(CCC[is.finite(CCC)]), 2))   # overall or grand median



# 2. pCCC(k)  -  (p. 8 Discussion)
summary(pCCCList[[3]])


# 3. Median CSMFAccuracy (middle left column p. 8, & Discussion)
median(csmfList[[3]]$CSMFAccuracy)   # Including Causadef
median(csmfList[[6]]$CSMFAccuracy)   # Not using Causadef


	
# 3.5. Regression coefficients (upper right column p. 8)
fitParams <- function(x,y) {
	fit       <- lm(y ~ x)
	alpha     <- coefficients(summary(lm(y ~ x)))[1,1]
	beta1     <- coefficients(summary(lm(y ~ x)))[2,1]
	r.squared <- summary(fit)$r.squared
	sigma     <- summary(fit)$sigma
	return(c(alpha=alpha, beta=beta1, sigma=sigma, r.squared=r.squared))
}

a = ddply(csmfList[[3]], .(cause), function(x) fitParams(x=x$obsFrac, y=x$predFrac))
a=cbind(CoD=a[,1], signif(a[,2:5], 3))
a[order(a$r.squared, decreasing=TRUE),]


# 4.  Full NxN matrix for the process
load('Data/retList.RData')
mat = pairsToMatrix(rowVals=retList$bayesDat$inferred, colVals=retList$bayesDat$gs_text, 
					square=TRUE, allVals=colnames(retList$probList[[1]]$QnBayes)) 
par(cex.main=0.8)
heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram='none', col=NULL, 
		margins=c(10,14), trace='none', density.info='none', key=FALSE,
		srtCol=45, cellnote=mat, notecex=0.7, notecol='black')
bringToTop();

# The following doesn't work.
mat2 <- 1/(mat + 1)
heatmap.2(mat2,  dendrogram='both', col=gray(32:0/32), 
		margins=c(10,14), trace='none', srtCol=45)
bringToTop();



