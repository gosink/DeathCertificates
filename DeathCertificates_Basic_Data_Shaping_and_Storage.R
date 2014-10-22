##############################################
#
# Start work with Abie Flaxman's (IHME) "Garbage Death Certificates".  This first
# piece of code mixes and matches data from various flat files and makes a reduced
# intelligble couple of data sets that I can use in all of the code moving forward.
#
# Details of the ICD codes:
# http://en.wikipedia.org/wiki/ICD-10
# http://en.wikipedia.org/wiki/ICD-10_Chapter_XX:_External_causes_of_morbidity_and_mortality
#
#   -John Gosink
#    7/29/14   - started project
#    9/10/14   - split off this set of methods as it's own file
#


# CMS32_DESC_LONG_DX.txt  ICD9 codes
# icd10cm_order_2015.txt  ICD10 codes


rm(list=ls())
setwd('C:/Users/John/Desktop/DeathCertificates/')
library(jsonlite)     # .JSON files
library(RJSONIO)
library(gdata)        # Excel files
library(plyr)
library(gtools)       # smartbind  (to bind dfs with different columns)
source('GC Methods.R')


# CREATE AND SAVE A EASY ACCESS DESCRIPTIONS OF THE ICD CODES.
# Comprehensive set of ICD-10 codes and descriptions from the CMS
icd10Vec <- readICDList(aFile='./Data/icd10cm_order_2015.txt')
save(icd10Vec, file='icd10Vec.RData')          # Saving the data at this step.


# CREATE AND SAVE A EASY ACCESS DESCRIPTIONS OF THE VARIOUS COLUMN HEADERS
datColDesc     <- read.csv('./Data/IHME_GC13_VA_DATA_CODEBOOK_Y2013M04D04.csv', header=TRUE, stringsAsFactors=FALSE)
datColDesc     <- datColDesc[nchar(datColDesc$variable) > 0,]
colDesc        <- datColDesc$question
names(colDesc) <- datColDesc$variable
save(colDesc, file='colDesc.RData')



# -----------------------------------------------
# READ IN THE BASIC DATA
# Now for the main data objects, start by reading in the Mexico death certificate data.   
# This supposedly has the data from 1589 Mexican subjects.   These are 'gold standard' 
# reviewed documents and 224 of them are noted as being garbage codes.
#
dat          <- read.csv('Data/mccd final_070311.csv', stringsAsFactors=FALSE);  
map_icd      <- fromJSON('./Data/map_icd.json', asText=FALSE)
dat$icd_name <- map_icd[dat$ICD]

# Now I have replicated Abie's things in the first part of:
# file:///C:/Users/John/Desktop/DeathCertificates/2014_01_30a_MEX_MCD_data.html
str(dat); summary(dat)      # 1589 obs. of  11 variables
sort(table(dat$icd_name))   # These are the International Classification of Diseases (224 are '_gc')


# Underlying causes of death.   This is actually a very large amount of information (530+ columns?)
# detailing an interview with one (or more?) surviving members of the deceased's family as well as
# a more thorough report from a pathologist(?).  A key result is the 11th field or so called "gs_text55"
# which has the updated/correct cause of death.   Compare/contrast this to "gs_text" in the first dat file.
#
# Here is the info from Bernardo:
'Abie, finally here is the info on the multiple cause of death DC coding for Mexico, with the dataset used in the 
validation study. I have saved it at:
J:\Project\VA\Archive\MEXICO MCCD
You will find 2 files. One contains only the info on the underlying cause of death and mapping to the VA 
reduced cause list. The other one has all the causes coded. Please note that in this one, the underlying 
cause of death is divided in 2 fields. The final digit of the coding is in the “desdoblado” column. Let me 
know, and if you like we can take a look at these together.'

# FOLD IN EXTENSIVE ADDITIONAL INFORMATION ABOUT EACH SUBJECT
underlyingInfo <- read.xls('./Data/Mexico_base de CD con 1638 registros.xls', sheet=1, header=TRUE, stringsAsFactors=FALSE)
underlyingInfo <- underlyingInfo[, colnames(underlyingInfo)!='X']
inx            <- is.finite(as.numeric(underlyingInfo$desdobla))
underlyingInfo[inx, 'causadef'] <- paste(underlyingInfo[inx, 'causadef'], underlyingInfo[inx, 'desdobla'], sep='')
head(underlyingInfo)

dat <- merge(dat, underlyingInfo, by.x='sid', by.y='Sid')
str(dat)

detailInfo <- NULL
for (aModule in dir('./Data/', pattern=glob2rx('IHME_GC13*csv'), full.names=TRUE)) {
#	print(aModule)
	temp       <- read.csv(aModule, header=TRUE, stringsAsFactors=FALSE)
	cat(aModule, dim(temp), '\n')
	detailInfo <- smartbind(detailInfo, temp, verbose=FALSE)
}

detailInfo <- subset(detailInfo, site=='Mexico')
dim(detailInfo)
detailInfo$sid  <- as.numeric(gsub('M-', '', detailInfo$sid))
detailInfo[1:20, 1:10]
detailInfo[sample(nrow(detailInfo),20), sample(ncol(detailInfo),10)]
length(unique(dat$sid))              # 1589 subjects...
length(unique(detailInfo$sid))       # 2031 subjects (More sids in the detailInfo)

# MERGE THE BASIC INFORMATION AND THE DETAILED INFORMATION TOGETHER
# After the next few steps dat is 1586 x 550.
dat <- merge(dat, detailInfo, by='sid')


# CLEAN UP AGE, GENDER AND OTHER THINGS..
dat$age         <- as.numeric(dat$g1_07a)
inx             <- is.na(dat$age)
dat[inx,'age']  <- as.numeric(dat[inx,'g1_06y']) - as.numeric(dat[inx,'g1_01y'])
plot(sort(dat$age)); grid()
dat$ageGroup   <- cut(dat$age, breaks=c(-2,2,16,50,200), labels=c('Infant', 'Child', 'Adult', 'Senior'))  
plot(table(dat$ageGroup))


dat$gender      <- dat$g1_05
dat[nchar(dat$gender)==0,'gender']  <- 'Unknown'

dat$education   <- dat$g1_09
dat[nchar(dat$education)==0,'education']  <- 'Unknown'
dat$education   <- factor(dat$education, levels=c('Unknown', 'No Schooling', 'Primary School', 
												'High School', 'College or Higher'))
table(dat$education)

# Tobacco use?   Probably cigarrette use is important with this dataset...
dat[sample(nrow(dat),20),c('a4_01', 'a4_02_1', 'a4_02_2', 'a4_02_3', 'a4_03', 'a4_04')]
plot(sort(dat$a4_04)); grid()
dat$cigsPerDay  <- dat$a4_04
plot(sort(dat$cigsPerDay)); grid()
dat$smokingGroup   <- cut(dat$cigsPerDay, breaks=c(-1,0,10,1000), labels=c('None', 'Light', 'Heavy'))
dat[which(dat$age < 10), 'smokingGroup'] = 'None'
table(dat$smokingGroup, useNA='always')

# Rename it to find it easier
dat$alcohol <- dat$a4_06
table(dat$alcohol, useNA='always')

# KEEP ONLY A SUBSET OF THE GIANT BLOCK OF DATA
# Reduce the columns to some manageable set.  At the top of this file I read in a tab delimited file
# that had explantions for each of the column headers.  I edited the file to include a column
# called 'discard'.  Zero means generally keep it.  1 is extra info.  2 is probably never needed by me.
# Now down to 1587 obs. of  67 variables when the following is executed.
str(datColDesc)
keepCols   <- c("sid", "age", "ageGroup", "gender", "education", "cigsPerDay", 'smokingGroup', 'alcohol',
				"source", "module.x", "gs_text", "va_code", "gs_code", "gs_assigned", 
				"gs_level.x", "ICD", "icd_name", "foliocer", "causadef", 
				"codcau11", "codcau21", "codcau31", "codcau41",  "codcau51", "site", 
				"gs_code34", "gs_text34", "va34", "gs_code46", "gs_text46", "va46", "gs_code55", 
				"gs_text55", "va55")
keepCols   <- c(keepCols, datColDesc$variable[datColDesc$discard==0])
dat        <- dat[, names(dat) %in% keepCols]
names(dat) <- gsub('\\.x|\\.y$', '', names(dat), perl=TRUE)


# SAVE AND STORE IT IN TWO DIFFERENT FORMATS
save(dat, file='dat.RData')          # Saving the data at this step.
write.csv(dat, file='dat.csv')



# ======== HOW TO GET THE DATA FROM NOW ON ==========
#load(file='dat.RData')          # Reloading the saved data at this step.
#load(file='icd10Vec.RData')
#load(file='colDesc.RData')




