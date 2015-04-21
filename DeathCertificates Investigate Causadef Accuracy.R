##############################################
#
# The "Assessing quality of medical death certification" paper
# in table 2 suggests that I can get 78% median accuracy 


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


with(dat, table(gs_text, causadef))
with(dat, heatmap.2(table(gs_text, causadef), trace='none', col=gray(32:0/32)))

