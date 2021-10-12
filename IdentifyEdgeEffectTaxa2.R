# IdentifyEdgeEffectTaxa2 (New Differential Abundance Analysis)
# October 11, 2021 (begun)

# The purpose of this script is to identify edge effect ASVs using differntial
# abundance analyses (DESeq2). As of October 11, 2021, THIS IS DIFFERENT and more
# up to date than the analyses in IdentifyEdgeEffectTaxa.R. Specifically, this
# script was made to incorporate feedback received on Sept. 28th in lab meeting.
# Here, I return to differential abundance analsyes between forest and patch 
# samples across the whole dataset (i.e. not divided by experimental unit/site).
# UNLIKE in 16SExploratoryDataAnalysisAug2021, here I apply a ubiquity filter. 
# To plot the ASV abundance, I take the median z-score across all 24 transects
# of that ASV's abundance and plot that. Finally, I also build plots that explore
# the ubiquity of ASVs at the site/EU level. 

# FUNCTIONS DEFINED IN THIS SCRIPT:
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVsampOccurrence determines the number of samples that each ASV occurs in 
ASVsampOccurrence <- function(physeq) { #input is a phyloseq object
  OTU <- otu_table(physeq)
  ASVmat <- as(OTU, "matrix") # make phyloseq ASV table into a non-phyloseq matrix
  sampleOccurence <- rowSums(ASVmat > 0) #sum up number of samples ASV occurs in 
  return(cbind(ASVmat, sampleOccurence)) #
}

# Proof of ASVsampOccurrence concept
# Make mock community to see if approach works
#Simulate non-negative integers simulating ASV "counts"
set.seed(19)
dat <- rpois(n=36, lambda = 1) 
# in cartoon example, as with real data above, ASVs are rows and samples are columns
datmat <- matrix(dat, nrow=9, ncol=4) 

index <- which(datmat > 0) #another way of doing it, R goes down columns
datmat[index] <- 1

# How many columns (samples) does each ASV (row) appear in?
ASVsampCount <- rowSums(datmat > 0)
ASVsampCount #looks correct!
cbind(datmat, ASVsampCount) 

# 3. ASVs_outta_ps gets the ASV table out of phyloseq 
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

# 4. Handy function from Pierre L for getting rid of annoying "V" columns in data:
# Source: https://stackoverflow.com/questions/32054368/use-first-row-data-as-column-names-in-r
header.true <- function(df) { # from Pierre L: https://stackoverflow.com/questions/32054368/use-first-row-data-as-column-names-in-r
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

# 5. zScore function computes the z-score for a given vector of numbers
# The function is broadly applicable, but ASVs should be rows for its usage in this script  
zScore <- function(dat) { # input is dataframe
  zScoreDf <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
  colnames(zScoreDf) <- colnames(dat)
  rownames(zScoreDf) <- rownames(dat)
  for (i in 1:nrow(dat)){
    ASVmean <- rowMeans(dat[i,]) #if you change rowMeans to mean, could use with matrices. Or could build in if/else statment
    ASVsd <- sd(dat[i,])
    for (j in 1:length(dat[i,])) {
      zScoreDf[i,j] <- (dat[i,j]-ASVmean)/ASVsd # subtract row mean from value, then divide by standard deviation for z score
    }
  }
  return(zScoreDf)
}

# Proof/testing out zScore function to make sure that it works
set.seed(19)
dat <- rpois(n=100, lambda = 1) 
# in cartoon example, as with real data above, ASVs are rows and samples are columns
datmat <- matrix(dat, nrow=20, ncol=5) 
datmat <- as.data.frame(datmat)
# testing it out below, it seems to work
test <- zScore(datmat)
dim(test) == dim(datmat) #yes
test[1,1] == (datmat[1,1] - rowMeans(datmat[1,]))/sd(datmat[1,])
(datmat[3,4] - rowMeans(datmat[3,]))/sd(datmat[3,]) == test[3,4]
(datmat[20,2] - rowMeans(datmat[20,]))/sd(datmat[20,]) == test[20,2]

###################################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

# Read in libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("grid")

# Load data
load("trimmedJustsoils.ps") #load phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset (from 
# 16SExploratoryDataAnalysisAug2021).

#################################################################################
# I. EXPLORING UBIQUITY
##################################################################################

# 1. What is the ubiquity of ASVs (looking at all EUs)?
ASVubiquity <- ASVsampOccurrence(trimmedJustsoils.ps)
# View(ASVubiquity) taxa/ASVs re rows, and samples are columns
dim(ASVubiquity)
abund <- rowSums(ASVubiquity)
ASVubiquityAbund <- cbind(ASVubiquity, abund)
# View(ASVubiquityAbund) #this shows that 59 of these ASVs do not appear at least 50 times,
# likely this is due to removal of the controls?

# Barplot showing ubiquity
#quartz()
barplot(table(ASVubiquity[,234]), ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 233)" )

# 2. How many EUs/sites does each ASV occur in?
pooledEU.ps <- merge_samples(trimmedJustsoils.ps, "EU") #merge samples by site
rownames(sample_data(pooledEU.ps)) #nice, this is by EU!
ASVubiquityEU <- ASVsampOccurrence(t(pooledEU.ps)) #not sure why I had to transpose this...
dim(ASVubiquityEU) #dimensions are correct and is working after transposing
#View(ASVubiquityEU)

# Barplot showing ubiquity by EU- shows that vast majority of ASVs, after dropping off rare ones, are present in more
# than one site!
barplot(table(ASVubiquityEU[,7]), ylab="Number of ASVs", xlab= "Number of EUs/sites Present in" )
