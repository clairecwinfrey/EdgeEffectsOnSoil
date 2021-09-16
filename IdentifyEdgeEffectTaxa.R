## Identifying possible ASVs by EU/Site
# (started) September 14, 2021

# The purpose of this script is to identify ASVs that occur in several samples
# per EU, and that show a difference between forest and patch. That way, I can
# begin narrowing down ASVs and then lower resoltution taxonomic groups that may
# show an edge effect.

# This script uses some objects made in 16SExploratoryDataAnalysisAug2021.R, but 
# DOES NOT rely on anything created in 16SEUEDASept2021.R. Thus, if I end up keeping
# the pipeline created here, this script can be used in place of 16SEUEDASept2021.R.

########################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

# # Read in libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")
library("mctoolsr")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("DESeq2") #for differential abundance analysis
library("grid")

# Load in objects made in 16SExploratoryDataAnalysisAugust2021
load(file = "EDA16SAug2021")
# Prune samples to separate by EU and check to make sure each looks right

EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU)
EU_53N_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53N")
unique(sample_data(EU_53N_Soils.ps)$EU)
EU_54S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_54S")
unique(sample_data(EU_54S_Soils.ps)$EU)
EU_8_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_8")
unique(sample_data(EU_8_Soils.ps)$EU)
EU_53S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53S")
unique(sample_data(EU_53S_Soils.ps)$EU)
EU_10_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_10")
unique(sample_data(EU_10_Soils.ps)$EU)

# Proof of concept
# Make mock community to see if approach works
#Simulate non-negative integers simulating ASV "counts"
set.seed(19)
dat <- rpois(n=36, lambda = 1) 
# in cartoon example, as with real data above, ASVs are rows and samples are columns
datmat <- matrix(dat, nrow=9, ncol=4) 

# How many columns (samples) does each ASV (row) appear in?
ASVsampCount <- rowSums(datmat > 0)
ASVsampCount #looks correct!
cbind(datmat, ASVsampCount) 

# Define function to get number of samples that each ASV occurs in 
ASVsampOccurrence <- function(physeq) { #input is a phyloseq object
  OTU <- otu_table(physeq)
  ASVmat <- as(OTU, "matrix") # make phyloseq ASV table into a non-phyloseq matrix
  sampleOccurence <- rowSums(ASVmat > 0) #sum up number of samples ASV occurs in 
  return(cbind(ASVmat, sampleOccurence)) #
}

#############
# EU 52
#############

ASVwithCount_52 <- ASVsampOccurrence(EU_52_Soils.ps)
# View(ASVwithCount_52) taxa/ASVs re rows, and samples are columns
dim(ASVwithCount_52)
barplot(table(ASVwithCount_52[,39]), ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 38)" )
length(which(ASVwithCount_52[,39] >= 4)) #7,514 ASVs are found in at least four samples!
length(which(ASVwithCount_52[,39] >= 10)) #3506 ASVs found in at least 10
length(which(ASVwithCount_52[,39] >= 20)) #1062
length(which(ASVwithCount_52[,39] >= 30)) #310
length(which(ASVwithCount_52[,39] >= 38)) #22 ASVs found in ALL samples!

#############
# EU 53N
#############
ASVwithCount_53N <- ASVsampOccurrence(EU_53N_Soils.ps)
dim(ASVwithCount_53N)
barplot(table(ASVwithCount_53N[,40]), main= "EU 53N", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)" )
length(which(ASVwithCount_53N[,39] >= 4)) #905
length(which(ASVwithCount_53N[,39] >= 10)) #222 ASVs found in at least 10
length(which(ASVwithCount_53N[,39] >= 20)) #75
length(which(ASVwithCount_53N[,39] >= 30)) #34
length(which(ASVwithCount_53N[,39] >= 39)) #25 found in all samples!

#############
# EU 54S
#############

ASVwithCount_54S <- ASVsampOccurrence(EU_54S_Soils.ps)
dim(ASVwithCount_54S)
barplot(table(ASVwithCount_54S[,39]), main= "EU 54S", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 38)" )
length(which(ASVwithCount_54S[,39] >= 4)) #7354
length(which(ASVwithCount_54S[,39] >= 10)) #3485
length(which(ASVwithCount_54S[,39] >= 20)) #1190
length(which(ASVwithCount_54S[,39] >= 30)) #390
length(which(ASVwithCount_54S[,39] >= 38)) #59 found in all samples!

#############
# EU 8
#############

ASVwithCount_8 <- ASVsampOccurrence(EU_8_Soils.ps)
dim(ASVwithCount_8)
barplot(table(ASVwithCount_8[,40]), main= "EU 8", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)" )
length(which(ASVwithCount_8[,39] >= 4)) #937
length(which(ASVwithCount_8[,39] >= 10)) #289
length(which(ASVwithCount_8[,39] >= 20)) #110
length(which(ASVwithCount_8[,39] >= 30)) #58
length(which(ASVwithCount_8[,39] >= 39)) #43 found in all samples!

#############
# EU 53S
#############

ASVwithCount_53S <- ASVsampOccurrence(EU_53S_Soils.ps)
dim(ASVwithCount_53S)
barplot(table(ASVwithCount_53S[,41]), main= "EU 53S", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 40)" )
length(which(ASVwithCount_53S[,39] >= 4)) #854
length(which(ASVwithCount_53S[,39] >= 10)) #299
length(which(ASVwithCount_53S[,39] >= 20)) #120
length(which(ASVwithCount_53S[,39] >= 30)) #65
length(which(ASVwithCount_53S[,39] >= 40)) #43 ASVs were found in ALL samples!

#############
# EU 10
#############

ASVwithCount_10 <- ASVsampOccurrence(EU_10_Soils.ps)
dim(ASVwithCount_10)
barplot(table(ASVwithCount_10[,40]), main= "EU 10", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)" )
length(which(ASVwithCount_10[,39] >= 4)) #895
length(which(ASVwithCount_10[,39] >= 10)) #296
length(which(ASVwithCount_10[,39] >= 20)) #133
length(which(ASVwithCount_10[,39] >= 30)) #80
length(which(ASVwithCount_10[,39] >= 39)) #50

