# UbiquityMedianSetup.R (Set-up/data wrangling before downstream)
# November 9, 2021 (begun)

# This script picks up where 16SExploratoryDataAnalysisAug2021 leaves off,
# applying a ubiquity filter to the samples (which have already been rarefied
# and where I have removed ASVs that don't occur at least 50 times  in dataset).
# The output of this script should be fed into all downstream analyses (represented 
# in scripts IdentifyEdgeEffectTaxa4 and DissimPatterns2.R)
###################################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

# Read in libraries
library("phyloseq")
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")
library("vegan")

# Load data
load("trimmedJustsoils.ps") #load phyloseq object that was made after rarefying,
# and keeping only ASVs that occurred at least 50 times across the (rarefied) 
# dataset (from 16SExploratoryDataAnalysisAug2021).

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
# from: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  if(taxa_are_rows(ASVTable)) {
    ASVTable <- t(ASVTable)
  }
  return(as.data.frame(ASVTable))
}

# 4. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

###############################################################################
#                   UBIQUITY FILTER
###############################################################################
# MOST OF NEXT 50 LINES BELOW COPIED FROM IDENTIFYEDGEEFFECTTAXA3.R
# In this section of the script, I take a look at the ASV ubiquity across samples
# and across EUs to determine what ubiquity filter to impose (I chose that the 
# ASVs need to occur in at least 45 samples). I use this to subset the ASV table.

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
#View(ASVubiquityEU) This shows how many times each ASV appears in each EU.

# Barplot showing ubiquity by EU- shows that vast majority of ASVs, after dropping off 
# rare ones, are present in more than one site!
barplot(table(ASVubiquityEU[,7]), ylab="Number of ASVs", xlab= "Number of EUs/sites Present in" )
# This shows that the majority of ASVs occur in all of the EUs

# Based on the figures above, I will create a subset of the ASV table, based on ASVs that 
# occur in at least 45 samples.

# Only use ASVs that occur in at least 45 samples
namesAll_45 <- names(which(ASVubiquity[,234] >= 45)) #this gives the names of the ASVs that occur at least
# 45 times
length(namesAll_45) #4480 
namesAll_45

# Remove all of the ASVs from the phyloseq object that don't occur at least 45 times
postUbiquity.ps <- prune_taxa(namesAll_45,trimmedJustsoils.ps) #4480 taxa as expected!
postUbiquity.ps

###############################################################################
#       MEDIAN ASV ABUNDANCE BY METER PER EU
###############################################################################
# This code gets median abundance for each ASV (median across 4 transects in each EU) at each meter
# along the transect
ASVsdf <- ASVs_outta_ps(postUbiquity.ps)
metaDf <- metadata_outta_ps(postUbiquity.ps)

# Make a new column that combines info from meter and EU
EUmeter <- rep(NA, nrow(metaDf)) #pre-allocate vector for new names, length is 268,800
for (i in 1:nrow(metaDf)){
  EUmeter[i] <- paste(c(metaDf[i,6], metaDf[i,9]), collapse="_")
} #this should work, see example below, but is taking forever!
EUmeter

medAbund <- ASVsdf %>% #First, add columns of interest
  mutate( #add columns
    EU= metaDf$EU,#add EU column
    Transect = metaDf$Transect, #add Transect column
    Meter = metaDf$Meter,# meter column
    EUmeter = EUmeter) %>% #EU and meter information
  pivot_longer(
    cols = c(1:ncol(ASVsdf)),
    names_to = "ASVname"
  ) %>% 
  dplyr::group_by(EUmeter, ASVname) %>% 
  dplyr::summarize(
    medianEUmeter = median(value), #get median abundance for each ASV, with samples grouped by EU and meter (i.e. EU Meter)
    n = n())  #n= how many samples at meter in that EU 
  
nrow(medAbund) == ncol(ASVsdf)*60 #so this is number of ASVs times number of new categories (i.e. 10 meters in transect * 6 EUs)

# To make sure that code above is working, manually get median of a few combinations 
abundCheck <- ASVsdf %>% #First, make it longer
  mutate( #add columns
    EU= metaDf$EU,#add EU column
    Transect = metaDf$Transect, #add Transect column
    Meter = metaDf$Meter,# meter column
    EUmeter = EUmeter) %>% #EU and meter information
  pivot_longer(
    cols = c(1:ncol(ASVsdf)),
    names_to = "ASVname"
  ) %>% 
  arrange(ASVname, EUmeter) #sort so that in numerical/alphabetical order
dim(abundCheck) #dim is 1043840 rows and 6 columns 
# View(abundCheck)

# Manually check a few to make sure that median function is working properly:
# ASV 1 across EU10, 10m 
EU_10_L10 <- abundCheck[1,6]
EU_10_T10 <- abundCheck[2,6]
EU_10_R10 <- abundCheck[3,6]
EU_10_B10 <- abundCheck[4,6]
median(as.numeric(c(EU_10_L10, EU_10_T10, EU_10_R10, EU_10_B10))) == medAbund[1,3]

# ASV 998 across EU8_90m (EU_8_R 90m is msising, I think dropped in rarefaction)
EU_8_T90 <- abundCheck[1043840,6]
EU_8_B90 <- abundCheck[1043839,6]
EU_8_L90 <- abundCheck[1043838,6]
median(as.numeric(c(EU_8_T90, EU_8_B90, EU_8_L90))) == medAbund[268800, 3]

# ASV 998 across EU 53S, meter 10
EU_53S_B10 <- abundCheck[1043724,6]
EU_53S_R10 <- abundCheck[1043725,6]
EU_53S_T10 <- abundCheck[1043726,6]
EU_53S_L10 <- abundCheck[1043727,6]
median(as.numeric(c(EU_53S_B10, EU_53S_R10, EU_53S_T10, EU_53S_L10))) == medAbund[138880,3]

###############################################################################
#                          MAKE NEW PHYLOSEQ OBJECT
###############################################################################
# This phyloseq object will incorporate all the change thus far that I've made
# to the dataset. 1) rarefied (now 233 samples), 2) only ASVs occuring at least
# 50 times across dataset, 3) only ASVs that occur in at least 45/233 samples,
# and 4) median abundances for each meter in each ASV (see above)

# NEED three files for phyloseq: OTU table, taxonomy table, and metadata
#### NEW OTU TABLE ####
# For the OTU table, need ASV IDs as rows and sample IDs as columns
# Make a new column that combines info from meter and EU

# Reformat so that ASV IDs are rows and sample IDs are columns
# 1.This will serve as the OTU table for the new phyloseq object 
ASVtabMedian <- medAbund %>% 
  ungroup %>% # necessary to remove columns
  select(-c(n)) %>% #remove n
  pivot_wider(names_from = ASVname, values_from= medianEUmeter) %>% # here, would be 4481 total columns (meter name and ASV name)
  column_to_rownames("EUmeter") %>% 
  t()
# View(ASVtabMedian) #correct format and values!

# 2. Taxonomy table -- ***GET THIS OUT OF PHYLOSEQ!!***
taxTabMedian <- tax_table(postUbiquity.ps) #since we did not drop any ASVs when coding above,
# this tax table should be the same as before 
unique(sort(rownames(taxTabMedian)) == sort(rownames(ASVtabMedian))) #Statement above is true!!!

# 3. Sample metadata
metaMedian <- metadata_outta_ps(postUbiquity.ps)
head(metaMedian)
