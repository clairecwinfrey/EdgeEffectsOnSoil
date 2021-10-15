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

# 5. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

# 6. zScore function computes the z-score for a given vector of numbers
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
#################################################################################

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

# Barplot showing ubiquity by EU- shows that vast majority of ASVs, after dropping off 
# rare ones, are present in more than one site!
barplot(table(ASVubiquityEU[,7]), ylab="Number of ASVs", xlab= "Number of EUs/sites Present in" )

#################################################################################
# II. DIFFERENTIAL ABUNDANCE ANALYSIS
#################################################################################
# In this section, I do a differential abundance analysis based on all samples put together,
# then use transect-specific Z-scores for the lines. 

# Only use ASVs that occur in at least 45 samples
namesAll_45 <- names(which(ASVubiquity[,234] >= 45)) #this gives the names of the ASVs that occur at least
# 45 times
length(namesAll_45) #4480 
namesAll_45

# Remove all of the ASVs from the phyloseq object that are not in at least ten times
allEUs_45times.ps <- prune_taxa(namesAll_45,trimmedJustsoils.ps) #4480 taxa as expected!

# Remove "edge" samples so that we can find samples that vary in edge versus forest withb new PS 
all_45timesNoEdge.ps <- subset_samples(allEUs_45times.ps, Habitat != "edge")
# Did we lose any ASVs that were only on the edge?
all_45timesNOEdge_count <- ASVsampOccurrence(all_45timesNoEdge.ps)
dim(all_45timesNOEdge_count) #4480 taxa as expected!
# length(which(NoEdge52_count[,35] == 0)) #no, none were found only on the edge that were not found elsewhere
# (but again, really rare taxa were trimmed upstream!)

# DIFFERENTIAL ABUNDANCE ANALYSIS
Deseq1_all_45times <- phyloseq_to_deseq2(all_45timesNoEdge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
set.seed(19) #is there any permutation or randomness involved below?
Deseqtested_all_45times <- DESeq(Deseq1_all_45times, test="Wald", fitType = "parametric")
DeSeq_res_all_45times <- results(Deseqtested_all_45times, cooksCutoff = FALSE)
alpha <- 0.001
sigtab_all_45times <- DeSeq_res_all_45times[which(DeSeq_res_all_45times$padj < alpha), ]
sigtab_all_45times <- cbind(as(sigtab_all_45times, "data.frame"), as(tax_table(all_45timesNoEdge.ps)[rownames(sigtab_all_45times), ], "matrix"))
head(sigtab_all_45times)
dim(sigtab_all_45times) #1696 ASVs out of the 7514 tested had a corrected p-value of less than 0.001; this was 4480 with 0.01 alpha
#View(sigtab_all_45times)

# How many samples are each of these ASVs found in and what is their overall abundance?
names_all_45times <- rownames(sigtab_all_45times) #pull out names of the ASVs
all_45timesNoEdgeASV <- ASVs_outta_ps(all_45timesNoEdge.ps)    #pull out ASV table
ASVtab_all_45times_da <- all_45timesNoEdgeASV[names_all_45times,] #get a smaller version of the ASV table that has only these ASVs
SampleCount_all_45times_da <- rowSums(ASVtab_all_45times_da > 0) #get number of soil samples that the ASV appears in
ASVtab_all_45times_da_counts <- cbind(ASVtab_all_45times_da, SampleCount_all_45times_da)
ASVabund <- rowSums(ASVtab_all_45times_da_counts[,1:210]) #get abundance of each ASV ACROSS all samples
ASVtab_all_45times_da_counts <- cbind(ASVtab_all_45times_da_counts, ASVabund)
# View(ASVtab_all_45times_da_counts) 

# BREAKING UP RESULTS BY TRANSECT:
# First, get the names of each ASV identified in the differential abundance analysis:
ASVnames_all_45times_da <- rownames(sigtab_all_45times)

# For the whole soils dataset (i.e. edges too), get ASV table with only the ASVs
# identified in diff abund analysis:
trimmedJustSoilsASVs <- ASVs_outta_ps(trimmedJustsoils.ps)
unique(rownames(trimmedJustSoilsASVs[ASVnames_all_45times_da,]) == ASVnames_all_45times_da) #the fact that
# these match means that we successfully pulled out only the ASVs of interest. 
trimmedJustSoilsDAASVs <- trimmedJustSoilsASVs[ASVnames_all_45times_da,]
dim(trimmedJustSoilsDAASVs) #1696 ASVs, as expected
trimmedJustSoilsDAASVs_t <- t(trimmedJustSoilsDAASVs) #tranpose so ASVs are columns
dim(trimmedJustSoilsDAASVs_t)

# Make a new dataframe with ASVs identified in diff abund and some metadata
trimmedJustSoilsMeta <- metadata_outta_ps(trimmedJustsoils.ps) #get metadata out of phyloseq
colnames(trimmedJustSoilsMeta)
# Check that sample order is the same in trimmedJustSoilsMeta_t and trimmedJustSoilsASVs:
unique(rownames(trimmedJustSoilsMeta) == rownames(trimmedJustSoilsDAASVs_t))
# Because the code above is true, we can cbind columns of interest from trimmedJustSoilsMeta to trimmedJustSoilsDAASVs_t
# Keep EU, transect, and meter from trimmedJustSoilsMeta
AllSoilsDAASVsandMeta <- cbind(trimmedJustSoilsDAASVs_t, trimmedJustSoilsMeta[,c(6,8,9)]) 
#View(AllSoilsDAASVsandMeta) #looks good!

# Separate by EU and by transect and then get z-scores for each ASV in each transect
#### EU 52 ####
## T ##
EU52_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "T") 
#View(EU52_T_ASVtab)
rownames(EU52_T_ASVtab) <- EU52_T_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU52_T_ASV_Zs <- zScore(EU52_T_ASVtab[,1:1696]) #don't select metadata columns
# View(EU52_T_ASV_Zs)

## B ##
EU52_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "B")
rownames(EU52_B_ASVtab) <- EU52_B_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU52_B_ASV_Zs <- zScore(EU52_B_ASVtab[,1:1696]) #don't select metadata columns

## L ##
EU52_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "L")
rownames(EU52_L_ASVtab) <- EU52_L_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU52_L_ASV_Zs <- zScore(EU52_L_ASVtab[,1:1696]) #don't select metadata columns

## R ##
EU52_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "R")
rownames(EU52_R_ASVtab) <- EU52_R_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU52_R_ASV_Zs <- zScore(EU52_R_ASVtab[,1:1696]) #don't select metadata columns

##############

#### EU 53N ####
## T ##
EU53N_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "T") 
#View(EU53N_T_ASVtab)
rownames(EU53N_T_ASVtab) <- EU53N_T_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53N_T_ASV_Zs <- zScore(EU53N_T_ASVtab[,1:1696]) #don't select metadata columns

## B ##
EU53N_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "B")
rownames(EU53N_B_ASVtab) <- EU53N_B_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53N_B_ASV_Zs <- zScore(EU53N_B_ASVtab[,1:1696]) #don't select metadata columns

## L ##
EU53N_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "L")
rownames(EU53N_L_ASVtab) <- EU53N_L_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53N_L_ASV_Zs <- zScore(EU53N_L_ASVtab[,1:1696]) #don't select metadata columns

## R ##
EU53N_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "R")
rownames(EU53N_R_ASVtab) <- EU53N_R_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53N_R_ASV_Zs <- zScore(EU53N_R_ASVtab[,1:1696]) #don't select metadata columns
##############

#### EU 54S ####
## T ##
EU54S_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "T") 
rownames(EU54S_T_ASVtab) <- EU54S_T_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU54S_T_ASV_Zs <- zScore(EU54S_T_ASVtab[,1:1696]) #don't select metadata columns

## B ##
EU54S_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "B")
rownames(EU54S_B_ASVtab) <- EU54S_B_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU54S_B_ASV_Zs <- zScore(EU54S_B_ASVtab[,1:1696]) #don't select metadata columns

## L ##
EU54S_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "L")
rownames(EU54S_L_ASVtab) <- EU54S_L_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU54S_L_ASV_Zs <- zScore(EU54S_L_ASVtab[,1:1696]) #don't select metadata columns

## R ##
EU54S_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "R")
rownames(EU54S_R_ASVtab) <- EU54S_R_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU54S_R_ASV_Zs <- zScore(EU54S_R_ASVtab[,1:1696]) #don't select metadata columns
##############

#### EU 8 ####
## T ##
EU8_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "T") 
rownames(EU8_T_ASVtab) <- EU8_T_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU8_T_ASV_Zs <- zScore(EU8_T_ASVtab[,1:1696]) #don't select metadata columns

## B ##
EU8_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "B")
rownames(EU8_B_ASVtab) <- EU8_B_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU8_B_ASV_Zs <- zScore(EU8_B_ASVtab[,1:1696]) #don't select metadata columns

## L ##
EU8_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "L")
rownames(EU8_L_ASVtab) <- EU8_L_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU8_L_ASV_Zs <- zScore(EU8_L_ASVtab[,1:1696]) #don't select metadata columns

## R ##
EU8_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "R")
rownames(EU8_R_ASVtab) <- EU8_R_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU8_R_ASV_Zs <- zScore(EU8_R_ASVtab[,1:1696]) #don't select metadata columns
##############

#### EU 53S ####
## T ##
EU53S_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "T") 
rownames(EU53S_T_ASVtab) <- EU53S_T_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53S_T_ASV_Zs <- zScore(EU53S_T_ASVtab[,1:1696]) #don't select metadata columns

## B ##
EU53S_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "B")
rownames(EU53S_B_ASVtab) <- EU53S_B_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53S_B_ASV_Zs <- zScore(EU53S_B_ASVtab[,1:1696]) #don't select metadata columns

## L ##
EU53S_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "L")
rownames(EU53S_L_ASVtab) <- EU53S_L_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53S_L_ASV_Zs <- zScore(EU53S_L_ASVtab[,1:1696]) #don't select metadata columns

## R ##
EU53S_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "R")
rownames(EU53S_R_ASVtab) <- EU53S_R_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU53S_R_ASV_Zs <- zScore(EU53S_R_ASVtab[,1:1696]) #don't select metadata columns
##############

#### EU 10 ####
## T ##
EU10_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "T") 
rownames(EU10_T_ASVtab) <- EU10_T_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU10_T_ASV_Zs <- zScore(EU10_T_ASVtab[,1:1696]) #don't select metadata columns

## B ##
EU10_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "B")
rownames(EU10_B_ASVtab) <- EU10_B_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU10_B_ASV_Zs <- zScore(EU10_B_ASVtab[,1:1696]) #don't select metadata columns

## L ##
EU10_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "L")
rownames(EU10_L_ASVtab) <- EU10_L_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU10_L_ASV_Zs <- zScore(EU10_L_ASVtab[,1:1696]) #don't select metadata columns

## R ##
EU10_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "R")
rownames(EU10_R_ASVtab) <- EU10_R_ASVtab$Meter #replace sample name with meter
# Get z-scores
EU10_R_ASV_Zs <- zScore(EU10_R_ASVtab[,1:1696]) #don't select metadata columns
##############

