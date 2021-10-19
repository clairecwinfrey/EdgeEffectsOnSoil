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

# 4. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

# 5. numOrder, from user technocrat at https://community.rstudio.com/t/re-arranging-columns-in-numerical-order/55207 
# re-orders columns in numerical order (least to greatest)
numOrder <- function(x){
  x %>% select(sort(colnames(df))) 
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

# Remove all of the ASVs from the phyloseq object that don't occur at least 45 times
allEUs_45times.ps <- prune_taxa(namesAll_45,trimmedJustsoils.ps) #4480 taxa as expected!

# Remove "edge" samples so that we can find samples that vary in edge versus forest with new PS 
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
ASVtab_all_45times_da_counts <- cbind(ASVtab_all_45times_da_counts, ASVabund) #note, minimum is in 38 samples, because we took out edges for DA
# View(ASVtab_all_45times_da_counts) 

#################################################################################
# III. Z-SCORE OF EACH ASV'S ABUNDANCE SEPARATED BY EU AND TRANSECT
#################################################################################
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
  filter(EU == "EU_52" & Transect == "T") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter)
#now, this is ordered by meter (but meter is last column)
#class(EU52_T_ASVtab)
# View(EU52_T_ASVtab)
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU52_T_ASV_Zs <- zScore(as.data.frame(t(EU52_T_ASVtab[,2:1697]))) 
colnames(EU52_T_ASV_Zs) <- EU52_T_ASVtab$Meter #these are in the correct order
# Test to make sure z-score is correct
testtest <- EU52_T_ASVtab[,2] - mean(EU52_T_ASVtab[,2])
testtesttestresults <- testtest/sd(EU52_T_ASVtab[,2])
EU52_T_ASV_Zs[1,] == testtesttestresults #TRUE, function is working


## B ##
EU52_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "B") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU52_B_ASV_Zs <- zScore(as.data.frame(t(EU52_B_ASVtab[,2:1697]))) 
colnames(EU52_B_ASV_Zs) <- EU52_B_ASVtab$Meter #these are in the correct order

## L ##
EU52_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "L") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU52_L_ASV_Zs <- zScore(as.data.frame(t(EU52_L_ASVtab[,2:1697]))) 
colnames(EU52_L_ASV_Zs) <- EU52_L_ASVtab$Meter #these are in the correct order

## R ##
EU52_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_52" & Transect == "R") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU52_R_ASV_Zs <- zScore(as.data.frame(t(EU52_R_ASVtab[,2:1697]))) 
colnames(EU52_R_ASV_Zs) <- EU52_R_ASVtab$Meter #these are in the correct order

#### EU 53N ####
## T ##
EU53N_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "T") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter)
#now, this is ordered by meter (but meter is last column)
#class(EU53N_T_ASVtab)
# View(EU53N_T_ASVtab)
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53N_T_ASV_Zs <- zScore(as.data.frame(t(EU53N_T_ASVtab[,2:1697]))) 
colnames(EU53N_T_ASV_Zs) <- EU53N_T_ASVtab$Meter #these are in the correct order
# Test to make sure z-score is correct
testtest <- EU53N_T_ASVtab[,2] - mean(EU53N_T_ASVtab[,2])
testtesttestresults <- testtest/sd(EU53N_T_ASVtab[,2])
EU53N_T_ASV_Zs[1,] == testtesttestresults #TRUE, function is working

## B ##
EU53N_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "B") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53N_B_ASV_Zs <- zScore(as.data.frame(t(EU53N_B_ASVtab[,2:1697]))) 
colnames(EU53N_B_ASV_Zs) <- EU53N_B_ASVtab$Meter #these are in the correct order

## L ##
EU53N_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "L") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53N_L_ASV_Zs <- zScore(as.data.frame(t(EU53N_L_ASVtab[,2:1697]))) 
colnames(EU53N_L_ASV_Zs) <- EU53N_L_ASVtab$Meter #these are in the correct order

## R ##
EU53N_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53N" & Transect == "R") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53N_R_ASV_Zs <- zScore(as.data.frame(t(EU53N_R_ASVtab[,2:1697]))) 
colnames(EU53N_R_ASV_Zs) <- EU53N_R_ASVtab$Meter #these are in the correct order

#### EU 54S ####
## T ##
EU54S_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "T") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter)
#now, this is ordered by meter (but meter is last column)
#class(EU54S_T_ASVtab)
# View(EU54S_T_ASVtab)
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU54S_T_ASV_Zs <- zScore(as.data.frame(t(EU54S_T_ASVtab[,2:1697]))) 
colnames(EU54S_T_ASV_Zs) <- EU54S_T_ASVtab$Meter #these are in the correct order
# Test to make sure z-score is correct
testtest <- EU54S_T_ASVtab[,2] - mean(EU54S_T_ASVtab[,2])
testtesttestresults <- testtest/sd(EU54S_T_ASVtab[,2])
EU54S_T_ASV_Zs[1,] == testtesttestresults #TRUE, function is working

## B ##
EU54S_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "B") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU54S_B_ASV_Zs <- zScore(as.data.frame(t(EU54S_B_ASVtab[,2:1697]))) 
colnames(EU54S_B_ASV_Zs) <- EU54S_B_ASVtab$Meter #these are in the correct order

## L ##
EU54S_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "L") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU54S_L_ASV_Zs <- zScore(as.data.frame(t(EU54S_L_ASVtab[,2:1697]))) 
colnames(EU54S_L_ASV_Zs) <- EU54S_L_ASVtab$Meter #these are in the correct order

## R ##
EU54S_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_54S" & Transect == "R") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU54S_R_ASV_Zs <- zScore(as.data.frame(t(EU54S_R_ASVtab[,2:1697]))) 
colnames(EU54S_R_ASV_Zs) <- EU54S_R_ASVtab$Meter #these are in the correct order

#### EU 8 ####
## T ##
EU8_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "T") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter)
#now, this is ordered by meter (but meter is last column)
#class(EU8_T_ASVtab)
# View(EU8_T_ASVtab)
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU8_T_ASV_Zs <- zScore(as.data.frame(t(EU8_T_ASVtab[,2:1697]))) 
colnames(EU8_T_ASV_Zs) <- EU8_T_ASVtab$Meter #these are in the correct order
# Test to make sure z-score is correct
testtest <- EU8_T_ASVtab[,2] - mean(EU8_T_ASVtab[,2])
testtesttestresults <- testtest/sd(EU8_T_ASVtab[,2])
EU8_T_ASV_Zs[1,] == testtesttestresults #TRUE, function is working

## B ##
EU8_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "B") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU8_B_ASV_Zs <- zScore(as.data.frame(t(EU8_B_ASVtab[,2:1697]))) 
colnames(EU8_B_ASV_Zs) <- EU8_B_ASVtab$Meter #these are in the correct order

## L ##
EU8_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "L") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU8_L_ASV_Zs <- zScore(as.data.frame(t(EU8_L_ASVtab[,2:1697]))) 
colnames(EU8_L_ASV_Zs) <- EU8_L_ASVtab$Meter #these are in the correct order

## R ##
EU8_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_8" & Transect == "R") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU8_R_ASV_Zs <- zScore(as.data.frame(t(EU8_R_ASVtab[,2:1697]))) 
colnames(EU8_R_ASV_Zs) <- EU8_R_ASVtab$Meter #these are in the correct order

#### EU 53S ####
## T ##
EU53S_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "T") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter)
#now, this is ordered by meter (but meter is last column)
#class(EU53S_T_ASVtab)
# View(EU53S_T_ASVtab)
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53S_T_ASV_Zs <- zScore(as.data.frame(t(EU53S_T_ASVtab[,2:1697]))) 
colnames(EU53S_T_ASV_Zs) <- EU53S_T_ASVtab$Meter #these are in the correct order
# Test to make sure z-score is correct
testtest <- EU53S_T_ASVtab[,2] - mean(EU53S_T_ASVtab[,2])
testtesttestresults <- testtest/sd(EU53S_T_ASVtab[,2])
EU53S_T_ASV_Zs[1,] == testtesttestresults #TRUE, function is working

## B ##
EU53S_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "B") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53S_B_ASV_Zs <- zScore(as.data.frame(t(EU53S_B_ASVtab[,2:1697]))) 
colnames(EU53S_B_ASV_Zs) <- EU53S_B_ASVtab$Meter #these are in the correct order

## L ##
EU53S_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "L") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53S_L_ASV_Zs <- zScore(as.data.frame(t(EU53S_L_ASVtab[,2:1697]))) 
colnames(EU53S_L_ASV_Zs) <- EU53S_L_ASVtab$Meter #these are in the correct order

## R ##
EU53S_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_53S" & Transect == "R") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU53S_R_ASV_Zs <- zScore(as.data.frame(t(EU53S_R_ASVtab[,2:1697]))) 
colnames(EU53S_R_ASV_Zs) <- EU53S_R_ASVtab$Meter #these are in the correct order

#### EU 10 ####
## T ##
EU10_T_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "T") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter)
#now, this is ordered by meter (but meter is last column)
#class(EU10_T_ASVtab)
# View(EU10_T_ASVtab)
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU10_T_ASV_Zs <- zScore(as.data.frame(t(EU10_T_ASVtab[,2:1697]))) 
colnames(EU10_T_ASV_Zs) <- EU10_T_ASVtab$Meter #these are in the correct order
# Test to make sure z-score is correct
testtest <- EU10_T_ASVtab[,2] - mean(EU10_T_ASVtab[,2])
testtesttestresults <- testtest/sd(EU10_T_ASVtab[,2])
EU10_T_ASV_Zs[1,] == testtesttestresults #TRUE, function is working

## B ##
EU10_B_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "B") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU10_B_ASV_Zs <- zScore(as.data.frame(t(EU10_B_ASVtab[,2:1697]))) 
colnames(EU10_B_ASV_Zs) <- EU10_B_ASVtab$Meter #these are in the correct order

## L ##
EU10_L_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "L") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU10_L_ASV_Zs <- zScore(as.data.frame(t(EU10_L_ASVtab[,2:1697]))) 
colnames(EU10_L_ASV_Zs) <- EU10_L_ASVtab$Meter #these are in the correct order

## R ##
EU10_R_ASVtab <- AllSoilsDAASVsandMeta %>% 
  filter(EU == "EU_10" & Transect == "R") %>% 
  tibble::rownames_to_column("SampleID") %>% #have to do this so that we can arrange by meter
  arrange(Meter) #order by meter
# Get z-scores. Don't select metadata columns, and flip so ASVs and other variables are rows and samples are columns
EU10_R_ASV_Zs <- zScore(as.data.frame(t(EU10_R_ASVtab[,2:1697]))) 
colnames(EU10_R_ASV_Zs) <- EU10_R_ASVtab$Meter #these are in the correct order

#################################################################################
# IV. GETTING HIGH, MEDIUM, AND LOW LOG2FOLDABUNDANCE DIFFERNTIALLY ABUNDANT ASVs
#################################################################################
# Get ASVs that have a high, medium, and low log2Fold Abundance for patch and separately
# for forest. Then pull out the Zscores from these ASVs using the made in part III above.

# First, manipulate "sigtab_all_45times" so that it can be a tibble:
# Make ASV names (i.e. rownames) into a column instead and make into a tibble
sigtab_all_45times.tb <- tibble::rownames_to_column(sigtab_all_45times, "ASV_name")

#### "PATCH" ASVs ####
# Get only ASVs with positive log fold change abundance and arrange in descending order log2foldchange
sigtab_all_45Pos <- sigtab_all_45times.tb %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange)) 
# 1. Top 10 highest logfold change (i.e. most log fold abundant in patch)
patchTop10 <- sigtab_all_45Pos[1:10,] # 10 most log fold abundant in patch

# 2. Bottom 10 starting at 0.0 (i.e. ASVs that are barely more significantly abundant in patch)
# Below, this gets the smallest log2FoldChanges
patchBottom10 <- sigtab_all_45Pos[(nrow(sigtab_all_45Pos) - 9):(nrow(sigtab_all_45Pos)),]

# 3. Middle 10 positive (i.e. ASVs that are enriched in patch to a medium degree)
# Find row of ASV that corresponds to the median, then take this + 5 and -4 for ten values
# Below, this code finds the row that has the log2FoldChange median value
medIndex1 <- with(sigtab_all_45Pos, which.min(log2FoldChange != quantile(log2FoldChange, .5, type = 1)))
# To check, get Log2FoldChange median and compare with value at row found above:
median(sigtab_all_45Pos$log2FoldChange) == sigtab_all_45Pos$log2FoldChange[medIndex1]
# This is true, so take four values less than median and five greater than for "middle 10")
patchMiddle10 <- sigtab_all_45Pos[((medIndex1 - 5):(medIndex1 + 4)),]

#######################
#### "FOREST" ASVs ####
# Most abundant in the forest
# Get only ASVs with negative log fold change abundance and sort most neg to least neg:
sigtab_all_45Neg <- sigtab_all_45times.tb %>% filter(log2FoldChange < 0) %>% arrange(log2FoldChange)  

# 4. Top 10 most negative logfold change (i.e. most log fold abundant in forest)
forestTop10 <- sigtab_all_45Neg[1:10,]

# 5. Bottom 10 starting at 0.0 and going more negative (i.e. ASVs that are barely more
# significantly abundant in forest
# Below, this gets the last 10 negative, i.e. the most positive of the negatives
forestBottom10 <- sigtab_all_45Neg[(nrow(sigtab_all_45Neg) - 9):(nrow(sigtab_all_45Neg)),]

# 6. Middle 10 negative (i.e. ASVs that are enriched in forest to medium degree)
# Find row of ASV that corresponds to the median, then take this + 5 and -4 for ten values
# Below, this code finds the row that has the log2FoldChange median value
medIndex2 <- with(sigtab_all_45Neg, which.min(log2FoldChange != quantile(log2FoldChange, .5, type = 1)))
# To check, get Log2FoldChange median and compare with value at row found above:
median(sigtab_all_45Neg$log2FoldChange) == sigtab_all_45Neg$log2FoldChange[medIndex2]
# This is true, so take four values less than median and five greater than for "middle 10")
forestMiddle10 <- sigtab_all_45Neg[((medIndex2 - 5):(medIndex2 + 4)),]

# Combine all the ASV names
HiMedLowASVs <- c(patchTop10$ASV_name, patchBottom10$ASV_name, patchMiddle10$ASV_name,
                  forestTop10$ASV_name, forestBottom10$ASV_name, forestMiddle10$ASV_name)
#################################################################################
# V. HIGH, MEDIUM, AND LOW DIFFERNTIALLY ABUNDANT ASVs BY EU
#################################################################################
# Get median z-scores (transects) by EU
# Dummy example for getting median of matrices. 
A <- matrix(c(2,4,3,5), 2)
B <- matrix(c(6,8,NA,9), 2)
C <- matrix(c(2,4,8,5), 2)

X <- list(A, B, C)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

apply(Y, c(1, 2), median, na.rm = TRUE) #this works!

#### EU 52 #####
A <- EU52_T_ASV_Zs[HiMedLowASVs,] #missing 70
# add 70 back in for A, as NAs
Adf <- as.data.frame(A)
Adf$"70" <- NA
Adf <- cbind(Adf$"10", Adf$"20", Adf$"30", Adf$"40", Adf$"50", Adf$"60", Adf$"70", Adf$"80", Adf$"90", Adf$"100")
rownames(Adf) <- rownames(A) #retrurn rownames
A <- as.matrix(Adf)
colnames(A) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100") #reset column names
A
B <- EU52_B_ASV_Zs[HiMedLowASVs,] #missing 90
# Add back in 90 m as NAs for B
Bdf <- as.data.frame(B)
Bdf$"90" <- NA
Bdf <- cbind(Bdf$"10", Bdf$"20", Bdf$"30", Bdf$"40", Bdf$"50", Bdf$"60", Bdf$"70", Bdf$"80", Bdf$"90", Bdf$"100")
rownames(Bdf) <- rownames(B) #return rownames
B <- as.matrix(Bdf)
colnames(B) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100") #reset column names
B
C <- EU52_L_ASV_Zs[HiMedLowASVs,]
D <- EU52_R_ASV_Zs[HiMedLowASVs,]

X <- list(A, B, C, D)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

EU52_topASVsMed <- apply(Y, c(1, 2), median, na.rm = TRUE)
rownames(EU52_topASVsMed) <- rownames(A) #rownames of each ASV
colnames(EU52_topASVsMed) <- colnames(A)
EU52_topASVsMed

#### EU 53N #####
A <- EU53N_T_ASV_Zs[HiMedLowASVs,] #missing 20
# add 20 back in for A, as NAs
Adf <- as.data.frame(A)
Adf$"20" <- NA
Adf <- cbind(Adf$"10", Adf$"20", Adf$"30", Adf$"40", Adf$"50", Adf$"60", Adf$"70", Adf$"80", Adf$"90", Adf$"100")
rownames(Adf) <- rownames(A) #retrurn rownames
A <- as.matrix(Adf)
colnames(A) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100") #reset column names
A
B <- EU53N_B_ASV_Zs[HiMedLowASVs,]
C <- EU53N_L_ASV_Zs[HiMedLowASVs,]
D <- EU53N_R_ASV_Zs[HiMedLowASVs,]

X <- list(A, B, C, D)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

EU53N_topASVsMed <- apply(Y, c(1, 2), median, na.rm = TRUE)
rownames(EU53N_topASVsMed) <- rownames(A) #rownames of each ASV
colnames(EU53N_topASVsMed) <- colnames(A)
EU53N_topASVsMed

#### EU 54S ####
A <- EU54S_T_ASV_Zs[HiMedLowASVs,]
B <- EU54S_B_ASV_Zs[HiMedLowASVs,]
C <- EU54S_L_ASV_Zs[HiMedLowASVs,]
D <- EU54S_R_ASV_Zs[HiMedLowASVs,]

X <- list(A, B, C, D)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

EU54S_topASVsMed <- apply(Y, c(1, 2), median, na.rm = TRUE)
rownames(EU54S_topASVsMed) <- rownames(A) #rownames of each ASV
colnames(EU54S_topASVsMed) <- colnames(A)
EU54S_topASVsMed

#### EU 8 #####
A <- EU8_T_ASV_Zs[HiMedLowASVs,]
B <- EU8_B_ASV_Zs[HiMedLowASVs,]
C <- EU8_L_ASV_Zs[HiMedLowASVs,]
D <- EU8_R_ASV_Zs[HiMedLowASVs,]

X <- list(A, B, C, D)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

EU8_topASVsMed <- apply(Y, c(1, 2), median, na.rm = TRUE)
rownames(EU8_topASVsMed) <- rownames(A) #rownames of each ASV
colnames(EU8_topASVsMed) <- colnames(A)
EU8_topASVsMed

#### EU 53S #####
A <- EU53S_T_ASV_Zs[HiMedLowASVs,]
B <- EU53S_B_ASV_Zs[HiMedLowASVs,]
C <- EU53S_L_ASV_Zs[HiMedLowASVs,]
D <- EU53S_R_ASV_Zs[HiMedLowASVs,]

X <- list(A, B, C, D)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

EU53S_topASVsMed <- apply(Y, c(1, 2), median, na.rm = TRUE)
rownames(EU53S_topASVsMed) <- rownames(A) #rownames of each ASV
colnames(EU53S_topASVsMed) <- colnames(A)
EU53S_topASVsMed

#### EU 10 #####
A <- EU10_T_ASV_Zs[HiMedLowASVs,]
B <- EU10_B_ASV_Zs[HiMedLowASVs,]
C <- EU10_L_ASV_Zs[HiMedLowASVs,]
D <- EU10_R_ASV_Zs[HiMedLowASVs,]

X <- list(A, B, C, D)
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

EU10_topASVsMed <- apply(Y, c(1, 2), median, na.rm = TRUE)
rownames(EU10_topASVsMed) <- rownames(A) #rownames of each ASV
colnames(EU10_topASVsMed) <- colnames(A)
EU10_topASVsMed
###########################################
# PLOTTING

#EU52_topASVsMed, EU53N_topASVsMed, EU54S_topASVsMed, EU8_topASVsMed, EU53S_topASVsMed, EU10_topASVsMed
quartz()
par(mfrow=c(2,3))
# this doesn't work! "'x' and 'y' must have same number of rows"--pivot_longer with more meters?
EU52_topASVsMed_plot <-  matplot(colnames(EU52_topASVsMed), EU52_topASVsMed, type = "l", xlab= "Meter",
                                 ylab= "Median Z-Score for ASV abundance", main= "Median ASV Z-Scores across EU 52")

# need to pivot longer and plot with ggplot2, I think
EU52_topASVsMedLonger <- as.data.frame(EU52_topASVsMed) %>% 
  rownames_to_column(var="ASV_name") %>% 
  pivot_longer(cols= "10":"100",
               names_to= "Meter", values_to = "Median_Zscore") 
# Make it so the meters are plotted in the correct order
EU52_topASVsMedLonger$Meter <- factor(EU52_topASVsMedLonger$Meter, levels = c("10", "20", "30", "40", "50","60", "70", "80", "90", "100"))

# why isn't y-axis showing up?
quartz()
ggplot(EU52_topASVsMedLonger, aes(x=Meter, y=Median_Zscore, color= ASV_name, group = ASV_name)) + 
  geom_line() + theme(axis.text.y = element_blank()) + ggtitle("Median ASV Z-Scores across EU 52") +
  theme(legend.position = "none") +
  scale_y_continuous("Median Z-score")
# add in this for scale_y_continuous(breaks=..., labels=...). ?

  
#################################################################################
# VI. HIGH, MEDIUM, AND LOW DIFFERNTIALLY ABUNDANT ASVs-- ACROSS ALL TRANSECTS
#################################################################################
