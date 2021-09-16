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
library("grid")

# Load in objects made in 16SExploratoryDataAnalysisAugust2021
load(file = "EDA16SAug2021")
# Prune samples to separate by EU and check to make sure each looks right

EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU) #only EU 52
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
all52ASVsnames <- names(which(ASVwithCount_52[,39] >= 38))
allASVs52
length(which(ASVwithCount_52[,39] == 0)) #532 taxa do not appear in this EU, but are found at the other sites

########### ASVs THAT OCCUR AT LEAST 4 TIMES ##########
names52_4 <- names(which(ASVwithCount_52[,39] >= 4)) #this gives the names of the ASVs that occur at least
# four times
length(names52_4)

# Remove all of the ASVs from the phyloseq object that are not in at least four times
EU_52_4times.ps <- prune_taxa(names52_4, EU_52_Soils.ps) #7514 taxa as expected!

# Remove "edge" samples so that we can find samples that vary in edge versus forest
EU_52_4timesNoEdge.ps <- subset_samples(EU_52_4times.ps, Habitat != "edge")
# Did we lose any ASVs that were only on the edge?
NoEdge52_count <- ASVsampOccurrence(EU_52_4timesNoEdge.ps)
dim(NoEdge52_count)
length(which(NoEdge52_count[,35] == 0)) #no, none were found only on the edge that were not found elsewhere
# (but again, really rare taxa were trimmed upstream!)

# DIFFERENTIAL ABUNDANCE ANALYSIS
Deseq1_52_4x <- phyloseq_to_deseq2(EU_52_4timesNoEdge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
Deseqtested_52_4x <- DESeq(Deseq1_52_4x, test="Wald", fitType = "parametric")
DeSeq_res_52_4x <- results(Deseqtested_52_4x, cooksCutoff = FALSE)
alpha <- 0.01
sigtab_52_4x <- DeSeq_res_52_4x[which(DeSeq_res_52_4x$padj < alpha), ]
sigtab_52_4x <- cbind(as(sigtab_52_4x, "data.frame"), as(tax_table(EU_52_4timesNoEdge.ps)[rownames(sigtab_52_4x), ], "matrix"))
head(sigtab_52_4x)
dim(sigtab_52_4x) #107 ASVs out of the 7514 tested had a corrected p-value of less than 0.01

# How many samples are each of these ASVs found in and what is their overall abundance?
names_52_4x <- rownames(sigtab_52_4x) #pull out names of the ASVs
ASVtab_52_4x_da <- otu_table(EU_52_4timesNoEdge.ps)[names_52_4x] #get a smaller version of the ASV table that has only these ASVs
SampleCount_52_4x_da <- rowSums(ASVtab_52_4x_da > 0) #get number of soil samples that the ASV appears in
ASVtab_52_4x_da_counts <- cbind(ASVtab_52_4x_da, SampleCount_52_4x_da)
ASVabund <- rowSums(ASVtab_52_4x_da_counts[,1:34]) #get abundance of each ASV ACROSS all samples
ASVtab_52_4x_da_counts <- cbind(ASVtab_52_4x_da_counts, ASVabund)

# What do these abundances look like in the forest versus the patch?
ASV_52_names <- rownames(sample_data(EU_52_4timesNoEdge.ps))
ASV_52_names == colnames(ASVtab_52_4x_da_counts[,1:34]) #because the names are the same (and in the same order)
# in the metadata and in the new ASV tab we made above, we can use this order to get the habitat of each sample
forest_index <- which(sample_data(EU_52_4timesNoEdge.ps)$Habitat== "forest") #get samples that are in forest
patch_index <- which(sample_data(EU_52_4timesNoEdge.ps)$Habitat== "patch") #get samples that are in the patch
# Use the forest and patch indices above to get average abundance in each
forest_sum <- rowSums(ASVtab_52_4x_da_counts[,forest_index]) #gets sum across each row of ONLY the forest samples
forest_meanAbund <- forest_sum/length(forest_index) #average number in each forest sample
patch_sum <- rowSums(ASVtab_52_4x_da_counts[,patch_index])
patch_meanAbund <- patch_sum/length(patch_index)

# Do they differ among transects?

# Combine it all together
ASVtab_52_4x_da_counts <- cbind(ASVtab_52_4x_da_counts, forest_meanAbund, patch_meanAbund)
head(ASVtab_52_4x_da_counts)
View(ASVtab_52_4x_da_counts)


################ ASVs THAT OCCUR AT LEAST 10 TIMES ###### (i.e. about 1/2 of samples)
names52_10 <- names(which(ASVwithCount_52[,39] >= 10)) #this gives the names of the ASVs that occur at least
# 10 times
length(names52_10) #3506 matches calculation higher up
names52_10

# Remove all of the ASVs from the phyloseq object that are not in at least four times
EU_52_10times.ps <- prune_taxa(names52_10, EU_52_Soils.ps) #3506 taxa as expected!

# Remove "edge" samples so that we can find samples that vary in edge versus forest
EU_52_10timesNoEdge.ps <- subset_samples(EU_52_10times.ps, Habitat != "edge")
# Did we lose any ASVs that were only on the edge?
NoEdge52_count <- ASVsampOccurrence(EU_52_10timesNoEdge.ps)
dim(NoEdge52_count)
length(which(NoEdge52_count[,35] == 0)) #no, none were found only on the edge that were not found elsewhere
# (but again, really rare taxa were trimmed upstream!)

# DIFFERENTIAL ABUNDANCE ANALYSIS
Deseq1_52_10x <- phyloseq_to_deseq2(EU_52_10timesNoEdge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
Deseqtested_52_10x <- DESeq(Deseq1_52_10x, test="Wald", fitType = "parametric")
DeSeq_res_52_10x <- results(Deseqtested_52_10x, cooksCutoff = FALSE)
alpha <- 0.01
sigtab_52_10x <- DeSeq_res_52_10x[which(DeSeq_res_52_10x$padj < alpha), ]
sigtab_52_10x <- cbind(as(sigtab_52_10x, "data.frame"), as(tax_table(EU_52_10timesNoEdge.ps)[rownames(sigtab_52_10x), ], "matrix"))
head(sigtab_52_10x)
dim(sigtab_52_10x) #107 ASVs out of the 7514 tested had a corrected p-value of less than 0.01
View(sigtab_52_10x)

# How many samples are each of these ASVs found in and what is their overall abundance?
names_52_10x <- rownames(sigtab_52_10x) #pull out names of the ASVs
ASVtab_52_10x_da <- otu_table(EU_52_10timesNoEdge.ps)[names_52_10x] #get a smaller version of the ASV table that has only these ASVs
SampleCount_52_10x_da <- rowSums(ASVtab_52_10x_da > 0) #get number of soil samples that the ASV appears in
ASVtab_52_10x_da_counts <- cbind(ASVtab_52_10x_da, SampleCount_52_10x_da)
ASVabund <- rowSums(ASVtab_52_10x_da_counts[,1:34]) #get abundance of each ASV ACROSS all samples
ASVtab_52_10x_da_counts <- cbind(ASVtab_52_10x_da_counts, ASVabund)

# What do these abundances look like in the forest versus the patch?
# Use the forest and patch indices created above (when looking at 4x) to get average abundance in each
forest_sum <- rowSums(ASVtab_52_10x_da_counts[,forest_index]) #gets sum across each row of ONLY the forest samples
forest_meanAbund <- forest_sum/length(forest_index) #average number in each forest sample
patch_sum <- rowSums(ASVtab_52_10x_da_counts[,patch_index])
patch_meanAbund <- patch_sum/length(patch_index)

# Do they differ among transects?

# Combine it all together
ASVtab_52_10x_da_counts <- cbind(ASVtab_52_10x_da_counts, forest_meanAbund, patch_meanAbund)
head(ASVtab_52_10x_da_counts)
#View(ASVtab_52_10x_da_counts)

# Only consider ASVs that appear at least 50 times across whole matrix
above_50 <- which(ASVtab_52_10x_da_counts[,36] >= 50)
ASVtab_52_10x_da_50 <- ASVtab_52_10x_da_counts[above_50,]
# View(ASVtab_52_10x_da_50)

# We want to plot these ASV abundances along the transect, but first we need to get separate matrices by transect
# Initial step is to separate out by transect by getting index for each transect
EU_52_Tindex <- which(sample_data(EU_52_10timesNoEdge.ps)$Transect=="T") # "top" transect
EU_52_Bindex <- which(sample_data(EU_52_10timesNoEdge.ps)$Transect=="B") # "bottom" transect
EU_52_Lindex <- which(sample_data(EU_52_10timesNoEdge.ps)$Transect=="L") # "left" transect
EU_52_Rindex <- which(sample_data(EU_52_10timesNoEdge.ps)$Transect=="R") # "right" transect

# Separate out ASVtab_52_10x_da_counts by transect (and keep columns on abundances and counts)
ASVtab_52_10x_T <- ASVtab_52_10x_da_50[,c(EU_52_Tindex,35:38)] 
# #replace sample number with Meter (this is fine b/c sample_data has 
# same sample order as ASVtab_52_10x_da_50)
colnames(ASVtab_52_10x_T)[1:8] <- sample_data(EU_52_10timesNoEdge.ps)$Meter[EU_52_Tindex]
ASVtab_52_10x_T # top transect 





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

