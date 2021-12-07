# IdentifyEdgeEffectTaxa4 (New Differential Abundance Analysis)
# November 9, 2021 (begun)

# This script identifies edge effect ASVs using differential abundance analyses
# (DESeq2). THIS IS DIFFERENT and more up to date than the earlier draft 
# analyses in eitherIdentifyEdgeEffectTaxa.R, IdentifyEdgeEffectTaxa2, or 
# IdentifyEdgeEffectsTaxa3 (which is kinda a scratch page of sorts). This uses
# the median ASV abundance across transects (i.e. "top", "bottom",
# "left" and "right") at each meter point along the transect for each EU, i.e.
# "medianEU.ps" calculated in UbiquityMedianSetup.R

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
load("medianEU.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45, and 4) finally getting the 
# median ASV abundance at each meter in each EU (median is four values from each
# transect)

# FUNCTIONS DEFINED IN THIS SCRIPT:
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVsampOccurrence determines the number of samples that each ASV occurs in 
# For proof that ASVsampOccurrence works, see UbiquityMedianSetup.R
ASVsampOccurrence <- function(physeq) { #input is a phyloseq object
  OTU <- otu_table(physeq)
  ASVmat <- as(OTU, "matrix") # make phyloseq ASV table into a non-phyloseq matrix
  sampleOccurence <- rowSums(ASVmat > 0) #sum up number of samples ASV occurs in 
  return(cbind(ASVmat, sampleOccurence)) #
}

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

#################################################################################
# I. DIFFERENTIAL ABUNDANCE ANALYSIS
#################################################################################
# In this section, I do a differential abundance analysis on patch versus forest 
# samples based on the final ASV table created in UbiquityMedianSetup.R, i.e., 
# that in "medianEU.ps". 

# Remove "edge" samples so that we can find samples that vary in edge versus forest with new PS 
medianEUNoEdge.ps <- subset_samples(medianEU.ps, Habitat != "edge") #CHECK, DOES THIS REMOVE ASVS THAT ARE NO LONGER PRESENT? I THINK NOT (SEE DISSIMILARITIESACROSSTRANSECTS.R FOR CODE?)
# Did we lose any ASVs that were only on the edge?
medianEUNoEdge_count <- ASVsampOccurrence(medianEUNoEdge.ps)
dim(medianEUNoEdge_count) #4480 taxa as expected! 
length(medianEUNoEdge_count[,35] == 0) #no, none were found only on the edge that were not found elsewhere
# (but again, really rare taxa were trimmed upstream!)

# DIFFERENTIAL ABUNDANCE ANALYSIS
Deseq1_medianEU <- phyloseq_to_deseq2(medianEUNoEdge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
set.seed(19) #is there any permutation or randomness involved below?
Deseqtested_medianEU <- DESeq(Deseq1_medianEU, test="Wald", fitType = "parametric") #SHOULD parametric be the right fitType?
DeSeq_res_medianEU <- results(Deseqtested_medianEU, cooksCutoff = FALSE)
alpha <- 0.001
sigtab_medianEU <- DeSeq_res_medianEU[which(DeSeq_res_medianEU$padj < alpha), ]
sigtab_medianEU <- cbind(as(sigtab_medianEU, "data.frame"), as(tax_table(medianEUNoEdge.ps)[rownames(sigtab_medianEU), ], "matrix"))
head(sigtab_medianEU)
dim(sigtab_medianEU) #678 ASVs out of the 4480 tested had a corrected p-value of less than 0.001; this was 1036 with 0.01 alpha
# (SEE BELOW)
#View(sigtab_medianEU)
class(sigtab_medianEU)

# If we make alpha less restrictive, i.e. 0.01
alpha01 <- 0.01
sigtab_medianEU01 <- DeSeq_res_medianEU[which(DeSeq_res_medianEU$padj < alpha01), ]
sigtab_medianEU01 <- cbind(as(sigtab_medianEU01, "data.frame"), as(tax_table(medianEUNoEdge.ps)[rownames(sigtab_medianEU01), ], "matrix"))
head(sigtab_medianEU01)
dim(sigtab_medianEU01) #1036 ASVs remaining... TOO high!

# For now, pull out top 15 and bottom 15 ASVs BY logfold change
sigtab_medEU_ordered <- sigtab_medianEU %>% 
  arrange(log2FoldChange)
topDAforest <- sigtab_medEU_ordered[c(1:15),]
topDApatch <- sigtab_medEU_ordered[(nrow(sigtab_medEU_ordered)-14):nrow(sigtab_medEU_ordered),]
dim(topDApatch) #15 long

top15_DeSeq <- rbind(topDAforest, topDApatch) #most negative are forest, most positive are patch
dim(top15_DeSeq) #30 rows, yaaaaaah

# Get ASV table of the ASVs found above
namesDStop15 <- rownames(top15_DeSeq) #got ASV names
medEUNoEdgeASV <- t(ASVs_outta_ps(medianEUNoEdge.ps))    #pull out OG ASV table and flip it
top15ASVtab_DS <- medEUNoEdgeASV[namesDStop15,] #get a smaller version of the ASV table that has only these ASVs
#View(top15ASVtab_DS) #this has the ASV abundances (median of course) for each sample for these top 30 ASVs

# Now transform abundances in the above dataframe to z-scores
top15zScores_DS <- zScore(as.data.frame(top15ASVtab_DS))
#View(top15zScores_DS)

# Add in the information from the diffAbund analysis, namely logFoldChange, and some taxonomic info
top15_da_all <- merge(top15zScores_DS, top15_DeSeq[,c(2,8:11)], by= "row.names")
# View(top15_da_all) niiiiiiiice

top15_da_all_longer <- top15_da_all %>% 
  as_tibble() %>% 
  pivot_longer(cols= EU_10_10:EU_8_90,
               names_to="EUmeter",
               values_to="z-score")


# How many samples are each of these ASVs found in and what is their overall abundance?
#namesDSmedEU <- rownames(sigtab_medianEU) #pull out names of the ASVs found in DESeq analysis
#medEUNoEdgeASV <- t(ASVs_outta_ps(medianEUNoEdge.ps))    #pull out OG ASV table and flip it
#medEUNoEdgeASV_DS <- medEUNoEdgeASV[namesDSmedEU,] #get a smaller version of the ASV table that has only these ASVs
#medEUNoEdgeASV_DS_count <- rowSums(medEUNoEdgeASV_DS > 0) #get number of soil samples that the ASV appears in
#medEUNoEdgeASV_DS_counts <- cbind(medEUNoEdgeASV_DS, medEUNoEdgeASV_DS_count)
#ASVabund <- rowSums(medEUNoEdgeASV_DS_counts[,1:ncol(medEUNoEdgeASV)]) #get abundance of each ASV ACROSS all samples
#medEUNoEdgeASV_DS_counts <- cbind(medEUNoEdgeASV_DS_counts, ASVabund) #minimum is in 12/60 meter samples! max is in 54/60
# View(medEUNoEdgeASV_DS_counts) 

# ###########################################
# PLOTTING

medEU_DS_ps <- prune_taxa(namesDSmedEU, medianEUNoEdge.ps) #only has ASVs from DS analysis (678 ASVs)
medEU_DS_ASVtab <- ASVs_outta_ps(medEU_DS_ps)
medEU_DS_ASVtabTaxTab <- taxtable_outta_ps(medEU_DS_ps)
medEU_DS_taxTab <- as(tax_table(medEU_DS_ps), "matrix") #get this out of phyloseq

# New 

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
#EU52_topASVsMedLonger$Meter <- factor(EU52_topASVsMedLonger$Meter, levels = c("10", "20", "30", "40", "50","60", "70", "80", "90", "100"))

#HiMedLowASVs <- c(patchTop10$ASV_name, patchBottom10$ASV_name, patchMiddle10$ASV_name,
#                  forestTop10$ASV_name, forestBottom10$ASV_name, forestMiddle10$ASV_name)


#patchTop10.tb <- EU52_topASVsMedLonger %>% 
 # filter(ASV_name == patchTop10$ASV_name) %>% 
#  mutate(ASVgroup = "patchTop10")
#patchBottom10.tb <- EU52_topASVsMedLonger %>% 
#  filter(ASV_name == patchBottom10$ASV_name) %>% 
#  mutate(ASVgroup = "patchBottom10")
#patchMiddle10.tb <- EU52_topASVsMedLonger %>% 
#  filter(ASV_name == patchMiddle10$ASV_name) %>% 
#  mutate(ASVgroup = "patchMiddle10")
#forestTop10.tb <- EU52_topASVsMedLonger %>% 
#  filter(ASV_name == forestTop10$ASV_name) %>% 
#  mutate(ASVgroup = "forestTop10")
#forestBottom10.tb <- EU52_topASVsMedLonger %>% 
#  filter(ASV_name == forestBottom10$ASV_name) %>% 
#  mutate(ASVgroup = "forestBottom10")
#forestMiddle10.tb <- EU52_topASVsMedLonger %>% 
#  filter(ASV_name == forestMiddle10$ASV_name) %>% 
#  mutate(ASVgroup = "forestMiddle10")
