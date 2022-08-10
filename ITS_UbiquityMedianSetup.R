# ITS_UbiquityMedianSetup.R (Set-up/data wrangling before downstream)
# July 31, 2022 (begun)

#########
# This script picks up where ITSExploratoryDataAnalysis.R leaves off,
# applying a ubiquity filter of appearing in at least 40 samples
# to the samples (which have already been rarefied
# and where I have removed ASVs that don't occur at least 50 times in dataset).
# Next, it formats the OTU table in ITS_postUbiquity.ps for FUNGuild 
# (see https://github.com/UMNFuN/FUNGuild). 
# I also get the find the median abundance of each ASV in each EU (median over all four transects,
# at each spot along the transect).

###################################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")       
library("tidyverse")       
library("vegan")
library("ggplot2")

# Load data
load("RobjectsSaved/trimmedITSjustsoils.ps") #load phyloseq object that was made after rarefying,
# and keeping only ASVs that occurred at least 50 times across the (rarefied) 
# dataset (from ITS_ExploratoryDataAnalysis.R).

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

# 4. # taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

###############################################################################
#                   UBIQUITY FILTER
###############################################################################
# MOST OF NEXT 50 LINES BELOW COPIED FROM IDENTIFYEDGEEFFECTTAXA3.R
# In this section of the script, I take a look at the ASV ubiquity across samples
# and across EUs to determine what ubiquity filter to impose (I chose that the 
# ASVs need to occur in at least 45 samples). I use this to subset the ASV table.

# 1. What is the ubiquity of ASVs (looking at all EUs)?
ITSASVubiquity <- ASVsampOccurrence(trimmedITSjustsoils.ps)
# View(ITSASVubiquity) taxa/ASVs are rows, and samples are columns
dim(ITSASVubiquity)
abund <- rowSums(ITSASVubiquity)
ITSASVubiquityAbund <- as.data.frame(cbind(ITSASVubiquity, abund))
# View(ITSASVubiquityAbund) #this shows that 59 of these ASVs do not appear at least 50 times,
# should be due to removal of the controls
dim(ITSASVubiquityAbund[which(ITSASVubiquityAbund$sampleOccurence >= 45),]) #Only 352 ASVs are found at least 45 times
dim(ITSASVubiquityAbund[which(ITSASVubiquityAbund$sampleOccurence >= 40),]) #413 are found at least 40 times
dim(ITSASVubiquityAbund[which(ITSASVubiquityAbund$sampleOccurence >= 35),]) #602 are found at least 40 times

# Barplot showing ubiquity
#quartz()
barplot(table(ITSASVubiquity[,234]), ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 233)", main= "ASV Ubiquity in ITS Soils" )

# 2. How many EUs/sites does each ASV occur in?
pooledEU.ps <- merge_samples(trimmedITSjustsoils.ps, "EU") #merge samples by site
rownames(sample_data(pooledEU.ps)) #nice, this is by EU!
ITSASVubiquityEU <- as.data.frame(ASVsampOccurrence(t(pooledEU.ps))) #not sure why I had to transpose this...
dim(ITSASVubiquityEU) #dimensions are correct and is working after transposing
# View(ITSASVubiquityEU) This shows how many times each ASV appears in each EU.
# View(ITSASVubiquityEU[which(ITSASVubiquityEU$sampleOccurence ==6),]) # 692 ASVs appeared in all of the EUs!

# Barplot showing ubiquity by EU- shows that vast majority of ASVs, after dropping off 
# rare ones, are present in more than one site!
barplot(table(ITSASVubiquityEU[,7]), ylab="Number of ASVs", xlab= "Number of EUs/sites Present in", main= "ASV Ubiquity in ITS Soils By Number of EUs" )
# This shows that the majority of ASVs occur in all of the EUs

# Only use ASVs that occur in at least 40 samples
ITSnamesAll_40 <- names(which(ITSASVubiquity[,234] >= 40)) #this gives the names of the ASVs that occur at least
# 45 times (i.e. in at least 45 samples)
length(ITSnamesAll_40) #413, matches above
ITSnamesAll_40

# Remove all of the ASVs from the phyloseq object that don't occur at least 40 times
ITS_postUbiquity.ps <- prune_taxa(ITSnamesAll_40,trimmedITSjustsoils.ps) #352 taxa as expected!
ITS_postUbiquity.ps

# NMDS for post-ubiquity samples (Bray-Curtis dissimilarity)
ITS_postUbiq.ord <- ordinate(ITS_postUbiquity.ps, "NMDS", "bray")

###############################################################################
#         FORMATTING ITS_postUbiquity.ps FOR FUNGuild 
###############################################################################
# format is tsv, where first row should be OTU_IOTU_ID(tab)sample1(tab)sample2(tab)sample3(tab)taxonomy(return)
# taxonomic levels should be separated by semicolons

postUbiqITS_ASVtab <- t(ASVs_outta_ps(ITS_postUbiquity.ps)) #get ASV table out of phyloseq, ASVs are rows
postUbiqITS_TaxTab <- taxtable_outta_ps(ITS_postUbiquity.ps) #get tax table out of phyloseq; ASVs are also rows!
rownames(postUbiqITS_TaxTab) == rownames(postUbiqITS_ASVtab)
FUNGuildTable1 <- merge(postUbiqITS_ASVtab, postUbiqITS_TaxTab, by= "row.names", all=TRUE)
colnames(FUNGuildTable1)[1] <- "OTU_ID"
# Collapse taxonomy into one column 
FUNGuildTable1 <- tidyr::unite(FUNGuildTable1, sep=";",col= "taxonomy", Kingdom:Species)
# View(FUNGuildTable1)
# Remove all of the pesky stuff before taxonomy
FUNGuildTable1$taxonomy <-gsub("k__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("p__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("c__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("o__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("f__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("g__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("s__","",as.character(FUNGuildTable1$taxonomy))
# colnames(FUNGuildTable1)

# Written to file Aug. 10, 2022
write.table(FUNGuildTable1, file = "FUNGuildTable.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
# ** Finally, remove any quotation marks that may appear in new file in a text editor (I used Atom)**

###############################################################################
#         APPROACH 1: MEDIAN ASV ABUNDANCE BY METER PER EU
###############################################################################
# This code gets median abundance for each ASV (median across 4 transects in each EU) at each meter
# along the transect
ASVsdf <- ASVs_outta_ps(ITS_postUbiquity.ps) #rows are samples, ASV abundance is columns 
metaDf <- metadata_outta_ps(ITS_postUbiquity.ps) #rows are samples, columns are all of the rest of the 
# metadata
metaDf[,9]
# Make a new column that combines info from meter and EU
EUmeter <- rep(NA, nrow(metaDf)) #pre-allocate vector for new names,
for (i in 1:nrow(metaDf)){
  EUmeter[i] <- paste(c(metaDf[i,6], metaDf[i,9]), collapse="_") #get EU and meter from sample data
} #this should work, see example below, but is taking forever!
length(EUmeter) #233, same as number of samples

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

nrow(medAbund) == ncol(ASVsdf)*60 #21120 is number of ASVs times number of new categories (i.e. 10 meters in transect * 6 EUs)

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
dim(abundCheck) #dim is 82016 rows and 6 columns 
# View(abundCheck)
# View(medAbund)

# Manually check one to make sure that median function is working properly:
# ASV 1 across EU10, 10m (only have T, B, and R)
EU_10_T10 <- abundCheck[1,6]
EU_10_R10 <- abundCheck[2,6]
EU_10_B10 <- abundCheck[3,6]
median(as.numeric(c(EU_10_L10, EU_10_T10, EU_10_R10, EU_10_B10))) == medAbund[1,3] #TRUE

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
ASVtabMedian <- medAbund[,1:3] %>% #remove n column 
  ungroup %>% # necessary to remove columns
  pivot_wider(names_from = ASVname, values_from= medianEUmeter) %>% 
  column_to_rownames("EUmeter") %>% 
  t()
# View(ASVtabMedian) #correct format and values!

# 2. Taxonomy table -- ***GET THIS OUT OF PHYLOSEQ!!***
taxTabMedian <- tax_table(ITS_postUbiquity.ps) #since we did not drop any ASVs when coding above,
# this tax table should be the same as before 
unique(sort(rownames(taxTabMedian)) == sort(rownames(ASVtabMedian))) #Statement above is true!!!
taxTabMedian <- as(taxTabMedian, "matrix") #get this out of phyloseq
class(taxTabMedian)

# 3. Sample metadata
postUbiquityMeta <- metadata_outta_ps(ITS_postUbiquity.ps)
head(postUbiquityMeta)
# Constructing this for now; can find a more elegant coding solution later
EUmeterNames <- colnames(ASVtabMedian)
EUmeterNames

metaDatMedian <- as.data.frame(matrix(nrow = length(EUmeterNames), ncol=3))
rownames(metaDatMedian) <- EUmeterNames
colnames(metaDatMedian) <- c("EU", "meter", "Habitat")
metaDatMedian[,1] <- c(rep("EU_10", 10), rep("EU_52", 10), rep("EU_53N", 10), rep("EU_53S", 10), rep("EU_54S", 10), rep("EU_8", 10))
metaDatMedian[,2] <- rep(c(10, 100, 20, 30, 40, 50, 60, 70, 80, 90),6)
metaDatMedian[,3] <- rep(c("patch", "forest", "patch", "patch", "patch", "edge", "forest", "forest", "forest", "forest"), 6)

# Transform all of the above to phyloseq objects
ASVtabMedianPS <- otu_table(ASVtabMedian, taxa_are_rows = TRUE)
taxTabMedianPS <- tax_table(taxTabMedian)
metaDatMedianPS <- sample_data(metaDatMedian)
ITS_medianEU.ps <- phyloseq(ASVtabMedianPS, taxTabMedianPS, metaDatMedianPS)
ITS_medianEU.ps #as expected, this has 4480 taxa and 60 samples!
colnames(otu_table(ITS_medianEU.ps))

#####################################################################
## FINALLY, SAVE THESE NEW PHYLOSEQ OBJECTS FOR EASY ACCESS #
# (saves are commented out at the end of this script because they do not 
# necessarily need to be re-saved everytime this script is run)

# 1. Phyloseq object after rarefying, removal of rares (rares= <50 reads across dataset 
# AND THEN ubiquity threshold removed < 40 samples across that dataset)
#save(ITS_postUbiquity.ps, file= "RobjectsSaved/ITS_postUbiquity.ps") #saved July 31, 2022
#save(ITS_medianEU.ps, file = "RobjectsSaved/ITS_medianEU.ps")


###################################################################
# December 14, 2021
# Updating taxonomic barplots!!

#####################
# TOP PHYLA
#####################

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)

ITS_medianEUNoEdge.ps <- subset_samples(ITS_medianEU.ps, Habitat != "edge")

sample_data(ITS_medianEUNoEdge.ps)$EUHabitat <- rep(NA, nrow(sample_data(ITS_medianEUNoEdge.ps))) #pre-allocate new category for habitat + EU
for (i in 1:nrow(sample_data(ITS_medianEUNoEdge.ps))) { #now make it!
  sample_data(ITS_medianEUNoEdge.ps)$EUHabitat[i] <- paste(c(sample_data(ITS_medianEUNoEdge.ps)$EU[i], sample_data(ITS_medianEUNoEdge.ps)$Habitat[i]), collapse="")
}

ITS_medianEU.ps.phylum.glom <-  tax_glom(ITS_medianEUNoEdge.ps, taxrank = "Phylum") 
tax_table(ITS_medianEU.ps.phylum.glom) #8 phyla
sample_data(ITS_medianEU.ps.phylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabunMed.phyla.0 <- transform_sample_counts(ITS_medianEU.ps.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(ITS_medianEU.ps.phylum.glom)) 

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunMed.phyla.1 <- merge_samples(relabunMed.phyla.0, group = c("EUHabitat"))
sample_data(relabunMed.phyla.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)
# Meter did an averaging thing; can just ignore it

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunMed.phyla.2 <- transform_sample_counts(relabunMed.phyla.1, function(x) x / sum(x))
sample_data(relabunMed.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabun.phyla.df <-psmelt(relabunMed.phyla.2)
dim(relabun.phyla.df) #
sum(relabun.phyla.df[,3]) 
colnames(relabun.phyla.df) 
relabun.phylatop99 <- relabun.phyla.df
relabun.phylatop99.5 <- relabun.phyla.df

relabun.phylatop99$Phylum[relabun.phylatop99$Abundance < 0.01] <- "< 1% abund."
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abund."

top_99p_phyla <- unique(relabun.phylatop99$Phylum)
top_99p_phyla

top_99.5p_phyla <- unique(relabun.phylatop99.5$Phylum)
top_99.5p_phyla

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Fungal phyla comprising at least 0.5% of total abundance (median and post-ubiquity)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# "Phyla comprising at least 1% of total abundance"
phylumPlot.99percent <- ggplot(data=relabun.phylatop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 8, face = "bold"))
quartz()
phylumPlot.99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Fungal phyla comprising at least 1% of total abundance (median and post-ubiquity)")

# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

top_99.5p_phyla <- relabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

###############################################################################
#         APPROACH 2: MEDIAN ASV ABUNDANCE BY METER ACROSS ALL EUS!
# THIS BELOW CAN PROBABLY BE DROPPED, AS IT WAS NOT USED DOWNSTREAM
###############################################################################
# This code gets median abundance for each ASV (median across 24 transects all EUs)
# at each meter.
# In other words, the EU is ignored, so the number of rows = number of ASVs in 
# ITS_postUbiquity.ps (4480) * points along transect (10). =

#ASVsdf <- ASVs_outta_ps(ITS_postUbiquity.ps) #rows are samples, ASV abundance is columns 
#metaDf <- metadata_outta_ps(ITS_postUbiquity.ps) #rows are samples, columns are all of the rest of the 
# metadata

#  I DONT THINK THAT THIS CALCULATED THE MEDIAN ABUDNANCES, JUST TOOK OTHER VARIABLE MEDIANS (SUM OF OTU TABLE, SEE MAN PAGE)
#medianMeter_noEUs.ps <- merge_samples(ITS_postUbiquity.ps, "Meter", fun= median)
#str(otu_table(medianMeter_noEUs.ps))
# Need to check if it actually calculated the median
#as.data.frame(otu_table(medianMeter_noEUs.ps))$ASV_1

# Make a new column that combines info from meter and EU
#EUmeter <- rep(NA, nrow(metaDf)) #pre-allocate vector for new names, length is 268,800
#for (i in 1:nrow(metaDf)){
 # EUmeter[i] <- paste(c(metaDf[i,6], metaDf[i,9]), collapse="_")
#} #this should work, see example below, but is taking forever!
#EUmeter

# NEED TO DOUBLE CHECK THIS MANUALLY
#medAbundnoEU <- ASVsdf %>% #First, add columns of interest
 # mutate( #add columns
 #   EU= metaDf$EU,#add EU column
 #   Transect = metaDf$Transect, #add Transect column
  #  Meter = metaDf$Meter,# meter column
  #  EUmeter = EUmeter) %>% #EU and meter information
# pivot_longer(
#    cols = c(1:ncol(ASVsdf)),
 #   names_to = "ASVname"
#  ) %>% 
#  dplyr::group_by(Meter, ASVname) %>% 
 # dplyr::summarize(
 #   medianMeter = median(value), #get median abundance for each ASV, with samples grouped by EU and meter (i.e. EU Meter)
 #   n = n())  #n= how many samples at meter in that EU 
# View(medAbundnoEU)

