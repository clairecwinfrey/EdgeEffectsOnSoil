# ITS ONLY Exploratory Data Analysis and Data Clean Up
# started May 13, 2022

# This script is to take a first look at ITS MiSeq data from the samples that I
# took in May 2021 to examine an edge effect across forest and savanna boun-
# daries at the Savannah River Site. Here, I look at sequence counts, rarefy,
# examine top classes and phyla, and look at ordinations. I repeat the ordinations
# and examination of top classes and phyla before and after removing rare species
# (here those that showed up LESS THAN 50 times across whole dataset).
# Then, I move forward with the rarefied, rare-removed data set and do preliminary
# investigations into why outlier samples might be outliers, how savannas/patches
# differ from forests, etc.

############# CLEAN UP!! #############

# FILES SAVED IN THIS SCRIPT (both in RObjectsSaved directory: (save code at end of script,
# but commented out so that new files are not accidentally re-saved/written over):
# 1. "EDA_ITSJuly2022"
# ITSrarefied.ps, samples_df, ASVsTrimmed, taxTrimmed, trimmedITSjustsoils.ps, outliersASVtax

# 2. "trimmedITSjustsoils.ps" 
# This is the phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset, including
# soils and controls (NOT biocrusts). No outliers were removed. 

# FUNCTIONS CREATED IN THIS SCRIPT:
# 1. Function to get the taxonomy table out of phyloseq:
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. Function to get ASV table out of phyloseq so that we can view it better
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

# 3. Function that gets OTU/ASV table out of phyloseq
# (from https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# 4. Function that gets average pH for each sample, based on the 2-3 pH readings taken
mean_pH_func <- function(pHs_mat) { #enter a dataframe that has 1-3 pHs per row
  # this function works for averaging 2-3 pHs (per row), but could easily be adjusted in
  # for loop to be more flexible.
  H_mat <- matrix(nrow=nrow(pHs_mat), ncol=ncol(pHs_mat))
  colnames(H_mat) <- colnames(pHs_mat)
  for(i in 1:nrow(pHs_mat)){ #get the H+ concentration for each pH reading
    H_mat[i,1] <- 10^(pHs_mat[i,1]*-1)
    H_mat[i,2] <- 10^(pHs_mat[i,2]*-1)
    H_mat[i,3] <- 10^(pHs_mat[i,3]*-1)
  }
  avgH <- rowMeans(H_mat, na.rm=TRUE) #get the average H+ concentration for each set of pHs
  avgpH <- -log10(avgH) #convert back to pH
  return(avgpH)
}


#######################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("ggplot2")      #graphics
library("readxl")       #necessary to import the data from Excel file
library("dplyr")        #filter and reformat data frames
library("tibble")       #Needed for converting column to row names
library("tidyverse")
library("mctoolsr")
library("vegan")
library("gridExtra")    #allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis

##################################################################################
# I. SET-UP, DATA CLEANING, RAREFACTION, AND FIRST TAXONOMIC & ORDINATION PLOTS
##################################################################################

seqtab_wTax_mctoolsr <- read.table("Raw_Data/ITS/03_tabletax/seqtab_wTax_mctoolsr.txt")
str(seqtab_wTax_mctoolsr)
head(seqtab_wTax_mctoolsr)
#View(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr)
rownames(seqtab_wTax_mctoolsr)

seqtab <- read.table("Raw_Data/ITS/03_tabletax/seqtab_final.txt", header=T)
dim(seqtab) #11,342 ASVs and 264 samples (soil samples, controls, etc.)
#View(seqtab)

tax_final.txt <- read.table("Raw_Data/ITS/03_tabletax/tax_final.txt")
str(tax_final.txt)
#View(tax_final.txt)

#### LOAD FILES IN FOR PHYLOSEQ #####
# 1. OTU table: ASVs/OTUs as rows, samples as columns. This is "seqtab"
otu_mat <- seqtab
seqtab$X #X is the ASV names
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X") #make it so that "X"(which is column with ASV names is the rownames columnn!)
colnames(otu_mat)
rownames(otu_mat)
#View(otu_mat)

######### CAN WORK ON THIS STEP MORE TO REMOVE phylum, class, order things before each taxon name ############
# 2. OTU taxonomy
# This is basically needs to be the same as the first and last columns of seqtab_wTax_mctoolsr
# So none of the guts that are essentially the ASV table
dim(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr) #samples
head(seqtab_wTax_mctoolsr)
head(seqtab_wTax_mctoolsr$V1) #THE ASV names!
head(seqtab_wTax_mctoolsr$V265) #the last and second to last columns are the taxonomy (i.e. V265 and V266)
# The last column (266) sonetimes has species info
head(seqtab_wTax_mctoolsr$V266) #species information only
step1 <- cbind(seqtab_wTax_mctoolsr$V1, seqtab_wTax_mctoolsr$V265) #bind together ASV IS names and taxonomy  
head(step1)
dim(step1)
step2 <- cbind(step1, seqtab_wTax_mctoolsr$V266)
head(step2) # columns are now "#ASV_ID", "taxonomy", and "Species"
colnames(step2) <- step2[1,] #make ASV_ID and taxonomy the column names
head(step2)
dim(step2)
step2 <- step2[-1,] #get rid of first row which is same as column names
head(step2)

# Make new columns based on the ";" separator
require(tidyr)
tax_sep <- separate(as.data.frame(step2), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                  "Order", "Family", "Genus", "ASV_ID"),
                    sep = ";")
head(tax_sep)
# Continuing tidying up taxonomy table
# View(tax_sep)
str(tax_sep)
all.equal(tax_sep$`#ASV_ID`, tax_sep$ASV_ID) #this is true, so can remove ASV_ID 
colnames(tax_sep)
tax_sep <- tax_sep[,-8] #remove duplicate ASV ID column
#View(tax_sep)
rownames(tax_sep)
tax_mat <- tax_sep %>% #make ASV_ID column the rownames to make it compatible with phyloseq
  tibble::column_to_rownames("#ASV_ID")
head(tax_mat)

# 3. Sample metadata-- this should be different for fungi than for 16S!!!!
### SRS_ITS_AllMetadataAug17.csv was double checked July 30, 2022 to make sure it was correct
metadata <- read.csv("Bioinformatics/SRS_ITS_AllMetadataAug17.csv") #this new csv has ALL sample metadata, but is missing some mean pH values
# not just for biological samples (i.e. controls too)
# View(metadata)
colnames(metadata)
rownames(metadata)
metadata$X <- NULL #dunno whatX is but it's not important 
samples_df <- metadata %>%
  tibble::column_to_rownames("ITS_MappingSampID")
str(samples_df)
rownames(samples_df)
colnames(samples_df)
samples_df$mean_pH <- mean_pH_func(samples_df[,3:5]) #get mean pH for each sample
# View(samples_df) #several of these were hand checked to make sure that function was working correctly

# Need to rename samples_df to match those in otu_mat; because they don't match
# each other except for the blanks and controls, those are the only ones that got
# phyloseq'd
#Remove Xs from the colnames on otu_mat
colnames(otu_mat)
newnames <- str_remove(colnames(otu_mat), "X")
colnames(otu_mat) <- newnames
length(rownames(samples_df)) == length(colnames(otu_mat)) #FALSE
# Where do they differ?
diffBWsamplesdf_otu_mat <- setdiff(rownames(samples_df), colnames(otu_mat)) 
diffBWsamplesdf_otu_mat #samplesdf has "ExtControlWater_3"  "ExtControlWater_5"  "ExtControlWater_7" 
# "ExtControlWater_8"  "ExtControlWater_16" "ExtControlWater_18", whereas otu_mat does not.
# THESE SAMPLES FAILED IN THE RUN, SO DROP THEM FROM SAMPLES_DF
Overlapsamplesdf_otu_mat <- intersect(rownames(samples_df), colnames(otu_mat)) #get the sample names in common
samples_df <- samples_df[Overlapsamplesdf_otu_mat,]
all.equal(sort(colnames(otu_mat)), sort(rownames(samples_df))) #Cool, names are now the same

# Transform otu table and taxonomy tables into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#View(otu_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
SRS_ITS_raw.ps <- phyloseq(OTU, TAX, samples)
SRS_ITS_raw.ps 
colnames(otu_table(SRS_ITS_raw.ps))

sample_names(SRS_ITS_raw.ps) # IMPORTANT: THESE ARE THE NAMES BASED ON THE MAPPING
# FILE USED FOR DEMULITPLEXING; I.E. THEY CORRESPOND TO ORDER LOADED INTO ROWS
# OF PLATES (*NOT NOT NOT*) DNA TUBE LABELS 

rank_names(SRS_ITS_raw.ps)
sample_variables(SRS_ITS_raw.ps)

##################################################################
# Check for anything non-fungal and ASVs non classified to at least phylum level
# and remove if present-- 
##################################################################

# Get taxonomy table out of phyloseq:
ITS_taxTable <- taxtable_outta_ps(SRS_ITS_raw.ps)
# View(ITS_taxTable) #raw has 11,342 ASVs

# Remove EVEERYTHING THAT ISN'T KINGDOM FUNGI:
tax_ITSnoNAs <- ITS_taxTable %>% 
  filter(Kingdom != "NA") %>% #Remove ASVs where kingdom is unknown ("NA" in first column) 
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in first column) 
# View(tax_ITSnoNAs) 
# Are any of these remaining a phylum other than fungi?
kingdomFung <- tax_ITSnoNAs %>% 
  filter(Kingdom == "k__Fungi")
dim(kingdomFung)[1] == dim(tax_ITSnoNAs)[1] #no, there were no kingdoms other than fungi

# We filtered out a total of 2854 taxa, for a resulting 8,488 ASVs.
dim(ITS_taxTable)[1] - dim(tax_ITSnoNAs)[1]

# We now have to trim the OTU table, because it still has the taxa in it that were filtered out
# above. It currently has 11,342 ASVs, whereas tax_ITSnoNAs has 8,488 ASVs. We'll want it this 
# way so that we re-make the phyloseq object and then rarefy correctly. 
dim(otu_table(SRS_ITS_raw.ps)) #11342 ASVs

ITS_ASVsToKeep <- rownames(tax_ITSnoNAs) #grab names of ASVs to keep (that were filtered out of tax
# table above)

ITS_noNAs.ps <- prune_taxa(ITS_ASVsToKeep, SRS_ITS_raw.ps) #new phyloseq object with just the stuff we want!
# Check to make sure new phyloseq object looks as expected
unique(rownames(tax_ITSnoNAs) == rownames(otu_table(ITS_noNAs.ps))) #yep! 

# TAKING A LOOK AT THE DATA
# How are number of reads distributed across the samples?
colnames(otu_table(ITS_noNAs.ps))
ITSseqspersample <- colSums(otu_table(ITS_noNAs.ps))
ITSSeqNumberPlot <-barplot(seqspersample, main="All ITS: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number")

# Just controls:
ITSseqsperctrl <- ITSseqspersample[244:length(ITSseqspersample)]
ITSSeqNumberPlotCtrls <- barplot(ITSseqsperctrl, main="ITS Controls: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.27)
ITSSeqNumberPlotCtrls #ExtBlank_3, ExtControlWater_1, ExtControlWater_2, ExtControlWater_22, and ExtControlWater_20 were pretty high

# Controls with easier axis names
simple <- c("ExtB1", "ExtB2", "ExtB3", "ExtW1", "ExtW10", "ExtW12", "ExtW13",
            "ExtW14", "ExtW15", "ExtW17", "ExtW19", "ExtW2", "ExtW20", "ExtW21",
            "ExtW22","ExtW4","ExtW9","PCRNTC1", "PCRNTC2", "PCRNTC3")
seqsperctrl_simple <- setNames(ITSseqsperctrl, simple)

#quartz()
barplot(seqsperctrl_simple, main="ITS Controls: Total Sequences Per Sample",
        ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.35)

######## RAREFACTION ##########
# What is the minimum for the non-control samples?
min(ITSseqspersample[1:243]) #0 - extraction blank 17! But all the PCR NTCs are 6 or lower too
max(ITSseqspersample[1:243]) #108353
mean(ITSseqspersample[1:243]) # 31289.93

# All samples?
min(ITSseqspersample) #0
max(ITSseqspersample) #108353
mean(ITSseqspersample) # 29079.19
sd(ITSseqspersample) #17177.62

# and the max and mins for the controls?
min(ITSseqspersample[244:length(ITSseqspersample)]) #0
max(ITSseqspersample[244:length(ITSseqspersample)]) #13171
mean(ITSseqspersample[244:length(ITSseqspersample)]) #2218.7
sd(ITSseqspersample[244:length(ITSseqspersample)]) #4071.68

# Plot this:
#quartz()
seqsPerExSamp <- ITSseqspersample[1:243]
ITSSeqNumberPlotExSamp <- barplot(sort(seqsPerExSamp), main="ITS Soils: Total Sequences Per Sample",
                               ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.4)

sort(seqsPerExSamp) #take a look in order from least to greatest
sort(seqspersample) #take a look in order from least to greatest (all samples and controls)

# We'll rarefy at 9717, which only costs us 6 samples (plus all but two of the controls)
# The "random rarefaction is made without replacement so that the variance of rarefied
# communities is rather related to rarefaction proportion than to the size of the sample". 
# (quoted bit from vegan's manual page)

set.seed(19)
ITSITSrarefied.ps <- rarefy_even_depth(ITS_noNAs.ps, sample.size = 9717, replace=FALSE, trimOTUs=TRUE)
# 258O ASVs were removed because they are no longer present in any sample after random subsampling

# Rarefaction curve:
samp.col <- c(rep("blue", 243), rep("grey", 20))

## Rarefaction curves look good! 
quartz()
rare.plot <- rarecurve(t(otu_table(ITS_noNAs.ps)), step = 3000, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")
rarecurve(t(otu_table(noeuksorNAs_ps)), step = 3000, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")

# How many are left?
colSums(otu_table(ITSITSrarefied.ps)) #All have 9,717
length(colSums(otu_table(ITSITSrarefied.ps))) #238 samples left
# Which samples are these:
sampleNames <- names(colSums(otu_table(ITSITSrarefied.ps)))
sampleNames

####################################################################
# TOP PHYLA AND CLASSES AND ORDINATIONS (BEFORE APPLYING RARITY OR UBIQUITY THRESHOLDS)
# The following lines is all done BEFORE removing rare taxa,
# (i.e., with rarefied dataset, but BEFORE removing those that do not occur at least 50 
# times in the dataset))

#####################
# TOP PHYLA
#####################

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
# TURN ASVs INTO PHYLUM LEVEL
ITSrarefied.phylum.glom <-  tax_glom(ITSrarefied.ps, taxrank = "Phylum") 
tax_table(ITSrarefied.phylum.glom) # 13 phyla

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
ITSrelabun.phyla.0 <- transform_sample_counts(ITSrarefied.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(ITSrelabun.phyla.0)) #ASVs are just a representative from each phylum

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
ITSrelabun.phyla.1 <- merge_samples(ITSrelabun.phyla.0, group = "EU")
sample_data(ITSrelabun.phyla.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
ITSrelabun.phyla.2 <- transform_sample_counts(ITSrelabun.phyla.1, function(x) x / sum(x))
sample_data(ITSrelabun.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
ITSrelabun.phyla.df <-psmelt(ITSrelabun.phyla.2)
dim(ITSrelabun.phyla.df) #117 ASVs by 13 metadata variables 
sum(ITSrelabun.phyla.df[,3])  #9
colnames(ITSrelabun.phyla.df) #metadata variables 
ITSrelabun.phylatop99 <- ITSrelabun.phyla.df
ITSrelabun.phylatop99.5 <- ITSrelabun.phyla.df

ITSrelabun.phylatop99$Phylum[ITSrelabun.phylatop99$Abundance < 0.01] <- "< 1% abund."
ITSrelabun.phylatop99.5$Phylum[ITSrelabun.phylatop99.5$Abundance < 0.005] <- "< .5% abund."

ITStop_99p_phyla <- unique(ITSrelabun.phylatop99$Phylum)
ITStop_99p_phyla #"p__Basidiomycota", p__Ascomycota", "p__Mortierellomycota" "p__Mucoromycota", "p__Glomeromycota", "< 1% abund."

ITStop_99.5p_phyla <- unique(ITSrelabun.phylatop99.5$Phylum)
ITStop_99.5p_phyla #"p__Basidiomycota", p__Ascomycota", "p__Mortierellomycota" "p__Mucoromycota", "p__Glomeromycota", "< 1% abund."

# Phyla comprising at least 0.5% of total abundance
ITSphylumPlot99.5percent <- ggplot(data=ITSrelabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
ITSphylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 0.5% of total abundance (before removing rare or non-ubiq.)")

# "Phyla comprising at least 1% of total abundance"
ITSphylumPlot.99percent <- ggplot(data=ITSrelabun.phylatop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
ITSphylumPlot.99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 1% of total abundance (before removing rare or non-ubiq.)")

# Get exact abundances of each phyla (top 99.5%):
colnames(ITSrelabun.phylatop99.5)
ITSrelabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

ITStop_99.5p_phyla <- ITSrelabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) %>% 
  View()

#####################
# TOP CLASSES:
#####################

# (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
# TURN ASVs INTO class LEVEL
ITSrarefied.class.glom <-  tax_glom(ITSrarefied.ps, taxrank = "Class") 
tax_table(ITSrarefied.class.glom) # good, this is only class (59 different phyla!)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
ITSrelabun.class.0 <- transform_sample_counts(ITSrarefied.class.glom, function(x) x / sum(x) )

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
ITSrelabun.class.1 <- merge_samples(ITSrelabun.class.0, group = "EU")
sample_data(ITSrelabun.class.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
ITSrelabun.class.2 <- transform_sample_counts(ITSrelabun.class.1, function(x) x / sum(x))
sample_data(ITSrelabun.class.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
ITSrelabun.class.df <-psmelt(ITSrelabun.class.2)
dim(ITSrelabun.class.df) #
sum(ITSrelabun.class.df[,3]) 
colnames(ITSrelabun.class.df) 
ITSrelabun.classtop99 <- ITSrelabun.class.df
ITSrelabun.classtop95 <- ITSrelabun.class.df

ITSrelabun.classtop99$Class[ITSrelabun.classtop99$Abundance < 0.01] <- "< 1% abund."
ITSrelabun.classtop95$Class[ITSrelabun.classtop95$Abundance < 0.05] <- "< 5% abund."

ITStop_99p_class <- unique(ITSrelabun.classtop99$Class)
ITStop_99p_class

ITStop_95p_class <- unique(ITSrelabun.classtop95$Class)
ITStop_95p_class

# Phyla comprising at least 5% of total abundance
ITSclassPlot.95pt <- ggplot(data=ITSrelabun.classtop95, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
ITSclassPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 5.5))  + ggtitle("Classes comprising at least 5% of total abundance (before removing rare or non-ubiq.)")

###########
# SOME ORDINATIONS

# Using a phyloseq tool for convenience!
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ITSord <- ordinate(ITSrarefied.ps, method = "NMDS", distance = "bray", trymax = 100)

# With no black outline around points
ITSrarefiedBrayNMDS <- phyloseq::plot_ordination(ITSrarefied.ps, ITSord, type= "samples", color= "EU")
quartz()
ITSrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on ITS Bray-Curtis Dissimilarities (before removing rare or non-ubiq.)")

# With black outline around points:
ITSrarefiedBrayNMDSoutlined <-ITSrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(aes(fill=EU),color="black",pch=21, size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities (before removing rare or non-ubiq.)")
quartz()
ITSrarefiedBrayNMDSoutlined

# Now, remove biocrust and controls from ordination:
ITSjustsoils.ps <- subset_samples(ITSrarefied.ps, Type != "BioCrust" & Type != "ExtContWater" & Type !="ExtBlank")
unique(sample_data(ITSjustsoils.ps)[,27]) #only soils
sample_names(ITSjustsoils.ps) #another check to show that we have only soils!

# Ordination plot of just soils:
set.seed(19)
ITSordSoils <- ordinate(ITSjustsoils.ps, method = "NMDS", distance = "bray", trymax = 100)
ITSsoilsrarefiedBrayNMDS <- phyloseq::plot_ordination(ITSjustsoils.ps, ITSordSoils, type= "samples", color= "EU")
quartz()
ITSsoilsrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(aes(fill=EU),color="black", pch=21, size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities (before removing rare or non-ubiq.)")

# Add in labels to figure out what weird samples are
ITSsoilsrarefiedBrayNMDS <- phyloseq::plot_ordination(ITSjustsoils.ps, ITSordSoils, type= "samples", color= "EU", label = "Sample.ID")
quartz()
ITSsoilsrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(aes(fill=EU),color="black", pch=21, size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities (before removing rare or non-ubiq.)")

# Visible outliers: (going clockwise from the top left on the ordination plot above):
outliers_ITS <- c("52D_R_10", "53ND_R_40", "53ND_B_20", "53ND_T_20", "52D_L_100")
outliers_16S <- c("53ND_B_20", " 10C_L_10", "53ND_R_40", "53SD_R_10", "52D_L_100", "52D_R_90", "53ND_L_80") #these are from 16SExploratoryDataAnalysisAug2021.R
intersect(outliers_ITS, outliers_16S) #about half of the outlier samples for 16S are also outliers in ITS. "53ND_R_40" "53ND_B_20" "52D_L_100"

##################################################################################
# II. REMOVAL OF "RARE" TAXA, NEW TAXONOMIC PLOTS AND ORDINATIONS
# (tried out different numbers of taxa, but settled on removing taxa that did not
# occur at least 50 times across the whole dataset to match 16S analysis
##################################################################################

###########################################################################
# REMOVE "RARE" TAXA AND THEN RE-RUN 1) sequences per sample, 
# 2) top phyla, 3) top classes, 4) ordinations
# Questions when comparing pre and post removal of rare taxa:
# 1) Do top phyla or classes change, or outliers, after this removal? If so, this is
# evidence of these shifts being driven by rare taxa:

ITSrarefiedNoBC.ps <- subset_samples(ITSrarefied.ps, Type != "BioCrust") #make a phyloseq object that does not have biocrusts, but
# that retains controls. 

rare_ASVtab <- ASVs_outta_ps(ITSrarefiedNoBC.ps) 
# Add column for abundance of ASVs across all (pre-rarefied) samples
rare_ASVtab$Abundance <- rowSums(rare_ASVtab)
dim(rare_ASVtab) #8,230 ASVs
# View(rare_ASVtab)

# Most are rare! 
quartz()
plot(rare_ASVtab$Abundance)
length(which(rare_ASVtab$Abundance <= 50)) #5,173 ASVs have 50 or fewer occurrences across the rarefied data set 
length(which(rare_ASVtab$Abundance <= 35)) #4,641 ASVs have 30 or fewer occurrences across the rarefied data set
length(which(rare_ASVtab$Abundance <= 10)) #2,639 ASVs have 10 or fewer occurrences across the rarefied data set

length(which(rare_ASVtab$Abundance >= 50)) #3,095 ASVs have 50 or more occurrences across the rarefied data set 

# We'll get rid of the ASVs that do not occur AT LEAST 50 times across our dataset:
keptASVsindex <- which(rare_ASVtab$Abundance >= 50) #gives row numbers to keep in the dataset 

ASVsTrimmed <- rare_ASVtab[keptASVsindex,] #keep only rows of the index, and all columns (i.e. all samples)
dim(ASVsTrimmed) #3095 ASVs across 236 samples, as expected 
quartz()
plot(ASVsTrimmed$Abundance, main= "Abundance of Each ASV (after removal of those <50 times)", xlab="ASV number", ylab= "ASV abundance") # tail is still long, because most ASVs are still rare
# relative to the most abundant ASVs

length(which(ASVsTrimmed$Abundance == 50)) # 38 ASVs were right at cut off!
range(ASVsTrimmed$Abundance) #50 65610
length(which(ASVsTrimmed$Abundance >= 1000)) #only 392 ASVs appear more than 1000 times
# across the data set, which is in average of about 4 times per sample!

# Get TAXONOMIC TABLE that matches the ASV table above using function we defined earlier
rare_taxTab <- taxtable_outta_ps(ITSrarefied.ps) 
# (first check that the ASVs (i.e. rownames) are in the same order for the two data frames so that we can use the 
# same index)
unique(rownames(rare_taxTab) == rownames(rare_ASVtab)) # all are the same, son we're good to move to the next step 
taxTrimmed <- rare_taxTab[keptASVsindex,] #keep only rows of the index made above for tax appearing at least 50 times
unique(rownames(taxTrimmed) == rownames(ASVsTrimmed)) #ALl true so this is the same!

# How many ASVs occur in every sample?
# Said another way: how many rows (i.e. ASVs) occur >= 1 time across all columns (i.e. samples)
# Or: are there any rows that contain no zeros?
# Inspiration from this post: https://stackoverflow.com/questions/9977686/how-to-remove-rows-with-any-zero-value

# Let's you find rows that contain zeros and then filter them out
rownames(ASVsTrimmed[!(apply(ASVsTrimmed, 1, function(y) any(y == 0))),])
# no ASVs are found in every sample. Surprising! But only 5 were found in every sample for 16S

# Merge trimmed ASV tables and taxonomy tables to get something that is formatted like
# "seqtab_wTax_mctoolsr"
colnames(ASVsTrimmed)
seqtab_wTax_trimmed <- cbind.data.frame(ASVsTrimmed[,-236], taxTrimmed) #-236 to remove "abundance" column
head(seqtab_wTax_trimmed)
View(seqtab_wTax_trimmed)

# Make Excel files of this!
# Saved to CSV July 31, 2022
write.csv(seqtab_wTax_trimmed, "SRSMay2021_ITS_seqtab_wtax_TRIMMED.csv")

############ RE-DO ANALYSES WITH DATASET WHERE THE RARE TAXA HAVE BEEN REMOVED (counting soils and controls
#### NEED TO REMAKE THESE INTO PHYLOSEQ OBJECTS FIRST (READ THROUGH THIS TOO)

class(taxTrimmed)
dim(ASVsTrimmed)
#View(ASVsTrimmed)

otu_mat <- as.matrix(ASVsTrimmed[,-236]) #remove "Abundance" column that is at the end
tax_mat <- as.matrix(taxTrimmed)
#View(otu_mat)
#View(tax_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
ITS_TrimmedSRS_16S.ps <- phyloseq(OTU, TAX, samples)
ITS_TrimmedSRS_16S.ps
dim(otu_table(ITS_TrimmedSRS_16S.ps))

#####################
# trimmedTop PHYLA
#####################

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_bartrimmedPlots.R)
# TURN ASVs INTO PHYLUM LEVEL
rarefiedTrimmed.phylum.glom <-  tax_glom(ITS_TrimmedSRS_16S.ps, taxrank = "Phylum") 
tax_table(rarefiedTrimmed.phylum.glom) # 10 different phyla!

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabunTrimmed.phyla.0 <- transform_sample_counts(rarefiedTrimmed.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabunTrimmed.phyla.0)) #same ASVs as above

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunTrimmed.phyla.1 <- merge_samples(relabunTrimmed.phyla.0, group = "EU")
sample_data(relabunTrimmed.phyla.1)  #shows by EU worked

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunTrimmed.phyla.2 <- transform_sample_counts(relabunTrimmed.phyla.1, function(x) x / sum(x))
sample_data(relabunTrimmed.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabunTrimmed.phyla.df <-psmelt(relabunTrimmed.phyla.2)
dim(relabunTrimmed.phyla.df) #
sum(relabunTrimmed.phyla.df[,3]) 
colnames(relabunTrimmed.phyla.df) 
relabunTrimmed.phylatrimmedTop99 <- relabunTrimmed.phyla.df
relabunTrimmed.phylatrimmedTop99.5 <- relabunTrimmed.phyla.df
relabunTrimmed.phylatrimmedTop95 <- relabunTrimmed.phyla.df

relabunTrimmed.phylatrimmedTop99$Phylum[relabunTrimmed.phylatrimmedTop99$Abundance < 0.01] <- "< 1% abund."
relabunTrimmed.phylatrimmedTop99.5$Phylum[relabunTrimmed.phylatrimmedTop99.5$Abundance < 0.005] <- "< .5% abund."
relabunTrimmed.phylatrimmedTop95$Phylum[relabunTrimmed.phylatrimmedTop95$Abundance < 0.05] <- "< 5% abund."

trimmedITStop_99p_phyla <- unique(relabunTrimmed.phylatrimmedTop99$Phylum)
trimmedITStop_99p_phyla

trimmedITStop_99.5p_phyla <- unique(relabunTrimmed.phylatrimmedTop99.5$Phylum)
trimmedITStop_99.5p_phyla

# Phyla comprising at least 0.5% of total abundance
phylumtrimmedPlot99.5percent <- ggplot(data=relabunTrimmed.phylatrimmedTop99.5, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
phylumtrimmedPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 0.5% of total abundance (ASVs <50 times removed)")

# "Phyla comprising at least 1% of total abundance"
phylumtrimmedPlot.99percent <- ggplot(data=relabunTrimmed.phylatrimmedTop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
phylumtrimmedPlot.99percent <- phylumtrimmedPlot.99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 1% of total abundance (ASVs <50 times removed)")
phylumtrimmedPlot.99percent

# Get exact abundances of each phyla (trimmedTop 99.5%):
colnames(relabunTrimmed.phylatrimmedTop99.5)
relabunTrimmed.phylatrimmedTop99.5[,2] #This is EU... now "Sample" because of the glomming!

trimmedITStop_99.5p_phyla <- relabunTrimmed.phylatrimmedTop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
trimmedITStop_99.5p_phyla

#####################
# trimmedTop CLASSES:
#####################

# (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_bartrimmedPlots.R)
# TURN ASVs INTO PHYLUM LEVEL
rarefiedTrimmed.class.glom <-  tax_glom(ITS_TrimmedSRS_16S.ps, taxrank = "Class") 
tax_table(rarefiedTrimmed.class.glom) # good, this is only class (38 different classes)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabunTrimmed.class.0 <- transform_sample_counts(rarefiedTrimmed.class.glom, function(x) x / sum(x) )
rownames(otu_table(relabunTrimmed.class.0)) # ASVs are just representative from each class

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunTrimmed.class.1 <- merge_samples(relabunTrimmed.class.0, group = "EU")
sample_data(relabunTrimmed.class.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunTrimmed.class.2 <- transform_sample_counts(relabunTrimmed.class.1, function(x) x / sum(x))
sample_data(relabunTrimmed.class.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabunTrimmed.class.df <-psmelt(relabunTrimmed.class.2)
dim(relabunTrimmed.class.df) #
sum(relabunTrimmed.class.df[,3]) 
colnames(relabunTrimmed.class.df) 
relabunTrimmed.classtrimmedTop99 <- relabunTrimmed.class.df
relabunTrimmed.classtrimmedTop95 <- relabunTrimmed.class.df

relabunTrimmed.classtrimmedTop99$Class[relabunTrimmed.classtrimmedTop99$Abundance < 0.01] <- "< 1% abund."
relabunTrimmed.classtrimmedTop95$Class[relabunTrimmed.classtrimmedTop95$Abundance < 0.05] <- "< 5% abund."

trimmedITStop_99p_class <- unique(relabunTrimmed.classtrimmedTop99$Class)
trimmedITStop_99p_class

trimmedITStop_95p_class <- unique(relabunTrimmed.classtrimmedTop95$Class)
trimmedITStop_95p_class

# Phyla comprising at least 5% of total abundance
classtrimmedPlot.95pt <- ggplot(data=relabunTrimmed.classtrimmedTop95, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
classtrimmedPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 5.5))  + ggtitle("Classes comprising at least 5% of total abundance (ASVs <50 times removed)")

# Make new phyloseq object by removing biocrust and controls before ordination:
trimmedITSjustsoils.ps <- subset_samples(ITS_TrimmedSRS_16S.ps, Type != "BioCrust" & Type != "ExtContWater" & Type != "ExtBlank")
unique(sample_data(trimmedITSjustsoils.ps)[,27]) #only soils
sample_names(trimmedITSjustsoils.ps) #another check to show that we have only soils!

# Plot ordination
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
trimOrd <- ordinate(trimmedITSjustsoils.ps, method = "NMDS", distance = "bray", trymax = 100) 

# With no black outline around points
trimmedBrayNMDS <- phyloseq::plot_ordination(trimmedITSjustsoils.ps, trimOrd, type= "samples", color= "EU")
quartz()
trimmedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(aes(fill=EU),color="black", pch=21, size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities (rare ASVs Removed)")

# Add in labels to figure out what weird samples are
labeledTrimmedBrayNMDS <- phyloseq::plot_ordination(trimmedITSjustsoils.ps, trimOrd, type= "samples", color= "EU", label = "Sample.ID")
quartz()
labeledTrimmedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities (rare ASVs Removed)")


##################################################################################
# III. EXPLORATION OF "OUTLIER" TAXA-- taxonomic barcharts, differential abundance analysis
##################################################################################

# Visible outliers: (going clockwise from the top left on the ordination plot "labeledTrimmedBrayNMDS" above):
outliersTrimmed <- c("52D_L_100", "52D_R_10", "53ND_T_20", "53ND_B_20", "53ND_R_10")
outliersTrimmed <- sort(outliersTrimmed) #sorted this because next lines of code will automatically sort
outliersTrimmed

# Do these outliers look weird?
# Get the rows for the outliers in sample_df
#View(samples_df)]
# Code below looks for the row numbers corresponding to each one of the outlier samples. It automatically sorts the data
outlier_index <- which(samples_df$Sample.ID=="52D_L_100" | samples_df$Sample.ID=="52D_R_10" | samples_df$Sample.ID=="53ND_T_20"
                       | samples_df$Sample.ID=="53ND_B_20" | samples_df$Sample.ID=="53ND_R_10")
outlier_index 
samples_df[outlier_index,1]==outliersTrimmed # Names using outlier_index match the outliers defined above!
outlierSampNumb <- rownames(samples_df[outlier_index,]) # now these are the actual names used in bioinformatics and in my plate schemes.
# They differ from outlier_index because the row numbers are not the same as the row names.
outlierSampNumb
# Can use "outlierSampNumb" with phyloseq's "prune_samples" function

# Just looking at ASVs
# Now, that we have the numbers corresponding to each sample, we can isolate the ASVs
ITS_TrimmedoutlierSamples <- get_taxa(trimmedITSjustsoils.ps, outlierSampNumb)
colnames(ITS_TrimmedoutlierSamples) == outlierSampNumb #Nice, this worked
str(ITS_TrimmedoutlierSamples)

# ALTERNATIVELY, CAN JUST USE PRUNE_SAMPLES IN PHYLOSEQ
sample_names(trimmedITSjustsoils.ps) #sample names are the numbers that we used in the PCR plate organization!
outlierSampNumb
# Also, these are the same as those in outlier_index
outliers.ps <- prune_samples(outlierSampNumb, trimmedITSjustsoils.ps)
sort(sample_data(outliers.ps)$Sample.ID) == sort(outliersTrimmed) #these are equal to the outliers we defined above!

outliersTax <- taxtable_outta_ps(outliers.ps)
outliersASVs <- ASVs_outta_ps(outliers.ps)
# View(outliersTax)

# Do the ASVs and the community abundances look weird for the outliers?
# Make an mctoolsr like table to view if we want!

outliersASVtax <- cbind.data.frame(outliersASVs, outliersTax)
#View(outliersASVtax)

outliersASVtax$Abundance <- rowSums(outliersASVtax[,1:5]) #What is the abundance of each ASV across all 5 samples? 
# Some of these outliers share a high amount of ASVs in the genus Archaeorhizomyces

# What do the top phyla and classes look like?

# Top Phyla
# TURN ASVs INTO PHYLUM LEVEL
outliers.phylum.glom <-  tax_glom(outliers.ps, taxrank = "Phylum") 
tax_table(outliers.phylum.glom) # only 10 phyla
sample_data(outliers.phylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
outliers.phylum.0 <- transform_sample_counts(outliers.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(outliers.phylum.0)) # I think that the ASVs are just representative from each class
sample_data(outliers.phylum.0)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
outliers.phylum.df <-psmelt(outliers.phylum.0)
dim(outliers.phylum.df ) 
sum(outliers.phylum.df[,3]) #nice this is 5, the number of samples!
colnames(outliers.phylum.df) 
outliers.phylumTop99 <- outliers.phylum.df
outliers.phylumTop95 <- outliers.phylum.df

outliers.phylumTop99$Phylum[outliers.phylumTop99$Abundance < 0.01] <- "< 1% abund."
outliers.phylumTop95$Phylum[outliers.phylumTop95$Abundance < 0.05] <- "< 5% abund."

outliers.phylumTop_99p <- unique(outliers.phylumTop99$Phylum)
outliers.phylumTop_99p #"p__Basidiomycota"     "p__Ascomycota"        "p__Mortierellomycota" "p__Mucoromycota"      "p__Glomeromycota"    
# [6] "< 1% abund."  
outliers.phylumTop_95p <- unique(outliers.phylumTop95$Phylum)
outliers.phylumTop_95p #"p__Basidiomycota"     "p__Ascomycota"        "p__Mortierellomycota" "p__Mucoromycota"      "< 5% abund."  

# Phyla comprising at least 1% of total abundance
outliers.phylumPlot.99pt <- ggplot(data=outliers.phylumTop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
outliers.phylumPlot.95pt <- outliers.phylumPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#999999", "#f781bf", "#a65628", "#ffff33", "#ff7f00", "#984ea3")) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 6))  + ggtitle("Outliers (after trimming): Phyla at least 1% of total abundance")
outliers.phylumPlot.95pt

# Compare this with soils aggregated by site/EU (trimmed, so aftr removing rare taxa):
# First, get the top phyla of just soil

# TURN ASVs INTO PHYLUM LEVEL
justsoilsphylum.glom <-  tax_glom(trimmedITSjustsoils.ps, taxrank = "Phylum") 
tax_table(justsoilsphylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
justsoils.phyla.0 <- transform_sample_counts(justsoilsphylum.glom, function(x) x / sum(x) )
rownames(otu_table(justsoils.phyla.0))

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
justsoils.phyla.1 <- merge_samples(justsoils.phyla.0, group = "EU")
sample_data(justsoils.phyla.1) 

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
justsoils.phyla.2 <- transform_sample_counts(justsoils.phyla.1, function(x) x / sum(x))
sample_data(justsoils.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
justsoils.phyla.df <-psmelt(justsoils.phyla.2)
dim(justsoils.phyla.df) 
justsoils.phyla.Top95 <- justsoils.phyla.df
justsoils.phyla.Top99 <- justsoils.phyla.df

justsoils.phyla.Top95$Phylum[justsoils.phyla.Top95$Abundance < 0.05] <- "< 5% abund."
justsoils.phyla.Top99$Phylum[justsoils.phyla.Top99$Abundance < 0.01] <- "< 1% abund."

# Plot 
# Phyla comprising at least 1% of total abundance 
# Colors are different for easier comparison with outliers in next section
justsoils.phylaPlot.95percent <- ggplot(data=justsoils.phyla.Top95, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
quartz()
justsoils.phylaPlot.95percent <- justsoils.phylaPlot.95percent + geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = c("#999999", "#f781bf", "#a65628", "#ffff33")) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 6))  + ggtitle("Soils (after trimming): Phyla at least 5% of total abundance")
justsoils.phylaPlot.95percent 

quartz()
grid.arrange(outliers.phylumPlot.95pt, justsoils.phylaPlot.95percent, nrow=1)


##################################################################################
# IV. COMPARING FOREST AND PATCHES WITH ORDINATIONS AND DIFFERENTIAL ABUNDANCE ANALYSIS
##################################################################################

#################
# Ordination of forest versus patch
#################

# Added "habitat" column to metadata in Excel, re-ran whole script
HabitatBrayNMDS <- phyloseq::plot_ordination(trimmedITSjustsoils.ps, trimOrd, type= "samples", color= "Habitat")
HabitatBrayNMDS <- HabitatBrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("(ITS) NMDS based on Bray-Curtis Dissimilarities (after trimming rare ASVs)")
quartz()
HabitatBrayNMDS 
# Cool, you can see that the forest and the patch separate out, with edge somewhat in between!

########################
# Does the Bray-Curtis distance b/w samples tend to increase with distance?
########################
ASVtab <- psotu2veg(trimmedITSjustsoils.ps)
#View(ASVtab)

ASVbrayDist <- vegdist(ASVtab, method = "bray")
ASVBrayDist.mat <- as.matrix(ASVbrayDist)
diag(ASVBrayDist.mat) <- NA

# Make data frame for indexing matrix
rownames(ASVBrayDist.mat) == rownames(sample_data(trimmedITSjustsoils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps <- colnames(ASVBrayDist.mat)
meter <- sample_data(trimmedITSjustsoils.ps)$Meter
index.df <- data.frame(samps, meter)

# Get row and column numbers of different meters
m10 <- which(index.df$meter == 10)
m20 <- which(index.df$meter == 20)
m30<- which(index.df$meter == 30)
m40 <- which(index.df$meter == 40)
m50 <- which(index.df$meter == 50)
m60 <- which(index.df$meter == 60)
m70 <- which(index.df$meter == 70)
m80<- which(index.df$meter == 80)
m90 <- which(index.df$meter == 90)
m100 <- which(index.df$meter == 100)

# Distances between patch at 10m and each point along transect
m10_comp10 <- ASVBrayDist.mat[m10, m10]
m10_comp20 <- ASVBrayDist.mat[m10, m20]
m10_comp30 <- ASVBrayDist.mat[m10, m30]
m10_comp40 <- ASVBrayDist.mat[m10, m40]
m10_comp50 <- ASVBrayDist.mat[m10, m50]
m10_comp60 <- ASVBrayDist.mat[m10, m60]
m10_comp70 <- ASVBrayDist.mat[m10, m70]
m10_comp80 <- ASVBrayDist.mat[m10, m80]
m10_comp90 <- ASVBrayDist.mat[m10, m90]
m10_comp100 <- ASVBrayDist.mat[m10, m100]

# Distances between forest at 100m and each point along transect
m100_comp100 <- ASVBrayDist.mat[m100, m100]
m100_comp90 <- ASVBrayDist.mat[m100, m90]
m100_comp80 <- ASVBrayDist.mat[m100, m80]
m100_comp70 <- ASVBrayDist.mat[m100, m70]
m100_comp60 <- ASVBrayDist.mat[m100, m60]
m100_comp50 <- ASVBrayDist.mat[m100, m50]
m100_comp40 <- ASVBrayDist.mat[m100, m40]
m100_comp30 <- ASVBrayDist.mat[m100, m30]
m100_comp20 <- ASVBrayDist.mat[m100, m20]
m100_comp10 <- ASVBrayDist.mat[m100, m10]

######## MAKE BOXPLOTS ######## 
bold_a <- expression(bold("Dissimilarity Relative to 10m (patch)"))
bold_b <- expression(bold("Dissimilarity Relative to 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box <- boxplot(list(m10_comp10, m10_comp20, m10_comp30, m10_comp40,
                          m10_comp50, m10_comp60, m10_comp70,
                          m10_comp80, m10_comp90, m10_comp100),
                     ylab = "Bray-Curtis Dissimilarity",
                     names = c("10", "20", "30", "40", "50",
                               "60", "70", "80", "90", "100"), cex.axis = 0.8,
                     xlab = "meter along transect",
                     cex.lab = 1,
                     ylim=c(0.0, 1.0))
mtext(text=bold_a, side=3, adj = -0.065, line = 2)
forest_box <- boxplot(list(m100_comp10, m100_comp20, m100_comp30, m100_comp40,
                           m100_comp50, m100_comp60, m100_comp70,
                           m100_comp80, m100_comp90, m100_comp100), 
                      ylab = "Bray-Curtis Dissimilarity",
                      names = c("10", "20", "30", "40", "50",
                                "60", "70", "80", "90", "100"), cex.axis = 0.8,
                      xlab = "meter along transect",
                      cex.lab = 1,
                      ylim=c(0.0, 1.0))
mtext(text=bold_b, side=3, adj = -0.065, line = 2)

##################################
# SAVE ALL OF THESE FOR EASY ACCESS
##################################

#Saved last July 31, 2022
save(ITSrarefied.ps, samples_df, ASVsTrimmed, taxTrimmed, trimmedITSjustsoils.ps, outliersASVtax, file = "RobjectsSaved/EDA_ITSJuly2022")
save(trimmedITSjustsoils.ps, file= "RobjectsSaved/trimmedITSjustsoils.ps") 
ITS_noNAsPreRarefaction_ASVTab <- ASVs_outta_ps(ITS_noNAs.ps) #this is the dataset pre-rarefaction, but after removing non-fungal ASVs andd ASVs not assigned at least to phylum.
ITS_noNAsPreRarefaction_TaxTable <- taxtable_outta_ps(ITS_noNAs.ps)

# The following three lines are each separate part of the phyloseq object before rarefaction (but after removing non-fungal ASVs andd ASVs not assigned at least to phylum.)
write.csv(ITS_noNAsPreRarefaction_ASVTab, file = "RobjectsSaved/ITS_noNAsPreRarefaction_ASVTab.csv")
write.csv(ITS_noNAsPreRarefaction_TaxTable, file= "RobjectsSaved/ITS_noNAsPreRarefaction_TaxTable.csv")
write.csv(samples_df, file= "RobjectsSaved/ITS_FinalMetadata.csv") 