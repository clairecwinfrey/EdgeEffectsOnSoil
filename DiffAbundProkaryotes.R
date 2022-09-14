# DiffAbundProkaryotes.R
# May 4, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect BACTERIAL AND ARCHAEAL (16S) ASVs using differential abundance analyses
# (DESeq2). This script:
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. postUbiquity.ps)).
# 2. Creates data frame with z-score (of ASV abundance) for each ASV within each EU 
# 3. Fits logistic curves to each one of the differentially abundant ASVs:
# -- i. gets Z-score (of ASV abundance) for each ASV within each EU 
# -- ii. for each ASV, then fits logistic models to THIS
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential abundance breakdown by forest, patch,
# and non-differentially abundant microbes 

# IdentifyEdgeEffectTaxa5.R investigated linear fit and logistic model fit of the
# the data, and found that logistic models were much better.

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("readxl")       # necessary to import the data from Excel file
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("grid")
library("drc") # might use for log fit thing

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: UbiquityMedianSetup.R

######
# FUNCTIONS DEFINED IN THIS SCRIPT (but often used first in other scripts):
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVs_outta_ps gets the ASV table out of phyloseq 
# from: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  if(taxa_are_rows(ASVTable)) {
    ASVTable <- t(ASVTable)
  }
  return(as.data.frame(ASVTable))
}

##########################################################################
# INDICATOR ANALYSIS AND TIDY DATAFRAME
##########################################################################

# Perform indicator species analysis on the post ubiquity dataset, and then
# make a dataframe which has a row corresponging to the abundance of each differentially
# abundant ASV in each place along the transect, as well as taxonomic info, and 
# whether or not that ASV was differentially abundant in patch or the forest.

# The next steps add 1 to all of the ASV abundance counts and re-make a phyloseq object.
# This was necessary because DESeq2 could not compute the geometric mean of the samples
# since every gene had at least one zero in it (and was thus thrown out of the analysis)
# (see recommendations and explanations here:
# https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564)
postUbiqASVs_16S <- ASVs_outta_ps(postUbiquity.ps)
postUbiqASVs_16SPlus1 <- t(postUbiqASVs_16S + 1) #adding one to every abundance count for differential abundance analysis;
# see explanation a few lines down. Invert so that I can make into a new phyloseq object.
# Make a new phyloseq object with these above
OTU = otu_table(postUbiqASVs_16SPlus1, taxa_are_rows = TRUE)
postUbiqASVs16S_Plus1.ps <- phyloseq(OTU, tax_table(postUbiquity.ps), sample_data(postUbiquity.ps))

# Remove edge samples to compare patch and matrix
NoEdgePlus1.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, Habitat != "edge") 
NoEdgePlus1.ps
step1 <- phyloseq::phyloseq_to_deseq2(NoEdgePlus1.ps, ~ Habitat) #set up DESeq2 dds object 
step2 <- DESeq2::DESeq(step1, test="Wald", fitType = "parametric") #differential expression analysis step;
#uses default Benjamini-Hochberg correction 
DeseqResults <- results(step2, cooksCutoff = FALSE) #make results object
DeseqResults <- DeseqResults[which(DeseqResults$padj < 0.001), ] #only get those ASVs below alpha level
DeseqResults <- cbind(as(DeseqResults, "data.frame"), as(tax_table(NoEdgePlus1.ps)[rownames(DeseqResults), ], "matrix")) #clean up format
DeseqResults$Habitat <- ifelse(DeseqResults$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
# Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
# View(DeseqResults) #2,169 diff abund ASVs)
sampDat <- sample_data(postUbiqASVs16S_Plus1.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ASVnamesDA <- rownames(DeseqResults) #2169 found
ASVsAll <- as.data.frame(t(ASVs_outta_ps(postUbiqASVs16S_Plus1.ps))) #get ASVs from original and transpose so that ASVs are rows 
# Create dataframe with everything of interest
diffAbunDat <- merge(DeseqResults, ASVsAll, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
diffAbunDat_tidy_PROKARYOTES <- diffAbunDat %>% 
  pivot_longer(cols= 16:ncol(diffAbunDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_tidy_PROKARYOTES)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
#View(diffAbunDat_tidy_PROKARYOTES) this has 505,377 rows, which is equal to 233 (number of samples) x 2,169 (number of diff abund ASVs)

save(diffAbunDat_tidy_PROKARYOTES, file= "RobjectsSaved/diffAbunDat_tidy_PROKARYOTES") #last saved September 13, 2022
##########

##########################################################################
# STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################

DeseqResultsMini <- DeseqResults[,c(8,14)]  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
postUbiqTaxTab <- taxtable_outta_ps(postUbiqASVs16S_Plus1.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == TRUE)) #2169
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE)) #3050 ASVs are NOT differentially abundant?
false_index <- which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE) #also 3050
# Now, construct a dataframe with the taxa that were not differentially abundant 
notDA_taxTab <- postUbiqTaxTab[false_index,] #get a taxonomy tab with ONLY the non-differentially abundant ASVs
notDA_taxTab$Habitat <- "AremainingASVs" #make a habitat column that labels these as NOT differentially abundant. A in front so that would be first
# in ggplot for ease.
# View(notDA_taxTab)
colnames(notDA_taxTab)
notDA_taxTabMini <- notDA_taxTab[,c(2,8)] #keep only phylum and habitat to match DeseqResultsMini
DAphylumAll <- rbind(DeseqResultsMini, notDA_taxTabMini) #this has ASV name, phylum, and whether diff abundant for ALL ASVs in postUbiquity analysis
# What are the phyla breakdown here?
# so effectively, we want to get numbers in each phyla in each group. Should just be able to plot this?

diffAbund_16S_stackedBarplotPhyla <- ggplot(DAphylumAll, aes(fill=Habitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant Bacterial and Archaeal ASVs")
# quartz()
diffAbund_16S_stackedBarplotPhyla
# Below saved August 1, 2022 so that it can be added to a 2 paneled plot with fungal plot!
# save(diffAbund_16S_stackedBarplotPhyla, file="RobjectsSaved/diffAbund_16S_stackedBarplotPhyla_plot")

# A few checks to make sure that the counting above is working as expected
# Acidobacteria
length(which(DAphylumAll$Phylum=="Acidobacteria")) #1043
acido_index <- which(DAphylumAll$Phylum=="Acidobacteria")
length(which(DAphylumAll[acido_index,]$Habitat=="forest")) #342 forest specialists within Acidobacteria
length(which(DAphylumAll[acido_index,]$Habitat=="patch")) #174 patch specialists within Acidobacteria
length(which(DAphylumAll[acido_index,]$Habitat=="AremainingASVs")) #527 remaining ASVs
(342+ 174 + 527) == length(which(DAphylumAll$Phylum=="Acidobacteria"))

#Chloroflexi
length(which(DAphylumAll$Phylum=="Chloroflexi")) #498
chloro_index <- which(DAphylumAll$Phylum=="Chloroflexi")
length(which(DAphylumAll[chloro_index,]$Habitat=="forest")) #17 forest specialists 
length(which(DAphylumAll[chloro_index,]$Habitat=="patch")) #252 patch specialists
length(which(DAphylumAll[chloro_index,]$Habitat=="AremainingASVs")) #229 remaining ASVs
(17+252+229) == length(which(DAphylumAll$Phylum=="Chloroflexi"))

##########################################################################
# PLOT OF Chloroflexi
##########################################################################

# Within the plot above, Chloroflexi are among the most interesting groups.
# Here, makes a plot of the different relative abundances of 

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
PostUbiq_16S_NOEdge.ps <- subset_samples(postUbiquity.ps, Habitat != "edge") #remove edge so we can compare patch versus forest

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
# relabunpostUbiqNOEdge.phyla.0 <- transform_sample_counts(PostUbiq_16S_NOEdge.ps.phylum.glom, function(x) x / sum(x) )
relabunpostUbiqNOEdge_16S.allASVs <- transform_sample_counts(PostUbiq_16S_NOEdge.ps, function(x) x / sum(x) )

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunpostUbiqNOEdge.allASVs <- merge_samples(relabunpostUbiqNOEdge_16S.allASVs, group = c("Habitat"))
sample_data(relabunpostUbiqNOEdge_16S.allASVs) #shows that we still have samples from each EU, biocrust, and extcontrol (water)
# Meter did an averaging thing; can just ignore it

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunpostUbiqNOEdge_16S.allASVs.2 <- transform_sample_counts(relabunpostUbiqNOEdge.allASVs, function(x) x / sum(x))
sample_data(relabunpostUbiqNOEdge_16S.allASVs.2)

# 
relabunpostUbiqNOEdge_16S.allASVs.df <-psmelt(relabunpostUbiqNOEdge_16S.allASVs.2)
head(relabunpostUbiqNOEdge_16S.allASVs.df) #
dim(relabunpostUbiqNOEdge_16S.allASVs.df) #10438, 33 because 5219 taxa

chloro_index <- which(relabunpostUbiqNOEdge_16S.allASVs.df$Phylum=="Chloroflexi")
chloroRelAbund <- relabunpostUbiqNOEdge_16S.allASVs.df[chloro_index, ]
#View(chloroRelAbund)
colnames(chloroRelAbund)[3] <- "Relative_abundance"
colnames(chloroRelAbund)[2] <- "Habitat_type"

# Make boxplot of relative abundances of each ASV in the Chloroflexi
chloroBoxPlot <- ggplot(chloroRelAbund, aes(x=Habitat_type, y=Relative_abundance, fill= Habitat_type)) +
  geom_boxplot() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  labs(title="Relative abundance of Chloroflexi ASVs", x="Habitat Type", y = "Relative abundance")
quartz()
chloroBoxPlot


