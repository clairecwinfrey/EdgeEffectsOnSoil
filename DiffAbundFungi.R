# DiffAbundFungi
# Aug 14, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect FUNGAL (ITS) ASVs using differential abundance 
# analyses (DESeq2) and makes relevant plots. It comes BEFORE fungalZscores.

# This script:
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. ITS_postUbiquity.ps)).
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential 
# abundance breakdown by forest, patch, and non-differentially abundant microbes 

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("reshape2")
library("growthcurver") #fits logistic curve to the abundance (z-score) data

# Load data
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: UbiquitypostUbiqSetup.R

######

# FUNCTIONS DEFINED IN THIS SCRIPT (but often used first in other scripts):
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVsampOccurrence determines the number of samples that each ASV occurs in 
# For proof that ASVsampOccurrence works, see UbiquitypostUbiqSetup.R
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

############################################################

##########################################################################
# 1.    INDICATOR ANALYSIS AND TIDY DATAFRAME
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
ITS_postUbiqASVs <- ASVs_outta_ps(ITS_postUbiquity.ps)
ITS_postUbiqASVsPlus1 <- t(ITS_postUbiqASVs + 1) #adding one to every abundance count for differential abundance analysis;
# see explanation a few lines down. Invert so that I can make into a new phyloseq object.
# Make a new phyloseq object with these above
OTU = otu_table(ITS_postUbiqASVsPlus1, taxa_are_rows = TRUE)
ITS_postUbiqASVsPlus1.ps <- phyloseq(OTU, tax_table(ITS_postUbiquity.ps), sample_data(ITS_postUbiquity.ps))

# Remove edge samples to compare patch and matrix
NoEdgePlus1.ps <- subset_samples(ITS_postUbiqASVsPlus1.ps, Habitat != "edge") 
step1 <- phyloseq::phyloseq_to_deseq2(NoEdgePlus1.ps, ~ Habitat) #set up DESeq2 dds object 
step2 <- DESeq2::DESeq(step1, test="Wald", fitType = "parametric") #differential expression analysis step
DeseqResults <- results(step2, cooksCutoff = FALSE) #make results object
DeseqResults <- DeseqResults[which(DeseqResults$padj < 0.001), ] #only get those ASVs below alpha level
DeseqResults <- cbind(as(DeseqResults, "data.frame"), as(tax_table(NoEdgePlus1.ps)[rownames(DeseqResults), ], "matrix")) #clean up format
DeseqResults$Habitat <- ifelse(DeseqResults$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
# Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
# View(DeseqResults) #287 diff abund ASVs)
sampDat <- sample_data(ITS_postUbiqASVsPlus1.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ASVnamesDA_FUNGI <- rownames(DeseqResults) #1696 found
ASVsAll <- as.data.frame(t(ASVs_outta_ps(ITS_postUbiqASVsPlus1.ps))) #get ASVs from original and transpose so that ASVs are rows 
# Create dataframe with everything of interest
diffAbunDat <- merge(DeseqResults, ASVsAll, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
diffAbunDat_tidy_FUNGI <- diffAbunDat %>% 
  pivot_longer(cols= 16:ncol(diffAbunDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_tidy_FUNGI)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
#View(diffAbunDat_tidy_FUNGI) this has 66,871 rows, which is equal to 233 (number of samples) x 287 (number of diff abund ASVs)
##########

#save(diffAbunDat_tidy_FUNGI, file= "RobjectsSaved/diffAbunDat_tidy_FUNGI") #saved Sept 9, 2022

# Re-make phyloseq object with just the differentially abundant taxa for ease of working with it in the future!
diffAbundFungiNames <- as.character(unique(diffAbunDat_tidy_FUNGI$ASV_name)) #added as character because at first it had "AsIs" attribute
fungiDiffAbund.ps <- phyloseq::prune_taxa(diffAbundFungiNames, ITS_postUbiquity.ps)
otu_table(fungiDiffAbund.ps)

# save(fungiDiffAbund.ps, file= "RobjectsSaved/fungiDiffAbund.ps") #saved September 21, 2022
##########################################################################
# 2. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################

DeseqResultsMini <- DeseqResults[,c(8,14)]  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
postUbiqTaxTab <- taxtable_outta_ps(ITS_postUbiqASVsPlus1.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == TRUE)) #287
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE)) #126 ASVs are NOT differentially abundant?
false_index <- which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE) #also 126
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

diffAbund_ITS_stackedBarplotPhyla <- ggplot(DAphylumAll, aes(fill=Habitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  scale_x_discrete(labels=c("p__Ascomycota" = "Ascomycota", "p__Basidiomycota" = "Basidiomycota",
                            "p__Calcarisporiellomycota" = "Calcarisporiellomycota",
                            "p__Glomeromycota" = "Glomeromycota", "p__Mortierellomycota"="Mortierellomycota",
                            "p__Mucoromycota" = "Mucoromycota", "p__Olpidiomycota"= "Olpidiomycota",
                            "p__Rozellomycota" = "Rozellomycota")) + 
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant Fungal ASVs") 
# quartz()
diffAbund_ITS_stackedBarplotPhyla

# A few checks to make sure that the counting above is working as expected
# Ascomycota
length(which(DAphylumAll$Phylum=="p__Ascomycota")) #249
asco_index <- which(DAphylumAll$Phylum=="p__Ascomycota")
length(which(DAphylumAll[asco_index,]$Habitat=="forest")) #56 forest specialists within Ascomycota
length(which(DAphylumAll[asco_index,]$Habitat=="patch")) #120 patch specialists within Ascomycota
length(which(DAphylumAll[asco_index,]$Habitat=="AremainingASVs")) #73 remaining ASVs within Ascomycota
(56+ 120 + 73) == length(which(DAphylumAll$Phylum=="p__Ascomycota"))
# this looks correct on the plot too!

# Construct two-paneled figure with 16S and ITS differentially abundant stacked barcharts side by side
# Load in previously made 16S figure (made in "EdgeEffectsbyASV_allSites.R")
load(file="RobjectsSaved/diffAbund_16S_stackedBarplotPhyla_plot")
# quartz()
grid.arrange(diffAbund_16S_stackedBarplotPhyla, diffAbund_ITS_stackedBarplotPhyla, nrow=2)

########################
# PLOT OF GLOMEROMYCOTA
########################

# Within the plot above, glomeromycota are among the most interesting groups.
# Here, makes a plot of the different relative abundances of 

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
ITS_postUbiqNOEdge.ps <- subset_samples(ITS_postUbiquity.ps, Habitat != "edge") #remove edge so we can compare patch versus forest

#ITS_postUbiqNOEdge.ps.phylum.glom <-  tax_glom(ITS_postUbiqNOEdge.ps, taxrank = "Phylum") 
#tax_table(ITS_postUbiqNOEdge.ps.phylum.glom) #8 phyla
#sample_data(ITS_postUbiqNOEdge.ps.phylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
# relabunpostUbiqNOEdge.phyla.0 <- transform_sample_counts(ITS_postUbiqNOEdge.ps.phylum.glom, function(x) x / sum(x) )
relabunpostUbiqNOEdge.allASVs <- transform_sample_counts(ITS_postUbiqNOEdge.ps, function(x) x / sum(x) )

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunpostUbiqNOEdge.allASVs <- merge_samples(relabunpostUbiqNOEdge.allASVs, group = c("Habitat"))
sample_data(relabunpostUbiqNOEdge.allASVs) #shows that we still have samples from each EU, biocrust, and extcontrol (water)
# Meter did an averaging thing; can just ignore it

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunpostUbiqNOEdge.allASVs.2 <- transform_sample_counts(relabunpostUbiqNOEdge.allASVs, function(x) x / sum(x))
sample_data(relabunpostUbiqNOEdge.allASVs.2)

# 
relabunpostUbiqNOEdge.allASVs.df <-psmelt(relabunpostUbiqNOEdge.allASVs.2)
head(relabunpostUbiqNOEdge.allASVs.df) #
dim(relabunpostUbiqNOEdge.allASVs.df) 

glomero_index <- which(relabunpostUbiqNOEdge.allASVs.df$Phylum=="p__Glomeromycota")
glomeroRelAbund <- relabunpostUbiqNOEdge.allASVs.df[glomero_index, ]
#View(glomeroRelAbund)
colnames(glomeroRelAbund)[3] <- "Relative_abundance"
colnames(glomeroRelAbund)[2] <- "Habitat_type"

# View(relabunpostUbiqNOEdge.phyla.df)

# Make boxplot of relative abundances of each ASV in the phylum
glomeroBoxPlot <- ggplot(glomeroRelAbund, aes(x=Habitat_type, y=Relative_abundance, fill= Habitat_type)) +
  geom_boxplot() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  labs(title="Relative abundance of Glomeromycota ASVs", x="Habitat Type", y = "Relative abundance")
quartz()
glomeroBoxPlot

# Get exact abundances of each phyla 
colnames(relabunpostUbiqNOEdge.phyla.df)
#relabunpostUbiqNOEdge.phyla.df[,2] #Now just all patch and forest!


# This plot is just to check that numbers look right!
relabunpostUbiqNOEdgePlot <- ggplot(data=relabunpostUbiqNOEdge.phyla.df, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
#quartz()
relabunpostUbiqNOEdgePlot + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Fungal phyla comprising at least 0.5% of total abundance (median and post-ubiquity)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#top_99.5p_phyla <- relabunpostUbiqNOEdge.phyla.df %>%
group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 



ITS_postUbiquity.ps
tax_table(ITS_postUbiquity.ps)

View(tax_table(ITS_postUbiquity.ps))


