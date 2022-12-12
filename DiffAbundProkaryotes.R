# DiffAbundProkaryotes.R
# May 4, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect BACTERIAL AND ARCHAEAL (16S) ASVs using differential abundance analyses
# (DESeq2). This script:

# PART 1: Looks at ALL EUs
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. postUbiquity.ps)).
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential abundance breakdown by forest, patch,
# and non-differentially abundant microbes 
# 3. Stacked barcharts of the abundance of non-specialists, forest specialists, and patch specialists across the transect
# 4. Comparisons and plots of ubiquity of non-specialists, forest specialists, and patch specialists 
# 5. What is the median abundance of each ASVs across samples?

# PART 2: Performs differential abundance analysis removing EU 53N from consideration, given how different its turnover from forest to patch is 
#c(supported by fis in DissAcrossTransectDiffAbundOnly.R and DissimilaritiesAcrosTransectsMedianAbund.R)

# (IdentifyEdgeEffectTaxa5.R investigated linear fit and logistic model fit of the
# the data, and found that logistic models were much better.)

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

# 3. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

#######################################################################################################
# PART 1: ALL EUS CONSIDERED
#######################################################################################################

##########################################################################
# 1. INDICATOR ANALYSIS AND TIDY DATAFRAME
##########################################################################

# Perform indicator species analysis on the post ubiquity dataset, and then
# make a dataframe which has a row corresponding to the abundance of each differentially
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
OTU <- otu_table(postUbiqASVs_16SPlus1, taxa_are_rows = TRUE)
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

#save(diffAbunDat_tidy_PROKARYOTES, file= "RobjectsSaved/diffAbunDat_tidy_PROKARYOTES") #last saved September 13, 2022
##########

# Re-make phyloseq object with just the differentially abundant taxa for ease of working with it in the future!
diffAbundProkaryotesNames <- as.character(unique(diffAbunDat_tidy_PROKARYOTES$ASV_name)) #added as character because at first it had "AsIs" attribute
prokaryotesDiffAbund.ps <- phyloseq::prune_taxa(diffAbundProkaryotesNames, postUbiquity.ps)
otu_table(prokaryotesDiffAbund.ps)

# save(prokaryotesDiffAbund.ps, file= "RobjectsSaved/prokaryotesDiffAbund.ps") #saved September 21, 2022

# Extend ASV table to have info on where each ASV specializes
ASVsPostUbiq <- t(ASVs_outta_ps(postUbiquity.ps))
DeseqResults_2 <- DeseqResults
colnames(DeseqResults_2)[14] <- "diffAbundHabitat" #rename this so that what habitat means is clearer
ASVsAllDiffAbund <- merge(ASVsPostUbiq, DeseqResults_2[,c(1:6,14)], by="row.names", all.x=TRUE)
colnames(ASVsAllDiffAbund)[1] <- "ASV_name"
# save(ASVsAllDiffAbund, file= "RobjectsSaved/ASVsAllDiffAbund") #saved December 11, 2022

# Make another ASV table that has taxonomy information too
taxPostUbiq <- taxtable_outta_ps(postUbiquity.ps)
taxPostUbiq <- taxPostUbiq %>% rownames_to_column(var="ASV_name") #change rownames to ASV nmae column for merging below
ASVsAllDiffAbund_tax <- merge(ASVsAllDiffAbund, taxPostUbiq, by="ASV_name", all.x=TRUE, all.y=TRUE)
ASVsAllDiffAbund_tax[which(is.na(ASVsAllDiffAbund_tax$diffAbundHabitat)),241] <- "non-specialist" #3050 are non-specialists
# The line below makes variable a factor so that the stacked barplot is in desired order
ASVsAllDiffAbund_tax$diffAbundHabitat <- factor(ASVsAllDiffAbund_tax$diffAbundHabitat, levels=c("non-specialist", "forest", "patch"))
# save(ASVsAllDiffAbund_tax, file= "RobjectsSaved/ASVsAllDiffAbund_tax") #saved December 11, 2022

##########################################################################
# 2. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################
diffAbund_16S_stackedBarplotPhyla <- ggplot(ASVsAllDiffAbund_tax, aes(fill=diffAbundHabitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey", "darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("non-specialist", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") 
# quartz()
diffAbund_16S_stackedBarplotPhyla
# Below saved December 11, 2022 so that it can be added to a 2 paneled plot with fungal plot!
# save(diffAbund_16S_stackedBarplotPhyla, file="RobjectsSaved/diffAbund_16S_stackedBarplotPhyla_plot")

# A few checks to make sure that the counting above is working as expected
# Acidobacteria
length(which(ASVsAllDiffAbund_tax$Phylum=="Acidobacteria")) #1043
acido_index <- which(ASVsAllDiffAbund_tax$Phylum=="Acidobacteria")
length(which(ASVsAllDiffAbund_tax[acido_index,]$diffAbundHabitat=="forest")) #342 forest specialists within Acidobacteria
length(which(ASVsAllDiffAbund_tax[acido_index,]$diffAbundHabitat=="patch")) #174 patch specialists within Acidobacteria
length(which(ASVsAllDiffAbund_tax[acido_index,]$diffAbundHabitat=="non-specialist")) #527 remaining ASVs
(342+ 174 + 527) == length(which(DAphylumAll$Phylum=="Acidobacteria"))

#Chloroflexi
length(which(ASVsAllDiffAbund_tax$Phylum=="Chloroflexi")) #498
chloro_index <- which(ASVsAllDiffAbund_tax$Phylum=="Chloroflexi")
length(which(ASVsAllDiffAbund_tax[chloro_index,]$diffAbundHabitat=="forest")) #17 forest specialists 
length(which(ASVsAllDiffAbund_tax[chloro_index,]$diffAbundHabitat=="patch")) #252 patch specialists
length(which(ASVsAllDiffAbund_tax[chloro_index,]$diffAbundHabitat=="non-specialist")) #229 remaining ASVs
(17+252+229) == length(which(ASVsAllDiffAbund_tax$Phylum=="Chloroflexi"))

########
# PLOT OF Chloroflexi
########

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
sample_data(relabunpostUbiqNOEdge_16S.allASVs) # Meter did an averaging thing; can just ignore it

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunpostUbiqNOEdge_16S.allASVs.2 <- transform_sample_counts(relabunpostUbiqNOEdge.allASVs, function(x) x / sum(x))
sample_data(relabunpostUbiqNOEdge_16S.allASVs.2)

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
#quartz()
chloroBoxPlot

##########################################################################
# 3. PLOTS FOR DIFFERENTIALLY ABUNDANT TAXA ALONG THE TRANSECT
##########################################################################
# This part makes a stacked barplot which has the percentage of reads that are forest, patch, or non-specialists at each point along the transect

# 1. First, need to get the median ASV abundance at each point along the transect, averaging across all EUs
# This code gets median abundance for each ASV (median across 24 transects all EUs)
# at each meter.
# In other words, the EU is ignored, so the number of rows = number of ASVs in 

# This code gets median abundance for each ASV at each meter along the transect (regardless of EU)
ASVsPostUbiq <- ASVs_outta_ps(postUbiquity.ps) #rows are samples, ASV abundance is columns 
postUbiqmetaDf <- metadata_outta_ps(postUbiquity.ps) #rows are samples, columns are all of the rest of the metadata
postUbiqmetaDf <- as.data.frame(as.matrix(postUbiqmetaDf)) #make format nicer
postUbiqmetaDf$Meter <- as.numeric(postUbiqmetaDf$Meter) #make meter numeric
unique(rownames(ASVsPostUbiq) == rownames(postUbiqmetaDf)) #this is true, which shows that we can use indices from
# metadata to isolate stuff from ASV table

# This is how to pull out rownames by meter
rownames(postUbiqmetaDf)[which(postUbiqmetaDf$Meter == 10)]
meterVec <- c(10,20,30,40,50,60,70,80,90,100)
meterRowNamesIndices <- vector("list", length(meterVec)) #this will have all of the rownames that correspond to samples in each meter
names(meterRowNamesIndices) <- paste(meterVec, "m", sep="_") 
for (j in 1:length(meterVec)){
  meterRowNamesIndices[[j]] <- rownames(postUbiqmetaDf)[which(postUbiqmetaDf$Meter == meterVec[j])]
}

# For each one of the ASVs, pull out all of the samples at each point
# first, make a dataframe to hold all the final stuff:
medianASVsByMeter <- as.data.frame(matrix(nrow=10, ncol=ncol(ASVsPostUbiq)))
colnames(medianASVsByMeter) <- colnames(ASVsPostUbiq)
rownames(medianASVsByMeter) <- names(meterRowNamesIndices)

for (k in 1:length(meterRowNamesIndices)){ #get mean ASV abundance at each meter
  medianASVsByMeter[k,] <- t(colMeans(ASVsPostUbiq[meterRowNamesIndices[[k]],])) #this pulls out all of the samples at each meter.
}
#View(medianASVsByMeter)

medianASVsByMeterTidy <-  medianASVsByMeter %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = ASV_1:ASV_9931, names_to = "ASV_name", values_to = "meanASVabundance")
# View(medianASVsByMeterTidy) #looks good!
colnames(medianASVsByMeterTidy)[1] <- "meter"
DAphylumAll_2 <- rownames_to_column(DAphylumAll) #grab this so that can get info on specialists
colnames(DAphylumAll_2)[1] <- "ASV_name"
medianASVsByMeterHabitat <- merge(medianASVsByMeterTidy, DAphylumAll_2, by="ASV_name")

#### GET RELATIVE ABUNDANCES #####
# Get the total ASV abundance within each meter
ASVmeterTotal <- as.data.frame(matrix(nrow=10, ncol=1))
rownames(ASVmeterTotal) <- paste(meterVec, "m", sep="_") 
colnames(ASVmeterTotal) <- "ASVmeterTotal"
for (j in 1:length(meterVec)){
  ASVmeterTotal[j,1] <- sum(medianASVsByMeterHabitat[which(medianASVsByMeterHabitat$meter ==  rownames(ASVmeterTotal)[j]),3]) #pull out meter by meter
}
ASVmeterTotal

# Finally get relative abundances by dividing each value in each row by the meter total
relAbundDfs <- vector("list", length(meterVec)) #this will have all of the rownames that correspond to samples in each meter
names(relAbundDfs) <- rownames(ASVmeterTotal)
for (h in 1:length(relAbundDfs)){
  relAbundDfs[[h]] <- medianASVsByMeterHabitat[which(medianASVsByMeterHabitat$meter ==  names(relAbundDfs)[[h]]),] #this pulls out just the data for each meter
  relAbundDfs[[h]]$relAbund <- (relAbundDfs[[h]]$meanASVabundance/ASVmeterTotal[h,1])*100
}
sum(relAbundDfs[[1]]$relAbund) #great, these add up to be 100!

### rbind all of the relAbunds[[all meters]] and plug this into ggplot2 below
relAbund_df <- do.call("rbind", relAbundDfs)
# View(relAbund_df)

# Plot it!
# Relative abundance
level_order <- names(meterRowNamesIndices) #set this to make in correct order from 10m to 100m
transectX <- c("40 m", "30 m","20 m", "10 m", "edge", "10 m", "20 m", "30 m", "40 m", "50 m")#for re-naming axis tick marks
relAbundTransectProks_plot <- ggplot(relAbund_df, aes(x = factor(meter, level = level_order), y = relAbund, fill = Habitat)) + 
  geom_bar(stat = "identity", position = "fill")  +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= NULL, labels=c("non-specialists", "forest specialists", "patch specialists")) +
  labs(y= "relative abundance", x = "meters from edge") + 
  ggtitle("Prokaryote Specialists Across the Transect") +
  scale_x_discrete(labels=transectX) + #change x-axis tick labels
  theme(axis.title=element_text(size=14))
#quartz()
relAbundTransectProks_plot

# Saved Dec 11, 2022
#save(relAbundTransectProks_plot, file= "RobjectsSaved/relAbundTransectProks_plot")

##########################################################################
# 4. COMPARISONS AND PLOTS OF UBIQUITY IN NON-SPECIALISTS VERSUS SPECIALISTS 
##########################################################################
ASVsAllDiffAbund_tax$sampOccurrence <- rowSums(ASVsAllDiffAbund_tax[,2:234] > 0) #just select samples and get occurrence across samples
# View(ASVsAllDiffAbund_tax)

# Plot the number of occurrences in each group:
# ggplot it!
ubiquityBySpecialistPlot <- ggplot(ASVsAllDiffAbund_tax, aes(x=diffAbundHabitat, y=sampOccurrence, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgreen", "darkgray", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("(Ubiquity) number of sites") +
  xlab("Habitat Specialist Group") +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))

quartz()
ubiquityBySpecialistPlot

# WITHOUT 100 M SAMPLES since there are more samples in the forest than the patch
ASVsAllDiffAbund_tax_no100m <- ASVsAllDiffAbund_tax
# find samples taken at 100 m
index100m <- which(sample_data(postUbiquity.ps)$Meter == 100) 
samples100m <- rownames(sample_data(postUbiquity.ps))[index100m]
samples100mCols <- which(colnames(ASVsAllDiffAbund_tax) %in% samples100m) #find the column names corresponding to these 100m samples
ASVsAllDiffAbund_tax_no100m <- ASVsAllDiffAbund_tax[,-samples100mCols] #make a new object without these columns
colnames(ASVsAllDiffAbund_tax_no100m)
ASVsAllDiffAbund_tax_no100m$sampOccurrence <- rowSums(ASVsAllDiffAbund_tax_no100m[,2:210] > 0) #just select samples and get occurrence across samples
# Note: lowest uumber of samples found in here is 31 samples (not 40) because 100 m samples were removed
#View(ASVsAllDiffAbund_tax_no100m)

# Test to see if I can do a one-way ANOVA
car::leveneTest(sampOccurrence ~ diffAbundHabitat, data = ASVsAllDiffAbund_tax_no100m) # shows that variances are NOT roughly homogenous

# Since assumptions aren't met, try pairwise Wilcox with Benjamini-Hochberg correction
prokUbiquity.wilcox <- pairwise.wilcox.test(ASVsAllDiffAbund_tax_no100m$sampOccurrence, ASVsAllDiffAbund_tax_no100m$diffAbundHabitat,
                                            p.adjust.method = "BH")
prokUbiquity.wilcox #shows patch > non-specialistand forest > non-specialist

# ggplot it!
prokUbiquityBySpecialistPlot_no100m <- ggplot(ASVsAllDiffAbund_tax_no100m, aes(x=diffAbundHabitat, y=sampOccurrence, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Ubiquity (number of sites)") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Prokaryotes: ubiquity of ASVs across specialist groups") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

quartz()
prokUbiquityBySpecialistPlot_no100m

# Saved Dec 11, 2022
#save(prokUbiquityBySpecialistPlot_no100m, file= "RobjectsSaved/prokUbiquityBySpecialistPlot_no100m")

##########################################################################
# 5. MEDIAN ABUNDANCE OF EACH ASV ACROSS SAMPLES
##########################################################################
colnames(ASVsAllDiffAbund_tax)
# Get median abundance across all samples
ASVsAllDiffAbund_tax$ASVmedAbund <- rowMedians(as.matrix(ASVsAllDiffAbund_tax[,2:234]))

# Make plot for median abundance 
prokMedASVabundPlot <- ggplot(ASVsAllDiffAbund_tax, aes(x=diffAbundHabitat, y=ASVmedAbund, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Median ASV abundance") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Prokaryotes: Median ASV abundance across samples") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

# quartz()
prokMedASVabundPlot

# MAKE PLOT OF MEDIAN ABUNDANCE IN PATCH OR FOREST SAMPLES
# The plot above shows that this is not super informative, so get it across samples, in the patch, on edge, and in forest
unique(colnames(ASVsAllDiffAbund_tax[2:234]) == rownames(postUbiqmetaDf)) #this is true, which shows that we can use indices from
# metadata to isolate stuff from ASV table
patchSampleIndex <- which(postUbiqmetaDf$Meter %in% c(10,20,30,40)) #isolate patch samples
forestSampleIndex <- which(postUbiqmetaDf$Meter %in% c(60,70,80,90,100)) #isolate forest samples
ASVsAllDiffAbund_tax$patchASVmedAbund <- rowMedians(as.matrix(ASVsAllDiffAbund_tax[,patchSampleIndex+1])) #add one because columns start at 2
ASVsAllDiffAbund_tax$forestASVmedAbund <- rowMedians(as.matrix(ASVsAllDiffAbund_tax[,forestSampleIndex+1])) #add one because columns start at 2
# Finally, make a new "longer" and "tidy" data frame
ASVsAllDiffAbund_taxAbunGroups <- ASVsAllDiffAbund_tax %>% 
  pivot_longer(cols="ASVmedAbund":"forestASVmedAbund", names_to= "habitatForMedAbundCalc", values_to="medAbund")
# View(ASVsAllDiffAbund_taxAbunGroups)

# grouped boxplot
prokMedAbundByHabitatPlot <- ggplot(ASVsAllDiffAbund_taxAbunGroups, aes(x=diffAbundHabitat, y=medAbund, color=habitatForMedAbundCalc)) + 
  #geom_boxplot() +
  #scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(size=1, alpha=0.9) +
  ylab("Median ASV abundance") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Prokaryotes: Median ASV abundance across samples") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

quartz()
prokMedAbundByHabitatPlot

#######################################################################################################
# PART 2: EU 53N EXCLUDED FROM ANALYSES
#######################################################################################################

##########################################################################
# I. INDICATOR ANALYSIS AND TIDY DATAFRAME
##########################################################################

# *REMOVE EU 53N*, then perform indicator species analysis on the post ubiquity dataset, and then
# make a dataframe which has a row corresponding to the abundance of each differentially
# abundant ASV in each place along the transect, as well as taxonomic info, and 
# whether or not that ASV was differentially abundant in patch or the forest.

# The next steps add 1 to all of the ASV abundance counts and re-make a phyloseq object.
# This was necessary because DESeq2 could not compute the geometric mean of the samples
# since every gene had at least one zero in it (and was thus thrown out of the analysis)
# (see recommendations and explanations here:
# https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564)

# Remove EU 53N from analyses
# postUbiq16Sno53N.ps <- subset_samples(postUbiquity.ps, EU != "EU_53N")
# unique(sample_data(postUbiq16Sno53N.ps)$EU) #this shows that EU 53N is no longer present!
# which(rowSums(otu_table(postUbiq16Sno53N.ps))==0) #this shows that we don't drop any ASVs by getting dropping EU 53N.
# In other words, there were no ASVs that were ONLY found in EU 53N, so we are good to move to the next step

# Get just the ASVs out of postUbiq16Sno53N
postUbiq16Sno53N_ASVs <- ASVs_outta_ps(postUbiq16Sno53N.ps)
postUbiq16Sno53N_ASVs_Plus1 <- t(postUbiq16Sno53N_ASVs + 1) #adding one to every abundance count for differential abundance analysis;
# see explanation a few lines down. Invert so that I can make into a new phyloseq object.
# Make a new phyloseq object with these above
OTU = otu_table(postUbiq16Sno53N_ASVs_Plus1, taxa_are_rows = TRUE)
postUbiq16Sno53N_Plus1.ps <- phyloseq(OTU, tax_table(postUbiq16Sno53N.ps), sample_data(postUbiq16Sno53N.ps))
# this new object has 5,219 taxa across 194 samples

# Remove edge samples to compare patch and matrix
NoEdgeNo53N_Plus1.ps <- subset_samples(postUbiq16Sno53N_Plus1.ps, Habitat != "edge") 
NoEdgeNo53N_Plus1.ps
unique(sample_data(NoEdgeNo53N_Plus1.ps)$EU) #no EU 53N
unique(sample_data(NoEdgeNo53N_Plus1.ps)$Habitat) #no edge
step1_no53N <- phyloseq::phyloseq_to_deseq2(NoEdgeNo53N_Plus1.ps, ~ Habitat) #set up DESeq2 dds object 
step2_no53N <- DESeq2::DESeq(step1_no53N, test="Wald", fitType = "parametric") #differential expression analysis step;
#uses default Benjamini-Hochberg correction 
DeseqResults_no53N <- results(step2_no53N, cooksCutoff = FALSE) #make results object
DeseqResults_no53N <- DeseqResults_no53N[which(DeseqResults_no53N$padj < 0.001), ] #only get those ASVs below alpha level
DeseqResults_no53N <- cbind(as(DeseqResults_no53N, "data.frame"), as(tax_table(NoEdgeNo53N_Plus1.ps)[rownames(DeseqResults_no53N), ], "matrix")) #clean up format
DeseqResults_no53N$Habitat <- ifelse(DeseqResults_no53N$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
# Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
# View(DeseqResults_no53N) #2,169 diff abund ASVs) 
sampDat <- sample_data(postUbiq16Sno53N_Plus1.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ASVnamesDA_no53N <- rownames(DeseqResults_no53N) #2169 found, same as without EU 53N
length(ASVnamesDA_no53N) #2169 found, same as without EU 53N
ASVsAll_no53N <- as.data.frame(t(ASVs_outta_ps(postUbiq16Sno53N_Plus1.ps))) #get ASVs from original and transpose so that ASVs are rows 
# Create dataframe with everything of interest
diffAbunDat_no53N <- merge(DeseqResults_no53N, ASVsAll_no53N, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
colnames(diffAbunDat_no53N)
diffAbunDat_no53N_tidy_PROKARYOTES <- diffAbunDat_no53N %>% 
  pivot_longer(cols= 16:ncol(diffAbunDat_no53N), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_no53N_tidy_PROKARYOTES)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
#View(diffAbunDat_no53N_tidy_PROKARYOTES) this has 420,786 rows, which is equal to 233 (number of samples) x 2,169 (number of diff abund ASVs)

#### EXPLORING HOW THE PRESENCE OF ABSENCE OF EU 53N MATTERS ####
length(intersect(ASVnamesDA_no53N, ASVnamesDA)) #this shows that 2,043 ASVs are shared, but 126 diff abund ASVs differ depending on whether EU 53N is present
with53N_index <- (ASVnamesDA %in% ASVnamesDA_no53N) #shows where ASVnamesDA is in ASVnamesDA_no53N
with53N_indexF <- which(with53N_index==FALSE)
onlyIfWith53N <- ASVnamesDA[with53N_indexF] #this shows those ASVs that are only differentially abundant if 53N is included

without53N_index <- (ASVnamesDA_no53N %in% ASVnamesDA) #shows where ASVnamesDA_no53N has the same as ASVnamesDA 
without53N_index_indexF <- which(without53N_index==FALSE)
onlyIf53Nexcluded <- ASVnamesDA_no53N[without53N_index_indexF] #this shows those ASVs that are only differentially abundant if 53N is EXCLUDED
intersect(onlyIfWith53N, onlyIf53Nexcluded) #shows that these two objects are distinct



#save(diffAbunDat_no53N_tidy_PROKARYOTES, file= "RobjectsSaved/diffAbunDat_no53N_PROKARYOTES") #last saved Oct. 1, 2022
##########

# Re-make phyloseq object with just the differentially abundant taxa for ease of working with it in the future!
diffAbund_no53N_ProkaryotesNames <- as.character(unique(diffAbunDat_no53N_tidy_PROKARYOTES$ASV_name)) #added as character because at first it had "AsIs" attribute
prokaryotes_no53N_DiffAbund.ps <- phyloseq::prune_taxa(diffAbund_no53N_ProkaryotesNames, postUbiq16Sno53N.ps)
otu_table(prokaryotes_no53N_DiffAbund.ps)
tax_table(prokaryotes_no53N_DiffAbund.ps)
sample_data(prokaryotes_no53N_DiffAbund.ps)

# save(prokaryotes_no53N_DiffAbund.ps, file= "RobjectsSaved/prokaryotes_no53N_DiffAbund.ps") #saved Oct. 1, 2022

# Extend ASV table to have info on where each ASV specializes
ASVsPostUbiq <- t(ASVs_outta_ps(postUbiq16Sno53N.ps))
DeseqResults_no53N2 <- DeseqResults_no53N
colnames(DeseqResults_no53N2)[14] <- "diffAbundHabitat" #rename this so that what habitat means is clearer
ASVsAllDiffAbund <- merge(ASVsPostUbiq, DeseqResults_no53N2[,c(1:6,14)], by="row.names", all.x=TRUE)
colnames(ASVsAllDiffAbund)[1] <- "ASV_name"
# save(ASVsAllDiffAbund, file= "RobjectsSaved/ASVsAllDiffAbund") #saved December 7, 2022

# Make another ASV table that has taxonomy information too
taxPostUbiq <- taxtable_outta_ps(postUbiq16Sno53N.ps)
taxPostUbiq <- taxPostUbiq %>% rownames_to_column(var="ASV_name") #change rownames to ASV nmae column for merging below
ASVsAllDiffAbund_tax <- merge(ASVsAllDiffAbund, taxPostUbiq, by="ASV_name", all.x=TRUE)
# save(ASVsAllDiffAbund_tax, file= "RobjectsSaved/ASVsAllDiffAbund_tax") #saved December 7, 2022
##########################################################################
# 2. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################

DeseqResults_no53NMini <- DeseqResults_no53N[,c(8,14)]  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
postUbiqTaxTab_no53N <- taxtable_outta_ps(postUbiq16Sno53N_Plus1.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(postUbiqTaxTab_no53N) %in% rownames(DeseqResults_no53NMini) == TRUE)) #2169
length(which(rownames(postUbiqTaxTab_no53N) %in% rownames(DeseqResults_no53NMini) == FALSE)) #3050 ASVs are NOT differentially abundant?
false_index_no53N <- which(rownames(postUbiqTaxTab_no53N) %in% rownames(DeseqResults_no53NMini) == FALSE) #also 3050
# Now, construct a dataframe with the taxa that were not differentially abundant 
notDA_taxTab_no53N <- postUbiqTaxTab_no53N[false_index_no53N,] #get a taxonomy tab with ONLY the non-differentially abundant ASVs
notDA_taxTab_no53N$Habitat <- "AremainingASVs" #make a habitat column that labels these as NOT differentially abundant. A in front so that would be first
# in ggplot for ease.
# View(notDA_taxTab)
colnames(notDA_taxTab_no53N)
notDA_taxTabMini_no53N <- notDA_taxTab_no53N[,c(2,8)] #keep only phylum and habitat to match DeseqResults_no53NMini
DAphylumAll_no53N <- rbind(DeseqResults_no53NMini, notDA_taxTabMini_no53N) #this has ASV name, phylum, and whether diff abundant for ALL ASVs in postUbiquity analysis
# What are the phyla breakdown here?
# so effectively, we want to get numbers in each phyla in each group. Should just be able to plot this?

diffAbund_no53N_16S_stackedBarplotPhyla <- ggplot(DAphylumAll_no53N, aes(fill=Habitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant Bacterial and Archaeal ASVs (53N excluded)")
# quartz()
diffAbund_no53N_16S_stackedBarplotPhyla
# Below saved Oct 1, 2022 so that it can be added to a 2 paneled plot with fungal plot!
# save(diffAbund_no53N_16S_stackedBarplotPhyla, file="RobjectsSaved/diffAbund_no53N_16S_stackedBarplotPhyla")

# A few checks to make sure that the counting above is working as expected
# Acidobacteria
length(which(DAphylumAll_no53N$Phylum=="Acidobacteria")) #1043
acido_index_no53N <- which(DAphylumAll_no53N$Phylum=="Acidobacteria")
length(which(DAphylumAll_no53N[acido_index_no53N,]$Habitat=="forest")) #331 forest specialists within Acidobacteria
length(which(DAphylumAll_no53N[acido_index_no53N,]$Habitat=="patch")) #186 patch specialists within Acidobacteria
length(which(DAphylumAll_no53N[acido_index_no53N,]$Habitat=="AremainingASVs")) #526 remaining ASVs
(342+ 174 + 527) == length(which(DAphylumAll_no53N$Phylum=="Acidobacteria"))

#Chloroflexi
length(which(DAphylumAll_no53N$Phylum=="Chloroflexi")) #498
chloro_index_no53N <- which(DAphylumAll_no53N$Phylum=="Chloroflexi")
length(which(DAphylumAll_no53N[chloro_index_no53N,]$Habitat=="forest")) #18 forest specialists 
length(which(DAphylumAll_no53N[chloro_index_no53N,]$Habitat=="patch")) #254 patch specialists
length(which(DAphylumAll_no53N[chloro_index_no53N,]$Habitat=="AremainingASVs")) #226 remaining ASVs
(18+254+226) == length(which(DAphylumAll_no53N$Phylum=="Chloroflexi"))

#View(DAphylumAll_no53N)
