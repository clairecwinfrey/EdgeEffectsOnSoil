# DiffAbundFungi
# Aug 14, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect FUNGAL (ITS) ASVs using differential abundance 
# analyses (DESeq2) and makes relevant plots. It comes BEFORE fungalZscores.

# PART 1: Looks at ALL EUs
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. ITS_postUbiquity.ps)).
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential abundance breakdown by forest, patch,
# and non-differentially abundant microbes 

# PART 2: Does the same thing as Part 1, *BUT* Removes EU 53N from consideration, given how different its turnover from forest to patch is 
# (supported by fis in DissAcrossTransectDiffAbundOnly.R and DissimilaritiesAcrosTransectsMedianAbund.R)


###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("gridExtra")    #allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library(vegan)

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

# 4. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

#######################################################################################################
# PART 1: ALL EUS CONSIDERED
#######################################################################################################

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
OTU <- otu_table(ITS_postUbiqASVsPlus1, taxa_are_rows = TRUE)
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
ASVnamesDA_FUNGI <- rownames(DeseqResults) #287 found
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

# Extend ASV table to have info on where each ASV specializes
ITS_ASVsPostUbiq <- t(ASVs_outta_ps(ITS_postUbiquity.ps))
DeseqResults_2 <- DeseqResults
colnames(DeseqResults_2)[14] <- "diffAbundHabitat" #rename this so that what habitat means is clearer
fungiASVsAllDiffAbund <- merge(ITS_ASVsPostUbiq, DeseqResults_2[,c(1:6,14)], by="row.names", all.x=TRUE)
colnames(fungiASVsAllDiffAbund)[1] <- "ASV_name"
# save(fungiASVsAllDiffAbund, file= "RobjectsSaved/fungiASVsAllDiffAbund") #saved December 11, 2022

# Make another ASV table that has taxonomy information too
fungiTaxPostUbiq <- taxtable_outta_ps(ITS_postUbiquity.ps)
fungiTaxPostUbiq <- fungiTaxPostUbiq %>% rownames_to_column(var="ASV_name") #change rownames to ASV nmae column for merging below
fungiASVsAllDiffAbund_tax <- merge(fungiASVsAllDiffAbund, fungiTaxPostUbiq, by="ASV_name", all.x=TRUE, all.y=TRUE)
fungiASVsAllDiffAbund_tax[which(is.na(fungiASVsAllDiffAbund_tax$diffAbundHabitat)),241] <- "non-specialist" #126 are non-specialists
# The line below makes variable a factor so that the stacked barplot is in desired order
fungiASVsAllDiffAbund_tax$diffAbundHabitat <- factor(fungiASVsAllDiffAbund_tax$diffAbundHabitat, levels=c("non-specialist", "forest", "patch"))
# save(fungiASVsAllDiffAbund_tax, file= "RobjectsSaved/fungiASVsAllDiffAbund_tax") #saved December 11, 2022


##########################################################################
# 2. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################
diffAbund_ITS_stackedBarplotPhyla <- ggplot(fungiASVsAllDiffAbund_tax, aes(fill=diffAbundHabitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), labels=c("non-specialist", "forest", "patch")) +
  scale_x_discrete(labels=c("p__Ascomycota" = "Ascomycota", "p__Basidiomycota" = "Basidiomycota",
                            "p__Calcarisporiellomycota" = "Calcarisporiellomycota",
                            "p__Glomeromycota" = "Glomeromycota", "p__Mortierellomycota"="Mortierellomycota",
                            "p__Mucoromycota" = "Mucoromycota", "p__Olpidiomycota"= "Olpidiomycota",
                            "p__Rozellomycota" = "Rozellomycota")) + 
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,400,by=25)) +
  ylab("number of ASVs in phylum") 
# quartz()
diffAbund_ITS_stackedBarplotPhyla

# A few checks to make sure that the counting above is working as expected
# Ascomycota
length(which(fungiASVsAllDiffAbund_tax$Phylum=="p__Ascomycota")) #249
asco_index <- which(fungiASVsAllDiffAbund_tax$Phylum=="p__Ascomycota")
length(which(fungiASVsAllDiffAbund_tax[asco_index,]$diffAbundHabitat=="forest")) #56 forest specialists within Ascomycota
length(which(fungiASVsAllDiffAbund_tax[asco_index,]$diffAbundHabitat=="patch")) #120 patch specialists within Ascomycota
length(which(fungiASVsAllDiffAbund_tax[asco_index,]$diffAbundHabitat=="non-specialist")) #73 remaining ASVs within Ascomycota
(56+ 120 + 73) == length(which(fungiASVsAllDiffAbund_tax$Phylum=="p__Ascomycota"))
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
sample_data(relabunpostUbiqNOEdge.allASVs) 
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


# View(tax_table(ITS_postUbiquity.ps))

##########################################################################
# 3. PLOTS FOR DIFFERENTIALLY ABUNDANT TAXA ALONG THE TRANSECT
##########################################################################
# This part makes a stacked barplot which has the percentage of reads that are forest, patch, or non-specialists at each point along the transect
# This code gets mean abundance for each ASV at each meter along the transect (regardless of EU))

# 1. First, need to get the mean ASV abundance at each point along the transect, averaging across all EUs

fungiASVsdf <- ASVs_outta_ps(ITS_postUbiquity.ps) #rows are samples, ASV abundance is columns 
fungimetaDf <- metadata_outta_ps(ITS_postUbiquity.ps) #rows are samples, columns are all of the rest of the metadata
fungimetaDf <- as.data.frame(as.matrix(fungimetaDf)) #make format nicer
fungimetaDf$Meter <- as.numeric(fungimetaDf$Meter) #make meter numeric
unique(rownames(fungiASVsdf) == rownames(fungimetaDf)) #this is true, which shows that we can use indices from
# metadata to isolate stuff from ASV table

# Get all the rownames in ASV table for samples corresponding to each meter
rownames(fungimetaDf)[which(fungimetaDf$Meter == 10)] #this is how to pull out rownames by meter
meterVec <- c(10,20,30,40,50,60,70,80,90,100)
fungiMeterRowNamesIndices <- vector("list", length(meterVec)) #this will have all of the rownames that correspond to samples in each meter
names(fungiMeterRowNamesIndices) <- paste(meterVec, "m", sep="_") 
for (j in 1:length(meterVec)){
  fungiMeterRowNamesIndices[[j]] <- rownames(fungimetaDf)[which(fungimetaDf$Meter == meterVec[j])]
}

# For each one of the ASVs, pull out all of the samples at each point
# first, make a dataframe to hold all the final stuff:
fungiMeanASVsByMeter <- as.data.frame(matrix(nrow=10, ncol=ncol(fungiASVsdf)))
colnames(fungiMeanASVsByMeter) <- colnames(fungiASVsdf)
rownames(fungiMeanASVsByMeter) <- names(fungiMeterRowNamesIndices)

# For loop to get mean ASV abundance at each meter
for (k in 1:length(fungiMeterRowNamesIndices)){ 
  fungiMeanASVsByMeter[k,] <- t(colMeans(fungiASVsdf[fungiMeterRowNamesIndices[[k]],])) #this pulls out all of the samples at each meter.
}
# View(fungiMeanASVsByMeter)
# Check that mean abundances are correct
# ASV 1 in 10m samples
tenMeterSamps <- rownames(fungimetaDf)[which(fungimetaDf$Meter==10)]
# for ASV 1
mean(fungiASVsdf[tenMeterSamps,1]) == fungiMeanASVsByMeter[1,1]
# ASV 6 in meter 100
hundredMeterSamps <- rownames(fungimetaDf)[which(fungimetaDf$Meter==100)]
mean(fungiASVsdf[hundredMeterSamps,6]) == fungiMeanASVsByMeter[10,6]

# Make this in long, "tidy" format, and then merge with object made earlier for specialist and taxonomic info
fungiMeanASVsByMeterTidy <- fungiMeanASVsByMeter %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = ASV_1:ASV_2512, names_to = "ASV_name", values_to = "meanASVabundance")
# View(fungiMeanASVsByMeterTidy) #looks good!
colnames(fungiMeanASVsByMeterTidy)[1] <- "meter"
fungiMeanASVsByMeterHabitat <- merge(fungiMeanASVsByMeterTidy, fungiASVsAllDiffAbund_tax[,c(1,241:247)], by="ASV_name") #add in which habitat specialist in & taxonomic info

#### GET RELATIVE ABUNDANCES #####
# Get the total ASV abundance within each meter
ASVmeterTotal <- as.data.frame(matrix(nrow=10, ncol=1))
rownames(ASVmeterTotal) <- paste(meterVec, "m", sep="_") 
colnames(ASVmeterTotal) <- "ASVmeterTotal"
for (j in 1:length(meterVec)){
  ASVmeterTotal[j,1] <- sum(fungiMeanASVsByMeterHabitat[which(fungiMeanASVsByMeterHabitat$meter ==  rownames(ASVmeterTotal)[j]),3]) #pull out meter by meter
}
ASVmeterTotal

# Finally get relative abundances by dividing each value in each row by the meter total
fungiRelAbundDfs <- vector("list", length(meterVec)) #this will have all of the rownames that correspond to samples in each meter
names(fungiRelAbundDfs) <- rownames(ASVmeterTotal)
for (h in 1:length(fungiRelAbundDfs)){
  fungiRelAbundDfs[[h]] <- fungiMeanASVsByMeterHabitat[which(fungiMeanASVsByMeterHabitat$meter ==  names(fungiRelAbundDfs)[[h]]),] #this pulls out just the data for each meter
  fungiRelAbundDfs[[h]]$relAbund <- (fungiRelAbundDfs[[h]]$meanASVabundance/ASVmeterTotal[h,1])*100
}
sum(fungiRelAbundDfs[[1]]$relAbund) #great, these add up to be 100!
sum(fungiRelAbundDfs[[2]]$relAbund)

# Checks:
# ASV 1, 20 meters
twentyMeterCheckData <- fungiMeanASVsByMeterHabitat[which(fungiMeanASVsByMeterHabitat$meter ==  names(fungiRelAbundDfs)[[2]]),]
twentyMeterCheckData[1,3] == fungiRelAbundDfs[[2]][1,3] #these are equal so this is working well so far
(twentyMeterCheckData[1,3]/ASVmeterTotal[2,1])*100 == fungiRelAbundDfs[[2]][1,11] #TRUE

# ASV 99, 100 meters
hundredMeterCheckData <- fungiMeanASVsByMeterHabitat[which(fungiMeanASVsByMeterHabitat$meter ==  names(fungiRelAbundDfs)[[10]]),]
which(hundredMeterCheckData$ASV_name == "ASV_99") #411th row
hundredMeterCheckData[411,3] == fungiRelAbundDfs[[10]][411,3] #these are equal
(hundredMeterCheckData[411,3]/ASVmeterTotal[10,1])*100 == fungiRelAbundDfs[[10]][411,11] #TRUE

### rbind all of the relAbunds[[all meters]] and plug this into ggplot2 below
fungiRelAbund_df <- do.call("rbind", fungiRelAbundDfs)
# View(fungiRelAbund_df)

# Plot it!
# Relative abundance
level_order <- names(fungiMeterRowNamesIndices) #set this to make in correct order from 10m to 100m
transectX <- c("40 m", "30 m","20 m", "10 m", "edge", "10 m", "20 m", "30 m", "40 m", "50 m")#for re-naming axis tick marks
relAbundTransectFungal_plot <- ggplot(fungiRelAbund_df, aes(x = factor(meter, level = level_order), y = relAbund, fill = diffAbundHabitat)) + 
  geom_bar(stat = "identity", position = "fill")  +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= NULL, labels=c("non-specialists", "forest specialists", "patch specialists")) +
  labs(y= "relative abundance", x = "meters from edge") + 
  ggtitle("Fungal Specialists Across the Transect") +
  scale_x_discrete(labels=transectX) + #change x-axis tick labels
  theme(axis.title=element_text(size=14))
# quartz()
relAbundTransectFungal_plot

# Saved Dec. 11, 2022
# save(relAbundTransectFungal_plot, file= "RobjectsSaved/relAbundTransectFungal_plot")

# Bring in relAbundTransectProks_plot from DiffAbundProkaryotes.R so that the plots can be plotted together
load("RobjectsSaved/relAbundTransectProks_plot")

# Plot prokaryotic and fungal plots together
require(gridExtra)
quartz()
grid.arrange(relAbundTransectProks_plot, relAbundTransectFungal_plot, ncol=1)

# Make another version without titles
relAbundTransectFungal_plotnoTitle <- relAbundTransectFungal_plot + ggtitle(element_blank()) #remove title
relAbundTransectProk_plotnoTitle <- relAbundTransectProks_plot + ggtitle(element_blank()) # remove title
grid.arrange(relAbundTransectProk_plotnoTitle, relAbundTransectFungal_plotnoTitle, ncol=1)

##########################################################################
# 4. COMPARISONS AND PLOTS OF UBIQUITY IN NON-SPECIALISTS VERSUS SPECIALISTS 
##########################################################################
colnames(fungiASVsAllDiffAbund_tax)
fungiASVsAllDiffAbund_tax$sampOccurrence <- rowSums(fungiASVsAllDiffAbund_tax[,2:234] > 0) #just select samples and get occurrence across samples
# View(fungiASVsAllDiffAbund_tax)

# Plot the number of occurrences in each group:
# ggplot it!
fungiUbiquityBySpecialistPlot <- ggplot(fungiASVsAllDiffAbund_tax, aes(x=diffAbundHabitat, y=sampOccurrence, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("(Ubiquity) number of sites") +
  xlab("Habitat Specialist Group") +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))

quartz()
fungiUbiquityBySpecialistPlot

# WITHOUT 100 M SAMPLES since there are more samples in the forest than the patch
fungiASVsAllDiffAbund_tax_no100m <- fungiASVsAllDiffAbund_tax
# find samples taken at 100 m
index100m <- which(sample_data(ITS_postUbiquity.ps)$Meter == 100) 
samples100m <- rownames(sample_data(ITS_postUbiquity.ps))[index100m]
samples100mCols <- which(colnames(fungiASVsAllDiffAbund_tax) %in% samples100m) #find the column names corresponding to these 100m samples
fungiASVsAllDiffAbund_tax_no100m <- fungiASVsAllDiffAbund_tax[,-samples100mCols] #make a new object without these columns
colnames(fungiASVsAllDiffAbund_tax_no100m)
fungiASVsAllDiffAbund_tax_no100m$sampOccurrence <- rowSums(fungiASVsAllDiffAbund_tax_no100m[,2:210] > 0) #just select samples and get occurrence across samples
# Note: lowest uumber of samples found in here is 32 samples (not 40) because 100 m samples were removed
# View(fungiASVsAllDiffAbund_tax_no100m)

# Test to see if I can do a one-way ANOVA
fungUbiq.anova <- aov(sampOccurrence ~ diffAbundHabitat, data = fungiASVsAllDiffAbund_tax_no100m)
car::leveneTest(sampOccurrence ~ diffAbundHabitat, data = fungiASVsAllDiffAbund_tax_no100m) # shows that variances are actually somewhat homogeneous
fungUbiq.anova.residual <- residuals(fungUbiq.anova) 
shapiro.test(fungUbiq.anova.residual) #shows that the residuals are NOT normal though
# Since assumptions aren't met, try pairwise Wilcox with Benjamini-Hochberg correction
fungiUbiquity.wilcox <- pairwise.wilcox.test(fungiASVsAllDiffAbund_tax_no100m$sampOccurrence, fungiASVsAllDiffAbund_tax_no100m$diffAbundHabitat,
                                            p.adjust.method = "BH")
fungiUbiquity.wilcox #shows patch > non-specialistand forest & > non-specialist

# ggplot it!
fungiUbiquityBySpecialistPlot_no100m <- ggplot(fungiASVsAllDiffAbund_tax_no100m, aes(x=diffAbundHabitat, y=sampOccurrence, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Ubiquity (number of sites)") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Fungi: ubiquity of ASVs across specialist groups") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

#quartz()
fungiUbiquityBySpecialistPlot_no100m

# Saved Dec 11, 2022
# save(fungiUbiquityBySpecialistPlot_no100m, file= "RobjectsSaved/fungiUbiquityBySpecialistPlot_no100m")

# Plot fungi and prokaryote plots side by side 
load(file= "RobjectsSaved/prokUbiquityBySpecialistPlot_no100m")
quartz()
grid.arrange(prokUbiquityBySpecialistPlot_no100m, fungiUbiquityBySpecialistPlot_no100m, ncol=2)

##########################################################################
# 5. MEDIAN ABUNDANCE OF EACH ASV ACROSS SAMPLES
##########################################################################
colnames(fungiASVsAllDiffAbund_tax)
# Get median abundance across all samples
fungiASVsAllDiffAbund_tax$ASVmedAbund <- rowMedians(as.matrix(fungiASVsAllDiffAbund_tax[,2:234]))

# Make plot for median abundance 
fungiMedASVabundPlot <- ggplot(fungiASVsAllDiffAbund_tax, aes(x=diffAbundHabitat, y=ASVmedAbund, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Median ASV abundance") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Fungi: Median ASV abundance across samples") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

# quartz()
fungiMedASVabundPlot
# save(fungiMedASVabundPlot, file= "RobjectsSaved/fungiMedASVabundPlot") #last saved Dec. 12, 2022

# MAKE PLOT OF MEDIAN ABUNDANCE IN PATCH OR FOREST SAMPLES
# The plot above shows that this is not super informative, so get it across samples, in the patch, on edge, and in forest
unique(colnames(fungiASVsAllDiffAbund_tax[2:234]) == rownames(fungimetaDf)) #this is true, which shows that we can use indices from
# metadata (plus 1!) to isolate stuff from ASV table
fungiPatchSampleIndex <- which(fungimetaDf$Meter %in% c(10,20,30,40)) #isolate patch samples
fungiForestSampleIndex <- which(fungimetaDf$Meter %in% c(60,70,80,90,100)) #isolate forest samples
fungiASVsAllDiffAbund_tax$patchASVmedAbund <- rowMedians(as.matrix(fungiASVsAllDiffAbund_tax[,fungiPatchSampleIndex+1])) #add one because columns start at 2
fungiASVsAllDiffAbund_tax$forestASVmedAbund <- rowMedians(as.matrix(fungiASVsAllDiffAbund_tax[,fungiForestSampleIndex+1])) #add one because columns start at 2
# Finally, make a new "longer" and "tidy" data frame
fungiASVsAllDiffAbund_taxAbunGroups <- fungiASVsAllDiffAbund_tax %>% 
  pivot_longer(cols="ASVmedAbund":"forestASVmedAbund", names_to= "habitatForMedAbundCalc", values_to="medAbund")
# View(fungiASVsAllDiffAbund_taxAbunGroups)

# grouped boxplot
fungiMedAbundByHabitatPlot <- ggplot(fungiASVsAllDiffAbund_taxAbunGroups, aes(x=diffAbundHabitat, y=medAbund, color=habitatForMedAbundCalc)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  ylab("Median ASV abundance") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Fungi: Median ASV abundance across samples") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

# quartz()
fungiMedAbundByHabitatPlot
# save(fungiMedAbundByHabitatPlot, file= "RobjectsSaved/fungiMedAbundByHabitatPlot") #last saved Dec. 12, 2022

##########
# MEAN ABUNDANCE 
##########
# This plot above is not very informative. Will take mean instead and then log 10 transform (mean so that there
# are values of zero
# Get mean abundance across all samples
colnames(fungiASVsAllDiffAbund_tax)
fungiASVsAllDiffAbund_tax$ASVmedAbund <- rowMedians(as.matrix(fungiASVsAllDiffAbund_tax[,2:234]))
fungiASVsAllDiffAbund_tax$ASVmeanAbund <- rowMeans(as.matrix(fungiASVsAllDiffAbund_tax[,2:234]))

# Make plot for mean abundance 
fungiMeanASVabundPlot <- ggplot(fungiASVsAllDiffAbund_tax, aes(x=diffAbundHabitat, y=ASVmeanAbund, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Mean ASV abundance") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Fungi: Mean ASV abundance across samples") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

# quartz()
fungiMeanASVabundPlot

log10fungi <- fungiMeanASVabundPlot  + scale_y_continuous(trans='log10') +
  ylab("log10 of Mean ASV abundance")
quartz()
log10fungi

log10(280)

##########
# MEAN ABUNDANCE LOG TRANSFORMED
##########

# Log transform ASV mean abundance values
fungiASVsAllDiffAbund_taxlog10 <- fungiASVsAllDiffAbund_tax #make a duplicate
fungiASVsAllDiffAbund_taxlog10$ASVmeanAbundLog10 <- log10(fungiASVsAllDiffAbund_taxlog10$ASVmeanAbund)
# View(fungiASVsAllDiffAbund_taxlog10)

# Make plot for these log transformed means
fungiLog10MeansASVabundPlot <- ggplot(fungiASVsAllDiffAbund_taxlog10 , aes(x=diffAbundHabitat, y=ASVmeanAbundLog10, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  geom_jitter(color="black", size=0.5, alpha=0.9) +
  ylab("Log10 of Mean ASV abundance") +
  xlab("Habitat specialist group") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  ggtitle("Fungi: Mean ASV abundance across samples") +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank())

# quartz()
fungiLog10MeansASVabundPlot

# save(fungiLog10MeansASVabundPlot, file= "RobjectsSaved/fungiLog10MeansASVabundPlot") #last saved Dec. 12, 2022


#######################################################################################################
# PART 2: EU 53N EXCLUDED FROM ANALYSES-- 
#######################################################################################################

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
postUbiqITSno53N.ps <- subset_samples(ITS_postUbiquity.ps, EU != "EU_53N")
unique(sample_data(postUbiqITSno53N.ps)$EU) #this shows that EU 53N is no longer present!
which(rowSums(otu_table(postUbiqITSno53N.ps))==0) #this shows that we don't drop any ASVs by getting dropping EU 53N.
# In other words, there were no ASVs that were ONLY found in EU 53N, so we are good to move to the next step

# Get just the ASVs out of postUbiqITSno53N.ps
postUbiqITSno53N_ASVs <- ASVs_outta_ps(postUbiqITSno53N.ps)
postUbiqITSno53N_ASVs_Plus1 <- t(postUbiqITSno53N_ASVs + 1) #adding one to every abundance count for differential abundance analysis;
# see explanation a few lines down. Invert so that I can make into a new phyloseq object.
# Make a new phyloseq object with these above
OTU = otu_table(postUbiqITSno53N_ASVs_Plus1, taxa_are_rows = TRUE)
postUbiqITSno53N_Plus1.ps <- phyloseq(OTU, tax_table(postUbiqITSno53N.ps), sample_data(postUbiqITSno53N.ps))
# this new object has 413 taxa across 193 samples

# Remove edge samples to compare patch and matrix
NoEdgeNo53N_ITS_Plus1.ps <- subset_samples(postUbiqITSno53N_Plus1.ps, Habitat != "edge") 
NoEdgeNo53N_ITS_Plus1.ps #413, 173
unique(sample_data(NoEdgeNo53N_ITS_Plus1.ps)$EU) #no EU 53N
unique(sample_data(NoEdgeNo53N_ITS_Plus1.ps)$Habitat) #no edge
step1_ITS_no53N <- phyloseq::phyloseq_to_deseq2(NoEdgeNo53N_ITS_Plus1.ps, ~ Habitat) #set up DESeq2 dds object 
step2_ITS_no53N <- DESeq2::DESeq(step1_ITS_no53N, test="Wald", fitType = "parametric") #differential expression analysis step;
#uses default Benjamini-Hochberg correction 
DeseqResults_ITS_no53N <- results(step2_ITS_no53N, cooksCutoff = FALSE) #make results object
DeseqResults_ITS_no53N <- DeseqResults_ITS_no53N[which(DeseqResults_ITS_no53N$padj < 0.001), ] #only get those ASVs below alpha level
DeseqResults_ITS_no53N <- cbind(as(DeseqResults_ITS_no53N, "data.frame"), as(tax_table(NoEdgeNo53N_ITS_Plus1.ps)[rownames(DeseqResults_ITS_no53N), ], "matrix")) #clean up format
DeseqResults_ITS_no53N$Habitat <- ifelse(DeseqResults_ITS_no53N$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for

# Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
# View(DeseqResults_ITS_no53N) #279 diff abund ASVs) 
sampDat <- sample_data(postUbiqITSno53N_Plus1.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
### CHANGE dataframe, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ITS_ASVnamesDA_no53N <- rownames(DeseqResults_ITS_no53N) 
length(ITS_ASVnamesDA_no53N) #279 found, compared to 287 as without EU 53N
ITS_ASVtabDA_no53N <- as.data.frame(t(ASVs_outta_ps(postUbiqITSno53N_Plus1.ps))) #get ASV table from original and transpose so that ASVs are rows 
dim(ITS_ASVtabDA_no53N) #193 samples, as wanted!
# Create dataframe with everything of interest
ITS_diffAbunDat_no53N <- merge(DeseqResults_ITS_no53N, ITS_ASVtabDA_no53N, by= "row.names") #grab ASV tab info from only those sa that are differentially abundant
str(ITS_diffAbunDat_no53N) #279 ASVs, 193 samples! 
# View(ITS_diffAbunDat_no53N)
diffAbunDat_no53N_tidy_FUNGI <- ITS_diffAbunDat_no53N %>% 
  pivot_longer(cols= 16:ncol(ITS_diffAbunDat_no53N), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_no53N_tidy_FUNGI)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
dim(diffAbunDat_no53N_tidy_FUNGI)
# View(diffAbunDat_no53N_tidy_FUNGI) this has 53,847 rows, which is equal to 193 (number of samples) x 279 (number of diff abund ASVs)

#### EXPLORING HOW THE PRESENCE OF ABSENCE OF EU 53N MATTERS ####
length(intersect(ITS_ASVnamesDA_no53N, ASVnamesDA_FUNGI)) #this shows that 266 ASVs are shared, out of 279 if 53N is excluded and out of 287 if it is included
with53N_index <- (ASVnamesDA_FUNGI %in% ITS_ASVnamesDA_no53N) #looks across all the names of diff abundant ASVs in ASVnamesDA_FUNGI and,
# along this vector, says whether or not that ASV is found in ITS_ASVnamesDA_no53N
with53N_indexF <- which(with53N_index==FALSE) #those ASVs that are in original analysis, but NOT when 53N is excluded
onlyIfWith53N <- ASVnamesDA_FUNGI[with53N_indexF] #this shows those ASVs that are only differentially abundant if 53N is included

without53N_index <- (ITS_ASVnamesDA_no53N %in% ASVnamesDA_FUNGI) #shows where, in ITS_ASVnamesDA_no53N, the ASVs are the same as those in ASVnamesDA_FUNGI
without53N_index_F <- which(without53N_index==FALSE) ##those ASVs that are only present when 53N is excluded (13 ASVs; and 13+266 =279, the # of ASVs found)
onlyIf53Nexcluded <- ITS_ASVnamesDA_no53N[without53N_index_F] #this shows those ASVs that are only differentially abundant if 53N is EXCLUDED
intersect(onlyIf53Nexcluded, onlyIfWith53N) #shows that these two objects are distinct, so one line of evidence showing I did it right

# save(diffAbunDat_no53N_tidy_FUNGI, file= "RobjectsSaved/diffAbunDat_no53N_tidy_FUNGI") #last saved Oct. 4, 2022
##########

# Re-make phyloseq object with just the differentially abundant taxa for ease of working with it in the future!
diffAbund_no53N_FUNGINames <- as.character(unique(diffAbunDat_no53N_tidy_FUNGI$ASV_name)) #added as character because at first it had "AsIs" attribute
fungi_no53N_DiffAbund.ps <- phyloseq::prune_taxa(diffAbund_no53N_FUNGINames, postUbiqITSno53N.ps)
dim(otu_table(fungi_no53N_DiffAbund.ps)) #279 ASVs ASVs across 193 samples
sample_data(fungi_no53N_DiffAbund.ps) #edges are there!

# save(fungi_no53N_DiffAbund.ps, file= "RobjectsSaved/fungi_no53N_DiffAbund.ps") #saved Oct. 4, 2022

##########################################################################
# II. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################
fungi_DeseqResults_no53NMini <- DeseqResults_ITS_no53N[,c(8,14)]  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
fungi_postUbiqTaxTabPLus1_no53N <- taxtable_outta_ps(postUbiqITSno53N_Plus1.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(fungi_postUbiqTaxTabPLus1_no53N) %in% rownames(fungi_DeseqResults_no53NMini) == TRUE)) #279! So these are the same 
length(which(rownames(fungi_postUbiqTaxTabPLus1_no53N) %in% rownames(fungi_DeseqResults_no53NMini) == FALSE)) #134 ASVs out of the 413 are not differentially abundant
false_index_no53N <- which(rownames(fungi_postUbiqTaxTabPLus1_no53N) %in% rownames(fungi_DeseqResults_no53NMini) == FALSE) #also 134

# Now, construct a dataframe with the taxa that were NOT differentially abundant 
notDA_taxTab_no53N_fungi <- fungi_postUbiqTaxTabPLus1_no53N[false_index_no53N,] #get a taxonomy tab with ONLY the non-differentially abundant ASVs
notDA_taxTab_no53N_fungi$Habitat <- "AremainingASVs" #make a habitat column that labels these as NOT differentially abundant. A in front so that would be first
# in ggplot for ease.
# View(notDA_taxTab_no53N_fungi)
colnames(notDA_taxTab_no53N_fungi)
notDA_taxTabMini_no53N_fungi <- notDA_taxTab_no53N_fungi[,c(2,8)] #keep only phylum and habitat to match DeseqResults_no53NMini
fungi_DAphylumAll_no53N <- rbind(fungi_DeseqResults_no53NMini, notDA_taxTabMini_no53N_fungi) #this has ASV name, phylum, and whether diff abundant for ALL ASVs in postUbiquity analysis
# What are the phyla breakdown here?
# so effectively, we want to get numbers in each phyla in each group. Should just be able to plot this?

diffAbund__no53N_ITS_stackedBarplotPhyla <- ggplot(fungi_DAphylumAll_no53N, aes(fill=Habitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  scale_x_discrete(labels=c("p__Ascomycota" = "Ascomycota", "p__Basidiomycota" = "Basidiomycota",
                            "p__Calcarisporiellomycota" = "Calcarisporiellomycota",
                            "p__Glomeromycota" = "Glomeromycota", "p__Mortierellomycota"="Mortierellomycota",
                            "p__Mucoromycota" = "Mucoromycota", "p__Olpidiomycota"= "Olpidiomycota",
                            "p__Rozellomycota" = "Rozellomycota")) + 
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,400,by=25)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant Fungal ASVs (53N excluded)") 
# quartz()
diffAbund__no53N_ITS_stackedBarplotPhyla


# Below saved Oct 4, 2022 so that it can be added to a 2 paneled plot with fungal plot!
# save(diffAbund__no53N_ITS_stackedBarplotPhyla, file="RobjectsSaved/diffAbund__no53N_ITS_stackedBarplotPhyla")

# MORE ABOUT DIFFERENCES WITH AND WITHOUT 53N --How do just the forest taxa overlap?
with53NforestNames_fungi <- rownames(DAphylumAll[which(DAphylumAll$Habitat=="forest"),]) #119 ASVs
with53NpatchNames_fungi <- rownames(DAphylumAll[which(DAphylumAll$Habitat=="patch"),]) #168 ASVs
with53NnotDANames_fungi <- rownames(DAphylumAll[which(DAphylumAll$Habitat=="AremainingASVs"),]) #126 ASVs
without53NforestNames_fungi <- rownames(fungi_DAphylumAll_no53N[which(fungi_DAphylumAll_no53N$Habitat=="forest"),]) # 114 ASVs! 
without53NpatchNames_fungi <- rownames(fungi_DAphylumAll_no53N[which(fungi_DAphylumAll_no53N$Habitat=="patch"),]) #165 ASVs
without53NnotDANames_fungi <- rownames(fungi_DAphylumAll_no53N[which(fungi_DAphylumAll_no53N$Habitat=="AremainingASVs"),]) #134
length(intersect(with53NforestNames_fungi, without53NforestNames_fungi)) #110, so 9 aren't shared (shown below too as length(onlyIfWith53N_forest))

with53N_indexForest <- (with53NforestNames_fungi %in% without53NforestNames_fungi) #shows the with53Nforest names that do or do not occur in without53NforestNames
with53N_indexForest_F <- which(with53N_indexForest==FALSE) #these are the rows of those that occur with53N_indexForest but NOT in without53NforestNames
onlyIfWith53N_forest <- with53NforestNames_fungi[with53N_indexForest_F] #this shows those ASVs that are only differentially abundant if 53N is included
length(onlyIfWith53N_forest) #9 unique ASVs if EU 53N is included

withOUT53N_indexForest <- (without53NforestNames_fungi %in% with53NforestNames_fungi) #shows the without53NforestNames names that do or do not in with53Nforest 
withOUT53N_indexForest_F <- which(withOUT53N_indexForest==FALSE) #these are the rows of those that occur in without53NforestNames but NOT in  with53N_indexForest 
onlyIfWithout53N_forest <- without53NforestNames_fungi[withOUT53N_indexForest_F] #this shows those ASVs that are only differentially abundant if 53N is EXCLUDED
length(onlyIfWithout53N_forest) #4. This makes sense because 110 were shared with other condition out of 114 differentially abundant ASVs




