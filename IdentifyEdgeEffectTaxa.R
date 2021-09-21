## Identifying possible ASVs by EU/Site
# (started) September 14, 2021

# The purpose of this script is to identify ASVs that occur in several samples
# per EU, and that show a difference between forest and patch. That way, I can
# begin narrowing down ASVs and then lower resoltution taxonomic groups that may
# show an edge effect.

# This script uses some objects made in 16SExploratoryDataAnalysisAug2021.R, but 
# DOES NOT rely on anything created in 16SEUEDASept2021.R. Thus, if I end up keeping
# the pipeline created here, this script can be used in place of 16SEUEDASept2021.R.

# Functions defined in this script:
# taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# ASVsampOccurrence determines the number of samples that each ASV occurs in 
ASVsampOccurrence <- function(physeq) { #input is a phyloseq object
  OTU <- otu_table(physeq)
  ASVmat <- as(OTU, "matrix") # make phyloseq ASV table into a non-phyloseq matrix
  sampleOccurence <- rowSums(ASVmat > 0) #sum up number of samples ASV occurs in 
  return(cbind(ASVmat, sampleOccurence)) #
}

# This function gets the ASV table out of phyloseq 
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

###################################################################################################

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

index <- which(datmat > 0) #another way of doing it, R goes down columns
datmat[index] <- 1

# How many columns (samples) does each ASV (row) appear in?
ASVsampCount <- rowSums(datmat > 0)
ASVsampCount #looks correct!
cbind(datmat, ASVsampCount) 


#############
# EU 52
#############

ASVwithCount_52 <- ASVsampOccurrence(EU_52_Soils.ps)
# View(ASVwithCount_52) taxa/ASVs re rows, and samples are columns
dim(ASVwithCount_52)
quartz()
barplot(table(ASVwithCount_52[,39]), ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 38)", main= "EU_52" )
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
# View(ASVtab_52_4x_da_counts)

########################################################################
################ ASVs THAT OCCUR AT LEAST 10 TIMES ###### (i.e. about 1/4 of samples)
########################################################################

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

######### BY WHOLE EU/SITE #######
# Here, I do a relative abundance analysis for the WHOLE site. In other words, I combine samples based
# on where they fall along the transect. In the next section, I leave the samples as individual samples
# and then look by transect. 

# Combine samples by meter (values across ASV table is SUM across all samples at the same meter point,
# and the rest of the metadata variables are the mean for each position along the meter)
EU_52_10_ByMeter.ps <- merge_samples(EU_52_10timesNoEdge.ps, group= "Meter")
EU_52_10_ByMeterASV <- as.data.frame(t(ASVs_outta_ps(EU_52_10_ByMeter))) #ASVs are rows, meter/samples are columns
dim(EU_52_10_ByMeterASV) #ASVs are rows, meter/samples are columns
# View(EU_52_10_ByMeterASV)
class(EU_52_10_ByMeterASV)

# Add back in habitat information (merging samples above finds mean of the metadata variables,
# and since habitat is categorical it goes to NAs)
sample_data(EU_52_10_ByMeter.ps)$Habitat <- c(rep("Patch", 4), rep("Forest",5))
sample_data(EU_52_10_ByMeter.ps)$Habitat 

# DIFFERENTIAL ABUNDANCE ANALYSIS
Deseq1_52_10x_ByMeter <- phyloseq_to_deseq2(EU_52_10_ByMeter.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
Deseqtested_52_10x_ByMeter <- DESeq(Deseq1_52_10x_ByMeter, test="Wald", fitType = "parametric")
DeSeq_res_52_10x_ByMeter <- results(Deseqtested_52_10x_ByMeter, cooksCutoff = FALSE)
alpha <- 0.01
sigtab_52_10x_ByMeter <- DeSeq_res_52_10x_ByMeter[which(DeSeq_res_52_10x_ByMeter$padj < alpha), ]
sigtab_52_10x_ByMeter <- cbind(as(sigtab_52_10x_ByMeter, "data.frame"), as(tax_table(EU_52_10_ByMeter.ps)[rownames(sigtab_52_10x_ByMeter), ], "matrix"))
head(sigtab_52_10x_ByMeter)
dim(sigtab_52_10x_ByMeter) #132 ASVs out of the 7514 tested had a corrected p-value of less than 0.01
#View(sigtab_52_10x_ByMeter)

# How many samples are each of these ASVs found in and what is their overall abundance?
names_52_10x_ByMeter <- rownames(sigtab_52_10x_ByMeter) #pull out names of the ASVs
ASVtab_52_10x_ByMeter_da <- EU_52_10_ByMeterASV[names_52_10x_ByMeter,] #get a smaller version of the ASV table that has only these ASVs
SampleCount_52_10x_ByMeter_da <- rowSums(ASVtab_52_10x_ByMeter_da > 0) #get number of soil samples that the ASV appears in
ASVtab_52_10x_ByMeter_da_counts <- cbind(ASVtab_52_10x_ByMeter_da, SampleCount_52_10x_ByMeter_da)
ASVabund <- rowSums(ASVtab_52_10x_ByMeter_da_counts[,1:9]) #get abundance of each ASV ACROSS all samples
ASVtab_52_10x_ByMeter_da_counts <- cbind(ASVtab_52_10x_ByMeter_da_counts, ASVabund)
#View(ASVtab_52_10x_ByMeter_da_counts)

# What do these abundances look like in the forest versus the patch?
# Use the forest and patch indices created above (when looking at 4x) to get average abundance in each
forest_sum <- rowSums(ASVtab_52_10x_ByMeter_da_counts[,5:9]) #gets sum across each row of ONLY the forest samples
forest_meanAbund <- forest_sum/5 #average number in each meter in the forest
patch_sum <- rowSums(ASVtab_52_10x_da_counts[,1:4])
patch_meanAbund <- patch_sum/4

# Combine it all together-- NON RELATIVE ABUNDANCE-- this is kinda like all the info!
ASVtab_52_10x_ByMeter_da_counts <- cbind(ASVtab_52_10x_ByMeter_da_counts, forest_meanAbund, patch_meanAbund)
head(ASVtab_52_10x_ByMeter_da_counts)
#View(ASVtab_52_10x_da_counts)

# Only consider ASVs that appear at least 50 times across whole matrix
above_50 <- which(ASVtab_52_10x_ByMeter_da_counts[,11] >= 50)
ASVtab_52_10x_ByMeter_50 <- ASVtab_52_10x_ByMeter_da_counts[above_50,]
# View(ASVtab_52_10x_ByMeter_50)

# Convert the counts (only ASVs that occur at least 50 times) to relative abundance
rowmax <- rep(NA, nrow(ASVtab_52_10x_ByMeter_50)) #make an empty vector to store rowmaxes
for (i in 1:nrow(ASVtab_52_10x_ByMeter_50)) {
  rowmax[i] <- max(ASVtab_52_10x_ByMeter_50[i, 1:9]) #get maximum abundance of each ASV
  relAbund_52_10x_ByMeter <- ASVtab_52_10x_ByMeter_50[,1:9]/rowmax
}
relAbund_52_10x_ByMeter
#View(relAbund_52_10x_ByMeter) # I manually calculated these and it looks good!

# Plot with Base R
relAbund_52_10x_ByMeter <- t(relAbund_52_10x_ByMeter) #now samples are x and ASVs are y
quartz()
relAbund_52_10x_ByMeter_plot <- matplot(rownames(relAbund_52_10x_ByMeter), relAbund_52_10x_ByMeter, type = "l", xlab= "Meter",
                     ylab= "Relative Abundance", main= "Changes in ASV abundance across EU 52- Meters Combined")

# Plot with ggplot and color by taxonomic group???
relAbund_52_10x_ByMeter_gg <- t(relAbund_52_10x_ByMeter) #switch it around again
ASVindex <- rownames(relAbund_52_10x_ByMeter_gg) #get all of these sample names 
taxTable52NoEdge <- taxtable_outta_ps(EU_52_10_ByMeter.ps) #pull out tax table
taxTable52NoEdge_10x_ByMeter <- taxTable52NoEdge[ASVindex,1:6] #get taxonomic information
dim(taxTable52NoEdge_10x_ByMeter) #108 ASVs 
relAbund_52_10x_ByMeter_gg <- t(merge(relAbund_52_10x_ByMeter_gg, taxTable52NoEdge_10x_ByMeter, by=0)) #by=0 makes it combine by rownames
# View(relAbund_52_10x_ByMeter_gg) #great, this looks as expected!
header.true <- function(df) { # from Pierre L: https://stackoverflow.com/questions/32054368/use-first-row-data-as-column-names-in-r
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
relAbund_52_10x_ByMeter_gg <- header.true(relAbund_52_10x_ByMeter_gg) #lose ASV names but that's okay
colnames(relAbund_52_10x_ByMeter_gg) <- rownames(taxTable52NoEdge_10x_ByMeter) #add them back in
relAbund_52_10x_ByMeter_gg <- as.data.frame(relAbund_52_10x_ByMeter_gg)
tax <- t(relAbund_52_10x_ByMeter_gg[10:15,]) #grab and save taxonomic information and invert it so I can make it into columns
relAbund_52_10x_ByMeter_gg <- relAbund_52_10x_ByMeter_gg[1:9,] #remove last rows that are taxonomic info
# help with ggplot from here: https://community.rstudio.com/t/how-plot-all-values-inside-a-data-frame-into-a-graph-using-ggplot-function/75021
relAbund_52_10x_ByMeter_gg <- rownames_to_column(relAbund_52_10x_ByMeter_gg, var="Meter")
relAbund_52_10x_ByMeter_gg <- relAbund_52_10x_ByMeter_gg %>% pivot_longer(cols= ASV_7:ASV_8741,
                              names_to= "ASV_name", values_to = "ASV_Rel_Abundance")
head(relAbund_52_10x_ByMeter_gg) #now each ASV in each sample has a value
tail(relAbund_52_10x_ByMeter_gg)
relAbund_52_10x_ByMeter_gg <- cbind.data.frame(relAbund_52_10x_ByMeter_gg, tax) #tax is now columns
# Make it so the meters are plotted in the correct order
relAbund_52_10x_ByMeter_gg$Meter <- factor(relAbund_52_10x_ByMeter_gg$Meter, levels = c("10", "20", "30", "40", "60", "70", "80", "90", "100"))
relAbund_52_10x_ByMeter_gg$ASV_Abundance <- as.numeric(relAbund_52_10x_ByMeter_gg$ASV_Abundance)
quartz()
ggplot(relAbund_52_10x_ByMeter_gg, aes(x=Meter, 
        y=ASV_Abundance), color= Phylum, group=ASV_name) +
    geom_point() +
    geom_line() +
    theme_classic()

quartz()
ggplot(relAbund_52_10x_ByMeter_gg, aes(Meter, ASV_Rel_Abundance, color = Phylum, group = ASV_name)) + 
  geom_line() + theme(axis.text.y = element_blank()) + ggtitle("Changes in ASV abundance across EU 52 (combined by meter)")




######### BY TRANSECT #########

# DIFFERENTIAL ABUNDANCE ANALYSIS
Deseq1_52_10x <- phyloseq_to_deseq2(EU_52_10timesNoEdge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
Deseqtested_52_10x <- DESeq(Deseq1_52_10x, test="Wald", fitType = "parametric")
DeSeq_res_52_10x <- results(Deseqtested_52_10x, cooksCutoff = FALSE)
alpha <- 0.01
sigtab_52_10x <- DeSeq_res_52_10x[which(DeSeq_res_52_10x$padj < alpha), ]
sigtab_52_10x <- cbind(as(sigtab_52_10x, "data.frame"), as(tax_table(EU_52_10timesNoEdge.ps)[rownames(sigtab_52_10x), ], "matrix"))
head(sigtab_52_10x)
dim(sigtab_52_10x) #132 ASVs out of the 7514 tested had a corrected p-value of less than 0.01

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
ASVtab_52_10x_T <- t(ASVtab_52_10x_T[,1:8])
rownames(ASVtab_52_10x_T) #mssing 70m
ASVtab_52_10x_T <- ASVtab_52_10x_T[c("10", "20", "30", "40", "60", "80", "90", "100"),] #rorder to make 10m - 100m
View(ASVtab_52_10x_T)

ASVtab_52_10x_B <- ASVtab_52_10x_da_50[,c(EU_52_Bindex,35:38)] 
colnames(ASVtab_52_10x_B)[1:8] <- sample_data(EU_52_10timesNoEdge.ps)$Meter[EU_52_Bindex]
rownames(ASVtab_52_10x_B) #missing 90m
ASVtab_52_10x_B <- t(ASVtab_52_10x_B[,1:8])
ASVtab_52_10x_B <- ASVtab_52_10x_B[c("10", "20", "30", "40", "60", "70", "80", "100"),] #rorder to make 10m - 100m


ASVtab_52_10x_L <- ASVtab_52_10x_da_50[,c(EU_52_Lindex,35:38)] 
colnames(ASVtab_52_10x_L)[1:9] <- sample_data(EU_52_10timesNoEdge.ps)$Meter[EU_52_Lindex]
ASVtab_52_10x_L <- t(ASVtab_52_10x_L[,1:9])
ASVtab_52_10x_L <- ASVtab_52_10x_L[c("10", "20", "30", "40", "60", "70", "80", "90", "100"),] #re-order to make 10m - 100m

ASVtab_52_10x_R <- ASVtab_52_10x_da_50[,c(EU_52_Rindex,35:38)] 
colnames(ASVtab_52_10x_R)[1:9] <- sample_data(EU_52_10timesNoEdge.ps)$Meter[EU_52_Rindex]
ASVtab_52_10x_R <- t(ASVtab_52_10x_R[,1:9])
ASVtab_52_10x_R <- ASVtab_52_10x_R[c("10", "20", "30", "40", "60", "70", "80", "90", "100"),] #re-order to make 10m - 100m


quartz()
par(mfrow= c(2,2))
T_52_plot <- matplot(rownames(ASVtab_52_10x_T), ASVtab_52_10x_T, type = "l", xlab= "Meter",
                    ylab= "Abundance", main= "Changes in ASV abundance across EU 52 T Transect")
B_52_plot <- matplot(rownames(ASVtab_52_10x_B), ASVtab_52_10x_B, type = "l", xlab= "Meter",
                     ylab= "Abundance", main= "Changes in ASV abundance across EU 52 B Transect")
L_52_plot <- matplot(rownames(ASVtab_52_10x_L), ASVtab_52_10x_L, type = "l", xlab= "Meter",
                     ylab= "Abundance", main= "Changes in ASV abundance across EU 52 L Transect")
R_52_plot <- matplot(rownames(ASVtab_52_10x_R), ASVtab_52_10x_R, type = "l", xlab= "Meter",
                     ylab= "Abundance", main= "Changes in ASV abundance across EU 52 R Transect")
            

# COMPARING BY METER AND BY TRANSECT 
sort(rownames(sigtab_52_10x)) == sort(rownames(sigtab_52_10x_ByMeter))
# By meter or by transect PULLED the same number of ASVs, but only some overlap. What are these?
ASVs_52_10x_MAndT <- Reduce(intersect, list(rownames(sigtab_52_10x), rownames(sigtab_52_10x_ByMeter)))
# Taxonomic info for these:
EU_52_10timesNoEdgeTax <- taxtable_outta_ps(EU_52_10timesNoEdge.ps)
ASVs_52_10x_MAndT_Tax <- EU_52_10timesNoEdgeTax[ASVs_52_10x_MAndT,]

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
barplot(table(ASVwithCount_10[,40]), main= "EU 10", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)")
length(which(ASVwithCount_10[,39] >= 4)) #895
length(which(ASVwithCount_10[,39] >= 10)) #296
length(which(ASVwithCount_10[,39] >= 20)) #133
length(which(ASVwithCount_10[,39] >= 30)) #80
length(which(ASVwithCount_10[,39] >= 39)) #50


###########
# Plot all ASV Sample Occurrence Plots Together 
quartz()
par(mfrow= c(2,3))
barplot(table(ASVwithCount_52[,39]), ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 38)", main= "EU_52" )
barplot(table(ASVwithCount_53N[,40]), main= "EU 53N", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)" )
barplot(table(ASVwithCount_54S[,39]), main= "EU 54S", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 38)" )
barplot(table(ASVwithCount_8[,40]), main= "EU 8", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)" )
barplot(table(ASVwithCount_53S[,41]), main= "EU 53S", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 40)" )
barplot(table(ASVwithCount_10[,40]), main= "EU 10", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)")
