# IdentifyEdgeEffectTaxa5 (New Differential Abundance Analysis)
# February 17, 2022 (begun)

# This script identifies edge effect ASVs using differential abundance analyses
# (DESeq2). This script performs differential abundance analyses on each
# individual EU, as per Noah's suggestions on Feb. 16 (see pic of white board
# and new schematic in PowerPoint presentation file). THIS DOES NOT USE THE 
# MEDIAN ABUNDANCES PER EU

###################################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

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
library("stringr") #for grep-like tools for data manipulation with character strings in data

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45
# R OBJECT MADE IN: UbiquityMedianSetup.R

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
    ASVmean <- rowMeans(dat[i,]) #if you change rowMeans to mean, could use with matrices. Or could build in if/else statement
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

###########################################################################
# DIFFERENTIAL ABUNDANCE ANALYSIS ON EACH EU
###########################################################################
 # First, split up data based on EU
# Prune samples to separate by EU and check to make sure each looks right
EU_52_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU) #only EU 52
EU_53N_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_53N")
unique(sample_data(EU_53N_Soils.ps)$EU)
EU_54S_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_54S")
unique(sample_data(EU_54S_Soils.ps)$EU)
EU_8_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_8")
unique(sample_data(EU_8_Soils.ps)$EU)
EU_53S_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_53S")
unique(sample_data(EU_53S_Soils.ps)$EU)
EU_10_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_10")
unique(sample_data(EU_10_Soils.ps)$EU)

sample_data(postUbiquity.ps)

####### EU 8 #######

ASVwithCount_8 <- ASVsampOccurrence(EU_8_Soils.ps)
dim(ASVwithCount_8)
barplot(table(ASVwithCount_8[,40]), main= "EU 8", ylab="Number of ASVs", xlab= "Number of Samples Present in (out of 39)" )
length(which(ASVwithCount_8[,40] >= 4)) #7275
length(which(ASVwithCount_8[,40] >= 10)) #3481
length(which(ASVwithCount_8[,40] >= 15)) #1894
length(which(ASVwithCount_8[,40] >= 20)) #1072
length(which(ASVwithCount_8[,40] >= 30)) #349
length(which(ASVwithCount_8[,40] >= 39)) #52 found in all samples!

# Remove "edge" samples so that we can find samples that vary in edge versus forest in the differential abundance analysis
EU_8NoEdge.ps <- subset_samples(EU_8_Soils.ps, Habitat != "edge")

# DIFFERENTIAL ABUNDANCE ANALYSIS # changed alpha to 0.001 to reduce how many....
Deseq1_EU8 <- phyloseq_to_deseq2(EU_8NoEdge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
Deseqtested_EU8 <- DESeq(Deseq1_EU8, test="Wald", fitType = "parametric")
DeseqResults_EU8 <- results(Deseqtested_EU8, cooksCutoff = FALSE)
alpha <- 0.001
sigtab_DeseqResults_EU8 <- DeseqResults_EU8[which(DeseqResults_EU8$padj < alpha), ]
sigtab_DeseqResults_EU8 <- cbind(as(sigtab_DeseqResults_EU8, "data.frame"), as(tax_table(EU_8NoEdge.ps)[rownames(sigtab_DeseqResults_EU8), ], "matrix"))
head(sigtab_DeseqResults_EU8)
dim(sigtab_DeseqResults_EU8) #1205/4480 had pvalue less than 0.001; 1617/4480 had p-value less than 0.01
View(sigtab_DeseqResults_EU8)

### SEPARATE BASED ON BEING IN FOREST OR PATCH ###
sigtab_DeseqResults_EU8$habitat <- ifelse(sigtab_DeseqResults_EU8$log2FoldChange<0, "Forest", "Patch")
#View(sigtab_DeseqResults_EU8)
forestDA_ASVs_EU8 <- sigtab_DeseqResults_EU8[which(sigtab_DeseqResults_EU8$habitat=="Forest"),] #324 ASVs
patchDA_ASVs_EU8 <- sigtab_DeseqResults_EU8[which(sigtab_DeseqResults_EU8$habitat=="Patch"),] #881 patch ASVs

# GET ABUNDANCES FOR EACH ASV ACROSS TRANSECT, WITH FOREST (60M - 100m) COMBINED 
forestNames_Deseq_EU8 <- rownames(forestDA_ASVs_EU8) #pull out names of the ASVs
EU8_ASVs <- as.data.frame(t(ASVs_outta_ps(EU_8_Soils.ps))) #transpose so that ASVs are rows
EU8_da_forestASVtab <- EU8_ASVs[forestNames_Deseq_EU8,] #get a smaller version of the ASV table that has only the diff abundance ASVs

# PREPARE DATAFRAME WITH DIFF ABUNDANCE RESULTS, SOME SAMPLE INFO, AND TAXONOMIC INFORMATION
sampDat <-  sample_data(EU_8_Soils.ps)[,c(1,8:9)] #get variables that I want out of sample_data
attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute
class(sampDat) #double check
sampDat <- tibble::rownames_to_column(sampDat, var="sampleNumber") #make sampleNumber a column
#merge differential abundance info with the ASV table for the forest
EU8_diffAbunDat <- merge(forestDA_ASVs_EU8, EU8_da_forestASVtab, by= "row.names") 
EU8_diffAbunDat_tidy <- EU8_diffAbunDat %>% 
  pivot_longer(cols= "105":"96", names_to= "sampleNumber", values_to= "ASVabundance") %>% #now has 12,636 rows because 324 ASVs * 39 samples 
  merge(sampDat, by="sampleNumber") #merge EU8_diffAbunData with the sampDat to get meter and such for 
#View(EU8_diffAbunDat_tidy)





# How many samples are each of these ASVs found in and what is their overall abundance? (Add edge samples back in in this chunk)
names_Deseq_EU8 <- rownames(sigtab_DeseqResults_EU8) #pull out names of the ASVs
# Get ASV table (out of phyloseq) for all points along transect:

EU8_ASVs <- as.data.frame(t(ASVs_outta_ps(EU_8_Soils.ps))) #transpose so that ASVs are rows
Deseq_EU8_ASVtab <- EU8_ASVs[names_Deseq_EU8,] #get a smaller version of the ASV table that has only the diff abundance ASVs
SampleCount_Deseq_EU8_ASVtab <- rowSums(Deseq_EU8_ASVtab > 0) #get number of soil samples that the ASV appears in
ASVtab_8_da_counts <- cbind(Deseq_EU8_ASVtab, SampleCount_Deseq_EU8_ASVtab)
ASVabund <- rowSums(ASVtab_8_da_counts[,1:39]) #get abundance of each ASV ACROSS all samples, but make sure to leave out sample count columm
ASVtab_8_da_counts <- cbind(ASVtab_8_da_counts, ASVabund)
# View(ASVtab_8_da_counts)

################################################
##### RELATIVE ABUNDANCE STUFF ######
# (CODING NOT FINISHED)


# What do these abundances look like in the forest, patch, and edge?
forest_sum <- rowSums(ASVtab_8_15x_ByMeter_da_counts[,6:10]) #gets sum across each row of ONLY the forest samples
forest_meanAbund <- forest_sum/5 #average number in each meter in the forest
patch_sum <- rowSums(ASVtab_8_15x_ByMeter_da_counts[,1:4])
patch_meanAbund <- patch_sum/4
edge_meanAbund <- ASVtab_8_15x_ByMeter_da_counts[,5] # 50m, and it is still basically a mean since this is across multiple edge samples (up to four
# I don't remember which samples were dropped here)

# Combine it all together-- NON RELATIVE ABUNDANCE-- this is kinda like all the info!
ASVtab_8_15x_ByMeter_all <- cbind(ASVtab_8_15x_ByMeter_da_counts, forest_meanAbund, patch_meanAbund, edge_meanAbund)
head(ASVtab_8_15x_ByMeter_all)
dim(ASVtab_8_15x_ByMeter_all) #still 209!!!
#View(ASVtab_8_15x_ByMeter_all)

# Convert the counts to relative abundance
rowmax <- rep(NA, nrow(ASVtab_8_15x_ByMeter_da_counts)) #make an empty vector to store rowmaxes
for (i in 1:nrow(ASVtab_8_15x_ByMeter_da_counts)) {
  rowmax[i] <- max(ASVtab_8_15x_ByMeter_da_counts[i, 1:10]) #get maximum abundance of each ASV
  relAbund_8_15x_ByMeter <- ASVtab_8_15x_ByMeter_da_counts[,1:10]/rowmax
}
relAbund_8_15x_ByMeter
#View(relAbund_8_10x_ByMeter) # I manually calculated a few of these and it looks good!

# Plot with Base R
relAbund_8_15x_ByMeter <- t(relAbund_8_15x_ByMeter) #now samples are x and ASVs are y
# quartz()
relAbund_8_15x_ByMeter_plot <- matplot(rownames(relAbund_8_15x_ByMeter), relAbund_8_15x_ByMeter, type = "l", xlab= "Meter",
                                       ylab= "Relative Abundance", main= "Changes in ASV abundance across EU 8")
xtick<-seq(0, 100, by=10) 
axis(side=1, at=xtick, labels = TRUE) #change x axis ticks 

# Plot with ggplot and color by taxonomic group
relAbund_8_15x_ByMeter_gg <- t(relAbund_8_15x_ByMeter) #switch it around again (now is ASVs as rows and meters as columns )
ASVindex <- rownames(relAbund_8_15x_ByMeter_gg) #get all of these ASV names 
taxTable8_15_ByMeter <- taxtable_outta_ps(EU_8_15_ByMeter.ps) #pull out tax table
taxTable8_15x_ByMeter <- taxTable8_15_ByMeter[ASVindex,1:6] #get taxonomic information
dim(taxTable8_15x_ByMeter) #209 ASVs , six taxonomic types (through genus)
relAbund_8_15x_ByMeter_gg <- t(merge(relAbund_8_15x_ByMeter_gg, taxTable8_15x_ByMeter, by=0)) #by=0 makes it combine by rownames
# View(relAbund_8_15x_ByMeter_gg) #rows are now meters and taxonomic info, columns are ASVs
relAbund_8_15x_ByMeter_gg <- header.true(relAbund_8_15x_ByMeter_gg) #lose ASV names but that's okay
colnames(relAbund_8_15x_ByMeter_gg) <- rownames(taxTable8_15x_ByMeter) #add ASV names back in
relAbund_8_15x_ByMeter_gg <- as.data.frame(relAbund_8_15x_ByMeter_gg) #makes values numeric again!
tax <- t(relAbund_8_15x_ByMeter_gg[11:16,]) #grab and save taxonomic information and invert it so I can make it into columns
relAbund_8_15x_ByMeter_gg <- relAbund_8_15x_ByMeter_gg[1:10,] #remove last rows that are taxonomic info
# help with ggplot from here: https://community.rstudio.com/t/how-plot-all-values-inside-a-data-frame-into-a-graph-using-ggplot-function/75021
relAbund_8_15x_ByMeter_gg <- rownames_to_column(relAbund_8_15x_ByMeter_gg, var="Meter")
relAbund_8_15x_ByMeter_gg <- relAbund_8_15x_ByMeter_gg %>% pivot_longer(cols= ASV_1:ASV_9313,
                                                                        names_to= "ASV_name", values_to = "ASV_Rel_Abundance")
head(relAbund_8_15x_ByMeter_gg) #now each ASV in each sample has a value
tail(relAbund_8_15x_ByMeter_gg) #no more taxonomic info on the end
relAbund_8_15x_ByMeter_gg <- cbind.data.frame(relAbund_8_15x_ByMeter_gg, tax) #tax is now columns
# Make it so the meters are plotted in the correct order
relAbund_8_15x_ByMeter_gg$Meter <- factor(relAbund_8_15x_ByMeter_gg$Meter, levels = c("10", "20", "30", "40", "50","60", "70", "80", "90", "100"))
# Finally, plot this!

quartz()
ggplot(relAbund_8_15x_ByMeter_gg, aes(Meter, ASV_Rel_Abundance, color = Phylum, group = ASV_name)) + 
  geom_line() + theme(axis.text.y = element_blank()) + ggtitle("Changes in ASV abundance across EU 8")

###################################


