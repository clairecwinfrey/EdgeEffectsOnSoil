# I6S_indicAnalysisEdgePatchMatrix.R
# started November 1, 2023

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("utilities")
library("Polychrome") #for color palette
library("indicspecies")

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# This object was created and saved in UbiquityMedianSetup.R

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
# 1. Code for function (coped from link above)
# New function: indSpecTable2
##############################
# indSpecTable2 take the output of the indicator species analysis,
# the taxonomy table for the samples, and the ASV table for the 
# samples. It returns a dataframe that has:
# 1) the ASV TaxonID,
# 2) the taxonomic name of the ASV, 3) the biserial correlation coefficient
# ("stat"), 4) the p-value for the indicator ASV, 5) the
# group of samples that the ASV is an indicator for, 6) and the relative
# abundance of each indicator ASV in each group of samples!

indSpecTable2 <- function(myMultipatt, taxTable, ASVtable, signLevel=0.05) { 
  # Arguments are: 
  # 1) myMultipatt: an object produced by multipatt function;
  # 2) taxTable: a taxonomy table that has the taxonIDs as row names and 
  # 7 columns for taxonomic levels (Kingdom, Phylum, Class...Species);
  # 3) signLevel: significance level for p-value cut-off
  # 4) ASV table where samples/sites are rows and TaxonIDs are columns
  # (this will be inverted in script)
  require(indicspecies)
  groups <- unique(myMultipatt$cluster) #gets groups from myMultipatt
  # The lines below pull out only the indicator ASVs that have a
  # p-value less than or equal to the significance level specified in the 
  # function call. Indicator ASVs are represented by TaxonID
  ASVstoKeep <- rownames(myMultipatt$sign)[which(myMultipatt$sign$p.value <= signLevel)]
  ### replace ASVdfs with ASVsToKeep
  ASVtable_t <- t(ASVtable) #make rows TaxonID 
  ASVtableIndex <- which(rownames(ASVtable_t) %in% ASVstoKeep)
  ASV_sigs <- ASVtable_t[ASVtableIndex,] #new ASV table only has significant ASVs
  # Below makes it so that individual samples are labeled according to clustering group
  colnames(ASV_sigs) <- myMultipatt$cluster
  # Combine by clustering groups
  combByCluster <- t(rowsum(t(ASV_sigs), group=colnames(ASV_sigs), na.rm =T)) 
  # Below, pre-allocate a dataframe that we'll fill in with the relative abundances
  # of each ASV in each group
  relAbundDf <- data.frame(matrix(nrow=nrow(combByCluster), ncol=(length(groups)))) 
  # The for loop below makes column names for relAbundDf based on names of groups used for clustering
  relAbundCol <- rep(NA,length(groups)) #pre-allocate vector to hold relative abundance columns
  for (i in 1:length(groups)){ #make column names for the relative abundance of each group in analysis
    relAbundCol[i] <- paste("relAbund in", colnames(combByCluster)[i])
  }
  colnames(relAbundDf) <- relAbundCol
  totalCounts <- colSums(combByCluster) #get total ASV counts in each clustering group
  # For loop calculates the relative abundance of each significant ASV in each grouping:
  for(j in 1:nrow(combByCluster)){  #loop over all the significant ASVs
    relAbundDf[j,] <- (combByCluster[j,]/totalCounts)*100
  }
  rownames(relAbundDf) <- rownames(combByCluster) #make row names TaxonID
  # Merge taxonomy table and relAbundDF by rows that they share
  mergedDf_1 <- merge(taxTable, relAbundDf, by=0)
  rownames(mergedDf_1) <- mergedDf_1[,1] #make the column "row.names" the actual rownames again
  mergedDf_1[,1] <- NULL #remove "row.names" column
  # Get the last 3 rows of myMultipatt$sign, which are index, stat, and p-value:
  accessSign <- myMultipatt$sign[,c((ncol(myMultipatt$sign)),(ncol(myMultipatt$sign) -1),(ncol(myMultipatt$sign) -2))]
  # Merge mergedDf_1 with accessSign to add index, stat, and p-value
  mergedDf_2 <- merge(mergedDf_1,accessSign, by=0)
  results <- mergedDf_2[with(mergedDf_2, order(index)),]
  rownames(results) <- results[,1] #make row names taxon ID
  results[,1] <- NULL #remove "row.names" column
  return(results)
}

# 2. PERFORM INDICATOR SPECIES ANALYSIS
# This will perform an indicator species analysis to determine which are edge, patch, and forest specialists

# Get the metadata for these samples
I6S_postUbiquity_meta <- as.data.frame(as.matrix(sample_data(postUbiquity.ps)))
I6S_postUbiquity_meta$Meter <- gsub(pattern= " ", replacement= "", I6S_postUbiquity_meta$Meter) #take out weird spaces in meter
unique(I6S_postUbiquity_meta$Habitat) #forest, patch, and edge 

# Based on the Bray-Curtis plots, I think that edge being 40, 50, and 60 is good. So add in a column to 
# do this 
I6S_postUbiquity_meta <- I6S_postUbiquity_meta %>%
  mutate(
    habitatEdge2 = case_when(
      Meter %in% c(10, 20, 30) ~ "savanna",
      Meter %in% c(40, 50, 60) ~ "edge",
      Meter %in% c(70, 80, 90, 100) ~ "forest",
      TRUE ~ NA_character_
    )
  )

#View(I6S_postUbiquity_meta)
cbind(I6S_postUbiquity_meta$Meter, I6S_postUbiquity_meta$habitatEdge2) #spot check shows that the code above was good!

# Get a string that has the habitat designation
I6S_habitatEdge2 <- I6S_postUbiquity_meta$habitatEdge2

# Get ASV and taxonomy tables out of phyloseq:
I6S_postUbiquity_ASVs <- ASVs_outta_ps(postUbiquity.ps)
# Double check that metadata and ASV tab still in same order, so that I can use the string of habitat designations that I made above
rownames(I6S_postUbiquity_ASVs) == rownames(I6S_postUbiquity_meta)
head(I6S_postUbiquity_ASVs) #ASVs are columms and samples are rows
dim(I6S_postUbiquity_ASVs)
I6S_postUbiquity_tax <- as.data.frame(phyloseq::tax_table(postUbiquity.ps), stringsAsFactors = F)
head(I6S_postUbiquity_tax) #looks good!

# Perform indicator analysis
set.seed(9)
I6S_habitatIndic <- multipatt(x=I6S_postUbiquity_ASVs, cluster=I6S_habitatEdge2, func="r.g", duleg = TRUE, control=how(nperm = 9999)) #duleg=true means that combinations
# are not considered

summary(I6S_habitatIndic) #shows all of the taxon IDs for ASVs associated with each group, as well as
# the  r.g. value (under stat, this is the correlation index that takes into account the
# variation within and among groups), and the p-value from a permutation test.

I6S_habitatIndic <- indSpecTable2(myMultipatt=I6S_habitatIndic, taxTable= I6S_postUbiquity_tax, ASVtable = I6S_postUbiquity_ASVs, signLevel=0.01)
# View(I6S_habitatIndic)

# saved Nov. 2, 2023 
# save(I6S_habitatIndic, file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_habitatIndicTableNov2")

# HERE, I JUST LOAD IT IN CASE THE CODE ABOVE HASN'T BEEN RUN (WHICH DOES TAKE A WHILE)
load(file = "~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/I6S_habitatIndicTableNov2") #file is I6S_habitatIndic

# remove p-values greater than 0.01. I did not keep this analysis in, but added comments throughout to describe what this would 
# result in at many steps
I6S_habitatIndic_0.01 <- I6S_habitatIndic %>% #I 
  filter(p.value < 0.01)

# Assign soil or phyllosphere based on index
for (j in 1:nrow(I6S_habitatIndic)){
  if (I6S_habitatIndic$index[j] == "1") {
    I6S_habitatIndic$group[j] <- "edge" 
  } 
  if (I6S_habitatIndic$index[j] == "2") {
    I6S_habitatIndic$group[j] <- "forest" 
  }
  if (I6S_habitatIndic$index[j] == "3") {
    I6S_habitatIndic$group[j] <- "savanna" 
  }
}

sum(I6S_habitatIndic$`relAbund in forest`)

I6S_edgeSpecialistsIndex <- which(I6S_habitatIndic$group == "edge") 
length(I6S_edgeSpecialistsIndex) #37 with p-value 0.01; 139 with p-value 0.05
I6S_edgeSpecialists <- rownames(I6S_habitatIndic)[I6S_edgeSpecialistsIndex]
I6S_forestSpecialistsIndex <- which(I6S_habitatIndic$group == "forest") 
length(I6S_forestSpecialistsIndex) # 738 with 0.01 p-value; 1070 with p-value 0.05
I6S_forestSpecialists <- rownames(I6S_habitatIndic)[I6S_forestSpecialistsIndex]
I6S_savSpecialistsIndex <- which(I6S_habitatIndic$group == "savanna") 
length(I6S_savSpecialistsIndex) #781 with p value of 0.01; 1095 with p-value 0.05
I6S_savSpecialists <- rownames(I6S_habitatIndic)[I6S_savSpecialistsIndex]

##########################################################################
# 3. PLOTS FOR DIFFERENTIALLY ABUNDANT TAXA ALONG THE TRANSECT
##########################################################################
# This part makes a stacked barplot which has the percentage of reads that are forest, patch, or non-specialists at each point along the transect
# This code gets mean abundance for each ASV at each meter along the transect (regardless of EU))

# 1. First, need to get the mean ASV abundance at each point along the transect, averaging across all EUs
bacterASVsdf <- ASVs_outta_ps(postUbiquity.ps) #rows are samples, ASV abundance is columns 
bactermetaDf <- as.data.frame(as.matrix(sample_data(postUbiquity.ps))) #rows are samples, columns are all of the rest of the metadata
bactermetaDf <- as.data.frame(as.matrix(bactermetaDf)) #make format nicer
bactermetaDf$Meter <- as.numeric(bactermetaDf$Meter) #make meter numeric
bacterTaxTab <- as.data.frame(as.matrix(tax_table(postUbiquity.ps)))
unique(rownames(bacterASVsdf) == rownames(bactermetaDf)) #this is true, which shows that we can use indices from
# metadata to isolate stuff from ASV table

# Get all the rownames in ASV table for samples corresponding to each meter
rownames(bactermetaDf)[which(bactermetaDf$Meter == 10)] #this is how to pull out rownames by meter
meterVec <- c(10,20,30,40,50,60,70,80,90,100)
bacterMeterRowNamesIndices <- vector("list", length(meterVec)) #this will have all of the rownames that correspond to samples in each meter
names(bacterMeterRowNamesIndices) <- paste(meterVec, "m", sep="_") 
for (j in 1:length(meterVec)){
  bacterMeterRowNamesIndices[[j]] <- rownames(bactermetaDf)[which(bactermetaDf$Meter == meterVec[j])]
}

# For each one of the ASVs, pull out all of the samples at each point
# first, make a dataframe to hold all the final stuff:
bacterMeanASVsByMeter <- as.data.frame(matrix(nrow=10, ncol=ncol(bacterASVsdf)))
colnames(bacterMeanASVsByMeter) <- colnames(bacterASVsdf)
rownames(bacterMeanASVsByMeter) <- names(bacterMeterRowNamesIndices)

# For loop to get mean ASV abundance at each meter
for (k in 1:length(bacterMeterRowNamesIndices)){ 
  bacterMeanASVsByMeter[k,] <- t(colMeans(bacterASVsdf[bacterMeterRowNamesIndices[[k]],])) #this pulls out all of the samples at each meter.
}
# View(bacterMeanASVsByMeter)
# Check that mean abundances are correct
# ASV 1 in 10m samples
tenMeterSamps <- rownames(bactermetaDf)[which(bactermetaDf$Meter==10)]
# for ASV 1
mean(bacterASVsdf[tenMeterSamps,1]) == bacterMeanASVsByMeter[1,1]
# ASV 6 in meter 100
hundredMeterSamps <- rownames(bactermetaDf)[which(bactermetaDf$Meter==100)]
mean(bacterASVsdf[hundredMeterSamps,6]) == bacterMeanASVsByMeter[10,6]

# Make this in long, "tidy" format, and then merge with object made earlier for specialist and taxonomic info
bacterMeanASVsByMeterTidy <- bacterMeanASVsByMeter %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = ASV_1:ASV_9931, names_to = "ASV_name", values_to = "meanASVabundance")
# View(bacterMeanASVsByMeterTidy) #looks good!
colnames(bacterMeanASVsByMeterTidy)[1] <- "meter"

# Add specialist category to taxonomy table
bacterTaxTab

bacterTaxTabWithSpecialists <- bacterTaxTab %>% #this object now has all bacterial species and their specialization types
  rownames_to_column(var="ASV_name") %>% 
  mutate(
    specialist = case_when(
      ASV_name %in% I6S_edgeSpecialists ~ "edge",
      ASV_name %in% I6S_forestSpecialists ~ "forest",
      ASV_name %in% I6S_savSpecialists ~ "savanna",
      TRUE ~ NA_character_
    )
  )
length(which(is.na(bacterTaxTabWithSpecialists$specialist)==TRUE)) #3663 with 0.01 threshold; 2915 are NAs with 0.05 threshold
all.equal((length(I6S_edgeSpecialists) + length(I6S_forestSpecialists) + length(I6S_savSpecialists) + 
             length(which(is.na(bacterTaxTabWithSpecialists$specialist)==TRUE))), nrow(bacterTaxTabWithSpecialists))
bacterTaxTabWithSpecialists$specialist[which(is.na(bacterTaxTabWithSpecialists$specialist)==TRUE)] <- "non-specialist"

# Annotate taxonomy based on specialist
colnames(bacterTaxTab)
bacterMeanASVsByMeterHabitat_edge <- merge(bacterMeanASVsByMeterTidy,bacterTaxTabWithSpecialists[,], by="ASV_name") #add in which habitat specialist in & taxonomic info
colnames(bacterMeanASVsByMeterHabitat_edge)[2] <- "meter"

#### GET RELATIVE ABUNDANCES #####
# Get the total ASV abundance within each meter
I6S_ASVmeterTotal_edge <- as.data.frame(matrix(nrow=10, ncol=1))
rownames(I6S_ASVmeterTotal_edge) <- paste(meterVec, "m", sep="_") 
colnames(I6S_ASVmeterTotal_edge) <- "ASVmeterTotal"
for (j in 1:length(meterVec)){
  I6S_ASVmeterTotal_edge[j,1] <- sum(bacterMeanASVsByMeterHabitat_edge[which(bacterMeanASVsByMeterHabitat_edge$meter == rownames(I6S_ASVmeterTotal_edge)[j]),3]) #pull out meter by meter
}
I6S_ASVmeterTotal_edge

# Finally get relative abundances by dividing each value in each row by the meter total
bacterRelAbundDfs_edge <- vector("list", length(meterVec)) #this will have all of the rownames that correspond to samples in each meter
names(bacterRelAbundDfs_edge) <- rownames(I6S_ASVmeterTotal_edge)
for (h in 1:length(bacterRelAbundDfs_edge)){
  bacterRelAbundDfs_edge[[h]] <- bacterMeanASVsByMeterHabitat_edge[which(bacterMeanASVsByMeterHabitat_edge$meter ==  names(bacterRelAbundDfs_edge)[[h]]),] #this pulls out just the data for each meter
  bacterRelAbundDfs_edge[[h]]$relAbund <- (bacterRelAbundDfs_edge[[h]]$meanASVabundance/I6S_ASVmeterTotal_edge[h,1])*100
}
sum(bacterRelAbundDfs_edge[[1]]$relAbund) #great, these add up to be 100!
sum(bacterRelAbundDfs_edge[[2]]$relAbund)

# Checks:
# ASV 1, 20 meters
twentyMeterCheckData <- bacterMeanASVsByMeterHabitat_edge[which(bacterMeanASVsByMeterHabitat_edge$meter ==  names(bacterRelAbundDfs_edge)[[2]]),]
twentyMeterCheckData[1,3] == bacterRelAbundDfs_edge[[2]][1,3] #these are equal so this is working well so far. This is getting the mean ASV abundance at 20m for ASV1
# Last check:
(twentyMeterCheckData[1,3]/I6S_ASVmeterTotal_edge[2,1]) #this divides the mean abundance at 20m for ASV1 by the total ASV abundance at that meter
(twentyMeterCheckData[1,3]/I6S_ASVmeterTotal_edge[2,1])*100 == bacterRelAbundDfs_edge[[2]][1,12] #TRUE, (twentyMeterCheckData[1,3]/I6S_ASVmeterTotal_edge[2,1]) gets 
colnames(bacterRelAbundDfs_edge[[2]])
# ASV 99, 100 meters
hundredMeterCheckData <- bacterMeanASVsByMeterHabitat_edge[which(bacterMeanASVsByMeterHabitat_edge$meter ==  names(bacterRelAbundDfs_edge)[[10]]),]
which(hundredMeterCheckData$ASV_name == "ASV_99") #411th row
hundredMeterCheckData[411,3] == bacterRelAbundDfs_edge[[10]][411,3] #these are equal
(hundredMeterCheckData[411,3]/I6S_ASVmeterTotal_edge[10,1])*100 == bacterRelAbundDfs_edge[[10]][411,12] #TRUE

### rbind all of the relAbunds[[all meters]] and plug this into ggplot2 below
bacterRelAbundEDGE_df <- do.call("rbind", bacterRelAbundDfs_edge)
# View(bacterRelAbundEDGE_df)
levels(bacterRelAbundEDGE_df$specialist)

# How many specialists are in each group (will add to plot below)
length(which(bacterTaxTabWithSpecialists$specialist == "savanna")) #1095
length(which(bacterTaxTabWithSpecialists$specialist == "forest")) #1070
length(which(bacterTaxTabWithSpecialists$specialist == "edge")) #139
length(which(bacterTaxTabWithSpecialists$specialist == "non-specialist")) #2,915

# Plot it!
# Relative abundance
level_order <- names(bacterMeterRowNamesIndices) #set this to make in correct order from 10m to 100m
relAbundTransectBacterEDGE_plot <- ggplot(bacterRelAbundEDGE_df, aes(x = factor(meter, level = level_order), y = relAbund, fill = specialist)) + 
  geom_bar(stat = "identity", position = "fill")  +
  scale_fill_manual(values=c("purple","darkgreen","darkgrey","goldenrod"), labels=c("edge (n=139)", "forested matrix (n=1,070)", "non-specialists (n=2,915)", "open patch (n= 1,095)")) +
  labs(y= "Relative abundance", x = "Meter on transect") + 
  scale_x_discrete(labels=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + #change x-axis tick labels
  guides(color = guide_legend(override.aes = list(size =7))) +
  theme(legend.title = element_text(size=14)) + #increase size of legend title
  theme(legend.text = element_text(size = 12)) + #increase size of legend font
  theme(axis.text=element_text(size=12), #increase font size of axis labels
        axis.title=element_text(size=12)) + #increase font size of axis title
  theme_bw()
# quartz()
relAbundTransectBacterEDGE_plot

# Saved November 2, 2023
# save(relAbundTransectBacterEDGE_plot, file= "RobjectsSaved/relAbundTransectBacterEDGE_plot")

#Bring in fungal version of this plot and edit it to match bacterial labeling 
load(file= "RobjectsSaved/relAbundTransectFungalEDGE_plot")
relAbundTransectFungalEDGE_plot <- relAbundTransectFungalEDGE_plot + scale_fill_manual(values=c("purple","darkgreen","darkgrey","goldenrod"), labels=c("edge", "forested matrix", "non-specialists", "open patch"))

# Plot fungal and bacterial side by side!
# quartz()
grid.arrange(relAbundTransectBacterEDGE_plot, relAbundTransectFungalEDGE_plot, nrow=2)
###############################
# TOP PHYLA PLOT (code from postUbiqGraphics_ITS.R)
###############################
# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
sampDat <- phyloseq::sample_data(I6S_postUbiquity_meta)
postUbiquity.ps <- phyloseq(sampDat, tax_table(postUbiquity.ps), otu_table(postUbiquity.ps))
# Combine taxa based on phylum
postUbiq.phylum.glom <-  tax_glom(postUbiquity.ps, taxrank = "Phylum") 
tax_table(postUbiq.phylum.glom) #8 phyla
sample_data(postUbiq.phylum.glom)

# Transform sample counts based on just glommed samples
relabun.phyla.0 <- transform_sample_counts(postUbiq.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0)) 
colSums(otu_table(relabun.phyla.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each phylum

###### TOP PHYLA BY HABITAT #####

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = c("habitatEdge2"))
sample_data(relabun.phyla.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
rowSums(otu_table(relabun.phyla.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5% of total abundance
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
relabun.phylatop99.5  <- relabun.phyla.df
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abundance"

# Olpidiomycota and Rozellomycota were phyla as well, but comprised less than .5% of abundance
unique(relabun.phylatop99.5$Phylum)
# Phyla comprising at least 0.5% of total abundance
# Get unique colors for these 7 groups
colors31 <- Polychrome::glasbey.colors(32)
swatch(colors31)
colorsForPhyla7 <- unname(colors31)[c(2:4,6:8)] #remove some unwanted colors
colorsForPhyla7[7] <- "gray48" #make last one gray

# Remove p_ in phylum names
relabun.phylatop99.5$Phylum <- gsub("p__","",as.character(relabun.phylatop99.5$Phylum))
uniquePhyla <- sort(unique(relabun.phylatop99.5$Phylum)) #7 unique phyla

# re-order phyla so that "< .5% abundance" is last
relabun.phylatop99.5$Phylum <- factor(relabun.phylatop99.5$Phylum, levels=uniquePhyla[c(2:7,1)])

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum))
phylumPlot99.5percent <- phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForPhyla7) +
  theme(axis.title.y = element_text(size = 16)) +
  ylab("Relative abundance") +
  theme(axis.text.y= element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=3)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
phylumPlot99.5percent
# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

######## STOPPED HERE OCTOBER 25, 2023 #####

top_99.5p_phyla <- relabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

# This gets relative abundance across all EUs
bacterRelabunFam_grouped <- relabun.phylatop99.5 %>% 
  group_by(Phylum) %>% 
  summarize(relAbundSummedAcrossEUs = sum(Abundance)) %>% 
  mutate(relAbund = relAbundSummedAcrossEUs/6)
#View(bacterRelabunFam_grouped)
sum(bacterRelabunFam_grouped$relAbund) #this adds up to 1!

###############################
# TOP PHYLA PLOT (code from postUbiqGraphics_ITS.R)
###############################
# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
sampDat <- phyloseq::sample_data(I6S_postUbiquity_meta)
postUbiquity.ps <- phyloseq(sampDat, tax_table(postUbiquity.ps), otu_table(postUbiquity.ps))
# Combine taxa based on phylum
postUbiq.phylum.glom <-  tax_glom(postUbiquity.ps, taxrank = "Phylum") 
tax_table(postUbiq.phylum.glom) #8 phyla
sample_data(postUbiq.phylum.glom)

# Transform sample counts based on just glommed samples
relabun.phyla.0 <- transform_sample_counts(postUbiq.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0)) 
colSums(otu_table(relabun.phyla.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each phylum

###### TOP PHYLA BY HABITAT #####

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = c("habitatEdge2"))
sample_data(relabun.phyla.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
rowSums(otu_table(relabun.phyla.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5% of total abundance
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
relabun.phylatop99.5  <- relabun.phyla.df
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abundance"

# Olpidiomycota and Rozellomycota were phyla as well, but comprised less than .5% of abundance
unique(relabun.phylatop99.5$Phylum)
# Phyla comprising at least 0.5% of total abundance
# Get unique colors for these 7 groups
colors31 <- Polychrome::glasbey.colors(32)
swatch(colors31)
colorsForPhyla7 <- unname(colors31)[c(2:4,6:8)] #remove some unwanted colors
colorsForPhyla7[7] <- "gray48" #make last one gray

# Remove p_ in phylum names
relabun.phylatop99.5$Phylum <- gsub("p__","",as.character(relabun.phylatop99.5$Phylum))
uniquePhyla <- sort(unique(relabun.phylatop99.5$Phylum)) #7 unique phyla

# re-order phyla so that "< .5% abundance" is last
relabun.phylatop99.5$Phylum <- factor(relabun.phylatop99.5$Phylum, levels=uniquePhyla[c(2:7,1)])

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum))
phylumPlot99.5percent <- phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForPhyla7) +
  theme(axis.title.y = element_text(size = 16)) +
  ylab("Relative abundance") +
  theme(axis.text.y= element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=3)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
phylumPlot99.5percent
# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

