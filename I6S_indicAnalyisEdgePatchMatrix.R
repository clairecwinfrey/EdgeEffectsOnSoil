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
bacterRelAbundEDGE_df$specialist <- factor(bacterRelAbundEDGE_df$specialist, levels=c("non-specialist","savanna", "edge", "forest"))
relAbundTransectBacterEDGE_plot <- ggplot(bacterRelAbundEDGE_df, aes(x = factor(meter, level = level_order), y = relAbund, fill = specialist)) + 
  geom_bar(stat = "identity", position = "fill")  +
  scale_fill_manual(values=c("darkgrey","goldenrod", "purple","darkgreen"), labels=c("non-specialists (n=2,915)", "open patch (n= 1,095)", "edge (n=139)", "forested matrix (n=1,070)")) +
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

# Saved November 18, 2023
# save(relAbundTransectBacterEDGE_plot, file= "RobjectsSaved/relAbundTransectBacterEDGE_plot")

#Bring in fungal version of this plot and edit it to match bacterial labeling 
load(file= "RobjectsSaved/relAbundTransectFungalEDGE_plot")

# Plot fungal and bacterial side by side!
# quartz()
grid.arrange(relAbundTransectBacterEDGE_plot, relAbundTransectFungalEDGE_plot, nrow=2)


###############################
# 4. MAKING TRANSECT PLOTS FOR EACH OF THE DIFFERENT EUS
###############################

#### GET RELATIVE ABUNDANCES #####
# GET MEAN ASV ABUNDANCE FOR EACH ASV WITHIN EACH EU
# 1. Need to get all of the rownames in the ASV table for samples correspondong to each meter, within each EU
I6S_meterIndices_byEU <- vector("list", length(unique(I6S_postUbiquity_meta$EU))) #this will a dataframe for each EU
names(I6S_meterIndices_byEU) <- unique(I6S_postUbiquity_meta$EU) #give this the names of the unique EUs

testEU_52_10mIndex <-intersect(which(bactermetaDf$EU == names(I6S_meterIndices_byEU)[1]), which(bactermetaDf$Meter == meterVec[1]))
bactermetaDf[testEU_52_10mIndex,] #this is correct!

# This for loop works in two nested lists. The list "I6S_meterIndices_byEU" contains 6 lists, one per EU. Within each of these EU lists,
# there is another list with 10 elements. Each of these 10 elements contains the sample names for the that EU, ordered 10-100 meters
for (i in 1:length(I6S_meterIndices_byEU)){
  # within each of these lists, separated by EU, make another list, one for each meter to hold sample indices for each meter
  I6S_meterIndices_byEU[[i]] <-vector("list", length(unique(meterVec)))
  for (j in 1:length(meterVec)){
    I6S_meterIndices_byEU[[i]][[j]] <-intersect(which(bactermetaDf$EU == names(I6S_meterIndices_byEU)[i]), which(bactermetaDf$Meter == meterVec[j]))
    names(I6S_meterIndices_byEU[[i]])[j] <- meterVec[j] #name each sub-list according to the meter it has indices for.
  }
}
# Check a few of these results:
bactermetaDf[I6S_meterIndices_byEU$EU_10$`100`,] #this shows that, indeed, these are all of the meter 100 samples within 10!
bactermetaDf[I6S_meterIndices_byEU$EU_53S$`80`,] #this shows that, indeed, these are all of the meter 80 samples within 53S!

# 2. This for loop gets mean ASV abundance at each meter in each EU (i.e., the mean abundance of each ASV across 4 transects)
# 2a. Pre-allocate lists and dataframes:
bacterMeanASVsByMeter_EUList <- vector("list", length=6) #make a list with 6 elements, one for each EU
names(bacterMeanASVsByMeter_EUList) <- unique(I6S_postUbiquity_meta$EU) #set the names equal to the names for each EU
# 2b. For loop to finish making preallocated lists (maybe not necessary...)
for (i in 1:length(bacterMeanASVsByMeter_EUList)){
  # within each of these lists, separated by EU, make a dataframe with 10 rows (one for each meter) and as many columns as ASVs
  bacterMeanASVsByMeter_EUList[[i]] <- as.data.frame(matrix(nrow= length(meterVec), ncol=ncol(bacterASVsdf)))
  colnames(bacterMeanASVsByMeter_EUList[[i]]) <- colnames(bacterASVsdf) #make the column names equal to the ASV names
  rownames(bacterMeanASVsByMeter_EUList[[i]]) <- paste0(meterVec, "m") #make the column names equal to the ASV names
}

# 2c. To do this next step, I need to go into each of the elements in the list "I6S_meterIndices_byEU", which correspond to each EU. The index i
# in the data frame above will correspond to the EU in I6S_meterIndices_byEU AND in bacterMeanASVsByMeter_EUList. The sublist in "I6S_meterIndices_byEU,
# goes with meterVec and meters. This index, k, will go with the meter sublists in I6S_meterIndices_byEU and the rows in each dataframe within
# bacterMeanASVsByMeter_EUList[[i]]. 
for (i in 1:length(bacterMeanASVsByMeter_EUList)){
  # within each of these lists, separated by EU, add have the mean ASV value across all 4 transects in each ASV (columns) with each meter as rows
  for (k in 1:length(meterVec)){ #loop over all 10 meters within each EU, i
    # here, i is list by EU, k is meter. So bacterMeanASVsByMeter_EUList is a list of 6 EUs, each with 10 lists within for each meter
    bacterMeanASVsByMeter_EUList[[i]][k,] <- colMeans(bacterASVsdf[which(rownames(bacterASVsdf) %in% rownames(bactermetaDf)[I6S_meterIndices_byEU[[i]][[k]]] == TRUE),])
  }
}

# 2d. Testing these results 
# 2d.1 -- this first part checks that the indexing is indeed working correctly, i.e. getting the correct meter/EU combo
testIndex <- which(rownames(bacterASVsdf) %in% rownames(bactermetaDf)[I6S_meterIndices_byEU$EU_53S$`80`] == TRUE) #
bacterASVsdf[which(rownames(bacterASVsdf) %in% rownames(bactermetaDf)[I6S_meterIndices_byEU$EU_53S$`80`] == TRUE),]
colMeans(bacterASVsdf[which(rownames(bacterASVsdf) %in% rownames(bactermetaDf)[I6S_meterIndices_byEU$EU_53S$`80`] == TRUE),])
rownames(bacterASVsdf[testIndex,]) #"154" "192" "200" "221"
bactermetaDf[rownames(bactermetaDf) %in% c("154","192","200","221") ==TRUE,] #yep, these are also 53SD, meter 80 
# 2d.2 -- this part checks that the mean ASVs are correct and accessible in my nested loops
colMeans(bacterASVsdf[which(rownames(bacterASVsdf) %in% rownames(bactermetaDf)[I6S_meterIndices_byEU$EU_53S$`80`] == TRUE),])
# All of these below are true, which means it's working
unname(bacterMeanASVsByMeter_EUList$EU_53S[8,] == colMeans(bacterASVsdf[which(rownames(bacterASVsdf) %in% rownames(bactermetaDf)[I6S_meterIndices_byEU$EU_53S$`80`] == TRUE),]))

# 3.Make this in long, "tidy" format,
bacterMeanASVsByMeterTIDY_EUList <- vector("list", length=6) #make a list with 6 elements, one for each EU
names(bacterMeanASVsByMeterTIDY_EUList) <- names(bacterMeanASVsByMeter_EUList)
for (i in 1:length(bacterMeanASVsByMeter_EUList)){ #looping over each EU's big dataframe
  bacterMeanASVsByMeterTIDY_EUList[[i]] <-  bacterMeanASVsByMeter_EUList[[i]] %>% 
    rownames_to_column() %>% 
    pivot_longer(cols = colnames(bacterMeanASVsByMeter_EUList[[i]])[1]:colnames(bacterMeanASVsByMeter_EUList[[i]])[ncol(bacterMeanASVsByMeter_EUList[[i]])], names_to = "ASV_name", values_to = "meanASVabundance")
  colnames(bacterMeanASVsByMeterTIDY_EUList[[i]])[1] <- "meter" #re-name rowname as meter
}

# 4. Merge with object made earlier for specialist and taxonomic information
colnames(bacterTaxTab)
bacterMeanASVsByMeterHabitatEdge_EUList <- vector("list", length=6) #make a list with 6 elements, one for each EU
names(bacterMeanASVsByMeterHabitatEdge_EUList) <- names(bacterMeanASVsByMeter_EUList)
for (i in 1:length(bacterMeanASVsByMeterHabitatEdge_EUList)){ #looping over each EU's big dataframe
  bacterMeanASVsByMeterHabitatEdge_EUList[[i]] <- merge(bacterMeanASVsByMeterTIDY_EUList[[i]],bacterTaxTabWithSpecialists[,], by="ASV_name") #add in which habitat specialist in & taxonomic info
}
names(bacterMeanASVsByMeterHabitatEdge_EUList)

# 5. Get the total ASV abundances within each meter within each EU (based on mean abundance at each)
bacterASVmeterTotal_edge_EUList <- vector("list", length=6) #make a list with 6 elements, one for each EU
names(bacterASVmeterTotal_edge_EUList) <- names(bacterMeanASVsByMeter_EUList)
# 5b. For loop to finish making preallocated lists (maybe not necessary...)
for (i in 1:length(bacterASVmeterTotal_edge_EUList)){
  # within each of these lists, separated by EU, make a dataframe with 10 rows (one for each meter) and only one column (for the single of interest)
  bacterASVmeterTotal_edge_EUList[[i]] <- as.data.frame(matrix(nrow=length(meterVec), ncol=1))
  colnames(bacterASVmeterTotal_edge_EUList[[i]]) <- "ASVmeterTotal" #make the column names equal to the ASV names
  rownames(bacterASVmeterTotal_edge_EUList[[i]]) <- paste0(meterVec, "m") #make the rownames names equal to the ASV names
  for (j in 1:length(meterVec)){
    bacterASVmeterTotal_edge_EUList[[i]][j,1] <- sum(bacterMeanASVsByMeterHabitatEdge_EUList[[i]][which(bacterMeanASVsByMeterHabitatEdge_EUList[[i]]$meter == rownames(bacterASVmeterTotal_edge_EUList[[i]])[j]),3]) #pull out meter by meter
  }
}
# 5c. Test a few
sumTest100_EU52 <- bacterMeanASVsByMeterHabitatEdge_EUList$EU_52 %>% 
  filter(meter == "100m")
sum(sumTest100_EU52$meanASVabundance) == bacterASVmeterTotal_edge_EUList$EU_52[10,] #yeah!
sumTest60_EU8 <- bacterMeanASVsByMeterHabitatEdge_EUList$EU_8 %>% 
  filter(meter == "60m")
sum(sumTest60_EU8$meanASVabundance) == bacterASVmeterTotal_edge_EUList$EU_8[6,] #yeah!

# 5d. Finally, make rownames a column for merging in the next step below
for (i in 1:length(bacterASVmeterTotal_edge_EUList)){
  bacterASVmeterTotal_edge_EUList[[i]] <- bacterASVmeterTotal_edge_EUList[[i]] %>% 
    rownames_to_column(var="meter")
}

# 6. Actually get the dang relative abundances by dividing each value in each row by the meter total
# For this last part, I need to take the meanASVabundance in each of the 6 bacterMeanASVsByMeterHabitatEdge_EUList, and divide it
# by the total for each meter, which is found in bacterASVmeterTotal_edge_EUList. Then this will be multiplied by 100. 
# In other words, in each in each bacterMeanASVsByMeterHabitatEdge_EUList, I can divide all of the values in the column called "meanASVabundance" by the 
# the bacterASVmeterTotal_edge_EUList that corresponds to its meter value.

# 6a. Merge each bacterMeanASVsByMeterHabitatEdge_EUList and each bacterASVmeterTotal_edge_EUList, organized by each EU, so that the total ASV abundance 
# is now in bacterMeanASVsByMeterHabitatEdge_EUList
for (j in 1:length(bacterMeanASVsByMeterHabitatEdge_EUList)){
  bacterMeanASVsByMeterHabitatEdge_EUList[[j]] <- merge(bacterMeanASVsByMeterHabitatEdge_EUList[[j]], bacterASVmeterTotal_edge_EUList[[j]], by="meter", all.x= TRUE)
}
# check a few of the merges
unique((bacterMeanASVsByMeterHabitatEdge_EUList$EU_52 %>% 
          filter(meter== "100m"))$ASVmeterTotal) == bacterASVmeterTotal_edge_EUList$EU_52[10,2] #nice this worked
unique((bacterMeanASVsByMeterHabitatEdge_EUList$EU_53S %>% 
          filter(meter== "70m"))$ASVmeterTotal) == bacterASVmeterTotal_edge_EUList$EU_53S[7,2] #nice this worked
# 6b. Divide each mean ASV abundance (within each EU), by the total ASV abundance in that EU at each meter and multiply by 100
bacterRelAbundDfs_edge_byEU <- bacterMeanASVsByMeterHabitatEdge_EUList #this will a dataframe for each EU for the final stuff, but since I want to build off of
# bacterMeanASVsByMeterHabitatEdge_EUList, I will duplicate it.
for (k in 1:length(bacterMeanASVsByMeterHabitatEdge_EUList)){
  # this next line makes it so that the ASV relative abundance is the last column (column 13) of this dataframe in each elements of the list
  bacterRelAbundDfs_edge_byEU[[k]][,13] <- (bacterMeanASVsByMeterHabitatEdge_EUList[[k]]$meanASVabundance/bacterMeanASVsByMeterHabitatEdge_EUList[[k]]$ASVmeterTotal)*100
  colnames(bacterRelAbundDfs_edge_byEU[[k]])[ncol(bacterRelAbundDfs_edge_byEU[[k]])] <- "ASVrelativeAbundance"
}
# 6c. Check a few to make sure that they add up to 100-- and they do!
sum((bacterRelAbundDfs_edge_byEU$EU_52 %>% 
       filter(meter == "10m"))$ASVrelativeAbundance)
sum((bacterRelAbundDfs_edge_byEU$EU_8 %>% 
       filter(meter == "70m"))$ASVrelativeAbundance)

# 7. Finally, perform a total global checks to make sure that all of this is right (one is fine since I do checks throughout!)
# 7a. EU_8, ASV 1. 
colnames(bacterRelAbundDfs_edge_byEU$EU_8)
# Pull out what I just calculated-- meter and ASV relative abundance
justEU8_ASV1rel <- bacterRelAbundDfs_edge_byEU$EU_8[which(bacterRelAbundDfs_edge_byEU$EU_8$ASV_name == "ASV_1"), c(1, 13)]
# Now calculate it from the get-go
EU_8_meter100names <- rownames(bactermetaDf %>% #this gets the sample names for EU 8 at meter 100
                                 filter(EU == "EU_8" & Meter == 100))
rownames(bactermetaDf[I6S_meterIndices_byEU$EU_8$`100`,]) == EU_8_meter100names #okay, so this is indexing as I wanted
# Get subset of the data
EU_8_meter100_ASVs <- bacterASVsdf[which(rownames(bacterASVsdf) %in% EU_8_meter100names == TRUE),] #
EU_8_meter100_ASVsMeans <- colMeans(EU_8_meter100_ASVs) #this gets the mean ASV abundance for each ASV across all 4 samples (i.e. across EU 8, meter 100 samples)
EU_8_meter100_ASVsSUM <- sum(EU_8_meter100_ASVsMeans)
sum(colSums(EU_8_meter100_ASVs)) == sum(EU_8_meter100_ASVs)
(EU_8_meter100_ASVsMeans[1]/EU_8_meter100_ASVsSUM)*100
sum((EU_8_meter100_ASVsMeans/EU_8_meter100_ASVsSUM)*100) #this adds up to 100, as expected!
justEU8_ASV1rel[1,2] == (EU_8_meter100_ASVsMeans[1]/EU_8_meter100_ASVsSUM)*100 #yay, this worked so my code worked!

# 8. FINALLY, PLOT IT
# Plot it!
# Relative abundance (I removed legend for better plotting of everything together, but can use legend for overall plot made above and add in
# in PowerPoint.)
level_order2 <- c("10m","20m","30m" ,"40m","50m","60m","70m","80m","90m", "100m" )  #set this to make in correct order from 10m to 100m
I6S_indicatorsByEUplotsList <- vector("list", length=6) #make a list with 6 elements, one for each EU
names(I6S_indicatorsByEUplotsList) <- names(bacterRelAbundDfs_edge_byEU) #name them like the EUs
for (m in 1:length(bacterRelAbundDfs_edge_byEU)){
  bacterRelAbundDfs_edge_byEU[[m]]$specialist <- factor(bacterRelAbundDfs_edge_byEU[[m]]$specialist, levels=c("non-specialist","savanna", "edge", "forest"))
  I6S_indicatorsByEUplotsList[[m]] <- ggplot(bacterRelAbundDfs_edge_byEU[[m]], aes(x = factor(meter, level = level_order2), y = ASVrelativeAbundance, fill = specialist)) + 
    geom_bar(stat = "identity", position = "fill")  +
    scale_fill_manual(values=c("darkgrey","goldenrod", "purple","darkgreen"), labels=c("non-specialists (n=2,915)", "open patch (n= 1,095)", "edge (n=139)", "forested matrix (n=1,070)")) +
    labs(y= "Relative abundance", x = "Meter on transect") + 
    scale_x_discrete(labels=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + #change x-axis tick labels
    guides(color = guide_legend(override.aes = list(size =7))) +
    theme(legend.title = element_text(size=8)) + #increase size of legend title
    theme(legend.text = element_text(size = 8)) + #increase size of legend font
    theme(axis.text=element_text(size=8), #increase font size of axis labels
          axis.title=element_text(size=8)) + #increase font size of axis title
    theme_bw() + theme(legend.position = "none") +
    ggtitle(paste0("plot for ", names(I6S_indicatorsByEUplotsList)[m]))
  
}

# quartz() 
grid.arrange(I6S_indicatorsByEUplotsList[[3]], I6S_indicatorsByEUplotsList[[5]], I6S_indicatorsByEUplotsList[[2]],
             I6S_indicatorsByEUplotsList[[1]], I6S_indicatorsByEUplotsList[[6]], I6S_indicatorsByEUplotsList[[4]],
             ncol=3)
# This above was plotted in the window to the right and saved as a JPEG, so that the weird horizontal bars wouldn't be there
# Saved November 18, 2023
# save(I6S_indicatorsByEUplotsList, file= "RobjectsSaved/I6S_indicatorsByEUplotsList")

##########################################################################
# 5. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################
head(bacterTaxTabWithSpecialists)
unique(bacterTaxTabWithSpecialists$specialist)
bacterTaxTabWithSpecialists$specialist <- factor(bacterTaxTabWithSpecialists$specialist, levels=c("non-specialist","savanna", "edge", "forest"))
diffAbund_16S_stackedBarplotPhyla <- ggplot(bacterTaxTabWithSpecialists, aes(fill=specialist, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","goldenrod", "purple","darkgreen"), labels=c("non-specialists (n=2,915)", "open patch (n= 1,095)", "edge (n=139)", "forested matrix (n=1,070)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  theme(legend.position = "none") + #remove legend since I'll change it in PP anyway
  theme(axis.text=element_text(size=12), #increase font size of axis labels
        axis.title=element_text(size=12)) + #increase font size of axis title
  theme(axis.title.x = element_blank()) + #remove x axis label
  theme(axis.title.y = element_blank()) + #remove y axis label
  scale_y_continuous(breaks=seq(0,1200,by=200), limits= c(0, 1200))

# quartz()
diffAbund_16S_stackedBarplotPhyla

# Below saved November 20, 2023 so that it can be added to a 2 paneled plot with fungal plot!
# save(diffAbund_16S_stackedBarplotPhyla, file="RobjectsSaved/I6S_indicStackedBarplotPhyla_plot")

# A few checks to make sure that the counting above is working as expected
# Acidobacteria
length(which(bacterTaxTabWithSpecialists$Phylum=="Acidobacteria")) #1043
acido_index <- which(bacterTaxTabWithSpecialists$Phylum=="Acidobacteria")
length(which(bacterTaxTabWithSpecialists[acido_index,]$specialist=="forest")) #327 forest specialists within Acidobacteria
length(which(bacterTaxTabWithSpecialists[acido_index,]$specialist=="savanna")) #159 patch specialists within Acidobacteria
length(which(bacterTaxTabWithSpecialists[acido_index,]$specialist=="edge")) #32 edge specialists within Acidobacteria
length(which(bacterTaxTabWithSpecialists[acido_index,]$specialist=="non-specialist")) #525 remaining ASVs
(327+ 159+ 32 + 525) == length(which(bacterTaxTabWithSpecialists$Phylum=="Acidobacteria"))

#Chloroflexi
length(which(bacterTaxTabWithSpecialists$Phylum=="Chloroflexi")) #498
chloro_index <- which(bacterTaxTabWithSpecialists$Phylum=="Chloroflexi")
length(which(bacterTaxTabWithSpecialists[chloro_index,]$specialist=="forest")) #10 forest specialists 
length(which(bacterTaxTabWithSpecialists[chloro_index,]$specialist=="savanna")) #244 patch specialists
length(which(bacterTaxTabWithSpecialists[chloro_index,]$specialist=="edge")) #7 patch specialists
length(which(bacterTaxTabWithSpecialists[chloro_index,]$specialist=="non-specialist")) #229 remaining ASVs
(10+244+7+237) == length(which(bacterTaxTabWithSpecialists$Phylum=="Chloroflexi"))

# Construct two-paneled figure with 16S and ITS differentially abundant stacked barcharts side by side
# Load in previously made 16S figure (made in "EdgeEffectsbyASV_allSites.R")
load(file="RobjectsSaved/ITS_indicStackedBarplotPhyla_plot")
# quartz()
grid.arrange(diffAbund_16S_stackedBarplotPhyla, diffAbund_ITS_stackedBarplotPhyla, nrow=2)

###############################
# 6. TOP PHYLA PLOT (code from postUbiqGraphics_ITS.R)
###############################
# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
sampDat <- phyloseq::sample_data(I6S_postUbiquity_meta) #this has the distinction of edge/forest/savanna with
# the edge as 40-60 meters
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
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = c("habitatEdge2")) #this has correct habitat stuff
sample_data(relabun.phyla.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
rowSums(otu_table(relabun.phyla.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5% of total abundance
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
relabun.phylatop99.5  <- relabun.phyla.df
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abundance"

unique(relabun.phylatop99.5$Phylum) #need 11, plus the gray for the < .5 abundance 
# Phyla comprising at least 0.5% of total abundance
# Get unique colors for these 11 groups
colors31 <- Polychrome::glasbey.colors(32)
swatch(colors31)
colorsForPhyla12 <- unname(colors31)[c(2:4,6:13)] #remove some unwanted colors
colorsForPhyla12[12] <- "gray48" #make last one gray

# Remove p_ in phylum names
relabun.phylatop99.5$Phylum <- gsub("p__","",as.character(relabun.phylatop99.5$Phylum))
uniquePhyla <- sort(unique(relabun.phylatop99.5$Phylum)) #12 unique phyla

# re-order phyla so that "< .5% abundance" is last
relabun.phylatop99.5$Phylum <- factor(relabun.phylatop99.5$Phylum, levels=uniquePhyla[c(2:11,12)])
habOrder <- c("savanna", "edge", "forest") #set correct order
# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x= factor(Sample, level = habOrder), y=Abundance, fill=Phylum))
phylumPlot99.5percent <- phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  scale_y_continuous(expand = c(0, 0)) + #this line removes weird white horizontal lines between colors of same color
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForPhyla12) +
  theme(axis.title.y = element_text(size = 16)) +
  ylab("Relative abundance") +
  theme(axis.text.y= element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=6)) + theme(legend.text = element_text(colour="black", size = 10)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
phylumPlot99.5percent
# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

#saved Nov. 16, 2023
# save(phylumPlot99.5percent, file= "/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Figures/prokaryotesPhylumHabitatEdgePlot99.5percent")

###############################
# 7. TOP FAMILY PLOT (code from postUbiqGraphics_ITS.R)
###############################

# Combine taxa based on family
postUbiq.family.glom <-  tax_glom(postUbiquity.ps, taxrank = "Family") 
sample_data(postUbiq.family.glom)

# Transform sample counts based on just glommed samples
relabun.family.0 <- transform_sample_counts(postUbiq.family.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.family.0)) 
colSums(otu_table(relabun.family.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each phylum

###### TOP FAMILIES BY HABITAT #####

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.family.1 <- merge_samples(relabun.family.0, group = c("habitatEdge2")) #this has correct habitat stuff
sample_data(relabun.family.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.family.2 <- transform_sample_counts(relabun.family.1, function(x) x / sum(x))
rowSums(otu_table(relabun.family.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5% of total abundance
relabun.family.df <-psmelt(relabun.family.2)
dim(relabun.family.df) #
relabun.familytop99  <- relabun.family.df
relabun.familytop99$Family[relabun.familytop99$Abundance < 0.01] <- "< 1% abundance"

unique(relabun.familytop99$Family) #need 19, plus the gray for the < 1 abundance 
# Phyla comprising at least 1% of total abundance
# Get unique colors for these 20 groups
swatch(colors31)
colorsForFams20 <- unname(colors31)[c(2:4,6:21)] #remove some unwanted colors
colorsForFams20[20] <- "gray48" #make last one gray

# Do some clean up
# Clean up family names
relabun.familytop99$Family[which(relabun.familytop99$Family == "Acidobacteriaceae_(Subgroup_1)")] <- "Acidobacteriaceae (s.g. 1)"
relabun.familytop99$Family[which(relabun.familytop99$Family == "Solibacteraceae_(Subgroup_3)")] <- "Solibacteraceae (s.g. 3)"
#relabun.familytop99$Family <- gsub("Acidobacteriaceae_(Subgroup_1)","Acidobacteriaceae (s.g. 1)",as.character(relabun.familytop99$Family))
#relabun.familytop99$Family <- gsub("Solibacteraceae_(Subgroup_3)","Solibacteraceae (s.g. 3)",as.character(relabun.familytop99$Family))
# make NAs "Unknown families"
relabun.familytop99$Family[which(relabun.familytop99$Family== "NA")] <- "Unknown families"
relabun.familytop99$Family[which(relabun.familytop99$Family== "Unknown_Family")] <- "Unknown families"
uniqueFamilies <- sort(unique(relabun.familytop99$Family)) #20 + 1 less than 1% abundance unique phyla

# re-order phyla so that "< .5% abundance" is last
relabun.familytop99$Family <- factor(relabun.familytop99$Family, levels=uniqueFamilies[c(2:16,18:20, 17, 1)])
habOrder <- c("savanna", "edge", "forest") #set correct order
# Families comprising at least 1% of total abundance
proksFamilyPlot99percent <- ggplot(data=relabun.familytop99, aes(x= factor(Sample, level = habOrder), y=Abundance, fill=Family))
proksFamilyPlot99percent <- proksFamilyPlot99percent + geom_bar(aes(), stat="identity", position="fill") +
  scale_y_continuous(expand = c(0, 0)) + #this line removes weird white horizontal lines between colors of same color
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForFams20) +
  theme(axis.title.y = element_text(size = 16)) +
  ylab("Relative abundance") +
  theme(axis.text.y= element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=10)) + theme(legend.text = element_text(colour="black", size = 10)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
proksFamilyPlot99percent

# Plot these together
# quartz()
grid.arrange(phylumPlot99.5percent, proksFamilyPlot99percent, nrow=1)

# Plot the plots without legends so that they are same formatting/size for paper. Will add
# back in legends from other versions in PowerPoint
phylumPlot99.5percent_noLeg <- phylumPlot99.5percent + theme(legend.position = "none")
phylumPlot99.5percent_noLeg
proksFamilyPlot99percent_noLeg <- proksFamilyPlot99percent + theme(legend.position = "none")
proksFamilyPlot99percent_noLeg
# quartz()
grid.arrange(phylumPlot99.5percent_noLeg, proksFamilyPlot99percent_noLeg, nrow=1)


#saved Nov. 16, 2023
# save(proksFamilyPlot99percent, file= "/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Figures/proksFamilyPlot99percent_plot")
