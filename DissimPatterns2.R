# Exploration of Dissimilarity Patterns 
# November 12, 2021

# This script calculates and plots the Bray-Curtis dissimilarities for the 
# meters across the transect. 
# ********** IMPORTANT ********** 
# 1. CURRENTLY, THIS TREATS TRANSECTS SEPARATELY, INSTEAD OF USING MEDIAN ABUNDANCE
# (FROM UBIQUITYMEDIANSETUP.R.). FURTHERMORE, THE INPUT HERE HAS NO UBIQUITY FILTER 
# (YET).

# Libraries
library("phyloseq")
library("ggplot2")      
library("dplyr")        
library("tibble")       
library("tidyr")
library("vegan")

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

load("trimmedJustsoils.ps") #load phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset (from 
# 16SExploratoryDataAnalysisAug2021).

# FUNCTIONS DEFINED IN THIS SCRIPT:
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

# 3. metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

# 4. TransectDissim obtains the Bray-Curtis dissimilarities between 1) 10 m and 
# all points on transect, 2) 100 m and all points on transect, and 3) 50 m
# (edge) and all points on transect. It specifically divides the EU in the 
# provided phyloseq object into its four component transects, calculates Bray
# Curtis dissimilarities, and then pulls out the comparisons mentioned above.
# It returns a list of 12 matrices (10 m comparisons, 100 m comparisons, and
# 50 m (edge) comparisons). #Note: for comparisons with calculating manually, i.e.
# proof that this function works, see "justDissimFunction.R" and DissimPatterns.R.
# (The latter does not directly compare, but manually calculates the Bray Curtis
# dissimilarities.)

transectDissim <- function(physeq){
  ASVsdf <- ASVs_outta_ps(physeq) #get ASV table using functioned defined earlier; samples are rows, ASVs = columns
  metaDf <- metadata_outta_ps(physeq) #get metadata using function defined earlier; samples are rows
  # Next three lines make a new df that has variables of interest and ASVs
  ASVsdf$Transect <- metaDf$Transect
  ASVsdf$Meter <- metaDf$Meter
  # Sample names (numbered) are rows. Columns are ASV names, except for last two which are Transect and Meter. 
  # Make a new column that combines info from meter and EU
  sampNames <- rep(NA, nrow(ASVsdf)) #pre-allocate vector for new names, length is 38 b/c 38 samples from EU 52
  for (i in 1:nrow(ASVsdf)){
    sampNames[i] <- paste(c(ASVsdf[i,(ncol(ASVsdf)-1)], ASVsdf[i,ncol(ASVsdf)]), collapse="_") #transect and meter
  } 
  ASVsdf$sampNames <- sampNames
  
  # The piped code below creates a list of 4 tibbles, corresponding to each transect
  transectLists <- ASVsdf %>%  #make multiple datasets based on transect
    dplyr::group_by(Transect) %>%  #group by EU
    dplyr::group_split()
  
  # Remove Meter and transect columns then change sampNames into the rownames 
  transectLists2 <- vector("list", 4) # preallocate list of vectors transformed transectLists, can drop once code is working
  # and just write over transectLists
  for (i in 1:length(transectLists)){
    transectLists2[[i]] <- transectLists[[i]] %>% 
      arrange(Meter) %>% #sort in ascending order
      select(-c(Meter, Transect)) %>% #drop meter and transect columns b/c info in sampNames
      column_to_rownames("sampNames") #make sampNames column in the rownames 
  }
  
  # The lines below get Bray-Curtis dissimilarities within transects (as created above)
  BCs <- vector("list", 4) #preallocate vector to store Bray-Curtis dissimilarities
  for (i in 1:length(transectLists2)){ #there are 4 tibbles within this, corresponding to each transect. 
    # they are in alphabetical order, starting with transect B
    # make a list of lists call Bray-Curtis with 4 subset lists
    # Get Bray-Curtis dissimilarities for a subset (not last two rows that are transect and meter)
    BCs[[i]] <- as.matrix(vegdist(transectLists2[[i]][,1:(ncol(transectLists2[[i]])-2)], method="bray"))
    diag(BCs[[i]]) <- NA
  }
  
  # The lines below pull out Bray Curtis distances within each transect for 10 m
  # and each other point along the transect, and 100 m and each point along transect
  ### AT SOME POINT, CAN ADD B/W EDGE AND ALL POINTS AND POOLED FOREST AND POOLED
  # FOREST AND ALL POINTS INTO CODE RIGHT HERE! WOULD NEED TO CHANGE PREALLOCATED LIST
  comps10m <- vector("list", 4) #preallocate space for all comps with 10 m for each transect
  for (i in 1:length(comps10m)){
    comps10m[[i]] <- BCs[[i]][1,]
  }
  comps100m <- vector("list", 4) #all comps with 10 m for each transect
  for (i in 1:length(comps100m)){
    comps100m[[i]] <- BCs[[i]][nrow(BCs[[i]]),]
  }
  comps50m <- vector("list", 4) #preallocate space for 50m (edge) comparisons
  index50m <- which(rownames(BCs[[1]])=="B_50") #get row number for edge
  # I use this above so that function can be flexible if samples are missing 
  for (i in 1:length(comps50m)){
    comps50m[[i]] <- BCs[[i]][index50m,]
  }
  #put it all together in list of vectors
  allComps <- c(comps10m, comps100m, comps50m)
  return(allComps)
}

###############################################################################
# Main part of script
###############################################################################

# To get dissimilarities, first split phyloseq object by EUs
EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
EU_53N_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53N")
EU_54S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_54S")
EU_8_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_8")
EU_53S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53S")
EU_10_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_10")

# Make a list with all of the EUs phyloseq objects 
listByEU <- list(EU_52_Soils.ps, EU_53N_Soils.ps, EU_54S_Soils.ps, 
                 EU_8_Soils.ps, EU_53S_Soils.ps, EU_10_Soils.ps)
# Set their names!
names(listByEU) <- c("EU_52", "EU_53N", "EU_54S", "EU_8", "EU_53S", "EU_10")
names(listByEU)

# Using function defined above, loop over all EUs to get all dissimilarity comparisons
allEUsComps <- vector("list", 6) #create an empty list for 6 elements, one per EU
for (i in 1:length(listByEU)){
  allEUsComps[[i]] <- transectDissim(listByEU[[i]])
}
names(allEUsComps) <- names(listByEU) #add names for each EU
names(allEUsComps) #correct!
length(allEUsComps$EU_52) #lengths look correct too
