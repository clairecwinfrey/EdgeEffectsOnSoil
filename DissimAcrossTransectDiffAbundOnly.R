# Bray-Curtis Dissimilarities - Using differentially abundant taxa ONLY  
# November 19, 2021

# This script calculates and plots the Bray-Curtis dissimilarities for the 
# meters across the transect.*IMPORTANT* This is only a subset of the data,
# specifically the differentially abundant 1) prokaryote and 2) fungal ASVs

# Libraries
library("phyloseq")
library("tidyverse")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Load data
load("RobjectsSaved/prokaryotesDiffAbund.ps") #load differentially abundant bacteria ps object made in DiffAbundProkaryotes.R
load("RobjectsSaved/diffAbunDat_tidy_PROKARYOTES") #made in DiffAbundProkaryotes.R
load("RobjectsSaved/prokaryotesASVmeans_AllEUs") #made in ProkaryotesEdgeEffectsZscores.R
load("RobjectsSaved/fungiDiffAbund.ps") #load differentially abundant fungal taxa os object made in DiffAbundFungi.R
load("RobjectsSaved/diffAbunDat_tidy_FUNGI") #made in DiffAbundFungi.R
load("RobjectsSaved/fungiASVmeans_AllEUs") #load differentially abundant ASVs' means calculated in FungiAbundanceZScores.R

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

# 3. pssd2veg takes the phyloseq metadata table and converts it to a dataframe
# from: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) { 
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# 4. 
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

# 5. TransectDissim3 obtains the Bray-Curtis dissimilarities between 1) 10 m and 
# all points on transect, 2) 100 m and all points on transect, and 3) 50 m
# (edge) and all points on transect. It specifically divides the EU in the 
# provided phyloseq object into its four component transects, calculates Bray
# Curtis dissimilarities, and then pulls out the comparisons mentioned above.
# It returns a dataframe with meter, Bray-Curtis value, transect, and whether
# the dissimilarity is compared to 10 m or to 100 m

# Note: for comparisons with calculating manually, i.e.proof that this
# function works, see "justDissimFunction.R" and DissimPatterns.R.
# (The latter does not directly compare, but manually calculates the Bray Curtis
# dissimilarities. The former also has a version of the function that results
# in only the comparison numbers (i.e. no Transect and meter as columns)

####### MAKEA FUNCTION THAT HAS IT FORMATTED FOR GGPLOT2!!!
transectDissim3 <- function(physeq) {
  ASVsdf <- ASVs_outta_ps(physeq) #get ASV table using functioned defined earlier; samples are rows, ASVs = columns
  metaDf <- metadata_outta_ps(physeq) #get metadata using function defined earlier; samples are rows
  # Next six lines make a new df that has variables of interest and ASVs
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
  
  # Use a for loop to arrange by meter and make sampNames row names (because code below doesn't work on list of lists)
  transectLists2 <- vector("list", 4) # preallocate list of vectors transformed transectLists, can drop once code is working
  # and just write over transectLists
  for (i in 1:length(transectLists)) {
    transectLists2[[i]] <- transectLists[[i]] %>% 
      arrange(Meter) %>% #sort in ascending order
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
  
  # Get 10 m and 100 meter comparisons, WITHIN TRANSECTS
  # The lines below pull out Bray Curtis distances within each transect for 10 m
  # and each other point along the transect, and 100 m and each point along transect
  ## Comparisons with 10 m
  comps10m <- vector("list", 4) #preallocate space for all comps with 10 m for each transect
  for (i in 1:length(comps10m)){
    #Below, preallocate matrix as each element in the list
    comps10m[[i]] <- matrix(data=NA, nrow= nrow(transectLists2[[i]]), ncol=3)
    for (j in 1:nrow(comps10m[[i]])) {
      #Make first column in comps10 "Meter" from transectLists2
      comps10m[[i]][j,1] <- transectLists2[[i]][j,ncol(transectLists2[[i]])] 
      #Make second column in comps10 B-C dissimilarity values with 10m
      comps10m[[i]][j,2] <- BCs[[i]][j,1] #so this is each row, adding in col 1 comps
      #Make third column in comps10 transect
      comps10m[[i]][j,3] <- transectLists2[[i]][j,(ncol(transectLists2[[i]])-1)]
    }
  }
  names(comps10m) <- c("B_10m", "L_10m", "R_10m", "T_10m")
  
  ## Comparisons with 100 m
  comps100m <- vector("list", 4) #pre-allocate space for all comps with 10 m for each transect
  for (i in 1:length(comps100m)){
    # Below, preallocate matrix as each element in the list
    comps100m[[i]] <- matrix(data=NA, nrow= nrow(transectLists2[[i]]), ncol=3)
    for (j in 1:nrow(comps100m[[i]])) {
      #Make first column in comps100 "Meter" from transectLists2
      comps100m[[i]][j,1] <- transectLists2[[i]][j,ncol(transectLists2[[i]])] 
      #Make second column in comps100 B-C dissimilarity values with 100m
      comps100m[[i]][j,2] <- (BCs[[i]][j,ncol(BCs[[i]])]) #last row in BCs is 100m 
      # (not missing for any transect in any EU)
      #Make third column in comps10 transect
      comps100m[[i]][j,3] <- transectLists2[[i]][j,(ncol(transectLists2[[i]])-1)]
    }
  }
  names(comps100m) <- c("B_100m", "L_100m", "R_100m", "T_100m")
  
  # Put it all together in a dataframe that is convenient for plotting
  # First get number of rows for use below
  nRow10m1 <- nrow(comps10m[[1]])
  nRow10m2 <- nrow(comps10m[[2]])
  nRow10m3 <- nrow(comps10m[[3]])
  nRow10m4 <- nrow(comps10m[[4]])
  nRow100m1 <- nrow(comps100m[[1]])
  nRow100m2 <- nrow(comps100m[[2]])
  nRow100m3 <- nrow(comps100m[[3]])
  nRow100m4 <- nrow(comps100m[[4]])
  
  # Put it all together in a dataframe that is convenient for plotting
  dissimResults <- as.data.frame(matrix(nrow= nrow(ASVsdf)*2, ncol= 4)) #make a dataframe
  dissimResults[1:nrow(ASVsdf),4] <- "10m" #first half of dataframe has 10m comparisons
  dissimResults[(nrow(ASVsdf)+1):(nrow(ASVsdf)+nrow(ASVsdf)),4] <- "100m" #last half were 100m comparisons
  dissimResults[1:nRow10m1,1:3] <- as.data.frame(comps10m[[1]])
  dissimResults[(nRow10m1 +1):(nRow10m1+nRow10m2), 1:3] <- as.data.frame(comps10m[[2]])
  dissimResults[(nRow10m1+nRow10m2 +1):(nRow10m1+nRow10m2+nRow10m3), 1:3] <- as.data.frame(comps10m[[3]])
  dissimResults[(nRow10m1+nRow10m2+nRow10m3+1):(nRow10m1+nRow10m2+nRow10m3+nRow10m4),1:3] <- as.data.frame(comps10m[[4]])
  dissimResults[(nRow10m1+nRow10m2+nRow10m3+nRow10m4+1):(nRow10m1+nRow10m2+nRow10m3+nRow10m4+nRow100m1),1:3] <- as.data.frame(comps100m[[1]])
  dissimResults[(nRow10m1+nRow10m2+nRow10m3+nRow10m4+nRow100m1+1):(nRow10m1+nRow10m2+nRow10m3+nRow10m4+nRow100m1+nRow100m2),1:3] <- as.data.frame(comps100m[[2]])
  dissimResults[(nRow10m1+nRow10m2+nRow10m3+nRow10m4+nRow100m1+nRow100m2+1):(nRow10m1+nRow10m2+nRow10m3+nRow10m4+nRow100m1+nRow100m2+nRow100m3),1:3] <- as.data.frame(comps100m[[3]])
  dissimResults[(nRow10m1+nRow10m2+nRow10m3+nRow10m4+nRow100m1+nRow100m2+nRow100m3+1):nrow(dissimResults),1:3] <- as.data.frame(comps100m[[4]])
  colnames(dissimResults) <- c("meter", "Bray-Curtis", "transect", "comparison")
  
  return(dissimResults)
}

###############################################################################
# Main part of script
###############################################################################

#############
# PROKARYOTES: 
#############

###################################
# 1. Getting dissimilarity values for all transects, all EUs

# First, get only the ASVs that are forest specialists
colnames(sample_data(prokaryotesDiffAbund.ps))
prok_forestNames <- as.character(diffAbunDat_tidy_PROKARYOTES$ASV_name[which(diffAbunDat_tidy_PROKARYOTES$Habitat=="forest")])
prokaryotes_forest_DiffAbund.ps <- prune_taxa(prok_forestNames, prokaryotesDiffAbund.ps)

# To get dissimilarities, first split phyloseq object by EUs and remove ASVs that don't occur in subset
EU_52_prok.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_52")
EU_52_prok.ps <- prune_taxa(taxa_sums(EU_52_prok.ps) > 0, EU_52_prok.ps) #remove non-occuring ASVs

EU_53N_prok.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_53N")
EU_53N_prok.ps <- prune_taxa(taxa_sums(EU_53N_prok.ps) > 0, EU_53N_prok.ps) #remove non-occuring ASVs

EU_54S_prok.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_54S")
EU_54S_prok.ps <- prune_taxa(taxa_sums(EU_54S_prok.ps) > 0, EU_54S_prok.ps) #remove non-occuring ASVs

EU_8_prok.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_8")
EU_8_prok.ps <- prune_taxa(taxa_sums(EU_8_prok.ps) > 0, EU_8_prok.ps) #remove non-occuring ASVs

EU_53S_prok.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_53S")
EU_53S_prok.ps <- prune_taxa(taxa_sums(EU_53S_prok.ps) > 0, EU_53S_prok.ps) #remove non-occuring ASVs

EU_10_prok.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_10")
EU_10_prok.ps <- prune_taxa(taxa_sums(EU_10_prok.ps) > 0, EU_10_prok.ps) #remove non-occuring ASVs

# Run function on each EU's phyloseq object and then make sure first two columns are numeric
prok_dissim_52 <- transectDissim3(EU_52_prok.ps)
prok_dissim_52[,1] <- as.numeric(prok_dissim_52[,1])
prok_dissim_52[,2] <- as.numeric(prok_dissim_52[,2])

prok_dissim_53N <- transectDissim3(EU_53N_prok.ps)
prok_dissim_53N[,1] <- as.numeric(prok_dissim_53N[,1])
prok_dissim_53N[,2] <- as.numeric(prok_dissim_53N[,2])

prok_dissim_54S <- transectDissim3(EU_54S_prok.ps)
prok_dissim_54S[,1] <- as.numeric(prok_dissim_54S[,1])
prok_dissim_54S[,2] <- as.numeric(prok_dissim_54S[,2])

prok_dissim_8 <- transectDissim3(EU_8_prok.ps)
prok_dissim_8[,1] <- as.numeric(prok_dissim_8[,1])
prok_dissim_8[,2] <- as.numeric(prok_dissim_8[,2])

prok_dissim_53S <- transectDissim3(EU_53S_prok.ps)
prok_dissim_53S[,1] <- as.numeric(prok_dissim_53S[,1])
prok_dissim_53S[,2] <- as.numeric(prok_dissim_53S[,2])

prok_dissim_10 <- transectDissim3(EU_10_prok.ps)
prok_dissim_10[,1] <- as.numeric(prok_dissim_10[,1])
prok_dissim_10[,2] <- as.numeric(prok_dissim_10[,2])

#############
# Plotting
# Function above spits this out as it needs to be to plot using ggplot2. 
# (I commented out lines comparing with 10m)

prok_ggEU_52 <- ggplot() + 
  #geom_line(data=prok_dissim_52[1:38,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=prok_dissim_52[39:76,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 52") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
#quartz()
prok_ggEU_52
prok_ggEU_53N <- ggplot() + 
  #geom_line(data=prok_dissim_53N[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=prok_dissim_53N[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 53N") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
prok_ggEU_53N
#quartz()
prok_ggEU_54S <- ggplot() + 
  #geom_line(data=prok_dissim_54S[1:38,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=prok_dissim_54S[39:76,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 54S") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
prok_ggEU_54S
#quartz()
prok_ggEU_8 <- ggplot() + 
  #geom_line(data=prok_dissim_8[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=prok_dissim_8[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 8") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
prok_ggEU_8
#quartz()
prok_ggEU_53S <- ggplot() + 
  #geom_line(data=prok_dissim_53S[1:40,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=prok_dissim_53S[41:80,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 53S") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
prok_ggEU_53S

#quartz()
prok_ggEU_10 <- ggplot() + 
  #geom_line(data=prok_dissim_10[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=prok_dissim_10[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 10") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
prok_ggEU_10 

quartz()
require(gridExtra)
grid.arrange(prok_ggEU_52, prok_ggEU_53N, prok_ggEU_54S, prok_ggEU_8, prok_ggEU_53S, prok_ggEU_10, ncol=3)

###################################
# 2. Getting dissimilarity values where MEAN is that used to calculate z-scores
# (i.e. each differentially abundant ASV at each meter)

prokaryotesASVmeans_AllEUs #made in ProkaryotesEdgeEffectsZscores.R
colnames(prokaryotesASVmeans_AllEUs)

# 1. CALCULATE BRAY-CURTIS DISSIMILARITIES ALONG TRANSECTS WITHIN EACH EU

# The lines below get the Bray-Curtis dissimilarities along the transect, comparing each point along
# the transect to 10 meters and to 100 meters 

# The piped code below creates a list of 6 tibbles, corresponding to each EU
prok_EULists <- prokaryotesASVmeans_AllEUs %>%  #make multiple datasets based on EU
  dplyr::group_by(EU) %>%  #group by EU
  dplyr::group_split() 
names(prok_EULists) <- sort(unique(prokaryotesASVmeans_AllEUs$EU))

# Use a for loop to make ASV_names row names (have to do this step after splitting by EUs so that there are no duplicated row names)
for (i in 1:length(prok_EULists)) {
  prok_EULists[[i]] <- prok_EULists[[i]] %>% 
    column_to_rownames("ASV_name")  #make sampNames column in the rownames
}

# For vegdist, rows need to be samples and columns should be species.
# 1) pull out only mean ASV abundances, 2) make ASV names the rownames,
# and 3) flip it all!
prok_EU_10_smaller <- as.data.frame(prok_EULists$EU_10[,1:10])
rownames(prok_EU_10_smaller) <- rownames(prok_EULists$EU_10)
prok_EU_10_smaller_t <- as.data.frame(t(prok_EU_10_smaller))
prok_EU_10_smaller_t$Meter <- c("10", "20", "30", "40", "50", "60",
                                "70", "80", "90", "100")

prok_EU_52_smaller <- as.data.frame(prok_EULists$EU_52[,1:10])
rownames(prok_EU_52_smaller) <- rownames(prok_EULists$EU_52)
prok_EU_52_smaller_t <- as.data.frame(t(prok_EU_52_smaller))
prok_EU_52_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                 "70", "80", "90", "100")

prok_EU_53N_smaller <- as.data.frame(prok_EULists$EU_53N[,1:10])
rownames(prok_EU_53N_smaller) <- rownames(prok_EULists$EU_53N)
prok_EU_53N_smaller_t <- as.data.frame(t(prok_EU_53N_smaller))
prok_EU_53N_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                 "70", "80", "90", "100")

prok_EU_53S_smaller <- as.data.frame(prok_EULists$EU_53S[,1:10])
rownames(prok_EU_53S_smaller) <- rownames(prok_EULists$EU_53S)
prok_EU_53S_smaller_t <- as.data.frame(t(prok_EU_53S_smaller))
prok_EU_53S_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                 "70", "80", "90", "100")

prok_EU_54S_smaller <- as.data.frame(prok_EULists$EU_54S[,1:10])
rownames(prok_EU_54S_smaller) <- rownames(prok_EULists$EU_54S)
prok_EU_54S_smaller_t <- as.data.frame(t(prok_EU_54S_smaller))
prok_EU_54S_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                 "70", "80", "90", "100")

prok_EU_8_smaller <- as.data.frame(prok_EULists$EU_8[,1:10])
rownames(prok_EU_8_smaller) <- rownames(prok_EULists$EU_8)
prok_EU_8_smaller_t <- as.data.frame(t(prok_EU_8_smaller))
prok_EU_8_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                 "70", "80", "90", "100")

# Make all of these into a list
prok_EUsLists <- list(prok_EU_10_smaller_t, prok_EU_52_smaller_t, prok_EU_53N_smaller_t, prok_EU_53S_smaller_t, 
                         prok_EU_54S_smaller_t, prok_EU_8_smaller_t)
length(prok_EUsLists) #6, as expected
names(prok_EUsLists) <- c("EU_10", "EU_52", "EU_53N", "EU_53S", "EU_54S", "EU_8")

# Calculate Bray-Curtis dissimilarities within EUs (as created above)
proks_meanASVs_BCs <- vector("list", 6) #pre-allocate vector to store Bray-Curtis dissimilarities
for (i in 1:length(prok_EUsLists)){ #there are 6 dataframes within this, corresponding to each EU. 
  # they are in numerical order, starting with EU 10
  # Get Bray-Curtis dissimilarities for a subset (not last three rows)
  proks_meanASVs_BCs[[i]] <- as.matrix(vegdist(prok_EUsLists[[i]][,1:ncol(prok_EUsLists[[i]])-1], method="bray")) #all but meter column
  diag(proks_meanASVs_BCs[[i]]) <- NA
  rownames(proks_meanASVs_BCs[[i]]) <- rownames(prok_EUsLists[[i]]) #keep rownames (for whatever reason it won't work automatically)
  colnames(proks_meanASVs_BCs[[i]]) <- rownames(prok_EUsLists[[i]]) #keep rownames (for whatever reason it won't work automatically)
}
names(proks_meanASVs_BCs) <- names(prok_EUsLists) #make names the same as EULists (works because for loop above runs through EU Lists)

# Get 10 m and 100 meter Bray-Curtis distances
## Comparisons with both 10m and 100m
prok_meanComps <- vector("list", 6) #pre-allocate space for all comps with 10 m for each EU
for (i in 1:length(prok_meanComps)){
  #Below, pre-allocate matrix as each element in the list
  prok_meanComps[[i]] <- matrix(data=NA, nrow= nrow(proks_meanASVs_BCs[[i]]), ncol=4)
  for (j in 1:nrow(prok_meanComps[[1]])) {
    #Make first column in comps10 "Meter" from EULists
    prok_meanComps[[i]][j,1] <- prok_EUsLists[[i]]$Meter[j]
    prok_meanComps[[i]][j,2] <- proks_meanASVs_BCs[[i]][j,1] #so this is each row, adding in col 1 comparisons (will have some NAs for comparison with 10m)
    #Make third column 100m comps
    prok_meanComps[[i]][j,3] <- proks_meanASVs_BCs[[i]][j,10]
    #Make fourth column in comps10 EU
    prok_meanComps[[i]][j,4] <- names(prok_EUsLists)[[i]]
  }
  colnames(prok_meanComps[[i]]) <- c("meter", "comp10m", "comp100m", "EU")
}
names(prok_meanComps) <- names(proks_meanASVs_BCs)

# Use a for loop to arrange data in a way compatible with ggplot2
prok_compsForPlot <- as.data.frame(matrix(data=NA, nrow= 120, ncol=4))
for (i in 1:length(prok_meanComps)) {
  prok_meanComps[[i]] <- as.data.frame(prok_meanComps[[i]]) %>% 
    pivot_longer(comp10m:comp100m,
                 names_to= "compType",
                 values_to = "BCdists")
}

# Combine all of the data!
prok_BCcompsPlot <- bind_rows(prok_meanComps$EU_10, prok_meanComps$EU_52, prok_meanComps$EU_53N, 
                              prok_meanComps$EU_53S, prok_meanComps$EU_54S, prok_meanComps$EU_8)
prok_BCcompsPlot$meter <- as.numeric(prok_BCcompsPlot$meter) 
prok_BCcompsPlot$BCdists <- as.numeric(prok_BCcompsPlot$BCdists)

### *NEED TO MAKE SURE THAT IT IS WORKING AS EXPECTED since the code was co-opted from DissimilaritiesAcrosTransectsMedians.R

#. PLOTTING USING DATAFRAME CREATED ABOVE-- BUT MAKE COMPARISONS ONLY WITH 100M INTO FOREST

# EU 10
# comparison with forest 100m only
# Matrix comparison lines ONLY (added to ESA 2022 presentation)
# get subset of data for only forest comparison
EU_10_prok_BCsforestOnly <- prok_BCcompsPlot[1:20,] %>% 
  filter(compType == "comp100m")

ggEU_10_prok_forestOnly <- ggplot(data=EU_10_prok_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 10 (prokaryotes)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")
# quartz()
ggEU_10_prok_forestOnly 

# EU 52
EU_52_prok_BCsforestOnly <- prok_BCcompsPlot[21:40,] %>% 
  filter(compType == "comp100m")

ggEU_52_prok_forestOnly <- ggplot(data=EU_52_prok_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 52 (prokaryotes)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")
# quartz()
ggEU_52_prok_forestOnly

# EU 53N
EU_53N_prok_BCsforestOnly <- prok_BCcompsPlot[41:60,] %>% 
  filter(compType == "comp100m")

ggEU_53N_prok_forestOnly <- ggplot(data=EU_53N_prok_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 53N (prokaryotes)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")

# EU 53S
EU_53S_prok_BCsforestOnly <- prok_BCcompsPlot[61:80,] %>% 
  filter(compType == "comp100m")

ggEU_53S_prok_forestOnly <- ggplot(data=EU_53S_prok_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 53S (prokaryotes)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")


# EU 54S
EU_54S_prok_BCsforestOnly <- prok_BCcompsPlot[81:100,] %>% 
  filter(compType == "comp100m")

ggEU_54S_prok_forestOnly <- ggplot(data=EU_54S_prok_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 54S (prokaryotes)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")

# EU 8
EU_8_prok_BCsforestOnly <- prok_BCcompsPlot[101:120,] %>% 
  filter(compType == "comp100m")

ggEU_8_prok_forestOnly <- ggplot(data=EU_8_prok_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 8 (prokaryotes)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")

##### PLOT ALL TOGETHER #####
quartz()
grid.arrange(ggEU_10_prok_forestOnly, ggEU_52_prok_forestOnly, ggEU_53N_prok_forestOnly, ggEU_53S_prok_forestOnly, ggEU_54S_prok_forestOnly, ggEU_8_prok_forestOnly, ncol=3)

###############################################
# FUNGI
################################################

###################################
# 1. Getting dissimilarity values for all transects, all EUs

# First, get only the ASVs that are forest specialists
fung_forestNames <- as.character(diffAbunDat_tidy_FUNGI$ASV_name[which(diffAbunDat_tidy_FUNGI$Habitat=="forest")])
prokaryotes_forest_DiffAbund.ps <- prune_taxa(fung_forestNames, fungiDiffAbund.ps)

# To get dissimilarities, first split phyloseq object by EUs and remove ASVs that don't occur in subset
EU_52_fung.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_52")
EU_52_fung.ps <- prune_taxa(taxa_sums(EU_52_fung.ps) > 0, EU_52_fung.ps) #remove non-occuring ASVs

EU_53N_fung.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_53N")
EU_53N_fung.ps <- prune_taxa(taxa_sums(EU_53N_fung.ps) > 0, EU_53N_fung.ps) #remove non-occuring ASVs

EU_54S_fung.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_54S")
EU_54S_fung.ps <- prune_taxa(taxa_sums(EU_54S_fung.ps) > 0, EU_54S_fung.ps) #remove non-occuring ASVs

EU_8_fung.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_8")
EU_8_fung.ps <- prune_taxa(taxa_sums(EU_8_fung.ps) > 0, EU_8_fung.ps) #remove non-occuring ASVs

EU_53S_fung.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_53S")
EU_53S_fung.ps <- prune_taxa(taxa_sums(EU_53S_fung.ps) > 0, EU_53S_fung.ps) #remove non-occuring ASVs

EU_10_fung.ps <- subset_samples(prokaryotes_forest_DiffAbund.ps, EU == "EU_10")
EU_10_fung.ps <- prune_taxa(taxa_sums(EU_10_fung.ps) > 0, EU_10_fung.ps) #remove non-occuring ASVs

# Run function on each EU's phyloseq object and then make sure first two columns are numeric
fung_dissim_52 <- transectDissim3(EU_52_fung.ps)
fung_dissim_52[,1] <- as.numeric(fung_dissim_52[,1])
fung_dissim_52[,2] <- as.numeric(fung_dissim_52[,2])

fung_dissim_53N <- transectDissim3(EU_53N_fung.ps)
fung_dissim_53N[,1] <- as.numeric(fung_dissim_53N[,1])
fung_dissim_53N[,2] <- as.numeric(fung_dissim_53N[,2])

fung_dissim_54S <- transectDissim3(EU_54S_fung.ps)
fung_dissim_54S[,1] <- as.numeric(fung_dissim_54S[,1])
fung_dissim_54S[,2] <- as.numeric(fung_dissim_54S[,2])

fung_dissim_8 <- transectDissim3(EU_8_fung.ps)
fung_dissim_8[,1] <- as.numeric(fung_dissim_8[,1])
fung_dissim_8[,2] <- as.numeric(fung_dissim_8[,2])

fung_dissim_53S <- transectDissim3(EU_53S_fung.ps)
fung_dissim_53S[,1] <- as.numeric(fung_dissim_53S[,1])
fung_dissim_53S[,2] <- as.numeric(fung_dissim_53S[,2])

fung_dissim_10 <- transectDissim3(EU_10_fung.ps)
fung_dissim_10[,1] <- as.numeric(fung_dissim_10[,1])
fung_dissim_10[,2] <- as.numeric(fung_dissim_10[,2])

#############
# Plotting
# Function above spits this out as it needs to be to plot using ggplot2. 
# I removed comparisons with 10m by commenting them out

fung_ggEU_52 <- ggplot() + 
  #geom_line(data=fung_dissim_52[1:38,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=fung_dissim_52[39:76,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 52") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
#quartz()
fung_ggEU_52
fung_ggEU_53N <- ggplot() + 
  #geom_line(data=fung_dissim_53N[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=fung_dissim_53N[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 53N") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
fung_ggEU_53N
#quartz()
fung_ggEU_54S <- ggplot() + 
  #geom_line(data=fung_dissim_54S[1:38,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=fung_dissim_54S[39:76,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 54S") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
fung_ggEU_54S
#quartz()
fung_ggEU_8 <- ggplot() + 
  #geom_line(data=fung_dissim_8[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=fung_dissim_8[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 8") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
fung_ggEU_8
#quartz()
fung_ggEU_53S <- ggplot() + 
  #geom_line(data=fung_dissim_53S[1:40,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=fung_dissim_53S[41:80,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 53S") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
fung_ggEU_53S
#quartz()
fung_ggEU_10 <- ggplot() + 
  #geom_line(data=fung_dissim_10[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=fung_dissim_10[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 10") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
fung_ggEU_10 

quartz()
require(gridExtra)
grid.arrange(fung_ggEU_52, fung_ggEU_53N, fung_ggEU_54S, fung_ggEU_8, fung_ggEU_53S, fung_ggEU_10, ncol=3)

# 2. Getting dissimilarity values where MEAN is that used to calculate z-scores
# (i.e. each differentially abundant ASV at each meter)

fungiASVmeans_AllEUs #made in fungiEdgeEffectsZscores.R
colnames(fungiASVmeans_AllEUs)

# 1. CALCULATE BRAY-CURTIS DISSIMILARITIES ALONG TRANSECTS WITHIN EACH EU

# The lines below get the Bray-Curtis dissimilarities along the transect, comparing each point along
# the transect to 10 meters and to 100 meters 

# The piped code below creates a list of 6 tibbles, corresponding to each EU
fung_EULists <- fungiASVmeans_AllEUs %>%  #make multiple datasets based on EU
  dplyr::group_by(EU) %>%  #group by EU
  dplyr::group_split() 
names(fung_EULists) <- sort(unique(fungiASVmeans_AllEUs$EU))

# Use a for loop to make ASV_names row names (have to do this step after splitting by EUs so that there are no duplicated row names)
for (i in 1:length(fung_EULists)) {
  fung_EULists[[i]] <- fung_EULists[[i]] %>% 
    column_to_rownames("ASV_name")  #make sampNames column in the rownames
}

# For vegdist, rows need to be samples and columns should be species.
# 1) pull out only mean ASV abundances, 2) make ASV names the rownames,
# and 3) flip it all!
fung_EU_10_smaller <- as.data.frame(fung_EULists$EU_10[,1:10])
rownames(fung_EU_10_smaller) <- rownames(fung_EULists$EU_10)
fung_EU_10_smaller_t <- as.data.frame(t(fung_EU_10_smaller))
fung_EU_10_smaller_t$Meter <- c("10", "20", "30", "40", "50", "60",
                                "70", "80", "90", "100")

fung_EU_52_smaller <- as.data.frame(fung_EULists$EU_52[,1:10])
rownames(fung_EU_52_smaller) <- rownames(fung_EULists$EU_52)
fung_EU_52_smaller_t <- as.data.frame(t(fung_EU_52_smaller))
fung_EU_52_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                 "70", "80", "90", "100")

fung_EU_53N_smaller <- as.data.frame(fung_EULists$EU_53N[,1:10])
rownames(fung_EU_53N_smaller) <- rownames(fung_EULists$EU_53N)
fung_EU_53N_smaller_t <- as.data.frame(t(fung_EU_53N_smaller))
fung_EU_53N_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                  "70", "80", "90", "100")

fung_EU_53S_smaller <- as.data.frame(fung_EULists$EU_53S[,1:10])
rownames(fung_EU_53S_smaller) <- rownames(fung_EULists$EU_53S)
fung_EU_53S_smaller_t <- as.data.frame(t(fung_EU_53S_smaller))
fung_EU_53S_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                  "70", "80", "90", "100")

fung_EU_54S_smaller <- as.data.frame(fung_EULists$EU_54S[,1:10])
rownames(fung_EU_54S_smaller) <- rownames(fung_EULists$EU_54S)
fung_EU_54S_smaller_t <- as.data.frame(t(fung_EU_54S_smaller))
fung_EU_54S_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                  "70", "80", "90", "100")

fung_EU_8_smaller <- as.data.frame(fung_EULists$EU_8[,1:10])
rownames(fung_EU_8_smaller) <- rownames(fung_EULists$EU_8)
fung_EU_8_smaller_t <- as.data.frame(t(fung_EU_8_smaller))
fung_EU_8_smaller_t$Meter <-  c("10", "20", "30", "40", "50", "60",
                                "70", "80", "90", "100")

# Make all of these into a list
fung_EUsLists <- list(fung_EU_10_smaller_t, fung_EU_52_smaller_t, fung_EU_53N_smaller_t, fung_EU_53S_smaller_t, 
                      fung_EU_54S_smaller_t, fung_EU_8_smaller_t)
length(fung_EUsLists) #6, as expected
names(fung_EUsLists) <- c("EU_10", "EU_52", "EU_53N", "EU_53S", "EU_54S", "EU_8")

# Calculate Bray-Curtis dissimilarities within EUs (as created above)
fungs_meanASVs_BCs <- vector("list", 6) #pre-allocate vector to store Bray-Curtis dissimilarities
for (i in 1:length(fung_EUsLists)){ #there are 6 dataframes within this, corresponding to each EU. 
  # they are in numerical order, starting with EU 10
  # Get Bray-Curtis dissimilarities for a subset (not last three rows)
  fungs_meanASVs_BCs[[i]] <- as.matrix(vegdist(fung_EUsLists[[i]][,1:ncol(fung_EUsLists[[i]])-1], method="bray")) #all but meter column
  diag(fungs_meanASVs_BCs[[i]]) <- NA
  rownames(fungs_meanASVs_BCs[[i]]) <- rownames(fung_EUsLists[[i]]) #keep rownames (for whatever reason it won't work automatically)
  colnames(fungs_meanASVs_BCs[[i]]) <- rownames(fung_EUsLists[[i]]) #keep rownames (for whatever reason it won't work automatically)
}
names(fungs_meanASVs_BCs) <- names(fung_EUsLists) #make names the same as EULists (works because for loop above runs through EU Lists)

# Get 10 m and 100 meter Bray-Curtis distances
## Comparisons with both 10m and 100m
fung_meanComps <- vector("list", 6) #pre-allocate space for all comps with 10 m for each EU
for (i in 1:length(fung_meanComps)){
  #Below, pre-allocate matrix as each element in the list
  fung_meanComps[[i]] <- matrix(data=NA, nrow= nrow(fungs_meanASVs_BCs[[i]]), ncol=4)
  for (j in 1:nrow(fung_meanComps[[1]])) {
    #Make first column in comps10 "Meter" from EULists
    fung_meanComps[[i]][j,1] <- fung_EUsLists[[i]]$Meter[j]
    fung_meanComps[[i]][j,2] <- fungs_meanASVs_BCs[[i]][j,1] #so this is each row, adding in col 1 comparisons (will have some NAs for comparison with 10m)
    #Make third column 100m comps
    fung_meanComps[[i]][j,3] <- fungs_meanASVs_BCs[[i]][j,10]
    #Make fourth column in comps10 EU
    fung_meanComps[[i]][j,4] <- names(fung_EUsLists)[[i]]
  }
  colnames(fung_meanComps[[i]]) <- c("meter", "comp10m", "comp100m", "EU")
}
names(fung_meanComps) <- names(fungs_meanASVs_BCs)

# Use a for loop to arrange data in a way compatible with ggplot2
fung_compsForPlot <- as.data.frame(matrix(data=NA, nrow= 120, ncol=4))
for (i in 1:length(fung_meanComps)) {
  fung_meanComps[[i]] <- as.data.frame(fung_meanComps[[i]]) %>% 
    pivot_longer(comp10m:comp100m,
                 names_to= "compType",
                 values_to = "BCdists")
}

# Combine all of the data!
fung_BCcompsPlot <- bind_rows(fung_meanComps$EU_10, fung_meanComps$EU_52, fung_meanComps$EU_53N, 
                              fung_meanComps$EU_53S, fung_meanComps$EU_54S, fung_meanComps$EU_8)
fung_BCcompsPlot$meter <- as.numeric(fung_BCcompsPlot$meter) 
fung_BCcompsPlot$BCdists <- as.numeric(fung_BCcompsPlot$BCdists)

### *NEED TO MAKE SURE THAT IT IS WORKING AS EXPECTED since the code was co-opted from DissimilaritiesAcrosTransectsMedians.R

#. PLOTTING USING DATAFRAME CREATED ABOVE-- BUT MAKE COMPARISONS ONLY WITH 100M INTO FOREST

# EU 10
# comparison with forest 100m only
# Matrix comparison lines ONLY (added to ESA 2022 presentation)
# get subset of data for only forest comparison
EU_10_fung_BCsforestOnly <- fung_BCcompsPlot[1:20,] %>% 
  filter(compType == "comp100m")

ggEU_10_fung_forestOnly <- ggplot(data=EU_10_fung_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 10 (fungi)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")
# quartz()
ggEU_10_fung_forestOnly 

# EU 52
EU_52_fung_BCsforestOnly <- fung_BCcompsPlot[21:40,] %>% 
  filter(compType == "comp100m")

ggEU_52_fung_forestOnly <- ggplot(data=EU_52_fung_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 52 (fungi)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")
# quartz()
ggEU_52_fung_forestOnly

# EU 53N
EU_53N_fung_BCsforestOnly <- fung_BCcompsPlot[41:60,] %>% 
  filter(compType == "comp100m")

ggEU_53N_fung_forestOnly <- ggplot(data=EU_53N_fung_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 53N (fungi)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")

# EU 53S
EU_53S_fung_BCsforestOnly <- fung_BCcompsPlot[61:80,] %>% 
  filter(compType == "comp100m")

ggEU_53S_fung_forestOnly <- ggplot(data=EU_53S_fung_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 53S (fungi)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")


# EU 54S
EU_54S_fung_BCsforestOnly <- fung_BCcompsPlot[81:100,] %>% 
  filter(compType == "comp100m")

ggEU_54S_fung_forestOnly <- ggplot(data=EU_54S_fung_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 54S (fungi)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")

# EU 8
EU_8_fung_BCsforestOnly <- fung_BCcompsPlot[101:120,] %>% 
  filter(compType == "comp100m")

ggEU_8_fung_forestOnly <- ggplot(data=EU_8_fung_BCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title= "EU 8 (fungi)", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")

##### PLOT ALL TOGETHER #####
quartz()
grid.arrange(ggEU_52_fung_forestOnly, ggEU_53N_fung_forestOnly, ggEU_54S_fung_forestOnly, 
             ggEU_8_fung_forestOnly, ggEU_53S_fung_forestOnly, ggEU_10_fung_forestOnly, ncol=3)