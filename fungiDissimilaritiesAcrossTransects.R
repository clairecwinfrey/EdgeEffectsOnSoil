# FUNGI- Bray-Curtis Dissimilarities within EUs Across Transects
# Exploration of Dissimilarity Patterns 
# Feb. 12, 2022

# This script calculates and plots the Bray-Curtis dissimilarities for the 
# meters across the transect. It uses the post-ubiquity (not median) abundances
# of each ASV for each transect. See description below, where ITS_postUbiquity.ps is
# loaded. In addition, it performs db-RDAs and forward model selection to test which
# environmental variables explain differences among samples

# Libraries
library("phyloseq")
library("ggplot2")      
library("dplyr")        
library("tibble")       
library("tidyr")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("ggord") #plotting dbRDAs

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch")

load("RobjectsSaved/ITS_postUbiquity.ps") #load ITS phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45. ITS_postUbiquity.ps was made and 
# saved in the script titled "ITS_UbiquityMedianSetup.R".

#### *** IMPORTANT: SHOULD LOAD IN ITS_postUbiquity.ps FOR THIS?! *****
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

# 4. TransectDissim3 obtains the Bray-Curtis dissimilarities between 1) 10 m and 
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
  metaDf <- pssd2veg(physeq) #get metadata using function defined earlier; samples are rows
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
# Getting dissimilarity values for all transects, all EUs
#############

# To get dissimilarities, first split phyloseq object by EUs and remove ASVs that don't occur in subset
EU_52_ITS_Soils.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_52")
EU_52_ITS_Soils.ps <- prune_taxa(taxa_sums(EU_52_ITS_Soils.ps) > 0, EU_52_ITS_Soils.ps) #remove non-occurring ASVs

EU_53N_ITS_Soils.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_53N")
EU_53N_ITS_Soils.ps <- prune_taxa(taxa_sums(EU_53N_ITS_Soils.ps) > 0, EU_53N_ITS_Soils.ps) #remove non-occurring  ASVs

EU_54S_ITS_Soils.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_54S")
EU_54S_ITS_Soils.ps <- prune_taxa(taxa_sums(EU_54S_ITS_Soils.ps) > 0, EU_54S_ITS_Soils.ps) #remove non-occurring ASVs

EU_8_ITS_Soils.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_8")
EU_8_ITS_Soils.ps <- prune_taxa(taxa_sums(EU_8_ITS_Soils.ps) > 0, EU_8_ITS_Soils.ps) #remove non-occurring ASVs

EU_53S_ITS_Soils.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_53S")
EU_53S_ITS_Soils.ps <- prune_taxa(taxa_sums(EU_53S_ITS_Soils.ps) > 0, EU_53S_ITS_Soils.ps) #remove non-occurring  ASVs

EU_10_ITS_Soils.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_10")
EU_10_ITS_Soils.ps <- prune_taxa(taxa_sums(EU_10_ITS_Soils.ps) > 0, EU_10_ITS_Soils.ps) #remove non-occurring  ASVs

# Run function on each EU's phyloseq object and then make sure first two columns are numeric
dissim_ITS_52 <- transectDissim3(EU_52_ITS_Soils.ps)
dissim_ITS_52[,1] <- as.numeric(dissim_ITS_52[,1])
dissim_ITS_52[,2] <- as.numeric(dissim_ITS_52[,2])

dissim_ITS_53N <- transectDissim3(EU_53N_ITS_Soils.ps)
dissim_ITS_53N[,1] <- as.numeric(dissim_ITS_53N[,1])
dissim_ITS_53N[,2] <- as.numeric(dissim_ITS_53N[,2])

dissim_ITS_54S <- transectDissim3(EU_54S_ITS_Soils.ps)
dissim_ITS_54S[,1] <- as.numeric(dissim_ITS_54S[,1])
dissim_ITS_54S[,2] <- as.numeric(dissim_ITS_54S[,2])

dissim_ITS_8 <- transectDissim3(EU_8_ITS_Soils.ps)
dissim_ITS_8[,1] <- as.numeric(dissim_ITS_8[,1])
dissim_ITS_8[,2] <- as.numeric(dissim_ITS_8[,2])

dissim_ITS_53S <- transectDissim3(EU_53S_ITS_Soils.ps)
dissim_ITS_53S[,1] <- as.numeric(dissim_ITS_53S[,1])
dissim_ITS_53S[,2] <- as.numeric(dissim_ITS_53S[,2])

dissim_ITS_10 <- transectDissim3(EU_10_ITS_Soils.ps)
dissim_ITS_10[,1] <- as.numeric(dissim_ITS_10[,1])
dissim_ITS_10[,2] <- as.numeric(dissim_ITS_10[,2])

#############
# Plotting
# Function above spits this out as it needs to be to plot using ggplot2. 

ggEU_52 <- ggplot() + 
  geom_line(data=dissim_ITS_52[1:38,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=dissim_ITS_52[39:76,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 52") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
ggEU_52
#quartz()
ggEU_53N <- ggplot() + 
  geom_line(data=dissim_ITS_53N[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=dissim_ITS_53N[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 53N") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
ggEU_53N
#quartz()
ggEU_54S <- ggplot() + 
  geom_line(data=dissim_ITS_54S[1:38,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=dissim_ITS_54S[39:76,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 54S") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
ggEU_54S
#quartz()
ggEU_8 <- ggplot() + 
  geom_line(data=dissim_ITS_8[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=dissim_ITS_8[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 8") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
ggEU_8
#quartz()
ggEU_53S <- ggplot() + 
  geom_line(data=dissim_ITS_53S[1:40,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=dissim_ITS_53S[41:80,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 53S") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
ggEU_53S

#quartz()
ggEU_10 <- ggplot() + 
  geom_line(data=dissim_ITS_10[1:39,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "goldenrod") +
  geom_line(data=dissim_ITS_10[40:78,], aes(x=meter, y=`Bray-Curtis`, group = transect), color = "darkgreen") +
  theme_bw() + ylim(0.3, 1.00) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("EU 10") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey")
ggEU_10 

#quartz()
require(gridExtra)
grid.arrange(ggEU_52, ggEU_53N, ggEU_54S, ggEU_8, ggEU_53S, ggEU_10, ncol=3)

###############################################
# DBRDA
################################################
# first, remove samples with missing pH values... since ordiR2Step won't work with them...
postUbiquitynoNA_ITS.ps <- subset_samples(ITS_postUbiquity.ps, Sample.ID != "53ND_T_60")
postUbiquitynoNA_ITS.ps <- subset_samples(postUbiquitynoNA_ITS.ps, Sample.ID != "53ND_L_100")

postUbiqnoNA_ITS_ASVs <- ASVs_outta_ps(postUbiquitynoNA_ITS.ps)

# Get Bray-Curtis dissimilarities
postUbiqnoNAs_ITS_BCdists <- vegdist(postUbiqnoNA_ITS_ASVs, method = "bray")

postUbiqNoNAsampsMeta <- pssd2veg(postUbiquitynoNA_ITS.ps)
#View(postUbiqNoNAsampsMeta)


## Make a db-RDA and test it
postUbiq_ITS_dbRDA.mod0 <- dbrda(postUbiqnoNAs_ITS_BCdists ~1, data= postUbiqNoNAsampsMeta)
postUbiq_ITS_dbRDA.full <- dbrda(postUbiqnoNAs_ITS_BCdists ~ mean_pH + Percent_Vegetation_Cover + 
                               meanDens + Habitat + Condition(EU), data= postUbiqNoNAsampsMeta)
set.seed(19)
postUbiq_ITS_dbRDAresults <- anova.cca(postUbiq_ITS_dbRDA.full, permutations = 9999)  #test significance of constraints
# Model: dbrda(formula = postUbiqnoNAs_ITS_BCdists ~ mean_pH + Percent_Vegetation_Cover + meanDens + Habitat + Condition(EU), data = postUbiqNoNAsampsMeta)
# Df SumOfSqs      F Pr(>F)    
# Model      5   11.452 8.9237  1e-04 ***
#   Residual 220   56.465   

set.seed(19)
postUbiq_ITS_dbRDAresultsByTerms <- anova.cca(postUbiq_ITS_dbRDA.full, by="margin", permutations = 9999) #test which variables are significant
# # Model: dbrda(formula = postUbiqnoNAs_ITS_BCdists ~ mean_pH + Percent_Vegetation_Cover + meanDens + Habitat + Condition(EU), data = postUbiqNoNAsampsMeta)
# Df SumOfSqs      F Pr(>F)    
# mean_pH                    1    1.319 5.1398 0.0001 ***
#   Percent_Vegetation_Cover   1    0.534 2.0802 0.0063 ** 
#   meanDens                   1    1.211 4.7165 0.0001 ***
#   Habitat                    2    1.365 2.6598 0.0001 ***
#   Residual                 220   56.465  

# SAVE RESULTS
#save(postUbiq_ITS_dbRDAresults, file="RobjectsSaved/postUbiq_ITS_dbRDAresults") #saved Feb. 12, 2023
#save(postUbiq_ITS_dbRDAresultsByTerms, file="RobjectsSaved/postUbiq_ITS_dbRDAresultsByTerms") #saved Feb. 12, 2023

# model selection
set.seed(93)
postUbiqITS_Mod_forsel <- ordiR2step(object=postUbiq_ITS_dbRDA.mod0, scope = postUbiq_ITS_dbRDA.full, permutations = 9999) #bumped down to 9999 because too many permutations kept messing up R 
postUbiqITS_Mod_forsel$anova
# Shows that canopy cover is the most important variable!!
# R2.adj Df    AIC      F Pr(>F)    
# + meanDens      0.10824  1 987.77 28.917  1e-04 ***
#   <All variables> 0.12993

# SAVE RESULTS
#save(postUbiqITS_Mod_forsel, file="RobjectsSaved/postUbiqITS_Mod_forsel") #saved Feb. 12, 2023

# canopy only
postUbiq_ITS_dbRDA.canopyOnly <- dbrda(postUbiqnoNAs_ITS_BCdists ~ 
                                    meanDens, data= postUbiqNoNAsampsMeta)
set.seed(93)
postUbiq_ITS_dbRDAcanopyOnlyresults <- anova.cca(postUbiq_ITS_dbRDA.canopyOnly, permutations = 9999)  #test significance of constraints

# canopy and habitat
postUbiq_ITS_dbRDA.canopyHabitatOnly <- dbrda(postUbiqnoNAs_ITS_BCdists ~ 
                                            meanDens + Habitat, data= postUbiqNoNAsampsMeta)
set.seed(93)
postUbiq_ITS_dbRDAcanopyHabitatOnly <- anova.cca(postUbiq_ITS_dbRDA.canopyOnly, permutations = 9999) 


cor(postUbiqNoNAsampsMeta$meanDens, postUbiqNoNAsampsMeta$mean_pH) #correlation is -0.318194


#quartz()
plot(postUbiq_ITS_dbRDA.canopyOnly)

quartz()
postUbiqCanHabPlot <- ggord(postUbiq_ITS_dbRDA.canopyOnly, postUbiqNoNAsampsMeta$Habitat, size =2)
postUbiqCanHabPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 20)) + theme(legend.position = "bottom") + theme(axis.text = element_text(size= 15)) + theme(axis.title = element_text(size = 17)) + theme(axis.title.x = element_text(vjust= -1.5)) +
  theme(axis.title.y = element_text(vjust = 1.5)) + theme(legend.box.margin = margin(10, -10, -10, -10))
