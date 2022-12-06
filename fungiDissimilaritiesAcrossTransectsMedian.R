# Bray-Curtis Dissimilarities - Plots and RDAs within EUs- Median Abundances 
# December 5, 2022

# This script calculates and plots the Bray-Curtis dissimilarities for the 
# meters across the transect. In addition, this script runs distance-based redundancy
# analyses (db_RDA) and tests their constraints with permutational ANOVAs. 
# THIS IS USES MEDIAN ABUNDANCES (see description of ITS_medianEU.ps object below)

# Libraries
library("phyloseq")
library("ggplot2")      
library("dplyr")        
library("tibble")       
library("tidyr")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Load data
load("RobjectsSaved/ITS_medianEU.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45, and 4) finally getting the 
# median ASV abundance at each meter in each EU (median is four values from each
# transect)
ITS_medianEU.ps

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

###############################################################################
# Main part of script

############################################################
# FOR MEDIAN ABUNDANCES (i.e. ITS_medianEU.ps)
################################################################
# 1. CALCULATE BRAY-CURTIS DISSIMILARITIES ALONG TRANSECTS WITHIN EACH EU

# The lines below get the Bray-Curtis dissimilarities along the transect, comparing each point along
# the transect to 10 meters and to 100 meters 
ASVsdf <- ASVs_outta_ps(ITS_medianEU.ps) #get ASV table using functioned defined earlier; samples are rows, ASVs = columns
metaDf <- metadata_outta_ps(ITS_medianEU.ps) #get metadata using function defined earlier; samples are rows
# Next six lines make a new df that has variables of interest and ASVs
ASVsdf$EU <- metaDf$EU
ASVsdf$Meter <- metaDf$meter
ASVsdf$sampNames <- rownames(ASVsdf) #preserve these as rownames so that they remain in next step

# The piped code below creates a list of 6 tibbles, corresponding to each EU
EULists <- ASVsdf %>%  #make multiple datasets based on EU
  dplyr::group_by(EU) %>%  #group by EU
  dplyr::group_split() 
names(EULists) <- sort(unique(ASVsdf$EU))

# Use a for loop to arrange by meter and make sampNames row names (because code below doesn't work on list of lists)
for (i in 1:length(EULists)) {
  EULists[[i]] <- EULists[[i]] %>% 
    arrange(Meter) %>% #sort in ascending order
    column_to_rownames("sampNames") #make sampNames column in the rownames 
}

# Calculate Bray-Curtis dissimilarities within EUs (as created above)
BCs <- vector("list", 6) #pre-allocate vector to store Bray-Curtis dissimilarities
for (i in 1:length(EULists)){ #there are 6 tibbles within this, corresponding to each EU. 
  # they are in numerical order, starting with EU 10
  # Get Bray-Curtis dissimilarities for a subset (not last three rows)
  BCs[[i]] <- as.matrix(vegdist(EULists[[i]][,1:(ncol(EULists[[i]])-3)], method="bray"))
  diag(BCs[[i]]) <- NA
  rownames(BCs[[i]]) <- rownames(EULists[[i]]) #keep rownames (for whatever reason it won't work automatically)
  colnames(BCs[[i]]) <- rownames(EULists[[i]]) #keep rownames (for whatever reason it won't work automatically)
}
names(BCs) <- names(EULists) #make names the same as EULists (works because for loop above runs through EU Lists)

# Get 10 m and 100 meter Bray-Curtis distances
## Comparisons with both 10m and 100m
comps <- vector("list", 6) #pre-allocate space for all comps with 10 m for each EU
for (i in 1:length(comps)){
  #Below, pre-allocate matrix as each element in the list
  comps[[i]] <- matrix(data=NA, nrow= nrow(BCs[[i]]), ncol=4)
  for (j in 1:nrow(comps[[1]])) {
    #Make first column in comps10 "Meter" from EULists
    comps[[i]][j,1] <- EULists[[i]]$Meter[j]
    comps[[i]][j,2] <- BCs[[i]][j,1] #so this is each row, adding in col 1 comparisons (will have some NAs for comparison with 10m)
    #Make third column 100m comps
    comps[[i]][j,3] <- BCs[[i]][j,10]
    #Make fourth column in comps10 EU
    comps[[i]][j,4] <- EULists[[i]]$EU[j]
  }
  colnames(comps[[i]]) <- c("meter", "comp10m", "comp100m", "EU")
}
names(comps) <- names(BCs)

# Use a for loop to arrange data in a way compatible with ggplot2
compsForPlot <- as.data.frame(matrix(data=NA, nrow= 120, ncol=4))
for (i in 1:length(comps)) {
  comps[[i]] <- as.data.frame(comps[[i]]) %>% 
    pivot_longer(comp10m:comp100m,
                 names_to= "compType",
                 values_to = "BCdists")
}

# Combine all of the data!
BCcompsPlot <- bind_rows(comps$EU_10, comps$EU_52, comps$EU_53N, comps$EU_53S, comps$EU_54S, comps$EU_8)
BCcompsPlot$meter <- as.numeric(BCcompsPlot$meter) 
BCcompsPlot$BCdists <- as.numeric(BCcompsPlot$BCdists)
#View(BCcompsPlot)
##############################################

#. 2 PLOTTING USING DATAFRAME CREATED ABOVE
# New labels for transect
transectX <- c("40", "30","20", "10", "edge", "10", "20", "30", "40", "50")#for re-naming axis tick marks

# EU 10
# Matrix AND Patch lines
ggEU_10 <- ggplot(data=BCcompsPlot[1:20,], aes(x=meter, y=BCdists, color= compType)) + 
  geom_line(size=2.7) +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 10", y = "Bray-Curtis dissimilarity", size= 10) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey") + 
  theme(text = element_text(size=15)) +
  scale_colour_manual(name = "compType", values = c("darkgreen", "goldenrod")) +
  theme(legend.position = "none")
# quartz()
ggEU_10 

# Matrix comparison lines ONLY (added to ESA 2022 presentation)
# get subset of data for only forest comparison
EU_10_MedBCsforestOnly <- BCcompsPlot[1:20,] %>% 
  filter(compType == "comp100m")

ggEU_10_forestOnly <- ggplot(data=EU_10_MedBCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 10", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")
# quartz()
ggEU_10_forestOnly 

# EU 52
ggEU_52 <- ggplot(data=BCcompsPlot[21:40,], aes(x=meter, y=BCdists, color= compType)) + 
  geom_line(size=2.7) +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 52", y = "Bray-Curtis dissimilarity", size= 10) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey") + 
  theme(text = element_text(size=15)) +
  scale_colour_manual(name = "compType", values = c("darkgreen", "goldenrod")) +
  theme(legend.position = "none")
# quartz()
ggEU_52 

# EU 53N
# quartz()
ggEU_53N <- ggplot(data=BCcompsPlot[41:60,], aes(x=meter, y=BCdists, color= compType)) + 
  geom_line(size=2.7) +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 53N", y = "Bray-Curtis dissimilarity", size= 10) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey") + 
  theme(text = element_text(size=15)) +
  scale_colour_manual(name = "compType", values = c("darkgreen", "goldenrod")) +
  theme(legend.position = "none")
# quartz()
ggEU_53N 

# EU 53S
ggEU_53S <- ggplot(data=BCcompsPlot[61:80,], aes(x=meter, y=BCdists, color= compType)) + 
  geom_line(size=2.7) +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 53S", y = "Bray-Curtis dissimilarity", size= 10) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey") + 
  theme(text = element_text(size=15)) +
  scale_colour_manual(name = "compType", values = c("darkgreen", "goldenrod")) +
  theme(legend.position = "none")
# quartz()
ggEU_53S 

# Matrix comparison lines ONLY 
# get subset of data for only forest comparison
EU_53S_MedBCsforestOnly <- BCcompsPlot[61:80,] %>% 
  filter(compType == "comp100m")

ggEU_53S_forestOnly <- ggplot(data=EU_53S_MedBCsforestOnly, aes(x=meter, y=BCdists)) + 
  geom_line(size=2.7, color= "darkgreen") +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 53S", y = "Bray-Curtis dissimilarity", size= 14) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=1.5) + 
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")
# quartz()
ggEU_53S_forestOnly 

# EU 54S
ggEU_54S <- ggplot(data=BCcompsPlot[81:100,], aes(x=meter, y=BCdists, color= compType)) + 
  geom_line(size=2.7) +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 54S", y = "Bray-Curtis dissimilarity", size= 10) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey") + 
  theme(text = element_text(size=15)) +
  scale_colour_manual(name = "compType", values = c("darkgreen", "goldenrod")) +
  theme(legend.position = "none")
# quartz()
ggEU_54S 

# EU 8
ggEU_8 <- ggplot(data=BCcompsPlot[101:120,], aes(x=meter, y=BCdists, color= compType)) + 
  geom_line(size=2.7) +
  theme_bw() + ylim(0.1, 1.00) +
  scale_x_continuous("meters from edge", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), labels=transectX) +
  labs(title= "EU 8", y = "Bray-Curtis dissimilarity", size= 10) + 
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey") + 
  theme(text = element_text(size=15)) +
  scale_colour_manual(name = "compType", values = c("darkgreen", "goldenrod")) +
  theme(legend.position = "none")
# quartz()
ggEU_8 

##### PLOT ALL TOGETHER #####
quartz()
grid.arrange(ggEU_52, ggEU_53N, ggEU_54S, ggEU_8, ggEU_53S, ggEU_10, ncol=3)

# PLOT FOREST COMPARISONS for EU 10 and EU 53 south side by side
quartz()
grid.arrange(ggEU_10_forestOnly, ggEU_53S_forestOnly, ncol=2)

#####################
# FOR load("RobjectsSaved/ITS_medianEU.ps")


################################################
# STATS!
################################################

ASVsdf <- ASVsdf %>% 
  arrange(EU, Meter)
dim(ASVsdf)

#View(ASVsdf[26,]) #cut out those with NAs for pH
#View(ASVsdf[30,]) #cut out those with NAs for pH
rownames(ASVsdf)

medBCdists <- vegdist(ASVsdf[c(1:25, 27:29, 31:60),1:413], method = "bray")
medBCdists

############### HERE ##################
load("pHCanVegMedEU") # bring in pHCanVegMedEU from EnvDataPlots.R
pHCanVegMedEU <- as.data.frame(pHCanVegMedEU[1:60,])
#View(pHCanVegMedEU)
pHCanVegMedEU$edgeDist <- abs(50-pHCanVegMedEU$Meter)

rownames(pHCanVegMedEU) <- rownames(ASVsdf)

set.seed(19)
medMod_dbRDA.mod0 <- dbrda(medBCdists ~1)
medMod_dbRDA <- dbrda(medBCdists ~ pHCanVegMedEU[c(1:25, 27:29, 31:60),3] + pHCanVegMedEU[c(1:25, 27:29, 31:60),4] + 
                        pHCanVegMedEU[c(1:25, 27:29, 31:60),5] + pHCanVegMedEU[c(1:25, 27:29, 31:60),6] +
                        Condition(pHCanVegMedEU[c(1:25, 27:29, 31:60),1]))
# col 3 is pH, col 4 is canopy cover, col 5 is vegetation cover, and col 6 is distance from edge
medMod_dbRDARes <- anova.cca(medMod_dbRDA, permutations = 9999)  #test significance of constraints
medMod_dbRDA$terms

#quartz()
plot(medMod_dbRDA)

# col 3 is pH, col 4 is canopy cover, col 5 is vegetation cover, and col 6 is distance from edge  
set.seed(93)
medMod_forsel <- ordiR2step(medMod_dbRDA.mod0, scope = medMod_dbRDA, permutations = 999)

medMod_forsel$anova #by far canopy cover best variable predictotr

colnames(pHCanVegMedEU[c(3,4)])
# Are median pHs and canopy cover correlated? Yes, yes they are!
cor(pHCanVegMedEU[c(1:25, 27:29, 31:60),4], pHCanVegMedEU[c(1:25, 27:29, 31:60),3])
