# OrdinationsMedian
# December 5, 2021

# This script creates ordinations (NMDS) to show how samples vary based on Experimental unit (EU)
# and whether they are in the patch, forest, or on the edge. Here, I use the median ASV
# abundances for each meter along the transect in each EU (see description of medianEU.ps below)

# Libraries
library("phyloseq")
library("ggplot2")      
library("dplyr")        
library("tibble")       
library("tidyr")
library("vegan")
library("gridExtra")    #allows you to make multiple plots on the same page with ggplot
library("grDevices") #has function chull for finding convex hulls

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

# Load data
load("medianEU.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45, and 4) finally getting the 
# median ASV abundance at each meter in each EU (median is four values from each
# transect)

################
# ORDINATIONS WITHIN EU

## Separate out data by each EU, remove non-occurring ASVs, make ordinations, then plot
EU_52_med.ps <- subset_samples(medianEU.ps, EU == "EU_52")
EU_52_med.ps <- prune_taxa(taxa_sums(EU_52_med.ps) > 0, EU_52_med.ps) #remove non-occurring ASVs, now #4142 taxa 
medOrd52 <- ordinate(EU_52_med.ps, method = "NMDS", distance = "bray", trymax = 100)
sample_data(EU_52_med.ps)$meter <- factor(sample_data(EU_52_med.ps)$meter) #make meter ordinal for plotting

# here meter is colored....
med52_BrayNMDS <- phyloseq::plot_ordination(EU_52_med.ps, medOrd52, type= "samples", color= "meter") +
  ggtitle("EU 52")
quartz()
med52_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_52") + 
  geom_label(label = sample_data(EU_52_med.ps)$meter)

# Make convex hulls around the outermost points
patchOrdPoints <- medOrd52$points[c(1,3:5),] #pulls out points for 10-40m
forestOrdPoints <- medOrd52$points[c(2,6:10),] #pulls out points for 60-100m
chullPatch <- grDevices::chull(patchOrdPoints[,1], patchOrdPoints[,2])

# Make polygons around the points
med52_BrayNMDS <- phyloseq::plot_ordination(EU_52_med.ps, medOrd52, type= "samples", color= "Habitat")
+ ggtitle("EU 52")
#quartz()
med52_BrayNMDS <- med52_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_52") + 
  geom_label(label = sample_data(EU_52_med.ps)$meter)

##########
# EU 53N
EU_53N_med.ps <- subset_samples(medianEU.ps, EU == "EU_53N")
EU_53N_med.ps <- prune_taxa(taxa_sums(EU_53N_med.ps) > 0, EU_53N_med.ps) #remove non-occurring ASVs
medOrd53N <- ordinate(EU_53N_med.ps, method = "NMDS", distance = "bray", trymax = 100)
sample_data(EU_53N_med.ps)$meter <- factor(sample_data(EU_53N_med.ps)$meter) #make meter ordinal for plotting
### WARNING MESSAGE: Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
  #stress is (nearly) zero: you may have insufficient data

med53N_BrayNMDS <- phyloseq::plot_ordination(EU_53N_med.ps, medOrd53N, type= "samples", color= "Habitat")
+ ggtitle("EU 53N")
#quartz()
med53N_BrayNMDS <- med53N_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_53N") + 
  geom_label(label = sample_data(EU_53N_med.ps)$meter)

##########
# EU 54S
EU_54S_med.ps <- subset_samples(medianEU.ps, EU == "EU_54S")
EU_54S_med.ps <- prune_taxa(taxa_sums(EU_54S_med.ps) > 0, EU_54S_med.ps) #remove non-occurring ASVs, now 4141 taxa
medOrd54S <- ordinate(EU_54S_med.ps, method = "NMDS", distance = "bray", trymax = 100)
sample_data(EU_54S_med.ps)$meter <- factor(sample_data(EU_54S_med.ps)$meter) #make meter ordinal for plotting
med54S_BrayNMDS <- phyloseq::plot_ordination(EU_54S_med.ps, medOrd54S, type= "samples", color= "Habitat") +
  ggtitle("EU 54S")
med54S_BrayNMDS <- med54S_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_54S") + 
  geom_label(label = sample_data(EU_54S_med.ps)$meter)
# quartz()

##########
# EU 8
EU_8_med.ps <- subset_samples(medianEU.ps, EU == "EU_8")
EU_8_med.ps <- prune_taxa(taxa_sums(EU_8_med.ps) > 0, EU_8_med.ps) #remove non-occurring ASVs
medOrd8 <- ordinate(EU_8_med.ps, method = "NMDS", distance = "bray", trymax = 100)
sample_data(EU_8_med.ps)$meter <- factor(sample_data(EU_8_med.ps)$meter) #make meter ordinal for plotting
med8_BrayNMDS <- phyloseq::plot_ordination(EU_8_med.ps, medOrd8, type= "samples", color= "Habitat") +
  ggtitle("EU 8")
#quartz()
med8_BrayNMDS <- med8_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_8") + 
  geom_label(label = sample_data(EU_8_med.ps)$meter)
### WARNING MESSAGE: Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
#stress is (nearly) zero: you may have insufficient data

##########
# EU 53S
EU_53S_med.ps <- subset_samples(medianEU.ps, EU == "EU_53S")
EU_53S_med.ps <- prune_taxa(taxa_sums(EU_53S_med.ps) > 0, EU_53S_med.ps) #remove non-occurring ASVs
medOrd53S <- ordinate(EU_53S_med.ps, method = "NMDS", distance = "bray", trymax = 100)
sample_data(EU_53S_med.ps)$meter <- factor(sample_data(EU_53S_med.ps)$meter) #make meter ordinal for plotting
med53S_BrayNMDS <- phyloseq::plot_ordination(EU_53S_med.ps, medOrd53S, type= "samples", color= "Habitat") +
  ggtitle("EU 53S")
#quartz()
med53S_BrayNMDS <- med53S_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_53S") + 
  geom_label(label = sample_data(EU_53S_med.ps)$meter)

##########
# EU 10
EU_10_med.ps <- subset_samples(medianEU.ps, EU == "EU_10")
EU_10_med.ps <- prune_taxa(taxa_sums(EU_10_med.ps) > 0, EU_10_med.ps) #remove non-occurring ASVs
medOrd10 <- ordinate(EU_10_med.ps, method = "NMDS", distance = "bray", trymax = 100)
sample_data(EU_10_med.ps)$meter <- factor(sample_data(EU_10_med.ps)$meter) #make meter ordinal for plotting
med10_BrayNMDS <- phyloseq::plot_ordination(EU_10_med.ps, medOrd10, type= "samples", color= "Habitat") +
  ggtitle("EU 10")
#quartz()
med10_BrayNMDS <- med10_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_10") + 
  geom_label(label = sample_data(EU_10_med.ps)$meter)
med10_BrayNMDS

### WARNING MESSAGE: Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
#stress is (nearly) zero: you may have insufficient data

##### PLOT ALL TOGETHER #####
quartz()
require(gridExtra)
grid.arrange(med52_BrayNMDS, med53N_BrayNMDS, med54S_BrayNMDS, med8_BrayNMDS, med53S_BrayNMDS, med10_BrayNMDS, ncol=3)

#################
# ALL DATA (ALL EUS) TOGETHER
# based on patch, edge, and forest
set.seed(19)
medOrdAll <- ordinate(medianEU.ps, method = "NMDS", distance = "bray", trymax = 100) 

# Differences among habitats (based on patch, edge, and forest)
medBrayNMDS <- phyloseq::plot_ordination(medianEU.ps, medOrdAll, type= "samples", color= "Habitat")
quartz()
medBrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Differences among EUs
medBrayNMDS_byEU <- phyloseq::plot_ordination(medianEU.ps, medOrdAll, type= "samples", color= "EU")
quartz()
medBrayNMDS_byEU + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")



