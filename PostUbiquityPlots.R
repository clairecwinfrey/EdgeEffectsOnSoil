# Post Ubiquity Plots (pre-median samples)
# December 5, 2021

# This script creates ordinations (NMDS) to show how samples vary based on Experimental unit (EU)
# and whether they are in the patch, forest, or on the edge. Here, I use the post ubiquity samples,
# but NOT the averaged samples based on taking the median at each transect meter for each EU. (See
# description of postUbiquity.ps below). In other words, this still has all the samples, 233, that 
# remained after rarefying).

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
load("postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45. PostUbiquity.ps was made and 
# saved in the script titled "UbiquityMedianSetup.R".

################
# ORDINATIONS WITHIN EU

## Separate out data by each EU, remove non-occurring ASVs, make ordinations, then plot
EU_52_ubiq.ps <- subset_samples(postUbiquity.ps, EU == "EU_52")
EU_52_ubiq.ps <- prune_taxa(taxa_sums(EU_52_ubiq.ps) > 0, EU_52_ubiq.ps) #remove non-occurring ASVs, now #4142 taxa 
ubiqOrd52 <- ordinate(EU_52_ubiq.ps, method = "NMDS", distance = "bray", trymax = 100)
#sample_data(EU_52_ubiq.ps)$meter <- factor(sample_data(EU_52_ubiq.ps)$meter) #make meter ordinal for plotting

# here meter is colored....
#ubiq52_BrayNMDS <- phyloseq::plot_ordination(EU_52_ubiq.ps, ubiqOrd52, type= "samples", color= "meter") +
#  ggtitle("EU 52")
#quartz()
#ubiq52_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_52") + 
#  geom_label(label = sample_data(EU_52_ubiq.ps)$meter)

# Make convex hulls around the outermost points
#patchOrdPoints <- ubiqOrd52$points[c(1,3:5),] #pulls out points for 10-40m
#forestOrdPoints <- ubiqOrd52$points[c(2,6:10),] #pulls out points for 60-100m
#chullPatch <- grDevices::chull(patchOrdPoints[,1], patchOrdPoints[,2])

# Make polygons around the points
ubiq52_BrayNMDS <- phyloseq::plot_ordination(EU_52_ubiq.ps, ubiqOrd52, type= "samples", color= "Habitat") +
  ggtitle("EU 52")
#quartz()
ubiq52_BrayNMDS <- ubiq52_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_52") + 
  geom_label(label = sample_data(EU_52_ubiq.ps)$Meter)
ubiq52_BrayNMDS

##########
# EU 53N
EU_53N_ubiq.ps <- subset_samples(postUbiquity.ps, EU == "EU_53N")
EU_53N_ubiq.ps <- prune_taxa(taxa_sums(EU_53N_ubiq.ps) > 0, EU_53N_ubiq.ps) #remove non-occurring ASVs
ubiqOrd53N <- ordinate(EU_53N_ubiq.ps, method = "NMDS", distance = "bray", trymax = 100)
#sample_data(EU_53N_ubiq.ps)$meter <- factor(sample_data(EU_53N_ubiq.ps)$meter) #make meter ordinal for plotting

ubiq53N_BrayNMDS <- phyloseq::plot_ordination(EU_53N_ubiq.ps, ubiqOrd53N, type= "samples", color= "Habitat") +
  ggtitle("EU 53N")

ubiq53N_BrayNMDS <- ubiq53N_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_53N") + 
  geom_label(label = sample_data(EU_53N_ubiq.ps)$Meter)
#quartz()
ubiq53N_BrayNMDS

##########
# EU 54S
EU_54S_ubiq.ps <- subset_samples(postUbiquity.ps, EU == "EU_54S")
EU_54S_ubiq.ps <- prune_taxa(taxa_sums(EU_54S_ubiq.ps) > 0, EU_54S_ubiq.ps) #remove non-occurring ASVs, now 4141 taxa
ubiqOrd54S <- ordinate(EU_54S_ubiq.ps, method = "NMDS", distance = "bray", trymax = 100)
#sample_data(EU_54S_ubiq.ps)$meter <- factor(sample_data(EU_54S_ubiq.ps)$meter) #make meter ordinal for plotting
ubiq54S_BrayNMDS <- phyloseq::plot_ordination(EU_54S_ubiq.ps, ubiqOrd54S, type= "samples", color= "Habitat") +
  ggtitle("EU 54S")
ubiq54S_BrayNMDS <- ubiq54S_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_54S") + 
  geom_label(label = sample_data(EU_54S_ubiq.ps)$Meter)
# quartz()
ubiq54S_BrayNMDS

##########
# EU 8
EU_8_ubiq.ps <- subset_samples(postUbiquity.ps, EU == "EU_8")
EU_8_ubiq.ps <- prune_taxa(taxa_sums(EU_8_ubiq.ps) > 0, EU_8_ubiq.ps) #remove non-occurring ASVs
ubiqOrd8 <- ordinate(EU_8_ubiq.ps, method = "NMDS", distance = "bray", trymax = 100)
#sample_data(EU_8_ubiq.ps)$meter <- factor(sample_data(EU_8_ubiq.ps)$meter) #make meter ordinal for plotting
ubiq8_BrayNMDS <- phyloseq::plot_ordination(EU_8_ubiq.ps, ubiqOrd8, type= "samples", color= "Habitat") +
  ggtitle("EU 8")
ubiq8_BrayNMDS <- ubiq8_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_8") + 
  geom_label(label = sample_data(EU_8_ubiq.ps)$Meter)
#quartz()
ubiq8_BrayNMDS

##########
# EU 53S
EU_53S_ubiq.ps <- subset_samples(postUbiquity.ps, EU == "EU_53S")
EU_53S_ubiq.ps <- prune_taxa(taxa_sums(EU_53S_ubiq.ps) > 0, EU_53S_ubiq.ps) #remove non-occurring ASVs
ubiqOrd53S <- ordinate(EU_53S_ubiq.ps, method = "NMDS", distance = "bray", trymax = 100)
### WARNING MESSAGE: Warning message:
# In metaMDS(veganifyOTU(physeq), distance, ...) :
#stress is (nearly) zero: you may have insufficient data
#sample_data(EU_53S_ubiq.ps)$meter <- factor(sample_data(EU_53S_ubiq.ps)$meter) #make meter ordinal for plotting
ubiq53S_BrayNMDS <- phyloseq::plot_ordination(EU_53S_ubiq.ps, ubiqOrd53S, type= "samples", color= "Habitat") +
  ggtitle("EU 53S")
ubiq53S_BrayNMDS <- ubiq53S_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_53S") + 
  geom_label(label = sample_data(EU_53S_ubiq.ps)$Meter)
#quartz()
ubiq53S_BrayNMDS 

##########
# EU 10
EU_10_ubiq.ps <- subset_samples(postUbiquity.ps, EU == "EU_10")
EU_10_ubiq.ps <- prune_taxa(taxa_sums(EU_10_ubiq.ps) > 0, EU_10_ubiq.ps) #remove non-occurring ASVs
ubiqOrd10 <- ordinate(EU_10_ubiq.ps, method = "NMDS", distance = "bray", trymax = 100)
#sample_data(EU_10_ubiq.ps)$meter <- factor(sample_data(EU_10_ubiq.ps)$meter) #make meter ordinal for plotting
ubiq10_BrayNMDS <- phyloseq::plot_ordination(EU_10_ubiq.ps, ubiqOrd10, type= "samples", color= "Habitat") +
  ggtitle("EU 10")
ubiq10_BrayNMDS <- ubiq10_BrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("EU_10") + 
  geom_label(label = sample_data(EU_10_ubiq.ps)$Meter)
#quartz()
ubiq10_BrayNMDS

##### PLOT ALL TOGETHER #####
quartz()
require(gridExtra)
grid.arrange(ubiq52_BrayNMDS, ubiq53N_BrayNMDS, ubiq54S_BrayNMDS, ubiq8_BrayNMDS, ubiq53S_BrayNMDS, ubiq10_BrayNMDS, ncol=3)

#################
# ALL DATA (ALL EUS) TOGETHER
ubiqOrdAll <- ordinate(postUbiquity.ps, method = "NMDS", distance = "bray", trymax = 1000) 

# Differences among habitats (based on patch, edge, and forest)
ubiqBrayNMDS <- phyloseq::plot_ordination(postUbiquity.ps, ubiqOrdAll, type= "samples", color= "Habitat")
quartz()
ubiqBrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Differences among EUs
ubiqBrayNMDS_byEU <- phyloseq::plot_ordination(postUbiquity.ps, ubiqOrdAll, type= "samples", color= "EU")
quartz()
ubiqBrayNMDS_byEU + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
