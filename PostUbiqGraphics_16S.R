# 16s- PostUbiquity Graphics 
# started August 2, 2022

# This script re-does ordinations and Bray-Curtis analyses with the postUbiquity analysis for prokaryote samples. In addition,
# it performs a PERMANOVA to test if samples differ based on habitat type.

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

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: EdgeEffectsAllSitesFUNGI.R

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

#################################
# ORDINATIONS
#################################
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ord_16S <- ordinate(postUbiquity.ps, method = "NMDS", distance = "bray", trymax = 100) 

# Ordination based on "habitat"/ plant community type
HabitatNMDS_16S_postUbiq <- phyloseq::plot_ordination(postUbiquity.ps, ord_16S, type= "samples", color= "Habitat")
HabitatNMDS_16S_postUbiq <- HabitatNMDS_16S_postUbiq +
  scale_color_manual(values=c("purple", "darkgreen", "goldenrod")) + #change color of points
  geom_point(size=3) + ggtitle("(Prokaryote) NMDS based on Bray-Curtis Dissimilarities") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 7))) + #change size of legend and title
#quartz()
HabitatNMDS_16S_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# Ordination based on EU (just to show that they are different!)
EU_NMDS_16S_postUbiq <- phyloseq::plot_ordination(postUbiquity.ps, ord_16S, type= "samples", color= "EU")
EU_NMDS_16S_postUbiq <- EU_NMDS_16S_postUbiq + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
EU_NMDS_16S_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# load ITS habitat ordination plot (as made in PostUbiqGraphics_ITS.R)
load(file="RobjectsSaved/diffAbund_ITS_stackedBarplotPhyla_plot")
#quartz()
grid.arrange(HabitatNMDS_16S_postUbiq, HabitatNMDS_ITS_postUbiq, ncol=2)

#########################################################################
# PERMANOVA TO TEST FOR THE EFFECT OF HABITAT TYPE ON THE DISTRIBUTIONS
#########################################################################

# First, need to get data out of phyloseq so that it can be processed in vegan
postUbiq_16SASVs <- ASVs_outta_ps(postUbiquity.ps) #get ASV table
str(postUbiq_16SASVs) #samples are rows, taxa are columns as they should be for vegan
# Getting rid of the phyloseq attribute for metadata is near impossible, so I have to hack it a bit
# by grabbing each column that I need (as vectors or strings), then combining these
sampNames <- rownames(sample_data(postUbiquity.ps))
rownames(postUbiq_16SASVs) == sampNames #since these match, we are good to go!
habitat <- sample_data(postUbiquity.ps)$Habitat
sampID <- sample_data(postUbiquity.ps)$Sample.ID
metaDatForVeg <- cbind(sampID, habitat)
rownames(metaDatForVeg) <- sampNames
metaDatForVeg <- as.data.frame(metaDatForVeg)

# Next get Bray-Curtis dissimilarities
postUbiq_16SBC <- vegdist(postUbiq_16SASVs, method= "bray")
set.seed(19) #set seed so that results are reproducible!
# Run the PERMANOVA
postUbiq_16SPermanova <- adonis(postUbiq_16SBC ~ habitat, data=metaDatForVeg, permutations= 99999)
# p < 0.001

