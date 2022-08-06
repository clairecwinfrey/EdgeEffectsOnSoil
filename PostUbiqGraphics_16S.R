# 16s- PostUbiquity Graphics 
# started August 2, 2022

# This script re-does ordinations and Bray-Curtis analyses with the postUbiquity analysis for prokaryote samples

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

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: EdgeEffectsAllSitesFUNGI.R

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
