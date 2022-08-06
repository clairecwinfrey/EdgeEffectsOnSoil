# ITS- PostUbiquity Graphics 
# started August 2, 2022

# This script re-does ordinations and Bray-Curtis analyses with the postUbiquity analysis for fungal samples

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
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: EdgeEffectsAllSitesFUNGI.R

#################################
# ORDINATIONS
#################################
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ord_ITS <- ordinate(ITS_postUbiquity.ps, method = "NMDS", distance = "bray", trymax = 100) 

# Ordination based on "habitat"/ plant community type (added to ESA presentation Aug 8)
HabitatNMDS_ITS_postUbiq <- phyloseq::plot_ordination(ITS_postUbiquity.ps, ord_ITS, type= "samples", color= "Habitat")
HabitatNMDS_ITS_postUbiq <- HabitatNMDS_ITS_postUbiq +
   scale_color_manual(values=c("purple", "darkgreen", "goldenrod")) + #change color of points
   geom_point(size=3) + ggtitle("(Fungal) NMDS based on Bray-Curtis Dissimilarities") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 7))) + #change size of legend and title
quartz()
HabitatNMDS_ITS_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# Save the plot made above (saved Aug. 6, 2022) so that we can make a two paneled figure with 16S stuff
save(HabitatNMDS_ITS_postUbiq, file="RobjectsSaved/diffAbund_ITS_stackedBarplotPhyla_plot")

# Ordination based on EU (just to show that they are different!)
EU_NMDS_ITS_postUbiq <- phyloseq::plot_ordination(ITS_postUbiquity.ps, ord_ITS, type= "samples", color= "EU")
EU_NMDS_ITS_postUbiq <- EU_NMDS_ITS_postUbiq + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
EU_NMDS_ITS_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

quartz()
grid.arrange(EU_NMDS_ITS_postUbiq, HabitatNMDS_ITS_postUbiq, ncol=2)
