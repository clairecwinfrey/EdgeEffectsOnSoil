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
library("utilities")

# Load data
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
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
ord_ITS <- ordinate(ITS_postUbiquity.ps, method = "NMDS", distance = "bray", trymax = 100) 

# Ordination based on "habitat"/ plant community type (added to ESA presentation Aug 8)
HabitatNMDS_ITS_postUbiq <- phyloseq::plot_ordination(ITS_postUbiquity.ps, ord_ITS, type= "samples", color= "Habitat")
HabitatNMDS_ITS_postUbiq <- HabitatNMDS_ITS_postUbiq +
   scale_color_manual(values=c("purple", "darkgreen", "goldenrod")) + #change color of points
   geom_point(size=3) + ggtitle("Fungi") +
  theme_bw() +
  theme(plot.title = element_text(size=22)) +
  guides(color = guide_legend(override.aes = list(size = 7))) + #change size of legend and title
  theme(axis.ticks = element_blank(), #remove x and y axis labels and tick marks
        axis.text = element_blank())

quartz()
HabitatNMDS_ITS_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# Save the plot made above (saved Aug. 6, 2022) so that we can make a two paneled figure with ITS stuff
#save(HabitatNMDS_ITS_postUbiq, file="RobjectsSaved/diffAbund_ITS_stackedBarplotPhyla_plot")

### MAKEA TWO PANELED PLOT WITH THIS AND FUNGAL ORDINATION PLOT SIDE BY SIDE #####
# LOAD 16S plot made in PostUbiqGraphics_16S.R
load(file="RobjectsSaved/HabitatNMDS_16S_postUbiq")
#quartz()
grid.arrange(HabitatNMDS_16S_postUbiq, HabitatNMDS_ITS_postUbiq, ncol=2)
# Note: reflected the fungal NMDS plot in PowerPoint to make patch and forest on the same side! 

# Ordination based on EU (just to show that they are different!)
EU_NMDS_ITS_postUbiq <- phyloseq::plot_ordination(ITS_postUbiquity.ps, ord_ITS, type= "samples", color= "EU")
EU_NMDS_ITS_postUbiq <- EU_NMDS_ITS_postUbiq + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
EU_NMDS_ITS_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!


#########################################################################
# PERMANOVA TO TEST FOR THE EFFECT OF HABITAT TYPE ON THE DISTRIBUTIONS
#########################################################################

# First, need to get data out of phyloseq so that it can be processed in vegan
postUbiq_ITSASVs <- ASVs_outta_ps(ITS_postUbiquity.ps) #get ASV table
str(postUbiq_ITSASVs) #samples are rows, taxa are columns as they should be for vegan
# Getting rid of the phyloseq attribute for metadata is near impossible, so I have to hack it a bit
# by grabbing each column that I need (as vectors or strings), then combining these
sampNames <- rownames(sample_data(ITS_postUbiquity.ps))
rownames(postUbiq_ITSASVs) == sampNames #since these match, we are good to go!
habitat <- sample_data(ITS_postUbiquity.ps)$Habitat
sampID <- sample_data(ITS_postUbiquity.ps)$Sample.ID
metaDatForVeg <- cbind(sampID, habitat)
rownames(metaDatForVeg) <- sampNames
metaDatForVeg <- as.data.frame(metaDatForVeg)

# Next get Bray-Curtis dissimilarities
postUbiq_ITSBC <- vegdist(postUbiq_ITSASVs, method= "bray")
set.seed(19) #set seed so that results are reproducible!
# Run the PERMANOVA
postUbiq_ITSPermanova <- adonis(postUbiq_ITSBC ~ habitat, data=metaDatForVeg, permutations= 99999)
# p < 0.001

# Save PERMANOVA results
# save(postUbiq_ITSPermanova, file= "RobjectsSaved/postUbiq_ITSPermanova") #last saved Aug 9, 2022
