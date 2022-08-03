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
load("RobjectsSaved/RobjectsSaved/ITS_postUbiqASVsPlus1_ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: EdgeEffectsAllSitesFUNGI.R