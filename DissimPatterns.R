# Exploration of Dissimilarity Patterns 
# October 7, 2021

# This script is to look at pattens of dissimilarity in our data. Specifically,
# Here I look at dissimilarities between: 1) 100 m (i.e. "most forest-y" and
# all other points, 2) 10 m (i.e. "most patch-y") and all other points, 3) the edge
# on the transect and all others along the transect in both directions, 4) the pooled
# forest samples versus all other points along the transect, and 5) the pooled patch
# samples versus all other points along the transect. I will likely use these to
# create GLMMs in the future.

# Importantly, analysis 1 and 2 have already mainly been done in, but here I break
# them up by EU

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

load("trimmedJustsoils.ps") #load phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset (from 
# 16SExploratoryDataAnalysisAug2021). 