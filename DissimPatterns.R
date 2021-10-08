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

# Libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

load("trimmedJustsoils.ps") #load phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset (from 
# 16SExploratoryDataAnalysisAug2021).

# Functions defined in this script:
# Function that gets OTU/ASV table out of phyloseq
# (from https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

#########################################

# Separate out samples by EU
EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU)
EU_53N_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53N")
unique(sample_data(EU_53N_Soils.ps)$EU)
EU_54S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_54S")
unique(sample_data(EU_54S_Soils.ps)$EU)
EU_8_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_8")
unique(sample_data(EU_8_Soils.ps)$EU)
EU_53S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53S")
unique(sample_data(EU_53S_Soils.ps)$EU)
EU_10_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_10")
unique(sample_data(EU_10_Soils.ps)$EU)


#########################################################################
# 1. and 2. Dissimilarity along transect from 10 m and 100 m 
# (this is only modified slightly from code in 16SEULevelEDASept2021.R)
#########################################################################

###########
# EU 52
##########
# Get ASV table
ASVtab_52 <- psotu2veg(EU_52_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_52 <- vegdist(ASVtab_52, method="bray")
BrayDist_52.mat <- as.matrix(BrayDist_52)
diag(BrayDist_52.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_52.mat) == rownames(sample_data(EU_52_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps52 <- colnames(BrayDist_52.mat) #this is the Mapping Sample ID of these
meter52 <- sample_data(EU_52_Soils.ps)$Meter
index52.df <- data.frame(samps52, meter52)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_52 <- which(index52.df$meter52 == 10)
m20_52 <- which(index52.df$meter52 == 20)
m30_52 <- which(index52.df$meter52 == 30)
m40_52 <- which(index52.df$meter52 == 40)
m50_52 <- which(index52.df$meter52 == 50)
m60_52 <- which(index52.df$meter52 == 60)
m70_52 <- which(index52.df$meter52 == 70)
m80_52 <- which(index52.df$meter52 == 80)
m90_52 <- which(index52.df$meter52 == 90)
m100_52 <- which(index52.df$meter52 == 100)

# B-C distances between patch at 10m and each point along transect
m10_52_comp10 <- BrayDist_52.mat[m10_52, m10_52] #add 10 and 10
m10_52_comp20 <- BrayDist_52.mat[m10_52, m20_52]
m10_52_comp30 <- BrayDist_52.mat[m10_52, m30_52]
m10_52_comp40 <- BrayDist_52.mat[m10_52, m40_52]
m10_52_comp50 <- BrayDist_52.mat[m10_52, m50_52]
m10_52_comp60 <- BrayDist_52.mat[m10_52, m60_52]
m10_52_comp70 <- BrayDist_52.mat[m10_52, m70_52]
m10_52_comp80 <- BrayDist_52.mat[m10_52, m80_52]
m10_52_comp90 <- BrayDist_52.mat[m10_52, m90_52]
m10_52_comp100 <- BrayDist_52.mat[m10_52, m100_52]


# B-C distances between forest at 100m and each point along transect
m100_52_comp100 <- BrayDist_52.mat[m100_52, m100_52] #add 100 and 100
m100_52_comp90 <- BrayDist_52.mat[m100_52, m90_52]
m100_52_comp80 <- BrayDist_52.mat[m100_52, m80_52]
m100_52_comp70 <- BrayDist_52.mat[m100_52, m70_52]
m100_52_comp60 <- BrayDist_52.mat[m100_52, m60_52]
m100_52_comp50 <- BrayDist_52.mat[m100_52, m50_52]
m100_52_comp40 <- BrayDist_52.mat[m100_52, m40_52]
m100_52_comp30 <- BrayDist_52.mat[m100_52, m30_52]
m100_52_comp20 <- BrayDist_52.mat[m100_52, m20_52]
m100_52_comp10 <- BrayDist_52.mat[m100_52, m10_52]

bold_eu52p <- expression(bold("EU 52: Dissimilarity from 10m (patch)"))
bold_eu52f <- expression(bold("EU 52: Dissimilarity from 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box52 <- boxplot(list(m10_52_comp10, m10_52_comp20, m10_52_comp30, m10_52_comp40,
                            m10_52_comp50, m10_52_comp60, m10_52_comp70,
                            m10_52_comp80, m10_52_comp90, m10_52_comp100),
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                       names = c("10", "20", "30", "40", "50",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu52p, side=3, adj = -0.065, line = 2)
forest_box52 <- boxplot(list(m100_52_comp10, m100_52_comp20, m100_52_comp30, m100_52_comp40,
                             m100_52_comp50, m100_52_comp60, m100_52_comp70,
                             m100_52_comp80, m100_52_comp90, m100_52_comp100), 
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu52f, side=3, adj = -0.065, line = 2)


###########
# EU 53N
##########

ASVtab_53N <- psotu2veg(EU_53N_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_53N <- vegdist(ASVtab_53N, method="bray")
BrayDist_53N.mat <- as.matrix(BrayDist_53N)
diag(BrayDist_53N.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53N.mat) == rownames(sample_data(EU_53N_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53N <- colnames(BrayDist_53N.mat) #this is the Mapping Sample ID of these
meter53N <- sample_data(EU_53N_Soils.ps)$Meter
index53N.df <- data.frame(samps53N, meter53N)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53N <- which(index53N.df$meter53N == 10)
m20_53N <- which(index53N.df$meter53N == 20)
m30_53N <- which(index53N.df$meter53N == 30)
m40_53N <- which(index53N.df$meter53N == 40)
m50_53N <- which(index53N.df$meter53N == 50)
m60_53N <- which(index53N.df$meter53N == 60)
m70_53N <- which(index53N.df$meter53N == 70)
m80_53N <- which(index53N.df$meter53N == 80)
m90_53N <- which(index53N.df$meter53N == 90)
m100_53N <- which(index53N.df$meter53N == 100)

# B-C distances between patch at 10m and each point along transect
m10_53N_comp10 <- BrayDist_53N.mat[m10_53N, m10_53N] #add 10 and 10
m10_53N_comp20 <- BrayDist_53N.mat[m10_53N, m20_53N]
m10_53N_comp30 <- BrayDist_53N.mat[m10_53N, m30_53N]
m10_53N_comp40 <- BrayDist_53N.mat[m10_53N, m40_53N]
m10_53N_comp50 <- BrayDist_53N.mat[m10_53N, m50_53N]
m10_53N_comp60 <- BrayDist_53N.mat[m10_53N, m60_53N]
m10_53N_comp70 <- BrayDist_53N.mat[m10_53N, m70_53N]
m10_53N_comp80 <- BrayDist_53N.mat[m10_53N, m80_53N]
m10_53N_comp90 <- BrayDist_53N.mat[m10_53N, m90_53N]
m10_53N_comp100 <- BrayDist_53N.mat[m10_53N, m100_53N]

# B-C distances between forest at 100m and each point along transect
m100_53N_comp100 <- BrayDist_53N.mat[m100_53N, m100_53N]
m100_53N_comp90 <- BrayDist_53N.mat[m100_53N, m90_53N]
m100_53N_comp80 <- BrayDist_53N.mat[m100_53N, m80_53N]
m100_53N_comp70 <- BrayDist_53N.mat[m100_53N, m70_53N]
m100_53N_comp60 <- BrayDist_53N.mat[m100_53N, m60_53N]
m100_53N_comp50 <- BrayDist_53N.mat[m100_53N, m50_53N]
m100_53N_comp40 <- BrayDist_53N.mat[m100_53N, m40_53N]
m100_53N_comp30 <- BrayDist_53N.mat[m100_53N, m30_53N]
m100_53N_comp20 <- BrayDist_53N.mat[m100_53N, m20_53N]
m100_53N_comp10 <- BrayDist_53N.mat[m100_53N, m10_53N]


bold_eu53Np <- expression(bold("EU 53N: Dissimilarity from 10m (patch)"))
bold_eu53Nf <- expression(bold("EU 53N: Dissimilarity from 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box53N <- boxplot(list(m10_53N_comp10, m10_53N_comp20, m10_53N_comp30, m10_53N_comp40,
                             m10_53N_comp50, m10_53N_comp60, m10_53N_comp70,
                             m10_53N_comp80, m10_53N_comp90, m10_53N_comp100),
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10m","20m", "30m", "40m", "50m",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu53Np, side=3, adj = -0.065, line = 2)
forest_box53N <- boxplot(list(m100_53N_comp10, m100_53N_comp20, m100_53N_comp30, m100_53N_comp40,
                              m100_53N_comp50, m100_53N_comp60, m100_53N_comp70,
                              m100_53N_comp80, m100_53N_comp90, m100_53N_comp100), 
                         ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                         names = c("10", "20", "30", "40", "50",
                                   "60", "70", "80", "90", "100"), cex.axis = 0.8,
                         cex.lab = 1,
                         ylim=c(0.0, 1.0))
mtext(text=bold_eu53Nf, side=3, adj = -0.065, line = 2)


###########
# EU 54S
##########

ASVtab_54S <- psotu2veg(EU_54S_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_54S <- vegdist(ASVtab_54S, method="bray")
BrayDist_54S.mat <- as.matrix(BrayDist_54S)
diag(BrayDist_54S.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_54S.mat) == rownames(sample_data(EU_54S_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps54S <- colnames(BrayDist_54S.mat) #this is the Mapping Sample ID of these
meter54S <- sample_data(EU_54S_Soils.ps)$Meter
index54S.df <- data.frame(samps54S, meter54S)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_54S <- which(index54S.df$meter54S == 10)
m20_54S <- which(index54S.df$meter54S == 20)
m30_54S <- which(index54S.df$meter54S == 30)
m40_54S <- which(index54S.df$meter54S == 40)
m50_54S <- which(index54S.df$meter54S == 50)
m60_54S <- which(index54S.df$meter54S == 60)
m70_54S <- which(index54S.df$meter54S == 70)
m80_54S <- which(index54S.df$meter54S == 80)
m90_54S <- which(index54S.df$meter54S == 90)
m100_54S <- which(index54S.df$meter54S == 100)

# B-C distances between patch at 10m and each point along transect
m10_54S_comp10 <- BrayDist_54S.mat[m10_54S, m10_54S]
m10_54S_comp20 <- BrayDist_54S.mat[m10_54S, m20_54S]
m10_54S_comp30 <- BrayDist_54S.mat[m10_54S, m30_54S]
m10_54S_comp40 <- BrayDist_54S.mat[m10_54S, m40_54S]
m10_54S_comp50 <- BrayDist_54S.mat[m10_54S, m50_54S]
m10_54S_comp60 <- BrayDist_54S.mat[m10_54S, m60_54S]
m10_54S_comp70 <- BrayDist_54S.mat[m10_54S, m70_54S]
m10_54S_comp80 <- BrayDist_54S.mat[m10_54S, m80_54S]
m10_54S_comp90 <- BrayDist_54S.mat[m10_54S, m90_54S]
m10_54S_comp100 <- BrayDist_54S.mat[m10_54S, m100_54S]

# B-C distances between forest at 100m and each point along transect
m100_54S_comp100 <- BrayDist_54S.mat[m100_54S, m100_54S]
m100_54S_comp90 <- BrayDist_54S.mat[m100_54S, m90_54S]
m100_54S_comp80 <- BrayDist_54S.mat[m100_54S, m80_54S]
m100_54S_comp70 <- BrayDist_54S.mat[m100_54S, m70_54S]
m100_54S_comp60 <- BrayDist_54S.mat[m100_54S, m60_54S]
m100_54S_comp50 <- BrayDist_54S.mat[m100_54S, m50_54S]
m100_54S_comp40 <- BrayDist_54S.mat[m100_54S, m40_54S]
m100_54S_comp30 <- BrayDist_54S.mat[m100_54S, m30_54S]
m100_54S_comp20 <- BrayDist_54S.mat[m100_54S, m20_54S]
m100_54S_comp10 <- BrayDist_54S.mat[m100_54S, m10_54S]


# Plot
bold_eu54Sp <- expression(bold("EU 54S: Dissimilarity from 10m (patch)"))
bold_eu54Sf <- expression(bold("EU 54S: Dissimilarity from 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box54S <- boxplot(list(m10_54S_comp10, m10_54S_comp20, m10_54S_comp30, m10_54S_comp40,
                             m10_54S_comp50, m10_54S_comp60, m10_54S_comp70,
                             m10_54S_comp80, m10_54S_comp90, m10_54S_comp100),
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10m","20m", "30m", "40m", "50m",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu54Sp, side=3, adj = -0.065, line = 2)
forest_box54S <- boxplot(list(m100_54S_comp10, m100_54S_comp20, m100_54S_comp30, m100_54S_comp40,
                              m100_54S_comp50, m100_54S_comp60, m100_54S_comp70,
                              m100_54S_comp80, m100_54S_comp90, m100_54S_comp100), 
                         ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                         names = c("10", "20", "30", "40", "50",
                                   "60", "70", "80", "90", "100"), cex.axis = 0.8,
                         cex.lab = 1,
                         ylim=c(0.0, 1.0))
mtext(text=bold_eu54Sf, side=3, adj = -0.065, line = 2)

###########
# EU 8
##########

ASVtab_8 <- psotu2veg(EU_8_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_8 <- vegdist(ASVtab_8, method="bray")
BrayDist_8.mat <- as.matrix(BrayDist_8)
diag(BrayDist_8.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_8.mat) == rownames(sample_data(EU_8_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps8 <- colnames(BrayDist_8.mat) #this is the Mapping Sample ID of these
meter8 <- sample_data(EU_8_Soils.ps)$Meter
index8.df <- data.frame(samps8, meter8)
sample_data(EU_8_Soils.ps)$Sample.ID

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_8 <- which(index8.df$meter8 == 10)
m20_8 <- which(index8.df$meter8 == 20)
m30_8 <- which(index8.df$meter8 == 30)
m40_8 <- which(index8.df$meter8 == 40)
m50_8 <- which(index8.df$meter8 == 50)
m60_8 <- which(index8.df$meter8 == 60)
m70_8 <- which(index8.df$meter8 == 70)
m80_8 <- which(index8.df$meter8 == 80)
m90_8 <- which(index8.df$meter8 == 90)
m100_8 <- which(index8.df$meter8 == 100)

# B-C distances between patch at 10m and each point along transect
m10_8_comp10 <- BrayDist_8.mat[m10_8, m10_8]
m10_8_comp20 <- BrayDist_8.mat[m10_8, m20_8]
m10_8_comp30 <- BrayDist_8.mat[m10_8, m30_8]
m10_8_comp40 <- BrayDist_8.mat[m10_8, m40_8]
m10_8_comp50 <- BrayDist_8.mat[m10_8, m50_8]
m10_8_comp60 <- BrayDist_8.mat[m10_8, m60_8]
m10_8_comp70 <- BrayDist_8.mat[m10_8, m70_8]
m10_8_comp80 <- BrayDist_8.mat[m10_8, m80_8]
m10_8_comp90 <- BrayDist_8.mat[m10_8, m90_8]
m10_8_comp100 <- BrayDist_8.mat[m10_8, m100_8]


# B-C distances between forest at 100m and each point along transect
m100_8_comp100 <- BrayDist_8.mat[m100_8, m100_8]
m100_8_comp90 <- BrayDist_8.mat[m100_8, m90_8]
m100_8_comp80 <- BrayDist_8.mat[m100_8, m80_8]
m100_8_comp70 <- BrayDist_8.mat[m100_8, m70_8]
m100_8_comp60 <- BrayDist_8.mat[m100_8, m60_8]
m100_8_comp50 <- BrayDist_8.mat[m100_8, m50_8]
m100_8_comp40 <- BrayDist_8.mat[m100_8, m40_8]
m100_8_comp30 <- BrayDist_8.mat[m100_8, m30_8]
m100_8_comp20 <- BrayDist_8.mat[m100_8, m20_8]
m100_8_comp10 <- BrayDist_8.mat[m100_8, m10_8]

# Plot
bold_eu8p <- expression(bold("EU 8: Dissimilarity from 10m (patch)"))
bold_eu8f <- expression(bold("EU 8: Dissimilarity from 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box8 <- boxplot(list(m10_8_comp10, m10_8_comp20, m10_8_comp30, m10_8_comp40,
                           m10_8_comp50, m10_8_comp60, m10_8_comp70,
                           m10_8_comp80, m10_8_comp90, m10_8_comp100),
                      col= "burlywood3",
                      ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                      names = c("10", "20", "30", "40", "50",
                                "60", "70", "80", "90", "100"), cex.axis = 0.8,
                      cex.lab = 1,
                      ylim=c(0.0, 1.0))
mtext(text=bold_eu8p, side=3, adj = -0.065, line = 2)
forest_box8 <- boxplot(list(m100_8_comp10, m100_8_comp20, m100_8_comp30, m100_8_comp40,
                            m100_8_comp50, m100_8_comp60, m100_8_comp70,
                            m100_8_comp80, m100_8_comp90, m100_8_comp100),
                       col= "darkgreen",
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                       names = c("10", "20", "30", "40", "50",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu8f, side=3, adj = -0.065, line = 2)


###########
# EU 53S
##########

ASVtab_53S <- psotu2veg(EU_53S_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_53S <- vegdist(ASVtab_53S, method="bray")
BrayDist_53S.mat <- as.matrix(BrayDist_53S)
diag(BrayDist_53S.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53S.mat) == rownames(sample_data(EU_53S_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53S <- colnames(BrayDist_53S.mat) #this is the Mapping Sample ID of these
meter53S <- sample_data(EU_53S_Soils.ps)$Meter
index53S.df <- data.frame(samps53S, meter53S)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53S <- which(index53S.df$meter53S == 10)
m20_53S <- which(index53S.df$meter53S == 20)
m30_53S <- which(index53S.df$meter53S == 30)
m40_53S <- which(index53S.df$meter53S == 40)
m50_53S <- which(index53S.df$meter53S == 50)
m60_53S <- which(index53S.df$meter53S == 60)
m70_53S <- which(index53S.df$meter53S == 70)
m80_53S <- which(index53S.df$meter53S == 80)
m90_53S <- which(index53S.df$meter53S == 90)
m100_53S <- which(index53S.df$meter53S == 100)

# B-C distances between patch at 10m and each point along transect
m10_53S_comp10 <- BrayDist_53S.mat[m10_53S, m10_53S]
m10_53S_comp20 <- BrayDist_53S.mat[m10_53S, m20_53S]
m10_53S_comp30 <- BrayDist_53S.mat[m10_53S, m30_53S]
m10_53S_comp40 <- BrayDist_53S.mat[m10_53S, m40_53S]
m10_53S_comp50 <- BrayDist_53S.mat[m10_53S, m50_53S]
m10_53S_comp60 <- BrayDist_53S.mat[m10_53S, m60_53S]
m10_53S_comp70 <- BrayDist_53S.mat[m10_53S, m70_53S]
m10_53S_comp80 <- BrayDist_53S.mat[m10_53S, m80_53S]
m10_53S_comp90 <- BrayDist_53S.mat[m10_53S, m90_53S]
m10_53S_comp100 <- BrayDist_53S.mat[m10_53S, m100_53S]


# B-C distances between forest at 100m and each point along transect
m100_53S_comp100 <- BrayDist_53S.mat[m100_53S, m100_53S]
m100_53S_comp90 <- BrayDist_53S.mat[m100_53S, m90_53S]
m100_53S_comp80 <- BrayDist_53S.mat[m100_53S, m80_53S]
m100_53S_comp70 <- BrayDist_53S.mat[m100_53S, m70_53S]
m100_53S_comp60 <- BrayDist_53S.mat[m100_53S, m60_53S]
m100_53S_comp50 <- BrayDist_53S.mat[m100_53S, m50_53S]
m100_53S_comp40 <- BrayDist_53S.mat[m100_53S, m40_53S]
m100_53S_comp30 <- BrayDist_53S.mat[m100_53S, m30_53S]
m100_53S_comp20 <- BrayDist_53S.mat[m100_53S, m20_53S]
m100_53S_comp10 <- BrayDist_53S.mat[m100_53S, m10_53S]

# Plot

bold_eu53Sp <- expression(bold("EU 53S: Dissimilarity from 10m (patch)"))
bold_eu53Sf <- expression(bold("EU 53S: Dissimilarity from 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box53S <- boxplot(list(m10_53S_comp10, m10_53S_comp20, m10_53S_comp30, m10_53S_comp40,
                             m10_53S_comp50, m10_53S_comp60, m10_53S_comp70,
                             m10_53S_comp80, m10_53S_comp90, m10_53S_comp100),
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu53Sp, side=3, adj = -0.065, line = 2)
forest_box53S <- boxplot(list(m100_53S_comp10, m100_53S_comp20, m100_53S_comp30, m100_53S_comp40,
                              m100_53S_comp50, m100_53S_comp60, m100_53S_comp70,
                              m100_53S_comp80, m100_53S_comp90, m10_53S_comp100), 
                         ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                         names = c("10", "20", "30", "40", "50",
                                   "60", "70", "80", "90", "100"), cex.axis = 0.8,
                         cex.lab = 1,
                         ylim=c(0.0, 1.0))
mtext(text=bold_eu53Sf, side=3, adj = -0.065, line = 2)

###########
# EU 10
##########
ASVtab_10 <- psotu2veg(EU_10_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_10 <- vegdist(ASVtab_10, method="bray")
BrayDist_10.mat <- as.matrix(BrayDist_10)
diag(BrayDist_10.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_10.mat) == rownames(sample_data(EU_10_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps10 <- colnames(BrayDist_10.mat) #this is the Mapping Sample ID of these
meter10 <- sample_data(EU_10_Soils.ps)$Meter
index10.df <- data.frame(samps10, meter10)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_10 <- which(index10.df$meter10 == 10)
m20_10 <- which(index10.df$meter10 == 20)
m30_10 <- which(index10.df$meter10 == 30)
m40_10 <- which(index10.df$meter10 == 40)
m50_10 <- which(index10.df$meter10 == 50)
m60_10 <- which(index10.df$meter10 == 60)
m70_10 <- which(index10.df$meter10 == 70)
m80_10 <- which(index10.df$meter10 == 80)
m90_10 <- which(index10.df$meter10 == 90)
m100_10 <- which(index10.df$meter10 == 100)

# B-C distances between patch at 10m and each point along transect
m10_10_comp10 <- BrayDist_10.mat[m10_10, m20_10]
m10_10_comp20 <- BrayDist_10.mat[m10_10, m20_10]
m10_10_comp30 <- BrayDist_10.mat[m10_10, m30_10]
m10_10_comp40 <- BrayDist_10.mat[m10_10, m40_10]
m10_10_comp50 <- BrayDist_10.mat[m10_10, m50_10]
m10_10_comp60 <- BrayDist_10.mat[m10_10, m60_10]
m10_10_comp70 <- BrayDist_10.mat[m10_10, m70_10]
m10_10_comp80 <- BrayDist_10.mat[m10_10, m80_10]
m10_10_comp90 <- BrayDist_10.mat[m10_10, m90_10]
m10_10_comp100 <- BrayDist_10.mat[m10_10, m100_10]

# B-C distances between forest at 100m and each point along transect
m100_10_comp100 <- BrayDist_10.mat[m100_10, m100_10]
m100_10_comp90 <- BrayDist_10.mat[m100_10, m90_10]
m100_10_comp80 <- BrayDist_10.mat[m100_10, m80_10]
m100_10_comp70 <- BrayDist_10.mat[m100_10, m70_10]
m100_10_comp60 <- BrayDist_10.mat[m100_10, m60_10]
m100_10_comp50 <- BrayDist_10.mat[m100_10, m50_10]
m100_10_comp40 <- BrayDist_10.mat[m100_10, m40_10]
m100_10_comp30 <- BrayDist_10.mat[m100_10, m30_10]
m100_10_comp20 <- BrayDist_10.mat[m100_10, m20_10]
m100_10_comp10 <- BrayDist_10.mat[m100_10, m10_10]

# Plot:
bold_eu10p <- expression(bold("EU 10: Dissimilarity from 10m (patch)"))
bold_eu10f <- expression(bold("EU 10: Dissimilarity from 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box10 <- boxplot(list(m10_10_comp10, m10_10_comp20, m10_10_comp30, m10_10_comp40,
                            m10_10_comp50, m10_10_comp60, m10_10_comp70,
                            m10_10_comp80, m10_10_comp90, m10_10_comp100),
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                       names = c("10", "20m", "30m", "40m", "50m",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu10p, side=3, adj = -0.065, line = 2)
forest_box10 <- boxplot(list(m100_10_comp10, m100_10_comp20, m100_10_comp30, m100_10_comp40,
                             m100_10_comp50, m100_10_comp60, m100_10_comp70,
                             m100_10_comp80, m100_10_comp90, m100_10_comp100), 
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu10f, side=3, adj = -0.065, line = 2)


######## MAKE BOXPLOTS ######## 
# Put everything together in one big ol' plot

quartz()
par(mfrow=c(3,4))
patch_box52 <- boxplot(list(m10_52_comp10, m10_52_comp20, m10_52_comp30, m10_52_comp40,
                            m10_52_comp50, m10_52_comp60, m10_52_comp70,
                            m10_52_comp80, m10_52_comp90, m10_52_comp100),
                       col= "burlywood3",
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                       names = c("10", "20", "30", "40", "50",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu52p, side=3, adj = -0.065, line = 2)
forest_box52 <- boxplot(list(m100_52_comp10, m100_52_comp20, m100_52_comp30, m100_52_comp40,
                             m100_52_comp50, m100_52_comp60, m100_52_comp70,
                             m100_52_comp80, m100_52_comp90, m100_52_comp100),
                        col= "darkgreen",
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu52f, side=3, adj = -0.065, line = 2)
patch_box53N <- boxplot(list(m10_53N_comp10, m10_53N_comp20, m10_53N_comp30, m10_53N_comp40,
                             m10_53N_comp50, m10_53N_comp60, m10_53N_comp70,
                             m10_53N_comp80, m10_53N_comp90, m10_53N_comp100),
                        col= "burlywood3",
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10m","20m", "30m", "40m", "50m",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu53Np, side=3, adj = -0.065, line = 2)
forest_box53N <- boxplot(list(m100_53N_comp10, m100_53N_comp20, m100_53N_comp30, m100_53N_comp40,
                              m100_53N_comp50, m100_53N_comp60, m100_53N_comp70,
                              m100_53N_comp80, m100_53N_comp90, m100_53N_comp100),
                         col= "darkgreen",
                         ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                         names = c("10", "20", "30", "40", "50",
                                   "60", "70", "80", "90", "100"), cex.axis = 0.8,
                         cex.lab = 1,
                         ylim=c(0.0, 1.0))
mtext(text=bold_eu53Nf, side=3, adj = -0.065, line = 2)
patch_box54S <- boxplot(list(m10_54S_comp10, m10_54S_comp20, m10_54S_comp30, m10_54S_comp40,
                             m10_54S_comp50, m10_54S_comp60, m10_54S_comp70,
                             m10_54S_comp80, m10_54S_comp90, m10_54S_comp100),
                        col= "burlywood3",
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10m","20m", "30m", "40m", "50m",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu54Sp, side=3, adj = -0.065, line = 2)
forest_box54S <- boxplot(list(m100_54S_comp10, m100_54S_comp20, m100_54S_comp30, m100_54S_comp40,
                              m100_54S_comp50, m100_54S_comp60, m100_54S_comp70,
                              m100_54S_comp80, m100_54S_comp90, m100_54S_comp100),
                         col= "darkgreen",
                         ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                         names = c("10", "20", "30", "40", "50",
                                   "60", "70", "80", "90", "100"), cex.axis = 0.8,
                         cex.lab = 1,
                         ylim=c(0.0, 1.0))
mtext(text=bold_eu54Sf, side=3, adj = -0.065, line = 2)
patch_box8 <- boxplot(list(m10_8_comp10, m10_8_comp20, m10_8_comp30, m10_8_comp40,
                           m10_8_comp50, m10_8_comp60, m10_8_comp70,
                           m10_8_comp80, m10_8_comp90, m10_8_comp100),
                      col= "burlywood3",
                      ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                      names = c("10", "20", "30", "40", "50",
                                "60", "70", "80", "90", "100"), cex.axis = 0.8,
                      cex.lab = 1,
                      ylim=c(0.0, 1.0))
mtext(text=bold_eu8p, side=3, adj = -0.065, line = 2)
forest_box8 <- boxplot(list(m100_8_comp10, m100_8_comp20, m100_8_comp30, m100_8_comp40,
                            m100_8_comp50, m100_8_comp60, m100_8_comp70,
                            m100_8_comp80, m100_8_comp90, m100_8_comp100),
                       col= "darkgreen",
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                       names = c("10", "20", "30", "40", "50",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu8f, side=3, adj = -0.065, line = 2)
patch_box53S <- boxplot(list(m10_53S_comp10, m10_53S_comp20, m10_53S_comp30, m10_53S_comp40,
                             m10_53S_comp50, m10_53S_comp60, m10_53S_comp70,
                             m10_53S_comp80, m10_53S_comp90, m10_53S_comp100),
                        col= "burlywood3",
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu53Sp, side=3, adj = -0.065, line = 2)
forest_box53S <- boxplot(list(m100_53S_comp10, m100_53S_comp20, m100_53S_comp30, m100_53S_comp40,
                              m100_53S_comp50, m100_53S_comp60, m100_53S_comp70,
                              m100_53S_comp80, m100_53S_comp90, m10_53S_comp100),
                         col= "darkgreen",
                         ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter", 
                         names = c("10", "20", "30", "40", "50",
                                   "60", "70", "80", "90", "100"), cex.axis = 0.8,
                         cex.lab = 1,
                         ylim=c(0.0, 1.0))
mtext(text=bold_eu53Sf, side=3, adj = -0.065, line = 2)
patch_box10 <- boxplot(list(m10_10_comp10, m10_10_comp20, m10_10_comp30, m10_10_comp40,
                            m10_10_comp50, m10_10_comp60, m10_10_comp70,
                            m10_10_comp80, m10_10_comp90, m10_10_comp100),
                       col= "burlywood3",
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter", 
                       names = c("10", "20m", "30m", "40m", "50m",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu10p, side=3, adj = -0.065, line = 2)
forest_box10 <- boxplot(list(m100_10_comp10, m100_10_comp20, m100_10_comp30, m100_10_comp40,
                             m100_10_comp50, m100_10_comp60, m100_10_comp70,
                             m100_10_comp80, m100_10_comp90, m100_10_comp100),
                        col= "darkgreen",
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), cex.axis = 0.8,
                        cex.lab = 1,
                        ylim=c(0.0, 1.0))
mtext(text=bold_eu10f, side=3, adj = -0.065, line = 2)

########################
# Does the Bray-Curtis distance b/w samples tend to increase with distance?
# this section copied from 16SExploratoryDataAnalysisAug2021)
########################
ASVtab <- psotu2veg(trimmedJustsoils.ps)
#View(ASVtab)

ASVbrayDist <- vegdist(ASVtab, method = "bray")
ASVBrayDist.mat <- as.matrix(ASVbrayDist)
diag(ASVBrayDist.mat) <- NA

# Make data frame for indexing matrix
rownames(ASVBrayDist.mat) == rownames(sample_data(trimmedJustsoils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps <- colnames(ASVBrayDist.mat)
meter <- sample_data(trimmedJustsoils.ps)$Meter
index.df <- data.frame(samps, meter)

# Get row and column numbers of different meters
m10 <- which(index.df$meter == 10)
m20 <- which(index.df$meter == 20)
m30<- which(index.df$meter == 30)
m40 <- which(index.df$meter == 40)
m50 <- which(index.df$meter == 50)
m60 <- which(index.df$meter == 60)
m70 <- which(index.df$meter == 70)
m80<- which(index.df$meter == 80)
m90 <- which(index.df$meter == 90)
m100 <- which(index.df$meter == 100)

# Distances between patch at 10m and each point along transect
m10_comp10 <- ASVBrayDist.mat[m10, m10]
m10_comp20 <- ASVBrayDist.mat[m10, m20]
m10_comp30 <- ASVBrayDist.mat[m10, m30]
m10_comp40 <- ASVBrayDist.mat[m10, m40]
m10_comp50 <- ASVBrayDist.mat[m10, m50]
m10_comp60 <- ASVBrayDist.mat[m10, m60]
m10_comp70 <- ASVBrayDist.mat[m10, m70]
m10_comp80 <- ASVBrayDist.mat[m10, m80]
m10_comp90 <- ASVBrayDist.mat[m10, m90]
m10_comp100 <- ASVBrayDist.mat[m10, m100]

# Distances between forest at 100m and each point along transect
m100_comp100 <- ASVBrayDist.mat[m100, m100]
m100_comp90 <- ASVBrayDist.mat[m100, m90]
m100_comp80 <- ASVBrayDist.mat[m100, m80]
m100_comp70 <- ASVBrayDist.mat[m100, m70]
m100_comp60 <- ASVBrayDist.mat[m100, m60]
m100_comp50 <- ASVBrayDist.mat[m100, m50]
m100_comp40 <- ASVBrayDist.mat[m100, m40]
m100_comp30 <- ASVBrayDist.mat[m100, m30]
m100_comp20 <- ASVBrayDist.mat[m100, m20]
m100_comp10 <- ASVBrayDist.mat[m100, m10]

######## MAKE BOXPLOTS ######## 
bold_a <- expression(bold("Dissimilarity Relative to 10m (patch)"))
bold_b <- expression(bold("Dissimilarity Relative to 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box <- boxplot(list(m10_comp10, m10_comp20, m10_comp30, m10_comp40,
                          m10_comp50, m10_comp60, m10_comp70,
                          m10_comp80, m10_comp90, m10_comp100),
                     ylab = "Bray-Curtis Dissimilarity",
                     names = c("10", "20", "30", "40", "50",
                               "60", "70", "80", "90", "100"), cex.axis = 0.8,
                     xlab = "meter along transect",
                     cex.lab = 1,
                     ylim=c(0.0, 1.0))
mtext(text=bold_a, side=3, adj = -0.065, line = 2)
forest_box <- boxplot(list(m100_comp10, m100_comp20, m100_comp30, m100_comp40,
                           m100_comp50, m100_comp60, m100_comp70,
                           m100_comp80, m100_comp90, m100_comp100), 
                      ylab = "Bray-Curtis Dissimilarity",
                      names = c("10", "20", "30", "40", "50",
                                "60", "70", "80", "90", "100"), cex.axis = 0.8,
                      xlab = "meter along transect",
                      cex.lab = 1,
                      ylim=c(0.0, 1.0))
mtext(text=bold_b, side=3, adj = -0.065, line = 2)

#########################################################################
# 3) the edge on the transect and all others along the transect in both directions
#########################################################################

###########
# EU 52
##########
# (using distances created in part 1)

# B-C distances between patch at 50m (edge) and each point along transect
m50_52_comp10 <- BrayDist_52.mat[m50_52, m10_52] #edge and 10m
m50_52_comp20 <- BrayDist_52.mat[m50_52, m20_52] #edge and 20m
m50_52_comp30 <- BrayDist_52.mat[m50_52, m30_52] #edge and 30m
m50_52_comp40 <- BrayDist_52.mat[m50_52, m40_52] #edge and 40m
m50_52_comp50 <- BrayDist_52.mat[m50_52, m50_52] #edge and self
m50_52_comp60 <- BrayDist_52.mat[m50_52, m60_52] #edge and 60m
m50_52_comp70 <- BrayDist_52.mat[m50_52, m70_52] #edge and 70m
m50_52_comp80 <- BrayDist_52.mat[m50_52, m80_52] #edge and 80m
m50_52_comp90 <- BrayDist_52.mat[m50_52, m90_52] #edge and 90m
m50_52_comp100 <- BrayDist_52.mat[m50_52, m100_52] #edge and 100m

bold_eu52edge <- expression(bold("EU 52: Dissimilarity from edge (50 m)"))
quartz()
edgecomp_box52 <- boxplot(list(m50_52_comp10, m50_52_comp20, m50_52_comp30, m50_52_comp40,
                            m50_52_comp50, m50_52_comp60, m50_52_comp70,
                            m50_52_comp80, m50_52_comp90, m50_52_comp100),
                       ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                       names = c("10", "20", "30", "40", "Edge (50 m)",
                                 "60", "70", "80", "90", "100"), cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu52edge, side=3, adj = -0.065, line = 2)

#########################################################################
# 4) The pooled forest samples versus all other points along the transect
# and 5) The pooled patch samples versus all other points along the transect
#########################################################################
EU_52_Soils.ps
sample_data(EU_52_Soils.ps)
pooledHabitat_52 <- merge_samples(EU_52_Soils.ps, "Habitat") #get samples pooled by habitat
pooledHabitat_52
pooledHabitat_52ASV <- psotu2veg(pooledHabitat_52)
rownames(pooledHabitat_52ASV)
pooledEdge_52ASV <- pooledHabitat_52ASV[1,] #pull out just pooled edge samples
pooledForest_52ASV <- pooledHabitat_52ASV[2,] #pull out pooled forest samples
pooledPatch_52ASV <- pooledHabitat_52ASV[3,] #pull out pooled patch samples

# Now to get dissimilarity matrices, I'll have to merge these pooled samples and with the
# non pooled samples of the other type

# Pooled forest versus all other samples along transect 
length(pooledForest_52ASV) #10,257 same as in whole EU_52_Soils.ps, even if some ASVs are zero
names(pooledForest_52ASV) #gives ASV names
summary(pooledForest_52ASV)
ASVtab_52 #from up above, this is ASV table of EU 52
ASVtab_52t <- t(ASVtab_52)
dim(ASVtab_52t)
#View(ASVtab_52t)
unique(names(pooledForest_52ASV) == rownames(ASVtab_52t)) 
# Because they match (line above), we can cbind to get new ASV table
pooledFvPtransect_52 <- t(cbind.data.frame(pooledForest_52ASV, ASVtab_52t))
#View(pooledFvPtransect_52) #needs to be samples as rows, ASVs as columns!

# Get Bray-Curtis dissimilarity matrix
pooledFvPtransect_52_Bray <- vegdist(pooledFvPtransect_52, method="bray")
pooledFvPtransect_52_Bray <- as.matrix(pooledFvPtransect_52_Bray)
diag(pooledFvPtransect_52_Bray) <- NA
View(pooledFvPtransect_52_Bray)
