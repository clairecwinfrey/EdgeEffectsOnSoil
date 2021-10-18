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

# from: https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_sep_variable.R
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
  # require(phyloseq)
  # require(plyr)
  
  ## Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }
  
  ## Extract sample meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  
  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  
  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }
  
  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  
  ## Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)
  
  ## Extract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })
  
  ## Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })
  
  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }
  
  return(res)
}

#################################################################
# SET UP FOR DOING 1 AND 2 BY TRANSECT AS OPPOSED TO ALL TOGETHER
# (analyses 1 and 2 are all together lower down)
#################################################################

# Separate out samples by EU and then by transect
# Then get B-C dissimilarities across these transects 
#### EU 52 #####
EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU)
EU_52_transSplit <- phyloseq_sep_variable(EU_52_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_52_T <- EU_52_transSplit$T
# Get just ASV table
ASVtab_52T <- psotu2veg(EU_52_T)
# View(ASVtab_52T) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_52T <- vegdist(ASVtab_52T, method="bray")
BrayDist_52T.mat <- as.matrix(BrayDist_52T)
diag(BrayDist_52T.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_52T.mat) == rownames(sample_data(EU_52_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps52T <- colnames(BrayDist_52T.mat) #this is the Mapping Sample ID of these
meter52T <- sample_data(EU_52_T)$Meter
index52T.df <- data.frame(samps52T, meter52T)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_52T <- which(index52T.df$meter52T == 10)
m20_52T <- which(index52T.df$meter52T == 20)
m30_52T <- which(index52T.df$meter52T == 30)
m40_52T <- which(index52T.df$meter52T == 40)
m50_52T <- which(index52T.df$meter52T == 50)
m60_52T <- which(index52T.df$meter52T == 60)
m70_52T <- which(index52T.df$meter52T == 70)
m80_52T <- which(index52T.df$meter52T == 80)
m90_52T <- which(index52T.df$meter52T == 90)
m100_52T <- which(index52T.df$meter52T == 100)

# B-C distances between patch at 10m and each point along transect
m10_52T_comp10 <- BrayDist_52T.mat[m10_52T, m10_52T] #add 10 and 10
m10_52T_comp20 <- BrayDist_52T.mat[m10_52T, m20_52T]
m10_52T_comp30 <- BrayDist_52T.mat[m10_52T, m30_52T]
m10_52T_comp40 <- BrayDist_52T.mat[m10_52T, m40_52T]
m10_52T_comp50 <- BrayDist_52T.mat[m10_52T, m50_52T]
m10_52T_comp60 <- BrayDist_52T.mat[m10_52T, m60_52T]
m10_52T_comp70 <- BrayDist_52T.mat[m10_52T, m70_52T]
m10_52T_comp80 <- BrayDist_52T.mat[m10_52T, m80_52T]
m10_52T_comp90 <- BrayDist_52T.mat[m10_52T, m90_52T]
m10_52T_comp100 <- BrayDist_52T.mat[m10_52T, m100_52T]

# B-C distances between forest at 100m and each point along transect
m100_52T_comp100 <- BrayDist_52T.mat[m100_52T, m100_52T] #add 100 and 100
m100_52T_comp90 <- BrayDist_52T.mat[m100_52T, m90_52T]
m100_52T_comp80 <- BrayDist_52T.mat[m100_52T, m80_52T]
m100_52T_comp70 <- BrayDist_52T.mat[m100_52T, m70_52T]
m100_52T_comp60 <- BrayDist_52T.mat[m100_52T, m60_52T]
m100_52T_comp50 <- BrayDist_52T.mat[m100_52T, m50_52T]
m100_52T_comp40 <- BrayDist_52T.mat[m100_52T, m40_52T]
m100_52T_comp30 <- BrayDist_52T.mat[m100_52T, m30_52T]
m100_52T_comp20 <- BrayDist_52T.mat[m100_52T, m20_52T]
m100_52T_comp10 <- BrayDist_52T.mat[m100_52T, m10_52T]


###############
### B ###
EU_52_B <- EU_52_transSplit$B
# Get just ASV table
ASVtab_52B <- psotu2veg(EU_52_B)
# View(ASVtab_52B) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_52B <- vegdist(ASVtab_52B, method="bray")
BrayDist_52B.mat <- as.matrix(BrayDist_52B)
diag(BrayDist_52B.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_52B.mat) == rownames(sample_data(EU_52_B))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps52B <- colnames(BrayDist_52B.mat) #this is the Mapping Sample ID of these
meter52B <- sample_data(EU_52_B)$Meter
index52B.df <- data.frame(samps52B, meter52B)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_52B <- which(index52B.df$meter52B == 10)
m20_52B <- which(index52B.df$meter52B == 20)
m30_52B <- which(index52B.df$meter52B == 30)
m40_52B <- which(index52B.df$meter52B == 40)
m50_52B <- which(index52B.df$meter52B == 50)
m60_52B <- which(index52B.df$meter52B == 60)
m70_52B <- which(index52B.df$meter52B == 70)
m80_52B <- which(index52B.df$meter52B == 80)
m90_52B <- which(index52B.df$meter52B == 90)
m100_52B <- which(index52B.df$meter52B == 100)

# B-C distances between patch at 10m and each point along transect
m10_52B_comp10 <- BrayDist_52B.mat[m10_52B, m10_52B] #add 10 and 10
m10_52B_comp20 <- BrayDist_52B.mat[m10_52B, m20_52B]
m10_52B_comp30 <- BrayDist_52B.mat[m10_52B, m30_52B]
m10_52B_comp40 <- BrayDist_52B.mat[m10_52B, m40_52B]
m10_52B_comp50 <- BrayDist_52B.mat[m10_52B, m50_52B]
m10_52B_comp60 <- BrayDist_52B.mat[m10_52B, m60_52B]
m10_52B_comp70 <- BrayDist_52B.mat[m10_52B, m70_52B]
m10_52B_comp80 <- BrayDist_52B.mat[m10_52B, m80_52B]
m10_52B_comp90 <- BrayDist_52B.mat[m10_52B, m90_52B]
m10_52B_comp100 <- BrayDist_52B.mat[m10_52B, m100_52B]

# B-C distances between forest at 100m and each point along transect
m100_52B_comp100 <- BrayDist_52B.mat[m100_52B, m100_52B] #add 100 and 100
m100_52B_comp90 <- BrayDist_52B.mat[m100_52B, m90_52B]
m100_52B_comp80 <- BrayDist_52B.mat[m100_52B, m80_52B]
m100_52B_comp70 <- BrayDist_52B.mat[m100_52B, m70_52B]
m100_52B_comp60 <- BrayDist_52B.mat[m100_52B, m60_52B]
m100_52B_comp50 <- BrayDist_52B.mat[m100_52B, m50_52B]
m100_52B_comp40 <- BrayDist_52B.mat[m100_52B, m40_52B]
m100_52B_comp30 <- BrayDist_52B.mat[m100_52B, m30_52B]
m100_52B_comp20 <- BrayDist_52B.mat[m100_52B, m20_52B]
m100_52B_comp10 <- BrayDist_52B.mat[m100_52B, m10_52B]
###############

### L ###
EU_52_L <- EU_52_transSplit$L
# Get just ASV table
ASVtab_52L <- psotu2veg(EU_52_L)
# View(ASVtab_52L) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_52L <- vegdist(ASVtab_52L, method="bray")
BrayDist_52L.mat <- as.matrix(BrayDist_52L)
diag(BrayDist_52L.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_52L.mat) == rownames(sample_data(EU_52_L))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps52L <- colnames(BrayDist_52L.mat) #this is the Mapping Sample ID of these
meter52L <- sample_data(EU_52_L)$Meter
index52L.df <- data.frame(samps52L, meter52L)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_52L <- which(index52L.df$meter52L == 10)
m20_52L <- which(index52L.df$meter52L == 20)
m30_52L <- which(index52L.df$meter52L == 30)
m40_52L <- which(index52L.df$meter52L == 40)
m50_52L <- which(index52L.df$meter52L == 50)
m60_52L <- which(index52L.df$meter52L == 60)
m70_52L <- which(index52L.df$meter52L == 70)
m80_52L <- which(index52L.df$meter52L == 80)
m90_52L <- which(index52L.df$meter52L == 90)
m100_52L <- which(index52L.df$meter52L == 100)

# B-C distances between patch at 10m and each point along transect
m10_52L_comp10 <- BrayDist_52L.mat[m10_52L, m10_52L] #add 10 and 10
m10_52L_comp20 <- BrayDist_52L.mat[m10_52L, m20_52L]
m10_52L_comp30 <- BrayDist_52L.mat[m10_52L, m30_52L]
m10_52L_comp40 <- BrayDist_52L.mat[m10_52L, m40_52L]
m10_52L_comp50 <- BrayDist_52L.mat[m10_52L, m50_52L]
m10_52L_comp60 <- BrayDist_52L.mat[m10_52L, m60_52L]
m10_52L_comp70 <- BrayDist_52L.mat[m10_52L, m70_52L]
m10_52L_comp80 <- BrayDist_52L.mat[m10_52L, m80_52L]
m10_52L_comp90 <- BrayDist_52L.mat[m10_52L, m90_52L]
m10_52L_comp100 <- BrayDist_52L.mat[m10_52L, m100_52L]

# B-C distances between forest at 100m and each point along transect
m100_52L_comp100 <- BrayDist_52L.mat[m100_52L, m100_52L] #add 100 and 100
m100_52L_comp90 <- BrayDist_52L.mat[m100_52L, m90_52L]
m100_52L_comp80 <- BrayDist_52L.mat[m100_52L, m80_52L]
m100_52L_comp70 <- BrayDist_52L.mat[m100_52L, m70_52L]
m100_52L_comp60 <- BrayDist_52L.mat[m100_52L, m60_52L]
m100_52L_comp50 <- BrayDist_52L.mat[m100_52L, m50_52L]
m100_52L_comp40 <- BrayDist_52L.mat[m100_52L, m40_52L]
m100_52L_comp30 <- BrayDist_52L.mat[m100_52L, m30_52L]
m100_52L_comp20 <- BrayDist_52L.mat[m100_52L, m20_52L]
m100_52L_comp10 <- BrayDist_52L.mat[m100_52L, m10_52L]
###############

### R ###
EU_52_R <- EU_52_transSplit$R
# Get just ASV table
ASVtab_52R <- psotu2veg(EU_52_R)
# View(ASVtab_52R) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_52R <- vegdist(ASVtab_52R, method="bray")
BrayDist_52R.mat <- as.matrix(BrayDist_52R)
diag(BrayDist_52R.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_52R.mat) == rownames(sample_data(EU_52_R))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps52R <- colnames(BrayDist_52R.mat) #this is the Mapping Sample ID of these
meter52R <- sample_data(EU_52_R)$Meter
index52R.df <- data.frame(samps52R, meter52R)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_52R <- which(index52R.df$meter52R == 10)
m20_52R <- which(index52R.df$meter52R == 20)
m30_52R <- which(index52R.df$meter52R == 30)
m40_52R <- which(index52R.df$meter52R == 40)
m50_52R <- which(index52R.df$meter52R == 50)
m60_52R <- which(index52R.df$meter52R == 60)
m70_52R <- which(index52R.df$meter52R == 70)
m80_52R <- which(index52R.df$meter52R == 80)
m90_52R <- which(index52R.df$meter52R == 90)
m100_52R <- which(index52R.df$meter52R == 100)

# B-C distances between patch at 10m and each point along transect
m10_52R_comp10 <- BrayDist_52R.mat[m10_52R, m10_52R] #add 10 and 10
m10_52R_comp20 <- BrayDist_52R.mat[m10_52R, m20_52R]
m10_52R_comp30 <- BrayDist_52R.mat[m10_52R, m30_52R]
m10_52R_comp40 <- BrayDist_52R.mat[m10_52R, m40_52R]
m10_52R_comp50 <- BrayDist_52R.mat[m10_52R, m50_52R]
m10_52R_comp60 <- BrayDist_52R.mat[m10_52R, m60_52R]
m10_52R_comp70 <- BrayDist_52R.mat[m10_52R, m70_52R]
m10_52R_comp80 <- BrayDist_52R.mat[m10_52R, m80_52R]
m10_52R_comp90 <- BrayDist_52R.mat[m10_52R, m90_52R]
m10_52R_comp100 <- BrayDist_52R.mat[m10_52R, m100_52R]

# B-C distances between forest at 100m and each point along transect
m100_52R_comp100 <- BrayDist_52R.mat[m100_52R, m100_52R] #add 100 and 100
m100_52R_comp90 <- BrayDist_52R.mat[m100_52R, m90_52R]
m100_52R_comp80 <- BrayDist_52R.mat[m100_52R, m80_52R]
m100_52R_comp70 <- BrayDist_52R.mat[m100_52R, m70_52R]
m100_52R_comp60 <- BrayDist_52R.mat[m100_52R, m60_52R]
m100_52R_comp50 <- BrayDist_52R.mat[m100_52R, m50_52R]
m100_52R_comp40 <- BrayDist_52R.mat[m100_52R, m40_52R]
m100_52R_comp30 <- BrayDist_52R.mat[m100_52R, m30_52R]
m100_52R_comp20 <- BrayDist_52R.mat[m100_52R, m20_52R]
m100_52R_comp10 <- BrayDist_52R.mat[m100_52R, m10_52R]
###############

# EU 53N
EU_53N_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53N")
unique(sample_data(EU_53N_Soils.ps)$EU)
EU_53N_transSplit <- phyloseq_sep_variable(EU_53N_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_53N_T <- EU_53N_transSplit$T
# Get just ASV table
ASVtab_53NT <- psotu2veg(EU_53N_T)
# View(ASVtab_53NT) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53NT <- vegdist(ASVtab_53NT, method="bray")
BrayDist_53NT.mat <- as.matrix(BrayDist_53NT)
diag(BrayDist_53NT.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53NT.mat) == rownames(sample_data(EU_53N_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53NT <- colnames(BrayDist_53NT.mat) #this is the Mapping Sample ID of these
meter53NT <- sample_data(EU_53N_T)$Meter
index53NT.df <- data.frame(samps53NT, meter53NT)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53NT <- which(index53NT.df$meter53NT == 10)
m20_53NT <- which(index53NT.df$meter53NT == 20)
m30_53NT <- which(index53NT.df$meter53NT == 30)
m40_53NT <- which(index53NT.df$meter53NT == 40)
m50_53NT <- which(index53NT.df$meter53NT == 50)
m60_53NT <- which(index53NT.df$meter53NT == 60)
m70_53NT <- which(index53NT.df$meter53NT == 70)
m80_53NT <- which(index53NT.df$meter53NT == 80)
m90_53NT <- which(index53NT.df$meter53NT == 90)
m100_53NT <- which(index53NT.df$meter53NT == 100)

# B-C distances between patch at 10m and each point along transect
m10_53NT_comp10 <- BrayDist_53NT.mat[m10_53NT, m10_53NT] #add 10 and 10
m10_53NT_comp20 <- BrayDist_53NT.mat[m10_53NT, m20_53NT]
m10_53NT_comp30 <- BrayDist_53NT.mat[m10_53NT, m30_53NT]
m10_53NT_comp40 <- BrayDist_53NT.mat[m10_53NT, m40_53NT]
m10_53NT_comp50 <- BrayDist_53NT.mat[m10_53NT, m50_53NT]
m10_53NT_comp60 <- BrayDist_53NT.mat[m10_53NT, m60_53NT]
m10_53NT_comp70 <- BrayDist_53NT.mat[m10_53NT, m70_53NT]
m10_53NT_comp80 <- BrayDist_53NT.mat[m10_53NT, m80_53NT]
m10_53NT_comp90 <- BrayDist_53NT.mat[m10_53NT, m90_53NT]
m10_53NT_comp100 <- BrayDist_53NT.mat[m10_53NT, m100_53NT]

# B-C distances between forest at 100m and each point along transect
m100_53NT_comp100 <- BrayDist_53NT.mat[m100_53NT, m100_53NT] #add 100 and 100
m100_53NT_comp90 <- BrayDist_53NT.mat[m100_53NT, m90_53NT]
m100_53NT_comp80 <- BrayDist_53NT.mat[m100_53NT, m80_53NT]
m100_53NT_comp70 <- BrayDist_53NT.mat[m100_53NT, m70_53NT]
m100_53NT_comp60 <- BrayDist_53NT.mat[m100_53NT, m60_53NT]
m100_53NT_comp50 <- BrayDist_53NT.mat[m100_53NT, m50_53NT]
m100_53NT_comp40 <- BrayDist_53NT.mat[m100_53NT, m40_53NT]
m100_53NT_comp30 <- BrayDist_53NT.mat[m100_53NT, m30_53NT]
m100_53NT_comp20 <- BrayDist_53NT.mat[m100_53NT, m20_53NT]
m100_53NT_comp10 <- BrayDist_53NT.mat[m100_53NT, m10_53NT]

###############
### B ###
EU_53N_B <- EU_53N_transSplit$B
# Get just ASV table
ASVtab_53NB <- psotu2veg(EU_53N_B)
# View(ASVtab_53NB) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53NB <- vegdist(ASVtab_53NB, method="bray")
BrayDist_53NB.mat <- as.matrix(BrayDist_53NB)
diag(BrayDist_53NB.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53NB.mat) == rownames(sample_data(EU_53N_B))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53NB <- colnames(BrayDist_53NB.mat) #this is the Mapping Sample ID of these
meter53NB <- sample_data(EU_53N_B)$Meter
index53NB.df <- data.frame(samps53NB, meter53NB)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53NB <- which(index53NB.df$meter53NB == 10)
m20_53NB <- which(index53NB.df$meter53NB == 20)
m30_53NB <- which(index53NB.df$meter53NB == 30)
m40_53NB <- which(index53NB.df$meter53NB == 40)
m50_53NB <- which(index53NB.df$meter53NB == 50)
m60_53NB <- which(index53NB.df$meter53NB == 60)
m70_53NB <- which(index53NB.df$meter53NB == 70)
m80_53NB <- which(index53NB.df$meter53NB == 80)
m90_53NB <- which(index53NB.df$meter53NB == 90)
m100_53NB <- which(index53NB.df$meter53NB == 100)

# B-C distances between patch at 10m and each point along transect
m10_53NB_comp10 <- BrayDist_53NB.mat[m10_53NB, m10_53NB] #add 10 and 10
m10_53NB_comp20 <- BrayDist_53NB.mat[m10_53NB, m20_53NB]
m10_53NB_comp30 <- BrayDist_53NB.mat[m10_53NB, m30_53NB]
m10_53NB_comp40 <- BrayDist_53NB.mat[m10_53NB, m40_53NB]
m10_53NB_comp50 <- BrayDist_53NB.mat[m10_53NB, m50_53NB]
m10_53NB_comp60 <- BrayDist_53NB.mat[m10_53NB, m60_53NB]
m10_53NB_comp70 <- BrayDist_53NB.mat[m10_53NB, m70_53NB]
m10_53NB_comp80 <- BrayDist_53NB.mat[m10_53NB, m80_53NB]
m10_53NB_comp90 <- BrayDist_53NB.mat[m10_53NB, m90_53NB]
m10_53NB_comp100 <- BrayDist_53NB.mat[m10_53NB, m100_53NB]

# B-C distances between forest at 100m and each point along transect
m100_53NB_comp100 <- BrayDist_53NB.mat[m100_53NB, m100_53NB] #add 100 and 100
m100_53NB_comp90 <- BrayDist_53NB.mat[m100_53NB, m90_53NB]
m100_53NB_comp80 <- BrayDist_53NB.mat[m100_53NB, m80_53NB]
m100_53NB_comp70 <- BrayDist_53NB.mat[m100_53NB, m70_53NB]
m100_53NB_comp60 <- BrayDist_53NB.mat[m100_53NB, m60_53NB]
m100_53NB_comp50 <- BrayDist_53NB.mat[m100_53NB, m50_53NB]
m100_53NB_comp40 <- BrayDist_53NB.mat[m100_53NB, m40_53NB]
m100_53NB_comp30 <- BrayDist_53NB.mat[m100_53NB, m30_53NB]
m100_53NB_comp20 <- BrayDist_53NB.mat[m100_53NB, m20_53NB]
m100_53NB_comp10 <- BrayDist_53NB.mat[m100_53NB, m10_53NB]
###############

### L ###
EU_53N_L <- EU_53N_transSplit$L
# Get just ASV table
ASVtab_53NL <- psotu2veg(EU_53N_L)
# View(ASVtab_53NL) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53NL <- vegdist(ASVtab_53NL, method="bray")
BrayDist_53NL.mat <- as.matrix(BrayDist_53NL)
diag(BrayDist_53NL.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53NL.mat) == rownames(sample_data(EU_53N_L))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53NL <- colnames(BrayDist_53NL.mat) #this is the Mapping Sample ID of these
meter53NL <- sample_data(EU_53N_L)$Meter
index53NL.df <- data.frame(samps53NL, meter53NL)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53NL <- which(index53NL.df$meter53NL == 10)
m20_53NL <- which(index53NL.df$meter53NL == 20)
m30_53NL <- which(index53NL.df$meter53NL == 30)
m40_53NL <- which(index53NL.df$meter53NL == 40)
m50_53NL <- which(index53NL.df$meter53NL == 50)
m60_53NL <- which(index53NL.df$meter53NL == 60)
m70_53NL <- which(index53NL.df$meter53NL == 70)
m80_53NL <- which(index53NL.df$meter53NL == 80)
m90_53NL <- which(index53NL.df$meter53NL == 90)
m100_53NL <- which(index53NL.df$meter53NL == 100)

# B-C distances between patch at 10m and each point along transect
m10_53NL_comp10 <- BrayDist_53NL.mat[m10_53NL, m10_53NL] #add 10 and 10
m10_53NL_comp20 <- BrayDist_53NL.mat[m10_53NL, m20_53NL]
m10_53NL_comp30 <- BrayDist_53NL.mat[m10_53NL, m30_53NL]
m10_53NL_comp40 <- BrayDist_53NL.mat[m10_53NL, m40_53NL]
m10_53NL_comp50 <- BrayDist_53NL.mat[m10_53NL, m50_53NL]
m10_53NL_comp60 <- BrayDist_53NL.mat[m10_53NL, m60_53NL]
m10_53NL_comp70 <- BrayDist_53NL.mat[m10_53NL, m70_53NL]
m10_53NL_comp80 <- BrayDist_53NL.mat[m10_53NL, m80_53NL]
m10_53NL_comp90 <- BrayDist_53NL.mat[m10_53NL, m90_53NL]
m10_53NL_comp100 <- BrayDist_53NL.mat[m10_53NL, m100_53NL]

# B-C distances between forest at 100m and each point along transect
m100_53NL_comp100 <- BrayDist_53NL.mat[m100_53NL, m100_53NL] #add 100 and 100
m100_53NL_comp90 <- BrayDist_53NL.mat[m100_53NL, m90_53NL]
m100_53NL_comp80 <- BrayDist_53NL.mat[m100_53NL, m80_53NL]
m100_53NL_comp70 <- BrayDist_53NL.mat[m100_53NL, m70_53NL]
m100_53NL_comp60 <- BrayDist_53NL.mat[m100_53NL, m60_53NL]
m100_53NL_comp50 <- BrayDist_53NL.mat[m100_53NL, m50_53NL]
m100_53NL_comp40 <- BrayDist_53NL.mat[m100_53NL, m40_53NL]
m100_53NL_comp30 <- BrayDist_53NL.mat[m100_53NL, m30_53NL]
m100_53NL_comp20 <- BrayDist_53NL.mat[m100_53NL, m20_53NL]
m100_53NL_comp10 <- BrayDist_53NL.mat[m100_53NL, m10_53NL]
###############

### R ###
EU_53N_R <- EU_53N_transSplit$R
# Get just ASV table
ASVtab_53NR <- psotu2veg(EU_53N_R)
# View(ASVtab_53NR) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53NR <- vegdist(ASVtab_53NR, method="bray")
BrayDist_53NR.mat <- as.matrix(BrayDist_53NR)
diag(BrayDist_53NR.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53NR.mat) == rownames(sample_data(EU_53N_R))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53NR <- colnames(BrayDist_53NR.mat) #this is the Mapping Sample ID of these
meter53NR <- sample_data(EU_53N_R)$Meter
index53NR.df <- data.frame(samps53NR, meter53NR)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53NR <- which(index53NR.df$meter53NR == 10)
m20_53NR <- which(index53NR.df$meter53NR == 20)
m30_53NR <- which(index53NR.df$meter53NR == 30)
m40_53NR <- which(index53NR.df$meter53NR == 40)
m50_53NR <- which(index53NR.df$meter53NR == 50)
m60_53NR <- which(index53NR.df$meter53NR == 60)
m70_53NR <- which(index53NR.df$meter53NR == 70)
m80_53NR <- which(index53NR.df$meter53NR == 80)
m90_53NR <- which(index53NR.df$meter53NR == 90)
m100_53NR <- which(index53NR.df$meter53NR == 100)

# B-C distances between patch at 10m and each point along transect
m10_53NR_comp10 <- BrayDist_53NR.mat[m10_53NR, m10_53NR] #add 10 and 10
m10_53NR_comp20 <- BrayDist_53NR.mat[m10_53NR, m20_53NR]
m10_53NR_comp30 <- BrayDist_53NR.mat[m10_53NR, m30_53NR]
m10_53NR_comp40 <- BrayDist_53NR.mat[m10_53NR, m40_53NR]
m10_53NR_comp50 <- BrayDist_53NR.mat[m10_53NR, m50_53NR]
m10_53NR_comp60 <- BrayDist_53NR.mat[m10_53NR, m60_53NR]
m10_53NR_comp70 <- BrayDist_53NR.mat[m10_53NR, m70_53NR]
m10_53NR_comp80 <- BrayDist_53NR.mat[m10_53NR, m80_53NR]
m10_53NR_comp90 <- BrayDist_53NR.mat[m10_53NR, m90_53NR]
m10_53NR_comp100 <- BrayDist_53NR.mat[m10_53NR, m100_53NR]

# B-C distances between forest at 100m and each point along transect
m100_53NR_comp100 <- BrayDist_53NR.mat[m100_53NR, m100_53NR] #add 100 and 100
m100_53NR_comp90 <- BrayDist_53NR.mat[m100_53NR, m90_53NR]
m100_53NR_comp80 <- BrayDist_53NR.mat[m100_53NR, m80_53NR]
m100_53NR_comp70 <- BrayDist_53NR.mat[m100_53NR, m70_53NR]
m100_53NR_comp60 <- BrayDist_53NR.mat[m100_53NR, m60_53NR]
m100_53NR_comp50 <- BrayDist_53NR.mat[m100_53NR, m50_53NR]
m100_53NR_comp40 <- BrayDist_53NR.mat[m100_53NR, m40_53NR]
m100_53NR_comp30 <- BrayDist_53NR.mat[m100_53NR, m30_53NR]
m100_53NR_comp20 <- BrayDist_53NR.mat[m100_53NR, m20_53NR]
m100_53NR_comp10 <- BrayDist_53NR.mat[m100_53NR, m10_53NR]
###############

# EU 54S
EU_54S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_54S")
unique(sample_data(EU_54S_Soils.ps)$EU)
EU_54S_transSplit <- phyloseq_sep_variable(EU_54S_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_54S_T <- EU_54S_transSplit$T
# Get just ASV table
ASVtab_54ST <- psotu2veg(EU_54S_T)
# View(ASVtab_54ST) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_54ST <- vegdist(ASVtab_54ST, method="bray")
BrayDist_54ST.mat <- as.matrix(BrayDist_54ST)
diag(BrayDist_54ST.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_54ST.mat) == rownames(sample_data(EU_54S_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps54ST <- colnames(BrayDist_54ST.mat) #this is the Mapping Sample ID of these
meter54ST <- sample_data(EU_54S_T)$Meter
index54ST.df <- data.frame(samps54ST, meter54ST)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_54ST <- which(index54ST.df$meter54ST == 10)
m20_54ST <- which(index54ST.df$meter54ST == 20)
m30_54ST <- which(index54ST.df$meter54ST == 30)
m40_54ST <- which(index54ST.df$meter54ST == 40)
m50_54ST <- which(index54ST.df$meter54ST == 50)
m60_54ST <- which(index54ST.df$meter54ST == 60)
m70_54ST <- which(index54ST.df$meter54ST == 70)
m80_54ST <- which(index54ST.df$meter54ST == 80)
m90_54ST <- which(index54ST.df$meter54ST == 90)
m100_54ST <- which(index54ST.df$meter54ST == 100)

# B-C distances between patch at 10m and each point along transect
m10_54ST_comp10 <- BrayDist_54ST.mat[m10_54ST, m10_54ST] #add 10 and 10
m10_54ST_comp20 <- BrayDist_54ST.mat[m10_54ST, m20_54ST]
m10_54ST_comp30 <- BrayDist_54ST.mat[m10_54ST, m30_54ST]
m10_54ST_comp40 <- BrayDist_54ST.mat[m10_54ST, m40_54ST]
m10_54ST_comp50 <- BrayDist_54ST.mat[m10_54ST, m50_54ST]
m10_54ST_comp60 <- BrayDist_54ST.mat[m10_54ST, m60_54ST]
m10_54ST_comp70 <- BrayDist_54ST.mat[m10_54ST, m70_54ST]
m10_54ST_comp80 <- BrayDist_54ST.mat[m10_54ST, m80_54ST]
m10_54ST_comp90 <- BrayDist_54ST.mat[m10_54ST, m90_54ST]
m10_54ST_comp100 <- BrayDist_54ST.mat[m10_54ST, m100_54ST]

# B-C distances between forest at 100m and each point along transect
m100_54ST_comp100 <- BrayDist_54ST.mat[m100_54ST, m100_54ST] #add 100 and 100
m100_54ST_comp90 <- BrayDist_54ST.mat[m100_54ST, m90_54ST]
m100_54ST_comp80 <- BrayDist_54ST.mat[m100_54ST, m80_54ST]
m100_54ST_comp70 <- BrayDist_54ST.mat[m100_54ST, m70_54ST]
m100_54ST_comp60 <- BrayDist_54ST.mat[m100_54ST, m60_54ST]
m100_54ST_comp50 <- BrayDist_54ST.mat[m100_54ST, m50_54ST]
m100_54ST_comp40 <- BrayDist_54ST.mat[m100_54ST, m40_54ST]
m100_54ST_comp30 <- BrayDist_54ST.mat[m100_54ST, m30_54ST]
m100_54ST_comp20 <- BrayDist_54ST.mat[m100_54ST, m20_54ST]
m100_54ST_comp10 <- BrayDist_54ST.mat[m100_54ST, m10_54ST]


###############
### B ###
EU_54S_B <- EU_54S_transSplit$B
# Get just ASV table
ASVtab_54SB <- psotu2veg(EU_54S_B)
# View(ASVtab_54SB) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_54SB <- vegdist(ASVtab_54SB, method="bray")
BrayDist_54SB.mat <- as.matrix(BrayDist_54SB)
diag(BrayDist_54SB.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_54SB.mat) == rownames(sample_data(EU_54S_B))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps54SB <- colnames(BrayDist_54SB.mat) #this is the Mapping Sample ID of these
meter54SB <- sample_data(EU_54S_B)$Meter
index54SB.df <- data.frame(samps54SB, meter54SB)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_54SB <- which(index54SB.df$meter54SB == 10)
m20_54SB <- which(index54SB.df$meter54SB == 20)
m30_54SB <- which(index54SB.df$meter54SB == 30)
m40_54SB <- which(index54SB.df$meter54SB == 40)
m50_54SB <- which(index54SB.df$meter54SB == 50)
m60_54SB <- which(index54SB.df$meter54SB == 60)
m70_54SB <- which(index54SB.df$meter54SB == 70)
m80_54SB <- which(index54SB.df$meter54SB == 80)
m90_54SB <- which(index54SB.df$meter54SB == 90)
m100_54SB <- which(index54SB.df$meter54SB == 100)

# B-C distances between patch at 10m and each point along transect
m10_54SB_comp10 <- BrayDist_54SB.mat[m10_54SB, m10_54SB] #add 10 and 10
m10_54SB_comp20 <- BrayDist_54SB.mat[m10_54SB, m20_54SB]
m10_54SB_comp30 <- BrayDist_54SB.mat[m10_54SB, m30_54SB]
m10_54SB_comp40 <- BrayDist_54SB.mat[m10_54SB, m40_54SB]
m10_54SB_comp50 <- BrayDist_54SB.mat[m10_54SB, m50_54SB]
m10_54SB_comp60 <- BrayDist_54SB.mat[m10_54SB, m60_54SB]
m10_54SB_comp70 <- BrayDist_54SB.mat[m10_54SB, m70_54SB]
m10_54SB_comp80 <- BrayDist_54SB.mat[m10_54SB, m80_54SB]
m10_54SB_comp90 <- BrayDist_54SB.mat[m10_54SB, m90_54SB]
m10_54SB_comp100 <- BrayDist_54SB.mat[m10_54SB, m100_54SB]

# B-C distances between forest at 100m and each point along transect
m100_54SB_comp100 <- BrayDist_54SB.mat[m100_54SB, m100_54SB] #add 100 and 100
m100_54SB_comp90 <- BrayDist_54SB.mat[m100_54SB, m90_54SB]
m100_54SB_comp80 <- BrayDist_54SB.mat[m100_54SB, m80_54SB]
m100_54SB_comp70 <- BrayDist_54SB.mat[m100_54SB, m70_54SB]
m100_54SB_comp60 <- BrayDist_54SB.mat[m100_54SB, m60_54SB]
m100_54SB_comp50 <- BrayDist_54SB.mat[m100_54SB, m50_54SB]
m100_54SB_comp40 <- BrayDist_54SB.mat[m100_54SB, m40_54SB]
m100_54SB_comp30 <- BrayDist_54SB.mat[m100_54SB, m30_54SB]
m100_54SB_comp20 <- BrayDist_54SB.mat[m100_54SB, m20_54SB]
m100_54SB_comp10 <- BrayDist_54SB.mat[m100_54SB, m10_54SB]
###############

### L ###
EU_54S_L <- EU_54S_transSplit$L
# Get just ASV table
ASVtab_54SL <- psotu2veg(EU_54S_L)
# View(ASVtab_54SL) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_54SL <- vegdist(ASVtab_54SL, method="bray")
BrayDist_54SL.mat <- as.matrix(BrayDist_54SL)
diag(BrayDist_54SL.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_54SL.mat) == rownames(sample_data(EU_54S_L))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps54SL <- colnames(BrayDist_54SL.mat) #this is the Mapping Sample ID of these
meter54SL <- sample_data(EU_54S_L)$Meter
index54SL.df <- data.frame(samps54SL, meter54SL)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_54SL <- which(index54SL.df$meter54SL == 10)
m20_54SL <- which(index54SL.df$meter54SL == 20)
m30_54SL <- which(index54SL.df$meter54SL == 30)
m40_54SL <- which(index54SL.df$meter54SL == 40)
m50_54SL <- which(index54SL.df$meter54SL == 50)
m60_54SL <- which(index54SL.df$meter54SL == 60)
m70_54SL <- which(index54SL.df$meter54SL == 70)
m80_54SL <- which(index54SL.df$meter54SL == 80)
m90_54SL <- which(index54SL.df$meter54SL == 90)
m100_54SL <- which(index54SL.df$meter54SL == 100)

# B-C distances between patch at 10m and each point along transect
m10_54SL_comp10 <- BrayDist_54SL.mat[m10_54SL, m10_54SL] #add 10 and 10
m10_54SL_comp20 <- BrayDist_54SL.mat[m10_54SL, m20_54SL]
m10_54SL_comp30 <- BrayDist_54SL.mat[m10_54SL, m30_54SL]
m10_54SL_comp40 <- BrayDist_54SL.mat[m10_54SL, m40_54SL]
m10_54SL_comp50 <- BrayDist_54SL.mat[m10_54SL, m50_54SL]
m10_54SL_comp60 <- BrayDist_54SL.mat[m10_54SL, m60_54SL]
m10_54SL_comp70 <- BrayDist_54SL.mat[m10_54SL, m70_54SL]
m10_54SL_comp80 <- BrayDist_54SL.mat[m10_54SL, m80_54SL]
m10_54SL_comp90 <- BrayDist_54SL.mat[m10_54SL, m90_54SL]
m10_54SL_comp100 <- BrayDist_54SL.mat[m10_54SL, m100_54SL]

# B-C distances between forest at 100m and each point along transect
m100_54SL_comp100 <- BrayDist_54SL.mat[m100_54SL, m100_54SL] #add 100 and 100
m100_54SL_comp90 <- BrayDist_54SL.mat[m100_54SL, m90_54SL]
m100_54SL_comp80 <- BrayDist_54SL.mat[m100_54SL, m80_54SL]
m100_54SL_comp70 <- BrayDist_54SL.mat[m100_54SL, m70_54SL]
m100_54SL_comp60 <- BrayDist_54SL.mat[m100_54SL, m60_54SL]
m100_54SL_comp50 <- BrayDist_54SL.mat[m100_54SL, m50_54SL]
m100_54SL_comp40 <- BrayDist_54SL.mat[m100_54SL, m40_54SL]
m100_54SL_comp30 <- BrayDist_54SL.mat[m100_54SL, m30_54SL]
m100_54SL_comp20 <- BrayDist_54SL.mat[m100_54SL, m20_54SL]
m100_54SL_comp10 <- BrayDist_54SL.mat[m100_54SL, m10_54SL]
###############

### R ###
EU_54S_R <- EU_54S_transSplit$R
# Get just ASV table
ASVtab_54SR <- psotu2veg(EU_54S_R)
# View(ASVtab_54SR) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_54SR <- vegdist(ASVtab_54SR, method="bray")
BrayDist_54SR.mat <- as.matrix(BrayDist_54SR)
diag(BrayDist_54SR.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_54SR.mat) == rownames(sample_data(EU_54S_R))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps54SR <- colnames(BrayDist_54SR.mat) #this is the Mapping Sample ID of these
meter54SR <- sample_data(EU_54S_R)$Meter
index54SR.df <- data.frame(samps54SR, meter54SR)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_54SR <- which(index54SR.df$meter54SR == 10)
m20_54SR <- which(index54SR.df$meter54SR == 20)
m30_54SR <- which(index54SR.df$meter54SR == 30)
m40_54SR <- which(index54SR.df$meter54SR == 40)
m50_54SR <- which(index54SR.df$meter54SR == 50)
m60_54SR <- which(index54SR.df$meter54SR == 60)
m70_54SR <- which(index54SR.df$meter54SR == 70)
m80_54SR <- which(index54SR.df$meter54SR == 80)
m90_54SR <- which(index54SR.df$meter54SR == 90)
m100_54SR <- which(index54SR.df$meter54SR == 100)

# B-C distances between patch at 10m and each point along transect
m10_54SR_comp10 <- BrayDist_54SR.mat[m10_54SR, m10_54SR] #add 10 and 10
m10_54SR_comp20 <- BrayDist_54SR.mat[m10_54SR, m20_54SR]
m10_54SR_comp30 <- BrayDist_54SR.mat[m10_54SR, m30_54SR]
m10_54SR_comp40 <- BrayDist_54SR.mat[m10_54SR, m40_54SR]
m10_54SR_comp50 <- BrayDist_54SR.mat[m10_54SR, m50_54SR]
m10_54SR_comp60 <- BrayDist_54SR.mat[m10_54SR, m60_54SR]
m10_54SR_comp70 <- BrayDist_54SR.mat[m10_54SR, m70_54SR]
m10_54SR_comp80 <- BrayDist_54SR.mat[m10_54SR, m80_54SR]
m10_54SR_comp90 <- BrayDist_54SR.mat[m10_54SR, m90_54SR]
m10_54SR_comp100 <- BrayDist_54SR.mat[m10_54SR, m100_54SR]

# B-C distances between forest at 100m and each point along transect
m100_54SR_comp100 <- BrayDist_54SR.mat[m100_54SR, m100_54SR] #add 100 and 100
m100_54SR_comp90 <- BrayDist_54SR.mat[m100_54SR, m90_54SR]
m100_54SR_comp80 <- BrayDist_54SR.mat[m100_54SR, m80_54SR]
m100_54SR_comp70 <- BrayDist_54SR.mat[m100_54SR, m70_54SR]
m100_54SR_comp60 <- BrayDist_54SR.mat[m100_54SR, m60_54SR]
m100_54SR_comp50 <- BrayDist_54SR.mat[m100_54SR, m50_54SR]
m100_54SR_comp40 <- BrayDist_54SR.mat[m100_54SR, m40_54SR]
m100_54SR_comp30 <- BrayDist_54SR.mat[m100_54SR, m30_54SR]
m100_54SR_comp20 <- BrayDist_54SR.mat[m100_54SR, m20_54SR]
m100_54SR_comp10 <- BrayDist_54SR.mat[m100_54SR, m10_54SR]
###############

# EU 8
EU_8_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_8")
unique(sample_data(EU_8_Soils.ps)$EU)
EU_8_transSplit <- phyloseq_sep_variable(EU_8_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_8_T <- EU_8_transSplit$T
# Get just ASV table
ASVtab_8T <- psotu2veg(EU_8_T)
# View(ASVtab_8T) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_8T <- vegdist(ASVtab_8T, method="bray")
BrayDist_8T.mat <- as.matrix(BrayDist_8T)
diag(BrayDist_8T.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_8T.mat) == rownames(sample_data(EU_8_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps8T <- colnames(BrayDist_8T.mat) #this is the Mapping Sample ID of these
meter8T <- sample_data(EU_8_T)$Meter
index8T.df <- data.frame(samps8T, meter8T)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_8T <- which(index8T.df$meter8T == 10)
m20_8T <- which(index8T.df$meter8T == 20)
m30_8T <- which(index8T.df$meter8T == 30)
m40_8T <- which(index8T.df$meter8T == 40)
m50_8T <- which(index8T.df$meter8T == 50)
m60_8T <- which(index8T.df$meter8T == 60)
m70_8T <- which(index8T.df$meter8T == 70)
m80_8T <- which(index8T.df$meter8T == 80)
m90_8T <- which(index8T.df$meter8T == 90)
m100_8T <- which(index8T.df$meter8T == 100)

# B-C distances between patch at 10m and each point along transect
m10_8T_comp10 <- BrayDist_8T.mat[m10_8T, m10_8T] #add 10 and 10
m10_8T_comp20 <- BrayDist_8T.mat[m10_8T, m20_8T]
m10_8T_comp30 <- BrayDist_8T.mat[m10_8T, m30_8T]
m10_8T_comp40 <- BrayDist_8T.mat[m10_8T, m40_8T]
m10_8T_comp50 <- BrayDist_8T.mat[m10_8T, m50_8T]
m10_8T_comp60 <- BrayDist_8T.mat[m10_8T, m60_8T]
m10_8T_comp70 <- BrayDist_8T.mat[m10_8T, m70_8T]
m10_8T_comp80 <- BrayDist_8T.mat[m10_8T, m80_8T]
m10_8T_comp90 <- BrayDist_8T.mat[m10_8T, m90_8T]
m10_8T_comp100 <- BrayDist_8T.mat[m10_8T, m100_8T]

# B-C distances between forest at 100m and each point along transect
m100_8T_comp100 <- BrayDist_8T.mat[m100_8T, m100_8T] #add 100 and 100
m100_8T_comp90 <- BrayDist_8T.mat[m100_8T, m90_8T]
m100_8T_comp80 <- BrayDist_8T.mat[m100_8T, m80_8T]
m100_8T_comp70 <- BrayDist_8T.mat[m100_8T, m70_8T]
m100_8T_comp60 <- BrayDist_8T.mat[m100_8T, m60_8T]
m100_8T_comp50 <- BrayDist_8T.mat[m100_8T, m50_8T]
m100_8T_comp40 <- BrayDist_8T.mat[m100_8T, m40_8T]
m100_8T_comp30 <- BrayDist_8T.mat[m100_8T, m30_8T]
m100_8T_comp20 <- BrayDist_8T.mat[m100_8T, m20_8T]
m100_8T_comp10 <- BrayDist_8T.mat[m100_8T, m10_8T]


###############
### B ###
EU_8_B <- EU_8_transSplit$B
# Get just ASV table
ASVtab_8B <- psotu2veg(EU_8_B)
# View(ASVtab_8B) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_8B <- vegdist(ASVtab_8B, method="bray")
BrayDist_8B.mat <- as.matrix(BrayDist_8B)
diag(BrayDist_8B.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_8B.mat) == rownames(sample_data(EU_8_B))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps8B <- colnames(BrayDist_8B.mat) #this is the Mapping Sample ID of these
meter8B <- sample_data(EU_8_B)$Meter
index8B.df <- data.frame(samps8B, meter8B)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_8B <- which(index8B.df$meter8B == 10)
m20_8B <- which(index8B.df$meter8B == 20)
m30_8B <- which(index8B.df$meter8B == 30)
m40_8B <- which(index8B.df$meter8B == 40)
m50_8B <- which(index8B.df$meter8B == 50)
m60_8B <- which(index8B.df$meter8B == 60)
m70_8B <- which(index8B.df$meter8B == 70)
m80_8B <- which(index8B.df$meter8B == 80)
m90_8B <- which(index8B.df$meter8B == 90)
m100_8B <- which(index8B.df$meter8B == 100)

# B-C distances between patch at 10m and each point along transect
m10_8B_comp10 <- BrayDist_8B.mat[m10_8B, m10_8B] #add 10 and 10
m10_8B_comp20 <- BrayDist_8B.mat[m10_8B, m20_8B]
m10_8B_comp30 <- BrayDist_8B.mat[m10_8B, m30_8B]
m10_8B_comp40 <- BrayDist_8B.mat[m10_8B, m40_8B]
m10_8B_comp50 <- BrayDist_8B.mat[m10_8B, m50_8B]
m10_8B_comp60 <- BrayDist_8B.mat[m10_8B, m60_8B]
m10_8B_comp70 <- BrayDist_8B.mat[m10_8B, m70_8B]
m10_8B_comp80 <- BrayDist_8B.mat[m10_8B, m80_8B]
m10_8B_comp90 <- BrayDist_8B.mat[m10_8B, m90_8B]
m10_8B_comp100 <- BrayDist_8B.mat[m10_8B, m100_8B]

# B-C distances between forest at 100m and each point along transect
m100_8B_comp100 <- BrayDist_8B.mat[m100_8B, m100_8B] #add 100 and 100
m100_8B_comp90 <- BrayDist_8B.mat[m100_8B, m90_8B]
m100_8B_comp80 <- BrayDist_8B.mat[m100_8B, m80_8B]
m100_8B_comp70 <- BrayDist_8B.mat[m100_8B, m70_8B]
m100_8B_comp60 <- BrayDist_8B.mat[m100_8B, m60_8B]
m100_8B_comp50 <- BrayDist_8B.mat[m100_8B, m50_8B]
m100_8B_comp40 <- BrayDist_8B.mat[m100_8B, m40_8B]
m100_8B_comp30 <- BrayDist_8B.mat[m100_8B, m30_8B]
m100_8B_comp20 <- BrayDist_8B.mat[m100_8B, m20_8B]
m100_8B_comp10 <- BrayDist_8B.mat[m100_8B, m10_8B]
###############

### L ###
EU_8_L <- EU_8_transSplit$L
# Get just ASV table
ASVtab_8L <- psotu2veg(EU_8_L)
# View(ASVtab_8L) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_8L <- vegdist(ASVtab_8L, method="bray")
BrayDist_8L.mat <- as.matrix(BrayDist_8L)
diag(BrayDist_8L.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_8L.mat) == rownames(sample_data(EU_8_L))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps8L <- colnames(BrayDist_8L.mat) #this is the Mapping Sample ID of these
meter8L <- sample_data(EU_8_L)$Meter
index8L.df <- data.frame(samps8L, meter8L)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_8L <- which(index8L.df$meter8L == 10)
m20_8L <- which(index8L.df$meter8L == 20)
m30_8L <- which(index8L.df$meter8L == 30)
m40_8L <- which(index8L.df$meter8L == 40)
m50_8L <- which(index8L.df$meter8L == 50)
m60_8L <- which(index8L.df$meter8L == 60)
m70_8L <- which(index8L.df$meter8L == 70)
m80_8L <- which(index8L.df$meter8L == 80)
m90_8L <- which(index8L.df$meter8L == 90)
m100_8L <- which(index8L.df$meter8L == 100)

# B-C distances between patch at 10m and each point along transect
m10_8L_comp10 <- BrayDist_8L.mat[m10_8L, m10_8L] #add 10 and 10
m10_8L_comp20 <- BrayDist_8L.mat[m10_8L, m20_8L]
m10_8L_comp30 <- BrayDist_8L.mat[m10_8L, m30_8L]
m10_8L_comp40 <- BrayDist_8L.mat[m10_8L, m40_8L]
m10_8L_comp50 <- BrayDist_8L.mat[m10_8L, m50_8L]
m10_8L_comp60 <- BrayDist_8L.mat[m10_8L, m60_8L]
m10_8L_comp70 <- BrayDist_8L.mat[m10_8L, m70_8L]
m10_8L_comp80 <- BrayDist_8L.mat[m10_8L, m80_8L]
m10_8L_comp90 <- BrayDist_8L.mat[m10_8L, m90_8L]
m10_8L_comp100 <- BrayDist_8L.mat[m10_8L, m100_8L]

# B-C distances between forest at 100m and each point along transect
m100_8L_comp100 <- BrayDist_8L.mat[m100_8L, m100_8L] #add 100 and 100
m100_8L_comp90 <- BrayDist_8L.mat[m100_8L, m90_8L]
m100_8L_comp80 <- BrayDist_8L.mat[m100_8L, m80_8L]
m100_8L_comp70 <- BrayDist_8L.mat[m100_8L, m70_8L]
m100_8L_comp60 <- BrayDist_8L.mat[m100_8L, m60_8L]
m100_8L_comp50 <- BrayDist_8L.mat[m100_8L, m50_8L]
m100_8L_comp40 <- BrayDist_8L.mat[m100_8L, m40_8L]
m100_8L_comp30 <- BrayDist_8L.mat[m100_8L, m30_8L]
m100_8L_comp20 <- BrayDist_8L.mat[m100_8L, m20_8L]
m100_8L_comp10 <- BrayDist_8L.mat[m100_8L, m10_8L]
###############

### R ###
EU_8_R <- EU_8_transSplit$R
# Get just ASV table
ASVtab_8R <- psotu2veg(EU_8_R)
# View(ASVtab_8R) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_8R <- vegdist(ASVtab_8R, method="bray")
BrayDist_8R.mat <- as.matrix(BrayDist_8R)
diag(BrayDist_8R.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_8R.mat) == rownames(sample_data(EU_8_R))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps8R <- colnames(BrayDist_8R.mat) #this is the Mapping Sample ID of these
meter8R <- sample_data(EU_8_R)$Meter
index8R.df <- data.frame(samps8R, meter8R)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_8R <- which(index8R.df$meter8R == 10)
m20_8R <- which(index8R.df$meter8R == 20)
m30_8R <- which(index8R.df$meter8R == 30)
m40_8R <- which(index8R.df$meter8R == 40)
m50_8R <- which(index8R.df$meter8R == 50)
m60_8R <- which(index8R.df$meter8R == 60)
m70_8R <- which(index8R.df$meter8R == 70)
m80_8R <- which(index8R.df$meter8R == 80)
m90_8R <- which(index8R.df$meter8R == 90)
m100_8R <- which(index8R.df$meter8R == 100)

# B-C distances between patch at 10m and each point along transect
m10_8R_comp10 <- BrayDist_8R.mat[m10_8R, m10_8R] #add 10 and 10
m10_8R_comp20 <- BrayDist_8R.mat[m10_8R, m20_8R]
m10_8R_comp30 <- BrayDist_8R.mat[m10_8R, m30_8R]
m10_8R_comp40 <- BrayDist_8R.mat[m10_8R, m40_8R]
m10_8R_comp50 <- BrayDist_8R.mat[m10_8R, m50_8R]
m10_8R_comp60 <- BrayDist_8R.mat[m10_8R, m60_8R]
m10_8R_comp70 <- BrayDist_8R.mat[m10_8R, m70_8R]
m10_8R_comp80 <- BrayDist_8R.mat[m10_8R, m80_8R]
m10_8R_comp90 <- BrayDist_8R.mat[m10_8R, m90_8R]
m10_8R_comp100 <- BrayDist_8R.mat[m10_8R, m100_8R]

# B-C distances between forest at 100m and each point along transect
m100_8R_comp100 <- BrayDist_8R.mat[m100_8R, m100_8R] #add 100 and 100
m100_8R_comp90 <- BrayDist_8R.mat[m100_8R, m90_8R]
m100_8R_comp80 <- BrayDist_8R.mat[m100_8R, m80_8R]
m100_8R_comp70 <- BrayDist_8R.mat[m100_8R, m70_8R]
m100_8R_comp60 <- BrayDist_8R.mat[m100_8R, m60_8R]
m100_8R_comp50 <- BrayDist_8R.mat[m100_8R, m50_8R]
m100_8R_comp40 <- BrayDist_8R.mat[m100_8R, m40_8R]
m100_8R_comp30 <- BrayDist_8R.mat[m100_8R, m30_8R]
m100_8R_comp20 <- BrayDist_8R.mat[m100_8R, m20_8R]
m100_8R_comp10 <- BrayDist_8R.mat[m100_8R, m10_8R]

# EU 53S
EU_53S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53S")
unique(sample_data(EU_53S_Soils.ps)$EU)
EU_53S_transSplit <- phyloseq_sep_variable(EU_53S_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_53S_T <- EU_53S_transSplit$T
# Get just ASV table
ASVtab_53ST <- psotu2veg(EU_53S_T)
# View(ASVtab_53ST) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53ST <- vegdist(ASVtab_53ST, method="bray")
BrayDist_53ST.mat <- as.matrix(BrayDist_53ST)
diag(BrayDist_53ST.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53ST.mat) == rownames(sample_data(EU_53S_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53ST <- colnames(BrayDist_53ST.mat) #this is the Mapping Sample ID of these
meter53ST <- sample_data(EU_53S_T)$Meter
index53ST.df <- data.frame(samps53ST, meter53ST)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53ST <- which(index53ST.df$meter53ST == 10)
m20_53ST <- which(index53ST.df$meter53ST == 20)
m30_53ST <- which(index53ST.df$meter53ST == 30)
m40_53ST <- which(index53ST.df$meter53ST == 40)
m50_53ST <- which(index53ST.df$meter53ST == 50)
m60_53ST <- which(index53ST.df$meter53ST == 60)
m70_53ST <- which(index53ST.df$meter53ST == 70)
m80_53ST <- which(index53ST.df$meter53ST == 80)
m90_53ST <- which(index53ST.df$meter53ST == 90)
m100_53ST <- which(index53ST.df$meter53ST == 100)

# B-C distances between patch at 10m and each point along transect
m10_53ST_comp10 <- BrayDist_53ST.mat[m10_53ST, m10_53ST] #add 10 and 10
m10_53ST_comp20 <- BrayDist_53ST.mat[m10_53ST, m20_53ST]
m10_53ST_comp30 <- BrayDist_53ST.mat[m10_53ST, m30_53ST]
m10_53ST_comp40 <- BrayDist_53ST.mat[m10_53ST, m40_53ST]
m10_53ST_comp50 <- BrayDist_53ST.mat[m10_53ST, m50_53ST]
m10_53ST_comp60 <- BrayDist_53ST.mat[m10_53ST, m60_53ST]
m10_53ST_comp70 <- BrayDist_53ST.mat[m10_53ST, m70_53ST]
m10_53ST_comp80 <- BrayDist_53ST.mat[m10_53ST, m80_53ST]
m10_53ST_comp90 <- BrayDist_53ST.mat[m10_53ST, m90_53ST]
m10_53ST_comp100 <- BrayDist_53ST.mat[m10_53ST, m100_53ST]

# B-C distances between forest at 100m and each point along transect
m100_53ST_comp100 <- BrayDist_53ST.mat[m100_53ST, m100_53ST] #add 100 and 100
m100_53ST_comp90 <- BrayDist_53ST.mat[m100_53ST, m90_53ST]
m100_53ST_comp80 <- BrayDist_53ST.mat[m100_53ST, m80_53ST]
m100_53ST_comp70 <- BrayDist_53ST.mat[m100_53ST, m70_53ST]
m100_53ST_comp60 <- BrayDist_53ST.mat[m100_53ST, m60_53ST]
m100_53ST_comp50 <- BrayDist_53ST.mat[m100_53ST, m50_53ST]
m100_53ST_comp40 <- BrayDist_53ST.mat[m100_53ST, m40_53ST]
m100_53ST_comp30 <- BrayDist_53ST.mat[m100_53ST, m30_53ST]
m100_53ST_comp20 <- BrayDist_53ST.mat[m100_53ST, m20_53ST]
m100_53ST_comp10 <- BrayDist_53ST.mat[m100_53ST, m10_53ST]


###############
### B ###
EU_53S_B <- EU_53S_transSplit$B
# Get just ASV table
ASVtab_53SB <- psotu2veg(EU_53S_B)
# View(ASVtab_53SB) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53SB <- vegdist(ASVtab_53SB, method="bray")
BrayDist_53SB.mat <- as.matrix(BrayDist_53SB)
diag(BrayDist_53SB.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53SB.mat) == rownames(sample_data(EU_53S_B))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53SB <- colnames(BrayDist_53SB.mat) #this is the Mapping Sample ID of these
meter53SB <- sample_data(EU_53S_B)$Meter
index53SB.df <- data.frame(samps53SB, meter53SB)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53SB <- which(index53SB.df$meter53SB == 10)
m20_53SB <- which(index53SB.df$meter53SB == 20)
m30_53SB <- which(index53SB.df$meter53SB == 30)
m40_53SB <- which(index53SB.df$meter53SB == 40)
m50_53SB <- which(index53SB.df$meter53SB == 50)
m60_53SB <- which(index53SB.df$meter53SB == 60)
m70_53SB <- which(index53SB.df$meter53SB == 70)
m80_53SB <- which(index53SB.df$meter53SB == 80)
m90_53SB <- which(index53SB.df$meter53SB == 90)
m100_53SB <- which(index53SB.df$meter53SB == 100)

# B-C distances between patch at 10m and each point along transect
m10_53SB_comp10 <- BrayDist_53SB.mat[m10_53SB, m10_53SB] #add 10 and 10
m10_53SB_comp20 <- BrayDist_53SB.mat[m10_53SB, m20_53SB]
m10_53SB_comp30 <- BrayDist_53SB.mat[m10_53SB, m30_53SB]
m10_53SB_comp40 <- BrayDist_53SB.mat[m10_53SB, m40_53SB]
m10_53SB_comp50 <- BrayDist_53SB.mat[m10_53SB, m50_53SB]
m10_53SB_comp60 <- BrayDist_53SB.mat[m10_53SB, m60_53SB]
m10_53SB_comp70 <- BrayDist_53SB.mat[m10_53SB, m70_53SB]
m10_53SB_comp80 <- BrayDist_53SB.mat[m10_53SB, m80_53SB]
m10_53SB_comp90 <- BrayDist_53SB.mat[m10_53SB, m90_53SB]
m10_53SB_comp100 <- BrayDist_53SB.mat[m10_53SB, m100_53SB]

# B-C distances between forest at 100m and each point along transect
m100_53SB_comp100 <- BrayDist_53SB.mat[m100_53SB, m100_53SB] #add 100 and 100
m100_53SB_comp90 <- BrayDist_53SB.mat[m100_53SB, m90_53SB]
m100_53SB_comp80 <- BrayDist_53SB.mat[m100_53SB, m80_53SB]
m100_53SB_comp70 <- BrayDist_53SB.mat[m100_53SB, m70_53SB]
m100_53SB_comp60 <- BrayDist_53SB.mat[m100_53SB, m60_53SB]
m100_53SB_comp50 <- BrayDist_53SB.mat[m100_53SB, m50_53SB]
m100_53SB_comp40 <- BrayDist_53SB.mat[m100_53SB, m40_53SB]
m100_53SB_comp30 <- BrayDist_53SB.mat[m100_53SB, m30_53SB]
m100_53SB_comp20 <- BrayDist_53SB.mat[m100_53SB, m20_53SB]
m100_53SB_comp10 <- BrayDist_53SB.mat[m100_53SB, m10_53SB]
###############

### L ###
EU_53S_L <- EU_53S_transSplit$L
# Get just ASV table
ASVtab_53SL <- psotu2veg(EU_53S_L)
# View(ASVtab_53SL) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53SL <- vegdist(ASVtab_53SL, method="bray")
BrayDist_53SL.mat <- as.matrix(BrayDist_53SL)
diag(BrayDist_53SL.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53SL.mat) == rownames(sample_data(EU_53S_L))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53SL <- colnames(BrayDist_53SL.mat) #this is the Mapping Sample ID of these
meter53SL <- sample_data(EU_53S_L)$Meter
index53SL.df <- data.frame(samps53SL, meter53SL)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53SL <- which(index53SL.df$meter53SL == 10)
m20_53SL <- which(index53SL.df$meter53SL == 20)
m30_53SL <- which(index53SL.df$meter53SL == 30)
m40_53SL <- which(index53SL.df$meter53SL == 40)
m50_53SL <- which(index53SL.df$meter53SL == 50)
m60_53SL <- which(index53SL.df$meter53SL == 60)
m70_53SL <- which(index53SL.df$meter53SL == 70)
m80_53SL <- which(index53SL.df$meter53SL == 80)
m90_53SL <- which(index53SL.df$meter53SL == 90)
m100_53SL <- which(index53SL.df$meter53SL == 100)

# B-C distances between patch at 10m and each point along transect
m10_53SL_comp10 <- BrayDist_53SL.mat[m10_53SL, m10_53SL] #add 10 and 10
m10_53SL_comp20 <- BrayDist_53SL.mat[m10_53SL, m20_53SL]
m10_53SL_comp30 <- BrayDist_53SL.mat[m10_53SL, m30_53SL]
m10_53SL_comp40 <- BrayDist_53SL.mat[m10_53SL, m40_53SL]
m10_53SL_comp50 <- BrayDist_53SL.mat[m10_53SL, m50_53SL]
m10_53SL_comp60 <- BrayDist_53SL.mat[m10_53SL, m60_53SL]
m10_53SL_comp70 <- BrayDist_53SL.mat[m10_53SL, m70_53SL]
m10_53SL_comp80 <- BrayDist_53SL.mat[m10_53SL, m80_53SL]
m10_53SL_comp90 <- BrayDist_53SL.mat[m10_53SL, m90_53SL]
m10_53SL_comp100 <- BrayDist_53SL.mat[m10_53SL, m100_53SL]

# B-C distances between forest at 100m and each point along transect
m100_53SL_comp100 <- BrayDist_53SL.mat[m100_53SL, m100_53SL] #add 100 and 100
m100_53SL_comp90 <- BrayDist_53SL.mat[m100_53SL, m90_53SL]
m100_53SL_comp80 <- BrayDist_53SL.mat[m100_53SL, m80_53SL]
m100_53SL_comp70 <- BrayDist_53SL.mat[m100_53SL, m70_53SL]
m100_53SL_comp60 <- BrayDist_53SL.mat[m100_53SL, m60_53SL]
m100_53SL_comp50 <- BrayDist_53SL.mat[m100_53SL, m50_53SL]
m100_53SL_comp40 <- BrayDist_53SL.mat[m100_53SL, m40_53SL]
m100_53SL_comp30 <- BrayDist_53SL.mat[m100_53SL, m30_53SL]
m100_53SL_comp20 <- BrayDist_53SL.mat[m100_53SL, m20_53SL]
m100_53SL_comp10 <- BrayDist_53SL.mat[m100_53SL, m10_53SL]
###############

### R ###
EU_53S_R <- EU_53S_transSplit$R
# Get just ASV table
ASVtab_53SR <- psotu2veg(EU_53S_R)
# View(ASVtab_53SR) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_53SR <- vegdist(ASVtab_53SR, method="bray")
BrayDist_53SR.mat <- as.matrix(BrayDist_53SR)
diag(BrayDist_53SR.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_53SR.mat) == rownames(sample_data(EU_53S_R))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps53SR <- colnames(BrayDist_53SR.mat) #this is the Mapping Sample ID of these
meter53SR <- sample_data(EU_53S_R)$Meter
index53SR.df <- data.frame(samps53SR, meter53SR)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_53SR <- which(index53SR.df$meter53SR == 10)
m20_53SR <- which(index53SR.df$meter53SR == 20)
m30_53SR <- which(index53SR.df$meter53SR == 30)
m40_53SR <- which(index53SR.df$meter53SR == 40)
m50_53SR <- which(index53SR.df$meter53SR == 50)
m60_53SR <- which(index53SR.df$meter53SR == 60)
m70_53SR <- which(index53SR.df$meter53SR == 70)
m80_53SR <- which(index53SR.df$meter53SR == 80)
m90_53SR <- which(index53SR.df$meter53SR == 90)
m100_53SR <- which(index53SR.df$meter53SR == 100)

# B-C distances between patch at 10m and each point along transect
m10_53SR_comp10 <- BrayDist_53SR.mat[m10_53SR, m10_53SR] #add 10 and 10
m10_53SR_comp20 <- BrayDist_53SR.mat[m10_53SR, m20_53SR]
m10_53SR_comp30 <- BrayDist_53SR.mat[m10_53SR, m30_53SR]
m10_53SR_comp40 <- BrayDist_53SR.mat[m10_53SR, m40_53SR]
m10_53SR_comp50 <- BrayDist_53SR.mat[m10_53SR, m50_53SR]
m10_53SR_comp60 <- BrayDist_53SR.mat[m10_53SR, m60_53SR]
m10_53SR_comp70 <- BrayDist_53SR.mat[m10_53SR, m70_53SR]
m10_53SR_comp80 <- BrayDist_53SR.mat[m10_53SR, m80_53SR]
m10_53SR_comp90 <- BrayDist_53SR.mat[m10_53SR, m90_53SR]
m10_53SR_comp100 <- BrayDist_53SR.mat[m10_53SR, m100_53SR]

# B-C distances between forest at 100m and each point along transect
m100_53SR_comp100 <- BrayDist_53SR.mat[m100_53SR, m100_53SR] #add 100 and 100
m100_53SR_comp90 <- BrayDist_53SR.mat[m100_53SR, m90_53SR]
m100_53SR_comp80 <- BrayDist_53SR.mat[m100_53SR, m80_53SR]
m100_53SR_comp70 <- BrayDist_53SR.mat[m100_53SR, m70_53SR]
m100_53SR_comp60 <- BrayDist_53SR.mat[m100_53SR, m60_53SR]
m100_53SR_comp50 <- BrayDist_53SR.mat[m100_53SR, m50_53SR]
m100_53SR_comp40 <- BrayDist_53SR.mat[m100_53SR, m40_53SR]
m100_53SR_comp30 <- BrayDist_53SR.mat[m100_53SR, m30_53SR]
m100_53SR_comp20 <- BrayDist_53SR.mat[m100_53SR, m20_53SR]
m100_53SR_comp10 <- BrayDist_53SR.mat[m100_53SR, m10_53SR]
###############

# EU 10
EU_10_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_10")
unique(sample_data(EU_10_Soils.ps)$EU)
EU_10_transSplit <- phyloseq_sep_variable(EU_10_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_10_T <- EU_10_transSplit$T
# Get just ASV table
ASVtab_10T <- psotu2veg(EU_10_T)
# View(ASVtab_10T) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_10T <- vegdist(ASVtab_10T, method="bray")
BrayDist_10T.mat <- as.matrix(BrayDist_10T)
diag(BrayDist_10T.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_10T.mat) == rownames(sample_data(EU_10_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps10T <- colnames(BrayDist_10T.mat) #this is the Mapping Sample ID of these
meter10T <- sample_data(EU_10_T)$Meter
index10T.df <- data.frame(samps10T, meter10T)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_10T <- which(index10T.df$meter10T == 10)
m20_10T <- which(index10T.df$meter10T == 20)
m30_10T <- which(index10T.df$meter10T == 30)
m40_10T <- which(index10T.df$meter10T == 40)
m50_10T <- which(index10T.df$meter10T == 50)
m60_10T <- which(index10T.df$meter10T == 60)
m70_10T <- which(index10T.df$meter10T == 70)
m80_10T <- which(index10T.df$meter10T == 80)
m90_10T <- which(index10T.df$meter10T == 90)
m100_10T <- which(index10T.df$meter10T == 100)

# B-C distances between patch at 10m and each point along transect
m10_10T_comp10 <- BrayDist_10T.mat[m10_10T, m10_10T] #add 10 and 10
m10_10T_comp20 <- BrayDist_10T.mat[m10_10T, m20_10T]
m10_10T_comp30 <- BrayDist_10T.mat[m10_10T, m30_10T]
m10_10T_comp40 <- BrayDist_10T.mat[m10_10T, m40_10T]
m10_10T_comp50 <- BrayDist_10T.mat[m10_10T, m50_10T]
m10_10T_comp60 <- BrayDist_10T.mat[m10_10T, m60_10T]
m10_10T_comp70 <- BrayDist_10T.mat[m10_10T, m70_10T]
m10_10T_comp80 <- BrayDist_10T.mat[m10_10T, m80_10T]
m10_10T_comp90 <- BrayDist_10T.mat[m10_10T, m90_10T]
m10_10T_comp100 <- BrayDist_10T.mat[m10_10T, m100_10T]

# B-C distances between forest at 100m and each point along transect
m100_10T_comp100 <- BrayDist_10T.mat[m100_10T, m100_10T] #add 100 and 100
m100_10T_comp90 <- BrayDist_10T.mat[m100_10T, m90_10T]
m100_10T_comp80 <- BrayDist_10T.mat[m100_10T, m80_10T]
m100_10T_comp70 <- BrayDist_10T.mat[m100_10T, m70_10T]
m100_10T_comp60 <- BrayDist_10T.mat[m100_10T, m60_10T]
m100_10T_comp50 <- BrayDist_10T.mat[m100_10T, m50_10T]
m100_10T_comp40 <- BrayDist_10T.mat[m100_10T, m40_10T]
m100_10T_comp30 <- BrayDist_10T.mat[m100_10T, m30_10T]
m100_10T_comp20 <- BrayDist_10T.mat[m100_10T, m20_10T]
m100_10T_comp10 <- BrayDist_10T.mat[m100_10T, m10_10T]


###############
### B ###
EU_10_B <- EU_10_transSplit$B
# Get just ASV table
ASVtab_10B <- psotu2veg(EU_10_B)
# View(ASVtab_10B) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_10B <- vegdist(ASVtab_10B, method="bray")
BrayDist_10B.mat <- as.matrix(BrayDist_10B)
diag(BrayDist_10B.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_10B.mat) == rownames(sample_data(EU_10_B))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps10B <- colnames(BrayDist_10B.mat) #this is the Mapping Sample ID of these
meter10B <- sample_data(EU_10_B)$Meter
index10B.df <- data.frame(samps10B, meter10B)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_10B <- which(index10B.df$meter10B == 10)
m20_10B <- which(index10B.df$meter10B == 20)
m30_10B <- which(index10B.df$meter10B == 30)
m40_10B <- which(index10B.df$meter10B == 40)
m50_10B <- which(index10B.df$meter10B == 50)
m60_10B <- which(index10B.df$meter10B == 60)
m70_10B <- which(index10B.df$meter10B == 70)
m80_10B <- which(index10B.df$meter10B == 80)
m90_10B <- which(index10B.df$meter10B == 90)
m100_10B <- which(index10B.df$meter10B == 100)

# B-C distances between patch at 10m and each point along transect
m10_10B_comp10 <- BrayDist_10B.mat[m10_10B, m10_10B] #add 10 and 10
m10_10B_comp20 <- BrayDist_10B.mat[m10_10B, m20_10B]
m10_10B_comp30 <- BrayDist_10B.mat[m10_10B, m30_10B]
m10_10B_comp40 <- BrayDist_10B.mat[m10_10B, m40_10B]
m10_10B_comp50 <- BrayDist_10B.mat[m10_10B, m50_10B]
m10_10B_comp60 <- BrayDist_10B.mat[m10_10B, m60_10B]
m10_10B_comp70 <- BrayDist_10B.mat[m10_10B, m70_10B]
m10_10B_comp80 <- BrayDist_10B.mat[m10_10B, m80_10B]
m10_10B_comp90 <- BrayDist_10B.mat[m10_10B, m90_10B]
m10_10B_comp100 <- BrayDist_10B.mat[m10_10B, m100_10B]

# B-C distances between forest at 100m and each point along transect
m100_10B_comp100 <- BrayDist_10B.mat[m100_10B, m100_10B] #add 100 and 100
m100_10B_comp90 <- BrayDist_10B.mat[m100_10B, m90_10B]
m100_10B_comp80 <- BrayDist_10B.mat[m100_10B, m80_10B]
m100_10B_comp70 <- BrayDist_10B.mat[m100_10B, m70_10B]
m100_10B_comp60 <- BrayDist_10B.mat[m100_10B, m60_10B]
m100_10B_comp50 <- BrayDist_10B.mat[m100_10B, m50_10B]
m100_10B_comp40 <- BrayDist_10B.mat[m100_10B, m40_10B]
m100_10B_comp30 <- BrayDist_10B.mat[m100_10B, m30_10B]
m100_10B_comp20 <- BrayDist_10B.mat[m100_10B, m20_10B]
m100_10B_comp10 <- BrayDist_10B.mat[m100_10B, m10_10B]
###############

### L ###
EU_10_L <- EU_10_transSplit$L
# Get just ASV table
ASVtab_10L <- psotu2veg(EU_10_L)
# View(ASVtab_10L) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_10L <- vegdist(ASVtab_10L, method="bray")
BrayDist_10L.mat <- as.matrix(BrayDist_10L)
diag(BrayDist_10L.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_10L.mat) == rownames(sample_data(EU_10_L))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps10L <- colnames(BrayDist_10L.mat) #this is the Mapping Sample ID of these
meter10L <- sample_data(EU_10_L)$Meter
index10L.df <- data.frame(samps10L, meter10L)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_10L <- which(index10L.df$meter10L == 10)
m20_10L <- which(index10L.df$meter10L == 20)
m30_10L <- which(index10L.df$meter10L == 30)
m40_10L <- which(index10L.df$meter10L == 40)
m50_10L <- which(index10L.df$meter10L == 50)
m60_10L <- which(index10L.df$meter10L == 60)
m70_10L <- which(index10L.df$meter10L == 70)
m80_10L <- which(index10L.df$meter10L == 80)
m90_10L <- which(index10L.df$meter10L == 90)
m100_10L <- which(index10L.df$meter10L == 100)

# B-C distances between patch at 10m and each point along transect
m10_10L_comp10 <- BrayDist_10L.mat[m10_10L, m10_10L] #add 10 and 10
m10_10L_comp20 <- BrayDist_10L.mat[m10_10L, m20_10L]
m10_10L_comp30 <- BrayDist_10L.mat[m10_10L, m30_10L]
m10_10L_comp40 <- BrayDist_10L.mat[m10_10L, m40_10L]
m10_10L_comp50 <- BrayDist_10L.mat[m10_10L, m50_10L]
m10_10L_comp60 <- BrayDist_10L.mat[m10_10L, m60_10L]
m10_10L_comp70 <- BrayDist_10L.mat[m10_10L, m70_10L]
m10_10L_comp80 <- BrayDist_10L.mat[m10_10L, m80_10L]
m10_10L_comp90 <- BrayDist_10L.mat[m10_10L, m90_10L]
m10_10L_comp100 <- BrayDist_10L.mat[m10_10L, m100_10L]

# B-C distances between forest at 100m and each point along transect
m100_10L_comp100 <- BrayDist_10L.mat[m100_10L, m100_10L] #add 100 and 100
m100_10L_comp90 <- BrayDist_10L.mat[m100_10L, m90_10L]
m100_10L_comp80 <- BrayDist_10L.mat[m100_10L, m80_10L]
m100_10L_comp70 <- BrayDist_10L.mat[m100_10L, m70_10L]
m100_10L_comp60 <- BrayDist_10L.mat[m100_10L, m60_10L]
m100_10L_comp50 <- BrayDist_10L.mat[m100_10L, m50_10L]
m100_10L_comp40 <- BrayDist_10L.mat[m100_10L, m40_10L]
m100_10L_comp30 <- BrayDist_10L.mat[m100_10L, m30_10L]
m100_10L_comp20 <- BrayDist_10L.mat[m100_10L, m20_10L]
m100_10L_comp10 <- BrayDist_10L.mat[m100_10L, m10_10L]
###############

### R ###
EU_10_R <- EU_10_transSplit$R
# Get just ASV table
ASVtab_10R <- psotu2veg(EU_10_R)
# View(ASVtab_10R) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_10R <- vegdist(ASVtab_10R, method="bray")
BrayDist_10R.mat <- as.matrix(BrayDist_10R)
diag(BrayDist_10R.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_10R.mat) == rownames(sample_data(EU_10_R))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps10R <- colnames(BrayDist_10R.mat) #this is the Mapping Sample ID of these
meter10R <- sample_data(EU_10_R)$Meter
index10R.df <- data.frame(samps10R, meter10R)

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_10R <- which(index10R.df$meter10R == 10)
m20_10R <- which(index10R.df$meter10R == 20)
m30_10R <- which(index10R.df$meter10R == 30)
m40_10R <- which(index10R.df$meter10R == 40)
m50_10R <- which(index10R.df$meter10R == 50)
m60_10R <- which(index10R.df$meter10R == 60)
m70_10R <- which(index10R.df$meter10R == 70)
m80_10R <- which(index10R.df$meter10R == 80)
m90_10R <- which(index10R.df$meter10R == 90)
m100_10R <- which(index10R.df$meter10R == 100)

# B-C distances between patch at 10m and each point along transect
m10_10R_comp10 <- BrayDist_10R.mat[m10_10R, m10_10R] #add 10 and 10
m10_10R_comp20 <- BrayDist_10R.mat[m10_10R, m20_10R]
m10_10R_comp30 <- BrayDist_10R.mat[m10_10R, m30_10R]
m10_10R_comp40 <- BrayDist_10R.mat[m10_10R, m40_10R]
m10_10R_comp50 <- BrayDist_10R.mat[m10_10R, m50_10R]
m10_10R_comp60 <- BrayDist_10R.mat[m10_10R, m60_10R]
m10_10R_comp70 <- BrayDist_10R.mat[m10_10R, m70_10R]
m10_10R_comp80 <- BrayDist_10R.mat[m10_10R, m80_10R]
m10_10R_comp90 <- BrayDist_10R.mat[m10_10R, m90_10R]
m10_10R_comp100 <- BrayDist_10R.mat[m10_10R, m100_10R]

# B-C distances between forest at 100m and each point along transect
m100_10R_comp100 <- BrayDist_10R.mat[m100_10R, m100_10R] #add 100 and 100
m100_10R_comp90 <- BrayDist_10R.mat[m100_10R, m90_10R]
m100_10R_comp80 <- BrayDist_10R.mat[m100_10R, m80_10R]
m100_10R_comp70 <- BrayDist_10R.mat[m100_10R, m70_10R]
m100_10R_comp60 <- BrayDist_10R.mat[m100_10R, m60_10R]
m100_10R_comp50 <- BrayDist_10R.mat[m100_10R, m50_10R]
m100_10R_comp40 <- BrayDist_10R.mat[m100_10R, m40_10R]
m100_10R_comp30 <- BrayDist_10R.mat[m100_10R, m30_10R]
m100_10R_comp20 <- BrayDist_10R.mat[m100_10R, m20_10R]
m100_10R_comp10 <- BrayDist_10R.mat[m100_10R, m10_10R]
###############

# Plot line graph of EU 10:
xMeters <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
# Dissimilarity from patch:
y_10T_comp10m <- c(m10_10T_comp10, m10_10T_comp20, m10_10T_comp30, m10_10T_comp40,
                   m10_10T_comp50, m10_10T_comp60, m10_10T_comp70, m10_10T_comp80,
                   m10_10T_comp90, m10_10T_comp100)
patch_10T <- cbind(xMeters, y_10T_comp10m)
y_10B_comp10m <- c(m10_10B_comp10, m10_10B_comp20, m10_10B_comp30, m10_10B_comp40,
                   m10_10B_comp50, m10_10B_comp60, m10_10B_comp70, m10_10B_comp80,
                   m10_10B_comp90, m10_10B_comp100)
patch_10B <- cbind(xMeters, y_10B_comp10m)
y_10L_comp10m <- c(m10_10L_comp10, m10_10L_comp20, m10_10L_comp30, m10_10L_comp40,
                   NA, m10_10L_comp60, m10_10L_comp70, m10_10L_comp80,
                   m10_10L_comp90, m10_10L_comp100)
patch_10L <- cbind(xMeters, y_10L_comp10m)

y_10R_comp10m <- c(m10_10R_comp10, m10_10R_comp20, m10_10R_comp30, m10_10R_comp40,
                   m10_10R_comp50, m10_10R_comp60, m10_10R_comp70, m10_10R_comp80,
                   m10_10R_comp90, m10_10R_comp100)
patch_10R <- cbind(xMeters, y_10R_comp10m)
## Dissimilarity from forest:
y_10T_comp100m <- c(m100_10T_comp10, m100_10T_comp20, m100_10T_comp30, m100_10T_comp40,
                    m100_10T_comp50, m100_10T_comp60, m100_10T_comp70, m100_10T_comp80,
                    m100_10T_comp90, m100_10T_comp100)
forest_10T <- cbind(xMeters, y_10T_comp100m)
y_10B_comp100m <- c(m100_10B_comp10, m100_10B_comp20, m100_10B_comp30, m100_10B_comp40,
                    m100_10B_comp50, m100_10B_comp60, m100_10B_comp70, m100_10B_comp80,
                    m100_10B_comp90, m100_10B_comp100)
forest_10B <- cbind(xMeters, y_10B_comp100m)
y_10L_comp100m <- c(m100_10L_comp10, m100_10L_comp20, m100_10L_comp30, m100_10L_comp40,
                    NA, m100_10L_comp60, m100_10L_comp70, m100_10L_comp80,
                    m100_10L_comp90, m100_10L_comp100)
forest_10L <- cbind(xMeters, y_10L_comp100m)

y_10R_comp100m <- c(m100_10R_comp10, m100_10R_comp20, m100_10R_comp30, m100_10R_comp40,
                    m100_10R_comp50, m100_10R_comp60, m100_10R_comp70, m100_10R_comp80,
                    m100_10R_comp90, m100_10R_comp100)
forest_10R <- cbind(xMeters, y_10R_comp100m)

PatchForest_10Trans_plot <- matplot(forest_10T[,1], forest_10T[,2], type = "l", xlab= "meter", col= "darkgreen",
                           ylab = "Bray-Curtis dissimilarity", main= "EU 10: Dissimilarity Relative to 10 m & 100 m",
                           ylim=c(0.2, 1.0), xlim = c(10, 100))
      matlines(forest_10B[,1], forest_10B[,2], col= "darkgreen")
      matlines(forest_10L[,1], forest_10L[,2], col= "darkgreen")
      matlines(forest_10R[,1], forest_10R[,2], col= "darkgreen")
      matlines(patch_10T[,1], patch_10T[,2], col= "goldenrod")
      matlines(patch_10B[,1], patch_10B[,2], col= "goldenrod")
      matlines(patch_10L[,1], patch_10L[,2], col= "goldenrod")
      matlines(patch_10R[,1], patch_10R[,2], col= "goldenrod")


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
                                 "60", "70", "80", "90", "100"), col= "goldenrod",
                       cex.axis = 0.8,
                       cex.lab = 1,
                       ylim=c(0.0, 1.0))
mtext(text=bold_eu10p, side=3, adj = -0.065, line = 2)
forest_box10 <- boxplot(list(m100_10_comp10, m100_10_comp20, m100_10_comp30, m100_10_comp40,
                             m100_10_comp50, m100_10_comp60, m100_10_comp70,
                             m100_10_comp80, m100_10_comp90, m100_10_comp100), 
                        ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                        names = c("10", "20", "30", "40", "50",
                                  "60", "70", "80", "90", "100"), col = "darkgreen", 
                        cex.axis = 0.8,
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
# (using distances created in part 1 (comparisons across all transects))

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

###########
# EU 10
##########
# (using distances created in part 1 (comparisons across all transects))

# B-C distances between patch at 50m (edge) and each point along transect
m50_10_comp10 <- BrayDist_10.mat[m50_10, m10_10] #edge and 10m
m50_10_comp20 <- BrayDist_10.mat[m50_10, m20_10] #edge and 20m
m50_10_comp30 <- BrayDist_10.mat[m50_10, m30_10] #edge and 30m
m50_10_comp40 <- BrayDist_10.mat[m50_10, m40_10] #edge and 40m
m50_10_comp50 <- BrayDist_10.mat[m50_10, m50_10] #edge and self
m50_10_comp60 <- BrayDist_10.mat[m50_10, m60_10] #edge and 60m
m50_10_comp70 <- BrayDist_10.mat[m50_10, m70_10] #edge and 70m
m50_10_comp80 <- BrayDist_10.mat[m50_10, m80_10] #edge and 80m
m50_10_comp90 <- BrayDist_10.mat[m50_10, m90_10] #edge and 90m
m50_10_comp100 <- BrayDist_10.mat[m50_10, m100_10] #edge and 100m

bold_eu10edge <- expression(bold("EU 10: Dissimilarity from edge (50 m)"))
quartz()
edgecomp_box10 <- boxplot(list(m50_10_comp10, m50_10_comp20, m50_10_comp30, m50_10_comp40,
                               m50_10_comp50, m50_10_comp60, m50_10_comp70,
                               m50_10_comp80, m50_10_comp90, m50_10_comp100),
                          ylab = "Bray-Curtis Dissimilarity", xlab = "Transect Meter",
                          names = c("10", "20", "30", "40", "Edge (50 m)",
                                    "60", "70", "80", "90", "100"), cex.axis = 0.8,
                          cex.lab = 1,
                          ylim=c(0.0, 1.0))
mtext(text=bold_eu10edge, side=3, adj = -0.065, line = 2)

#########################################################################
# 4) The pooled forest samples versus all other points along the transect
# and 5) The pooled patch samples versus all other points along the transect
#########################################################################

##########
# EU 52
##########

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

##########
# EU 10
##########

######
EU_10_Soils.ps
sample_data(EU_10_Soils.ps)
pooledHabitat_10 <- merge_samples(EU_10_Soils.ps, "Habitat") #get samples pooled by habitat
pooledHabitat_10
pooledHabitat_10ASV <- psotu2veg(pooledHabitat_10)
rownames(pooledHabitat_10ASV)
pooledEdge_10ASV <- pooledHabitat_10ASV[1,] #pull out just pooled edge samples
pooledForest_10ASV <- pooledHabitat_10ASV[2,] #pull out pooled forest samples
pooledPatch_10ASV <- pooledHabitat_10ASV[3,] #pull out pooled patch samples

# Now to get dissimilarity matrices, I'll have to merge these pooled samples and with the
# non pooled samples of the other type

# Pooled forest versus all other samples along transect 
length(pooledForest_10ASV) #10,257 same as in whole EU_10_Soils.ps, even if some ASVs are zero
names(pooledForest_10ASV) #gives ASV names
summary(pooledForest_10ASV)
ASVtab_10 #from up above, this is ASV table of EU 10
ASVtab_10t <- t(ASVtab_10)
dim(ASVtab_10t)
#View(ASVtab_10t)
unique(names(pooledForest_10ASV) == rownames(ASVtab_10t)) 
# Because they match (line above), we can cbind to get new ASV table
pooledFvPtransect_10 <- t(cbind.data.frame(pooledForest_10ASV, ASVtab_10t))
#View(pooledFvPtransect_10) #needs to be samples as rows, ASVs as columns!

# Get Bray-Curtis dissimilarity matrix
pooledFvPtransect_10_Bray <- vegdist(pooledFvPtransect_10, method="bray")
pooledFvPtransect_10_Bray <- as.matrix(pooledFvPtransect_10_Bray)
diag(pooledFvPtransect_10_Bray) <- NA
#View(pooledFvPtransect_10_Bray) # this has pooled forest samples versus all other samples (i.e. patch and edge)

# Need to convert samples to meter?

# Pull out B-Curtis dissimilarities between each sample and pooled forest

