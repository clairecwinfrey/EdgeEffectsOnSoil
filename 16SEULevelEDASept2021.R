## Site level 16S analysis
# (started) September 7, 2021

########################

# Load in objects made in 16SExploratoryDataAnalysisAugust2021

# Prune samples to separate by EU and check to make sure eaach looks right

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

# Define function to get ASV table out of phyloseq format
# from https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# For each site, get ASV table, Bray-Curtis dissimilarities

###########
# EU 52
##########

# Get ASV table
ASVtab_52 <- psotu2veg(EU_52_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_52 <- vegdist(ASVtab_52, method="bray")
BrayDist_52.mat <- as.matrix(BrayDist_52)
diag(BrayDist_52.mat) <- NA
BrayDist_52.mat[lower.tri(BrayDist_52.mat)] <- NA #because symmetrical, can make lower triangle zero

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
m100_52_comp90 <- BrayDist_52.mat[m100_52, m90_52]
m100_52_comp80 <- BrayDist_52.mat[m100_52, m80_52]
m100_52_comp70 <- BrayDist_52.mat[m100_52, m70_52]
m100_52_comp60 <- BrayDist_52.mat[m100_52, m60_52]
m100_52_comp50 <- BrayDist_52.mat[m100_52, m50_52]
m100_52_comp40 <- BrayDist_52.mat[m100_52, m40_52]
m100_52_comp30 <- BrayDist_52.mat[m100_52, m30_52]
m100_52_comp20 <- BrayDist_52.mat[m100_52, m20_52]
m100_52_comp10 <- BrayDist_52.mat[m100_52, m10_52]

###########
# EU 53N
##########

ASVtab_53N <- psotu2veg(EU_53N_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_53N <- vegdist(ASVtab_53N, method="bray")
BrayDist_53N.mat <- as.matrix(BrayDist_53N)
diag(BrayDist_53N.mat) <- NA
BrayDist_53N.mat[lower.tri(BrayDist_53N.mat)] <- NA #because symmetrical, can make lower triangle zero

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
m100_53N_comp90 <- BrayDist_53N.mat[m100_53N, m90_53N]
m100_53N_comp80 <- BrayDist_53N.mat[m100_53N, m80_53N]
m100_53N_comp70 <- BrayDist_53N.mat[m100_53N, m70_53N]
m100_53N_comp60 <- BrayDist_53N.mat[m100_53N, m60_53N]
m100_53N_comp50 <- BrayDist_53N.mat[m100_53N, m50_53N]
m100_53N_comp40 <- BrayDist_53N.mat[m100_53N, m40_53N]
m100_53N_comp30 <- BrayDist_53N.mat[m100_53N, m30_53N]
m100_53N_comp20 <- BrayDist_53N.mat[m100_53N, m20_53N]
m100_53N_comp10 <- BrayDist_53N.mat[m100_53N, m10_53N]

###########
# EU 54S
##########

ASVtab_54S <- psotu2veg(EU_54S_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_53N <- vegdist(ASVtab_53N, method="bray")
BrayDist_53N.mat <- as.matrix(BrayDist_53N)
diag(BrayDist_53N.mat) <- NA
BrayDist_53N.mat[lower.tri(BrayDist_53N.mat)] <- NA #because symmetrical, can make lower triangle zero

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
m100_53N_comp90 <- BrayDist_53N.mat[m100_53N, m90_53N]
m100_53N_comp80 <- BrayDist_53N.mat[m100_53N, m80_53N]
m100_53N_comp70 <- BrayDist_53N.mat[m100_53N, m70_53N]
m100_53N_comp60 <- BrayDist_53N.mat[m100_53N, m60_53N]
m100_53N_comp50 <- BrayDist_53N.mat[m100_53N, m50_53N]
m100_53N_comp40 <- BrayDist_53N.mat[m100_53N, m40_53N]
m100_53N_comp30 <- BrayDist_53N.mat[m100_53N, m30_53N]
m100_53N_comp20 <- BrayDist_53N.mat[m100_53N, m20_53N]
m100_53N_comp10 <- BrayDist_53N.mat[m100_53N, m10_53N]

###########
# EU 8
##########

ASVtab_8 <- psotu2veg(EU_8_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_8 <- vegdist(ASVtab_8, method="bray")
BrayDist_8.mat <- as.matrix(BrayDist_8)
diag(BrayDist_8.mat) <- NA
BrayDist_8.mat[lower.tri(BrayDist_8.mat)] <- NA #because symmetrical, can make lower triangle zero

# Make data frame for indexing matrix
colnames(BrayDist_8.mat) == rownames(sample_data(EU_8_Soils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps8 <- colnames(BrayDist_8.mat) #this is the Mapping Sample ID of these
meter8 <- sample_data(EU_8_Soils.ps)$Meter
index8.df <- data.frame(samps8, meter8)

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
m100_8_comp90 <- BrayDist_8.mat[m100_8, m90_8]
m100_8_comp80 <- BrayDist_8.mat[m100_8, m80_8]
m100_8_comp70 <- BrayDist_8.mat[m100_8, m70_8]
m100_8_comp60 <- BrayDist_8.mat[m100_8, m60_8]
m100_8_comp50 <- BrayDist_8.mat[m100_8, m50_8]
m100_8_comp40 <- BrayDist_8.mat[m100_8, m40_8]
m100_8_comp30 <- BrayDist_8.mat[m100_8, m30_8]
m100_8_comp20 <- BrayDist_8.mat[m100_8, m20_8]
m100_8_comp10 <- BrayDist_8.mat[m100_8, m10_8]

###########
# EU 53S
##########

ASVtab_53S <- psotu2veg(EU_53S_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_53S <- vegdist(ASVtab_53S, method="bray")
BrayDist_53S.mat <- as.matrix(BrayDist_53S)
diag(BrayDist_53S.mat) <- NA
BrayDist_53S.mat[lower.tri(BrayDist_53S.mat)] <- NA #because symmetrical, can make lower triangle zero

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
m100_53S_comp90 <- BrayDist_53S.mat[m100_53S, m90_53S]
m100_53S_comp80 <- BrayDist_53S.mat[m100_53S, m80_53S]
m100_53S_comp70 <- BrayDist_53S.mat[m100_53S, m70_53S]
m100_53S_comp60 <- BrayDist_53S.mat[m100_53S, m60_53S]
m100_53S_comp50 <- BrayDist_53S.mat[m100_53S, m50_53S]
m100_53S_comp40 <- BrayDist_53S.mat[m100_53S, m40_53S]
m100_53S_comp30 <- BrayDist_53S.mat[m100_53S, m30_53S]
m100_53S_comp20 <- BrayDist_53S.mat[m100_53S, m20_53S]
m100_53S_comp10 <- BrayDist_53S.mat[m100_53S, m10_53S]

###########
# EU 10
##########
ASVtab_10 <- psotu2veg(EU_10_Soils.ps)

# Get Bray-Curtis dissimilarities
BrayDist_10 <- vegdist(ASVtab_10, method="bray")
BrayDist_10.mat <- as.matrix(BrayDist_10)
diag(BrayDist_10.mat) <- NA
BrayDist_10.mat[lower.tri(BrayDist_10.mat)] <- NA #because symmetrical, can make lower triangle zero

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
m100_10_comp90 <- BrayDist_10.mat[m100_10, m90_10]
m100_10_comp80 <- BrayDist_10.mat[m100_10, m80_10]
m100_10_comp70 <- BrayDist_10.mat[m100_10, m70_10]
m100_10_comp60 <- BrayDist_10.mat[m100_10, m60_10]
m100_10_comp50 <- BrayDist_10.mat[m100_10, m50_10]
m100_10_comp40 <- BrayDist_10.mat[m100_10, m40_10]
m100_10_comp30 <- BrayDist_10.mat[m100_10, m30_10]
m100_10_comp20 <- BrayDist_10.mat[m100_10, m20_10]
m100_10_comp10 <- BrayDist_10.mat[m100_10, m10_10]


######## MAKE BOXPLOTS ######## 
# So far, this is just for EU 10. Still need to label/distinguish transects and take out outliers
bold_eu10p <- expression(bold("EU 10: Dissimilarity Relative to 10m (patch)"))
bold_eu10f <- expression(bold("EU 10: Dissimilarity Relative to 100m (forest)"))
quartz()
par(mfrow=c(1,2))
patch_box10 <- boxplot(list(m10_10_comp20, m10_10_comp30, m10_10_comp40,
                          m10_10_comp50, m10_10_comp60, m10_10_comp70,
                          m10_10_comp80, m10_10_comp90, m10_10_comp100),
                     ylab = "Bray-Curtis Dissimilarity",
                     names = c("20m", "30m", "40m", "50m",
                               "60m", "70m", "80m", "90m", "100m"), cex.axis = 0.8,
                     cex.lab = 1,
                     ylim=c(0.0, 1.0))
mtext(text=bold_eu10p, side=3, adj = -0.065, line = 2)
forest_box10 <- boxplot(list(m100_10_comp10, m100_10_comp20, m100_10_comp30, m100_10_comp40,
                           m100_10_comp50, m100_10_comp60, m100_10_comp70,
                           m100_10_comp80, m100_10_comp90), 
                      ylab = "Bray-Curtis Dissimilarity",
                      names = c("10m", "20m", "30m", "40m", "50m",
                                "60m", "70m", "80m", "90m"), cex.axis = 0.8,
                      cex.lab = 1,
                      ylim=c(0.0, 1.0))
mtext(text=bold_eu10f, side=3, adj = -0.065, line = 2)

