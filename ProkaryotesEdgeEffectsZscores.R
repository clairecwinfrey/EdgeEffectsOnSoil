# ProkaryotesEdgeEffects- Z-scores
# September 13, 2022

################################
# ~~~~DESCRIPTION~~~~
################################
# This script identifies gets z-scores for the abundance of differentially abundant PROKARYOTE
# (bacterial and archaeal) taxa
# (the differentially abundant taxa were identified in DiffAbundProkaryotes.R)

# This script:
# 2. Creates data frame with z-score (of ASV abundance) for each ASV within each EU calculated as:
# --- (1) Mean abundance of each ASV within a given meter on the transect (e.g. mean abundance 
#          of ASV_1 across all transects at meter 10 in EU 10)
# --  (2) Subtracts the mean of that ASV across all points in that EU (e.g. ASV_1's mean value across all 
#           points in EU_10) from the value in (1)
# --  (3) Finally, the difference of (1) - (2) is divided by the standard deviation of the ASVs abundance 
#           across all points in that EU.    
# 3. Fits logistic curves to each one of the differentially abundant ASVs using the Z-scores calculated 
#           above
#IMPORTANTLY, although this is the final z-score script with much of the messier, step-by-step scratch work
# taken out, additional scratch work showing that these calculations  and code work can be found at: _________________________

################################
# SET UP
################################

setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")
load("RobjectsSaved/diffAbunDat_tidy_PROKARYOTES") #from DiffAbundProkaryotes.R (and the largely obsolete 
# EdgeEdgebyASV_allSites.R)

library("tidyverse")

##########################################################################
# 1.  Z-SCORES OF EACH ASV WITHIN EACH EU 
##########################################################################
# 1. Format diffAbunDat_tidy_prokary nd a pre-allocated list for analysis.
# Pull out only forest specialists
diffAbunDat_tidy_prokary_forest <- diffAbunDat_tidy_PROKARYOTES %>% 
  filter(Habitat == "forest")
# View(diffAbunDat_tidy_prokary_forest)
head(diffAbunDat_tidy_prokary_forest)

# 1.1 Split up by EU
da_tidy_prokary_forest_ByEU <- diffAbunDat_tidy_prokary_forest %>% 
  group_split(EU)

# Rename them all according to EU
unique(da_tidy_prokary_forest_ByEU[[1]]$EU)
unique(da_tidy_prokary_forest_ByEU[[2]]$EU)
unique(da_tidy_prokary_forest_ByEU[[3]]$EU)
unique(da_tidy_prokary_forest_ByEU[[4]]$EU)
unique(da_tidy_prokary_forest_ByEU[[5]]$EU)
unique(da_tidy_prokary_forest_ByEU[[6]]$EU)
names(da_tidy_prokary_forest_ByEU) <- c("EU_10", "EU_52", "EU_53N", "EU_53S", "EU_54S", "EU_8")

# How big should each dataframe be?
length(unique(diffAbunDat_tidy_prokary_forest$ASV_name)) #1121

#########################################
# EU 10 
#########################################
# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_10_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_prokary_forest_ByEU$EU_10$ASV_name)), ncol=20)) 
colnames(EU_10_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_10_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_10 <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_10$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_10 <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_10$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_prokary.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_prokary_forest_ByEU$EU_10$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_10[[i]] <- da_tidy_prokary_forest_ByEU$EU_10 %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_prokary_forest_ByEU$EU_10$ASV_name)[i])
  uniqueMeters_10[[i]] <- unique(sort(ASVinfo_10[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_10_ASVs_df$EU <- unique(da_tidy_prokary_forest_ByEU$EU_10$EU)
  EU_10_ASVs_df$ASV_name[i] <- unique(da_tidy_prokary_forest_ByEU$EU_10$ASV_name)[i]
  EU_10_ASVs_df$ASV_EUmean[i] <- mean(ASVinfo_10[[i]]$ASVabundance)
  EU_10_ASVs_df$ASV_EUsd[i] <- sd(ASVinfo_10[[i]]$ASVabundance)
  EU_10_ASVs_df[i,15:20] <- unique(ASVinfo_10[[i]][,10:15])
  for (j in 1:length(uniqueMeters_10[[1]])) { #loop over all of the unique meters
    EU_10_ASVs_df[i,j] <- mean(filter(ASVinfo_10[[i]], Meter == uniqueMeters_10[[i]][j])$ASVabundance)
  }
}
#View(EU_10_ASVs_df)
# This for loop below gets the z-scores of each ASV within the given EU. 
# First, pre-allocate something to store z-scores 
zScoresEU_10 <- matrix(data= NA, nrow=nrow(EU_10_ASVs_df), ncol=10)
colnames(zScoresEU_10) <- c("z_10", "z_20", "z_30", "z_40", "z_50", "z_60", "z_70", "z_80", "z_90", "z_100")
for (k in 1:nrow(EU_10_ASVs_df)){ #this for loop goes over each row
  for (m in 1:10){ #this for loop goes for each meter in each row (of which there are ten)
    zScoresEU_10[k,m] <-  (EU_10_ASVs_df[k,m] - EU_10_ASVs_df$ASV_EUmean[k])/EU_10_ASVs_df$ASV_EUsd[k]
  }
}
zScoresEU_10_all <- cbind(zScoresEU_10, EU_10_ASVs_df[,11:20])
dim(zScoresEU_10_all)
# View(zScoresEU_10_all)

####### Double check a few ! #######

# 1. Double check first set of nested for loops:
# i. ASV 87, meter 10
# View(da_tidy_prokary_forest_ByEU$EU_10)
info_3754_10_10a <- filter(da_tidy_prokary_forest_ByEU$EU_10, ASV_name=="ASV_3754" & Meter == 10)
EU_10_ASVs_df$mean_10[which(EU_10_ASVs_df$ASV_name=="ASV_3754")] ==  mean(info_3754_10_10a$ASVabundance)
info_3754_10 <- filter(da_tidy_prokary_forest_ByEU$EU_10, ASV_name=="ASV_3754")
mean(info_3754_10$ASVabundance) == EU_10_ASVs_df[which(EU_10_ASVs_df$ASV_name=="ASV_3754"),12] #TRUE
sd(info_3754_10$ASVabundance) == EU_10_ASVs_df$ASV_EUsd[which(EU_10_ASVs_df$ASV_name=="ASV_3754")] #TRUE
info_3754_tax <- unique(da_tidy_prokary_forest_ByEU$EU_10[which(da_tidy_prokary_forest_ByEU$EU_10$ASV_name=="ASV_3754"),10:15]) #pull out all taxonomic info
info_3754_tax_2 <-  info_3754_10_10a %>% filter(ASV_name=="ASV_3754")
info_3754_tax == unique(info_3754_tax_2[,10:15]) #yes, this is the same!

# ii. ASV 7, meter 80
info_444_10_80a <- filter(da_tidy_prokary_forest_ByEU$EU_10, ASV_name=="ASV_444" & Meter == 80)
EU_10_ASVs_df$mean_80[which(EU_10_ASVs_df$ASV_name=="ASV_444")] ==  mean(info_444_10_80a$ASVabundance)
info_444_10 <- filter(da_tidy_prokary_forest_ByEU$EU_10, ASV_name=="ASV_444")
mean(info_444_10$ASVabundance) == EU_10_ASVs_df[which(EU_10_ASVs_df$ASV_name=="ASV_444"),12] #TRUE
sd(info_444_10$ASVabundance) == EU_10_ASVs_df$ASV_EUsd[which(EU_10_ASVs_df$ASV_name=="ASV_444")] #TRUE
info_444_tax <- unique(da_tidy_prokary_forest_ByEU$EU_10[which(da_tidy_prokary_forest_ByEU$EU_10$ASV_name=="ASV_444"),10:15]) #pull out all taxonomic info
info_444_tax_2 <-  info_444_10 %>% filter(ASV_name=="ASV_444")
info_444_tax == unique(info_444_tax_2[,10:15]) #yes, this is the same!

# 2. Double check second set of nested for loops:
# CHECKS (simple coding since I already checked that EU_10_ASVs_df was working)
zScoresEU_10_all[1,1] == (EU_10_ASVs_df[1,1] - EU_10_ASVs_df$ASV_EUmean[1])/EU_10_ASVs_df$ASV_EUsd[1]
zScoresEU_10_all[3,10] == (EU_10_ASVs_df[3,10] - EU_10_ASVs_df$ASV_EUmean[3])/EU_10_ASVs_df$ASV_EUsd[3]
zScoresEU_10_all[15,10] == (EU_10_ASVs_df[15,10] - EU_10_ASVs_df$ASV_EUmean[15])/EU_10_ASVs_df$ASV_EUsd[15]


#########################################
# EU 52
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_52_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_prokary_forest_ByEU$EU_52$ASV_name)), ncol=20)) 
colnames(EU_52_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_52_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_52 <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_52$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_52 <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_52$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_prokary.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_prokary_forest_ByEU$EU_52$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_52[[i]] <- da_tidy_prokary_forest_ByEU$EU_52 %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_prokary_forest_ByEU$EU_52$ASV_name)[i])
  uniqueMeters_52[[i]] <- unique(sort(ASVinfo_52[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_52_ASVs_df$EU <- unique(da_tidy_prokary_forest_ByEU$EU_52$EU)
  EU_52_ASVs_df$ASV_name[i] <- unique(da_tidy_prokary_forest_ByEU$EU_52$ASV_name)[i]
  EU_52_ASVs_df$ASV_EUmean[i] <- mean(ASVinfo_52[[i]]$ASVabundance)
  EU_52_ASVs_df$ASV_EUsd[i] <- sd(ASVinfo_52[[i]]$ASVabundance)
  EU_52_ASVs_df[i,15:20] <- unique(ASVinfo_52[[i]][,10:15])
  for (j in 1:length(uniqueMeters_52[[1]])) { #loop over all of the unique meters
    EU_52_ASVs_df[i,j] <- mean(filter(ASVinfo_52[[i]], Meter == uniqueMeters_52[[i]][j])$ASVabundance)
  }
}

# This for loop below gets the z-scores of each ASV within the given EU. 
# First, pre-allocate something to store z-scores 
zScoresEU_52 <- matrix(data= NA, nrow=nrow(EU_52_ASVs_df), ncol=10)
colnames(zScoresEU_52) <- c("z_10", "z_20", "z_30", "z_40", "z_50", "z_60", "z_70", "z_80", "z_90", "z_100")
for (k in 1:nrow(EU_52_ASVs_df)){ #this for loop goes over each row
  for (m in 1:10){ #this for loop goes for each meter in each row (of which there are ten)
    zScoresEU_52[k,m] <-  (EU_52_ASVs_df[k,m] - EU_52_ASVs_df$ASV_EUmean[k])/EU_52_ASVs_df$ASV_EUsd[k]
  }
}
zScoresEU_52_all <- cbind(zScoresEU_52, EU_52_ASVs_df[,11:20])
# View(zScoresEU_52_all)
dim(zScoresEU_52_all)

#########################################
# EU 53N
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_53N_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_prokary_forest_ByEU$EU_53N$ASV_name)), ncol=20)) 
colnames(EU_53N_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_53N_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_53N <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_53N$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_53N <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_53N$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_prokary.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_prokary_forest_ByEU$EU_53N$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_53N[[i]] <- da_tidy_prokary_forest_ByEU$EU_53N %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_prokary_forest_ByEU$EU_53N$ASV_name)[i])
  uniqueMeters_53N[[i]] <- unique(sort(ASVinfo_53N[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_53N_ASVs_df$EU <- unique(da_tidy_prokary_forest_ByEU$EU_53N$EU)
  EU_53N_ASVs_df$ASV_name[i] <- unique(da_tidy_prokary_forest_ByEU$EU_53N$ASV_name)[i]
  EU_53N_ASVs_df$ASV_EUmean[i] <- mean(ASVinfo_53N[[i]]$ASVabundance)
  EU_53N_ASVs_df$ASV_EUsd[i] <- sd(ASVinfo_53N[[i]]$ASVabundance)
  EU_53N_ASVs_df[i,15:20] <- unique(ASVinfo_53N[[i]][,10:15])
  for (j in 1:length(uniqueMeters_53N[[1]])) { #loop over all of the unique meters
    EU_53N_ASVs_df[i,j] <- mean(filter(ASVinfo_53N[[i]], Meter == uniqueMeters_53N[[i]][j])$ASVabundance)
  }
}

# This for loop below gets the z-scores of each ASV within the given EU. 
# First, pre-allocate something to store z-scores 
zScoresEU_53N <- matrix(data= NA, nrow=nrow(EU_53N_ASVs_df), ncol=10)
colnames(zScoresEU_53N) <- c("z_10", "z_20", "z_30", "z_40", "z_50", "z_60", "z_70", "z_80", "z_90", "z_100")
for (k in 1:nrow(EU_53N_ASVs_df)){ #this for loop goes over each row
  for (m in 1:10){ #this for loop goes for each meter in each row (of which there are ten)
    zScoresEU_53N[k,m] <-  (EU_53N_ASVs_df[k,m] - EU_53N_ASVs_df$ASV_EUmean[k])/EU_53N_ASVs_df$ASV_EUsd[k]
  }
}
zScoresEU_53N_all <- cbind(zScoresEU_53N, EU_53N_ASVs_df[,11:20])
# View(zScoresEU_53N_all)
dim(zScoresEU_53N_all)

#########################################
# EU 53S
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_53S_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_prokary_forest_ByEU$EU_53S$ASV_name)), ncol=20)) 
colnames(EU_53S_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_53S_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_53S <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_53S$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_53S <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_53S$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_prokary.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_prokary_forest_ByEU$EU_53S$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_53S[[i]] <- da_tidy_prokary_forest_ByEU$EU_53S %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_prokary_forest_ByEU$EU_53S$ASV_name)[i])
  uniqueMeters_53S[[i]] <- unique(sort(ASVinfo_53S[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_53S_ASVs_df$EU <- unique(da_tidy_prokary_forest_ByEU$EU_53S$EU)
  EU_53S_ASVs_df$ASV_name[i] <- unique(da_tidy_prokary_forest_ByEU$EU_53S$ASV_name)[i]
  EU_53S_ASVs_df$ASV_EUmean[i] <- mean(ASVinfo_53S[[i]]$ASVabundance)
  EU_53S_ASVs_df$ASV_EUsd[i] <- sd(ASVinfo_53S[[i]]$ASVabundance)
  EU_53S_ASVs_df[i,15:20] <- unique(ASVinfo_53S[[i]][,10:15])
  for (j in 1:length(uniqueMeters_53S[[1]])) { #loop over all of the unique meters
    EU_53S_ASVs_df[i,j] <- mean(filter(ASVinfo_53S[[i]], Meter == uniqueMeters_53S[[i]][j])$ASVabundance)
  }
}

# This for loop below gets the z-scores of each ASV within the given EU. 
# First, pre-allocate something to store z-scores 
zScoresEU_53S <- matrix(data= NA, nrow=nrow(EU_53S_ASVs_df), ncol=10)
colnames(zScoresEU_53S) <- c("z_10", "z_20", "z_30", "z_40", "z_50", "z_60", "z_70", "z_80", "z_90", "z_100")
for (k in 1:nrow(EU_53S_ASVs_df)){ #this for loop goes over each row
  for (m in 1:10){ #this for loop goes for each meter in each row (of which there are ten)
    zScoresEU_53S[k,m] <-  (EU_53S_ASVs_df[k,m] - EU_53S_ASVs_df$ASV_EUmean[k])/EU_53S_ASVs_df$ASV_EUsd[k]
  }
}
zScoresEU_53S_all <- cbind(zScoresEU_53S, EU_53S_ASVs_df[,11:20])
# View(zScoresEU_53S_all)
dim(zScoresEU_53S_all)

#########################################
# EU 54S
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_54S_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_prokary_forest_ByEU$EU_54S$ASV_name)), ncol=20)) 
colnames(EU_54S_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_54S_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_54S <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_54S$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_54S <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_54S$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_prokary.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_prokary_forest_ByEU$EU_54S$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_54S[[i]] <- da_tidy_prokary_forest_ByEU$EU_54S %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_prokary_forest_ByEU$EU_54S$ASV_name)[i])
  uniqueMeters_54S[[i]] <- unique(sort(ASVinfo_54S[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_54S_ASVs_df$EU <- unique(da_tidy_prokary_forest_ByEU$EU_54S$EU)
  EU_54S_ASVs_df$ASV_name[i] <- unique(da_tidy_prokary_forest_ByEU$EU_54S$ASV_name)[i]
  EU_54S_ASVs_df$ASV_EUmean[i] <- mean(ASVinfo_54S[[i]]$ASVabundance)
  EU_54S_ASVs_df$ASV_EUsd[i] <- sd(ASVinfo_54S[[i]]$ASVabundance)
  EU_54S_ASVs_df[i,15:20] <- unique(ASVinfo_54S[[i]][,10:15])
  for (j in 1:length(uniqueMeters_54S[[1]])) { #loop over all of the unique meters
    EU_54S_ASVs_df[i,j] <- mean(filter(ASVinfo_54S[[i]], Meter == uniqueMeters_54S[[i]][j])$ASVabundance)
  }
}

# This for loop below gets the z-scores of each ASV within the given EU. 
# First, pre-allocate something to store z-scores 
zScoresEU_54S <- matrix(data= NA, nrow=nrow(EU_54S_ASVs_df), ncol=10)
colnames(zScoresEU_54S) <- c("z_10", "z_20", "z_30", "z_40", "z_50", "z_60", "z_70", "z_80", "z_90", "z_100")
for (k in 1:nrow(EU_54S_ASVs_df)){ #this for loop goes over each row
  for (m in 1:10){ #this for loop goes for each meter in each row (of which there are ten)
    zScoresEU_54S[k,m] <-  (EU_54S_ASVs_df[k,m] - EU_54S_ASVs_df$ASV_EUmean[k])/EU_54S_ASVs_df$ASV_EUsd[k]
  }
}
zScoresEU_54S_all <- cbind(zScoresEU_54S, EU_54S_ASVs_df[,11:20])
# View(zScoresEU_54S_all)
dim(zScoresEU_54S_all)

#########################################
# EU 8
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_8_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_prokary_forest_ByEU$EU_8$ASV_name)), ncol=20)) 
colnames(EU_8_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_8_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_8 <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_8$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_8 <- vector("list", length(unique(da_tidy_prokary_forest_ByEU$EU_8$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_prokary.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_prokary_forest_ByEU$EU_8$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_8[[i]] <- da_tidy_prokary_forest_ByEU$EU_8 %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_prokary_forest_ByEU$EU_8$ASV_name)[i])
  uniqueMeters_8[[i]] <- unique(sort(ASVinfo_8[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_8_ASVs_df$EU <- unique(da_tidy_prokary_forest_ByEU$EU_8$EU)
  EU_8_ASVs_df$ASV_name[i] <- unique(da_tidy_prokary_forest_ByEU$EU_8$ASV_name)[i]
  EU_8_ASVs_df$ASV_EUmean[i] <- mean(ASVinfo_8[[i]]$ASVabundance)
  EU_8_ASVs_df$ASV_EUsd[i] <- sd(ASVinfo_8[[i]]$ASVabundance)
  EU_8_ASVs_df[i,15:20] <- unique(ASVinfo_8[[i]][,10:15])
  for (j in 1:length(uniqueMeters_8[[1]])) { #loop over all of the unique meters
    EU_8_ASVs_df[i,j] <- mean(filter(ASVinfo_8[[i]], Meter == uniqueMeters_8[[i]][j])$ASVabundance)
  }
}

# This for loop below gets the z-scores of each ASV within the given EU. 
# First, pre-allocate something to store z-scores 
zScoresEU_8 <- matrix(data= NA, nrow=nrow(EU_8_ASVs_df), ncol=10)
colnames(zScoresEU_8) <- c("z_10", "z_20", "z_30", "z_40", "z_50", "z_60", "z_70", "z_80", "z_90", "z_100")
for (k in 1:nrow(EU_8_ASVs_df)){ #this for loop goes over each row
  for (m in 1:10){ #this for loop goes for each meter in each row (of which there are ten)
    zScoresEU_8[k,m] <-  (EU_8_ASVs_df[k,m] - EU_8_ASVs_df$ASV_EUmean[k])/EU_8_ASVs_df$ASV_EUsd[k]
  }
}
zScoresEU_8_all <- cbind(zScoresEU_8, EU_8_ASVs_df[,11:20])
# View(zScoresEU_8_all)
dim(zScoresEU_8_all)

#########################################
# 2. MERGE ALL TOGETHER AND DOUBLE CHECK SOME FALUES
#########################################

# rbind all of these from above together for just one dataframe
zScores_allEUs <- rbind(zScoresEU_10_all, zScoresEU_52_all, zScoresEU_53N_all, zScoresEU_53S_all, zScoresEU_54S_all, zScoresEU_8_all)
# View(zScores_allEUs)
nrow(zScores_allEUs) == 1121*6 #as long as expected.

# 1. Double check a few:
# i. a few meters in EU 53N, ASV 1025 -- ALL CORRECT
zScoresEU_53N_1025 <- filter(zScores_allEUs, ASV_name=="ASV_1025" & EU == "EU_53N")
info_53N_1025 <- filter(da_tidy_prokary_forest_ByEU$EU_53N, ASV_name =="ASV_1025")
mean_53N_1025 <- mean(info_53N_1025$ASVabundance)
st_53N_1025 <- sd(info_53N_1025$ASVabundance)
mean_53N_1025_50m <- mean(info_53N_1025$ASVabundance[which(info_53N_1025$Meter==50)])
zScoresEU_53N_1025$z_50 == (mean_53N_1025_50m - mean_53N_1025)/st_53N_1025 #TRUE!
mean_53N_1025_70m <- mean(info_53N_1025$ASVabundance[which(info_53N_1025$Meter==70)])
zScoresEU_53N_1025$z_70 == (mean_53N_1025_70m - mean_53N_1025)/st_53N_1025 #TRUE!
mean_53N_1025_100m <- mean(info_53N_1025$ASVabundance[which(info_53N_1025$Meter==100)])
zScoresEU_53N_1025$z_100 == (mean_53N_1025_100m - mean_53N_1025)/st_53N_1025 #TRUE!

# ii. a few meters in EU 8, ASV 4320-- all correct!
zScoresEU_8_4320 <- filter(zScores_allEUs, ASV_name=="ASV_4320" & EU == "EU_8")
info_8_4320 <- filter(da_tidy_prokary_forest_ByEU$EU_8, ASV_name =="ASV_4320")
mean_8_4320 <- mean(info_8_4320$ASVabundance)
st_8_4320 <- sd(info_8_4320$ASVabundance)
mean_EU_8_4320_20m <- mean(info_8_4320$ASVabundance[which(info_8_4320$Meter==20)])
zScoresEU_8_4320$z_20 == (mean_EU_8_4320_20m - mean_8_4320)/st_8_4320 #TRUE!
mean_EU_8_4320_40m <- mean(info_8_4320$ASVabundance[which(info_8_4320$Meter==40)])
zScoresEU_8_4320$z_40 == (mean_EU_8_4320_40m - mean_8_4320)/st_8_4320 #TRUE!

# iii. a few meters in EU 10, ASV 99 - all correct!
zScoresEU_10_99 <- filter(zScores_allEUs, ASV_name=="ASV_99" & EU == "EU_10")
info_10_99 <- filter(da_tidy_prokary_forest_ByEU$EU_10, ASV_name =="ASV_99")
mean_10_99 <- mean(info_10_99$ASVabundance)
st_10_99 <- sd(info_10_99$ASVabundance)
mean_EU_10_99_60m <- mean(info_10_99$ASVabundance[which(info_10_99$Meter==60)])
zScoresEU_10_99$z_60 == (mean_EU_10_99_60m - mean_10_99)/st_10_99 #These are both NaNs, so this is working!
mean_EU_10_99_80m <- mean(info_10_99$ASVabundance[which(info_10_99$Meter==80)])
zScoresEU_10_99$z_80 == (mean_EU_10_99_80m - mean_10_99)/st_10_99 #These are both NaNs, so this is working!


##########################################################################
# 3.  FITTING LOGISTIC CURVES
##########################################################################

# What is minimum z-score (so we know what to add to make positive)
for (i in 1:10){
  print(min(zScores_allEUs[,i], na.rm=TRUE)) #-1.442465, which occurs at meter 10. So we should be able to add 1.5 or 2.
}
