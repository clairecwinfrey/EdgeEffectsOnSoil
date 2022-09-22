# FungalEdgeEffects- Z-scores
# Begun around August 15, 2022

################################
# ~~~~DESCRIPTION~~~~
################################
# This script identifies gets z-scores for the abundance of differentially abundant fungal taxa
# (the differentially abundant taxa were identified in DiffAbundFungi.R)

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
load("RobjectsSaved/diffAbunDat_tidy_FUNGI") #saved Aug 14, 2022 in EdgeEffectsByASVByEU_FUNGI.R (and DiffAbunFungi.R, I think)

library("tidyverse")
library("growthcurver")

##########################################################################
# 1.  Z-SCORES OF EACH ASV WITHIN EACH EU 
##########################################################################
# 1. Format diffAbunDat_tidy_FUNGI nd a pre-allocated list for analysis.
# Pull out only forest specialists
diffAbunDat_tidy_FUNGI_forest <- diffAbunDat_tidy_FUNGI %>% 
  filter(Habitat == "forest")
# Make a vector of all of the 119 forest-specializing ASVs
ASVnamesDA_FUNGI_forest <- unique(diffAbunDat_tidy_FUNGI_forest$ASV_name) 
# View(diffAbunDat_tidy_FUNGI_forest)
head(diffAbunDat_tidy_FUNGI_forest)

# 1.1 Split up by EU
da_tidy_FUNGI_forest_ByEU <- diffAbunDat_tidy_FUNGI_forest %>% 
  group_split(EU)

# Rename them all according to EU
unique(da_tidy_FUNGI_forest_ByEU[[1]]$EU)
unique(da_tidy_FUNGI_forest_ByEU[[2]]$EU)
unique(da_tidy_FUNGI_forest_ByEU[[3]]$EU)
unique(da_tidy_FUNGI_forest_ByEU[[4]]$EU)
unique(da_tidy_FUNGI_forest_ByEU[[5]]$EU)
unique(da_tidy_FUNGI_forest_ByEU[[6]]$EU)
names(da_tidy_FUNGI_forest_ByEU) <- c("EU_10", "EU_52", "EU_53N", "EU_53S", "EU_54S", "EU_8")

#########################################
# EU 10 
#########################################
# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_10_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name)), ncol=20)) 
colnames(EU_10_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_10_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_10 <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_10 <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_FUNGI.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_10[[i]] <- da_tidy_FUNGI_forest_ByEU$EU_10 %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name)[i])
  uniqueMeters_10[[i]] <- unique(sort(ASVinfo_10[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_10_ASVs_df$EU <- unique(da_tidy_FUNGI_forest_ByEU$EU_10$EU)
  EU_10_ASVs_df$ASV_name[i] <- unique(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name)[i]
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
# View(zScoresEU_10_all)

####### Double check a few ! #######

# 1. Double check first set of nested for loops:
# i. ASV 87, meter 10
# View(da_tidy_FUNGI_forest_ByEU$EU_10)
info_87_10_10a <- filter(da_tidy_FUNGI_forest_ByEU$EU_10, ASV_name=="ASV_87" & Meter == 10)
EU_10_ASVs_df[1,1] ==  mean(info_87_10_10a$ASVabundance)
info_87_10_10b <- filter(da_tidy_FUNGI_forest_ByEU$EU_10, ASV_name=="ASV_87")
mean(info_87_10_10b$ASVabundance) == EU_10_ASVs_df[1,12]
sd(info_87_10_10b$ASVabundance) == EU_10_ASVs_df[1,13]
info_87_tax <- unique(da_tidy_FUNGI_forest_ByEU$EU_10[which(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name=="ASV_87"),10:15]) #pull out all taxonomic info
info_87_tax_2 <-  unique(EU_10_ASVs_df %>% filter(ASV_name=="ASV_87") %>%  
                         select(Phylum:Species))
info_87_tax == info_87_tax_2 #yes, this is the same!

# ii. ASV 7, meter 80
info_7_10_80a <- filter(da_tidy_FUNGI_forest_ByEU$EU_10, ASV_name=="ASV_7" & Meter == 80)
EU_10_ASVs_df[7,8] ==  mean(info_7_10_80a$ASVabundance)
info_7_10_80b <- filter(da_tidy_FUNGI_forest_ByEU$EU_10, ASV_name=="ASV_7")
mean(info_7_10_80b$ASVabundance) == EU_10_ASVs_df[7,12]
sd(info_7_10_80b$ASVabundance) == EU_10_ASVs_df[7,13]

# iii. ASV 1215 (ie. second ASV), meter 50
info_1215_10_50a <- filter(da_tidy_FUNGI_forest_ByEU$EU_10, ASV_name=="ASV_1215" & Meter == 50)
EU_10_ASVs_df[5,5] == mean(info_1215_10_50a$ASVabundance)
info_1215_tax <- unique(da_tidy_FUNGI_forest_ByEU$EU_10[which(da_tidy_FUNGI_forest_ByEU$EU_10$ASV_name=="ASV_1215"),10:15]) #pull out all taxonomic info
info_1215_tax_2 <-  EU_10_ASVs_df %>% filter(ASV_name=="ASV_1215") %>%  
                                select(Phylum:Species)
info_1215_tax == info_1215_tax_2 #yes, this is the same!

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
EU_52_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_FUNGI_forest_ByEU$EU_52$ASV_name)), ncol=20)) 
colnames(EU_52_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_52_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_52 <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_52$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_52 <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_52$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_FUNGI.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_FUNGI_forest_ByEU$EU_52$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_52[[i]] <- da_tidy_FUNGI_forest_ByEU$EU_52 %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_FUNGI_forest_ByEU$EU_52$ASV_name)[i])
  uniqueMeters_52[[i]] <- unique(sort(ASVinfo_52[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_52_ASVs_df$EU <- unique(da_tidy_FUNGI_forest_ByEU$EU_52$EU)
  EU_52_ASVs_df$ASV_name[i] <- unique(da_tidy_FUNGI_forest_ByEU$EU_52$ASV_name)[i]
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

#########################################
# EU 53N
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_53N_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_FUNGI_forest_ByEU$EU_53N$ASV_name)), ncol=20)) 
colnames(EU_53N_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_53N_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_53N <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_53N$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_53N <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_53N$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_FUNGI.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_FUNGI_forest_ByEU$EU_53N$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_53N[[i]] <- da_tidy_FUNGI_forest_ByEU$EU_53N %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_FUNGI_forest_ByEU$EU_53N$ASV_name)[i])
  uniqueMeters_53N[[i]] <- unique(sort(ASVinfo_53N[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_53N_ASVs_df$EU <- unique(da_tidy_FUNGI_forest_ByEU$EU_53N$EU)
  EU_53N_ASVs_df$ASV_name[i] <- unique(da_tidy_FUNGI_forest_ByEU$EU_53N$ASV_name)[i]
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

#########################################
# EU 53S
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_53S_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_FUNGI_forest_ByEU$EU_53S$ASV_name)), ncol=20)) 
colnames(EU_53S_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_53S_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_53S <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_53S$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_53S <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_53S$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_FUNGI.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_FUNGI_forest_ByEU$EU_53S$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_53S[[i]] <- da_tidy_FUNGI_forest_ByEU$EU_53S %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_FUNGI_forest_ByEU$EU_53S$ASV_name)[i])
  uniqueMeters_53S[[i]] <- unique(sort(ASVinfo_53S[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_53S_ASVs_df$EU <- unique(da_tidy_FUNGI_forest_ByEU$EU_53S$EU)
  EU_53S_ASVs_df$ASV_name[i] <- unique(da_tidy_FUNGI_forest_ByEU$EU_53S$ASV_name)[i]
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

#########################################
# EU 54S
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_54S_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_FUNGI_forest_ByEU$EU_54S$ASV_name)), ncol=20)) 
colnames(EU_54S_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_54S_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_54S <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_54S$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_54S <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_54S$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_FUNGI.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_FUNGI_forest_ByEU$EU_54S$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_54S[[i]] <- da_tidy_FUNGI_forest_ByEU$EU_54S %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_FUNGI_forest_ByEU$EU_54S$ASV_name)[i])
  uniqueMeters_54S[[i]] <- unique(sort(ASVinfo_54S[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_54S_ASVs_df$EU <- unique(da_tidy_FUNGI_forest_ByEU$EU_54S$EU)
  EU_54S_ASVs_df$ASV_name[i] <- unique(da_tidy_FUNGI_forest_ByEU$EU_54S$ASV_name)[i]
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

#########################################
# EU 8
#########################################

# 1. Pre-allocate a dataframes
# First is for holding the mean abundance of each ASV at each point along the transect
# columns in each dataframe are EU, all points on transect, ASV EU mean, and ASV EU sd. 15 rows corresponding to the 119 ASVs
EU_8_ASVs_df <- as.data.frame(matrix(nrow=length(unique(da_tidy_FUNGI_forest_ByEU$EU_8$ASV_name)), ncol=20)) 
colnames(EU_8_ASVs_df) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name", "Phylum","Class", "Order", "Family", "Genus", "Species")
EU_8_ASVs_df 
# preallocate thing to store ASV info for each ASV
ASVinfo_8 <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_8$ASV_name))) #119 different storage spots!
# pre-allocate thing to store unique meters!
uniqueMeters_8 <- vector("list", length(unique(da_tidy_FUNGI_forest_ByEU$EU_8$ASV_name)))

#2. Run 2 nested for loops to get data sorted by ASV
# These two for loops below (indexed with i and j): (1) get mean and sd of each ASV's abundance across all samples within 
# the specificed EU and (2): get the mean abundance of each ASV at each meter in the EU. There are a lot of ones because
# in EdgeEffectByASVbyEU_FUNGI.R, I added one to each ASV abundance to make it work with DESeq. 
for (i in 1:length(unique(da_tidy_FUNGI_forest_ByEU$EU_8$ASV_name))){ #this is [[1]] since I am just trying to do EU 10 right now. It goes over all 119 ASVs
  # 1 will be replaced with an index once I do all EUs at once!
  ASVinfo_8[[i]] <- da_tidy_FUNGI_forest_ByEU$EU_8 %>% #this pulls out just the info that we need for each ASV in each EU
    filter(ASV_name == unique(da_tidy_FUNGI_forest_ByEU$EU_8$ASV_name)[i])
  uniqueMeters_8[[i]] <- unique(sort(ASVinfo_8[[i]]$Meter)) #get all the unique meters (this way so that in future cases, 
  # can account for ASV not being present at a given meter). Sort so that it goes in numeric order
  EU_8_ASVs_df$EU <- unique(da_tidy_FUNGI_forest_ByEU$EU_8$EU)
  EU_8_ASVs_df$ASV_name[i] <- unique(da_tidy_FUNGI_forest_ByEU$EU_8$ASV_name)[i]
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

#########################################
# 2. MERGE ALL TOGETHER AND DOUBLE CHECK SOME VALUES
#########################################

# rbind all of these from above together for just one dataframe
fungi_zScores_allEUs <- rbind(zScoresEU_10_all, zScoresEU_52_all, zScoresEU_53N_all, zScoresEU_53S_all, zScoresEU_54S_all, zScoresEU_8_all)
# View(fungi_zScores_allEUs)

# 1. Double check a few:
# i. a few meters in EU 53N, ASV 414
zScoresEU_53N_414 <- filter(zScores_allEUs, ASV_name=="ASV_414" & EU == "EU_53N")
info_53N_414 <- filter(da_tidy_FUNGI_forest_ByEU$EU_53N, ASV_name =="ASV_414")
mean_53N_414 <- mean(info_53N_414$ASVabundance)
st_53N_414 <- sd(info_53N_414$ASVabundance)
mean_53N_414_50m <- mean(info_53N_414$ASVabundance[which(info_53N_414$Meter==50)])
zScoresEU_53N_414$z_50 == (mean_53N_414_50m - mean_53N_414)/st_53N_414 #TRUE!
mean_53N_414_70m <- mean(info_53N_414$ASVabundance[which(info_53N_414$Meter==70)])
zScoresEU_53N_414$z_70 == (mean_53N_414_70m - mean_53N_414)/st_53N_414 #TRUE!
mean_53N_414_100m <- mean(info_53N_414$ASVabundance[which(info_53N_414$Meter==100)])
zScoresEU_53N_414$z_100 == (mean_53N_414_100m - mean_53N_414)/st_53N_414 #TRUE!

# ii. a few meters in EU 8, ASV 14
zScoresEU_8_14 <- filter(zScores_allEUs, ASV_name=="ASV_14" & EU == "EU_8")
info_8_14 <- filter(da_tidy_FUNGI_forest_ByEU$EU_8, ASV_name =="ASV_14")
mean_8_14 <- mean(info_8_14$ASVabundance)
st_8_14 <- sd(info_8_14$ASVabundance)
mean_EU_8_14_20m <- mean(info_8_14$ASVabundance[which(info_8_14$Meter==20)])
zScoresEU_8_14$z_40 == (mean_EU_8_14_20m - mean_8_14)/st_8_14 #TRUE!
mean_EU_8_14_40m <- mean(info_8_14$ASVabundance[which(info_8_14$Meter==40)])
zScoresEU_8_14$z_40 == (mean_EU_8_14_40m - mean_8_14)/st_8_14 #TRUE!

# iii. a few meters in EU 10, ASV 646
zScoresEU_10_646 <- filter(zScores_allEUs, ASV_name=="ASV_646" & EU == "EU_10")
info_10_646 <- filter(da_tidy_FUNGI_forest_ByEU$EU_10, ASV_name =="ASV_646")
mean_10_646 <- mean(info_10_646$ASVabundance)
st_10_646 <- sd(info_10_646$ASVabundance)
mean_EU_10_646_60m <- mean(info_10_646$ASVabundance[which(info_10_646$Meter==60)])
zScoresEU_10_646$z_60 == (mean_EU_10_646_60m - mean_10_646)/st_10_646 #These are both NaNs, so this is working!
mean_EU_10_646_80m <- mean(info_10_646$ASVabundance[which(info_10_646$Meter==80)])
zScoresEU_10_646$z_80 == (mean_EU_10_646_80m - mean_10_646)/st_10_646 #These are both NaNs, so this is working!

# save(fungi_zScores_allEUs, file="RobjectsSaved/fungi_zScores_allEUs") #save it all (last saved Sept 13, 2022)
# Save all of the ASV means used to calculate z-scores

fungiASVmeans_AllEUs <- rbind(EU_10_ASVs_df, EU_52_ASVs_df, EU_53N_ASVs_df, EU_53S_ASVs_df, EU_54S_ASVs_df, EU_8_ASVs_df)
# save(fungiASVmeans_AllEUs, file="RobjectsSaved/fungiASVmeans_AllEUs") #save it all (last saved Sept 21, 2022)

##########################################################################
# 3.  FITTING LOGISTIC CURVES
##########################################################################
# This section is largely based on the explanation for the growthcurver package,
# found here: https://rpubs.com/angelov/growthcurver ######

# Make this in long form so that it works below:
zScores_allEUs_FUNGI_longer <- zScores_allEUs %>% 
  pivot_longer(cols=z_10:z_100, names_to="Meter", values_to = "abundZ_score")
#  weirdly, adding this code here makes the dataframe a list so we'll do it below instead
# lapply(gsub, pattern = "z_", replacement = "", fixed = TRUE) #remove "z_"s so that meter is numeric" 

zScores_allEUs_FUNGI_longer[] <- lapply(zScores_allEUs_FUNGI_longer, gsub, pattern = "z_", replacement = "", fixed = TRUE) #remove "z_"s so that meter is numeric" %>% 
# Make various meters numeric
zScores_allEUs_FUNGI_longer$ASV_EUmean <- as.numeric(zScores_allEUs_FUNGI_longer$ASV_EUmean)
zScores_allEUs_FUNGI_longer$ASV_EUsd <- as.numeric(zScores_allEUs_FUNGI_longer$ASV_EUsd)
zScores_allEUs_FUNGI_longer$Meter <- as.numeric(zScores_allEUs_FUNGI_longer$Meter)
zScores_allEUs_FUNGI_longer$abundZ_score <- as.numeric(zScores_allEUs_FUNGI_longer$abundZ_score)

head(zScores_allEUs_FUNGI_longer) #looks good!

save(zScores_allEUs_FUNGI_longer, file= "RobjectsSaved/zScores_allEUs_FUNGI_longer") #saved September 13, 2022

# What is minimum z-score? (to use the function below, everything needs to be positive)
for (i in 1:10){
  print(min(zScores_allEUs[,i], na.rm=TRUE))
}
#minimum Z-score is -1.104, at meter 30. So adding 1.5 to everything should make everything positive
# (and since the minimum Z-score for the prokaryotic data is -1.442465, which occurs at meter 10, 
# this should work for both)

# First, because these are growth plots, add a 1 to all of the z-scores for correct fit
zScores_allEUs_FUNGI_longerPlus1.5 <- zScores_allEUs_FUNGI_longer
zScores_allEUs_FUNGI_longerPlus1.5$abundZ_score <- zScores_allEUs_FUNGI_longer$abundZ_score + 1.5

# Plot them (raw, i.e. not plus 1)
z_plot1 <- ggplot(zScores_allEUs_FUNGI_longer, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=0.7) +
  theme_bw() + ggtitle("Fungi: All ASV abundance Z_scores (raw)")
quartz()
z_plot1

# Plot plus 1
z_plus1_plot <- ggplot(zScores_allEUs_FUNGI_longerPlus1.5, aes(x = Meter, y = abundZ_score, color=EU)) + geom_point(alpha=0.7) +
  theme_bw() + ggtitle("Fungi: All ASV abundance Z_scores (plus 1.5)")
quartz()
z_plus1_plot

# Making the model (just on one ASV for now!)
ASV_87_df <- zScores_allEUs_FUNGI_longerPlus1.5 %>% filter(ASV_name=="ASV_87")
#View(ASV_87_df)

modelASV87plus1mod <- growthcurver::SummarizeGrowth(data_t=ASV_87_df$Meter, data_n=ASV_87_df$abundZ_score, bg_correct = "none")
modelASV87plus1mod$vals #gives all of the values
predict(modelASV87plus1mod$model) # gives you the predicted abundance values (according to the model)
str(modelASV87plus1mod)
modelASV87plus1mod$vals$r #this is growth rate constant NOT r-squared
modelASV87plus1mod$vals$sigma #0.4782246 -- the smaller the better!
# sigma is a measure of the goodnesss of fit of the parameters of the logistic equation for the data; 
# it is the residual standard error from the nonlinear regression model. Smaller sigma values indicate
# a better fit of the logistic curve to the data than larger values.
#quartz()
plot(predict(modelASV87plus1mod$model))
#quartz()
plot(modelASV87plus1mod$model$m$fitted())
modelASV87plus1mod$model$m$fitted()

### Plotting ##
# Base R
plot(modelASV87plus1mod) #ugly!

# ggplot
unique(ASV_87_df[,c(5:9)]) #get taxonomic information for below
modelASV87plus1mod$vals$t_mid #inflection point = -21.50041
ASV87plus1mod_plot1 <- ggplot(ASV_87_df, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=1.3) + theme_bw()
ASV87plus1mod_plot1
# Adding predicted values
ASV_87_df.predicted <- data.frame(Meter = ASV_87_df$Meter, pred.Zabund = modelASV87plus1mod$model$m$fitted())
ASV87plus1mod_plot2 <- ASV87plus1mod_plot1 + geom_line(data=ASV_87_df.predicted, aes(y=pred.Zabund), color="red", size= 2) + ggtitle("Cortinarius sp. (Basidiomycota)") +
  #geom_point(x=modelASV87plus1mod$vals$t_mid, y = 0.78, color= "red", size=6) + #Here I just guessed y based on how it looked! but inflection point is off the page!
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=2) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) #make it so all meters show on x-axis!
#quartz()
ASV87plus1mod_plot2


ASV_87_df.predicted[c(1:10),] == ASV_87_df.predicted[c(11:20),] #this data frame is longer than it needs to be, but that's okay!
ASV_87_df.predicted[which(ASV_87_df.predicted$Meter==10),2]  # At 10 meters, abundance prediction (+ 1.5) is 1.266708 on the y-axis
# Depth = 80% of distance between edge and abundance at end point, or 
modelASV87plus1mod$vals$t_mid - 10 

# Making the model for all ASVs

#Get the summary metrics for the entire plate of sample data provided
#with the Growthcurver package

#First, load the example data provided with Growthcurver. Note that there is
#a column named "time" -- this is necessary for Growthcurver to know which
#column contains the time measurements. In this dataset, the repeated
#measurements from a single well in a plate are given in a column of data.

myPlate <- growthdata
names(myPlate)
head(myPlate)

#Next, do the analysis for all the columns.
summary_plate <- SummarizeGrowthByPlate(plate = myPlate)

#The output is a data frame that contains the information on the best
#fit for each column of data.
head(summary_plate) 
?SummarizeGrowthByPlate


head(zScores_allEUs_FUNGI_longerPlus1.5)
zScores_allEUs_FUNGI_plate <- zScores_allEUs_FUNGI_longerPlus1.5 %>% 
  pivot_wider(names_from = ASV_name, values_from=abundZ_score) #make dataframe
  # "wider" so that there are columns for every ASV_name (analogous to samples in myPlate example) 
colnames(zScores_allEUs_FUNGI_plate)

zScores_allEUs_FUNGI_plate2 <- zScores_allEUs_FUNGI_plate[, c(10, 11:129, 1:9)] #re-order so meter is first
#View(zScores_allEUs_FUNGI_plate2)
colnames(zScores_allEUs_FUNGI_plate2)[1] <- "time" #make time to be able to use the SummarizeGrowthByPlate function
zScores_allEUs_FUNGI_plate_shorter <- zScores_allEUs_FUNGI_plate[,1:120]
View(zScores_allEUs_FUNGI_plate_shorter)

# Attempt #1 - DOES NOT WORK!
#SummarizeGrowthByPlate(zScores_allEUs_FUNGI_plate_shorter, bg_correct="none", plot_file= "RobjectsSaved/fungiDiffAbundGrowthCurver")

# Another way?
zScores_allEUs_FUNGI_plate_shorter2 <- zScores_allEUs_FUNGI_longerPlus1.5[,c(1,4,11:12)] %>% 
pivot_wider(names_from = ASV_name, values_from=abundZ_score) #make dataframe
colnames(zScores_allEUs_FUNGI_plate_shorter2)
zScores_allEUs_FUNGI_plate_shorter2 <- zScores_allEUs_FUNGI_plate_shorter2[, c(2, 3:121, 1)] #re-order so meter is first
colnames(zScores_allEUs_FUNGI_plate_shorter2)[1] <- "time" #make time to be able to use the SummarizeGrowthByPlate function
zScores_allEUs_FUNGI_plate_shorter2_noEU <- zScores_allEUs_FUNGI_plate_shorter2[,1:120]
zScores_allEUs_FUNGI_plate_shorter2_noEU <- zScores_allEUs_FUNGI_plate_shorter2_noEU %>% 
  mutate_all(function(x) ifelse(is.nan(x), NA, x)) #make all NaNs into NAs for function below
head(zScores_allEUs_FUNGI_plate_shorter2_noEU)
#View(zScores_allEUs_FUNGI_plate_shorter2_noEU)

# Attempt #2 -- THIS WORKED!!!!!!!!!!
fungLogFits <- SummarizeGrowthByPlate(zScores_allEUs_FUNGI_plate_shorter2_noEU, bg_correct="none", plot_file= "RobjectsSaved/fungiDiffAbundGrowthCurver")
str(fungLogFits)
# 13/119 could not be fit, probably because of NAs/NaNs. Will try to run these separately on their own later!

colnames(fungLogFits)[1] <- "ASV_name" #re-name this ASV_name

# Get taxonomic information by merging with another data frame! 
fungLogFitsTaxa <- merge(fungLogFits, unique(zScores_allEUs_FUNGI_longerPlus1.5[,4:9]), by="ASV_name", all.x= FALSE, all.y= FALSE)
# View(fungLogFitsTaxa)

# Now, make a plot showing how taxonomy is affected by the "edge" == t_mid???!
head(fungLogFitsTaxa)
fungLogFits_edgePlot <- ggplot(fungLogFitsTaxa, aes(x=t_mid, fill=Phylum)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  #xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs")
# quartz()
  fungLogFits_edgePlot

fungLogFits_edgePlot <- ggplot(DAphylumAll, aes(fill=Habitat, x=Phylum)) +
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  scale_x_discrete(labels=c("p__Ascomycota" = "Ascomycota", "p__Basidiomycota" = "Basidiomycota",
                            "p__Calcarisporiellomycota" = "Calcarisporiellomycota",
                            "p__Glomeromycota" = "Glomeromycota", "p__Mortierellomycota"="Mortierellomycota",
                            "p__Mucoromycota" = "Mucoromycota", "p__Olpidiomycota"= "Olpidiomycota",
                            "p__Rozellomycota" = "Rozellomycota")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant Fungal ASVs")
# quartz()
diffAbund_ITS_stackedBarplotPhyla


fungLogFitsTaxaLess200 <- fungLogFitsTaxa[which(fungLogFitsTaxa$t_mid <= 200),]
dim(fungLogFitsTaxaLess200)
fungLogFitsTaxaEdgeOver200m <- fungLogFitsTaxa[which(fungLogFitsTaxa$t_mid >= 200),] #15 are over 100 meters
View(fungLogFitsTaxaEdgeOver200m)

questionable_index1 <- which(fungLogFitsTaxa$note=="questionable fit"|fungLogFitsTaxa$note== "cannot fit data") #62 long
questionableFit <- fungLogFitsTaxa[questionable_index1, ]
nrow(questionableFit)  #62
which(questionableFit$t_mid >200)
length(which(questionableFit$t_mid < 0)) #46

questionable_index <- which(fungLogFitsTaxaLess200$note=="questionable fit"|fungLogFitsTaxaLess200$note== "cannot fit data") #62 long
fungLogFitsTaxaLess200onlyGood <- fungLogFitsTaxaLess200[-questionable_index,]
dim(fungLogFitsTaxaLess200onlyGood) #only 37 taxa! :/

cannotFit_index <- which(fungLogFitsTaxaLess200$note== "cannot fit data") #16 long
fungLogFitsTaxaLess200_onlyFit <- fungLogFitsTaxaLess200[-cannotFit_index,]

fungLogFits_edgePlotLess200onlyGood <- ggplot(fungLogFitsTaxaLess200onlyGood, aes(x=t_mid, fill=Phylum, color=Phylum)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=37)")
# quartz()
fungLogFits_edgePlotLess200onlyGood


fungLogFits_edgePlotLess200onlyFitPhylum <- ggplot(fungLogFitsTaxaLess200_onlyFit, aes(x=t_mid, fill=Phylum, color=Phylum)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=83)")
# quartz()
fungLogFits_edgePlotLess200onlyFitPhylum

fungLogFits_edgePlotLess200onlyFitClass <- ggplot(fungLogFitsTaxaLess200_onlyFit, aes(x=t_mid, fill=Class, color=Class)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=83)")
# quartz()
fungLogFits_edgePlotLess200onlyFitClass



fungLogFits_edgePlotLess200onlyGood_Class <- ggplot(fungLogFitsTaxaLess200onlyGood, aes(x=t_mid, fill=Class, color=Class)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=37)")
# quartz()
fungLogFits_edgePlotLess200onlyGood_Class

fungLogFits_edgePlotLess200onlyGood_Order <- ggplot(fungLogFitsTaxaLess200onlyGood, aes(x=t_mid, fill=Order, color=Order)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=37)")
# quartz()
fungLogFits_edgePlotLess200onlyGood_Order

fungLogFits_edgePlotLess200onlyGood_Family <- ggplot(fungLogFitsTaxaLess200onlyGood, aes(x=t_mid, fill=Family, color=Family)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=37)")
# quartz()
fungLogFits_edgePlotLess200onlyGood_Family


# Bring in FUNGuild info! (from FUNGuildExploration.R)
load("RobjectsSaved/fungResultsDA")
fungLogFitsTaxaFungGuildLess200onlyFit <- merge(fungLogFitsTaxaLess200_onlyFit, fungResultsDA[,c(1, 238:242)], by="ASV_name", all.x= TRUE, all.y= FALSE)
# View(fungLogFitsTaxaFungGuildLess200onlyFit)

fungLogFitsFungGuildLess200onlyFit_trophic_plot <- ggplot(fungLogFitsTaxaFungGuildLess200onlyFit, aes(x=t_mid, fill=Trophic.Mode, color=Trophic.Mode)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in trophic mode") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=83)")
# quartz()
fungLogFitsFungGuildLess200onlyFit_trophic_plot

fungLogFitsFungGuildLess200onlyFit_guild_plot <- ggplot(fungLogFitsTaxaFungGuildLess200onlyFit, aes(x=t_mid, fill=Guild, color=Guild)) +
  geom_histogram(binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  #scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in guild") +
  xlab("edge (as inflection point)") +
  ggtitle("Functional Edges of Differentially Abundant Fungal ASVs (n=83)")
# quartz()
fungLogFitsFungGuildLess200onlyFit_guild_plot



