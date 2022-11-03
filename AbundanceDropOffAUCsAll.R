# AbundanceDropOff.R
# Started November 2, 2022

# Description: This script looks at how the abundance of different ASVs changes across the transect. Data is processed
# for both prokarytoes and fungi here (but separately). Here, I fit a flexible spline to the abundance data from 
# 10-100 meters (but 90 and 100 meters are averaged so that there are 40 m in patch and forest), get the
# area under this curve (AUC), and then find where along the transect (i.e. which meter) the AUC is 25% that of the
# maximum. 

################################
options(scipen = 999) #make it so that annoying scientific notation doesn't pop up

library(tidyverse)
library(splines)
library(gridExtra)
library(phyloseq)
library(purrr)
library(npreg) #for ss function
library(flux) #for auc function

# setwd:
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")
load("RobjectsSaved/diffAbunDat_no53N_tidy_FUNGI") #made in DiffAbundFungi
load("RobjectsSaved/diffAbunDat_no53N_PROKARYOTES") #object is diffAbunDat_no53N_tidy_PROKARYOTES

#View(diffAbunDat_no53N_tidy_FUNGI)

################################
# MAIN PART OF SCRIPT

########################################
##### 1. SETTING UP DATA FOR ANALYSES
########################################
######## FUNGI ########
### i. Isolate only ASVs that were differentially abundant and higher in the forest 
# View(diffAbunDat_no53N_tidy_FUNGI) #number of rows is equal to the 5 EUs times 10 points along the transect times 114 differentially abundant fungal taxa
length(unique(diffAbunDat_no53N_tidy_FUNGI$ASV_name))*193 == nrow(diffAbunDat_no53N_tidy_FUNGI) #number of rows is equal to 193 samples times 279 differentially abundant taxa
diffAbunDat_no53N_tidy_FUNGI_forest <- diffAbunDat_no53N_tidy_FUNGI %>% filter(Habitat== "forest")
# View(diffAbunDat_no53N_tidy_FUNGI_forest) # ASVabundance is abundance at each site

### ii. Make dataframes for each of the ASVs in each EU (so 5 EUs x 114 ASVs = 570 dataframes)
diffAbun_no53N_tidy_FUNGI_forest_Lists <- diffAbunDat_no53N_tidy_FUNGI_forest %>% 
  group_split(EU, ASV_name)
length(diffAbun_no53N_tidy_FUNGI_forest_Lists) #5 EUs x 114 ASVs = 570 dataframes
# Investigate how this broke it up to then check names:
fungiOrderedASVs <- vector(length=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
fungiOrderedEUs <- vector(length=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
for (i in 1:length(diffAbun_no53N_tidy_FUNGI_forest_Lists)){
  fungiOrderedASVs[i] <- print(diffAbun_no53N_tidy_FUNGI_forest_Lists[[i]]$ASV_name)
  fungiOrderedEUs[i] <- print(diffAbun_no53N_tidy_FUNGI_forest_Lists[[i]]$EU)
}
fungiOrderedASVs #order is 1, 10, 1060, 1069, 1081, 109, 110, 1127, 113, 115, etc. 
fungiOrderedEUs #order is 10, 52, 53S, 54S, 8

# Rename names of items in list based on combining EU and ASV information
# Make a vector for naming the lists
fungiMergedNames <- vector(length=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
for (i in 1:length(fungiOrderedASVs)){
  fungiMergedNames[i] <- paste(fungiOrderedEUs[i], fungiOrderedASVs[i], sep="_")
}
names(diffAbun_no53N_tidy_FUNGI_forest_Lists) <- fungiMergedNames #give names
# Check a few (they all look as expected!)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_10_ASV_1203$ASV_name)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_10_ASV_1203$EU)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_53S_ASV_591$ASV_name)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_53S_ASV_591$EU)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_8_ASV_97$ASV_name)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_8_ASV_97$EU)

######## PROKARYOTES ########
### i. Isolate only ASVs that were differentially abundant and higher in the forest 
length(unique(diffAbunDat_no53N_tidy_PROKARYOTES$ASV_name))*194 == nrow(diffAbunDat_no53N_tidy_PROKARYOTES) #number of rows is equal 194 samples x 2169 differentially abundant taxa
diffAbunDat_no53N_tidy_PROKARYOTES_forest <- diffAbunDat_no53N_tidy_PROKARYOTES %>% filter(Habitat== "forest")
# View(diffAbunDat_no53N_tidy_PROKARYOTES_forest) # ASVabundance is abundance at each site
length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))

### ii. Make dataframes for each of the ASVs in each EU (so 5 EUs x 1099 ASVs = 5495 dataframes in list)
diffAbun_no53N_tidy_PROKARYOTES_forest_Lists <- diffAbunDat_no53N_tidy_PROKARYOTES_forest %>% 
  group_split(EU, ASV_name)
length(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists) #5 EUs x 1099 ASVs = 5495 dataframes
# Investigate how this broke it up to then check names:
proOrderedASVs <- vector(length=length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))*5)
proOrderedEUs <- vector(length=length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))*5)
for (i in 1:length(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists)){
  proOrderedASVs[i] <- print(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[i]]$ASV_name)
  proOrderedEUs[i] <- print(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[i]]$EU)
}
proOrderedASVs #order is 1002, 1006, 1060, 1009, etc.
unique(proOrderedEUs) #order is 10, 52, 53S, 54S, 8

# Rename names of items in list based on combining EU and ASV information
# Make a vector for naming the lists
proMergedNames <- c()
for (i in 1:length(proOrderedASVs)){
  proMergedNames[i] <- paste(proOrderedEUs[i], proOrderedASVs[i], sep="_")
}
names(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists) <- proMergedNames #give names
# Check a few (they all look as expected!)
unique(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists$EU_10_ASV_4758$ASV_name) #4758
unique(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists$EU_52_ASV_5330$EU) #EU 52
unique(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists$EU_52_ASV_5330$ASV_name) #"ASV_5330"
unique(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists$EU_52_ASV_523$EU) #EU 52
unique(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists$EU_54S_ASV_836$ASV_name) #836

########################################
##### 2. REPLACE METERS 90 AND 100 WITH MEAN OF 90 AND 100
########################################
######## FUNGI ########
##### i. Get mean abundance of 90 and 100 meters for each ASV within each EU, for each transect
### i.i first, get mean abundance from 90-100 meters in forest.
# Pre-allocate dataframes
fungi_forestMeansTransects <- as.data.frame(matrix(nrow=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5, ncol=7)) #570 rows (115 ASVs x 5 EUs)
colnames(fungi_forestMeansTransects) <- c("EU", "ASV_name", "ASVplusEU", "B_meanForestAbund", "L_meanForestAbund", "R_meanForestAbund", "T_meanForestAbund")
df <- as.data.frame(matrix(nrow=12, ncol=5)) #12 rows, meters 90 and 100 x (up to) 4 transects
colnames(df) <- c("ASV_name", "EU", "Transect", "Meter", "ASV_abundance")
subsetDats <- rep(list(df), length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5) #need one for each ASV and EU combo, i.e. 570

# Make for loop that gets the mean from 90-100 in all of the elements of the list "diffAbun_no53N_tidy_FUNGI_forest_Lists"
for (j in 1:length(diffAbun_no53N_tidy_FUNGI_forest_Lists)){ #gets ASV name, EU, transect, and meter for each EU/ASV combo, but only for meters 90-100
  # The line below fills in as many rows as are needed, hence using the part with %in% twice
  subsetDats[[j]][1:nrow(diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]][diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]]$Meter %in% c(90, 100),]),] <- 
    diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]][diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]]$Meter %in% c(90, 100),c(2,19:21, 17)]
  fungi_forestMeansTransects[j,1] <- unique(na.omit(subsetDats[[j]]$EU))
  fungi_forestMeansTransects[j,2] <- unique(na.omit(subsetDats[[j]]$ASV_name))
  fungi_forestMeansTransects[j,3] <- paste(fungi_forestMeansTransects[j,1], fungi_forestMeansTransects[j,2], sep="_") #make this name match those of "diffAbun_no53N_tidy_FUNGI_forest_Lists"
  fungi_forestMeansTransects[j,4] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "B"),]$ASV_abundance)
  fungi_forestMeansTransects[j,5] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "L"),]$ASV_abundance)
  fungi_forestMeansTransects[j,6] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "R"),]$ASV_abundance)
  fungi_forestMeansTransects[j,7] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "T"),]$ASV_abundance)
}

# Checking a few! 
mean(subsetDats[[1]][which(subsetDats[[1]]$Transect == "T"),]$ASV_abundance) == fungi_forestMeansTransects$T_meanForestAbund[1]
mean(subsetDats[[333]][which(subsetDats[[333]]$Transect == "B"),]$ASV_abundance) == fungi_forestMeansTransects$B_meanForestAbund[333]
mean(subsetDats[[509]][which(subsetDats[[509]]$Transect == "R"),]$ASV_abundance) == fungi_forestMeansTransects$R_meanForestAbund[509]

#############
### i.ii. Put mean abundance from 90-100 for each transect in place of separate 90 and 100 values
# Make all of these dataframes (not tibbles) so that they are compatible with some of the code below
diffAbun_FUNGI_comboForest_dfLists <- vector("list", length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5) #make a empty list to store them in
for (i in 1:length(diffAbun_FUNGI_comboForest_dfLists)){
  diffAbun_FUNGI_comboForest_dfLists[[i]] <- as.data.frame(diffAbun_no53N_tidy_FUNGI_forest_Lists[[i]])
  print(class(diffAbun_FUNGI_comboForest_dfLists[[i]])) #make sure that they are all data frames
}

# Make names equal to the original objects
names(diffAbun_FUNGI_comboForest_dfLists) <- names(diffAbun_no53N_tidy_FUNGI_forest_Lists)

# Replace 90-100 meter values in diffAbun_no53N_tidy_FUNGI_forest_Lists_forestCombo with fungi_forestMeans values. Will call it meter 90, for plotting purposes
# Preallocate a 114*5 long list of vector for each ASV in each EU to hold the indices for 90 and 100m
indicestoKeep <- vector("list", length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
# Outer for loop removes meters 90 and 100; inner for loop replaces meter 80 with mean of each transect across meters 80, 90, and 100
for (k in 1:length(diffAbun_FUNGI_comboForest_dfLists)){ 
  indicestoKeep[[k]] <- which(diffAbun_FUNGI_comboForest_dfLists[[k]]$Meter %in% c(10, 20, 30, 40, 50, 60, 70, 80, 90)==TRUE)
  diffAbun_FUNGI_comboForest_dfLists[[k]] <- diffAbun_FUNGI_comboForest_dfLists[[k]][indicestoKeep[[k]],] #keep only meters 10-90
  for (w in 1:4){ #loop over the four different transects
    diffAbun_FUNGI_comboForest_dfLists[[k]][which(diffAbun_FUNGI_comboForest_dfLists[[k]]$Transect %in% 
                                                    unique(sort(diffAbun_FUNGI_comboForest_dfLists[[k]]$Transect))[w] & 
                                                    diffAbun_FUNGI_comboForest_dfLists[[k]]$Meter == 90),17] <- #make 90 get the mean value for 90=100
      fungi_forestMeansTransects[k,3+w] 
  }
}

# Check a few:
# ASV 1, EU 10, meter = 80
check_10_1_90.df <- diffAbun_FUNGI_comboForest_dfLists[[1]][which(diffAbun_FUNGI_comboForest_dfLists[[1]]$Meter == 90),]
check_10_1_90.df[which(check_10_1_90.df$Transect=="B"),17] == fungi_forestMeansTransects[1,4] #okay, so this is getting assigned as expected!

# ASV 75, EU 53S, meter = 80
check_53S_75_90.df <- diffAbun_FUNGI_comboForest_dfLists[[326]][which(diffAbun_FUNGI_comboForest_dfLists[[326]]$Meter == 90),]
check_53S_75_90.df[which(check_53S_75_90.df$Transect=="L"),17] == fungi_forestMeansTransects[326,3+2] #getting assigned as expected!

# ASV 272, EU 53S, meter = 30 -- check to see if other meters are still looking as they should
check_53S_272_30.df <- diffAbun_FUNGI_comboForest_dfLists[[278]][which(diffAbun_FUNGI_comboForest_dfLists[[278]]$Meter == 30),]
check_53S_272_30.df$ASVabundance == diffAbun_no53N_tidy_FUNGI_forest_Lists[[278]] %>% filter(Meter == 30) %>% select(ASVabundance) #Yes!

# ASV 36, EU 53S, meter = 30 -- check to see if other meters are still looking as they should
check_54S_36_60.df <- diffAbun_FUNGI_comboForest_dfLists[[404]][which(diffAbun_FUNGI_comboForest_dfLists[[404]]$Meter == 60),]
check_54S_36_60.df$ASVabundance == diffAbun_no53N_tidy_FUNGI_forest_Lists[[404]] %>% filter(Meter == 60) %>% select(ASVabundance) #Yes!


######## PROKARYOTES ########
##### i. Get mean abundance of 90 and 100 meters for each ASV within each EU, for each transect
### i.i first, get mean abundance from 90-100 meters in forest.
# Pre-allocate dataframes
pro_forestMeansTransects <- as.data.frame(matrix(nrow=length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))*5, ncol=7)) #570 rows (115 ASVs x 5 EUs)
colnames(pro_forestMeansTransects) <- c("EU", "ASV_name", "ASVplusEU", "B_meanForestAbund", "L_meanForestAbund", "R_meanForestAbund", "T_meanForestAbund")
pro_df <- as.data.frame(matrix(nrow=12, ncol=5)) #12 rows, meters 90 and 100 x (up to) 4 transects
colnames(pro_df) <- c("ASV_name", "EU", "Transect", "Meter", "ASV_abundance")
proSubsetDats <- rep(list(pro_df), length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))*5) #need one for each ASV and EU combo, i.e. 570

# Make for loop that gets the mean from 90-100 in all of the elements of the list "diffAbun_no53N_tidy_FUNGI_forest_Lists"
for (j in 1:length(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists)){ #gets ASV name, EU, transect, and meter for each EU/ASV combo, but only for meters 90-100
  # The line below fills in as many rows as are needed, hence using the part with %in% twice
  proSubsetDats[[j]][1:nrow(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[j]][diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[j]]$Meter %in% c(90, 100),]),] <- 
    diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[j]][diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[j]]$Meter %in% c(90, 100),c(2,19:21, 17)]
  pro_forestMeansTransects[j,1] <- unique(na.omit(proSubsetDats[[j]]$EU))
  pro_forestMeansTransects[j,2] <- unique(na.omit(proSubsetDats[[j]]$ASV_name))
  pro_forestMeansTransects[j,3] <- paste(pro_forestMeansTransects[j,1], pro_forestMeansTransects[j,2], sep="_") #make this name match those of "diffAbun_no53N_tidy_FUNGI_forest_Lists"
  pro_forestMeansTransects[j,4] <- mean(proSubsetDats[[j]][which(proSubsetDats[[j]]$Transect == "B"),]$ASV_abundance)
  pro_forestMeansTransects[j,5] <- mean(proSubsetDats[[j]][which(proSubsetDats[[j]]$Transect == "L"),]$ASV_abundance)
  pro_forestMeansTransects[j,6] <- mean(proSubsetDats[[j]][which(proSubsetDats[[j]]$Transect == "R"),]$ASV_abundance)
  pro_forestMeansTransects[j,7] <- mean(proSubsetDats[[j]][which(proSubsetDats[[j]]$Transect == "T"),]$ASV_abundance)
}

# Checking a few! 
mean(proSubsetDats[[1]][which(proSubsetDats[[1]]$Transect == "T"),]$ASV_abundance) == pro_forestMeansTransects$T_meanForestAbund[1]
mean(proSubsetDats[[333]][which(proSubsetDats[[333]]$Transect == "B"),]$ASV_abundance) == pro_forestMeansTransects$B_meanForestAbund[333]
mean(proSubsetDats[[509]][which(proSubsetDats[[509]]$Transect == "R"),]$ASV_abundance) == pro_forestMeansTransects$R_meanForestAbund[509]

#############
### i. ii. Put mean abundance from 90-100 for each transect in place of separate 90 and 100 values
# Make all of these dataframes (not tibbles) so that they are compatible with some of the code below
diffAbun_PROKARYOTES_comboForest_dfLists <- vector("list", length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))*5) #make a empty list to store them in
for (i in 1:length(diffAbun_PROKARYOTES_comboForest_dfLists)){
  diffAbun_PROKARYOTES_comboForest_dfLists[[i]] <- as.data.frame(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[i]])
  print(class(diffAbun_PROKARYOTES_comboForest_dfLists[[i]])) #make sure that they are all data frames
}

# Make names equal to the original objects
names(diffAbun_PROKARYOTES_comboForest_dfLists) <- names(diffAbun_no53N_tidy_PROKARYOTES_forest_Lists)

# Replace 90-100 meter values in diffAbun_no53N_tidy_FUNGI_forest_Lists_forestCombo with fungi_forestMeans values. Will call it meter 90, for plotting purposes
# Preallocate a 114*5 long list of vector for each ASV in each EU to hold the indices for 90 and 100m
proIndicestoKeep <- vector("list", length(unique(diffAbunDat_no53N_tidy_PROKARYOTES_forest$ASV_name))*5)
# Outer for loop removes meters 90 and 100; inner for loop replaces meter 80 with mean of each transect across meters 80, 90, and 100
for (k in 1:length(diffAbun_PROKARYOTES_comboForest_dfLists)){ 
  proIndicestoKeep[[k]] <- which(diffAbun_PROKARYOTES_comboForest_dfLists[[k]]$Meter %in% c(10, 20, 30, 40, 50, 60, 70, 80, 90)==TRUE)
  diffAbun_PROKARYOTES_comboForest_dfLists[[k]] <- diffAbun_PROKARYOTES_comboForest_dfLists[[k]][proIndicestoKeep[[k]],] #keep only meters 10-90
  for (w in 1:4){ #loop over the four different transects
    diffAbun_PROKARYOTES_comboForest_dfLists[[k]][which(diffAbun_PROKARYOTES_comboForest_dfLists[[k]]$Transect %in% 
                                                    unique(sort(diffAbun_PROKARYOTES_comboForest_dfLists[[k]]$Transect))[w] & 
                                                    diffAbun_PROKARYOTES_comboForest_dfLists[[k]]$Meter == 90),17] <- #make 90 get the mean value for 90=100
      pro_forestMeansTransects[k,3+w] 
  }
}

# Check a few:
# ASV 1002, EU 10, meter = 80
check_10_1002_90.df <- diffAbun_PROKARYOTES_comboForest_dfLists[[1]][which(diffAbun_PROKARYOTES_comboForest_dfLists[[1]]$Meter == 90),]
check_10_1002_90.df[which(check_10_1002_90.df$Transect=="B"),17] == pro_forestMeansTransects[1,4] #okay, so this is getting assigned as expected!

# ASV 1987, EU 10, meter = 80
check_10_1987_90.df <- diffAbun_PROKARYOTES_comboForest_dfLists[[326]][which(diffAbun_PROKARYOTES_comboForest_dfLists[[326]]$Meter == 90),]
check_10_1987_90.df[which(check_10_1987_90.df$Transect=="L"),17] == pro_forestMeansTransects[326,3+2] #getting assigned as expected!

# ASV 3653, EU 54S, meter = 30 -- check to see if other meters are still looking as they should
check_54S_3653_30.df <- diffAbun_PROKARYOTES_comboForest_dfLists[[4000]][which(diffAbun_PROKARYOTES_comboForest_dfLists[[4000]]$Meter == 30),]
check_54S_3653_30.df$ASVabundance == diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[4000]] %>% filter(Meter == 30) %>% select(ASVabundance) #Yes!

# ASV 1568, EU 53S, meter = 30 -- check to see if other meters are still looking as they should
check_53S_1568_60.df <- diffAbun_PROKARYOTES_comboForest_dfLists[[2399]][which(diffAbun_PROKARYOTES_comboForest_dfLists[[2399]]$Meter == 60),]
check_53S_1568_60.df$ASVabundance == diffAbun_no53N_tidy_PROKARYOTES_forest_Lists[[2399]] %>% filter(Meter == 60) %>% select(ASVabundance) #Yes!

########################################
##### 3. SPLINES AND AUC CURVES
########################################
######## FUNGI ########
#############
### i. Plot the data compiled above (with averages at 90, no splines yet) 
diffAbun_FUNGI_comboForest_plots <- diffAbun_FUNGI_comboForest_dfLists %>% 
  map(~ggplot(.x, aes(x=Meter, y=ASVabundance)) + #plot everything in a list!
        geom_point(size=3, aes(color= Transect)) +
        theme_bw(base_size=16))

# Take a look at a few of them
diffAbun_FUNGI_comboForest_plots[[10]]
diffAbun_FUNGI_comboForest_plots[[533]]
diffAbun_FUNGI_comboForest_plots[[1]]
diffAbun_FUNGI_comboForest_plots[[444]]
diffAbun_FUNGI_comboForest_plots[[337]]
diffAbun_FUNGI_comboForest_plots[[339]]
diffAbun_FUNGI_comboForest_plots[[248]]
diffAbun_FUNGI_comboForest_plots[[187]]
diffAbun_FUNGI_comboForest_plots[[107]]

# Give all of these plots the same names as the (original) item that it came from
names(diffAbun_FUNGI_comboForest_plots) <- names(diffAbun_no53N_tidy_FUNGI_forest_Lists)
diffAbun_FUNGI_comboForest_plots[[1]]
diffAbun_FUNGI_comboForest_plots[[500]]

#############
# ii. Make splines for everything using ss function with small lambda to make error as small as possible  
diffAbundFungi_ssSplines <- vector("list", length(diffAbun_FUNGI_comboForest_dfLists)) #pre-allocate list
fungi_ssSplines <- as.data.frame(matrix(nrow=length(diffAbundFungi_ssSplines), ncol=4)) #570 rows (115 ASVs x 5 EUs) #pre-allocate
colnames(fungi_ssSplines) <- c("EU_ASV", "EU", "ASV_name", "adjRsq")
for (u in 1:length(diffAbundFungi_ssSplines)) {
  tryCatch({ #add an error handling function that does nothing so that the for loop will continue
    # even though some splines can't be fit (and thus have no adjusted r-squared)
    diffAbundFungi_ssSplines[[u]] <- ss(x= diffAbun_FUNGI_comboForest_dfLists[[u]]$Meter, y= diffAbun_FUNGI_comboForest_dfLists[[u]]$ASVabundance, all.knots= TRUE, lambda = 1e-15)
    fungi_ssSplines[u,1] <- names(diffAbun_FUNGI_comboForest_dfLists)[[u]] #make names the same as that ASV/EU combo
    fungi_ssSplines[u,2] <- unique(diffAbun_FUNGI_comboForest_dfLists[[u]]$EU)
    fungi_ssSplines[u,3] <- unique(diffAbun_FUNGI_comboForest_dfLists[[u]]$ASV_name)
    fungi_ssSplines[u,4] <- summary(diffAbundFungi_ssSplines[[u]])$adj.r.squared #make this the adj R -squared for the "model"
  }, error=function(e){})
}
# warning message:  In sqrt(sse/(n - df)) : NaNs produced 0 means that there were some NaNs in some step of the model fitting;
# warings message: cor(object$data$y, fitted) : the standard deviation is zero means that for some of these fits, 
# In the cases where there was no variation, it is likely that that ASV was not found in that given EU.

names(diffAbundFungi_ssSplines) <- names(diffAbun_FUNGI_comboForest_dfLists) #rename spline model lists

# Take a look at a random one
summary(diffAbundFungi_ssSplines[[50]])
fungi_ssSplines[50,] #oh, very bad r-squared

diffAbundFungi_ssSplines[[100]]
plot(diffAbundFungi_ssSplines[[50]])
points(x= diffAbundFungi_ssSplines[[50]]$data[,1], y= diffAbundFungi_ssSplines[[50]]$data[,2], pch=19)
points(x= diffAbun_FUNGI_comboForest_dfLists[[50]]$Meter, y= diffAbun_FUNGI_comboForest_dfLists[[50]]$ASVabundance, col=as.factor(diffAbun_FUNGI_comboForest_dfLists[[50]]$Transect), pch=19)


# I can't save all of these easily in a plot, and I can't make an object that saves the points added to the plots...
plotSplines <- function(ModelListElement){ #input is item in the list
  plot(ModelListElement)
  points(x= ModelListElement$data[,1], y=ModelListElement$data[,2], col=as.factor(ModelListElement$Transect), pch=19)
  main= paste(unique(ModelListElement$EU), unique(ModelListElement$ASV_name), sep="_")
  xlab= "Meter"
  ylab="ASV abundance"
}

plotSplines(diffAbundFungi_ssSplines[[100]])
# Okay, so to make this ggplot-able, would need to make a dataframe for each ASV/EU combination, all of which would be contained in a list.
# These datafarames would have 3 columns- x $data[,1] and $data[,2] from each element in diffAbundFungi_ssSplines, and diffAbun_FUNGI_comboForest_dfLists[[i]]$Transect
# combine this with

# shows where these data are in the ss output object
diffAbundFungi_ssSplines[[50]]$data[,2] == diffAbun_FUNGI_comboForest_dfLists[[50]]$ASVabundance

# iii. Explore the adjusted r-squared values
fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .5),] #These have a really high r-squared!
nrow(fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .5),])/nrow(fungi_ssSplines)*100 #only 7, or 1.754386% have an adj r squared greater than .5
fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .5),]
nrow(fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .3),]) #51
nrow(fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .3),])/nrow(fungi_ssSplines)*100 # About 9% have an adjusted r-squared above .3. This is equivalent to 51.
nrow(fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .25),]) #77
nrow(fungi_ssSplines[which(fungi_ssSplines$adjRsq >= .25),])/nrow(fungi_ssSplines)*100 # 13.5% have an adjusted r-squared above .25. This is equivalent to 77.

# A few example plots:
# 1. adjusted r-squared is 0.2781415
fungi_ssSplines$adjRsq[566]
plot(diffAbundFungi_ssSplines[[566]])
points(x= diffAbun_FUNGI_comboForest_dfLists[[566]]$Meter, y= diffAbun_FUNGI_comboForest_dfLists[[566]]$ASVabundance, col=as.factor(diffAbun_FUNGI_comboForest_dfLists[[566]]$Transect), pch=19)

# 2. r-squared is 0.4899863
fungi_ssSplines$adjRsq[533]
plot(diffAbundFungi_ssSplines[[533]])
points(x= diffAbun_FUNGI_comboForest_dfLists[[533]]$Meter, y= diffAbun_FUNGI_comboForest_dfLists[[533]]$ASVabundance, col=as.factor(diffAbun_FUNGI_comboForest_dfLists[[533]]$Transect), pch=19)
# THIS IS WEIRD, IT IS HIGHER IN THE PATCH!!!
# i. First plot the raw data and see how it looks-- OKAY higher in some but not all patches
diffAbun_FUNGI_comboForest_plots$EU_10_ASV_486 #one high in patch (low abundances)
diffAbun_FUNGI_comboForest_plots$EU_53S_ASV_486 #many higher in forest
diffAbun_FUNGI_comboForest_plots$EU_8_ASV_486 #higher in patch (many points but low abundance)
diffAbun_FUNGI_comboForest_plots$EU_54S_ASV_486 #many higher in forest (many points, high abundance)
diffAbun_FUNGI_comboForest_plots$EU_52_ASV_486 #a few higher in forest (low abundance)

# 3. 
fungi_ssSplines[which(fungi_ssSplines$EU_ASV == "EU_54S_ASV_486"),4] #adj r-squared is 0.624358
plot(diffAbundFungi_ssSplines$EU_54S_ASV_486)
points(x= diffAbun_FUNGI_comboForest_dfLists$EU_54S_ASV_486$Meter, y= diffAbun_FUNGI_comboForest_dfLists$EU_54S_ASV_486$ASVabundance, col=as.factor(diffAbun_FUNGI_comboForest_dfLists$EU_54S_ASV_486$Transect), pch=19)

# 4. 
fungi_ssSplines[which(fungi_ssSplines$EU_ASV == "EU_8_ASV_17"),4] #adj r-squared is 0.8833391
plot(diffAbundFungi_ssSplines$EU_8_ASV_17)
points(x=diffAbun_FUNGI_comboForest_dfLists$EU_8_ASV_17$Meter, y= diffAbun_FUNGI_comboForest_dfLists$EU_8_ASV_17$ASVabundance, col=as.factor(diffAbun_FUNGI_comboForest_dfLists$EU_8_ASV_17$Transect), pch=19)

# vii. Histogram of all the r-squared
hist(fungi_ssSplines$adjRsq)

# 4. Filter out ASV/EU combos that have an adjusted r-squared value less than 0.25 and those that do not have at 
# least 100 reads of the ASV present (across all samples in that EU)
# i. Apply adjusted p-value threshold
fungi_ToKeepRsqIndex <- fungi_ssSplines[which(fungi_ssSplines$adjRsq >= 0.25),1] #get only those ASV/EU combos with an adj. r-squared at least 0.25
length(fungi_ToKeepRsqIndex) #77

fungi_ASV_EUsListsToKeep1 <- which(names(diffAbun_FUNGI_comboForest_dfLists) %in% fungi_ToKeepRsqIndex)
fungi_listSubset1 <- (diffAbun_FUNGI_comboForest_dfLists)[fungi_ASV_EUsListsToKeep1] #get a subset of the lists that pass the .25 criterion
length(fungi_listSubset1) #77, as expected
names(fungi_listSubset1) == fungi_ToKeepRsqIndex #yes, this worked!

# ii. Apply 100 read number threshold (across meters 10-90 for everything)
fungi_sumsDF <- as.data.frame(matrix(nrow=length(fungi_listSubset1), ncol=2))
colnames(fungi_sumsDF) <- c("ASV_EUname", "totalASVabundance")
for (i in 1:nrow(fungi_sumsDF)){
  fungi_sumsDF[i,1] <- names(fungi_listSubset1)[[i]] #get names for each
  fungi_sumsDF[i,2] <- sum(fungi_listSubset1[[i]]$ASVabundance) #get ASV abundance for each 
}
fungi_ASV_EUsListsToKeep2 <- fungi_sumsDF[which(fungi_sumsDF$totalASVabundance >= 100),1] #names of the ASVs/EUs that have at least 100 reads
length(fungi_ASV_EUsListsToKeep2) #65 passed this filter

fungi_goodSplinesIndex <- which(fungi_ssSplines$EU_ASV %in% fungi_ASV_EUsListsToKeep2 == TRUE)
fungi_goodSplines <- fungi_ssSplines[fungi_goodSplinesIndex,]
nrow(fungi_goodSplines) #65 as expected!

# Select only those in ASV_EUsListsToKeep2, i.e. a subset of the lists that pass the .25 criterion AND the abundance threshold
fungalGoodSplinesAllDatList <- (diffAbun_FUNGI_comboForest_dfLists)[fungi_ASV_EUsListsToKeep2] #65
names(fungalGoodSplinesAllDatList)

# iii. Get AUC values
# First, pre-allocate a list of dataframes to hold everything
df <- as.data.frame(matrix(nrow=7, ncol=8)) #12 rows (meters 20-90)
colnames(df) <- c("ASV_name", "EU", "EU_ASV", "Meter", "AUC", "numberPoints", "i_index")
aucList_fungi <- rep(list(df), length(fungalGoodSplinesAllDatList)) #need one for each good ASV and EU combo, i.e. 65
fungi_removalIndices <- vector("list", length(fungalGoodSplinesAllDatList))

for (g in 1:length(aucList_fungi)){ #loop over all (remaining) combinations of ASV and EU
  for (i in 1:8){ #loop over 8, i.e. the number of meters 
    aucList_fungi[[g]][(i), 1] <- unique(diffAbun_FUNGI_comboForest_dfLists[[g]]$ASV_name) #add in ASV information
    aucList_fungi[[g]][(i), 2] <- unique(diffAbun_FUNGI_comboForest_dfLists[[g]]$EU) #add in EU information
    aucList_fungi[[g]][(i), 3] <- paste(unique(diffAbun_FUNGI_comboForest_dfLists[[g]]$EU), unique(diffAbun_FUNGI_comboForest_dfLists[[g]]$ASV_name), sep="_") #add in EU/ASV information
    aucList_fungi[[g]][(i), 4] <- (i*10)
    aucList_fungi[[g]][(i+1), 5] <- auc(sort(unique(diffAbundFungi_ssSplines[[g]]$data[,1]))[1:(i+1)], diffAbundFungi_ssSplines[[g]]$y[1:(i+1)])
    aucList_fungi[[g]][(i), 6] <- length(1:i+1)
    aucList_fungi[[g]][(i+1), 7] <- i
  }
  fungi_removalIndices[[g]] <- c(which(is.na(aucList_fungi[[g]]$i_index)), which(is.na(aucList_fungi[[g]]$Meter))) #clean up by removing rows with NAs 
  aucList_fungi[[g]] <- aucList_fungi[[g]][-fungi_removalIndices[[g]],]
}

# This is an example of how to find the meter value where the AUC is closest to 25% that of the maximum
aucList_fungi[[1]][which.min(abs(aucList_fungi[[1]]$AUC-max(aucList_fungi[[1]]$AUC)*0.25)),4]

# For loop that finds the meter value where the AUC is closest to 25% that of the maximum
aucList_fungi2 <- aucList_fungi
for (k in 1:length(aucList_fungi2)) {
  aucList_fungi2[[k]] <- aucList_fungi[[k]] %>% 
    select(ASV_name:AUC) %>% 
    mutate(MaxAUC = max(aucList_fungi[[k]]$AUC)) %>% 
    mutate(quarterMaxAUC = MaxAUC*0.25) %>% 
    mutate(quarterMeter = aucList_fungi[[k]][which.min(abs(AUC-quarterMaxAUC)),4])
}

aucList_fungi2[[56]]

fungalGoodSplinesListShorter <- fungalGoodSplinesList #make another list to shorten and clean, below
# Now, make the fungalGoodSplinesList smaller and tidier
for (j in 1:length(fungalGoodSplinesAllDatList)) {
  fungalGoodSplinesListShorter[[j]] <- fungalGoodSplinesAllDatList[[j]] %>% 
    select(c(SampleNumberID, ASV_name, Kingdom:Species, ASVabundance, EU, Transect, Meter)) %>% 
    mutate(adjR2 = fungi_goodSplines$adjRsq[j]) %>% #give them their r-squared values
    mutate(AUCmeter = unique(aucList_fungi2[[j]]$quarterMeter))
}

#View(fungalGoodSplinesListShorter[[55]])
length(fungalGoodSplinesListShorter) #65

# Make a dataframe that is conducive for plotting in ggplot2
fungalGoodSplinesForPlot <- as.data.frame(matrix(nrow=length(fungalGoodSplinesListShorter), ncol= ncol(fungalGoodSplinesListShorter[[1]])))
dim(fungalGoodSplinesForPlot)
colnames(fungalGoodSplinesForPlot) <- colnames(fungalGoodSplinesListShorter[[1]])
for (w in 1:nrow(fungalGoodSplinesForPlot)){
  fungalGoodSplinesForPlot[w,] <- fungalGoodSplinesListShorter[[w]][1,]
}
#View(fungalGoodSplinesForPlot)

# This shows that there are 45 unique ASVs here
length(unique(fungalGoodSplinesForPlot$ASV_name)) #45 unique ASVs

# ggplot it!
fungalAUCmetersPlot1 <- ggplot(fungalGoodSplinesForPlot, aes(fill=Phylum, x=AUCmeter)) + 
  geom_bar(position="stack", stat="count")
quartz()
fungalAUCmetersPlot1

######## PROKARYOTES ########
#############
### i. Plot the data compiled above (with averages at 90, no splines yet) 
