# FungalAbundanceDropOff.R
# Started October 6, 2022

# This script looks at how the abundance of different ASVs changes across the transect. To do so, it uses area 
# under the curve with splines to evaluate where the abundance drops off. Specifically, within each ASV 
# within each EU and for each transect , it gets the mean abundance from 80 meters to 100 meters, designating 
# this the "forest" abundance (i.e. ASV 1, EU 10, meter 80 would have 4 values-- the mean ASV abundance of transect 
# B at 80, 90, and 100 m, L's mean across 80,90, and 100, R's mean across 80,90, and 100, and T's mean across 80,90,
# and 100 m . The rest of the abundances along the transect are left as is. Then a spline is fit to these points. 
# The area under the curve is then calculated for the area under the spline line from 10 meters to the averaged forest
# abundance. I then find the point along the x-axis where the area under the curve is half of this "maximum" area.

################################
# SET UP
################################
options(scipen = 999) #make it so that annoying scientific notation doesn't pop up

library("tidyverse")
library("splines")
library("gridExtra")
library(phyloseq)
library(purrr)

# setwd("/Volumes/Memorex\ USB/") #for working on a different computer only:
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40. All EUs, including EU 53N
# R OBJECT MADE IN: UbiquitypostUbiqSetup.R
load("RobjectsSaved/diffAbunDat_no53N_tidy_FUNGI") #made in DiffAbundFungi

#View(diffAbunDat_no53N_tidy_FUNGI)

################################
# MAIN PART OF SCRIPT

####################
##### 1. ISOLATE ONLY ASVS THAT WERE MORE ABUNDANT IN THE FOREST #####
####################
nrow(diffAbunDat_no53N_tidy_FUNGI) #number of rows is equal to the 5 EUs times 10 points along the transect times 114 differentially abundant fungal taxa
colnames(diffAbunDat_no53N_tidy_FUNGI)
diffAbunDat_no53N_tidy_FUNGI_forest <- diffAbunDat_no53N_tidy_FUNGI %>% filter(Habitat== "forest")
# View(diffAbunDat_no53N_tidy_FUNGI_forest) # ASVabundance is abundance at each site

####################
##### 2. MAKE DATAFRAMES FOR EACH OF THE ASVS IN EACH EU, (so 5 EUs x 114 ASVs = 570 dataframes) #####
####################
diffAbun_no53N_tidy_FUNGI_forest_Lists <- diffAbunDat_no53N_tidy_FUNGI_forest %>% 
  group_split(EU, ASV_name)
length(diffAbun_no53N_tidy_FUNGI_forest_Lists) #5 EUs x 114 ASVs = 570 dataframes
# Investigate how this broke it up to then check names:
orderedASVs <- vector(length=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
orderedEUs <- vector(length=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
for (i in 1:length(diffAbun_no53N_tidy_FUNGI_forest_Lists)){
  orderedASVs[i] <- print(diffAbun_no53N_tidy_FUNGI_forest_Lists[[i]]$ASV_name)
  orderedEUs[i] <- print(diffAbun_no53N_tidy_FUNGI_forest_Lists[[i]]$EU)
}
orderedASVs #order is 1, 10, 1060, 1069, 1081, 109, 110, 1127, 113, 115
orderedEUs #order is 10, 52., 53S, 54S, 8

# Rename names of items in list based on combining EU and ASV information
# Make a vector for naming the lists
mergedNames <- vector(length=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)
for (i in 1:length(orderedASVs)){
  mergedNames[i] <- paste(orderedEUs[i], orderedASVs[i], sep="_")
}
names(diffAbun_no53N_tidy_FUNGI_forest_Lists) <- mergedNames #give names
# Check a few (they all look as expected!)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_10_ASV_1203$ASV_name)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_10_ASV_1203$EU)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_53S_ASV_591$ASV_name)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_53S_ASV_591$EU)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_8_ASV_97$ASV_name)
unique(diffAbun_no53N_tidy_FUNGI_forest_Lists$EU_8_ASV_97$EU)

####################
##### 3. FOR EACH ASV IN EACH EU, GET MEAN ABUNDANCE FROM 80 -100 METERS TO REPRESENT "FOREST ABUNDANCE" #####
# (i.e. the mean of abundances at 80m, 90m, and 100m (up to 12 points))
####################
# DID NOT END UP USING THE METHOD IN 3, BUT KEEPING ANYWAY JUST IN CASE

diffAbun_no53N_tidy_FUNGI_forest_Lists

# i. pre-allocate dataframes
fungi_forestMeans <- as.data.frame(matrix(nrow=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5, ncol=5)) #570 rows (115 ASVs x 5 EUs)
colnames(fungi_forestMeans) <- c("EU", "ASV_name", "ASVplusEU", "meanForestAbund", "halfMeanForestAbund")
df <- as.data.frame(matrix(nrow=12, ncol=5)) #12 rows, meters 80, 90, and 100 x (up to) 4 transects
colnames(df) <- c("ASV_name", "EU", "Transect", "Meter", "ASV_abundance")
subsetDats <- rep(list(df), length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5) #need one for each ASV and EU combo, i.e. 570

# ii. Make for loop that gets the mean from 80-100 in all of the elements of the list "diffAbun_no53N_tidy_FUNGI_forest_Lists"
for (j in 1:length(diffAbun_no53N_tidy_FUNGI_forest_Lists)){ #gets ASV name, EU, transect, and meter for each EU/ASV combo, but only for meters 80-100
  # The line below fills in as many rows as are needed, hence using the part with %in% twice
  subsetDats[[j]][1:nrow(diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]][diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]]$Meter %in% c(80, 90, 100),]),] <- 
    diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]][diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]]$Meter %in% c(80, 90, 100),c(2,19:21, 17)]
  fungi_forestMeans[j,1] <- unique(na.omit(subsetDats[[j]]$EU))
  fungi_forestMeans[j,2] <- unique(na.omit(subsetDats[[j]]$ASV_name))
  fungi_forestMeans[j,3] <- paste(fungi_forestMeans[j,1], fungi_forestMeans[j,2], sep="_") #make this name match those of "diffAbun_no53N_tidy_FUNGI_forest_Lists"
  fungi_forestMeans[j,4] <- mean(subsetDats[[j]]$ASV_abundance, na.rm = TRUE)
  fungi_forestMeans[j,5] <- (mean(subsetDats[[j]]$ASV_abundance, na.rm = TRUE))*0.5
}

# Check a couple from above to make sure for loop is working correctly (using different method than %in%):
# EU 10, ASV 1
abund80_10_1 <- diffAbun_no53N_tidy_FUNGI_forest_Lists[[1]]$ASVabundance[which(diffAbun_no53N_tidy_FUNGI_forest_Lists[[1]]$Meter == 80)] #all abundances at 80m
abund90_10_1 <- diffAbun_no53N_tidy_FUNGI_forest_Lists[[1]]$ASVabundance[which(diffAbun_no53N_tidy_FUNGI_forest_Lists[[1]]$Meter == 90)] #all abundances at 90m
abund100_10_1 <- diffAbun_no53N_tidy_FUNGI_forest_Lists[[1]]$ASVabundance[which(diffAbun_no53N_tidy_FUNGI_forest_Lists[[1]]$Meter == 100)] #all abundances at 100m
(sum(abund80_10_1, abund90_10_1, abund100_10_1)/(length(abund80_10_1) + length(abund90_10_1) + length(abund100_10_1))) == fungi_forestMeans[1,]$meanForestAbund #TRUE

# EU 52, ASV 59
abund80_52_59 <- diffAbun_no53N_tidy_FUNGI_forest_Lists[[200]]$ASVabundance[which(diffAbun_no53N_tidy_FUNGI_forest_Lists[[200]]$Meter == 80)] #all abundances at 80m
abund90_52_59 <- diffAbun_no53N_tidy_FUNGI_forest_Lists[[200]]$ASVabundance[which(diffAbun_no53N_tidy_FUNGI_forest_Lists[[200]]$Meter == 90)] #all abundances at 90m
abund100_52_59 <- diffAbun_no53N_tidy_FUNGI_forest_Lists[[200]]$ASVabundance[which(diffAbun_no53N_tidy_FUNGI_forest_Lists[[200]]$Meter == 100)] #all abundances at 100m
(sum(abund80_52_59, abund90_52_59, abund100_52_59)/(length(abund80_52_59) + length(abund90_52_59) + length(abund100_52_59))) == fungi_forestMeans[200,]$meanForestAbund #TRUE

####################
##### 4. (WITHIN EACH TRANSECT) and FOR EACH ASV IN EACH EU, GET MEAN ABUNDANCE FROM 80 -100 METERS TO REPRESENT "FOREST ABUNDANCE" #####
# (i.e. 4 values-- the mean of abundances at 80m, 90m, and 100m for T, B, L, and R transects)
####################
# diffAbun_no53N_tidy_FUNGI_forest_Lists

# i. pre-allocate dataframes
fungi_forestMeansTransects <- as.data.frame(matrix(nrow=length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5, ncol=7)) #570 rows (115 ASVs x 5 EUs)
colnames(fungi_forestMeansTransects) <- c("EU", "ASV_name", "ASVplusEU", "B_meanForestAbund", "L_meanForestAbund", "R_meanForestAbund", "T_meanForestAbund")
df <- as.data.frame(matrix(nrow=12, ncol=5)) #12 rows, meters 80, 90, and 100 x (up to) 4 transects
colnames(df) <- c("ASV_name", "EU", "Transect", "Meter", "ASV_abundance")
subsetDats <- rep(list(df), length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5) #need one for each ASV and EU combo, i.e. 570

# ii. Make for loop that gets the mean from 80-100 in all of the elements of the list "diffAbun_no53N_tidy_FUNGI_forest_Lists"
for (j in 1:length(diffAbun_no53N_tidy_FUNGI_forest_Lists)){ #gets ASV name, EU, transect, and meter for each EU/ASV combo, but only for meters 80-100
  # The line below fills in as many rows as are needed, hence using the part with %in% twice
  subsetDats[[j]][1:nrow(diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]][diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]]$Meter %in% c(80, 90, 100),]),] <- 
    diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]][diffAbun_no53N_tidy_FUNGI_forest_Lists[[j]]$Meter %in% c(80, 90, 100),c(2,19:21, 17)]
  fungi_forestMeansTransects[j,1] <- unique(na.omit(subsetDats[[j]]$EU))
  fungi_forestMeansTransects[j,2] <- unique(na.omit(subsetDats[[j]]$ASV_name))
  fungi_forestMeansTransects[j,3] <- paste(fungi_forestMeans[j,1], fungi_forestMeans[j,2], sep="_") #make this name match those of "diffAbun_no53N_tidy_FUNGI_forest_Lists"
  fungi_forestMeansTransects[j,4] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "B"),]$ASV_abundance)
  fungi_forestMeansTransects[j,5] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "L"),]$ASV_abundance)
  fungi_forestMeansTransects[j,6] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "R"),]$ASV_abundance)
  fungi_forestMeansTransects[j,7] <- mean(subsetDats[[j]][which(subsetDats[[j]]$Transect == "T"),]$ASV_abundance)
}

# Checking a few! (subsetDats is correct based on the checking I did in #4)
mean(subsetDats[[1]][which(subsetDats[[1]]$Transect == "T"),]$ASV_abundance) == fungi_forestMeansTransects$T_meanForestAbund[1]
mean(subsetDats[[333]][which(subsetDats[[333]]$Transect == "B"),]$ASV_abundance) == fungi_forestMeansTransects$B_meanForestAbund[333]
mean(subsetDats[[509]][which(subsetDats[[509]]$Transect == "R"),]$ASV_abundance) == fungi_forestMeansTransects$R_meanForestAbund[509]

####################
##### 5. PUT MEAN ABUNDANCE AT 80 TO 100 METERS ("FOREST ABUNDANCE") FOR EACH TRANSECT IN PLACE OF 80-100m #####
####################
# i. Make all of these dataframes (not tibbles) so that they are compatible with some of the code below
diffAbun_FUNGI_comboForest_dfLists <- vector("list", length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5) #make a empty list to store them in
for (i in 1:length(diffAbun_FUNGI_comboForest_dfLists)){
  diffAbun_FUNGI_comboForest_dfLists[[i]] <- as.data.frame(diffAbun_no53N_tidy_FUNGI_forest_Lists[[i]])
  print(class(diffAbun_FUNGI_comboForest_dfLists[[i]])) #make sure that they are all data frames
}

# ii. Replace 80-100 meter values in diffAbun_no53N_tidy_FUNGI_forest_Lists_forestCombo with fungi_forestMeans values. Will make it meter 80, for plotting purposes
head(diffAbun_FUNGI_comboForest_dfLists[[1]])
# Preallocate a 114*5 long list of vector for each ASV in each EU to hold the indices for 90 and 100m
indicestoKeep <- vector("list", length(unique(diffAbunDat_no53N_tidy_FUNGI_forest$ASV_name))*5)

# Outer for loop removes meters 90 and 100; inner for loop replaces meter 80 with mean of each transect across meters 80, 90, and 100
for (k in 1:length(diffAbun_FUNGI_comboForest_dfLists)){ 
  indicestoKeep[[k]] <- which(diffAbun_FUNGI_comboForest_dfLists[[k]]$Meter %in% c(10, 20, 30, 40, 50, 60, 70, 80)==TRUE)
  diffAbun_FUNGI_comboForest_dfLists[[k]] <- diffAbun_FUNGI_comboForest_dfLists[[k]][indicestoKeep[[k]],] #keep only meters 10-80
  for (w in 1:4){ #loop over the four different transects
    diffAbun_FUNGI_comboForest_dfLists[[k]][which(diffAbun_FUNGI_comboForest_dfLists[[k]]$Transect %in% 
                                                                                   unique(sort(diffAbun_FUNGI_comboForest_dfLists[[k]]$Transect))[w] & 
                                                    diffAbun_FUNGI_comboForest_dfLists[[k]]$Meter == 80),17] <-
      fungi_forestMeansTransects[k,3+w] 
  }
}

# iii. Now, check a few of these:
# ASV 1, EU 10, meter = 80
check_10_1_80.df <- diffAbun_FUNGI_comboForest_dfLists[[1]][which(diffAbun_FUNGI_comboForest_dfLists[[1]]$Meter == 80),]
check_10_1_80.df[which(check_10_1_80.df$Transect=="B"),17] == fungi_forestMeansTransects[1,4] #okay, so this is getting assigned as expected!

# ASV 75, EU 53S, meter = 80
check_53S_75_80.df <- diffAbun_FUNGI_comboForest_dfLists[[326]][which(diffAbun_FUNGI_comboForest_dfLists[[326]]$Meter == 80),]
check_53S_75_80.df[which(check_53S_75_80.df$Transect=="L"),17] == fungi_forestMeansTransects[326,3+2] #okay, so this is getting assigned as expected!
 
# ASV 272, EU 53S, meter = 30 -- check to see if other meters are still looking as they should
check_53S_272_30.df <- diffAbun_FUNGI_comboForest_dfLists[[278]][which(diffAbun_FUNGI_comboForest_dfLists[[278]]$Meter == 30),]
check_53S_272_30.df$ASVabundance == diffAbun_no53N_tidy_FUNGI_forest_Lists[[278]] %>% filter(Meter == 30) %>% select(ASVabundance) #Yes!

# ASV 36, EU 53S, meter = 30 -- check to see if other meters are still looking as they should
check_54S_36_60.df <- diffAbun_FUNGI_comboForest_dfLists[[404]][which(diffAbun_FUNGI_comboForest_dfLists[[404]]$Meter == 60),]
check_54S_36_60.df$ASVabundance == diffAbun_no53N_tidy_FUNGI_forest_Lists[[404]] %>% filter(Meter == 60) %>% select(ASVabundance) #Yes!

####################
##### 6. MAKE SPLINES FOR THE DATA IN EACH ONE OF THESE DATAFRAMES ABOVE #####
####################
diffAbun_FUNGI_comboForest_dfLists

# i. First plot the raw data and see how it looks 
diffAbun_FUNGI_comboForest_plots <- diffAbun_FUNGI_comboForest_dfLists %>% 
  map(~ggplot(.x, aes(x=Meter, y=ASVabundance)) + #plot everything in a list!
        geom_point(size=3, aes(color= Transect)) +
        theme_bw(16))

# Give all of these plots the same names as the (original) item that it came from
names(diffAbun_FUNGI_comboForest_plots) <- names(diffAbun_no53N_tidy_FUNGI_forest_Lists)
diffAbun_FUNGI_comboForest_plots[[1]]
diffAbun_FUNGI_comboForest_plots[[500]]
diffAbun_FUNGI_comboForest_plots[[145]]
length(diffAbun_FUNGI_comboForest_plots)

# THIS IS WHERE I STOPPED!
# ii. Generate a model with splines for each species
diffAbun_FUNGI_comboForest_SPLINEs <- diffAbun_FUNGI_comboForest_dfLists %>% 
  map(~lm(ASVabundance ~ bs(Meter), data= .x)) #bs is a B-spline basis for polynomial splines

# Plot model predictions with geom_line
# This makes a new series of xs (meters) to get spline lines on, from 10 m (farthest point in patch) to 80 m
# (which represents farthest in the forest, but is average across 80,90, 100 m
xDatForPreds <- data.frame(Meter=seq(10, 80, length=71)) 

fungalDA_subset_model_preds <- map_df(diffAbun_FUNGI_comboForest_SPLINEs, ~cbind(xDatForPreds, ASVabundance=predict(.x, xDatForPreds))) 
View(fungalDA_subset_model_preds)


# Plot these predictions
fungalSpline_subset_plot <- ggplot(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset, aes(Meter, abundZ_score, colour=ASV_name)) +
  # Show raw data
  geom_point() +
  # Compare internal geom_smooth calculations with our models
  geom_smooth(size=1.5, linetype="11", se=TRUE, formula=y ~ bs(x), method="lm") +
  # Plot model predictions
  geom_line(data=fungalDA_subset_model_preds) +
  theme_classic()






  
  
  
  map(~ggplot(.x, aes(x=Meter, y=abundZ_score)) + #plot everything in a list!
        geom_point(size=3, aes(color= EU)) +
        theme_bw(16) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE))


# iii. Method 3 with just a few of the ASVs (first 10 ASVs, across all EUs, I think, i.e. 60 here)
zScores_no53N_FUNGI_Plus1.5_forestCombo_subset <- zScores_no53N_FUNGI_Plus1.5_forestCombo[zScores_no53N_FUNGI_Plus1.5_forestCombo$ASV_name %in% 
                                                                                            unique(zScores_no53N_FUNGI_Plus1.5_forestCombo$ASV_name)[1:10],]
# Generate a model with splines for each species
fungalDAbyASV_models_subset <- zScores_no53N_FUNGI_Plus1.5_forestCombo_subset %>% 
  split(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset$ASV_name) %>% 
  map(~lm(abundZ_score ~ bs(Meter), data= .x)) #bs is a B-spline basis for polynomial splines

# Plot model predictions with geom_line
# length = 20 makes it so 20 x values, I think. Max Meter is 80, i.e. whole forest
new_dat_subset <- data.frame(Meter=seq(min(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset$Meter), max(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset$Meter), length=20)) 
fungalDA_subset_model_preds <- map_df(fungalDAbyASV_models_subset, ~cbind(new_dat_subset, abundZ_score=predict(.x, new_dat_subset)),
                                      .id="ASV_name") 

# Plot these predictions
fungalSpline_subset_plot <- ggplot(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset, aes(Meter, abundZ_score, colour=ASV_name)) +
  # Show raw data
  geom_point() +
  # Compare internal geom_smooth calculations with our models
  geom_smooth(size=1.5, linetype="11", se=TRUE, formula=y ~ bs(x), method="lm") +
  # Plot model predictions
  geom_line(data=fungalDA_subset_model_preds) +
  theme_classic()























# iv. Divide up zScores_no53N_FUNGI_Plus1.5_forestCombo into 114 lists, one per ASV, then plot
zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV <- zScores_no53N_FUNGI_Plus1.5_forestCombo %>% 
  group_split(ASV_name)
names(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV) <- sort(unique(zScores_no53N_FUNGI_Plus1.5_forestCombo$ASV_name)) #rename to make ASV names
# View(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV[[1]]) #each one has 40 rows, i.e. 5 EUs x values for 10-70, plus 80 (i.e. combined forest))
zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_109

# v. Plot them all together, in a list of plots using map function (purrr package from tidyverse)
# 1. This gets a line using a linear model, i.e. lm
zScores_no53N_FUNGI_forestCombo_plots1 <- zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV %>% 
  map(~ggplot(.x, aes(x=Meter, y=abundZ_score)) + #plot everything in a list!
        geom_point(size=3, aes(color= EU)) +
        theme_bw(16) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE))

names(zScores_no53N_FUNGI_forestCombo_plots1) <- names(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV)
zScores_no53N_FUNGI_forestCombo_plots1$ASV_819
#View(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819)

# 2. This gets one line for everything using "loess", i.e. local regression fitting
zScores_no53N_FUNGI_forestCombo_plots2 <- zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV %>% 
  map(~ggplot(.x, aes(x=Meter, y=abundZ_score)) + #plot everything in a list!
        geom_point(size=3, aes(color= EU)) +
        theme_bw(16) +
        geom_smooth(method = "loess"))

names(zScores_no53N_FUNGI_forestCombo_plots2) <- names(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV)
zScores_no53N_FUNGI_forestCombo_plots2$ASV_819

# 3. A polynomial interpolation (looks the same as lm above)
zScores_no53N_FUNGI_forestCombo_plots3 <- zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV %>% 
  map(~ggplot(.x, aes(x=Meter, y=abundZ_score)) + #plot everything in a list!
        geom_point(size=3, aes(color= EU)) +
        theme_bw(16) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE))

names(zScores_no53N_FUNGI_forestCombo_plots3) <- names(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV)
zScores_no53N_FUNGI_forestCombo_plots3$ASV_819

# 4. Splines. Don't know yet how to do it for every ASV at once, so here's just ASV_819
# i. Method 1:
spline.819 <- as.data.frame(spline(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819$Meter, zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819$abundZ_score))
p819 <- ggplot(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819, aes(Meter, abundZ_score)) +
  geom_point(size=3, aes(color= EU))
p819 + geom_line(data = spline.819, aes(x = x, y = y))







################### THIS WAY USES Z-SCORES AND IS NO LONGER THE APPROACH I WANT TO TAKE ########

GI_forestCombo_plots3$ASV_819

# 4. Splines. Don't know yet how to do it for every ASV at once, so here's just ASV_819
# i. Method 1:
spline.819 <- as.data.frame(spline(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819$Meter, zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819$abundZ_score))
p819 <- ggplot(zScores_no53N_FUNGI_Plus1.5_forestCombo_ByASV$ASV_819, aes(Meter, abundZ_score)) +
  geom_point(size=3, aes(color= EU))
p819 + geom_line(data = spline.819, aes(x = x, y = y))

# ii. Method 2 (with help from https://community.rstudio.com/t/cannot-plot-spline-regression-with-ggplot/31792/2)
# Generate a model with splines for each species
fungalDAbyASV_models <- zScores_no53N_FUNGI_Plus1.5_forestCombo %>% 
  split(zScores_no53N_FUNGI_Plus1.5_forestCombo$ASV_name) %>% 
  map(~lm(abundZ_score ~ bs(Meter), data= .x)) #bs is a B-spline basis for polynomial splines

# Plot model predictions with geom_line
# length = 20 makes it so 20 x values, I think. Max Meter is 80, i.e. whole forest
new_dat <- data.frame(Meter=seq(min(zScores_no53N_FUNGI_Plus1.5_forestCombo$Meter), max(zScores_no53N_FUNGI_Plus1.5_forestCombo$Meter), length=20)) 
fungalDA_model_preds <- map_df(fungalDAbyASV_models, ~cbind(new_dat, abundZ_score=predict(.x, new_dat)),
                               .id="ASV_name") 

# Plot all of the predictions
fungalSpline_plot <- ggplot(zScores_no53N_FUNGI_Plus1.5_forestCombo, aes(Meter, abundZ_score, colour=ASV_name)) +
  # Show raw data
  geom_point() +
  # Compare internal geom_smooth calculations with our models
  geom_smooth(size=1.5, linetype="11", se=FALSE, formula=y ~ bs(x), method="lm") +
  # Plot model predictions
  geom_line(data=fungalDA_model_preds) +
  theme_classic()
quartz()
fungalSpline_plot

# iii. Method 3 with just a few of the ASVs
zScores_no53N_FUNGI_Plus1.5_forestCombo_subset <- zScores_no53N_FUNGI_Plus1.5_forestCombo[zScores_no53N_FUNGI_Plus1.5_forestCombo$ASV_name %in% 
                                                                                            unique(zScores_no53N_FUNGI_Plus1.5_forestCombo$ASV_name)[1:10],]
# Generate a model with splines for each species
fungalDAbyASV_models_subset <- zScores_no53N_FUNGI_Plus1.5_forestCombo_subset %>% 
  split(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset$ASV_name) %>% 
  map(~lm(abundZ_score ~ bs(Meter), data= .x)) #bs is a B-spline basis for polynomial splines

# Plot model predictions with geom_line
# length = 20 makes it so 20 x values, I think. Max Meter is 80, i.e. whole forest
new_dat_subset <- data.frame(Meter=seq(min(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset$Meter), max(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset$Meter), length=20)) 
fungalDA_subset_model_preds <- map_df(fungalDAbyASV_models_subset, ~cbind(new_dat_subset, abundZ_score=predict(.x, new_dat_subset)),
                                      .id="ASV_name") 

# Plot these predictions
fungalSpline_subset_plot <- ggplot(zScores_no53N_FUNGI_Plus1.5_forestCombo_subset, aes(Meter, abundZ_score, colour=ASV_name)) +
  # Show raw data
  geom_point() +
  # Compare internal geom_smooth calculations with our models
  geom_smooth(size=1.5, linetype="11", se=TRUE, formula=y ~ bs(x), method="lm") +
  # Plot model predictions
  geom_line(data=fungalDA_subset_model_preds) +
  theme_classic()
