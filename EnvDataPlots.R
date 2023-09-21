# Environmental and physicochemical data and plots
# December 8, 2021

# This script creates plots for soil physiochemical and environmental data.

setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot


# FUNCTIONS DEFINED IN THIS SCRIPT:
#1.  Write function to get mean pH
Mean_pH_func <- function(pHs_mat) { #enter a dataframe that has the pHs, as made above
  # this function works for averaging 2-3 pHs (per row), but could easily be adjusted in
  # for loop to be more flexible.
  H_mat <- matrix(nrow=nrow(pHs_mat), ncol=ncol(pHs_mat))
  colnames(H_mat) <- colnames(pHs_mat)
  for(i in 1:nrow(pHs_mat)){ #get the H+ concentration for each pH reading
    H_mat[i,1] <- 10^(pHs_mat[i,1]*-1)
    H_mat[i,2] <- 10^(pHs_mat[i,2]*-1)
    H_mat[i,3] <- 10^(pHs_mat[i,3]*-1)
    H_mat[i,4] <- 10^(pHs_mat[i,4]*-1) #just added this for maybe 4 transects??
  }
  avgH <- rowMeans(H_mat, na.rm=TRUE) #get the average H+ concentration for each set of pHs
  avgpH <- -log10(avgH) #convert back to pH
  return(avgpH)
}

# 2. This function takes a vector of different pHs as input and returns a single mean pH (so it's
# less complicated than Mean_pH_func)
mean_pH_vecFunc <- function(pHs_vector){  
  meanHydroniumConc <- mean(10^(pHs_vector*-1)) #gets mean hydronium concentrations of the pHs in vector
  mean_pH <- -log10(meanHydroniumConc) #converts this back to pH by taking negative log of hydronium concentrations
  return(mean_pH)
} 

# 3. metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

############################################################

load("RobjectsSaved/medianEU.ps") #ps object with median abundances 
#load("Bioinformatics/trimmedJustsoils.ps") #this has the veg, canopy, pH that I need to average across transects
allMetaDat <- read.csv(file="SRS_allMetadata_Dec8.csv") #created in DataCleaning_pH_Jul23
#View(allMetaDat)

#### Get average pH for each EU at each point along the transect, by taking average of the 4 transects
head(allMetaDat)
class(allMetaDat)
# Get mean pH of all of the mean pHs for each EU at a given point along the transect
meanpH_EU <- allMetaDat %>% 
  group_by(EU, Meter) %>% 
  summarize(meanpHEu = mean_pH_vecFunc(mean_pH)) 
# View(meanpH_EU) #missing 53N_60 and 53N_100 because of an NA at 1/4 transects.
# Mssing T-60 and L_100 in 53N
# So, I'll have to add it back in manually (row 26 and row 30)
meanpH_EU <- meanpH_EU[-61,] #remove weird all NA final row (these were controls so don't have pH)
meanpH_EU[26,3] <- mean_pH_vecFunc(allMetaDat$mean_pH[c(86,96,106)]) #median of values in (meanpHs of) 53N_B_60, 53N_L_60, and 53N_R_60
meanpH_EU[30,3] <- mean_pH_vecFunc(allMetaDat$mean_pH[c(90,110,120)]) #median of values in (meanpHs of) 53N_B_100, 53N_R_100, and 53N_T_100
#View(meanpH_EU)

# These lines below show that getting pH is working as expected
test <- intersect(which(allMetaDat$EU==10), which(allMetaDat$Meter == 100))
EU10_100_all <- allMetaDat[test,]
EU_10_100pHs <- c(EU10_100_all$pH_1, EU10_100_all$pH_2)

EU_10_100_meanpH_step1 <-10^(EU_10_100pHs*-1)
EU_10_100_meanpH_step2 <- mean(EU_10_100_meanpH_step1)
EU_10_100_meanpH_final <- -log10(EU_10_100_meanpH_step2)
mean_pH_vecFunc(EU_10_100pHs) ==   EU_10_100_meanpH_final #TRUE

# Get average canopy cover
avgCanopyEU <- allMetaDat %>% 
  group_by(EU, Meter) %>% 
  summarize(avgCanopyEu = mean(meanDens))
#View(avgCanopyEU)
avgCanopyEU <- avgCanopyEU[-61,] #remove weird last row of NAs

# Get average vegetation cover
avgVegEU <- allMetaDat %>% 
  group_by(EU, Meter) %>% 
  summarize(avgVegCov = mean(Percent_Vegetation_Cover))
avgVegEU <- avgVegEU[-61,] #remove weird last row of NAs

# Before merge, check that EU and meter are in same order (should be because dplyr!)
unique(meanpH_EU[,1:2] == avgCanopyEU[,1:2]) #all true
unique(avgCanopyEU[,1:2] == avgVegEU[,1:2]) #all true 

pHCanVegMedEU <- cbind(meanpH_EU, avgCanopyEU[,3], avgVegEU[,3])
#View(pHCanVegMedEU)

######################
# Plotting
 
##########
# 4 plots showing pH, vegetation cover, canopy cover, and temperature across transects

# Colors for legend
EUcolors <- c("EU 10" = "#d73027", "EU 52" = "#4575b4", "EU 53N" = "#fc8d59",
            "EU 53S" = "#91bfdb", "EU 54S" = "#fee090", "EU 8" = "#e0f3f8")

# pH
pHEUMed_plot <- ggplot(pHCanVegMedEU, aes(x=Meter, y=meanpHEu, group = EU, color = EU)) + 
  geom_line(size=2.7) +
  scale_colour_manual(values=c(`10`="#d73027",`52`="#4575b4",`53N`="#fc8d59",`53S`="#91bfdb",`54S`="#fee090",`8`="#e0f3f8"))+
  theme_bw() +
  scale_x_continuous("Meter on transect", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("pH across transect") + 
  geom_vline(xintercept = 50, linetype= "dashed",  color= "black", size=1.2) +
  labs(y= "Mean pH", size=10) +
  theme(text = element_text(size=15)) +
  theme(legend.position = "none") #remove legend
pHEUMed_plot

# Canopy cover
CanEU_plot <- ggplot(pHCanVegMedEU, aes(x=Meter, y=avgCanopyEu, group = EU, color = EU)) + 
  geom_line(size=2.7) +
  scale_colour_manual(values=c(`10`="#d73027",`52`="#4575b4",`53N`="#fc8d59",`53S`="#91bfdb",`54S`="#fee090",`8`="#e0f3f8"))+
  theme_bw() +
  scale_x_continuous("Meter on transect", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("Canopy cover across transect") + 
  geom_vline(xintercept = 50, linetype= "dashed",  color= "black", size=1.2) +
  labs(y= "Mean canopy cover (%)", size=10) +
  theme(text = element_text(size=15)) +
  theme(legend.position = "none") #remove legend
CanEU_plot

# Vegetation cover
vegEU_plot <- ggplot(pHCanVegMedEU, aes(x=Meter, y=avgVegCov, group = EU, color = EU)) + 
  geom_line(size=2.7) +
  scale_colour_manual(values=c(`10`="#d73027",`52`="#4575b4",`53N`="#fc8d59",`53S`="#91bfdb",`54S`="#fee090",`8`="#e0f3f8"))+
  theme_bw() +
  scale_x_continuous("Meter on transect", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("Vegetation cover across transect") + 
  geom_vline(xintercept = 50, linetype= "dashed",  color= "black", size=1.2) +
  labs(y= "Mean vegetation cover (%)", size=10) +
  theme(text = element_text(size=15)) +
  theme(legend.position = "none") #remove legend
  
# quartz()
vegEU_plot

# This plot has the legend included, and is used to make legend for whole figure
vegEU_plotwithLegend <- ggplot(pHCanVegMedEU, aes(x=Meter, y=avgVegCov, group = EU, color = EU)) + 
  geom_line(size=2.7) +
  scale_colour_manual(values=c(`10`="#d73027",`52`="#4575b4",`53N`="#fc8d59",`53S`="#91bfdb",`54S`="#fee090",`8`="#e0f3f8"))+
  theme_bw() +
  scale_x_continuous("Meter on transect", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("Vegetation cover across transect") + 
  geom_vline(xintercept = 50, linetype= "dashed",  color= "black", size=1.2) +
  labs(y= "Mean vegetation cover (%)", size=10) +
  theme(text = element_text(size=15)) +
  theme(legend.position="bottom",
        legend.key.size = unit(3, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=13)) #change legend text font size
# quartz()
vegEU_plotwithLegend


##### PLOT ALL TOGETHER #####
# quartz()
grid.arrange(pHEUMed_plot, CanEU_plot, vegEU_plot, ncol=3)

# Add habitat column:
pHCanVegMedEU$habitat <- rep(NA, nrow(pHCanVegMedEU))
pHCanVegMedEU$habitat[which(pHCanVegMedEU$Meter %in% c(10,20,30,40))] <- "patch"
pHCanVegMedEU$habitat[which(pHCanVegMedEU$Meter %in% c(60,70,80,90,100))] <- "matrix"
pHCanVegMedEU$habitat[which(pHCanVegMedEU$Meter %in% 50)] <- "edge"

############################################################
############### CANOPY COVER "COLOR BARS" ###################
############################################################

# Aug 9, 2022 (based on tutorial here: https://ggplot2.tidyverse.org/reference/guide_colourbar.html)

#### EU 10
df <- expand.grid(X1 = 1:10, X2 = 1:10)
df$value <- df$X1 * df$X2

pHCanVegMedEU$EU[1:10]
pHCanVegMedEU[1:10,4] #just EU_10 canopy cover data

# What is the mean difference between forest and patch canopy cover
EU10_avgPatchCan <- mean(pHCanVegMedEU$avgCanopyEu[1:4]) #mean patch
EU10_avgForestCan <- mean(pHCanVegMedEU$avgCanopyEu[6:10]) #mean forest
canContrast_EU10 <- EU10_avgForestCan - EU10_avgPatchCan
canContrast_EU10

canEU10_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
dim(canEU10_df)
canCoverVec10 <- pHCanVegMedEU$avgCanopyEu[1:10] #add  just the canopy cover data for EU 10
canEU10_df$value <- rep(canCoverVec10, 10) #add canopy cover data
# View(canEU10_df) #Now, each row in this 10 x 10 grid should be the same, and should be the mean canopy
# cover across transect 10
canEU10_df$value #since only the first row is what we are using for the color bar, we will add in some 
# high and low numbers in the first row so that color bar goes from 0 to 100 percent canopy cover. So any row but the first one can be used as the color bar
canEU10_df$value[1:10] <- c(rep(0,5), rep(100,5))

canCoverEU_10_grid <- ggplot(canEU10_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 10") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
#quartz() # Horizontal color bar legend
canCoverEU_10_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                        direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()
# Vertical color bar legend
#quartz() 
canCoverEU_10_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                   direction= "vertical", barwidth = 1.5, barheight = 6)) + 
  theme_bw()

###########################################################
### 53S
canEU53S_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
dim(canEU53S_df)
pHCanVegMedEU$EU[31:40] #53S

# What is the mean difference between forest and patch canopy cover
EU53S_avgPatchCan <- mean(pHCanVegMedEU$avgCanopyEu[31:34]) #mean patch
EU53S_avgForestCan <- mean(pHCanVegMedEU$avgCanopyEu[36:40]) #mean forest
canContrast_EU53S <- EU53S_avgForestCan - EU53S_avgPatchCan
canContrast_EU53S

canCoverVec53S <- pHCanVegMedEU$avgCanopyEu[31:40] #add  just the canopy cover data for EU 53S

canEU53S_df$value <- rep(canCoverVec53S, 10) #add canopy cover data
# View(canEU53S_df) #Now, each row in this 10 x 10 grid should be the same, and should be the mean canopy
# cover across transect 53S
canEU53S_df$value #since only the first row is what we are using for the color bar, we will add in some 
# high and low numbers in the first row so that color bar goes from 0 to 100 percent canopy cover. So any row but the first one can be used as the color bar
canEU53S_df$value[1:10] <- c(rep(0,5), rep(100,5))

canCoverEU_53S_grid <- ggplot(canEU53S_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 53S") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
#quartz() # Horizontal color bar legend
canCoverEU_53S_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                                    direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()

# Vertical color bar legend
#quartz() 
canCoverEU_53S_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                    direction= "vertical", barwidth = 1.5, barheight = 6)) + 
  theme_bw()

###########################################################
### 54S
canEU54S_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
dim(canEU54S_df)
pHCanVegMedEU$EU[41:50] #54S

# What is the mean difference between forest and patch canopy cover
EU54S_avgPatchCan <- mean(pHCanVegMedEU$avgCanopyEu[41:44]) #mean patch
EU54S_avgForestCan <- mean(pHCanVegMedEU$avgCanopyEu[46:50]) #mean forest
canContrast_EU54S <- EU54S_avgForestCan - EU54S_avgPatchCan

canCoverVec54S <- pHCanVegMedEU$avgCanopyEu[41:50] #add  just the canopy cover data for EU 54S
canEU54S_df$value <- rep(canCoverVec54S, 10) #add canopy cover data
# View(canEU54S_df) #Now, each row in this 10 x 10 grid should be the same, and should be the mean canopy
# cover across transect 54S
canEU54S_df$value #since only the first row is what we are using for the color bar, we will add in some 
# high and low numbers in the first row so that color bar goes from 0 to 100 percent canopy cover. So any row but the first one can be used as the color bar
canEU54S_df$value[1:10] <- c(rep(0,5), rep(100,5))

canCoverEU_54S_grid <- ggplot(canEU54S_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 54S") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
#quartz() # Horizontal color bar legend
canCoverEU_54S_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                                    direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()
# Vertical color bar legend
canCoverEU_54S_gridVert <- canCoverEU_54S_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                    direction= "vertical", barwidth = 1.5, barheight = 6, ticks.colour = NA, frame.color = NULL)) + 
  theme_bw()
#quartz() 
canCoverEU_54S_gridVert

###########################################################
### 53N
canEU53N_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
dim(canEU53N_df)
pHCanVegMedEU$EU[21:30] #53N

# What is the mean difference between forest and patch canopy cover
EU53N_avgPatchCan <- mean(pHCanVegMedEU$avgCanopyEu[21:24]) #mean patch
EU53N_avgForestCan <- mean(pHCanVegMedEU$avgCanopyEu[26:30]) #mean forest
canContrast_EU53N <- EU53N_avgForestCan - EU53N_avgPatchCan

canCoverVec53N <- pHCanVegMedEU$avgCanopyEu[21:30] #add  just the canopy cover data for EU 53N
canEU53N_df$value <- rep(canCoverVec53N, 10) #add canopy cover data
# View(canEU53N_df) #Now, each row in this 10 x 10 grid should be the same, and should be the mean canopy
# cover across transect 53N
canEU53N_df$value #since only the first row is what we are using for the color bar, we will add in some 
# high and low numbers in the first row so that color bar goes from 0 to 100 percent canopy cover. So any row but the first one can be used as the color bar
canEU53N_df$value[1:10] <- c(rep(0,5), rep(100,5))


canCoverEU_53N_grid <- ggplot(canEU53N_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 53N") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
#quartz() # Horizontal color bar legend
canCoverEU_53N_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                                    direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()


# Vertical color bar legend
#quartz() 
canCoverEU_53N_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                    direction= "vertical", barwidth = 1.5, barheight = 6)) + 
  theme_bw()

###########################################################
### 52
canEU52_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
dim(canEU52_df)
pHCanVegMedEU$EU[11:20] #52

# What is the mean difference between forest and patch canopy cover?
EU52_avgPatchCan <- mean(pHCanVegMedEU$avgCanopyEu[11:14]) #mean patch
EU52_avgForestCan <- mean(pHCanVegMedEU$avgCanopyEu[16:20]) #mean forest
canContrast_EU52 <- EU52_avgForestCan - EU52_avgPatchCan
canContrast_EU52

canCoverVec52 <- pHCanVegMedEU$avgCanopyEu[11:20] #add  just the canopy cover data for EU 52
canEU52_df$value <- rep(canCoverVec52, 10) #add canopy cover data
# View(canEU52_df) #Now, each row in this 10 x 10 grid should be the same, and should be the mean canopy
# cover across transect 52
canEU52_df$value #since only the first row is what we are using for the color bar, we will add in some 
# high and low numbers in the first row so that color bar goes from 0 to 100 percent canopy cover. So any row but the first one can be used as the color bar
canEU52_df$value[1:10] <- c(rep(0,5), rep(100,5))

canCoverEU_52_grid <- ggplot(canEU52_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 52") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
#quartz() # Horizontal color bar legend
canCoverEU_52_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                                   direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()
# Vertical color bar legend
#quartz() 
canCoverEU_52_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                   direction= "vertical", barwidth = 1.5, barheight = 6)) + 
  theme_bw()

###########################################################
### 8
canEU8_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
pHCanVegMedEU$EU[51:60] #8

# What is the mean difference between forest and patch canopy cover?
EU8_avgPatchCan <- mean(pHCanVegMedEU$avgCanopyEu[51:54]) #mean patch
EU8_avgForestCan <- mean(pHCanVegMedEU$avgCanopyEu[56:60]) #mean forest
canContrast_EU8 <- EU8_avgForestCan - EU8_avgPatchCan
canContrast_EU8

canCoverVec8 <- pHCanVegMedEU$avgCanopyEu[51:60] #add  just the canopy cover data for EU 8
canEU8_df$value <- rep(canCoverVec8, 10) #add canopy cover data
# View(canEU8_df) 
canEU8_df$value #since only the first row is what we are using for the color bar, we will add in some 
# high and low numbers in the first row so that color bar goes from 0 to 100 percent canopy cover. So any row but the first one can be used as the color bar
canEU8_df$value[1:10] <- c(rep(0,5), rep(100,5))

canCoverEU_8_grid <- ggplot(canEU8_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 8") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
#quartz() # Horizontal color bar legend
canCoverEU_8_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                                   direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()
# Vertical color bar legend
#quartz() 
canCoverEU_8_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                   direction= "vertical", barwidth = 1.5, barheight = 6)) + 
  theme_bw()

meanCanopyDiff <- as.data.frame(matrix(nrow=6, ncol=4))
colnames(meanCanopyDiff) <- c("EU", "meanCanopyDiff", "meanPatchCan", "meanForestCan")
meanCanopyDiff[1,] <- c("EU_10", canContrast_EU10, EU10_avgPatchCan, EU10_avgForestCan)
meanCanopyDiff[2,] <- c("EU_53S", canContrast_EU53S, EU53S_avgPatchCan, EU53S_avgForestCan)
meanCanopyDiff[3,] <- c("EU_54S", canContrast_EU54S, EU54S_avgPatchCan, EU54S_avgForestCan)
meanCanopyDiff[4,] <- c("EU_53N", canContrast_EU53N, EU53N_avgPatchCan, EU53N_avgForestCan)
meanCanopyDiff[5,] <- c("EU_52", canContrast_EU52, EU52_avgPatchCan, EU52_avgForestCan)
meanCanopyDiff[6,] <- c("EU_8", canContrast_EU8, EU8_avgPatchCan, EU8_avgForestCan)
#View(meanCanopyDiff)

# save pH, canopy cover, and vegetation data for transects
# save(pHCanVegMedEU, file="pHCanVegMedEU")

############################################################
#                  TRANSECT SOIL DATA 
############################################################
# Bring in soil end transect data. NOTE: Although these data also have
# pH, I won't use it since I also took it.
soilDat <- read.csv("LabWorkAndSoilAnalyses/cleanedUpSoilAnalysisDataFromCSU_F4-F51.csv")
head(soilDat)

# Add EU and meter columns to soilDat
for (i in 1:nrow(soilDat)){
  if (grepl(pattern= "100", x=soilDat[i,1])){
    soilDat$Meter[i] <- 100
  } else {
    soilDat$Meter[i] <- 10
  }
  # add in an additional if statement to account for "90" pattern
  if (grepl(pattern= "90", x=soilDat[i,1])){
    soilDat$Meter[i] <- 100
  }
  if (grepl(pattern= "52D", x=soilDat[i,1])){
    soilDat$EU[i] <- "52"
  }
  if (grepl(pattern= "53SD", x=soilDat[i,1])){
    soilDat$EU[i] <- "53S"
  }
  if (grepl(pattern= "53N", x=soilDat[i,1])){
    soilDat$EU[i] <- "53N"
  }
  if (grepl(pattern= "10C", x=soilDat[i,1])){
    soilDat$EU[i] <- "10"
  }
  if (grepl(pattern= "54S", x=soilDat[i,1])){
    soilDat$EU[i] <- "54S"
  }
  if (grepl(pattern= "8D", x=soilDat[i,1])){
    soilDat$EU[i] <- "8"
  }
}

#View(soilDat)

head(pHCanVegMedEU)
head(soilDat)
class(soilDat$Meter)
soilDat$Meter <- as.character(soilDat$Meter)

# How do some of these characteristics vary by EU?
aggregate(percentOM ~ EU, data= soilDat, mean)
summary(lm(percentOM  ~ -1 + EU, data= soilDat))
confint(lm(percentOM  ~ -1 + EU, data= soilDat))
aggregate(NO3.N ~ EU, data= soilDat, mean) #this is same across the board!
aggregate(P ~ EU, data= soilDat, mean)
summary(lm(P  ~ -1 + EU, data= soilDat))
confint(lm(P  ~ -1 + EU, data= soilDat)) #while all these confidence intervals overlap, EU52 clearly has more P overall
aggregate(K ~ EU, data= soilDat, mean)
summary(lm(K  ~ -1 + EU, data= soilDat))
confint(lm(K  ~ -1 + EU, data= soilDat)) #while all these confidence intervals overlap, EU52 clearly has more P overall

# How do some of these characteristics vary by meter? And stat them!
aggregate(percentOM ~ Meter, data= soilDat, mean)
t.test(percentOM ~ Meter, data= soilDat) #greater in 100 than in 10. t = -2.5289, df = 38.217, p-value = 0.01569
aggregate(NO3.N ~ Meter, data= soilDat, mean)
# Nitrate-nitrogen has the same value for both, so t.test doesn't make sense
aggregate(P ~ Meter, data= soilDat, mean)
t.test(P ~ Meter, data= soilDat) #t = 1.5492, df = 22.84, p-value = 0.1351. NOT significant, but KINDA close (with 10 more than 100) 
# mean for 10 is 19.53478, mean for meter 100 is 12.17600)
aggregate(K ~ Meter, data= soilDat, mean)
t.test(K ~ Meter, data= soilDat) # No difference among 10 and 100, t = -0.036609, df = 45.137, p-value = 0.971

# How do these characteristics vary by EU and Meter (i.e. 10 or 100)?
# Test of aggregation:
soilDat[soilDat$EU == "52",c(1,5)][c(1,4,5,8),] #just meter ten
mean(soilDat[soilDat$EU == "52",c(1,5)][c(1,4,5,8),]$percentOM) #this is for the 10 meters of EU 52
# shows that aggregate is working to clump meter and EU correctly 
aggregate(percentOM  ~ EU + Meter, data= soilDat, mean) 
aggregate(NO3.N  ~ EU + Meter, data= soilDat, mean) 
# They are all the same!
aggregate(P  ~ EU + Meter, data= soilDat, mean) 
aggregate(K  ~ EU + Meter, data= soilDat, mean) 


