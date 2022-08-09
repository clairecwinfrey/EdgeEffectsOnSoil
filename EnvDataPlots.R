# Environemtnal and physicochemical data and plots
# December 8, 2021

# This script creates plots for soil physiochemical and environmental data.

setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

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

# 2. metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

############################################################

load("Bioinformatics/medianEU.ps") #ps object with median abundances 
#load("Bioinformatics/trimmedJustsoils.ps") #this has the veg, canopy, pH that I need to average across transects
allMetaDat <- read.csv(file="SRS_allMetadata_Dec8.csv") #created in DataCleaning_pH_Jul23
View(allMetaDat)

#### Get average pH for each EU at each point along the transect, by taking average of the 4 transects
head(allMetaDat)
class(allMetaDat)

# Get average pH
medianpH_EU <- allMetaDat %>% 
  group_by(EU, Meter) %>% 
  summarize(medpHEu = median(mean_pH))
# View(medianpH_EU) #missing 53N_60 and 53N_100 because of an NA at 1/4 transects.
# Mssing T-60 and L_100 in 53N
# So, I'll have to add it back in manually (row 26 and row 30)
medianpH_EU <- medianpH_EU[-61,] #remove weird all NA final row
medianpH_EU[26,3] <- median(allMetaDat$mean_pH[c(86,96,106)]) #median of values in (meanpHs of) 53N_B_60, 53N_L_60, and 53N_R_60
medianpH_EU[30,3] <- median(allMetaDat$mean_pH[c(90,110,120)]) #median of values in (meanpHs of) 53N_B_100, 53N_R_100, and 53N_T_100

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
unique(medianpH_EU[,1:2] == avgCanopyEU[,1:2]) #all true
unique(avgCanopyEU[,1:2] == avgVegEU[,1:2]) #all true 

pHCanVegMedEU <- cbind(medianpH_EU, avgCanopyEU[,3], avgVegEU[,3])
View(pHCanVegMedEU)

######################
# Plotting

# A simple, base R plot for canopy cover and  vegetation cover (EU 8)
quartz()
plot(envEU8$Meter, 
     envEU8$avgCanopyEu, type= "l",
cex.lab=1.3, 
xlab="Meter", 
ylab="Percent Cover", 
xaxp = c(10, 100, 9.9), 
ylim=c(0,100), 
cex=1.4,
lwd=3)
lines(envEU8$Meter, 
      envEU8$avgVegCov, 
      lwd=3,
      lty="dashed")
legend("topleft", 
       legend= c("Canopy Cover", "Vegetation Cover"), 
       lty= c("solid", "dashed"), 
       bty= "o", 
       lwd =2, 
       cex=.8, 
       inset= c(0.1, 0.1))
 
##########
# 4 plots showing pH, vegetation cover, canopy cover, and temperature across transects

# Colors for legend
EUcolors <- c("EU 10" = "#d73027", "EU 52" = "#4575b4", "EU 53N" = "#fc8d59",
            "EU 53S" = "#91bfdb", "EU 54S" = "#fee090", "EU 8" = "#e0f3f8")

# pH
pHEUMed_plot <- ggplot() + 
  geom_line(data=pHCanVegMedEU[1:10,], aes(x=Meter, y=medpHEu, group = EU), color = "#d73027", size= 1.5) +
  geom_line(data=pHCanVegMedEU[11:20,], aes(x=Meter, y=medpHEu, group = EU), color = "#4575b4", size= 1.5) +
  geom_line(data=pHCanVegMedEU[21:30,], aes(x=Meter, y=medpHEu, group = EU), color = "#fc8d59", size= 1.5) +
  geom_line(data=pHCanVegMedEU[31:40,], aes(x=Meter, y=medpHEu, group = EU), color = "#91bfdb", size= 1.5) +
  geom_line(data=pHCanVegMedEU[41:50,], aes(x=Meter, y=medpHEu, group = EU), color = "#fee090", size= 1.5) +
  geom_line(data=pHCanVegMedEU[51:60,], aes(x=Meter, y=medpHEu, group = EU), color = "#e0f3f8", size= 1.5) +
  theme_bw() + ylim(3, 7) +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("Mean pH across transect") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey") +
  labs(y= "Mean pH") + theme_bw()
quartz()
pHEUMed_plot

# Canopy cover
CanEU_plot <- ggplot() + 
  geom_line(data=pHCanVegMedEU[1:10,], aes(x=Meter, y=avgCanopyEu), color = "#d73027", size= 1.5) +
  geom_line(data=pHCanVegMedEU[11:20,], aes(x=Meter, y=avgCanopyEu), color = "#4575b4", size= 1.5) +
  geom_line(data=pHCanVegMedEU[21:30,], aes(x=Meter, y=avgCanopyEu), color = "#fc8d59", size= 1.5) +
  geom_line(data=pHCanVegMedEU[31:40,], aes(x=Meter, y=avgCanopyEu), color = "#91bfdb", size= 1.5) +
  geom_line(data=pHCanVegMedEU[41:50,], aes(x=Meter, y=avgCanopyEu), color = "#fee090", size= 1.5) +
  geom_line(data=pHCanVegMedEU[51:60,], aes(x=Meter, y=avgCanopyEu), color = "#e0f3f8", size= 1.5) +
  theme_bw() + 
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("Mean canopy cover across transect") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey") +
  labs(y= "Mean canopy cover (%)") + theme_bw()
quartz()
CanEU_plot

# Vegetation cover
vegEU_plot <- ggplot() + 
  geom_line(data=pHCanVegMedEU[1:10,], aes(x=Meter, y=avgVegCov, group = EU), color = "#d73027", size= 1.5) +
  geom_line(data=pHCanVegMedEU[11:20,], aes(x=Meter, y=avgVegCov, group = EU), color = "#4575b4", size= 1.5) +
  geom_line(data=pHCanVegMedEU[21:30,], aes(x=Meter, y=avgVegCov, group = EU), color = "#fc8d59", size= 1.5) +
  geom_line(data=pHCanVegMedEU[31:40,], aes(x=Meter, y=avgVegCov, group = EU), color = "#91bfdb", size= 1.5) +
  geom_line(data=pHCanVegMedEU[41:50,], aes(x=Meter, y=avgVegCov, group = EU), color = "#fee090", size= 1.5) +
  geom_line(data=pHCanVegMedEU[51:60,], aes(x=Meter, y=avgVegCov, group = EU), color = "#e0f3f8", size= 1.5) +
  theme_bw() +
  scale_x_continuous("Transect Meter", breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  ggtitle("Mean vegetation cover across transect") + geom_vline(xintercept = 50, linetype= "dashed", color= "grey") +
  labs(y= "Mean vegetation cover (%)") + theme_bw()
quartz()
vegEU_plot

##### PLOT ALL TOGETHER #####
quartz()
require(gridExtra)
grid.arrange(pHEUMed_plot, CanEU_plot, vegEU_plot, ncol=3)

####### CANOPY COVER "COLOR BARS" #########
# Aug 9, 2022 (based on tutorial here: https://ggplot2.tidyverse.org/reference/guide_colourbar.html)
########################
df <- expand.grid(X1 = 1:10, X2 = 1:10)
df$value <- df$X1 * df$X2

pHCanVegMedEU[1:10,4] #just EU_10 canopy cover data
pHCanVegMedEU[31:40,4] #just EU_53S canopy cover data

canEU10_df <- expand.grid(X = 1:10, Y = 1:10) #make a 10 x 10 grid,
# where X is the x-axis, and Y the y-axis
dim(canEU10_df)
canCoverVec10 <- pHCanVegMedEU$avgCanopyEu[1:10] #add  just the canopy cover data for EU 10
canEU10_df$value <- rep(canCoverVec10, 10) #add canopy cover data
# View(canEU10_df) #Now, each row in this 10 x 10 grid should be the same, and should be the mean canopy
# cover across transect 10

canCoverEU_10_grid <- ggplot(canEU10_df, aes(X, Y)) + geom_tile(aes(fill = value)) + 
  ggtitle("Mean Canopy Cover EU 10") +
  scale_fill_gradient(low = "white", high = "darkgreen")
# Horizontal color bar legend
quartz() # Horizontal color bar legend
canCoverEU_10_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "top", nbin = 100, title= "Canopy cover",
                                        direction= "horizontal", barwidth = 6, barheight = 1.5)) + 
  theme_bw()
# Vertical color bar legend
quartz() 
canCoverEU_10_grid + guides(fill = guide_colourbar(ticks=FALSE, label.position = "left", nbin = 100, title= "Canopy cover",
                                                   direction= "vertical", barwidth = 1.5, barheight = 6)) + 
  theme_bw()


# save pH, canopy cover, and vegetation data for transects
# save(pHCanVegMedEU, file="pHCanVegMedEU")


