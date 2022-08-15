# EdgeEffectbyASV_BYEU -- ITS/FUNGI
# Aug 14, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect FUNGAL (ITS) ASVs using differential abundance 
# analyses (DESeq2). IT IS NEARLY IDENTICAL TO EDGEEFFECTSALLSITESFUNGI.R. The major
# change is fewer scratch work and the Z-score is calculated differently (see 
# explanation below).

# This script:
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. ITS_postUbiquity.ps)).
# 2. Creates data frame with z-score (of ASV abundance) for each ASV within each EU calculated as:
# --- (1) Mean abundance of each ASV within a given meter on the transect (e.g. mean abundance 
#          of ASV_1 across all transects at meter 10 in EU 10)
# --  (2) Subtracts the mean of that ASV across all points in that EU (e.g. ASV_1's mean value across all 
#           points in EU_10) from the value in (1)
# --  (3) Finally, the difference of (1) - (2) is divided by the standard deviation of the ASVs abundance 
#           across all points in that EU.    
# 3. Fits logistic curves to each one of the differentially abundant ASVs using the Z-scores calculated 
#           above
# 4. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential 
#          abundance breakdown by forest, patch, and non-differentially abundant microbes 


# IdentifyEdgeEffectTaxa5.R investigated linear fit and logistiC model fit of the
# the data, and found that logistic models were much better (for bacterial/archaeal data)

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("readxl")       # necessary to import the data from Excel file
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("grid")
library("stringr") #for grep-like tools for data manipulation with character strings in data

# Load data
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: UbiquitypostUbiqSetup.R

######

# FUNCTIONS DEFINED IN THIS SCRIPT (but often used first in other scripts):
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVsampOccurrence determines the number of samples that each ASV occurs in 
# For proof that ASVsampOccurrence works, see UbiquitypostUbiqSetup.R
ASVsampOccurrence <- function(physeq) { #input is a phyloseq object
  OTU <- otu_table(physeq)
  ASVmat <- as(OTU, "matrix") # make phyloseq ASV table into a non-phyloseq matrix
  sampleOccurence <- rowSums(ASVmat > 0) #sum up number of samples ASV occurs in 
  return(cbind(ASVmat, sampleOccurence)) #
}

# 3. ASVs_outta_ps gets the ASV table out of phyloseq 
# from: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  if(taxa_are_rows(ASVTable)) {
    ASVTable <- t(ASVTable)
  }
  return(as.data.frame(ASVTable))
}

# THIS FUNCTION DOES NOT SEEM TO ACTUALLY STRIP THE ATTRIBUTE :(
# 4. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

# 5. # log.fit function creates a logistic model and plots ASV values
# To write this, I adapted the original function that Julian sent to me on Jan 18, based on the link and
# email explanation he sent me on Dec 11, 2021. The link that he got the equation from was:
# https://stats.stackexchange.com/questions/47802/whats-the-most-pain-free-way-to-fit-logistic-growth-curves-in-r
log.fit <- function(y, x){ 
  # y = a/ (1 + bc^-x)
  log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
  # phi1 = Asym = a numeric parameter representing the asymptote
  # phi2 = xmid = a numeric parameter representing the x value AT THE INFLECTION point of the curve
  # scal = "a numeric scale parameter on the input axis. These values are created by nls function
  # by creating "initial estimates"
  
  # Here, I think that the C parameter is kinda like the magnitude or the asymptotde. Is basically the limiting value or where the 
  # function rises up to and eventually levels off
  C <- summary(log.ss)$coef[1] #this gives the estimate for phi1 
  # a 
  A <- exp((summary(log.ss)$coef[2]) * (1/summary(log.ss)$coef[3])) 
  # phi2 = xmid = a numeric parameter representing the x value AT THE INFLECTION point of the curve
  #k
  K <- (1 / summary(log.ss)$coef[3])
  #scal = "a numeric scale parameter on the input axis. This is the reciprocal of that.
  
  plot(y ~ x, main = "Logistic Function", xlab= "distance (m)", ylab= "ASV abundance")
  lines(0:max(x), predict(log.ss, data.frame(x=0:max(x))), col="red")
  
  r1 <- sum((x - mean(x))^2)
  r2 <- sum(residuals(log.ss)^2)
  
  r_sq <- (r1 - r2) / r1 #check, is this right?
  
  out <- data.frame(cbind(c(C=C, a=A, k=K, R.value=sqrt(r_sq))))
  names(out)[1] <- "Logistic Curve"
  
  return(out)
}

# 6. 
log.fitdiffAbundFunct <- function(y, x, ASVnames){ #one important difference here is that ASV names name the ASVs as it plot it
  # y = a/ (1 + bc^-x)
  log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
  # phi1 = Asym = a numeric parameter representing the asymptote
  # phi2 = xmid = a numeric parameter representing the x value AT THE INFLECTION point of the curve
  # scal = "a numeric scale parameter on the input axis. These values are created by nls function
  # by creating "initial estimates"
  
  # Here, I think that the C parameter is kinda like the magnitude or the asymptotde. Is basically the limiting value or where the 
  # function rises up to and eventually levels off
  C <- summary(log.ss)$coef[1] #this gives the estimate for phi1 
  # a 
  A <- exp((summary(log.ss)$coef[2]) * (1/summary(log.ss)$coef[3])) 
  # phi2 = xmid = a numeric parameter representing the x value AT THE INFLECTION point of the curve
  #k
  K <- (1 / summary(log.ss)$coef[3])
  #scal = "a numeric scale parameter on the input axis. This is the reciprocal of that.
  for (j in 1:length(ASVnames)){
    plot(y ~ x, main = paste("Logistic Function for", ASVnames[j]), xlab= "distance (m)", ylab= "ASV abundance (Z-score)")
    lines(0:max(x), predict(log.ss, data.frame(x=0:max(x))), col="red")
  }
  r1 <- sum((x - mean(x))^2)
  r2 <- sum(residuals(log.ss)^2)
  
  r_sq <- (r1 - r2) / r1 #check, is this right?
  
  out <- data.frame(cbind(c(C=C, a=A, k=K, R.value=sqrt(r_sq))))
  names(out)[1] <- "Logistic Curve"
  
  return(out)
}

############################################################

##########################################################################
# 1.    INDICATOR ANALYSIS AND TIDY DATAFRAME
##########################################################################
# Perform indicator species analysis on the post ubiquity dataset, and then
# make a dataframe which has a row corresponging to the abundance of each differentially
# abundant ASV in each place along the transect, as well as taxonomic info, and 
# whether or not that ASV was differentially abundant in patch or the forest.

# The next steps add 1 to all of the ASV abundance counts and re-make a phyloseq object.
# This was necessary because DESeq2 could not compute the geometric mean of the samples
# since every gene had at least one zero in it (and was thus thrown out of the analysis)
# (see recommendations and explanations here:
# https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564)
ITS_postUbiqASVs <- ASVs_outta_ps(ITS_postUbiquity.ps)
ITS_postUbiqASVsPlus1 <- t(ITS_postUbiqASVs + 1) #adding one to every abundance count for differential abundance analysis;
# see explanation a few lines down. Invert so that I can make into a new phyloseq object.
# Make a new phyloseq object with these above
OTU = otu_table(ITS_postUbiqASVsPlus1, taxa_are_rows = TRUE)
ITS_postUbiqASVsPlus1.ps <- phyloseq(OTU, tax_table(ITS_postUbiquity.ps), sample_data(ITS_postUbiquity.ps))

# Remove edge samples to compare patch and matrix
NoEdgePlus1.ps <- subset_samples(ITS_postUbiqASVsPlus1.ps, Habitat != "edge") 
step1 <- phyloseq::phyloseq_to_deseq2(NoEdgePlus1.ps, ~ Habitat) #set up DESeq2 dds object 
step2 <- DESeq2::DESeq(step1, test="Wald", fitType = "parametric") #differential expression analysis step
DeseqResults <- results(step2, cooksCutoff = FALSE) #make results object
DeseqResults <- DeseqResults[which(DeseqResults$padj < 0.001), ] #only get those ASVs below alpha level
DeseqResults <- cbind(as(DeseqResults, "data.frame"), as(tax_table(NoEdgePlus1.ps)[rownames(DeseqResults), ], "matrix")) #clean up format
DeseqResults$Habitat <- ifelse(DeseqResults$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
# Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
# View(DeseqResults) #287 diff abund ASVs)
sampDat <- sample_data(ITS_postUbiqASVsPlus1.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ASVnamesDA_FUNGI <- rownames(DeseqResults) #1696 found
ASVsAll <- as.data.frame(t(ASVs_outta_ps(ITS_postUbiqASVsPlus1.ps))) #get ASVs from original and transpose so that ASVs are rows 
# Create dataframe with everything of interest
diffAbunDat <- merge(DeseqResults, ASVsAll, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
diffAbunDat_tidy_FUNGI <- diffAbunDat %>% 
  pivot_longer(cols= 16:ncol(diffAbunDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_tidy_FUNGI)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
#View(diffAbunDat_tidy_FUNGI) this has 66,871 rows, which is equal to 233 (number of samples) x 287 (number of diff abund ASVs)
##########

#save(diffAbunDat_tidy_FUNGI, file= "RobjectsSaved/diffAbunDat_tidy_FUNGI") #saved Aug 14, 2022

##########################################################################
# 2.  Z-SCORES OF EACH ASV WITHIN EACH EU 
##########################################################################
# 1. Format diffAbunDat_tidy_FUNGI nd a pre-allocated list for analysis.
# Pull out only forest specialists
diffAbunDat_tidy_FUNGI_forest <- diffAbunDat_tidy_FUNGI %>% 
  filter(Habitat == "forest")
# Make a vector of all of the 119 forest-specializing ASVs
ASVnamesDA_FUNGI_forest <- unique(diffAbunDat_tidy_FUNGI_forest$ASV_name) 
# View(diffAbunDat_tidy_FUNGI_forest)

ASVmeterZAbunds_fungi <- vector("list", length(ASVnamesDA_FUNGI_forest)) #pre-allocate space for each ASV's abundance for each ASV


############ TRYING OUT A FEW METHODS OF Z-SCORES ############
# test with ASV_1 in EU_10
# THIS MAKES A DATAFRAME THAT I CAN USE FOR ALL VALUES OF ASV_1 ACROSS ALL EUs
ASV1meansAcrossEUs <- as.data.frame(matrix(nrow=6, ncol=13)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
colnames(ASV1meansAcrossEUs) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd")
#names(new) <- paste(ASVnamesDA_FUNGI_forest[1],"acrossEUs", sep= "")


######### EU_10 #######
# First, isolate 1 ASV in EU_10
ASV_1_EU_10_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_10") #at this point, have only ASV_1 in EU_10
head(ASV_1_EU_10_tidyTest_1)
head(ASV1meansAcrossEUs)

# Get ASV mean and stdev across all points on the transect
ASV1meansAcrossEUs$EU[1] <- unique(ASV_1_EU_10_tidyTest_1$EU)
ASV1meansAcrossEUs$ASV_EUmean[1] <- mean(ASV_1_EU_10_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU 10
ASV1meansAcrossEUs$ASV_EUsd[1] <- sd(ASV_1_EU_10_tidyTest_1$ASVabundance)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV1_EU10 <- sort(unique(ASV_1_EU_10_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV1_EU10)){ #loop over number of unique meters (i.e. 10)
  ASV1meansAcrossEUs[1,which(colnames(ASV1meansAcrossEUs)==
                               paste("mean_",uniqueMetersASV1_EU10[i], sep=""))] <- mean(ASV_1_EU_10_tidyTest_1
                                [which(ASV_1_EU_10_tidyTest_1$Meter==uniqueMetersASV1_EU10[i]), 17]) #mean for meter 10
}

############### Repeat this for all of the other EUs ########################
unique(diffAbunDat_tidy_FUNGI$EU)
######### EU_52 #########
ASV_1_EU_52_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_52") 
head(ASV_1_EU_52_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV1meansAcrossEUs$EU[2] <- unique(ASV_1_EU_52_tidyTest_1$EU)
ASV1meansAcrossEUs$ASV_EUmean[2] <- mean(ASV_1_EU_52_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV1meansAcrossEUs$ASV_EUsd[2] <- sd(ASV_1_EU_52_tidyTest_1$ASVabundance)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV1_EU52 <- sort(unique(ASV_1_EU_52_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV1_EU52)){ #loop over number of unique meters (i.e. 10)
  ASV1meansAcrossEUs[2,which(colnames(ASV1meansAcrossEUs)==
                               paste("mean_",uniqueMetersASV1_EU52[i], sep=""))] <- mean(ASV_1_EU_52_tidyTest_1
                                                                                [which(ASV_1_EU_52_tidyTest_1$Meter==uniqueMetersASV1_EU52[i]), 17]) #mean for meter 10
}
  
###########  EU_53N ########## 
ASV_1_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_53N") 
head(ASV_1_EU_53N_tidyTest_1)
ASV_1_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_53N") 
head(ASV_1_EU_53N_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV1meansAcrossEUs$EU[3] <- unique(ASV_1_EU_53N_tidyTest_1$EU)
ASV1meansAcrossEUs$ASV_EUmean[3] <- mean(ASV_1_EU_53N_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV1meansAcrossEUs$ASV_EUsd[3] <- sd(ASV_1_EU_53N_tidyTest_1$ASVabundance)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV1_EU53N <- sort(unique(ASV_1_EU_53N_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV1_EU53N)){ #loop over number of unique meters (i.e. 10)
  ASV1meansAcrossEUs[3,which(colnames(ASV1meansAcrossEUs)==
                               paste("mean_",uniqueMetersASV1_EU53N[i], sep=""))] <- mean(ASV_1_EU_53N_tidyTest_1
                                                                                          [which(ASV_1_EU_53N_tidyTest_1$Meter==uniqueMetersASV1_EU53N[i]), 17]) #mean for meter 10
}

###########  EU_54S ########## 
ASV_1_EU_54S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_54S") 
head(ASV_1_EU_54S_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV1meansAcrossEUs$EU[4] <- unique(ASV_1_EU_54S_tidyTest_1$EU)
ASV1meansAcrossEUs$ASV_EUmean[4] <- mean(ASV_1_EU_54S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV1meansAcrossEUs$ASV_EUsd[4] <- sd(ASV_1_EU_54S_tidyTest_1$ASVabundance)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV1_EU54S <- sort(unique(ASV_1_EU_54S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV1_EU54S)){ #loop over number of unique meters (i.e. 10)
  ASV1meansAcrossEUs[4,which(colnames(ASV1meansAcrossEUs)==
                               paste("mean_",uniqueMetersASV1_EU54S[i], sep=""))] <- mean(ASV_1_EU_54S_tidyTest_1
                                                                                          [which(ASV_1_EU_54S_tidyTest_1$Meter==uniqueMetersASV1_EU54S[i]), 17]) #mean for meter 10
}

##########  EU_8 ########## 
ASV_1_EU_8_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_8") 
head(ASV_1_EU_8_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV1meansAcrossEUs$EU[5] <- unique(ASV_1_EU_8_tidyTest_1$EU)
ASV1meansAcrossEUs$ASV_EUmean[5] <- mean(ASV_1_EU_8_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV1meansAcrossEUs$ASV_EUsd[5] <- sd(ASV_1_EU_8_tidyTest_1$ASVabundance)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV1_EU8 <- sort(unique(ASV_1_EU_8_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV1_EU8)){ #loop over number of unique meters (i.e. 10)
  ASV1meansAcrossEUs[5,which(colnames(ASV1meansAcrossEUs)==
                               paste("mean_",uniqueMetersASV1_EU8[i], sep=""))] <- mean(ASV_1_EU_8_tidyTest_1
                                                                                        [which(ASV_1_EU_8_tidyTest_1$Meter==uniqueMetersASV1_EU8[i]), 17]) #mean for meter 10
}

########## EU_53S ###########
ASV_1_EU_53S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[1]) %>%
  filter(EU == "EU_53S") 
head(ASV_1_EU_53S_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV1meansAcrossEUs$EU[6] <- unique(ASV_1_EU_53S_tidyTest_1$EU)
ASV1meansAcrossEUs$ASV_EUmean[6] <- mean(ASV_1_EU_53S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV1meansAcrossEUs$ASV_EUsd[6] <- sd(ASV_1_EU_53S_tidyTest_1$ASVabundance)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV1_EU53S <- sort(unique(ASV_1_EU_53S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV1_EU53S)){ #loop over number of unique meters (i.e. 10)
  ASV1meansAcrossEUs[6,which(colnames(ASV1meansAcrossEUs)==
                               paste("mean_",uniqueMetersASV1_EU53S[i], sep=""))] <- mean(ASV_1_EU_53S_tidyTest_1
                                                                                          [which(ASV_1_EU_53S_tidyTest_1$Meter==uniqueMetersASV1_EU53S[i]), 17]) #mean for meter 10
}

############## CHECKING FIRST FOR LOOP THAT CREATED ASV1meansAcrossEUs ###############
###### Making sure the for loop above works ######
index <- which(ASV_1_EU_10_tidyTest_1$Meter==uniqueMetersASV1_EU10[1]) #these are the rows
ASV_1_EU_10_tidyTest_1[index,] #this shows that it does in fact pull out all the transects for meter 10
ASV_1_EU_10_tidyTest_1[index, 17] #gives abundances. 
mean(ASV_1_EU_10_tidyTest_1[index, 17]) #okay! this is what we want!
# So, putting it all together:
mean(ASV_1_EU_10_tidyTest_1[which(ASV_1_EU_10_tidyTest_1$Meter==uniqueMetersASV1_EU10[1]), 17])

# Does indexing for column name work?
ASV1meansAcrossEUs[,2] # Meter==uniqueMetersASV1_EU10[i] #we want the column to be whatever is the name of the meter
paste(uniqueMetersASV1_EU10[1])
whichColumn <- paste("mean_",uniqueMetersASV1_EU10[1], sep="")
ASV1meansAcrossEUs[1,]
colnames(ASV1meansAcrossEUs)
which(colnames(ASV1meansAcrossEUs)=="mean_50")
ASV1meansAcrossEUs[1,which(colnames(ASV1meansAcrossEUs)==paste("mean_",uniqueMeters[9], sep=""))]
ASV1meansAcrossEUs[1,9]

# Check a few values:
check1 <- ASV_1_EU_10_tidyTest_1 %>%
  filter(Meter == "10")
checkAvg <- (28+1+468)/3
checkAvg ==  ASV1meansAcrossEUs[1,1]
mean(ASV_1_EU_10_tidyTest_1$ASVabundance) == ASV1meansAcrossEUs[1,12]
######################################################################

############ Z-SCORES (STILL FOR ASV1 ACROSS ALL EUS) ###########

# THIS CREATES A SEPARATE DATAFRAME FOR THE Z-SCORE CALCULATIONS 
ZscoresEUs_ASV1 <- as.data.frame(matrix(nrow=6, ncol=13)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
colnames(ZscoresEUs_ASV1) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100", "EU", "ASV_EUmean", "ASV_EUsd")
ZscoresEUs_ASV1[,11:13] <- ASV1meansAcrossEUs[,11:13]

# Now do vectorized calculation to take all of the mean ASV values (at each meter in each EU), minus mean in that EU, all divided by stdev
ZscoresEUs_ASV1[,1:10] <- ((ASV1meansAcrossEUs[,1:10]- ASV1meansAcrossEUs[,12])/ASV1meansAcrossEUs[,13])
# Check a few
ZscoresEUs_ASV1[2,2] == ((ASV1meansAcrossEUs[2,2] - ASV1meansAcrossEUs[2,12])/ASV1meansAcrossEUs[2,13])
ZscoresEUs_ASV1[5,10] == ((ASV1meansAcrossEUs[5,10] - ASV1meansAcrossEUs[5,12])/ASV1meansAcrossEUs[5,13])

# Finally, need to make this in long form so that it works below:
ZscoresEUs_ASV1_longer <- ZscoresEUs_ASV1 %>% 
  pivot_longer(cols=`10`:`100`, names_to="Meter", values_to = "abundZ_score") 
ZscoresEUs_ASV1_longer$Meter <- as.numeric(ZscoresEUs_ASV1_longer$Meter)

# View(ZscoresEUs_ASV1_longer) #now there are 60 values, just as we wanted!
##########################################################################
# 3.  FITTING LOGISTIC CURVES
##########################################################################

# Trying first with just ASV 1
logFit_ASV1_f <- log.fitdiffAbundFunct(y= ZscoresEUs_ASV1_longer$abundZ_score, x= ZscoresEUs_ASV1_longer$Meter, ASVnames="ASV_1")
logFit_ASV1_f

#

##########################################################################
# 4. STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################

DeseqResultsMini <- DeseqResults[,c(8,14)]  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
postUbiqTaxTab <- taxtable_outta_ps(ITS_postUbiqASVsPlus1.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == TRUE)) #287
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE)) #126 ASVs are NOT differentially abundant?
false_index <- which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE) #also 126
# Now, construct a dataframe with the taxa that were not differentially abundant 
notDA_taxTab <- postUbiqTaxTab[false_index,] #get a taxonomy tab with ONLY the non-differentially abundant ASVs
notDA_taxTab$Habitat <- "AremainingASVs" #make a habitat column that labels these as NOT differentially abundant. A in front so that would be first
# in ggplot for ease.
# View(notDA_taxTab)
colnames(notDA_taxTab)
notDA_taxTabMini <- notDA_taxTab[,c(2,8)] #keep only phylum and habitat to match DeseqResultsMini
DAphylumAll <- rbind(DeseqResultsMini, notDA_taxTabMini) #this has ASV name, phylum, and whether diff abundant for ALL ASVs in postUbiquity analysis
# What are the phyla breakdown here?
# so effectively, we want to get numbers in each phyla in each group. Should just be able to plot this?

diffAbund_ITS_stackedBarplotPhyla <- ggplot(DAphylumAll, aes(fill=Habitat, x=Phylum)) + 
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

# A few checks to make sure that the counting above is working as expected
# Ascomycota
length(which(DAphylumAll$Phylum=="p__Ascomycota")) #249
asco_index <- which(DAphylumAll$Phylum=="p__Ascomycota")
length(which(DAphylumAll[asco_index,]$Habitat=="forest")) #56 forest specialists within Ascomycota
length(which(DAphylumAll[asco_index,]$Habitat=="patch")) #120 patch specialists within Ascomycota
length(which(DAphylumAll[asco_index,]$Habitat=="AremainingASVs")) #73 remaining ASVs within Ascomycota
(56+ 120 + 73) == length(which(DAphylumAll$Phylum=="p__Ascomycota"))
# this looks correct on the plot too!

# Construct two-paneled figure with 16S and ITS differentially abundant stacked barcharts side by side
# Load in previously made 16S figure (made in "EdgeEffectsbyASV_allSites.R")
load(file="RobjectsSaved/diffAbund_16S_stackedBarplotPhyla_plot")
# quartz()
grid.arrange(diffAbund_16S_stackedBarplotPhyla, diffAbund_ITS_stackedBarplotPhyla, nrow=2)

########################
# PLOT OF GLOMEROMYCOTA
########################

# Within the plot above, glomeromycota are among the most interesting groups.
# Here, makes a plot of the different relative abundances of 

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
ITS_postUbiqNOEdge.ps <- subset_samples(ITS_postUbiquity.ps, Habitat != "edge") #remove edge so we can compare patch versus forest

#ITS_postUbiqNOEdge.ps.phylum.glom <-  tax_glom(ITS_postUbiqNOEdge.ps, taxrank = "Phylum") 
#tax_table(ITS_postUbiqNOEdge.ps.phylum.glom) #8 phyla
#sample_data(ITS_postUbiqNOEdge.ps.phylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
# relabunpostUbiqNOEdge.phyla.0 <- transform_sample_counts(ITS_postUbiqNOEdge.ps.phylum.glom, function(x) x / sum(x) )
relabunpostUbiqNOEdge.allASVs <- transform_sample_counts(ITS_postUbiqNOEdge.ps, function(x) x / sum(x) )

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunpostUbiqNOEdge.allASVs <- merge_samples(relabunpostUbiqNOEdge.allASVs, group = c("Habitat"))
sample_data(relabunpostUbiqNOEdge.allASVs) #shows that we still have samples from each EU, biocrust, and extcontrol (water)
# Meter did an averaging thing; can just ignore it

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunpostUbiqNOEdge.allASVs.2 <- transform_sample_counts(relabunpostUbiqNOEdge.allASVs, function(x) x / sum(x))
sample_data(relabunpostUbiqNOEdge.allASVs.2)

# 
relabunpostUbiqNOEdge.allASVs.df <-psmelt(relabunpostUbiqNOEdge.allASVs.2)
head(relabunpostUbiqNOEdge.allASVs.df) #
dim(relabunpostUbiqNOEdge.allASVs.df) 

glomero_index <- which(relabunpostUbiqNOEdge.allASVs.df$Phylum=="p__Glomeromycota")
glomeroRelAbund <- relabunpostUbiqNOEdge.allASVs.df[glomero_index, ]
#View(glomeroRelAbund)
colnames(glomeroRelAbund)[3] <- "Relative_abundance"
colnames(glomeroRelAbund)[2] <- "Habitat_type"

# View(relabunpostUbiqNOEdge.phyla.df)

# Make boxplot of relative abundances of each ASV in the phylum
glomeroBoxPlot <- ggplot(glomeroRelAbund, aes(x=Habitat_type, y=Relative_abundance, fill= Habitat_type)) +
  geom_boxplot() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  labs(title="Relative abundance of Glomeromycota ASVs", x="Habitat Type", y = "Relative abundance")
quartz()
glomeroBoxPlot

# Get exact abundances of each phyla 
#colnames(relabunpostUbiqNOEdge.phyla.df)
#relabunpostUbiqNOEdge.phyla.df[,2] #Now just all patch and forest!


# This plot is just to check that numbers look right!
#relabunpostUbiqNOEdgePlot <- ggplot(data=relabunpostUbiqNOEdge.phyla.df, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
#quartz()
relabunpostUbiqNOEdgePlot + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Fungal phyla comprising at least 0.5% of total abundance (median and post-ubiquity)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#top_99.5p_phyla <- relabunpostUbiqNOEdge.phyla.df %>%
group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
