# DiffAbundFungi
# Aug 14, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect FUNGAL (ITS) ASVs using differential abundance 
# analyses (DESeq2) and makes relevant plots. It comes BEFORE fungalZscores.

# This script:
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. ITS_postUbiquity.ps)).
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential 
# abundance breakdown by forest, patch, and non-differentially abundant microbes 

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("reshape2")
library("growthcurver") #fits logistic curve to the abundance (z-score) data

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
head(diffAbunDat_tidy_FUNGI_forest)

ASVmeterZAbunds_fungi <- vector("list", length(ASVnamesDA_FUNGI_forest)) #pre-allocate space for each ASV's abundance for each ASV


############ TRYING OUT A FEW METHODS OF Z-SCORES ############
# test with ASV_1 in EU_10
# PREALLOCATE DATAFRAME THAT I CAN USE FOR ALL VALUES OF ASV_1 ACROSS ALL EUs
ASV1meansAcrossEUs <- as.data.frame(matrix(nrow=6, ncol=14)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
head(ASV1meansAcrossEUs)
colnames(ASV1meansAcrossEUs) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name")
#names(new) <- paste(ASVnamesDA_FUNGI_forest[1],"acrossEUs", sep= "")
head(ASV1meansAcrossEUs)

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
ASV1meansAcrossEUs$ASV_name[1] <- unique(ASV_1_EU_10_tidyTest_1$ASV_name)
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
ASV1meansAcrossEUs$ASV_name[2] <- unique(ASV_1_EU_52_tidyTest_1$ASV_name)
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
ASV1meansAcrossEUs$ASV_name[3] <- unique(ASV_1_EU_53N_tidyTest_1$ASV_name)
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
ASV1meansAcrossEUs$ASV_name[4] <- unique(ASV_1_EU_54S_tidyTest_1$ASV_name)

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
ASV1meansAcrossEUs$ASV_name[5] <- unique(ASV_1_EU_8_tidyTest_1$ASV_name)

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
ASV1meansAcrossEUs$ASV_name[6] <- unique(ASV_1_EU_53S_tidyTest_1$ASV_name)
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
ASV1meansAcrossEUs[1,which(colnames(ASV1meansAcrossEUs)==paste("mean_",uniqueMetersASV1_EU10[9], sep=""))]
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

##### Trying a new following explanation here: https://rpubs.com/angelov/growthcurver ######
# First, because these are growth plots, add a 1 to all of the z-scores for correct fit
ZscoresEUs_ASV1_longerPlus1 <- ZscoresEUs_ASV1_longer
ZscoresEUs_ASV1_longerPlus1$abundZ_score <- ZscoresEUs_ASV1_longer$abundZ_score + 1
ggplot(ZscoresEUs_ASV1_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=0.7) +
  theme_bw()

modelASV1plus1.f <- growthcurver::SummarizeGrowth(data_t=ZscoresEUs_ASV1_longerPlus1$Meter, data_n=ZscoresEUs_ASV1_longerPlus1$abundZ_score, bg_correct = "none")
modelASV1plus1.f$vals #gives all of the values
predict(modelASV1plus1.f$model) # gives you the predicted abundance values (according to the model)
str(modelASV1plus1.f)
modelASV1plus1.f$vals$r #this is growth rate constant NOT r-squared
modelASV1plus1.f$vals$sigma #0.4752275 -- the smaller the better!
# sigma is a measure of the goodnesss of fit of the parameters of the logistic equation for the data; 
# it is the residual standard error from the nonlinear regression model. Smaller sigma values indicate
# a better fit of the logistic curve to the data than larger values.
plot(predict(modelASV1plus1.f$model) )
plot(modelASV1plus1.f$model$m$fitted())
modelASV1plus1.f$model$m$fitted()

### Plotting ##
# Base R
plot(modelASV1plus1.f) #ugly!

# ggplot
modelASV1plus1.f$vals$t_mid #inflection point
modelASV1plus1.f$vals$sigma #0.4772041
p1 <- ggplot(ZscoresEUs_ASV1_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=1.3) + theme_bw()
p1


# Adding predicted values
df.predicted <- data.frame(Meter = ZscoresEUs_ASV1_longerPlus1$Meter, pred.Zabund = modelASV1plus1.f$model$m$fitted())
p2 <- p1 + geom_line(data=df.predicted, aes(y=pred.Zabund), color="red", size= 4) + ggtitle("Mortierella humilis (fungi)") +
  geom_point(x=modelASV1plus1.f$vals$t_mid, y = 0.78, color= "red", size=6) + #Here I just guessed a y based on how it looked!
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=2) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) #make it so all meters show on x-axis!
quartz()
p2

diffAbunDat_tidy_FUNGI[1,] #getting info for ASV_1 to add to ggtitle above

df.predicted[c(1:10),] == df.predicted[c(11:20),] #this data frame is longer than it needs to be, but that's okay!
df.predicted[which(df.predicted$Meter==10),2]  # At 10 meters, abundance prediction (+ 1) is 0.4173783 on the y-axis
# Depth = 80% of distance between edge and abundance at end point, or 
modelASV1plus1.f$vals$t_mid - 10 #26.79592
26.79592*.2 # 5.359184 is about at barrier

################# TRYING WITH ASV_416!

# test with ASV_416 in EU_10
# PREALLOCATE DATAFRAME THAT I CAN USE FOR ALL VALUES OF ASV_1 ACROSS ALL EUs
ASV_416_meansAcrossEUs <- as.data.frame(matrix(nrow=6, ncol=14)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
head(ASV_416_meansAcrossEUs)
colnames(ASV_416_meansAcrossEUs) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name")
#names(new) <- paste(ASVnamesDA_FUNGI_forest[1],"acrossEUs", sep= "")
head(ASV_416_meansAcrossEUs)

######### EU_10 #######
# First, isolate 1 ASV in EU_10
ASV_416_EU_10_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[2]) %>% #second is ASV 416
  filter(EU == "EU_10") #at this point, have only ASV_416 in EU_10
head(ASV_416_EU_10_tidyTest_1)
head(ASV_416_meansAcrossEUs)

# Get ASV mean and stdev across all points on the transect
ASV_416_meansAcrossEUs$EU[1] <- unique(ASV_416_EU_10_tidyTest_1$EU)
ASV_416_meansAcrossEUs$ASV_EUmean[1] <- mean(ASV_416_EU_10_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU 10
ASV_416_meansAcrossEUs$ASV_EUsd[1] <- sd(ASV_416_EU_10_tidyTest_1$ASVabundance)
ASV_416_meansAcrossEUs$ASV_name[1] <- unique(ASV_416_EU_10_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_416_EU10 <- sort(unique(ASV_416_EU_10_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_416_EU10)){ #loop over number of unique meters (i.e. 10)
  ASV_416_meansAcrossEUs[1,which(colnames(ASV_416_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_416_EU10[i], sep=""))] <- mean(ASV_416_EU_10_tidyTest_1
                                                                                                [which(ASV_416_EU_10_tidyTest_1$Meter==uniqueMetersASV_416_EU10[i]), 17]) #mean for meter 10
}

############### Repeat this for all of the other EUs ########################
unique(diffAbunDat_tidy_FUNGI$EU)
######### EU_52 #########
ASV_416_EU_52_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[2]) %>%
  filter(EU == "EU_52") 
head(ASV_416_EU_52_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_416_meansAcrossEUs$EU[2] <- unique(ASV_416_EU_52_tidyTest_1$EU)
ASV_416_meansAcrossEUs$ASV_EUmean[2] <- mean(ASV_416_EU_52_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_416_meansAcrossEUs$ASV_EUsd[2] <- sd(ASV_416_EU_52_tidyTest_1$ASVabundance)
ASV_416_meansAcrossEUs$ASV_name[2] <- unique(ASV_416_EU_52_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_416_EU52 <- sort(unique(ASV_416_EU_52_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_416_EU52)){ #loop over number of unique meters (i.e. 10)
  ASV_416_meansAcrossEUs[2,which(colnames(ASV_416_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_416_EU52[i], sep=""))] <- mean(ASV_416_EU_52_tidyTest_1
                                                                                                [which(ASV_416_EU_52_tidyTest_1$Meter==uniqueMetersASV_416_EU52[i]), 17]) #mean for meter 10
}

###########  EU_53N ########## 
ASV_416_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[2]) %>%
  filter(EU == "EU_53N") 
head(ASV_416_EU_53N_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_416_meansAcrossEUs$EU[3] <- unique(ASV_416_EU_53N_tidyTest_1$EU)
ASV_416_meansAcrossEUs$ASV_EUmean[3] <- mean(ASV_416_EU_53N_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_416_meansAcrossEUs$ASV_EUsd[3] <- sd(ASV_416_EU_53N_tidyTest_1$ASVabundance)
ASV_416_meansAcrossEUs$ASV_name[3] <- unique(ASV_416_EU_53N_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_416_EU53N <- sort(unique(ASV_416_EU_53N_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_416_EU53N)){ #loop over number of unique meters (i.e. 10)
  ASV_416_meansAcrossEUs[3,which(colnames(ASV_416_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_416_EU53N[i], sep=""))] <- mean(ASV_416_EU_53N_tidyTest_1
                                                                                                 [which(ASV_416_EU_53N_tidyTest_1$Meter==uniqueMetersASV_416_EU53N[i]), 17]) #mean for meter 10
}

###########  EU_54S ########## 
ASV_416_EU_54S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[2]) %>%
  filter(EU == "EU_54S") 
head(ASV_416_EU_54S_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_416_meansAcrossEUs$EU[4] <- unique(ASV_416_EU_54S_tidyTest_1$EU)
ASV_416_meansAcrossEUs$ASV_EUmean[4] <- mean(ASV_416_EU_54S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_416_meansAcrossEUs$ASV_EUsd[4] <- sd(ASV_416_EU_54S_tidyTest_1$ASVabundance)
ASV_416_meansAcrossEUs$ASV_name[4] <- unique(ASV_416_EU_54S_tidyTest_1$ASV_name)

# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_416_EU54S <- sort(unique(ASV_416_EU_54S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_416_EU54S)){ #loop over number of unique meters (i.e. 10)
  ASV_416_meansAcrossEUs[4,which(colnames(ASV_416_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_416_EU54S[i], sep=""))] <- mean(ASV_416_EU_54S_tidyTest_1
                                                                                                 [which(ASV_416_EU_54S_tidyTest_1$Meter==uniqueMetersASV_416_EU54S[i]), 17]) #mean for meter 10
}

##########  EU_8 ########## 
ASV_416_EU_8_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[2]) %>%
  filter(EU == "EU_8") 
head(ASV_416_EU_8_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_416_meansAcrossEUs$EU[5] <- unique(ASV_416_EU_8_tidyTest_1$EU)
ASV_416_meansAcrossEUs$ASV_EUmean[5] <- mean(ASV_416_EU_8_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_416_meansAcrossEUs$ASV_EUsd[5] <- sd(ASV_416_EU_8_tidyTest_1$ASVabundance)
ASV_416_meansAcrossEUs$ASV_name[5] <- unique(ASV_416_EU_8_tidyTest_1$ASV_name)

# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_416_EU8 <- sort(unique(ASV_416_EU_8_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_416_EU8)){ #loop over number of unique meters (i.e. 10)
  ASV_416_meansAcrossEUs[5,which(colnames(ASV_416_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_416_EU8[i], sep=""))] <- mean(ASV_416_EU_8_tidyTest_1
                                                                                               [which(ASV_416_EU_8_tidyTest_1$Meter==uniqueMetersASV_416_EU8[i]), 17]) #mean for meter 10
}

########## EU_53S ###########
ASV_416_EU_53S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[2]) %>%
  filter(EU == "EU_53S") 
head(ASV_416_EU_53S_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_416_meansAcrossEUs$EU[6] <- unique(ASV_416_EU_53S_tidyTest_1$EU)
ASV_416_meansAcrossEUs$ASV_EUmean[6] <- mean(ASV_416_EU_53S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_416_meansAcrossEUs$ASV_EUsd[6] <- sd(ASV_416_EU_53S_tidyTest_1$ASVabundance)
ASV_416_meansAcrossEUs$ASV_name[6] <- unique(ASV_416_EU_53S_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_416_EU53S <- sort(unique(ASV_416_EU_53S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_416_EU53S)){ #loop over number of unique meters (i.e. 10)
  ASV_416_meansAcrossEUs[6,which(colnames(ASV_416_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_416_EU53S[i], sep=""))] <- mean(ASV_416_EU_53S_tidyTest_1
                                                                                                 [which(ASV_416_EU_53S_tidyTest_1$Meter==uniqueMetersASV_416_EU53S[i]), 17]) #mean for meter 10
}

############## CHECKING FIRST FOR LOOP THAT CREATED ASV_416_meansAcrossEUs ###############
###### Making sure the for loop above works ######
index <- which(ASV_416_EU_10_tidyTest_1$Meter==uniqueMetersASV_416_EU10[1]) #these are the rows
ASV_416_EU_10_tidyTest_1[index,] #this shows that it does in fact pull out all the transects for meter 10
ASV_416_EU_10_tidyTest_1[index, 17] #gives abundances. 
mean(ASV_416_EU_10_tidyTest_1[index, 17]) #okay! this is what we want!
# So, putting it all together:
mean(ASV_416_EU_10_tidyTest_1[which(ASV_416_EU_10_tidyTest_1$Meter==uniqueMetersASV_416_EU10[1]), 17])

# Does indexing for column name work?
ASV_416_meansAcrossEUs[,2] # Meter==uniqueMetersASV_416_EU10[i] #we want the column to be whatever is the name of the meter
paste(uniqueMetersASV_416_EU10[1])
whichColumn <- paste("mean_",uniqueMetersASV_416_EU10[1], sep="")
ASV_416_meansAcrossEUs[1,]
colnames(ASV_416_meansAcrossEUs)
which(colnames(ASV_416_meansAcrossEUs)=="mean_50")
ASV_416_meansAcrossEUs[1,which(colnames(ASV_416_meansAcrossEUs)==paste("mean_",uniqueMetersASV_416_EU10[9], sep=""))]
ASV_416_meansAcrossEUs[1,9]

# Check a few values:
check1 <- ASV_416_EU_10_tidyTest_1 %>%
  filter(Meter == "10")
checkAvg <- (1+1+1)/3
checkAvg ==  ASV_416_meansAcrossEUs[1,1]
mean(ASV_416_EU_10_tidyTest_1$ASVabundance) == ASV_416_meansAcrossEUs[1,12]
######################################################################

######################################################################

############ Z-SCORES (STILL FOR ASV416 ACROSS ALL EUS) ###########

# THIS CREATES A SEPARATE DATAFRAME FOR THE Z-SCORE CALCULATIONS 
ZscoresEUs_ASV416 <- as.data.frame(matrix(nrow=6, ncol=13)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
colnames(ZscoresEUs_ASV416) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100", "EU", "ASV_EUmean", "ASV_EUsd")
ZscoresEUs_ASV416[,11:13] <- ASV_416_meansAcrossEUs[,11:13]

# Now do vectorized calculation to take all of the mean ASV values (at each meter in each EU), minus mean in that EU, all divided by stdev
ZscoresEUs_ASV416[,1:10] <- ((ASV_416_meansAcrossEUs[,1:10]- ASV_416_meansAcrossEUs[,12])/ASV_416_meansAcrossEUs[,13])
# Check a few
ZscoresEUs_ASV416[2,2] == ((ASV_416_meansAcrossEUs[2,2] - ASV_416_meansAcrossEUs[2,12])/ASV_416_meansAcrossEUs[2,13])
ZscoresEUs_ASV416[5,10] == ((ASV_416_meansAcrossEUs[5,10] - ASV_416_meansAcrossEUs[5,12])/ASV_416_meansAcrossEUs[5,13])

# Finally, need to make this in long form so that it works below:
ZscoresEUs_ASV416_longer <- ZscoresEUs_ASV416 %>% 
  pivot_longer(cols=`10`:`100`, names_to="Meter", values_to = "abundZ_score") 
ZscoresEUs_ASV416_longer$Meter <- as.numeric(ZscoresEUs_ASV416_longer$Meter)

# View(ZscoresEUs_ASV416_longer) #now there are 60 values, just as we wanted!
##########################################################################
# 3.  FITTING LOGISTIC CURVES
##########################################################################

# Trying first with just ASV 416-- doesn't work!
logFit_ASV416_f <- log.fitdiffAbundFunct(y= ZscoresEUs_ASV416_longer$abundZ_score, x= ZscoresEUs_ASV416_longer$Meter, ASVnames="ASV_416")
logFit_ASV416_f #Error in qr.solve(QR.B, cc) : singular matrix 'a' in solve

##### Trying a new following explanation here: https://rpubs.com/angelov/growthcurver ######
# First, because these are growth plots, add a 1 to all of the z-scores for correct fit
ZscoresEUs_ASV416_longerPlus1 <- ZscoresEUs_ASV416_longer
ZscoresEUs_ASV416_longerPlus1$abundZ_score <- ZscoresEUs_ASV416_longer$abundZ_score + 1
ggplot(ZscoresEUs_ASV416_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=0.7) +
  theme_bw()

modelASV416plus1.f <- growthcurver::SummarizeGrowth(data_t=ZscoresEUs_ASV416_longerPlus1$Meter, data_n=ZscoresEUs_ASV416_longerPlus1$abundZ_score, bg_correct = "none")
modelASV416plus1.f$vals #gives all of the values
predict(modelASV416plus1.f$model) # gives you the predicted abundance values (according to the model)
str(modelASV416plus1.f)
modelASV416plus1.f$vals$r #this is growth rate constant NOT r-squared
modelASV416plus1.f$vals$sigma #0.4752275 -- the smaller the better!
# sigma is a measure of the goodnesss of fit of the parameters of the logistic equation for the data; 
# it is the residual standard error from the nonlinear regression model. Smaller sigma values indicate
# a better fit of the logistic curve to the data than larger values.
plot(predict(modelASV416plus1.f$model) )
plot(modelASV416plus1.f$model$m$fitted())
modelASV416plus1.f$model$m$fitted()

### Plotting ##
# Base R
plot(modelASV416plus1.f) #ugly!

# ggplot
ASV_416_EU_10_tidyTest_1 #getting info for ASV_416 to add to ggtitle
modelASV416plus1.f$vals$t_mid #inflection point
modelASV416plus1.f$vals$sigma #0.4772041
p1 <- ggplot(ZscoresEUs_ASV416_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=1.3) + theme_bw()
p1
# Adding predicted values
df.predicted <- data.frame(Meter = ZscoresEUs_ASV416_longerPlus1$Meter, pred.Zabund = modelASV416plus1.f$model$m$fitted())
p2 <- p1 + geom_line(data=df.predicted, aes(y=pred.Zabund), color="red", size= 4) + ggtitle("Mortierella sp. (Mortierellomycota)") +
  #geom_point(x=modelASV416plus1.f$vals$t_mid, y = 0.78, color= "red", size=6) + #Here I just guessed a y based on how it looked!
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=2) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) #make it so all meters show on x-axis!
quartz()
p2



df.predicted[c(1:10),] == df.predicted[c(11:20),] #this data frame is longer than it needs to be, but that's okay!
df.predicted[which(df.predicted$Meter==10),2]  # At 10 meters, abundance prediction (+ 1) is 0.4173783 on the y-axis
# Depth = 80% of distance between edge and abundance at end point, or 
modelASV416plus1.f$vals$t_mid - 10 #26.83405




diffAbunDat_tidy_FUNGI[3,]
######################################################
################# ################# TRYING WITH ASV_792!

# test with ASV_792 in EU_10
# PREALLOCATE DATAFRAME THAT I CAN USE FOR ALL VALUES OF ASV_1 ACROSS ALL EUs
ASV_792_meansAcrossEUs <- as.data.frame(matrix(nrow=6, ncol=14)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
head(ASV_792_meansAcrossEUs)
colnames(ASV_792_meansAcrossEUs) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name")
#names(new) <- paste(ASVnamesDA_FUNGI_forest[3],"acrossEUs", sep= "")
head(ASV_792_meansAcrossEUs)

######### EU_10 #######
# First, isolate 1 ASV in EU_10
ASV_792_EU_10_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>% #second is ASV 792
  filter(EU == "EU_10") #at this point, have only ASV_792 in EU_10
head(ASV_792_EU_10_tidyTest_1)
head(ASV_792_meansAcrossEUs)

# Get ASV mean and stdev across all points on the transect
ASV_792_meansAcrossEUs$EU[1] <- unique(ASV_792_EU_10_tidyTest_1$EU)
ASV_792_meansAcrossEUs$ASV_EUmean[1] <- mean(ASV_792_EU_10_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU 10
ASV_792_meansAcrossEUs$ASV_EUsd[1] <- sd(ASV_792_EU_10_tidyTest_1$ASVabundance)
ASV_792_meansAcrossEUs$ASV_name[1] <- unique(ASV_792_EU_10_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_792_EU10 <- sort(unique(ASV_792_EU_10_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_792_EU10)){ #loop over number of unique meters (i.e. 10)
  ASV_792_meansAcrossEUs[1,which(colnames(ASV_792_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_792_EU10[i], sep=""))] <- mean(ASV_792_EU_10_tidyTest_1
                                                                                                [which(ASV_792_EU_10_tidyTest_1$Meter==uniqueMetersASV_792_EU10[i]), 17]) #mean for meter 10
}

############### Repeat this for all of the other EUs ########################
unique(diffAbunDat_tidy_FUNGI$EU)
######### EU_52 #########
ASV_792_EU_52_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_52") 
head(ASV_792_EU_52_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_792_meansAcrossEUs$EU[2] <- unique(ASV_792_EU_52_tidyTest_1$EU)
ASV_792_meansAcrossEUs$ASV_EUmean[2] <- mean(ASV_792_EU_52_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_792_meansAcrossEUs$ASV_EUsd[2] <- sd(ASV_792_EU_52_tidyTest_1$ASVabundance)
ASV_792_meansAcrossEUs$ASV_name[2] <- unique(ASV_792_EU_52_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_792_EU52 <- sort(unique(ASV_792_EU_52_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_792_EU52)){ #loop over number of unique meters (i.e. 10)
  ASV_792_meansAcrossEUs[2,which(colnames(ASV_792_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_792_EU52[i], sep=""))] <- mean(ASV_792_EU_52_tidyTest_1
                                                                                                [which(ASV_792_EU_52_tidyTest_1$Meter==uniqueMetersASV_792_EU52[i]), 17]) #mean for meter 10
}

###########  EU_53N ########## 
ASV_792_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_53N") 
head(ASV_792_EU_53N_tidyTest_1)
ASV_792_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_53N") 
head(ASV_792_EU_53N_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_792_meansAcrossEUs$EU[3] <- unique(ASV_792_EU_53N_tidyTest_1$EU)
ASV_792_meansAcrossEUs$ASV_EUmean[3] <- mean(ASV_792_EU_53N_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_792_meansAcrossEUs$ASV_EUsd[3] <- sd(ASV_792_EU_53N_tidyTest_1$ASVabundance)
ASV_792_meansAcrossEUs$ASV_name[3] <- unique(ASV_792_EU_53N_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_792_EU53N <- sort(unique(ASV_792_EU_53N_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_792_EU53N)){ #loop over number of unique meters (i.e. 10)
  ASV_792_meansAcrossEUs[3,which(colnames(ASV_792_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_792_EU53N[i], sep=""))] <- mean(ASV_792_EU_53N_tidyTest_1
                                                                                                 [which(ASV_792_EU_53N_tidyTest_1$Meter==uniqueMetersASV_792_EU53N[i]), 17]) #mean for meter 10
}

###########  EU_54S ########## 
ASV_792_EU_54S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_54S") 
head(ASV_792_EU_54S_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_792_meansAcrossEUs$EU[4] <- unique(ASV_792_EU_54S_tidyTest_1$EU)
ASV_792_meansAcrossEUs$ASV_EUmean[4] <- mean(ASV_792_EU_54S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_792_meansAcrossEUs$ASV_EUsd[4] <- sd(ASV_792_EU_54S_tidyTest_1$ASVabundance)
ASV_792_meansAcrossEUs$ASV_name[4] <- unique(ASV_792_EU_54S_tidyTest_1$ASV_name)

# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_792_EU54S <- sort(unique(ASV_792_EU_54S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_792_EU54S)){ #loop over number of unique meters (i.e. 10)
  ASV_792_meansAcrossEUs[4,which(colnames(ASV_792_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_792_EU54S[i], sep=""))] <- mean(ASV_792_EU_54S_tidyTest_1
                                                                                                 [which(ASV_792_EU_54S_tidyTest_1$Meter==uniqueMetersASV_792_EU54S[i]), 17]) #mean for meter 10
}

##########  EU_8 ########## 
ASV_792_EU_8_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_8") 
head(ASV_792_EU_8_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_792_meansAcrossEUs$EU[5] <- unique(ASV_792_EU_8_tidyTest_1$EU)
ASV_792_meansAcrossEUs$ASV_EUmean[5] <- mean(ASV_792_EU_8_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_792_meansAcrossEUs$ASV_EUsd[5] <- sd(ASV_792_EU_8_tidyTest_1$ASVabundance)
ASV_792_meansAcrossEUs$ASV_name[5] <- unique(ASV_792_EU_8_tidyTest_1$ASV_name)

# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_792_EU8 <- sort(unique(ASV_792_EU_8_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_792_EU8)){ #loop over number of unique meters (i.e. 10)
  ASV_792_meansAcrossEUs[5,which(colnames(ASV_792_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_792_EU8[i], sep=""))] <- mean(ASV_792_EU_8_tidyTest_1
                                                                                               [which(ASV_792_EU_8_tidyTest_1$Meter==uniqueMetersASV_792_EU8[i]), 17]) #mean for meter 10
}

########## EU_53S ###########
ASV_792_EU_53S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_53S") 
head(ASV_792_EU_53S_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_792_meansAcrossEUs$EU[6] <- unique(ASV_792_EU_53S_tidyTest_1$EU)
ASV_792_meansAcrossEUs$ASV_EUmean[6] <- mean(ASV_792_EU_53S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_792_meansAcrossEUs$ASV_EUsd[6] <- sd(ASV_792_EU_53S_tidyTest_1$ASVabundance)
ASV_792_meansAcrossEUs$ASV_name[6] <- unique(ASV_792_EU_53S_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_792_EU53S <- sort(unique(ASV_792_EU_53S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_792_EU53S)){ #loop over number of unique meters (i.e. 10)
  ASV_792_meansAcrossEUs[6,which(colnames(ASV_792_meansAcrossEUs)==
                                   paste("mean_",uniqueMetersASV_792_EU53S[i], sep=""))] <- mean(ASV_792_EU_53S_tidyTest_1
                                                                                                 [which(ASV_792_EU_53S_tidyTest_1$Meter==uniqueMetersASV_792_EU53S[i]), 17]) #mean for meter 10
}

############## CHECKING FIRST FOR LOOP THAT CREATED ASV_792_meansAcrossEUs ###############
###### Making sure the for loop above works ######
index <- which(ASV_792_EU_10_tidyTest_1$Meter==uniqueMetersASV_792_EU10[1]) #these are the rows
ASV_792_EU_10_tidyTest_1[index,] #this shows that it does in fact pull out all the transects for meter 10
ASV_792_EU_10_tidyTest_1[index, 17] #gives abundances. 
mean(ASV_792_EU_10_tidyTest_1[index, 17]) #okay! this is what we want!
# So, putting it all together:
mean(ASV_792_EU_10_tidyTest_1[which(ASV_792_EU_10_tidyTest_1$Meter==uniqueMetersASV_792_EU10[1]), 17])

# Does indexing for column name work?
ASV_792_meansAcrossEUs[,2] # Meter==uniqueMetersASV_792_EU10[i] #we want the column to be whatever is the name of the meter
paste(uniqueMetersASV_792_EU10[1])
whichColumn <- paste("mean_",uniqueMetersASV_792_EU10[1], sep="")
ASV_792_meansAcrossEUs[1,]
colnames(ASV_792_meansAcrossEUs)
which(colnames(ASV_792_meansAcrossEUs)=="mean_50")
ASV_792_meansAcrossEUs[1,which(colnames(ASV_792_meansAcrossEUs)==paste("mean_",uniqueMetersASV_792_EU10[9], sep=""))]
ASV_792_meansAcrossEUs[1,9]

# Check a few values:
check1 <- ASV_792_EU_10_tidyTest_1 %>%
  filter(Meter == "10")
checkAvg <- (1+1+1)/3
checkAvg ==  ASV_792_meansAcrossEUs[1,1]
mean(ASV_792_EU_10_tidyTest_1$ASVabundance) == ASV_792_meansAcrossEUs[1,12]
######################################################################

######################################################################

############ Z-SCORES (STILL FOR ASV792 ACROSS ALL EUS) ###########

# THIS CREATES A SEPARATE DATAFRAME FOR THE Z-SCORE CALCULATIONS 
ZscoresEUs_ASV792 <- as.data.frame(matrix(nrow=6, ncol=13)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
colnames(ZscoresEUs_ASV792) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100", "EU", "ASV_EUmean", "ASV_EUsd")
ZscoresEUs_ASV792[,11:13] <- ASV_792_meansAcrossEUs[,11:13]

# Now do vectorized calculation to take all of the mean ASV values (at each meter in each EU), minus mean in that EU, all divided by stdev
ZscoresEUs_ASV792[,1:10] <- ((ASV_792_meansAcrossEUs[,1:10]- ASV_792_meansAcrossEUs[,12])/ASV_792_meansAcrossEUs[,13])
# Check a few
ZscoresEUs_ASV792[2,2] == ((ASV_792_meansAcrossEUs[2,2] - ASV_792_meansAcrossEUs[2,12])/ASV_792_meansAcrossEUs[2,13])
ZscoresEUs_ASV792[5,10] == ((ASV_792_meansAcrossEUs[5,10] - ASV_792_meansAcrossEUs[5,12])/ASV_792_meansAcrossEUs[5,13])

# Finally, need to make this in long form so that it works below:
ZscoresEUs_ASV792_longer <- ZscoresEUs_ASV792 %>% 
  pivot_longer(cols=`10`:`100`, names_to="Meter", values_to = "abundZ_score") 
ZscoresEUs_ASV792_longer$Meter <- as.numeric(ZscoresEUs_ASV792_longer$Meter)

# View(ZscoresEUs_ASV792_longer) #now there are 60 values, just as we wanted!
##########################################################################
# 3.  FITTING LOGISTIC CURVES
##########################################################################

# Trying first with just ASV 792-- doesn't work!
logFit_ASV792_f <- log.fitdiffAbundFunct(y= ZscoresEUs_ASV792_longer$abundZ_score, x= ZscoresEUs_ASV792_longer$Meter, ASVnames="ASV_792")
logFit_ASV792_f #Error in qr.solve(QR.B, cc) : singular matrix 'a' in solve

##### Trying a new following explanation here: https://rpubs.com/angelov/growthcurver ######
# First, because these are growth plots, add a 1 to all of the z-scores for correct fit
ZscoresEUs_ASV792_longerPlus1 <- ZscoresEUs_ASV792_longer
ZscoresEUs_ASV792_longerPlus1$abundZ_score <- ZscoresEUs_ASV792_longer$abundZ_score + 1
ggplot(ZscoresEUs_ASV792_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=0.7) +
  theme_bw()

modelASV792plus1.f <- growthcurver::SummarizeGrowth(data_t=ZscoresEUs_ASV792_longerPlus1$Meter, data_n=ZscoresEUs_ASV792_longerPlus1$abundZ_score, bg_correct = "none")
modelASV792plus1.f$vals #gives all of the values
predict(modelASV792plus1.f$model) # gives you the predicted abundance values (according to the model)
str(modelASV792plus1.f)
modelASV792plus1.f$vals$r #this is growth rate constant NOT r-squared
modelASV792plus1.f$vals$sigma #0.4752275 -- the smaller the better!
# sigma is a measure of the goodnesss of fit of the parameters of the logistic equation for the data; 
# it is the residual standard error from the nonlinear regression model. Smaller sigma values indicate
# a better fit of the logistic curve to the data than larger values.
plot(predict(modelASV792plus1.f$model) )
plot(modelASV792plus1.f$model$m$fitted())
modelASV792plus1.f$model$m$fitted()

### Plotting ##
# Base R
plot(modelASV792plus1.f) #ugly!

# ggplot
diffAbunDat_tidy_FUNGI[2,] #getting info for ASV_792 to add to ggtitle
modelASV792plus1.f$vals$t_mid #inflection point
modelASV792plus1.f$vals$sigma #0.4772041
p1 <- ggplot(ZscoresEUs_ASV792_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=1.3) + theme_bw()
p1
# Adding predicted values
ASV_792_EU_10_tidyTest_1 #look for name for below
df.predicted <- data.frame(Meter = ZscoresEUs_ASV792_longerPlus1$Meter, pred.Zabund = modelASV792plus1.f$model$m$fitted())
p2 <- p1 + geom_line(data=df.predicted, aes(y=pred.Zabund), color="red", size= 4) + ggtitle("Phialocephala fortini (Ascomycota)") +
  #geom_point(x=modelASV792plus1.f$vals$t_mid, y = 0.78, color= "red", size=6) + #Here I just guessed a y based on how it looked!
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=2) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) #make it so all meters show on x-axis!
quartz()
p2



df.predicted[c(1:10),] == df.predicted[c(11:20),] #this data frame is longer than it needs to be, but that's okay!
df.predicted[which(df.predicted$Meter==10),2]  # At 10 meters, abundance prediction (+ 1) is 0.4173783 on the y-axis
# Depth = 80% of distance between edge and abundance at end point, or 
modelASV792plus1.f$vals$t_mid - 10 

###################################################
#################  TRYING WITH ASV_23! #################

# test with ASV_23 in EU_10
# PREALLOCATE DATAFRAME THAT I CAN USE FOR ALL VALUES OF ASV_1 ACROSS ALL EUs
ASV_23_meansAcrossEUs <- as.data.frame(matrix(nrow=6, ncol=14)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
head(ASV_23_meansAcrossEUs)
colnames(ASV_23_meansAcrossEUs) <- c("mean_10", "mean_20", "mean_30", "mean_40", "mean_50", "mean_60", "mean_70", "mean_80", "mean_90", "mean_100", "EU", "ASV_EUmean", "ASV_EUsd", "ASV_name")
#names(new) <- paste(ASVnamesDA_FUNGI_forest[3],"acrossEUs", sep= "")
head(ASV_23_meansAcrossEUs)

######### EU_10 #######
# First, isolate 1 ASV in EU_10
ASV_23_EU_10_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[4]) %>% #second is ASV 23
  filter(EU == "EU_10") #at this point, have only ASV_23 in EU_10
head(ASV_23_EU_10_tidyTest_1)
head(ASV_23_meansAcrossEUs)

# Get ASV mean and stdev across all points on the transect
ASV_23_meansAcrossEUs$EU[1] <- unique(ASV_23_EU_10_tidyTest_1$EU)
ASV_23_meansAcrossEUs$ASV_EUmean[1] <- mean(ASV_23_EU_10_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU 10
ASV_23_meansAcrossEUs$ASV_EUsd[1] <- sd(ASV_23_EU_10_tidyTest_1$ASVabundance)
ASV_23_meansAcrossEUs$ASV_name[1] <- unique(ASV_23_EU_10_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_23_EU10 <- sort(unique(ASV_23_EU_10_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_23_EU10)){ #loop over number of unique meters (i.e. 10)
  ASV_23_meansAcrossEUs[1,which(colnames(ASV_23_meansAcrossEUs)==
                                  paste("mean_",uniqueMetersASV_23_EU10[i], sep=""))] <- mean(ASV_23_EU_10_tidyTest_1
                                                                                              [which(ASV_23_EU_10_tidyTest_1$Meter==uniqueMetersASV_23_EU10[i]), 17]) #mean for meter 10
}

############### Repeat this for all of the other EUs ########################
unique(diffAbunDat_tidy_FUNGI$EU)
######### EU_52 #########
ASV_23_EU_52_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[4]) %>%
  filter(EU == "EU_52") 
head(ASV_23_EU_52_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_23_meansAcrossEUs$EU[2] <- unique(ASV_23_EU_52_tidyTest_1$EU)
ASV_23_meansAcrossEUs$ASV_EUmean[2] <- mean(ASV_23_EU_52_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_23_meansAcrossEUs$ASV_EUsd[2] <- sd(ASV_23_EU_52_tidyTest_1$ASVabundance)
ASV_23_meansAcrossEUs$ASV_name[2] <- unique(ASV_23_EU_52_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_23_EU52 <- sort(unique(ASV_23_EU_52_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_23_EU52)){ #loop over number of unique meters (i.e. 10)
  ASV_23_meansAcrossEUs[2,which(colnames(ASV_23_meansAcrossEUs)==
                                  paste("mean_",uniqueMetersASV_23_EU52[i], sep=""))] <- mean(ASV_23_EU_52_tidyTest_1
                                                                                              [which(ASV_23_EU_52_tidyTest_1$Meter==uniqueMetersASV_23_EU52[i]), 17]) #mean for meter 10
}

###########  EU_53N ########## 
ASV_23_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[4]) %>%
  filter(EU == "EU_53N") 
head(ASV_23_EU_53N_tidyTest_1)
ASV_23_EU_53N_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[3]) %>%
  filter(EU == "EU_53N") 
head(ASV_23_EU_53N_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_23_meansAcrossEUs$EU[3] <- unique(ASV_23_EU_53N_tidyTest_1$EU)
ASV_23_meansAcrossEUs$ASV_EUmean[3] <- mean(ASV_23_EU_53N_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_23_meansAcrossEUs$ASV_EUsd[3] <- sd(ASV_23_EU_53N_tidyTest_1$ASVabundance)
ASV_23_meansAcrossEUs$ASV_name[3] <- unique(ASV_23_EU_53N_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_23_EU53N <- sort(unique(ASV_23_EU_53N_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_23_EU53N)){ #loop over number of unique meters (i.e. 10)
  ASV_23_meansAcrossEUs[3,which(colnames(ASV_23_meansAcrossEUs)==
                                  paste("mean_",uniqueMetersASV_23_EU53N[i], sep=""))] <- mean(ASV_23_EU_53N_tidyTest_1
                                                                                               [which(ASV_23_EU_53N_tidyTest_1$Meter==uniqueMetersASV_23_EU53N[i]), 17]) #mean for meter 10
}

###########  EU_54S ########## 
ASV_23_EU_54S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[4]) %>%
  filter(EU == "EU_54S") 
head(ASV_23_EU_54S_tidyTest_1)
# Get ASV mean and stdev across all points on the transect
ASV_23_meansAcrossEUs$EU[4] <- unique(ASV_23_EU_54S_tidyTest_1$EU)
ASV_23_meansAcrossEUs$ASV_EUmean[4] <- mean(ASV_23_EU_54S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_23_meansAcrossEUs$ASV_EUsd[4] <- sd(ASV_23_EU_54S_tidyTest_1$ASVabundance)
ASV_23_meansAcrossEUs$ASV_name[4] <- unique(ASV_23_EU_54S_tidyTest_1$ASV_name)

# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_23_EU54S <- sort(unique(ASV_23_EU_54S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_23_EU54S)){ #loop over number of unique meters (i.e. 10)
  ASV_23_meansAcrossEUs[4,which(colnames(ASV_23_meansAcrossEUs)==
                                  paste("mean_",uniqueMetersASV_23_EU54S[i], sep=""))] <- mean(ASV_23_EU_54S_tidyTest_1
                                                                                               [which(ASV_23_EU_54S_tidyTest_1$Meter==uniqueMetersASV_23_EU54S[i]), 17]) #mean for meter 10
}

##########  EU_8 ########## 
ASV_23_EU_8_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[4]) %>%
  filter(EU == "EU_8") 
head(ASV_23_EU_8_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_23_meansAcrossEUs$EU[5] <- unique(ASV_23_EU_8_tidyTest_1$EU)
ASV_23_meansAcrossEUs$ASV_EUmean[5] <- mean(ASV_23_EU_8_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_23_meansAcrossEUs$ASV_EUsd[5] <- sd(ASV_23_EU_8_tidyTest_1$ASVabundance)
ASV_23_meansAcrossEUs$ASV_name[5] <- unique(ASV_23_EU_8_tidyTest_1$ASV_name)

# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_23_EU8 <- sort(unique(ASV_23_EU_8_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_23_EU8)){ #loop over number of unique meters (i.e. 10)
  ASV_23_meansAcrossEUs[5,which(colnames(ASV_23_meansAcrossEUs)==
                                  paste("mean_",uniqueMetersASV_23_EU8[i], sep=""))] <- mean(ASV_23_EU_8_tidyTest_1
                                                                                             [which(ASV_23_EU_8_tidyTest_1$Meter==uniqueMetersASV_23_EU8[i]), 17]) #mean for meter 10
}

########## EU_53S ###########
ASV_23_EU_53S_tidyTest_1 <- diffAbunDat_tidy_FUNGI %>% 
  filter(ASV_name == ASVnamesDA_FUNGI_forest[4]) %>%
  filter(EU == "EU_53S") 
head(ASV_23_EU_53S_tidyTest_1)

# Get ASV mean and stdev across all points on the transect
ASV_23_meansAcrossEUs$EU[6] <- unique(ASV_23_EU_53S_tidyTest_1$EU)
ASV_23_meansAcrossEUs$ASV_EUmean[6] <- mean(ASV_23_EU_53S_tidyTest_1$ASVabundance) #get ASV mean ACROSS all samples in EU
ASV_23_meansAcrossEUs$ASV_EUsd[6] <- sd(ASV_23_EU_53S_tidyTest_1$ASVabundance)
ASV_23_meansAcrossEUs$ASV_name[6] <- unique(ASV_23_EU_53S_tidyTest_1$ASV_name)
# Now, need to get mean ASV abundance at each point along the transect (i.e. average of up to 4 transects)
uniqueMetersASV_23_EU53S <- sort(unique(ASV_23_EU_53S_tidyTest_1$Meter)) #get all the unique meters (this way so that in future cases, 
# can account for ASV not being present at a given meter). Sort so that it goes in numeric order
for (i in 1:length(uniqueMetersASV_23_EU53S)){ #loop over number of unique meters (i.e. 10)
  ASV_23_meansAcrossEUs[6,which(colnames(ASV_23_meansAcrossEUs)==
                                  paste("mean_",uniqueMetersASV_23_EU53S[i], sep=""))] <- mean(ASV_23_EU_53S_tidyTest_1
                                                                                               [which(ASV_23_EU_53S_tidyTest_1$Meter==uniqueMetersASV_23_EU53S[i]), 17]) #mean for meter 10
}

############## CHECKING FIRST FOR LOOP THAT CREATED ASV_23_meansAcrossEUs ###############
###### Making sure the for loop above works ######
index <- which(ASV_23_EU_10_tidyTest_1$Meter==uniqueMetersASV_23_EU10[1]) #these are the rows
ASV_23_EU_10_tidyTest_1[index,] #this shows that it does in fact pull out all the transects for meter 10
ASV_23_EU_10_tidyTest_1[index, 17] #gives abundances. 
mean(ASV_23_EU_10_tidyTest_1[index, 17]) #okay! this is what we want!
# So, putting it all together:
mean(ASV_23_EU_10_tidyTest_1[which(ASV_23_EU_10_tidyTest_1$Meter==uniqueMetersASV_23_EU10[1]), 17])

# Does indexing for column name work?
ASV_23_meansAcrossEUs[,2] # Meter==uniqueMetersASV_23_EU10[i] #we want the column to be whatever is the name of the meter
paste(uniqueMetersASV_23_EU10[1])
whichColumn <- paste("mean_",uniqueMetersASV_23_EU10[1], sep="")
ASV_23_meansAcrossEUs[1,]
colnames(ASV_23_meansAcrossEUs)
which(colnames(ASV_23_meansAcrossEUs)=="mean_50")
ASV_23_meansAcrossEUs[1,which(colnames(ASV_23_meansAcrossEUs)==paste("mean_",uniqueMetersASV_23_EU10[9], sep=""))]
ASV_23_meansAcrossEUs[1,9]

# Check a few values:
check1 <- ASV_23_EU_10_tidyTest_1 %>%
  filter(Meter == "10")
checkAvg <- (1+1+1)/3
checkAvg ==  ASV_23_meansAcrossEUs[1,1]
mean(ASV_23_EU_10_tidyTest_1$ASVabundance) == ASV_23_meansAcrossEUs[1,12]
######################################################################

######################################################################

############ Z-SCORES (STILL FOR ASV23 ACROSS ALL EUS) ###########

# THIS CREATES A SEPARATE DATAFRAME FOR THE Z-SCORE CALCULATIONS 
ZscoresEUs_ASV23 <- as.data.frame(matrix(nrow=6, ncol=13)) #has 6 rows for 6 EUs and columns for EU, all points on transect, ASV EU mean, and ASV EU sd
colnames(ZscoresEUs_ASV23) <- c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100", "EU", "ASV_EUmean", "ASV_EUsd")
ZscoresEUs_ASV23[,11:13] <- ASV_23_meansAcrossEUs[,11:13]

# Now do vectorized calculation to take all of the mean ASV values (at each meter in each EU), minus mean in that EU, all divided by stdev
ZscoresEUs_ASV23[,1:10] <- ((ASV_23_meansAcrossEUs[,1:10]- ASV_23_meansAcrossEUs[,12])/ASV_23_meansAcrossEUs[,13])
# Check a few
ZscoresEUs_ASV23[2,2] == ((ASV_23_meansAcrossEUs[2,2] - ASV_23_meansAcrossEUs[2,12])/ASV_23_meansAcrossEUs[2,13])
ZscoresEUs_ASV23[5,10] == ((ASV_23_meansAcrossEUs[5,10] - ASV_23_meansAcrossEUs[5,12])/ASV_23_meansAcrossEUs[5,13])

# Finally, need to make this in long form so that it works below:
ZscoresEUs_ASV23_longer <- ZscoresEUs_ASV23 %>% 
  pivot_longer(cols=`10`:`100`, names_to="Meter", values_to = "abundZ_score") 
ZscoresEUs_ASV23_longer$Meter <- as.numeric(ZscoresEUs_ASV23_longer$Meter)

# View(ZscoresEUs_ASV23_longer) #now there are 60 values, just as we wanted!
##########################################################################
# 3.  FITTING LOGISTIC CURVES
##########################################################################

# Trying first with just ASV 23-- doesn't work!
logFit_ASV23_f <- log.fitdiffAbundFunct(y= ZscoresEUs_ASV23_longer$abundZ_score, x= ZscoresEUs_ASV23_longer$Meter, ASVnames="ASV_23")
logFit_ASV23_f #Error in qr.solve(QR.B, cc) : singular matrix 'a' in solve

##### Trying a new following explanation here: https://rpubs.com/angelov/growthcurver ######
# First, because these are growth plots, add a 1 to all of the z-scores for correct fit
ZscoresEUs_ASV23_longerPlus1 <- ZscoresEUs_ASV23_longer
ZscoresEUs_ASV23_longerPlus1$abundZ_score <- ZscoresEUs_ASV23_longer$abundZ_score + 1
ggplot(ZscoresEUs_ASV23_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=0.7) +
  theme_bw()

modelASV23plus1.f <- growthcurver::SummarizeGrowth(data_t=ZscoresEUs_ASV23_longerPlus1$Meter, data_n=ZscoresEUs_ASV23_longerPlus1$abundZ_score, bg_correct = "none")
modelASV23plus1.f$vals #gives all of the values
predict(modelASV23plus1.f$model) # gives you the predicted abundance values (according to the model)
str(modelASV23plus1.f)
modelASV23plus1.f$vals$r #this is growth rate constant NOT r-squared
modelASV23plus1.f$vals$sigma #0.4644079 -- the smaller the better!
# sigma is a measure of the goodnesss of fit of the parameters of the logistic equation for the data; 
# it is the residual standard error from the nonlinear regression model. Smaller sigma values indicate
# a better fit of the logistic curve to the data than larger values.
plot(predict(modelASV23plus1.f$model) )
plot(modelASV23plus1.f$model$m$fitted())
modelASV23plus1.f$model$m$fitted()

### Plotting ##
# Base R
plot(modelASV23plus1.f) #ugly!

# ggplot
diffAbunDat_tidy_FUNGI[2,] #getting info for ASV_23 to add to ggtitle
modelASV23plus1.f$vals$t_mid #inflection point
modelASV23plus1.f$vals$sigma #0.4772041
p1 <- ggplot(ZscoresEUs_ASV23_longerPlus1, aes(x = Meter, y = abundZ_score)) + geom_point(alpha=1.3) + theme_bw()
p1
# Adding predicted values
ASV_23_EU_10_tidyTest_1 #look for name for below
df.predicted <- data.frame(Meter = ZscoresEUs_ASV23_longerPlus1$Meter, pred.Zabund = modelASV23plus1.f$model$m$fitted())
p2 <- p1 + geom_line(data=df.predicted, aes(y=pred.Zabund), color="red", size= 4) + ggtitle("Russula subsulphurea (Basidiomycota)") +
  #geom_point(x=modelASV23plus1.f$vals$t_mid, y = 0.78, color= "red", size=6) + #Here I just guessed a y based on how it looked!
  geom_vline(xintercept = 50, linetype= "dashed", color= "darkgrey", size=2) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) #make it so all meters show on x-axis!
quartz()
p2



df.predicted[c(1:10),] == df.predicted[c(11:20),] #this data frame is longer than it needs to be, but that's okay!
df.predicted[which(df.predicted$Meter==10),2]  # At 10 meters, abundance prediction (+ 1) is 0.4173783 on the y-axis
# Depth = 80% of distance between edge and abundance at end point, or 
modelASV23plus1.f$vals$t_mid - 10 


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
