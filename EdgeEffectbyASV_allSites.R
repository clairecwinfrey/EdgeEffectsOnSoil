# EdgeEdgebyASV_allSites.R
# May 4, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect ASVs using differential abundance analyses
# (DESeq2). This script:
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. postUbiquity.ps)).
# 2. Creates data frame with z-score (of ASV abundance) for each ASV within each EU 
# 3. Fits logistic curves to each one of the differentially abundant ASVs:
# -- i. gets Z-score (of ASV abundance) for each ASV within each EU 
# -- ii. for each ASV, then fits logistic models to THIS
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential abundance breakdown by forest, patch,
# and non-differentially abundant microbes 

# IdentifyEdgeEffectTaxa5.R investigated linear fit and logisti model fit of the
# the data, and found that logistic models were much better.

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("grid")
library("stringr") #for grep-like tools for data manipulation with character strings in data
library("drc") # might use for log fit thing

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <45
# R OBJECT MADE IN: UbiquityMedianSetup.R

######

# FUNCTIONS DEFINED IN THIS SCRIPT (but often used first in other scripts):
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVsampOccurrence determines the number of samples that each ASV occurs in 
# For proof that ASVsampOccurrence works, see UbiquityMedianSetup.R
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
    plot(y ~ x, main = paste("Logistic Function for", ASVnames[j]), xlab= "distance (m)", ylab= "ASV abundance")
    lines(0:max(x), predict(log.ss, data.frame(x=0:max(x))), col="red")
  }
  r1 <- sum((x - mean(x))^2)
  r2 <- sum(residuals(log.ss)^2)
  
  r_sq <- (r1 - r2) / r1 #check, is this right?
  
  out <- data.frame(cbind(c(C=C, a=A, k=K, R.value=sqrt(r_sq))))
  names(out)[1] <- "Logistic Curve"
  
  return(out)
}

# 7. zScore function computes the z-score for a given vector of numbers
# The function is broadly applicable, but ASVs should be rows for its usage in this script  
zScore <- function(dat) { # input is dataframe
  zScoreDf <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
  colnames(zScoreDf) <- colnames(dat)
  rownames(zScoreDf) <- rownames(dat)
  for (i in 1:nrow(dat)){
    ASVmean <- rowMeans(dat[i,]) #if you change rowMeans to mean, could use with matrices. Or could build in if/else statement
    ASVsd <- sd(dat[i,])
    for (j in 1:length(dat[i,])) {
      zScoreDf[i,j] <- (dat[i,j]-ASVmean)/ASVsd # subtract row mean from value, then divide by standard deviation for z score
    }
  }
  return(zScoreDf)
}

# Proof/testing out zScore function to make sure that it works
set.seed(19)
dat <- rpois(n=100, lambda = 1) 
# in cartoon example, as with real data above, ASVs are rows and samples are columns
datmat <- matrix(dat, nrow=20, ncol=5) 
datmat <- as.data.frame(datmat)
# testing it out below, it seems to work
test <- zScore(datmat)
dim(test) == dim(datmat) #yes
test[1,1] == (datmat[1,1] - rowMeans(datmat[1,]))/sd(datmat[1,])
(datmat[3,4] - rowMeans(datmat[3,]))/sd(datmat[3,]) == test[3,4]
(datmat[20,2] - rowMeans(datmat[20,]))/sd(datmat[20,]) == test[20,2]
############################################################


############################################################
# FUNCTION FOR GETTING LOGISTICAL FIT OF EACH INDICATOR ASV
############################################################
#### REPLACE POSTUBIQUIRY THING WITH PHYSEQ TO RESTORE FUNCTION########
# This function
diffAbundfunct2 <- function(physeq, alpha){ 
  #physeq= phyloseq object (single EU represented),alpha is alpha significance level for DESeq analysis
  NoEdge.ps <- subset_samples(EU_52_Soils.ps, Habitat != "edge") #remove edge samples
  step1 <- phyloseq::phyloseq_to_deseq2(NoEdge.ps, ~ Habitat) #set up DESeq2 dds object 
  step2 <- DESeq2::DESeq(step1, test="Wald", fitType = "parametric") #differential expression analysis step;
  #uses default Benjamini-Hochberg correction 
  DeseqResults <- results(step2, cooksCutoff = FALSE) #make results object
  DeseqResults <- DeseqResults[which(DeseqResults$padj < alpha), ] #only get those ASVs below alpha level
  DeseqResults <- cbind(as(DeseqResults, "data.frame"), as(tax_table(NoEdge.ps)[rownames(DeseqResults), ], "matrix")) #clean up format
  DeseqResults$Habitat <- ifelse(DeseqResults$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
  forestDA_ASVs <- DeseqResults[which(DeseqResults$Habitat=="forest"),] #just forest ASVs
  patchDA_ASVs <- DeseqResults[which(DeseqResults$Habitat=="patch"),] #just patch ASVs
  # Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
  sampDat <- sample_data(EU_52_Soils.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
  attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
  sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
  # Next 8 lines get abundances for each ASV across the transect for FOREST
  forestNames <- rownames(forestDA_ASVs) #pull out ASV names
  ASVsAll <- as.data.frame(t(ASVs_outta_ps(EU_52_Soils.ps))) #get ASVs from original and transpose so that ASVs are rows
  forestASVtab <- ASVsAll[forestNames,] #get a smaller version of the ASV table that has only the diff abundance ASVs
  # Create dataframe with everything of interest
  forest_diffAbunDat <- merge(forestDA_ASVs, forestASVtab, by= "row.names") #FOREST samples merging based on shared ASV (row) names
  forest_diffAbunDat_tidy <- forest_diffAbunDat %>% 
    pivot_longer(cols= 16:ncol(forest_diffAbunDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
    merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
  # REPEAT THE ABOVE FOR PATCH
  patchNames <- rownames(patchDA_ASVs) #pull out ASV names
  patchASVtab <- ASVsAll[patchNames,] #get a smaller version of the ASV table that has only the diff abundance ASVs
  patch_diffAbundDat <- merge(patchDA_ASVs, patchASVtab, by= "row.names") #PATCH samples merging based on shared ASV (row) names
  patch_diffAbunDat_tidy <- patch_diffAbundDat %>% 
    pivot_longer(cols= 16:ncol(patch_diffAbundDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of patch ASVs x sample number
    merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
  # Next two lines print out how many differentially abundant ASVs were found in the patch and the forest
  print(paste("analysis found", nrow(forestDA_ASVs), "differentially abundant ASVs in forest")) 
  print(paste("analysis found", nrow(patchDA_ASVs), "differentially abundant ASVs in patch"))
  diffAbund_df <- vector("list", 2)
  diffAbund_df[[1]] <- forest_diffAbunDat_tidy
  diffAbund_df[[2]] <- patch_diffAbunDat_tidy
  names(diffAbund_df) <- c("forestASVs", "patchASVs")
  # The rest of the function from here on gets the logistic fit of each curve
  # First for forest!
  forestASVmeterAbunds <- vector("list", length(forestNames)) #pre-allocate space for each ASV's abundance for each ASV
  forestlogFitList <- vector("list", length(forestNames)) #pre-allocate space for each ASV's logistic fit info
  names(forestASVmeterAbunds) <- forestNames
  names(forestlogFitList) <- paste(forestNames, "_forest")
  for (i in 1:length(forestASVmeterAbunds)){
    tryCatch({
      forestASVmeterAbunds[[i]] <- forest_diffAbunDat_tidy[which(forest_diffAbunDat_tidy$Row.names==forestNames[i]),c(10,13, 17:21)] #keeps phylum, family, ASV abundance, sample ID, EU, transect, and meter
      forestASVmeterAbunds[[i]] <- forestASVmeterAbunds[[i]] %>% 
        dplyr::arrange(Meter) #arrange by meter, in descending order
      forestlogFitList[[i]] <- log.fitdiffAbundFunct(x=forestASVmeterAbunds[[i]]$Meter, y=forestASVmeterAbunds[[i]]$ASVabundance, ASVnames=forestNames) #errors are ignored so that those that fit can be fit
    }, error=function(e){})
  }
  # Now for patch!
  patchASVmeterAbunds <- vector("list", length(patchNames)) #pre-allocate space for each ASV's abundance for each ASV
  patchlogFitList <- vector("list", length(patchNames)) #pre-allocate space for each ASV's logistic fit info
  names(patchASVmeterAbunds) <- patchNames
  names(patchlogFitList) <- paste(patchNames, "_patch")
  for (i in 1:length(patchASVmeterAbunds)){
    tryCatch({
      patchASVmeterAbunds[[i]] <- patch_diffAbunDat_tidy[which(patch_diffAbunDat_tidy$Row.names==patchNames[i]),c(10,13, 17:21)] #keeps phylum, family, ASV abundance, sample ID, EU, transect, and meter
      patchASVmeterAbunds[[i]] <- patchASVmeterAbunds[[i]] %>% 
        dplyr::arrange(Meter) #arrange by meter, in descending order
      patchlogFitList[[i]] <- log.fitdiffAbundFunct(x=patchASVmeterAbunds[[i]]$Meter, y=patchASVmeterAbunds[[i]]$ASVabundance, ASVnames=patchNames) #errors are ignored so that those that fit can be fit
    }, error=function(e){})
  }
  #what percentage did not get a logistic fit for one reason or another?
  forestSuccess <- round((1-sum(sapply(forestlogFitList, is.null))/length(forestlogFitList))*100)
  print(paste(forestSuccess,"% of forest ASVs were fit to a logistic model"))
  patchSuccess <- round((1-sum(sapply(patchlogFitList, is.null))/length(patchlogFitList))*100)
  print(paste(patchSuccess,"% of patch ASVs were fit to a logistic model"))
  logFitsList <- vector("list", 2) #make space for the results for forest and patch
  logFitsList <- c(forestlogFitList, patchlogFitList)
  # Now to clean up the results to make it in a tidy dataframe
  logFits.df <- as.data.frame(matrix(nrow = length(logFitsList), ncol = 5)) #pre-allocate matrix to hold everything in. 
  # rows will correspond to diff abund ASVs and cols variables and output
  colnames(logFits.df) <- c("forestOrPatch", "asymptote", "inflectionPoint", "k", "R-value")
  rownames(logFits.df) <- names(logFitsList)
  logFits.df[,1] <- c(rep("forest", length(forestlogFitList)), rep("patch", length(patchlogFitList)))
  for (i in 1:length(logFitsList)) {
    tryCatch({ #put in tryCatch so that it would still fill in values for the ones that DID work
      logFits.df[i,2] <- logFitsList[[i]][1,]
      logFits.df[i,3] <- logFitsList[[i]][2,]
      logFits.df[i,4] <- logFitsList[[i]][3,]
      logFits.df[i,5] <- logFitsList[[i]][4,]
    }, error=function(e){})
  }
  return(logFits.df)
}

try1 <- diffAbundfunct2(postUbiquity.ps, alpha=0.001)
#1695
541 + 1155
View(try1)

##########################################################################
# INDICATOR ANALYSIS AND TIDY DATAFRAME
##########################################################################

# Perform indicator species analysis on the post ubiquity dataset, and then
# make a dataframe which has a row corresponging to the abundance of each differentially
# abundant ASV in each place along the transect, as well as taxonomic info, and 
# whether or not that ASV was differentially abundant in patch or the forest.
  NoEdge.ps <- subset_samples(postUbiquity.ps, Habitat != "edge") #remove edge samples
  step1 <- phyloseq::phyloseq_to_deseq2(NoEdge.ps, ~ Habitat) #set up DESeq2 dds object 
  step2 <- DESeq2::DESeq(step1, test="Wald", fitType = "parametric") #differential expression analysis step;
  #uses default Benjamini-Hochberg correction 
  DeseqResults <- results(step2, cooksCutoff = FALSE) #make results object
  DeseqResults <- DeseqResults[which(DeseqResults$padj < 0.001), ] #only get those ASVs below alpha level
  DeseqResults <- cbind(as(DeseqResults, "data.frame"), as(tax_table(NoEdge.ps)[rownames(DeseqResults), ], "matrix")) #clean up format
  DeseqResults$Habitat <- ifelse(DeseqResults$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
  # Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
  sampDat <- sample_data(postUbiquity.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
  attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
  sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
  ### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
  ASVnamesDA <- rownames(DeseqResults) #1696 found
  ASVsAll <- as.data.frame(t(ASVs_outta_ps(postUbiquity.ps))) #get ASVs from original and transpose so that ASVs are rows (4,480 ASVs in OG)
  # Create dataframe with everything of interest
  diffAbunDat <- merge(DeseqResults, ASVsAll, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
  diffAbunDat_tidy <- diffAbunDat %>% 
    pivot_longer(cols= 16:ncol(diffAbunDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
    merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_tidy)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
 #View(diffAbunDat_tidy) this has 395,168 rows, which is equal to 233 (number of samples) x 1696 (number of diff abund ASVs)
 ##########

# THIS BELOW DOESN'T WORK SO I'M GREYING OUT FOR NOW
#diffAbunDat_wide <- diffAbunDat %>% 
#  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
#colnames(diffAbunDat_wide)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
  

##########################################################################
# GET Z-SCORES OF EACH ASV WITHIN EU
##########################################################################
# 1. Samples based on raw (i.e. not median) abundance
# Prune samples to separate by EU and check to make sure each looks right
EU_52_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU) #only EU 52
EU_53N_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_53N")
unique(sample_data(EU_53N_Soils.ps)$EU)
EU_54S_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_54S")
unique(sample_data(EU_54S_Soils.ps)$EU)
EU_8_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_8")
unique(sample_data(EU_8_Soils.ps)$EU)
EU_53S_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_53S")
unique(sample_data(EU_53S_Soils.ps)$EU)
EU_10_Soils.ps <- subset_samples(postUbiquity.ps, EU == "EU_10")
unique(sample_data(EU_10_Soils.ps)$EU)

# 2. Make list of all of these ASV tables by EU:
ASVtabByEU_List <- vector("list", 6) #make a list to hold the ASV tabs from each 
# Fill each element in the list in by each EU (all currently have 4480 taxa, i.e. the whole dataset)
ASVtabByEU_List[[1]] <- EU_52_Soils.ps #4480 taxa
ASVtabByEU_List[[2]] <- EU_53N_Soils.ps
ASVtabByEU_List[[3]] <- EU_54S_Soils.ps
ASVtabByEU_List[[4]] <- EU_8_Soils.ps
ASVtabByEU_List[[5]] <- EU_53S_Soils.ps
ASVtabByEU_List[[6]] <- EU_10_Soils.ps

# 3. Remove ASVs that are not differentially abundant
for (j in 1:length(ASVtabByEU_List)) {
  otu_table(ASVtabByEU_List[[j]]) <- otu_table(ASVtabByEU_List[[j]])[ASVnamesDA] #pull out only diff abundant ASVs
}

# 4. Remove ASVs that, for any given EU, are not present 
# pre-allocate zero-index list thing
zeroIndexList <- vector("list", 6) #pre-allocate a new list of 6 vectors, one per EU
for (k in 1:length(ASVtabByEU_List)) {
  # next line finds each ASV in each EU, that is not present at least once and makes an index for it
  zeroIndexList[[k]] <- which(rowSums(otu_table(ASVtabByEU_List[[k]]))==0) 
  # line below takes out all of these ASVs that do not appear in each EU.
  otu_table(ASVtabByEU_List[[k]]) <- otu_table(ASVtabByEU_List[[k]])[-zeroIndexList[[k]]] 
}
# Check a few to see if code above worked
which(rowSums(otu_table(ASVtabByEU_List[[1]]))==0) 
which(rowSums(otu_table(ASVtabByEU_List[[6]]))==0)

################
#### Scratch work for (4) this above #####
#which(rowSums(otu_table(ASVtabByEU_List[[1]]))==0)
#otu_table(ASVtabByEU_List[[1]])
#zero_index_52 <- which(rowSums(otu_table(EU_52_Soils.ps))==0)
#names(zero_index_52)
#dim(otu_table(EU_52_Soils.ps)[zero_index_52]) # 13 ASVs across 38 samples
#otu_table(EU_52_Soils.ps) <- otu_table(EU_52_Soils.ps)[-zero_index_52] # 4467 ASVs
#dim(otu_table(EU_52_Soils.ps))
#otu_table(EU_52_Soils.ps)[[1]]
#otu_table(ASVtabByEU_List[[1]]) <- otu_table(ASVtabByEU_List[[1]])[-zero_index_52] 
################

# 5. Rename elements of list by EU
names(ASVtabByEU_List)[[1]] <- "EU_52_z"
names(ASVtabByEU_List)[[2]] <- "EU_53N_z"
names(ASVtabByEU_List)[[3]] <- "EU_54S_z"
names(ASVtabByEU_List)[[4]] <- "EU_8_z"
names(ASVtabByEU_List)[[5]]<- "EU_53S_z"
names(ASVtabByEU_List)[[6]] <- "EU_10_z"

ASVtabByEU_List[[1]] #shows that these are still phyloseq objects

# 6. Get each element in the list out of phyloseq and get z-scores
for (k in 1:length(ASVtabByEU_List)){ #ASVtabByEU_List is a list of phyloseq objects that has only the differentially abundant ASVs
  # in each EU
  ASVtabByEU_List[[k]] <- t(ASVs_outta_ps(ASVtabByEU_List[[k]])) #get ASV table out of phyloseq and invert
 
  ASVtabByEU_List[[k]] <- zScore(as.data.frame(ASVtabByEU_List[[k]])) #get z-score for each ASV (where mean and stdev are over all in that EU)
}

################
# This area of "scratch work" is to see if it matters to replace 0s with 0.01 in the z-score step above
# In a nutshell, only 118/64896 observations are zero, or 0.18%, so I think not.
#ASVzs0 <- ASVtabByEU_List[[2]] #ad
#ASVzs0.01 <- ASVtabByEU_List[[2]] #made with replacing all 0s with 0.01 with this line in the for loop:
# # ASVtabByEU_List[[k]][ASVtabByEU_List[[k]] == 0] <- 0.01
#ASVzs0 <- ASVtabByEU_List[[2]]
#length(ASVtabByEU_List[[2]])
#ASVzs0.01 <- ASVtabByEU_List[[2]]
#length(which(ASVzs0==0))
#length(which(ASVzs0.01==0))
#ASVzs0.01[(which(ASVzs0==0))]
#ASVzs0[(which(ASVzs0==0))]
#dim(ASVzs0.01)
#dim(ASVzs0)
#all.equal(ASVzs0, ASVzs0.01)
################


# Does it change much if we do NOT replace the 0 with 0.01? (before I did this, I re-ran the code up to the previous for loop)
ASVtabByEU_List2 <- vector("list", 6) #pre-allocate a new list of 6 vectors, one per EU
for (k in 1:length(ASVtabByEU_List2)){ #ASVtabByEU_List is a list of phyloseq objects that has only the differentially abundant ASVs
  # in each EU
  ASVtabByEU_List2[[k]] <- t(ASVs_outta_ps(ASVtabByEU_List[[k]])) #get ASV table out of phyloseq and invert
  ASVtabByEU_List2[[k]] <- zScore(as.data.frame(ASVtabByEU_List[[k]])) #get z-score for each ASV (where mean and stdev are over all in that EU)
}

head(ASVtabByEU_List2[[2]])

####### scratch work for this above:
### if you use this scratch work, ust re-create ASVtabByEU_List objects
ASVtabByEU_List[[1]] <- ASVs_outta_ps(ASVtabByEU_List[[1]])
class(ASVtabByEU_List[[1]])
tiny <- ASVtabByEU_List[[1]][1:12, 1:12] #try on a smaller subset
tiny[tiny == 0] <- 0.01 #it works!
tiny2 <- zScore(as.data.frame(tiny))

ASVtabByEU_List[[1]][ASVtabByEU_List[[1]] == 0] <- 0.01
############################

# 7. Now, merge all these elements of the list by rowname, so that the there are as many rows as ASVs
# and as many columns as samples (i.e. 233)
## will have to periodically reset rownames to col one, so that can merge by them
merge_1 <- merge(ASVtabByEU_List[[1]], ASVtabByEU_List[[2]], by= "row.names", all=TRUE) #merge first two EUs, all = TRUE merges even those ASVs that don't occur in both
dim(ASVtabByEU_List[[1]])
dim(ASVtabByEU_List[[2]])
dim(merge_1) #1694, 78
#View(merge_1) # now there are some NAs, but that should be okay for next step, I think
merge_1 <- merge_1 %>% 
  column_to_rownames(var="Row.names")
head(merge_1)
class(merge_1)
merge_2 <- merge(merge_1, ASVtabByEU_List[[3]], by= "row.names",  all=TRUE)#merge first 2 EUs with EU 3
dim(merge_2) #1696, 116
merge_2 <- merge_2 %>% 
  column_to_rownames(var="Row.names")
head(merge_2)
merge_3 <- merge(merge_2, ASVtabByEU_List[[4]], by= "row.names", all=TRUE) #merge first 3 EUs with EU 4
dim(merge_3) #1696, 155
merge_3 <- merge_3 %>% 
  column_to_rownames(var="Row.names")
head(merge_3)
# View(merge_3)
merge_4 <- merge(merge_3, ASVtabByEU_List[[5]], by= "row.names", all=TRUE) #merge first 4 EUs with EU 5
dim(merge_4) #1696, 195
merge_4 <- merge_4 %>% 
  column_to_rownames(var="Row.names")
head(merge_4)
merge_5 <- merge(merge_4, ASVtabByEU_List[[6]], by= "row.names", all=TRUE)#merge first 5 EUs with EU 6
dim(merge_5) #1696, 234
merge_5 <- merge_5 %>% 
  column_to_rownames(var="Row.names")
head(merge_5)
abundZscores_allEUs <- merge_5
#View(abundZscores_allEUs)
sum(apply(abundZscores_allEUs,2,is.nan)) #no NaNs 
length(which(is.na(abundZscores_allEUs)))  #1980 (what does this mean?). Are there that many NAs?! 
dim(abundZscores_allEUs) #1696  233
1696*233 #395168 total values

abundZscores_allEUs[385,]

##########################################################################
# LOGISTIC FIT FOR EACH ASV
##########################################################################

# 1. Take abundZscores_allEUs and get relevant metadata so that can do logistic regression thing
head(DeseqResults) #made earlier in this script
sampDat # also made earlier in this script, has sample.ID, EU, Transect, Meter, SampleNumberID
### this repeats a lot of the script from earlier
### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ASVnamesDA <- rownames(DeseqResults) #1696 found
#ASVsAll <- as.data.frame(t(ASVs_outta_ps(postUbiquity.ps))) #get ASVs from original and transpose so that ASVs are rows (4,480 ASVs in OG)
# Create dataframe with everything of interest
diffAbun_ZDat <- merge(DeseqResults, abundZscores_allEUs, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
#View(diffAbun_ZDat)
diffAbunDat_Z_tidy <- diffAbun_ZDat %>% 
  pivot_longer(cols= 16:ncol(diffAbun_ZDat), names_to= "SampleNumberID", values_to= "ZabundASV") %>% #has # of rows equal to # of diff abund ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_Z_tidy)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".

# 2. Get logistic fits for all of these ASVs!
ASVmeterAbunds <- vector("list", length(ASVnamesDA)) #pre-allocate space for each ASV's abundance for each ASV
logFitList <- vector("list", length(ASVnamesDA)) #pre-allocate space for each ASV's logistic fit info
names(ASVmeterAbunds) <- ASVnamesDA
names(logFitList) <- ASVnamesDA
for (i in 1:length(ASVmeterAbunds)){
  tryCatch({
    ASVmeterAbunds[[i]] <- diffAbunDat_Z_tidy[which(diffAbunDat_Z_tidy$ASV_name==ASVnamesDA[i]),c(10,13, 17:21)] #keeps phylum, family, ASV abundance, sample ID, EU, transect, and meter
    ASVmeterAbunds[[i]] <- ASVmeterAbunds[[i]] %>% 
      dplyr::arrange(Meter) #arrange by meter, in descending order
  }, error=function(e){})
}    

# greying out for now for testing purposes
#    logFitList[[i]] <- log.fitdiffAbundFunct(x=ASVmeterAbunds[[i]]$Meter, y=ASVmeterAbunds[[i]]$ASVabundance, ASVnames=ASVnamesDA) #errors are ignored so that those that fit can be fit
# }, error=function(e){})
# }

# does this work?
log.fitdiffAbundFunct(y= ASVmeterAbunds[[1]]$ZabundASV, x= ASVmeterAbunds[[1]]$Meter, ASVnames=names(ASVmeterAbunds)[[1]]) 
  #Error in qr.solve(QR.B, cc) : singular matrix 'a' in solve
# don't know what this issue is, but this is it!!

log.fitdiffAbundFunct(x=ASVmeterAbunds[[1]]$Meter, y=ASVmeterAbunds[[i]]$ASVabundance, ASVnames=ASVnamesDA)
diffAbunDat_Z_tidy[which(diffAbunDat_Z_tidy$ASV_name==ASVnamesDA[1]),c(10,13, 17:21)]

#### BUT!!!! WE MAY BE ABLE TO PLOT IT ANYWAY!!
# THIS DIDNT WORK!!!
plotVec <- rep(NA, length(ASVmeterAbunds)) #pre-allocate
for (i in 1:length(length(ASVmeterAbunds))){
  plotVec[i] <- plot(ASVmeterAbunds[[i]]$ZabundASV ~ ASVmeterAbunds[[i]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[i]]), xlab= "distance (m)", ylab= "ASV abundance")
}

for (i in 1:length(length(ASVmeterAbunds))){
  plot(ASVmeterAbunds[[i]]$ZabundASV ~ ASVmeterAbunds[[i]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[i]]), xlab= "distance (m)", ylab= "ASV abundance")
}


for (j in 1:length(length(ASVmeterAbunds))){
  plotVec[i] <- plot(ASVmeterAbunds[[i]]$ZabundASV ~ ASVmeterAbunds[[i]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[i]]), xlab= "distance (m)", ylab= "ASV abundance")
  print(plot(plotVec[i]))
}

quartz()
plot(ASVmeterAbunds[[1]]$ZabundASV ~ ASVmeterAbunds[[1]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[1]]), xlab= "distance (m)", ylab= "ASV abundance")

quartz()
plot(ASVmeterAbunds[[2]]$ZabundASV ~ ASVmeterAbunds[[2]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[2]]), xlab= "distance (m)", ylab= "ASV abundance")

quartz()
plot(ASVmeterAbunds[[100]]$ZabundASV ~ ASVmeterAbunds[[100]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[100]]), xlab= "distance (m)", ylab= "ASV abundance")

quartz()
plot(ASVmeterAbunds[[243]]$ZabundASV ~ ASVmeterAbunds[[243]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[243]]), xlab= "distance (m)", ylab= "ASV abundance")

quartz()
plot(ASVmeterAbunds[[1000]]$ZabundASV ~ ASVmeterAbunds[[1000]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[1000]]), xlab= "distance (m)", ylab= "ASV abundance")

quartz()
plot(ASVmeterAbunds[[681]]$ZabundASV ~ ASVmeterAbunds[[681]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[681]]), xlab= "distance (m)", ylab= "ASV abundance")


which(names(ASVmeterAbunds)=="ASV_1280")

my_plot <- recordPlot()
plot.new()  
my_plot

length(ASVmeterAbunds)
##########################################################################
# STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################

DeseqResultsMini <- DeseqResults[,c(8,14)] %>%  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
postUbiqTaxTab <- taxtable_outta_ps(postUbiquity.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == TRUE)) #1696 
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE)) #2784 #which ASVs are NOT differentially abundant?
false_index <- which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE) #also 2,784
notDA_taxTab <- postUbiqTaxTab[false_index,] #get a taxonomy tab with ONLY the non-differentially abundant ASVs
notDA_taxTab$Habitat <- "AremainingASVs" #make a habitat column that labels these as NOT differentially abundant. A in front so that would be first
# in ggplot for ease.
# View(notDA_taxTab)
colnames(notDA_taxTab)
notDA_taxTabMini <- notDA_taxTab[,c(2,8)] #keep only phylum and habitat to match DeseqResultsMini
DAphylumAll <- rbind(DeseqResultsMini, notDA_taxTabMini) #this has ASV name, phylum, and whether diff abundant for ALL ASVs in postUbiquity analysis
# What are the phyla breakdown here?
# so effectively, we want to get numbers in each phyla in each group. Should just be able to plot this?
quartz()
ggplot(DAphylumAll, aes(fill=Habitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant ASVs")

# A few checks to make sure that the counting above is working as expected
# Acidobacteria
length(which(DAphylumAll$Phylum=="Acidobacteria")) #930
acido_index <- which(DAphylumAll$Phylum=="Acidobacteria")
length(which(DAphylumAll[acido_index,]$Habitat=="forest")) #178 forest specialists within Acidobacteria
length(which(DAphylumAll[acido_index,]$Habitat=="patch")) #207 patch specialists within Acidobacteria
length(which(DAphylumAll[acido_index,]$Habitat=="AremainingASVs")) #545 remaining ASVs
178+ 207 + 545 # =930

# Firmicutes
length(which(DAphylumAll$Phylum=="Firmicutes")) #18
firmi_index <- which(DAphylumAll$Phylum=="Firmicutes")
length(which(DAphylumAll[firmi_index,]$Habitat=="forest")) #0 forest specialists 
length(which(DAphylumAll[firmi_index,]$Habitat=="patch")) #2 patch specialists
length(which(DAphylumAll[firmi_index,]$Habitat=="AremainingASVs")) #16 remaining ASVs

#Chloroflexi
length(which(DAphylumAll$Phylum=="Chloroflexi")) #18
chloro_index <- which(DAphylumAll$Phylum=="Chloroflexi")
length(which(DAphylumAll[chloro_index,]$Habitat=="forest")) #4 forest specialists 
length(which(DAphylumAll[chloro_index,]$Habitat=="patch")) #238 patch specialists
length(which(DAphylumAll[chloro_index,]$Habitat=="AremainingASVs")) #180 remaining ASVs


################################################################################################
# DISCARD ALL OF THIS (PROBABLY)
################################################################################################
#####################
# The rest of this is to be discarded (likely): 

# REMOVE ASVS THAT DO NOT OCCUR IN EU
#########
# Proof to show that this above works:
zero_index_52 <- which(rowSums(otu_table(EU_52_Soils.ps))==0)
names(zero_index_52)
dim(otu_table(EU_52_Soils.ps)[zero_index_52]) # 13 ASVs across 38 samples
otu_table(EU_52_Soils.ps) <- otu_table(EU_52_Soils.ps)[-zero_index_52] # 4467 ASVs
dim(otu_table(EU_52_Soils.ps))
#########


colnames(diffAbunDat_tidy)
diffAbunDat_tidy_wider <- diffAbunDat_tidy %>% 
  group_by(EU, ASV_name) %>% 
  pivot_wider(names_from= Transect, values_from = ASVabundance)

# Write a for loop that makes a new column for each ASV at a certain EU_transect_meter
# This way, there is a unique identifier for each ASV within each sample (beyond just combination
# of info in columns)
for (i in 1:nrow(diffAbunDat_tidy)){
  diffAbunDat_tidy$EU_TransectMeter[i] <- paste(diffAbunDat_tidy$ASV_name[i], "_", diffAbunDat_tidy$EU[i], diffAbunDat_tidy$Transect[i], "_", diffAbunDat_tidy$Meter[i], sep="")
}
head(diffAbunDat_tidy)

small <- diffAbunDat_tidy[1:100, c(2, 10, 17, 19:ncol(diffAbunDat_tidy))]
View(small)
for (i in 1:nrow(small)){
  small$ASV_nameEU[i] <- paste(small$ASV_name[i], "_", small$EU[i], sep="")
}
ncol(small)
small[,9] <- small[,1]
colnames(small)[9] <- "ASV_name"
View(small)
small[,1] <- small[,8]
colnames(small)[1] <- "ASV_nameEU"
View(small)
small[,8] <- NULL
View(small)

small2 <- small %>% 
  pivot_wider(names_from= EU_TransectMeter, values_from = ASVabundance) 
View(small2)

small3 <- small %>% 
  group_by(ASV_nameEU) %>% 
  pivot_wider(names_from= EU_TransectMeter, values_from = ASVabundance) 
View(small3)

small$ASVabundance


small_try <- small %>% 
  group_by(EU, ASV_name) %>%
  pivot_wider(names_from= EU_TransectMeter, values_from = ASVabundance)
View(small_try)

# Now make it wider so that 
diffAbunDat_tidy_wider <- diffAbunDat_tidy %>% 
  pivot_wider(names_from= EU_TransectMeter, values_from = ASVabundance)

diffAbunDat_tidy3 <- diffAbunDat_tidy2 %>% 
  pivot_wider(names_from= EU_TransectMeter, values_from = ASVabundance)

diffAbunDat_tidy_wider <- diffAbunDat_tidy %>% 
  group_by(EU, ASV_name) %>% 
  pivot_wider(names_from=EU_TransectMeter, values_from = ASVabundance)
head(diffAbunDat_tidy_wider)


View(diffAbunDat_tidy_wider)

#create data frame
df <- data.frame(player=rep(c('A', 'B'), each=4),
                 year=rep(c(1, 1, 2, 2), times=2),
                 stat=rep(c('points', 'assists'), times=4),
                 amount=c(14, 6, 18, 7, 22, 9, 38, 4))

#view data frame
df

df2 <- df %>% 
  pivot_wider(names_from=stat, values_from=amount)
