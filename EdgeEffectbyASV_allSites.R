# EdgeEdgebyASV_allSites.R
# May 4, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect ASVs using differential abundance analyses
# (DESeq2). This script:
# 1. Performs differential abundance analyses on all EUs together (working on the 
# dataset after removing v rare taxa and applying a ubiquity threshold 
# (i.e. postUbiquity.ps)).
# 2. Make a stacked barchart that shows number of differentially abundant ASVs in each  differential abundance breakdown by forest, patch,
# and non-differentially abundant microbes 
# 3. Fits logistic curves to each one of the differentially abundant ASVs:
# -- i. gets Z-score (of ASV abundance) for each ASV within each EU 
# -- ii. for each ASV, then fits logistic models through THIS

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

log.fitdiffAbundFunct <- function(y, x, ASVnames){ 
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

# This function
diffAbundfunct2 <- function(physeq, alpha){ 
  #physeq= phyloseq object (single EU represented),alpha is alpha significance level for DESeq analysis
  NoEdge.ps <- subset_samples(physeq, Habitat != "edge") #remove edge samples
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
  sampDat <- sample_data(physeq)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
  attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
  sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
  # Next 8 lines get abundances for each ASV across the transect for FOREST
  forestNames <- rownames(forestDA_ASVs) #pull out ASV names
  ASVsAll <- as.data.frame(t(ASVs_outta_ps(physeq))) #get ASVs from original and transpose so that ASVs are rows
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

##########################################################################
# GET Z-SCORES OF EACH ASV WITHIN EU
##########################################################################
colnames(diffAbunDat_tidy)
diffAbunDat_tidy_wider <- diffAbunDat_tidy %>% 
  group_by(EU, ASV_name) %>% 
  pivot_wider(names_from= Transect, values_from = ASVabundance)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
