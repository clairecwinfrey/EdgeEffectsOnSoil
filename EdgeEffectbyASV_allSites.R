# EdgeEdgebyASV_allSites.R -- 16S
# May 4, 2022 (begun)

# ~~~~CLEAN UP DESCRIPTION~~~~
# This script identifies edge effect BACTERIAL AND ARCHAEAL (16S) ASVs using differential abundance analyses
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

# IdentifyEdgeEffectTaxa5.R investigated linear fit and logistic model fit of the
# the data, and found that logistic models were much better.

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
library("drc") # might use for log fit thing

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
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

# The next steps add 1 to all of the ASV abundance counts and re-make a phyloseq object.
# This was necessary because DESeq2 could not compute the geometric mean of the samples
# since every gene had at least one zero in it (and was thus thrown out of the analysis)
# (see recommendations and explanations here:
# https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564)
postUbiqASVs_16S <- ASVs_outta_ps(postUbiquity.ps)
postUbiqASVs_16SPlus1 <- t(postUbiqASVs_16S + 1) #adding one to every abundance count for differential abundance analysis;
# see explanation a few lines down. Invert so that I can make into a new phyloseq object.
# Make a new phyloseq object with these above
OTU = otu_table(postUbiqASVs_16SPlus1, taxa_are_rows = TRUE)
postUbiqASVs16S_Plus1.ps <- phyloseq(OTU, tax_table(postUbiquity.ps), sample_data(postUbiquity.ps))

# Remove edge samples to compare patch and matrix
NoEdgePlus1.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, Habitat != "edge") 
NoEdgePlus1.ps
  step1 <- phyloseq::phyloseq_to_deseq2(NoEdgePlus1.ps, ~ Habitat) #set up DESeq2 dds object 
  step2 <- DESeq2::DESeq(step1, test="Wald", fitType = "parametric") #differential expression analysis step;
  #uses default Benjamini-Hochberg correction 
  DeseqResults <- results(step2, cooksCutoff = FALSE) #make results object
  DeseqResults <- DeseqResults[which(DeseqResults$padj < 0.001), ] #only get those ASVs below alpha level
  DeseqResults <- cbind(as(DeseqResults, "data.frame"), as(tax_table(NoEdgePlus1.ps)[rownames(DeseqResults), ], "matrix")) #clean up format
  DeseqResults$Habitat <- ifelse(DeseqResults$log2FoldChange<0, "forest", "patch") #make new column specifying which ecosystem each ASV is for
  # Next few lines prepare dataframe with diff abundance results, some sample info, and taxonomic info 
  # View(DeseqResults) #2,169 diff abund ASVs)
  sampDat <- sample_data(postUbiqASVs16S_Plus1.ps)[,c(1,6,8:9)] #Sample.ID, EU, Transect, Meter
  attr(sampDat, "class") <- "data.frame" #remove phyloseq attribute, make dataframe
  sampDat <- tibble::rownames_to_column(sampDat, var="SampleNumberID") #make mapping file ID a column instead of rownames
  ### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
  ASVnamesDA <- rownames(DeseqResults) #2169 found
  ASVsAll <- as.data.frame(t(ASVs_outta_ps(postUbiqASVs16S_Plus1.ps))) #get ASVs from original and transpose so that ASVs are rows 
  # Create dataframe with everything of interest
  diffAbunDat <- merge(DeseqResults, ASVsAll, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
  diffAbunDat_tidy <- diffAbunDat %>% 
    pivot_longer(cols= 16:ncol(diffAbunDat), names_to= "SampleNumberID", values_to= "ASVabundance") %>% #has # of rows equal to # of forest ASVs x sample number
    merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_tidy)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
 #View(diffAbunDat_tidy) this has 505,377 rows, which is equal to 233 (number of samples) x 2,169 (number of diff abund ASVs)
 ##########

##########################################################################
# STACKED BARCHART OF DIFFERENTIAL ABUNDANCE PHYLA NUMBER IN EACH CATEGORY
##########################################################################

DeseqResultsMini <- DeseqResults[,c(8,14)]  #diff abund analysis: just ASV name (as rownames), phylum, and habitat 
postUbiqTaxTab <- taxtable_outta_ps(postUbiqASVs16S_Plus1.ps) #get full taxonomy table from post ubiquity dataset
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == TRUE)) #2169
length(which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE)) #3050 ASVs are NOT differentially abundant?
false_index <- which(rownames(postUbiqTaxTab) %in% rownames(DeseqResultsMini) == FALSE) #also 3050
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

diffAbund_16S_stackedBarplotPhyla <- ggplot(DAphylumAll, aes(fill=Habitat, x=Phylum)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  scale_y_continuous(breaks=seq(0,1000,by=100)) +
  ylab("number of ASVs in phylum") +
  ggtitle("Differentially Abundant Bacterial and Archaeal ASVs")
# quartz()
diffAbund_16S_stackedBarplotPhyla
# Below saved August 1, 2022 so that it can be added to a 2 paneled plot with fungal plot!
# save(diffAbund_16S_stackedBarplotPhyla, file="RobjectsSaved/diffAbund_16S_stackedBarplotPhyla_plot")

# A few checks to make sure that the counting above is working as expected
# Acidobacteria
length(which(DAphylumAll$Phylum=="Acidobacteria")) #1043
acido_index <- which(DAphylumAll$Phylum=="Acidobacteria")
length(which(DAphylumAll[acido_index,]$Habitat=="forest")) #342 forest specialists within Acidobacteria
length(which(DAphylumAll[acido_index,]$Habitat=="patch")) #174 patch specialists within Acidobacteria
length(which(DAphylumAll[acido_index,]$Habitat=="AremainingASVs")) #527 remaining ASVs
(342+ 174 + 527) == length(which(DAphylumAll$Phylum=="Acidobacteria"))

#Chloroflexi
length(which(DAphylumAll$Phylum=="Chloroflexi")) #498
chloro_index <- which(DAphylumAll$Phylum=="Chloroflexi")
length(which(DAphylumAll[chloro_index,]$Habitat=="forest")) #17 forest specialists 
length(which(DAphylumAll[chloro_index,]$Habitat=="patch")) #252 patch specialists
length(which(DAphylumAll[chloro_index,]$Habitat=="AremainingASVs")) #229 remaining ASVs
(17+252+229) == length(which(DAphylumAll$Phylum=="Chloroflexi"))

##########################################################################
# PLOT OF Chloroflexi
##########################################################################

# Within the plot above, Chloroflexi are among the most interesting groups.
# Here, makes a plot of the different relative abundances of 

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
PostUbiq_16S_NOEdge.ps <- subset_samples(postUbiquity.ps, Habitat != "edge") #remove edge so we can compare patch versus forest

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
# relabunpostUbiqNOEdge.phyla.0 <- transform_sample_counts(PostUbiq_16S_NOEdge.ps.phylum.glom, function(x) x / sum(x) )
relabunpostUbiqNOEdge_16S.allASVs <- transform_sample_counts(PostUbiq_16S_NOEdge.ps, function(x) x / sum(x) )

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunpostUbiqNOEdge.allASVs <- merge_samples(relabunpostUbiqNOEdge_16S.allASVs, group = c("Habitat"))
sample_data(relabunpostUbiqNOEdge_16S.allASVs) #shows that we still have samples from each EU, biocrust, and extcontrol (water)
# Meter did an averaging thing; can just ignore it

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunpostUbiqNOEdge_16S.allASVs.2 <- transform_sample_counts(relabunpostUbiqNOEdge.allASVs, function(x) x / sum(x))
sample_data(relabunpostUbiqNOEdge_16S.allASVs.2)

# 
relabunpostUbiqNOEdge_16S.allASVs.df <-psmelt(relabunpostUbiqNOEdge_16S.allASVs.2)
head(relabunpostUbiqNOEdge_16S.allASVs.df) #
dim(relabunpostUbiqNOEdge_16S.allASVs.df) #10438, 33 because 5219 taxa

chloro_index <- which(relabunpostUbiqNOEdge_16S.allASVs.df$Phylum=="Chloroflexi")
chloroRelAbund <- relabunpostUbiqNOEdge_16S.allASVs.df[chloro_index, ]
#View(chloroRelAbund)
colnames(chloroRelAbund)[3] <- "Relative_abundance"
colnames(chloroRelAbund)[2] <- "Habitat_type"

# Make boxplot of relative abundances of each ASV in the Chloroflexi
chloroBoxPlot <- ggplot(chloroRelAbund, aes(x=Habitat_type, y=Relative_abundance, fill= Habitat_type)) +
  geom_boxplot() +
  scale_fill_manual(values=c("darkgreen", "goldenrod")) +
  labs(title="Relative abundance of Chloroflexi ASVs", x="Habitat Type", y = "Relative abundance")
quartz()
chloroBoxPlot


##########################################################################
# GET Z-SCORES OF EACH ASV WITHIN EU
##########################################################################
# 1. Samples based on raw (i.e. not median) abundance
# Prune samples to separate by EU and check to make sure each looks right
# Each has 5219 ASVs (those left after applying ubiquity threshold of 40). However, not all are found in each EU
EU_52_Soils.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, EU == "EU_52")
unique(sample_data(EU_52_Soils.ps)$EU) #only EU 52
EU_53N_Soils.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, EU == "EU_53N")
unique(sample_data(EU_53N_Soils.ps)$EU)
EU_54S_Soils.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, EU == "EU_54S")
unique(sample_data(EU_54S_Soils.ps)$EU)
EU_8_Soils.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, EU == "EU_8")
unique(sample_data(EU_8_Soils.ps)$EU)
EU_53S_Soils.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, EU == "EU_53S")
unique(sample_data(EU_53S_Soils.ps)$EU)
EU_10_Soils.ps <- subset_samples(postUbiqASVs16S_Plus1.ps, EU == "EU_10")
unique(sample_data(EU_10_Soils.ps)$EU)

# 2. Make list of all of these ASV tables by EU:
ASVtabByEU_List <- vector("list", 6) #make a list to hold the ASV tabs from each 
# Fill each element in the list in by each EU (all currently have 5219 taxa, i.e. the whole dataset)
# ASVs are columns, samples are rows; all have 4480 taxa
ASVtabByEU_List[[1]] <- ASVs_outta_ps(EU_52_Soils.ps) #5219 taxa
dim(ASVtabByEU_List[[1]])
ASVtabByEU_List[[2]] <- ASVs_outta_ps(EU_53N_Soils.ps)
dim(ASVtabByEU_List[[2]])
ASVtabByEU_List[[3]] <- ASVs_outta_ps(EU_54S_Soils.ps)
dim(ASVtabByEU_List[[3]])
ASVtabByEU_List[[4]] <- ASVs_outta_ps(EU_8_Soils.ps)
dim(ASVtabByEU_List[[4]])
ASVtabByEU_List[[5]] <- ASVs_outta_ps(EU_53S_Soils.ps)
dim(ASVtabByEU_List[[5]])
ASVtabByEU_List[[6]] <- ASVs_outta_ps(EU_10_Soils.ps)
dim(ASVtabByEU_List[[6]])

# 3. Remove ASVs that are not differentially abundant (calculated across all EUs)
for (j in 1:length(ASVtabByEU_List)) {
  ASVtabByEU_List[[j]] <- ASVtabByEU_List[[j]][ASVnamesDA] #pull out only diff abundant ASVs
  print(dim(ASVtabByEU_List[[j]]))
}
# shows that all now have only 2,169 ASVs, or the differentially abundant ones

# How many NAs, if any, exist at this step?
length(which(is.na(ASVtabByEU_List))) #i think that this shows me what I want, but just in case:
which(is.na(ASVtabByEU_List[[1]]))
which(is.na(ASVtabByEU_List[[2]]))
which(is.na(ASVtabByEU_List[[3]]))
which(is.na(ASVtabByEU_List[[4]]))
which(is.na(ASVtabByEU_List[[5]]))
which(is.na(ASVtabByEU_List[[6]]))
# No NAs anywhere 

######### REMOVING ASVS THAT ARE ZERO IN A GIVEN EU... Removed and greyed out since we 
# looking at all the EUs together, so ASVs that are zero in a given EU are valid. So commenting out for now##
# 4. Remove ASVs that, for any given EU, are not present 
# pre-allocate zero-index list thing
#zeroIndexList <- vector("list", 6) #pre-allocate a new list of 6 vectors, one per EU
#for (k in 1:length(ASVtabByEU_List)) {
  # next line finds each ASV in each EU, that is not present at least once and makes an index for it
#  zeroIndexList[[k]] <- which(rowSums(otu_table(ASVtabByEU_List[[k]]))==0) 
  # line below takes out all of these ASVs that do not appear in each EU.
#  otu_table(ASVtabByEU_List[[k]]) <- otu_table(ASVtabByEU_List[[k]])[-zeroIndexList[[k]]] 
#}
# Check a few to see if code above worked
# These two checks below pull out the rows numbers that are zeros, and the "column" name above it is the ASV name
#which(rowSums(otu_table(ASVtabByEU_List[[1]]))==0)  #yep, no more zeros!
#which(rowSums(otu_table(ASVtabByEU_List[[6]]))==0) #yep, no more zeros!

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

# 4. Rename elements of list by EU
names(ASVtabByEU_List)[[1]] <- "EU_52_z"
names(ASVtabByEU_List)[[2]] <- "EU_53N_z"
names(ASVtabByEU_List)[[3]] <- "EU_54S_z"
names(ASVtabByEU_List)[[4]] <- "EU_8_z"
names(ASVtabByEU_List)[[5]]<- "EU_53S_z"
names(ASVtabByEU_List)[[6]] <- "EU_10_z"

t(ASVtabByEU_List[[1]])  #data.frame objects

# 5. Get each element in the list out of phyloseq and get z-scores
for (k in 1:length(ASVtabByEU_List)){ #ASVtabByEU_List is a list of phyloseq objects that has only the differentially abundant ASVs
  # in each EU
  ASVtabByEU_List[[k]] <- t(ASVtabByEU_List[[k]]) #make ASVs rows and samples columns to work with z-score function
  ASVtabByEU_List[[k]] <- zScore(as.data.frame(ASVtabByEU_List[[k]])) #get z-score for each ASV (where mean and stdev are over all in that EU)
}

##### MORE SCRATCH WORK ########
# This section tests the z-score code above, showing that where the NaNs appear is where a given ASV does not show up at 
# all in an EU. In other words, we get NaNs because if an ASV does not appear in an EU, then all samples within that EU
# are zero, so the standard deviation is zero (meaning all numbers are equal). Dividing by a stdev of zero, as occurs in
# the Z-score calculation, results in a NaN in R. These are real!
# 
# So this is where the NAs/NaNs are coming in
# This is the Z-score function. Let's figure out where it's putting NAs or NaNs in
EU_52_df <- as.data.frame(t(ASVs_outta_ps(EU_52_Soils.ps))) #ASVs are rows, samples are columns
which(is.na(EU_52_df)) #none are NAs at this stage
# Now, pulling apart z score function to see what's going on
  zScoreDftest <- matrix(NA, nrow=nrow(EU_52_df), ncol=ncol(EU_52_df)) #pre-allocate
  colnames(zScoreDftest) <- colnames(EU_52_df)
  rownames(zScoreDftest) <- rownames(EU_52_df)
  for (i in 1:nrow(EU_52_df)){
    ASVmeanTest <- rowMeans(EU_52_df[i,]) #if you change rowMeans to mean, could use with matrices. Or could build in if/else statement
    ASVsdTest <- sd(EU_52_df[i,])
    for (j in 1:length(EU_52_df[i,])) {
      zScoreDftest[i,j] <- (EU_52_df[i,j]-ASVmeanTest)/ASVsdTest # subtract row mean from value, then divide by standard deviation for z score
    }
  }
which(is.nan(zScoreDftest)) == which(is.na(zScoreDftest)) #these are the same. So is.na() finds NAs or NaNs?
length(which(is.nan(zScoreDftest))) #1064 NaNs why?
EU52_NAsIndex <- which(is.nan(zScoreDftest))
EU52_NAsIndexRowCols <- which(is.nan(zScoreDftest), arr.ind = TRUE)
zScoreDftest[EU52_NAsIndexRowCols] #all these are NaNs
EU_52_df.mat <- as.matrix(EU_52_df) #make matrix to use with index
EU_52_df.mat[EU52_NAsIndex] #so all the areas of zeros are where NaNs came in 
# Are there any places were zeros did NOT result in NaNs? YES!
length(which(EU_52_df.mat==0)) #yes #129033
# Maybe it is just those areas where NONE of the samples in this dataset had at least one instance of this ASV.
names(which(rowSums(EU_52_df.mat)==0)) #these ASVs are not present in any sample in this EU
# do these ASVs match up with those giving NaNs above?
unique(rownames(EU52_NAsIndexRowCols)) == names(which(rowSums(EU_52_df.mat)==0))
unique(rownames(EU52_NAsIndexRowCols)) 
################################################


################
####### Scratch work to show that the z-score bit works above:
### if you use this scratch work, just re-create ASVtabByEU_List objects
#ASVtabByEU_List[[1]] <- ASVs_outta_ps(ASVtabByEU_List[[1]])
#class(ASVtabByEU_List[[1]])
#tiny <- ASVtabByEU_List[[1]][1:12, 1:12] #try on a smaller subset
#tiny[tiny == 0] <- 0.01 #it works!
#tiny2 <- zScore(as.data.frame(tiny))

# ASVtabByEU_List[[1]][ASVtabByEU_List[[1]] == 0] <- 0.01

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
# *end scratch work*
############################

#View(ASVtabByEU_List[[2]]) End result is a dataframe with ASV number as row and sample # as column

# 6. # Now I will combine all of the data across EUs to make one dataframe, ( just doing this now b/c z-scores needed to be calculated
# within each EU earlier). In other words, I now merge all these elements of the list by rowname, so that the there are as 
# many rows as ASVs and as many columns as samples (i.e. 233)
## will have to periodically reset rownames to col one, so that can merge by them
merge_1 <- merge(ASVtabByEU_List[[1]], ASVtabByEU_List[[2]], by= "row.names", all=TRUE) #merge first two EUs, all = TRUE merges even those ASVs that don't occur in both
dim(ASVtabByEU_List[[1]])
dim(ASVtabByEU_List[[2]])
dim(merge_1) #1694, 78
#View(merge_1) # now there are some NAs
merge_1 <- merge_1 %>% 
  column_to_rownames(var="Row.names")
head(merge_1)
class(merge_1)
merge_2 <- merge(merge_1, ASVtabByEU_List[[3]], by= "row.names",  all=TRUE)#merge first 2 EUs with EU 3
dim(merge_2) 
merge_2 <- merge_2 %>% 
  column_to_rownames(var="Row.names")
head(merge_2)
merge_3 <- merge(merge_2, ASVtabByEU_List[[4]], by= "row.names", all=TRUE) #merge first 3 EUs with EU 4
dim(merge_3) 
merge_3 <- merge_3 %>% 
  column_to_rownames(var="Row.names")
head(merge_3)
# View(merge_3)
merge_4 <- merge(merge_3, ASVtabByEU_List[[5]], by= "row.names", all=TRUE) #merge first 4 EUs with EU 5
dim(merge_4) #2169 195
merge_4 <- merge_4 %>% 
  column_to_rownames(var="Row.names")
head(merge_4)
merge_5 <- merge(merge_4, ASVtabByEU_List[[6]], by= "row.names", all=TRUE)#merge first 5 EUs with EU 6
dim(merge_5) #2169  234
merge_5 <- merge_5 %>% #
  column_to_rownames(var="Row.names")
head(merge_5) 
dim(merge_5) #now has 233 columns because what was first column is now the rownames
#View(merge_5) 

abundZscores_allEUs <- merge_5
#View(abundZscores_allEUs)

# How many NAs are there in the dataframe? NAs should be places where a particular sample
# did not have that ASV present
length(which(is.na(abundZscores_allEUs)))  #there are 3807 NAs across dataset 
index <- which(is.na(abundZscores_allEUs))
length(index)
dim(abundZscores_allEUs) #2169  233
# another way of looking for NAs
NAs <- sapply(abundZscores_allEUs, function(x) sum(is.na(x)))
sum(NAs) #3807
# Which ASVs have NAs (that is, do not show up in every EU)
abundZscores_allEUs_NAsIndexRowCols <- which(is.na(abundZscores_allEUs), arr.ind = TRUE)
NonOverlappingASVs <- unique(rownames(abundZscores_allEUs_NAsIndexRowCols))
NonOverlappingASVs #these ASVs have some NaNs associated with them

##########################################################################
# LOGISTIC FIT FOR EACH ASV
##########################################################################

# 1. Take abundZscores_allEUs and get relevant metadata so that can do logistic regression thing
head(DeseqResults) #made earlier in this script
sampDat # also made earlier in this script, has sample.ID, EU, Transect, Meter, SampleNumberID
### this repeats a lot of the script from earlier
### HERE IS WHERE WE CHANGE THINGS, SINCE WE DON'T NEED OR WANT SEPARATE DFS FOR FOREST AND PATCH
ASVnamesDA <- rownames(DeseqResults) #2169 found
#ASVsAll <- as.data.frame(t(ASVs_outta_ps(postUbiquity.ps))) #get ASVs from original and transpose so that ASVs are rows 
# Create dataframe with everything of interest
diffAbun_ZDat <- merge(DeseqResults, abundZscores_allEUs, by= "row.names") #grab ASV tab info from only those samples that are differentially abundant
#View(diffAbun_ZDat)
diffAbunDat_Z_tidy <- diffAbun_ZDat %>% 
  pivot_longer(cols= 16:ncol(diffAbun_ZDat), names_to= "SampleNumberID", values_to= "ZabundASV") %>% #has # of rows equal to # of diff abund ASVs x sample number
  merge(sampDat, by="SampleNumberID") #merge with the sampDat to get metadata variables of interest.
colnames(diffAbunDat_Z_tidy)[2] <- "ASV_name" #rename "Row.names" column to be "ASV_name".
# View(diffAbunDat_Z_tidy)

# 2. Get logistic fits for all of these ASVs!
ASVmeterAbunds <- vector("list", length(ASVnamesDA)) #pre-allocate space for each ASV's abundance for each ASV
logFitList <- vector("list", length(ASVnamesDA)) #pre-allocate space for each ASV's logistic fit info
names(ASVmeterAbunds) <- ASVnamesDA
names(logFitList) <- ASVnamesDA
# ASVmeterAbunds for loop has Phylum, family, ZabundASV, Sample.ID (e.g. 53SD_B_10) EU, Transect, and Meter , with a separate dataframe for each ASV
for (i in 1:length(ASVmeterAbunds)){
  tryCatch({
    # the next line loops over all of the differentially abundant ASVs 
    ASVmeterAbunds[[i]] <- diffAbunDat_Z_tidy[which(diffAbunDat_Z_tidy$ASV_name==ASVnamesDA[i]),c(10,13, 17:21)] #keeps phylum, family, ASV abundance (Zscore), sample ID, EU, transect, and meter
    ASVmeterAbunds[[i]] <- ASVmeterAbunds[[i]] %>% 
      dplyr::arrange(Meter) #arrange by meter, in increasing order (i.e. 10, 20, 30, etc. Places 100m where it goes!)
  }, error=function(e){})
}    


# greying out for now for testing purposes
#    logFitList[[i]] <- log.fitdiffAbundFunct(x=ASVmeterAbunds[[i]]$Meter, y=ASVmeterAbunds[[i]]$ASVabundance, ASVnames=ASVnamesDA) #errors are ignored so that those that fit can be fit
# }, error=function(e){})
# }

# does this work? #YES
log.fitdiffAbundFunct(y= ASVmeterAbunds[[1]]$ZabundASV, x= ASVmeterAbunds[[1]]$Meter, ASVnames=names(ASVmeterAbunds)[[1]]) 
  #Error in qr.solve(QR.B, cc) : singular matrix 'a' in solve
# don't know what this issue is, but this is it!!


log.fitdiffAbundFunct(x=ASVmeterAbunds[[1]]$Meter, y=ASVmeterAbunds[[i]]$ASVabundance, ASVnames=ASVnamesDA) # Error in if (any(nEQ <- vNms != make.names(vNms))) vNms[nEQ] <- paste0("`",  : missing value where TRUE/FALSE needed
diffAbunDat_Z_tidy[which(diffAbunDat_Z_tidy$ASV_name==ASVnamesDA[1]),c(10,13, 17:21)]

#### BUT!!!! WE MAY BE ABLE TO PLOT IT ANYWAY!!
# THIS did not Work
DA_ASVsplotList <- vector("list", length(ASVnamesDA)) #pre-allocate list to store all of these plots?
for (i in 1:length(ASVnamesDA)){ #loop over all of the indicator ASVs
  DA_ASVsplotList[[i]] <- plot(ASVmeterAbunds[[i]]$ZabundASV ~ ASVmeterAbunds[[i]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[i]]), xlab= "distance (m)", ylab= "ASV abundance")
}

plot(ASVmeterAbunds[[1]]$ZabundASV ~ ASVmeterAbunds[[1]]$Meter, main = paste("Logistic Function for", names(ASVmeterAbunds)[[1]]), xlab= "distance (m)", ylab= "ASV abundance")


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

