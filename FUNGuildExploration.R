# FUNGuildExploration.R
# Aug 11, 2022

###################################################################################################
# This script formats data for and then explores FUNGuild Output
###################################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")       
library("tidyverse")       
library("vegan")
library("ggplot2")

# Load data
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after rarefying,
# and keeping only ASVs that occurred at least 50 times across the (rarefied) 
# dataset (from ITS_UbiquityMedianSetup.R).
load(file= "RobjectsSaved/fungiTaxTabWithSpecialists_Nov20_2023") #indicator species analysis results, called fungiTaxTabWithSpecialists, from ITS_indicAnalysisEdgePatchMatrix.R

# FUNCTIONS DEFINED IN THIS SCRIPT:
# 1. taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. ASVs_outta_ps gets the ASV table out of phyloseq 
# from: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  if(taxa_are_rows(ASVTable)) {
    ASVTable <- t(ASVTable)
  }
  return(as.data.frame(ASVTable))
}

# 3. # taxtable_outta_ps takes the phyloseq ASV table and converts it to a dataframe
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

###############################################################################
#         FORMATTING ITS_postUbiquity.ps FOR FUNGuild 
###############################################################################
# format is tsv, where first row should be OTU_IOTU_ID(tab)sample1(tab)sample2(tab)sample3(tab)taxonomy(return)
# taxonomic levels should be separated by semicolons

postUbiqITS_ASVtab <- t(ASVs_outta_ps(ITS_postUbiquity.ps)) #get ASV table out of phyloseq, ASVs are rows
postUbiqITS_TaxTab <- taxtable_outta_ps(ITS_postUbiquity.ps) #get tax table out of phyloseq; ASVs are also rows!
rownames(postUbiqITS_TaxTab) == rownames(postUbiqITS_ASVtab)
FUNGuildTable1 <- merge(postUbiqITS_ASVtab, postUbiqITS_TaxTab, by= "row.names", all=TRUE)
colnames(FUNGuildTable1)[1] <- "OTU_ID"
# Collapse taxonomy into one column 
FUNGuildTable1 <- tidyr::unite(FUNGuildTable1, sep=";",col= "taxonomy", Kingdom:Species)
# View(FUNGuildTable1)
# Remove all of the pesky stuff before taxonomy
FUNGuildTable1$taxonomy <-gsub("k__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("p__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("c__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("o__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("f__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("g__","",as.character(FUNGuildTable1$taxonomy))
FUNGuildTable1$taxonomy <-gsub("s__","",as.character(FUNGuildTable1$taxonomy))
# colnames(FUNGuildTable1)

# Written to file Aug. 10, 2022
# write.table(FUNGuildTable1, file = "FUNGuildTable.txt", sep = "\t",
         #   row.names = TRUE, col.names = TRUE)
# ** Finally, remove any quotation marks that may appear in new file in a text editor (I used Atom) and make
# sure that column names and data line up.

# SEE THE FILE FUNGuildStepByStep for how I ran FUNGuild

###############################################################################
#         Exploring FUNGuild Results
###############################################################################

# load file:
rawFungResults <- read.delim(file= "~/Desktop/CU_Research/SoilEdgeEffectsResearch/FUNGuildTable.guilds_matched.txt", header = T)
# View(rawFungResults) #has X in front of all the sample names, but this should be okay!
# Out of the 413 ASVs that were fed in, 257 had a hit!
colnames(rawFungResults)[1] <- "ASV_name"

# How many indicator fungal taxa could we match?
257/413*100 #62.2276% were matched

# Get only those ASVs with probable or highly probable results (i.e. NO possibly)
highProbIndex <- which(rawFungResults$Confidence.Ranking == "Highly Probable") #80 highly probable
probIndex <- which(rawFungResults$Confidence.Ranking == "Probable") #151 probable
fungAssignedIndex <- c(highProbIndex, probIndex)
fungResults <- rawFungResults[fungAssignedIndex,]
# View(fungResults)   #231 out of 257 were probable or highly probable! 

# How many indicator fungal taxa were matched as highly probable or probable?
231/413*100 #55.9322% were matched

#Join with indicator species analysis results
fungResults_indic <- left_join(x=fungResults, y=fungiTaxTabWithSpecialists, by="ASV_name")
dim(fungResults_indic) #231 rows (equal to ASVs that were matched with FungGuild and 25r columns (i.e., samples
# and other information))
#View(fungResults_indic)

# save(fungResults_indic, file="RobjectsSaved/fungResults_indic") #saved Nov. 20, 2023

# How many ASVs are remaining in each specialist category? These will be added to plot below
length(which(fungResults_indic$specialist== "edge")) #12
length(which(fungResults_indic$specialist== "savanna")) #64
length(which(fungResults_indic$specialist== "forest")) #35
length(which(fungResults_indic$specialist== "non-specialist")) #120


# Make a stacked bar plot of trophic modes in different habitat types
FUNGuild_stackedBarplotTrophic <- ggplot(fungResults_indic, aes(fill=specialist, x=Trophic.Mode)) +
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","goldenrod", "purple","darkgreen"), labels=c("non-specialists (n=120)", "open patch (n= 164)", "edge (n=12)", "forested matrix (n=35)")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  ylab("number of ASVs") +
  xlab("Trophic mode") +
  ggtitle("Trophic Modes from FUNGuild")

# quartz()
FUNGuild_stackedBarplotTrophic
# 
# # Make a stacked bar plot of different in different habitat types
FUNGuild_stackedBarplotGuild <- ggplot(fungResults_indic, aes(fill=specialist, x=Guild)) +
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","goldenrod", "purple","darkgreen"), labels=c("non-specialists (n=120)", "open patch (n= 164)", "edge (n=12)", "forested matrix (n=35)")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  ylab("number of ASVs") +
  xlab("Guild") +
  ggtitle("Guilds from FUNGuild")
# quartz()
FUNGuild_stackedBarplotGuild
# 
# 
# # Make a stacked barplot with only the EMF and AMF represented by habitat type
  # 1. Get indices of both
# Get only arbuscular mycorrhizal
arbMycoIndex <- which(fungResults_indic$Guild == "Arbuscular Mycorrhizal")
length(unique(fungResults_indic$Guild))
guilds <- unique(fungResults_indic$Guild)
# 
length(arbMycoIndex) #19
arbMycoIndex
# Get only EMF
EMFIndex <- which(fungResults_indic$Guild == "Ectomycorrhizal")
length(EMFIndex) #48
AM_EMFindex <- c(arbMycoIndex, EMFIndex)
# 
# What percentage of the fungal reads (that were matched with FUNGuild) were ECM or AM?
colnames(fungResults_indic)
fungResults_indic$ASVabund <- rowSums(fungResults_indic[,2:234]) #get how many times each ASV appears across samples and make this a new row
allSummedFUNGuild <- sum(fungResults_indic$ASVabund) #The abundance of all of these added up is 838202
AMEMPabund <- sum(fungResults_indic$ASVabund[AM_EMFindex]) #The abundance of JUST the ECM/AM is 390251
390251/838202*100 #46.55811% of all of the reads where a guild could be assigned were ECM or AM fungi

# 2. Plot!
FUNGuild_stackedBarplotAM_EMF <- ggplot(fungResults_indic[AM_EMFindex,], aes(fill=specialist, x=Guild)) +
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","goldenrod", "purple","darkgreen"), labels=c("non-specialists (n=120)", "open patch (n= 164)", "edge (n=12)", "forested matrix (n=35)")) +
  scale_x_discrete(labels=c("AM", "ECM")) +
  theme_bw() +
  theme(legend.title= element_blank(), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15)) +
  theme(axis.title.x = element_blank()) + #remove x axis label
  theme(axis.title.y = element_blank()) #remove y axis label

# quartz()
FUNGuild_stackedBarplotAM_EMF

#View(fungResultsDA)

# saved Nov. 22, 2022
#save(fungResultsDA, file= "RobjectsSaved/fungResultsDA")















####### THE REST OF THE SCRIPT BELOW WAS UTILIZED WHEN I WAS STILL DOING DIFFERENTIAL ABUNDANCE
# # This first for loop makes it so that ASVs are either labeled differentially abundant or not
# diffAbunASVs <- unique(diffAbunDat_tidy$ASV_name)
# length(diffAbunASVs)
# for (i in 1:nrow(fungResults)) {
#   if (fungResults$ASV_name[i] %in% diffAbunASVs) { #so this gives the positions in fungResults$OTU_ID
#     # where one of the ASVs in diffAbunASVs is 
#     fungResults$DiffAbund[i] <- "DiffAbund"
#   } else {
#     fungResults$DiffAbund[i] <- "notDiffAbund"
#   }
# }
# 
# length(which(fungResults$ASV_name %in% diffAbunASVs == T)) #159 out of the 287 differentially abundant taxa are probable or highly probable
# # 231/413 (all in the post ubiquity dataset were differentially abundant)
# 231/413*100
# 
# 
# # Join with differentially abundant data!
# fungResultsDA <- left_join(x=fungResults, y=diffAbunDat_tidy[,c(2,8,10:16)], by="ASV_name")
# fungResultsDA <- fungResultsDA %>% dplyr::distinct() #remove duplicate rows
# #View(fungResultsDA)
# NAindex <- which(is.na(fungResultsDA$Habitat))
# fungResultsDA$Habitat[NAindex] <- "notDiffAbund" #make those not diff abund in patch or forest labeled as such!
# 
# unique(fungResultsDA$Guild)
# unique(fungResultsDA$Trophic.Mode)
# #View(fungResultsDA)
# 
# #save(fungResultsDA, file="RobjectsSaved/fungResultsDA") #saved Sept 14, 2022
# 
# # Make a stacked bar plot of trophic modes in different habitat types
# FUNGuild_stackedBarplotTrophic <- ggplot(fungResultsDA, aes(fill=Habitat, x=Trophic.Mode)) + 
#   geom_bar(position="stack", stat="count") +
#   scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
#   theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
#   ylab("number of ASVs") +
#   xlab("Trophic mode") +
#   ggtitle("Trophic Modes from FUNGuild") 
# 
# # quartz()
# FUNGuild_stackedBarplotTrophic
# 
# # Make a stacked bar plot of different in different habitat types
# FUNGuild_stackedBarplotGuild <- ggplot(fungResultsDA, aes(fill=Habitat, x=Guild)) + 
#   geom_bar(position="stack", stat="count") +
#   scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
#   theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
#   ylab("number of ASVs") +
#   xlab("Guild") +
#   ggtitle("Guilds from FUNGuild") 
# # quartz()
# FUNGuild_stackedBarplotGuild
# 
# # quartz()
# FUNGuild_stackedBarplotGuild
# # View(fungResultsDA)
# 
# # Make a stacked barplot with only the EMF and AMF represented by habitat type (added to ESA 2022 presentation)
#   # 1. Get indices of both
# # Get only arbuscular mycorrhizal 
# arbMycoIndex <- which(fungResultsDA$Guild == "Arbuscular Mycorrhizal")
# length(unique(fungResultsDA$Guild))
# guilds <- unique(fungResultsDA$Guild)
# 
# length(arbMycoIndex) #19
# arbMycoIndex
# # Get only EMF
# EMFIndex <- which(fungResultsDA$Guild == "Ectomycorrhizal")
# length(EMFIndex) #48
# AM_EMFindex <- c(arbMycoIndex, EMFIndex)
# 
# # What percentage of the fungal reads (that were matched with FUNGuild) were ECM or AM?
# colnames(fungResultsDA)
# fungResultsDA$ASVabund <- rowSums(fungResultsDA[,2:234]) #get how many times each ASV appears across samples and make this a new row
# allSummedFUNGuild <- sum(fungResultsDA$ASVabund) #The abundance of all of these added up is 838202
# AMEMPabund <- sum(fungResultsDA$ASVabund[AM_EMFindex]) #The abundance of JUST the ECM/AM is 390251
# 390251/838202*100 #46.55811% of all of the reads where a guild could be assigned were ECM or AM fungi
# 
# 
#    # 2. Plot!
# FUNGuild_stackedBarplotAM_EMF <- ggplot(fungResultsDA[AM_EMFindex,], aes(fill=Habitat, x=Guild)) + 
#   geom_bar(position="stack", stat="count") +
#   scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), labels=c("non-specialists", "forest", "savanna")) +
#   scale_x_discrete(labels=c("AM", "ECM")) +
#   theme_bw() + 
#   theme(legend.title= element_blank(), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15)) +
#   theme(axis.title.x = element_blank()) + #remove x axis label
#   theme(axis.title.y = element_blank()) #remove y axis label
# 
# # quartz()
# FUNGuild_stackedBarplotAM_EMF
# 
# #View(fungResultsDA)
# 
# # saved Nov. 22, 2022
# #save(fungResultsDA, file= "RobjectsSaved/fungResultsDA")

