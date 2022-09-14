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
load("RobjectsSaved/diffAbunDat_tidy") #differential abundance results made in EdgeEffectsAllSitesFUNGI.R

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
write.table(FUNGuildTable1, file = "FUNGuildTable.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
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

# Get only those ASVs with probable or highly probable results (i.e. NO possibly)
highProbIndex <- which(rawFungResults$Confidence.Ranking == "Highly Probable") #80 highly probable
probIndex <- which(rawFungResults$Confidence.Ranking == "Probable") #151 probable
fungAssignedIndex <- c(highProbIndex, probIndex)
fungResults <- rawFungResults[fungAssignedIndex,]
# View(fungResults)   #231 out of 257 were probable or highly probable! 

# This first for loop makes it so that ASVs are either labeled differentially abundant or not
diffAbunASVs <- unique(diffAbunDat_tidy$ASV_name)
length(diffAbunASVs)
for (i in 1:nrow(fungResults)) {
  if (fungResults$ASV_name[i] %in% diffAbunASVs) { #so this gives the positions in fungResults$OTU_ID
    # where one of the ASVs in diffAbunASVs is 
    fungResults$DiffAbund[i] <- "DiffAbund"
  } else {
    fungResults$DiffAbund[i] <- "notDiffAbund"
  }
}

length(which(fungResults$ASV_name %in% diffAbunASVs == T)) #159 out of the 287 differentially abundant taxa are probable or highly probable
# 231/413 (all in the post ubiquity dataset were differentially abundant)
231/413*100


# Join with differentially abundant data!
fungResultsDA <- left_join(x=fungResults, y=diffAbunDat_tidy[,c(2,8,10:16)], by="ASV_name")
fungResultsDA <- fungResultsDA %>% dplyr::distinct() #remove duplicate rows
#View(fungResultsDA)
NAindex <- which(is.na(fungResultsDA$Habitat))
fungResultsDA$Habitat[NAindex] <- "notDiffAbund" #make those not diff abund in patch or forest labeled as such!

unique(fungResultsDA$Guild)
unique(fungResultsDA$Trophic.Mode)
View(fungResultsDA)

save(fungResultsDA, file="RobjectsSaved/fungResultsDA") #saved Sept 14, 2022

# Make a stacked bar plot of trophic modes in different habitat types
FUNGuild_stackedBarplotTrophic <- ggplot(fungResultsDA, aes(fill=Habitat, x=Trophic.Mode)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  ylab("number of ASVs") +
  xlab("Trophic mode") +
  ggtitle("Trophic Modes from FUNGuild") 

#quartz()
FUNGuild_stackedBarplotTrophic

# Make a stacked bar plot of different in different habitat types
FUNGuild_stackedBarplotGuild <- ggplot(fungResultsDA, aes(fill=Habitat, x=Guild)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  ylab("number of ASVs") +
  xlab("Guild") +
  ggtitle("Guilds from FUNGuild") 
quartz()
FUNGuild_stackedBarplotGuild

# quartz()
FUNGuild_stackedBarplotGuild
# View(fungResultsDA)

# Make a stacked barplot with only the EMF and AMF represented by habitat type (added to ESA 2022 presentation)
  # 1. Get indices of both
# Get only arbuscular mycorrhizal 
arbMycoIndex <- which(fungResultsDA$Guild == "Arbuscular Mycorrhizal")
length(arbMycoIndex) #19
arbMycoIndex
# Get only EMF
EMFIndex <- which(fungResultsDA$Guild == "Ectomycorrhizal")
length(EMFIndex) #48
AM_EMFindex <- c(arbMycoIndex, EMFIndex)

   # 2. Plot!
FUNGuild_stackedBarplotAM_EMF <- ggplot(fungResultsDA[AM_EMFindex,], aes(fill=Habitat, x=Guild)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_manual(values=c("darkgrey","darkgreen","goldenrod"), name= "Differentially abundant in:", labels=c("not differentially abundant", "forest", "patch")) +
  theme(axis.text.x = element_text(angle = 90), legend.title= element_blank()) +
  ylab("number of ASVs") +
  xlab("Guild") +
  ggtitle("Guilds from FUNGuild") 

# quartz()
FUNGuild_stackedBarplotAM_EMF
