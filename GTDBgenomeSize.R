# Do ASV names match in my stuff and GTDB_matches_SRS_4ntMismatch?

###################################################################################################
# SET UP
###################################################################################################
# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("tidyverse")    
library("vegan")
library("gridExtra")    # allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis
library("grid")


######
# FUNCTIONS DEFINED IN THIS SCRIPT (but often used first in other scripts):
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

# 3. # metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}


# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: UbiquityMedianSetup.R

load("RobjectsSaved/ASVsAllDiffAbund_tax_no100mlog10")

# Make a ggplot with genome size (couldn't do so in server (see file SoilSporersWithMadin.R))
load("RobjectsSaved/genomeSizesDataFrame")
colnames(sizesDf)[1] <- "ASV_name" #re-name for merging further down

GTDB_matches_SRS_4ntMismatch <- read.csv("GTDB_matches_SRS_4ntMismatch.csv")
head(GTDB_matches_SRS_4ntMismatch)
dim(GTDB_matches_SRS_4ntMismatch)

postUbiqTax <- taxtable_outta_ps(postUbiquity.ps)
postUbiqTax2 <- postUbiqTax %>% 
  rownames_to_column(var= "ASV") #make ASVs a column called ASV so that it can be merged with GTDB dataframe by ASV
dim(postUbiqTax2)
dim(GTDB_matches)

GenomesinPostUbiqIndex <- which(GTDB_matches_SRS_4ntMismatch$ASV %in% rownames(postUbiqTax)== TRUE)
GTDB_matches <- GTDB_matches_SRS_4ntMismatch[GenomesinPostUbiqIndex,]
colnames(GTDB_matches) <- c("X", "ASV", "Genome", "Identity", "Length", "Mismatch", "xx",      
                            "rr", "rf", "vv", "bg", "b", "hh", "gKingdom", "gPhylum",
                            "gClass", "gOrder", "gFamily", "gGenus" ) #rename so that they are distinct

postUbiqGTDBmerged <- merge(postUbiqTax2, GTDB_matches, by="ASV", all.y= FALSE)
#View(postUbiqGTDBmerged)
#View(postUbiqGTDBmerged[which(postUbiqGTDBmerged$Mismatch <= 1),]) #length = 345
unique(postUbiqGTDBmerged$gGenus == postUbiqGTDBmerged$Genus) #these are all true (or NA), so we should be able to use it!!!
unique(postUbiqGTDBmerged$gFamily == postUbiqGTDBmerged$Family) #these are all true, so we should be able to use it!!!
unique(postUbiqGTDBmerged$gPhylum == postUbiqGTDBmerged$Phylum) #these are all true, so we should be able to use it!!!

which(rownames(tax_table(postUbiquity.ps))=="ASV_7")
tax_table(postUbiquity.ps)[5,]
which(rownames(tax_table(postUbiquity.ps))=="ASV_989")
tax_table(postUbiquity.ps)[943,]

##########################################################################
# PLOT GENOME SIZES PER SPECIALIST GROUP
##########################################################################
genSizeBoxPlot <- ggplot(sizesDf, aes(x=diffAbundHabitat, y=genome_size, fill=diffAbundHabitat)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("darkgreen", "darkgray", "goldenrod")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  theme(plot.title = element_text(size=18)) +
  theme(legend.key.size = unit(1.7, 'cm')) +
  theme(legend.title=element_blank()) +
  ylab("Genome size") +
  xlab("Specialist group")

quartz()
genSizeBoxPlot

##########################################################################
# ARE TAXA WITH LARGER GENOMES MORE UBIQITOUS?
##########################################################################
# Make a plot with genome size on x, ubiquity on the y.
head(sizesDf)
head(ASVsAllDiffAbund_tax_no100mlog10)
ASVsAllDiffAbund_tax_no100mlog10GenSize <- merge(ASVsAllDiffAbund_tax_no100mlog10, sizesDf, 
                                                 by=c("ASV_name", "Kingdom", "Phylum", "Class", 
                                                      "Order", "Family", "Genus", "Species"), all.x=FALSE)
head(ASVsAllDiffAbund_tax_no100mlog10GenSize)

genSizeUbiqPlot <- ggplot(ASVsAllDiffAbund_tax_no100mlog10GenSize, aes(y=genome_size, x=sampOccurrence, color=diffAbundHabitat.x)) + 
  geom_point() +
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c("darkgray", "darkgreen", "goldenrod")) +
  ylab("Genome Size") +
  xlab("Ubiquity (number of sites occupied)") +
  ggtitle("Relationship between ASV genome size and ubiquity") +
  theme(plot.title = element_text(size=18)) +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  theme_bw()

quartz()
genSizeUbiqPlot

##########################################################################
# WHAT PHYLA ARE REPRESENTED IN EACH OF THE MATCHED GENOMES?
##########################################################################
# Plot it!
groupsOrder <- c("non-specialists", "forest specialists", "patch specialists")
matchedGenomesByPhyla <- ggplot(sizesDf, aes(fill = Phylum, x = diffAbundHabitat)) + 
  geom_bar(position = "stack", stat="count")  +
  scale_fill_manual(values = c("#f781bf", "#a65628", "#b2df8a", "#ffff33", "#ff7f00", "#984ea3", "#a6cee3", "#4daf4a", "#377eb8", "#e41a1c")) +
  labs(y= "number of matched genomes", x= "Specialist group") + 
  ggtitle("Taxonomy of matched genomes") +
  scale_x_discrete(labels=groupsOrder) + #change x-axis tick labels
  theme(axis.title=element_text(size=14)) +
  theme(axis.text.x=element_text(size=12))
# quartz()
matchedGenomesByPhyla

