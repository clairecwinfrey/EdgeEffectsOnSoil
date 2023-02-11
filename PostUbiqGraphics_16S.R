# 16s- PostUbiquity Graphics 
# started August 2, 2022

# This script re-does ordinations and Bray-Curtis analyses with the postUbiquity analysis for prokaryote samples. 
# In addition, it creates stacked barplots to show the top phyla and families (using postUbiquity.ps, see description below)
# Finally, it performs a PERMANOVA to test if samples differ based on habitat type.

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
library("utilities")

# Load data
load("RobjectsSaved/postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# This object was created and saved in UbiquityMedianSetup.R


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

###############################
# TOP PHYLA PLOT
###############################
# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)

# Combine taxa based on phylum
postUbiq.phylum.glom <-  tax_glom(postUbiquity.ps, taxrank = "Phylum") 
tax_table(postUbiq.phylum.glom) #24 phyla
sample_data(postUbiq.phylum.glom)

# Transform sample counts based on just glommed samples
relabun.phyla.0 <- transform_sample_counts(postUbiq.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0)) 
colSums(otu_table(relabun.phyla.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each phylum

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = c("EU"))
sample_data(relabun.phyla.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
rowSums(otu_table(relabun.phyla.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5 or 1% of total abundance
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
# Make copies 
relabun.phylatop99 <- relabun.phyla.df
relabun.phylatop99.5 <- relabun.phyla.df

relabun.phylatop99$Phylum[relabun.phylatop99$Abundance < 0.01] <- "< 1% abund."
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abund."

top_99p_phyla <- unique(relabun.phylatop99$Phylum)
top_99p_phyla

top_99.5p_phyla <- unique(relabun.phylatop99.5$Phylum)
top_99.5p_phyla

# Surprised that Gemmatimonadetes not above >1%! But Gemmatimonadetes in top 16 and
# comprises at least 0.5% of the total abundance

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum))
phylumPlot99.5percent <- phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.y= element_text(size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=14)) +
  guides(fill=guide_legend(nrow=5)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

#quartz()
phylumPlot99.5percent
# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

top_99.5p_phyla <- relabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()
#View(top_99.5p_phyla)

# This gets the top phyla's relative abundance in each group.
# It shows that the top 4 are Acidobcteria, Proteobacteria, Verrucomicrobia,
# and Planctomycetes, nd Chloroflexi and Actinos are comaprable 
relabun.phylatop99.5_grouped <- relabun.phylatop99.5 %>% 
  group_by(Phylum, Sample) %>% 
  summarize(relAbund = sum(Abundance))
# View(relabun.phylatop99.5_grouped)

# This gets relative abundnce across all EUs
relabun.phyla_grouped2 <- relabun.phylatop99.5 %>% 
  group_by(Phylum) %>% 
  summarize(relAbundSummedAcrossEUs = sum(Abundance)) %>% 
  mutate(relAbund = relAbundSummedAcrossEUs/6)
#View(relabun.phyla_grouped2)

sum(relabun.phyla_grouped2$relAbund) #this adds up to 1!

###############################
# TOP FAMILIES PLOT
###############################
# Combine taxa based on family
postUbiq.family.glom <-  tax_glom(postUbiquity.ps, taxrank = "Family") 
tax_table(postUbiq.family.glom) # 208 families
sample_data(postUbiq.family.glom)

# Transform sample counts based on just glommed samples
relabun.fam.0 <- transform_sample_counts(postUbiq.family.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.fam.0)) 
colSums(otu_table(relabun.fam.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each family

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.fam.1 <- merge_samples(relabun.fam.0, group = c("EU"))
sample_data(relabun.fam.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.fam.2 <- transform_sample_counts(relabun.fam.1, function(x) x / sum(x))
rowSums(otu_table(relabun.fam.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least 1 or 5% of total abundance
relabun.fam.df <-psmelt(relabun.fam.2)
dim(relabun.fam.df) #
# Make copies 
relabun.famtop99 <- relabun.fam.df
relabun.famtop95 <- relabun.fam.df

relabun.famtop99$Family[relabun.famtop99$Abundance < 0.01] <- "< 1% abund."
relabun.famtop95$Family[relabun.famtop95$Abundance < 0.05] <- "< 5% abund."

top_99p_fam <- unique(relabun.famtop99$Family)
top_99p_fam 

top_95p_fam <- unique(relabun.famtop95$Family)
top_95p_fam #9

# PLOT FOR FAMILIES COMPRISING AT LEAST 1% of abundance
# For aesthetic/fitting reasons, remove subgroup info for two of the families from legend
relabun.famtop99$Family <- gsub("_(Subgroup_3)","",as.character(relabun.famtop99$Family),fixed = T) #remove subgroup from Solibacteraceae
relabun.famtop99$Family <- gsub("_(Subgroup_1)","",as.character(relabun.famtop99$Family),fixed = T) #remove subgroup from Acidobacteriaceae

familyPlot99percent <- ggplot(data=relabun.famtop99, aes(x=Sample, y=Abundance, fill=Family))

familyPlot99percent <- familyPlot99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.y= element_text(size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=14)) +
  guides(fill=guide_legend(nrow=8)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

#quartz()
familyPlot99percent
# Get exact abundances of each family (top 99%) in each sample:
top_99_fam <- relabun.famtop99 %>%
  group_by(Sample, Family) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

# This gets the top families's relative abundance in each group.
# It shows that the top 4 are Acidobcteria, Proteobacteria, Verrucomicrobia,
# and Planctomycetes, nd Chloroflexi and Actinos are comaprable 
relabun.famtop99_grouped <- relabun.famtop99 %>% 
  group_by(Family, Sample) %>% 
  summarize(relAbund = sum(Abundance))
# View(relabun.famtop99_grouped)

# This gets relative abundance across all EUs
relabun.famtop99_grouped2 <- relabun.famtop99 %>% 
  group_by(Family) %>% 
  summarize(relAbundSummedAcrossEUs = sum(Abundance)) %>% 
  mutate(relAbund = relAbundSummedAcrossEUs/6)
#View(relabun.famtop99_grouped2)

sum(relabun.famtop99_grouped2$relAbund) #this adds up to 1!
#################################
# ORDINATIONS
#################################
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ord_16S <- ordinate(postUbiquity.ps, method = "NMDS", distance = "bray", trymax = 100) 
# stress = 0.1105877 

# Ordination based on "habitat"/ plant community type 
HabitatNMDS_16S_postUbiq <- phyloseq::plot_ordination(postUbiquity.ps, ord_16S, type= "samples", color= "Habitat")
HabitatNMDS_16S_postUbiq <- HabitatNMDS_16S_postUbiq +
  scale_color_manual(values=c("purple", "darkgreen", "goldenrod"), #change color of points
                     breaks= c("edge", "forest", "patch"),
                     labels= c("edge", "forest", "savanna")) + #rename categories
  geom_point(size=3) + ggtitle("Bacteria and Archaea") +
  theme_bw() +
  theme(plot.title = element_text(size=18)) +
  guides(color = guide_legend(override.aes = list(size =7))) +
  theme(legend.title = element_text(size=14)) + #increase size of legend title
  theme(legend.text = element_text(size = 12)) + #increase size of legend font
  theme(axis.text=element_text(size=12), #increase font size of axis labels
        axis.title=element_text(size=12)) #increase font size of axis title

# + #change size of legend and title
#theme(axis.ticks = element_blank(), #remove x and y axis labels and tick marks
#     axis.text = element_blank())

#quartz()
HabitatNMDS_16S_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# save plot
save(HabitatNMDS_16S_postUbiq, file="RobjectsSaved/HabitatNMDS_16S_postUbiq") #saved Jan. 2, 2023

# Ordination based on EU (just to show that they are different!)
EU_NMDS_16S_postUbiq <- phyloseq::plot_ordination(postUbiquity.ps, ord_16S, type= "samples", color= "EU")
EU_NMDS_16S_postUbiq <- EU_NMDS_16S_postUbiq + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
EU_NMDS_16S_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# load ITS habitat ordination plot (as made in PostUbiqGraphics_ITS.R)
load(file="RobjectsSaved/diffAbund_ITS_stackedBarplotPhyla_plot")
#quartz()
grid.arrange(HabitatNMDS_16S_postUbiq, HabitatNMDS_ITS_postUbiq, ncol=2)

#########################################################################
# PERMANOVA TO TEST FOR THE EFFECT OF HABITAT TYPE ON THE DISTRIBUTIONS
#########################################################################

# First, need to get data out of phyloseq so that it can be processed in vegan
postUbiq_16SASVs <- ASVs_outta_ps(postUbiquity.ps) #get ASV table
str(postUbiq_16SASVs) #samples are rows, taxa are columns as they should be for vegan
# Getting rid of the phyloseq attribute for metadata is near impossible, so I have to hack it a bit
# by grabbing each column that I need (as vectors or strings), then combining these
sampNames <- rownames(sample_data(postUbiquity.ps))
rownames(postUbiq_16SASVs) == sampNames #since these match, we are good to go!
habitat <- sample_data(postUbiquity.ps)$Habitat
sampID <- sample_data(postUbiquity.ps)$Sample.ID
metaDatForVeg <- cbind(sampID, habitat)
rownames(metaDatForVeg) <- sampNames
metaDatForVeg <- as.data.frame(metaDatForVeg)

# Next get Bray-Curtis dissimilarities
postUbiq_16SBC <- vegdist(postUbiq_16SASVs, method= "bray")
set.seed(19) #set seed so that results are reproducible!
# Run the PERMANOVA
postUbiq_16SPermanova <- adonis(postUbiq_16SBC ~ habitat, data=metaDatForVeg, permutations= 99999)
# p < 0.001

# Save PERMANOVA results
# save(postUbiq_16SPermanova, file= "RobjectsSaved/postUbiq_16SPermanova") #last saved Aug 9, 2022

# WHAT ABOUT A DB-RDA-- SAME AS THE PERMANOVA
mod1 <- dbrda(postUbiq_16SBC ~ habitat, data=metaDatForVeg)

dbRDAprokHabitatResults <- anova.cca(mod1, permutations= 99999)
# this shows the exact same thing as the PERMANOVA

# Save DBRDA results
# save(dbRDAprokHabitatResults, file= "RobjectsSaved/dbRDAprokHabitatResults") #last saved Dec. 14, 2022
