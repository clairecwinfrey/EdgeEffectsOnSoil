# ITS- PostUbiquity Graphics 
# started August 2, 2022

# This script re-does ordinations and Bray-Curtis analyses with the postUbiquity analysis for fungal samples.
# In addition, it creates stacked barcharts of the relative abundance of top phyla and families, and performs PERMANOVAs 
# and db-RDAs nd makes visualizations to explore how sample segregate by habitat type (db-RDAs performed to test both 
# between habitats across sites and between habitats within sites)

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
library("Polychrome") #for color palette

# Load data
load("RobjectsSaved/ITS_postUbiquity.ps") #load phyloseq object that was made after 1) rarefying,
# the 2) keeping only ASVs that occurred at least 50 times across the (rarefied), 
# then 3) filtering out ASVs with a ubiquity of <40
# R OBJECT MADE IN: ITS_UbiquityMedianSetup.R

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
postUbiq.phylum.glom <-  tax_glom(ITS_postUbiquity.ps, taxrank = "Phylum") 
tax_table(postUbiq.phylum.glom) #8 phyla
sample_data(postUbiq.phylum.glom)

# Transform sample counts based on just glommed samples
relabun.phyla.0 <- transform_sample_counts(postUbiq.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0)) 
colSums(otu_table(relabun.phyla.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each phylum

###### TOP PHYLA BY EU #####

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = c("EU"))
sample_data(relabun.phyla.1) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
rowSums(otu_table(relabun.phyla.2)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5% of total abundance
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
relabun.phylatop99.5  <- relabun.phyla.df
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abundance"

# Olpidiomycota and Rozellomycota were phyla as well, but comprised less than .5% of abundance

# Phyla comprising at least 0.5% of total abundance
# Get unique colors for these 13 groups
colors31 <- Polychrome::glasbey.colors(32)
swatch(colors31)
colorsForPhyla7 <- unname(colors31)[c(2:4,6:8)] #remove some unwanted colors
colorsForPhyla7[7] <- "gray48" #make last one gray

# Remove p_ in phylum names
relabun.phylatop99.5$Phylum <- gsub("p__","",as.character(relabun.phylatop99.5$Phylum))
uniquePhyla <- sort(unique(relabun.phylatop99.5$Phylum)) #7 unique phyla

# re-order phyla so that "< .5% abundance" is last
relabun.phylatop99.5$Phylum <- factor(relabun.phylatop99.5$Phylum, levels=uniquePhyla[c(2:7,1)])

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum))
phylumPlot99.5percent <- phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForPhyla7) +
  theme(axis.title.y = element_text(size = 16)) +
  ylab("Relative abundance") +
  theme(axis.text.y= element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=3)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
phylumPlot99.5percent
# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

top_99.5p_phyla <- relabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

# This gets relative abundance across all EUs
fungiRelabunFam_grouped <- relabun.phylatop99.5 %>% 
  group_by(Phylum) %>% 
  summarize(relAbundSummedAcrossEUs = sum(Abundance)) %>% 
  mutate(relAbund = relAbundSummedAcrossEUs/6)
#View(fungiRelabunFam_grouped)
sum(fungiRelabunFam_grouped$relAbund) #this adds up to 1!

###### TOP PHYLA BY HABITAT #####

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.phyla.1_habitat <- merge_samples(relabun.phyla.0, group = c("Habitat"))
sample_data(relabun.phyla.1_habitat) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.phyla.2_habitat <- transform_sample_counts(relabun.phyla.1_habitat, function(x) x / sum(x))
rowSums(otu_table(relabun.phyla.2_habitat)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least .5% of total abundance
relabun.phyla.df_habitat <-psmelt(relabun.phyla.2_habitat)
dim(relabun.phyla.df_habitat) #
relabun.phylatop99.5_habitat  <- relabun.phyla.df_habitat
relabun.phylatop99.5_habitat$Phylum[relabun.phylatop99.5_habitat$Abundance < 0.005] <- "< .5% abundance"

# Phyla comprising at least 0.5% of total abundance
sort(unique(relabun.phylatop99.5_habitat$Phylum))

# Remove p_ in phylum names
relabun.phylatop99.5_habitat$Phylum <- gsub("p__","",as.character(relabun.phylatop99.5_habitat$Phylum))
uniquePhyla_habitat <- sort(unique(relabun.phylatop99.5_habitat$Phylum)) #7 unique phyla

# re-order phyla so that "< .5% abundance" is last
relabun.phylatop99.5_habitat$Phylum <- factor(relabun.phylatop99.5_habitat$Phylum, levels=uniquePhyla_habitat[c(2:7,1)])

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent_habitat <- ggplot(data=relabun.phylatop99.5_habitat, aes(x=Sample, y=Abundance, fill=Phylum))
phylumPlot99.5percent_habitat <- phylumPlot99.5percent_habitat + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForPhyla7) +
  theme(axis.title.y = element_text(size = 16)) +
  ylab("Relative abundance") +
  scale_x_discrete(labels=c("edge", "forest", "savanna")) +
  theme(axis.text.y= element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=3)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
phylumPlot99.5percent_habitat

###############################
# TOP FAMILIES PLOT
###############################
# Combine taxa based on family
postUbiq.family.glom <-  tax_glom(ITS_postUbiquity.ps, taxrank = "Family") 
tax_table(postUbiq.family.glom) # 107 families total
sample_data(postUbiq.family.glom)

# Transform sample counts based on just glommed samples
relabun.fam.0 <- transform_sample_counts(postUbiq.family.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.fam.0)) 
colSums(otu_table(relabun.fam.0)) #right now these all sum to one, which shows that right now, it is relative
# abundance by sample and that the code is working as expected.
# ASVs are just representative from each family

###### TOP FAMILIES BY EU #####

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

relabun.famtop99$Family[relabun.famtop99$Abundance < 0.01] <- "< 1% abundance"
relabun.famtop95$Family[relabun.famtop95$Abundance < 0.05] <- "< 5% abundance"

top_99p_fam <- unique(relabun.famtop99$Family)
top_99p_fam #28 families

top_95p_fam <- unique(relabun.famtop95$Family)
top_95p_fam #10 families

# Remove f__ in family names
relabun.famtop99$Family <- gsub("f__","",as.character(relabun.famtop99$Family))
relabun.famtop95$Family <- gsub("f__","",as.character(relabun.famtop95$Family))

# Clean up Mucoromycotina_fam_Incertae_sedis
relabun.famtop99$Family <- gsub("_fam_Incertae_sedis",", inc. sed.",as.character(relabun.famtop99$Family))
relabun.famtop95$Family <- gsub("_fam_Incertae_sedis",", inc. sed.",as.character(relabun.famtop95$Family))
relabun.famtop99$Family <- gsub("NA","Unclassified",as.character(relabun.famtop99$Family))

# Family comprising at least 1% of total abundance
uniqueFams <- sort(unique(relabun.famtop99$Family)) #29 unique families
# Get unique colors for these 13 groups
colors31 <- Polychrome::glasbey.colors(32)
swatch(colors31)
colorsForFams29 <- unname(colors31)[c(2:30)] #remove some unwanted colors
# switch around a few colors
colorsForFams29[5] <- colorsForFams29[6]
colorsForFams29[6] <- unname(colors31)[5]
colorsForFams29[25] <- unname(colors31)[30]
colorsForFams29[29] <- "gray48" #make last one gray

# re-order families so that "Unclassified" and "< 1% abundance" are last
relabun.famtop99$Family <- factor(relabun.famtop99$Family, levels=uniqueFams[c(2:28,29,1)])

# PLOT FOR FAMILIES COMPRISING AT LEAST 1% of abundance
familyPlot99percent <- ggplot(data=relabun.famtop99, aes(x=Sample, y=Abundance, fill=Family))

familyPlot99percent <- familyPlot99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForFams29) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.y= element_text(size=16)) +
  ylab("Relative abundance") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=10)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
familyPlot99percent
# Get exact abundances of each family (top 99%) in each sample:
colnames(relabun.famtop99)
relabun.famtop99[,2] #This is EU... now "Sample" because of the glomming!

top_99_fam <- relabun.famtop99 %>%
  group_by(Sample, Family) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

# This gets relative abundance across all EUs
fungiTopFamsGrouped <- relabun.famtop99 %>% 
  group_by(Family) %>% 
  summarize(relAbundSummedAcrossEUs = sum(Abundance)) %>% 
  mutate(relAbund = relAbundSummedAcrossEUs/6)
#View(fungiTopFamsGrouped)
sum(fungiTopFamsGrouped$relAbund) #this adds up to 1!

###### TOP FAMILIES BY FOREST/SAVANNA #####

# Merge samples so that we only have combined abundances for EU (i.e. experimental replicate)
relabun.fam.1_habitat <- merge_samples(relabun.fam.0, group = c("Habitat"))
sample_data(relabun.fam.1_habitat) #this just confirms that samples were combined by EU, other variables are averaged but can be ignored

# Convert to proportions again after merging samples by EU.
relabun.fam.2_habitat <- transform_sample_counts(relabun.fam.1_habitat, function(x) x / sum(x))
rowSums(otu_table(relabun.fam.2_habitat)) #these show that all of the ASVs now sum to one, which is exactly what we want!

# Get taxa that that are at least 1 or 5% of total abundance
relabun.fam.df_habitat <-psmelt(relabun.fam.2_habitat)
dim(relabun.fam.df_habitat) #
# Make copies 
relabun.famtop99_habitat <- relabun.fam.df_habitat
relabun.famtop95_habitat <- relabun.fam.df_habitat

relabun.famtop99_habitat$Family[relabun.famtop99_habitat$Abundance < 0.01] <- "< 1% abundance"
relabun.famtop95_habitat$Family[relabun.famtop95_habitat$Abundance < 0.05] <- "< 5% abundance"

top_99p_fam_habitat <- unique(relabun.famtop99_habitat$Family)
top_99p_fam_habitat #26 families

top_95p_fam_habitat <- unique(relabun.famtop95_habitat$Family)
top_95p_fam_habitat #8 families

# Remove f__ in family names
relabun.famtop99_habitat$Family <- gsub("f__","",as.character(relabun.famtop99_habitat$Family))

# Clean up NAs
relabun.famtop99_habitat$Family <- gsub("NA","Unclassified",as.character(relabun.famtop99_habitat$Family))

# Family comprising at least 1% of total abundance
top_99p_fam_habitat <- sort(unique(relabun.famtop99_habitat$Family))
# Get unique colors for these 26 groups
swatch(colors31)
colorsForFams26 <- unname(colors31)[c(2:26)] #remove some unwanted colors
# switch around a few colors
colorsForFams26[26] <- "gray48" #make last one gray

# re-order families so that "Unclassified" and "< 1% abundance" are last
relabun.famtop99_habitat$Family <- factor(relabun.famtop99_habitat$Family, levels=top_99p_fam_habitat[c(2:25,26,1)])

# Plot it
familyPlot99percent_Habitat <- ggplot(data=relabun.famtop99_habitat, aes(x=Sample, y=Abundance, fill=Family))

familyPlot99percent_Habitat <- familyPlot99percent_Habitat + geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(values = colorsForFams26) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.y= element_text(size=16)) +
  ylab("Relative abundance") +
  scale_x_discrete(labels=c("edge", "forest", "savanna")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x= element_text(size=16)) +
  guides(fill=guide_legend(nrow=10)) + theme(legend.text = element_text(colour="black", size = 14)) +
  theme(legend.title= element_blank()) #remove legend title

# quartz()
familyPlot99percent_Habitat

##### CALCULATE TOP FAMILIES IN EACH HABITAT ####
colnames(relabun.fam.df_habitat)
relabun.fam.df_habitat[,2] #"Sample" here is now habitat type, because of the glomming done above
length(unique(relabun.fam.df_habitat$Family)) #82

topFamsHabitat_ITS <- relabun.fam.df_habitat %>%
    select(Sample, Abundance, Phylum, Class, Order, Family) %>% 
    group_by(Sample, Family) %>%
    arrange(-Abundance)

# View(topFamsHabitat_ITS)

# Which are the top in each habitat type? This went into the discussion.
topFamsForest_ITS <- topFamsHabitat_ITS[which(topFamsHabitat_ITS$Sample == "forest"),]
topFamsPatch_ITS <- topFamsHabitat_ITS[which(topFamsHabitat_ITS$Sample == "patch"),]
sum(topFamsPatch_ITS$Abundance)
sum(topFamsForest_ITS$Abundance)

# This gets relative abundance across ALL habitats (i.e. Habitat distinction is ignored)
fungiTopFamsGrouped_Habs <- relabun.fam.df_habitat %>% 
  group_by(Family) %>% 
  summarize(relAbundSummedAcrossHabs = sum(Abundance)) %>% 
  mutate(relAbund = relAbundSummedAcrossHabs/3)
# View(fungiTopFamsGrouped_Habs)
sum(fungiTopFamsGrouped_Habs$relAbund) #this adds up to 1!

##### CALCULATE TOP PHYLA IN EACH HABITAT ####
colnames(relabun.phyla.df_habitat)
relabun.phyla.df_habitat[,2] #"Sample" here is now habitat type, because of the glomming done above
length(unique(relabun.phyla.df_habitat$Phylum)) #8

topPhylaHabitat_ITS <- relabun.phyla.df_habitat %>%
  select(Sample, Abundance, Phylum) %>% 
  group_by(Sample, Phylum) %>%
  arrange(-Abundance)
topPhylaHabitat_ITS

# Which are the top in each habitat type? This went into the discussion.
topPhylaForest_ITS <- topPhylaHabitat_ITS[which(topPhylaHabitat_ITS$Sample == "forest"),]
topPhylaPatch_ITS <- topPhylaHabitat_ITS[which(topPhylaHabitat_ITS$Sample == "patch"),]
sum(topPhylaForest_ITS$Abundance)
sum(topPhylaPatch_ITS$Abundance)
topPhylaForest_ITS
topPhylaPatch_ITS

# PLOT TOP PHYLUM AND TOP FAMILIES BY EU PLOT SIDE BY SIDE
# quartz()
grid.arrange(phylumPlot99.5percent, familyPlot99percent, ncol=2)

# PLOT TOP PHYLA AND TOP FAMILIES BY HABITAT SIDE BY SIDE
# quartz()
grid.arrange(phylumPlot99.5percent_habitat, familyPlot99percent_Habitat, ncol=2)

# Saved May 12, 2023
# save(phylumPlot99.5percent_habitat, file= "RobjectsSaved/fungiPhylumPlot99.5percent_habitat")

#################################
# ORDINATIONS
#################################
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ord_ITS <- ordinate(ITS_postUbiquity.ps, method = "NMDS", distance = "bray", trymax = 100) 
# stress= 0.1484425 

# Ordination based on "habitat"/ plant community type 
HabitatNMDS_ITS_postUbiq <- phyloseq::plot_ordination(ITS_postUbiquity.ps, ord_ITS, type= "samples", color= "Habitat")
HabitatNMDS_ITS_postUbiq <- HabitatNMDS_ITS_postUbiq +
  scale_color_manual(values=c("purple", "darkgreen", "goldenrod"), #change color of points
                      breaks= c("edge", "forest", "patch"),
                      labels= c("edge", "forest", "savanna")) + #rename categories
  geom_point(size=3) + ggtitle("Fungi") +
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

# quartz()
HabitatNMDS_ITS_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!

# Save the plot made above (saved Jan.2, 2022) so that we can make a two paneled figure with ITS stuff
#save(HabitatNMDS_ITS_postUbiq, file="RobjectsSaved/HabitatNMDS_ITS_postUbiq")

### MAKEA TWO PANELED PLOT WITH THIS AND FUNGAL ORDINATION PLOT SIDE BY SIDE #####
# LOAD 16S plot made in PostUbiqGraphics_16S.R
load(file="RobjectsSaved/HabitatNMDS_16S_postUbiq")
## quartz()
grid.arrange(HabitatNMDS_16S_postUbiq, HabitatNMDS_ITS_postUbiq, ncol=2)
# Note: reflected the fungal NMDS plot in PowerPoint to make patch and forest on the same side! 

# Ordination based on EU (just to show that they are different!)
EU_NMDS_ITS_postUbiq <- phyloseq::plot_ordination(ITS_postUbiquity.ps, ord_ITS, type= "samples", color= "EU")
EU_NMDS_ITS_postUbiq <- EU_NMDS_ITS_postUbiq + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
EU_NMDS_ITS_postUbiq #cool, you can see that the forest and the patch separate out, with edge somewhat in between!


#########################################################################
# PERMANOVA AND DBRDAS TO TEST FOR THE EFFECT OF HABITAT TYPE ON THE DISTRIBUTIONS
#########################################################################

# First, need to get data out of phyloseq so that it can be processed in vegan
postUbiq_ITSASVs <- ASVs_outta_ps(ITS_postUbiquity.ps) #get ASV table
str(postUbiq_ITSASVs) #samples are rows, taxa are columns as they should be for vegan
# Getting rid of the phyloseq attribute for metadata is near impossible, so I have to hack it a bit
# by grabbing each column that I need (as vectors or strings), then combining these
sampNames <- rownames(sample_data(ITS_postUbiquity.ps))
rownames(postUbiq_ITSASVs) == sampNames #since these match, we are good to go!
habitat <- sample_data(ITS_postUbiquity.ps)$Habitat
sampID <- sample_data(ITS_postUbiquity.ps)$Sample.ID
metaDatForVeg <- cbind(sampID, habitat)
rownames(metaDatForVeg) <- sampNames
metaDatForVeg <- as.data.frame(metaDatForVeg)

# Next get Bray-Curtis dissimilarities
postUbiq_ITSBC <- vegdist(postUbiq_ITSASVs, method= "bray")
set.seed(19) #set seed so that results are reproducible!
# Run the PERMANOVA
postUbiq_ITSPermanova <- adonis(postUbiq_ITSBC ~ habitat, data=metaDatForVeg, permutations= 9999)
# p < 0.001

# Save PERMANOVA results
# save(postUbiq_ITSPermanova, file= "RobjectsSaved/postUbiq_ITSPermanova") #last saved Aug 9, 2022

# WHAT ABOUT A DB-RDA-- SAME AS THE PERMANOVA
mod1Fungi <- dbrda(postUbiq_ITSBC ~ habitat, data=metaDatForVeg)

dbRDAfungiHabitatResults <- anova.cca(mod1Fungi, permutations= 9999)
# this shows the exact same thing as the PERMANOVA

# dbRDAfungiHabitatResults MODEL OUTPUT:
# Model: dbrda(formula = postUbiq_ITSBC ~ habitat, data = metaDatForVeg)
# Df SumOfSqs      F Pr(>F)    
# Model      2    7.901 12.483  1e-04 ***
#   Residual 230   72.786   

# Save DBRDA results
# save(dbRDAfungiHabitatResults, file= "RobjectsSaved/dbRDAfungiHabitatResults") #last saved June 4, 2022

#########################################################################
# DOES EVERY EU DIFFER WITH RESPECT TO HABITAT TYPE?
#########################################################################
## Separate out data by each EU, remove non-occurring ASVs, make ordinations, then perform dbRDAs

##########
# EU 52
EU_52_ITSubiq.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_52")
EU_52_ITSubiq.ps <- prune_taxa(taxa_sums(EU_52_ITSubiq.ps) > 0, EU_52_ITSubiq.ps) #remove non-occurring ASVs, now 410 taxa
EU_52_ITSASVs <- ASVs_outta_ps(EU_52_ITSubiq.ps) #get ASV table
EU_52_ITSBC <- vegdist(EU_52_ITSASVs, method= "bray") #get B-C dissimilarities
EU52_metaDat <- as.data.frame(as.matrix(sample_data(EU_52_ITSubiq.ps)))
EU_52_ITSmod <- dbrda(EU_52_ITSBC ~ Habitat, data=EU52_metaDat) #make dbRDA model
set.seed(19)
EU_52_ITSdbRDAHabitatResults <- anova.cca(EU_52_ITSmod, permutations= 9999) #test significance of constraints
# 3.2491  1e-04 ***

##########
# EU 53N
EU_53N_ITSubiq.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_53N")
EU_53N_ITSubiq.ps <- prune_taxa(taxa_sums(EU_53N_ITSubiq.ps) > 0, EU_53N_ITSubiq.ps) #remove non-occurring ASVs, 407 ASVs
EU_53N_ITSASVs <- ASVs_outta_ps(EU_53N_ITSubiq.ps) #get ASV table
EU_53N_ITSBC <- vegdist(EU_53N_ITSASVs, method= "bray") #get B-C dissimilarities
EU_53N_metaDat <- as.data.frame(as.matrix(sample_data(EU_53N_ITSubiq.ps)))
EU_53N_ITSmod <- dbrda(EU_53N_ITSBC ~ Habitat, data=EU_53N_metaDat) #make dbRDA model
set.seed(19)
EU_53N_ITSdbRDAHabitatResults <- anova.cca(EU_53N_ITSmod, permutations= 9999) #test significance of constraints
# 2.176  2e-04 ***

##########
# EU 54S
EU_54S_ITSubiq.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_54S")
EU_54S_ITSubiq.ps <- prune_taxa(taxa_sums(EU_54S_ITSubiq.ps) > 0, EU_54S_ITSubiq.ps) #remove non-occurring ASVs, now 406
EU_54S_ITSASVs <- ASVs_outta_ps(EU_54S_ITSubiq.ps) #get ASV table
EU_54S_ITSBC <- vegdist(EU_54S_ITSASVs, method= "bray") #get B-C dissimilarities
EU_54S_metaDat <- as.data.frame(as.matrix(sample_data(EU_54S_ITSubiq.ps)))
EU_54S_ITSmod <- dbrda(EU_54S_ITSBC ~ Habitat, data=EU_54S_metaDat) #make dbRDA model
set.seed(19)
EU_54S_ITSdbRDAHabitatResults <- anova.cca(EU_54S_ITSmod, permutations= 9999) #test significance of constraints
# 3.9145  1e-04 ***

##########
# EU 8
EU_8_ITSubiq.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_8")
EU_8_ITSubiq.ps <- prune_taxa(taxa_sums(EU_8_ITSubiq.ps) > 0, EU_8_ITSubiq.ps) #remove non-occurring ASVs, now 404
EU_8_ITSASVs <- ASVs_outta_ps(EU_8_ITSubiq.ps) #get ASV table
EU_8_ITSBC <- vegdist(EU_8_ITSASVs, method= "bray") #get B-C dissimilarities
EU_8_metaDat <- as.data.frame(as.matrix(sample_data(EU_8_ITSubiq.ps)))
EU_8_ITSmod <- dbrda(EU_8_ITSBC ~ Habitat, data=EU_8_metaDat) #make dbRDA model
set.seed(19)
EU_8_ITSdbRDAHabitatResults <- anova.cca(EU_8_ITSmod, permutations= 9999) #test significance of constraints
# 13.1  1e-04 ***

##########
# EU 53S
EU_53S_ITSubiq.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_53S")
EU_53S_ITSubiq.ps <- prune_taxa(taxa_sums(EU_53S_ITSubiq.ps) > 0, EU_53S_ITSubiq.ps) #remove non-occurring ASVs, now 392
EU_53S_ITSASVs <- ASVs_outta_ps(EU_53S_ITSubiq.ps) #get ASV table
EU_53S_ITSBC <- vegdist(EU_53S_ITSASVs, method= "bray") #get B-C dissimilarities
EU_53S_metaDat <- as.data.frame(as.matrix(sample_data(EU_53S_ITSubiq.ps)))
EU_53S_ITSmod <- dbrda(EU_53S_ITSBC ~ Habitat, data=EU_53S_metaDat) #make dbRDA model
set.seed(19)
EU_53S_ITSdbRDAHabitatResults <- anova.cca(EU_53S_ITSmod, permutations= 9999) #test significance of constraints
# 3.68  1e-04 ***

##########
# EU 10
EU_10_ITSubiq.ps <- subset_samples(ITS_postUbiquity.ps, EU == "EU_10")
EU_10_ITSubiq.ps <- prune_taxa(taxa_sums(EU_10_ITSubiq.ps) > 0, EU_10_ITSubiq.ps) #remove non-occurring ASVs, now 392
EU_10_ITSASVs <- ASVs_outta_ps(EU_10_ITSubiq.ps) #get ASV table
EU_10_ITSBC <- vegdist(EU_10_ITSASVs, method= "bray") #get B-C dissimilarities
EU_10_metaDat <- as.data.frame(as.matrix(sample_data(EU_10_ITSubiq.ps)))
EU_10_ITSmod <- dbrda(EU_10_ITSBC ~ Habitat, data=EU_10_metaDat) #make dbRDA model
set.seed(19)
EU_10_ITSdbRDAHabitatResults <- anova.cca(EU_10_ITSmod, permutations= 9999) #test significance of constraints
# 7.7019  1e-04 ***