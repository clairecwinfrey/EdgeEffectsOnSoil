# Exploratory Data Analysis
# Aug 17, 2021

# This script is to take a first look at 16S MiSeq data from the samples that I
# took in May 2021 to examine an edge effect across forest and savanna boun-
# daries at the Savannah River Site. Here, I look at sequence counts, rarefy,
# examine top classes and phyla, and look at ordinations. I repeat the ordinations
# and examination of top classes and phyla before and after removing rare species
# (here those that showed up LESS THAN 50 times across whole dataset).
# Then, I move forward with the rarefied, rare-removed data set and do preliminary
# investigations into why outlier samples might be outliers, how savannas/patches
# differ from forests, etc.

# Read in libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")
library("mctoolsr")
library("vegan")

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")
list.files()
seqtab_wTax_mctoolsr <- read.table("seqtab_wTax_mctoolsr.txt")
str(seqtab_wTax_mctoolsr)
head(seqtab_wTax_mctoolsr)
#View(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr)
rownames(seqtab_wTax_mctoolsr)

seqtab <- read.table("seqtab_final.txt", header=T)
dim(seqtab)
#View(seqtab)

tax_final.txt <- read.table("tax_final.txt")
str(tax_final.txt)
#View(tax_final.txt)

#### LOAD FILES IN FOR PHYLOSEQ #####
# 1. OTU table: ASVs/OTUs as rows, samples as columns. This is "seqtab"
otu_mat <- seqtab
seqtab$X #X is the ASV names
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X") #make it so that "X"(which is column with ASV names is the rownames columnn!)
colnames(otu_mat)
rownames(otu_mat)
#View(otu_mat)

# 2. OTU taxonomy
# This is basically the same as the first and last columns of seqtab_wTax_mctoolsr
# So none of the guts that are essentially the ASV table
dim(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr) #samples
head(seqtab_wTax_mctoolsr$V1) #THE ASV names!
head(seqtab_wTax_mctoolsr$V271) #the taxonomy!
step1 <- cbind(seqtab_wTax_mctoolsr$V1, seqtab_wTax_mctoolsr$V271)
head(step1)
colnames(step1) <- step1[1,] #make ASV_ID and taxonomy the column names
head(step1)
dim(step1)
step1 <- step1[-1,] #get rid of first row which is same as column names
str(step1)
#View(step1)
# as.data.frame(step1) #This didn't work so I foreced it below

# Make new columns based on the ";" separator
require(tidyr)
tax_sep <- separate(as.data.frame(step1), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                  "Order", "Family", "Genus", "Species"),
                    sep = ";")

# View(tax_sep)
str(tax_sep)
all.equal(tax_sep$`#ASV_ID`, tax_sep$Species) #Yep
# tax_final.txt shows that no species in this data set... A choice made in bioinformatics script?
colnames(tax_sep)
tax_sep$Species <- NA
#View(tax_sep)

# IDEM for tax sep
tax_mat <- tax_sep %>%
  tibble::column_to_rownames("#ASV_ID")


# 3. Sample metadata
metadata <- read.csv("SRS_AllMetadataAug17.csv") #this new csv has ALL sample metadata
# not just for biological samples (i.e. controls too)
# View(metadata)
colnames(metadata)
rownames(metadata)
metadata$X <- NULL #dunno whatX is but it's not important 
samples_df <- metadata %>%
  tibble::column_to_rownames("MappingSampID")
str(samples_df)
rownames(samples_df)

# Need to rename samples_df to match those in otu_mat; because they don't match
# each other except for the blanks and controls, those are the only ones that got
# phyloseq'd
length(rownames(samples_df)) == length(colnames(otu_mat))
#Will remove x's from the colnames on otu_mat
# Go into atom search and remove Xs
newnames <- c("1","10","100","101",
              "102","103","104","105",
              "106","107","108","109",
              "11","110","111","112",
              "113","114","115","116",
              "117","118","119","12",
              "120","121","122","123",
              "124","125","126","127",
              "128","129","13","130",
              "131","132","133","134",
              "135","136","137","138",
              "139","14","140","141",
              "142","143","144","145",
              "146","147","148","149",
              "15","150","151","152",
              "153","154","155","156",
              "157","158","159","16",
              "160","161","162","163",
              "164","165","166","167",
              "168","169","17","170",
              "171","172","173","174",
              "175","176","177","178",
              "179","18","180","181",
              "182","183","184","185",
              "186","187","188","189",
              "19","190","191","192",
              "193","194","195","196",
              "197","198","199","2",
              "20","200","201","202",
              "203","204","205","206",
              "207","208","209","21",
              "210","211","212","213",
              "214","215","216","217",
              "218","219","22","220",
              "221","222","223","224",
              "225","226","227","228",
              "229","23","230","231",
              "232","233","234","235",
              "236","237","238","239",
              "24","240","241","242",
              "243","25","26","27",
              "28","29","3","30",
              "31","32","33","34",
              "35","36","37","38",
              "39","4","40","41",
              "42","43","44","45",
              "46","47","48","49",
              "5","50","51","52",
              "53","54","55","56",
              "57","58","59","6",
              "60","61","62","63",
              "64","65","66","67",
              "68","69","7","70",
              "71","72","73","74",
              "75","76","77","78",
              "79","8","80","81",
              "82","83","84","85",
              "86","87","88","89",
              "9","90","91","92",
              "93","94","95","96",
              "97","98","99","ExtBlank_1",
              "ExtBlank_2","ExtBlank_3","ExtControlWater_1","ExtControlWater_10",
              "ExtControlWater_11","ExtControlWater_12","ExtControlWater_13","ExtControlWater_14",
              "ExtControlWater_15","ExtControlWater_16","ExtControlWater_17","ExtControlWater_18",
              "ExtControlWater_19","ExtControlWater_2","ExtControlWater_20","ExtControlWater_21",
              "ExtControlWater_22","ExtControlWater_3","ExtControlWater_4","ExtControlWater_5",
              "ExtControlWater_6","ExtControlWater_7","ExtControlWater_8","ExtControlWater_9",
              "PCR_NTC")

colnames(otu_mat) <- newnames
all.equal(sort(colnames(otu_mat)), sort(rownames(samples_df))) #Cool, names are now the same

# Transform otu table and taxonomy tables into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#View(otu_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
SRS_16S.ps <- phyloseq(OTU, TAX, samples)
SRS_16S.ps 
colnames(otu_table(SRS_16S.ps))

sample_names(SRS_16S.ps) # IMPORTANT: THESE ARE THE NAMES BASED ON THE MAPPING
# FILE USED FOR DEMULITPLEXING; I.E. THEY CORRESPOND TO ORDER LOADED INTO ROWS
# OF PLATES (*NOT NOT NOT*) DNA TUBE LABELS 

rank_names(SRS_16S.ps)
sample_variables(SRS_16S.ps)


# FILTER OUT MITOCHONDRIA, CHLOROPLASTS,a and KINGDOM EUKARYOTA
# chloroplast is an order and mitochondria is a family 

# Make a function to get the taxonomy table out of phyloseq:
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# Get taxonomy table out of phyloseq:
SRS_taxTable <- taxtable_outta_ps(SRS_16S.ps)
#View(SRS_taxTable)

# These lines below are technically not needed, since I found a way to do it with phyloseq
# tax_noeuksorNAs is no longer used, but serves as a good check to make sure that phyloseq
# is doing what I want it to. 
# Remove chloroplasts, mitochondria, and all eukaryotes:
tax_noeuksorNAs <- SRS_taxTable %>% 
  filter(Order != "Chloroplast") %>% 
  filter(Family != "Mitochondria") %>% 
  filter(Kingdom != "Eukaryota") %>% 
  filter(Kingdom != "NA") %>% #Remove ASVs where kingdom is unknown ("NA" in first column) 
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in first column) 
# View(tax_noeuksorNAs) 
dim(SRS_taxTable)[1] - dim(tax_noeuksorNAs)[1]
# We filtered out a total of 3255 taxa.

# We now have to trim the OTU table, because it still has the taxa in it that were filtered out
# above. It currently has 38255 taxa, whereas tax_noeuksorNAs has 38255 taxa. We'll want it this 
# way so that we re-make the phyloseq object and then rarefy correctly. Could have done this all in
# same step as above, but oh well!

# Find ASVs where the order is Chloroplast
chloros <- SRS_taxTable %>% 
  filter(Order == "Chloroplast")
dim(chloros) #149 chloroplast ASVs
chloro_names <- rownames(chloros)
# Find ASVs where the family is Mitochondria
mitos <- SRS_taxTable %>% 
  filter(Family == "Mitochondria")
dim(mitos) #959 mitochondria ASVs
mito_names <- rownames(mitos)
# Find ASVs where kingdom is Eukaryota
kingdomEuks <- SRS_taxTable %>% 
  filter(Kingdom == "Eukaryota")
dim(kingdomEuks) #482
euks_names <- rownames(kingdomEuks)
# Find ASVs where kingdom is NAS
kingdomNAs <- SRS_taxTable %>% 
  filter(Kingdom == "NA")
dim(kingdomNAs) #490
kNAs_names <- rownames(kingdomNAs)
# Find ASVs where phylum is NAS
PhylumNAs <- SRS_taxTable %>% 
  filter(Phylum == "NA")
dim(PhylumNAs) #2147
pNAs_names <- rownames(PhylumNAs)

149 + 959 + 482 + 490 + 2147 #4227, this is more than what we removed above because some kingdom that were NAs had NA phyla too
bad_ASVs <- c(chloro_names, mito_names, euks_names, kNAs_names, pNAs_names)
length(bad_ASVs) == 149 + 959 + 482 + 490 + 2147 

# Inspiration for using phyloseq here (Joey711's code: https://github.com/joey711/phyloseq/issues/652
# Remove the ASVs listed above:
all_Taxa <- taxa_names(SRS_16S.ps) #get all tax names in original, uncleaned dataset
ASVstoKeep <- all_Taxa[!(all_Taxa %in% bad_ASVs)]
length(ASVstoKeep) #36,175, matches what we created above!
noeuksorNAs_ps <- prune_taxa(ASVstoKeep, SRS_16S.ps) #new phyloseq object with just the stuff we want!

# TAKING A LOOK AT THE DATA
# How are number of reads distributed across the samples?
otu_table(noeuksorNAs_ps)
seqspersample <- colSums(otu_table(noeuksorNAs_ps))
SeqNumberPlot <-barplot(seqspersample, main="All: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number")

# Just controls:
seqsperctrl <- seqspersample[244:length(seqspersample)]
SeqNumberPlotCtrls <- barplot(seqsperctrl, main="Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.27)
SeqNumberPlotCtrls

# Controls with easier axis names
simple <- c("ExtB1", "ExtB2", "ExtB3", "ExtW1", "ExtW10", "ExtW11", "ExtW12",
"ExtW13", "ExtW14", "ExtW15", "ExtW16", "ExtW17", "ExtW18", "ExtW19",
"ExtW2", "ExtW20", "ExtW21","ExtW22", "ExtW3","ExtW4",
"ExtW5","ExtW6","ExtW7","ExtW8", "ExtW9", "PCRNTC")
seqsperctrl_simple <- setNames(seqsperctrl, simple)

quartz()
barplot(seqsperctrl_simple, main="Controls: Total Sequences Per Sample",
                                    ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.35)

# Let's rarefy:
# What is the minimum for the non-control samples?
min(seqspersample[1:243]) #1224
max(seqspersample[1:243]) #157870
mean(seqspersample[1:243]) # 26730.58

min(seqspersample) #77
max(seqspersample) #157870
mean(seqspersample) # 24645.93
sd(seqspersample)

min(seqspersample[244:length(seqspersample)]) #77
max(seqspersample[244:length(seqspersample)]) #157870
mean(seqspersample[244:length(seqspersample)]) # 3015.692
sd(seqspersample[244:length(seqspersample)])

# Plot this:
seqsPerExSamp <- seqspersample[1:243]
SeqNumberPlotExSamp <- barplot(seqsPerExSamp, main="Soils: Total Sequences Per Sample",
                                    ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.4)

sort(seqsPerExSamp) #take a look in order from least to greatest
sort(seqspersample) #take a look in order from least to greatest (all samples and controls)

# We'll rarefy at 14973, which only costs us 7 samples (plus all but two of the controls)
# The "random rarefaction is made without replacement so that the variance of rarefied
# communities is rather related to rarefaction proportion than to the size of the sample". (quoted bit from 
# vegan's manual page)

set.seed(19)
rarefied.ps <- rarefy_even_depth(noeuksorNAs_ps, sample.size = 14973, replace=FALSE, trimOTUs=TRUE)
# 1121OTUs were removed because they are no longer present in any sample after random subsampling

# Rarefaction curve:
samp.col = c(rep("blue", 243), rep("grey", 26))

rare.plot <- rarecurve(t(otu_table(noeuksorNAs_ps)), step = 3000, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")
rarecurve(t(otu_table(noeuksorNAs_ps)), step = 3000, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")

# How many are left?
colSums(otu_table(rarefied.ps)) #All have 14973
length(colSums(otu_table(rarefied.ps))) #238 samples left (2 less than when I rarefied w/o removing phylum= NA)
# Which samples are these:
sampleNames <- names(colSums(otu_table(rarefied.ps)))
sampleNames

####################################################################
# The following lines through about line 513 is all done BEFORE removing rare taxa)

#####################
# TOP PHYLA
#####################

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
# TURN ASVs INTO PHYLUM LEVEL
rarefied.phylum.glom <-  tax_glom(rarefied.ps, taxrank = "Phylum") 
tax_table(rarefied.phylum.glom) # good, this is only phyla (42 different phyla!)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.phyla.0 <- transform_sample_counts(rarefied.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0)) #weirdly, this is ASV_.... Is this right? 
# I think that the ASVs are just representative from each phylum

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = "EU")
sample_data(relabun.phyla.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
sample_data(relabun.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
sum(relabun.phyla.df[,3]) 
colnames(relabun.phyla.df) 
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
quartz()
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
 theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 0.5% of total abundance")

# "Phyla comprising at least 1% of total abundance"
quartz()
phylumPlot.99percent <- ggplot(data=relabun.phylatop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
phylumPlot.99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 1% of total abundance")

## Below is junk from earlier script... keep for easier manipulation later!
# scale_fill_manual(values = c("#4575b4", "#d73027", "#fc8d59", "#fee090", "#91bfdb", "grey"), 
               #   name= "Phylum", breaks= c("D_1__Firmicutes", "D_1__Proteobacteria", "D_1__Bacteroidetes", "D_1__Actinobacteria", "D_1__Fusobacteria", "< 1% abund."), 
                #  labels =c("Firmicutes", "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Fusobacteria", "< 1% abund.")) +

# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

top_99.5p_phyla <- relabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

#####################
# TOP CLASSES:
#####################

# (adopted from my code at:
   # https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
   # TURN ASVs INTO PHYLUM LEVEL
rarefied.class.glom <-  tax_glom(rarefied.ps, taxrank = "Class") 
tax_table(rarefied.class.glom) # good, this is only class (42 different phyla!)
length(unique(tax_table(rarefied.class.glom))) #854 classes... although some NAs in there too

 # TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.class.0 <- transform_sample_counts(rarefied.class.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.class.0)) # I think that the ASVs are just representative from each class
 
 # MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabun.class.1 <- merge_samples(relabun.class.0, group = "EU")
sample_data(relabun.class.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)
 
 # CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.class.2 <- transform_sample_counts(relabun.class.1, function(x) x / sum(x))
sample_data(relabun.class.2)
 
 # NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabun.class.df <-psmelt(relabun.class.2)
dim(relabun.class.df) #
sum(relabun.class.df[,3]) 
colnames(relabun.class.df) 
relabun.classtop99 <- relabun.class.df
relabun.classtop95 <- relabun.class.df

relabun.classtop99$Class[relabun.classtop99$Abundance < 0.01] <- "< 1% abund."
relabun.classtop95$Class[relabun.classtop95$Abundance < 0.05] <- "< 5% abund."

top_99p_class <- unique(relabun.classtop99$Class)
top_99p_class

top_95p_class <- unique(relabun.classtop95$Class)
top_95p_class

# Phyla comprising at least 5% of total abundance
quartz()
classPlot.95pt <- ggplot(data=relabun.classtop95, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
classPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 5.5))  + ggtitle("Classes comprising at least 5% of total abundance")


# Ignore the lines below... not currently using them
# Bray-Curtis dissimilarities based on square-root transformed data
# Code from mctoolsR
#braydist <- vegan::vegdist(t(otu_table(rarefied.ps)), method = "bray")
#summary(braydist)
# bray.mat <- as.matrix(braydist)
# View(bray.mat) looks as expected
# Plot ordination:

# ord <- calc_ordination(braydist, 'nmds')

# Using a phyloseq tool for convenience!
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ord <- ordinate(rarefied.ps, method = "NMDS", distance = "bray", trymax = 100)

# With no black outline around points
quartz()
rarefiedBrayNMDS <- phyloseq::plot_ordination(rarefied.ps, ord, type= "samples", color= "EU")
rarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# With black outline around points:
quartz()
rarefiedBrayNMDSoutlined <-rarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(aes(fill=EU),color="black",pch=21, size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
rarefiedBrayNMDSoutlined

# With transect included as shape:
# This plot is not very informative!
quartz()
rarefiedBrayNMDStran <- rrarefiedBrayNMDS <- phyloseq::plot_ordination(rarefied.ps, ord, type= "samples", color= "EU", shape= "Transect")
rarefiedBrayNMDStran + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Now, remove biocrust and controls from ordination:
justsoils.ps <- subset_samples(rarefied.ps, Type != "BioCrust" & Type != "ExtContWater")
unique(sample_data(justsoils.ps)[,27]) #only soils
sample_names(justsoils.ps) #another check to show that we have only soils!

# Ordination plot of just soils:
set.seed(19)
ordSoils <- ordinate(justsoils.ps, method = "NMDS", distance = "bray", trymax = 100)
quartz()
soilsrarefiedBrayNMDS <- phyloseq::plot_ordination(justsoils.ps, ordSoils, type= "samples", color= "EU")
soilsrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Add in labels to figure out what weird samples are
quartz()
soilsrarefiedBrayNMDS <- phyloseq::plot_ordination(justsoils.ps, ordSoils, type= "samples", color= "EU", label = "Sample.ID")
soilsrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Visible outliers: (going clockwise from the top left on the ordination plot above):
# 53ND_B_20, 10C_L_10, 53ND_R_40, 53SD_R_10, 52D_R_10, 52D_L_100, (Maybe!) 52D_R_90, 53ND_L_80
outliers <- c("53ND_B_20", " 10C_L_10", "53ND_R_40", "53SD_R_10", "52D_L_100", "52D_R_90", "53ND_L_80")

###########################################################################
# REMOVE "RARE" TAXA AND THEN RE-RUN 1) sequences per sample, 
# 2) top phyla, 3) top classes, 4) ordinations
# Questions when comparing pre and post removal of rare taxa:
# 1) Do top phyla or classes change, or outliers, after this removal? If so, this is
# evidence of these shifts being driven by rare taxa:

# Use function to get ASV table out of phyloseq so that we can view it better
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

rare_ASVtab <- ASVs_outta_ps(rarefied.ps) 
# Add column for abundance of ASVs across all (pre-rarefied) samples
rare_ASVtab$Abundance <- rowSums(rare_ASVtab)

# Most are rare! 
plot(rare_ASVtab$Abundance)
length(which(rare_ASVtab$Abundance <= 50)) #23,775 ASVs have 50 or fewer occurrences across the rarefied data set 
length(which(rare_ASVtab$Abundance <= 35)) #20,744 ASVs have 30 or fewer occurrences across the rarefied data set
length(which(rare_ASVtab$Abundance <= 10)) #9,818 ASVs have 10 or fewer occurrences across the rarefied data set

length(which(rare_ASVtab$Abundance >= 50)) #10,257 ASVs have 50 or more occurrences across the rarefied data set 

# We'll get rid of the ASVs that do not occur AT LEAST 50 times across our dataset:
keptASVsindex <- which(rare_ASVtab$Abundance >= 50) #gives row numbers to keep in the dataset 
max(keptASVsindex) # last row to keep is 15719 (but not all before are kept probably b/c of rarfying)
rare_ASVtab$Abundance[15719] #has exactly fifty 

ASVsTrimmed <- rare_ASVtab[keptASVsindex,] #keep only rows of the index, and all columns (i.e. all samples)
dim(ASVsTrimmed) #10257 ASVs across 239 samples, as expected 
plot(ASVsTrimmed$Abundance) # tail is still long, because most ASVs are still rare
# relative to the most abundant ASVs

length(which(ASVsTrimmed$Abundance == 50)) # 153 ASVs were right at cut off!
range(ASVsTrimmed$Abundance) #50 75046
length(which(ASVsTrimmed$Abundance >= 1000)) #only 517 ASVs appear more than 1000 times
# across the data set, which is in average of about 4 times per sample!

# How many ASVs occur in every sample?
# Said another way: how many rows (i.e. ASVs) occur >= 1 time across all columns (i.e. samples)
# Or: are there any rows that contain no zeros?
# Inspiration from this post: https://stackoverflow.com/questions/9977686/how-to-remove-rows-with-any-zero-value

# Let's you find rows that contain zeros and then filter them out
rownames(ASVsTrimmed[!(apply(ASVsTrimmed, 1, function(y) any(y == 0))),])
# "ASV_1"  "ASV_3"  "ASV_4"  "ASV_12" "ASV_17" are found in every sample 
# This corresponds to rows:
ASVsinAllindex <- c(1, 3, 4, 10, 15)
rownames(ASVsTrimmed)[ASVsinAllindex] # yes, this matches above
# What are these ASVs?
rarefied_taxa[ASVsinAllindex,]

