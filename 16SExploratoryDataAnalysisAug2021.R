# Exploratory Data Analysis
# Aug 17, 2021

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyr")

setwd("/Users/clairewinfrey/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")
list.files()
seqtab_wTax_mctoolsr <- read.table("seqtab_wTax_mctoolsr.txt", header=T)
str(seqtab_wTax_mctoolsr)
head(seqtab_wTax_mctoolsr)
View(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr)
rownames(seqtab_wTax_mctoolsr)

seqtab <- read.table("seqtab_final.txt", header=T)
dim(seqtab)
#View(seqtab)

tax_final.txt <- read.table("tax_final.txt")
str(tax_final.txt)
#View(tax_final.txt)


library("mctoolsr")
### USING A HYBRID APPROACH OF http://leffj.github.io/mctools
# Going off of Matt's suggestion-- first column ASV ID, then headers sample ID and last column taxonomy all together



#### LOAD FILES IN FOR PHYLOSEQ #####
# 1. OTU table: ASVs/OTUs as rows, samples as columns. This is "seqtab"
otu_mat <- seqtab
seqtab$X
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X") #make it so that "X"(which is column with ASV names is the rownames columnn!)
colnames(otu_mat)
#View(otu_mat)

# 2. OTU taxonomy
# This is basically the same as the first and last columns of seqtab_wTax_mctoolsr
# So none of the guts that are essentially the ASV table
dim(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr)
head(seqtab_wTax_mctoolsr$V1) #THE ASV names!
head(seqtab_wTax_mctoolsr$V271) #the taxonomy!
step1 <- cbind(seqtab_wTax_mctoolsr$V1, seqtab_wTax_mctoolsr$V271)
str(step1)
colnames(step1) <- step1[1,]
step1 <- step1[-1,]
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
View(tax_sep)

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
View(SRS_taxTable)


# Remove chloroplasts, mitochondria, and all eukaryotes:
SRS_taxTable_noeuks <- SRS_taxTable %>% 
  filter(Order != "Chloroplast") %>% 
  filter(Family != "Mitochondria") %>% 
  filter(Kingdom != "Eukaryota")
# View(SRS_taxTable_noeuks) 36,665 rows is exactly what we expect after the removal of all of these taxa!




### DID NOT DECIDE TO DROP ANY SINGLE OR DOUBLETONS FOR NOW 
# Filter out blanks and such for now
SRS_16S_soilonly <- SRS_16S.ps %>%
  subset_samples(Type == "soil")

SRS_16S_soilonly 
