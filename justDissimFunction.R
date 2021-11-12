# justDissimFunction.R
# In this script, I show that the function "transectDissim" is calculating
# values as expected. Here, I explicitly compare the results of the output of the
# "transectDisim" function (see DissimPatterns2.R) with a manual calculation
# of some of these comparison values (full manual calculations are in
# DissimPatterns.R which is no longer a part of the workflow for this
# project, as its code is replaced with that in DissimPatterns2.R).

# The "transectDissim" function here matches that found in DissimPatterns2.R"

# Libraries
library("phyloseq")
library("ggplot2")      
library("dplyr")        
library("tibble")       
library("tidyr")
library("vegan")

setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch/Bioinformatics")

load("trimmedJustsoils.ps") #load phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset (from 
# 16SExploratoryDataAnalysisAug2021).

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

# 3. metadata_outta_ps takes the phyloseq metadata table and converts it to a dataframe
metadata_outta_ps <- function(physeq){ #input is a phyloseq object
  metaDat <- sample_data(physeq)
  return(as.data.frame(metaDat))
}

# 4. TransectDissim obtains the Bray-Curtis dissimilarities between 1) 10 m and 
# all points on transect, 2) 100 m and all points on transect, and 3) 50 m
# (edge) and all points on transect. It specifically divides the EU in the 
# provided phyloseq object into its four component transects, calculates Bray
# Curtis dissimilarities, and then pulls out the comparisons mentioned above.
# It returns a list of 12 matrices (10 m comparisons, 100 m comparisons, and
# 50 m (edge) comparisons). #Note: for comparisons with calculating manually, i.e.
# proof that this function works, see "justDissimFunction.R" and DissimPatterns.R.
# (The latter does not directly compare, but manually calculates the Bray Curtis
# dissimilarities.)

transectDissim <- function(physeq){
  ASVsdf <- ASVs_outta_ps(physeq) #get ASV table using functioned defined earlier; samples are rows, ASVs = columns
  metaDf <- metadata_outta_ps(physeq) #get metadata using function defined earlier; samples are rows
  # Next three lines make a new df that has variables of interest and ASVs
  ASVsdf$Transect <- metaDf$Transect
  ASVsdf$Meter <- metaDf$Meter
  # Sample names (numbered) are rows. Columns are ASV names, except for last two which are Transect and Meter. 
  # Make a new column that combines info from meter and EU
  sampNames <- rep(NA, nrow(ASVsdf)) #pre-allocate vector for new names, length is 38 b/c 38 samples from EU 52
  for (i in 1:nrow(ASVsdf)){
    sampNames[i] <- paste(c(ASVsdf[i,(ncol(ASVsdf)-1)], ASVsdf[i,ncol(ASVsdf)]), collapse="_") #transect and meter
  } 
  ASVsdf$sampNames <- sampNames
  
  # The piped code below creates a list of 4 tibbles, corresponding to each transect
  transectLists <- ASVsdf %>%  #make multiple datasets based on transect
    dplyr::group_by(Transect) %>%  #group by EU
    dplyr::group_split()

  # Remove Meter and transect columns then change sampNames into the rownames 
  transectLists2 <- vector("list", 4) # preallocate list of vectors transformed transectLists, can drop once code is working
  # and just write over transectLists
  for (i in 1:length(transectLists)){
    transectLists2[[i]] <- transectLists[[i]] %>% 
      arrange(Meter) %>% #sort in ascending order
      select(-c(Meter, Transect)) %>% #drop meter and transect columns b/c info in sampNames
      column_to_rownames("sampNames") #make sampNames column in the rownames 
  }
  
  # The lines below get Bray-Curtis dissimilarities within transects (as created above)
  BCs <- vector("list", 4) #preallocate vector to store Bray-Curtis dissimilarities
  for (i in 1:length(transectLists2)){ #there are 4 tibbles within this, corresponding to each transect. 
    # they are in alphabetical order, starting with transect B
    # make a list of lists call Bray-Curtis with 4 subset lists
    # Get Bray-Curtis dissimilarities for a subset (not last two rows that are transect and meter)
    BCs[[i]] <- as.matrix(vegdist(transectLists2[[i]][,1:(ncol(transectLists2[[i]])-2)], method="bray"))
    diag(BCs[[i]]) <- NA
  }
  
  # The lines below pull out Bray Curtis distances within each transect for 10 m
  # and each other point along the transect, and 100 m and each point along transect
  ### AT SOME POINT, CAN ADD B/W EDGE AND ALL POINTS AND POOLED FOREST AND POOLED
  # FOREST AND ALL POINTS INTO CODE RIGHT HERE! WOULD NEED TO CHANGE PREALLOCATED LIST
  comps10m <- vector("list", 4) #preallocate space for all comps with 10 m for each transect
  for (i in 1:length(comps10m)){
    comps10m[[i]] <- BCs[[i]][1,]
  }
  comps100m <- vector("list", 4) #all comps with 10 m for each transect
  for (i in 1:length(comps100m)){
    comps100m[[i]] <- BCs[[i]][nrow(BCs[[i]]),]
  }
  comps50m <- vector("list", 4) #preallocate space for 50m (edge) comparisons
  index50m <- which(rownames(BCs[[1]])=="B_50") #get row number for edge
  # I use this above so that function can be flexible if samples are missing 
  for (i in 1:length(comps50m)){
    comps50m[[i]] <- BCs[[i]][index50m,]
  }
  #put it all together in list of vectors
  allComps <- c(comps10m, comps100m, comps50m)
  return(allComps)
}

# 4. phylosq_sep_variable is needed to run the code that manually calculates Bray-Curtis dissimilarities
# from: https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_sep_variable.R
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
  # require(phyloseq)
  # require(plyr)
  
  ## Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }
  
  ## Extract sample meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  
  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  
  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }
  
  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  
  ## Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)
  
  ## Extract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })
  
  ## Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })
  
  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }
  
  return(res)
}

###############################################################################
# USING transectDissim function 
###############################################################################

# To get dissimilarities, first split phyloseq object by EUs
EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
EU_53N_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53N")
EU_54S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_54S")
EU_8_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_8")
EU_53S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53S")
EU_10_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_10")

# Make a list with all of the EUs phyloseq objects 
listByEU <- list(EU_52_Soils.ps, EU_53N_Soils.ps, EU_54S_Soils.ps, 
                 EU_8_Soils.ps, EU_53S_Soils.ps, EU_10_Soils.ps)
# Set their names!
names(listByEU) <- c("EU_52", "EU_53N", "EU_54S", "EU_8", "EU_53S", "EU_10")
names(listByEU)

# Using function defined above, loop over all EUs to get all dissimilarity comparisons
allEUsComps <- vector("list", 6) #create an empty list for 6 elements, one per EU
for (i in 1:length(listByEU)){
  allEUsComps[[i]] <- transectDissim(listByEU[[i]])
}
names(allEUsComps) <- names(listByEU) #add names for each EU
names(allEUsComps) #correct!
length(allEUsComps$EU_52) #lengths look correct too

###############################################################################
# MANUAL CALCULATION & COMPARISON WITH transectDissim OUTPUT
###############################################################################

# The lines below pull out Bray Curtis distances within each transect for 10 m
# and each other point along the transect, and 100 m and each point along transect,
# MANUALLY, to double check results of transectDissim.
# I use EU 52 transect T as an example:

# Separate out samples by EU and then by transect
# Then get B-C dissimilarities across these transects 
#### EU 52 #####
unique(sample_data(EU_52_Soils.ps)$EU) #all are EU 52, as expected!
EU_52_transSplit <- phyloseq_sep_variable(EU_52_Soils.ps, variable= "Transect", drop_zeroes = T)
### T ###
EU_52_T <- EU_52_transSplit$T
# Get just ASV table
ASVtab_52T <- ASVs_outta_ps(EU_52_T)
# View(ASVtab_52T) #samples are rows, ASVs are columns

# Get Bray-Curtis dissimilarities
BrayDist_52T <- vegdist(ASVtab_52T, method="bray")
BrayDist_52T.mat <- as.matrix(BrayDist_52T)
diag(BrayDist_52T.mat) <- NA

# Make data frame for indexing matrix
colnames(BrayDist_52T.mat) == rownames(sample_data(EU_52_T))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps52T <- colnames(BrayDist_52T.mat) #this is the Mapping Sample ID of these
meter52T <- sample_data(EU_52_T)$Meter
index52T.df <- data.frame(samps52T, meter52T) #columns line up and match up with metadata

# Get row numbers for each meter. Because these the rows (and columns) in the BC matrix are the same as the rows
# in the index dataframe, this is how we'll get the rows and columns to use to subset the BC matrix
m10_52T <- which(index52T.df$meter52T == 10)
m20_52T <- which(index52T.df$meter52T == 20)
m30_52T <- which(index52T.df$meter52T == 30)
m40_52T <- which(index52T.df$meter52T == 40)
m50_52T <- which(index52T.df$meter52T == 50)
m60_52T <- which(index52T.df$meter52T == 60)
m70_52T <- which(index52T.df$meter52T == 70)
m80_52T <- which(index52T.df$meter52T == 80)
m90_52T <- which(index52T.df$meter52T == 90)
m100_52T <- which(index52T.df$meter52T == 100)

# B-C distances between patch at 10m and each point along transect
m10_52T_comp10 <- BrayDist_52T.mat[m10_52T, m10_52T] #should be NA
m10_52T_comp20 <- BrayDist_52T.mat[m10_52T, m20_52T]
m10_52T_comp30 <- BrayDist_52T.mat[m10_52T, m30_52T]
m10_52T_comp40 <- BrayDist_52T.mat[m10_52T, m40_52T]
m10_52T_comp50 <- BrayDist_52T.mat[m10_52T, m50_52T]
m10_52T_comp60 <- BrayDist_52T.mat[m10_52T, m60_52T]
m10_52T_comp70 <- BrayDist_52T.mat[m10_52T, m70_52T]
m10_52T_comp80 <- BrayDist_52T.mat[m10_52T, m80_52T]
m10_52T_comp90 <- BrayDist_52T.mat[m10_52T, m90_52T]
m10_52T_comp100 <- BrayDist_52T.mat[m10_52T, m100_52T]

EU52T_10mcomp_manual <- c(m10_52T_comp10, m10_52T_comp20, m10_52T_comp30, m10_52T_comp40, 
                   m10_52T_comp50, m10_52T_comp60, m10_52T_comp70, m10_52T_comp80,
                   m10_52T_comp90, m10_52T_comp100)

# B-C distances between forest at 100m and each point along transect
m100_52T_comp100 <- BrayDist_52T.mat[m100_52T, m100_52T] #add 100 and 100
m100_52T_comp90 <- BrayDist_52T.mat[m100_52T, m90_52T]
m100_52T_comp80 <- BrayDist_52T.mat[m100_52T, m80_52T]
m100_52T_comp70 <- BrayDist_52T.mat[m100_52T, m70_52T]
m100_52T_comp60 <- BrayDist_52T.mat[m100_52T, m60_52T]
m100_52T_comp50 <- BrayDist_52T.mat[m100_52T, m50_52T]
m100_52T_comp40 <- BrayDist_52T.mat[m100_52T, m40_52T]
m100_52T_comp30 <- BrayDist_52T.mat[m100_52T, m30_52T]
m100_52T_comp20 <- BrayDist_52T.mat[m100_52T, m20_52T]
m100_52T_comp10 <- BrayDist_52T.mat[m100_52T, m10_52T]

EU52T_100mcomp_manual <- c(m100_52T_comp10, m100_52T_comp20, m100_52T_comp30,
                  m100_52T_comp40, m100_52T_comp50, m100_52T_comp60, m100_52T_comp70,
                  m100_52T_comp80, m100_52T_comp90, m100_52T_comp100)

# Are the manual values the same as my function calculates?
EU_52T_10m_function <- unlist(allEUsComps$EU_52[4]) #10 m comparisons
names(EU_52T_10m_function) <- NULL #remove names so we can compare
EU_52T_100m_function <- unlist(allEUsComps$EU_52[8]) #100 m comparisons
names(EU_52T_100m_function) <- NULL

print(paste(EU_52T_10m_function, EU52T_10mcomp_manual)) #yes, these are the same!
print(paste(EU_52T_100m_function, EU52T_100mcomp_manual)) #yes, these are the same!