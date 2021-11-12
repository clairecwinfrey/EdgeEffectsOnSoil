
# Pull out individual EUs from total dataset
EU_52_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_52")
EU_53N_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53N")
EU_54S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_54S")
EU_8_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_8")
EU_53S_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_53S")
EU_10_Soils.ps <- subset_samples(trimmedJustsoils.ps, EU == "EU_10")
listByEU <- list(EU_52_Soils.ps, EU_53N_Soils.ps, EU_54S_Soils.ps, 
                 EU_8_Soils.ps, EU_53S_Soils.ps, EU_10_Soils.ps)

# TransectDissim obtains the Bray-Curtis dissimilarities between 1) 10 m and 
# all points on transect, 2) 100 m and all points on transect, and 3) 50 m
# (edge) and all points on transect. It specifically divides the EU in the 
# provided phyloseq object into its four component transects, calculates Bray
# Curtis dissimilarities, and then pulls out the comparisons mentioned above.
# It returns a list of 12 matrices (10 m comparisons, 100 m comparisons, and
# 50 m (edge) comparisons).

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

comps52 <- transectDissim(EU_52_Soils.ps)

allEUsComps<- vector("list", 4)
for (i in 1:length(listByEU)){
  allEUsComps[[i]] <- transectDissim(listByEU[[i]])
}
allEUsComps
