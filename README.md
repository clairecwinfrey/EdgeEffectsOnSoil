# EdgeEffectsOnSoil
This repository contains the scripts used to process 16S and ITS data from soils collected at the Savannah River Site (SRS) in May 2021. Scripts were written by Claire Winfrey, with substantial input on analyses from Julian Resasco and Noah Fierer. Below is a description of each script, listed in the order that scripts should be run to replicate results. Analyses on 16S and ITS data were done in separate scripts, in parallel.

**16S**
1. SRS_16S_Bioinformatics.R- This script demultiplexes raw 16S rDNA reads from the MiSeq, uses dada2 to trim, denoise, dereplicate, and merge reads, and assigns bacterial/archaeal taxonomy (using SILVA version 132).
2. 16SExploratoryDataAnalysisAug2021.R - This script takes a first look at 16S MiSeq data. I look at sequence counts, rarefy, examine top classes and phyla, and look at ordinations bsaed on different EUs, different plant communities, etc. I repeat the ordinations and examination of top classes and phyla before and after removing rare species. I also investigate samples that seem to be outliers in ordination space. Finally, this script creates the cleaned up dataset used in all downstream scripts and analyses.
3. UbiquityMedianSetup.R- This script applies a ubiquity filter to the "trimmed" dataset made in 16SExploratoryDataAnalysisAug2021.R.
4. EdgeEdgebyASV_allSites.R - This script performs a differential abundance analysis using DESeq2 package to find patch and forest specialist ASVs for ONLY bacterial and archaeal data. Importantly, this script uses the post ubiquity dataset created in UbiquityMedianSetup.R, but with 1 count added to each of the ASV abundance counts (necessary for the ITS DESeq analysis). Then, it transforms the abundance of each of these ASVs into a z-score (where the mean and standard deviation for calculation of that ASV's z-score are drawn from within a given EU but at each point along the transect). Then, using these z-scores, the script finds the logistic fit of the curve of each ASV as it goes over the transect. Finally, this script makes a stacked barplot showing the number of ASVs in each phyla of the postUbiquity dataset, and whether it is an indicator species for the patch, the matrix, or neither.

**ITS**
1. SRS_ITS_Bioinformatics.R- This script demultiplexes raw ITS rDNA reads from the MiSeq, uses dada2 to trim, denoise, dereplicate, and merge reads, and assigns bacterial/archaeal taxonomy (using UNITE database, Oct. 5, 2021 release).
2. ITS_ExploratoryDataAnalysis.R - This script explores the ITS sequencing data, roughly doing the same thing for ITS data as 16SExploratoryDataAnalysisAug2021.R. However, unlike in the 16S script, no differential abundance analysis is done at this stage to explore patch versus matrix.
3. ITS_UbiquityMedianSetup.R- this script applies a ubiquity filter to the "trimmed" dataset made in ITS_ExploratoryDataAnalysis.R.
4. EdgeEffectsAllSitesFUNGI.R-- see description above for EdgeEdgebyASV_allSites.R, but this script uses the ITS data.
5. PostUbiqGraphics.R - makes NMDS ordinations and Bray-Curtis dissimilarities along the transect for dataset created in ITS_UbiquityMedianSetup.R (i.e. after removing rare ASVs and not very ubiquitous ones)


**Scripts with additional analyses that probably won't make it into paper**
