READ ME file for Waraniak JM, Batchelor S, Wagner T, and Keagy J. Landscape transcriptomics analysis detects thermal stress responses and potential adaptive variation in wild brook trout (Salvelinus fontinalis) during successive heatwaves. 

This READ ME file describes the datasets and R code necessary to run all analyses described in the above manuscript.

- - - - - - - - - - - - - - - - - 
Note: Download all files into the same directory and set this directory as your working directory when running R code.

01_TempData.R
Format: R Script
R code to process raw data from HOBO temperature loggers and field records of sampled brook trout to create files that are used as inputs by other R scripts.
Required input files: 2022_heatwave_temps.xlsx, 2022_trout_sampling_data.xlsx
Output files used in other scripts: BrookTrout_Heatwaves_DailyTemps.txt, trout_enviro_data.txt

02_ReadCountProcessing.R
Format: R Script
R code to process transcript count output files from SALMON bioinformatics pipeline. Filters rare transcripts, converts raw counts to log-transformed counts per million, and creates files to be used as inputs by other R scripts.
Required input files: SALMON quant.sf files stored in folder named "trimmed" (available from NCBI GEO Project Accession No. GSE273935), trout_enviro_data.txt, BrookTrout_Sexes.txt, BrookTrout_LengthWeight.txt
Output files used in other scripts: salmon_logcounts_environ.txt

03_BrookTroutHeatwave_Map.R
Format: R Script
R code to create Figure 1 as seen in Manuscript.
Required input files: NHDFlowline.shp from National Hydrology Database (mid-Atlantic Region), BrookTrout_Heatwaves_DailyTemps.txt, 2022_trout_sampling_data.xlsx

04-1_Salmon_WGCNA_NetworkConstruction_Auto.R
Format: R Script
R code to run the network construction segment of the WGCNA analysis non-interactively. This procedure requires a large amount of memory (~200GB) so it is recommended to run on a cluster. This R script can be submitted as a batch job for a Slurm or PBS job manager.
Required input files: salmon_logcounts_environ.txt
Output files used in other scripts: Salmon_WGCNA_netwrkConstruction_auto.RData, salmonTOM-block.1.RData

04_WGCNA_Analyses.R
Format: R Script
R code to run the weighted gene correlations network analysis and generalized additive mixed models of WGCNA modules. Creates Figures 2 & 3 in manuscript.
Required input files: salmon_logcounts_environ.txt, BrookTroutHeatwave_netwrkConstruction_auto.RData, salmonTOM-block.1.RData, Phylofish_GOData.txt,  BrookTrout_Heatwaves_DailyTemps.txt, ME_gam_heatmap.txt
Output files used in other scripts: WGCNA.RData

05_DE_analysis.R
Format: R Script
R code to run the ImpulseDE2 analysis to identify differentially expressed genes and generalized additive mixed models on DE genes to see how they relate to environmental covariates. Creates Figure 5 in the manuscript. 
Required input files: salmon_logcounts_environ.txt, SALMON quant.sf files stored in folder named "trimmed" (available from NCBI GEO Project Accession No. GSE273935), WGCNA.RData

06_Trajectory_Analysis.R
Format: R Script
R code to run the trajectory analysis and produce Figure 4 in the manuscript. 
Required input files: salmon_logcounts_environ.txt, BrookTrout_Heatwaves_DailyTemps.txt, Phylofish_GOData.txt

07_RADseq_analysis.R
Format: R Script
R code to process SNP data output from Stacks, run PCA on SNP data, prepare input files and process output files for fineRADstructure, and run redundancy analysis variance partitioning. 
Required input files: radiator_data_20240506@0951.rad, populations.haps_chunks.out, populations.haps_chunks.mcmc.xml, populations.haps_chunks.mcmcTree.xml, salmon_logcounts_environ.txt

2022_trout_sampling_data.xlsx
Format: Microsoft Excel Workbook
Sheet: README - contains detailed descriptions of all columns in other sheets
Sheet: HW site and sampling event Info - contains time and location information for each sampling event
	* Column A: stream name
	* Column B: date of the sampling event (MM/DD/YYYY)
	* Column C: number of brook trout caught during sampling
	* Column D: number of brown trout caught during sampling
	* Column E: starting time of sampling event
	* Column F: completion time of sampling event
	* Column G: in-stream temperature as measured by handheld thermometer at start of sampling event (ºC)
	* Column H: latitude and longitude values in decimal degrees, separated by a comma
Sheet: HW trout field data - field data collected for each individual sampled during sampling events
	* Column A: date of sampling event (MM/DD/YYYY)
	* Column B: stream name
	* Column C: unique ID code for each sampled individual
	* Column D: species of the sampled individual (Brook or Brown trout)
	* Column E: time when electrofishing run that captured individual started
	* Column F: time when electrofishing run that captured individual ended
	* Column G: time when processing of individual started (weight/length was measured, tissue samples were collected)
	* Column H: total length of individual (mm)
	* Column I: mass of individual (g)
	* Column J: In-stream temperature as measured by handheld thermometer at start of electrofishing run that captured individual (ºC)
	* Column K: additional comments
	* Column L: whether or not a gill tissue sample was taken from a captured individual (Y/N)
	* Column M: whether or not a fin tissue sample was taken from a captured individual (Y/N)
Sheet: RNA QC - Quality control measurements to quantify amount and integrity of sequenceable RNA from each gill tissue sample
	* Column A: unique ID code for each sampled individual
	* Column B: date of RNA/DNA extraction (MM/DD/YYYY)
	* Column C: first Qubit estimate of RNA concentration (ng/µl)
	* Column D: second Qubit estimate of RNA concentration (ng/µl)
	* Column E: Nanodrop estimate of RNA concentration (ng/µl)
	* Column F: Ratio of 260 to 280 nm wavelengths from Nanodrop
	* Column G: Ratio of 260 to 230 nm wavelength from Nanodrop
	* Column H: date of first test using sex ID markers (MM/DD/YYYY)
	* Column I: result of first sex ID marker test (M/F)
	* Column J: date of follow-up test using sex ID markers (MM/DD/YYYY)
	* Column K: result of follow-up sex ID marker test (M/F) 
	* Column L: date of TapeStation test to determine mRNA quality and concentration (MM/DD/YYYY)
	* Column M: TapeStation RINe score showing quality of mRNA in sample (closer to 10 is better)
	* Column N: TapeStation estimate of mRNA concentration in sample (ng/µl)

2022heatwave_temps.xlsx
Format: Microsoft Excel Workbook
Sheet: README - contains detailed descriptions of all columns in other sheets
Sheet: logger info - contains information on each of the HOBO loggers
	* Column A: logger ID
	* Column B: stream name
	* Column C: whether the temperature is measued in water or air
	* Column D: latitude/longitude coordinates of where the logger was placed (decimal degrees)
	* Column E: frequency of temperature recording by each logger (minutes)
	* Column F: date the logger was set up in stream (MM/DD/YYYY)
	* Column G: notes on where the logger was placed
	* Column H: additional comments
Sheet: All temp data - data collected by all temperature loggers 
	* Column A: stream name
	* Column B: logger ID
	* Column C: timestamp of when temperature record was taken (MM/DD/YYYY HH:MM)
	* Column D: temperature (ºC)

BrookTrout_LengthWeight.txt
Format: Tab-Delimited Text File
Contains records for all 215 individuals captured during field sampling events and a header row for column names. The first column is individual ID, the second column is the recorded total length of that individual in millimeters, and the third row is the recorded mass of that individual in grams.

BrookTrout_Sexes.txt
Format: Tab-Delimited Text File
Contains records from sex identification tests for 119 individuals from which gill tissue samples were collected and DNA/RNA was extracted and a header row for column names. The first column is individual ID and the second column is the determined sex (M for male or F for female).

ME_gam_heatmap.txt
Format: Tab-Delimited Text File
Contains results from the generalized additive mixed models conducted on 33 WGCNA modules (in WGCNA_Analyses.R) and a header row for column names. The first column is the module color name. All other columns are the p-value from the GAMM determining the statistical significance of the explanatory variables tested. P-values > 0.05 have been changed to 1. 

Phylofish_GOData.txt
Format: Tab-Delimited Text File
Simplified version of data downloaded from BioMart search for Gene Ontology terms associated with transcripts from the brook trout transcriptome available on Phylofish (as of 01/22/2024). This data file has already been filtered to only include transcripts that passed the abundance filters in our dataset. The first column is a seven digit GO ID code, the second column is the associated GO term, the third column is the ontology category (B - biological process, C - cellular component, M - molecular function), and the fourth column is the associated brook trout transcript name. Note that GO IDs can be repeated, but each entry is associated with a unique transcript.

radiator_data_20240506@0951.rad
Format: Genomic Data Structure
Filtered SNP genotypes for 47 individuals genotyped by RAD-seq. 

