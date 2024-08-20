#read in libraries
library(WGCNA)
library(data.table)
salmon_data <- fread("PATH_TO_YOUR_DIRECTORY/salmon_logcounts_environ.txt", drop=1) #read in dataset with environmental and expression data
enviro_data <- salmon_data[,1:16] #create a separate data frame with just environmental data
express_data <- salmon_data[,-(1:16)] #create a separate data frame with just expression data
rownames(express_data) <- salmon_data$Sample_ID #name the expression data rows with sample IDs so future analyses of the expression data can be aligned to environmental data

enableWGCNAThreads() #allow multiple threads for faster computing
net <- blockwiseModules(express_data, maxBlockSize = 35000, power = 4, TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, pamRespectsDendro = F, saveTOMs = T, saveTOMFileBase = "salmonTOM", verbose=5, nthreads=8) #create the TOM and save for loading into WGCNA_Analyses.R
moduleLabels <- net$colors #save an object with the color names as the labels for the modules
moduleColors <- labels2colors(net$colors) #create a copy of the color names reformatted so they can print colors
MEs <- net$MEs #save module eigengene values to a dataframe
geneTree <- net$dendrograms[[1]] #save the dendrogram of transcripts
save(MEs, moduleLabels, moduleColors, geneTree, file="PATH_TO_YOUR_DIRECTORY/BrookTroutHeatwave_netwrkConstruction_auto.RData") #save these results to an .RData file 
