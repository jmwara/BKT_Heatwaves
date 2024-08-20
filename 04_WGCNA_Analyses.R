#read in the libraries
library(data.table)
library(MuMIn)
library(WGCNA)
library(dplyr)
library(corrplot)
library(mgcv)
#BiocManager::install("GOfuncR")
library(GOfuncR)
library(ggplot2)
library(gridExtra)
#BiocManager::install("ggtree")
library(ggtree)
library(ggpubr)
library(tidyr)

salmon_data <- fread("salmon_logcounts_environ.txt", drop=1) #read in dataset with environmental and expression data
enviro_data <- salmon_data[,1:16] #create a separate data frame with just environmental data
express_data <- salmon_data[,-(1:16)] #create a separate data frame with just expression data
rownames(express_data) <- salmon_data$Sample_ID #name the expression data rows with sample IDs so future analyses of the expression data can be aligned to environmental data

sample_dendro <- hclust(dist(express_data), method="average") #create a dendrogram based on the euclidean distance among the samples based on the expression data
sample_dendro$labels <- rownames(express_data) #add the sample names to the dendrogram 
enviro_data_norm <- stdize(enviro_data[,4:12]) #center & standardize the quantitative environmental data so it is easier to see patterns in variation 
traitColors <- numbers2colors(enviro_data_norm) #assign colors to the standardized environmental data, the low values are blue, average values are white, and high values are red
plotDendroAndColors(sample_dendro, traitColors, groupLabels = names(enviro_data[,4:12])) #plot the dendogram with a heat map of the environmental values underneath the dendrogram

disableWGCNAThreads() #turn off multi-threaded WGCNA analyses (just for the pickSoftThreshold, for some reason it was running into problems on the laptop)
powers <- c(c(1:10), seq(from=12, to=20, by=2)) #create a vector of powers to test (1-10 and then 12-20 evens only)
sft_signed <- pickSoftThreshold(express_data, powerVector = powers, verbose = 5, networkType = "signed") #run the soft threshold selector using the power values specified

plot(sft_signed$fitIndices$Power, -sign(sft_signed$fitIndices$slope)*sft_signed$fitIndices$SFT.R.sq, xlab="Soft Threshold (Power)", ylab="Scale Free Topology Model Fit, signed R^2", main=paste("Scale Independence")) #plot the power vs. the slope*r squared 
plot(sft_signed$fitIndices$Power, sft_signed$fitIndices$mean.k., xlab="Soft Threshold (Power)", ylab="Mean Connectivity", type="n", main=paste("Mean Connectivity")) #plot the 
text(sft_signed$fitIndices$Power, sft_signed$fitIndices$mean.k., labels = powers, cex=0.9, col="red")

load(file="BrookTroutHeatwave_netwrkConstruction_auto.RData") #read in the results of the one block network construction run on the cluster
table(moduleColors) #print the module name and the numbers of genes in each module
plotDendroAndColors(geneTree, moduleColors, dendroLabels = F, hang=0.03, addGuide = T, guideHang = 0.05) #plot the gene expression dendrogram and the module colors from the one block network construction

enableWGCNAThreads() #allow the R session to use multiple threads to construct WGCNA networks
load("salmonTOM-block.1.RData") #load in the TOM from the cluster session (too computationally intensive for a laptop to do in one block)
salmon_signed_TOM <- TOM #rename the TOM
#Test the effect of cut height on the modules detected
salmon_signed_network_cH0.995 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with detect cut height of 0.995
salmon_signed_network_cH0.996 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.996) #run the module detection with detect cut height of 0.996
salmon_signed_network_cH0.997 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.997) #run the module detection with detect cut height of 0.997
salmon_signed_network_cH0.998 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.998) #run the module detection with detect cut height of 0.998
salmon_signed_network_cH0.999 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.999) #run the module detection with detect cut height of 0.999
plotDendroAndColors(salmon_signed_network_cH0.995$dendrograms[[1]], cbind(salmon_signed_network_cH0.995$colors, matchLabels(labels2colors(salmon_signed_network_cH0.996$colors), salmon_signed_network_cH0.995$colors), matchLabels(labels2colors(salmon_signed_network_cH0.997$colors), salmon_signed_network_cH0.995$colors), matchLabels(labels2colors(salmon_signed_network_cH0.998$colors), salmon_signed_network_cH0.995$colors), matchLabels(labels2colors(salmon_signed_network_cH0.999$colors), salmon_signed_network_cH0.995$colors)), c("detectCutHeight=0.995", "detectCutHeight=0.996", "detectCutHeight=0.997", "detectCutHeight=0.998", "detectCutHeight=0.998"), dendroLabels = F, hang=0.03, addGuide = T, guideHang=0.05) #plot the gene dendrogram and module assignments for the various tests of different merge cut heights
detectCutHeight.df <- data.frame(ch995=salmon_signed_network_cH0.995$colors, ch996=salmon_signed_network_cH0.996$colors, ch997=salmon_signed_network_cH0.997$colors, ch998=salmon_signed_network_cH0.998$colors, ch999=salmon_signed_network_cH0.999$colors) #create a data.frame of the module assignments from each of the merge cut height tests
ch995.consistentgenes <- sum(detectCutHeight.df$ch995==detectCutHeight.df$ch996 & detectCutHeight.df$ch995==detectCutHeight.df$ch997 & detectCutHeight.df$ch995==detectCutHeight.df$ch998 & detectCutHeight.df$ch995==detectCutHeight.df$ch999) #calculate the average number of genes assigned to the same module under the conditions of each test. In this case, all genes were assigned to the same modules for all the values tested

#Test the effect of the deep split parameter on the modules detected
salmon_signed_network_cH0.995_ds0 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 0, detectCutHeight = 0.995) #run the module detection with deep split of 0
salmon_signed_network_cH0.995_ds1 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 1, detectCutHeight = 0.995) #run the module detection with deep split of 1
salmon_signed_network_cH0.995_ds3 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 3, detectCutHeight = 0.995) #run the module detection with deep split of 3
salmon_signed_network_cH0.995_ds4 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.15, deepSplit = 4, detectCutHeight = 0.995) #run the module detection with deep split of 4
plotDendroAndColors(salmon_signed_network_cH0.995$dendrograms[[1]], cbind(salmon_signed_network_cH0.995_ds0$colors, matchLabels(labels2colors(salmon_signed_network_cH0.995_ds1$colors), salmon_signed_network_cH0.995_ds0$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995$colors), salmon_signed_network_cH0.995_ds0$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995_ds3$colors), salmon_signed_network_cH0.995_ds0$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995_ds4$colors), salmon_signed_network_cH0.995_ds0$colors)), c("deepSplit=0", "deepSplit=1", "deepSplit=2", "deepSplit=3", "deepSplit=4"), dendroLabels = F, hang=0.03, addGuide = T, guideHang=0.05) #plot the gene dendrogram and module assignments for the various tests of different deep split values
deepSplit.df <- data.frame(ds0=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds0$colors), salmon_signed_network_cH0.995_ds4$colors), ds1=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds1$colors), salmon_signed_network_cH0.995_ds4$colors), ds2=matchLabels(labels2colors(salmon_signed_network_cH0.995$colors), salmon_signed_network_cH0.995_ds4$colors), ds3=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds3$colors), salmon_signed_network_cH0.995_ds4$colors), ds4=salmon_signed_network_cH0.995_ds4$colors) #create a data.frame of the module assignments from each of the deep split tests
ds0.consistentgenes <- (sum(deepSplit.df$ds0==deepSplit.df$ds1) + sum(deepSplit.df$ds0==deepSplit.df$ds2) + sum(deepSplit.df$ds0==deepSplit.df$ds3) + sum(deepSplit.df$ds0==deepSplit.df$ds4))/4 #calculate the average number of genes sharing the same module assignment in the deep split=0 test compared to the other tests
ds1.consistentgenes <- (sum(deepSplit.df$ds1==deepSplit.df$ds0) + sum(deepSplit.df$ds1==deepSplit.df$ds2) + sum(deepSplit.df$ds1==deepSplit.df$ds3) + sum(deepSplit.df$ds1==deepSplit.df$ds4))/4 #calculate the average number of genes sharing the same module assignment in the deep split=1 test compared to the other tests
ds2.consistentgenes <- (sum(deepSplit.df$ds2==deepSplit.df$ds0) + sum(deepSplit.df$ds2==deepSplit.df$ds1) + sum(deepSplit.df$ds2==deepSplit.df$ds3) + sum(deepSplit.df$ds2==deepSplit.df$ds4))/4 #calculate the average number of genes sharing the same module assignment in the deep split=2 test compared to the other tests
ds3.consistentgenes <- (sum(deepSplit.df$ds3==deepSplit.df$ds0) + sum(deepSplit.df$ds3==deepSplit.df$ds1) + sum(deepSplit.df$ds3==deepSplit.df$ds2) + sum(deepSplit.df$ds3==deepSplit.df$ds4))/4 #calculate the average number of genes sharing the same module assignment in the deep split=3 test compared to the other tests
ds4.consistentgenes <- (sum(deepSplit.df$ds4==deepSplit.df$ds0) + sum(deepSplit.df$ds4==deepSplit.df$ds1) + sum(deepSplit.df$ds4==deepSplit.df$ds2) + sum(deepSplit.df$ds4==deepSplit.df$ds3))/4 #calculate the average number of genes sharing the same module assignment in the deep split=4 test compared to the other tests

#Test the effect of the merge cut height parameter on the modules detected
salmon_signed_network_cH0.995_ds2_mcH0.075 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.075, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a merge cut height of 0.075
salmon_signed_network_cH0.995_ds2_mcH0.1 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.1, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a merge cut height of 0.1
salmon_signed_network_cH0.995_ds2_mcH0.2 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.2, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a merge cut height of 0.2
salmon_signed_network_cH0.995_ds2_mcH0.25 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a merge cut height of 0.25
plotDendroAndColors(salmon_signed_network_cH0.995$dendrograms[[1]], cbind(salmon_signed_network_cH0.995_ds2_mcH0.075$colors, matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.1$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.2$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.25$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors)), c("mergeCutHeight=0.075", "mergeCutHeight=0.1", "mergeCutHeight=0.15", "mergeCutHeight=0.2", "mergeCutHeight=0.25"), dendroLabels = F, hang=0.03, addGuide = T, guideHang=0.05) #plot the gene dendrogram and module assignments for the various tests of different merge cut heights
mergeCutHeight.df <- data.frame(mcH075=salmon_signed_network_cH0.995_ds2_mcH0.075$colors, mcH01=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.1$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors), mcH015=matchLabels(labels2colors(salmon_signed_network_cH0.995$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors), mcH02=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.2$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors), mcH025=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.25$colors), salmon_signed_network_cH0.995_ds2_mcH0.075$colors)) #create a data.frame of the module assignments from each of the merge cut height tests
mcH075.consistentgenes <- (sum(mergeCutHeight.df$mcH075==mergeCutHeight.df$mcH01) + sum(mergeCutHeight.df$mcH075==mergeCutHeight.df$mcH015) + sum(mergeCutHeight.df$mcH075==mergeCutHeight.df$mcH02) + sum(mergeCutHeight.df$mcH075==mergeCutHeight.df$mcH025))/4 #calculate the average number of genes sharing the same module assignment in the merge cut height=0.075 test compared to the other tests
mcH01.consistentgenes <- (sum(mergeCutHeight.df$mcH01==mergeCutHeight.df$mcH075) + sum(mergeCutHeight.df$mcH01==mergeCutHeight.df$mcH015) + sum(mergeCutHeight.df$mcH01==mergeCutHeight.df$mcH02) + sum(mergeCutHeight.df$mcH01==mergeCutHeight.df$mcH025))/4 #calculate the average number of genes sharing the same module assignment in the merge cut height=0.1 test compared to the other tests
mcH015.consistentgenes <- (sum(mergeCutHeight.df$mcH015==mergeCutHeight.df$mcH075) + sum(mergeCutHeight.df$mcH015==mergeCutHeight.df$mcH01) + sum(mergeCutHeight.df$mcH015==mergeCutHeight.df$mcH02) + sum(mergeCutHeight.df$mcH015==mergeCutHeight.df$mcH025))/4 #calculate the average number of genes sharing the same module assignment in the merge cut height=0.15 test compared to the other tests
mcH02.consistentgenes <- (sum(mergeCutHeight.df$mcH02==mergeCutHeight.df$mcH075) + sum(mergeCutHeight.df$mcH02==mergeCutHeight.df$mcH01) + sum(mergeCutHeight.df$mcH02==mergeCutHeight.df$mcH015) + sum(mergeCutHeight.df$mcH02==mergeCutHeight.df$mcH025))/4 #calculate the average number of genes sharing the same module assignment in the merge cut height=0.2 test compared to the other tests
mcH025.consistentgenes <- (sum(mergeCutHeight.df$mcH025==mergeCutHeight.df$mcH075) + sum(mergeCutHeight.df$mcH025==mergeCutHeight.df$mcH01) + sum(mergeCutHeight.df$mcH025==mergeCutHeight.df$mcH015) + sum(mergeCutHeight.df$mcH025==mergeCutHeight.df$mcH02))/4 #calculate the average number of genes sharing the same module assignment in the merge cut height=0.25 test compared to the other tests

salmon_signed_network_cH0.995_ds2_mcH0.15_mms20 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 20, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a minimum module size of 20
salmon_signed_network_cH0.995_ds2_mcH0.15_mms40 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 40, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a minimum module size of 40
salmon_signed_network_cH0.995_ds2_mcH0.15_mms50 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 50, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.995) #run the module detection with a minimum module size of 50
plotDendroAndColors(salmon_signed_network_cH0.995$dendrograms[[1]], cbind(salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors, matchLabels(labels2colors(salmon_signed_network_cH0.995$colors), salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors), salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors), matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.15_mms50$colors), salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors)), c("minModuleSize=20", "minModuleSize=30", "minModuleSize=40", "minModuleSize=50"), dendroLabels = F, hang=0.03, addGuide = T, guideHang=0.05) #plot the gene dendrogram and module assignments for the various tests of different minimum module sizes
minModuleSize.df <- data.frame(mms20=salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors, mms30=matchLabels(labels2colors(salmon_signed_network_cH0.995$colors), salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors), mms40=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors), salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors), mms50=matchLabels(labels2colors(salmon_signed_network_cH0.995_ds2_mcH0.15_mms50$colors), salmon_signed_network_cH0.995_ds2_mcH0.15_mms20$colors)) #create a data.frame of the module assignments from each of the minimum module size tests
mms20.consistentgenes <- (sum(minModuleSize.df$mms20==minModuleSize.df$mms30) + sum(minModuleSize.df$mms20==minModuleSize.df$mms40) + sum(minModuleSize.df$mms20==minModuleSize.df$mms50))/3 #calculate the average number of genes sharing the same module assignment in the minimum module size=20 test compared to the other tests
mms30.consistentgenes <- (sum(minModuleSize.df$mms30==minModuleSize.df$mms20) + sum(minModuleSize.df$mms30==minModuleSize.df$mms40) + sum(minModuleSize.df$mms30==minModuleSize.df$mms50))/3 #calculate the average number of genes sharing the same module assignment in the minimum module size=30 test compared to the other tests
mms40.consistentgenes <- (sum(minModuleSize.df$mms40==minModuleSize.df$mms20) + sum(minModuleSize.df$mms40==minModuleSize.df$mms30) + sum(minModuleSize.df$mms40==minModuleSize.df$mms50))/3 #calculate the average number of genes sharing the same module assignment in the minimum module size=40 test compared to the other tests
mms50.consistentgenes <- (sum(minModuleSize.df$mms50==minModuleSize.df$mms20) + sum(minModuleSize.df$mms50==minModuleSize.df$mms30) + sum(minModuleSize.df$mms50==minModuleSize.df$mms40))/3 #calculate the average number of genes sharing the same module assignment in the minimum module size=50 test compared to the other tests

salmon_signed_MEs <- salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$MEs #copy the module eigengenes from the module assignment test with the most consistent parameters
WGCNA.plot <- plotDendroAndColors(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$dendrograms[[1]], salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors, "Modules", dendroLabels = F, hang=0.03, addGuide = T, guideHang=0.05) #plot the results from the chosen parameter set
jpeg("Figure2.jpeg", width=12, height=8, unit="in", res=300) #create an empty 12x8 jpeg 
print(plotDendroAndColors(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$dendrograms[[1]], salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors, "Modules", dendroLabels = F, hang=0.03, addGuide = T, guideHang=0.05)) #print the dendrogram figure to the jpeg
dev.off() #turn device off to stop printing items to the jpeg file

salmon_eigengene.df <- cbind(enviro_data, salmon_signed_MEs) #combine the field data and the module eigengene vales
salmon_eigengene.df$Site <- as.factor(salmon_eigengene.df$Site) #make sure Site is set as a factor variable
st_salmon_eigengene.df <- stdize(salmon_eigengene.df, omit=c("Sample_ID", "Date", "Site", "Sample_Period", "Sex")) #rescale all the eigengene values to have a mean of 0 and standard deviation of 1
start_date <- as.Date("2022-07-15") #set a day 0 for the study period
for(i in 1:nrow(st_salmon_eigengene.df)){
  st_salmon_eigengene.df$Date2[i] <- difftime(as.Date(salmon_eigengene.df$Date[i]), start_date, units = "days") #calculate the number of days from the start of the study period for each sampling event
}
st_salmon_eigengene.df$Sex <- as.factor(salmon_eigengene.df$Sex) #make sure the Sex variable is in factor format
st_salmon_eigengene.df$Site <- as.factor(st_salmon_eigengene.df$Site) #make sure the Site variable is in factor format

corrplot(cor(salmon_signed_MEs, enviro_data[,4:12])) #print a correlation plot showing direction and magnitude of correlations between module eigengenes and the environmental variables

plotEigengeneNetworks(salmon_signed_MEs, setLabels = colnames(salmon_signed_MEs)) #plot the dendrogram of modules 
save.image("WGCNA.RData")

load("WGCNA.RData") #if you've run the above code already, you can skip to this step 
Phylofish_GO <- fread("Phylofish_GOData.txt") #read in the Gene Ontology custom data frame from the Phylofish BioMart search (the Phylofish_GOData.txt file supplied in this repository has alread gone through next 3 lines to filter entries, otherwise it would have been to large to include in the repository)
Phylofish_GO_salmon <- Phylofish_GO[Phylofish_GO$Name %in% colnames(express_data),] #create a re-ordered and filtered GO data frame that matches the transcript order of the gene expression data frame
Phylofish_GO_salmon <- Phylofish_GO_salmon[Phylofish_GO_salmon$`Go name` != "",] #remove entries with no GO data from the Gene Ontology dataset (some entries had significant alignments to known genes, but had no GO data)
Phylofish_GO_salmon_input <- data.frame(gene=Phylofish_GO_salmon$Name, go_ID=Phylofish_GO_salmon$Code) #create a data frame with just gene IDs and GO codes to use in the GO enrichment analysis

geneModuleMembership <- as.data.frame(cor(express_data, salmon_signed_MEs, use = "p")) #calculate correlation coefficients between transcript count data and module eigengenes (module membership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples = nrow(express_data))) #calculate the p-values for the correlation coefficients
names(geneModuleMembership) <- paste("MM", substring(colnames(salmon_signed_MEs), 3), sep="") #set column names of the geneModuleMembership data frame (remove the MEs from in front of the module color names and add an MM) 
names(MMPvalue) <- paste("p.MM", substring(colnames(salmon_signed_MEs), 3), sep="") #set column names for the MMPvalue data frame (remove the MEs from in front of the module color names and add a p.MM)

daily_temps.df <- read.delim("~/Desktop/Brook_Trout_Heatwave/Brook_Trout_Heatwave_SampleSelect/BrookTrout_Heatwaves_DailyTemps.txt") #read in temperature data
daily_temps.df <- daily_temps.df[daily_temps.df$Date > "2022-07-14" & daily_temps.df$Date < "2022-08-23",] #filter temperature data to only include observations during the study period

#GAMs
darkorange.gam <- gam(z.MEdarkorange ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(darkorange.gam) #check the summary for the full gam model
darkorange.gam <- gam(z.MEdarkorange ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(darkorange.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="darkorange") #find the number of genes in the specified module
darkorange_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMdarkorange)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
darkorange_GO <- go_enrich(darkorange_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
darkorange_GO.refined <- refine(darkorange_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

white.gam <- gam(z.MEwhite ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(white.gam) #check the summary for the full gam model
white.gam <- gam(z.MEwhite ~ s(z.Length) + s(Date2) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(white.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="white") #find the number of genes in the specified module
white_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMwhite)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
white_GO <- go_enrich(white_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
white_GO.refined <- refine(white_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

lightgreen.gam <- gam(z.MElightgreen ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(lightgreen.gam) #check the summary for the full gam model
lightgreen.gam <- gam(z.MElightgreen ~ s(Date2) + + s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(lightgreen.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="lightgreen") #find the number of genes in the specified module
lightgreen_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMlightgreen)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
lightgreen_GO <- go_enrich(lightgreen_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
lightgreen_GO.refined <- refine(lightgreen_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

darkred.gam <- gam(z.MEdarkred ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(darkred.gam) #check the summary for the full gam model
darkred.gam <- gam(z.MEdarkred ~ s(Date2) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(darkred.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="darkred") #find the number of genes in the specified module
darkred_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMdarkred)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
darkred_GO <- go_enrich(darkred_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
darkred_GO.refined <- refine(darkred_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

orange.gam <- gam(z.MEorange ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(orange.gam) #check the summary for the full gam model
orange.gam <- gam(z.MEorange ~ s(Date2) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(orange.gam) #check the summary for the full gam model
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="orange") #find the number of genes in the specified module
orange_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMorange)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
orange_GO <- go_enrich(orange_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
orange_GO.refined <- refine(orange_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

pink.gam <- gam(z.MEpink ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(pink.gam) #check the summary for the full gam model
pink.gam <- gam(z.MEpink ~ s(z.Length) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(pink.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="pink") #find the number of genes in the specified module
pink_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMpink)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
pink_GO <- go_enrich(pink_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
pink_GO.refined <- refine(pink_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

yellow.gam <- gam(z.MEyellow ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(yellow.gam) #check the summary for the full gam model
yellow.gam <- gam(z.MEyellow ~ s(Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(yellow.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="yellow") #find the number of genes in the specified module
yellow_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMyellow)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
yellow_GO <- go_enrich(yellow_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
yellow_GO.refined <- refine(yellow_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

tan.gam <- gam(z.MEtan ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(tan.gam) #check the summary for the full gam model
tan.gam <- gam(z.MEtan ~ s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(tan.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="tan") #find the number of genes in the specified module
tan_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMtan)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
tan_GO <- go_enrich(tan_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
tan_GO.refined <- refine(tan_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

brown.gam <- gam(z.MEbrown ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(brown.gam) #check the summary for the full gam model
brown.gam <- gam(z.MEbrown ~ s(z.Length), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(brown.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="brown") #find the number of genes in the specified module
brown_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMbrown)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
brown_GO <- go_enrich(brown_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
brown_GO.refined <- refine(brown_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

turquoise.gam <- gam(z.MEturquoise ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(turquoise.gam) #check the summary for the full gam model
turquoise.gam <- gam(z.MEturquoise ~ s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(turquoise.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="turquoise") #find the number of genes in the specified module
turquoise_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMturquoise)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
turquoise_GO <- go_enrich(turquoise_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
turquoise_GO.refined <- refine(turquoise_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

greenyellow.gam <- gam(z.MEgreenyellow ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(greenyellow.gam) #check the summary for the full gam model
greenyellow.gam <- gam(z.MEgreenyellow ~ s(Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(greenyellow.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="greenyellow") #find the number of genes in the specified module
greenyellow_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMgreenyellow)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
greenyellow_GO <- go_enrich(greenyellow_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
greenyellow_GO.refined <- refine(greenyellow_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

blue.gam <- gam(z.MEblue ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(blue.gam) #check the summary for the full gam model
blue.gam <- gam(z.MEblue ~ s(Date2) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(blue.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="blue") #find the number of genes in the specified module
blue_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMblue)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
blue_GO <- go_enrich(blue_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
blue_GO.refined <- refine(blue_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

darkgreen.gam <- gam(z.MEdarkgreen ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(darkgreen.gam) #check the summary for the full gam model
darkgreen.gam <- gam(z.MEdarkgreen ~ s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(darkgreen.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="darkgreen") #find the number of genes in the specified module
darkgreen_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMdarkgreen)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
darkgreen_GO <- go_enrich(darkgreen_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
darkgreen_GO.refined <- refine(darkgreen_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

magenta.gam <- gam(z.MEmagenta ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(magenta.gam) #check the summary for the full gam model
magenta.gam <- gam(z.MEmagenta ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(magenta.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="magenta") #find the number of genes in the specified module
magenta_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMmagenta)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
magenta_GO <- go_enrich(magenta_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
magenta_GO.refined <- refine(magenta_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

lightcyan.gam <- gam(z.MElightcyan ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(lightcyan.gam) #check the summary for the full gam model
lightcyan.gam <- gam(z.MElightcyan ~ s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(lightcyan.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="lightcyan") #find the number of genes in the specified module
lightcyan_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMlightcyan)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
lightcyan_GO <- go_enrich(lightcyan_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
lightcyan_GO.refined <- refine(lightcyan_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

midnightblue.gam <- gam(z.MEmidnightblue ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(midnightblue.gam) #check the summary for the full gam model
midnightblue.gam <- gam(z.MEmidnightblue ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(midnightblue.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="midnightblue") #find the number of genes in the specified module
midnightblue_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMmidnightblue)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
midnightblue_GO <- go_enrich(midnightblue_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
midnightblue_GO.refined <- refine(midnightblue_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

darkturquoise.gam <- gam(z.MEdarkturquoise ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(darkturquoise.gam) #check the summary for the full gam model
darkturquoise.gam <- gam(z.MEdarkturquoise ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(darkturquoise.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="darkturquoise") #find the number of genes in the specified module
darkturquoise_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMdarkturquoise)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
darkturquoise_GO <- go_enrich(darkturquoise_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
darkturquoise_GO.refined <- refine(darkturquoise_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

cyan.gam <- gam(z.MEcyan ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(cyan.gam) #check the summary for the full gam model
cyan.gam <- gam(z.MEcyan ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(cyan.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="cyan") #find the number of genes in the specified module
cyan_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMcyan)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
cyan_GO <- go_enrich(cyan_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
cyan_GO.refined <- refine(cyan_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

green.gam <- gam(z.MEgreen ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(green.gam) #check the summary for the full gam model
green.gam <- gam(z.MEgreen ~ s(Date2) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(green.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="green") #find the number of genes in the specified module
green_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMgreen)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
green_GO <- go_enrich(green_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
green_GO.refined <- refine(green_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

skyblue.gam <- gam(z.MEskyblue ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(skyblue.gam) #check the summary for the full gam model
skyblue.gam <- gam(z.MEskyblue ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(skyblue.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="skyblue") #find the number of genes in the specified module
skyblue_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMskyblue)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
skyblue_GO <- go_enrich(skyblue_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
skyblue_GO.refined <- refine(skyblue_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

black.gam <- gam(z.MEblack ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(black.gam) #check the summary for the full gam model
black.gam <- gam(z.MEblack ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(black.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="black") #find the number of genes in the specified module
black_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMblack)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
black_GO <- go_enrich(black_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
black_GO.refined <- refine(black_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

red.gam <- gam(z.MEred ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(red.gam) #check the summary for the full gam model
red.gam <- gam(z.MEred ~ s(Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(red.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="red") #find the number of genes in the specified module
red_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMred)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
red_GO <- go_enrich(red_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
red_GO.refined <- refine(red_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

steelblue.gam <- gam(z.MEsteelblue ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(steelblue.gam) #check the summary for the full gam model
steelblue.gam <- gam(z.MEsteelblue ~ s(Date2) + s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(steelblue.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="steelblue") #find the number of genes in the specified module
steelblue_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMsteelblue)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
steelblue_GO <- go_enrich(steelblue_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
steelblue_GO.refined <- refine(steelblue_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

royalblue.gam <- gam(z.MEroyalblue ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(royalblue.gam) #check the summary for the full gam model
royalblue.gam <- gam(z.MEroyalblue ~ s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(royalblue.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="royalblue") #find the number of genes in the specified module
royalblue_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMroyalblue)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
royalblue_GO <- go_enrich(royalblue_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
royalblue_GO.refined <- refine(royalblue_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

violet.gam <- gam(z.MEviolet ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(violet.gam) #check the summary for the full gam model
violet.gam <- gam(z.MEviolet ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(violet.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="violet") #find the number of genes in the specified module
violet_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMviolet)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
violet_GO <- go_enrich(violet_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
violet_GO.refined <- refine(violet_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

paleturquoise.gam <- gam(z.MEpaleturquoise ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(paleturquoise.gam) #check the summary for the full gam model
paleturquoise.gam <- gam(z.MEpaleturquoise ~ s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(paleturquoise.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="paleturquoise") #find the number of genes in the specified module
paleturquoise_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMpaleturquoise)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
paleturquoise_GO <- go_enrich(paleturquoise_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
paleturquoise_GO.refined <- refine(paleturquoise_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

saddlebrown.gam <- gam(z.MEsaddlebrown ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(saddlebrown.gam) #check the summary for the full gam model
saddlebrown.gam <- gam(z.MEsaddlebrown ~ s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(saddlebrown.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="saddlebrown") #find the number of genes in the specified module
saddlebrown_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMsaddlebrown)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
saddlebrown_GO <- go_enrich(saddlebrown_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
saddlebrown_GO.refined <- refine(saddlebrown_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

grey60.gam <- gam(z.MEgrey60 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(grey60.gam) #check the summary for the full gam model
grey60.gam <- gam(z.MEgrey60 ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(grey60.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="grey60") #find the number of genes in the specified module
grey60_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMgrey60)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
grey60_GO <- go_enrich(grey60_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
grey60_GO.refined <- refine(grey60_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

purple.gam <- gam(z.MEpurple ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(purple.gam) #check the summary for the full gam model
purple.gam <- gam(z.MEpurple ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(purple.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="purple") #find the number of genes in the specified module
purple_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMpurple)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
purple_GO <- go_enrich(purple_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
purple_GO.refined <- refine(purple_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

salmon.gam <- gam(z.MEsalmon ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(salmon.gam) #check the summary for the full gam model
salmon.gam <- gam(z.MEsalmon ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(salmon.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="salmon") #find the number of genes in the specified module
salmon_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMsalmon)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
salmon_GO <- go_enrich(salmon_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
salmon_GO.refined <- refine(salmon_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

darkgrey.gam <- gam(z.MEdarkgrey ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(darkgrey.gam) #check the summary for the full gam model
darkgrey.gam <- gam(z.MEdarkgrey ~ s(Date2) + s(z.DailyAvgTemp, by=Site), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(darkgrey.gam) #recheck the summary
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="darkgrey") #find the number of genes in the specified module
darkgrey_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMdarkgrey)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
darkgrey_GO <- go_enrich(darkgrey_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
darkgrey_GO.refined <- refine(darkgrey_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

lightyellow.gam <- gam(z.MElightyellow ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(lightyellow.gam) #check the summary for the full gam model
lightyellow.gam <- gam(z.MElightyellow ~ s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(lightyellow.gam) #check the summary for the full gam model
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="lightyellow") #find the number of genes in the specified module
lightyellow_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMlightyellow)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
lightyellow_GO <- go_enrich(lightyellow_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
lightyellow_GO.refined <- refine(lightyellow_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

grey.gam <- gam(z.MEgrey ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(grey.gam) #check the summary for the full gam model
grey.gam <- gam(z.MEgrey ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_eigengene.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(grey.gam) #check the summary for the full gam model
sum(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors=="grey") #find the number of genes in the specified module
grey_GO_input <- data.frame(gene=row.names(geneModuleMembership), candidate=-log(MMPvalue$p.MMgrey)) #prepare the inputs for the GO analysis with the gene names and -log p-values of module membership for the specified module
grey_GO <- go_enrich(grey_GO_input, annotations = Phylofish_GO_salmon_input, test="wilcoxon", silent = F) #run the GO analysis with a Wilcoxon rank-sum test using the module membership p-value as the score and the custom GO database from Phylofish for the annotations
grey_GO.refined <- refine(grey_GO, fwer=0.05, annotations = Phylofish_GO_salmon_input) #refine the GO results, removing genes from significant child categories and re-tests the significance of GO categories to favor more specific results in the GO analysis

newdata.df <- data.frame(z.DailyAvgTemp=rep(rep(rep(c(seq(min(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="Big Poe",]$z.DailyAvgTemp), max(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="Big Poe",]$z.DailyAvgTemp), length.out=250), seq(min(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="East Branch",]$z.DailyAvgTemp), max(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="East Branch",]$z.DailyAvgTemp), length.out=250), seq(min(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="Shaver's Creek",]$z.DailyAvgTemp), max(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="Shaver's Creek",]$z.DailyAvgTemp), length.out=250), seq(min(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="Standing Stone",]$z.DailyAvgTemp), max(st_salmon_eigengene.df[st_salmon_eigengene.df$Site=="Standing Stone",]$z.DailyAvgTemp), length.out=250)), 35), 100), 2), Site=rep(rep(rep(c(rep("Big Poe", 250), rep("East Branch", 250), rep("Shaver's Creek", 250), rep("Standing Stone", 250)), 35), 100), 2), Date2=rep(rep(sort(rep(1:35, 1000)), 100), 2), z.Length=rep(sort(rep(seq(min(st_salmon_eigengene.df$z.Length), max(st_salmon_eigengene.df$z.Length), length.out=100), 35000)), 2), Sex=c(rep("M", 3500000), rep("F", 3500000))) #create a new data frame with a range of variables to plug into GAMs to generate plots of the predicted eigengene values across various ranges of the 
newdata.df$DailyAvgTemp <- rep(rep(rep(c(seq(min(salmon_eigengene.df[salmon_eigengene.df$Site=="Big Poe",]$DailyAvgTemp), max(salmon_eigengene.df[salmon_eigengene.df$Site=="Big Poe",]$DailyAvgTemp), length.out=250), seq(min(salmon_eigengene.df[salmon_eigengene.df$Site=="East Branch",]$DailyAvgTemp), max(salmon_eigengene.df[salmon_eigengene.df$Site=="East Branch",]$DailyAvgTemp), length.out=250), seq(min(salmon_eigengene.df[salmon_eigengene.df$Site=="Shaver's Creek",]$DailyAvgTemp), max(salmon_eigengene.df[salmon_eigengene.df$Site=="Shaver's Creek",]$DailyAvgTemp), length.out=250), seq(min(salmon_eigengene.df[salmon_eigengene.df$Site=="Standing Stone",]$DailyAvgTemp), max(salmon_eigengene.df[salmon_eigengene.df$Site=="Standing Stone",]$DailyAvgTemp), length.out=250)), 35), 100), 2) #add a column to the model input data with the real values of temperature aligned with the standardized values so they can be used on the plot
newdata.df$Length <- rep(sort(rep(seq(min(salmon_eigengene.df$Length), max(salmon_eigengene.df$Length), length.out=100), 35000)), 2) #add a column to the model input data with the real values of length aligned with the standardized values so they can be used on the plot
newdata1.df <- newdata.df[newdata.df$Sex=="M", ] #pull out input options of males only so calculating predicted values goes faster for models without the Sex parameter
newdata2.df <- newdata.df[newdata.df$Length==107,] #pull out input options of one length only so calculating predicted values goes faster for models without the Length parameter
newdata3.df <- newdata.df[newdata.df$Length==107 & newdata.df$Sex=="M",] #pull out input options of males of one size only so calculating predicted values goes faster for models without the Sex and Length parameters
newdata4.df <- newdata.df[newdata.df$Date2==17,] #pull out input options of fish from one date only so calculating predicted values goes faster for models without the Date parameter

blue_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(blue.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEblue) + mean(salmon_eigengene.df$MEblue), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=blue_plotdata.df[blue_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=blue_plotdata.df[blue_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=blue_plotdata.df[blue_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=blue_plotdata.df[blue_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEblue, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Blue Module Eigengene") #edit the axis labels

cyan_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(cyan.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEcyan) + mean(salmon_eigengene.df$MEcyan), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
cyan.plot <- ggplot()+ 
  geom_line(data=cyan_plotdata.df[cyan_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value for each site
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEcyan, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(italic("cyan")))#edit the axis labels

darkgreen_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(darkgreen.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEdarkgreen) + mean(salmon_eigengene.df$MEdarkgreen), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=darkgreen_plotdata.df[darkgreen_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value for each site
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkgreen, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Dark Green Module Eigengene") #edit the axis labels

darkgrey_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(darkgrey.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEdarkgrey) + mean(salmon_eigengene.df$MEdarkgrey), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=darkgrey_plotdata.df[darkgrey_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=darkgrey_plotdata.df[darkgrey_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=darkgrey_plotdata.df[darkgrey_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=darkgrey_plotdata.df[darkgrey_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkgrey, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size=14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Dark Grey Module Eigengene") #edit the axis labels

darkorange_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(darkorange.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEdarkorange) + mean(salmon_eigengene.df$MEdarkorange), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkorange, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Dark Orange Module Eigengene") #edit the axis labels

darkorange.plot <- ggplot()+ 
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkorange, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  labs(x=" ", y=" ", title=expression(italic("darkorange"))) #edit the axis labels

darkorange_plotdata.df <- tibble(z.DailyAvgTemp=newdata4.df$z.DailyAvgTemp, DailyAvgTemp=newdata4.df$DailyAvgTemp, z.preds=predict.gam(darkorange.gam, newdata=newdata4.df), preds=z.preds * sd(salmon_eigengene.df$MEdarkorange) + mean(salmon_eigengene.df$MEdarkorange), Date2=newdata4.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata4.df$Site, Length=newdata4.df$Length, z.Length=newdata4.df$z.Length, Sex=newdata4.df$Sex) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Length==107 & darkorange_plotdata.df$Sex=="M",], aes(x=DailyAvgTemp, y=preds, color="Min Length"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Length==211 & darkorange_plotdata.df$Sex=="M",], aes(x=DailyAvgTemp, y=preds, color="Max Length"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Length > 158.75 & darkorange_plotdata.df$Length < 159.75 & darkorange_plotdata.df$Sex=="M",], aes(x=DailyAvgTemp, y=preds, color="Median Length"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  scale_color_manual(name="Length", values=c("Min Length"="blue", "Max Length"="red", "Median Length"="palegoldenrod"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkorange, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Dark Orange Module Eigengene") #edit the axis labels
ggplot()+ 
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Length > 158.75 & darkorange_plotdata.df$Length < 159.75 & darkorange_plotdata.df$Sex=="M",], aes(x=DailyAvgTemp, y=preds, color="Male"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=darkorange_plotdata.df[darkorange_plotdata.df$Length > 158.75 & darkorange_plotdata.df$Length < 159.75 & darkorange_plotdata.df$Sex=="F",], aes(x=DailyAvgTemp, y=preds, color="Female"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  scale_color_manual(name="Sex", values=c("Male"="lightblue", "Female"="pink"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkorange, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Dark Orange Module Eigengene") #edit the axis labels

darkred_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(darkred.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEdarkred) + mean(salmon_eigengene.df$MEdarkred), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=darkred_plotdata.df[darkred_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=darkred_plotdata.df[darkred_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=darkred_plotdata.df[darkred_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=darkred_plotdata.df[darkred_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkred, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Dark Red Module Eigengene") #edit the axis labels

darkturquoise_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(darkturquoise.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEdarkturquoise) + mean(salmon_eigengene.df$MEdarkturquoise), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
darkturquoise.plot <- ggplot()+ 
  geom_line(data=darkturquoise_plotdata.df[darkturquoise_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEdarkturquoise, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(italic("darkturquoise"))) #edit the axis labels

grey_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(grey.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEgrey) + mean(salmon_eigengene.df$MEgrey), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=grey_plotdata.df[grey_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=grey_plotdata.df[grey_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=grey_plotdata.df[grey_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=grey_plotdata.df[grey_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEgrey, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Grey Module Eigengene") #edit the axis labels

grey60_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(grey60.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEgrey60) + mean(salmon_eigengene.df$MEgrey60), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
grey60.plot <- ggplot()+ 
  geom_line(data=grey60_plotdata.df[grey60_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEgrey60, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(italic("grey60"))) #edit the axis labels

lightcyan_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(lightcyan.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MElightcyan) + mean(salmon_eigengene.df$MElightcyan), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=lightcyan_plotdata.df[lightcyan_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=lightcyan_plotdata.df[lightcyan_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=lightcyan_plotdata.df[lightcyan_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=lightcyan_plotdata.df[lightcyan_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MElightcyan, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Light Cyan Module Eigengene") #edit the axis labels

lightgreen_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(lightgreen.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MElightgreen) + mean(salmon_eigengene.df$MElightgreen), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=lightgreen_plotdata.df[lightgreen_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=lightgreen_plotdata.df[lightgreen_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=lightgreen_plotdata.df[lightgreen_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=lightgreen_plotdata.df[lightgreen_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MElightgreen, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Light Green Module Eigengene") #edit the axis labels
lightgreen.plot <- ggplot()+ 
  geom_line(data=lightgreen_plotdata.df[lightgreen_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MElightgreen, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  labs(x=" ", y=" ", title=expression(italic("lightgreen"))) #edit the axis labels

lightyellow_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(lightyellow.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MElightyellow) + mean(salmon_eigengene.df$MElightyellow), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=lightyellow_plotdata.df[lightyellow_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MElightyellow, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Light Yellow Module Eigengene") #edit the axis labels

midnightblue_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(midnightblue.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEmidnightblue) + mean(salmon_eigengene.df$MEmidnightblue), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=midnightblue_plotdata.df[midnightblue_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEmidnightblue, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Midnight Blue Module Eigengene")+ #edit the axis labels
  ylim(-0.2, 0.2)

orange_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(orange.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEorange) + mean(salmon_eigengene.df$MEorange), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=orange_plotdata.df[orange_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=orange_plotdata.df[orange_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=orange_plotdata.df[orange_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=orange_plotdata.df[orange_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEorange, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Orange Module Eigengene") #edit the axis labels

paleturquoise_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(paleturquoise.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEpaleturquoise) + mean(salmon_eigengene.df$MEpaleturquoise), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=paleturquoise_plotdata.df[paleturquoise_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEpaleturquoise, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Pale Turquoise Module Eigengene") #edit the axis labels

pink_plotdata.df <- tibble(z.DailyAvgTemp=newdata4.df$z.DailyAvgTemp, DailyAvgTemp=newdata4.df$DailyAvgTemp, z.preds=predict.gam(pink.gam, newdata=newdata4.df), preds=z.preds * sd(salmon_eigengene.df$MEpink) + mean(salmon_eigengene.df$MEpink), Date2=newdata4.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata4.df$Site, Length=newdata4.df$Length, z.Length=newdata4.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=pink_plotdata.df, aes(x=Length, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=Length, y=MEpink, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Length (mm)", y="Pink Module Eigengene") #edit the axis labels

purple_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(purple.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEpurple) + mean(salmon_eigengene.df$MEpurple), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=purple_plotdata.df[purple_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=purple_plotdata.df[purple_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=purple_plotdata.df[purple_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=purple_plotdata.df[purple_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEpurple, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Purple Module Eigengene") #edit the axis labels
purple.plot <- ggplot()+ 
  geom_line(data=purple_plotdata.df[purple_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEpurple, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  labs(x=" ", y=" ", title=expression(italic("purple"))) #edit the axis labels

royalblue_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(royalblue.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEroyalblue) + mean(salmon_eigengene.df$MEroyalblue), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=royalblue_plotdata.df[royalblue_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEroyalblue, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Royal Blue Module Eigengene") #edit the axis labels

saddlebrown_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(saddlebrown.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEsaddlebrown) + mean(salmon_eigengene.df$MEsaddlebrown), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=saddlebrown_plotdata.df[saddlebrown_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=saddlebrown_plotdata.df[saddlebrown_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=saddlebrown_plotdata.df[saddlebrown_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=saddlebrown_plotdata.df[saddlebrown_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEsaddlebrown, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Saddlebrown Module Eigengene") #edit the axis labels

salmon_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(salmon.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEsalmon) + mean(salmon_eigengene.df$MEsalmon), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
salmon.plot <- ggplot()+ 
  geom_line(data=salmon_plotdata.df[salmon_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEsalmon, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(italic("salmon"))) #edit the axis labels

steelblue_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(steelblue.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEsteelblue) + mean(salmon_eigengene.df$MEsteelblue), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=steelblue_plotdata.df[steelblue_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=steelblue_plotdata.df[steelblue_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=steelblue_plotdata.df[steelblue_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=steelblue_plotdata.df[steelblue_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEsteelblue, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Steel Blue Module Eigengene") #edit the axis labels

tan_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(tan.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEtan) + mean(salmon_eigengene.df$MEtan), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
tan.plot <- ggplot()+ 
  geom_line(data=tan_plotdata.df[tan_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEtan, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 16)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Tan Module Eigengene") #edit the axis labels

violet_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(violet.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEviolet) + mean(salmon_eigengene.df$MEviolet), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=violet_plotdata.df[violet_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEviolet, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Violet Module Eigengene") #edit the axis labels

white_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(white.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_eigengene.df$MEwhite) + mean(salmon_eigengene.df$MEwhite), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of module eigengenes based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=white_plotdata.df[white_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=white_plotdata.df[white_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=white_plotdata.df[white_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=white_plotdata.df[white_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_eigengene.df, aes(x=DailyAvgTemp, y=MEwhite, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic(base_size = 14)+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="Steel Blue Module Eigengene") #edit the axis labels

stream.leg <- get_legend(tan.plot)

#Figures for Paper
eigengene.plots <- ggarrange(cyan.plot, darkorange.plot, darkturquoise.plot, grey60.plot, lightgreen.plot, purple.plot, salmon.plot, stream.leg, nrow = 2, ncol=4) #create multiplot with eigengenes that are mentioned in the text, with a separate panel for the legend for all plots
eigengene.plots <- annotate_figure(eigengene.plots, bottom = text_grob("24-h Average Temperature (ºC)", size=16), left=text_grob("Module Eigengene Value", size=16, rot=90)) #add shared axis labels for all plots
ggsave("Figure3.jpeg", eigengene.plots, width=16, height=8, units="in", dpi=300) #save the multiplot as figure 3

#WGCNA dendro plots
salmon_signed_eigencolors <- salmon_signed_MEs #save a copy of the eigengene data
colnames(salmon_signed_eigencolors) <- substring(colnames(salmon_signed_MEs), first=3, last=40L) #rename the eigengene columns to drop the ME in front of the color names
plotEigengeneNetworks(salmon_signed_eigencolors, setLabels = colnames(salmon_signed_eigencolors), plotHeatmaps = T, colorLabels = TRUE, signed = T, coloredBarplot = T, barplotMeans = F, barplotErrors = F, plotPreservation = "standard") #plot the WGCNA dendrogram
eigen_tree <- hclust(dist(t(salmon_signed_eigencolors))) #create a dendrogram of the eigengenes for each module 
eigentree.plot <- ggtree(eigen_tree)+ #plot the eigengene dendrogram
  geom_tiplab(align=T, geom="text")+ #add tip labels aligned to the end nodes of the dendrogram
  geom_tiplab(aes(label=("XXXX")), linesize = 0, align=T, geom="label", color=colnames(salmon_signed_eigencolors), fill=colnames(salmon_signed_eigencolors), offset = 0.3)+ #set the tip labels to show up as blocks of color 
  theme(plot.margin=unit(c(0.2, 0, 2.8, 0.5), "cm")) #adjust plot margins

gam.heat.df <- read.table("ME_gam_heatmap.txt", header = T) #read in table of compiled P-values for all environmental variables in the GAMMs (ns p-values were made to equal 1 for easier plotting)
gam.heat.df$Eigengene <- factor(gam.heat.df$Eigengene, levels=unique(gam.heat.df$Eigengene)) #make sure the module variable in the p-value data frame is a factor in the same order as the dendrogram
gam_heat.df <- gam.heat.df %>% as_tibble() %>% gather(key="Variable", value="p-val", Temperature:Temp.Date, factor_key = T) #create a long version of the p-value data frame
eigen_heatplot.plot <- ggplot(data=gam_heat.df, aes(x=Variable, y=fct_rev(Eigengene), fill=-log(`p-val`)))+ #create a new ggplot object from the long p-value data frame with GAMM variables on the x axis, modules on the y axis, and colored by log transformed p-values
  geom_tile()+ #create a colored tile geom
  scale_fill_gradient(low="white", high="red")+ #set the high p-values to white and low p-values to red 
  theme_classic(base_size = 14)+ #increase the base font size and change the plot theme 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) #remove the module names and y-axis (it will be aligned with the dendrogram)

ggsave(filename="Figure2-2.jpeg", grid.arrange(eigentree.plot, eigen_heatplot.plot, nrow=1), width=13.5, height = 8.5, units="in", dpi=300) #print a jpeg of the GAMM variable p-value heatmap with module dendrogram


