library(devtools)
#install_github("YosefLab/ImpulseDE2")
library(ImpulseDE2)
library(data.table)
library(stringr)
library(dplyr)
library(GOfuncR)
library(ggplot2)
library(MuMIn)
library(mgcv)
library(gridExtra)
library(WGCNA)
library(ggpubr)
library(poppr)
library(ggtree)
library(dendroextras)
library(tidyr)
library(ggarrow)

salmon_data <- fread("~/Desktop/BrookTrout_ReadCounts/BrookTrout_ReadCounts/salmon_logcounts_environ.txt", drop=1) #read in dataset with environmental and expression data
enviro_data <- salmon_data[,1:16] #create a separate data frame with just environmental data
express_data <- salmon_data[,-(1:16)] #create a separate data frame with just expression data
rownames(express_data) <- salmon_data$Sample_ID #name the expression data rows with sample IDs so future analyses of the expression data can be aligned to environmental data
salmon_files <- list.files("~/Desktop/BrookTrout_ReadCounts/salmon_counts/trimmed", pattern = "quant") #create a list of the file names with salmon read counts
for(i in 1:length(salmon_files)){ #create a for loop to read in the read count files and automatically name them with the sample ID from the file name
  readcount <- read.delim(str_c("~/Desktop/BrookTrout_ReadCounts/salmon_counts/trimmed/", salmon_files[i], "/quant.sf"))
  assign(str_c(str_split_1(salmon_files[i], "_")[1],"_salmon.rc"), readcount)
}

salmon_count.list <- objects(pattern = "_salmon.rc") #create a list with all of the salmon read count objects in the R environment
salmon_counts <- NULL #create a null object to turn into a table with the salmon read counts from all of the samples
for(i in 1:length(salmon_count.list)){ #create a for loop to create a list with tables for each sample with just the gene names and read counts
  salmon_counts[i] <- list(get(salmon_count.list[i]) %>% dplyr::select("Name", "NumReads"))
}
salmon_counts <- salmon_counts %>% purrr::reduce(left_join, by="Name") #match up the read counts by gene name across all the samples and put them into a single data frame
colnames(salmon_counts) <- c("Name", lapply(salmon_files, function(x) str_split_1(x, pattern= "_")[1])) #rename the columns of the data frame with sample IDs
sample_annot <- data.frame(Sample=enviro_data$Sample_ID, Condition="case", Time=as.numeric(enviro_data$Sample_Period), TimeCateg=as.character(enviro_data$Sample_Period), Batch=as.factor(enviro_data$Site)) #arrange data for ImpulseDE design into a data frame; includes sample ID, a Condition column (only one condition, case), Time a numeric column for the sample period, TimeCateg a character-class version of the Time column, and Batch a column with the 4 stream names as factors
salmon_counts_filtered <- salmon_counts[salmon_counts$Name %in% colnames(express_data),] #filter the count data so that it only includes genes that passed the occurrence filter
salmon_counts_filtered <- t(salmon_counts_filtered) #transpose the count data so it can be converted to a matrix
salmon_counts_filtered <- salmon_counts_filtered[-c(1, 56), ] #remove the gene name and outlier sample EB044 columns from the transposed read count data
salmon_counts_filtered.int <- matrix(as.integer(salmon_counts_filtered), nrow(salmon_counts_filtered), ncol(salmon_counts_filtered)) #create a numeric matrix of the read count data, rounded to nearest integer for ImpulseDE
salmon_counts_filtered.int <- t(salmon_counts_filtered.int) #re-transpose the integer count data matrix so it can be added to the sample_annot ImpulseDE input data frame
colnames(salmon_counts_filtered.int) <- sample_annot$Sample #add sample ids as the column names for the integer count data matrix
row.names(salmon_counts_filtered.int) <- salmon_counts$Name[salmon_counts$Name %in% colnames(express_data)] #add the gene names as the row names to the integer count data matrix
rownames(sample_annot) <- sample_annot$Sample #copy sample ids in the sample annot ImpulseDE data frame to the row names (otherwise results will lack the gene IDs)

sample_annot2 <- data.frame(Sample=enviro_data$Sample_ID, Condition="case", Time=as.numeric(enviro_data$Sample_Period), TimeCateg=as.character(enviro_data$Sample_Period), Batch=as.factor(enviro_data$Site)) #create a copy of the sample_annot input data frame
for(i in 1:nrow(sample_annot2)){ #replace the times during the second heatwave (5,6,7,8) with the corresponding times (1,2,3,4) during the first heatwave
  sample_annot2$Time[i] <- if(sample_annot2$Time[i] <= 4) sample_annot2$Time[i] else
    if(sample_annot2$Time[i]==5) 1 else
      if(sample_annot2$Time[i]==6) 2 else
        if(sample_annot2$Time[i]==7) 3 else 4
}
sample_annot2$TimeCateg <- as.character(sample_annot2$Time) #update the TimeCateg column in the sample_annot copy to reflect the changes in the Time column
row.names(sample_annot2) <- sample_annot2$Sample #copy sample ids in the sample annot ImpulseDE data frame to the row names (otherwise results will lack the gene IDs)

heatwavept.riverbatch.iDE2test <- runImpulseDE2(matCountData = salmon_counts_filtered.int, dfAnnotation = sample_annot2, boolCaseCtrl = F, vecConfounders = c("Batch"), scaNProc = 1, scaQThres = 0.05) #run ImpulseDE2 with streams as confounders and corresponding time points from the two heatwaves as replicates 
heatwave.riverbatch.impulse.plots <- plotGenes(scaNTopIDs = 43, objectImpulseDE2 = heatwavept.riverbatch.iDE2test, boolCaseCtrl = F)
print(heatwave.riverbatch.impulse.plots[[16]])

heatwave_DEgenes.vec <- heatwavept.riverbatch.iDE2test@vecDEGenes #get a vector of the significantly differentially expressed gene names

#Alternatively, to avoid running ImpulseDE you can set the list of DE genes in the manuscript here
heatwave_DEgenes.vec <- c("QSF_AN32B.3.8", "QSF_CDC37.1.3", "QSF_CIRBP.24.24", "QSF_CIRBP.7.24", "QSF_CSF3R.3.3", "QSF_DH12A.2.2", "QSF_ES1.1.3", "QSF_HMGB2.9.11", "QSF_HS90A.1.1", "QSF_HS90B.1.1", "QSF_HSP47.2.3", "QSF_HSP47.3.3", "QSF_HSP70B.1.1", "QSF_HSPA4L.1.1", "QSF_KHDR1.2.2", "QSF_LOC100136615.4.4", "QSF_LOC100194703.44.106", "QSF_LOC100691603.1.1", "QSF_LOC100696101.2.2", "QSF_LOC100696517.1.1", "QSF_LOC100699058.1.1", "QSF_LOC101166160.1.1", "QSF_LOC101166550.1.1", "QSF_LOC562935.1.1", "QSF_P2RX4.82.83", "QSF_P4HA1.6.11", "QSF_P4HA2.1.3","QSF_PCY2.5.29", "QSF_RO32.1.2", "QSF_SARNP.8.8", "QSF_SFRS6.2.5", "QSF_SRSF2.6.8", "QSF_TCPG.1.3", "QSF_THOC4.1.2", "QSF_TIA1.4.5", "QSF_TRA2B.3.7", "QSF_TRAD1.20.35", "QSF_TRI39.3.8", "QSF_TSC22D3.2.2", "QSF_UGT1A5.1.1", "QSF_USP9X.2.4", "QSF_VRK3.145.924", "QSF_ZGC_162825.1.2")
salmon_DEgenes.df <- subset(salmon_data, select=colnames(salmon_data) %in% heatwave_DEgenes.vec) #subset count data for the 43 differentially expressed genes
salmon_DEgenes.df <- cbind(enviro_data, salmon_DEgenes.df) #create data frame that combines field data and differentially expressed transcript count data for each sample
salmon_DEgenes.df$Site <- as.factor(salmon_DEgenes.df$Site) #make sure the Site variable is a factor
st_salmon_DEgenes.df <- stdize(salmon_DEgenes.df, omit=c("Sample_ID", "Date", "Site", "Sample_Period", "Sex")) #standardize all of the quantitative variables to have a mean of 0 and standard deviation of 1
start_date <- as.Date("2022-07-15") #set the start date
for(i in 1:nrow(st_salmon_DEgenes.df)){ #create a second Date variable that expresses date as the number of since the start date set above (the day before sampling started)
  st_salmon_DEgenes.df$Date2[i] <- difftime(as.Date(salmon_DEgenes.df$Date[i]), start_date, units = "days")
}
st_salmon_DEgenes.df$Sex <- as.factor(salmon_DEgenes.df$Sex) #make sure the sex column is a factor

#LOC100699058.1.1
LOC100699058.1.1.gam <- gam(z.QSF_LOC100699058.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC100699058.1.1.gam) #check the summary for the full gam model
LOC100699058.1.1.gam <- gam(z.QSF_LOC100699058.1.1 ~ s(z.Length) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #re-run the model with just the significant model terms
summary(LOC100699058.1.1.gam) #recheck the summary

#HS90A.1.1
HS90A.1.1.gam <- gam(z.QSF_HS90A.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HS90A.1.1.gam) #check the summary for the full gam model
HS90A.1.1.gam <- gam(z.QSF_HS90A.1.1 ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HS90A.1.1.gam) #recheck the summary

#HSPA4L.1.1
HSPA4L.1.1.gam <- gam(z.QSF_HSPA4L.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HSPA4L.1.1.gam) #check the summary for the full gam model
HSPA4L.1.1.gam <- gam(z.QSF_HSPA4L.1.1 ~ s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HSPA4L.1.1.gam) #recheck the summary

#TRA2B.3.7
TRA2B.3.7.gam <- gam(z.QSF_TRA2B.3.7 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(TRA2B.3.7.gam) #check the summary for the full gam model
TRA2B.3.7.gam <- gam(z.QSF_TRA2B.3.7 ~ s(Date2) + s(z.DailyAvgTemp) + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(TRA2B.3.7.gam) #recheck the summary

#TCPG.1.3
TCPG.1.3.gam <- gam(z.QSF_TCPG.1.3 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(TCPG.1.3.gam) #check the summary for the full gam model
TCPG.1.3.gam <- gam(z.QSF_TCPG.1.3 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run gam with just the significant terms
summary(TCPG.1.3.gam) #recheck the summary

#ES1.1.3
ES1.1.3.gam <- gam(z.QSF_ES1.1.3 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(ES1.1.3.gam) #check the summary for the full gam model
ES1.1.3.gam <- gam(z.QSF_ES1.1.3 ~ s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(ES1.1.3.gam) #recheck the summary

#AN32B.3.8
AN32B.3.8.gam <- gam(z.QSF_AN32B.3.8 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(AN32B.3.8.gam) #check the summary for the full gam model
AN32B.3.8.gam <- gam(z.QSF_AN32B.3.8 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(AN32B.3.8.gam) #recheck the summary

#UGT1A5.1.1
UGT1A5.1.1.gam <- gam(z.QSF_UGT1A5.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(UGT1A5.1.1.gam) #check the summary for the full gam model
UGT1A5.1.1.gam <- gam(z.QSF_UGT1A5.1.1 ~ s(Date2) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(UGT1A5.1.1.gam) #recheck the summary

#LOC100696517.1.1
LOC100696517.1.1.gam <- gam(z.QSF_LOC100696517.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC100696517.1.1.gam) #check the summary for the full gam model
LOC100696517.1.1.gam <- gam(z.QSF_LOC100696517.1.1 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the signficant terms
summary(LOC100696517.1.1.gam) #recheck the summary

#HS90B.1.1
HS90B.1.1.gam <- gam(z.QSF_HS90B.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HS90B.1.1.gam) #check the summary for the full gam model
HS90B.1.1.gam <- gam(z.QSF_HS90B.1.1 ~ s(Date2) + s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HS90B.1.1.gam) #recheck the summary

#PCY2.5.29
PCY2.5.29.gam <- gam(z.QSF_PCY2.5.29 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(PCY2.5.29.gam) #check the summary for the full gam model
PCY2.5.29.gam <- gam(z.QSF_PCY2.5.29 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(PCY2.5.29.gam) #recheck the summary

#P4HA2.1.3
P4HA2.1.3.gam <- gam(z.QSF_P4HA2.1.3 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(P4HA2.1.3.gam) #check the summary for the full gam model
P4HA2.1.3.gam <- gam(z.QSF_P4HA2.1.3 ~ s(z.Length) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(P4HA2.1.3.gam) #recheck the summary

#VRK3.145.924
VRK3.145.924.gam <- gam(z.QSF_VRK3.145.924 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(VRK3.145.924.gam) #check the summary for the full gam model
VRK3.145.924.gam <- gam(z.QSF_VRK3.145.924 ~ s(z.Length), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(VRK3.145.924.gam) #recheck the summary

#CIRBP.7.24
CIRBP.7.24.gam <- gam(z.QSF_CIRBP.7.24 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(CIRBP.7.24.gam) #check the summary for the full gam model
CIRBP.7.24.gam <- gam(z.QSF_CIRBP.7.24 ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(CIRBP.7.24.gam) #recheck the summary

#ZGC_162825.1.2
ZGC_162825.1.2.gam <- gam(z.QSF_ZGC_162825.1.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(ZGC_162825.1.2.gam) #check the summary for the full gam model
ZGC_162825.1.2.gam <- gam(z.QSF_ZGC_162825.1.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(ZGC_162825.1.2.gam) #recheck the summary

#LOC101166550.1.1
LOC101166550.1.1.gam <- gam(z.QSF_LOC101166550.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC101166550.1.1.gam) #check the summary for the full gam model
LOC101166550.1.1.gam <- gam(z.QSF_LOC101166550.1.1 ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC101166550.1.1.gam) #recheck the summary

#USP9X.2.4
USP9X.2.4.gam <- gam(z.QSF_USP9X.2.4 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(USP9X.2.4.gam) #check the summary for the full gam model
USP9X.2.4.gam <- gam(z.QSF_USP9X.2.4 ~ s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(USP9X.2.4.gam) #recheck the summary

#TRI39.3.8
TRI39.3.8.gam <- gam(z.QSF_TRI39.3.8 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(TRI39.3.8.gam) #check the summary for the full gam model
TRI39.3.8.gam <- gam(z.QSF_TRI39.3.8 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(TRI39.3.8.gam) #recheck the summary

#TSC22D3.2.2
TSC22D3.2.2.gam <- gam(z.QSF_TSC22D3.2.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(TSC22D3.2.2.gam) #check the summary for the full gam model
TSC22D3.2.2.gam <- gam(z.QSF_TSC22D3.2.2 ~ s(Date2) + s(z.DailyAvgTemp, by=Site) + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(TSC22D3.2.2.gam) #recheck the summary

#SFRS6.2.5
SFRS6.2.5.gam <- gam(z.QSF_SFRS6.2.5 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(SFRS6.2.5.gam) #check the summary for the full gam model
SFRS6.2.5.gam <- gam(z.QSF_SFRS6.2.5 ~ s(Date2) + s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(SFRS6.2.5.gam) #recheck the summary

#LOC100194703.44.106
LOC100194703.44.106.gam <- gam(z.QSF_LOC100194703.44.106 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC100194703.44.106.gam) #check the summary for the full gam model
LOC100194703.44.106.gam <- gam(z.QSF_LOC100194703.44.106 ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC100194703.44.106.gam) #recheck the summary

#DH12A.2.2
DH12A.2.2.gam <- gam(z.QSF_DH12A.2.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(DH12A.2.2.gam) #check the summary for the full gam model
DH12A.2.2.gam <- gam(z.QSF_DH12A.2.2 ~ s(Date2) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(DH12A.2.2.gam) #recheck the summary

#TRAD1.20.35
TRAD1.20.35.gam <- gam(z.QSF_TRAD1.20.35 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(TRAD1.20.35.gam) #check the summary for the full gam model
TRAD1.20.35.gam <- gam(z.QSF_TRAD1.20.35 ~ s(Date2) + s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(TRAD1.20.35.gam) #recheck the summary

#CSF3R.3.3
CSF3R.3.3.gam <- gam(z.QSF_CSF3R.3.3 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(CSF3R.3.3.gam) #check the summary for the full gam model
CSF3R.3.3.gam <- gam(z.QSF_CSF3R.3.3 ~ s(Date2) + s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(CSF3R.3.3.gam) #recheck the summary

#THOC4.1.2
THOC4.1.2.gam <- gam(z.QSF_THOC4.1.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(THOC4.1.2.gam) #check the summary for the full gam model
THOC4.1.2.gam <- gam(z.QSF_THOC4.1.2 ~ s(Date2) + s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(THOC4.1.2.gam) #recheck the summary

#HSP47.3.3
HSP47.3.3.gam <- gam(z.QSF_HSP47.3.3 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HSP47.3.3.gam) #check the summary for the full gam model
HSP47.3.3.gam <- gam(z.QSF_HSP47.3.3 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HSP47.3.3.gam) #recheck the summary

#HSP47.2.3
HSP47.2.3.gam <- gam(z.QSF_HSP47.2.3 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HSP47.2.3.gam) #check the summary for the full gam model
HSP47.2.3.gam <- gam(z.QSF_HSP47.2.3 ~ factor(Sex) + s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HSP47.2.3.gam) #recheck the summary

#P4HA1.6.11
P4HA1.6.11.gam <- gam(z.QSF_P4HA1.6.11 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(P4HA1.6.11.gam) #check the summary for the full gam model
P4HA1.6.11.gam <- gam(z.QSF_P4HA1.6.11 ~ s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(P4HA1.6.11.gam) #recheck the summary

#LOC100136615.4.4
LOC100136615.4.4.gam <- gam(z.QSF_LOC100136615.4.4 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC100136615.4.4.gam) #check the summary for the full gam model
LOC100136615.4.4.gam <- gam(z.QSF_LOC100136615.4.4 ~ s(z.Length) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC100136615.4.4.gam) #recheck the summary

#LOC100696101.2.2
LOC100696101.2.2.gam <- gam(z.QSF_LOC100696101.2.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC100696101.2.2.gam) #check the summary for the full gam model
LOC100696101.2.2.gam <- gam(z.QSF_LOC100696101.2.2 ~ s(z.DailyAvgTemp) + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC100696101.2.2.gam) #recheck the summary

#TIA1.4.5
TIA1.4.5.gam <- gam(z.QSF_TIA1.4.5 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(TIA1.4.5.gam) #check the summary for the full gam model
TIA1.4.5.gam <- gam(z.QSF_TIA1.4.5 ~ s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(TIA1.4.5.gam) #recheck the summary

#LOC101166160.1.1
LOC101166160.1.1.gam <- gam(z.QSF_LOC101166160.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC101166160.1.1.gam) #check the summary for the full gam model
LOC101166160.1.1.gam <- gam(z.QSF_LOC101166160.1.1 ~ s(z.Length) + s(Date2) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC101166160.1.1.gam) #recheck the summary

#LOC100691603.1.1
LOC100691603.1.1.gam <- gam(z.QSF_LOC100691603.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC100691603.1.1.gam) #check the summary for the full gam model
LOC100691603.1.1.gam <- gam(z.QSF_LOC100691603.1.1 ~ s(Date2) + s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC100691603.1.1.gam) #recheck the summary

#CDC37.1.3
CDC37.1.3.gam <- gam(z.QSF_CDC37.1.3 ~ s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model (failed to converge with the Sex variable, so it was removed)
summary(CDC37.1.3.gam) #check the summary for the full gam model
CDC37.1.3.gam <- gam(z.QSF_CDC37.1.3 ~ s(Date2) + s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(CDC37.1.3.gam) #recheck the summary

#LOC562935.1.1
LOC562935.1.1.gam <- gam(z.QSF_LOC562935.1.1 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(LOC562935.1.1.gam) #check the summary for the full gam model
LOC562935.1.1.gam <- gam(z.QSF_LOC562935.1.1 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(LOC562935.1.1.gam) #recheck the summary

#SARNP.8.8
SARNP.8.8.gam <- gam(z.QSF_SARNP.8.8 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(SARNP.8.8.gam) #check the summary for the full gam model
SARNP.8.8.gam <- gam(z.QSF_SARNP.8.8 ~ s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(SARNP.8.8.gam) #recheck the summary

#RO32.1.2
RO32.1.2.gam <- gam(z.QSF_RO32.1.2 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(RO32.1.2.gam) #check the summary for the full gam model
RO32.1.2.gam <- gam(z.QSF_RO32.1.2 ~ s(z.Length) + s(z.DailyAvgTemp), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(RO32.1.2.gam) #recheck the summary

#SRSF2.6.8
SRSF2.6.8.gam <- gam(z.QSF_SRSF2.6.8 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(SRSF2.6.8.gam) #check the summary for the full gam model
SRSF2.6.8.gam <- gam(z.QSF_SRSF2.6.8 ~ factor(Sex) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(SRSF2.6.8.gam) #recheck the summary

#HMGB2.9.11
HMGB2.9.11.gam <- gam(z.QSF_HMGB2.9.11 ~ factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HMGB2.9.11.gam) #check the summary for the full gam model
HMGB2.9.11.gam <- gam(z.QSF_HMGB2.9.11 ~ s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HMGB2.9.11.gam) #recheck the summary

#KHDR1.2.2
KHDR1.2.2.gam <- gam(z.QSF_KHDR1.2.2 ~  s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(KHDR1.2.2.gam) #check the summary for the full gam model
KHDR1.2.2.gam <- gam(z.QSF_KHDR1.2.2 ~  s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(KHDR1.2.2.gam) #recheck the summary

#CIRBP.24.24
CIRBP.24.24.gam <- gam(z.QSF_CIRBP.24.24 ~  factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(CIRBP.24.24.gam) #check the summary for the full gam model
CIRBP.24.24.gam <- gam(z.QSF_CIRBP.24.24 ~  s(z.DailyAvgTemp) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(CIRBP.24.24.gam) #recheck the summary

#HSP70B.1.1
HSP70B.1.1.gam <- gam(z.QSF_HSP70B.1.1 ~  factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(HSP70B.1.1.gam) #check the summary for the full gam model
HSP70B.1.1.gam <- gam(z.QSF_HSP70B.1.1 ~  s(z.Length) + s(z.DailyAvgTemp, by=Site), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(HSP70B.1.1.gam) #recheck the summary

#P2RX4.82.83
P2RX4.82.83.gam <- gam(z.QSF_P2RX4.82.83 ~  factor(Sex) + s(z.Length) + s(Date2) + s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re") + ti(z.DailyAvgTemp, Date2), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the full gam model
summary(P2RX4.82.83.gam) #check the summary for the full gam model
P2RX4.82.83.gam <- gam(z.QSF_P2RX4.82.83 ~  s(z.DailyAvgTemp) + s(z.DailyAvgTemp, by=Site) + s(Site, bs="re"), data=st_salmon_DEgenes.df, na.action = "na.fail", gamma=2, select=T) #run the gam with just the significant terms
summary(P2RX4.82.83.gam) #recheck the summary


#GAM Plots
newdata.df <- data.frame(z.DailyAvgTemp=rep(rep(rep(c(seq(min(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="Big Poe",]$z.DailyAvgTemp), max(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="Big Poe",]$z.DailyAvgTemp), length.out=250), seq(min(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="East Branch",]$z.DailyAvgTemp), max(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="East Branch",]$z.DailyAvgTemp), length.out=250), seq(min(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="Shaver's Creek",]$z.DailyAvgTemp), max(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="Shaver's Creek",]$z.DailyAvgTemp), length.out=250), seq(min(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="Standing Stone",]$z.DailyAvgTemp), max(st_salmon_DEgenes.df[st_salmon_DEgenes.df$Site=="Standing Stone",]$z.DailyAvgTemp), length.out=250)), 35), 100), 2), Site=rep(rep(rep(c(rep("Big Poe", 250), rep("East Branch", 250), rep("Shaver's Creek", 250), rep("Standing Stone", 250)), 35), 100), 2), Date2=rep(rep(sort(rep(1:35, 1000)), 100), 2), z.Length=rep(sort(rep(seq(min(st_salmon_DEgenes.df$z.Length), max(st_salmon_DEgenes.df$z.Length), length.out=100), 35000)), 2), Sex=c(rep("M", 3500000), rep("F", 3500000))) #create a new data frame with a range of variables to plug into GAMs to generate plots of the predicted eigengene values across various ranges of the 
newdata.df$DailyAvgTemp <- rep(rep(rep(c(seq(min(salmon_DEgenes.df[salmon_DEgenes.df$Site=="Big Poe",]$DailyAvgTemp), max(salmon_DEgenes.df[salmon_DEgenes.df$Site=="Big Poe",]$DailyAvgTemp), length.out=250), seq(min(salmon_DEgenes.df[salmon_DEgenes.df$Site=="East Branch",]$DailyAvgTemp), max(salmon_DEgenes.df[salmon_DEgenes.df$Site=="East Branch",]$DailyAvgTemp), length.out=250), seq(min(salmon_DEgenes.df[salmon_DEgenes.df$Site=="Shaver's Creek",]$DailyAvgTemp), max(salmon_DEgenes.df[salmon_DEgenes.df$Site=="Shaver's Creek",]$DailyAvgTemp), length.out=250), seq(min(salmon_DEgenes.df[salmon_DEgenes.df$Site=="Standing Stone",]$DailyAvgTemp), max(salmon_DEgenes.df[salmon_DEgenes.df$Site=="Standing Stone",]$DailyAvgTemp), length.out=250)), 35), 100), 2) #add a column to the model input data with the real values of temperature aligned with the standardized values so they can be used on the plot
newdata.df$Length <- rep(sort(rep(seq(min(salmon_DEgenes.df$Length), max(salmon_DEgenes.df$Length), length.out=100), 35000)), 2) #add a column to the model input data with the real values of length aligned with the standardized values so they can be used on the plot
newdata1.df <- newdata.df[newdata.df$Sex=="M", ] #pull out input options of males only so calculating predicted values goes faster for models without the Sex parameter
newdata2.df <- newdata.df[newdata.df$Length==107,] #pull out input options of one length only so calculating predicted values goes faster for models without the Length parameter
newdata3.df <- newdata.df[newdata.df$Length==107 & newdata.df$Sex=="M",] #pull out input options of males of one size only so calculating predicted values goes faster for models without the Sex and Length parameters
newdata4.df <- newdata.df[newdata.df$Date2==17,] #pull out input options of fish from one date only so calculating predicted values goes faster for models without the Date parameter
study_predict <- data.frame(z.Length=rep(0, 144), Sex=rep("M", 144), Date=rep(seq(from=as.Date("2022-07-15"), to=as.Date("2022-08-19"), by=1), 4), Date2=rep(0:35, 4), Site=c(rep("Big Poe", 36), rep("East Branch", 36), rep("Shaver's Creek", 36), rep("Standing Stone", 36))) #create a data frame with one point for each day during the study period in each site (using average total length and male for sex)
daily_temps.df <- read.delim("~/Desktop/Brook_Trout_Heatwave/Brook_Trout_Heatwave_SampleSelect/BrookTrout_Heatwaves_DailyTemps.txt") #read in daily stream temperature data 
daily_temps.df$Date <- as.Date(daily_temps.df$Date) #convert the stream temperature data column Date from character to Date format
study_predict <- left_join(study_predict, daily_temps.df, by=c("Date", "Site")) #add the temperature readings to the predict data frame we made above
study_predict$z.DailyAvgTemp <- (study_predict$Temp.avg - mean(enviro_data$DailyAvgTemp))/sd(enviro_data$DailyAvgTemp) #create a new variable transforming the daily average temperature into the standardized values used in the GAMMs 
express_means <- salmon_DEgenes.df %>% group_by(Date, Site) %>% summarise(across(starts_with("QSF"), ~ mean(.x)))

#AN32B.3.8
AN32B.3.8_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(AN32B.3.8.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_AN32B.3.8) + mean(salmon_DEgenes.df$QSF_AN32B.3.8), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
cirbp.plot <- ggplot()+ 
  geom_line(data=AN32B.3.8_plotdata.df[AN32B.3.8_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_AN32B.3.8, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y=" ", title=expression(paste("AN32B.3.8 (", italic("cirbp"), ")"))) #edit the axis labels

AN32B.3.8_backcast.df <- tibble(study_predict, z.preds=predict.gam(AN32B.3.8.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_AN32B.3.8) + mean(salmon_DEgenes.df$QSF_AN32B.3.8))
cirbp_predict1.plot <- ggplot()+
  #geom_line(data=AN32B.3.8_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_AN32B.3.8, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_point(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site), size=4)+
  geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4, color=Site), alpha=1, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ . *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("AN32B.3.8 (", italic("cirbp"), ")"))) #edit the axis labels
  
cirbp_predict2.plot <- ggplot()+
  geom_line(data=AN32B.3.8_backcast.df, aes(x=Date, y=preds, color=Site), linewidth=2)+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_AN32B.3.8, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  #geom_point(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site), size=4)+
  #geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4), alpha=0.6, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ . *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("AN32B.3.8 (", italic("cirbp"), ")"))) #edit the axis labels

#CDC37.1.3
CDC37.1.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(CDC37.1.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_CDC37.1.3) + mean(salmon_DEgenes.df$QSF_CDC37.1.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=CDC37.1.3_plotdata.df[CDC37.1.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=CDC37.1.3_plotdata.df[CDC37.1.3_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=CDC37.1.3_plotdata.df[CDC37.1.3_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=CDC37.1.3_plotdata.df[CDC37.1.3_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_CDC37.1.3, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="CDC37.1.3 expression (log reads per million)") #edit the axis labels
cdc37.plot <- ggplot()+ 
  geom_line(data=CDC37.1.3_plotdata.df[CDC37.1.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_CDC37.1.3, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title = expression(paste("CDC37.1.3 (", italic("cdc37"),")"))) #edit the axis labels

CDC37.1.3_backcast.df <- tibble(study_predict, z.preds=predict.gam(CDC37.1.3.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_CDC37.1.3) + mean(salmon_DEgenes.df$QSF_CDC37.1.3))
cdc37_predict.plot <- ggplot()+
  geom_line(data=CDC37.1.3_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_CDC37.1.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/3.5), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ . *3.5, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("CDC37.1.3 (", italic("cdc37"), ")"))) #edit the axis labels

#CIRBP.24.24
CIRBP.24.24_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(CIRBP.24.24.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_CIRBP.24.24) + mean(salmon_DEgenes.df$QSF_CIRBP.24.24), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
cirbpb1.plot <- ggplot()+ 
  geom_line(data=CIRBP.24.24_plotdata.df[CIRBP.24.24_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_CIRBP.24.24, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title = expression(paste("CIRBP.24.24 (", italic("cirbpb"),")"))) #edit the axis labels

CIRBP.24.24_backcast.df <- tibble(study_predict, z.preds=predict.gam(CIRBP.24.24.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_CIRBP.24.24) + mean(salmon_DEgenes.df$QSF_CIRBP.24.24))
cirbp1_predict.plot <- ggplot()+
  geom_line(data=CIRBP.24.24_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_CIRBP.24.24, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 1), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -1) *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("CIRBP.24.24 (", italic("cirbp"), ")"))) #edit the axis labels


#CIRBP.7.24
CIRBP.7.24_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(CIRBP.7.24.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_CIRBP.7.24) + mean(salmon_DEgenes.df$QSF_CIRBP.7.24), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=CIRBP.7.24_plotdata.df[CIRBP.7.24_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=CIRBP.7.24_plotdata.df[CIRBP.7.24_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=CIRBP.7.24_plotdata.df[CIRBP.7.24_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=CIRBP.7.24_plotdata.df[CIRBP.7.24_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_CIRBP.7.24, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="CIRBP.7.24 expression (log reads per million)") #edit the axis labels
cirbpb.2.plot <- ggplot()+ 
  geom_line(data=CIRBP.7.24_plotdata.df[CIRBP.7.24_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_CIRBP.7.24, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title = expression(paste("CIRBP.7.24 (", italic("cirbpb"),")"))) #edit the axis labels

CIRBP.7.24_backcast.df <- tibble(study_predict, z.preds=predict.gam(CIRBP.7.24.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_CIRBP.7.24) + mean(salmon_DEgenes.df$QSF_CIRBP.7.24))
cirbp2_predict.plot <- ggplot()+
  geom_line(data=CIRBP.7.24_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_CIRBP.7.24, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 2), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -2) *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("CIRBP.7.24 (", italic("cirbp"), ")"))) #edit the axis labels

#CSF3R.3.3
CSF3R.3.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(CSF3R.3.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_CSF3R.3.3) + mean(salmon_DEgenes.df$QSF_CSF3R.3.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=CSF3R.3.3_plotdata.df[CSF3R.3.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=CSF3R.3.3_plotdata.df[CSF3R.3.3_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=CSF3R.3.3_plotdata.df[CSF3R.3.3_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=CSF3R.3.3_plotdata.df[CSF3R.3.3_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_CSF3R.3.3, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="CSF3R.3.3 expression (log reads per million)") #edit the axis labels

CSF3R.3.3_backcast.df <- tibble(study_predict, z.preds=predict.gam(CSF3R.3.3.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_CSF3R.3.3) + mean(salmon_DEgenes.df$QSF_CSF3R.3.3))
csf3r_predict.plot <- ggplot()+
  geom_line(data=CSF3R.3.3_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_CSF3R.3.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 1), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -1) *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("CSF3R.3.3 (", italic("csf3r"), ")"))) #edit the axis labels

#DH12A.2.2
DH12A.2.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(DH12A.2.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_DH12A.2.2) + mean(salmon_DEgenes.df$QSF_DH12A.2.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=DH12A.2.2_plotdata.df[DH12A.2.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=DH12A.2.2_plotdata.df[DH12A.2.2_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=DH12A.2.2_plotdata.df[DH12A.2.2_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=DH12A.2.2_plotdata.df[DH12A.2.2_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_DH12A.2.2, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="DH12A.2.2 expression (log reads per million)") #edit the axis labels

DH12A.2.2_backcast.df <- tibble(study_predict, z.preds=predict.gam(DH12A.2.2.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_DH12A.2.2) + mean(salmon_DEgenes.df$QSF_DH12A.2.2))
hsd20b2_predict.plot <- ggplot()+
  geom_line(data=DH12A.2.2_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_DH12A.2.2, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 1), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -1) *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("DH12A.2.2 (", italic("hsd20b2"), ")"))) #edit the axis labels

#ES1.1.3
ES1.1.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(ES1.1.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_ES1.1.3) + mean(salmon_DEgenes.df$QSF_ES1.1.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
gatd3.plot <- ggplot()+ 
  geom_line(data=ES1.1.3_plotdata.df[ES1.1.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_ES1.1.3, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title = expression(paste("ES1.1.3 (", italic("gatd3"),")"))) #edit the axis labels

ES1.1.3_backcast.df <- tibble(study_predict, z.preds=predict.gam(ES1.1.3.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_ES1.1.3) + mean(salmon_DEgenes.df$QSF_ES1.1.3))
gatd3_predict.plot <- ggplot()+
  geom_line(data=ES1.1.3_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_ES1.1.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4.5), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ . *4.5, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("ES1.1.3 (", italic("gatd3"), ")"))) #edit the axis labels

#HMGB2.9.11
HMGB2.9.11_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HMGB2.9.11.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HMGB2.9.11) + mean(salmon_DEgenes.df$QSF_HMGB2.9.11), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=HMGB2.9.11_plotdata.df[HMGB2.9.11_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HMGB2.9.11, color=Site), alpha=0.3)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="HMGB2.9.11 expression (log reads per million)") #edit the axis labels

HMGB2.9.11_backcast.df <- tibble(study_predict, z.preds=predict.gam(HMGB2.9.11.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HMGB2.9.11) + mean(salmon_DEgenes.df$QSF_HMGB2.9.11))
hmgb2b_predict.plot <- ggplot()+
  geom_line(data=HMGB2.9.11_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HMGB2.9.11, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/6 +3.5), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -3.5) *6, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HMGB2.9.11 (", italic("hmgb2b"), ")"))) #edit the axis labels

#HS90A.1.1
HS90A.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HS90A.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HS90A.1.1) + mean(salmon_DEgenes.df$QSF_HS90A.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=HS90A.1.1_plotdata.df[HS90A.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=HS90A.1.1_plotdata.df[HS90A.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=HS90A.1.1_plotdata.df[HS90A.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=HS90A.1.1_plotdata.df[HS90A.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HS90A.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="HS90A.1.1 expression (log reads per million)") #edit the axis labels
hsp90ab1.1.plot <- ggplot()+ 
  geom_line(data=HS90A.1.1_plotdata.df[HS90A.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HS90A.1.1, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("HSP90A.1.1 (", italic("hsp90ab1"),")"))) #edit the axis labels

HS90A.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(HS90A.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HS90A.1.1) + mean(salmon_DEgenes.df$QSF_HS90A.1.1))
hsp90ab1_predict.plot <- ggplot()+
  geom_line(data=HS90A.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HS90A.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/3 - 2), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. + 2) *3, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HS90A.1.1 (", italic("hsp90ab1"), ")"))) #edit the axis labels

#HS90B.1.1
HS90B.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HS90B.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HS90B.1.1) + mean(salmon_DEgenes.df$QSF_HS90B.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=HS90B.1.1_plotdata.df[HS90B.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=HS90B.1.1_plotdata.df[HS90B.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=HS90B.1.1_plotdata.df[HS90B.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=HS90B.1.1_plotdata.df[HS90B.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HS90B.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="HS90B.1.1 expression (log reads per million)") #edit the axis labels
hsp90ab1.2.plot <- ggplot()+ 
  geom_line(data=HS90B.1.1_plotdata.df[HS90B.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HS90B.1.1, color=Site), alpha=0.3, show.legend=F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("HSP90B.1.1 (", italic("hsp90ab1"),")"))) #edit the axis labels

HS90B.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(HS90B.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HS90B.1.1) + mean(salmon_DEgenes.df$QSF_HS90B.1.1))
hsp90ab1.2_predict.plot <- ggplot()+
  geom_line(data=HS90B.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HS90B.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/2.5 + 2.5), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -2.5) *2.5, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HS90B.1.1 (", italic("hsp90ab1"), ")"))) #edit the axis labels

#HSP47.2.3
HSP47.2.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HSP47.2.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP47.2.3) + mean(salmon_DEgenes.df$QSF_HSP47.2.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
serpinh1.1.plot <- ggplot()+ 
  geom_line(data=HSP47.2.3_plotdata.df[HSP47.2.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HSP47.2.3, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("HSP47.2.3 (", italic("serpinh1b"),")"))) #edit the axis labels

HSP47.2.3_backcast.df <- tibble(study_predict, z.preds=predict.gam(HSP47.2.3.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP47.2.3) + mean(salmon_DEgenes.df$QSF_HSP47.2.3))
serpinh1.1_predict.plot <- ggplot()+
  geom_line(data=HSP47.2.3_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HSP47.2.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/3-4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. +4) *3, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HSP47.2.3 (", italic("serpinh1b"), ")"))) #edit the axis labels

#HSP47.3.3
HSP47.3.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HSP47.3.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP47.3.3) + mean(salmon_DEgenes.df$QSF_HSP47.3.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
serpinh1b.2.plot <- ggplot()+ 
  geom_line(data=HSP47.3.3_plotdata.df[HSP47.3.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HSP47.3.3, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("HSP47.3.3 (", italic("serpinh1b"),")")))+ #edit the axis labels
  ylim(2, 8.5)

HSP47.3.3_backcast.df <- tibble(study_predict, z.preds=predict.gam(HSP47.3.3.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP47.3.3) + mean(salmon_DEgenes.df$QSF_HSP47.3.3))
serpinh1.2_predict.plot <- ggplot()+
  geom_line(data=HSP47.3.3_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HSP47.3.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/2-4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. +4) *2, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HSP47.3.3 (", italic("serpinh1b"), ")"))) #edit the axis labels

#HSP70B.1.1
HSP70B.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HSP70B.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP70B.1.1) + mean(salmon_DEgenes.df$QSF_HSP70B.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
hsp70.plot <- ggplot()+ 
  geom_line(data=HSP70B.1.1_plotdata.df[HSP70B.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HSP70B.1.1, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("HSP70B.1.1 (", italic("hsp70"),")"))) #edit the axis labels

HSP70B.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(HSP70B.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP70B.1.1) + mean(salmon_DEgenes.df$QSF_HSP70B.1.1))
hsp70_predict.plot <- ggplot()+
  geom_line(data=HSP70B.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HSP70B.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/2-4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. +4) *2, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HSP70B.1.1 (", italic("hsp70"), ")"))) #edit the axis labels

#HSPA4L.1.1
HSPA4L.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HSPA4L.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSPA4L.1.1) + mean(salmon_DEgenes.df$QSF_HSPA4L.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
hspa4l.plot <- ggplot()+ 
  geom_line(data=HSPA4L.1.1_plotdata.df[HSPA4L.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HSPA4L.1.1, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("HSPA4L.1.1 (", italic("hspa4l"),")")))+ #edit the axis labels
  ylim(-0.5, 3.1)

HSPA4L.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(HSPA4L.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSPA4L.1.1) + mean(salmon_DEgenes.df$QSF_HSPA4L.1.1))
hspa4l_predict.plot <- ggplot()+
  geom_line(data=HSPA4L.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HSPA4L.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4-3), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. +3) *4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("HSPA4L.1.1 (", italic("hspa4l"), ")"))) #edit the axis labels

#KHDR1.2.2
KHDR1.2.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(KHDR1.2.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_KHDR1.2.2) + mean(salmon_DEgenes.df$QSF_KHDR1.2.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
khdrbs1.plot <- ggplot()+ 
  geom_line(data=KHDR1.2.2_plotdata.df[KHDR1.2.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_KHDR1.2.2, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("KHDR1.2.2 (", italic("khdrbs1"),")"))) #edit the axis labels

KHDR1.2.2_backcast.df <- tibble(study_predict, z.preds=predict.gam(KHDR1.2.2.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_KHDR1.2.2) + mean(salmon_DEgenes.df$QSF_KHDR1.2.2))
khdrbs1_predict.plot <- ggplot()+
  geom_line(data=KHDR1.2.2_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_KHDR1.2.2, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/6 + 4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -4)*6, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("KHDR1.2.2 (", italic("khdrbs1"), ")"))) #edit the axis labels

#LOC100136615.4.4
LOC100136615.4.4_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC100136615.4.4.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100136615.4.4) + mean(salmon_DEgenes.df$QSF_LOC100136615.4.4), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC100136615.4.4_plotdata.df[LOC100136615.4.4_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100136615.4.4, color=Site), alpha=0.3)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC100136615.4.4 expression (log reads per million)") #edit the axis labels

LOC100136615.4.4_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC100136615.4.4.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100136615.4.4) + mean(salmon_DEgenes.df$QSF_LOC100136615.4.4))
LOC129840284_predict.plot <- ggplot()+
  geom_line(data=LOC100136615.4.4_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC100136615.4.4, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 1), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -1)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC100136615.4.4 (", italic("uncharacterized LOC129840284"), ")"))) #edit the axis labels

#LOC100194703.44.106
LOC100194703.44.106_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC100194703.44.106.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100194703.44.106) + mean(salmon_DEgenes.df$QSF_LOC100194703.44.106), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC100194703.44.106_plotdata.df[LOC100194703.44.106_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=LOC100194703.44.106_plotdata.df[LOC100194703.44.106_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=LOC100194703.44.106_plotdata.df[LOC100194703.44.106_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=LOC100194703.44.106_plotdata.df[LOC100194703.44.106_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100194703.44.106, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC100194703.44.106 expression (log reads per million)") #edit the axis labels

LOC100194703.44.106_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC100194703.44.106.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100194703.44.106) + mean(salmon_DEgenes.df$QSF_LOC100194703.44.106))
grin2ab_predict.plot <- ggplot()+
  geom_line(data=LOC100194703.44.106_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC100194703.44.106, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 -3), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. +3)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC100194703.44.106 (", italic("grin2ab"), ")"))) #edit the axis labels

#LOC100691603.1.1
LOC100691603.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC100691603.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100691603.1.1) + mean(salmon_DEgenes.df$QSF_LOC100691603.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC100691603.1.1_plotdata.df[LOC100691603.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=LOC100691603.1.1_plotdata.df[LOC100691603.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=LOC100691603.1.1_plotdata.df[LOC100691603.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=LOC100691603.1.1_plotdata.df[LOC100691603.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100691603.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC100691603.1.1 expression (log reads per million)") #edit the axis labels
dnaja.plot <- ggplot()+ 
  geom_line(data=LOC100691603.1.1_plotdata.df[LOC100691603.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100691603.1.1, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("LOC100691603.1.1 (", italic("dnaja"),")"))) #edit the axis labels

LOC100691603.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC100691603.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100691603.1.1) + mean(salmon_DEgenes.df$QSF_LOC100691603.1.1))
dnaja_predict.plot <- ggplot()+
  geom_line(data=LOC100691603.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC100691603.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 -1), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. +1)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC100691603.1.1 (", italic("dnaja"), ")"))) #edit the axis labels


#LOC100696101.2.2
LOC100696101.2.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC100696101.2.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100696101.2.2) + mean(salmon_DEgenes.df$QSF_LOC100696101.2.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC100696101.2.2_plotdata.df[LOC100696101.2.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=LOC100696101.2.2_plotdata.df[LOC100696101.2.2_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=LOC100696101.2.2_plotdata.df[LOC100696101.2.2_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=LOC100696101.2.2_plotdata.df[LOC100696101.2.2_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100696101.2.2, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC100696101.2.2 expression (log reads per million)") #edit the axis labels

LOC100696101.2.2_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC100696101.2.2.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100696101.2.2) + mean(salmon_DEgenes.df$QSF_LOC100696101.2.2))
arrdc2_predict.plot <- ggplot()+
  geom_line(data=LOC100696101.2.2_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC100696101.2.2, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 +1), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -1)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC100696101.2.2 (", italic("arrdc2"), ")"))) #edit the axis labels

#LOC100696517.1.1
LOC100696517.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC100696517.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100696517.1.1) + mean(salmon_DEgenes.df$QSF_LOC100696517.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC100696517.1.1_plotdata.df[LOC100696517.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=LOC100696517.1.1_plotdata.df[LOC100696517.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=LOC100696517.1.1_plotdata.df[LOC100696517.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=LOC100696517.1.1_plotdata.df[LOC100696517.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100696517.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC100696517.1.1 expression (log reads per million)") #edit the axis labels
ddit4.plot <- ggplot()+ 
  geom_line(data=LOC100696517.1.1_plotdata.df[LOC100696517.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100696517.1.1, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("LOC100696517.1.1 (", italic("ddit4"),")"))) #edit the axis labels

LOC100696517.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC100696517.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100696517.1.1) + mean(salmon_DEgenes.df$QSF_LOC100696517.1.1))
ddit4_predict.plot <- ggplot()+
  geom_line(data=LOC100696517.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC100696517.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. )*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC100696517.1.1 (", italic("ddit4"), ")"))) #edit the axis labels

#LOC100699058.1.1
LOC100699058.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC100699058.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100699058.1.1) + mean(salmon_DEgenes.df$QSF_LOC100699058.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC100699058.1.1_plotdata.df[LOC100699058.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC100699058.1.1, color=Site), alpha=0.3)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC100699058.1.1 expression (log reads per million)") #edit the axis labels

LOC100699058.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC100699058.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC100699058.1.1) + mean(salmon_DEgenes.df$QSF_LOC100699058.1.1))
plod1a_predict.plot <- ggplot()+
  geom_line(data=LOC100699058.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC100699058.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. )*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC100699058.1.1 (", italic("plod1a"), ")"))) #edit the axis labels

#LOC101166160.1.1
LOC101166160.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC101166160.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC101166160.1.1) + mean(salmon_DEgenes.df$QSF_LOC101166160.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC101166160.1.1_plotdata.df[LOC101166160.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=LOC101166160.1.1_plotdata.df[LOC101166160.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=LOC101166160.1.1_plotdata.df[LOC101166160.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=LOC101166160.1.1_plotdata.df[LOC101166160.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC101166160.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC101166160.1.1 expression (log reads per million)") #edit the axis labels

LOC101166160.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC101166160.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC101166160.1.1) + mean(salmon_DEgenes.df$QSF_LOC101166160.1.1))
LOC129861452_predict.plot <- ggplot()+
  geom_line(data=LOC101166160.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC101166160.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. )*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC101166160.1.1 (", italic("uncharacterized LOC129861452"), ")"))) #edit the axis labels

#LOC101166550.1.1
LOC101166550.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC101166550.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC101166550.1.1) + mean(salmon_DEgenes.df$QSF_LOC101166550.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=LOC101166550.1.1_plotdata.df[LOC101166550.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=LOC101166550.1.1_plotdata.df[LOC101166550.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=LOC101166550.1.1_plotdata.df[LOC101166550.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=LOC101166550.1.1_plotdata.df[LOC101166550.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC101166550.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="LOC101166550.1.1 expression (log reads per million)") #edit the axis labels

LOC101166550.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC101166550.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC101166550.1.1) + mean(salmon_DEgenes.df$QSF_LOC101166550.1.1))
sich73_predict.plot <- ggplot()+
  geom_line(data=LOC101166550.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC101166550.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 +2), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -2)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC101166550.1.1 (", italic("si:ch73-233m11.2"), ")"))) #edit the axis labels

#LOC562935.1.1
LOC562935.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(LOC562935.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC562935.1.1) + mean(salmon_DEgenes.df$QSF_LOC562935.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
hspa8b.plot <- ggplot()+ 
  geom_line(data=LOC562935.1.1_plotdata.df[LOC562935.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_LOC562935.1.1, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("LOC562935.1.1 (", italic("hspa8b"),")")))+ #edit the axis labels
  ylim(6.6, 8.2)

LOC562935.1.1_backcast.df <- tibble(study_predict, z.preds=predict.gam(LOC562935.1.1.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_LOC562935.1.1) + mean(salmon_DEgenes.df$QSF_LOC562935.1.1))
hspa8b_predict.plot <- ggplot()+
  geom_line(data=LOC562935.1.1_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_LOC562935.1.1, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 +2), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -2)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("LOC562935 (", italic("hspa8b"), ")"))) #edit the axis labels

#P2RX4.82.83
P2RX4.82.83_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(P2RX4.82.83.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_P2RX4.82.83) + mean(salmon_DEgenes.df$QSF_P2RX4.82.83), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ube2b.plot <- ggplot()+ 
  geom_line(data=P2RX4.82.83_plotdata.df[P2RX4.82.83_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_P2RX4.82.83, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("P2RX4.82.83 (", italic("ube2b"),")"))) #edit the axis labels

#P4HA1.6.11
P4HA1.6.11_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(P4HA1.6.11.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_P4HA1.6.11) + mean(salmon_DEgenes.df$QSF_P4HA1.6.11), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
p4ha1b.plot <- ggplot()+ 
  geom_line(data=P4HA1.6.11_plotdata.df[P4HA1.6.11_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_P4HA1.6.11, color=Site), alpha=0.3, show.legend=F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("P4HA1.6.11 (", italic("p4ha1b"),")"))) #edit the axis labels

#P4HA2.1.3
P4HA2.1.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(P4HA2.1.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_P4HA2.1.3) + mean(salmon_DEgenes.df$QSF_P4HA2.1.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=P4HA2.1.3_plotdata.df[P4HA2.1.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_P4HA2.1.3, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="P4HA2.1.3 expression (log reads per million)") #edit the axis labels

#PCY2.5.29
PCY2.5.29_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(PCY2.5.29.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_PCY2.5.29) + mean(salmon_DEgenes.df$QSF_PCY2.5.29), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=PCY2.5.29_plotdata.df[PCY2.5.29_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_PCY2.5.29, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="PCY2.5.29 expression (log reads per million)") #edit the axis labels

#RO32.1.2
RO32.1.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(RO32.1.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_RO32.1.2) + mean(salmon_DEgenes.df$QSF_RO32.1.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
hnrnpa3.plot <- ggplot()+ 
  geom_line(data=RO32.1.2_plotdata.df[RO32.1.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_RO32.1.2, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("R032.1.2 (", italic("hnrnpa3"),")"))) #edit the axis labels

#SARNP.8.8
SARNP.8.8_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(SARNP.8.8.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_SARNP.8.8) + mean(salmon_DEgenes.df$QSF_SARNP.8.8), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=SARNP.8.8_plotdata.df[SARNP.8.8_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=SARNP.8.8_plotdata.df[SARNP.8.8_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=SARNP.8.8_plotdata.df[SARNP.8.8_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=SARNP.8.8_plotdata.df[SARNP.8.8_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_SARNP.8.8, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="SARNP.8.8 expression (log reads per million)") #edit the axis labels

#SFRS6.2.5
SFRS6.2.5_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(SFRS6.2.5.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_SFRS6.2.5) + mean(salmon_DEgenes.df$QSF_SFRS6.2.5), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=SFRS6.2.5_plotdata.df[SFRS6.2.5_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=SFRS6.2.5_plotdata.df[SFRS6.2.5_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=SFRS6.2.5_plotdata.df[SFRS6.2.5_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=SFRS6.2.5_plotdata.df[SFRS6.2.5_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_SFRS6.2.5, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="SFRS6.2.5 expression (log reads per million)") #edit the axis labels
sfrs6.plot <- ggplot()+ 
  geom_line(data=SFRS6.2.5_plotdata.df[SFRS6.2.5_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend=F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_SFRS6.2.5, color=Site), alpha=0.3, show.legend=F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("SFRS6.2.5 (", italic("sfrs6"),")"))) #edit the axis labels

#SRSF2.6.8
SRSF2.6.8_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(SRSF2.6.8.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_SRSF2.6.8) + mean(salmon_DEgenes.df$QSF_SRSF2.6.8), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
srsf2.plot <- ggplot()+ 
  geom_line(data=SRSF2.6.8_plotdata.df[SRSF2.6.8_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_SRSF2.6.8, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("SRSF2.6.8 (", italic("srsf2"),")"))) #edit the axis labels

#TCPG.1.3
TCPG.1.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TCPG.1.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TCPG.1.3) + mean(salmon_DEgenes.df$QSF_TCPG.1.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
TCPG.plot <- ggplot()+ 
  geom_line(data=TCPG.1.3_plotdata.df[TCPG.1.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TCPG.1.3, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="TCPG.1.3 expression (log reads per million)") #edit the axis labels
stream.leg <- get_legend(TCPG.plot)

#THOC4.1.2
THOC4.1.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(THOC4.1.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_THOC4.1.2) + mean(salmon_DEgenes.df$QSF_THOC4.1.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=THOC4.1.2_plotdata.df[THOC4.1.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=THOC4.1.2_plotdata.df[THOC4.1.2_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=THOC4.1.2_plotdata.df[THOC4.1.2_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=THOC4.1.2_plotdata.df[THOC4.1.2_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_THOC4.1.2, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="THOC4.1.2 expression (log reads per million)") #edit the axis labels

#TIA1.4.5
TIA1.4.5_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TIA1.4.5.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TIA1.4.5) + mean(salmon_DEgenes.df$QSF_TIA1.4.5), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
tia1.plot <- ggplot()+ 
  geom_line(data=TIA1.4.5_plotdata.df[TIA1.4.5_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TIA1.4.5, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("TIA1.4.5 (", italic("tia1"),")"))) #edit the axis labels

#TRA2B.3.7
TRA2B.3.7_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TRA2B.3.7.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TRA2B.3.7) + mean(salmon_DEgenes.df$QSF_TRA2B.3.7), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting

ggplot()+ 
  geom_line(data=TRA2B.3.7_plotdata.df[TRA2B.3.7_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=TRA2B.3.7_plotdata.df[TRA2B.3.7_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=TRA2B.3.7_plotdata.df[TRA2B.3.7_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=TRA2B.3.7_plotdata.df[TRA2B.3.7_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TRA2B.3.7, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="TRA2B.3.7 expression (log reads per million)") #edit the axis labels

TRA2B.3.7_backcast.df <- tibble(study_predict, z.preds=predict.gam(TRA2B.3.7.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_TRA2B.3.7) + mean(salmon_DEgenes.df$QSF_TRA2B.3.7))
tra2b_predict.plot <- ggplot()+
  geom_line(data=TRA2B.3.7_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_TRA2B.3.7, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 0.5), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -0.5)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("TRA2B.3.7 (", italic("tra2b"), ")"))) #edit the axis labels

#TRAD1.20.35
TRAD1.20.35_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TRAD1.20.35.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TRAD1.20.35) + mean(salmon_DEgenes.df$QSF_TRAD1.20.35), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=TRAD1.20.35_plotdata.df[TRAD1.20.35_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=TRAD1.20.35_plotdata.df[TRAD1.20.35_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=TRAD1.20.35_plotdata.df[TRAD1.20.35_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=TRAD1.20.35_plotdata.df[TRAD1.20.35_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TRAD1.20.35, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="TRAD1.20.35 expression (log reads per million)") #edit the axis labels

#TRI39.3.8
TRI39.3.8_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TRI39.3.8.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TRI39.3.8) + mean(salmon_DEgenes.df$QSF_TRI39.3.8), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
trim25.plot <- ggplot()+ 
  geom_line(data=TRI39.3.8_plotdata.df[TRI39.3.8_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TRI39.3.8, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("TR139.3.8 (", italic("E3 ubiquitin/ISG15 ligase TRIM25-like"),")")))+ #edit the axis labels
  ylim(4, 6.5)

#TSC22D3.2.2
TSC22D3.2.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TSC22D3.2.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TSC22D3.2.2) + mean(salmon_DEgenes.df$QSF_TSC22D3.2.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=TSC22D3.2.2_plotdata.df[TSC22D3.2.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=TSC22D3.2.2_plotdata.df[TSC22D3.2.2_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=TSC22D3.2.2_plotdata.df[TSC22D3.2.2_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=TSC22D3.2.2_plotdata.df[TSC22D3.2.2_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TSC22D3.2.2, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="TSC22D3.2.2 expression (log reads per million)") #edit the axis labels

#UGT1A5.1.1
UGT1A5.1.1_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(UGT1A5.1.1.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_UGT1A5.1.1) + mean(salmon_DEgenes.df$QSF_UGT1A5.1.1), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=UGT1A5.1.1_plotdata.df[UGT1A5.1.1_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=UGT1A5.1.1_plotdata.df[UGT1A5.1.1_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=UGT1A5.1.1_plotdata.df[UGT1A5.1.1_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=UGT1A5.1.1_plotdata.df[UGT1A5.1.1_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_UGT1A5.1.1, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="UGT1A5.1.1 expression (log reads per million)") #edit the axis labels

#USP9X.2.4
USP9X.2.4_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(USP9X.2.4.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_USP9X.2.4) + mean(salmon_DEgenes.df$QSF_USP9X.2.4), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=USP9X.2.4_plotdata.df[USP9X.2.4_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=USP9X.2.4_plotdata.df[USP9X.2.4_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=USP9X.2.4_plotdata.df[USP9X.2.4_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=USP9X.2.4_plotdata.df[USP9X.2.4_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_USP9X.2.4, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="USP9X.2.4 expression (log reads per million)") #edit the axis labels
usp9.plot <- ggplot()+ 
  geom_line(data=USP9X.2.4_plotdata.df[USP9X.2.4_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_USP9X.2.4, color=Site), alpha=0.3, show.legend = F)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x=" ", y=" ", title=expression(paste("A - USP9X.2.4 (", italic("usp9"),")"))) #edit the axis labels

USP9X.2.4_backcast.df <- tibble(study_predict, z.preds=predict.gam(USP9X.2.4.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_USP9X.2.4) + mean(salmon_DEgenes.df$QSF_USP9X.2.4))
usp9_predict.plot <- ggplot()+
  geom_line(data=USP9X.2.4_backcast.df, aes(x=Date, y=preds, color=Site), show.legend = F)+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_USP9X.2.4, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 +2), alpha=0.5, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (. -2)*4, name="Temperature"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="Expression (Counts Per Million)", title=expression(paste("B - USP9X.2.4 (", italic("usp9"), ")"))) #edit the axis labels


#ZGC_162825.1.2
ZGC_162825.1.2_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(ZGC_162825.1.2.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_ZGC_162825.1.2) + mean(salmon_DEgenes.df$QSF_ZGC_162825.1.2), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
ggplot()+ 
  geom_line(data=ZGC_162825.1.2_plotdata.df[ZGC_162825.1.2_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first heatwave
  geom_line(data=ZGC_162825.1.2_plotdata.df[ZGC_162825.1.2_plotdata.df$Date=="2022-08-08",], aes(x=DailyAvgTemp, y=preds, color="Peak Heatwave 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second heatwave
  geom_line(data=ZGC_162825.1.2_plotdata.df[ZGC_162825.1.2_plotdata.df$Date=="2022-07-30",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 1"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the first recovery period
  geom_line(data=ZGC_162825.1.2_plotdata.df[ZGC_162825.1.2_plotdata.df$Date=="2022-08-18",], aes(x=DailyAvgTemp, y=preds, color="Recovery Period 2"))+ #add the model line for the relationship between temperature and the module eigengene value during the peak of the second recovery period
  scale_color_manual(name="Date", values=c("Peak Heatwave 1"="firebrick1", "Peak Heatwave 2"="firebrick4", "Recovery Period 1"="slateblue1", "Recovery Period 2"="slateblue4"))+ #set the colors to use for each of the modeled lines
  ggnewscale::new_scale_colour()+ #set a break so we can use a new discrete color scale
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_ZGC_162825.1.2, color=Site), alpha=0.3)+ #add the actual eigengene values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Daily Average Temperature (ºC)", y="ZGC_162825.1.2 expression (log reads per million)") #edit the axis labels

#Plots for manuscript
#AN32B.3.8
AN32B.3.8_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(AN32B.3.8.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_AN32B.3.8) + mean(salmon_DEgenes.df$QSF_AN32B.3.8), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
cirbp.plot <- ggplot()+ 
  geom_line(data=AN32B.3.8_plotdata.df[AN32B.3.8_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_AN32B.3.8, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="", title=expression(paste("A - AN32B.3.8 (", italic("cirbp"), ")"))) #edit the axis labels

AN32B.3.8_backcast.df <- tibble(study_predict, z.preds=predict.gam(AN32B.3.8.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_AN32B.3.8) + mean(salmon_DEgenes.df$QSF_AN32B.3.8))
cirbp_predict1.plot <- ggplot()+
  #geom_line(data=AN32B.3.8_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_AN32B.3.8, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_point(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site), size=4)+
  geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4, color=Site), alpha=1, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ . *4, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="")
  
cirbp_predict2.plot <- ggplot()+
  geom_line(data=AN32B.3.8_backcast.df, aes(x=Date, y=preds, color=Site), linewidth=2)+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_AN32B.3.8, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  #geom_point(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site), size=4)+
  #geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_AN32B.3.8, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4), alpha=0.6, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ . *4, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="")
  
#HSP47.2.3
HSP47.2.3_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(HSP47.2.3.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP47.2.3) + mean(salmon_DEgenes.df$QSF_HSP47.2.3), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
serpinh1.plot <- ggplot()+ 
  geom_line(data=HSP47.2.3_plotdata.df[HSP47.2.3_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_HSP47.2.3, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="", title=expression(paste("B - HSP47.2.3 (", italic("serpinh1b"), ")"))) #edit the axis labels

HSP47.2.3_backcast.df <- tibble(study_predict, z.preds=predict.gam(HSP47.2.3.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_HSP47.2.3) + mean(salmon_DEgenes.df$QSF_HSP47.2.3))
serpinh1_predict1.plot <- ggplot()+
  #geom_line(data=HSP47.2.3_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HSP47.2.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_point(data=express_means, aes(x=Date, y=QSF_HSP47.2.3, color=Site), size=4)+
  geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_HSP47.2.3, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/3 - 4, color=Site), alpha=1, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (.+4) *3, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="")

serpinh1_predict2.plot <- ggplot()+
  geom_line(data=HSP47.2.3_backcast.df, aes(x=Date, y=preds, color=Site), linewidth=2)+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_HSP47.2.3, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  #geom_point(data=express_means, aes(x=Date, y=QSF_HSP47.2.3, color=Site), size=4)+
  #geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_HSP47.2.3, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/3 - 4), alpha=0.6, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (.+4) *3, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="")

#TRA2B.3.7
TRA2B.3.7_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(TRA2B.3.7.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_TRA2B.3.7) + mean(salmon_DEgenes.df$QSF_TRA2B.3.7), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
tra2b.plot <- ggplot()+ 
  geom_line(data=TRA2B.3.7_plotdata.df[TRA2B.3.7_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_TRA2B.3.7, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="", title=expression(paste("C - TRA2B.3.7 (", italic("tra2b"), ")"))) #edit the axis labels

TRA2B.3.7_backcast.df <- tibble(study_predict, z.preds=predict.gam(TRA2B.3.7.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_TRA2B.3.7) + mean(salmon_DEgenes.df$QSF_TRA2B.3.7))
tra2b_predict1.plot <- ggplot()+
  #geom_line(data=TRA2B.3.7_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_TRA2B.3.7, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_point(data=express_means, aes(x=Date, y=QSF_TRA2B.3.7, color=Site), size=4)+
  geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_TRA2B.3.7, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4, color=Site), alpha=1, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (.) *4, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="")

tra2b_predict2.plot <- ggplot()+
  geom_line(data=TRA2B.3.7_backcast.df, aes(x=Date, y=preds, color=Site), linewidth=2)+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_TRA2B.3.7, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  #geom_point(data=express_means, aes(x=Date, y=QSF_TRA2B.3.7, color=Site), size=4)+
  #geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_TRA2B.3.7, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4), alpha=0.6, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (.) *4, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="", y="")

#USP9X.2.4
USP9X.2.4_plotdata.df <- tibble(z.DailyAvgTemp=newdata3.df$z.DailyAvgTemp, DailyAvgTemp=newdata3.df$DailyAvgTemp, z.preds=predict.gam(USP9X.2.4.gam, newdata=newdata3.df), preds=z.preds * sd(salmon_DEgenes.df$QSF_USP9X.2.4) + mean(salmon_DEgenes.df$QSF_USP9X.2.4), Date2=newdata3.df$Date2, Date=as.Date(Date2, origin="2022-07-15"), Site=newdata3.df$Site, Length=newdata3.df$Length, z.Length=newdata3.df$z.Length) #generate predicted values of DE genes expression based on the best GAM and combine them with all of the explanatory variable needed for plotting
usp9.plot <- ggplot()+ 
  geom_line(data=USP9X.2.4_plotdata.df[USP9X.2.4_plotdata.df$Date=="2022-07-22",], aes(x=DailyAvgTemp, y=preds, color=Site), show.legend = F)+ #add the model line for the relationship between temperature and the DE gene value during the peak of the first heatwave
  geom_point(data=salmon_DEgenes.df, aes(x=DailyAvgTemp, y=QSF_USP9X.2.4, color=Site), alpha=0.3, show.legend = F)+ #add the actual gene expression values for each individual
  theme_classic()+ #change the plot theme to make it look neater
  #facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Temperature (ºC)", y="", title=expression(paste("D - USP9X.2.4 (", italic("usp9"), ")"))) #edit the axis labels

USP9X.2.4_backcast.df <- tibble(study_predict, z.preds=predict.gam(USP9X.2.4.gam, study_predict), preds=z.preds * sd(salmon_DEgenes.df$QSF_USP9X.2.4) + mean(salmon_DEgenes.df$QSF_USP9X.2.4))
usp9_predict1.plot <- ggplot()+
  #geom_line(data=USP9X.2.4_backcast.df, aes(x=Date, y=preds, color=Site))+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_USP9X.2.4, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  geom_point(data=express_means, aes(x=Date, y=QSF_USP9X.2.4, color=Site), size=4)+
  geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_USP9X.2.4, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 2, color=Site), alpha=1, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (.-2) *4, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+ #facet the plot by Site so the differences can be seen more easily
  labs(x="Date", y="")

usp9_predict2.plot <- ggplot()+
  geom_line(data=USP9X.2.4_backcast.df, aes(x=Date, y=preds, color=Site), linewidth=2)+
  geom_point(data=salmon_DEgenes.df, aes(x=Date, y=QSF_USP9X.2.4, color=Site), alpha=0.3, show.legend = F)+#add the actual gene expression values for each individual
  #geom_point(data=express_means, aes(x=Date, y=QSF_USP9X.2.4, color=Site), size=4)+
  #geom_arrow_chain(data=express_means, aes(x=Date, y=QSF_USP9X.2.4, color=Site))+
  geom_line(data=daily_temps.df, aes(x=Date, y=Temp.avg/4 + 2), alpha=0.6, linetype=2)+
  scale_y_continuous(sec.axis=sec_axis(transform = ~ (.-2) *4, name="Temperature (ºC)"))+
  xlim(as.Date("2022-07-15"), as.Date("2022-08-19"))+
  theme_classic()+ #change the plot theme to make it look neater
  facet_wrap(~Site, nrow=1)+
  labs(x="Date", y="")

#Figure 6
Fig6_v1.plots <- ggarrange(cirbp.plot,cirbp_predict1.plot, serpinh1.plot, serpinh1_predict1.plot, tra2b.plot, tra2b_predict1.plot,usp9.plot, usp9_predict1.plot, widths=c(1,4), heights=c(1,1,1,1.1), ncol=2, nrow=4, legend ="right", common.legend = T)
Fig6_v1.plots <- annotate_figure(Fig6_v1.plots, left=text_grob("Expression (log reads per million)", size=14, rot=90))
ggsave("Figure6.jpeg", Fig6_v1.plots, width=12, height=9, units="in", dpi=300)
Fig6_v2.plots <- ggarrange(cirbp.plot,cirbp_predict2.plot, serpinh1.plot, serpinh1_predict2.plot, tra2b.plot, tra2b_predict2.plot,usp9.plot, usp9_predict2.plot, widths=c(1,4), heights=c(1,1,1,1.1), ncol=2, nrow=4, legend ="right", common.legend = T)
Fig6_v2.plots <- annotate_figure(Fig6_v2.plots, left=text_grob("Expression (log reads per million)", size=14, rot=90))
ggsave("Figure6v2.jpeg", Fig6_v2.plots, width=12, height=9, units="in", dpi=300)

#Figure S1
DEgene.plots <- ggarrange(dnaja.plot, hspa4l.plot, gatd3.plot, hnrnpa3.plot, hsp90ab1.2.plot, stream.leg, cdc37.plot, serpinh1.1.plot, cirbpb1.plot, sfrs6.plot, p4ha1b.plot, NULL, hsp70.plot, hsp90ab1.1.plot, srsf2.plot, usp9.plot, ddit4.plot, ube2b.plot, trim25.plot, tia1.plot, khdrbs1.plot, cirbpb.2.plot, serpinh1b.2.plot, hspa8b.plot, nrow=4, ncol=6, common.legend=F) #list the DE gene graphs to be added to Figure 5, organized into 4 rows & 6 columns
DEgene.plots <- annotate_figure(DEgene.plots, left=text_grob("Expression (log reads per million)", size=18, rot=90), bottom=text_grob("24-h Average Temperature (ºC)", size=18)) #add shared axis labels 
ggsave("FigureS1.jpeg", DEgene.plots, width=18, height=12, units="in", dpi=300) #save a jpeg of Figure 5

#Which WGCNA modules do DE genes belong to?
load("WGCNA.RData") #load in the wgcna results 
salmon_signed_network_cH0.995_ds2_mcH0.15_mms40 <- blockwiseModules(express_data, maxBlockSize = 35000, loadTom=salmon_signed_TOM, power=8, networkType = "signed", minModuleSize = 40, mergeCutHeight = 0.15, deepSplit = 2, detectCutHeight = 0.995) #create the modules from the WGCNA dendrogram
DEgenes_colors.df <- data.frame(colors=as.factor(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors[names(salmon_signed_network_cH0.995_ds2_mcH0.15_mms40$colors)%in% heatwave_DEgenes.vec])) #match up genes in the WGCNA to the genes listed in the DE gene vector, add the module color names to the DE gene vector
DEgenes_colors.df %>% group_by(colors) %>% summarise(length(colors)) #show how many DE genes were assigned to each module

#DE Genes tree & GAM heatmap
GAM_results.df <- read.table("BKT_DEgenes_GAMresults.txt", header = T)
DEgenes.tree <- hclust(dist(cor(salmon_DEgenes.df[,17:59], salmon_DEgenes.df[,17:59])))
GAM_results.df$Transcript <- reorder(GAM_results.df$Transcript, -c(31, 10, 41, 39, 21, 22, 9, 33, 11, 12, 2, 3, 18, 17, 42, 20, 30, 15, 25, 27, 5, 43, 32, 19, 6, 1, 4, 7, 13, 34, 16, 36, 37, 38, 35, 29, 40, 8, 26, 23, 14, 28, 24))
GAM_results.df$GeneName <- reorder(GAM_results.df$GeneName, -c(31, 10, 41, 39, 21, 22, 9, 33, 11, 12, 2, 3, 18, 17, 42, 20, 30, 15, 25, 27, 5, 43, 32, 19, 6, 1, 4, 7, 13, 34, 16, 36, 37, 38, 35, 29, 40, 8, 26, 23, 14, 28, 24))
tree_labels <- data.frame(DEgenes.tree$labels, newlabel=label_pad(str_sub(DEgenes.tree$labels, start=5L, end=-1L)))
DEgenes_tree.ggplot <- ggtree(DEgenes.tree, layout = "rectangular", ladderize = T) %<+% tree_labels + xlim(NA,5)+
  theme(plot.margin = unit(c(0, -0.3, 0.5, 0.5), "cm"))+
  geom_hilight(mapping=aes(subset=node %in% c(45,46), fill=c("Up-Regulated", "Down-Regulated")), alpha=0.5, align="right")+
  scale_fill_manual(values=c("Up-Regulated"="red", "Down-Regulated"="blue"), guide=guide_legend(position="bottom", title="Expression Profile During Heatwave", theme(legend.title.position = "top")))+
  geom_tiplab(aes(label=newlabel), align = T, family="mono", fontface=2)

GAM_results_long.df <- GAM_results.df %>% as_tibble() %>% gather(key="Variable", value="p-val", Temp:Temp.Date, factor_key = T)
DEgenes_GAMheat.ggplot <- ggplot(data=GAM_results_long.df, aes(x=Variable, y=Transcript:GeneName, fill=-log(`p-val`+ 0.0000000000000001)))+
  geom_tile()+ #create a colored tile geom
  scale_fill_gradient(low="white", high="red", limits=c(0, 40), guide=guide_colorbar(title="-log(P-Value)"))+ #set the high p-values to white and low p-values to red 
  ggnewscale::new_scale_fill()+
  geom_tile(aes(x="GAMM R^2", y=Transcript:GeneName, fill=Rsquared))+
  scale_fill_gradient(low="white", high="blue", limits=c(0,1), guide=guide_colorbar(title="R^2"))+ #set the low R^2 to white and high R^2 to blue 
  theme_classic(base_size = 14)+ #increase the base font size and change the plot theme 
  theme(axis.text.y=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) #remove the module names and y-axis (it will be aligned with the dendrogram)

DEgenes_GAMMplots.plot <- grid.arrange(DEgenes_tree.ggplot, DEgenes_GAMheat.ggplot, ncol=2, widths=c(0.9,1.6))
ggsave("Figure_4.jpeg", plot=DEgenes_GAMMplots.plot, width=12, height=9, units = "in", dpi=300)
  

  
