#read in libraries
library(dplyr)
library(tidyr)
library(purrr)
#BiocManager::install("edgeR", force=T)
library(edgeR)
library(stringr)
library(ggplot2)
library(ade4)

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
salmon_counts_cpm <- cpm(salmon_counts[,-1]) #standardize the raw gene counts to counts per million reads (ignore first column with gene names)
occurrence_filter_test <- data.frame(N=1:115, genes0=0, genes05=0, genes1=0, genes2=0) #create an empty data frame to test how different minimum transcript count cutoffs and sample occurence cutoffs would affect the final dataset 
for(i in 1:115){
  occurrence_filter_test$genes05[i] <- sum(rowSums(salmon_counts_cpm > 0.5) >= i) #calculate how many transcripts have at least 0.5 reads per million across N samples
  occurrence_filter_test$genes1[i] <- sum(rowSums(salmon_counts_cpm > 1) >= i) #calculate how many transcripts have at least 1 read per million across N samples
  occurrence_filter_test$genes2[i] <- sum(rowSums(salmon_counts_cpm > 2) >= i) #calculate how many transcripts have at least 2 reads per million across N samples
  occurrence_filter_test$genes0[i] <- sum(rowSums(salmon_counts_cpm > 0) >= i) #calculate how many transcripts have at least 1 read across N samples
}
ggplot(data=occurrence_filter_test)+ #plot results of occurrence filter tests
  geom_point(aes(x=N, y=genes0), pch=1, color="blue")+ #points and line for any transcripts
  geom_line(aes(x=N, y=genes0), color="blue")+
  geom_point(aes(x=N, y=genes05), pch=1, color="purple")+ #points and line for transcritps with cpm > 0.5
  geom_line(aes(x=N, y=genes05), color="purple")+ 
  geom_point(aes(x=N, y=genes1), pch=1, color="red")+ #points and line for transcripts with cpm > 1
  geom_line(aes(x=N, y=genes1), color="red")+
  geom_point(aes(x=N, y=genes2), pch=1, color="orange")+ #points and line for transctipts with cpm > 2
  geom_line(aes(x=N, y=genes2), color="orange")+
  theme_classic()+ #change the plot theme
  labs(x="Sample Threshold", y="Transcripts Retained") #add axis labels

cpm05 <- rowSums(salmon_counts_cpm > 1) >= 60 #identify the genes with at least 0.5 reads per million across at least 2 samples
salmon_counts_cpm05 <- cbind(salmon_counts$Name, salmon_counts_cpm)[cpm05,] #get rid of genes in the count per million data frame with very low expression levels that do not meet the above criteria
salmon_counts_05 <- cbind(salmon_counts)[cpm05,] #get rid of genes in the raw read count data frame with very low expression levels that do not meet the above criteria
salmon_counts.dge <- DGEList(salmon_counts_05) #convert to Digital Gene Expression list format
plotMDS(salmon_counts.dge) #create NMDS plot of salmon count data
salmon_counts_filtered.dge <- salmon_counts.dge[,-55] #remove outlier sample EB044
plotMDS(salmon_counts_filtered.dge) #replot NMDS without EB044 

trout_enviro_data.df <- read.delim("~/Desktop/Brook_Trout_Heatwave/Brook_Trout_Heatwave_SampleSelect/trout_enviro_data.txt", sep="\t") #read in individual trout environmental variable data frame
trout_enviro_data.df$Sample_ID <- as.factor(trout_enviro_data.df$Sample_ID) #set sample id to factor class
trout_enviro_data.df$Date <- as.Date(trout_enviro_data.df$Date) #set sample dates to Date class
trout_enviro_data.df$Site <- as.factor(trout_enviro_data.df$Site) #set sample site to factor class
trout_sexes.df <- read.delim("~/Desktop/BrookTrout_ReadCounts/BrookTrout_ReadCounts/BrookTrout_Sexes.txt", header = T) #read in data frame with fish sex information
trout_enviro_data_filtered.df <- right_join(trout_enviro_data.df, data.frame(Sample_ID=rownames(salmon_counts.dge$samples)), by="Sample_ID") #combine environmental information with read counts matched up by sample ID
trout_enviro_data_filtered.df <- left_join(trout_enviro_data_filtered.df, trout_sexes.df, by="Sample_ID") #add sex data to environmental/read count data frame, matched up by sample ID
trout_size.df <- read.delim("~/Desktop/BrookTrout_ReadCounts/BrookTrout_ReadCounts/BrookTrout_LengthWeight.txt", header=T) #read in data frame with length and weight data
trout_enviro_data_filtered.df <- left_join(trout_enviro_data_filtered.df, trout_size.df, join_by(Sample_ID==Sample.ID), relationship = "one-to-one") #add length and weight data to environmental/read count/sex data frame matched up by sample ID
trout_enviro_data_filtered.df <- trout_enviro_data_filtered.df[-55,] #remove outlier sample EB044

salmon_logcounts <- cpm(salmon_counts_filtered.dge, log = T) #convert count data to log scale
salmon.pca <- dudi.pca(salmon_logcounts, scannf=F, nf=5, center=T, scale=T) #run PCA on the transcript data frame
salmon.pca$eig[1]/sum(salmon.pca$eig) #calculate the variance explained by PC1
salmon.pca$eig[2]/sum(salmon.pca$eig) #calcualte the variance explained by PC2
salmon.pca$eig[3]/sum(salmon.pca$eig) #calculate the variance explained by PC3
salmon.pca$eig[4]/sum(salmon.pca$eig) #calculate the variance explained by PC4

ggplot(data=salmon.pca$c1, aes(x=CS1, y=CS2, color=trout_enviro_data_filtered.df$Site))+ #plot the first two principal component axes, colored by sample site
  geom_point()+ #add sample points
  stat_ellipse()+ #add 80% CI ellipses around each sample site group
  theme_classic()+ #change the plot theme
  labs(x="PC1 (84.7%)", y="PC2 (1.1%)", color="Sample Site") #set axis and legend titles

ggplot(data=salmon.pca$c1, aes(x=CS3, y=CS4, color=trout_enviro_data_filtered.df$Site))+ #plot the third and fourth principal component axes, colored by sample site
  geom_point()+ #add sample points
  stat_ellipse()+ #add 80% CI ellipses around each sample site group
  theme_classic()+ #change the plot theme
  labs(x="PC3 (0.9%)", y="PC2 (0.7%)", color="Sample Site") #set axis and legend titles

trout_enviro_data.df <- left_join(data.frame(Sample_ID=(colnames(salmon_counts[,-c(1,55)]))), trout_enviro_data_filtered.df, by="Sample_ID") #add sample ID names to the environmental data data frame
salmon_logcounts_t <- t(salmon_logcounts) #transpose the log-transformed read counts data frame (so samples are in rows, transcripts in columns)
colnames(salmon_logcounts_t) <- salmon_counts.dge$genes$Name #add transcript names as column names to the read count data frame
salmon_logcounts_environ.df <- right_join(trout_enviro_data.df, data.frame(Sample_ID=rownames(salmon_logcounts_t), salmon_logcounts_t), by="Sample_ID") #combine the read count and environmental variable data frames matched up by sample
write.table(salmon_logcounts_environ.df, file="~/Desktop/BrookTrout_ReadCounts/BrookTrout_ReadCounts/salmon_logcounts_environ.txt") #write data frame as tab-delimited text file


