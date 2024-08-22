#devtools::install_github("thierrygosselin/radiator")
library(radiator)
#BiocManager::install("SNPRelate")
library(vegan)
library(ade4)
library(adegenet)
library(ggplot2)
library(hierfstat)
library(data.table)
library(MuMIn)
library(rdacca.hp)

#snps.data <- read_vcf("populations.snps.vcf") #read in vcf file from Stacks
#snps_filtered.data <- filter_rad(snps.data) #run the radiator filters, see manuscript for values
snps_filtered.data <- read_rad("radiator_data_20240506@0951.rad") #read in filtered SNP dataset
snps_filtered.genind <- genomic_converter(snps_filtered.data, output = "genind") #convert filtered SNP dataset to genind format
snps_filtered.genind.sum <- summary(snps_filtered.genind$genind) #create summary object of the genind (mainly to look at NA percentage)
snps_filtered.table <- tab(snps_filtered.genind$genind, freq=T, NA.method="mean") #convert filtered SNP genind to tabular format, replace missing genotypes with the average genotype for that population
pca.snps <- dudi.pca(snps_filtered.table, center=T, scale=F, scannf=F, nf=30) #run PCA on the SNP table
s.label(pca.snps$li) #plot a quick view of the PCA 

pca.snps.df <- pca.snps$li #create a data frame based on the principal component scores
pca.snps.df$Site <- as.factor(c(rep("Big_Poe", 10), rep("East_Branch", 11), rep("Shavers_Creek", 16), rep("Standing_Stone", 10))) #add populations to the principal components data frame

RADseq_pca.plot <- ggplot(data=pca.snps.df, aes(x=Axis1, y=Axis2, color=Site, fill=Site))+ #create a ggplot with PCs 1 & 2 along the x & y axes, and color and fill geoms by population
  geom_point()+ #add points to the SNP PCA plot
  stat_ellipse(geom="polygon", alpha=0.5, level=0.8, color=NA)+ #add transparent filled 80% CI ellipses with no outlines to the SNP PCA plot
  theme_classic()+ #change the plot theme to make it look cleaner
  geom_vline(xintercept = 0)+ #add a bolder line along the y axis to make it stand out more
  geom_hline(yintercept = 0)+ #add a bolder line along the x axis to make it stand out more
  xlab("PC1 (10.3%)")+ #add x axis label
  ylab("PC2 (7.4%)") #add y axis label
ggsave("RADseq_pca.jpeg", plot=RADseq_pca.plot, width=8, height=8, units = "in", dpi=300) #save a jpeg of the SNP PCA plot

ggplot(data=pca.snps.df, aes(x=Axis3, y=Axis4, color=Site, fill=Site))+ #create a ggplot with PCs 3 & 4 along the x & y axes, and color and fill geoms by population
  geom_point()+ #add points to the SNP PCA plot
  stat_ellipse(geom="polygon", alpha=0.5, level=0.8, color=NA)+ #add transparent filled 80% CI ellipses with no outlines to the SNP PCA plot
  theme_classic()+ #change the plot theme to make it look cleaner
  geom_vline(xintercept = 0)+ #add a bolder line along the y axis to make it stand out more
  geom_hline(yintercept = 0) #add a bolder line along the x axis to make it stand out more

pop_snps.df <- genind2hierfstat(snps_filtered.genind$genind, pop=as.factor(c(rep("Big_Poe", 10), rep("East_Branch", 11), rep("Shavers_Creek", 16), rep("Standing_Stone", 10)))) #convert the filtered SNP genind object to hierfstat format, make sure population names get added correctly
snps.Fst <- genet.dist(pop_snps.df, method="Fst") #calculate pairwise Fst between all populations
snps.bootvc <- boot.vc(levels=as.factor(pca.snps.df$Site), loci=pop_snps.df[,-1]) #calculate confidence intervals for F-statistic and heterozygosity estimates
snps.bootvc$ci #print overall F-statistic and heterozygosity confidence intervals

snps.genpop <- genind2genpop(snps_filtered.genind$genind, pop=as.factor(c(rep("Big_Poe", 10), rep("East_Branch", 11), rep("Shavers_Creek", 16), rep("Standing_Stone", 10)))) #convert the filtered SNP genind object to genpop format
snps_hclust <- hclust(dist.genpop(snps.genpop, method=1), method="average") #calculate UPGMA tree for populations
plot(snps_hclust) #view a quick plot the the UPGMA tree

salmon_data <- fread("salmon_logcounts_environ.txt", drop=1) #read in dataset with environmental and expression data
rad_salmon_data <- salmon_data[salmon_data$Sample_ID %in% rownames(pca.snps.df), ] #filter the transcriptomic data to only retain samples that also have SNP genomic data
enviro_data <- rad_salmon_data[,1:16] #create a separate data frame with just environmental data
express_data <- data.frame(rad_salmon_data$QSF_AN32B.3.8, rad_salmon_data$QSF_CDC37.1.3, rad_salmon_data$QSF_CIRBP.24.24, rad_salmon_data$QSF_CIRBP.7.24, rad_salmon_data$QSF_CSF3R.3.3, rad_salmon_data$QSF_DH12A.2.2, rad_salmon_data$QSF_ES1.1.3, rad_salmon_data$QSF_HMGB2.9.11, rad_salmon_data$QSF_HS90A.1.1, rad_salmon_data$QSF_HS90B.1.1, rad_salmon_data$QSF_HSP47.2.3, rad_salmon_data$QSF_HSP47.3.3, rad_salmon_data$QSF_HSP70B.1.1, rad_salmon_data$QSF_HSPA4L.1.1, rad_salmon_data$QSF_KHDR1.2.2, rad_salmon_data$QSF_LOC100136615.4.4, rad_salmon_data$QSF_LOC100194703.44.106, rad_salmon_data$QSF_LOC100691603.1.1, rad_salmon_data$QSF_LOC100696101.2.2, rad_salmon_data$QSF_LOC100696517.1.1, rad_salmon_data$QSF_LOC100699058.1.1, rad_salmon_data$QSF_LOC101166160.1.1, rad_salmon_data$QSF_LOC101166550.1.1, rad_salmon_data$QSF_LOC562935.1.1, rad_salmon_data$QSF_P2RX4.82.83, rad_salmon_data$QSF_P4HA1.6.11, rad_salmon_data$QSF_P4HA2.1.3, rad_salmon_data$QSF_PCY2.5.29, rad_salmon_data$QSF_RO32.1.2, rad_salmon_data$QSF_SARNP.8.8, rad_salmon_data$QSF_SFRS6.2.5, rad_salmon_data$QSF_SRSF2.6.8, rad_salmon_data$QSF_TCPG.1.3, rad_salmon_data$QSF_THOC4.1.2, rad_salmon_data$QSF_TIA1.4.5, rad_salmon_data$QSF_TRA2B.3.7, rad_salmon_data$QSF_TRAD1.20.35, rad_salmon_data$QSF_TRI39.3.8, rad_salmon_data$QSF_TSC22D3.2.2, rad_salmon_data$QSF_UGT1A5.1.1, rad_salmon_data$QSF_USP9X.2.4, rad_salmon_data$QSF_VRK3.145.924, rad_salmon_data$QSF_ZGC_162825.1.2) #create a separate data frame with just expression data for the differentially expressed genes identified by ImpulseDE2
rownames(express_data) <- rad_salmon_data$Sample_ID #name the expression data rows with sample IDs so future analyses of the expression data can be aligned to environmental data
#add principal components from the SNP PCA to the environmental data frame
enviro_data$PC1 <- pca.snps$li$Axis1
enviro_data$PC2 <- pca.snps$li$Axis2
enviro_data$PC3 <- pca.snps$li$Axis3
enviro_data$PC4 <- pca.snps$li$Axis4
enviro_data$PC5 <- pca.snps$li$Axis5
enviro_data$PC6 <- pca.snps$li$Axis6

start_date <- as.Date("2022-07-15") #set a start date a few days before the beginning of the sampling period
enviro_data$Date2 <- difftime(as.Date(enviro_data$Date), start_date, units = "days") #create a second Date object calculating the number of days from the beginning of the first heatwave

norm_express_data <- stdize(express_data) #standardize the expression data so all transcripts have a mean of 0 and standard deviation of 1
expression_full.rda <- rda(norm_express_data ~ Date2 + DailyAvgTemp + as.factor(Sex) + as.factor(Site) + PC1 + PC2 + PC3 +PC4 +PC5 + PC6, data=enviro_data) #run an RDA on the expression data with date, temperature, sample stream, and the first 6 SNP principal components 
anova.cca(expression_full.rda) #run an anova to check the significance of the variance explained by the constrained factors of the RDA
vif.cca(expression_full.rda) #check variance inflation factors of the contrained variables of the RDA

enviro_data_test.df <- data.frame(Date=enviro_data$Date2, DailyAvgTemp=enviro_data$DailyAvgTemp, Sex=as.factor(enviro_data$Sex), Site=as.factor(enviro_data$Site), PC1=enviro_data$PC1, PC2=enviro_data$PC2, PC3=enviro_data$PC3, PC4=enviro_data$PC4, PC5=enviro_data$PC5, PC6=enviro_data$PC6) #create an input file with the explanatory variables for our RDA model
varpart.allvars <- permu.hp(dv=norm_express_data, iv=enviro_data_test.df, permutations=500) #run a permuted variance partitioning analysis for dealing with highly collinear explanatory variables
varpart.obs <- rdacca.hp(dv=norm_express_data, iv=enviro_data_test.df, var.part = T) #calculate RDA axes and model significance using permutation tests for highly collinear explanatory variables

pdf(file="BTHeatwave_RDA.pdf",height=10,width=14) #create an empty pdf
plot(expression_full.rda, scaling=3) #add rda plot to the pdf
dev.off() #stop printing objects to the pdf
