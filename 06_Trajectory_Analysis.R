library(devtools)
#install_github("mlcollyer/RRPP")
library(RRPP)
library(geomorph)
library(data.table)
library(ggplot2)
library(ggarrow)
library(MuMIn)
library(dplyr)
library(WGCNA)
library(GOfuncR)
library(ggpubr)

salmon_data <- fread("salmon_logcounts_environ.txt", drop=1) #read in dataset with environmental and expression data
enviro_data <- salmon_data[-30,1:16] #create a separate data frame with just environmental data; remove an outlier BP075
express_data <- salmon_data[-30,-(1:16)] #create a separate data frame with just expression data; remove an outlier BP075 
rownames(express_data) <- salmon_data$Sample_ID[-30] #name the expression data rows with sample IDs so future analyses of the expression data can be aligned to environmental data
start_date <- as.Date("2022-07-15") #set the start date
for(i in 1:nrow(enviro_data)){ #create a second Date variable that expresses date as the number of since the start date set above (the day before sampling started)
  enviro_data$Date2[i] <- difftime(as.Date(enviro_data$Date[i]), start_date, units = "days")
}

temp_traj.df <- geomorph.data.frame(expression=as.matrix(stdize(express_data)), Site=as.factor(enviro_data$Site), Date=enviro_data$Date2, Time=as.factor(enviro_data$Sample_Period), Temp=enviro_data$DailyAvgTemp) #create a geomorph data frame with the transcript count data as the morphological data, then with site, sample period, and daily average water temperature as additional variables
temp.fit <- lm.rrpp(expression ~ Temp*Site*Time, data=temp_traj.df, int.first=T, iter=499, seed=8392, print.progress=T) #fit an RRPP model with a three-way interaction between temperature, site, and sample period
anova.lm.rrpp(temp.fit)

null.fit <- lm.rrpp(expression ~ Site, data=temp_traj.df, iter=499, print.progress = T) #fit an RRPP model with just the site effects so we can pull out the consistent effects of the heatwave in the trajectory analysis
temp.traj <- trajectory.analysis(temp.fit, fit.null=null.fit, groups=temp_traj.df$Site, traj.pts = temp_traj.df$Time, print.progress=T) #run a trajectory analysis on the RRPP model with separate trajectories calculated for each sample site and the sample periods as the trajectory points

summary(temp.traj, attribute="MD", show.trajectories=T) #show the results for the magnitude differences test of the trajectory analysis
summary(temp.traj, attribute="TC", angle.type="deg") #show the results for the orientation test of the trajectory analysis
summary(temp.traj, attribute="SD") #show the results for the shape differences test of the trajectory analysis

temp_traj.plot <- plot(temp.traj, pch=as.numeric(enviro_data$Sample_Period) + 10, bg=as.factor(enviro_data$Sample_Period), cex=0.7, col="grey") #create a default plot object from the trajectory analysis to make it easier to extract PC points and trajectory points
add.trajectories(temp_traj.plot) #add the trajectory points to the default plot

trajectory_pc_plot.df <- data.frame(Sample_ID=enviro_data$Sample_ID, Site=enviro_data$Site, Time=as.integer(enviro_data$Sample_Period), PC1=temp_traj.plot$pc.points[,1], PC2=temp_traj.plot$pc.points[,2], PC3=temp_traj.plot$pc.points[,3], PC4=temp_traj.plot$pc.points[,4]) #create a data frame containing the coordinates for the first four principal component axes used in the trajectory analysis
trajectory_traj_plot.df <- data.frame(Site=c(rep("Big Poe", 8), rep("East Branch", 8), rep("Shaver's Creek", 8), rep("Standing Stone", 8)), Sample.Period=rep(c(1:8), 4), Traj1=c(temp.traj$trajectories$obs$`Big Poe`[,1], temp.traj$trajectories$obs$`East Branch`[,1], temp.traj$trajectories$obs$`Shaver's Creek`[,1],  temp.traj$trajectories$obs$`Standing Stone`[,1]), Traj2=c(temp.traj$trajectories$obs$`Big Poe`[,2], temp.traj$trajectories$obs$`East Branch`[,2], temp.traj$trajectories$obs$`Shaver's Creek`[,2],  temp.traj$trajectories$obs$`Standing Stone`[,2]), Traj3=c(temp.traj$trajectories$obs$`Big Poe`[,3], temp.traj$trajectories$obs$`East Branch`[,3], temp.traj$trajectories$obs$`Shaver's Creek`[,3],  temp.traj$trajectories$obs$`Standing Stone`[,3]), Traj4=c(temp.traj$trajectories$obs$`Big Poe`[,4], temp.traj$trajectories$obs$`East Branch`[,4], temp.traj$trajectories$obs$`Shaver's Creek`[,4],  temp.traj$trajectories$obs$`Standing Stone`[,4])) #create a data frame containing the trajectory points on the first four PC axes from the trajectory analysis

full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df, aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.5)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe",], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch",], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek",], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone",], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_point(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period==1,], aes(x=Traj1, y=Traj2), shape=1, size=5, stroke=2, show.legend = F)+
  geom_text(data=trajectory_traj_plot.df, aes(x=Traj1, y=Traj2, label=Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"))+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits 
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 
ggsave("full_trajectory.jpeg", full_trajectory.plot, width=8, height=7, units="in", dpi=300)

traj12.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df[trajectory_pc_plot.df$Time < 3,], aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"), guide="none")+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj23.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df[trajectory_pc_plot.df$Time < 4,], aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period==1:2,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & (trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & (trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & (trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & (trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period==1 | trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3,], aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"), guide="none")+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj34.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df[trajectory_pc_plot.df$Time < 5,], aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & (trajectory_traj_plot.df$Sample.Period==1 | trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & (trajectory_traj_plot.df$Sample.Period==1 | trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & (trajectory_traj_plot.df$Sample.Period==1 | trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & (trajectory_traj_plot.df$Sample.Period==1 | trajectory_traj_plot.df$Sample.Period==2 | trajectory_traj_plot.df$Sample.Period==3),], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period==3:4,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period==3:4,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period==3:4,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period==3:4,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period==1:4,], aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"), guide="none")+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj45.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df[trajectory_pc_plot.df$Time < 6,], aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period < 5,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period < 5,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period < 5,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period < 5,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period >= 4 & trajectory_traj_plot.df$Sample.Period <= 5,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period >= 4 & trajectory_traj_plot.df$Sample.Period <= 5,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period >= 4 & trajectory_traj_plot.df$Sample.Period <= 5,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period >= 4 & trajectory_traj_plot.df$Sample.Period <= 5,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period < 6,], aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"), guide="none")+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj56.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df[trajectory_pc_plot.df$Time < 7,], aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period < 6,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period < 6,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period < 6,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period < 6,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period >= 5 & trajectory_traj_plot.df$Sample.Period <= 6,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period >= 5 & trajectory_traj_plot.df$Sample.Period <= 6,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period >= 5 & trajectory_traj_plot.df$Sample.Period <= 6,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period >= 5 & trajectory_traj_plot.df$Sample.Period <= 6,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period < 7,], aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"), guide="none")+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj67.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df[trajectory_pc_plot.df$Time < 8,], aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period < 7,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period < 7,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period < 7,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period < 7,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period >= 6 & trajectory_traj_plot.df$Sample.Period <= 7,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period >= 6 & trajectory_traj_plot.df$Sample.Period <= 7,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period >= 6 & trajectory_traj_plot.df$Sample.Period <= 7,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period >= 6 & trajectory_traj_plot.df$Sample.Period <= 7,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Sample.Period < 8,], aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"), guide="none")+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj78.plot <- full_trajectory.plot <- ggplot(data=trajectory_pc_plot.df, aes(x=PC1, y=PC2, color=Site))+ #create a ggplot with the PC points for axis 1 and 2 of the trajectory analysis, color the points by site
  geom_point(alpha=0.2)+ #add the pc points to the plot, make them slightly transparent
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period < 8,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period < 8,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period < 8,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period < 8,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7, alpha=0.2, linetype = 2)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Big Poe" & trajectory_traj_plot.df$Sample.Period >= 7 & trajectory_traj_plot.df$Sample.Period <= 8,], aes(x=Traj1, y=Traj2, color="Big Poe"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Big Poe site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="East Branch" & trajectory_traj_plot.df$Sample.Period >= 7 & trajectory_traj_plot.df$Sample.Period <= 8,], aes(x=Traj1, y=Traj2, color="East Branch"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the East Branch site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Shaver's Creek" & trajectory_traj_plot.df$Sample.Period >= 7 & trajectory_traj_plot.df$Sample.Period <= 8,], aes(x=Traj1, y=Traj2, color="Shaver's Creek"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Shaver's Creek site
  geom_arrow_chain(data=trajectory_traj_plot.df[trajectory_traj_plot.df$Site=="Standing Stone" & trajectory_traj_plot.df$Sample.Period >= 7 & trajectory_traj_plot.df$Sample.Period <= 8,], aes(x=Traj1, y=Traj2, color="Standing Stone"), linewidth=0.7)+ #add arrows showing the direction and order of the trajectories for the Standing Stone site
  geom_text(data=trajectory_traj_plot.df, aes(x=Traj1, y=Traj2, label = Sample.Period, color=Site), size=5)+ #add the trajectory points to the plot, make them slightly larger than the pc points and solid black, vary the point shape by sample site
  scale_color_manual(values=c("Big Poe"= "#F8766D", "East Branch"="#7CAE00", "Shaver's Creek"="#00BFC4", "Standing Stone"="#C77CFF"))+ #color the arrows by site and add them to the legend with the correct labels
  theme_classic()+ #change the plot theme to make it look neater
  lims(x=c(-110, 160), y=c(-170, 170))+ #set axis limits
  labs(x="PC1 (14.85%)\n <- Ribosomal Components | No Enriched GO Terms ->", y="PC2 (10.86%)\n <- Oxygen Carrier Activity | Immune Response ->") #add axis labels 

traj.panels <- ggarrange(traj12.plot, traj23.plot, traj34.plot, traj45.plot, traj56.plot, traj67.plot, traj78.plot, ncol=7 , nrow=1, widths = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1.1))
ggsave("traj_panels.jpeg", traj.panels, width=24, height=4, units="in", dpi=300) #save a jpeg of Figure 4

daily_temps.df <- read.delim("BrookTrout_Heatwaves_DailyTemps.txt") #read in daily temperature data frame
daily_temps.df$Date <- as.Date.character(daily_temps.df$Date) #convert Date column in daily temps data frame to R Date format
trajectory_traj_plot.df$Date <- c("2022-07-20", "2022-07-23", "2022-07-27", "2022-08-01", "2022-08-04", "2022-08-08", "2022-08-11", "2022-08-18", "2022-07-19", "2022-07-22", "2022-07-26", "2022-08-02", "2022-08-05", "2022-08-09", "2022-08-12", "2022-08-19",  "2022-07-19", "2022-07-22", "2022-07-26", "2022-08-02", "2022-08-05", "2022-08-09", "2022-08-12", "2022-08-19", "2022-07-20", "2022-07-23", "2022-07-27", "2022-08-01", "2022-08-04", "2022-08-08", "2022-08-11", "2022-08-18") #add a sampling dates column to the trajectory point data frame
trajectory_traj_plot.df$Date <- as.Date.character(trajectory_traj_plot.df$Date) #convert the dates column in the trajectory point data frame to R date format

Traj1_temp.plot <- ggplot()+ #create new ggplot object
  geom_path(data=daily_temps.df, aes(x=Date, y=Temp.avg, color=Site), linetype=2)+ #plot path of daily average temperatures
  geom_point(data=trajectory_traj_plot.df, aes(x=Date, y=(Traj1/12) + 15), size=2)+ #plot points showing trajectory point value on PC1 during each sampling date (rescaled to fit on same plot as temperature)
  geom_arrow_chain(data=trajectory_traj_plot.df, aes(x=Date, y=(Traj1/12) + 15, color=Site))+ #plot arrows linking trajectory points
  theme_classic()+ #change the plot theme
  scale_y_continuous(sec.axis=sec_axis(transform=~(.-15) * 12, name="Trajectory PC1"))+ #add a secondary y-axis to display trajectory PC1 values
  xlim(c(as.Date.character("2022-07-15"), as.Date.character("2022-08-20")))+ #standardize the x axis to focus on the sampling period
  ylab("Daily Average Temperature (ºC)")+ #add a y-axis label
  facet_wrap(~Site, nrow=1) #create a row of panels showing each stream individually

Traj2_temp.plot <- ggplot()+ #create new ggplot object
  geom_path(data=daily_temps.df, aes(x=Date, y=Temp.avg, color=Site), linetype=2)+ #plot path of daily average temperatures
  geom_point(data=trajectory_traj_plot.df, aes(x=Date, y=(Traj2/20) + 15), size=2)+ #plot points showing trajectory point value on PC2 during each sampling date (rescaled to fit on same plot as temperature)
  geom_arrow_chain(data=trajectory_traj_plot.df, aes(x=Date, y=(Traj2/20) + 15, color=Site))+ #plot arrows linking trajectory points
  theme_classic()+ #change the plot theme
  scale_y_continuous(sec.axis=sec_axis(transform=~(.-15) * 20, name="Trajectory PC2"))+ #add a secondary y-axis to display trajectory PC2 values
  xlim(c(as.Date.character("2022-07-15"), as.Date.character("2022-08-20")))+ #standardize the x axis to focus on the sampling period
  ylab("Daily Average Temperature (ºC)")+ #add a y-axis label
  facet_wrap(~Site, nrow=1) #create a row of panels showing each stream individually

traj_temp.plots <- ggarrange(Traj1_temp.plot, Traj2_temp.plot, nrow=2)
ggsave("traj_temp.jpeg", traj_temp.plots, width=16, height=5, unit="in", dpi=300)

express_site.day <- salmon_data %>% group_by(Site, as.factor(Sample_Period)) %>% summarize( across(starts_with("QSF"), mean)) #get average PC scores (per sampling event) for each transcript
temp_traj.loadings <- cor(temp.traj$LS.means$obs, express_site.day[,-c(1:2)]) #calculate the loadings of transcripts on the trajectory PC axes
temp_traj.loadings_pval <- as.data.frame(corPvalueStudent(as.matrix(t(temp_traj.loadings)), nSamples = nrow(express_data))) #calculate p-values for significance of loading values

traj1_siggenes <- rownames(temp_traj.loadings_pval[temp_traj.loadings_pval$PC1 < 0.05/32670,]) #create a list of genes that are significantly related to PC1 (Bonferroni corrected)
traj2_siggenes <- rownames(temp_traj.loadings_pval[temp_traj.loadings_pval$PC2 < 0.05/32670,]) #create a list of genes that are significantly related to PC2 (Bonferroni corrected)
traj1_siggenes_up <- colnames(temp_traj.loadings[,temp_traj.loadings[1,] > 0])[colnames(temp_traj.loadings[,temp_traj.loadings[1,] > 0]) %in% traj1_siggenes] #split the list of genes related to PC1 into positively correlated
traj1_siggenes_down <- colnames(temp_traj.loadings[,temp_traj.loadings[1,] < 0])[colnames(temp_traj.loadings[,temp_traj.loadings[1,] < 0]) %in% traj1_siggenes] #and negatively correlated
traj2_siggenes_up <- colnames(temp_traj.loadings[,temp_traj.loadings[1,] > 0])[colnames(temp_traj.loadings[,temp_traj.loadings[1,] > 0]) %in% traj2_siggenes] #split the list of genes related to PC2 into positively correlated
traj2_siggenes_down <- colnames(temp_traj.loadings[,temp_traj.loadings[1,] < 0])[colnames(temp_traj.loadings[,temp_traj.loadings[1,] < 0]) %in% traj2_siggenes] #and negatively correlated

Phylofish_GO <- fread("Phylofish_GOData.txt") #read in Phylofish gene ontology data frame
Phylofish_GO_salmon <- Phylofish_GO[Phylofish_GO$Name %in% colnames(express_data),] #filter the GO data frame to only include transcripts used in the analysis
Phylofish_GO_salmon <- Phylofish_GO_salmon[Phylofish_GO_salmon$`Go name` != "",] #remove records with no associated GO terms
Phylofish_GO_salmon_input <- data.frame(gene=Phylofish_GO_salmon$Name, go_ID=Phylofish_GO_salmon$Code) #create a simplified data frame of the GO terms to input into the go_enrich function

traj1_up_go_input <- data.frame(gene=rownames(temp_traj.loadings_pval), candidate=as.numeric(rownames(temp_traj.loadings_pval) %in% traj1_siggenes_up)) #create input data frame for enrichment analysis where all transcripts are given a value of 1 if they are significantly positively associated with PC1 or 0 if not
traj1_up.GO <- go_enrich(traj1_up_go_input, annotations = Phylofish_GO_salmon_input, test="hyper", silent = T, n_randsets = 1000) #run the GO analysis, using the list of genes significantly positively related to PC1 as the candidate gene list in a hypergeometric test against the full gene list
traj1_up.GO.ref <- refine(traj1_up.GO) #refine the results of the GO analysis (Note: will throw up error b/c PC1 had no GO terms significantly enriched in the positive direction)

traj1_down_go_input <- data.frame(gene=rownames(temp_traj.loadings_pval), candidate=as.numeric(rownames(temp_traj.loadings_pval) %in% traj1_siggenes_down)) #create input data frame for enrichment analysis where all transcripts are given a value of 1 if they are significantly negatively associated with PC1 or 0 if not
traj1_down.GO <- go_enrich(traj1_down_go_input, annotations = Phylofish_GO_salmon_input, test="hyper", silent = T, n_randsets = 1000) #run the GO analysis, using the list of genes significantly negatively related to PC1 as the candidate gene list in a hypergeometric test against the full gene list
traj1_down.GO.ref <- refine(traj1_down.GO, annotations = Phylofish_GO_salmon_input) #refine the results of the GO analysis

traj2_up_go_input <- data.frame(gene=rownames(temp_traj.loadings_pval), candidate=as.numeric(rownames(temp_traj.loadings_pval) %in% traj2_siggenes_up)) #create input data frame for enrichment analysis where all transcripts are given a value of 1 if they are significantly positively associated with PC2 or 0 if not
traj2_up.GO <- go_enrich(traj2_up_go_input, annotations = Phylofish_GO_salmon_input, test="hyper", silent = T, n_randsets = 1000) #run the GO analysis, using the list of genes significantly positively related to PC2 as the candidate gene list in a hypergeometric test against the full gene list
traj2_up.GO.ref <- refine(traj2_up.GO, annotations = Phylofish_GO_salmon_input) #refine the results of the GO analysis

traj2_down_go_input <- data.frame(gene=rownames(temp_traj.loadings_pval), candidate=as.numeric(rownames(temp_traj.loadings_pval) %in% traj2_siggenes_down)) #create input data frame for enrichment analysis where all transcripts are given a value of 1 if they are significantly negatively associated with PC2 or 0 if not
traj2_down.GO <- go_enrich(traj2_down_go_input, annotations = Phylofish_GO_salmon_input, test="hyper", silent = T, n_randsets = 1000) #run the GO analysis, using the list of genes significantly negatively related to PC2 as the candidate gene list in a hypergeometric test against the full gene list
traj2_down.GO.ref <- refine(traj2_down.GO, annotations = Phylofish_GO_salmon_input) #refine the results of the GO analysis

