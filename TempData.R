#load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)

temp_data <- read_excel("../2022heatwave temps.xlsx", sheet="All temp data") #read in excel sheet with all data from all temperature loggers in all four streams
colnames(temp_data) <- c("Site", "Logger", "Date", "Temp") #rename the columns so they are easier to work with
temp_data_water.df <- temp_data[temp_data$Logger!="SCA_Air" & temp_data$Logger!="EB_Air" & temp_data$Logger!="BP_Air" & temp_data$Logger!="SS_Air",] #create a data frame with only the water temperature 
temp_data_water_avg.df <- temp_data_water.df %>% group_by(Date, Site) %>% summarise(Temp = mean(Temp)) #average temperatures collected at the same time point across loggers from the same stream 
temp_data_water_avg.df <- temp_data_water_avg.df[temp_data_water_avg.df$Date < as.POSIXct("2022-09-01"),] #only keep temperature records from before September 1

ggplot(data=temp_data_water_avg.df, aes(x=Date, y=Temp, color=Site))+ #create ggplot with date on the x-axis, temperature on the y-axis, and separate sites using color
  geom_line()+ #add a line object connecting temperature record points through time
  theme_classic() #change the plot theme to make it look neater

temp_data_water_dailyavg.df <- temp_data_water_avg.df %>% group_by(Date=as.Date(Date), Site) %>% summarise(Temp=mean(Temp)) #calculate the average daily temperatures for each stream (using time point averages across loggers)
temp_data_water_dailymax.df <- temp_data_water_avg.df %>% group_by(Date=as.Date(Date), Site) %>% summarise(Temp=max(Temp, na.rm = T)) #calculate the maximum daily temperatures for each stream (using time point averages across loggers)
temp_data_water_dailyrange.df <- temp_data_water_avg.df %>% group_by(Date=as.Date(Date), Site) %>% summarise(Temp=max(Temp)-min(Temp)) #calculate the minimum daily temperatures for each stream (using time point averages across loggers)
temp_data_water_dailysd.df <- temp_data_water_avg.df %>% group_by(Date=as.Date(Date), Site) %>% summarise(Temp=sd(Temp)) #calculate the standard deviation in daily temperatures for each stream (using time point averages across loggers)

ggplot(data=temp_data_water_dailyavg.df, aes(x=Date, y=Temp, color=Site))+ #create ggplot with date on the x-axis, temperature on the y-axis, and separate sites using color
  geom_line()+ #add a line object connecting daily average temperatures through time
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Avg Daily Temp (ºC)") #edit the y-axis label

ggplot(data=temp_data_water_dailymax.df, aes(x=Date, y=Temp, color=Site))+ #create ggplot with date on the x-axis, temperature on the y-axis, and separate sites using color
  geom_line()+ #add a line object connecting daily maximum temperatures through time
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Maximum Daily Temp (ºC)") #edit the y-axis label

ggplot(data=temp_data_water_dailyrange.df, aes(x=Date, y=Temp, color=Site))+ #create ggplot with date on the x-axis, temperature on the y-axis, and separate sites using color
  geom_line()+ #add a line object connecting daily temperature range through time
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Daily Temp Range (ºC)") #edit the y-axis label

ggplot(data=temp_data_water_dailysd.df, aes(x=Date, y=Temp, color=Site))+ #create ggplot with date on the x-axis, temperature on the y-axis, and separate sites using color
  geom_line()+ #add a line object connecting daily temperature standard deviation through time
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Daily Temp Standard Dev (ºC)") #edit the y-axis label

daily_temps.df <- left_join(temp_data_water_dailyavg.df, temp_data_water_dailymax.df, by=c("Date", "Site"), suffix=c(".avg", ".max")) #create a data frame that combines the average and maximum temperatures matched to each time point recorded from all four streams
write.table(daily_temps.df, file="BrookTrout_Heatwaves_DailyTemps.txt", quote=F, row.names = F, sep="\t") #save the data frame as a tab-delimited text file

trout_field_data.df <- read_excel("../2022 trout sampling data.xlsx", sheet="HW trout field data") #read in data sheet for individual fish with information on where and when they were sampled
trout_RNA_QC.df <- read_excel("../2022 trout sampling data.xlsx", sheet="RNA QC") #read in data sheet for individual fish with RNA quality control data and sexing assay results
trout_sample_data.df <- left_join(trout_RNA_QC.df, trout_field_data.df, by="Sample ID") #merge the field data information with the RNA QC data for all fish that have RNA extracted
trout_sample_data.df <- trout_sample_data.df[trout_sample_data.df$Species=="Brook",] #only keep brook trout samples
date(trout_sample_data.df$`Bucket Begin`) <- trout_sample_data.df$Date #replace the default date (Jan 1, 1899) that R adds to times without dates in the "Bucket Begin" column with the sampling date from the "Date" column
date(trout_sample_data.df$`Bucket End`) <- trout_sample_data.df$Date #replace the default date (Jan 1, 1899) that R adds to times without dates in the "Bucket End" column with the sampling date from the "Date" column

BP_sampletemp.plot <- ggplot(data=temp_data_water_dailyavg.df[temp_data_water_dailyavg.df$Site=="Big Poe",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis, and set color to use for the line
  geom_line(color="#F8766D")+ #add a line object connecting daily average temperatures through time
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="Big Poe",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Avg Daily Temp (ºC)", title="Big Poe")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

EB_sampletemp.plot <- ggplot(data=temp_data_water_dailyavg.df[temp_data_water_dailyavg.df$Site=="East Branch",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_line(color="#7CAE00")+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="East Branch",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Avg Daily Temp (ºC)", title = "East Branch")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

SCA_sampletemp.plot <- ggplot(data=temp_data_water_dailyavg.df[temp_data_water_dailyavg.df$Site=="Shaver's Creek",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis, and set color to use for the line
  geom_line(color="#00BFC4")+ #add a line object connecting daily average temperatures through time
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="Shaver's Creek Above",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Avg Daily Temp (ºC)", title="Shaver's Creek")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

SS_sampletemp.plot <- ggplot(data=temp_data_water_dailyavg.df[temp_data_water_dailyavg.df$Site=="Standing Stone",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_line(color="#C77CFF")+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="Standing Stone",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Avg Daily Temp (ºC)", title = "Standing Stone")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

stacked_meantemp.plot <- grid.arrange(BP_sampletemp.plot, EB_sampletemp.plot, SCA_sampletemp.plot, SS_sampletemp.plot, nrow=4) #stack the temperature plots of the four rivers with sampling events indicated
ggsave("AvgTemp_SampleDates.jpg", stacked_meantemp.plot, width=6, height=14, units="in")

BP_hightemp.plot <- ggplot(data=temp_data_water_dailymax.df[temp_data_water_dailymax.df$Site=="Big Poe",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis, and set color to use for the line
  geom_line(color="#F8766D")+ #add a line object connecting daily average temperatures through time
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="Big Poe",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Daily High Temp (ºC)", title="Big Poe")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(15, 26) #add y-axis limits so they will be consistent across different plots

EB_hightemp.plot <- ggplot(data=temp_data_water_dailymax.df[temp_data_water_dailymax.df$Site=="East Branch",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_line(color="#7CAE00")+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="East Branch",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Daily High Temp (ºC)", title = "East Branch")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(15, 26) #add y-axis limits so they will be consistent across different plots

SCA_hightemp.plot <- ggplot(data=temp_data_water_dailymax.df[temp_data_water_dailymax.df$Site=="Shaver's Creek",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis, and set color to use for the line
  geom_line(color="#00BFC4")+ #add a line object connecting daily average temperatures through time
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="Shaver's Creek Above",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Daily High Temp (ºC)", title="Shaver's Creek")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(15, 26) #add y-axis limits so they will be consistent across different plots

SS_hightemp.plot <- ggplot(data=temp_data_water_dailymax.df[temp_data_water_dailymax.df$Site=="Standing Stone",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_line(color="#C77CFF")+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.Date(trout_sample_data.df[trout_sample_data.df$Site=="Standing Stone",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Daily High Temp (ºC)", title = "Standing Stone")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(15, 26) #add y-axis limits so they will be consistent across different plots

stacked_hightemp.plot <- grid.arrange(BP_hightemp.plot, EB_hightemp.plot, SCA_hightemp.plot, SS_hightemp.plot, nrow=4) #stack the temperature plots of the four rivers with sampling events indicated
ggsave("HighTemp_SampleDates.jpg", stacked_hightemp.plot, width=6, height=14, units="in")

sample_event_tempvars.df <- data.frame(Sample_Time=unique(trout_sample_data.df$`Bucket Begin`)) #create a data frame with one entry for each unique sampling event
datesite.df <- data.frame(Sample_Time=trout_sample_data.df$`Bucket Begin`, Site=trout_sample_data.df$Site, Sample_Temp=trout_sample_data.df$`Temperature in Stream (C)`) #create a data frame with the sampling date, site, and water temperature measured at the time of sampling for each fish
sample_event_tempvars.df <- left_join(sample_event_tempvars.df, datesite.df, by="Sample_Time", multiple="any") #match each unique sampling event to the correct site & temperature
sample_event_tempvars.df$Site <- sub("Shaver's Creek Above", "Shaver's Creek", sample_event_tempvars.df$Site) #edit the site name for Shaver's Creek so it matches the site name in the data frame with temperature logger data

#calculate the average water temperature in the 24h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$OneDayMean[i] <- mean(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-1440), ]$Temp, na.rm = T)
}

#calculate the average water temperature in the 72h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$ThreeDayMean[i] <- mean(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440*3)), ]$Temp, na.rm = T)
}

#calculate the highest recorded water temperature in the 24h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$OneDayMax[i] <- max(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440)), ]$Temp, na.rm = T)
}

#calculate the highest recorded water temperature in the 72h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$ThreeDayMax[i] <- max(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440*3)), ]$Temp, na.rm = T)
}

#calculate the change in water temperature from the same (or closest to the same) time of day 24h prior to the water temperature during each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$OneDayChange[i] <- mean(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440)), ]$Temp, na.rm = T) - mean(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= minutes(-1440) & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-2880), ]$Temp, na.rm = T)
}

#calculate the change in water temperature from the same (or closest to the same) time of day 72h prior to the water temperature during each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$ThreeDayChange[i] <- sample_event_tempvars.df$Sample_Temp[i] - temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date==min(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440*3)), ]$Date),]$Temp
}

#calculate the standard deviation in water temperature from the 24h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$OneDayVar[i] <- sd(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440)), ]$Temp, na.rm = T)
}

#calculate the standard deviation in water temperature from the 72h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$ThreeDayVar[i] <- sd(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440*3)), ]$Temp, na.rm = T)
}

#calculate the number of temperature logger readings above 21ºC (brook trout thermal limit) from the 24h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$OneDayStress[i] <- sum(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440)), ]$Temp > 21, na.rm = T)
}

#calculate the number of temperature logger readings above 21ºC (brook trout thermal limit) from the 72h prior to each sampling event
for(i in 1:nrow(sample_event_tempvars.df)){
  sample_event_tempvars.df$ThreeDayStress[i] <- sum(temp_data_water_avg.df[temp_data_water_avg.df$Site==sample_event_tempvars.df$Site[i] & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] <= 0 & temp_data_water_avg.df$Date - sample_event_tempvars.df$Sample_Time[i] >= minutes(-(1440*3)), ]$Temp > 21, na.rm = T)
}

sample_event_tempvars.df <- sample_event_tempvars.df[-c(3,4,7,8,15,17,19,21,23,25,26,30,31,33,35,36,38,42,48,50,52,53,55,56,58,59,62,63),] #remove records of sampled brook trout that were not selected for sequencing

#Under/over temperature graphs centered on 20º C, roughly the temp where brook trout begin to be thermally stressed
watertemp_underover_stress_SS.plot <- ggplot(data=temp_data_water_avg.df[temp_data_water_avg.df$Site=="Standing Stone",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_ribbon(aes(ymin=20, ymax=pmax(Temp, 20), x=Date), fill="red", alpha=0.5)+
  geom_ribbon(aes(ymin=pmin(20, Temp), ymax=20, x=Date), fill="blue", alpha=0.5)+
  #geom_line(color="#C77CFF", linewidth=1)+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.POSIXct(trout_sample_data.df[trout_sample_data.df$Site=="Standing Stone",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Temp (ºC)", title = "Standing Stone")+ #edit the y-axis label
  xlim(as.POSIXct("2022-07-15"), as.POSIXct("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots
watertemp_underover_stress_SS.plot

watertemp_underover_stress_SCA.plot <- ggplot(data=temp_data_water_avg.df[temp_data_water_avg.df$Site=="Shaver's Creek",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_ribbon(aes(ymin=20, ymax=pmax(Temp, 20), x=Date), fill="red", alpha=0.5)+
  geom_ribbon(aes(ymin=pmin(20, Temp), ymax=20, x=Date), fill="blue", alpha=0.5)+
  #geom_line(color="#C77CFF", linewidth=1)+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.POSIXct(trout_sample_data.df[trout_sample_data.df$Site=="Shaver's Creek Above",]$`Bucket Begin`), linetype=2, alpha=0.6)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Temp (ºC)", title = "Shaver's Creek")+ #edit the y-axis label
  xlim(as.POSIXct("2022-07-15"), as.POSIXct("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots
watertemp_underover_stress_SCA.plot

watertemp_underover_stress_BP.plot <- ggplot(data=temp_data_water_avg.df[temp_data_water_avg.df$Site=="Big Poe",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_ribbon(aes(ymin=20, ymax=pmax(Temp, 20), x=Date), fill="red", alpha=0.5)+
  geom_ribbon(aes(ymin=pmin(20, Temp), ymax=20, x=Date), fill="blue", alpha=0.5)+
  #geom_line(color="#C77CFF", linewidth=1)+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.POSIXct(trout_sample_data.df[trout_sample_data.df$Site=="Big Poe",]$`Bucket Begin`), linetype=2, alpha=0.25)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Temp (ºC)", title = "Big Poe")+ #edit the y-axis label
  xlim(as.POSIXct("2022-07-15"), as.POSIXct("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots
watertemp_underover_stress_BP.plot

watertemp_underover_stress_EB.plot <- ggplot(data=temp_data_water_avg.df[temp_data_water_avg.df$Site=="East Branch",], aes(x=Date, y=Temp))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_ribbon(aes(ymin=20, ymax=pmax(Temp, 20), x=Date), fill="red", alpha=0.5)+
  geom_ribbon(aes(ymin=pmin(20, Temp), ymax=20, x=Date), fill="blue", alpha=0.5)+
  #geom_line(color="#C77CFF", linewidth=1)+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_vline(xintercept=as.POSIXct(trout_sample_data.df[trout_sample_data.df$Site=="East Branch",]$`Bucket Begin`), linetype=2, alpha=0.25)+
  theme_classic()+ #change the plot theme to make it look neater
  labs(y="Temp (ºC)", title = "East Branch")+ #edit the y-axis label
  xlim(as.POSIXct("2022-07-15"), as.POSIXct("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots
watertemp_underover_stress_EB.plot

watertemp_underover_stress_plots.plot <- grid.arrange(watertemp_underover_stress_BP.plot, watertemp_underover_stress_EB.plot, watertemp_underover_stress_SCA.plot, watertemp_underover_stress_SS.plot, nrow=4)
ggsave("WaterTemp_UnderOver_Stress_Plots.jpg", watertemp_underover_stress_plots.plot, width = 6, height = 12, units = "in")

