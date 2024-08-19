library(ggmap)
library(ggplot2)
library(cowplot)
library(sp)
library(sf)
library(raster)
library(grid)
library(ggspatial)
library(readxl)
library(tidyverse)
library(gridExtra)

Big_Poe <- c("Big Poe", 40.82, -77.48) #set lat/long for Big Poe sample site
East_Branch <- c("East Branch", 40.662, -77.726) #set lat/long for East Branch sample site
Shavers_Creek <- c("Shaver's Creek", 40.685, -77.90) #set lat/long for Shaver's Creek sample site
Standing_Stone <- c("Standing Stone", 40.712, -77.70) #set lat/long for Standing Stone sample site
BKT_sites <- rbind(Big_Poe, East_Branch, Shavers_Creek, Standing_Stone) #combine lat/long for four sites into single matrix
BKT_sites <- as.data.frame(BKT_sites) #change sample site lat/long matrix to a data frame
colnames(BKT_sites) <- c("Site", "lat", "lon") #set column names for sample site lat/long data frame
BKT_sites.spdf <- SpatialPointsDataFrame(cbind(as.numeric(BKT_sites$lon), as.numeric(BKT_sites$lat)), data=BKT_sites) #create a spatial points data frame from sample site lat/long data frame
proj4string(BKT_sites.spdf) <- CRS("+init=epsg:4326") #give the spatial data frame the standard WGS 1984 crs
BKT_sites.sf <- st_as_sf(BKT_sites.spdf, coords=1:2) #convert the spatial data frame to a sf object for plotting in ggplot
BKT_sites.sf$lon <- as.numeric(BKT_sites.sf$lon) #make sure longitude values are in numeric format 
BKT_sites.sf$lat <- as.numeric(BKT_sites.sf$lat) #make sure latitude values are in numeric format 

register_stadiamaps("f179183b-b95a-4b00-b32a-7f13310e90e8") #add Stadia Maps api key
PA_basemap <- get_stadiamap(bbox=c(left=-78.0, bottom=40.54, right=-77.4, top=40.87), maptype = "stamen_terrain_background", zoom=12, source="stadia") #download terrain map of central pennsylvania & surrounding area from stamen

PA_basemap_lighter <- matrix(adjustcolor(PA_basemap, alpha.f=0.5), nrow=nrow(PA_basemap)) #create a new matrix with adjused colors of the Pennsylvania terrain map so they appear lighter and elements to be added later will stand out more
attributes(PA_basemap_lighter) <- attributes(PA_basemap) #apply all the attributes from the original Pennsylvania terrain map to the new map with lighter colors

PA_basemap.ggplot <- ggmap(PA_basemap_lighter)+ #save a plot of the light colored Pennsylvania terrain map that can be used as a background for later maps
  theme_half_open(font_size=10)+ #change the map theme to get rid of the background and make the axes look prettier
  scale_y_continuous(expand=expansion(mult=c(0,0)))+ #get rid of extra white space between basemap and y axis
  scale_x_continuous(expand=expansion(mult=c(0,0)))+ #get rid of extra white space between basemap and x axis
  labs(x="Longitude", y="Latitude")+ #add axis labels to the plot
  coord_sf(crs=4326)
  
#Add sample site locations to the study area basemap
PA_basemap.ggplot +
  geom_sf(data=BKT_sites.sf, aes(color=Site))+ #add sample point coordinates to the map and color them by original Region designations from the data sheets
  coord_sf(xlim=c(-78, -77.4), ylim=c(40.54, 40.87)) #add a coordinate system to the plot, assume the input crs is 4326 for any elements without an available crs

#EXTRACT HUC8 LEVEL SUBBASINS FOR EACH SAMPLE POINT
NHD_flowline.shp <- st_read(dsn="~/Desktop/FHC_raw_microsat_data/FHC_microsat_analysis/FHC_microsat_popgen/NHDPlusMA/NHDPlus02/NHDSnapshot/Hydrography/", layer="NHDFlowline") #read in shapefile with National Hydrology Database flowlines for the mid-Atlantic region
four_rivers.shp <- NHD_flowline.shp[NHD_flowline.shp$GNIS_NAME=="Big Poe Creek" | NHD_flowline.shp$GNIS_NAME=="East Branch Standing Stone Creek" | NHD_flowline.shp$GNIS_NAME=="Shaver Creek" | NHD_flowline.shp$GNIS_NAME=="Standing Stone Creek",] #pull out the flowlines of the study streams from the mid-Atlantic NHD flowlines shapefile
four_rivers.shp <- st_zm(four_rivers.shp) #format the study streams shapefile so it can be plotted with ggplot
secondary_rivers.shp <- NHD_flowline.shp[NHD_flowline.shp$GNIS_NAME=="Penns Creek" | NHD_flowline.shp$GNIS_NAME=="Little Poe Creek" | NHD_flowline.shp$GNIS_NAME=="Laurel Run" | NHD_flowline.shp$GNIS_NAME=="Greenlee Run" | NHD_flowline.shp$GNIS_NAME=="Croyle Run" | NHD_flowline.shp$GNIS_NAME=="Black Lick Run" | NHD_flowline.shp$GNIS_NAME=="Ross Run" | NHD_flowline.shp$GNIS_NAME=="Detweiler Run" | NHD_flowline.shp$GNIS_NAME=="Boal Gap Run" | NHD_flowline.shp$GNIS_NAME=="Sinking Creek" | NHD_flowline.shp$GNIS_NAME=="Galbraith Gap Run" | NHD_flowline.shp$GNIS_NAME=="Spring Creek" | NHD_flowline.shp$GNIS_NAME=="Cedar Run" | NHD_flowline.shp$GNIS_NAME=="Potter Run" | NHD_flowline.shp$GNIS_NAME=="Britton Run" | NHD_flowline.shp$GNIS_NAME=="Laurel Creek" | NHD_flowline.shp$GNIS_NAME=="Lingle Creek" | NHD_flowline.shp$GNIS_NAME=="Tea Creek" | NHD_flowline.shp$GNIS_NAME=="Honey Creek" | NHD_flowline.shp$GNIS_NAME=="Havice Creek" | NHD_flowline.shp$GNIS_NAME=="Treaster Run" | NHD_flowline.shp$GNIS_NAME=="Strodes Run" | NHD_flowline.shp$GNIS_NAME=="Buck Run" | NHD_flowline.shp$GNIS_NAME=="Kishacoquillas Creek" | NHD_flowline.shp$GNIS_NAME=="Meadow Creek" | NHD_flowline.shp$GNIS_NAME=="Panther Run" | NHD_flowline.shp$GNIS_NAME=="Pine Swamp Run" | NHD_flowline.shp$GNIS_NAME=="Delaware Creek" | NHD_flowline.shp$GNIS_NAME=="Weikert Run" | NHD_flowline.shp$GNIS_NAME=="Croyle Run" | NHD_flowline.shp$GNIS_NAME=="Swift Run" | NHD_flowline.shp$GNIS_NAME=="Rock Run" | NHD_flowline.shp$GNIS_NAME=="Muddy Creek" | NHD_flowline.shp$GNIS_NAME=="Kettle Run" | NHD_flowline.shp$GNIS_NAME=="Lost Creek" | NHD_flowline.shp$GNIS_NAME=="Cocolamus Creek" | NHD_flowline.shp$GNIS_NAME=="Stony Run" | NHD_flowline.shp$GNIS_NAME=="South Branch Middle Creek" | NHD_flowline.shp$GNIS_NAME=="Ulsh Gap Run" | NHD_flowline.shp$GNIS_NAME=="Jacks Creek" | NHD_flowline.shp$GNIS_NAME=="Mowry Run" | NHD_flowline.shp$GNIS_NAME=="Shindle Run" | NHD_flowline.shp$GNIS_NAME=="Gregory Run" | NHD_flowline.shp$GNIS_NAME=="Henrys Run" | NHD_flowline.shp$GNIS_NAME=="Armond Run" | NHD_flowline.shp$GNIS_NAME=="Herod Run" | NHD_flowline.shp$GNIS_NAME=="Soft Run" | NHD_flowline.shp$GNIS_NAME=="Little Kishacoquillas Creek" | NHD_flowline.shp$GNIS_NAME=="Spruce Creek" | NHD_flowline.shp$GNIS_NAME=="Beaver Branch" | NHD_flowline.shp$GNIS_NAME=="Halfmoon Creek" | NHD_flowline.shp$GNIS_NAME=="Warriors Mark Run" | NHD_flowline.shp$GNIS_NAME=="Slab Cabin Run" | NHD_flowline.shp$GNIS_NAME=="Roaring Run" | NHD_flowline.shp$GNIS_NAME=="Shingletown Branch" | NHD_flowline.shp$GNIS_NAME=="Buffalo Run" | NHD_flowline.shp$GNIS_NAME=="Logan Branch",] #pull out the flowline shapefiles of the other streams in the watershed to plot
secondary_rivers.shp <- st_zm(secondary_rivers.shp) #format the shapefile of the other streams in the watershed so they can be plotted with ggplot

#READ IN TEMPERATURE AND SAMPLE EVENT DATE & TIME DATA
daily_temps.df <- read.delim("~/Desktop/Brook_Trout_Heatwave/Brook_Trout_Heatwave_SampleSelect/BrookTrout_Heatwaves_DailyTemps.txt") #read in the temperature logger data
daily_temps.df$Date <- as.Date(daily_temps.df$Date) #convert the date column in the temperature data frame to Date format
trout_sample_data.df <- read_excel("~/Desktop/Brook_Trout_Heatwave/2022 trout sampling data.xlsx", sheet="HW trout field data") #read in the trout sampling data to get sampling dates
trout_sample_data.df$Date <- as.Date(trout_sample_data.df$Date) #convert the date column in the trout data to Date format
trout_sample_data.df <- trout_sample_data.df[trout_sample_data.df$Date != as.Date("2022-08-16"), ] #remove the August 16 sampling events since those samples were not sequenced
trout_sample_data.df <- trout_sample_data.df[trout_sample_data.df$Date != as.Date("2022-08-15"), ] #remove the August 15 sampling events since those samples were not sequenced
trout_sample_data.df$Site[trout_sample_data.df$Site=="Shaver's Creek Above"] <- "Shaver's Creek" #edit the name for Shaver's Creek in the trout data frame so it matches the other data frames
trout_sample_data.df <- left_join(trout_sample_data.df, daily_temps.df, by=c("Date", "Site")) #combine the trout sample data with the daily temperature data
  
#PLOT SAMPLE SITES RECLASSIFIED WITH HUC8 LEVEL SUBBASINS AND RIVER SHAPEFILES
PA_basemap_rivers.ggplot <- PA_basemap.ggplot + #save a new map using the basemap as a starting point and adding in the following elements
  geom_sf(data=secondary_rivers.shp, aes(geometry=geometry), color="deepskyblue3", alpha=0.5, inherit.aes = F)+ #add flowlines of other rivers in the watershed, make them blue and semi-transparent so they blend in with the background a bit
  geom_sf(data=four_rivers.shp, aes(geometry=geometry, color=GNIS_NAME), linewidth=2, inherit.aes = F)+ #add the sample stream flowlines, color them with the colors used for samples sites throughout the manuscript, and make the lines thicker so they are even more visible 
  scale_color_discrete(na.translate=F)+ #don't include NAs in the color legend
  geom_sf(data=BKT_sites.sf, aes(geometry=geometry), color="black", size=2)+ #add sample point coordinates to the map and color them by original Region designations from the data sheets
  geom_sf_text(data=BKT_sites.sf, aes(geometry=geometry, label=Site), nudge_x = -0.02, nudge_y = 0.005)+ #add sample site name labels to the map, nudge them slightly up and to the left so they are not printed directly on top of the sample site points
  coord_sf(xlim=c(-78, -77.4), ylim=c(40.57, 40.87), crs=st_crs(3857), default_crs = st_crs(4326))+ #add a coordinate system to the plot, assume the input crs is 4326 for any elements without an available crs
  annotate("text", x=-77.88, y=40.79, label="State College")+ #add a label over the town of State College on the map
  annotation_north_arrow(height=unit(2.5, "cm"), width=unit(2.5, "cm"), pad_x=unit(16.5, "cm"), pad_y=unit(1.5, "cm"), style=north_arrow_minimal())+ #add a north arrow
  annotation_scale(pad_x=unit(12, "cm"), pad_y=unit(0.5, "cm"))+ #add a scale bar
  theme_classic(base_size = 16)+
  labs(color="Sampled Streams", title = "A") #edit the color legend title
  
inset.ggplot <- ggplot()+ #create an empty inset plot to show the extent of the study area shown in the main map
  borders(database = "world", fill="grey80", colour=NA, xlim=c(-85, -65), ylim=c(30, 50))+ #plot shapefile for world landmasses, fill with light grey and remove border lines
  borders(database = "state", fill="grey80", colour="grey70")+ #add state shapefiles to show US state borders in the inset map
  borders(database = "lakes", fill="grey60", colour=NA, xlim=c(-85, -65), ylim=c(30, 50))+ #get shapefile for world lakes, fill with dark grey and remove border lines
  geom_rect(aes(xmin=-78, xmax=-77.4, ymin=40.57, ymax=40.87), color="red", fill=NA)+ #add a bounding box showing the extent of the main map
  theme(panel.background=element_rect(fill="grey60"), panel.grid = element_line(color="grey70"), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_rect(fill=NA), plot.margin = margin(0,0,0,0, "cm"))+ #remove the axes and change the colors of the long/lat grid in the extent plot
  coord_sf(xlim=c(-85, -70), ylim=c(35, 45), crs=st_crs(3857), default_crs = st_crs(4326)) #change the crs of the extent plot so the projection matches the projection of the Pennsylvania basemap

BP_sampletemp.plot <- ggplot(data=daily_temps.df[daily_temps.df$Site=="Big Poe",], aes(x=Date, y=Temp.avg))+ #create ggplot with date on the x-axis, temperature on the y-axis, and set color to use for the line
  geom_line(color="#F8766D")+ #add a line object connecting daily average temperatures through time
  geom_point(data=trout_sample_data.df[trout_sample_data.df$Site=="Big Poe",], aes(x=Date, y=Temp.avg, fill=Temp.avg), color="#F8766D", size=3, pch=21)+ #add points for sample dates in this stream
  scale_fill_gradient(low="blue", high="red", limits=c(14.5, 21.5) , guide=F)+ #fill the sample date points with a color based on the sample date temperature
  geom_hline(yintercept = 20, linetype=2, alpha=0.6)+ #add a dotted horizontal line at 20ºC to indicate acute stress threshold
  theme_classic(base_size=16)+ #change the plot theme to make it look neater
  labs(y="Average Daily Temp (ºC)", title="B", subtitle = "Big Poe")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

EB_sampletemp.plot <- ggplot(data=daily_temps.df[daily_temps.df$Site=="East Branch",], aes(x=Date, y=Temp.avg))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_line(color="#7CAE00")+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_point(data=trout_sample_data.df[trout_sample_data.df$Site=="East Branch",], aes(x=Date, y=Temp.avg, fill=Temp.avg), color="#7CAE00", size=3, pch=21)+ #add points for sample dates in this stream
  scale_fill_gradient(low="blue", high="red", limits=c(14.5, 21.5) , guide=F)+ #fill the sample date points with a color based on the sample date temperature
  geom_hline(yintercept = 20, linetype=2, alpha=0.6)+ #add a dotted horizontal line at 20ºC to indicate acute stress threshold
  theme_classic(base_size=16)+ #change the plot theme to make it look neater
  labs(y=" ", title = "C", subtitle = "East Branch")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

SCA_sampletemp.plot <- ggplot(data=daily_temps.df[daily_temps.df$Site=="Shaver's Creek",], aes(x=Date, y=Temp.avg))+ #create ggplot with date on the x-axis, temperature on the y-axis, and set color to use for the line
  geom_line(color="#00BFC4")+ #add a line object connecting daily average temperatures through time
  geom_point(data=trout_sample_data.df[trout_sample_data.df$Site=="Shaver's Creek",], aes(x=Date, y=Temp.avg, fill=Temp.avg), color="#00BFC4", size=3, pch=21)+ #add points for sample dates in this stream
  scale_fill_gradient(low="blue", high="red", limits=c(14.5, 21.5) , guide=F)+ #fill the sample date points with a color based on the sample date temperature
  geom_hline(yintercept = 20, linetype=2, alpha=0.6)+ #add a dotted horizontal line at 20ºC to indicate acute stress threshold
  theme_classic(base_size=16)+ #change the plot theme to make it look neater
  labs(y=" ", title="D", subtitle="Shaver's Creek")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

SS_sampletemp.plot <- ggplot(data=daily_temps.df[daily_temps.df$Site=="Standing Stone",], aes(x=Date, y=Temp.avg))+ #create ggplot with date on the x-axis, temperature on the y-axis
  geom_line(color="#C77CFF")+ #add a line object connecting daily average temperatures through time, pick a color so it matches the stream color from the previous set of plots
  geom_point(data=trout_sample_data.df[trout_sample_data.df$Site=="Standing Stone",], aes(x=Date, y=Temp.avg, fill=Temp.avg), color="#C77CFF", size=3, pch=21)+ #add points for sample dates in this stream
  scale_fill_gradient(low="blue", high="red", limits=c(14.5, 21.5) , guide=F)+ #fill the sample date points with a color based on the sample date temperature
  geom_hline(yintercept = 20, linetype=2, alpha=0.6)+ #add a dotted horizontal line at 20ºC to indicate acute stress threshold
  theme_classic(base_size=16)+ #change the plot theme to make it look neater
  labs(y=" ", title = "E", subtitle="Standing Stone")+ #edit the y-axis label
  xlim(as.Date("2022-07-15"), as.Date("2022-08-22"))+ #set x-axis limits so it will be consistent across different plots and to remove plot area outside of sampling period
  ylim(14, 22) #add y-axis limits so they will be consistent across different plots

map_layout <- rbind(c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(2,3,4,5)) #set up multiplot layout (4 rows, 4 columns), map should take up all four columns of first 3 rows, temperature plots take up one column each in the bottom row
stacked_meantemp.plot <- grid.arrange(PA_basemap_rivers.ggplot, BP_sampletemp.plot, EB_sampletemp.plot, SCA_sampletemp.plot, SS_sampletemp.plot, layout_matrix=map_layout) #save an object with the map and temperature plots of the four rivers with sampling events indicated

jpeg("Figure1.jpeg", width=14, height=12, units="in", res=300) #open an 8x8" 300dpi jpeg to save the sample site plot
print(grid.arrange(PA_basemap_rivers.ggplot, BP_sampletemp.plot, EB_sampletemp.plot, SCA_sampletemp.plot, SS_sampletemp.plot, layout_matrix=map_layout)) #add the main map and temperature plots
print(inset.ggplot, vp=viewport(0.8, 0.8, width=0.2, height=0.2)) #add the extent inset plot in the upper right hand corner of the jpeg
dev.off() #stop writing to the jpeg

