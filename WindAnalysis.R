#Code Works!!!

#Libraries ----
library(data.table)
library(dplyr)
library(ggplot2)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) 


#Wind Correction ----
#Function conducts transformation to correct for wind in the observed concentrations.
WindStandard <- function(Concentration,Wind) {
    Rise = 0.02
    Depth = 0.2
    WaveHeight <- ifelse(Wind == 0, 0, ifelse(Wind > 0 & Wind < 2, 0.1, ifelse(Wind >= 2 & Wind < 3.5, 0.25,ifelse(Wind >= 3.5 & Wind < 5, 0.8,ifelse(Wind >= 5 & Wind < 8.5, 1.25,ifelse(Wind >= 8.5 & Wind < 11, 2.25,NA)))))) #Based on beaufort scale
    k = 0.4
    Uw <- ((0.0012 *1.22* Wind^2)/1030)^(1/2) #Validated with Kukulka pretty much exactly the same value. This is correct now in m/s. #0.00012 * FakeWind from Reiser #Frictional Velocity of water from Pugh, #Wind Stress (https://marine.rutgers.edu/dmcs/ms501/2004/Notes/Wilkin20041014.htm)
    Ao <- 1.5 * Uw * k * WaveHeight
    DimCon <- Concentration / (1-exp(-1*Depth*Rise*(1/Ao)))
    print(DimCon)
}
#A source that has wind speed.
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form
#There are two variables listed, one for the x and one for the y coordinate of velocity. 

#Also this source allows download
#https://rda.ucar.edu/#!lfd?nb=y&b=proj&v=ECMWF%20ERA%2015%20Reanalysis

#This is data we are currently using.
#https://rda.ucar.edu/datasets/ds627.1/

#May want to look into this one to correct older datasets.
#https://rda.ucar.edu/datasets/ds626.0/

#Need to get the monthly wind data since the other wind data doesn't cover the whole dataset. 
#https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2
#Catalog https://www.ncei.noaa.gov/thredds/catalog.html
#https://www.ncei.noaa.gov/thredds/catalog/cfs_reanl_ts/201008/catalog.html?dataset=cfs_reanl_ts/201008/wnd10m.l.gdas.201008.grb2
#https://www.ncei.noaa.gov/thredds/catalog/model-cfs-allfiles/cfsv2_analysis_timeseries/2012/201212/catalog.html
#https://www.ncei.noaa.gov/thredds/catalog/cfs_reanl_ts/catalog.html

#Grab Wind Data from URLs----

#URL <- "https://www.ncei.noaa.gov/thredds/fileServer/uv/daily/2000s/uv" #Partial URL to Grab the data from
#TrawlData <- read.csv("TrawlData/ProcessedFiles/FinalMerge.csv") #Lat Long Points and date
#Date <- sort(unique(paste(TrawlData$Year,ifelse(TrawlData$Month < 10, paste("0", as.character(TrawlData$Month), sep=""), as.character(TrawlData$Month)),ifelse(TrawlData$Day < 10, paste("0", as.character(TrawlData$Day), sep=""), as.character(TrawlData$Day)), sep = ""))) #The rest of the URL
#DateTwoThou <- as.numeric(Date)[as.numeric(Date) > 20150915] #Cropping the rest of the URL
#DateTwoThou <- DateTwoThou[!is.na(DateTwoThou)]
#DateTwoThou <- DateTwoThou[2:2554] #cropping the url more
#DateTwoThou <- DateTwoThou[2:10]
#Folder <- "WindMonthly/" #The folder the downloaded data will go to
#for(Day in DateTwoThou ){
#  download.file(url = paste(URL, Day, "rt.nc", sep = ""), destfile = paste(Folder, Day, ".nc", sep = ""), mode = "wb" )
#}

###ScrapeWindDataFrom All Survey Points ----

Location <- read.csv("TrawlData/ProcessedFiles/FinalMerge.csv")

correctedsources <- c("C. Moore", "F.Galgani", "H.Carson", "J. Reisser", "M. Eriksen", "Cozar 2015", "Cozar 2017")

LocationDF <- Location %>%
  mutate(Mesh.Size == as.numeric(Mesh.Size)) %>%
  dplyr::filter(Mesh.Size < 4000 & Mesh.Size > 0) %>%
  dplyr::filter(is.na(TopLim)| TopLim > 100000) %>%
  #dplyr::select(Latitude, Longitude, Month, Day, Year, X..km2, Microplastics...km.2, Macroplastics....km.2, Source, Notes, v10Mps) %>%
  distinct() %>%
  mutate(Corrected = ifelse(Source %in% correctedsources, "Yes", "No")) %>%
  rename(Lon = Longitude, Lat = Latitude, km2 = X..km2, Micro = Microplastics...km.2, Macro = Macroplastics....km.2, MeasuredWind = v10Mps)

# Extract all wind data from historic record ----
# Takes a long time, could be optimized. 
for(row in 1:nrow(LocationDF)) {
  if(LocationDF[row, "Year"] < 1979) next
  wind <- readGDAL(paste("WindMonthly/ei.moda.an.sfc.regn128sc.", LocationDF[row, "Year"], ifelse(LocationDF[row, "Month"] > 9, LocationDF[row, "Month"], paste("0", LocationDF[row, "Month"], sep = "")),  "0100.cowger407132", sep = ""))
  wind <- raster(wind)
  LocationDF[row, "wind"] <- raster::extract(wind, Locations[row,], method='simple')
}

#Tested to make sure it didn't matter if the projection was there. It doesn't.
#wind <- readGDAL(paste("G:/My Drive/GrayLab/Plastics/Articles Publish/Active/OceanTimeTrend/WindMonthly/ei.moda.an.sfc.regn128sc.", LocationDF[100, "Year"], ifelse(LocationDF[100, "Month"] > 9, LocationDF[100, "Month"], paste("0", LocationDF[100, "Month"], sep = "")),  "0100.cowger407132", sep = ""))
#try<-  extract(wind, Locations[100,], method='simple')
#Try<- extract(wind, LocationTest[100,], method='simple')

# Conduct correction with the new wind data ----
LocationDF$CorrectedConcentration <- ifelse(LocationDF$Corrected == "Yes", LocationDF$km2, WindStandard(LocationDF$km2, LocationDF$wind))

# Write the corrected dataset ----
write.csv(LocationDF, "TrawlData/ProcessedFiles/LocationsCleanedWithWind.csv")
