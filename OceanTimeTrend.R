
#Libraries ----
library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)
library(dismo)
library(rgdal)
library(mgcv)
library(ggplot2)
library(dplyr)
library(leaflet)
library(viridis)
library(gridExtra)
library(purrr)
library(tidyr)
library(magick)
library(plotly)
library(NADA)
library(EnvStats)


#Functions----

BootMedian <- function(data) {
  B <- 10000
  ave <- numeric(B)
  n = length(data)
  
  set.seed(34345)
  for (i in 1:B) {
    boot <- sample(1:n, size=n, replace = TRUE)
    ave[i] <- mean(data[boot])
  }
  return(quantile(ave, c(0.025, 0.5, 0.975), na.rm = T))
}

getmode <- function(v) {
  uniqv <- unique(v[!is.nan(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmean <- function(v) {
  mean(v, na.rm = T)
}

getgeomean <- function(x) {
  10^mean(log10(x + 1))
}


#General Data Cleanup ----
VanSeb <- as.matrix(read.csv("ModelGrid/vansebillemodel_abundance.csv", header = F))
Ar <- as.matrix(read.csv("ModelGrid/vansebillemodel_area.csv", header = F))
AreaModel <- as.vector(Ar)
OceanLand <- as.matrix(read.csv("ModelGrid/mask.csv", header = F))

LocationWindCorrected <- read.csv("TrawlData/ProcessedFiles/LocationsCleanedWithWind.csv")
corrected <- unique(LocationWindCorrected$Source[LocationWindCorrected$Corrected == "Yes"])

#with wind
Locations <- LocationWindCorrected %>%
  dplyr::filter(Corrected == "No") %>%
  dplyr::mutate(km2 = CorrectedConcentration) %>%
  dplyr::filter(!is.na(km2)) %>%
  dplyr::mutate(Latcopy = Lat) %>%
  filter(Mesh.Size > 50)

Locations$LonErik<- ifelse(Locations$Lon < 0, Locations$Lon + 360, Locations$Lon)
coordinates(Locations) <- ~LonErik + Lat

#Create a raster of Eriks model----
Van <- raster(VanSeb, xmn = 0, xmx= 360, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(Van)
plot(Locations, add = T) #Check to make sure the datasets are aligned properly.

OceanLandDF <- as.data.frame(as.vector(OceanLand))
VanDF <- as.data.frame(as.vector(VanSeb))
OceanLandDF$Van <- VanDF$`as.vector(VanSeb)`
Area <- raster(Ar, xmn = 0, xmx= 360, ymn=-90, ymx=90)

OceanBasinVector <- as.numeric(as.character(OceanLandDF$`as.vector(OceanLand)`))
OceanBasinVector <- ifelse(is.nan(OceanBasinVector), NA, OceanBasinVector)
dim(OceanBasinVector)<- c(181,361)
OceanBasin <- raster(OceanBasinVector, xmn = 0, xmx= 360, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

OceanLandDFScaled <- OceanLandDF %>%
  mutate(Basin = as.factor(`as.vector(OceanLand)`)) %>%
  #dplyr::group_by(Basin) %>% ----
  mutate(scaled = (OceanLandDF$Van-mean(Van, na.rm = T))/sd(Van, na.rm = T))

OceanScaled <- OceanLandDFScaled$scaled#-mean(log(OceanLandDF$Van+1), na.rm = T))/sd(log(OceanLandDF$Van+1), na.rm = T)
dim(OceanScaled)<- c(181,361)
OceanScaled <- raster(OceanScaled, xmn = 0, xmx= 360, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

plot(OceanScaled)
plot(OceanBasin)
plot(Locations, add = T) #Check to make sure the datasets are aligned properly.

#Extract the ocean scaled and basin values. ----
Locations$Basin <- raster::extract(OceanBasin, Locations)
Locations$OceanScaled <- raster::extract(OceanScaled, Locations)

LocationFram <- Locations@data
LocationFram <- LocationFram %>%
  mutate(Day = ifelse(is.na(Day), 15, Day))

LocationFram$Date <- as.numeric(as.Date(with(LocationFram, paste(Year, Month, Day, sep="-")), "%Y-%m-%d")) - as.numeric(as.Date("1979-01-15")) #Cycle through the date creating new raster each time. 
LocationFram$DateFormatted <- (as.Date(LocationFram$Date,origin = "1979-01-15"))

LocationFram <- filter(LocationFram, !is.na(Date)) %>%
  filter(!is.na(Basin)) %>%
  mutate(Basin = as.factor(as.character(Basin)))

LocationFramCountYear <- LocationFram %>%
  mutate(YearRange = cut(Year, breaks = seq(1970, 2020, by = 5))) %>%
  group_by(YearRange, Basin) %>%
  dplyr::summarise(datapoints = n())%>%
  ungroup() %>%
  pivot_wider(names_from = Basin, values_from = datapoints)

names(LocationFramCountYear) <- c("YearRange", "South Atlantic", "North Pacific", "North Atlantic", "South Pacific",  "Indian", "Mediterranean")

write.csv(LocationFramCountYear, "AnnualDataPoints.csv")

#Model Time Trend with Ocean, Basin, and Date using residuals ----

LocationFramBasinSet <- LocationFram %>%
  mutate(Basin = ifelse(Basin == 1, "North Pacific", ifelse(Basin == 2, "South Pacific", ifelse(Basin == 3, "North Atlantic", ifelse(Basin == 4, "South Atlantic", ifelse(Basin == 5, "Indian", ifelse(Basin == 6, "Mediterranean", NA))))))) %>%
  mutate(censored = km2 == 0) 

hist(log(LocationFramBasinSet$km2 + 1))
summary(LocationFramBasinSet$censored)

fit = cenros(LocationFramBasinSet$km2, LocationFramBasinSet$censored)

hist(log(fit$modeled))

set.seed(211)

fittedvalues <- sample(fit$modeled[fit$censored], length(fit$modeled[fit$censored]), replace = F)

LocationFramBasinSet <- LocationFramBasinSet %>%
  mutate(fitkm2 = ifelse(censored, fittedvalues, km2)) %>%
  mutate(fitlogkm2 = log(fitkm2))

hist(log(LocationFramBasinSet$fitkm2))

NorthPacific <- LocationFramBasinSet %>%
  dplyr::filter(Basin == "North Pacific")

Dataset <- LocationFramBasinSet

hist(Dataset$Date)
hist(Dataset$fitlogkm2)
hist(Dataset$OceanScaled)

summary(Dataset)

#Dataset %>%
#  filter(Mesh.Size < 300)

#time to concentration
time <- gam(km2 ~ as.factor(Basin) + s(Date), data = Dataset, family = tw())
summary(time)
plot(time)

time_bias <- gam(OceanScaled ~ s(Date), data = Dataset, family = tw())
summary(time_bias)
plot(time_bias)

ggplot(Dataset, aes(x = DateFormatted, y = OceanScaled)) + geom_smooth() + geom_rug(sides="b", alpha = 0.1) + facet_wrap(.~Basin, scales = "free") + theme_classic()


hist(log(Dataset$km2+1))
#Using this currently
gamsmoothresiduals = gam(fitlogkm2 ~ OceanScaled + as.factor(Basin), data=Dataset)
summary(gamsmoothresiduals)
AIC(gamsmoothresiduals)
plot(gamsmoothresiduals)
qq.gam(gamsmoothresiduals, type = "response")


Dataset$PredictedVals <- predict.gam(gamsmoothresiduals, Dataset, type = "response")
#Dataset$PredictedVals <- predict.lm(gamsmoothresiduals, Dataset, type = "response")

Dataset$Residuals <- Dataset$fitlogkm2 - Dataset$PredictedVals

hist(Dataset$Residuals)
#Dataset$Residualslog <- log(Dataset$Residuals)

#Dataset$Residuals2 <- residuals.gam(gamsmoothresiduals, type = "deviance")

#Visualize difference in residuals
#ggplot() + geom_point(aes(x = residuals.gam(gamsmoothresiduals, type = "deviance"), y = residuals.gam(gamsmoothresiduals, type = "response"))) + theme_classic() + scale_y_log10()

gamsmoothresidualsdate = gam(Residuals ~ s(Date), data=Dataset) 
summary(gamsmoothresidualsdate)
AIC(gamsmoothresidualsdate)
plot(gamsmoothresidualsdate)
qq.gam(gamsmoothresidualsdate, type = "response")

dftimetrend <- data.frame(Date = seq(0,14752))
dftimetrend$DateFormatted <- (as.Date(dftimetrend$Date,origin = "1979-01-15"))
dftimetrend$prediction <- predict.gam(gamsmoothresidualsdate, dftimetrend, type = "response")
predictionse <- predict.gam(gamsmoothresidualsdate, dftimetrend, type = "response", se.fit = T)
dftimetrend$lower <- dftimetrend$prediction - (2* predictionse$se.fit)
dftimetrend$upper <- dftimetrend$prediction + (2* predictionse$se.fit)

#dftimetrend$se <- gamsmoothresidualsdate$
mean_count_concentration_km2 = mean(Dataset$fitkm2)
mean_mass_concentration_g_km2 = mean(Dataset$fitkm2) *  1.36 * 10^-2 
total_ocean_area_km2 = 361900000

dftimetrend$averageconcentration_count_km2 = mean_count_concentration_km2 * exp(dftimetrend$prediction)
dftimetrend$averageconcentrationlower_count_km2 = mean_count_concentration_km2 * exp(dftimetrend$lower)
dftimetrend$averageconcentrationupper_count_km2 = mean_count_concentration_km2 * exp(dftimetrend$upper)

dftimetrend$averagemassconcentration_g_km2 =  mean_mass_concentration_g_km2 * exp(dftimetrend$prediction)
dftimetrend$averagemassconcentration_g_km2_lower =  mean_mass_concentration_g_km2 * exp(dftimetrend$lower)
dftimetrend$averagemassconcentration_g_km2_upper =  mean_mass_concentration_g_km2 * exp(dftimetrend$upper)

dftimetrend$oceantotalcount =  total_ocean_area_km2 * dftimetrend$averageconcentration_count_km2 
dftimetrend$oceantotalcountlower = total_ocean_area_km2 * dftimetrend$averageconcentrationlower_count_km2
dftimetrend$oceantotalcountupper = total_ocean_area_km2* dftimetrend$averageconcentrationupper_count_km2

dftimetrend$oceantotalmassmetrictonnes =  total_ocean_area_km2 * dftimetrend$averagemassconcentration_g_km2 * 10^-6
dftimetrend$oceantotalmassmetrictonnes_lower =  total_ocean_area_km2 * dftimetrend$averagemassconcentration_g_km2_lower * 10^-6
dftimetrend$oceantotalmassmetrictonnes_upper =  total_ocean_area_km2 * dftimetrend$averagemassconcentration_g_km2_upper * 10^-6

#mean(c(1, 10, 100)*40)
#ggplot(dftimetrend) + geom_line(aes(x = Date, y = averageconcentration)) + scale_y_log10()

ggplot(dftimetrend) + geom_line(aes(x = DateFormatted, y = oceantotalmassmetrictonnes)) + theme_gray_etal()

ggplot(dftimetrend) + 
  geom_line(aes(x = DateFormatted, y = exp(prediction)), size = 1) + 
  geom_line(aes(x = DateFormatted, y = exp(lower)), size = 0.5) + 
  geom_line(aes(x = DateFormatted, y = exp(upper)), size = 0.5) + 
  theme_gray_etal() + 
  scale_y_log10(breaks = seq(1, 13, by = 2), labels = seq(1, 13, by = 2)) + 
  labs(x = "Year", y = "Multiple of Average Conditions")

evennums <- c(2,4,6,8)

ggplot(dftimetrend) + 
  geom_line(aes(x = DateFormatted, y = oceantotalcount), size = 1) + 
  geom_line(aes(x = DateFormatted, y = oceantotalcountlower), size = 0.1) + 
  geom_line(aes(x = DateFormatted, y = oceantotalcountupper), size = 0.1) + 
  theme_gray_etal() + 
  scale_y_log10(breaks = c(8 * 10^12, evennums * 10^13, evennums * 10^14), labels = c(8 , evennums * 10, evennums * 10^2), limits = c(8*10^12, 8*10^14)) +
  #scale_y_log10(breaks = c( 5*10^12, 10^13, 5*10^13, 10^14, 5*10^14), labels = c(5*10^12, 10^13, 5*10^13, 10^14, 5*10^14), limits = c(5*10^12, 5*10^14)) + 
  labs(x = "Year", y = "Trillions of Particles or Tens of Thousands of Metric Tonnes")

write.csv(dftimetrend, "timetrend.csv")

yearly_mean <- dftimetrend %>%
    mutate(year = format(DateFormatted, format = "%Y")) %>%
    group_by(year) %>%
    summarise(mean_total_count_trillion_particles = mean(oceantotalcount)/10^12, mean_total_mass_million_metric_tonnes = mean(oceantotalmassmetrictonnes)/10^6)

write.csv(yearly_mean, "yearly_mean.csv")

#Map with residuals ----
OceanRankedModel <- data.frame(OceanScaled = as.vector(OceanScaled), Area = as.vector(Area), Basin = as.vector(OceanBasin)) %>%
  mutate(Basin = ifelse(Basin == 1, "North Pacific", ifelse(Basin == 2, "South Pacific", ifelse(Basin == 3, "North Atlantic", ifelse(Basin == 4, "South Atlantic", ifelse(Basin == 5, "Indian", ifelse(Basin == 6, "Mediterranean", NA)))))))

OceanRankedModel$Prediction_Grid <- predict.gam(gamsmoothresiduals, OceanRankedModel, type = "response")

Days = as.numeric(seq(as.Date("1980/1/1"), as.Date("2020/1/1"), "year")) - as.numeric(as.Date("1979-01-15"))

OceanRankedModel$Date <- as.numeric(as.Date("2019/1/1")) - as.numeric(as.Date("1979-01-15"))
OceanRankedModel$Residual <- predict.gam(gamsmoothresidualsdate, OceanRankedModel, type="response")
OceanRankedModel$predicted_concentration_num_km2 <- exp(OceanRankedModel$Prediction + OceanRankedModel$Residual)

OceanPrediction <- OceanRankedModel$predicted_concentration_num_km2
dim(OceanPrediction)<- c(361, 181) #When this gets flipped, the map looks right. 
OceanPrediction <- t(OceanPrediction)
#OceanPrediction <- t(OceanPrediction)
OceanPredictionRaster <- raster(OceanPrediction, xmn = 1, xmx= 361, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#writeRaster(OceanPredictionRaster, paste(WD,"/", Y, ".tif", sep=""), format="GTiff", overwrite=TRUE)
plot(OceanPredictionRaster)

raster <- as(OceanPredictionRaster, "SpatialPixelsDataFrame")
raster <- as.data.frame(raster)
colnames(raster) <- c("value", "x", "y")
raster$rank <- rank(raster$value)/length(raster$value)
max(raster$rank)

Plot <- ggplot(data=raster) + 
  geom_tile(aes(x=x,y=y,fill=rank)) + 
  coord_equal() +
  scale_fill_distiller(palette = "Spectral", breaks = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(raster$value, c(0, 0.25, 0.5, 0.75, 1)), 0) ) +
  #scale_fill_brewer(type = "div", palette = "RdYlBu",  drop = F) +
  theme_gray_etal() + 
  theme(panel.background = element_rect(fill = "darkslategray")) + 
  labs(x = "Longitude", y = "Latitude", fill = "Concentration (log10(#/km2))", title = "2019") #+ annotation_custom(ggplotGrob(YearPlot), xmin = 50, ymin = 30, xmax = 120, ymax = 70)

