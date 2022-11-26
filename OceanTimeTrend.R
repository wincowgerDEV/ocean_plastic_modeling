
#Libraries ----
library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(rgeos)
library(rgdal)
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(NADA)

#Functions ----
theme_gray_etal<- function(base_size = 12, bgcolor = NA){
    half_line <- base_size/2
    theme(
        line = element_line(colour = "black", size = rel(1.5), 
                            linetype = 1, lineend = "butt"), 
        rect = element_rect(fill = NA, colour = "black",
                            size = 0.5, linetype = 1),
        text = element_text(face = "plain",
                            colour = "black", size = base_size,
                            lineheight = 0.9,  hjust = 0.5,
                            vjust = 0.5, angle = 0, 
                            margin = margin(), debug = FALSE), 
        
        axis.line = element_blank(), 
        axis.text = element_text(size = rel(1.5), colour = "grey10"),
        axis.text.x = element_text(margin = margin(t = half_line/2), 
                                   vjust = 1), 
        axis.text.y = element_text(margin = margin(r = half_line/2),
                                   hjust = 1),
        axis.ticks = element_line(colour = "black", size=1), 
        axis.ticks.length = unit(half_line*0.75, "pt"), 
        axis.title = element_text(size = rel(1.5), colour = "black"),
        axis.title.x = element_text(margin = margin(t = half_line*5,
                                                    b = half_line)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(r = half_line*5,
                                                    l = half_line)),
        
        legend.background = element_rect(colour = NA), 
        legend.key = element_rect(colour = NA),
        legend.key.size = unit(2, "lines"), 
        legend.key.height = NULL,
        legend.key.width = NULL, 
        legend.text = element_text(size = rel(1)),
        legend.text.align = NULL,
        legend.title = element_text(size = rel(1)), 
        legend.title.align = NULL, 
        legend.position = "right", 
        legend.direction = NULL,
        legend.justification = "center", 
        legend.box = NULL, 
        
        panel.background = element_rect(fill=bgcolor,colour = "black", size = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(half_line, "pt"), panel.margin.x = NULL, 
        panel.spacing.y = NULL, panel.ontop = FALSE, 
        
        #Facet Labels
        strip.background = element_blank(),
        strip.text = element_text(face="bold",colour = "black", size = rel(1.5)),
        strip.text.x = element_text(margin = margin(t = half_line,
                                                    b = half_line)), 
        strip.text.y = element_text(angle = 0, 
                                    margin = margin(l = half_line, 
                                                    r = half_line)),
        strip.switch.pad.grid = unit(5, "lines"),
        strip.switch.pad.wrap = unit(5, "lines"), 
        
        
        plot.background = element_rect(colour = rgb(119,136,153, max = 255)), 
        plot.title = element_text(size = rel(1.5), 
                                  margin = margin(b = half_line * 1.2)),
        plot.margin = margin(4*half_line, 4*half_line, 4*half_line, 4*half_line),
        complete = TRUE)
}

#General Data Cleanup ----
VanSeb <- as.matrix(read.csv("ModelGrid/vansebillemodel_abundance.csv", header = F))
Ar <- as.matrix(read.csv("ModelGrid/vansebillemodel_area.csv", header = F))
AreaModel <- as.vector(Ar)
OceanLand <- as.matrix(read.csv("ModelGrid/mask.csv", header = F))

LocationWindCorrected <- read.csv("TrawlData/ProcessedFiles/LocationsCleanedWithWind.csv")
corrected <- unique(LocationWindCorrected$Source[LocationWindCorrected$Corrected == "Yes"])

Locations_full_test <- LocationWindCorrected
Locations_full_test$LonErik <- ifelse(Locations_full_test$Lon < 0, Locations_full_test$Lon + 360, Locations_full_test$Lon)
Locations_full_test$LonErik <- ifelse(Locations_full_test$Lon < 0, Locations_full_test$Lon + 360, Locations_full_test$Lon)

coordinates(Locations_full_test) <- ~LonErik + Lat

#with wind

change_correction <- c("M. Eriksen", "F.Galgani", "H.Carson", "J. Reisser", "C. Moore")
change_wind <- c("M. Eriksen", 
                 "Wong et al. 1974", 
                 "F.Galgani", 
                 "Shaw & Mapes  1979",
                 "Carpenter 1972",
                 "Shaw 1975",
                 "H.Carson", 
                 "J. Reisser", 
                 "C. Moore")

# infer wind values 
set.seed(3984)
#LocationWindCorrected$wind[is.na(LocationWindCorrected$wind)] <- sample(LocationWindCorrected$wind[!is.na(LocationWindCorrected$wind)], sum(is.na(LocationWindCorrected$wind)), replace = T)

Locations <- LocationWindCorrected %>%
  mutate(Corrected = ifelse(Source %in% change_correction, "No", Corrected)) %>%
  mutate(CorrectedConcentration = ifelse(Source %in% change_wind, WindStandard(km2, wind), CorrectedConcentration)) %>%
  dplyr::filter(Corrected == "No") %>%
  dplyr::mutate(km2 = CorrectedConcentration) %>%
  dplyr::filter(!is.na(km2)) %>%
  dplyr::mutate(Latcopy = Lat) %>%
  filter(Mesh.Size > 50)

anti_locations <- LocationWindCorrected %>%
    anti_join(Locations %>% select(X.1))

Locations$LonErik<- ifelse(Locations$Lon < 0, Locations$Lon + 360, Locations$Lon)
Locations$LatErik <- Locations$Lat
coordinates(Locations) <- ~LonErik + LatErik

#Create a raster of Van Sebille model----
Van <- raster(VanSeb, xmn = 0, xmx= 360, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(Van)
#plot(Locations_full_test, add = T)
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

LocationFram$Date <- as.numeric(as.Date(with(LocationFram, paste(Year, Month, Day, sep="-")), "%Y-%m-%d")) - as.numeric(as.Date("1979-1-15")) #Cycle through the date creating new raster each time. 
LocationFram$DateFormatted <- (as.Date(LocationFram$Date,origin = "1979-1-15"))

LocationFram <- filter(LocationFram, !is.na(Date)) %>%
  filter(!is.na(Basin)) %>%
  mutate(Basin = as.factor(as.character(Basin)))

LocationFramCountYear <- LocationFram %>%
  mutate(YearRange = cut(Year, breaks = seq(1980, 2020, by = 5))) %>%
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

# Predict the nondetects ----
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

#Model Using just a GAM ----
time <- gam(fitkm2 ~ Basin + OceanScaled + s(Date), data = Dataset, family = tw())
summary(time)
plot(time)
qq.gam(time, type = "response")

Basin = as.character(OceanLandDFScaled$Basin)
Basins = ifelse(Basin == 1, "North Pacific", ifelse(Basin == 2, "South Pacific", ifelse(Basin == 3, "North Atlantic", ifelse(Basin == 4, "South Atlantic", ifelse(Basin == 5, "Indian", ifelse(Basin == 6, "Mediterranean", NA))))))

global_basins <- expand.grid(paste(Basins, OceanLandDFScaled$scaled, sep = "_"), Date = seq(0,14752, by = 365)) %>%
                            mutate(Basin = gsub("_.*", "", Var1), OceanScaled = as.numeric(gsub(".*_", "", Var1))) %>%
                            select(-Var1)

global_basins$PredictedVals <- predict.gam(time, global_basins, type = "response")

global_basins_2 <- global_basins %>%
    group_by(Date) %>%
    summarise(sum = sum(PredictedVals, na.rm = T))
#This approach definitely has some issues, very large numbers. Likely due to poor qq plot fits. 

summary_uncertainties <- Dataset %>%
    group_by(Year, Basin) %>%
    summarize(count = n(), cv = sd(fitkm2)/mean(fitkm2))

write.csv(summary_uncertainties, "uncertainties.csv")

ggplot(summary_uncertainties) + geom_point(aes(x = Year, y = count, color = Basin, shape = Basin), size = 5, alpha = 0.75) + scale_color_viridis_d() + scale_y_log10() + theme_gray_etal()
ggplot(summary_uncertainties) + geom_point(aes(x = Year, y = cv, color = Basin, shape = Basin), size = 5, alpha = 0.75) + scale_y_continuous(limits = c(0, 15)) + scale_color_viridis_d()  + theme_gray_etal()
#ggplot(summary_uncertainties) + geom_histogram(aes(x = cv, fill = Basin, shape = Basin), size = 5, alpha = 0.75) + scale_y_continuous(limits = c(0, 15)) + scale_color_viridis_d()  + theme_gray_etal()

summary(summary_uncertainties$cv)

ggplot(Dataset, aes(x = DateFormatted, y = fitkm2)) + geom_smooth() + geom_rug(sides="b", alpha = 0.1) + facet_wrap(.~Basin, scales = "free") + theme_classic()

gamsmoothresiduals = gam(fitlogkm2 ~ OceanScaled + as.factor(Basin), data=Dataset)
summary(gamsmoothresiduals)
AIC(gamsmoothresiduals)
plot(gamsmoothresiduals)
qq.gam(gamsmoothresiduals, type = "response")

Dataset$PredictedVals <- predict.gam(gamsmoothresiduals, Dataset, type = "response")

Dataset$Residuals <- Dataset$fitlogkm2 - Dataset$PredictedVals

hist(Dataset$Residuals)

write.csv(Dataset, "TrawlData/ProcessedFiles/model_dataset.csv")

#Visualize difference in residuals
#ggplot() + geom_point(aes(x = residuals.gam(gamsmoothresiduals, type = "deviance"), y = residuals.gam(gamsmoothresiduals, type = "response"))) + theme_classic() + scale_y_log10()

gamsmoothresidualsdate = gam(Residuals ~ s(Date), data=Dataset) 
summary(gamsmoothresidualsdate)
AIC(gamsmoothresidualsdate)
plot(gamsmoothresidualsdate)
qq.gam(gamsmoothresidualsdate, type = "response")

dftimetrend <- data.frame(Date = seq(0,14752))
dftimetrend$DateFormatted <- (as.Date(dftimetrend$Date,origin = "1979-1-15"))
dftimetrend$prediction <- predict.gam(gamsmoothresidualsdate, dftimetrend, type = "response")
predictionse <- predict.gam(gamsmoothresidualsdate, dftimetrend, type = "response", se.fit = T)
dftimetrend$lower <- dftimetrend$prediction - (2* predictionse$se.fit)
dftimetrend$upper <- dftimetrend$prediction + (2* predictionse$se.fit)

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
  scale_y_log10(breaks = c(8 * 10^12, evennums * 10^13, evennums * 10^14), labels = c(8 , evennums * 10, evennums * 10^2)) +
  #scale_y_log10(breaks = c( 5*10^12, 10^13, 5*10^13, 10^14, 5*10^14), labels = c(5*10^12, 10^13, 5*10^13, 10^14, 5*10^14), limits = c(5*10^12, 5*10^14)) + 
  labs(x = "Year", y = "Trillions of Particles or Tens of Thousands of Metric Tonnes")

write.csv(dftimetrend, "timetrend.csv")

yearly_mean <- dftimetrend %>%
    mutate(year = format(DateFormatted, format = "%Y")) %>%
    group_by(year) %>%
    summarise(mean_total_count_trillion_particles = mean(oceantotalcount)/10^12, mean_total_mass_million_metric_tonnes = mean(oceantotalmassmetrictonnes)/10^6)

write.csv(yearly_mean, "yearly_mean.csv")


