# Code and Raw Data for Manuscript 

"A growing smog of more than 145 trillion plastic particles is afloat in the world’s oceans – urgent solutions required."

## Code
- OceanTimeTrend.R: Creates the ocean plastic time trend model. 
- WindAnalysis.R: Conducts wind correction for the ocean plastic observations. 

## Full ocean observation dataset for modeling is [HERE](https://github.com/wincowgerDEV/ocean_plastic_modeling/blob/main/TrawlData/ProcessedFiles/LocationsCleanedWithWind.csv)
Metadata:
- Day: Day of the month
- Month: Month of the year
- Year: Year
- Latitude: WGS 84 decimal degrees
- Longitude: WGS 84 decimal degrees
- X..km2: Total count concentration per square kilometer
- Total.Wieght..g..km.2: total mass in grams per square kilometer. 
- Source: The source the data is from 
- Mesh.Size: net size in microns
- Notes: Any notes I made about the data. 
- CutOff: cutoff particle size in microns between microplastic and macroplastic
- TopLim: upper particle size limit in microns. 
- Microplastics...km.2: count concentration of microplastics per square kilometer
- Macroplastics....km.2: count concentration of macroplastics per square kilometer
- Microplastics.Weight..g...km.2: mass concentration of microplastics grams per square kilometer
- Macroplastics.Weight..g...km.2: mass concentration of microplastics grams per square kilometer
- Ship: ship name
- TowAreaM2: area the net covered if available in square meters. 
- CountPlasticFragments: total count of plastic fragments
- CountPlasticPellets: total count of plastic pellets
- CountTotalPlastics: total count of all plastics
- FragmentDensityNoSqKm: concentration of plastic fragments in number per square kilometer. 
- PelletDensityNoSqKm: concentration of plastic pellets in number per square kilometer. 
- v10Mps: wind speed in meters per second if provided. 

## [TrawlData/RawFiles](https://github.com/wincowgerDEV/ocean_plastic_modeling/tree/main/TrawlData/RawFiles) holds the raw data from all datasets that were merged to build the dataset used in this study. 

## [WindMonthly](https://github.com/wincowgerDEV/ocean_plastic_modeling/tree/main/WindMonthly) holds the monthly wind raster grid used to assess historic wind conditions. 

## [ModelGrid](https://github.com/wincowgerDEV/ocean_plastic_modeling/tree/main/ModelGrid) holds the oceanographic model used as a reference for expected concentrations. 

