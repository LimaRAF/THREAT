#### Loading packages ####
require(raster)
require(sf)
require(foreach)
require(doSNOW)
require(snow)

### Defining the general path to read/save files (tifs and shape files) ###
setwd("~/Renato")

## EOO shapefiles
EOO.poly <- readRDS("EOO.poly.rds")

## ESA land cover
hab.map <- raster::stack("hab_map.tif")

## LandUse classes
hab.class <- readRDS("hab.class.rds")

## Sourcing the function
source("EOO.habitat.R")

## Running the function
toto <- EOO.habitat(EOO.poly[1:5], hab.map, hab.class, years = 23, 
                    parallel = TRUE, NbeCores = 7)
saveRDS(toto, "EOO_habitat_loss.rds")
