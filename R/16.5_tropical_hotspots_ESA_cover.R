########################################################
#### GETTING ESA 2018 LCC MAP FOR TROPICAL HOTSPOTS ####
########################################################
## Getting packages, general paths and files ##
require(raster)
require(rgdal)
require(rgeos)
require(sp)

## Defining the general path to read/save files (tifs and shape files)
path = "E://ownCloud//W_GIS" # Renato's path

## Hotspots shapefiles
hotspots <- readOGR(dsn=paste(path,"//WO_Global_Biodiversity_Hotspots//hotspots_2016_1",sep=""),
                    layer="hotspots_2016_1")

trop.hotspots <- c("Atlantic Forest", "Caribbean Islands", 
                   "Coastal Forests of Eastern Africa", "East Melanesian Islands",
                   "Eastern Afromontane", "Guinean Forests of West Africa",
                   "Indo-Burma", "Madagascar and the Indian Ocean Islands",
                   "Mesoamerica", "New Caledonia", "Philippines", 
                   "Sundaland", "Tropical Andes",
                   "Tumbes-Choco-Magdalena", "Wallacea", "Western Ghats and Sri Lanka")
hotspots <- hotspots[hotspots@data$NAME %in% trop.hotspots,]
plot(hotspots)
  
## Tnc regions
# tnc <- readOGR(dsn=paste(path,"//Am_Lat_Biorregions_TNC",sep=""),layer="tnc_terr_ecoregions")
# tnc <- tnc[tnc@data$WWF_REALM2 %in% c("Afrotropic", "Neotropic"),]
# tnc <- tnc[tnc@data$ECO_NAME %in% 
#              c("Amazon-Orinoco-Southern Caribbean Mangroves",
#                "Central Congolian Lowland Forests",
#                "ChocÃ³-DariÃ©n Moist Forests"),]

#### ESA 2018 Land use map ####
require(ncdf4)
lu <- raster(paste(path,
                  "//WO_ESA_Land_Cover_map//v.2.1.1//C3S-LC-L4-LCCS-Map-300m-P1Y-2018-v2.1.1.nc",sep=""))
pj = crs(lu)

# clipping global map to hotspot's extents
for(i in 1:length(trop.hotspots)) {
  hot.i <- hotspots[hotspots@data$NAME %in% trop.hotspots[i],]
  proj4string(hot.i) = pj 
  hot.i1 <- extent(hot.i) - c(1,-1,1,-1)
  clip <- crop(lu, hot.i1)
  nome.arquivo <- paste(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//ESA-LCCS-Map-300m-2018-v2.1.1_",
                        gsub(" ","_",trop.hotspots[i]),
                        ".tif",sep="")
  writeRaster(clip, filename= nome.arquivo,
              format="GTiff",overwrite=TRUE)
}

# clipping hotspot's extents to hotspot's limits
list.rasters <- list.files(paste0(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//"),
                           full.names = TRUE)
for(i in 1:length(list.rasters)) {
  raster.i <- raster(list.rasters[i])
  pj.i <- crs(raster.i)
  
  hot.i <- gsub("ESA-LCCS-Map-300m-2018-v2.1.1_|\\.tif", "", 
                sapply(strsplit(list.rasters[i],"//"), tail, 1))
  hot.i <- gsub("_", " ", hot.i)
  hot.i1 <- hotspots[hotspots@data$NAME %in% hot.i,]
  proj4string(hot.i1) <- pj.i
  
  clip1 <- crop(raster.i, hot.i1)
  nome.arquivo <- paste(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//ESA-LCCS-Map-300m-2018-v2.1.1_",
                        gsub(" ","_", hot.i),
                        "_clipped.tif",sep="")
  writeRaster(clip1, filename= nome.arquivo,
              format="GTiff",overwrite=TRUE)
}

# masking hotspot's limits
list.rasters <- list.files(paste0(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//"),
                           full.names = TRUE)
list.rasters1 <- list.rasters[grepl("clipped", list.rasters)] 
for(i in 1:length(list.rasters1)) {
  raster.i <- raster(list.rasters1[i])
  pj.i <- crs(raster.i)
  
  hot.i <- gsub("ESA-LCCS-Map-300m-2018-v2.1.1_|_clipped\\.tif", "", 
                sapply(strsplit(list.rasters1[i],"//"), tail, 1))
  hot.i <- gsub("_", " ", hot.i)
  hot.i1 <- hotspots[hotspots@data$NAME %in% hot.i,]
  
  if ("outer limit" %in% hot.i1@data$Type) # removendo contornos não terrestres 
    hot.i1 <- hot.i1[!hot.i1@data$Type %in% "outer limit",]
  
  proj4string(hot.i1) <- pj.i
  
  clip2 <- mask(raster.i, hot.i1)
  gc()
  nome.arquivo <- paste(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//ESA-LCCS-Map-300m-2018-v2.1.1_",
                        gsub(" ","_", hot.i),
                        "_masked.tif",sep="")
  writeRaster(clip2, filename= nome.arquivo,
              format="GTiff",overwrite=TRUE)
}

#### GETTING THE VEGETATION COVER PER HOTSPOT ####
list.rasters <- list.files(paste0(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//"),
                           full.names = TRUE)
list.rasters2 <- list.rasters[grepl("masked", list.rasters)] 
result <- vector("list", length(list.rasters2))
nomes <- gsub("ESA-LCCS-Map-300m-2018-v2.1.1_|_masked\\.tif", "", 
              sapply(strsplit(list.rasters2,"//"), tail, 1))
names(result) <- nomes
classes = c(10, 11, 12, 20, 30, 40, 50, 60, 61, 62, 70,71, 72, 80, 81, 82, 90, 
            100, 110, 120, 121,122, 130, 140, 150, 151,152, 153, 160, 170, 180, 
            190, 200, 201,202, 210, 220) 
for(i in 1:length(list.rasters2)) {
  raster.i <- raster(list.rasters2[i])
  
  tmp <- table(factor(raster.i[], levels = classes))
  cellStats(!is.na(raster.i), sum)
  sum(area(raster.i, na.rm = TRUE)[])
  pixs <- sum(tmp)
  area <- cellStats(area(raster.i, na.rm = TRUE), sum, na.rm = TRUE)
  tmp1 <- round(100*tmp/pixs, 6)
  gc()
  result[[i]] <- c(area_km2 = area, n.pixels = pixs, tmp1)
}

# Loading the ESA Land Use legends and replacing the class names
tmp <- read.csv(paste(path,"//WO_ESA_Land_Cover_map//v2.0.7//ESACCI-LC-Legend.csv",sep=""), 
                as.is=TRUE)
for(i in 1:length(result)) names(result[[i]])[-c(1,2)] <- 
                              tmp$LegendTreeCo[match(names(result[[i]])[-c(1,2)], tmp$NB_LAB)] 
result1 <- do.call(rbind, result)
result2 <- cbind.data.frame(hotspot.region = gsub("_"," ",row.names(result1)),
                            area_km2 = result1[,1],
                            n.pixel = result1[,2],
                            t(rowsum(t(result1[,-c(1,2)]), colnames(result1[,-c(1,2)]))))

# Saving the results
write.csv(result2, "data//esa_2018_lc_per_hotspot.csv")
