########################################################
#### GETTING ESA 2018 LCC MAP FOR TROPICAL HOTSPOTS ####
########################################################
rm(list = ls())
gc()

## Getting packages, general paths and files ##
require(raster)
require(rgdal)
require(rgeos)
require(cleangeo)
require(sp)
require(ncdf4)

## Hotspots shapefiles
hotspots <- readOGR(dsn="data/data-raw/hotspots_2016_1", layer="hotspots_2016_1")

trop.hotspots <- c("Atlantic Forest", "Caribbean Islands", 
                   "Coastal Forests of Eastern Africa", #"East Melanesian Islands",
                   "Eastern Afromontane", "Guinean Forests of West Africa",
                   "Indo-Burma", "Madagascar and the Indian Ocean Islands",
                   "Mesoamerica", "New Caledonia", "Philippines", 
                   "Sundaland", "Tropical Andes",
                   "Tumbes-Choco-Magdalena", "Wallacea", "Western Ghats and Sri Lanka")
hotspots <- hotspots[hotspots@data$NAME %in% trop.hotspots,]
hotspots.no.lim <- hotspots[!hotspots@data$Type %in% "outer limit",]

## Tnc regions
tnc <- readRDS("data/data-raw/ecoregions_2017/Ecoregions2017.rds")
tnc <- tnc[tnc@data$REALM %in% c("Afrotropic", "Neotropic", "Australasia"),]

## New guinea
ext.ng <- extent(tnc)
ext.ng[1] <- 128; ext.ng[2] <- 153
ext.ng[3] <- -14; ext.ng[4] <- 1
ng <- crop(tnc, ext.ng)
exclude.ng <- "Cape York|Arnhem|Kimberly|Victoria|Carpentaria|Banda Sea|Seram Rain Forests|Halmahera Rain Forests|Solomon Islands Rain Forests"
ng <- ng[ !grepl(exclude.ng, ng@data$ECO_NAME, ignore.case = TRUE),]
# plot(ng, col = 1:6)
# plot(ext.ng, add = TRUE)
ng.lim <- aggregate(ng)
ng.lim1 <- rgeos::gBuffer(ng.lim, byid=TRUE, width=0) #correcting possible overlapping polygons
ng.lim1 <- cleangeo::clgeo_Clean(ng.lim1) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
ng.lim1 <- SpatialPolygonsDataFrame(ng.lim1, data.frame(NAME = "New Guinea", 
                                                        Type = "hotspot area"))
#inspecting: ok!
# plot(ng.lim1, col = 1:6)
# plot(hotspots, add = TRUE,
#      border = "blue", lwd = 2)

## Amazon
am.wwf <- readOGR(dsn="data/data-raw/amazon_wwf_limits/data/shape"
                            ,layer="amazonia")
ext.am <- extent(am.wwf)
ext.am[1] <- -80; ext.am[2] <- -44
ext.am[3] <- -19; ext.am[4] <- 11
am <- crop(tnc, ext.am)
exclude.am <- "Cerrado|MaranhÃ£o BabaÃ§u|Sechura|Atacama|Andean|Chaco|Yungas|Puna|ChocÃ³|Llanos|Isthmian|Pantanal|Cordillera|Montane|Xeric|Magdalena-Urab|Catatumbo|Ecuadorian Dry Forests|Alto ParanÃ¡ Atlantic Forests|Bahia Interior Forests|Lake|Chiquitano|MaraÃ±Ã³n"
am <- am[ !grepl(exclude.am, am@data$ECO_NAME, ignore.case = TRUE),]

get_these <- !is.na(over(am, am.wwf)$ID)
am <- am[get_these, ]
am.all <- disaggregate(am)
get_these <- !is.na(over(am.all, am.wwf)$ID)
am.all <- am.all[get_these, ]
am.lim <- aggregate(am.all)
am.lim1 <- rgeos::gBuffer(am.lim, byid=TRUE, width=0) #correcting possible overlapping polygons
am.lim1 <- cleangeo::clgeo_Clean(am.lim1) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
am.lim1 <- SpatialPolygonsDataFrame(am.lim1, data.frame(NAME = "Amazon", 
                                                        Type = "hotspot area"))
#inspecting: ok!
# par(mar=c(1,1,1,1))
# plot(am.lim, col = 1:6)
# plot(hotspots[hotspots@data$NAME %in% "Tropical Andes",], add = TRUE,
#      border = "blue", lwd = 2)
# plot(am.wwf, add = TRUE, border = "red", lwd = 2)

## Central Africa
ext.ca <- extent(tnc)
ext.ca[1] <- 4; ext.ca[2] <- 40
ext.ca[3] <- -20; ext.ca[4] <- 13
ca <- crop(tnc, ext.ca)
# ca <- ca[grepl("Forests|Mangroves", ca@data$WWF_MHTNAM),]
include.ca <- c("Central African mangroves",
                "Central Congolian lowland forests",
                "Congolian coastal forests",
                "Eastern Congolian swamp forests",
                "Northeast Congolian lowland forests",
                "Northern Congolian Forest-Savanna",
                "Northwest Congolian lowland forests",
                "Southern Congolian forest-savanna",
                "Atlantic Equatorial Coastal Forests",
                "Western Congolian forest-savanna",
                "Western Congolian swamp forests")
ca1 <- ca[ca@data$ECO_NAME %in% include.ca,]
ca.lim <- gDifference(ca1, hotspots.no.lim[hotspots.no.lim@data$NAME %in%
                                  c("Eastern Afromontane",
                                    "Guinean Forests of West Africa"),],
            byid = FALSE)
ca.lim1 <- rgeos::gBuffer(ca.lim, byid=TRUE, width=0) #correcting possible overlapping polygons
ca.lim1 <- cleangeo::clgeo_Clean(ca.lim1) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
ca.all <- disaggregate(ca.lim1)
get_these <- gArea(ca.all, byid = TRUE) > 200
ca.all <- ca.all[get_these, ]
ca.lim <- aggregate(ca.all)
ca.lim <- SpatialPolygonsDataFrame(ca.lim, data.frame(NAME = "Central Africa", 
                                                        Type = "hotspot area"))

#Comparing with the shp send by gilles from the Sosef et al. 2019 paper
ca.sosef <- readOGR(dsn="data/data-raw/central_africa",
                            layer="cental_african_forest_sosef2019")
#Inspecting
plot(ca.sosef, border = "white")
plot(ca.lim, add = TRUE, border = "red", lwd = 2)
plot(hotspots.no.lim, add = TRUE, border = "blue", lwd = 2)
plot(ca.sosef, add = TRUE)

## binding all maps
trop.hotspots <- c(trop.hotspots, "New Guinea", "Central Africa", "Amazon")
hotspots <- rbind(hotspots,
                       ng.lim1, ca.lim, am.lim1)

## Saving the final shapefiles
saveRDS(hotspots, "data/hotspots_limits.rds")


###############################
#### ESA 2018 Land use map ####
###############################
# Clipping and masking ESA Land Cover data to the biodiversity hotspots
# Note: Due to the size of the rasters and the processing time 
# only the final results were saved within the project (file 'esa_2018_lc_per_hotspot.csv')

path = "E://ownCloud//W_GIS" # External HD path
# 
# lu <- raster(paste(path,
#                   "//WO_ESA_Land_Cover_map//v.2.1.1//C3S-LC-L4-LCCS-Map-300m-P1Y-2018-v2.1.1.nc",sep=""))
# pj = crs(lu)
# 
# # clipping global map to hotspot's extents
# for(i in 1:length(trop.hotspots)) {
# # for(i in 16:18) {
#   hot.i <- hotspots[hotspots@data$NAME %in% trop.hotspots[i],]
#   proj4string(hot.i) = pj 
#   hot.i1 <- extent(hot.i) - c(1,-1,1,-1)
#   clip <- crop(lu, hot.i1)
#   nome.arquivo <- paste(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//ESA-LCCS-Map-300m-2018-v2.1.1_",
#                         gsub(" ","_",trop.hotspots[i]),
#                         ".tif",sep="")
#   writeRaster(clip, filename= nome.arquivo,
#               format="GTiff",overwrite=TRUE)
# }
# 
# # clipping hotspot's extents to hotspot's limits
# list.rasters <- list.files(paste0(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//"),
#                            full.names = TRUE)
# for(i in 1:length(list.rasters)) {
# # for(i in c(1,10,35)) {
#     raster.i <- raster(list.rasters[i])
#   pj.i <- crs(raster.i)
#   
#   hot.i <- gsub("ESA-LCCS-Map-300m-2018-v2.1.1_|\\.tif", "", 
#                 sapply(strsplit(list.rasters[i],"//"), tail, 1))
#   hot.i <- gsub("_", " ", hot.i)
#   hot.i1 <- hotspots[hotspots@data$NAME %in% hot.i,]
#   proj4string(hot.i1) <- pj.i
#   
#   clip1 <- crop(raster.i, hot.i1)
#   nome.arquivo <- paste(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//ESA-LCCS-Map-300m-2018-v2.1.1_",
#                         gsub(" ","_", hot.i),
#                         "_clipped.tif",sep="")
#   writeRaster(clip1, filename= nome.arquivo,
#               format="GTiff",overwrite=TRUE)
# }
# 
# # masking hotspot's limits
# list.rasters <- list.files(paste0(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//"),
#                            full.names = TRUE)
# list.rasters1 <- list.rasters[grepl("clipped", list.rasters)] 
# for(i in 1:length(list.rasters1)) {
# # for(i in c(1,4,13)) {
#     raster.i <- raster(list.rasters1[i])
#   pj.i <- crs(raster.i)
#   
#   hot.i <- gsub("ESA-LCCS-Map-300m-2018-v2.1.1_|_clipped\\.tif", "", 
#                 sapply(strsplit(list.rasters1[i],"//"), tail, 1))
#   hot.i <- gsub("_", " ", hot.i)
#   hot.i1 <- hotspots[hotspots@data$NAME %in% hot.i,]
#   
#   if ("outer limit" %in% hot.i1@data$Type) # removendo contornos não terrestres 
#     hot.i1 <- hot.i1[!hot.i1@data$Type %in% "outer limit",]
#   
#   proj4string(hot.i1) <- pj.i
#   
#   clip2 <- mask(raster.i, hot.i1)
#   gc()
#   nome.arquivo <- paste(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//ESA-LCCS-Map-300m-2018-v2.1.1_",
#                         gsub(" ","_", hot.i),
#                         "_masked.tif",sep="")
#   writeRaster(clip2, filename= nome.arquivo,
#               format="GTiff",overwrite=TRUE)
# }
# 
# 
# #### GETTING THE VEGETATION COVER PER HOTSPOT ####
# list.rasters <- list.files(paste0(path,"//WO_ESA_Land_Cover_map//v.2.1.1//hotspots//"),
#                            full.names = TRUE)
# list.rasters2 <- list.rasters[grepl("masked", list.rasters)] 
# result <- vector("list", length(list.rasters2))
# nomes <- gsub("ESA-LCCS-Map-300m-2018-v2.1.1_|_masked\\.tif", "", 
#               sapply(strsplit(list.rasters2,"//"), tail, 1))
# names(result) <- nomes
# classes = c(10, 11, 12, 20, 30, 40, 50, 60, 61, 62, 70,71, 72, 80, 81, 82, 90, 
#             100, 110, 120, 121,122, 130, 140, 150, 151,152, 153, 160, 170, 180, 
#             190, 200, 201,202, 210, 220) 
# for(i in 1:length(list.rasters2)) {
#   cat(i, "\n")
#   raster.i <- raster(list.rasters2[i])
#   
#   tmp <- table(factor(raster.i[], levels = classes))
#   cellStats(!is.na(raster.i), sum)
#   sum(area(raster.i, na.rm = TRUE)[])
#   pixs <- sum(tmp)
#   area <- cellStats(area(raster.i, na.rm = TRUE), sum, na.rm = TRUE)
#   tmp1 <- round(100*tmp/pixs, 6)
#   gc()
#   result[[i]] <- c(area_km2 = area, n.pixels = pixs, tmp1)
# }
# 
# # Loading the ESA Land Use legends and replacing the class names
# tmp <- read.csv("data/ESACCI-LC-Legend.csv", as.is=TRUE)
# for(i in 1:length(result)) names(result[[i]])[-c(1,2)] <- 
#                               tmp$LegendTreeCo[match(names(result[[i]])[-c(1,2)], tmp$NB_LAB)] 
# result1 <- do.call(rbind, result)
# result2 <- cbind.data.frame(hotspot.region = gsub("_"," ",row.names(result1)),
#                             area_km2 = result1[,1],
#                             n.pixel = result1[,2],
#                             t(rowsum(t(result1[,-c(1,2)]), colnames(result1[,-c(1,2)]))))
# 
# # Saving the results
# write.csv(result2, "data//esa_2018_lc_per_hotspot.csv")
