################################################
#### GETTING DESCRIPTIONS WITHIN SPECIE EOO ####
################################################
rm(list = ls())
require(ConR)
require(raster)
require(sf)
require(lwgeom)
source("C://Users//renato//Documents//raflima//R_packages//ConR//R//EOO.habitat.R")
path <- "E://ownCloud//W_GIS" # Renato's path

#### GETTING SPECIES EOOs ####
EOO.poly <- readRDS(file = "data/spp.convex.hull.polys_sf_uncropped.rds")

#### ESA land-use cover: forest cover as habitat ####
## Getting forest cover classes fro ESA LC map
toto = read.csv(paste(path,"//WO_ESA_Land_Cover_map//v2.0.7//ESACCI-LC-Legend.csv",sep=""), as.is=TRUE)
hab.class <- toto$NB_LAB[grepl("ForestCover", toto$LegendTreeCoSimp)]
#hab.class <- c("50", "60", "61", "62", "70", "71", "72", "80", "81", "82", "90")

## Getting ESA cropped rasters for the AF, aggregating and saving
# fcov <- raster::stack(paste(path,"//WO_ESA_Land_Cover_map//v2.0.7//ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7_AF_50km_buffer_mask.tif",sep=""))
# names(fcov) = paste("ESA",c(1992,1995,2000,2005,2010,2015), sep=".")
# # aggregating the raster (making the resolution coarser)
# tmp <- fcov[[c(1,2,3,6)]]
# tmp <- raster::aggregate(tmp, fact = 3, fun = modal, na.rm = TRUE)
# tmp1 <- tmp[[c(3,4)]]
# writeRaster(tmp1, filename="data/ESA_Land_Cover_map_2000_2015_AF_1km.tif",
#             format="GTiff",overwrite=TRUE)
# rm(fcov, tmp, tmp1)

## Loading aggregated rasters and extracting the habitat loss/quality
ano1 = 2000
hab.map <- raster::stack(paste0("data/ESA_Land_Cover_map_",ano1,"_2015_AF_1km.tif"))
names(hab.map) <- paste("ESA",c(ano1,2015), sep=".")
anos <- 2015 - ano1 

toto.1 <- EOO.habitat(EOO.poly[1:1000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
toto.2 <- EOO.habitat(EOO.poly[1001:2000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
toto.3 <- EOO.habitat(EOO.poly[2001:3000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
toto.4 <- EOO.habitat(EOO.poly[3001:4000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
toto.5 <- EOO.habitat(EOO.poly[4001:4862,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
toto <- rbind.data.frame(toto.1, toto.2, toto.3, toto.4, toto.5,
                          stringsAsFactors = FALSE)
toto0<- merge(data.frame(tax = EOO.poly$tax, stringsAsFactors = FALSE), 
              toto, by = "tax", all = TRUE, sort = FALSE)
toto0<-toto0[order(toto0$tax),] 
saveRDS(toto0, 
        paste0("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo Extincao na MA/data analysis/data/EOO_hab_loss_",
               ano1,"_2015.rds"))

#Any rate (habitat loss > 0.1%)
table(toto0$rate.loss>0.1)
table(toto0$rate.loss>toto0$recover)

#### PROTECTED AREAS ####
strict.ucs <- readRDS("data/StrictUCsNeotrop_simplified_clipped.rds")
types <- vapply(sf::st_geometry(strict.ucs), function(x) { class(x)[2]}, "")
polys <- strict.ucs[ grepl("*POLYGON", types), ]
hab.map <- polys
rm(polys, strict.ucs)

toto1 <- EOO.habitat(EOO.poly, hab.map, hab.class = NULL, years = NULL, ID_shape = "WDPAID",
                    parallel = TRUE, NbeCores = 6)
saveRDS(toto1, "data/EOO_StrictUCs.rds")


#### HUMAN INFLUENCES #### 
# data between 1995 - 2004 (0 low and 64 is high human influence)
hab.map <- raster::raster(paste(path,"//Am_Lat_Global_Human_Influence_Index_v2//hii_neotrop.tif",sep=""))

#t0 <- Sys.time()
toto2.1 <- EOO.habitat(EOO.poly[1:1000,], hab.map, hab.class = NULL, years = NULL,
                       parallel = TRUE, NbeCores = 6)
toto2.2 <- EOO.habitat(EOO.poly[1001:2000,], hab.map, hab.class = NULL, years = NULL,
                       parallel = TRUE, NbeCores = 6)
toto2.3 <- EOO.habitat(EOO.poly[2001:3000,], hab.map, hab.class = NULL, years = NULL,
                       parallel = TRUE, NbeCores = 6)
toto2.4 <- EOO.habitat(EOO.poly[3001:4000,], hab.map, hab.class = NULL, years = NULL,
                       parallel = TRUE, NbeCores = 6)
toto2.5 <- EOO.habitat(EOO.poly[4001:length(EOO.poly$geometry),], hab.map, hab.class = NULL, years = NULL,
                       parallel = TRUE, NbeCores = 6)
toto2 <- rbind.data.frame(toto2.1, toto2.2, toto2.3, toto2.4, toto2.5,
                          stringsAsFactors = FALSE)
#Sys.time() - t0
saveRDS(toto2, "data/EOO_HII.rds")


#### EOO inside the Atlantic Forest limits ####
af <- rgdal::readOGR(dsn=paste(path,"//AF_limites_milton",sep=""),layer="merge_limites_MA11428_TNC_ARG_BR_PAR")
af <- rgeos::gBuffer(af, byid=TRUE, width=0) #correcting possible overlapping polygons
hab.map <- cleangeo::clgeo_Clean(af) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
hab.map <- rgeos::gBuffer(hab.map, byid=TRUE, width=0.05) #adding some margin (~5km buffer) for possible projection and spatial resolution differences
#Again for safety...
hab.map <- rgeos::gBuffer(hab.map, byid=TRUE, width=0) #correcting possible overlapping polygons
hab.map <- cleangeo::clgeo_Clean(hab.map) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
rm(af)

toto3 <- EOO.habitat(EOO.poly, hab.map, hab.class = NULL, years = NULL, ID_shape = "OBJECTID",
                     parallel = TRUE, NbeCores = 6)
toto3$prop.EOO[!is.na(toto3$prop.EOO) & toto3$prop.EOO>100] <- 100
saveRDS(toto3, "data/EOO_AtlanticForest.rds")

#### Agricultural Land Use in Brazil ####
#Note: 1940 until 2014
#The distribution and intensity of overall agricultural land use in Brazil from 1940 to 1995:
#Reference: Dias et al. 2016 Global change biology 22(8): 2887-2903.
#new explicitly spatialized database of agriculture areas in Brazil which includes 
#croplands (total between 1940 and 2012), pasturelands (natural and planted between 1940 and 2012). 
#We reconstructed the agricultural historical patterns by combining agricultural census data and remote sensing data 
#for the whole of Brazil at 30" spatial resolution (approximately 1 km x 1 km).

lu <- raster::stack(paste(path,"//BR_LandUse_1940_2014//LUCULT19402014.nc",sep=""), varname = "landuse")
years = names(lu)[c(1,6,11,16,21,26,31,36,41,46,51,53,56,61,66,71,75)]
lu <- lu[[years]]
pj <- crs(lu)
grid = spTransform(hex.grid1, pj)
