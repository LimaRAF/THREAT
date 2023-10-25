################################################
#### GETTING DESCRIPTIONS WITHIN SPECIE EOO ####
################################################
rm(list = ls())
gc()

#### LOADING THE REQUIRED PACKAGES AND FUNCTIONS ####
require(ConR)
require(raster)
require(sf)
require(lwgeom)
require(rgdal)
require(rgeos)
require(cleangeo)
source("R//99_functions.R")

#### GETTING SPECIES EOOs ####
EOO.poly <- readRDS(file = "data/spp.convex.hull.polys_sf_uncropped.rds")
EOO.poly <- sf::st_make_valid(EOO.poly)

#### ESA land-use cover: forest cover as habitat ####
## Getting forest cover classes fro ESA LC map
toto <- read.csv("data/ESACCI-LC-Legend.csv", as.is=TRUE)
hab.class <- toto$NB_LAB[grepl("ForestCover", toto$LegendTreeCoSimp)]

## Getting ESA cropped rasters for the AF, aggregating and saving
## DONT RUN: ALREADY CALLING THE SAVED OBJECT BELOW
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
hab.map <- raster::stack(paste0("data/ESA_Land_Cover_map_", ano1, "_2015_AF_1km.tif"))
names(hab.map) <- paste("ESA",c(ano1,2015), sep=".")
anos <- 2015 - ano1 

sf::sf_use_s2(FALSE)
toto.1 <- EOO.habitat(EOO.poly[1:1000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
saveRDS(toto.1, "data/tmp1.eoo.rds")
toto.2 <- EOO.habitat(EOO.poly[1001:2000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
saveRDS(toto.2, "data/tmp2.eoo.rds")
toto.3 <- EOO.habitat(EOO.poly[2001:3000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
saveRDS(toto.3, "data/tmp3.eoo.rds")
toto.4 <- EOO.habitat(EOO.poly[3001:4000,], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
saveRDS(toto.4, "data/tmp4.eoo.rds")
toto.5 <- EOO.habitat(EOO.poly[4001:dim(EOO.poly)[1],], hab.map, hab.class = hab.class, years = anos,
                      parallel = TRUE, NbeCores = 5)
saveRDS(toto.5, "data/tmp5.eoo.rds")
sf::sf_use_s2(TRUE)
toto <- rbind.data.frame(toto.1, toto.2, toto.3, toto.4, toto.5,
                          stringsAsFactors = FALSE)

toto0 <- merge(data.frame(tax = EOO.poly$tax, stringsAsFactors = FALSE), 
               toto, by = "tax", all = TRUE, sort = FALSE)
toto0 <- toto0[order(toto0$tax),] 
saveRDS(toto0, 
        paste0("data/EOO_hab_loss_",
               ano1,"_2015_uncropped.rds"))

#Any rate (habitat loss > 0.1%)
table(toto0$rate.loss>0.1)
table(toto0$rate.loss>toto0$recover)

#### PROTECTED AREAS ####
strict.ucs <- readRDS("data/StrictUCsNeotrop_simplified_clipped.rds")
types <- vapply(sf::st_geometry(strict.ucs), function(x) { class(x)[2]}, "")
polys <- strict.ucs[ grepl("*POLYGON", types), ]
hab.map <- sf::st_make_valid(polys)
rm(polys, strict.ucs)

toto1 <- EOO.habitat(EOO.poly, hab.map, hab.class = NULL, years = NULL, #ID_shape = "WDPAID",
                     parallel = TRUE, NbeCores = 6)
# saveRDS(toto1, "data/EOO_StrictUCs.rds")
saveRDS(toto1, "data/EOO_StrictUCs_cropped.rds")


#### HUMAN INFLUENCES #### NO LONGER BEING CALCULATED
# # data between 1995 - 2004 (0 low and 64 is high human influence)
# hab.map <- raster::raster(paste(path,"//Am_Lat_Global_Human_Influence_Index_v2//hii_neotrop.tif",sep=""))
# 
# #t0 <- Sys.time()
# toto2.1 <- EOO.habitat(EOO.poly[1:1000,], hab.map, hab.class = NULL, years = NULL,
#                        parallel = TRUE, NbeCores = 6)
# saveRDS(toto2.1, "data/tmp2.1.eoo.rds")
# toto2.2 <- EOO.habitat(EOO.poly[1001:2000,], hab.map, hab.class = NULL, years = NULL,
#                        parallel = TRUE, NbeCores = 6)
# saveRDS(toto2.2, "data/tmp2.2.eoo.rds")
# toto2.3 <- EOO.habitat(EOO.poly[2001:3000,], hab.map, hab.class = NULL, years = NULL,
#                        parallel = TRUE, NbeCores = 6)
# saveRDS(toto2.3, "data/tmp2.3.eoo.rds")
# # Problematic species: Myriocarpa cordifolia, Ocotea revolutifolia
# ids <- 3001:4000
# prob.ids <- c(which(EOO.poly$tax %in% c("Myriocarpa cordifolia", "Ocotea revolutifolia")))
# ids <- ids[!ids %in% prob.ids]
# toto2.4 <- EOO.habitat(EOO.poly[ids,], hab.map, hab.class = NULL, years = NULL,
#                        parallel = TRUE, NbeCores = 6)
# saveRDS(toto2.4, "data/tmp2.4.eoo.rds")
# # toto2.4.1 <- EOO.habitat(EOO.poly[3001:3305,], hab.map, hab.class = NULL, years = NULL,
# #                        parallel = TRUE, NbeCores = 6)
# # toto2.4.2 <- EOO.habitat(EOO.poly[3307:4000,], hab.map, hab.class = NULL, years = NULL,
# #                        parallel = TRUE, NbeCores = 6)
# toto2.5 <- EOO.habitat(EOO.poly[4001:length(EOO.poly$geometry),], hab.map, hab.class = NULL, years = NULL,
#                        parallel = TRUE, NbeCores = 6)
# saveRDS(toto2.5, "data/tmp2.5.eoo.rds")
# # toto2 <- rbind.data.frame(toto2.1, toto2.2, toto2.3, toto2.4, toto2.5,
# #                           stringsAsFactors = FALSE)
# tmp <- toto2.1[1:length(prob.ids),]
# tmp$tax <- EOO.poly$tax[prob.ids]
# tmp[,2:dim(tmp)[2]] <- NA
# toto2 <- rbind.data.frame(toto2.1, toto2.2, toto2.3, toto2.4, tmp, toto2.5,
#                           stringsAsFactors = FALSE)
# toto2 <- toto2[order(toto2$tax),] 
# 
# #Sys.time() - t0
# # saveRDS(toto2, "data/EOO_HII.rds")
# saveRDS(toto2, "data/EOO_HII_cropped.rds")

#### EOO inside the Atlantic Forest limits ####
af <- rgdal::readOGR(dsn= "data/AF_limits/",
                     layer="merge_limites_MA11428_TNC_ARG_BR_PAR")
af <- rgeos::gBuffer(af, byid=TRUE, width=0) #correcting possible overlapping polygons
hab.map <- cleangeo::clgeo_Clean(af) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
hab.map <- rgeos::gBuffer(hab.map, byid=TRUE, width=0.05) #adding some margin (~5km buffer) for possible projection and spatial resolution differences
#Again for safety...
hab.map <- rgeos::gBuffer(hab.map, byid=TRUE, width=0) #correcting possible overlapping polygons
hab.map <- cleangeo::clgeo_Clean(hab.map) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
rm(af)

sf::sf_use_s2(FALSE)
toto3 <- EOO.habitat(EOO.poly, hab.map, hab.class = NULL, years = NULL, #ID_shape = "OBJECTID",
                     parallel = TRUE, NbeCores = 6)
sf::sf_use_s2(TRUE)
toto3$prop.EOO[!is.na(toto3$prop.EOO) & toto3$prop.EOO>100] <- 100
# saveRDS(toto3, "data/EOO_AtlanticForest.rds")
saveRDS(toto3, "data/EOO_AtlanticForest_cropped.rds")
