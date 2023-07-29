########################################################
#### A DATABASE AND GAZETTEER OF CONSERVATION UNITS ####
########################################################
rm(list = ls())
gc()

## Loading packages
require(wdpar)
require(lwgeom)
require(sf)

#### DOWNLOADING, CLIPPING AND SIMPLIFYING WDPA FOR THE NEOTROPICS ####
# Downloading
raw_pa_data <- wdpa_fetch("global", 
                          download_dir = file.path(here::here(), 
                                                   "data", "data-raw"))
# Loading
raw_pa_data <- wdpa_read(file.path(here::here(), 
                                   "data", "data-raw", "WDPA_Jul2023_Public.zip"))

# Removing non-neotropical PAs
neotrop.countries <- data.frame(PARENT_ISO3 = raw_pa_data$PARENT_ISO3)
neotrop.countries$PARENT <- 
  countrycode::countrycode(neotrop.countries$PARENT_ISO3, "iso3c", "country.name")
tmp <- list.files("E://ownCloud//W_GIS//Am_Lat_ADM_GADM_v3.6", pattern = ".rds")
tmp <- as.character(unlist(sapply(tmp, function(x) strsplit(x, "_")[[1]][2])))
neotrop.countries$neotrop.check <- neotrop.countries$PARENT_ISO3 %in% tmp
neotrop.countries$neotrop.check[neotrop.countries$neotrop.check %in% FALSE & 
                                  grepl("US-FL", raw_pa_data$SUB_LOC)] <- TRUE
raw_pa_data <- raw_pa_data[neotrop.countries$neotrop.check,]

#Removing Marine-only PAs
tmp <- raw_pa_data[, "MARINE", drop = TRUE]
tmp <- !tmp %in% "2" 
raw_pa_data <- raw_pa_data[tmp,]

#Removing PAs not currently implemented
tmp <- raw_pa_data[, "STATUS", drop = TRUE]
tmp <- tmp %in% c("Designated", "Established", "Inscribed") 
raw_pa_data <- raw_pa_data[tmp,]

#Filtering only PA of strict protecion (IUCN categories I, II and III)
tmp <- raw_pa_data[, c("DESIG", "DESIG_ENG", "IUCN_CAT"), drop = TRUE]
tmp$string <- paste(tmp$DESIG,tmp$DESIG_ENG, sep="_")
#write.csv(tmp[!duplicated(tmp1),], "PA_designation.csv", row.names = FALSE)
tmp1 <- read.csv("data/PA_designation.csv", as.is = TRUE, na.string = c(NA,""," "))
tmp1$string <- paste(tmp1$DESIG,tmp1$DESIG_ENG, sep="_")
tmp2 <- dplyr::left_join(tmp, tmp1, by= "string")
tmp3 <- tmp2$IUCN_CAT_THREAT %in% c("Ia", "Ib", "II", "III")
raw_pa_data <- raw_pa_data[tmp3,]

## Filtering uncecessary columns
raw_pa_data1 <- raw_pa_data[, c(1,4,7,8,13,17:21,26,27, 29)]
# saveRDS(raw_pa_data1, "data/data-raw/StrictUCsNeotrop.rds")

# Cleaning the database -------------------------------------------------------
raw_pa_data1 <- readRDS("data/data-raw/StrictUCsNeotrop.rds")

# Repair invalid geometry
raw_pa_data1 <- sf::st_make_valid(raw_pa_data1)

# getting the points
GEOMETRY_TYPE <- attributes(st_geometry(raw_pa_data1$geometry))$classes
table(GEOMETRY_TYPE == "MULTIPOINT" & raw_pa_data1$REP_AREA<0.01) #No empty areas for points!

# geometries are wrapped to the dateline
raw_pa_data1 <- st_wrap_dateline(raw_pa_data1, 
                                 options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
                 quiet = TRUE) 

# Reproject data to coordinate system specified in argument to crs
raw_pa_data1 <- st_transform(raw_pa_data1, 5641)  

# Fix any invalid geometries that have manifested
raw_pa_data1 <- sf::st_make_valid(raw_pa_data1)

# # Buffer areas represented as point localities to circular areas
# raio <- raw_pa_data1$REP_AREA[GEOMETRY_TYPE == "MULTIPOINT"] # area in km2
# raio <- sqrt(raio/pi)
# buffs <- st_buffer(raw_pa_data1[GEOMETRY_TYPE == "MULTIPOINT",], dist = raio*1000 )
# plot(raw_pa_data1[GEOMETRY_TYPE == "MULTIPOINT",1], pch=19)
# plot(buffs[,1])
# raw_pa_data1[GEOMETRY_TYPE == "MULTIPOINT",] <- buffs

# Snap the geometries to a grid to fix any remaining geometry issues 
raw_pa_data2 <- st_snap_to_grid(raw_pa_data1, 1)

# Simplify the protected area geometries to reduce computational burden
raw_pa_data3 <- sf::st_simplify(raw_pa_data2, dTolerance=10, preserveTopology = TRUE)

# Fix any invalid geometries that have manifested
pa_data_clean <- sf::st_make_valid(raw_pa_data3)

# Removing slivers (geometies wih less than 0.1 square meters)
tmp <- pa_data_clean$REP_AREA>1e-7
pa_data_clean <- pa_data_clean[tmp,]
#saveRDS(pa_data_clean, "data/StrictUCsNeotrop_simplified_proj.rds")

#Reprojecting the map
pa_data_clean1 <- st_transform(pa_data_clean, 4326)

#Importing Neotropical bondaries
bond_data <- readRDS("data/data-raw/gadm36_Neotrop_0_sp_simplified.rds")
bond_data <- st_as_sf(bond_data)
pa_data_clean2 <- st_transform(pa_data_clean1, st_crs(bond_data))
plot(st_geometry(pa_data_clean2))
plot(st_geometry(bond_data), add = TRUE)
#saveRDS(pa_data_clean2, "data/StrictUCsNeotrop_simplified.rds")

## Clip UCs with coast line
pa_data_clean3 <- st_intersection(pa_data_clean2, bond_data)
saveRDS(pa_data_clean3, "data/StrictUCsNeotrop_simplified_clipped.rds")
