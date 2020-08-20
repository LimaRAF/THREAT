###################################################
#### ASSESSING SPECIES EEO AND AOO - CRITERIA B ###
###################################################
rm(list=ls())

### AFAZERES ###
# LER O GUI DO CNCFLORA PARA AS AVALIA??ES: QUAIS COORDENADAS USAR??


#### LOADING PACKAGES ###
devtools::install_github("gdauby/ConR", ref = "master", force = TRUE)
#devtools::install_github("gdauby/ConR@devel", force = TRUE)
library("ConR")
require(data.table)
source("./R/suggestions_for_ConR.r")
source("C://Users//renato//Documents//raflima//R_packages//ConR//R//EOO.sensitivity.R")
source("C://Users//renato//Documents//raflima//R_packages//ConR//R//over.valid.poly.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//subpop.comp.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//proj_crs.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//coord.check.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//subpop.estimation.R")

#### LOADING THE NEOTROPICS MAP ###
neotrop.simp <- readRDS("data/Contour_Neotrop_simplified_tol_005_no_small.rds")

#### LOADING THREAT OCCURRENCE DATA (HERBARIUM + TREECO) ###
oc.data <- readRDS("data/threat_occ_data_final.rds")

#Putting data in the ConR format
MyData <- oc.data[, c("ddlat","ddlon",
                      "tax","higher.tax.rank",
                      "coly","vouchers",
                      "detBy","dety",
                      "tax.check.final","UC",
                      "dist.eoo","tax.conf","source")]
rm(oc.data)

# For testing...
#ids <- grepl("Myrtaceae", MyData$higher.tax.rank)
#ids = grep("Psidium ubatubense",MyData$tax)
#MyData <- MyData[ids,] 

#### OBTAINING THE NAME AND NUMBER OF OCCURRENCES (TOTAL AND SPATIALLY-UNIQUE) PER SPECIES ####
resultado <- readRDS("data/assess_iucn_spp.rds")
names(resultado)[which(names(resultado) %in% c("family.correct1","species.correct2"))] <- c("family","tax")


####################################
#### EXTENT OF OCCURRENCE (EOO) ####
####################################

## Convex Hull method - took 2.5 and 23.5 min in mine machine without saving the
#species maps, using the very light and the ok-resolution Neotropical map, respectively
system.time(
EOO.hull <- EOO.computing(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                                 method = "convex.hull",
                                 method.less.than3 = "not comp",
                                 export_shp = TRUE,
                                 exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                                 write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                                 write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                                 parallel = TRUE, NbeCores = 6) # run on parallel? How many cores?
)

#extracting the EOO from the object
EOO <- do.call(rbind.data.frame, EOO.hull[grepl("EOO", names(EOO.hull))])
sp_names <- unique(MyData$tax[grepl("high", MyData$tax.check.final)])
dimnames(EOO) <- list(sp_names, "EOO")

shps <- EOO.hull[grepl("spatial.polygon", names(EOO.hull))]
for (i in 1:length(shps))
  slot(slot(shps[[i]], "polygons")[[1]], "ID") <- sp_names[!is.na(EOO$EOO)][i]
shps <- do.call(rbind, shps)
shps_df <- SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
shps_df$tax <- as.character(shps_df$tax) 

#inpecting
par(mar=c(1,0,1,0))
sp <- "Abarema limae" # restrictd range
sp <- "Anadenanthera peregrina" # more widespread
plot(shps_df[sp,])
points(MyData[MyData$tax %in% sp, 2:1], pch=21)
points(MyData[grepl("high", MyData$tax.check.final) & MyData$tax %in% sp, 2:1],
       pch = 19)
plot(neotrop.simp, add=TRUE, border = "grey")

#saving species EOO (convex hulls)
shps_df_sf <- sf::st_as_sf(shps_df)
#saveRDS(shps_df, "data/spp.convex.hull.polys.rds")
saveRDS(shps_df_sf, "data/spp.convex.hull.polys_sf.rds")
rm(EOO.hull, EOO, shps, shps_df, shps_df_sf); gc()

## Convex Hull method for each class of confidence level 
# took 8.25 and 73.3 minutes in my machine with the very ligth and with the ok-resolution Neotropical maps, respectively!
system.time(
EOO.hull.sensitivity <- EOO.sensitivity (MyData[, c(1:3,9)],
                                         levels.order = c("low", "medium", "high"),
                                         occ.based = FALSE,
                                         exclude.area = TRUE,
                                         country_map = neotrop.simp,
                                         method.range = "convex.hull",
                                         write_results=FALSE,
                                         parallel = TRUE,
                                         NbeCores = 7,
                                         show_progress = TRUE)
)
saveRDS(EOO.hull.sensitivity, "data/EOO.convex.hull.rds")
rm(EOO.hull.sensitivity); gc()

## Alpha-Hull method
#### GILLES I AM GETTING AN ERROR HERE: ####
#"Evaluation error: TopologyException: Input geom 1 is invalid: Self-intersection at or near point -68.779132918930316 -14.59956504251698 at -68.779132918930316 -14.59956504251698."
# 
# EOO.alpha <- EOO.computing(MyData[grepl("high", MyData$tax.check.final),], 
#                                   method.range = "alpha.hull", alpha = 1, # maybe 2 looks like a better compromise
#                                   exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
#                                   write_results=FALSE, # If TRUE, a csv fiel is created in the working directory
#                                   parallel = TRUE, NbeCores = 7) # run on parallel? How many cores?
# EOO.alpha1 <- EOO.computing(MyData[grepl("high", MyData$tax.check.final),], 
#                             method.range = "alpha.hull", alpha = 2, # maybe 2 looks like a better compromise
#                             exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
#                             write_results=FALSE, # If TRUE, a csv fiel is created in the working directory
#                             parallel = TRUE, NbeCores = 7) # run on parallel? How many cores?
# rm(EOO.alpha,EOO.alpha1,EOO.alpha.sensitivity); gc()


########################
#### SUBPOPULATIONS ####
########################
#This function applies the method called circular buffer method for estimating the number of subpopulation (Rivers et al., 2010).
#From Rivers et al. 2010: The ideal radius of the buffer is debatable; however when dispersal
#characteristics of the species are unknown then a sliding scale, such as the 1/10th maximum inter-point distance, is the preferred choice, as it is species-specific and not sensitive to collection density
radius = subpop.radius(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                       factor.div = 10, quant.max = 0.95)
radius1 = subpop.radius(MyData[, c(1:3)], 
                       factor.div = 10, quant.max = 0.95)

#Getting family-specific radius for missing estimates (only high confidence level)
tmp <- merge(radius, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp <- tmp[order(tmp$tax),]
fam.radius <- aggregate(as.double(tmp$radius), list(tmp$family), median, na.rm = TRUE)
radius.new <- merge(tmp, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new$radius[is.na(radius.new$radius)] <- 
  radius.new$x[is.na(radius.new$radius)]
radius.new <- radius.new[order(radius.new$tax), ]
table(radius$tax == radius$tax)

#Getting family-specific radius for missing estimates (any confidence level)
tmp1 <- merge(radius1, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp1 <- tmp1[order(tmp1$tax),]
fam.radius <- aggregate(as.double(tmp1$radius), list(tmp1$family), median, na.rm = TRUE)
radius.new1 <- merge(tmp1, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new1$radius[is.na(radius.new1$radius)] <- 
  radius.new1$x[is.na(radius.new1$radius)]
radius.new1 <- radius.new1[order(radius.new1$tax), ]
table(radius1$tax == radius1$tax)

#Getting the subpopulations
SUB <- my.subpop.comp(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                      Resol_sub_pop= radius.new[,c("tax", "radius")])
SUB1 <- my.subpop.comp(MyData[, c(1:3)], 
                      Resol_sub_pop= radius.new1[,c("tax", "radius")])

#extracting the numb. of subpopulations from the object (high confidence level)
tmp <- sapply(SUB, function(x) x[grepl("Number of subpopulation", names(x))])
sp_names <- gsub(".Number of subpopulation","",names(tmp))
n.sub.pop <- do.call(rbind.data.frame, tmp)
dimnames(n.sub.pop) <- list(sp_names, "Number.subpops")

#### CHECK HERE: EXTRACTION OF POLYGONS NOT WORKING (IS IT NECESSARY??) ####
tmp <- sapply(SUB, function(x) x[grepl("subpop.poly", names(x))])
sp_names <- gsub(".subpop.poly","",names(tmp))
for (i in 1:length(tmp))
  slot(slot(tmp[[i]], "polygons")[[1]], "ID") <- sp_names[i]
shps <- do.call(rbind, tmp)
shps_df <- SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
shps_df$tax <- as.character(shps_df$tax) 

#extracting the numb. of subpopulations from the object (any confidence level)
tmp <- sapply(SUB1, function(x) x[grepl("Number of subpopulation", names(x))])
sp_names <- gsub(".Number of subpopulation","",names(tmp))
n.sub.pop1 <- do.call(rbind.data.frame, tmp)
dimnames(n.sub.pop1) <- list(sp_names, "Number.subpops")

#Merging and saving
tmp1 <- merge(n.sub.pop1, n.sub.pop, by="row.names", all = TRUE, sort = FALSE)
names(tmp1) <- c("tax","Number.subpops.level.1", "Number.subpops.level.2")
tmp1 <- tmp1[order(tmp1$tax),]
saveRDS(tmp1, "data/number.subpops.rds")
rm(tmp,tmp1,radius,radius1,radius.new,radius.new1,SUB,SUB1,n.sub.pop,n.sub.pop1); gc()

  
#################################
#### AREA OF OCCUPANCY (AOO) ####
#################################
#only high confidence level
system.time(
AOO <- AOO.computing(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3)], 
                          Cell_size_AOO = 2, # Size of the grid cell in kilometers
                    nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                    parallel = TRUE, NbeCores = 7, 
                    show_progress= TRUE, export_shp=FALSE)
)
#any confidence level
system.time(
AOO1 = AOO.computing(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon), c(1:3)], 
                     Cell_size_AOO = 2, # Size of the grid cell in kilometers
                    nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                    parallel = TRUE, NbeCores = 7, 
                    show_progress= TRUE, export_shp=FALSE)
)
#Saving
tmp2 <- merge(AOO1, AOO, by="row.names", all = TRUE, sort = FALSE)
names(tmp2) <- c("tax","AOO.level.1", "AOO.level.2")
tmp2 <- tmp2[order(tmp2$tax),]
saveRDS(tmp2, "data/AOO.rds")
rm(AOO, AOO1, tmp2); gc()


#############################
#### NUMBER OF LOCATIONS ####
#############################

## Loading an pre-processing the Protected Areas Map ###
strict.ucs <- readRDS("data/StrictUCsNeotrop_simplified_clipped.rds")
types <- vapply(sf::st_geometry(strict.ucs), function(x) { class(x)[2]}, "")
polys <- strict.ucs[ grepl("*POLYGON", types), ]
polys.sp <- as(polys, "Spatial")
geom.col <- strict.ucs[ grepl("GEOMETRYCOLLECTION", types), ]
polys.sp1<- list("vector", length(geom.col))
for (i in 1:length(geom.col)){
  tmp <- sf::st_cast(geom.col[i,], "GEOMETRY")
  #tmp1 <- sf::st_cast(tmp)
  tmp1 <- tmp
  types <- vapply(sf::st_geometry(tmp1), function(x) { class(x)[2]}, "")
  tmp2 <- tmp1[ grepl("*POLYGON", types), ]
  tmp3 <- sf::st_cast(tmp2, "MULTIPOLYGON")
  polys.sp1[[i]] <- tmp3
}
polys.sp1 <- do.call(rbind, polys.sp1)
polys.sp1 <- as(polys.sp1, "Spatial")
strict.ucs.spdf <- rbind(polys.sp, polys.sp1)
rm(polys, polys.sp, polys.sp1,geom.col,strict.ucs)

system.time(
locs <- my.locations.comp(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3),], 
                                method = "fixed_grid",
                                nbe_rep = 0, # number of raster with random starting position
                                Cell_size_locations = 2, #grid size in kilometres used for estimating the number of location
                                Rel_cell_size=0.05,
                                protec.areas = strict.ucs.spdf, #a SpatialPolygonsDataFrame, shapefile with protected areas.
                                ID_shape_PA= "NAME", # field name of protec.areas with ID
                                method_protected_area="no_more_than_one", 
                                parallel=TRUE, NbeCores=7)
)

table(names(locs[[3]]) == names(locs[[4]]))
saveRDS(cbind.data.frame(locs[[3]],locs[4]), "data/tmp.rds")

system.time(
locs1 <- my.locations.comp(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon), c(1:3),], 
                            method = "fixed_grid",
                            nbe_rep = 0, # number of raster with random starting position
                            Cell_size_locations = 2, #grid size in kilometres used for estimating the number of location
                            Rel_cell_size=0.05,
                            protec.areas = strict.ucs.spdf, #a SpatialPolygonsDataFrame, shapefile with protected areas.
                            ID_shape_PA= "NAME", # field name of protec.areas with ID
                            method_protected_area="no_more_than_one", 
                            parallel=TRUE, NbeCores=7)
)
table(names(locs1[[3]]) == names(locs1[[4]]))
saveRDS(cbind.data.frame(locs1[[3]],locs1[4]), "data/temporary/tmp1.rds")

##Combining the results
tmp <- cbind.data.frame(locs[[3]],locs[4])
names(tmp) <- c("LocUCs.level2", "LocOutUCs.level2")
tmp1 <- cbind.data.frame(locs1[[3]],locs1[4])
names(tmp1) <- c("LocUCs.level1", "LocOutUCs.level1")
localities <- merge(tmp1, tmp, by= "row.names", all = TRUE, sort = FALSE)
names(localities)[1] <- "tax"
localities$tax <- gsub("_", " ", localities$tax)
localities <- localities[order(localities$tax),]
saveRDS(localities, "data/number_localities.rds")


#### MERGING ALL RESULTS AND SAVING CRITERIA B PARAMETERS ####
EOO <- readRDS("data/EOO.convex.hull.rds")
AOO <- readRDS("data/AOO.rds")
subpops <- readRDS("data/number.subpops.rds")
localities <- readRDS("data/number_localities.rds")
table(EOO$Species == AOO$tax)
subpops <- subpops[order(subpops$tax),]
table(EOO$Species == subpops$tax)
table(EOO$Species == localities$tax)

critB_high <- cbind.data.frame(EOO[, c("Species", "Occs.level.3", "EOO.level.3")],
                               AOO = as.double(AOO[,"AOO.level.2"]),
                               Nbe_subPop = as.double(subpops[,"Number.subpops.level.2"]),
                               Nbe_loc = as.double(localities[,"LocOutUCs.level2"]),
                               Nbe_loc_PA = as.double(localities[,"LocUCs.level2"]),
                               stringsAsFactors = FALSE)
critB_low <- cbind.data.frame(EOO[, c("Species", "Occs.level.1", "EOO.level.1")],
                               AOO = as.double(AOO[,"AOO.level.1"]),
                               Nbe_subPop = as.double(subpops[,"Number.subpops.level.1"]),
                               Nbe_loc = as.double(localities[,"LocOutUCs.level1"]),
                               Nbe_loc_PA = as.double(localities[,"LocUCs.level1"]),
                               stringsAsFactors = FALSE)
names(critB_high)[1:3] <- names(critB_low)[1:3] <- c("tax", "Nbe_occs", "EOO")
critB_high <- critB_high[,c("tax","EOO","AOO","Nbe_subPop","Nbe_loc","Nbe_loc_PA","Nbe_occs")]
critB_low <- critB_low[,c("tax","EOO","AOO","Nbe_subPop","Nbe_loc","Nbe_loc_PA","Nbe_occs")]
saveRDS(critB_high, "data/criteriaB_metrics_high_confidence.rds")
saveRDS(critB_low, "data/criteriaB_metrics_low_confidence.rds")


#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA B ####
#############################
critB_high <- readRDS("data/criteriaB_metrics_high_confidence.rds")
critB_low <- readRDS("data/criteriaB_metrics_low_confidence.rds")

args(criteria.B)



#### GILLES EXTRA CODES ####

library(tidyverse)

# all_source <- 
#   list.files("D:/MonDossierR/ConR_orig/R/", full.names = T)
# for (i in 1:length(all_source)) source(all_source[i])

oc.data_tb <- 
  oc.data %>% 
  as_tibble() %>% 
  dplyr::select(ddlat, ddlon, tax)

## list of species with at least one georeferenced occurrence
oc.data_tb %>% 
  filter(!is.na(ddlat), !is.na(ddlon)) %>% 
  distinct(tax)

## list of species with no georeferenced occurrences
oc.data_tb %>% 
  filter(is.na(ddlat) | is.na(ddlon)) %>% 
  distinct(tax)

## no paralleling 2224.4 secondes
system.time(results_cb <- 
              criterion_B(x = oc.data_tb))
## paralleling
system.time(results_cb <- 
              criterion_B(x = oc.data_tb, parallel = T))
dim(results_cb)

## number and proportion of species in each category
results_cb %>% 
  as_tibble() %>% 
  group_by(category) %>% 
  count() %>% 
  mutate(prop = n/nrow(results_cb)*100)
system.time(
  results_cb <-
    EOO.computing(
      XY = oc.data_tb,
      parallel = T,
      exclude.area = F,
      country_map = neotrop.simp
    )
)

