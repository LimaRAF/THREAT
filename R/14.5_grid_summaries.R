#####################################################
#####################################################
#### OBTAINING THE RED LIST INDEX PER GRID CELL #####
#####################################################
#####################################################
rm(list=ls())

## Loading packages
require(data.table)
require(sp)

# -----------------------------------------------------------------------------
####################################
#### ADAPTATIVE RESOLUTION GRID #### 
####################################


## LOADING THE DATA  ##
## Threat occurrence data ##
occs <- readRDS("data/threat_occ_data_final.rds")
occs <- occs[,c("tax","ddlat","ddlon")]
#removing possible missing coordinates
occs$ddlat = as.numeric(occs$ddlat)
occs$ddlon = as.numeric(occs$ddlon)
occs <- occs[!is.na(ddlat) | !is.na(ddlon),]

## Threat assessments
all.crit <- readRDS("data/all.criteria.rds")
# Getting NT for each criteria
subcriteria <- c("category_A", "category_B","category_C","category_D")
subcriteria1 <- list(c("A2"), c("B1","B2"),c("C1", "C2"), "D")
for(i in 1:length(subcriteria)) {
  all.crit[,subcriteria[i]] <- 
    ConR::near.threatened(all.crit[,subcriteria[i]],
                          all.crit$EOO, all.crit$AOO, all.crit$declineB,
                          all.crit$reduction_A12, all.crit$pop.size, all.crit$pop.size.low,
                          all.crit$nbe_loc_total, all.crit$sever.frag, ext.fluct = NULL,
                          subpop = all.crit$Nbe_subPop, subcriteria = subcriteria1[[i]])
}
# Filtering and merging the two datasets
colunas <- c("species", "endemic", "cat.reg.clean",
             "category_A", "category_B","category_C","category_D")
occs <- merge(occs, all.crit[,colunas], by.x = "tax", by.y = "species")

## Threat category as species names
occs$category <- occs$cat.reg.clean
# occs$category[occs$category %in% "NA"] <- "LC"
col.order <- c("tax", "ddlat", "ddlon", "endemic","category",
               "category_A", "category_B", "category_C", "category_D")
occs <- occs[ , .SD, .SDcols = col.order]
names(occs) <- c("Name", "Lat", "Long", "Endemism", "Category",
                 "Category_A", "Category_B", "Category_C", "Category_D")
occs$Ordem <- 1:dim(occs)[1]


##Transforming the occurrences into spatial points
occs.sp = SpatialPointsDataFrame(coords = cbind(occs$Long, occs$Lat), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs"),
                                 data = occs[,c("Name", "Endemism", "Category", 
                                                "Category_A", "Category_B", 
                                                "Category_C", "Category_D",
                                                "Ordem")])

## Getting only the occurrences within the AF buffer
path <- "E://ownCloud//W_GIS" # Renato's path
af <- rgdal::readOGR(dsn=paste(path,"//AF_limites_milton",sep=""),layer="merge_limites_MA11428_TNC_ARG_BR_PAR")
af <- rgeos::gBuffer(af, byid=TRUE, width=0) #correcting possible overlapping polygons
af <- cleangeo::clgeo_Clean(af) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
af <- raster::aggregate(af, dissolve=TRUE)
#creating the Atlantic Forest 0.25 degree (~25 km) buffer
af.buf <- rgeos::gBuffer(af, byid=TRUE, width = 0.25)
af.buf <- raster::disaggregate(af.buf)
af.buf <- cleangeo::clgeo_Clean(af.buf)
#Removing buffer inner-holes
toto <- SpatialPolygons(list(Polygons(list(af.buf@polygons[[1]]@Polygons[[1]]),ID=1)))
crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")  # geographical, datum WGS84
proj4string(toto) <- crs.geo
af.buf <- rbind(af.buf[2:9], toto)
af.buf <- raster::aggregate(af.buf)
# #cropping
# occs.sp1 <- raster::crop(occs.sp, af.buf)
# #calculating the distance between points and all grid cells
# # plot(grid1)
#   
# ## Saving data for the Infomap Bioregions application
# occs1 <- as.data.frame(occs.sp1)
# names(occs1)[names(occs1) %in% c("coords.x1", "coords.x2")] <- c("Long", "Lat")
# # write.csv(occs1, "data/threat_occ_data_all.csv", sep = ",")
# 
# ## Crossing the saved shapefiles/grids from Infomap Bioregions application with the AF map
# nomes <- c("100_max500","200_max500","200_max500_cellmax2","200_max1000","200_max1000_cellmax2")
# for(i in 1:length(nomes)) {
#   nome.i <- paste0("./data/adaptative_resolution_AF_min", nomes[i])
#   shp <- rgdal::readOGR(dsn=nome.i, layer = "threat_occ_data_all")
#   af1 <- spTransform(af,  raster::crs(shp))
#   shp1 <- shp[(over(shp, af1) %in% 1),]
#   shp1@data$plotOrder <- shp1@plotOrder
#   shp1@data$ID <- NA
#   for(i in 1:length(shp1)) { shp1@data$ID[i] <- as.character(slot(slot(shp1, "polygons")[[i]],"ID")) } 
#   shp1@data <- shp1@data[, c("plotOrder","ID")]
#   saveRDS(shp1, paste0(nome.i,"/clean_grid.rds"))
# }


## Loading the shapefile from Infomap Bioregions application
grid0 <- readRDS("./data/adaptative_resolution_AF_min200_max500_cellmax2/clean_grid.rds")

## removing overllaping parts of the grid cells
ids.polys <- grid0$ID
grid_sf <- sf::st_as_sf(grid0)
grid_inter <- sf::st_intersection(grid_sf)
# x_inter <- x_inter[x_inter$n.overlaps == 1,]
grid_inter <- grid_inter[rownames(grid_inter) %in% ids.polys,]
grid <- as(grid_inter, "Spatial")

#Checking the result
# par(mfrow = c(1,2), mar = c(0.5,0.5,0.5,0.5))
# plot(grid0); plot(af.buf, add=TRUE)
# plot(grid); plot(af.buf, add=TRUE)

## Crossing the occurrences with the grid
occs.over <- over(occs.sp, grid) 
occs1 <- cbind.data.frame(occs, occs.over[,c("plotOrder", "ID")])

## Re-assing grid ID for occurrences outside the grid
# plot(grid)
# plot(occs.sp1[is.na(occs1$ID),], add= TRUE, col =2, pch = 19, cex = 0.5)
tmp <- occs.sp[is.na(occs1$ID),]
#cropping for the AF 25 km buffer
tmp <- raster::crop(tmp, af.buf)

#calculating the distance between points and all grid cells
# plot(grid1)
# plot(tmp, add= TRUE)
tmp.proj <- spTransform(tmp,  CRS("+init=epsg:5641"))
grid.proj <- spTransform(grid,  CRS("+init=epsg:5641"))
tmp1 <- rgeos::gDistance(tmp.proj, grid.proj, byid=TRUE)
result = cbind.data.frame(tax = as.character(tmp$Name), 
                          cell.ID = NA_character_, stringsAsFactors = FALSE)
for(i in 1:dim(tmp1)[2]) {
  IDi = which.min(tmp1[,i])
  dist1 = tmp1[which.min(tmp1[,i]),i]
  # if(dist1 > 0 & dist1 < 10000) result$cell.ID[i] = as.character(paste("ID", IDi, sep = ""))
  if(dist1 > 0 & dist1 < 10000) result$cell.ID[i] = as.character(IDi)
}
occs1$ID[match(tmp$Ordem, occs1$Ordem)] <- result$cell.ID

## Removing points falling outside the grid limits
occs1 <- occs1[!is.na(occs1$ID), ]
# plot(grid)
# points(occs1$Long, occs1$Lat, pch=19, col=2, cex = 0.3)
tail(sort(table(occs1$ID, useNA = "always")))

## Removing grid cells without any occurrence? not for now
# cells.ids <- sort(unique(occs1$ID))
# grid1 <- grid[grid$ID %in% cells.ids,]

## Clipping the grid to get the cell area
af.buf <- spTransform(af.buf,  CRS("+init=epsg:5641"))
grid.proj1 <- raster::crop(grid.proj, af.buf)

#### Summary of the sampling effort per grid cell ####
## Getting the new (clipped) polygon centroids
source("R//polygon_centroid.R")
grid.proj2 <- centroid(grid.proj1)
grid.center <- polygonCenter(grid.proj)
grid.center1 <- polygonCenter(grid.proj1)

## Getting the summaries per grid cell
grid.result = cbind.data.frame(order = 1:length(grid.proj), 
                               ID = as.character(sapply(slot(grid.proj, "polygons"), function(x) slot(x, "ID"))),
                               plotOrder = slot(grid.proj, "plotOrder"),
                               area_ha = as.double(rgeos::gArea(grid.proj, byid = TRUE)/10000),
                               # lat = sp::coordinates(grid.proj)[,2],
                               # long = sp::coordinates(grid.proj)[,1],
                               lat = sp::coordinates(rgeos::gCentroid(grid.proj, byid = TRUE))[,2],
                               long = sp::coordinates(rgeos::gCentroid(grid.proj, byid = TRUE))[,1],
                               lat1 = grid.center[,2],
                               long1 = grid.center[,1]
)
grid.result$ID = as.character(grid.result$ID)

grid.result.clip = cbind.data.frame(ID = as.character(sapply(slot(grid.proj1, "polygons"), function(x) slot(x, "ID"))),
                               area_ha.clip = as.double(rgeos::gArea(grid.proj1, byid = TRUE)/10000),
                               # lat.clip = sp::coordinates(rgeos::gCentroid(grid.proj1, byid = TRUE))[,2],
                               # long.clip = sp::coordinates(rgeos::gCentroid(grid.proj1, byid = TRUE))[,1],
                               lat.clip = grid.center1[,2],
                               long.clip = grid.center1[,1]
)
grid.result.clip$ID = as.character(grid.result.clip$ID)
grid.result <- dplyr::left_join(grid.result, grid.result.clip)

##Inspecting the results
# plot(grid.proj, xlim=c(3874270, 4100313), ylim = c(7371300, 8070668), border = "grey")
# plot(grid.proj1, add = TRUE)
# points(grid.result[,c("long","lat")], pch = 19, cex=0.6, col = "red")
# points(grid.result[,c("long1","lat1")], pch = 19, cex=0.6, col = "blue")
# points(grid.result[,c("long2","lat2")], pch = 19, cex=0.6, col = "darkgreen")
# text(grid.result[,c("long")],grid.result[,c("lat")], grid.result[,c("plotOrder")], cex=0.6, col = "red")
# text(grid.result[,c("long1")],grid.result[,c("lat1")], grid.result[,c("plotOrder")], cex=0.6, col = "blue")
# text(grid.result[,c("long.clip")],grid.result[,c("lat.clip")], grid.result[,c("plotOrder")], cex=0.6, col = "darkgreen")
# text(grid.result[,c("long.clip1")],grid.result[,c("lat.clip1")], grid.result[,c("plotOrder")], cex=0.6, col = "darkorange")


### Number of occurrences, number of species and sampling coverage per site ###
tmp = iNEXT::DataInfo(as.data.frame.matrix(table(occs1$Name, occs1$ID))) # getting metrics
tmp$SC[tmp$n==1] = 0.0001 # correcting values of SampCover for sites with N<10 and S ~~ N
tmp$SC[tmp$n==2&tmp$S.obs==2] = 0.0005; tmp$SC[tmp$n==2&tmp$S.obs==1] = 0.001
tmp$SC[tmp$n==3&tmp$S.obs==3] = 0.001; tmp$SC[tmp$n==3&tmp$S.obs==2] = 0.005; tmp$SC[tmp$n==3&tmp$S.obs==1] = 0.01
tmp$SC[tmp$n==4&tmp$S.obs==4] = 0.01; tmp$SC[tmp$n==4&tmp$S.obs==3] = 0.025; tmp$SC[tmp$n==4&tmp$S.obs==2] = 0.05; tmp$SC[tmp$n==4&tmp$S.obs==1] = 0.1
tmp$SC[tmp$n==5&tmp$S.obs==5] = 0.1; tmp$SC[tmp$n==5&tmp$S.obs==4] = 0.125; tmp$SC[tmp$n==5&tmp$S.obs==3] = 0.15; tmp$SC[tmp$n==5&tmp$S.obs==2] = 0.175; tmp$SC[tmp$n==5&tmp$S.obs==1] = 0.2
plot(tmp$SC ~ log(tmp$n))
names(tmp)[1:4] = c("ID","N.total","S.total","SampCover")
grid.result <- merge(grid.result, tmp[,c("ID","N.total","S.total","SampCover")], 
                     by.x= "ID", by.y= "ID", all.x=TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

### Number of occurrences, number of species and sampling coverage per site ###
endemics <- occs1$Endemism %in% "endemic"
tmp = iNEXT::DataInfo(as.data.frame.matrix(table(occs1$Name[endemics], occs1$ID[endemics]))) # getting metrics
tmp$SC[tmp$n==1] = 0.0001 # correcting values of SampCover for sites with N<10 and S ~~ N
tmp$SC[tmp$n==2&tmp$S.obs==2] = 0.0005; tmp$SC[tmp$n==2&tmp$S.obs==1] = 0.001
tmp$SC[tmp$n==3&tmp$S.obs==3] = 0.001; tmp$SC[tmp$n==3&tmp$S.obs==2] = 0.005; tmp$SC[tmp$n==3&tmp$S.obs==1] = 0.01
tmp$SC[tmp$n==4&tmp$S.obs==4] = 0.01; tmp$SC[tmp$n==4&tmp$S.obs==3] = 0.025; tmp$SC[tmp$n==4&tmp$S.obs==2] = 0.05; tmp$SC[tmp$n==4&tmp$S.obs==1] = 0.1
tmp$SC[tmp$n==5&tmp$S.obs==5] = 0.1; tmp$SC[tmp$n==5&tmp$S.obs==4] = 0.125; tmp$SC[tmp$n==5&tmp$S.obs==3] = 0.15; tmp$SC[tmp$n==5&tmp$S.obs==2] = 0.175; tmp$SC[tmp$n==5&tmp$S.obs==1] = 0.2
plot(tmp$SC ~ log(tmp$n))
names(tmp)[1:4] = c("ID","N.total.end","S.total.end","SampCover.end")
grid.result <- merge(grid.result, tmp[,c("ID","N.total.end","S.total.end","SampCover.end")], 
                     by.x= "ID", by.y= "ID", all.x=TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

### Calculating the RLI per grid cell
grid.data <- occs1
grid.data$Category.new <- grid.data$Category
grid.data$Category.new[grid.data$Category %in% "NA"] <- "LC"
resultado <- vector("list", dim(grid.result)[1])
names(resultado) <- grid.result$ID
resultado1 <- resultado
boot = 4999
boot.log = FALSE
for (i in 1: length(resultado)) {
  nome.i <- names(resultado)[i]

  ## ALL POPULATIONS
  grid.data.i <- grid.data[grid.data$ID %in% nome.i,]  
  if (dim(grid.data.i)[1] > 0) { # All records
    #RLI
    rli.i <- red::rli(grid.data.i$Category, boot = boot.log, runs = boot)
    rli.i.A <- red::rli(grid.data.i$Category_A, boot = boot.log, runs = boot)
    rli.i.B <- red::rli(grid.data.i$Category_B, boot = boot.log, runs = boot)
    rli.i.C <- red::rli(grid.data.i$Category_C, boot = boot.log, runs = boot)
    rli.i.D <- red::rli(grid.data.i$Category_D, boot = boot.log, runs = boot)
    rli.i.new <- red::rli(grid.data.i$Category.new, boot = boot.log, runs = boot)
    rli.i <- c(rli.i, rli.i.A, rli.i.B, rli.i.C, rli.i.D, rli.i.new)
    
    #Proportion of threatened
    threat.cats <- c("EN","VU","CR")
    rm_these <- !grid.data.i$Category %in% "NA"
    prop.i <- dim(grid.data.i[grid.data.i$Category %in% threat.cats & rm_these,])[1]/
                  dim(grid.data.i[rm_these,])[1]
    prop.i.A <- dim(grid.data.i[grid.data.i$Category_A %in% threat.cats & rm_these,])[1]/
                  dim(grid.data.i[rm_these,])[1]
    prop.i.B <- dim(grid.data.i[grid.data.i$Category_B %in% threat.cats & rm_these,])[1]/
                  dim(grid.data.i[rm_these,])[1]
    prop.i.C <- dim(grid.data.i[grid.data.i$Category_C %in% threat.cats & rm_these,])[1]/
                  dim(grid.data.i[rm_these,])[1]
    prop.i.D <- dim(grid.data.i[grid.data.i$Category_D %in% threat.cats & rm_these,])[1]/
                  dim(grid.data.i[rm_these,])[1]
    prop.i.new <- dim(grid.data.i[grid.data.i$Category %in% threat.cats,])[1]/
                    dim(grid.data.i)[1]
    prop.i <- c(prop.i, prop.i.A, prop.i.B, prop.i.C, prop.i.D, prop.i.new)
  } else {
    # rli.i <- c("LowCL" = NA,  "Median" = NA, "UpCL" = NA)
    rli.i <- c("RLI" = NA, "RLI_catA" = NA, "RLI_catB" = NA, 
               "RLI_catC" = NA, "RLI_catD" = NA, "RLI_new" = NA)
    prop.i <- c("Threat" = NA, "Threat_catA" = NA, "Threat_catB" = NA, 
               "Threat_catC" = NA, "Threat_catD" = NA, "Threat_new" = NA)
  }
  
  grid.data.i <- grid.data.i[!duplicated(grid.data.i$Name),]  
  if (dim(grid.data.i)[1] > 0) { # One record per species
    #RLI
    rli.i.pa <- red::rli(grid.data.i$Category, boot = boot.log, runs = boot)
    rli.i.A <- red::rli(grid.data.i$Category_A, boot = boot.log, runs = boot)
    rli.i.B <- red::rli(grid.data.i$Category_B, boot = boot.log, runs = boot)
    rli.i.C <- red::rli(grid.data.i$Category_C, boot = boot.log, runs = boot)
    rli.i.D <- red::rli(grid.data.i$Category_D, boot = boot.log, runs = boot)
    rli.i.new <- red::rli(grid.data.i$Category.new, boot = boot.log, runs = boot)
    rli.i.pa <- c(rli.i.pa, rli.i.A, rli.i.B, rli.i.C, rli.i.D, rli.i.new)
    
    #Proportion of threatened
    threat.cats <- c("EN","VU","CR")
    rm_these <- !grid.data.i$Category %in% "NA"
    prop.i.pa <- dim(grid.data.i[grid.data.i$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.A <- dim(grid.data.i[grid.data.i$Category_A %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.B <- dim(grid.data.i[grid.data.i$Category_B %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.C <- dim(grid.data.i[grid.data.i$Category_C %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.D <- dim(grid.data.i[grid.data.i$Category_D %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.new <- dim(grid.data.i[grid.data.i$Category %in% threat.cats,])[1]/
      dim(grid.data.i)[1]
    prop.i.pa <- c(prop.i.pa, prop.i.A, prop.i.B, prop.i.C, prop.i.D, prop.i.new)
    
  } else {
    # rli.i <- c("LowCL" = NA,  "Median" = NA, "UpCL" = NA)
    rli.i.pa <- c("RLI" = NA, "RLI_catA" = NA, "RLI_catB" = NA, 
               "RLI_catC" = NA, "RLI_catD" = NA, "RLI_new" = NA)
    prop.i.pa <- c("Threat" = NA, "Threat_catA" = NA, "Threat_catB" = NA, 
                "Threat_catC" = NA, "Threat_catD" = NA, "Threat_new" = NA)
  }
  
  ## ONLY ENDEMICS
  grid.data.ii <- grid.data.i[grid.data.i$Endemism %in% "endemic",]
  if (dim(grid.data.ii)[1] > 0) { # All records
    #RLI
    rli.i.end <- red::rli(grid.data.ii$Category, boot = boot.log, runs = boot)
    rli.i.end.A <- red::rli(grid.data.ii$Category_A, boot = boot.log, runs = boot)
    rli.i.end.B <- red::rli(grid.data.ii$Category_B, boot = boot.log, runs = boot)
    rli.i.end.C <- red::rli(grid.data.ii$Category_C, boot = boot.log, runs = boot)
    rli.i.end.D <- red::rli(grid.data.ii$Category_D, boot = boot.log, runs = boot)
    rli.i.end.new <- red::rli(grid.data.ii$Category.new, boot = boot.log, runs = boot)
    rli.i.end <- c(rli.i.end, rli.i.end.A, rli.i.end.B, rli.i.end.C, rli.i.end.D, rli.i.end.new)
    
    #Proportion of threatened
    threat.cats <- c("EN","VU","CR")
    rm_these <- !grid.data.ii$Category %in% "NA"
    prop.i.end <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.A <- dim(grid.data.ii[grid.data.ii$Category_A %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.B <- dim(grid.data.ii[grid.data.ii$Category_B %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.C <- dim(grid.data.ii[grid.data.ii$Category_C %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.D <- dim(grid.data.ii[grid.data.ii$Category_D %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.new <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats,])[1]/
      dim(grid.data.ii)[1]
    prop.i.end <- c(prop.i.end, prop.i.end.A, prop.i.end.B, prop.i.end.C, prop.i.end.D, prop.i.end.new)
  } else {
    # rli.i.end <- c("LowCL" = NA,  "Median" = NA, "UpCL" = NA)
    rli.i.end <- c("RLI" = NA, "RLI_catA" = NA, "RLI_catB" = NA, 
               "RLI_catC" = NA, "RLI_catD" = NA, "RLI_new" = NA)
    prop.i.end <- c("Threat" = NA, "Threat_catA" = NA, "Threat_catB" = NA, 
                "Threat_catC" = NA, "Threat_catD" = NA, "Threat_new" = NA)
  }
  
  grid.data.ii <- grid.data.ii[!duplicated(grid.data.ii$Name),]  
  if (dim(grid.data.i)[1] > 0) { # One record per species
    #RLI
    rli.i.end.pa <- red::rli(grid.data.ii$Category, boot = boot.log, runs = boot)
    rli.i.end.A <- red::rli(grid.data.ii$Category_A, boot = boot.log, runs = boot)
    rli.i.end.B <- red::rli(grid.data.ii$Category_B, boot = boot.log, runs = boot)
    rli.i.end.C <- red::rli(grid.data.ii$Category_C, boot = boot.log, runs = boot)
    rli.i.end.D <- red::rli(grid.data.ii$Category_D, boot = boot.log, runs = boot)
    rli.i.end.new <- red::rli(grid.data.ii$Category.new, boot = boot.log, runs = boot)
    rli.i.end.pa <- c(rli.i.end.pa, rli.i.end.A, rli.i.end.B, rli.i.end.C, rli.i.end.D, rli.i.end.new)
    
    #Proportion of threatened
    threat.cats <- c("EN","VU","CR")
    prop.i.end.pa <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.A <- dim(grid.data.ii[grid.data.ii$Category_A %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.B <- dim(grid.data.ii[grid.data.ii$Category_B %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.C <- dim(grid.data.ii[grid.data.ii$Category_C %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.D <- dim(grid.data.ii[grid.data.ii$Category_D %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.new <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats,])[1]/
      dim(grid.data.ii)[1]
    prop.i.end.pa <- c(prop.i.end.pa, prop.i.end.A, prop.i.end.B, prop.i.end.C, prop.i.end.D, prop.i.end.new)
    
  } else {
    rli.i.end.pa <- c("RLI" = NA, "RLI_catA" = NA, "RLI_catB" = NA, 
                  "RLI_catC" = NA, "RLI_catD" = NA, "RLI_new" = NA)
    prop.i.end.pa <- c("Threat" = NA, "Threat_catA" = NA, "Threat_catB" = NA, 
                    "Threat_catC" = NA, "Threat_catD" = NA, "Threat_new" = NA)
  }
  
  res <- c(rli.i, rli.i.end, rli.i.pa, rli.i.end.pa)
  names(res) <- c("RLI", "RLI_catA", "RLI_catB", "RLI_catC", "RLI_catD", "RLI_new",
                  "RLI.end", "RLI_catA.end", "RLI_catB.end", "RLI_catC.end", "RLI_catD.end", "RLI_new.end",
                  "RLI.pa", "RLI_catA.pa", "RLI_catB.pa", "RLI_catC.pa", "RLI_catD.pa", "RLI_new.pa",
                  "RLI.end.pa", "RLI_catA.end.pa", "RLI_catB.end.pa", "RLI_catC.end.pa", "RLI_catD.end.pa", "RLI_new.end.pa")
  resultado[[i]] <- res
  
  res1 <- c(prop.i, prop.i.end, prop.i.pa, prop.i.end.pa)
  names(res1) <- c("Threat", "Threat_catA", "Threat_catB", "Threat_catC", "Threat_catD", "Threat_new",
                  "Threat.end", "Threat_catA.end", "Threat_catB.end", "Threat_catC.end", "Threat_catD.end", "Threat_new.end",
                  "Threat.pa", "Threat_catA.pa", "Threat_catB.pa", "Threat_catC.pa", "Threat_catD.pa", "Threat_new.pa",
                  "Threat.end.pa", "Threat_catA.end.pa", "Threat_catB.end.pa", "Threat_catC.end.pa", "Threat_catD.end.pa", "Threat_new.end.pa")
  resultado1[[i]] <- res1
  
  
}
#RLI
all <- do.call(rbind.data.frame, resultado)
names(all) <- names(resultado[[1]])
all$ID <- names(resultado)
head(all)
apply(all[,-c(25)], 2, summary)
#Proportion of threatened
all1 <- do.call(rbind.data.frame, resultado1)
names(all1) <- names(resultado1[[1]])
all1$ID <- names(resultado1)
head(all1)
apply(all1[,-c(25)], 2, summary)
#merging
all <- merge(all, all1, 
                     by= "ID", all = TRUE, sort = FALSE)
grid.result <- merge(grid.result, all, 
                     by= "ID", all.x = TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]


#### Final edits of the grid summary ####
## NO NEED - EFFECT WERE ONLY MARGINAL ##

# ## Is there an effect of cell area on N.total and RLI? Yes, for all variables!
# toto = grid.result[grid.result$area_ha<500000,]
# toto = grid.result[!is.na(grid.result$area_ha),]
# par(mfrow=c(1,1))
# plot(log(toto$N.total) ~ toto$area_ha)
# abline(lm(log(toto$N.total) ~ toto$area_ha))
# plot(log(toto$S.total) ~ toto$area_ha)
# abline(lm(log(toto$S.total) ~ toto$area_ha))
# plot(toto$RLI ~ toto$area_ha)
# abline(lm(toto$RLI ~ toto$area_ha))
# plot(log(toto$N.total) ~ toto$area_ha.clip)
# abline(lm(log(toto$N.total) ~ toto$area_ha.clip))
# plot(toto$RLI ~ toto$area_ha.clip)
# abline(lm(toto$RLI ~ toto$area_ha.clip))
# 
# 
# ### Defining the relationship between cell area and sampling intensity and correcting the diversity estimates ###
# #function that perform the corrections
# # dados = grid.result
# # lim.area = 500000
# # nome.N = "N.total"
# # nome.est = "RLI"
# 
# correct.estimates = function(dados,lim.area,nome.N,nome.est,...) {
#   # data for cells < lim.area
#   toto = dados[dados$area_ha < lim.area,]
#   N = log(toto[,nome.N])
#   ids = !is.na(N) & !is.infinite(N)
#   area = toto$area_ha[ids]
#   N = N[ids]
#   # fitting and comparing the models: N ~ area
#   mod.N = lm(N ~ area)
#   mod.N1 = lm(N ~ area + I(area^2))
#   aic = bbmle::AICtab(mod.N,mod.N1,sort=FALSE)
#   if(aic$dAIC[1]>log(8)) mod.final = mod.N1 else mod.final = mod.N
#   par(mfrow=c(1,3), las=1,mar=c(4,4,1,0.5),mgp=c(2,0.5,0),tcl=-.3,...)
#   plot(N ~ area)
#   curve(predict(mod.final, newdata= data.frame(area= x)),add=TRUE, lwd=2, col=2)
#   # calculating the corrected sampling effort per grid cell
#   if (aic$dAIC[1] > log(8)) { 
#     N.est = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*area + coef(mod.final)[3]*area^2),2)
#     N.max = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*1000000 + coef(mod.final)[3]*1000000^2),2)
#   } else {
#     N.est = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*area),2)
#     N.max = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*1000000),2)
#   }
#   prop = (N.max - N.est)/N.est
#   decay = try(nls(prop ~ a*exp(b*area), start=list(a=-0.2,b=0.00001)),TRUE)
#   if (class(decay) == "try-error") decay = try(nls(prop~a*exp(b*area),start=list(a=-0.1,b=0.000001)),TRUE)
#   plot(prop ~ area);   curve(coef(decay)[1]*exp(coef(decay)[2]*x), lwd=2, col=2, add=T)
#   prop[prop < 0.0] = coef(decay)[1] * exp(coef(decay)[2] * area[prop < 0.0]) * 0.5
#   N.correct = toto[,nome.N][ids] + round(prop*toto[,nome.N][ids],0)
#   
#   # data for cells < 250000
#   toto = dados[dados$area_ha > lim.area,]
#   ids = !is.na(toto[,nome.N]) & !is.na(toto[,nome.est]) & !is.infinite(log(toto[,nome.N])) & !is.infinite(log(toto[,nome.est]))
#   N1 = log(toto[,nome.N])[ids]
#   est = toto[,nome.est][ids]
#   w = toto$SampCover[ids]
#   # fitting and comparing the models: N ~ area
#   mod.est  = lm(est ~ N1, weights = w^3)
#   mod.est1 = lm(est ~ N1 + I(N1^2), weights = w^3)
#   aic = bbmle::AICtab(mod.est,mod.est1,sort=FALSE)
#   if (aic$dAIC[1]>log(8)) mod.final1 = mod.est1 else mod.final1 = mod.est
#   plot(est ~ N1)
#   curve(predict(mod.final1, newdata= data.frame(N1= x)),add=TRUE, lwd=2, col=2)
#   # calculating the corrected sampling effort per grid cell
#   if (aic$dAIC[1]>log(8)) { 
#     est.pred = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*N + coef(mod.final1)[3]*N^2),4)
#     est.est = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*log(N.correct) + coef(mod.final1)[3]*log(N.correct)^2),4)
#   } else {
#     est.pred = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*N),4)
#     est.est = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*log(N.correct)),4)
#   }
#   prop1 = (est.est - est.pred)/est.est
#   decay = try(nls(prop1~a*exp(b*area),start=list(a=0.001,b=0.000001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=5,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=1,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=0.1,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=-0.2,b=-0.00001)),TRUE)
#   # plot(prop1 ~ area);   curve(coef(decay)[1]*exp(coef(decay)[2]*x), lwd=2, col=2, add=T)
#   # prop1[prop1<0.0] = coef(decay)[1]*exp(coef(decay)[2]*area[prop1<0.0])*0.1
#   prop1[prop1>0.0] = coef(decay)[1]*exp(coef(decay)[2]*area[prop1>0.0])*0.1
#   
#   # plot(prop1 ~ area)
#   est.correct = dados[,nome.est]
#   est.correct[dados$area_ha<lim.area & !is.na(dados[,nome.N])] = est.correct[dados$area_ha<lim.area & !is.na(dados[,nome.N])] +
#     prop1 * est.correct[dados$area_ha<lim.area & !is.na(dados[,nome.N])]
#   plot(est.correct ~ dados[,nome.est], xlab="estimate",ylab="corrected est.")
#   legend("topleft", nome.est, cex=1.2, bty="n");abline(0,1,lty=2)
#   par(mfrow=c(1,1))
#   return(est.correct)
# }
# 
# #generating the corrected estimates for each group
# grid.result$Median.RLI.correct = correct.estimates(grid.result,215000,"N.total","Median.RLI")
# grid.result$Median.RLI.end.correct = correct.estimates(grid.result,215000,"N.total","Median.RLI.end")


#### SAME CALCULATION BUT ONLY FOR THE SPECIES CATEGORIZED AS CR ####
cr_spp <- occs1$Category %in% "CR"
### Number of occurrences, number of species and sampling coverage per site ###
tmp = iNEXT::DataInfo(as.data.frame.matrix(table(occs1$Name[cr_spp], 
                                                 occs1$ID[cr_spp]))) # getting metrics
tmp$SC[tmp$n==1] = 0.0001 # correcting values of SampCover for sites with N<10 and S ~~ N
tmp$SC[tmp$n==2&tmp$S.obs==2] = 0.0005; tmp$SC[tmp$n==2&tmp$S.obs==1] = 0.001
tmp$SC[tmp$n==3&tmp$S.obs==3] = 0.001; tmp$SC[tmp$n==3&tmp$S.obs==2] = 0.005; tmp$SC[tmp$n==3&tmp$S.obs==1] = 0.01
tmp$SC[tmp$n==4&tmp$S.obs==4] = 0.01; tmp$SC[tmp$n==4&tmp$S.obs==3] = 0.025; tmp$SC[tmp$n==4&tmp$S.obs==2] = 0.05; tmp$SC[tmp$n==4&tmp$S.obs==1] = 0.1
tmp$SC[tmp$n==5&tmp$S.obs==5] = 0.1; tmp$SC[tmp$n==5&tmp$S.obs==4] = 0.125; tmp$SC[tmp$n==5&tmp$S.obs==3] = 0.15; tmp$SC[tmp$n==5&tmp$S.obs==2] = 0.175; tmp$SC[tmp$n==5&tmp$S.obs==1] = 0.2
plot(tmp$SC ~ log(tmp$n))
cols.cr <- c("ID","N.total_CR","S.total_CR","SampCover_CR")
names(tmp)[1:4] <- cols.cr
grid.result <- merge(grid.result, tmp[,cols.cr], 
                     by.x= "ID", by.y= "ID", all.x=TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

### Number of occurrences, number of species and sampling coverage per site ###
endemics <- occs1$Endemism %in% "endemic" & occs1$Category %in% "CR"
tmp = iNEXT::DataInfo(as.data.frame.matrix(table(occs1$Name[endemics], occs1$ID[endemics]))) # getting metrics
tmp$SC[tmp$n==1] = 0.0001 # correcting values of SampCover for sites with N<10 and S ~~ N
tmp$SC[tmp$n==2&tmp$S.obs==2] = 0.0005; tmp$SC[tmp$n==2&tmp$S.obs==1] = 0.001
tmp$SC[tmp$n==3&tmp$S.obs==3] = 0.001; tmp$SC[tmp$n==3&tmp$S.obs==2] = 0.005; tmp$SC[tmp$n==3&tmp$S.obs==1] = 0.01
tmp$SC[tmp$n==4&tmp$S.obs==4] = 0.01; tmp$SC[tmp$n==4&tmp$S.obs==3] = 0.025; tmp$SC[tmp$n==4&tmp$S.obs==2] = 0.05; tmp$SC[tmp$n==4&tmp$S.obs==1] = 0.1
tmp$SC[tmp$n==5&tmp$S.obs==5] = 0.1; tmp$SC[tmp$n==5&tmp$S.obs==4] = 0.125; tmp$SC[tmp$n==5&tmp$S.obs==3] = 0.15; tmp$SC[tmp$n==5&tmp$S.obs==2] = 0.175; tmp$SC[tmp$n==5&tmp$S.obs==1] = 0.2
tmp$SC[tmp$n==6&tmp$S.obs==6] = 0.125; tmp$SC[tmp$n==6&tmp$S.obs==5] = 0.15; tmp$SC[tmp$n==6&tmp$S.obs==4] = 0.175; tmp$SC[tmp$n==6&tmp$S.obs==3] = 0.2; tmp$SC[tmp$n==6&tmp$S.obs==2] = 0.225; tmp$SC[tmp$n==6&tmp$S.obs==1] = 0.25
tmp$SC[tmp$n==7&tmp$S.obs==7] = 0.15; tmp$SC[tmp$n==7&tmp$S.obs==6] = 0.175; tmp$SC[tmp$n==7&tmp$S.obs==5] = 0.2; tmp$SC[tmp$n==7&tmp$S.obs==4] = 0.225; tmp$SC[tmp$n==7&tmp$S.obs==3] = 0.25; tmp$SC[tmp$n==7&tmp$S.obs==2] = 0.275; tmp$SC[tmp$n==7&tmp$S.obs==1] = 0.3
tmp$SC[tmp$n==8&tmp$S.obs==8] = 0.175; tmp$SC[tmp$n==8&tmp$S.obs==7] = 0.2; tmp$SC[tmp$n==8&tmp$S.obs==6] = 0.225; tmp$SC[tmp$n==8&tmp$S.obs==5] = 0.25; tmp$SC[tmp$n==8&tmp$S.obs==4] = 0.275; tmp$SC[tmp$n==8&tmp$S.obs==3] = 0.3; tmp$SC[tmp$n==8&tmp$S.obs==2] = 0.325; tmp$SC[tmp$n==8&tmp$S.obs==1] = 0.35
tmp$SC[tmp$n==9&tmp$S.obs==9] = 0.2; tmp$SC[tmp$n==9&tmp$S.obs==8] = 0.225; tmp$SC[tmp$n==9&tmp$S.obs==7] = 0.25; tmp$SC[tmp$n==9&tmp$S.obs==6] = 0.275; tmp$SC[tmp$n==9&tmp$S.obs==5] = 0.3; tmp$SC[tmp$n==9&tmp$S.obs==4] = 0.325; tmp$SC[tmp$n==9&tmp$S.obs==3] = 0.35; tmp$SC[tmp$n==9&tmp$S.obs==2] = 0.375; tmp$SC[tmp$n==9&tmp$S.obs==1] = 0.4
tmp$SC[tmp$n==10&tmp$S.obs==3] = 0.225
tmp$SC[tmp$n==12&tmp$S.obs==6] = 0.325; tmp$SC[tmp$n==12&tmp$S.obs==5] = 0.35; tmp$SC[tmp$n==12&tmp$S.obs==4] = 0.375
tmp$SC[tmp$n==15&tmp$S.obs==4] = 0.425
tmp$SC[tmp$n==20&tmp$S.obs==6] = 0.475
tmp$SC[tmp$n==24&tmp$S.obs==8] = 0.500
tmp$SC[tmp$n==25&tmp$S.obs==8] = 0.515
# tmp$SC[tmp$n==33&tmp$S.obs==7] = 0.715
plot(tmp$SC ~ log(tmp$n))
cols.cr <- c("ID","N.total.end_CR","S.total.end_CR","SampCover.end_CR")
names(tmp)[1:4] = cols.cr
grid.result <- merge(grid.result, tmp[,cols.cr], 
                     by.x= "ID", by.y= "ID", all.x=TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]


### Calculating the % of Threat per grid cell
grid.data <- occs1
grid.data$Category.new <- grid.data$Category
grid.data$Category.new[grid.data$Category %in% "NA"] <- "LC"
resultado <- vector("list", dim(grid.result)[1])
names(resultado) <- grid.result$ID
resultado1 <- resultado
# boot = 4999
# boot.log = FALSE
for (i in 1: length(resultado)) {
  nome.i <- names(resultado)[i]
  
  ## ALL POPULATIONS
  grid.data.i <- grid.data[grid.data$ID %in% nome.i,]  
  if (dim(grid.data.i)[1] > 0) { # All records

    #Proportion of threatened
    threat.cats <- c("CR")
    rm_these <- !grid.data.i$Category %in% "NA"
    prop.i <- dim(grid.data.i[grid.data.i$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.new <- dim(grid.data.i[grid.data.i$Category %in% threat.cats,])[1]/
      dim(grid.data.i)[1]
    prop.i <- c(prop.i, prop.i.new)
  } else {
    prop.i <- c("Threat_CR" = NA, "Threat_new_CR" = NA)
  }
  
  grid.data.i <- grid.data.i[!duplicated(grid.data.i$Name),]  
  if (dim(grid.data.i)[1] > 0) { # One record per species

    #Proportion of threatened
    threat.cats <- c("CR")
    rm_these <- !grid.data.i$Category %in% "NA"
    prop.i.pa <- dim(grid.data.i[grid.data.i$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.i[rm_these,])[1]
    prop.i.new <- dim(grid.data.i[grid.data.i$Category %in% threat.cats,])[1]/
      dim(grid.data.i)[1]
    prop.i.pa <- c(prop.i.pa, prop.i.new)
    
  } else {
    prop.i.pa <- c("Threat_CR" = NA, "Threat_new_CR" = NA)
  }
  
  ## ONLY ENDEMICS
  grid.data.ii <- grid.data.i[grid.data.i$Endemism %in% "endemic",]
  if (dim(grid.data.ii)[1] > 0) { # All records

    #Proportion of threatened
    threat.cats <- c("CR")
    rm_these <- !grid.data.ii$Category %in% "NA"
    prop.i.end <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.new <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats,])[1]/
      dim(grid.data.ii)[1]
    prop.i.end <- c(prop.i.end, prop.i.end.new)
  } else {
    prop.i.end <- c("Threat_CR" = NA, "Threat_new_CR" = NA)
  }
  
  grid.data.ii <- grid.data.ii[!duplicated(grid.data.ii$Name),]  
  if (dim(grid.data.i)[1] > 0) { # One record per species

    #Proportion of threatened
    threat.cats <- c("CR")
    prop.i.end.pa <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats & rm_these,])[1]/
      dim(grid.data.ii[rm_these,])[1]
    prop.i.end.new <- dim(grid.data.ii[grid.data.ii$Category %in% threat.cats,])[1]/
      dim(grid.data.ii)[1]
    prop.i.end.pa <- c(prop.i.end.pa, prop.i.end.new)
    
  } else {
    prop.i.end.pa <- c("Threat_CR" = NA, "Threat_new_CR" = NA)
  }
  
  res1 <- c(prop.i, prop.i.end, prop.i.pa, prop.i.end.pa)
  names(res1) <- c("Threat_CR", "Threat_new_CR",
                   "Threat.end_CR", "Threat_new.end_CR",
                   "Threat.pa_CR", "Threat_new.pa_CR",
                   "Threat.end.pa_CR", "Threat_new.end.pa_CR")
  resultado1[[i]] <- res1
}
#Proportion of threatened
all1 <- do.call(rbind.data.frame, resultado1)
names(all1) <- names(resultado1[[1]])
all1$ID <- names(resultado1)
head(all1)
apply(all1[,-c(9)], 2, summary)
#merging
grid.result <- merge(grid.result, all1, 
                     by= "ID", all.x = TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

## Saving the grid results
saveRDS(grid.result, "data/grid.results_adaptive_resol.rds")



# -----------------------------------------------------------------------------
########################
#### HEXAGONAL GRID #### 
########################
rm(list=ls())
gc()

## LOADING THE DATA  ##
## Threat occurrence data ##
occs <- readRDS("data/threat_occ_data_final.rds")
occs <- occs[,c("tax","ddlat","ddlon","tax.check2","tax.check.final")]
#removing possible missing coordinates
occs$ddlat <- as.numeric(occs$ddlat)
occs$ddlon <- as.numeric(occs$ddlon)
occs <- occs[!is.na(ddlat) | !is.na(ddlon),]

## Threat assessments
all.crit <- readRDS("data/all.criteria.rds")
# Getting NT for each criteria
subcriteria <- c("category_A", "category_B","category_C","category_D")
subcriteria1 <- list(c("A2"), c("B1","B2"),c("C1", "C2"), "D")
for(i in 1:length(subcriteria)) {
  all.crit[,subcriteria[i]] <- 
    ConR::near.threatened(all.crit[,subcriteria[i]],
                          all.crit$EOO, all.crit$AOO, all.crit$declineB,
                          all.crit$reduction_A12, all.crit$pop.size, all.crit$pop.size.low,
                          all.crit$nbe_loc_total, all.crit$sever.frag, ext.fluct = NULL,
                          subpop = all.crit$Nbe_subPop, subcriteria = subcriteria1[[i]])
}
# Filtering and merging the two datasets
colunas <- c("species", "endemic", "cat.reg.clean",
             "category_A", "category_B","category_C","category_D")
occs <- merge(occs, all.crit[, colunas], by.x = "tax", by.y = "species")

##Transforming the occurrences into spatial points
occs.sp <- SpatialPointsDataFrame(coords = cbind(occs$ddlon, occs$ddlat),
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs"),
                                 data = occs[,c("tax","tax.check2","tax.check.final", "endemic","cat.reg.clean",
                                                "category_A", "category_B","category_C","category_D")])
occs.sp$ordem <- 1:dim(occs.sp)[1]
occs.sp <- spTransform(occs.sp, CRS("+init=epsg:5641"))

## Getting the grid and editing it
grid <- readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//af_hex_grid_50km.rds")

## Making sure grid IDs are the same internally
grid1 <- grid
tmp <- SpatialPoints(coords = coordinates(grid), proj4string = raster::crs(grid))
tmp1 <- over(tmp, grid)
for(i in 1:length(grid1)) { slot(slot(grid1, "polygons")[[i]],"ID") = as.character(paste("ID",i,sep=""))} 

## Projecting the Grid
grid1 <- spTransform(grid1, CRS("+init=epsg:5641"))
# plot(grid1)

## Overlaying the occurences against the grid
occs.sp@data$cell.ID = paste("ID", over(occs.sp, grid1),sep="")

## Re-assing grid ID for points falling near the coast
## Or finding the grid cell for points with projection issues
tmp <- occs.sp[occs.sp@data$cell.ID %in% "IDNA",]

#First, lets remove points outside the AF extent
ext <- raster::extent(grid1) + c(-50000,50000,-50000,50000)
tmp <- raster::crop(tmp, ext)

#Second, lets remove points too far away from the AF grid
path <- "E://ownCloud//W_GIS" # Renato's path
am.lat <- rgdal::readOGR(dsn=paste(path,"//Am_Lat_ADM_ArcGis",sep=""),layer="LatinAmerica")
#s.am <- readRDS(paste(path,"//Am_Lat_ADM_GADM_v3.6//gadm36_SouthAm_0_sp.rds",sep=""))
brasil = readRDS(paste(path,"//Am_Lat_ADM_GADM_v3.6//gadm36_BRA_1_sp.rds",sep=""))
af <- rgdal::readOGR(dsn=paste(path,"//AF_limites_milton",sep=""),layer="merge_limites_MA11428_TNC_ARG_BR_PAR")
af <- rgeos::gBuffer(af, byid=TRUE, width=0) #correcting possible overlapping polygons
af <- cleangeo::clgeo_Clean(af) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
af <- raster::aggregate(af, dissolve=TRUE)
#creating the Atlantic Forest 0.1 degree (~10 km) buffer
af.buf <- rgeos::gBuffer(af, byid=TRUE, width = 0.10)
af.buf <- raster::disaggregate(af.buf)
af.buf <- cleangeo::clgeo_Clean(af.buf)
#Removing buffer inner-holes
toto <- SpatialPolygons(list(Polygons(list(af.buf@polygons[[1]]@Polygons[[1]]),ID=1)))
crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")  # geographical, datum WGS84
proj4string(toto) <- crs.geo
af.buf <- rbind(af.buf[2:9], toto)
af.buf <- raster::aggregate(af.buf)
#projecting
am.lat.proj= spTransform(am.lat,  CRS("+init=epsg:5641"))
brasil.proj= spTransform(brasil,  CRS("+init=epsg:5641"))
af.buf.proj = spTransform(af.buf, CRS("+init=epsg:5641"))
af.proj = spTransform(af, CRS("+init=epsg:5641"))

# Cropping
# plot(af.buf.proj)
# plot(tmp, add= TRUE)
tmp <- raster::crop(tmp, af.buf.proj)
#calculating the distance between points and all grid cells
# plot(grid1)
# plot(tmp, add= TRUE)
tmp1 <- rgeos::gDistance(tmp, grid1, byid=TRUE)
result = cbind.data.frame(tax = as.character(tmp$tax), 
                          cell.ID = "IDNA", stringsAsFactors = FALSE)
for(i in 1:dim(tmp1)[2]){
  IDi = which.min(tmp1[,i])
  dist1 = tmp1[which.min(tmp1[,i]),i]
  if(dist1>0 & dist1<5000) result$cell.ID[i] = as.character(paste("ID", IDi, sep = ""))
}
occs.sp@data$cell.ID[match(tmp$ordem,occs.sp@data$ordem)] <- result$cell.ID

## Removing points falling outside the grid limits
occs.sp <- occs.sp[!occs.sp@data$cell.ID %in% "IDNA", ]

#### SUMMARIES REPRESENTATIVENESS OF THE SAMPLING EFFORT PER GRID CELL FOR EACH SPECIES GROUP ####
source("R//polygon_centroid.R")
grid.center <- polygonCenter(grid1)

## Getting the summaries per grid cell
grid.result = cbind.data.frame(order = 1:length(grid1), 
                               ID = as.character(sapply(slot(grid1, "polygons"), function(x) slot(x, "ID"))),
                               plotOrder = slot(grid1, "plotOrder"),
                               area_ha = as.double(rgeos::gArea(grid1, byid = TRUE)/10000),
                               # lat = sp::coordinates(grid1)[,2],
                               # long = sp::coordinates(grid1)[,1],
                               lat = sp::coordinates(rgeos::gCentroid(grid1, byid = TRUE))[,2],
                               long = sp::coordinates(rgeos::gCentroid(grid1, byid = TRUE))[,1],
                               lat1 = grid.center[,2], 
                               long1 = grid.center[,1])
grid.result$ID = as.character(grid.result$ID)

### Number of occurrences, number of species and sampling coverage per site ###
tmp = iNEXT::DataInfo(as.data.frame.matrix(table(occs.sp$tax, occs.sp$cell.ID))) # getting metrics
tmp$SC[tmp$n==1] = 0.0001 # correcting values of SampCover for sites with N<10 and S ~~ N
tmp$SC[tmp$n==2&tmp$S.obs==2] = 0.0005; tmp$SC[tmp$n==2&tmp$S.obs==1] = 0.001
tmp$SC[tmp$n==3&tmp$S.obs==3] = 0.001; tmp$SC[tmp$n==3&tmp$S.obs==2] = 0.005; tmp$SC[tmp$n==3&tmp$S.obs==1] = 0.01
tmp$SC[tmp$n==4&tmp$S.obs==4] = 0.01; tmp$SC[tmp$n==4&tmp$S.obs==3] = 0.025; tmp$SC[tmp$n==4&tmp$S.obs==2] = 0.05; tmp$SC[tmp$n==4&tmp$S.obs==1] = 0.1
tmp$SC[tmp$n==5&tmp$S.obs==5] = 0.1; tmp$SC[tmp$n==5&tmp$S.obs==4] = 0.125; tmp$SC[tmp$n==5&tmp$S.obs==3] = 0.15; tmp$SC[tmp$n==5&tmp$S.obs==2] = 0.175; tmp$SC[tmp$n==5&tmp$S.obs==1] = 0.2
tmp$SC[tmp$n==6&tmp$S.obs==6] = 0.125; tmp$SC[tmp$n==6&tmp$S.obs==5] = 0.15; tmp$SC[tmp$n==6&tmp$S.obs==4] = 0.175; tmp$SC[tmp$n==6&tmp$S.obs==3] = 0.2; tmp$SC[tmp$n==6&tmp$S.obs==2] = 0.225; tmp$SC[tmp$n==6&tmp$S.obs==1] = 0.25
tmp$SC[tmp$n==7&tmp$S.obs==7] = 0.15; tmp$SC[tmp$n==7&tmp$S.obs==6] = 0.175; tmp$SC[tmp$n==7&tmp$S.obs==5] = 0.2; tmp$SC[tmp$n==7&tmp$S.obs==4] = 0.225; tmp$SC[tmp$n==7&tmp$S.obs==3] = 0.25; tmp$SC[tmp$n==7&tmp$S.obs==2] = 0.275; tmp$SC[tmp$n==7&tmp$S.obs==1] = 0.3
tmp$SC[tmp$n==8&tmp$S.obs==8] = 0.175; tmp$SC[tmp$n==8&tmp$S.obs==7] = 0.2; tmp$SC[tmp$n==8&tmp$S.obs==6] = 0.225; tmp$SC[tmp$n==8&tmp$S.obs==5] = 0.25; tmp$SC[tmp$n==8&tmp$S.obs==4] = 0.275; tmp$SC[tmp$n==8&tmp$S.obs==3] = 0.3; tmp$SC[tmp$n==8&tmp$S.obs==2] = 0.325; tmp$SC[tmp$n==8&tmp$S.obs==1] = 0.35
tmp$SC[tmp$n==9&tmp$S.obs==9] = 0.2; tmp$SC[tmp$n==9&tmp$S.obs==8] = 0.225; tmp$SC[tmp$n==9&tmp$S.obs==7] = 0.25; tmp$SC[tmp$n==9&tmp$S.obs==6] = 0.275; tmp$SC[tmp$n==9&tmp$S.obs==5] = 0.3; tmp$SC[tmp$n==9&tmp$S.obs==4] = 0.325; tmp$SC[tmp$n==9&tmp$S.obs==3] = 0.35; tmp$SC[tmp$n==9&tmp$S.obs==2] = 0.375; tmp$SC[tmp$n==9&tmp$S.obs==1] = 0.4
tmp$SC[tmp$n==12&tmp$S.obs==5]= 0.35
#plot(tmp$SC ~ log(tmp$n))
samp.cov.cutoff = quantile(tmp$SC[tmp$n>=10], prob=c(0.05,0.1,0.25,0.5,0.75,0.8,0.85,0.9,0.95,0.99))
names(tmp)[1:4] = c("ID","N.total","S.total","SampCover")
grid.result <- merge(grid.result, tmp[,c("ID","N.total","S.total","SampCover")], 
                     by.x= "ID", by.y= "ID", all.x=TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

### Number of occurrences, number of species and sampling coverage per site ###
endemics <- occs.sp$endemic %in% "endemic"
tmp = iNEXT::DataInfo(as.data.frame.matrix(table(occs.sp$tax[endemics], occs.sp$cell.ID[endemics]))) # getting metrics
tmp$SC[tmp$n==1] = 0.0001 # correcting values of SampCover for sites with N<10 and S ~~ N
tmp$SC[tmp$n==2&tmp$S.obs==2] = 0.0005; tmp$SC[tmp$n==2&tmp$S.obs==1] = 0.001
tmp$SC[tmp$n==3&tmp$S.obs==3] = 0.001; tmp$SC[tmp$n==3&tmp$S.obs==2] = 0.005; tmp$SC[tmp$n==3&tmp$S.obs==1] = 0.01
tmp$SC[tmp$n==4&tmp$S.obs==4] = 0.01; tmp$SC[tmp$n==4&tmp$S.obs==3] = 0.025; tmp$SC[tmp$n==4&tmp$S.obs==2] = 0.05; tmp$SC[tmp$n==4&tmp$S.obs==1] = 0.1
tmp$SC[tmp$n==5&tmp$S.obs==5] = 0.1; tmp$SC[tmp$n==5&tmp$S.obs==4] = 0.125; tmp$SC[tmp$n==5&tmp$S.obs==3] = 0.15; tmp$SC[tmp$n==5&tmp$S.obs==2] = 0.175; tmp$SC[tmp$n==5&tmp$S.obs==1] = 0.2
tmp$SC[tmp$n==6&tmp$S.obs==6] = 0.125; tmp$SC[tmp$n==6&tmp$S.obs==5] = 0.15; tmp$SC[tmp$n==6&tmp$S.obs==4] = 0.175; tmp$SC[tmp$n==6&tmp$S.obs==3] = 0.2; tmp$SC[tmp$n==6&tmp$S.obs==2] = 0.225; tmp$SC[tmp$n==6&tmp$S.obs==1] = 0.25
tmp$SC[tmp$n==7&tmp$S.obs==7] = 0.15; tmp$SC[tmp$n==7&tmp$S.obs==6] = 0.175; tmp$SC[tmp$n==7&tmp$S.obs==5] = 0.2; tmp$SC[tmp$n==7&tmp$S.obs==4] = 0.225; tmp$SC[tmp$n==7&tmp$S.obs==3] = 0.25; tmp$SC[tmp$n==7&tmp$S.obs==2] = 0.275; tmp$SC[tmp$n==7&tmp$S.obs==1] = 0.3
tmp$SC[tmp$n==8&tmp$S.obs==8] = 0.175; tmp$SC[tmp$n==8&tmp$S.obs==7] = 0.2; tmp$SC[tmp$n==8&tmp$S.obs==6] = 0.225; tmp$SC[tmp$n==8&tmp$S.obs==5] = 0.25; tmp$SC[tmp$n==8&tmp$S.obs==4] = 0.275; tmp$SC[tmp$n==8&tmp$S.obs==3] = 0.3; tmp$SC[tmp$n==8&tmp$S.obs==2] = 0.325; tmp$SC[tmp$n==8&tmp$S.obs==1] = 0.35
tmp$SC[tmp$n==9&tmp$S.obs==9] = 0.2; tmp$SC[tmp$n==9&tmp$S.obs==8] = 0.225; tmp$SC[tmp$n==9&tmp$S.obs==7] = 0.25; tmp$SC[tmp$n==9&tmp$S.obs==6] = 0.275; tmp$SC[tmp$n==9&tmp$S.obs==5] = 0.3; tmp$SC[tmp$n==9&tmp$S.obs==4] = 0.325; tmp$SC[tmp$n==9&tmp$S.obs==3] = 0.35; tmp$SC[tmp$n==9&tmp$S.obs==2] = 0.375; tmp$SC[tmp$n==9&tmp$S.obs==1] = 0.4
tmp$SC[tmp$n==12&tmp$S.obs==5]= 0.35
#plot(tmp$SC ~ log(tmp$n))
names(tmp)[1:4] = c("ID","N.total.end","S.total.end","SampCover.end")
grid.result <- merge(grid.result, tmp[,c("ID","N.total.end","S.total.end","SampCover.end")], 
                     by.x= "ID", by.y= "ID", all.x=TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

### Calculating the RLI per grid cell
grid.data <- occs.sp@data
grid.data$cat.reg.clean.new <- grid.data$cat.reg.clean
grid.data$cat.reg.clean.new[grid.data$cat.reg.clean %in% "NA"] <- "LC"
resultado <- vector("list", dim(grid.result)[1])
names(resultado) <- grid.result$ID

boot = 499
boot.log = FALSE
for (i in 1: length(resultado)) {
  nome.i <- names(resultado)[i]
  grid.data.i <- grid.data[grid.data$cell.ID == nome.i,]  
  if (dim(grid.data.i)[1] > 0) {
    rli.i <- red::rli(grid.data.i$cat.reg.clean, boot = boot.log, runs = boot)
    rli.i.A <- red::rli(grid.data.i$category_A, boot = boot.log, runs = boot)
    rli.i.B <- red::rli(grid.data.i$category_B, boot = boot.log, runs = boot)
    rli.i.C <- red::rli(grid.data.i$category_C, boot = boot.log, runs = boot)
    rli.i.D <- red::rli(grid.data.i$category_D, boot = boot.log, runs = boot)
    rli.i.new <- red::rli(grid.data.i$cat.reg.clean.new, boot = boot.log, runs = boot)
    rli.i <- c(rli.i, rli.i.A, rli.i.B, rli.i.C, rli.i.D, rli.i.new)
  } else {
    # rli.i <- c("LowCL" = NA,  "Median" = NA, "UpCL" = NA)
    rli.i <- c("RLI" = NA, "RLI_catA" = NA, "RLI_catB" = NA, 
               "RLI_catC" = NA, "RLI_catD" = NA, "RLI_new" = NA)
  }
  
  grid.data.ii <- grid.data.i[grid.data.i$endemic %in% "endemic",]
  if (dim(grid.data.ii)[1] > 0) {
    rli.i.end <- red::rli(grid.data.ii$cat.reg.clean, boot = boot.log, runs = boot)
    rli.i.end.A <- red::rli(grid.data.ii$category_A, boot = boot.log, runs = boot)
    rli.i.end.B <- red::rli(grid.data.ii$category_B, boot = boot.log, runs = boot)
    rli.i.end.C <- red::rli(grid.data.ii$category_C, boot = boot.log, runs = boot)
    rli.i.end.D <- red::rli(grid.data.ii$category_D, boot = boot.log, runs = boot)
    rli.i.end.new <- red::rli(grid.data.ii$cat.reg.clean.new, boot = boot.log, runs = boot)
    rli.i.end <- c(rli.i.end, rli.i.end.A, rli.i.end.B, rli.i.end.C, rli.i.end.D, rli.i.end.new)
  } else {
    # rli.i.end <- c("LowCL" = NA,  "Median" = NA, "UpCL" = NA)
    rli.i.end <- c("RLI" = NA, "RLI_catA" = NA, "RLI_catB" = NA, 
                   "RLI_catC" = NA, "RLI_catD" = NA, "RLI_new" = NA)
    
  }
  res <- c(rli.i, rli.i.end)
  names(res) <- c("RLI", "RLI_catA", "RLI_catB", "RLI_catC", "RLI_catD", "RLI_new",
                  "RLI.end", "RLI_catA.end", "RLI_catB.end", "RLI_catC.end", "RLI_catD.end", "RLI_new.end")
  resultado[[i]] <- res
}
all <- do.call(rbind.data.frame, resultado)
names(all) <- names(resultado[[1]])
all$ID <- names(resultado)
head(all)
grid.result <- merge(grid.result, all, 
                     by= "ID", all.x = TRUE, sort = FALSE)
grid.result <- grid.result[order(grid.result$order),]

#########################################
#### Final edits of the grid summary ####
#########################################

## Is there an effect of cell area on N.total, RLI and RLI.end? Yes, for all variables!
toto = grid.result[grid.result$area_ha<21500,]
par(mfrow=c(1,1))
plot(log(toto$N.total) ~ toto$area_ha)
abline(lm(log(toto$N.total) ~ toto$area_ha))
plot(log(toto$S.total) ~ toto$area_ha)
abline(lm(log(toto$S.total) ~ toto$area_ha))
plot(toto$RLI ~ toto$area_ha)
abline(lm(toto$RLI ~ toto$area_ha))

### Defining the relationship between cell area and sampling intensity and correcting the diversity estimates ###
#function that perform the corrections
# dados = grid.result
# lim.area = 215000
# nome.N = "N.total"
# nome.est = "Median.RLI"

# correct.estimates = function(dados,lim.area,nome.N,nome.est,...) {
#   # data for cells < 250000
#   toto = dados[dados$area_ha<lim.area,]
#   N = log(toto[,nome.N])
#   ids = !is.na(N) & !is.infinite(N)
#   area = toto$area_ha[ids]
#   N = N[ids]
#   # fitting and comparing the models: N ~ area
#   mod.N = lm(N ~ area)
#   mod.N1 = lm(N ~ area + I(area^2))
#   aic = bbmle::AICtab(mod.N,mod.N1,sort=FALSE)
#   if(aic$dAIC[1]>log(8)) mod.final = mod.N1 else mod.final = mod.N
#   par(mfrow=c(1,3), las=1,mar=c(4,4,1,0.5),mgp=c(2,0.5,0),tcl=-.3,...)
#   plot(N ~ area)
#   curve(predict(mod.final, newdata= data.frame(area= x)),add=TRUE, lwd=2, col=2)
#   # calculating the corrected sampling effort per grid cell
#   if (aic$dAIC[1] > log(8)) { 
#     N.est = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*area + coef(mod.final)[3]*area^2),2)
#     N.max = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*216506 + coef(mod.final)[3]*216506^2),2)
#   } else {
#     N.est = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*area),2)
#     N.max = round(exp(coef(mod.final)[1] + coef(mod.final)[2]*216506),2)
#   }
#   prop = (N.max - N.est)/N.est
#   decay = try(nls(prop ~ a*exp(b*area), start=list(a=25,b=-0.0001)),TRUE)
#   if (class(decay) == "try-error") decay = try(nls(prop~a*exp(b*area),start=list(a=10,b=-0.0001)),TRUE)
#   #plot(prop ~ area);   curve(coef(decay)[1]*exp(coef(decay)[2]*x), lwd=2, col=2, add=T)
#   prop[prop < 0.0] = coef(decay)[1] * exp(coef(decay)[2] * area[prop < 0.0]) * 0.5
#   N.correct = toto[,nome.N][ids] + round(prop*toto[,nome.N][ids],0)
#   
#   # data for cells < 250000
#   toto = dados[dados$area_ha > 170000,]
#   ids = !is.na(toto[,nome.N]) & !is.na(toto[,nome.est]) & !is.infinite(log(toto[,nome.N])) & !is.infinite(log(toto[,nome.est]))
#   N1 = log(toto[,nome.N])[ids]
#   est = toto[,nome.est][ids]
#   w = toto$SampCover[ids]
#   # fitting and comparing the models: N ~ area
#   mod.est  = lm(est ~ N1, weights = w^3)
#   mod.est1 = lm(est ~ N1 + I(N1^2), weights = w^3)
#   aic = bbmle::AICtab(mod.est,mod.est1,sort=FALSE)
#   if (aic$dAIC[1]>log(8)) mod.final1 = mod.est1 else mod.final1 = mod.est
#   plot(est ~ N1)
#   curve(predict(mod.final1, newdata= data.frame(N1= x)),add=TRUE, lwd=2, col=2)
#   # calculating the corrected sampling effort per grid cell
#   if (aic$dAIC[1]>log(8)) { 
#     est.pred = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*N + coef(mod.final1)[3]*N^2),4)
#     est.est = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*log(N.correct) + coef(mod.final1)[3]*log(N.correct)^2),4)
#   } else {
#     est.pred = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*N),4)
#     est.est = round(exp(coef(mod.final1)[1] + coef(mod.final1)[2]*log(N.correct)),4)
#   }
#   prop1 = (est.est - est.pred)/est.est
#   decay = try(nls(prop1~a*exp(b*area),start=list(a=25,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=5,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=1,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=0.1,b=-0.0001)),TRUE)
#   if(class(decay) == "try-error") decay = try(nls(prop1~a*exp(b*area),start=list(a=-0.2,b=-0.00001)),TRUE)
#   # plot(prop1 ~ area);   curve(coef(decay)[1]*exp(coef(decay)[2]*x), lwd=2, col=2, add=T)
#   # prop1[prop1<0.0] = coef(decay)[1]*exp(coef(decay)[2]*area[prop1<0.0])*0.1
#   prop1[prop1>0.0] = coef(decay)[1]*exp(coef(decay)[2]*area[prop1>0.0])*0.1
#   
#   # plot(prop1 ~ area)
#   est.correct = dados[,nome.est]
#   est.correct[dados$area_ha<lim.area & !is.na(dados[,nome.N])] = est.correct[dados$area_ha<lim.area & !is.na(dados[,nome.N])] +
#     prop1 * est.correct[dados$area_ha<lim.area & !is.na(dados[,nome.N])]
#   plot(est.correct ~ dados[,nome.est], xlab="estimate",ylab="corrected est.")
#   legend("topleft", nome.est, cex=1.2, bty="n");abline(0,1,lty=2)
#   par(mfrow=c(1,1))
#   return(est.correct)
# }
# 
# #generating the corrected estimates for each group
# grid.result$Median.RLI.correct = correct.estimates(grid.result,215000,"N.total","Median.RLI")
# grid.result$Median.RLI.end.correct = correct.estimates(grid.result,215000,"N.total","Median.RLI.end")

## Saving the grid results
saveRDS(list(grid.result, samp.cov.cutoff), "data/grid.results_50km.rds")

