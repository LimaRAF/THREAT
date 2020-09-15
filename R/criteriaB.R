##########################################################
##########################################################
#### ASSESSING SPECIES GEOGRAPHIC RANGE - CRITERION B ####
##########################################################
##########################################################
rm(list=ls())

### AFAZERES ###
# LER O GUI DO CNCFLORA PARA AS AVALIA??ES: QUAIS COORDENADAS USAR??


#### LOADING PACKAGES ###
devtools::install_github("gdauby/ConR", ref = "master", force = TRUE)
#devtools::install_github("gdauby/ConR@devel", force = TRUE)
detach("package:ConR", unload=TRUE)
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
                 repos = NULL, 
                 type = "source")
library("ConR")
require(data.table)
require(raster)
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

## GILLES COMMENT: Are you sure you want to crop the EOO (exclude area)? IUCN guidelines suggest to not do it for criterion B
## RENATO UPDATE: EOO is now computed without excluding the sea areas of EOO

system.time(
  EOO.hull <- EOO.computing(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                            method = "convex.hull",
                            method.less.than3 = "not comp",
                            export_shp = TRUE,
                            exclude.area = FALSE, country_map = NULL, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
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
saveRDS(shps_df_sf, "data/spp.convex.hull.polys_sf_uncropped.rds")
rm(EOO.hull, EOO, shps, shps_df, shps_df_sf); gc()

## Convex Hull method for each class of confidence level 
# took 8.25 and 73.3 minutes in my machine with the very ligth and with the ok-resolution Neotropical maps, respectively!
system.time(
  EOO.hull.sensitivity <- EOO.sensitivity (MyData[, c(1:3,9)],
                                           levels.order = c("low", "medium", "high"),
                                           occ.based = FALSE,
                                           exclude.area = FALSE,
                                           country_map = NULL,
                                           method.range = "convex.hull",
                                           write_results=FALSE,
                                           parallel = TRUE,
                                           NbeCores = 7,
                                           show_progress = TRUE)
)
saveRDS(EOO.hull.sensitivity, "data/EOO.convex.hull_uncropped.rds")
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
                       factor.div = 10, quant.max = 0.9)
radius1 = subpop.radius(MyData[, c(1:3)], 
                        factor.div = 10, quant.max = 0.9)

#Getting family-specific radius for missing estimates (only high confidence level)
tmp <- merge(radius, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp <- tmp[order(tmp$tax),]
fam.radius <- aggregate(as.double(tmp$radius), list(tmp$family), median, na.rm = TRUE)
radius.new <- merge(tmp, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new$radius[is.na(radius.new$radius)] <- 
  radius.new$x[is.na(radius.new$radius)]
radius.new$radius[as.double(radius.new$radius) < 4] <- 4
radius.new <- radius.new[order(radius.new$tax), ]
table(radius$tax == radius$tax)
hist(as.double(radius.new$radius), nclass = 40)

#Getting family-specific radius for missing estimates (any confidence level)
tmp1 <- merge(radius1, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp1 <- tmp1[order(tmp1$tax),]
fam.radius <- aggregate(as.double(tmp1$radius), list(tmp1$family), median, na.rm = TRUE)
radius.new1 <- merge(tmp1, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new1$radius[is.na(radius.new1$radius)] <- 
  radius.new1$x[is.na(radius.new1$radius)]
radius.new1$radius[as.double(radius.new1$radius) < 4] <- 4
radius.new1 <- radius.new1[order(radius.new1$tax), ]
table(radius1$tax == radius1$tax)
hist(as.double(radius.new1$radius), nclass = 40)

##Saving the information on subpopulation radius
tmp2 <- merge(radius.new1[,-4], radius.new[,-c(1,4)], by="tax", all.x= TRUE, sort = FALSE)
names(tmp2)[3:4] <- c("est.radius.level.1", "est.radius.level.2")
tmp2 <- tmp2[,c("family","tax","est.radius.level.1", "est.radius.level.2")]
tmp2 <- tmp2[order(tmp2$tax),]
saveRDS(tmp2, "data/est.radius.subpops.rds")

# tmp2 <- readRDS("data/est.radius.subpops.rds")
# radius.new <- tmp2[,c("tax", "est.radius.level.2")]
# radius.new1 <- tmp2[,c("tax", "est.radius.level.1")]
# names(radius.new) <- names(radius.new1) <- c("tax", "radius")
# radius.new$radius <- as.double(radius.new$radius)
# radius.new1$radius <- as.double(radius.new1$radius)

#Getting the subpopulations
SUB <- subpop.comp(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                   Resol_sub_pop= radius.new[,c("tax", "radius")],
                   parallel = TRUE, NbeCores = 5)
SUB1 <- subpop.comp(MyData[, c(1:3)],
                    Resol_sub_pop= radius.new1[,c("tax", "radius")],
                    parallel = TRUE, NbeCores = 5)

# #extracting the numb. of subpopulations from the object (high confidence level)
# tmp <- sapply(SUB, function(x) x[grepl("Number of subpopulation", names(x))])
# sp_names <- gsub(".Number of subpopulation","",names(tmp))
# n.sub.pop <- do.call(rbind.data.frame, tmp)
# dimnames(n.sub.pop) <- list(sp_names, "Number.subpops")
# 
# #### CHECK HERE: EXTRACTION OF POLYGONS NOT WORKING (IS IT NECESSARY??) ####
# tmp <- sapply(SUB, function(x) x[grepl("subpop.poly", names(x))])
# sp_names <- gsub(".subpop.poly","",names(tmp))
# for (i in 1:length(tmp))
#   slot(slot(tmp[[i]], "polygons")[[1]], "ID") <- sp_names[i]
# shps <- do.call(rbind, tmp)
# shps_df <- SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
# shps_df$tax <- as.character(shps_df$tax) 
# 
# #extracting the numb. of subpopulations from the object (any confidence level)
# tmp <- sapply(SUB1, function(x) x[grepl("Number of subpopulation", names(x))])
# sp_names <- gsub(".Number of subpopulation","",names(tmp))
# n.sub.pop1 <- do.call(rbind.data.frame, tmp)
# dimnames(n.sub.pop1) <- list(sp_names, "Number.subpops")

#Merging and saving
tmp1 <- merge(SUB1, SUB, by="row.names", all = TRUE, sort = FALSE)
names(tmp1) <- c("tax","Number.subpops.level.1", "Number.subpops.level.2")
tmp1 <- tmp1[order(tmp1$tax),]
saveRDS(tmp1, "data/number.subpops.rds")

### FIGURE SV ###
#Plotting the estimated radius against species mode of dispersal and seed mass
hab <- read.csv("data/threat_habitats.csv", as.is = TRUE)
tmp3 <- merge(radius1, radius, by="tax", all = TRUE, sort = FALSE)
tmp3 <- tmp3[order(tmp3$tax),]
tmp3 <- tmp3[!is.na(tmp3$radius.y),]
PopData <- merge(tmp3, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE)
PopData <- merge(PopData, resultado[,c("tax", "total.occs", "non.dup.occs")], by= "tax", all.x = TRUE)
PopData <- PopData[!is.na(PopData$taxon_id),]
PopData <- PopData[order(PopData$tax),]

sumario  <- summary(as.double(PopData$radius.y))
sumario1 <- summary(as.double(PopData$radius.y)[PopData$non.dup.occs>30])

jpeg(filename = "figures/Figure_SV.jpg", width = 3000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(1,2))
par(mar=c(4,3,0.75,0.5), mgp=c(2,0.35,0),tcl=-0.2,las=1)
plot(log10(as.double(radius.y)) ~ log10(non.dup.occs),
     data = PopData, #[PopData$non.dup.occs>30,],
     main = "", ylab = "Estimated radius (km)", xlab = "Number of occurrences",
     xaxt = "n", yaxt = "n")
axis(1, at=log10(c(5,15,30,50,100,1000)), 
     labels = c(5,15,30,50,100,1000))
axis(2, at=log10(c(1,5,50,200,500)), 
     labels = c(1,5,50,200,500))

par(mar=c(4,3,0.75,0.5), mgp=c(2,0.35,0),tcl=-0.2,las=1)
hist(as.double(PopData$radius.y), nclass=40,
     main = "", xlab = "Estimated radius (km)")
hist(as.double(PopData$radius.y)[PopData$non.dup.occs>30], nclass=40,
     add=TRUE, col="grey")
legend("topright", c(paste0("All species: ",
                     paste(round(sumario[c(2,3,5)],1), collapse = "-")," km"),
                     paste0("Only >30 obs.: ",
                     paste(round(sumario1[c(2,3,5)],1), collapse = "-")," km")),
       fill=c("white","grey"),  bty= "n")

# legend("topleft",legend=do.call(expression, list(
#   bquote({Marginal~italic(R)^2}~"="~.(round(attributes(sumario)$rsqs[[1]]*100,1))~"%"),
#   bquote({Conditional~italic(R)^2}~"="~.(round(attributes(sumario)$rsqs[[2]]*100,1))~"%"),
#   bquote("Number of obs. = "~.(attributes(sumario)$n)),
#   bquote(italic(chi)^2~"="~.(round(anova.best$Chisq[2],1)))
#   #bquote(italic(F)~"- test="~.(round(summary(toto)$fstatistic[1],2)))
# )),
# bty="n",horiz=F,cex=1.1,adj=c(0.1,0)
# )

## Not shown
# plot(log(as.double(radius.y)) ~ log(SeedMass_g),
#      data = PopData[PopData$non.dup.occs>30,])
# abline(lm(log(as.double(radius.y)) ~ log(SeedMass_g),
#           data = PopData[PopData$non.dup.occs>30,]),
#        lwd =2, col=2)
# 
# boxplot(log(as.double(radius.y)) ~ dispersal.syndrome,
#         data = PopData[PopData$dispersal.syndrome %in% c("anemochoric", "autochoric", "zoochoric"),],
#         notch =TRUE, varwidth = TRUE)
# anova(lm(log(as.double(radius.y)) ~ dispersal.syndrome,
#    data = PopData[PopData$dispersal.syndrome %in% c("anemochoric", "autochoric", "zoochoric"),]))

dev.off()
rm(tmp,tmp1,tmp2,radius,radius1,radius.new,radius.new1,SUB,SUB1,n.sub.pop,n.sub.pop1); gc()

#################################
#### AREA OF OCCUPANCY (AOO) ####
#################################
##Setting raster options
raster::removeTmpFiles(0.5)
raster::rasterOptions(tmpdir = "E://raster_TEMP//", 
              tmptime = 1.5, 
              #timer = TRUE,
              #tolerance = 0.5,
              chunksize = 1e+08,
              maxmemory = 1e+09)

#only high confidence level
system.time(
  AOO <- AOO.computing(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3)], 
                       Cell_size_AOO = 2, # Size of the grid cell in kilometers
                       nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                       parallel = TRUE, NbeCores = 5, 
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
# geom.col <- strict.ucs[ grepl("GEOMETRYCOLLECTION", types), ]
# polys.sp1<- list("vector", length(geom.col))
# for (i in 1:length(geom.col)){
#   tmp <- sf::st_cast(geom.col[i,], "GEOMETRY")
#   tmp1 <- sf::st_cast(tmp)
#   types <- vapply(sf::st_geometry(tmp), function(x) { class(x)[2]}, "")
#   tmp2 <- tmp1[ grepl("*POLYGON", types), ]
#   tmp3 <- sf::st_cast(tmp2, "MULTIPOLYGON")
#   polys.sp1[[i]] <- tmp3
# }
# polys.sp1 <- do.call(rbind, polys.sp1)
# polys.sp1 <- as(polys.sp1, "Spatial")
# strict.ucs.spdf <- rbind(polys.sp, polys.sp1)
strict.ucs.spdf <- polys.sp
rm(polys, polys.sp, polys.sp1,geom.col,strict.ucs)

system.time(
  locs <- my.locations.comp(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3),], 
                            method = "fixed_grid",
                            nbe_rep = 30, # number of raster with random starting position
                            Cell_size_locations = 10, #grid size in kilometres used for estimating the number of location
                            Rel_cell_size = 0.05, 
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
                             nbe_rep = 30, # number of raster with random starting position
                             Cell_size_locations = 10, #grid size in kilometres used for estimating the number of location
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


##############################
#### SEVERE FRAGMENTATION ####
##############################
detach("package:ConR", unload=TRUE)
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
                 repos = NULL, 
                 type = "source")
library("ConR")
?get.patches

##Setting raster options
raster::removeTmpFiles(0.5)
raster::rasterOptions(tmpdir = "E://raster_TEMP//", 
                      tmptime = 1.5, 
                      #timer = TRUE,
                      #tolerance = 0.5,
                      chunksize = 1e+08,
                      maxmemory = 1e+09)

##Getting the population radius for each species and merging
tmp2 <- readRDS("data/est.radius.subpops.rds")
radius.new <- tmp2[,c("tax", "est.radius.level.2")]
radius.new1 <- tmp2[,c("tax", "est.radius.level.1")]
names(radius.new) <- names(radius.new1) <- c("tax", "radius")
radius.new$radius <- as.double(radius.new$radius)
radius.new1$radius <- as.double(radius.new1$radius)
toto <- MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3),]
toto <- merge(toto, radius.new, by = "tax", all.x = TRUE, sort = FALSE)
toto <- toto[, c("ddlat","ddlon","tax","radius")]
toto1 <- MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon), c(1:3),]
toto1 <- merge(toto1, radius.new1, by = "tax", all.x = TRUE, sort = FALSE)
toto1 <- toto1[, c("ddlat","ddlon","tax","radius")]


##Getting the list of species data
list_data <- ConR:::coord.check(XY = toto,
                     listing = TRUE)
list_data1 <- ConR:::coord.check(XY = toto1,
                                listing = TRUE)

cl <- snow::makeSOCKcluster(5)
doSNOW::registerDoSNOW(cl)
`%d%` <- foreach::`%dopar%`
output <-
  foreach::foreach(
    x = 1:length(list_data),
    .combine = 'c',
    .options.snow = NULL
  ) %d% {
    res <- ConR:::get.patches(list_data[[x]],
                       cell_size = 2,
                       nbe_rep = 0,
                       proj_type = "cea",
                       Resol_sub_pop = as.double(unique(list_data[[x]]$radius)),
                       subpop_poly = NULL,
                       dist_isolated = as.double(unique(list_data[[x]]$radius)) * 1
                       #dist_isolated = 100
    )    
    res
  }

snow::stopCluster(cl)

ConR:::get.patches(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3),],
            cell_size = 2,
            nbe_rep = 0,
            proj_type = "cea",
            Resol_sub_pop = 10,
            subpop_poly = NULL,
            dist_isolated
)

  
#################################################################################################################################################################################################H    
#################################################################################################################################################################################################H    
##############################################################  
#### MERGING ALL RESULTS AND SAVING CRITERIA B PARAMETERS ####
##############################################################

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
#Loading packages
devtools::install_github("gdauby/ConR", ref = "devel", force = TRUE)
library(ConR)
library(dplyr)

#Getting the estimates for criterion B
critB_high <- readRDS("data/criteriaB_metrics_high_confidence.rds")
critB_high$tax <- as.character(critB_high$tax)
critB_low <- readRDS("data/criteriaB_metrics_low_confidence.rds")
critB_low$tax <- as.character(critB_low$tax)

#Getting estimates of species continuing decline
#GET MISSING DECLINE FROM AOO DECLINE (Criteria A)
critC <- readRDS("data/criterionC_all_prop_mature.rds")
est.decline <- sapply(strsplit(critC$cont.decline,"\\|"), tail, 1)  
est.decline <- gsub("\\(|\\)|[0-9]", "", est.decline)
est.decline <- gsub(" -", "", est.decline)
table(est.decline, critC$any.decline, useNA = "always")
critC$decline <- critC$any.decline
critC$decline[!critC$decline %in% "Decreasing" & est.decline %in% "Decreasing"] <- "Decreasing" 
  
critB_high <- merge(critB_high, critC[,c("species","decline")],
                    by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_high <- critB_high[order(critB_high$tax),]
critB_low <- merge(critB_low, critC[,c("species","decline")],
                    by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_low <- critB_low[order(critB_low$tax),]

##For assuming, that all other species (rarer) are decreasing as well
critB_high$decline[is.na(critB_high$decline)] <- "Decreasing"
critB_low$decline[is.na(critB_low$decline)] <- "Decreasing"

#Combining the info on number of localities and % in PAs
critB_high <- 
  critB_high %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)
critB_low <- 
  critB_low %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)

#Perfoming the assessments
results_Cb_high <- cat_criterion_b(EOO = critB_high$EOO,
                              AOO = critB_high$AOO,
                              locations = critB_high$nbe_loc_total,
                              #protected = critB_high$protected, 
                              decline = critB_high$decline,
                              protected.threshold = 100 ## Gilles, maybe set a minimun number of occurrences as well? Most of 100%PA species had <3 occurrences (90% of the time)
)
1 - table(results_Cb_high$ranks_B12a, useNA = "always")/dim(critB_high)[1] #20.75%

results_Cb_low <- cat_criterion_b(EOO = critB_low$EOO,
                                   AOO = critB_low$AOO,
                                   locations = critB_low$nbe_loc_total,
                                   #protected = critB_low$protected, 
                                   decline = critB_low$decline,
                                   protected.threshold = 100 ## Gilles, maybe set a minimun number of occurrences as well? Most of 100%PA species had <3 occurrences (90% of the time)
)
1 - table(results_Cb_low$ranks_B12a, useNA = "always")/dim(critB_low)[1] #20.75%

#Saving the results
results_Cb_high <- do.call(cbind.data.frame, c(results_Cb_high, stringsAsFactors = FALSE))
critB_high.all <- critB_high[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","decline")]
critB_high.all$category_B <- critB_high.all$category_B_code <- NA_character_
critB_high.all[!is.na(critB_high.all$AOO), c("category_B", "category_B_code")] <-
  results_Cb_high
critB_high.all[is.na(critB_high.all$AOO),]

# critB_high <- # was getting an error because of the ~10 species without classification
#   critB_high %>% 
#   mutate(category_B = results_Cb_high$ranks_B12a) %>%
#   mutate(category_B_code = results_Cb_high$cat_codes)
# critB_high %>% 
#   group_by(cb_category) %>% 
#   summarise(count = n())

results_Cb_low <- do.call(cbind.data.frame, c(results_Cb_low, stringsAsFactors = FALSE))
critB_low.all <- critB_low[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","decline")]
critB_low.all$category_B <- critB_low.all$category_B_code <- NA_character_
critB_low.all[!is.na(critB_low.all$AOO), c("category_B", "category_B_code")] <-
  results_Cb_low
critB_low.all[critB_low.all$tax %in% critB_high.all$tax[is.na(critB_high.all$AOO)],]

saveRDS(as.data.frame(critB_high.all), "data/criterionB_all_high_tax_confidence.rds")
saveRDS(as.data.frame(critB_low.all), "data/criterionB_all_low_tax_confidence.rds")

#### FIGURE: LOW VS. HIGH CONFIDENCE LEVEL ####
#### CHECK HERE: this figure should compare the optimum scheme with the high tax right? ####

jpeg(filename = "figures/Figure_SU.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

mat <- as.matrix(table(paste0(critB_high.all$category_B,"_hi"), paste0(critB_low.all$category_B,"_lo")))
mat <- mat[c(1,2,5,3), c(3,4,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
grid.col = c(CR_lo = "red", EN_lo = "darkorange", VU_lo = "gold", LCorNT_lo = "khaki",
             CR_hi = "red", EN_hi = "darkorange", VU_hi = "gold", LCorNT_hi = "khaki")
col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 15] = mat[mat < 15]*2
mat[mat > 0 & mat < 5] = 10

#plotting the diagram
par(mfrow=c(1,1))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible) = FALSE
chordDiagram(mat, big.gap = 20, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
             grid.col = grid.col, col = col_mat,
             self.link = 1, link.visible = visible,
             #h=0.9,
             #w=1,
             #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
             link.lwd = 4
             #h.ratio = 0.7
             #reduce_to_mid_line = FALSE,
             #w2=0.5,
             #rou=0.2
             #point1 = rep(0,16)
)
#Putting legends on
sec.ind <- c("CR","EN","VU","LC or NT","LC or NT","VU","EN","CR")
for(si in get.all.sector.index()) {
  lab <- sec.ind[which(si == get.all.sector.index())]
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  if(si == "VU_lo") {
    circos.text(mean(xlim), mean(ylim), labels = "VU", 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = FALSE, col = "black")
    
  } else {
    circos.text(mean(xlim), mean(ylim), labels = lab, 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = TRUE, col = "black")
  }  
}
legend("topleft","High confidence", bty="n", cex=1.2)
legend("topright","Any confidence", bty="n", cex=1.2)
dev.off()


#### GILLES EXTRA CODES ####

library(tidyverse)
library(ConR)

# all_source <- 
#   list.files("D:/MonDossierR/ConR_orig/R/", full.names = T)
# for (i in 1:length(all_source)) source(all_source[i])


## list of species with at least one georeferenced occurrence
MyData %>% 
  filter(!is.na(ddlat), !is.na(ddlon)) %>% 
  distinct(tax)

## list of species with no georeferenced occurrences
MyData %>% 
  filter(is.na(ddlat) | is.na(ddlon)) %>% 
  distinct(tax)

## no paralleling 2224.4 secondes
system.time(results_cb <- 
              criterion_B(x = MyData))
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

