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
# devtools::install_github("gdauby/ConR@devel", force = TRUE)
# detach("package:ConR", unload=TRUE)
# install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
#                  repos = NULL, 
#                  type = "source")
library("ConR")
require(data.table)
require(raster)
require(stars)
require(sf)

source("./R/suggestions_for_ConR.r")
source("C://Users//renato//Documents//raflima//R_packages//Backups//ConR//R1//EOO.sensitivity.R")
source("C://Users//renato//Documents//raflima//R_packages//Backups//ConR//R1//over.valid.poly.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//subpop.comp.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//proj_crs.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//coord.check.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//subpop.estimation.R")
# source("C://Users//renato//Documents//raflima//R_packages//ConR//R//AOO.estimation.R")

#### LOADING THE NEOTROPICS MAP ###
neotrop.simp <- readRDS("data/Contour_Neotrop_simplified_tol_005_no_small.rds")

#### LOADING THREAT OCCURRENCE DATA (HERBARIUM + TREECO) ###
oc.data <- readRDS("data/threat_occ_data_final.rds")

#Putting data in the ConR format
MyData <- oc.data[, c("ddlat","ddlon",
                      "tax","higher.tax.rank",
                      "coly","vouchers",
                      "detBy","dety",
                      "tax.check2","tax.check.final","UC",
                      "dist.eoo","tax.conf","source")]
MyData$tax.check2 <- MyData$tax.check2 %in% "TRUE" # the super,hyper high confidence level
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
## RENATO UPDATE 2: EOO is now also computed excluding the sea areas of EOO 
#for obtaining the protected areas, forest cover and HII within species EOO

system.time(
  EOO.hull <- ConR:::EOO.computing(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                            method = "convex.hull",
                            method.less.than3 = "not comp",
                            export_shp = TRUE,
                            # exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                            exclude.area = FALSE, country_map = NULL, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                            write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                            write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                            parallel = TRUE, NbeCores = 6) # run on parallel? How many cores?
) #Took 40 secs! 26 minutes if we exclude the sea areas

#extracting the EOO from the object
EOO <- do.call(rbind.data.frame, EOO.hull[grepl("EOO", names(EOO.hull))])
sp_names <- unique(MyData$tax[grepl("high", MyData$tax.check.final)])
dimnames(EOO) <- list(sp_names, "EOO")
# EOO <- EOO.hull[[1]]
# EOO1 <- EOO.hull1[[1]]

shps <- EOO.hull[grepl("spatial.polygon", names(EOO.hull))]
for (i in 1:length(shps))
  slot(slot(shps[[i]], "polygons")[[1]], "ID") <- sp_names[!is.na(EOO$EOO)][i]
shps <- do.call(rbind, shps)
shps_df <- SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
#shps_df <-EOO.hull[[2]]
shps_df$tax <- as.character(shps_df$tax) 

#inpecting
par(mar=c(1,0,1,0))
sp <- "Abarema limae" # restrictd range
sp <- "Anadenanthera peregrina" # more widespread
# plot(sf::st_geometry(shps_df[shps_df$taxa %in% sp,1]))
plot(shps_df[sp,])
points(MyData[MyData$tax %in% sp, 2:1], pch=21)
points(MyData[grepl("high", MyData$tax.check.final) & MyData$tax %in% sp, 2:1],
       pch = 19)
plot(neotrop.simp, add=TRUE, border = "grey")

#saving species EOO (convex hulls)
shps_df_sf <- sf::st_as_sf(shps_df)
# saveRDS(shps_df, "data/spp.convex.hull.polys.rds")
# saveRDS(shps_df_sf, "data/spp.convex.hull.polys_sf_cropped.rds")
saveRDS(shps_df_sf, "data/spp.convex.hull.polys_sf_uncropped.rds")
rm(EOO.hull, EOO, shps, shps_df, shps_df_sf); gc()

## Convex Hull method for each class of confidence level 
# took 8.25 and 73.3 minutes in my machine with the very ligth and with the ok-resolution Neotropical maps, respectively!
# took 96 secs without the map (7 cores)
system.time(
  EOO.hull.sensitivity <- EOO.sensitivity (MyData[, c(1:3,10)],
                                           levels.order = c("low", "medium", "high"), # high here is the optimum scheme described in the methods
                                           occ.based = FALSE,
                                           exclude.area = FALSE,
                                           country_map = NULL,
                                           method.range = "convex.hull",
                                           write_results=FALSE,
                                           parallel = TRUE,
                                           NbeCores = 7,
                                           show_progress = TRUE)
)
EOO.hull.sensitivity1 <- EOO.sensitivity (MyData[, c(1:3,9)],
                                                  levels.order = c(FALSE, TRUE),
                                                  occ.based = FALSE,
                                                  exclude.area = FALSE,
                                                  country_map = NULL,
                                                  method.range = "convex.hull",
                                                  write_results=FALSE,
                                                  parallel = TRUE,
                                                  NbeCores = 7,
                                                  show_progress = TRUE)

## Merging the info (all occurrences, optimum, and high) and re-organizing
EOO.hull.sensitivity1 <- EOO.hull.sensitivity1[,c("Species","Occs.level.2","EOO.level.2","EOO.increase.1")]
names(EOO.hull.sensitivity1) <- c("Species","Occs.level.4","EOO.level.4","EOO.increase.1.1")
EOO.hull.sensitivity <- merge(EOO.hull.sensitivity, EOO.hull.sensitivity1, by="Species", all = TRUE, sort = FALSE)
EOO.hull.sensitivity <- EOO.hull.sensitivity[,c("Species","Occs.level.1","Occs.level.2","Occs.level.3","Occs.level.4",
                                                "EOO.level.1","EOO.level.2","EOO.level.3","EOO.level.4",
                                                "EOO.increase.1","EOO.increase.2","EOO.increase.1.1")]
EOO.hull.sensitivity <- EOO.hull.sensitivity[order(EOO.hull.sensitivity$Species),]
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
detach("package:ConR", unload=TRUE)
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local
                 repos = NULL,
                 type = "source")
library("ConR")


#This function applies the method called circular buffer method for estimating the number of subpopulation (Rivers et al., 2010).
#From Rivers et al. 2010: The ideal radius of the buffer is debatable; however when dispersal
#characteristics of the species are unknown then a sliding scale, such as the 1/10th maximum inter-point distance, is the preferred choice, as it is species-specific and not sensitive to collection density
radius <- subpop.radius(MyData[grepl("high", MyData$tax.check.final), c(1:3)], 
                       factor.div = 10, quant.max = 0.9)
radius1 <- subpop.radius(MyData[, c(1:3)], 
                        factor.div = 10, quant.max = 0.9)
radius2 <- subpop.radius(MyData[MyData$tax.check2, c(1:3)], 
                        factor.div = 10, quant.max = 0.9)


#Getting family-specific radius for missing estimates (only optimum confidence level) 
tmp <- merge(radius, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp <- tmp[order(tmp$tax),]
tmp$genus <- as.character(sapply(strsplit(tmp$tax, " "), function(x) x[1]))
gen.radius <- aggregate(as.double(tmp$radius), list(tmp$genus), median, na.rm = TRUE)
fam.radius <- aggregate(as.double(tmp$radius), list(tmp$family), median, na.rm = TRUE)
radius.new <- merge(tmp, gen.radius, by.x = 'genus', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new <- merge(radius.new, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new$radius.final <- as.double(radius.new$radius)
radius.new$radius.final[is.na(radius.new$radius.final)] <- 
  radius.new$x.x[is.na(radius.new$radius.final)]
radius.new$radius.final[is.na(radius.new$radius.final)] <- 
  radius.new$x.y[is.na(radius.new$radius.final)]
radius.new$radius.final[radius.new$radius.final < 4] <- 4
radius.new <- radius.new[order(radius.new$tax), ]
table(radius$tax == tmp$tax)
#hist(as.double(radius.new$radius.final), nclass = 40)

#Getting family-specific radius for missing estimates (any confidence level)
tmp1 <- merge(radius1, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp1 <- tmp1[order(tmp1$tax),]
tmp1$genus <- as.character(sapply(strsplit(tmp1$tax, " "), function(x) x[1]))
gen.radius <- aggregate(as.double(tmp1$radius), list(tmp1$genus), median, na.rm = TRUE)
fam.radius <- aggregate(as.double(tmp1$radius), list(tmp1$family), median, na.rm = TRUE)
radius.new1 <- merge(tmp1, gen.radius, by.x = 'genus', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new1 <- merge(radius.new1, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new1$radius.final <- as.double(radius.new1$radius)
radius.new1$radius.final[is.na(radius.new1$radius.final)] <- 
  radius.new1$x.x[is.na(radius.new1$radius.final)]
radius.new1$radius.final[is.na(radius.new1$radius.final)] <- 
  radius.new1$x.y[is.na(radius.new1$radius.final)]
radius.new1$radius.final[radius.new1$radius.final < 4] <- 4
radius.new1 <- radius.new1[order(radius.new1$tax), ]
table(radius1$tax == tmp1$tax)
#hist(as.double(radius.new1$radius.final), nclass = 40)

#Getting family-specific radius for missing estimates (only high confidence level)
tmp2 <- merge(radius2, resultado[,c("family","tax")], by = 'tax', all.x = TRUE, sort=FALSE)
tmp2 <- tmp2[order(tmp2$tax),]
tmp2$genus <- as.character(sapply(strsplit(tmp2$tax, " "), function(x) x[1]))
gen.radius <- aggregate(as.double(tmp2$radius), list(tmp2$genus), median, na.rm = TRUE)
fam.radius <- aggregate(as.double(tmp2$radius), list(tmp2$family), median, na.rm = TRUE)
radius.new2 <- merge(tmp2, gen.radius, by.x = 'genus', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new2 <- merge(radius.new2, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new2$radius.final <- as.double(radius.new2$radius)
radius.new2$radius.final[is.na(radius.new2$radius.final)] <- 
  radius.new2$x.x[is.na(radius.new2$radius.final)]
radius.new2$radius.final[is.na(radius.new2$radius.final)] <- 
  radius.new2$x.y[is.na(radius.new2$radius.final)]
radius.new2$radius.final[is.na(radius.new2$radius.final)] <- 
  median(radius.new2$x.y, na.rm = TRUE)
radius.new2$radius.final[radius.new2$radius.final < 4] <- 4
radius.new2 <- radius.new2[order(radius.new2$tax), ]
table(radius2$tax == tmp2$tax)
#hist(as.double(radius.new2$radius.final), nclass = 40)


##Saving the information on subpopulation radius
tmp2 <- merge(radius.new1[,c(3,4)], radius.new[,c(3,4)], by="tax", all.x= TRUE, sort = FALSE)
tmp2 <- merge(tmp2, radius.new2[,c(3,4)], by="tax", all.x= TRUE, sort = FALSE)
names(tmp2) <- c("tax","est.radius.level.1", "est.radius.level.2", "est.radius.level.3")
tmp2 <- tmp2[order(tmp2$tax),]
saveRDS(tmp2, "data/est.radius.subpops.rds")

tmp2 <- merge(radius.new1[,c(3,7)], radius.new[,c(3,7)], by="tax", all.x= TRUE, sort = FALSE)
tmp2 <- merge(tmp2, radius.new2[,c(3,7)], by="tax", all.x= TRUE, sort = FALSE)
names(tmp2) <- c("tax","final.radius.level.1", "final.radius.level.2", "final.radius.level.3")
tmp2 <- tmp2[order(tmp2$tax),]
saveRDS(tmp2, "data/final.radius.subpops.rds")

#Getting the subpopulations
radius.new <- radius.new[,c("tax", "radius.final")]
radius.new1 <- radius.new1[,c("tax", "radius.final")]
radius.new2 <- radius.new2[,c("tax", "radius.final")]
names(radius.new) <- names(radius.new1) <- names(radius.new2) <- 
  c("tax", "radius")

SUB <- subpop.comp(MyData[grepl("high", MyData$tax.check.final), c(1:3)], # optimum tax. confidence
                   Resol_sub_pop= radius.new[,c("tax", "radius")],
                   parallel = TRUE, NbeCores = 5)
SUB1 <- subpop.comp(MyData[, c(1:3)],
                    Resol_sub_pop= radius.new1[,c("tax", "radius")], # any confidence
                    parallel = TRUE, NbeCores = 5)
SUB2 <- subpop.comp(MyData[MyData$tax.check2, c(1:3)],
                    Resol_sub_pop= radius.new2[,c("tax", "radius")], # only high confidence
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
tmp1 <- merge(tmp1, SUB2, by.x="Row.names", by.y="row.names", all = TRUE, sort = FALSE)
names(tmp1) <- c("tax","Number.subpops.level.1", "Number.subpops.level.2", "Number.subpops.level.3")
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
PopData <- PopData[!is.na(PopData$internal_taxon_id),]
PopData <- PopData[order(PopData$tax),]

sumario  <- summary(as.double(PopData$radius.y))
sumario1 <- summary(as.double(PopData$radius.y)[PopData$non.dup.occs>30])

jpeg(filename = "figures/Figure_SV.jpg", width = 3000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(1,2))
par(mar=c(4,3,0.75,0.5), mgp=c(2,0.35,0),tcl=-0.2,las=1)
plot(log10(as.double(radius.y)) ~ log10(non.dup.occs),
     data = PopData[PopData$habito %in% "tree",], #[PopData$non.dup.occs>30,],
     main = "", ylab = "Estimated radius (km)", xlab = "Number of occurrences",
     xaxt = "n", yaxt = "n", ylim = log10(c(0.05,800)),
     col = adjustcolor("red", alpha.f = 0.4), pch = 19)
points(log10(as.double(radius.y)) ~ log10(non.dup.occs),
       data = PopData[PopData$habito %in% "shrub",],  
       col = adjustcolor("blue", alpha.f = 0.4), pch = 15)
axis(1, at=log10(c(5,15,30,50,100,1000)), 
     labels = c(5,15,30,50,100,1000))
axis(2, at=log10(c(1,5,50,200,500)), 
     labels = c(1,5,50,200,500))
legend("bottomright", c("Trees", "Shrubs"), 
       pch = c(19,15), col=c("red","blue"), cex=1.1, bty="n")


par(mar=c(4,3,0.75,0.5), mgp=c(2,0.35,0),tcl=-0.2,las=1)
hist(as.double(PopData$radius.y), nclass=40,
     main = "", xlab = "Estimated radius (km)", col = "white")
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
detach("package:ConR", unload=TRUE)
devtools::install_github("gdauby/ConR", ref = "master", force = TRUE)
library("ConR")

##Setting raster options
raster::removeTmpFiles(0.5)
raster::rasterOptions(tmpdir = "E://raster_TEMP//", 
              tmptime = 1.5, 
              #timer = TRUE,
              #tolerance = 0.5,
              #chunksize = 1e+08,
              maxmemory = 1e+09)

#only high confidence level
system.time( # Took using v1.3
  AOO <- AOO.computing(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3)], 
                       Cell_size_AOO = 2, # Size of the grid cell in kilometers
                       nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                       parallel = TRUE, NbeCores = 5, 
                       show_progress= TRUE, export_shp=FALSE)
)
#any confidence level
system.time(
  AOO1 <- AOO.computing(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon), c(1:3)], 
                       Cell_size_AOO = 2, # Size of the grid cell in kilometers
                       nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                       parallel = TRUE, NbeCores = 6, 
                       show_progress= TRUE, export_shp=FALSE)
)
system.time(
  AOO2 <- AOO.computing(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & MyData$tax.check2, c(1:3)], 
                       Cell_size_AOO = 2, # Size of the grid cell in kilometers
                       nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                       parallel = TRUE, NbeCores = 6, 
                       show_progress= TRUE, export_shp=FALSE)
)

#Saving
tmp2 <- merge(AOO1, AOO, by="row.names", all = TRUE, sort = FALSE)
tmp2 <- merge(tmp2, AOO2, by.x="Row.names", by.y="row.names", all = TRUE, sort = FALSE)
names(tmp2) <- c("tax","AOO.level.1", "AOO.level.2", "AOO.level.3")
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
saveRDS(cbind.data.frame(locs1[[3]],locs1[4]), "data/tmp1.rds")

system.time(
  locs2 <- my.locations.comp(MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & MyData$tax.check2, c(1:3),], 
                             method = "fixed_grid",
                             nbe_rep = 30, # number of raster with random starting position
                             Cell_size_locations = 10, #grid size in kilometres used for estimating the number of location
                             Rel_cell_size=0.05,
                             protec.areas = strict.ucs.spdf, #a SpatialPolygonsDataFrame, shapefile with protected areas.
                             ID_shape_PA= "NAME", # field name of protec.areas with ID
                             method_protected_area="no_more_than_one", 
                             parallel=TRUE, NbeCores=7)
)
table(names(locs2[[3]]) == names(locs2[[4]]))
saveRDS(cbind.data.frame(locs2[[3]],locs2[4]), "data/tmp2.rds")

##Combining the results
tmp <- cbind.data.frame(locs[[3]],locs[4])
names(tmp) <- c("LocUCs.level2", "LocOutUCs.level2")
tmp1 <- cbind.data.frame(locs1[[3]],locs1[4])
names(tmp1) <- c("LocUCs.level1", "LocOutUCs.level1")
tmp2 <- cbind.data.frame(locs2[[3]],locs2[4])
names(tmp2) <- c("LocUCs.level3", "LocOutUCs.level3")


localities <- merge(tmp1, tmp, by= "row.names", all = TRUE, sort = FALSE)
localities <- merge(localities, tmp2, by.x= "Row.names", by.y= "row.names", all = TRUE, sort = FALSE)
names(localities)[1] <- "tax"
localities$tax <- gsub("_", " ", localities$tax)
localities <- localities[order(localities$tax),]
saveRDS(localities, "data/number_localities.rds")


##############################
#### SEVERE FRAGMENTATION ####
##############################

## Setting a temporary file directory in the external HD
write("TMPDIR = 'E:/Rtemp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

#Loading packages
require(data.table)
require(raster)
require(stars)
require(sf)
detach("package:ConR", unload=TRUE)
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
                 repos = NULL, 
                 type = "source")
library("ConR")

##Setting raster options
raster::showTmpFiles()
raster::removeTmpFiles(h = 4)
raster:::tmpDir()
raster::rasterOptions(tmpdir = 'E://Rtemp//raster')
raster::rasterOptions(tmptime = 4,
                      chunksize = 1e+08,
                      maxmemory = 1e+09)

##Getting the population radius for each species
tmp2 <- readRDS("data/final.radius.subpops.rds")
radius.new <- tmp2[,c("tax", "final.radius.level.2")]
radius.new1 <- tmp2[,c("tax", "final.radius.level.1")]
radius.new2 <- tmp2[,c("tax", "final.radius.level.3")]
names(radius.new) <- names(radius.new1) <- names(radius.new2) <-c("tax", "radius")
radius.new$radius <- as.double(radius.new$radius)
radius.new1$radius <- as.double(radius.new1$radius)
radius.new2$radius <- as.double(radius.new2$radius)

##Getting the AOO for each species
AOO <- readRDS("data/AOO.rds")

##Creating the new objects for analysis
toto <- MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & grepl("high", MyData$tax.check.final), c(1:3),]
toto <- merge(toto, radius.new, by = "tax", all.x = TRUE, sort = FALSE)
toto <- merge(toto, AOO[,c("tax", "AOO.level.2")], by = "tax", all.x = TRUE, sort = FALSE)
toto <- toto[, c("ddlat","ddlon","tax","radius","AOO.level.2")]
names(toto)[5] <- "AOO"

toto1 <- MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon), c(1:3),]
toto1 <- merge(toto1, radius.new1, by = "tax", all.x = TRUE, sort = FALSE)
toto1 <- merge(toto1, AOO[,c("tax", "AOO.level.1")], by = "tax", all.x = TRUE, sort = FALSE)
toto1 <- toto1[, c("ddlat","ddlon","tax","radius","AOO.level.1")]
names(toto1)[5] <- "AOO"

toto2 <- MyData[!is.na(MyData$ddlat) & !is.na(MyData$ddlon) & MyData$tax.check2, c(1:3),]
toto2 <- merge(toto2, radius.new2, by = "tax", all.x = TRUE, sort = FALSE)
toto2 <- merge(toto2, AOO[,c("tax", "AOO.level.3")], by = "tax", all.x = TRUE, sort = FALSE)
toto2 <- toto2[, c("ddlat","ddlon","tax","radius","AOO.level.3")]
names(toto2)[5] <- "AOO"

dados <- list(toto, toto1, toto2)
rm(toto, toto1, toto2)
resultado <- vector("list", length(dados))
fator <- 1
dist <- 100

require(parallel)
require(doParallel)
for(i in 1:length(dados)) {

#i=1

  ##Getting the list of species data
  list_data <- ConR:::coord.check(XY = dados[[i]],
                                  listing = TRUE)

  cl <- makeCluster(detectCores()-4)
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(raster))
  clusterEvalQ(cl, library(ConR))
  clusterEvalQ(cl, library(stars))

  clusterExport(cl, list("list_data", "fator", "dist"),envir=environment())
  output <- parLapply(cl, 1:length(list_data), fun= function(x) {
  #output <- parLapply(cl, 1:5, fun= function(x) {
  
  # output <- vector("list", length(list_data)) 
  # for(x in 1:length(list_data)) {
    
    res <- try(ConR:::get.patches(
      list_data[[x]],
      cell_size = 2,
      nbe_rep = 0,
      AOO = as.double(unique(list_data[[x]]$AOO)),
      Resol_sub_pop = as.double(unique(list_data[[x]]$radius)),
      subpop_poly = NULL,
      #dist_isolated = as.double(unique(list_data[[x]]$radius)) * fator,
      dist_isolated = dist,
      proj_type = "cea",
      export_shp = FALSE
    ),TRUE)
      
      if(class(res) == "try-error") {
        res <- NULL
        res <- "ERRO"
        names(res) <- unique(list_data[[x]]$tax)
      }
  return(res)
  })
      
  #   cat(x,"\n")
  #   output[[x]] <-res 
  # }    

  stopImplicitCluster()
  stopCluster(cl)

  output1 <- do.call("c", output)
  df1 <- data.frame(tax = names(output1), frag = as.character(output1), stringsAsFactors = FALSE) 
  resultado[[i]] <- df1
  #saveRDS(df1, paste0("data/tmp",i,".rds"))
} 
resultados <- merge(resultado[[2]], resultado[[1]], by= "tax", all.x = TRUE, sort = FALSE)
resultados <- merge(resultados, resultado[[3]], by= "tax", all.x = TRUE, sort = FALSE)
names(resultados)[2:4] <- c("frag.level.1","frag.level.2","frag.level.3") 
resultados <- resultados[order(resultados$tax),]
saveRDS(resultados, paste0("data/severe_frag_dist_",dist,"km.rds"))
rm(dados)

##Assessing the impact of dist_isolated
dist25 <- readRDS("data/severe_frag_dist_25km.rds")[,c("tax", "frag.level.2")]
dist50 <- readRDS("data/severe_frag_dist_50km.rds")[,c("tax", "frag.level.2")]
dist75 <- readRDS("data/severe_frag_dist_75km.rds")[,c("tax", "frag.level.2")]
dist100 <- readRDS("data/severe_frag_dist_100km.rds")[,c("tax", "frag.level.2")]
dist_rad <- readRDS("data/severe_frag_dist_radius.rds")[,c("tax", "frag.level.2")]
AOO <- readRDS("data/AOO.rds")
subpops <- readRDS("data/number.subpops.rds")
dist0 <- merge(AOO[,c("tax", "AOO.level.2")], subpops[,c("tax", "Number.subpops.level.2")],
               by = "tax", all = TRUE, sort = FALSE) 
dist0$frag.level.2 <- 100 * dist0$Number.subpops.level.2/ (dist0$AOO.level.2 / 4)
  
jpeg(filename = "figures/Figure_TT.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

par(mar=c(3,3,1,1), mgp=c(1.5,0.25,0),tcl=-0.2,las=1)
plot(0:6, c(0.8,16,33, 50, 66, 83,100), 
     col= "white", xaxt = "n", cex.lab = 1.2,
     xlab = "Large dispersal distance", ylab="Fragmentation level (%)")#, log = "y")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5,5.5), labels = c("0 km", "25 km","50 km", "75 km", "100 km", "radius"))
boxplot(as.double(dist0$frag.level.2), add=TRUE, notch = TRUE, at=0.5, yaxt = "n")
boxplot(as.double(dist25$frag.level.2), add=TRUE, notch = TRUE, at=1.5, yaxt = "n")
boxplot(as.double(dist50$frag.level.2), add=TRUE, notch = TRUE, at=2.5, yaxt = "n")
boxplot(as.double(dist75$frag.level.2), add=TRUE, notch = TRUE, at=3.5, yaxt = "n")
boxplot(as.double(dist100$frag.level.2), add=TRUE, notch = TRUE, at=4.5, yaxt = "n")
boxplot(as.double(dist_rad$frag.level.2), add=TRUE, notch = TRUE, at=5.5, yaxt = "n")
abline(h = 50, lty =3)
legenda <- as.character(c(round((100*table(as.double(dist0$frag.level.2)>50)/dim(dist0)[1])[2],1),
             round((100*table(as.double(dist25$frag.level.2)>50)/dim(dist25)[1])[2],1),
             round((100*table(as.double(dist50$frag.level.2)>50)/dim(dist50)[1])[2],1),
             round((100*table(as.double(dist75$frag.level.2)>50)/dim(dist75)[1])[2],1),
             round((100*table(as.double(dist100$frag.level.2)>50)/dim(dist100)[1])[2],1),
             round((100*table(as.double(dist_rad$frag.level.2)>50)/dim(dist_rad)[1])[2],1)))
xpos <- c(0.5,1.5,2.5,3.5,4.5,5.5) + 0.35
for(i in 1:length(xpos)) text(xpos[i], 53.5, paste0(legenda[i], "%"))
dev.off()

#################################################################################################################################################################################################H    
#################################################################################################################################################################################################H    
##############################################################  
#### MERGING ALL RESULTS AND SAVING CRITERIA B PARAMETERS ####
##############################################################
rm(list=ls())

#Reading the saved metrics
EOO <- readRDS("data/EOO.convex.hull_uncropped.rds")
AOO <- readRDS("data/AOO.rds")
subpops <- readRDS("data/number.subpops.rds")
localities <- readRDS("data/number_localities.rds")
sev.frag <- merge(AOO, subpops,
               by = "tax", all = TRUE, sort = FALSE) 
sev.frag$frag.level.1 <- 100 * sev.frag$Number.subpops.level.1/ (sev.frag$AOO.level.1 / 4)
sev.frag$frag.level.2 <- 100 * sev.frag$Number.subpops.level.2/ (sev.frag$AOO.level.2 / 4)
sev.frag$frag.level.3 <- 100 * sev.frag$Number.subpops.level.3/ (sev.frag$AOO.level.3 / 4)

#Checking the species order
table(EOO$Species == AOO$tax)
table(EOO$Species == subpops$tax)
table(EOO$Species == localities$tax)
table(EOO$Species == sev.frag$tax)

#Building the data frame for each level of taxonomic confidence
critB_opt <- cbind.data.frame(EOO[, c("Species", "Occs.level.3", "EOO.level.3")],
                               AOO = as.double(AOO[,"AOO.level.2"]),
                               Nbe_subPop = as.double(subpops[,"Number.subpops.level.2"]),
                               Nbe_loc = as.double(localities[,"LocOutUCs.level2"]),
                               Nbe_loc_PA = as.double(localities[,"LocUCs.level2"]),
                               Frag.level = as.double(sev.frag[,"frag.level.2"]),
                               stringsAsFactors = FALSE)
critB_low <- cbind.data.frame(EOO[, c("Species", "Occs.level.1", "EOO.level.1")],
                               AOO = as.double(AOO[,"AOO.level.1"]),
                               Nbe_subPop = as.double(subpops[,"Number.subpops.level.1"]),
                               Nbe_loc = as.double(localities[,"LocOutUCs.level1"]),
                               Nbe_loc_PA = as.double(localities[,"LocUCs.level1"]),
                               Frag.level = as.double(sev.frag[,"frag.level.1"]),
                               stringsAsFactors = FALSE)
critB_high <- cbind.data.frame(EOO[, c("Species", "Occs.level.4", "EOO.level.4")],
                               AOO = as.double(AOO[,"AOO.level.3"]),
                               Nbe_subPop = as.double(subpops[,"Number.subpops.level.3"]),
                               Nbe_loc = as.double(localities[,"LocOutUCs.level3"]),
                               Nbe_loc_PA = as.double(localities[,"LocUCs.level3"]),
                               Frag.level = as.double(sev.frag[,"frag.level.3"]),
                               stringsAsFactors = FALSE)


names(critB_opt)[1:3] <- names(critB_low)[1:3] <- names(critB_high)[1:3] <- c("tax", "Nbe_occs", "EOO")
critB_opt <- critB_opt[,c("tax","EOO","AOO","Nbe_subPop","Nbe_loc","Nbe_loc_PA","Nbe_occs")]
critB_low <- critB_low[,c("tax","EOO","AOO","Nbe_subPop","Nbe_loc","Nbe_loc_PA","Nbe_occs")]
critB_high <- critB_high[,c("tax","EOO","AOO","Nbe_subPop","Nbe_loc","Nbe_loc_PA","Nbe_occs")]
saveRDS(critB_opt, "data/criteriaB_metrics_optim_confidence.rds")
saveRDS(critB_low, "data/criteriaB_metrics_low_confidence.rds")
saveRDS(critB_high, "data/criteriaB_metrics_high_confidence.rds")


#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA B ####
#############################
#Loading packages
detach("package:ConR", unload=TRUE)
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local
                 repos = NULL,
                 type = "source")
library("ConR")
library(dplyr)

#Getting the estimates for criterion B
critB_opt <- readRDS("data/criteriaB_metrics_optim_confidence.rds")
critB_opt$tax <- as.character(critB_opt$tax)
critB_low <- readRDS("data/criteriaB_metrics_low_confidence.rds")
critB_low$tax <- as.character(critB_low$tax)
critB_high <- readRDS("data/criteriaB_metrics_high_confidence.rds")
critB_high$tax <- as.character(critB_high$tax)

#Getting estimates of species continuing decline based on the index of abundance (Criterion C)
critC <- readRDS("data/criterionC_all_prop_mature.rds")
est.decline <- sapply(strsplit(critC$cont.decline,"\\|"), tail, 1)  
est.decline <- gsub("\\(|\\)|[0-9]", "", est.decline)
est.decline <- gsub(" -", "", est.decline)
table(est.decline, critC$any.decline, useNA = "always")
critC$declineC <- critC$any.decline
critC$declineC[!critC$declineC %in% "Decreasing" & est.decline %in% "Decreasing"] <- "Decreasing" 
critC$declineC[critC$declineC %in% "Increasing" | critC$declineC %in% "Stable"] <- "Not Decreasing" 

#Getting estimates of species continuing decline based on habitat loss (Criterion B)
eoo.decline <- readRDS("data/EOO_hab_loss_2000_2015.rds")
eoo.decline$declineB <- eoo.decline$rel.loss
eoo.decline$declineB[!is.na(eoo.decline$declineB) & eoo.decline$declineB >= 1] <- 1
eoo.decline$declineB[!is.na(eoo.decline$declineB) & as.double(eoo.decline$declineB) < 1 & as.double(eoo.decline$declineB) >= 0] <- 0
eoo.decline$declineB[is.nan(eoo.decline$declineB)] <- NA
eoo.decline$declineB[eoo.decline$declineB %in% 1] <- "Decreasing"
eoo.decline$declineB[eoo.decline$declineB %in% 0] <- "Not Decreasing"
table(eoo.decline$declineB, useNA = "always")

#Correcting the decrease for pioneer species
hab <- read.csv("data/threat_habitats.csv", as.is = TRUE)
tmp <- merge(eoo.decline, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
tmp <- tmp[order(tmp$tax),]
table(tmp$tax == eoo.decline$tax)
eoo.decline$net.loss <- eoo.decline$recover - eoo.decline$loss 
hist(eoo.decline$net.loss, nclass=80)
eoo.decline$declineB[eoo.decline$net.loss >=0.1 & tmp$ecol.group %in% "pioneer"] <- "Not Decreasing"
table(eoo.decline$declineB, useNA = "always")

#Merging any decline info with the criterion B metrics
critB_opt <- merge(critB_opt, critC[,c("species","declineC")],
                    by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_opt <- merge(critB_opt, eoo.decline[,c("tax","declineB")],
                   by = "tax", all.x = TRUE, sort = FALSE)
critB_opt <- critB_opt[order(critB_opt$tax),]

critB_low <- merge(critB_low, critC[,c("species","declineC")],
                   by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_low <- merge(critB_low, eoo.decline[,c("tax","declineB")],
                   by = "tax", all.x = TRUE, sort = FALSE)
critB_low <- critB_low[order(critB_low$tax),]

critB_high <- merge(critB_high, critC[,c("species","declineC")],
                   by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_high <- merge(critB_high, eoo.decline[,c("tax","declineB")],
                   by = "tax", all.x = TRUE, sort = FALSE)
critB_high <- critB_high[order(critB_high$tax),]

#Comparing the decline from criterion C and B and correcting if necessary (priority for abudance decline over EOO decline) For assuming, that all other species (rarer) are decreasing as well
hab <- read.csv("data/threat_habitats.csv", as.is = TRUE)
tmp <- merge(critB_opt, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
table(critB_opt$tax == tmp$tax)
critB_opt$declineB[critB_opt$declineB %in% "Not Decreasing" & critB_opt$declineC %in% "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
#table(critB_opt$declineB, critB_opt$declineC, useNA = "always")
tmp <- merge(critB_low, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
table(critB_low$tax == tmp$tax)
critB_low$declineB[critB_low$declineB %in% "Not Decreasing" & critB_low$declineC %in% "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
#table(critB_low$declineB, critB_low$declineC, useNA = "always")
tmp <- merge(critB_high, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
table(critB_high$tax == tmp$tax)
critB_high$declineB[critB_high$declineB %in% "Not Decreasing" & critB_high$declineC %in% "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
#table(critB_high$declineB, critB_high$declineC, useNA = "always")

#Combining the info on number of localities and % in PAs
require(dplyr)
critB_opt <- 
  critB_opt %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)
critB_low <- 
  critB_low %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)
critB_high <- 
  critB_high %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)

#Creating the severely fragmented colum
critB_opt$sever.frag <- 100 * critB_opt$Nbe_subPop/ (critB_opt$AOO / 4) > 50
critB_low$sever.frag <- 100 * critB_low$Nbe_subPop/ (critB_low$AOO / 4) > 50
critB_high$sever.frag <- 100 * critB_high$Nbe_subPop/ (critB_high$AOO / 4) > 50


#### PERFORMING THE ASSESSMENTES OF CRITERION B ####
results_Cb_opt <- cat_criterion_b(EOO = critB_opt$EOO,
                              AOO = critB_opt$AOO,
                              locations = critB_opt$nbe_loc_total,
                              sever.frag = critB_opt$sever.frag,
                              #protected = critB_opt$protected, 
                              decline = critB_opt$declineB,
                              protected.threshold = 100 
)
sum((100 * table(results_Cb_opt$ranks_B, useNA = "always")/dim(critB_opt)[1])[c(1,2,4)]) #17.25%; now 15.65%

results_Cb_low <- cat_criterion_b(EOO = critB_low$EOO,
                                   AOO = critB_low$AOO,
                                   locations = critB_low$nbe_loc_total,
                                   sever.frag = critB_low$sever.frag,
                                   #protected = critB_low$protected, 
                                   decline = critB_low$declineB,
                                   protected.threshold = 100
)
sum((100 * table(results_Cb_low$ranks_B, useNA = "always")/dim(critB_low)[1])[c(1,2,4)]) #14.68%; now 13.57%

results_Cb_high <- cat_criterion_b(EOO = critB_high$EOO,
                                   AOO = critB_high$AOO,
                                   locations = critB_high$nbe_loc_total,
                                   sever.frag = critB_high$sever.frag,
                                   #protected = critB_high$protected, 
                                   decline = critB_high$declineB,
                                   protected.threshold = 100 
)
sum((100 * table(results_Cb_high$ranks_B, useNA = "always")/dim(critB_high)[1])[c(1,2,4)]) #26.95%; now 24.99%


#Saving the results
results_Cb_opt <- do.call(cbind.data.frame, c(results_Cb_opt, stringsAsFactors = FALSE))
critB_opt.all <- critB_opt[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_opt.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_opt.all[, c("category_B", "category_B_code","category_B1","category_B2")] <-
  results_Cb_opt
critB_opt.all[is.na(critB_opt.all$AOO),]

results_Cb_low <- do.call(cbind.data.frame, c(results_Cb_low, stringsAsFactors = FALSE))
critB_low.all <- critB_low[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_low.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_low.all[, c("category_B", "category_B_code","category_B1","category_B2")] <-
  results_Cb_low
critB_low.all[is.na(critB_low.all$AOO),]

results_Cb_high <- do.call(cbind.data.frame, c(results_Cb_high, stringsAsFactors = FALSE))
critB_high.all <- critB_high[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_high.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_high.all[, c("category_B", "category_B_code","category_B1","category_B2")] <-
  results_Cb_high
critB_high.all[is.na(critB_high.all$AOO),]

#Saving
saveRDS(as.data.frame(critB_opt.all), "data/criterionB_all_optim_tax_confidence.rds")
saveRDS(as.data.frame(critB_low.all), "data/criterionB_all_low_tax_confidence.rds")
saveRDS(as.data.frame(critB_high.all), "data/criterionB_all_high_tax_confidence.rds")



#### FIGURE: OPTIMUM VS. HIGH CONFIDENCE LEVEL ####
critB_opt.all <- readRDS("data/criterionB_all_optim_tax_confidence.rds")
critB_high.all <- readRDS("data/criterionB_all_high_tax_confidence.rds")
res <- readRDS("data/tax_conf_effect_on_RLI.rds")
require(circlize)

jpeg(filename = "figures/Figure_SU_new.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(1,2))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)

## OPTIMUM VS. HIGH
mat <- as.matrix(table(paste0(critB_opt.all$category_B,"_opt"), 
                       paste0(critB_high.all$category_B,"_hi")))
mat <- mat[c(1,2,5,3), c(3,5,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
grid.col = c(CR_hi = "red", EN_hi = "darkorange", VU_hi = "gold", LCorNT_hi = "khaki",
             CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 15] = mat[mat < 15]*2
mat[mat > 0 & mat < 5] = 10

#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible) = FALSE
chordDiagram(mat, big.gap = 15, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
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
  if(si == "VU_hi") {
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
legend("topright","Optim. confidence", bty="n", cex=1.2)
legend("topleft",legend=expression(bold("A")),
       bty="n",horiz=F,cex=1.5,x.intersp=-1,y.intersp=-0.3)

## ANY VS. HIGH (very similar to opt vs. high: not plotting)
# mat <- as.matrix(table(paste0(critB_low.all$category_B,"_lo"), paste0(critB_high.all$category_B,"_hi")))
# mat <- mat[c(1,2,4,3), c(3,5,2,1)]
# colnames(mat) <- gsub(" ", "", colnames(mat))
# rownames(mat) <- gsub(" ", "", rownames(mat))
# 
# #Defining the colors of tracks and links
# grid.col = c(CR_hi = "red", EN_hi = "darkorange", VU_hi = "gold", LCorNT_hi = "khaki",
#              CR_lo = "red", EN_lo = "darkorange", VU_lo = "gold", LCorNT_lo = "khaki")
# col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
# col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
# col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
# #col_mat[mat < 5] = "#00000000"
# mat[mat < 15] = mat[mat < 15]*2
# mat[mat > 0 & mat < 5] = 10
# 
# #plotting the diagram
# circos.clear()
# circos.par(start.degree = 90)
# visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
# #diag(visible) = FALSE
# lava::revdiag(visible) = FALSE
# chordDiagram(mat, big.gap = 15, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
#              grid.col = grid.col, col = col_mat,
#              self.link = 1, link.visible = visible,
#              #h=0.9,
#              #w=1,
#              #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
#              link.lwd = 4
#              #h.ratio = 0.7
#              #reduce_to_mid_line = FALSE,
#              #w2=0.5,
#              #rou=0.2
#              #point1 = rep(0,16)
# )
# #Putting legends on
# sec.ind <- c("CR","EN","VU","LC or NT","LC or NT","VU","EN","CR")
# for(si in get.all.sector.index()) {
#   lab <- sec.ind[which(si == get.all.sector.index())]
#   xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#   ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#   if(si == "VU_hi") {
#     circos.text(mean(xlim), mean(ylim), labels = "VU", 
#                 sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
#                 facing = "bending", niceFacing = FALSE, col = "black")
#     
#   } else {
#     circos.text(mean(xlim), mean(ylim), labels = lab, 
#                 sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
#                 facing = "bending", niceFacing = TRUE, col = "black")
#   }  
# }
# legend("topleft","High confidence", bty="n", cex=1.2)
# legend("topright","Any confidence", bty="n", cex=1.2)

## Assessing the confidence level cutoff 
par(mar=c(3,3.5,1.5,0.5))
cut <- seq(0.05,0.95, by=0.05)
ids <- cut < 0.85
par(mar=c(3,3.5,0.75,0.5), mgp=c(2,0.25,0), tcl=-0.2, las=1)
plot(res[,2][ids] ~ cut[ids], ylim = c(0.85, .92),
     xlab = "Tax. confidence level", ylab = "Red List Index", 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2, pch=19)
axis(1, at=cut, cex.axis = 1, labels = cut * 100)
lines(res[,1][ids] ~ cut[ids] , lty = 2)
lines(res[,3][ids] ~ cut[ids] , lty = 2)
# abline(h = opt.rli[2])
# abline(h = high.rli[2], col = 2)
points(res[,5][ids] ~ cut[ids], col = "red", pch=19)
lines(res[,4][ids] ~ cut[ids], col = "red", lty = 2)
lines(res[,6][ids] ~ cut[ids], col = "red", lty = 2)
# abline(v=c(0.6,0.75))
legend("topleft", expression(bold(B)), bty="n", cex=1.3,
       x.intersp=-0.7, y.intersp=0.1)
dev.off()


#### GILLES' EXTRA CODES ####

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

