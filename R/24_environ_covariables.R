#####################################################H
#### GETTING ALTITUDE DATA FOR THREAT OCCURRENCES ####
#####################################################H
rm(list = ls())

### Loading packages
require(raster)
require(rgdal)
require(dplyr)
require(parallel)
require(doParallel) 	

#### EXTRACTING ALL ENVIRONMENTAL CO-VARIABLES FROM TREECO ####

## Loading plot coordinates
occs <- readRDS("data/threat_species_by_country.rds")
occs1 <- occs
occs1[ , ordem := .I]
occs1[ , latitude.work1 := as.numeric(latitude.work1)]
occs1[ , longitude.work1 := as.numeric(longitude.work1)]
occs2 <- occs1[!is.na(longitude.work1), ]
occs2 <- occs2[!is.na(latitude.work1), ]
coords1 <- as.data.frame(occs2)

## Getting data frame with spatial coordinates (points)
coordinates(coords1) <- c("longitude.work1", "latitude.work1")  # set spatial coordinates

## Defining geographical projection
crs.geo <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs")  # geographical, datum WGS84
proj4string(coords1) <- crs.geo  # define projection system of our data
surveys1 = coords1

##file to receive the environmental data
# dados <- surveys1$ordem
dados <- data.frame(dados = surveys1$ordem)


##DEM - ALTITUDE##
#A altimetria foi baixada em metros
path0 = "E://ownCloud//W_GIS//Am_Lat_MDT//dem"
myfiles <- list.files(path = path0, pattern = "tif",  full.names = TRUE)
myfiles <- myfiles[grep("tif.aux.xml|tif.ovr",myfiles,invert=TRUE)]

#extracting the distance to the cosatline for each survey (paralellized)
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(rgdal))
clusterExport(cl, list("surveys1", "myfiles"),envir=environment())
tmp2 = parLapply(cl,1:length(myfiles), fun= function(i) {
  # loading the corresponding landscape from the path folder
  dcl <- try(raster(myfiles[i]), TRUE)
  # getting the coordinate of the survey
  ext <- round(extent(dcl), 0)
  locs <- surveys1[surveys1@coords[,"latitude.work1"] >= ext@ymin & 
                     surveys1@coords[,"latitude.work1"] <= ext@ymax,]
  locs <- locs[locs@coords[,"longitude.work1"] >= ext@xmin & 
                 locs@coords[,"longitude.work1"] <= ext@xmax,]
  # checking if there are coordinates within tile ranges
  if(length(locs) == 0) { 
    #res = cbind.data.frame(NA,NA,NA)
    #names(res) = c("ordem","value","value1") 
  } else { 
    # homogeneizing projections
    pj = crs(dcl)
    locs1 = try(spTransform(locs, pj),TRUE)
    # extracting the raster values
    cc.ex = extract(dcl, locs1)
    # extracting the raster values from nearby
    # cc.prox = extract(dcl,locs1, fun=median, buffer=30*sqrt(2)+(30*sqrt(2)/2))
    #cat(i,"\n")	
    # saving the results
    # res = cbind.data.frame(locs1$ordem,cc.ex,cc.prox)
    res = cbind.data.frame(locs1$ordem,cc.ex)
    # names(res) = c("ordem","value","value1")
    names(res) = c("ordem","value")
    #if(is.null(land)) land = res else land = rbind.data.frame(land,res)
    res
  }
})
stopImplicitCluster()
stopCluster(cl)
alt = do.call("rbind.data.frame", tmp2)
alt = merge(dados, alt, by.x="dados", by.y="ordem",all.x=TRUE) 
alt = alt[match(dados$dados, alt$dados, nomatch=0),]
table(alt$dados == dados$dados)
dados$ALT <- alt$value
dados$ALT[!is.na(dados$ALT) & dados$ALT < 0] <- 0

#########################
#### SAVING ENV DATA ####
#########################
names(dados)[1] = "ordem"
dados1 <- merge(occs1, dados, by="ordem")
data.table::setkeyv(dados1, "ordem")
dados1[ , ordem := NULL]
saveRDS(dados1, file = "data/threat_species_by_country_environmental.rds", compress = "gzip")

