############################
#### ESTIMATING POP SIZE ###
############################

### THING TO DO ###
# Growth-form specific density corrections for studies using DBH cutoffs different than 5.0cm
# Exclude non Atlantic Forest sites - OK!
require(gstat)
require(rgdal)
require(rgeos)
require(sp)
require(raster)
require(plotrix)

### READING DATA

#Reading abundance data
abund = read.csv("D://Documentos//Pos Doc//Databases//Species abundances//abundances.csv",as.is=T)
abund = abund[1:(dim(abund)[1]-4),]

#Reading plot data
refs = read.csv("D://Documentos//Pos Doc//Databases//References//references.csv",as.is=T)

### PREPARING VEG TABLE

#Removing non-tree and dead individuals
abund = abund[grep("Liana|liana",abund$obs,inv=T),] ## removing the lianas
abund = abund[grep("Bamboo|bamboo",abund$obs,inv=T),] ## removing bamboos
abund = abund[grep("Mortas", abund$species.correct, inv=T),] ## removing dead individuals
abund = abund[grep("Correcao", abund$species.correct, inv=T),] ## removing dead individuals

#Considering species determined as cf. and removing varieties or subspecies
abund$species.correct = gsub(" cf\\. "," ",abund$species.correct) # removing the cf. from species names
abund$species.correct = sapply(strsplit(abund$species.correct, " var\\."),function(x) x[1]) # removing identifications under species level
abund$species.correct = sapply(strsplit(abund$species.correct, " subsp\\."),function(x) x[1]) # removing identifications under species level

#Removing records not determined at species level
abund1 = abund[grep(" aff\\. |Indet\\.|sp\\.|Indets\\.", abund$species.correct, inv=T),]

#Filtering only SiteCodes with abundance data
refs = refs[refs$SiteCode %in% unique(abund1$SiteCode),] 

#Getting SiteCode list
sites = refs[grep("Atlantic|Pampa",refs$domain),]$SiteCode
length(sites)

#Filtering only Atlantic Forest sites 
abund2 = abund1[abund1$SiteCode %in% sites,]
refs1 = refs[refs$SiteCode %in% sites,]

#Verifying if all sites were included
tmp = as.character(refs1[refs1$SiteCode %in% sites,]$ordem) 
tmp[!tmp %in% unique(unlist(strsplit(abund2$ordem,"\\|")))] #OK

#Obtaining the vegetation table per SiteCode
tmp = xtabs(abund2$N ~ abund2$SiteCode + abund2$species.correct) #contigency table
vegtab = as.data.frame.matrix(tmp) #transforming table into data.frame
vegtab = vegtab[grep("_",row.names(vegtab),inv=T),] # removing abundance data for DBH cutoff < 5cm
#sites = row.names(vegtab) # removing SiteCode names for DBH cutoff < 5cm
spec = colnames(vegtab) # number of species
n_spec = length(spec) # number of species


### PREPARING SITE TABLE
#Getting total sampling effort for each SiteCode
effort = aggregate(as.double(refs1$effort_ha),list(refs1$SiteCode),sum,na.rm=T)
effort[effort$Group.1 == "SPibiu1",]$x = 2.5

#Getting dbh cutoff criteria for each SiteCode (original, results and abundance digitalized)
dbh = aggregate(as.character(refs1$dbh_cutoff1),list(refs1$SiteCode),function(x) unique(as.character(x)))
dbh_res = aggregate(as.character(refs1$dbh_results),list(refs1$SiteCode),function(x) unique(as.character(x)))
dbh_dig = aggregate(as.character(refs1$status_diagnostico),list(refs1$SiteCode),function(x) unique(as.character(x)))

tmp = merge(dbh,dbh_res,by="Group.1") #merging dbh cutoffs
dbh = merge(tmp,dbh_dig,by="Group.1") #merging dbh cutoffs

#Replacing the original dbh by the dbh that actually entered the abundance dataset 
dbh[!dbh$x.y %in% "" & !dbh$x %in% c("ok NAB digitado! DAP>10","ok NAB digitado! Checado DAP>10"),]$x.x = "DBH>=4.8-5.0cm"
dbh[dbh$x %in% c("ok NAB digitado! DAP>10","ok NAB digitado! Checado DAP>10"),]$x.x = "DBH>=10.0cm"
dbh[grep("DBH>=4.8-5.0cm",dbh$x.x),]$x.x = "DBH>=4.8-5.0cm"
dbh[grep("DGH>=2.9cm",dbh$x.x),]$x.x = "DGH>=3.0-3.2cm"
dbh[grep("DBH30>=5.0cm",dbh$x.x),]$x.x = "DGH>=4.8-5.0cm"
dbh[grep("DBH>=9.5-9.6cm",dbh$x.x),]$x.x = "DBH>=10.0cm"

##DBH cutoff correction based on individual sites
#DBH>=2.5-2.6cm => DBH>=5.0cm = *0.55 
#DGH>=3.0-3.2cm => DBH>=5.0cm = *0.50
#DBH>=3.0-3.2cm => DBH>=5.0cm = *0.65
#DGH>=4.8-5.0cm => DBH>=5.0cm = *0.70
#DBH>=10.0cm => DBH>=5.0cm = *1.5

##DBH cutoff correction based on individual sites
refs1[! refs1$dbh_result %in% "",]$dbh_cutoff1 = refs1[! refs1$dbh_result %in% "",]$dbh_results
refs1[grep("DBH>=4.8-5.0cm",refs1$dbh_cutoff1),]$dbh_cutoff1 = "DBH>=4.8-5.0cm"
refs1[grep("DGH>=2.9cm",refs1$dbh_cutoff1),]$dbh_cutoff1 = "DGH>=3.0-3.2cm"
refs1[grep("DBH30>=5.0cm",refs1$dbh_cutoff1),]$dbh_cutoff1 = "DGH>=4.8-5.0cm"
refs1[grep("DBH>=9.5-9.6cm",refs1$dbh_cutoff1),]$dbh_cutoff1 = "DBH>=10.0cm"
DA = as.data.frame.matrix(aggregate(as.double(refs1$DA),list(refs1$dbh_cutoff1),mean,na.rm=T))
DA$ratio = round(DA[DA$Group.1 == "DBH>=4.8-5.0cm",]$x/DA$x,2)

#Removing info (duplicated lines) for multiple-plot SiteCodes
refs1 = refs1[!duplicated(refs1$SiteCode),] 

#Getting coordinates for each survey
refs1[is.na(refs1$lat_correct),]$lat_correct = refs1[is.na(refs1$lat_correct),]$lat
refs1[is.na(refs1$long_correct),]$long_correct = refs1[is.na(refs1$long_correct),]$long

#Merging SiteCode coordinates with the contigency table 
data_env = merge(refs1[,c("SiteCode","long_correct","lat_correct"),],vegtab,by.x="SiteCode",by.y="row.names")
data_env$long_correct = as.double(data_env$long_correct) 
data_env$lat_correct = as.double(data_env$lat_correct)

#Obtaining the same database but accounting for sampling effort (density = trees/ha)
data_env_dens = data_env
table(effort$Group.1 == data_env$SiteCode) #checking order of SiteCode
data_env_dens[,4:dim(data_env)[2]] = data_env_dens[,4:dim(data_env)[2]]/effort$x #dividing abundance by samplling effort

#Obtaining the same database but accounting for sampling effort (density = trees/ha) and dbh cutoff criteria (corrections for DBH>=5.0cm)
data_env_dens_dbh = data_env_dens
table(dbh$Group.1 == data_env_dens_dbh$SiteCode) #checking order of SiteCode
data_env_dens_dbh[dbh$x.x == "DBH>=10.0cm",4:dim(data_env)[2]] = data_env_dens_dbh[dbh$x.x == "DBH>=10.0cm",4:dim(data_env)[2]]*2 #DBH>=10.0cm
data_env_dens_dbh[dbh$x.x == "DBH>=2.5-2.6cm",4:dim(data_env)[2]] = data_env_dens_dbh[dbh$x.x == "DBH>=2.5-2.6cm",4:dim(data_env)[2]]*0.45 #DBH>=2.5-2.6cm
data_env_dens_dbh[dbh$x.x == "DBH>=3.0-3.2cm",4:dim(data_env)[2]] = data_env_dens_dbh[dbh$x.x == "DBH>=3.0-3.2cm",4:dim(data_env)[2]]*0.65 #DBH>=3.0-3.2cm
data_env_dens_dbh[dbh$x.x == "DGH>=3.0-3.2cm",4:dim(data_env)[2]] = data_env_dens_dbh[dbh$x.x == "DGH>=3.0-3.2cm",4:dim(data_env)[2]]*0.50 #DGH>=3.0-3.2cm
data_env_dens_dbh[dbh$x.x == "DGH>=4.8-5.0cm",4:dim(data_env)[2]] = data_env_dens_dbh[dbh$x.x == "DGH>=4.8-5.0cm",4:dim(data_env)[2]]*0.75 #DGH>=4.8-5.0cm


### PREPARING THE ATLANTIC FOREST GRID FOR PREDICTIONS
## Defining the general path to read/save files (tifs and shape files)
path = "C://Users//RenatoLima//ownCloud//W_GIS"

## Reading the shapefiles with Brazil and its states limits
af.lim = readOGR(dsn=paste(path,"//AF_Biorregions_TNC",sep=""),layer="Floresta_Atlantica_TNC_2")
#af.lim = readOGR(dsn=paste(path,"//AF_limits_lei_11428",sep=""),layer="mata_atlantica11428")
#af.lim = readOGR(dsn=paste(path,"//AF_limits_extended",sep=""),layer="merge_limites_MA_buffer20km_ajustado")
brasil <- readOGR(dsn=paste(path,"//BR_ADM_limits",sep=""),layer=paste("BRA_adm",0,sep=""))
estados <- readOGR(dsn=paste(path,"//BR_ADM_limits",sep=""),layer=paste("BRA_adm",1,sep=""))

## Buffering around the AF limits shapefile
#af.lim1 = gBuffer(af.lim,width = 0.1,byid=F) ##NEED TO PROJECT SHAPEFILE

### Define the grid extent:
bb = bbox(af.lim)

x.range <- as.numeric(c(bb[1,1],bb[1,2]))  # min/max longitude of the interpolation area
y.range <- as.numeric(c(bb[2,1],bb[2,2]))  # min/max latitude of the interpolation area

#Create a data frame from all combinations of the supplied vectors or factors. See the description of the return value for precise details of the way this is done. Set spatial coordinates to create a Spatial object. Assign gridded structure:
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1), 
	y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE
proj4string(grd) <- crs(af.lim)


#filtering only the grid inside the AF limits
grd1 = over(grd,af.lim)
data_to_pred = grd[!is.na(grd1$ECO_CODE),]
#data_to_pred = grd[!is.na(grd1$OBJECTID_1),]
data_to_pred = as.data.frame(data_to_pred)
names(data_to_pred) = c("Longitude","Latitude")

i=1249

t1 = Sys.time()
for (i in 1:n_spec){
  taxon_data = cbind.data.frame(data_env_dens_dbh$long_correct,data_env_dens_dbh$lat_correct,data_env_dens_dbh[,spec[i]])

  # adjust names to conform with data to predict
  names(taxon_data) = c("Longitude","Latitude","taxon")

  nmx = 100; mxd = 2
  taxon.gstat = gstat(id = "taxon", formula = taxon ~ 1, 
                      locations = ~ Latitude + Longitude,
                      data = taxon_data, nmax = nmx, maxdist = mxd,
                      set = list(idp = 2))
  #taxon.pred = predict.gstat(taxon.gstat, data_to_pred) ##ASK HANS FOR FUNCTION 'predict.gstat'
  taxon.pred = predict(taxon.gstat, data_to_pred)$taxon.pred	
  taxon.pred[is.na(taxon.pred)] = 0	
	
  MapSum[,i] = taxon.pred
  cat("iteration ",i, " completed", "\n")
}
Sys.time() - t1

#loess fit to predicted data
sp=0.15
z.loess = loess(taxon.pred ~ data_to_pred$Longitude * data_to_pred$Latitude,
                  span = 0.25, degree = 2, se = T,
                  normalize = TRUE, family = "gaussian",
                  surface = "direct") #!surface is direct to be able to extrapolate
#calculate the predicted values for the AF grid
grid.z.predict = predict(z.loess, data_to_pred, se = T)

#replace all fits < co by zero
grid.z.predict$fit[grid.z.predict$fit < 0] = 0

    grid.col   = vector(length = length(data_to_pred$Longitude))
    grid.min   = min(grid.z.predict$fit, na.rm = TRUE)
    grid.min1  = min(taxon.pred, na.rm = TRUE)
    grid.max   = max(grid.z.predict$fit, na.rm = TRUE)
    grid.max1  = max(taxon.pred, na.rm = TRUE)
    grid.range = grid.max - grid.min
    grid.range1= grid.max1 - grid.min1    

    n.colors = 256
    grid.pal = colorRampPalette(c("white", "black"))(n.colors) ## (n)

    grid.col   = 1 -(grid.z.predict$fit - grid.min)/grid.range
    grid.col1  = 1 -(taxon.pred - grid.min)/grid.range

    grid.col   = grid.pal[1+round((n.colors-1)*(grid.z.predict$fit - grid.min)/grid.range)]
    grid.col1  = grid.pal[1+round((n.colors-1)*(taxon.pred - grid.min1)/grid.range1)]

    taxon_data.xy = taxon_data 	 
    coordinates(taxon_data.xy) <- ~ Longitude + Latitude
    proj4string(taxon_data.xy) <- crs(af.lim)
    cex = rescale(taxon_data.xy$taxon, c(0.5,3))
    cex[is.infinite(log(taxon_data.xy$taxon))] = NA


par(mfrow=c(1,3))
par(mai=c(0.275,0.375,0.325,0.025),mgp=c(2,0.2,0),tcl=-0.2,las=1)
    plot(data_to_pred$Longitude,data_to_pred$Latitude, 
         main = paste(spec[i]," - observations",sep=""),
         xlab = "Longitude", ylab = "Latitude",
         xlim = c(-50, -45), ylim = c(-33,-3), asp = 30/30,
         xaxp = c(-60, -35, 10), yaxp = c(-33, -3, 4),
         pch = 22, cex = 0.1,
         col = "white", bg = "white")
    plot(af.lim,add=T,border=grey(0.5))	
    plot(brasil,add=T,) 
    plot(estados,add=T)
    points(taxon_data.xy,col="red",pch=19,cex=cex)
    legend ("bottomright", c(paste("mean= ",round(mean(taxon_data$taxon),2),"; sd= ",round(sd(taxon_data$taxon),2),sep=""),
		paste("median= ",round(median(taxon_data$taxon),1),sep=""),
		paste("max= ",round(max(taxon_data$taxon),1),sep="")),bty="n",cex=1.2)

    plot(data_to_pred$Longitude,data_to_pred$Latitude, 
         main = paste(spec[i]," - loess (span= ",sp,")",sep=""),
         xlab = "Longitude", ylab = "Latitude",
         xlim = c(-50, -45), ylim = c(-33,-3), asp = 30/30,
         xaxp = c(-60, -35, 10), yaxp = c(-33, -3, 4),
         pch = 22, cex = 0.1,
         col = grid.col, bg = grid.col)
    plot(af.lim,add=T,border=grey(0.5))	
    plot(brasil,add=T,) 
    plot(estados,add=T)
    legend ("bottomright", c(paste("mean= ",round(mean(grid.z.predict$fit),2),"; sd= ",round(sd(grid.z.predict$fit),2),sep=""),
		paste("median= ",round(median(grid.z.predict$fit),1),sep=""),
		paste("max= ",round(max(grid.z.predict$fit),1),sep="")),bty="n",cex=1.2)


    plot(data_to_pred$Longitude,data_to_pred$Latitude, 
         main = paste(spec[i]," - gstat (nmax = ",nmx,", maxd = ",mxd,")",sep=""),
         xlab = "Longitude", ylab = "Latitude",
         xlim = c(-50, -45), ylim = c(-33,-3), asp = 30/30,
         xaxp = c(-60, -35, 10), yaxp = c(-33, -3, 4),
         pch = 22, cex = 0.1,
         col = grid.col1, bg = grid.col1)
    plot(af.lim,add=T,border=grey(0.5))	
    plot(brasil,add=T,) 
    plot(estados,add=T)
    legend ("bottomright", c(paste("mean= ",round(mean(taxon.pred),2),"; sd= ",round(sd(taxon.pred),2),sep=""),
		paste("median= ",round(median(taxon.pred),1),sep=""),
		paste("max= ",round(max(taxon.pred),1),sep="")),bty="n",cex=1.2)
