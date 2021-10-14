####################################
####################################
#### MAPS OF THREATENED SPECIES ####
####################################
####################################
rm(list=ls())
gc()

## Loading packages
require(data.table)
require(sp)
require(gstat)

## Getting only the necessau shapefiles
path <- "E://ownCloud//W_GIS" # Renato's path
am.lat <- rgdal::readOGR(dsn=paste(path,"//Am_Lat_ADM_ArcGis",sep=""),layer="LatinAmerica")
brasil = readRDS(paste(path,"//Am_Lat_ADM_GADM_v3.6//gadm36_BRA_1_sp.rds",sep=""))
af <- rgdal::readOGR(dsn=paste(path,"//AF_limites_milton",sep=""),layer="merge_limites_MA11428_TNC_ARG_BR_PAR")
af <- rgeos::gBuffer(af, byid=TRUE, width=0) #correcting possible overlapping polygons
af <- cleangeo::clgeo_Clean(af) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
af <- raster::aggregate(af, dissolve=TRUE)
#projecting
am.lat.proj= spTransform(am.lat,  CRS("+init=epsg:5641"))
brasil.proj= spTransform(brasil,  CRS("+init=epsg:5641"))
af.proj = spTransform(af, CRS("+init=epsg:5641"))



####################################################################################################################################################################################H
####################################################################################################################################################################################H
#### FIGURES ####
#################

# -------------------------------------------------------------------------
##################################
#### ADAPTIVE RESOLUTION GRID ####
##################################
## reading the results per grid cell
grid.result = readRDS("data/grid.results_adaptive_resol.rds")#[[1]]
# samp.cov.cutoff = readRDS("data/grid.results_adaptive_resol.rds")[[2]]
samp.cov.cutoff = quantile(grid.result$SampCover[!is.na(grid.result$N.total) & 
                                                   grid.result$N.total >= 10], 
                           prob=c(0.05,0.1,0.25,0.5,0.75,0.8,0.85,0.9,0.95,0.99))


## Defining the sampling coverage cut-off for analysis
#For rarefaction, a temptative minimum sampling is:
plot(RLI ~ N.total, grid.result, log="x"); abline(v=c(30,50,100))
plot(RLI.end ~ N.total, grid.result, log="xy"); abline(v=c(30,50,100))
plot(RLI ~ S.total, grid.result, log="x"); abline(v=c(30,50,100, 250))
plot(RLI.end ~ S.total, grid.result, log="xy"); abline(v=c(30,50,100, 250))

plot(Threat ~ N.total, grid.result, log="x"); abline(v=c(30,50,100))
plot(Threat.end ~ N.total, grid.result, log="xy"); abline(v=c(30,50,100))
plot(Threat ~ S.total, grid.result, log="x"); abline(v=c(30,50,100, 250))
plot(Threat.end ~ S.total, grid.result, log="xy"); abline(v=c(30,50,100, 250))
plot(Threat.pa ~ S.total, grid.result, log="x"); abline(v=c(30,50,100, 250))
plot(Threat.end.pa ~ S.total, grid.result, log="x"); abline(v=c(30,50,100, 250))


#What does it mean to remove cells above a given quantile of the sampling coverage distribution?
result = NULL; result1 = NULL
for (i in 1:length(samp.cov.cutoff)) {
  toto = samp.cov.cutoff[i]
  toto1= summary(grid.result$N.total[grid.result$SampCover < toto])
  toto2= summary(grid.result$N.total[grid.result$SampCover >= toto])
  toto3= c(toto1[c(3,4,1,6)],toto2[c(3,4,1,6)]) 
  toto4= summary(grid.result$S.total[grid.result$SampCover < toto])
  toto5= summary(grid.result$S.total[grid.result$SampCover >= toto])
  toto6= c(toto4[c(3,4,1,6)],toto5[c(3,4,1,6)]) 
  if(is.null(result)) result = toto3 else result = rbind.data.frame(result, toto3)
  if(is.null(result1)) result1 = toto6 else result1 = rbind.data.frame(result1, toto6)
}
row.names(result) = row.names(result1) = names(samp.cov.cutoff)
colnames(result) = colnames(result1) = c("med.excl","mean.excl","min.excl","max.excl","med.incl","mean.incl","min.incl","max.incl")
result
result1
cut = samp.cov.cutoff[3] # using the 50% quantile (mean) of the sampling coverage distribution
head(sort(grid.result$N.total[grid.result$SampCover>cut]), 50) ## REMOVE ALSO CELLS WITH LESS THAN < 30 OCCURRENCES

#Inspecting other possibilities of data filtering
ns <- c(50,100,250,500)
spp <- c(30,50,100,250)
grid.result$select = ((grid.result$SampCover >= samp.cov.cutoff[5]  & grid.result$N.total >= ns[1]) | 
                        (grid.result$SampCover >= samp.cov.cutoff[4] & grid.result$SampCover < samp.cov.cutoff[5] & grid.result$N.total >= ns[2]) | 
                        (grid.result$SampCover >= samp.cov.cutoff[3] & grid.result$SampCover < samp.cov.cutoff[4] & grid.result$N.total >= ns[3]))
grid.result$select1 = ((grid.result$SampCover >= samp.cov.cutoff[5]  & grid.result$S.total >= spp[1]) | 
                         (grid.result$SampCover >= samp.cov.cutoff[4] & grid.result$SampCover < samp.cov.cutoff[5] & grid.result$S.total >= spp[2]) | 
                         (grid.result$SampCover >= samp.cov.cutoff[3] & grid.result$SampCover < samp.cov.cutoff[4] & grid.result$S.total >= spp[3]))
grid.result$select[grid.result$N.total < 250] <- FALSE
grid.result$select1[grid.result$N.total < 250] <- FALSE

#Result by N.total
par(mfrow = c(1,2))
plot(grid.result$SampCover ~ grid.result$N.total, log="x")
abline(v= ns, col=c(1,2,3,4)); abline(h=samp.cov.cutoff[3:5], col=c(3,2,1))
points(SampCover ~ N.total, subset(grid.result,select), col=2, pch=21)
summary(grid.result$N.total[grid.result$select])
head(sort(grid.result$N.total[grid.result$select]),50)
summary(grid.result$SampCover[grid.result$select])
#Result by S.total
plot(grid.result$SampCover ~ grid.result$S.total, log="x")
abline(v= spp, col=c(1,2,3,4)); abline(h=samp.cov.cutoff[3:5], col=c(3,2,1))
points(SampCover ~ S.total, subset(grid.result,select1), col=2, pch=21)
summary(grid.result$S.total[grid.result$select1])
head(sort(grid.result$S.total[grid.result$select1]),50)
summary(grid.result$SampCover[grid.result$select1])

## For the SOM
summary(grid.result$N.total[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$S.total[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$RLI.pa[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$Threat.pa[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$Threat.end.pa[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$SampCover[!is.na(grid.result$select) & grid.result$select])

############################################################H
#### CENTRES OF HIGH CONCENTRATION OF THREATENED SPECIES ####
############################################################H
## Getting the grid and editing it
grid0 <- readRDS("./data/adaptative_resolution_AF_min200_max500_cellmax2/clean_grid.rds")

## removing overllaping parts of the grid cells
ids.polys <- grid0$ID
grid_sf <- sf::st_as_sf(grid0)
grid_inter <- sf::st_intersection(grid_sf)
# x_inter <- x_inter[x_inter$n.overlaps == 1,]
grid_inter <- grid_inter[rownames(grid_inter) %in% ids.polys,]
grid <- as(grid_inter, "Spatial")

## Making sure grid IDs are the same internally
grid1 <- grid
tmp <- SpatialPoints(coords = coordinates(grid), proj4string = raster::crs(grid))
tmp1 <- over(tmp, grid)
for(i in 1:length(grid1)) { slot(slot(grid1, "polygons")[[i]],"ID") = as.character(paste("ID",i,sep=""))} 

## Projecting the Grid
grid1 <- spTransform(grid1, CRS("+init=epsg:5641"))

# cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen", "forestgreen")
cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
colfunc <- colorRampPalette(cores)
grid.result$RLI.col = as.character(cut(grid.result$RLI, 
                                       # breaks = quantile(grid.result$RLI, na.rm=TRUE, prob = seq(0,1,0.05)), 
                                       breaks = seq(0,1,0.05), 
                                       labels = colfunc(20), include.lowest = TRUE))
grid.result$RLI.end.col = as.character(cut(grid.result$RLI.end, 
                                           # breaks = quantile(grid.result$RLI.end, na.rm=TRUE, prob = seq(0,1,0.05)), 
                                           breaks = seq(0,1,0.05), 
                                           labels = colfunc(20), include.lowest = TRUE))


#### EXPLORING THE OBSERVED RESULTS ####
par(mfrow=c(1,2))
plot(grid1, border=grey(0.9)) ## ENDEMIC SPECIES RICHNESS (SR)
legend("bottomright", "All species", bty='n',cex=2)
plot(grid1[!is.na(grid.result$RLI),],
     col=grid.result$RLI.col[!is.na(grid.result$RLI)],
     border=grid.result$RLI.col[!is.na(grid.result$RLI)],
     #border=grey(0.9),
     add=TRUE)
plot(af.proj, add=TRUE)
plot(grid1,border=grey(0.9)) ## OCASIONAL SPECIES RICHNESS (SR)
legend("bottomright", "Endemics",bty='n',cex=2)
plot(grid1[!is.na(grid.result$RLI.end),],
     col=grid.result$RLI.end.col[!is.na(grid.result$RLI.end)],
     border=grid.result$RLI.end.col[!is.na(grid.result$RLI.end)],
     #border=grey(0.9),
     add=TRUE)
plot(af.proj, add=TRUE)

#### CREATING THE GRID FOR PREDICTION ###
# Define the grid extent
bb = bbox(grid1)
bb[1] = bb[1] - 20000
bb[3] = bb[3] - 20000
bb[4] = bb[4] + 20000
x.range <- as.numeric(c(bb[1,1],bb[1,2]))  # min/max longitude of the interpolation area
y.range <- as.numeric(c(bb[2,1],bb[2,2]))  # min/max latitude of the interpolation area
rec_grid <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 5000), 
                        y = seq(from = y.range[1], to = y.range[2], by = 5000))  # expand points to grid
coordinates(rec_grid) <- ~x + y
#gridded(rec_grid) <- TRUE
proj4string(rec_grid) <- raster::crs(grid1)
rec_grid = as(rec_grid, "SpatialPoints")

#Selecting the grid cells intersecting the map
grid1.clip <- raster::crop(grid1, af.proj)
tmp = sp::over(rec_grid, grid1.clip)
# tmp = sp::over(rec_grid, grid1)
rec_grid1 = rec_grid[!is.na(tmp$ID),]
#rec_grid1 = rec_grid[gIntersects(rec_grid, grid1, byid = TRUE),]

#Creating the grid for interpolation
x.out = as.double(coordinates(rec_grid1)[,1])
y.out = as.double(coordinates(rec_grid1)[,2])
grid.out = cbind.data.frame(x.out, y.out)
names(grid.out)[1:2] = c("x","y")
coordinates(grid.out) = ~x+y
proj4string(grid.out) <- raster::crs(grid1)
#gridded(grid.out) <- TRUE

##Preparing data for interpolation
grid2 = grid1[!is.na(grid.result$select) & grid.result$select,]
x.data = as.double(coordinates(grid2)[,1])
y.data = as.double(coordinates(grid2)[,2])
dados = grid.result[!is.na(grid.result$select) & grid.result$select,]

dados$RLI_catCD <- apply(dados[,c("RLI_catC","RLI_catD")], 1, min)
dados$RLI_catCD.end <- apply(dados[,c("RLI_catC.end","RLI_catD.end")], 1, min)
dados$RLI_catCD.pa <- apply(dados[,c("RLI_catC.pa","RLI_catD.pa")], 1, min)

dados$Threat_catCD <- apply(dados[,c("Threat_catC","Threat_catD")], 1, max)
dados$Threat_catCD.end <- apply(dados[,c("Threat_catC.end","Threat_catD.end")], 1, max)
dados$Threat_catCD.pa <- apply(dados[,c("Threat_catC.pa","Threat_catD.pa")], 1, max)


grid2 = cbind.data.frame(x.data,y.data,dados)
names(grid2)[1:2] = c("x","y")
coordinates(grid2) = ~x+y
proj4string(grid2) <- raster::crs(grid1)

cor.af = "black" # "darkgreen" # 
cor.br = "grey"
cor.am = "black" 
cutoff = "25%"
am.lat.proj1 = raster::crop(am.lat.proj, bb)
brasil.proj1 = raster::crop(brasil.proj, bb)

#vars = c("N.total","S.total")
# vars = c("RLI","RLI.pa")
# vars = c("RLI_new","RLI_new.pa")
# vars = c("RLI", "RLI_catA", "RLI_catB", "RLI_catCD",
#          "RLI.pa", "RLI_catA.pa", "RLI_catB.pa", "RLI_catCD.pa")
# vars = c("Threat", "RLI.pa", "Threat.end", "RLI.end.pa")
# 
# nomes = c("A - Number of occurrences","B - Number of species")
# nomes = rep(c("A - All records","B - Presence only"),4)
# nomes = rep(c("A - All criteria","B - Criterion A","C - Criterion B", "D - Criteria C and D"),2)
# nomes = rep(c("A - All populations, % Threat",
#               "B - All populations, Red List Index",
#               "C - Endemic, % Threat", "D - Endemics, Red List Index"),2)
#Figure 3
vars = c("RLI.pa","RLI.end.pa")
nomes = rep(c("A - All populations","B - Endemics"),4)
figura <- "Figure3_"

#Figure SQ
vars = c("Threat.pa","Threat.end.pa")
nomes = rep(c("A - All populations","B - Endemics"),4)
figura <- "FigureSQ_"

#Figure # abortando...
# vars = c("Threat_catA.pa","Threat_catD.pa", "Threat_catC.pa","Threat_catD.pa")
# nomes = rep(c("A - Criterion A","B - Criterion B", "C - Criterion C", "D - Criterion D"),2)
# figura <- "FigureWW_"



for(i in 1:length(vars)) {
  var = vars[i]
  nome.fig = paste("figures/",figura, LETTERS[i],"_adpt_",var,".jpg",sep="")
  #nome.fig = paste("Figure2_",var,"_BW.jpg",sep="")
  nome.shp = gsub("\\.","_",paste("shp_centre_",var,sep=""))
  
  if (var == "N.total") {
    dados1 = grid.result[!is.na(grid.result[,var]),] 
    dados2 = log(dados1[,var]+1)
    
    grid2.1 = grid1
    x.data.1 = as.double(coordinates(grid2.1)[,1])
    y.data.1 = as.double(coordinates(grid2.1)[,2])
    grid2.1 = cbind.data.frame(x.data.1, y.data.1, grid.result)
    names(grid2.1)[1:2] = c("x","y")
    coordinates(grid2.1) = ~x+y
    proj4string(grid2.1) <- crs(grid1)
    grid3 = grid2.1[!is.na(grid2.1@data[,var]),]
  } else { 
    dados2 = dados[!is.na(dados[,var]), var]
    # if(grepl("CWE",var)) {
    #   dados2 = log(dados1[,var])
    #   infs = min(dados2[!is.infinite(dados2)]) - 0.05
    #   dados2[is.infinite(dados2)] = infs
    # } else {
    #   dados2 = log(dados1[,var]+1)
    # }
    grid3 = grid2[!is.na(grid2@data[,var]),]
  }
  #Inverse distance weighting
  #idw <- idw(dados2 ~ 1, locations= grid3, newdata= grid.out, idp =1, nmax=30)
  #idw.output = cbind.data.frame(coordinates(grid.out), idw$var1.pred)
  #names(idw.output)<-c("x","y","z")
  #par(mfrow=c(1,1), mar=c(2,2,0,0))
  # plot(af.proj, border="white")
  # image(idw.output, col = colfunc(21), add=TRUE,  breaks = c(0,quantile(dados2, na.rm=TRUE, prob = seq(0,1,0.05)) + seq(0,0.02,0.001)))
  # plot(af.proj, add=TRUE, border="grey")
  # idw.output1 = xtabs(z~x+y, idw.output)
  # cl <- contourLines(as.double(attributes(idw.output1)$dimnames$x),
  #             as.double(attributes(idw.output1)$dimnames$y),
  #             idw.output1, levels = quantile(dados2, prob= cutoff), nlevels=1)
  # shp <- ContourLines2SLDF(cl, proj4string=crs(grid1))
  # shp.idw = disaggregate(shp)
  # shp.idw1= shp.idw[gLength(shp.idw, byid = TRUE)>=mean(gLength(shp.idw, byid = TRUE)),]
  #plot(shp.idw, add=TRUE, lwd=2)
  
  ## Kriging
  # Fitting the semivariogram
  semivariog <- variogram(dados2 ~ 1, locations= grid3, data= grid3)
  
  #using fit.variogram
  rm(fvl, fve, fvs, fvg)
  #the data looks like it might be an exponential shape, so we will try that first with the values estimated from the empirical 
  fvl <- fit.variogram(semivariog, vgm(1, "Lin", 600000, 10))
  plot(semivariog, fvl)
  fve <- fit.variogram(semivariog, vgm(1, "Exp", 600000, 10000))
  plot(semivariog, fve)
  fvs <- fit.variogram(semivariog, vgm(1, "Sph", 600000, 10000))
  plot(semivariog, fvs)
  fvg <- fit.variogram(semivariog, vgm(1, "Gau", 600000, 10))
  plot(semivariog, fvg)
  #fvm <- fit.variogram(semivariog, vgm(1, "Mat", 600000, 10))
  #plot(semivariog, fvm)
  
  # Chosing the best model with the best fit for kriging
  SSErrs = c(fve = attr(fve, "SSErr"), 
             fvs  = attr(fvs, "SSErr"), 
             fvl  = attr(fvg, "SSErr"),
             fvg  = attr(fvg, "SSErr"))
  #fvm  = attr(fvm, "SSErr"))
  tmp = list(fve = fve, fvs = fvs, fvl = fvl, fvg = fvg)
  best.m = tmp[[head(names(SSErrs)[which(SSErrs == min(SSErrs))],1)]]
  
  # Kriging
  krig <- gstat::krige(formula= dados2 ~ 1, locations= grid3, newdata= grid.out, model= best.m, nmax = 25)
  #krig1<-krige(formula= dados2 ~ x+y, locations= grid3, newdata= grid.out, model= fit.variog.exp, nmax = 20)
  krig.output = cbind.data.frame(coordinates(grid.out), krig$var1.pred)
  names(krig.output)<-c("x","y","z")
  krig.output0 <- raster::rasterFromXYZ(krig.output)
  
  jpeg(filename = nome.fig, width = 3000*1.5, height = 3000*2, units = "px", pointsize = 12,
       res = 600, family = "sans", type="cairo", bg="white")
  
  #Definig the breaks for colors 
  qts = quantile(dados2, na.rm=TRUE, prob = seq(0,1,0.05))
  # qts = quantile(krig$var1.pred, na.rm=TRUE, prob = seq(0,1,0.05))
  zeros = sum(qts==0)
  uns = sum(qts==1)
  
  # if(grepl("CWE", var)) zeros = sum(qts==infs)
  cres = c(rep("#FFFFFF",zeros), colfunc(20-zeros))
  brks = qts
  if (zeros>0) {
    brks[brks==0] = seq(0,(zeros-1)*0.00001,0.00001)
    if(grepl("CWE", var)) brks[brks==infs] = seq(infs,infs+(zeros-1)*0.00001,0.00001)
  }
  if (uns>0) {
    brks[brks==1] = seq(1 - (uns-1)*0.00001, 1, 0.00001)
  }
  #ploting
  par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
  plot(af.proj, border="white")
  # raster::image(krig.output, col = cres, add=TRUE,  breaks = brks)
  if (grepl("Threat", var))
    cres <- rev(cres)
  # raster::image(krig.output0, col = cres, add=TRUE,  breaks = brks)
  
  # fixed breaks for scale
  brks1 <- seq(0.65, 1, 0.02)
  cores1 <- c("darkred", "red", "gold", "green")
  colfunc1 <- colorRampPalette(cores1)
  cres1 <- colfunc(length(brks1)-1)
  if (grepl("Threat", var))
    cres1 <- rev(cres1)
  raster::image(krig.output0, col = cres1, add=TRUE,  breaks = brks1)
  # krig.output1 = xtabs(z~x+y, krig.output)
  # cl <- contourLines(as.double(attributes(krig.output1)$dimnames$x),
  #                    as.double(attributes(krig.output1)$dimnames$y),
  #                    krig.output1, 
  #                    # levels = qts[which(names(qts) %in% c("75%","80%","85%","90%","95%"))])
  #                    levels = qts[which(names(qts) %in% c("5%","10%","15%","20%","25%"))])
  # shp <- maptools::ContourLines2SLDF(cl, proj4string = raster::crs(grid1))
  # shp@data$quantile = head(c("5%","10%","15%","20%","25%"), length(shp)) 
  # if(grepl("CWE",var)) {
  #   shp@data[,var] = NA
  #   shp@data[,var] = exp(as.double(as.character(shp@data$level)))
  # } else {
  #   shp@data[,var] = NA
  #   shp@data[,var] = as.double(as.character(shp@data$level))
  # }  
  # shp.krig.cut = sp::disaggregate(shp[which(shp@data$level == qts[which(names(qts)==cutoff)]),])
  # shp.krig.cut1 = shp.krig.cut[gLength(shp.krig.cut, byid = TRUE) >= mean(gLength(shp.krig.cut, byid = TRUE)),]
  # shp.krig.cut1 = shp.krig.cut[rgeos::gLength(shp.krig.cut, byid = TRUE) >=2*pi*25000,]
  plot(brasil.proj1, add=TRUE, lty="1414", border= cor.br)
  plot(af.proj, add=TRUE, border= cor.af)
  plot(am.lat.proj1, add=TRUE, border = cor.am)
  nome = nomes[i]
  legend("topleft", nome, cex=1.5, bty="n",inset=c(0.02,0.03))
  # plot(shp.krig.cut1, add=TRUE, lwd=2, col=1)
  par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
  at.leg <- seq(0,100,20)
  at.leg1 <- seq(60,100,10)
  vals.leg <- paste0(at.leg,"%")
  fields::image.plot( legend.only=TRUE, add=TRUE,
                      col = cres1,
                      legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
                      legend.line = 1, legend.cex = 0.5,
                      # zlim=c(2.5,97.5),
                      zlim=c(60,100),
                      #breaks = seq(0,100,20),
                      axis.args=list(at= at.leg1,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
                                     # labels = round(qts[names(qts) %in% vals.leg],2)))
                                     labels= at.leg1/100))
  if (grepl("Threat", var))  
    legend(4500000,6024198,expression(bold("Prop. Threatened Species")),cex=1.1,bty="n",y.intersp=-0.01)
  if (grepl("RLI", var))  
    legend(4900000,6024198,expression(bold("Red List Index")),cex=1.1,bty="n",y.intersp=-0.01)
  #colorbar.plot(5549396, 7757332, strip=cbind(0:100), horizontal=F, 
  #	    strip.width = 0.04, strip.length = 0.6,
  #      col=cres)
  #axis(2, at=seq(0,100,20), labels = TRUE, las=1)
  #plot(shp.idw1, add=TRUE, lwd=2, col=6)
  #plot(shp.krig1, add=TRUE, lwd=2, col=1)
  
  #Saving the shapefile
  # shp@data = shp@data[,2:3]
  #writeOGR(obj=shp, dsn='.', layer=nome.shp, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## Piecewise linear interpolation
  # IL = interp(grid3, z= var, xo=x.data, yo=y.data,  
  #             nx=100, ny= 100)
  #plot(af.proj, border="white")
  # image(IL, add=TRUE)
  # x = coordinates(IL)[,1]; y = coordinates(IL)[,2]; z = IL@data[,var]
  # obj = xtabs(z ~ x + y)
  # cl <- contourLines(as.double(attributes(obj)$dimnames$x),
  #                    as.double(attributes(obj)$dimnames$y),
  #                    obj, levels = quantile(log(dados+1), prob= cutoff), nlevels=1)
  # shp <- ContourLines2SLDF(cl, proj4string=crs(grid1))
  # shp.lin = disaggregate(shp)
  # shp.lin1 = shp.lin[gLength(shp.lin, byid = TRUE)>mean(gLength(shp.lin, byid = TRUE)),]
  # plot(shp.lin1, add=TRUE, lwd=2, col=5)
  
  ## Bilinear interpolation
  #obj = xtabs(dados1~x.data + y.data)
  #obj = list(as.double(attributes(obj)$dimnames$x.data),
  #           as.double(attributes(obj)$dimnames$y.data),
  #           obj)
  #interp.surface.grid(obj, coordinates(grid.out))
  dev.off()
}




# -------------------------------------------------------------------------
##############################
#### 50 KM HEXAGONAL GRID ####
##############################
## reading the results per grid cell
grid.result = readRDS("data/grid.results_50km.rds")[[1]]
samp.cov.cutoff = readRDS("data/grid.results_50km.rds")[[2]]

## Defining the sampling coverage cut-off for analysis
#For rarefaction, a temptative minimum sampling is:
plot(RLI ~ N.total, grid.result, log="x"); abline(v=c(30,50,100))
plot(RLI.end ~ N.total, grid.result, log="xy"); abline(v=c(30,50,100))
plot(RLI ~ S.total, grid.result, log="x"); abline(v=c(30,50,100))
plot(RLI.end ~ S.total, grid.result, log="xy"); abline(v=c(30,50,100))

#What does it mean to remove cells above a given quantile of the sampling coverage distribution?
result = NULL; result1 = NULL
for (i in 1:length(samp.cov.cutoff)) {
  toto = samp.cov.cutoff[i]
  toto1= summary(grid.result$N.total[grid.result$SampCover < toto])
  toto2= summary(grid.result$N.total[grid.result$SampCover >= toto])
  toto3= c(toto1[c(3,4,1,6)],toto2[c(3,4,1,6)]) 
  toto4= summary(grid.result$S.total[grid.result$SampCover < toto])
  toto5= summary(grid.result$S.total[grid.result$SampCover >= toto])
  toto6= c(toto4[c(3,4,1,6)],toto5[c(3,4,1,6)]) 
  if(is.null(result)) result = toto3 else result = rbind.data.frame(result, toto3)
  if(is.null(result1)) result1 = toto6 else result1 = rbind.data.frame(result1, toto6)
}
row.names(result) = row.names(result1) = names(samp.cov.cutoff)
colnames(result) = colnames(result1) = c("med.excl","mean.excl","min.excl","max.excl","med.incl","mean.incl","min.incl","max.incl")
result
result1
cut = samp.cov.cutoff[4] # using the 50% quantile (mean) of the sampling coverage distribution
head(sort(grid.result$N.total[grid.result$SampCover>cut]), 50) ## REMOVE ALSO CELLS WITH LESS THAN < 30 OCCURRENCES

#Inspecting other possibilities of data filtering
ns <- c(50,100,250)
spp <- c(30,50,100)
grid.result$select = ((grid.result$SampCover >= samp.cov.cutoff[5]  & grid.result$N.total >= ns[1]) | 
                        (grid.result$SampCover >= samp.cov.cutoff[4] & grid.result$SampCover < samp.cov.cutoff[5] & grid.result$N.total >= ns[2]) | 
                        (grid.result$SampCover >= samp.cov.cutoff[3] & grid.result$SampCover < samp.cov.cutoff[4] & grid.result$N.total >= ns[3]))
grid.result$select1 = ((grid.result$SampCover >= samp.cov.cutoff[5]  & grid.result$S.total >= spp[1]) | 
                        (grid.result$SampCover >= samp.cov.cutoff[4] & grid.result$SampCover < samp.cov.cutoff[5] & grid.result$S.total >= spp[2]) | 
                        (grid.result$SampCover >= samp.cov.cutoff[3] & grid.result$SampCover < samp.cov.cutoff[4] & grid.result$S.total >= spp[3]))
#Result by N.total
par(mfrow = c(1,2))
plot(grid.result$SampCover ~ grid.result$N.total, log="x")
abline(v= ns, col=c(1,2,3)); abline(h=samp.cov.cutoff[3:5], col=c(3,2,1))
points(SampCover ~ N.total, subset(grid.result,select), col=2, pch=21)
summary(grid.result$N.total[grid.result$select])
  head(sort(grid.result$N.total[grid.result$select]),50)
summary(grid.result$SampCover[grid.result$select])
#Result by S.total
plot(grid.result$SampCover ~ grid.result$S.total, log="x")
abline(v= spp, col=c(1,2,3)); abline(h=samp.cov.cutoff[3:5], col=c(3,2,1))
points(SampCover ~ S.total, subset(grid.result,select1), col=2, pch=21)
summary(grid.result$S.total[grid.result$select1])
head(sort(grid.result$S.total[grid.result$select1]),50)
summary(grid.result$SampCover[grid.result$select1])


############################################################H
#### CENTRES OF HIGH CONCENTRATION OF THREATENED SPECIES ####
############################################################H
## Getting the grid and editing it
grid <- readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//af_hex_grid_50km.rds")

## Making sure grid IDs are the same internally
grid1 <- grid
tmp <- SpatialPoints(coords = coordinates(grid), proj4string = raster::crs(grid))
tmp1 <- over(tmp, grid)
for(i in 1:length(grid1)) { slot(slot(grid1, "polygons")[[i]],"ID") = as.character(paste("ID",i,sep=""))} 

## Projecting the Grid
grid1 <- spTransform(grid1, CRS("+init=epsg:5641"))

# cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen", "forestgreen")
cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
colfunc <- colorRampPalette(cores)
grid.result$RLI.col = as.character(cut(grid.result$RLI, 
                                       # breaks = quantile(grid.result$RLI, na.rm=TRUE, prob = seq(0,1,0.05)), 
                                       breaks = seq(0,1,0.05), 
                                       labels = colfunc(20), include.lowest = TRUE))
grid.result$RLI.end.col = as.character(cut(grid.result$RLI.end, 
                                         # breaks = quantile(grid.result$RLI.end, na.rm=TRUE, prob = seq(0,1,0.05)), 
                                         breaks = seq(0,1,0.05), 
                                         labels = colfunc(20), include.lowest = TRUE))


#### EXPLORING THE OBSERVED RESULTS ####
par(mfrow=c(1,2))
plot(grid1, border=grey(0.9)) ## ENDEMIC SPECIES RICHNESS (SR)
legend("bottomright", "All species", bty='n',cex=2)
plot(grid1[!is.na(grid.result$RLI)],
     col=grid.result$RLI.col[!is.na(grid.result$RLI)],
     border=grid.result$RLI.col[!is.na(grid.result$RLI)],
     #border=grey(0.9),
     add=TRUE)
plot(af.proj, add=TRUE)
plot(grid1,border=grey(0.9)) ## OCASIONAL SPECIES RICHNESS (SR)
legend("bottomright", "Endemics",bty='n',cex=2)
plot(grid1[!is.na(grid.result$RLI.end)],
     col=grid.result$RLI.end.col[!is.na(grid.result$RLI.end)],
     border=grid.result$RLI.end.col[!is.na(grid.result$RLI.end)],
     #border=grey(0.9),
     add=TRUE)
plot(af.proj, add=TRUE)

#### CREATING THE GRID FOR PREDICTION ###
# Define the grid extent
bb = bbox(grid1)
bb[1] = bb[1] - 20000
bb[3] = bb[3] - 20000
bb[4] = bb[4] + 20000
x.range <- as.numeric(c(bb[1,1],bb[1,2]))  # min/max longitude of the interpolation area
y.range <- as.numeric(c(bb[2,1],bb[2,2]))  # min/max latitude of the interpolation area
rec_grid <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 5000), 
                        y = seq(from = y.range[1], to = y.range[2], by = 5000))  # expand points to grid
coordinates(rec_grid) <- ~x + y
#gridded(rec_grid) <- TRUE
proj4string(rec_grid) <- raster::crs(grid1)
rec_grid = as(rec_grid, "SpatialPoints")

#Selecting the grid cells intersecting the map
tmp = sp::over(rec_grid, grid1)
rec_grid1 = rec_grid[!is.na(tmp),]
#rec_grid1 = rec_grid[gIntersects(rec_grid, grid1, byid = TRUE),]

#Creating the grid for interpolation
x.out = as.double(coordinates(rec_grid1)[,1])
y.out = as.double(coordinates(rec_grid1)[,2])
grid.out = cbind.data.frame(x.out, y.out)
names(grid.out)[1:2] = c("x","y")
coordinates(grid.out) = ~x+y
proj4string(grid.out) <- raster::crs(grid1)
#gridded(grid.out) <- TRUE

##Preparing data for interpolation
grid2 = grid1[!is.na(grid.result$select1) & grid.result$select1,]
x.data = as.double(coordinates(grid2)[,1])
y.data = as.double(coordinates(grid2)[,2])
dados = grid.result[!is.na(grid.result$select1) & grid.result$select1,]

dados$RLI_catCD <- apply(dados[,c("RLI_catC","RLI_catD")], 1, min)
dados$RLI_catCD.end <- apply(dados[,c("RLI_catC.end","RLI_catD.end")], 1, min)

grid2 = cbind.data.frame(x.data,y.data,dados)
names(grid2)[1:2] = c("x","y")
coordinates(grid2) = ~x+y
proj4string(grid2) <- raster::crs(grid1)

cor.af = "black" # "darkgreen" # 
cor.br = "grey"
cor.am = "black" 
cutoff = "25%"
#vars = c("N.total","S.total")
vars = c("RLI","RLI.end")
vars = c("RLI.end", "RLI_catA.end", "RLI_catB.end", "RLI_catCD.end")


#nomes = c("A - Number of occurrences","B - Number of species")
nomes = rep(c("A - All species","B - Endemics only"),4)
nomes = rep(c("A - All criteria","B - Criterion A","C - Criterion B", "D - Criteria C and D"),2)
am.lat.proj1 = raster::crop(am.lat.proj, bb)
brasil.proj1 = raster::crop(brasil.proj, bb)

for(i in 1:length(vars)) {
  var = vars[i]
  nome.fig = paste("figures/Figure3_hex_",var,".jpg",sep="")
  #nome.fig = paste("Figure2_",var,"_BW.jpg",sep="")
  nome.shp = gsub("\\.","_",paste("shp_centre_",var,sep=""))
  
  if (var == "N.total") {
    dados1 = grid.result[!is.na(grid.result[,var]),] 
    dados2 = log(dados1[,var]+1)
    
    grid2.1 = grid1
    x.data.1 = as.double(coordinates(grid2.1)[,1])
    y.data.1 = as.double(coordinates(grid2.1)[,2])
    grid2.1 = cbind.data.frame(x.data.1, y.data.1, grid.result)
    names(grid2.1)[1:2] = c("x","y")
    coordinates(grid2.1) = ~x+y
    proj4string(grid2.1) <- crs(grid1)
    grid3 = grid2.1[!is.na(grid2.1@data[,var]),]
  } else { 
    dados2 = dados[!is.na(dados[,var]), var]
    # if(grepl("CWE",var)) {
    #   dados2 = log(dados1[,var])
    #   infs = min(dados2[!is.infinite(dados2)]) - 0.05
    #   dados2[is.infinite(dados2)] = infs
    # } else {
    #   dados2 = log(dados1[,var]+1)
    # }
    grid3 = grid2[!is.na(grid2@data[,var]),]
  }
  #Inverse distance weighting
  #idw <- idw(dados2 ~ 1, locations= grid3, newdata= grid.out, idp =1, nmax=30)
  #idw.output = cbind.data.frame(coordinates(grid.out), idw$var1.pred)
  #names(idw.output)<-c("x","y","z")
  #par(mfrow=c(1,1), mar=c(2,2,0,0))
  # plot(af.proj, border="white")
  # image(idw.output, col = colfunc(21), add=TRUE,  breaks = c(0,quantile(dados2, na.rm=TRUE, prob = seq(0,1,0.05)) + seq(0,0.02,0.001)))
  # plot(af.proj, add=TRUE, border="grey")
  # idw.output1 = xtabs(z~x+y, idw.output)
  # cl <- contourLines(as.double(attributes(idw.output1)$dimnames$x),
  #             as.double(attributes(idw.output1)$dimnames$y),
  #             idw.output1, levels = quantile(dados2, prob= cutoff), nlevels=1)
  # shp <- ContourLines2SLDF(cl, proj4string=crs(grid1))
  # shp.idw = disaggregate(shp)
  # shp.idw1= shp.idw[gLength(shp.idw, byid = TRUE)>=mean(gLength(shp.idw, byid = TRUE)),]
  #plot(shp.idw, add=TRUE, lwd=2)
  
  ## Kriging
  # Fitting the semivariogram
  semivariog <- variogram(dados2 ~ 1, locations= grid3, data= grid3)
  
  #using fit.variogram
  rm(fvl, fve, fvs, fvg)
  #the data looks like it might be an exponential shape, so we will try that first with the values estimated from the empirical 
  fvl <- fit.variogram(semivariog, vgm(1, "Lin", 600000, 10))
  plot(semivariog, fvl)
  fve <- fit.variogram(semivariog, vgm(1, "Exp", 600000, 10000))
  plot(semivariog, fve)
  fvs <- fit.variogram(semivariog, vgm(1, "Sph", 600000, 10000))
  plot(semivariog, fvs)
  fvg <- fit.variogram(semivariog, vgm(1, "Gau", 600000, 10))
  plot(semivariog, fvg)
  #fvm <- fit.variogram(semivariog, vgm(1, "Mat", 600000, 10))
  #plot(semivariog, fvm)
  
  # Chosing the best model with the best fit for kriging
  SSErrs = c(fve = attr(fve, "SSErr"), 
             fvs  = attr(fvs, "SSErr"), 
             fvl  = attr(fvg, "SSErr"),
             fvg  = attr(fvg, "SSErr"))
  #fvm  = attr(fvm, "SSErr"))
  tmp = list(fve = fve, fvs = fvs, fvl = fvl, fvg = fvg)
  best.m = tmp[[head(names(SSErrs)[which(SSErrs == min(SSErrs))],1)]]
  
  # Kriging
  krig <- gstat::krige(formula= dados2 ~ 1, locations= grid3, newdata= grid.out, model= best.m, nmax = 25)
  #krig1<-krige(formula= dados2 ~ x+y, locations= grid3, newdata= grid.out, model= fit.variog.exp, nmax = 20)
  krig.output = cbind.data.frame(coordinates(grid.out), krig$var1.pred)
  names(krig.output)<-c("x","y","z")
  krig.output0 <- raster::rasterFromXYZ(krig.output)

  jpeg(filename = nome.fig, width = 3000*1.5, height = 3000*2, units = "px", pointsize = 12,
       res = 600, family = "sans", type="cairo", bg="white")
  
  #Defining the breaks for colors 
  qts = quantile(krig$var1.pred, na.rm=TRUE, prob = seq(0,1,0.05))
  # qts = quantile(dados2, na.rm=TRUE, prob = seq(0,1,0.05))
  zeros = sum(qts==0)
  # if(grepl("CWE", var)) zeros = sum(qts==infs)
  cres = c(rep("#FFFFFF",zeros), colfunc(20-zeros))
  brks = qts
  if (zeros>0) {
    brks[brks==0] = seq(0,(zeros-1)*0.00001,0.00001)
    if(grepl("CWE", var)) brks[brks==infs] = seq(infs,infs+(zeros-1)*0.00001,0.00001)
  }
  #ploting
  par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
  plot(af.proj, border="white")
  # raster::image(krig.output, col = cres, add=TRUE,  breaks = brks)
  raster::image(krig.output0, col = cres, add=TRUE,  breaks = brks)
  
  #Defining the breaks for colors 
  brks1 <- seq(0, 1, 0.025)
  # cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
  cores1 <- c("darkred", "red", "gold", "green")
  colfunc1 <- colorRampPalette(cores1)
  cres1 <- colfunc1(length(brks1)-1)
  raster::image(krig.output0, col = cres1, add=TRUE,  breaks = brks1)
  # krig.output1 = xtabs(z~x+y, krig.output)
  # cl <- contourLines(as.double(attributes(krig.output1)$dimnames$x),
  #                    as.double(attributes(krig.output1)$dimnames$y),
  #                    krig.output1, 
  #                    # levels = qts[which(names(qts) %in% c("75%","80%","85%","90%","95%"))])
  #                    levels = qts[which(names(qts) %in% c("5%","10%","15%","20%","25%"))])
  # shp <- maptools::ContourLines2SLDF(cl, proj4string = raster::crs(grid1))
  # shp@data$quantile = head(c("5%","10%","15%","20%","25%"), length(shp)) 
  # if(grepl("CWE",var)) {
  #   shp@data[,var] = NA
  #   shp@data[,var] = exp(as.double(as.character(shp@data$level)))
  # } else {
  #   shp@data[,var] = NA
  #   shp@data[,var] = as.double(as.character(shp@data$level))
  # }  
  # shp.krig.cut = sp::disaggregate(shp[which(shp@data$level == qts[which(names(qts)==cutoff)]),])
  # shp.krig.cut1 = shp.krig.cut[gLength(shp.krig.cut, byid = TRUE) >= mean(gLength(shp.krig.cut, byid = TRUE)),]
  # shp.krig.cut1 = shp.krig.cut[rgeos::gLength(shp.krig.cut, byid = TRUE) >=2*pi*25000,]
  plot(brasil.proj1, add=TRUE, lty="1414", border= cor.br)
  plot(af.proj, add=TRUE, border= cor.af)
  plot(am.lat.proj1, add=TRUE, border = cor.am)
  nome = nomes[i]
  legend("topleft", nome, cex=1.5, bty="n",inset=c(0.02,0.03))
  # plot(shp.krig.cut1, add=TRUE, lwd=2, col=1)
  par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
  at.leg <- seq(0,100,10)
  at.leg1 <- seq(10,100,20)
  vals.leg <- paste0(at.leg,"%")
  fields::image.plot( legend.only=TRUE, add=TRUE,
                      # col = cres,
                      col = cres1,
                      legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
                      legend.line = 1, legend.cex = 0.5,
                      # zlim=c(2.5,97.5),
                      zlim=c(5,95),
                      #breaks = seq(0,100,20),
                      axis.args=list(at= at.leg1,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
                             # labels = round(qts[names(qts) %in% vals.leg],2)))
                             labels= at.leg1/100))
  legend(4900000,6024198,expression(bold("Red List Index")),cex=1.1,bty="n",y.intersp=-0.01)
  #colorbar.plot(5549396, 7757332, strip=cbind(0:100), horizontal=F, 
  #	    strip.width = 0.04, strip.length = 0.6,
  #      col=cres)
  #axis(2, at=seq(0,100,20), labels = TRUE, las=1)
  #plot(shp.idw1, add=TRUE, lwd=2, col=6)
  #plot(shp.krig1, add=TRUE, lwd=2, col=1)
  
  #Saving the shapefile
  # shp@data = shp@data[,2:3]
  #writeOGR(obj=shp, dsn='.', layer=nome.shp, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  ## Piecewise linear interpolation
  # IL = interp(grid3, z= var, xo=x.data, yo=y.data,  
  #             nx=100, ny= 100)
  #plot(af.proj, border="white")
  # image(IL, add=TRUE)
  # x = coordinates(IL)[,1]; y = coordinates(IL)[,2]; z = IL@data[,var]
  # obj = xtabs(z ~ x + y)
  # cl <- contourLines(as.double(attributes(obj)$dimnames$x),
  #                    as.double(attributes(obj)$dimnames$y),
  #                    obj, levels = quantile(log(dados+1), prob= cutoff), nlevels=1)
  # shp <- ContourLines2SLDF(cl, proj4string=crs(grid1))
  # shp.lin = disaggregate(shp)
  # shp.lin1 = shp.lin[gLength(shp.lin, byid = TRUE)>mean(gLength(shp.lin, byid = TRUE)),]
  # plot(shp.lin1, add=TRUE, lwd=2, col=5)
  
  ## Bilinear interpolation
  #obj = xtabs(dados1~x.data + y.data)
  #obj = list(as.double(attributes(obj)$dimnames$x.data),
  #           as.double(attributes(obj)$dimnames$y.data),
  #           obj)
  #interp.surface.grid(obj, coordinates(grid.out))
  dev.off()
}




