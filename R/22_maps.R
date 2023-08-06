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
require(rgdal)
require(rgeos)
require(cleangeo)
require(raster)
require(fields)
require(sp)
require(gstat)

## Getting only the necessau shapefiles
am.lat <- rgdal::readOGR(dsn="data/Am_Lat_ADM", layer="LatinAmerica")
brasil <- readRDS("data/Am_Lat_ADM/gadm36_BRA_1_sp.rds")

af <- rgdal::readOGR(dsn="data/AF_limits/", layer="merge_limites_MA11428_TNC_ARG_BR_PAR")
af <- rgeos::gBuffer(af, byid=TRUE, width=0) #correcting possible overlapping polygons
af <- cleangeo::clgeo_Clean(af) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
af <- raster::aggregate(af, dissolve=TRUE)
#projecting
am.lat.proj= spTransform(am.lat,  CRS("+init=epsg:5641"))
brasil.proj= spTransform(brasil,  CRS("+init=epsg:5641"))
af.proj = spTransform(af, CRS("+init=epsg:5641"))



# -------------------------------------------------------------------------
##################################
#### ADAPTIVE RESOLUTION GRID ####
##################################
## reading the results per grid cell
# toto <- readRDS("data/grid.results_adaptive_resol_old.rds")
grid.result <- readRDS("data/grid.results_adaptive_resol.rds")#[[1]]
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
plot(RLI_erro_pad ~ N.total, grid.result, log="x"); abline(v=c(30,50,100))

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
# grid.result$select[is.na(grid.result$select)] <- FALSE
# grid.result$select1[is.na(grid.result$select1)] <- FALSE

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
summary(grid.result$SampCover.end[grid.result$select1])

## For the SOM
table(grid.result$select1, useNA = "always")
summary(grid.result$N.total[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$S.total[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$RLI.pa[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$Threat.pa[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$Threat.end.pa[!is.na(grid.result$select) & grid.result$select])
summary(grid.result$SampCover[!is.na(grid.result$select) & grid.result$select])

############################################################H
#### CENTRES OF HIGH CONCENTRATION OF THREATENED SPECIES ####
############################################################H
# ## Getting the grid and editing it
# grid0 <- readRDS("./data/adaptative_resolution_AF_min200_max500_cellmax2/clean_grid.rds")
# 
# ## removing overllaping parts of the grid cells
# ids.polys <- grid0$ID
# grid_sf <- sf::st_as_sf(grid0)
# grid_inter <- sf::st_intersection(grid_sf)
# # x_inter <- x_inter[x_inter$n.overlaps == 1,]
# grid_inter <- grid_inter[rownames(grid_inter) %in% ids.polys,]
# grid <- as(grid_inter, "Spatial")
# 
# ## Making sure grid IDs are the same internally
# grid1 <- grid
# tmp <- SpatialPoints(coords = coordinates(grid), proj4string = raster::crs(grid))
# tmp1 <- over(tmp, grid)
# for(i in 1:length(grid1)) { slot(slot(grid1, "polygons")[[i]],"ID") = as.character(paste("ID",i,sep=""))} 

## Projecting the Grid
grid <- readRDS("./data/adaptative_resolution_AF_min200_max500_cellmax2/clean_grid_final.rds")
grid1 <- spTransform(grid, CRS("+init=epsg:5641"))

#### EXPLORING THE OBSERVED RESULTS ####
## Defining the color scales
# cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen", "forestgreen")
brks <- seq(0,1,0.05)
cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
colfunc <- colorRampPalette(cores)

brks1 <- seq(0, 0.35, 0.025)
cores1 <- c("darkred","red", "gold", "green")
colfunc1 <- colorRampPalette(cores1)
cres1 <- rev(colfunc1(length(brks1)-1))

# brks.threat <- seq(0.28, 1, 0.0375)
brks.threat <- seq(0, 1, 0.05)
cores.threat <- c("darkred", "red", "gold", "green")
colfunc.threat <- colorRampPalette(cores.threat)
cres1.threat <- rev(colfunc.threat(length(brks.threat)-1))

cores2 <- c("darkred","red", "lightpink2","cornflowerblue", "blue", "blue4")
colfunc2 <- colorRampPalette(cores2)
brks2 <- seq(0,0.21,0.01)
grid.result$RLI.col = as.character(cut(grid.result$RLI, 
         # breaks = quantile(grid.result$RLI, na.rm=TRUE, prob = seq(0,1,0.05)), 
         breaks = brks, 
         labels = colfunc(20), include.lowest = TRUE))
grid.result$RLI.end.col = as.character(cut(grid.result$RLI.end, 
             # breaks = quantile(grid.result$RLI.end, na.rm=TRUE, prob = seq(0,1,0.05)), 
             breaks = brks, 
             labels = colfunc(20), include.lowest = TRUE))
grid.result$Threat.col = as.character(cut(grid.result$Threat.pa, 
            breaks = brks,  
            labels = rev(colfunc(length(brks)-1)), include.lowest = TRUE))
grid.result$Threat.end.col = as.character(cut(grid.result$Threat.end.pa, 
                breaks = brks,  
                labels = rev(colfunc(length(brks)-1)), include.lowest = TRUE))
brks.N <- seq(0.25, 4.25, 0.2)
grid.result$N.col = as.character(cut(log10(grid.result$N.total), 
       breaks = brks.N,  
       labels = rev(colfunc(length(brks.N)-1)), include.lowest = TRUE))
brks.S <- seq(0.25, 3.6, 0.2)
grid.result$S.col = as.character(cut(log10(grid.result$S.total), 
       breaks = brks.S,  
       labels = rev(colfunc(length(brks.S)-1)), include.lowest = TRUE))
grid.result$RLI.error.col = as.character(cut(grid.result$RLI_erro_pad, 
               # breaks = quantile(grid.result$RLI, na.rm=TRUE, prob = seq(0,1,0.05)), 
               breaks = brks2, 
               labels = rev(colfunc2(length(brks2)-1)), include.lowest = TRUE))
grid.result$RLI.error.end.col = as.character(cut(grid.result$RLI_erro_pad.end, 
                   # breaks = quantile(grid.result$RLI.end, na.rm=TRUE, prob = seq(0,1,0.05)), 
                   breaks = brks2, 
                   labels = rev(colfunc2(length(brks2)-1)), include.lowest = TRUE))

# Grids with no data -> "White" ()
grid.result$Threat.col[is.na(grid.result$select)] <- "#FFFFFF"
#grid.result$Threat.col[is.na(grid.result$Threat.col)] <- "#FFFFFF"
grid.result$Threat.end.col[is.na(grid.result$select)] <- "#FFFFFF"
#grid.result$Threat.end.col[is.na(grid.result$Threat.end.col)] <- "#FFFFFF"

# Grids with insuficient data -> "Light Grey" ("#D3D3D3") or "Grey" ("#808080") or "Dark Grey" ("#A9A9A9")
grid.result$Threat.col[!is.na(grid.result$select) & 
                         !grid.result$select] <- "#D3D3D3"
grid.result$Threat.end.col[!is.na(grid.result$select) & 
                             !grid.result$select] <- "#D3D3D3" 
cor.af = "black" # "darkgreen" # 
cor.br = "grey"
cor.am = "black" 
          
#Plot extent plus cropping
bb = bbox(grid1)
bb[1] = bb[1] - 20000
bb[3] = bb[3] - 20000
bb[4] = bb[4] + 20000
am.lat.proj1 = raster::crop(am.lat.proj, bb)
brasil.proj1 = raster::crop(brasil.proj, bb)
        
        
###################
#### FIGURE SR ####
###################
at.leg <- seq(00,100,20)
at.leg1 <- seq(10,100,10)
vals.leg <- paste0(at.leg,"%")
        
#### REVER ESSAS FIGURA ####
        
#PANEL A
jpeg(filename = "figures/Figure_SRA.jpg", width = 3000*1.5, height = 3000*2, 
             units = "px", pointsize = 12,
             res = 600, family = "sans", type="cairo", bg="white")
        
par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
plot(af.proj, border="white")
plot(grid1, border=grey(0.75), add = TRUE)
plot(brasil.proj1, add=TRUE, 
  #lty="1414", 
 border= cor.br)
plot(am.lat.proj1, add=TRUE, border = cor.am)
plot(grid1, col = adjustcolor(grid.result$Threat.col,
  alpha.f = 0.75,
  offset = c(0.1, 0.1, 0.1, 0.15)), # <- "more white"
         #transform = diag(c(1,1,1,0.5))),
         #border = "white",
         add=TRUE)
# plot(af.proj, add=TRUE, border= cor.af)
legend("topleft", "A - All populations", cex=1.5, bty="n",inset=c(0.02,0.03))
        
#at.leg <- c(0,1,2,3,4)
par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
fields::image.plot( legend.only=TRUE, add=TRUE,
                            col = adjustcolor(cres1.threat, alpha.f = 0.75),
                            legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
                            legend.line = 1, legend.cex = 0.5,
                            zlim=c(0,100),
                            axis.args=list(at= at.leg,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
             labels= at.leg/100))
legend(4500000,6024198,expression(bold("Prop. Threatened Species")),
 cex=1.1,bty="n",y.intersp=-0.01)
dev.off()
        
#PANEL B
jpeg(filename = "figures/Figure_SRB.jpg", width = 3000*1.5, height = 3000*2, 
 units = "px", pointsize = 12,
 res = 600, family = "sans", type="cairo", bg="white")
        
par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
plot(af.proj, border="white")
plot(grid1, border=grey(0.75), add = TRUE)
plot(brasil.proj1, add=TRUE, 
             #lty="1414", 
             border= cor.br)
plot(am.lat.proj1, add=TRUE, border = cor.am)
        
rm.cells1 <- !grid.result$Threat.end.col %in% c("#FFFFFF", "#D3D3D3")
plot(grid1[rm.cells1, ],
             col = adjustcolor(grid.result$Threat.end.col[rm.cells1],
 alpha.f = 0.75,
 offset = c(0.1, 0.1, 0.1, 0.15)), # <- "more white"
             # transform = diag(c(1,1,1,0.5))),
             border = "black",
             add=TRUE)
legend("topleft", "B - Endemics", cex=1.5, bty="n",inset=c(0.02,0.03))
        
par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
fields::image.plot( legend.only=TRUE, add=TRUE,
                      col = adjustcolor(cres1, alpha.f = 0.75),
                      legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
                      legend.line = 1, legend.cex = 0.5,
                      zlim=c(10,100),
                        axis.args=list(at= at.leg,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
                                           labels= at.leg/100))
legend(4500000,6024198,expression(bold("Prop. Threatened Species")),
               cex=1.1,bty="n",y.intersp=-0.01)
dev.off()
        
        
###################
#### FIGURE UU ####
###################
# LOG10 of the number of occurrences and species per grid cell
        
#PANEL A
jpeg(filename = "figures/Figure_UU_A.jpg", width = 3000*1.5, height = 3000*2, 
  units = "px", pointsize = 12,
  res = 600, family = "sans", type="cairo", bg="white")
        
par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
plot(af.proj, border="white")
plot(grid1, border=grey(0.75), add = TRUE)
plot(brasil.proj1, add=TRUE, 
    #lty="1414", 
    border= cor.br)
plot(am.lat.proj1, add=TRUE, border = cor.am)
        
# rm.cells <- !grid.result$Threat.col %in% c("#FFFFFF", "#D3D3D3")
plot(grid1, col = adjustcolor(grid.result$N.col,
        alpha.f = 0.75,
        offset = c(0.1, 0.1, 0.1, 0.15)), # <- "more white"
        # transform = diag(c(1,1,1,0.5))),
        border = "black",
        add=TRUE)
# plot(af.proj, add=TRUE, border= cor.af)
legend("topleft", "A - Number of occurrences", cex=1.5, bty="n",inset=c(0.02,0.03))
        
cres1.N <- rev(colfunc(length(brks.N)-1))
at.leg <- c(0,1,2,3,4)
par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
fields::image.plot( legend.only=TRUE, add=TRUE,
  col = adjustcolor(cres1.N, alpha.f = 0.75),
  legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
  legend.line = 1, legend.cex = 0.5,
  zlim=c(0,5),
  axis.args=list(at= at.leg,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
    labels= c(1,10,100,1000,10000)))
# legend(4500000,6024198,expression(bold("Prop. Threatened Species")),
#        cex=1.1,bty="n",y.intersp=-0.01)
dev.off()
        
#PANEL B
jpeg(filename = "figures/Figure_UU_B.jpg", width = 3000*1.5, height = 3000*2, 
 units = "px", pointsize = 12,
 res = 600, family = "sans", type="cairo", bg="white")
        
par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
plot(af.proj, border="white")
plot(grid1, border=grey(0.75), add = TRUE)
plot(brasil.proj1, add=TRUE, 
    #lty="1414", 
    border= cor.br)
plot(am.lat.proj1, add=TRUE, border = cor.am)
        
plot(grid1, col = adjustcolor(grid.result$S.col,
      alpha.f = 0.75,
      offset = c(0.1, 0.1, 0.1, 0.15)), # <- "more white"
      # transform = diag(c(1,1,1,0.5))),
      border = "black",
      add=TRUE)
legend("topleft", "B - Number of species", cex=1.5, bty="n",inset=c(0.02,0.03))
        
cres1.S <- rev(colfunc(length(brks.S)-1))
at.leg <- c(0,1,2,3)
        
par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
fields::image.plot( legend.only=TRUE, add=TRUE,
        col = adjustcolor(cres1.S, alpha.f = 0.75),
        legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
        legend.line = 1, legend.cex = 0.5,
        zlim=c(0,4),
        axis.args=list(at= at.leg,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
                                           labels= c(1,10,100,1000)))
# legend(4500000,6024198,expression(bold("Prop. Threatened Species")),
# cex=1.1,bty="n",y.intersp=-0.01)
dev.off()
        
# -------------------------------------------------------------------------
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
#  "RLI.pa", "RLI_catA.pa", "RLI_catB.pa", "RLI_catCD.pa")
# vars = c("Threat", "RLI.pa", "Threat.end", "RLI.end.pa")
# 
# nomes = c("A - Number of occurrences","B - Number of species")
# nomes = rep(c("A - All records","B - Presence only"),4)
# nomes = rep(c("A - All criteria","B - Criterion A","C - Criterion B", "D - Criteria C and D"),2)
# nomes = rep(c("A - All populations, % Threat",
#       "B - All populations, Red List Index",
#               "C - Endemic, % Threat", "D - Endemics, Red List Index"),2)

#Figure 3
vars = c("RLI.pa","RLI.end.pa")
nomes = rep(c("A - All species","B - Endemics"),4)
figura <- "Figure3_"
              
#Figure SQ
vars = c("Threat.pa","Threat.end.pa")
nomes = rep(c("A - All species","B - Endemics"),4)
figura <- "FigureSQ_"
              
#Figure SS (Aqui só os painéis C e D - A e B são gerados abaixo no item DISTRIBUTIONS AND INTERPOLATION FOR CR SPECIES)
vars = c("Threat.pa_CR", "Threat.end.pa_CR")
nomes = rep(c("C - Predicted, all species", "D - Predicted, endemics"),4)
figura <- "Figure_SS_"
              
#Figure VV
vars = c("RLI_erro_pad", "RLI_erro_pad.end")
nomes = rep(c("A - All species","B - Endemics"),4)
figura <- "Figure_VV_"
              
              
for(i in 1:length(vars)) {
  var = vars[i]
  nome.fig = paste("figures/",figura, LETTERS[i],"_adpt_",var,".jpg",sep="")
  #nome.fig = paste("Figure2_",var,"_BW.jpg",sep="")
  nome.shp = gsub("\\.","_",paste("shp_centre_",var,sep=""))
  
  if (var == "N.total") {
    dados1 = grid.result[!is.na(grid.result[,var]),] 
    dados2 = log10(dados1[,var]+1)
    
    grid2.1 = grid1
    x.data.1 = as.double(coordinates(grid2.1)[,1])
    y.data.1 = as.double(coordinates(grid2.1)[,2])
    grid2.1 = cbind.data.frame(x.data.1, y.data.1, grid.result)
    names(grid2.1)[1:2] = c("x","y")
    coordinates(grid2.1) = ~x+y
    proj4string(grid2.1) <- raster::crs(grid1)
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
  
  ## Kriging
  # Fitting the semivariogram
  semivariog <- gstat::variogram(dados2 ~ 1, locations= grid3, data= grid3)
  
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
  if (grepl("RLI", var) & grepl("erro", var))
    krig.output = cbind.data.frame(coordinates(grid.out), krig$var1.var)
  
  names(krig.output)<-c("x","y","z")
  krig.output0 <- raster::rasterFromXYZ(krig.output)
  krig.output2 <- raster::crop(krig.output0, af.proj)
  krig.output3 <- raster::mask(krig.output2, af.proj)
  
  jpeg(filename = nome.fig, width = 3000*1.5, height = 3000*2, units = "px", pointsize = 12,
       res = 600, family = "sans", type="cairo", bg="white")
  # pdf(file = "figures/Figure2.pdf", width = 5, height = 7, pointsize = 12,
  #     family = "sans", bg="white")

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
  
  #Plot settings
  # raster::image(krig.output, col = cres, add=TRUE,  breaks = brks)
  if (grepl("Threat", var))
    cres <- rev(cres)
  # raster::image(krig.output0, col = cres, add=TRUE,  breaks = brks)
  
  # fixed breaks for scale
  brks1 <- seq(0.65, 1, 0.02)
  if (grepl("Threat", var) & grepl("_CR", var))
    brks1 <- seq(0, 0.45, 0.02)
  
  cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
  colfunc <- colorRampPalette(cores)
  
  cores1 <- c("darkred", "red", "gold", "green")
  colfunc1 <- colorRampPalette(cores1)
  cres1 <- colfunc(length(brks1)-1)
  if (grepl("Threat", var))
    cres1 <- rev(cres1)
  
  #ploting
  par(mfrow=c(1,1), las=1, mar=c(0,0,0,0), mgp=c(1.7,0.4,0), tcl=-.3)
  plot(af.proj, border="white")
  # raster::image(krig.output3, col = cres1, add=TRUE) #,  breaks = brks1)
  
  if (grepl("RLI", var) & !grepl("erro", var)) {
    # cores.plot <- cres1
    quebras <- seq(0.3, 1, 0.025)
    # quebras <- seq(0.30, 0.75, 0.025)
    cores <- c("darkred", "red", "darkorange", "gold", "yellow")
    colfunc <- colorRampPalette(cores)
    cores.plot <- colfunc(length(quebras)-1)
    raster::image(krig.output3, col = cores.plot, add=TRUE) #,  breaks = brks1)
    zlim <- c(30,75)
    breaks <- seq(0,100,20)
    at.leg <- seq(30,70,10)
    at.leg.lab <- at.leg/100
  }  
  
  if (grepl("Threat", var)) {
    # cores.plot <- cres1
    quebras <- seq(0,1,0.05)
    cores <- c("darkred", "red", "darkorange", "gold", "yellow")
    colfunc <- colorRampPalette(cores)
    cores.plot <- rev(colfunc(length(quebras)-1))
    
    raster::image(krig.output3, col = cores.plot, add=TRUE) #,  breaks = brks1)
    zlim <- c(35,100)
    breaks <- seq(0,100,20)
    at.leg <- seq(0,100,20)
    at.leg.lab <- at.leg/100
  }
  
  if (grepl("Threat", var) & grepl("CR", var)) {
    # cores.plot <- cres1
    quebras <- seq(0, 0.35, 0.025)
    cores <- c("darkred","red", "gold", "green")
    colfunc <- colorRampPalette(cores)
    cores.plot <- rev(colfunc(length(quebras)-1))
    
    raster::image(krig.output3, col = cores.plot, add=TRUE) #,  breaks = brks1)
    zlim <- c(0,55)
    breaks <- seq(0,55,10)
    at.leg <- seq(0,55,10)
    at.leg.lab <- at.leg/100
  }
  
  # at.leg <- seq(0,100,20)
  # at.leg1 <- seq(60,100,10)
  # at.leg1 <- seq(0,60,10)
  # vals.leg <- paste0(at.leg,"%")
  
  if (grepl("RLI", var) & grepl("erro", var)) {
    # cores.plot <- cres1
    quebras <- seq(7e-05, 0.00018, 0.075e-04)
    cores <- c("darkred","red3","red", "lightpink2","cornflowerblue", "cyan","blue1")
    colfunc <- colorRampPalette(cores)
    cores.plot <- rev(colfunc(length(quebras)-1))
    
    raster::image(krig.output3, col = cores.plot, add=TRUE) #,  breaks = brks1)
    zlim <- c(7e-05,0.00018)
    breaks <- seq(0,100,20)
    at.leg <- seq(7e-05,0.00018,0.15e-04)
    at.leg.lab <- at.leg/100
    
  }    
  
  plot(brasil.proj1, add=TRUE, lty="1414", border= cor.br)
  plot(af.proj, add=TRUE, border= cor.af)
  plot(am.lat.proj1, add=TRUE, border = cor.am)
  nome = nomes[i]
  legend("topleft", nome, cex=1.5, bty="n",inset=c(0.02,0.03))
  # plot(shp.krig.cut1, add=TRUE, lwd=2, col=1)
  
  par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
  fields::image.plot( legend.only=TRUE, add=TRUE,
                      col = cores.plot,
                      legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
                      legend.line = 1, legend.cex = 0.5,
                      zlim = zlim,
                      # zlim=c(2.5,97.5),
                      # zlim=c(60,100),
                      # zlim=c(0,55),
                      #breaks = breaks,
                      #breaks = seq(0,100,20),
                      axis.args=list(at= at.leg,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
                                     # labels = round(qts[names(qts) %in% vals.leg],2)))
                                     labels= at.leg.lab))
  if (grepl("Threat", var))  
    legend(4500000,6024198,expression(bold("Prop. Threatened Species")),cex=1.1,bty="n",y.intersp=-0.01)
  if (grepl("RLI", var) & !grepl("erro", var))  
    legend(4900000,6024198,expression(bold("Red List Index")),cex=1.1,bty="n",y.intersp=-0.01)
  if (grepl("RLI", var) & grepl("erro", var))  
    legend(4800000,6024198, expression(bold("Prediction error")),cex=1.1,bty="n",y.intersp=-0.01)
  #colorbar.plot(5549396, 7757332, strip=cbind(0:100), horizontal=F, 
  #	    strip.width = 0.04, strip.length = 0.6,
  #      col=cres)
  #axis(2, at=seq(0,100,20), labels = TRUE, las=1)
  #plot(shp.idw1, add=TRUE, lwd=2, col=6)
  #plot(shp.krig1, add=TRUE, lwd=2, col=1)
  
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
grid <- readRDS("data/af_hex_grid_50km.rds")

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


#### FIGURE 3 ####
cor.af = "black" # "darkgreen" # 
cor.br = "grey"
cor.am = "black" 
cutoff = "25%"
#vars = c("N.total","S.total")
vars = c("RLI","RLI.end")
# vars = c("RLI.end", "RLI_catA.end", "RLI_catB.end", "RLI_catCD.end")


#nomes = c("A - Number of occurrences","B - Number of species")
nomes = rep(c("A - All species","B - Endemics only"),4)
# nomes = rep(c("A - All criteria","B - Criterion A","C - Criterion B", "D - Criteria C and D"),2)
am.lat.proj1 = raster::crop(am.lat.proj, bb)
brasil.proj1 = raster::crop(brasil.proj, bb)
  
for(i in 1:length(vars)) {
  var = vars[i]
 #nome.fig = paste("figures/Figure3_hex_",var,".jpg",sep="")
 nome.fig = paste("Figure2_",var,"_BW.jpg",sep="")
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
 
 ## Kriging
 # Fitting the semivariogram
 semivariog <- gstat::variogram(dados2 ~ 1, locations= grid3, data= grid3)
 
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
 
 dev.off()
}
                    
# -------------------------------------------------------------------------
##################################################
#### SPECIES WITH RECORDS OLDER THAN 50 YEARS ####
##################################################
#### FIGURE SP  ####
old <- readRDS("data/old_records_xy.rds")

# creating the spatial object and projecting
old.xy <- as.data.frame(old[, .SD, .SDcols = c("ddlon", "ddlat")])
coordinates(old.xy) <- ~ddlon + ddlat
proj4string(old.xy) <- raster::crs(af)
old.xy <- spTransform(old.xy, raster::crs(af.proj))

#Defining the panel extents
all <- sul <- centro <- norte <- raster::extent(grid1)
all[2] <- 6250000
sul[1:2] <- c(3800000, 4450000)
sul[3:4] <- c(6600000, 7100000)

centro[1:2] <- c(4511378, 5570000)
centro[3:4] <- c(7200000, 8600000)

norte[1:2] <- c(5168602, 6250000)
norte[3:4] <- c(8966880, 9750000)

plots.ext <- list(all, norte, centro, sul)
names(plots.ext) <- c("A- Atlantic Forest", "B - Northern", "C - Central",
 "D - Southern")

# colors and symbols
cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
cores.old <- rep("darkorange", dim(old)[1])
cores.old[old$category %in% "CR"] <- "red"
cores.old[old$category %in% "CR_PE"] <- "purple"
cores.old <- adjustcolor(cores.old, alpha.f = 0.7)
    
# pchs.old <- rep(21, dim(old)[1])
# pchs.old[old$only.type] <- 24

# Cities (large ones and those close to the occurrences)
cidades.all <- rnaturalearth::ne_download(scale = 10, 
  type = 'populated_places_simple', category = 'cultural')
cidades.all <- cidades.all[cidades.all@data$iso_a2 %in% c("BR","AR","PY"), ]
cidades.extra <- c("Tauá", "Baturité", "João Pessoa", "Aracaju", 
   "Juazeiro do Norte", "Caaruaru", "Campina Grande",
   "Ilhéus", "Porto Seguro", "Ubaitaba", "Almenara",
   "São Mateus", "Linhares", "Governador Valadares", "Cachoeiro de Itapemirim",
   "Caratinga", "Ponte Nova", "Sete Lagoas", "Juiz de Fora", "Ipatinga",
   "Barbacena", "Petrópolis", "Cabo Frio", "Nova Friburgo", "Barra Mansa",
   "Santos", "Santo André", "Taubaté","Itanhaem",
   "Palmas", "Brusque", "Itajaí", "Caçador", "Lajes")
cidades.all <- cidades.all[cidades.all@data$name %in% cidades.extra,]
cidades <- rnaturalearth::ne_download(scale = 50, 
  type = 'populated_places_simple', category = 'cultural')
cidades <- cidades[cidades@data$iso_a2 %in% c("BR","AR","PY"), ]
cidades <- rbind(cidades, cidades.all)
cidades.rm <- c("Sorocaba", "Goiânia", "Alvorada", "Feira de Santana",
"Ciudad del Este")
cidades <- cidades[!cidades@data$name %in% cidades.rm,]
# cidades <- cidades[cidades@data$adm0name %in% "Brazil", ]
cidades@data <- cidades@data[, c("featurecla","adm0name","adm1name","name",
 "pop_max","latitude","longitude")]
# Positioning the city names in the legend
cidades@data$pos <- 4
pos_1 <- c("Governador Valadares", "Belo Horizonte", "Vitória da Conquista",
   "Rio de Janeiro", "Campina Grande")
cidades@data$pos[cidades@data$name %in% pos_1] <- 1
# pos_2 <- c("Campina Grande")
# cidades@data$pos[cidades@data$name %in% pos_2] <- 2
pos_3 <- c("Sete Lagoas", "Petrópolis", "São Paulo")
cidades@data$pos[cidades@data$name %in% pos_3] <- 3
# Projecting the cities data
cidades.proj <- spTransform(cidades, raster::crs(af.proj))

# Defining the image size
imag.size <- list(c(1.5, 2), c(2, 1.5), c(1.5, 2), c(2, 1.5))

for(i in 1:length(plots.ext)) {
  
  nome.fig <- paste0("figures/Figure_SP",
                 substring(names(plots.ext)[i], 1,1),".jpg")
                      xlim.i <- plots.ext[[i]][1:2]
                          ylim.i <- plots.ext[[i]][3:4]
                          fatores <- imag.size[[i]]
                          
jpeg(filename = nome.fig, width = 3000*fatores[1], height = 3000*fatores[2], 
 units = "px", pointsize = 12,
 res = 600, family = "sans", type="cairo", bg="white")
                          
#Defining the breaks for colors 
if (i == 1) {
  par(mfrow=c(1,1), las=1, mar=c(0,0,0,0), mgp=c(1.7,0.4,0), tcl=-.3)
  plot(af.proj, border="white", xlim = xlim.i, ylim = ylim.i)
  plot(brasil.proj1, add=TRUE, 
       #lty="1414", 
       border= cor.br)
  plot(af.proj, add=TRUE, border= cor.af)
  plot(am.lat.proj1, add=TRUE, border = cor.am, lwd = 2)
  plot(old.xy, add = TRUE,
       pch = 21, bg = cores.old, cex = 1)
  
  plot(plots.ext[[2]], add = TRUE)  
  plot(plots.ext[[3]], add = TRUE)  
  plot(plots.ext[[4]], add = TRUE)  
  dif.x <- 50000
  dif.y <- 50000
  text(plots.ext[[2]][2] - dif.x, 
       plots.ext[[2]][3] + dif.y, expression(bold(B)), cex = 1.2)
  text(plots.ext[[3]][2] - dif.x, 
       plots.ext[[3]][3] + dif.y, expression(bold(C)), cex = 1.2)
  text(plots.ext[[4]][2] - dif.x, 
       plots.ext[[4]][3] + dif.y, expression(bold(D)), cex = 1.2)
  # text(6200000, 9030000, expression(bold(B)), cex = 1.2)
  # text(5471266, 7222421, expression(bold(C)), cex = 1.2)
  # text(3777789, 6806204, expression(bold(D)), cex = 1.2)
  
  legend("topleft", names(plots.ext)[i], cex=1.5, bty="n",inset=c(0.0,0.05))
}

if (i != 1) {
  par(mfrow=c(1,1), las=0, mar=c(3,3,1,2), mgp=c(1.7,0.4,0), tcl=-.3)
  af.proj.crop <- raster::crop(af.proj, plots.ext[[i]])
  plot(af.proj.crop, border="white", xlim = xlim.i, ylim = ylim.i,  
       axes = TRUE)
  # plot(af.proj, border="white", xlim = xlim.i, ylim = ylim.i, axes = TRUE)
  plot(brasil.proj1, add=TRUE, 
       #lty="1414", 
       border= cor.br)
  plot(af.proj, add=TRUE, border= cor.af, lwd = 2)
  # plot(am.lat.proj1, add=TRUE, border = cor.am, lwd = 2)
  plot(cidades.proj, add=TRUE, pch = 15, cex = 1.3, col = "black")
  plot(old.xy, add = TRUE,
       pch = 21, bg = cores.old, cex = 1.5)
  text(sp::coordinates(cidades.proj),cidades.proj@data$name, 
       pos = cidades.proj@data$pos)
  legend("bottomright", c("CR_PE","CR","EN"),
         pch=21, 
         pt.bg = adjustcolor(c("purple", "red", "darkorange"), alpha.f = 0.8), 
         bty= "n", cex=1.3)
  legend("topleft", names(plots.ext)[i], cex=1.5, bty="n",inset=c(-0.02,0.0))
  }
dev.off()
}


# -------------------------------------------------------------------------
##################
#### FIGURE 4 ####
##################
#The Red List Index for other tropical forests

#Getting the world map
# world0 <- rnaturalearth::ne_download(scale = 110, type = 'map_units', category = 'cultural')
world0 <- rnaturalearth::ne_download(scale = 110, type = 'land', category = 'physical')
world0 <- rgeos::gBuffer(world0, byid=TRUE, width=0) #correcting possible overlapping polygons
world0 <- cleangeo::clgeo_Clean(world0) # fixing possible problems with the new aggregated polygon (e.g. orphaned holes)
ext.wo <- extent(world0)
ext.wo[1] <- -130
ext.wo[3] <- -57
ext.wo[4] <- 33
wo1 <- raster::crop(world0, ext.wo)
plot(wo1)
# Projecting the world map
crs.robin <- sp::CRS("+proj=robin +datum=WGS84 +no_defs")
world.proj <- sp::spTransform(wo1, crs.robin) #Robinson

#Getting the hostposts and other tropical regions
hot <- readRDS("data/hotspots_limits.rds")
hot.proj <- sp::spTransform(hot, crs.robin) #Robinson
hot.proj.outer <- hot.proj[hot.proj@data$Type %in% "outer limit",]
hot.proj <- hot.proj[hot.proj@data$Type %in% "hotspot area",]

#How much of the earth land surface the hostspots represent?
100*rgeos::gArea(hot.proj)/rgeos::gArea(world.proj)


#Defining the legend colors
brks1 <- seq(0.3, 1, 0.025)
# cores <- c("darkred", "red", "darkorange", "gold", "yellowgreen")
# cores1 <- c("darkred", "red", "gold", "green")
cores1 <- c("darkred", "red", "darkorange", "gold", "yellow", "yellowgreen", "chartreuse3")
colfunc1 <- colorRampPalette(cores1)
cres1 <- colfunc1(length(brks1)-1)

#Defining the region colors
info <- readRDS("data/hotspots_results.rds")
info <- info[match(hot.proj$NAME, info$hotspot.region), ]
info$RLI.col = as.character(cut(info$median_rli, 
    # breaks = quantile(grid.result$RLI, na.rm=TRUE, prob = seq(0,1,0.05)), 
        breaks = brks1, 
        labels = colfunc1(length(brks1)-1), include.lowest = TRUE))

#Legend regions
coords <- sp::coordinates(rgeos::gCentroid(hot.proj, byid = TRUE))
coords <- cbind(coords, nome = hot.proj$NAME)
table.order <- c("Atlantic Forest", "Caribbean Islands", 
                 "Coastal Forests of Eastern Africa", "Eastern Afromontane", 
                 "Guinean Forests of West Africa", "Indo-Burma", "Madagascar and the Indian Ocean Islands", 
                 "Mesoamerica",  "New Caledonia", "Philippines", "Sundaland", "Tropical Andes", 
                 "Tumbes-Choco-Magdalena", "Wallacea", "Western Ghats and Sri Lanka", "Amazon", 
                 "Central Africa", "New Guinea")
coords <- coords[match(table.order, coords[,3]), ]
coords <- cbind(coords, code = 1:dim(coords)[1])
# Adapting the legend coordinates for better visualization
x.coords <- y.coords <- rep(0, dim(coords)[1])
names(x.coords) <- names(y.coords) <- coords[,3]
x.coords[c(1, 2, 3, 4, 7, 8, 10, 11, 12, 13, 14, 15)] <- 
  c(900000,1600000,400000,500000,500000,-200000,850000,-900000,-1000000,-500000,
    -500000,-600000)
coords[,1] <- as.double(coords[,1]) + x.coords
y.coords[c(3, 4, 5, 8, 9, 11, 13, 14, 18)] <- 
  c(700000,-100000,-500000,-500000,-500000,-800000,-200000,-1000000,600000)
coords[,2] <- as.double(coords[,2]) + y.coords

#Plotting the figure itself
jpeg(filename = "figures/Figure4.jpg", width = 5400*2, height = 5400*1.1, units = "px", pointsize = 12,
     res = 600, family = "sans", type="cairo", bg="white")
pdf(file = "figures/Figure4.pdf", width = 14, height = 6, pointsize = 12,
    family = "sans", bg="white")
par(mfrow=c(1,1), las=1,mar=c(0.5,0,0.5,0.5),mgp=c(1.7,0.4,0),tcl=-.3)
## Tropical forest area
sp::plot(world.proj, border = "black", col = "lightgrey", # main = names(cores)[i],
         xlim=c(-10000000,16000000), ylim = c(-5700000, 1300000))
sp::plot(hot.proj, border = "black", col = info$RLI.col, add = TRUE)
sp::plot(hot.proj.outer, border = "black", add = TRUE)

# legend("topleft", expression(bold("A)")), cex=1.5, bty="n", x.intersp =  -1, y.intersp = 1.5)
at.leg <- seq(0.3,1,0.1)
fields::image.plot( legend.only=TRUE, col = cres1, add=TRUE,legend.shrink = 1, legend.width = 2, legend.mar = 1.2, 
                    legend.line = 0, legend.cex = 0.5, 
                    zlim= c(0.3,1),
                    smallplot = c(0.1,0.875,0.135,0.16),horizontal = TRUE,	 
                    axis.args=list(at = at.leg,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=0.9, 
           labels = at.leg))
legend(10300000,-5700000,expression(bold("Red List Index")),
       cex=1.2, bty="n", y.intersp=-0.01)
text(coords[,1:2], coords[,4], cex = 1.1, col=1, font = 2)
dev.off()


# -------------------------------------------------------------------------
########################################################
#### DISTRIBUTIONS AND INTERPOLATION FOR CR SPECIES ####
########################################################
#### FIGURE SS  ####
## PANELS C AND D TOGETHER WITH THE KRIGING FIGURES ABOVE

colunas <- c("Threat.pa", "Threat.end.pa", "Threat.pa_CR", "Threat.end.pa_CR")
apply(grid.result[, colunas], 2 , summary)

#### OBSERVED OCCURRENCESS ####
## Defining the color scales
brks1 <- seq(0, 0.35, 0.025)
cores1 <- c("darkred","red", "gold", "green")
colfunc1 <- colorRampPalette(cores1)
cres1 <- rev(colfunc1(length(brks1)-1))
grid.result$Threat.col_CR = as.character(cut(grid.result$Threat.pa_CR, 
                     breaks = brks1,  
                     labels = rev(colfunc1(length(brks1)-1)), include.lowest = TRUE))
grid.result$Threat.end.col_CR = as.character(cut(grid.result$Threat.end.pa_CR, 
 breaks = brks1,  
 labels = rev(colfunc1(length(brks1)-1)), include.lowest = TRUE))

# Grids with no data -> "White" ()
grid.result$Threat.col_CR[is.na(grid.result$select)] <- "#FFFFFF"
grid.result$Threat.end.col_CR[is.na(grid.result$select)] <- "#FFFFFF"
grid.result$Threat.end.col_CR[is.na(grid.result$Threat.end.col_CR)] <- "#FFFFFF"
      
# Grids with insuficient overall data -> "Light Grey" ("#D3D3D3") or "Grey" ("#808080") or "Dark Grey" ("#A9A9A9")
grid.result$Threat.col_CR[!is.na(grid.result$select) & 
                            !grid.result$select] <- "#D3D3D3"
grid.result$Threat.end.col_CR[!is.na(grid.result$select) &
                                !grid.result$select] <- "#D3D3D3" 

# Grids with suficient data but insuficient CR daya -> "Light Grey" ("#D3D3D3")
grid.result$Threat.col_CR[!is.na(grid.result$select) &
                            grid.result$select &
                            grid.result$SampCover_CR < 0.25] <- "#D3D3D3"
  
grid.result$Threat.end.col_CR[!is.na(grid.result$select) &
                                !grid.result$select &
                                grid.result$SampCover.end_CR < 0.15] <- "#D3D3D3"
  
# Setting the colors of the limits
cor.af = "black" # "darkgreen" # 
cor.br = grey(0.25)
cor.am = "black" 

#Plot extent plus cropping
bb = bbox(grid1)
bb[1] = bb[1] - 20000
bb[3] = bb[3] - 20000
bb[4] = bb[4] + 20000
am.lat.proj1 = raster::crop(am.lat.proj, bb)
brasil.proj1 = raster::crop(brasil.proj, bb)

## First panels ##
zlim <- c(35,100)
breaks <- seq(0,100,20)
at.leg <- seq(0,100,20)
at.leg.lab <- at.leg/100

# at.leg <- seq(0,100,20)
at.leg1 <- seq(0,55,10)
vals.leg <- paste0(at.leg,"%")

#PANEL A
jpeg(filename = "figures/Figure_SSA.jpg", width = 3000*1.5, height = 3000*2, 
 units = "px", pointsize = 12,
 res = 600, family = "sans", type="cairo", bg="white")

par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
plot(af.proj, border="white")
plot(grid1, border=grey(0.75), add = TRUE)
plot(brasil.proj1, add=TRUE, 
 #lty="1414", 
 border= cor.br)
plot(am.lat.proj1, add=TRUE, border = cor.am)

rm.cells <- !grid.result$Threat.col_CR %in% c("#FFFFFF", "#D3D3D3")
plot(grid1[rm.cells, ],
 col = adjustcolor(grid.result$Threat.col_CR[rm.cells],
   alpha.f = 0.6,
   offset = c(0.1, 0.1, 0.1, 0.15)), # <- "more white"
 # transform = diag(c(1,1,1,0.5))),
 border = "black",
 add=TRUE)
# plot(af.proj, add=TRUE, border= cor.af)
legend("topleft", "A - Observed, all species", cex=1.5, bty="n",inset=c(0.02,0.03))

par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
fields::image.plot( legend.only=TRUE, add=TRUE,
col = adjustcolor(cres1, alpha.f = 0.75),
legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
legend.line = 1, legend.cex = 0.5,
zlim=c(0,55),
axis.args=list(at= at.leg1,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
 labels= at.leg1/100))
legend(4500000,6024198,expression(bold("Prop. Threatened Species")),
   cex=1.1,bty="n",y.intersp=-0.01)
dev.off()

#PANEL B
jpeg(filename = "figures/Figure_SSB.jpg", width = 3000*1.5, height = 3000*2, 
 units = "px", pointsize = 12,
 res = 600, family = "sans", type="cairo", bg="white")

par(mfrow=c(1,1), las=1,mar=c(0,0,0,0),mgp=c(1.7,0.4,0),tcl=-.3)
plot(af.proj, border="white")
plot(grid1, border=grey(0.75), add = TRUE)
plot(brasil.proj1, add=TRUE, 
 #lty="1414", 
 border= cor.br)
plot(am.lat.proj1, add=TRUE, border = cor.am)

rm.cells1 <- !grid.result$Threat.end.col_CR %in% c("#FFFFFF", "#D3D3D3")
plot(grid1[rm.cells1, ],
 col = adjustcolor(grid.result$Threat.end.col_CR[rm.cells1],
   alpha.f = 0.6,
   offset = c(0.1, 0.1, 0.1, 0.15)), # <- "more white"
 # transform = diag(c(1,1,1,0.5))),
 border = "black",
 add=TRUE)
legend("topleft", "B - Observed, endemics", cex=1.5, bty="n",inset=c(0.02,0.03))

par(fig = c(0.5, 0.85, 0.015, 0.4), mar=c(0,0,0,0), new=TRUE)
fields::image.plot( legend.only=TRUE, add=TRUE,
    col = adjustcolor(cres1, alpha.f = 0.75),
    legend.shrink = 1, legend.width = 4, legend.mar = 1.2,
    legend.line = 1, legend.cex = 0.5,
    zlim=c(0,55),
    axis.args=list(at= at.leg1,  mgp=c(1.5,0.6,0), tcl=-.3, cex.axis=1, 
         labels= at.leg1/100))
legend(4500000,6024198,expression(bold("Prop. Threatened Species")),
           cex=1.1,bty="n",y.intersp=-0.01)
dev.off()
    