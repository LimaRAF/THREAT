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
# require(soiltexture)
# require(taxize)
# require(flora)
# require(gdata)
# require(rgeos)
# require(maptools)
# require(geosphere)

#### EXTRACTING ALL ENVIRONMENTAL CO-VARIABLES FROM TREECO ####

## Defining the general path to read/save files (tifs and shape files)
path = "E://ownCloud//W_GIS" 

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


# ##########################
# ### CLIMATIC VARIABLES ###
# ##########################
# 
# ##RELATIVE HUMIDITY - Brazil 100x100m (Alvares et al. 2013)
# ur <- raster(paste(path,"//BR_relative_humidity//UR_annual_BR1.tif",sep=""))
# pj = crs(ur) #defining and standardizing the geographical projections
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(ur,surveys2) # Relative humidity
# tmp1= extract(ur,surveys2,method="bilinear") # Relative humidity
# dados = cbind.data.frame(dados,umidade.rel=tmp,umidade.rel1=tmp1)
# rm(ur)
# 
# ##TEMPERATURE - Brazil 100x100m (Alvares et al. 2013)
# meses = c("_tmed_ann_fim1.tif",paste("_tmed_",tolower(month.abb),"_fim1.tif",sep=""))
# for(i in 0:12) {
#   tmp=c() # Precipitation
#   filename = paste(path,"//BR_mean_air_temp//",i,meses[i+1],sep="")
#   T = raster(filename)
#   pj = crs(T) 
#   if(is.na(pj)) {
#     pj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#     proj4string(T) <- pj } else {}
#   surveys2 = spTransform(surveys1, pj)
#   
#   tmp= extract(T,surveys2) # Annual temp
#   tmp1= extract(T,surveys2,method="bilinear")
#   tmp[is.na(tmp)] = tmp1[is.na(tmp)]
#   dados=cbind.data.frame(dados,tmp)
# }
# names(dados)[(dim(dados)[2]-12):dim(dados)[2]] = c("temp",paste("tmed",month.abb,sep="."))
# 
# ##PRECIPITATION - Brazil 100x100m (Alvares et al. 2013)
# meses = c("_precip_Ann_BR1.tif",paste("_precip_",tolower(month.abb),"1.tif",sep=""))
# for(i in 0:12) {
#   tmp=c() # Precipitation
#   filename = paste(path,"//BR_rainfall//",i,meses[i+1],sep="")
#   prec = raster(filename)
#   pj = crs(prec) 
#   if(is.na(pj)) {
#     pj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#     proj4string(prec) <- pj } else {}
#   surveys2 = spTransform(surveys1, pj)
#   
#   tmp= extract(prec,surveys2) # Annual ppt
#   tmp1= extract(prec,surveys2,method="bilinear")
#   tmp[is.na(tmp)] = tmp1[is.na(tmp)]
#   dados=cbind.data.frame(dados,tmp)
# }
# names(dados)[(dim(dados)[2]-12):dim(dados)[2]] = c("ppt",paste("ppt",month.abb,sep="."))
# dados$ppt.Jan = round(dados$ppt.Jan,0)
# 
# ## NUMBER OF DRY MONTHS, WALSH'S SCORES AND SAZONALITY INDEX ###
# #DRY MONTHS (<100m and <60mm)
# dados$dry.months100 = apply(dados[,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))],1,function(x) sum(x<100,na.rm=T))
# dados$dry.months60 = apply(dados[,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))],1,function(x) sum(x<60,na.rm=T))
# 
# #PPT DRIEST QUARTER (3 MONTHS) - equivalent to "ppt.3Driest" from BIOCLIM
# tmp = apply(dados[,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))],1,function(x) which(x==min(x)))
# for(i in 1:length(tmp)) {
#   tmp1 = dados[i,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))]
#   tmp2 = as.double(tmp[[i]])
#   if(length(tmp2)==1) {
#     opt1 = min(c(	try(sum(tmp1[[tmp2]],tmp1[[tmp2+1]],tmp1[[tmp2+2]]),TRUE),
#                   try(sum(tmp1[[tmp2-1]],tmp1[[tmp2]],tmp1[[tmp2+1]]),TRUE),
#                   try(sum(tmp1[[tmp2-2]],tmp1[[tmp2-1]],tmp1[[tmp2]]),TRUE)),na.rm=T) } else {}
#   if(length(tmp2)==2) {
#     opt1 = min(c(	try(sum(tmp1[[tmp2[1]]],tmp1[[tmp2[2]]],tmp1[[min(tmp2)-1]]),TRUE),
#                   try(sum(tmp1[[tmp2[1]]],tmp1[[tmp2[2]]],tmp1[[max(tmp2)+1]]),TRUE),
#                   try(sum(tmp1[[tmp2[1]]],tmp1[[tmp2[2]]],tmp1[[min(tmp2)+1]]),TRUE),
#                   try(sum(tmp1[[tmp2[1]]],tmp1[[tmp2[2]]],tmp1[[max(tmp2)-1]]),TRUE)),na.rm=T) } else {}
#   if(length(tmp2)==3) {
#     tmp3 = min(tmp2):max(tmp2)
#     opt1 = min(c(	sum(tmp1[[tmp3[1]]],tmp1[[tmp3[2]]],tmp1[[tmp3[3]]]),
#                   sum(tmp1[[tmp3[length(tmp3)-2]]],tmp1[[tmp3[length(tmp3)-1]]],tmp1[[tmp3[length(tmp3)]]])),na.rm=T)  } else {}
#   opt2 = sum(tmp1[,c("ppt.Jun","ppt.Jul","ppt.Aug")],na.rm=T)
#   opt3 = sum(tmp1[,c("ppt.Nov","ppt.Dec","ppt.Jan")],na.rm=T)
#   opt4 = sum(tmp1[,c("ppt.Dec","ppt.Jan","ppt.Feb")],na.rm=T)
#   opt5 = min(c(opt2,opt3,opt4),na.rm=T)
#   #cat(i,"\n")
#   if(as.double(opt1)<=as.double(opt5)) tmp[[i]] = as.double(opt1) else tmp[[i]] = as.double(opt5)  
# }
# dados$ppt.3Driest = round(as.double(unlist(tmp)),0)
# 
# ## WALSH'S SCORES or WALSH'S PERHUMIDITY INDEX
# tmp = dados[,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))]
# tmp[!is.na(tmp)&tmp<=50] = -2
# tmp[!is.na(tmp)&tmp>200] = 2
# tmp[!is.na(tmp)&tmp>100] = 1
# tmp[!is.na(tmp)&tmp>50] = -1
# 
# tmp1= cbind.data.frame(tmp,tmp[,1])
# res = NULL
# for(j in 1:dim(tmp1)[1]) {
#   rslt=rep(NA,12)
#   for (i in 1:11) {
#     foc = tmp1[j,i]
#     prox = tmp1[j,i+1]
#     rslt[1+i] = foc<0&prox>0
#     rslt[1] = tmp1[j,12]<0&tmp1[j,13]>0
#     ad = sum(rslt,na.rm=T)*0.5
#     res[j] = ad
#   }
# }
# dados$walsh.score = apply(tmp,1,sum,na.rm=T) + res
# 
# ## WALSH'S SEASONALITY INDEX
# tmp = dados[,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))]
# tmp1 = apply(tmp,1,sum,na.rm=T)
# dados$walsh.seasonality.index = NA
# 
# for(i in 1:dim(tmp)[1]) {
#   dados$walsh.seasonality.index[i] =sum(abs(tmp[i,] - tmp1[i]/12))/tmp1[i]
# }
# 
# ## PPT SEASONALITY
# tmp = dados[,c(grep("ppt.Jan",names(dados)):grep("ppt.Dec",names(dados)))]
# tmp1 = apply(tmp,1,sd)
# tmp2 = apply(tmp,1,mean)
# dados$ppt.Seasonality = tmp1/tmp2
# 
# ### CLIMATIC WATER DEFICIT/STRESS ###
# # CWD: grain size is 0.041666666666667 degree, or about 1 cell every 5 km
# # CWD= sum(PPT - PET)
# cwd <- raster(paste(path,"//WO_CWD//CWD.bil",sep=""))
# 
# tmp= extract(cwd,surveys2)
# tmp1= extract(cwd,surveys2,method="bilinear") 
# tmp[is.na(tmp)] = tmp1[is.na(tmp)] 
# dados = cbind.data.frame(dados,CWD=tmp,CWD1=tmp1)
# rm(cwd)
# 
# ##Koeppen
# koep <- readOGR(dsn=paste(path,"//BR_koppen",sep=""),layer="BR_koppen")
# pj = crs(koep) # Homogenizing projections
# surveys2 = spTransform(surveys1, pj)
# tmp = over(surveys2,koep, fn=NULL)
# dados = cbind.data.frame(dados,koeppen=as.character(tmp$gridcode))
# dados$koeppen = as.character(dados$koeppen)
# 
# dados[dados$koeppen %in% 1,]$koeppen = "Cwa" 
# dados[dados$koeppen %in% 2,]$koeppen = "Am"
# dados[dados$koeppen %in% 3,]$koeppen = "Af" #ok
# dados[dados$koeppen %in% 4,]$koeppen = "Cfa"
# dados[dados$koeppen %in% 5,]$koeppen = "Cwb"
# #dados[dados$koeppen %in% 6,]$koeppen = "Csb" #nested within Csa (altitudes >1000m)
# #dados[dados$koeppen %in% 7,]$koeppen = "Csa" #nested in As/BSh (altitudes 800-1000m)
# dados[dados$koeppen %in% 8,]$koeppen = "Cfb"
# #dados[dados$koeppen %in% 9,]$koeppen = "BSh"
# #dados[dados$koeppen %in% 10,]$koeppen = "As"
# #dados[dados$koeppen %in% 11,]$koeppen = "Cwc" 
# dados[dados$koeppen %in% 12,]$koeppen = "Aw"
# rm(koep)
# 
# ### ENVIRONMENTAL STREESS ###
# # E= 1.e-3 x (0.178xTseasonality - 0.938xCWD - 6.61xPPTseasonality)
# # Tseasonality= bioclim 4; PPTseasonality= bioclim 15
# e <- raster(paste(path,"//WO_Environmental_Stress_Chave_etal_2014//E.tif",sep=""))
# pj = crs(e) 
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(e,surveys2) # annual evapot
# tmp1= extract(e,surveys2,method="bilinear") # annual evapot
# tmp[is.na(tmp)] = round(tmp1[is.na(tmp)],7) 
# dados = cbind.data.frame(dados,EnvStress=tmp,EnvStress1=tmp1)
# 
# ### BIOCLIM ###
# #getting the list of files available
# files = list.files(paste(path,"//WO_BioClim_v2.0",sep=""),pattern="LA.tif",full.name=TRUE)
# #extracting the climate variable month-by-month
# for(i in 1:length(files)) {
#   tmp=c() # Bioclim
#   bioc = raster(files[i])
#   surveys2 = spTransform(surveys1, crs(bioc))
#   tmp= extract(bioc,surveys2)
#   tmp[is.na(tmp)] = extract(bioc,surveys2[is.na(tmp),],method="bilinear")
#   dados=cbind.data.frame(dados,tmp)
# }
# names(dados)[(dim(dados)[2]-18):dim(dados)[2]] = 
#   c("Mean.Temp","Diurnal.Range","Isothermality","Temp.Seasonality",
#     "Max.Temp.Warmest","Min.Temp.Coldest","Temp.Annual.Range",
#     "Mean.Temp.3Wettest","Mean.Temp.3Driest","Mean.Temp.3Warmest",
#     "Mean.Temp.3Coldest","Annual.ppt","ppt.Wettest","ppt.Driest",
#     "ppt.Seasonality.wclim","ppt.3Wettest","ppt.3Driest.wclim","ppt.3Warmest",
#     "ppt.3Coldest")
# 
# ## CLOUD COVER (%) AND FROST DAY FREQUENCY (days) - CRU version 4.02 (1901-2017), 0.5 degress resolution
# #getting the list of files available
# cloud = raster(paste(path,"//WO_CRU_version4//cru_ts4.02.1901.2017.cloud_mean.tif",sep=""))
# pj = crs(cloud)
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(cloud, surveys2) # cloud cover
# tmp1= extract(cloud, surveys2, method="bilinear") # cloud cover
# tmp[is.na(tmp)] = tmp1[is.na(tmp)]
# dados = cbind.data.frame(dados,CloudCover=tmp,CloudCover1=tmp1)
# 
# #getting the list of files available
# frost = raster(paste(path,"//WO_CRU_version4//cru_ts4.02.1901.2017.frost_mean.tif",sep=""))
# pj = crs(frost)
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(frost, surveys2) # cloud cover
# tmp1= extract(frost, surveys2, method="bilinear") # cloud cover
# tmp[is.na(tmp)] = tmp1[is.na(tmp)]
# dados = cbind.data.frame(dados,FrostDays=tmp,FrostDays1=tmp1)
# 
# #getting the list of files available
# pet = raster(paste(path,"//WO_CRU_version4//cru_ts4.02.1901.2017.pet_mean.tif",sep=""))
# pj = crs(pet)
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(pet, surveys2) # cloud cover
# tmp1= extract(pet, surveys2, method="bilinear") # cloud cover
# tmp[is.na(tmp)] = tmp1[is.na(tmp)]
# dados = cbind.data.frame(dados,PET=tmp,PET1=tmp1)
# 
# ### ACTUAL EVAPOTRANSPIRATION ### based on WorldClim data (1950-2000) - resolution of 30 secs (mm 
# #calculated based on PET, vegetation type and soil characteristics
# aet <- raster(paste(path,"//WO_mean_annual_AET//aet_yr//aet_yr_AMSUL.tif",sep=""))
# pj = crs(aet) 
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(aet,surveys2) # annual evapot
# tmp1= extract(aet,surveys2,method="bilinear") # annual evapot
# tmp[is.na(tmp)] = round(tmp1[is.na(tmp)],0) 
# dados = cbind.data.frame(dados,AET=tmp,AET1=tmp1)
# 
# ### SOIL WATER STRESS ### based on WorldClim data (1950-2000) - resolution of 30 secs (mm) 
# #calculated based on Effective PPT (Eppt), Actual Evapo (AET) and Runoff (R):  dSWC = Eppt - AET - R
# #extracting the climate variable month-by-month
# #res=NULL
# for(i in 1:12) {
#   tmp=c() # Minimum temperature
#   #getting the list of files available
#   file = list.files(paste(path,"//WO_soil_water_stress//swc_fr_",i,"//",sep=""),pattern="w001001.adf",full.name=TRUE)
#   swc = raster(file)
#   surveys2 = spTransform(surveys1, crs(swc))
#   tmp= extract(swc,surveys2)
#   tmp[is.na(tmp)] = extract(swc,surveys2[is.na(tmp),],method="bilinear")
#   #if(is.null(res)) res=tmp else res=cbind.data.frame(res,tmp)
#   dados=cbind.data.frame(dados,tmp)
# }
# names(dados)[(dim(dados)[2]-11):dim(dados)[2]] = paste("SWC",month.abb,sep=".")
# dados$SWC.mean = apply(dados[,(dim(dados)[2]-11):dim(dados)[2]],1,mean)
# 
# rm(aet,bioc,cloud,frost,prec,swc,T,pet,e)
# 
# #################
# ### SOIL DATA ###
# #################
# 
# ### SOIL TYPES - IBGE 2017 - formely was Embrapa 2011
# soil.br <- readRDS(paste(path,"//BR_solos_IBGE//Pedologia_area_Brasil_treeco.rds",sep=""))
# pj = crs(soil.br) # Homogenizing projections
# surveys2 = spTransform(surveys1, pj)
# tmp = over(surveys2,soil.br, fn=NULL)
# dados = cbind.data.frame(dados,symbol=as.character(tmp$LEGENDA),soil0=as.character(tmp$COMPONENTE),soil1=as.character(tmp$COMPONENT1),soil2=as.character(tmp$COMPONENT2),soil3=as.character(tmp$COMPONENT3),
#                          Order=as.character(tmp$Order), Suborder=as.character(tmp$Suborder),	
#                          Group=as.character(tmp$Group), Subgroup=as.character(tmp$Subgroup))
# rm(soil.br)
# dados$soil0 = as.character(dados$soil0); dados$soil0 = gsub('\\\n','',dados$soil0) 
# dados$soil1 = as.character(dados$soil1); dados$soil1 = gsub('\\\n','',dados$soil1) 
# dados$soil2 = as.character(dados$soil2); dados$soil2 = gsub('\\\n','',dados$soil2) 
# dados$soil3 = as.character(dados$soil3); dados$soil3 = gsub('\\\n','',dados$soil3) 
# 
# ### DATA FROM THE Harmonized World Soil Database ### - v 1.2
# #Reading and extracting the values
# hwsd <- raster(paste(path,"//WO_soil_HWSD_121//hwsd_LA.tif",sep=""))
# crs(hwsd) = crs.geo
# tmp = extract(hwsd,surveys1)
# 
# #Loading the data file and the metadata codes and descriptions
# dado = read.csv(paste(path,"//WO_soil_HWSD_121//HWSD_DATA.csv",sep=""),sep=",",na.string=c(""," ",NA))
# meta = read.csv(paste(path,"//WO_soil_HWSD_121//HWSD_Metadata.csv",sep=""),sep=",",na.string=c(""," ",NA))
# fao = read.csv(paste(path,"//WO_soil_HWSD_121//D_SYMBOL90.csv",sep=""),sep=",",na.string=c(""," ",NA))
# 
# #Merging and re-classifying information
# tmp = left_join(data.frame(ordem=surveys1@data$ordem,ID=tmp),dado,by="ID")
# tmp[,"ISSOIL"] = gsub("0","0_no",tmp[,"ISSOIL"]); tmp[,"ISSOIL"] = gsub("1","1_yes",tmp[,"ISSOIL"]) # is soil? 
# tmp = merge(tmp,fao[,c("CODE","VALUE")],by.x="SU_CODE90",by.y="CODE",all.x=TRUE) # FAO 1990 soil types
# tmp1 = strsplit(as.character(meta[meta$FIELD %in% "DRAINAGE",]$CODES[1]),";")[[1]]
# tmp1 = cbind.data.frame(code=sapply(strsplit(tmp1,"="),function(x) x[1]),DRAINAGE=gsub("=","_",tmp1))
# tmp = merge(tmp,tmp1[,c("code","DRAINAGE")],by.x="DRAINAGE",by.y="code",all.x=TRUE) # Drainage
# tmp1 = strsplit(as.character(meta[meta$FIELD %in% "T_USDA_TEX_CLASS",]$CODES[1]),";")[[1]]
# tmp1 = cbind.data.frame(code=sapply(strsplit(tmp1,"="),function(x) x[1]),T_usda_texture=gsub("=","_",tmp1))
# tmp = merge(tmp,tmp1[,c("code","T_usda_texture")],by.x="T_USDA_TEX_CLASS",by.y="code",all.x=TRUE) # Soil texture
# 
# #Filtering and renaming variables
# tmp = tmp[,c("ordem","ISSOIL","VALUE","DRAINAGE","T_GRAVEL","T_SAND",            
#              "T_SILT","T_CLAY","T_usda_texture","T_REF_BULK_DENSITY","T_BULK_DENSITY",    
#              "T_OC","T_PH_H2O","T_CEC_CLAY","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE",             
#              "S_GRAVEL","S_SAND","S_SILT","S_CLAY","S_USDA_TEX_CLASS","S_PH_H2O","S_CEC_SOIL","S_TEB")]
# names(tmp) = gsub("T_","top\\.",names(tmp))
# names(tmp) = gsub("S_","sub\\.",names(tmp))
# tmp1 = tmp[match(dados$dados,tmp$ordem,nomatch=0),]
# #table(tmp1$ordem == dados$dados)
# dados = cbind.data.frame(dados,tmp1[,-(1)])
# 
# ### Soil Quality - Harmonized World Soil Database - 30 segs (~800-1000m)
# #Classes are:
# #1: No or slight limitations
# #2: Moderate limitations
# #3: Sever limitations
# #4: Very severe limitations
# #5: Mainly non-soil
# #6: Permafrost area
# #7: Water bodies
# #Only classes 1 to 4 correspond to soil limitations for plant growth (corn). 
# #Class 1 is generally rated between 80 and 100% of the growth potential, class 2 between 60 and 80%, class 3 between 40 and 60%, and class 4 less than 40%.
# 
# #Reading and extracting the values
# cats = c("nutrients","nutr.retention","rooting","oxygen","salts","toxicity","mechanization")
# res = NULL
# for(i in 1:length(cats)) {
#   sq <- raster(paste(path,"//WO_soil_HWSD_121//soil_quality//sq",i,"_LA.tif",sep=""))
#   tmp = extract(sq,surveys1)
#   #tmp1 = extract(sq,surveys1[tmp %in% 0,],buffer=9000)
#   #tmp1 = unlist(sapply(tmp1,function(x) modal(x[x>0])))
#   #tmp[tmp==0] = tmp1
#   if(is.null(res)) res= tmp else res=cbind.data.frame(res,tmp)
# }
# names(res) = cats
# res1 = apply(res[,c(1:5)],1,sum,na.rm=TRUE) #calculating soil quality (no toxicity and mechanization)
# #res1[res1>=25] = NA #Excluding water bodies
# #res1[res1<5] = NA #Excluding sites with no data
# dados$soil.quality = res1
# #res[res>4] = NA
# names(res) = paste(names(res),"hwsd",sep=".")
# dados = cbind.data.frame(dados,res)
# 
# #removing objects
# rm(hwsd,dado,meta,fao,sq)
# 
# 
# ### SOIL APTITUDE AND MEAN PROPERTIES FROM TREECO ###
# soil.br1 <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Site conditions//soil_aptitude_treeco.csv",as.is=TRUE)
# soil.br1 <- soil.br1[,c("main_soil","drainage","depth.emb","swa","fert.apt","al.tox","hydra","Sand","pH_H2O","OM","CEC_pH7","SumBases","SatBases","SatAlum","CationRetention","FieldCapacity")]
# #soil.br1$main_soil
# 
# surveys1$main_soil = as.character(dados$Subgroup)
# surveys1$main_soil[is.na(surveys1$main_soil)|surveys1$main_soil=="NA"] = as.character(dados$Group[is.na(surveys1$main_soil)|surveys1$main_soil=="NA"])
# surveys1$main_soil[is.na(surveys1$main_soil)|surveys1$main_soil=="NA"] = as.character(dados$Suborder[is.na(surveys1$main_soil)|surveys1$main_soil=="NA"])
# tmp = left_join (as.data.frame(surveys1),soil.br1,by="main_soil")
# 
# dados = cbind.data.frame(dados,tmp[,c("drainage","depth.emb","swa","fert.apt","al.tox","Sand","pH_H2O","OM","CEC_pH7","SumBases","SatBases","SatAlum","CationRetention","FieldCapacity")])
# rm(soil.br1,tmp1)
# 
# # extracting the aptitude values of each soil properties
# dados$drainage1 = 	as.numeric(substr(dados$drainage,1,3)) 
# dados$depth.emb1=	as.numeric(substr(dados$depth.emb,1,3))
# dados$swa1=			as.numeric(substr(dados$swa,1,3))
# dados$fert.apt1=		as.numeric(substr(dados$fert.apt,1,3))
# dados$al.tox1=		as.numeric(substr(dados$al.tox,1,3))
# 
# #Calculating the semi-quantitative soil fitness variables
# dados$soil.aptitude = dados$drainage1 + dados$depth.emb1 + dados$fert.apt1 + dados$al.tox1
# 
# 
# #######################
# #### LANDFORM DATA ####
# #######################
# 
# ##ASPECT##
# #Orienta??o de vertente calculada utilizando a fun??o "aspect" GDAL-QGIS
# #0-360 (o Norte est? entre 1 e 359; o Sul ? 180); pixels sem declividade/orienta??o ficaram com valor 0
# #data from NASA Shuttle Radar Topography Mission (year 2000) - Georeferenced to WGS84 
# #resolution of 1 arc degree (30.8m at equator)
# 
# ##Note: Aspect as either degrees (more intuitive) or radians (often used for calculations). 
# #In particular, aspect is frequently converted to metrics of 'northness' and 'eastness', each ranging
# #-1 - 1. This is because an aspect of 0 and 360 is the same (both face north), so northness and eastness
# #break this down into two other components. Northness is calculated as the cosine of aspect (in radians),
# #and eastness is calculated as the sine of aspect (again, in radians). 1 represents East or North, and -1
# #represents West or South, respectively.
# 
# #some useful functions for calculations
# rad2deg <- function(rad) {(rad * 180) / (pi)}
# deg2rad <- function(deg) {(deg * pi) / (180)}
# 
# #getting the list of tiles
# path0 = "E://ownCloud//W_GIS//Am_Lat_MDT//aspect"
# myfiles <- list.files(path = path0, pattern = "tif",  full.names = TRUE)
# myfiles <- myfiles[grep("tif.aux.xml|tif.ovr",myfiles,invert=TRUE)]
# 
# #extracting the distance to the cosatline for each survey (paralellized)
# cl <- makeCluster(detectCores()-1)
# clusterEvalQ(cl, library(sp))
# clusterEvalQ(cl, library(raster))
# clusterEvalQ(cl, library(rgdal))
# clusterExport(cl, list("surveys1","myfiles","rad2deg","deg2rad"),envir=environment())
# tmp2 = parLapply(cl,1:length(myfiles), fun= function(i) {
#   # loading the corresponding landscape from the path folder
#   dcl <- try(raster(myfiles[i]),TRUE)
#   # getting the coordinate of the survey
#   ext = round(extent(dcl),0)
#   locs= surveys1[surveys1@coords[,"Latitude"]>=ext@ymin & surveys1@coords[,"Latitude"]<=ext@ymax,]
#   locs= locs[locs@coords[,"Longitude"]>=ext@xmin & locs@coords[,"Longitude"]<=ext@xmax,]
#   # checking if there are coordinates within tile ranges
#   if(length(locs)==0) { 
#     #res = cbind.data.frame(NA,NA,NA)
#     #names(res) = c("ordem","value","value1") 
#   } else { 
#     # homogeneizing projections
#     pj = crs(dcl)
#     locs1 = try(spTransform(locs, pj),TRUE)
#     # extracting the raster values (exact pixel) 
#     cc.ex = extract(dcl,locs1)
#     # extracting the raster values from nearby
#     cc = extract(dcl,locs1, fun=NULL, buffer=30*sqrt(2)+(30*sqrt(2)/2))
#     # manipulating the values to get the median exposition		
#     cc = sapply(cc, function(x) x[!x%in%0])
#     cc.rad = lapply(cc,deg2rad) #converting degrees to radians
#     cc.cos = sapply(cc.rad,cos) #converting radians to cosine (best for 'northness': 1=North; -1=South)
#     cc.sin = sapply(cc.rad,sin) #converting radians to cosine (best for 'eastness': 1=East; -1=West)
#     cc.med = sapply(cc.cos,median,na.rm=TRUE) #calculating the median cosine
#     cc.med1= sapply(cc.sin,median,na.rm=TRUE) #calculating the median sine
#     angle.med = rad2deg(atan2(cc.med1,cc.med)) #re-converting cosine and sine to degrees
#     angle.med[!is.na(angle.med)&angle.med<0] = 360 + angle.med[!is.na(angle.med)&angle.med<0] #fixing negative angles (degrees>180)
#     cc.prox = angle.med
#     #cat(i,"\n")	
#     # saving the results
#     res = cbind.data.frame(locs1$ordem,cc.ex,cc.prox)
#     names(res) = c("ordem","value","value1")
#     #if(is.null(land)) land = res else land = rbind.data.frame(land,res)
#     res
#   }
# })
# stopImplicitCluster()
# stopCluster(cl)
# expo = do.call("rbind.data.frame",tmp2)
# expo = merge(dados,expo,by.x="dados",by.y="ordem",all.x=TRUE) 
# expo = expo[match(dados$dados,expo$dados,nomatch=0),]
# table(expo$dados==dados$dados)
# dados = cbind.data.frame(dados,EXPO=expo$value,EXPO1=expo$value1)

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
dados = cbind.data.frame(dados, ALT = alt$value)   

# ##SLOPE##
# #-Declividade - utilizando a fun??o "slope" GDAL-QGIS
# #0-90 (declividade em graus)
# 
# path0 = "E://ownCloud//W_GIS//Am_Lat_MDT//slope"
# myfiles <- list.files(path = path0, pattern = "tif", full.names = TRUE)
# myfiles <- myfiles[grep("tif.aux.xml|tif.ovr",myfiles,invert=TRUE)]
# 
# #extracting the distance to the cosatline for each survey (paralellized)
# cl <- makeCluster(detectCores()-1)
# clusterEvalQ(cl, library(sp))
# clusterEvalQ(cl, library(raster))
# clusterEvalQ(cl, library(rgdal))
# clusterExport(cl, list("surveys1","myfiles"),envir=environment())
# tmp2 = parLapply(cl,1:length(myfiles), fun= function(i) {
#   # loading the corresponding landscape from the path folder
#   dcl <- try(raster(myfiles[i]),TRUE)
#   # getting the coordinate of the survey
#   ext = round(extent(dcl),0)
#   locs= surveys1[surveys1@coords[,"Latitude"]>=ext@ymin & surveys1@coords[,"Latitude"]<=ext@ymax,]
#   locs= locs[locs@coords[,"Longitude"]>=ext@xmin & locs@coords[,"Longitude"]<=ext@xmax,]
#   # checking if there are coordinates within tile ranges
#   if(length(locs)==0) { 
#     #res = cbind.data.frame(NA,NA,NA)
#     #names(res) = c("ordem","value","value1") 
#   } else { 
#     # homogeneizing projections
#     pj = crs(dcl)
#     locs1 = try(spTransform(locs, pj),TRUE)
#     # extracting the raster values
#     cc.ex = extract(dcl,locs1)
#     # extracting the raster values from nearby
#     cc.prox = extract(dcl,locs1, fun=median, buffer=30*sqrt(2)+(30*sqrt(2)/2))
#     #cat(i,"\n")	
#     # saving the results
#     res = cbind.data.frame(locs1$ordem,cc.ex,cc.prox)
#     names(res) = c("ordem","value","value1")
#     #if(is.null(land)) land = res else land = rbind.data.frame(land,res)
#     res
#   }
# })
# stopImplicitCluster()
# stopCluster(cl)
# decliv = do.call("rbind.data.frame",tmp2)
# decliv = merge(dados,decliv,by.x="dados",by.y="ordem",all.x=TRUE) 
# decliv = decliv[match(dados$dados,decliv$dados,nomatch=0),]
# table(decliv$dados==dados$dados)
# dados = cbind.data.frame(dados,DECLIV=decliv$value,DECLIV1=decliv$value1)
# 
# #Fazendo ultimas corre??es
# dados[!is.na(dados$EXPO)&dados$DECLIV %in% 0,]$EXPO = 0 #deveria ser um NULL mas para evitar erros a exposi??o ficou 0
# 
# 
# ### POTENTIAL ANNUAL DIRECT INCIDENT RADIATION ###
# #From McCune, B. & Keon, D. 2002. Equations for potential annual direct incident radiation and heat load. Journal of vegetation science, 13(4), 603-606.
# 
# #sites at the northern hemisphere
# nort = surveys1$lat1>=0
# #folded aspect for the southern hemisphere
# tmp = 180 - abs(dados$EXPO - 180) #folding "exact" slope aspect
# tmp[!nort] = 180 - tmp[!nort] #converting to south hemisphere
# tmp[dados$DECLIV %in% 0] = 0
# 
# tmp1 = 180 - abs(dados$EXPO1 - 180) #folding "median" slope aspect
# tmp1[!nort] = 180 - tmp1[!nort] #converting to south hemisphere
# tmp1[dados$DECLIV1 %in% 0] = 0
# 
# #getting latitude and slope aspect and declivity (in radians)
# tmp = cbind.data.frame(lat=-surveys1$Latitude,
#                        slope=dados$DECLIV,
#                        aspect=dados$EXPO,
#                        slope.1=dados$DECLIV1,
#                        aspect.1=dados$EXPO1,
#                        lat1=TT.deg2rad(-surveys1$Latitude),
#                        slope1=TT.deg2rad(dados$DECLIV),
#                        aspect1=TT.deg2rad(tmp),
#                        slope2=TT.deg2rad(dados$DECLIV1),
#                        aspect2=TT.deg2rad(tmp1)) 
# 
# #Getting the annual potential direct radiation (MJ/cm2/yr)
# #Eq.1 is for all slopes (<90degrees) and Eq.2 is for slopes <60 degrees 
# table(tmp$slope>=60)
# 
# ##DirRad_Eq2
# #tmp$Eq1.exact = exp(-1.467+1.582*cos(tmp$lat1)*cos(tmp$slope1)-1.5*cos(tmp$aspect1)* sin(tmp$slope1)*sin(tmp$lat1)-0.262*sin(tmp$lat1)*sin(tmp$slope1)+0.607* sin(tmp$aspect1)*sin(tmp$slope1))
# tmp$Eq2.exact = exp(-1.236+1.35*cos(tmp$lat1)*cos(tmp$slope1)-1.376*cos(tmp$aspect1)* sin(tmp$slope1)*sin(tmp$lat1)-0.331*sin(tmp$lat1)*sin(tmp$slope1)+0.375* sin(tmp$aspect1)*sin(tmp$slope1))
# #tmp$Eq1.median = exp(-1.467+1.582*cos(tmp$lat1)*cos(tmp$slope2)-1.5*cos(tmp$aspect2)* sin(tmp$slope2)*sin(tmp$lat1)-0.262*sin(tmp$lat1)*sin(tmp$slope2)+0.607* sin(tmp$aspect2)*sin(tmp$slope2))
# tmp$Eq2.median = exp(-1.236+1.35*cos(tmp$lat1)*cos(tmp$slope2)-1.376*cos(tmp$aspect2)* sin(tmp$slope2)*sin(tmp$lat1)-0.331*sin(tmp$lat1)*sin(tmp$slope2)+0.375* sin(tmp$aspect2)*sin(tmp$slope2))
# 
# #Saving just the two versions of Eq.2, since there is only two sites with slope>=60
# dados = cbind.data.frame(dados,DirRad.ex=tmp$Eq2.exact,DirRad.med=tmp$Eq2.median)
# 
# #removing the objects that will not be used anymore
# rm(dcl,i,land,locs,locs1,myfiles,res,tmp,slp,alti,expo)
# rm(ad,aet,bioc,cwd,dcl,dest,e,ext,filename,foc,i,j,land,legenda,locs,locs1,meses,myfiles,opt1,opt2,
#    opt3,opt4,opt5,tmp,tmp2,tmp3,ur,prec,prox,res,result,rslt,T)
# 
# 
# ### BIORREGIONS ### - TNC
# bioreg <- readOGR(dsn=paste(path,"//Am_Lat_Biorregions_TNC",sep=""),layer="tnc_terr_ecoregions")
# 
# pj = crs(bioreg) # Homogenizing projections
# surveys2 = spTransform(surveys1, pj)
# tmp = over(surveys2,bioreg, fn=NULL)
# dados = cbind.data.frame(dados,bioreg=tmp$ECO_NAME,veget.tnc=tmp$WWF_MHTNAM)
# rm(bioreg)
# 
# 
# ### HUMAN INFLUENCES ### between 1995 - 2004 (0 low and 65 is high)
# hum <- raster(paste(path,"//Am_Lat_Global_Human_Influence_Index_v2//hii_s_amer//w001001.adf",sep=""))
# pj = crs(hum)
# surveys2 = spTransform(surveys1, pj)
# tmp= extract(hum,surveys2)
# tmp1= extract(hum,surveys2,method="bilinear") 
# dados = cbind.data.frame(dados,HumanInfluence=tmp,HumanInfluence1=tmp1)
# rm(hum)

#########################
#### SAVING ENV DATA ####
#########################
names(dados)[1] = "ordem"
dados1 <- merge(occs1, dados, by="ordem")
data.table::setkeyv(dados1, "ordem")
dados1[ , ordem := NULL]
saveRDS(dados1, file = "data/threat_species_by_country_environmental.rds", compress = "gzip")

