###################################################
#### ASSESSING SPECIES EEO AND AOO - CRITERIA B ###
###################################################

### AFAZERES ###
# PASSAR OS ESTUDOS DO CRITERIO B PARA C?
# LER O GUI DO CNCFLORA PARA AS AVALIA??ES: QUAIS COORDENADAS USAR??

# TESTAR ConR USANDO OCORRENCIAS EM UNIDADES DE CONSERVACAO
# TESTAR FUNCAO ADAPTADA USANDO OCORRENCIAS EM UNIDADES DE CONSERVACAO
# Olhar com mais cuidado a lista de species assigned para mim pelo GTA
renato.gta = read.csv("renato_assignments_GTA.csv", as.is= TRUE)

# ADICIONAR UM PASSO? ALL PLOT OCCURRENCES WITHIN SPECIES EOO PRIOR TO OBTAIN AOO?
# GET PREVIOUS ASSESSMENTS: any assessments with more than 10 years old should be re-assessed 

#### LOADING PACKAGES ###
# library(data.table)
# library(rgeos)
# library(rgdal)
# library(redlistr)
# require(maptools)
# require(readxl)
# require(alphahull)
# require(spatstat)
require(conR)
source("suggestions_for_ConR.r")


#### LOADING THE NEOTROPICS MAP ###
neotrop = readRDS("E:/ownCloud/W_GIS/Am_Lat_ADM_GADM_v3.6/gadm36_Neotrop_0_sp_simplified.rds")
neotrop <- gBuffer(neotrop, byid=TRUE, width=0)

## projecting and simplifying the neotropical contours
neotrop.simp <- gSimplify(neotrop,tol=0.05)
neotrop.simp <- gBuffer(neotrop.simp, byid=TRUE, width=0)
sum(gIsValid(neotrop.simp, byid=TRUE)==FALSE)
neotrop.simp.proj <- spTransform(neotrop.simp, CRS("+init=epsg:5641"))  # define projection system of our data



##################################################
#### CALCULATING THE METRICS FOR EACH SPECIES ####
##################################################
## Reading herbarium data
oc.data <- readRDS("threat_occ_data.rds")

#Putting data in the format demanded by the package
MyData <- cbind.data.frame(ddlat = as.double(oc.data$latitude.work1),
                           ddlon = as.double(oc.data$longitude.work1),
                           tax = as.character(oc.data$species.correct2),
                           higher.tax.rank = oc.data$family.correct1,
                           coly = as.double(oc.data$ano),
                           vouchers = oc.data$dup.ID1, # change for dup.ID1
                           detBy = oc.data$determinador.name,
                           dety = oc.data$ano.det,
                           tax.check2 = oc.data$tax.check2,
                           UC = oc.data$UC,
                           stringsAsFactors = FALSE)
#ids <- grepl("Myrtaceae", MyData$higher.tax.rank)
#ids = grep("Psidium ubatubense",MyData$tax)
#MyData <- MyData[ids,] 

#### OBTAINING THE NAME AND NUMBER OF OCCURRENCES (TOTAL AND SPATIALLY-UNIQUE) PER SPECIES ####
resultado <- readRDS("assess_iucn_spp.rds")


#### CLASSIFYING SPECIES ACCORDING TO THE NUMBER OF OCCURRENCES AT EACH TAX. CONFIDENCE LEVEL ####
tmp <- as.data.frame.matrix(table(MyData$tax, MyData$tax.check2))
tmp$all = apply(tmp[,c("cannot_check","FALSE","TRUE","TRUE_OTHER","TRUE_TBC")], 1, sum)
tmp$false.true.true1 = apply(tmp[,c("FALSE","TRUE","TRUE_OTHER","TRUE_TBC")], 1, sum)
tmp$true.true1 = apply(tmp[,c("TRUE","TRUE_OTHER","TRUE_TBC")], 1, sum)
tmp$true.true2 = apply(tmp[,c("TRUE","TRUE_TBC")], 1, sum)
tmp$true = tmp[,c("TRUE")]
tmp <- tmp[,6:10]
tmp$class = NA
tmp$class[tmp$true >=75] <- "true_only"
tmp$class[is.na(tmp$class) & tmp$TRUE_TBC > 0] <- "true_tbc"
# tmp$class[is.na(tmp$class) & tmp$TRUE_TBC > 0 & tmp$true.true2 <75] <- "true_tbc"
# tmp$class[is.na(tmp$class) & tmp$TRUE_TBC > 0 & tmp$true.true2 >=75] <- "true_tbc_sample"
tmp$class[is.na(tmp$class) & tmp$true>= 30] <- "true_other"
tmp$class[is.na(tmp$class) & tmp$true>= 15] <- "true_false"
tmp$class[is.na(tmp$class) & tmp$true< 15] <- "all"


tmp$class[is.na(tmp$class) & tmp$true>= 30 & tmp$TRUE_OTHER<30] <- "true_other"
tmp$class[is.na(tmp$class) & tmp$true.true1 >= 75 & tmp$true >= 30] <- "true_other"
tmp$class[is.na(tmp$class) & tmp$true.true1 >= 75 & tmp$true < 30] <- "true_other"

tmp$class[is.na(tmp$class) & tmp$true > 30 & (tmp$TRUE_OTHER<50)] <- "true_other"
tmp$class[is.na(tmp$class) & tmp$true > 30 & (tmp$TRUE_OTHER>=50)] <- "true_other_sample"
tmp$class[is.na(tmp$class) & tmp$all < 10] <- "all"

resultado <- cbind.data.frame(resultado, tmp, stringsAsFactors = FALSE)
#table(resultado$tax == names(table(MyData$tax[ids])))
resultado$non.dup.occurs.tax <- table(MyData$tax[grepl("TRUE", MyData$tax.check2)])

## Classifying the number of spatially unique occurrences
resultado$group = findInterval(resultado$non.dup.occurs, c(0,3,10,15,75,5000))

resultado$group = NA
resultado$group[resultado$non.dup.occurs<3] = "less.than.three"
resultado$group[is.na(resultado$group) & resultado$non.dup.occurs>=3 & resultado$non.dup.occurs<10] = "three.to.ten"
resultado$group[is.na(resultado$group) & resultado$non.dup.occurs>=10 & resultado$non.dup.occurs<15] = "ten.to.fiften"
resultado$group[is.na(resultado$group) & resultado$non.dup.occurs>=15] = "more.than.fiften"
table(resultado$group)

####################################
#### EXTENT OF OCCURRENCE (EOO) ####
####################################

#### REMOVING MISSING COORDINATES ####
## Doing it only for distribution analyses 
# oc.data$longitude.work1[is.na(oc.data$longitude.work1)] = 
#   oc.data$longitude.work[is.na(oc.data$longitude.work1)]
# oc.data$latitude.work1[is.na(oc.data$latitude.work1)] = 
#   oc.data$latitude.work[is.na(oc.data$latitude.work1)]
# oc.data <- oc.data[!is.na(latitude.work1) | !is.na(longitude.work1),] 
# oc.data <- oc.data[!latitude.work1 %in% "no_coord" | !longitude.work  %in% "no_coord",] 
# range(oc.data$latitude.work1)
# range(oc.data$longitude.work1)

## Convex Hull method
EOO.hull <- ConR:::EOO.computing(MyData[ids & grepl("TRUE", MyData$tax.check2),], export_shp = FALSE,
                             exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                             write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                             write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                             parallel = FALSE, NbeCores = 6) # run on parallel? How many cores?
#plot(EOO.hull$spatial.polygon_2, col="red")
#plot(neotrop.simp, add=T)

## Alpha-Hull method
EOO.alpha <- ConR:::EOO.computing(MyData[ids & grepl("TRUE", MyData$tax.check2),], export_shp = FALSE,
                          method.range =  "alpha.hull", alpha = 2, # maybe 2 is a better choice
                          exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                          write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                          write_results=FALSE, file.name = "EOO.alpha", # If TRUE, a csv fiel is created in the working directory
                          parallel = TRUE, NbeCores = 6) # run on parallel? How many cores?

#### SUBPOPULATIONS ####
#This function applies the method called circular buffer method for estimating the number of subpopulation (Rivers et al., 2010).
#From Rivers et al. 2010: The ideal radius of the buffer is debatable; however when dispersal
#characteristics of the species are unknown then a sliding scale, such as the 1/10th maximum inter-point distance, is the preferred choice, as it is species-specific and not sensitive to collection density
radius = subpop.radius(MyData[ids,], factor.div = 10, quant.max = 0.95)
SUB <- my.subpop.comp(MyData[ids,], Resol_sub_pop=30)
SUB1 <- my.subpop.comp(MyData[ids,], Resol_sub_pop= radius)
sapply(SUB, function(x) x[[1]])
sapply(SUB1, function(x) x[[1]])

#### AREA OF OCCUPANCY (AOO) ####
AOO = AOO.computing(MyData[ids,], Cell_size_AOO = 2, # Size of the grid cell in kilometers
                    nbe.rep.rast.AOO = 2, # number of raster with random starting position for estimating the AOO
                    parallel = FALSE, NbeCores = 6, 
                    show_progress= TRUE, export_shp=FALSE)

#### NUMBER OF LOCATIONS ####
locs = my.locations.comp(MyData[ids,], 
                                method = "fixed_grid",
                                nbe_rep = 0, # number of raster with random starting position
                                Cell_size_locations = 2, #grid size in kilometres used for estimating the number of location
                                Rel_cell_size=0.05,
                                protec.areas = NULL, #a SpatialPolygonsDataFrame, shapefile with protected areas.
                                ID_shape_PA= "WDPA_PID", # field name of protec.areas with ID
                                method_protected_area="no_more_than_one", 
                                parallel=FALSE, NbeCores=2)

#### IUCN ASSESSMENTS ####
# putting the results in the good format for evaluation
resultado.hull = rbind.data.frame(EEO = EOO.hull$EOO, 
                              AOO = as.double(AOO),
                              Nbe_unique_occ. = resultado$non.dup.occurs,
                              Nbe_subPop = as.double(sapply(SUB, function(x) x[[1]])),  
                              Nbe_loc = as.double(locs[[2]])
                              #Nbe_loc_PA = 
                              )
colnames(resultado.hull) = resultado$tax
rownames(resultado.hull) = c("EOO", "AOO", "Nbe_unique_occ.", "Nbe_subPop", "Nbe_loc")
DATA = resultado.hull

# obtaining the results using 
df = data.frame(ID=1:length(neotrop.simp)) 
neotrop.simp.spdf = SpatialPolygonsDataFrame(neotrop.simp, df)
IUCN = IUCN.eval(MyData[ids,], country_map = neotrop.simp.spdf, Cell_size_AOO = 2, Cell_size_locations = 2, 
                 Resol_sub_pop = 30, method_locations = "fixed_grid", Rel_cell_size = 0.05, 
                 DrawMap = FALSE, add.legend = TRUE, 
                 file_name = NULL, export_shp = FALSE, write_shp = FALSE, 
                 write_results = FALSE, protec.areas = NULL, map_pdf = FALSE, draw.poly.EOO=TRUE,
                 exclude.area = TRUE, method_protected_area = "no_more_than_one", 
                 ID_shape_PA = "WDPA_PID", 
                 buff_width = 0.1, SubPop=TRUE, alpha=1, buff.alpha=0.1, 
                 method.range="convex.hull", nbe.rep.rast.AOO=0,
                 showWarnings=TRUE, write_file_option="excel", 
                 parallel=TRUE, NbeCores=6)

criteria.B(resultado.hull, poly_borders = NULL, 
                       protec.areas = NULL, 
                       file_name = NULL, add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                       EOO.threshold = c(20000, 5000, 100), AOO.threshold = c(2000, 500, 10), Loc.threshold = c(10, 5, 1))


