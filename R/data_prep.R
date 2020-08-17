#####################################################
#####################################################
#### PREPARING DATA FROM THREAT IUCN ASSESSMENTS ####
#####################################################
#####################################################
#rm(list=ls())

#### LOADING PACKAGES ####
library(data.table)
library(rgeos)
library(dplyr)
source("R/suggestions_for_ConR.r")

########################################H
#### SIMPLIFIED NEOTROPICAL COUNTOUR ####
########################################H
neotrop <- readRDS("E:/ownCloud/W_GIS/Am_Lat_ADM_GADM_v3.6/gadm36_Neotrop_0_sp.rds")
neotrop <- gBuffer(neotrop, byid=TRUE, width=0)

#projecting and simplifying the neotropical contours
neotrop.simp <- gSimplify(neotrop, tol=0.0001, topologyPreserve = TRUE)
neotrop.simp <- gBuffer(neotrop.simp, byid=TRUE, width=0)
neotrop.simp1 <- gSimplify(neotrop, tol=0.001, topologyPreserve = TRUE)
neotrop.simp1 <- gBuffer(neotrop.simp1, byid=TRUE, width=0)
neotrop.simp2 <- gSimplify(neotrop, tol=0.005, topologyPreserve = TRUE)
neotrop.simp2 <- gBuffer(neotrop.simp2, byid=TRUE, width=0)
neotrop.simp3 <- gSimplify(neotrop, tol=0.01, topologyPreserve = TRUE)
neotrop.simp3 <- gBuffer(neotrop.simp3, byid=TRUE, width=0)

#how many polygons left?
length(neotrop)
length(neotrop.simp)
length(neotrop.simp1)
length(neotrop.simp2)
length(neotrop.simp3)

# any bad polygons remaining?
sum(gIsValid(neotrop.simp, byid=TRUE)==FALSE) #no!
sum(gIsValid(neotrop.simp1, byid=TRUE)==FALSE) #no!
sum(gIsValid(neotrop.simp2, byid=TRUE)==FALSE) #no!
sum(gIsValid(neotrop.simp3, byid=TRUE)==FALSE) #no!

#removing very small polygons
length(neotrop.simp2)
area <- gArea(neotrop.simp2, byid=TRUE) 
neotrop.simp4 <- neotrop.simp2[area > 0.004*0.004,] #removing polygons smaller than ~16 ha/0.16 km2
length(neotrop.simp4)
 

#Inspecting more closely...
plot(neotrop.simp1, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=1)
plot(neotrop.simp2, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=2)
plot(neotrop.simp3, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=3)
plot(neotrop.simp4, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=4)

#Saving
#saveRDS(neotrop.simp, file = "data//Contour_Neotrop_simplified_very_large.rds")
saveRDS(neotrop.simp1, file = "data//Contour_Neotrop_simplified_tol_001.rds")
saveRDS(neotrop.simp2, file = "data//Contour_Neotrop_simplified_tol_005.rds")
saveRDS(neotrop.simp3, file = "data//Contour_Neotrop_simplified_tol_01.rds")
saveRDS(neotrop.simp4, file = "data//Contour_Neotrop_simplified_tol_005_no_small.rds")


#####################################################################################################################################################################h
#####################################################################################################################################################################h
##################################H
#### HERBARIUM OCCURRENCE DATA ####
##################################H

#### LOADING OCCURRENCE DATA ###
## Loading the paths to occurrence data 
paths = dir("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//occurrence_data",full.names=TRUE)
paths = paths[grepl('merged_outliers.csv',paths)]

## Loading the different collections
lista = vector("list",length(paths))
for (i in 1:length(paths)) {
  path.csv = paths[i]
  #tmp = read.csv(path.csv,as.is=TRUE,na.strings=c(""," ","NA"))
  tmp = fread(path.csv,na.strings=c(""," ","NA"))
  lista[[i]] = tmp
}
## Merging the different collections
oc.data = rbindlist(lista) #merging the unicate and duplicated tables
rm(lista, tmp, paths, i)
gc()

#################################H
#### MANAGING OCCURRENCE DATA ####
#################################H

#### REMOVING SPECIES NOT IN THE ATLANTIC FOREST CHECKLIST ###
## Creating a vector with the names of specimens determined at species and infra-specific levels
oc.data[,species.correct2 := sapply(strsplit(species.correct1," "), function(x) paste(x[1], x[2],sep=" ")),]
setkeyv(oc.data, "species.correct2") ## setting 'species.correct2' as key to the data.table (makes computations faster)

spp.af = readRDS("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/spp.af.rds")
#Removing Unresolved, Uncertain application and etc
spp.af <- spp.af[!spp.af$TreeCo_status %in% c("remove","remove?","remove? Unresolved","Incorrect","Missapplied","Illegitimate","Illegitimate name","Uncertain application","Unresolved","not in the AF"),]
#names with possible problems 
#sort(table(spp.af$TreeCo_status))
#names confirmed for the Atlantic Forest
spp.af.true <- spp.af[af.check2 %in% "TRUE", species.correct2]
#Adding some missing names (don't know why they were excluded; probably due to missing valid determinations or coordinates)
miss.spp = c("Agonandra brasiliensis","Aiouea bracteata","Aspidosperma brasiliense",
             "Aspidosperma quirandy","Brasiliopuntia schulzii","Cochlospermum vitifolium","Derris leucanthus","Ficus pakkensis","Handroanthus selachidentatus",
             "Handroanthus spongiosus","Ilex cognata","Luetzelburgia purpurea","Mezilaurus synandra","Mimosa caesalpiniifolia",
             "Ocotea grandifructa","Palicourea subcuspidata","Palicourea didymocarpos","Persea punctata","Piptadenia killipii","Piranhea securinega",
             "Quillaja lancifolia","Quillaja lanceolata","Quillaja sellowiana",
             "Tovomita longifolia","Trischidium decipiens","Zanthoxylum unifoliolatum")
oc.data <- oc.data[species.correct2 %in% c(spp.af.true, miss.spp),] # occurrences inside the AF
#oc.data[,uniqueN(species.correct2)] # 5183 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### SOME LAST MINUTE, NEW SYNONYMS TO BE REPLACED ###
oc.data[species.correct2 %in% "Brasiliopuntia schulzii", status := "ok"]
oc.data[species.correct2 %in% "Brasiliopuntia schulzii", species.correct2 := "Brasiliopuntia brasiliensis"]
oc.data[species.correct2 %in% "Derris leucanthus", status := "ok"]
oc.data[species.correct2 %in% "Derris leucanthus", species.correct2 := "Muellera campestris"]
spp2rep <- c("Palicourea subcuspidata","Palicourea didymocarpa","Palicourea didymocarpos")
oc.data[species.correct2 %in% spp2rep, status := "ok"]
oc.data[species.correct2 %in% spp2rep, species.correct2 := "Palicourea didymocarpa"]
spp2rep <- c("Quillaja brasiliensis","Quillaja lancifolia","Quillaja lanceolata","Quillaja sellowiana")
oc.data[species.correct2 %in% spp2rep, species.correct2 := "Quillaja lancifolia"]
oc.data[loc.correct %like% "brazil" & species.correct2 %in% "Ixora grandifolia", species.correct2 := "Ixora muelleri"]
oc.data[species.correct2 %in% "Misanteca duartei", species.correct2 := "Licaria guianensis"]
oc.data[species.correct2 %like% "Zanthoxylum" & numTombo %in% c("jpb_54464","ase_6935","ase_29975","ase_30740","ase_35130"), species.correct2 := "Zanthoxylum unifoliolatum"]
oc.data[species.correct2 %in% c("Calyptranthes grandifolia","Calyptranthes brasiliensis"), species.correct2 := "Calyptranthes brasiliensis"]
oc.data[species.correct2 %in% c("Ossaea loligomorpha","Miconia loligomorpha"), species.correct2 := "Leandra loligomorpha"]
oc.data[species.correct2 %in% c("Plinia brachybotrya","Plinia pseudodichasiantha"), species.correct2 := "Plinia pseudodichasiantha"]
#oc.data[,uniqueN(species.correct2)] # 5178 species

#### REMOVING SPECIES THAT SHOULD NOT BE IN THE LIST (EXOTICS, CULTIVATED OR NOT IN THE AF) ###
rm.spp <- c("Annona calophylla","Bauhinia galpinii","Bunchosia glandulifera","Cereus repandus","Eugenia acapulcensis",      
            "Garcinia longifolia","Grabowskia boerhaaviifolia","Humiriastrum cuspidatum","Ixora grandifolia","Maytenus macrocarpa","Miconia pausana",           
            "Mimosa caesalpiniifolia","Mimosa rufescens","Myrceugenia obtusa","Myrsine ligustrina","Ocotea megaphylla", 
            "Prunus oleifolia","Randia obovata","Rauvolfia tetraphylla","Schinus fasciculata","Senegalia fiebrigii",
            "Sesbania macroptera","Sloanea multiflora","Solanum confusum","Solanum lanceifolium","Tarenaya spinosa",
            "Vachellia macracantha","Zanthoxylum schreberi")
oc.data <- oc.data[!species.correct2 %in% rm.spp,] # occurrences inside the AF

### Removing all woody bamboos - No abudance data ###
oc.data <- oc.data[!family.correct1 %in% "Poaceae",]
#oc.data[,uniqueN(species.correct2)] # 5101 species

#### REMOVING DUPLICATES ###
gc()
toto1 <- dim(oc.data)[1]
## Ordering the dataset to more-easily remove duplicates
#first by dup.ID, then decreasing dup.numb, then decreasing dup.prop, then decreasing source name (preferenc to splink, than jabot, then splink_new, then gbif)
oc.data[, source1 := as.double(as.character(factor(source, levels = c("gbif","sndb","jabot","splink","splink_new"), labels = c(5,4,2,1,3)))), ]
setorder(oc.data,  dup.ID, -dup.numb, -dup.prop, source1)               ## re-ordering the data.table, using setorder()
oc.data[,source1:=NULL] #removing the extra column created for ranking
#creating the unique accession number and multiple-acession number for each specimen
oc.data[, dup.ID1 := dup.ID]
oc.data[is.na(dup.ID1), dup.ID1 := numTombo, by = numTombo]
#removing the duplicates
oc.data <- unique(oc.data, by = "dup.ID1")
#removing the extra column created for ranking
#oc.data[,dup.ID1:=NULL] 
toto1 - dim(oc.data)[1]; 100*(toto1 - dim(oc.data)[1])/toto1 ## 730,297 (23.68%) removed
# oc.data[,uniqueN(species.correct2)] # 5101 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING DATA NOT GEOGRAPHICALLY VALIDATED ###
## But keeping records for species without coordinates but know to be in the Atlantic Forest
gc()
toto2 = dim(oc.data)[1]
oc.data <- oc.data[geo.check1 %like% "ok_county|ok_locality" | 
                     (geo.check1 %like% "ok_state" & af.check2 == TRUE)]
#oc.data[,uniqueN(species.correct2)] # 5101 species
#table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING MISSING COORDINATES ###
## Doing it after, only before distribution analyses 
oc.data$longitude.work1[is.na(oc.data$longitude.work1)] <-
        oc.data$longitude.work[is.na(oc.data$longitude.work1)]
oc.data$latitude.work1[is.na(oc.data$latitude.work1)] <-
        oc.data$latitude.work[is.na(oc.data$latitude.work1)]
oc.data$longitude.work1[oc.data$longitude.work1 %in% "no_coord" | oc.data$latitude.work1 %in% "no_coord"] <- NA
oc.data$latitude.work1[oc.data$longitude.work1 %in% "no_coord" | oc.data$latitude.work1 %in% "no_coord"] <- NA
# oc.data <- oc.data[!is.na(latitude.work1) | !is.na(longitude.work1),]
# range(as.double(oc.data$latitude.work1), na.rm = TRUE)
# range(as.double(oc.data$longitude.work1), na.rm = TRUE)
toto2 - dim(oc.data)[1]; 100*(toto2 - dim(oc.data)[1])/toto2 ## 961,869 (40.87% of the non-duplicated) removed
100*(toto2 - dim(oc.data)[1])/toto1 ## extra 31.19% in respect to all records
#oc.data[,uniqueN(species.correct2)] # 5101 species
#table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING SPATIAL DUPLICATES ###
gc()
toto3 = dim(oc.data)[1]
setkeyv(oc.data, "species.correct2") ## setting 'species.correct2' as key to the data.table (makes computations faster)
oc.data[ , coord.string := as.character(paste(latitude.work1, longitude.work1, sep="_")),]
oc.data[ , spat.dups := duplicated(coord.string, incomparables = "NA_NA"), by = species.correct2]
oc.data <- oc.data[spat.dups == FALSE,,]
oc.data[,coord.string:=NULL] 
oc.data[,spat.dups:=NULL] 
toto3 - dim(oc.data)[1]; 100*(toto3 - dim(oc.data)[1])/toto2 ## 572,540 (24.33% of remaining records) removed
100*(toto3 - dim(oc.data)[1])/toto1 ## extra 18.57% in respect to all records
# oc.data[,uniqueN(species.correct2)] # 5101 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING SPATIAL OUTLIERS  ###
## GETTING SPECIES INFORMATION: TBC, geographical range and cultivation #### Getting species information and generating the vectors for the groups of species
spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
tmp = unique(oc.data[,species.correct2,])
spp = spp[spp$SpeciesKTSA %in% tmp,]
uso = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species Uses//plant_uses.csv", na= c(""," ",NA))
uso = uso[as.character(uso$Name_submitted) %in% tmp,]

## Preliminary editing of the list
spp$nomes = spp$species.reflora
#replacing names not found in ReFlora by the TNRS/Hans taxonomy
spp$nomes[is.na(spp$nomes)] = spp$Accepted_name[is.na(spp$nomes)]
#replacing names not found in ReFlora and TNRS by the Plant List
spp$nomes[is.na(spp$nomes)] = spp$NewTaxon.TPL[is.na(spp$nomes)]

## Getting the list of species per group
#tbc = spp$nomes[spp$TBC %in% "TBC"]
cult= as.character(uso$Name_submitted[as.character(uso$group_renato) %in% "cultivated"])

## REMOVING TRUE OUTLIERS ##
#1- remover true outliers para (todas?) as especies: classic Mahalanobis distance >3 & robust (Mahalanobis) distance>16
#table(oc.data$true.out, useNA = "always")
gc()
toto4 = dim(oc.data)[1]
oc.data <- oc.data[is.na(true.out) |true.out %in% FALSE] # removing 3117 true outliers
# oc.data[,uniqueN(species.correct2)] # 5101 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

## REMOVING PROBABLE OUTLIERS ##
gc()
#2- remover probable outliers para cultivated species:  classic Mahalanobis distance >2.5 & robust (Mahalanobis) distance>11.5
tmp <- data.table(oc.data[species.correct2 %in% cult])
#table(tmp$probable.out, useNA = "always")
tmp1 <- data.table(tmp[probable.out %in% TRUE])
oc.data <- oc.data[!numTombo %in% tmp1$numTombo]
toto4 - dim(oc.data)[1]; 100*(toto4 - dim(oc.data)[1])/toto3 ## 3170 (0.23% of remaining records) removed
100*(toto4 - dim(oc.data)[1])/toto1 ## extra 0.1% in respect to all records
# oc.data[,uniqueN(species.correct2)] # 5101 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#3- remover non-core specimens for heavily cultivated native species (e.g. Araucaria angustifolia): rob.out.99 == FALSE
## REMOVE OCCURRENCES OUTSIDE THE NEOTROPICS, EXCEPT FOR PANTROPICAL SPECIES?
## REMOVE OCCURRENCE DATA FROM WELL-KNOW BOTANICAL GARDENS GROWING SOUTH AMERICA SPECIES?

#### REMOVING NON-NEOTROPICAL OCCURRENCES ###
gc()
toto5 = dim(oc.data)[1]
# Finding possible unassigned non-neotropical occurrences
non.neotrop <- table(oc.data[neotrop.check == FALSE | is.na(neotrop.check), unique(country)])
non.neotrop <- data.frame(code = names(non.neotrop), n.occ = as.double(non.neotrop), county = names(non.neotrop),
                          stringsAsFactors = FALSE)
non.neotrop$county[nchar(non.neotrop$county) == 2] <- 
  countrycode::countrycode(non.neotrop$county[nchar(non.neotrop$county) == 2], "iso2c", "country.name")
non.neotrop <- non.neotrop[!non.neotrop$code %in% c("Brazil","Brasil","BRASIL","BR","BB","BO","BS","BZ",
                                                    "CL","CO","CR","CU","DM","DO","GF","GY","JM","MQ",
                                                    "Mexico","MX","NI","PA","PE","PR","SR","SV",
                                                    "US","VE","VI","ZZ"),]
oc.data[country %in% non.neotrop$code & is.na(neotrop.check), neotrop.check := FALSE]
yes.neotrop <- table(oc.data[is.na(neotrop.check), unique(country)])
yes.neotrop <- data.frame(code = names(yes.neotrop), n.occ = as.double(yes.neotrop), county = names(yes.neotrop),
                          stringsAsFactors = FALSE)
yes.neotrop$county[nchar(yes.neotrop$county) == 2] <- 
  countrycode::countrycode(yes.neotrop$county[nchar(yes.neotrop$county) == 2], "iso2c", "country.name")
yes.neotrop <- yes.neotrop[!yes.neotrop$code %in% c("ZZ"),]
oc.data[country %in% yes.neotrop$code & is.na(neotrop.check), neotrop.check := TRUE]
#Getting the number of non-neotropical occurrences per species
extra.neotrop.spp = table(oc.data[neotrop.check == FALSE & af.check2 %in% c("cannot_check"), species.correct2])
oc.data <- oc.data[neotrop.check == TRUE | is.na(neotrop.check) |(neotrop.check == FALSE & af.check2 %in% c("TRUE","FALSE"))]
toto5 - dim(oc.data)[1]; 100*(toto5 - dim(oc.data)[1])/toto4 ## 485 (0.11% of remaining records) removed
100*(toto5 - dim(oc.data)[1])/toto1 ## extra 0.02% in respect to all records
#oc.data[,uniqueN(species.correct2)] # 5101 species
#table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### GETTING SPECIES INFORMATION: TBC, geographical range and cultivation ### 
# Getting species information and generating the vectors for the groups of species
spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
spp = spp[spp$SpeciesKTSA %in% oc.data$species.correct2,]
# uso = read_excel("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species Uses//plant_uses.xlsx", na= c(""," ",NA))
# uso = uso[uso$Name_submitted %in% oc.data$species.correct2,]
# 
# ## Preliminary editing of the list
spp$nomes = spp$species.reflora
#replacing names not found in ReFlora by the TNRS/Hans taxonomy
spp$nomes[is.na(spp$nomes)] = spp$Accepted_name[is.na(spp$nomes)]
#replacing names not found in ReFlora and TNRS by the Plant List
spp$nomes[is.na(spp$nomes)] = spp$NewTaxon.TPL[is.na(spp$nomes)]

## Getting the list of species per group
tbc = spp$nomes[spp$TBC %in% "TBC"]
# cult= uso$Name_submitted[uso$group_renato == "cultivated"]

#### SEPARATING DATA IDENTIFIED AND NOT IDENTIFIED BY TAXONOMISTS ###
## FLAGGING:
#valid/TRUE: type specimen, identification validated by the family specialist, or collector is the specialist of the same family 
#probably_valid/TRUE_...: identifications made by plant taxonomists or TBC (taxa with low taxonomic complexity)
#invalid/FALSE: determiner is somebody else (probably not a taxonomist) 
#cannot_check: empty or anonymous identifiedBy

#Creating a vector with the names of specimens determined at species and infra-specific levels
oc.data[,tax.check2 := as.character(tax.check1),]

#Validating all occurrences identified by other-families taxonomists
tmp = data.table(oc.data[tax.check %in% "TRUE" & is.na(typeStatus),])
tmp = table(tmp$determinador.name)
#### CHECK HERE: what should be the best cut-off to set the most frequent taxonomists? Set to the 1 quartile? Leaving the 20 that was fixed before... ####
taxonomists = names(tmp[tmp >= 20]) # quantile(tmp, prob = seq(0.2,0.5,0.05)) 
taxonomists = taxonomists[!taxonomists %in% c("Anonimo","Anonymous","Semdeterminador")]
oc.data[determinador.name %in% taxonomists & tax.check2 %in% "FALSE", tax.check2 := "TRUE_OTHER",]

#Validating all occurrences of TBC
oc.data[species.correct1 %in% tbc & tax.check2 %in% c("FALSE","cannot_check"), tax.check2 := "TRUE_TBC",]
#oc.data[,uniqueN(species.correct2)] # 5101 species
#table(is.na(oc.data$latitude.work1))

#some last minute validation to avoid the removal of species from the checklist
oc.data[species.correct2 %in% "Agonandra brasiliensis" & determinador.name1 %like% c("Hiepko|Heipko|Marquete|Groppo"), tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Aiouea bracteata" & determinador.name1 %like% c("Baitello|Lorea"), tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Quillaja lancifolia" & determinador.name1 %in% c("Luebert, F."), tax.check2 := "TRUE" ]
oc.data[species.correct2 %like% "Handroanthus" & determinador.name1  %like% c("Santo, F.S.E.|Espirito Santo, F.S.|Santo, F.E.|Silva-Castro|ilva, M.M."), tax.check2 := "TRUE" ]

#### TRYING TO OBTAIN INFORMATION ON CONSERVATION UNITS FROM LOCALITY INFO ###
locals <- oc.data$locality
#Editing the localities
unwanted_array = list('Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Ç'='C', 'È'='E', 'É'='E',
                      'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                      'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 'ß'='S', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c',
                      'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                      'ö'='o', 'ø'='o', 'ü'='u', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )
locals <- chartr(paste(names(unwanted_array), collapse=''), paste(unwanted_array, collapse=''), locals)
locals <- gsub('\\.\\.',".", locals)
locals <- gsub(' - | -|- ',"-", locals)
locals <- gsub('^-- ',"", locals)
locals <- gsub(' --$',"", locals)
locals <- gsub('^\\.|^,',"", locals)
locals <- gsub('\\.$|,$',"", locals)
locals <- gsub('^-\\. ','', locals)
locals <- gsub('-\\.$','', locals)
locals <- gsub('^-','', locals)
locals <- gsub('-$','', locals)
locals <- gsub("     "," ", locals)
locals <- gsub("    "," ", locals)
locals <- gsub("   "," ", locals)
locals <- gsub("  "," ", locals)
locals <- stringr::str_trim(locals)
oc.data$locality1 <- locals
# locals <- table(locals)
# tmp <- data.frame(locality = names(locals), records = as.double(locals))
# tmp <- tmp[order(tmp$records, decreasing = TRUE),]
# tmp <- tmp[!duplicated(tmp$locality),]
# PI <- c("Parque Estadual|Parque Nacional|Parque Municipal|PARNA |P\\.E\\.|P\\.N\\.|Estacao Ecologica|Estacao Biologica|E\\.B\\.|Reserva Biologica|REBIO |R\\.E\\.B\\.I\\.O\\.|Monumento Natural|Reserva da Vida Silvestre|RVS |Santuario Ecologico|Refugio Biologico|Reserva Ecologica|Reserva Estadual|Parque do Estado|CPCN Pro-Mata|Ilha do Cardoso|Ilha Grande|Parque Natural Municipal|PArque Ecologico|Parque do Sabia|Nucleo Picinguaba/PESM")
# US <- c("Area de Relevante Interesse Ecologico|ARIE |Floresta Nacional|FLONA |Estacao Experimental|ESEX |Jardim Botanico|Zoologico|Zoobotanico|Horto Florestal|Refugio Biologico|Reserva Extrativista|RESEX |Reserva de Fauna|Reserva da Fauna|Reserva de Desenvolvimento Sustentavel|RDS |Reserva Particular do Patrimonio Natural|RPPN |R\\.P\\.P\\.N\\.|Reserva Florestal|Reserva Municipal|Serra da Cantareira|Fazenda de Ensino e Pesquisa|AC, Fazenda Santa Elisa|Fazenda Agua Limpa|Serra do Japi|Cidade Universitaria Armando de Salles Oliveira")
# tmp$PI <- grepl(PI, tmp$locality, ignore.case = TRUE) & !grepl("ntorno d|roximo ao|caminho d", tmp$locality, ignore.case = TRUE)
# tmp$US <- grepl(US, tmp$locality, ignore.case = TRUE) & !grepl("ntorno d|roximo ao|caminho d", tmp$locality, ignore.case = TRUE)
# tmp$APA <- grepl("Area de Protecao Ambiental|APA |A\\.P\\.A\\.", tmp$locality) & !grepl("ntorno do|roximo ao|caminho do", tmp$locality, ignore.case = TRUE)
# write.csv(tmp, "UCs.csv")
ucs <- read.csv("UCs.csv", as.is = TRUE)
oc.data <- dplyr::left_join(oc.data, ucs, by= "locality1")

#### SAVING THE ESSENTIAL DATA FOR THE ASSESSMENTS ###

## Filtering only the necessary columns
oc.data1 <- oc.data[,c("latitude.work1","longitude.work1","species.correct2","family.correct1",
         "coletor.name","coletor.number","ano","numTombo","dup.ID1",
         "typeStatus","determinador.name1","ano.det1","tax.check2",
         "geo.check1","origin.coord1","af.check2","UC"),]
saveRDS(oc.data1, file = "data/threat_occ_data.rds", compress = "gzip")

#### OBTAINING THE NAME AND NUMBER OF OCCURRENCES (TOTAL AND SPATIALLY-UNIQUE) PER SPECIES ###
resultado <- oc.data[!duplicated(oc.data$species.correct2), c("family.correct1", "species.correct2")]
#table(resultado$species.correct2 == names(table(oc.data$species.correct2)))
resultado$total.occs <- as.double(table(oc.data$species.correct2))
tmp <- merge(resultado, 
             data.frame(species.correct2 = names(extra.neotrop.spp), extra.neotrop.occurs = as.double(extra.neotrop.spp)), 
                        by = "species.correct2", all.x = TRUE, sort = FALSE)
tmp <- tmp[order(tmp$species.correct2),]
#table(resultado$species.correct2 == tmp$species.correct2)
resultado$non.neotrop.occs <- tmp$extra.neotrop.occurs
MyData <- cbind.data.frame(ddlat = as.double(oc.data$latitude.work1),
                           ddlon = as.double(oc.data$longitude.work1),
                           tax = as.character(oc.data$species.correct2))
MyData <- MyData[!is.na(MyData$ddlat) | !is.na(MyData$ddlon),]
tmp <- unique.coords(MyData)
#table(resultado$species.correct2 == tmp$tax)
resultado$non.dup.occs <- tmp$non.dup.occurs
resultado$total.occs <- as.double(apply(resultado[,c("total.occs","non.neotrop.occs")], 1, sum, na.rm = TRUE))

#####################################################################################################################################################################H
#####################################################################################################################################################################H
#######################H
#### INVENTORY DATA ####
#######################H
rm(list=ls()[-which(ls() == "resultado")])

##Loading functions
path = "C:/Users/renato/Documents/raflima/Pos Doc/Databases/TreeCo Database Management"
source(paste(path,"references_data_prep.R",sep="/"))
#source(paste(path,"abundances_data_prep.R",sep="/"))

## EDITING TREECO SITE AND ABUNDANCE DATA ##
#Loading the data
surveys = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//References//references.csv",as.is=TRUE, encoding = "UTF-8")
names(surveys)[1] = "refID"
trees = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species abundances//abundances.csv",as.is=TRUE, na.string=c(""," ",NA))

# Running the pre-processing of the dataset
trees <- trees[!is.na(trees$species.correct),]
trees <- trees[!grepl("_",trees$SiteCode),]
sites <- unlist(strsplit(unique(trees$ordem),"\\|"))
surveys1 = survey_data_prep(surveys,linhas=sites)$data_frame

## READING AND EDITING THE ABUNDANCE AND PLOT DATA ##
#trees1 <- abundance_data_prep(trees, surveys1, linhas= sites)$abundances_per_ordem
trees1 <- aggregate(trees$N, list(trees$ordem, trees$species.correct), sum, na.rm=TRUE)
names(trees1)[1:3] = c("ordem","species.correct", "N")
trees1$RecordID <- aggregate(trees$RecordID, list(trees$ordem, trees$species.correct), function(x) x[1])$x
tmp <- aggregate(trees$family, list(trees$species.correct), unique)
names(tmp) = c("species.correct","family")
trees1 <- left_join(trees1, tmp, by = "species.correct")
trees1 <- trees1[,c("RecordID","ordem","species.correct","family","N")]

##Binding other important info
abund = trees[,c("RecordID","taxon.rank","tax_check","notes","accession","coletor","number","record","speciesLink","status.spLink","DetBy","DetDate"),]
trees1 <- left_join(trees1, abund, by = "RecordID")

#Filtering only the important info
nomes = c("refID","effort_ha","domain","forest_type","forest_subtype",
          "Unidade_de_conservacao","Unidade_de_conservacao1",
          "country","state","county","locality","locality1",
          "lat1","long1","confiabilidade","Reference","year","year_data",
          "ordem","SiteCode","status_diagnostico","SppListVerification",
          "dbh_cutoff1","dbh","dbh1")
surveys.final = surveys1[,nomes]

## Merging abundance and survey data
#Filtering the AF surveys with abundance data
surveys.final$ordem = as.character(surveys.final$ordem)
cwm = trees1
dados = cwm[cwm$ordem %in% surveys.final$ordem,]
dados$ordem = as.character(dados$ordem)

#Getting the multiple-site inventories
tmp = as.character(unique(cwm[!cwm$ordem %in% surveys.final$ordem,]$ordem))
tmp = tmp[grepl("\\|",tmp)]
tmp1 = strsplit(tmp,"\\|") 
tmp2 = sapply(tmp1,function(x) any(x %in% surveys.final$ordem))
tmp3 = tmp[tmp2]
dados1 = cwm[as.character(cwm$ordem) %in% tmp3,]

##Getting environmental-info for multiple-site inventories
tmp = as.character(dados1$ordem)
tmp1 = strsplit(tmp,"\\|") 
sites = surveys.final[surveys.final$ordem %in% unlist(tmp1),]
sites2 = sites

#create a final table with unique SiteCodes. 
#if factors: collapse levels
#if numerics: weighted mean by effort_ha
#if effort_ha: sum effort_ha
ordem = as.vector(dados1$ordem)

#columns indices that are factors (+ include ordem columns of which idx is 1)
factors.idx = c(1, which(sapply(1:ncol(sites2), function(x) is.character(sites2[,x]))==T)) 
factors.idx = factors.idx[!duplicated(factors.idx)]

#check if these columns are actually factors
colnames(sites2)[factors.idx]

sites3 = NULL
for(i in 1:length(ordem)) {
  tmp1 = strsplit(ordem[i],"\\|")[[1]] 
  tmp = sites2[which(sites2$ordem %in% tmp1),,drop=F]
  if(nrow(tmp)==1) {
    out = tmp
  } else {
    out = droplevels(sites2[1,,drop=F])
    #collapse factors
    out[,factors.idx] = sapply(factors.idx, 
                               function(y) 
                                 paste(unique(as.vector(tmp[,y])), collapse="|"))
    #average weighted by effort_ha
    out[,-factors.idx] = apply(tmp[,-factors.idx], 2, 
                               function(y) 
                                 weighted.mean(y, tmp$effort_ha, na.rm=T))
    #modify effort_ha
    out$effort_ha = sum(tmp$effort_ha)
  }
  if(is.null(sites3)) sites3 = out else sites3 = rbind.data.frame(sites3,out)	
}
surveys2 = sites3[!duplicated(sites3$ordem),]
rm(sites,sites2)

##MERGING CWM WITH ENVIRONMENTAL INFO
##Single site plots
dados2 = merge(dados,surveys.final, by="ordem",all.x=TRUE)

##Multiple-site plots
dados3 = merge(dados1,surveys2,by="ordem",all.x=TRUE)

##Combining all the data
trees.final = rbind.data.frame(dados2, dados3)
names(trees.final) = gsub("\\.x", "", names(trees.final))
trees.final = trees.final[order(trees.final$RecordID),]

##Editing collection numbers
tmp = trees.final[grepl('[0-9]', trees.final$accession),]
access = tmp$accession
access = as.character(unlist(sapply(strsplit(access,";"), head, n=1)))
access = as.character(unlist(sapply(strsplit(access,","), head, n=1)))
access = as.character(unlist(sapply(strsplit(access,"-"), head, n=1)))
access = as.character(unlist(sapply(strsplit(access,"_"), head, n=1)))
access = gsub('[A-Z][0-9]',' ',access)
access = gsub('  ',' ',access)
access = stringr::str_trim(access)
access = strsplit(access," ")
numTombo = paste(tolower(as.character(sapply(access,head,1))),
             as.character(sapply(access,tail,1)),sep="_")
trees.final$numTombo <- NA
trees.final$numTombo[grepl('[0-9]',trees.final$accession)] <- numTombo

## FILTERING AND SAVING DATA FOR SPECIFIC GROUPS
## Removing records no identified at species level
trees.final <- trees.final[trees.final$taxon.rank %in% c("species","variety","subspecies"),]

## Removing records not in the AF checklist
trees.final$species.correct2 <- sapply(strsplit(trees.final$species.correct," "), function(x) paste(x[1], x[2],sep=" "))
spp.af = readRDS("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/spp.af.rds")
spp.af.true <- spp.af[af.check2 %in% "TRUE", species.correct2]
#Adding some missing names (don't know why they were excluded; probably due to missing valid determinations or coordinates)
miss.spp = c("Agonandra brasiliensis","Aiouea bracteata","Aspidosperma brasiliense",
             "Aspidosperma quirandy","Cochlospermum vitifolium","Ficus pakkensis","Handroanthus selachidentatus",
             "Handroanthus spongiosus","Ilex cognata","Luetzelburgia purpurea","Mezilaurus synandra","Mimosa caesalpiniifolia",
             "Ocotea grandifructa","Persea punctata","Piptadenia killipii","Piranhea securinega","Quillaja lancifolia",
             "Tovomita longifolia","Trischidium decipiens")
trees.final <- trees.final[trees.final$species.correct2 %in% c(spp.af.true, miss.spp),]

## Removing records already in the herbarium data
tombos <- readRDS("data/threat_occ_data.rds")
tombos <- tombos$dup.ID1
tombos <- as.character(unlist(strsplit(tombos, "\\|")))
trees.final <- trees.final[!trees.final$numTombo %in% tombos,]

## Removing sites with poor coordinates ## 
trees.final <- trees.final[trees.final$confiabilidade %in% c("precisa","exata","boa","boa|precisa","precisa|exata","precisa|exata|boa"),]

## Removing missing coordinates ##
trees.final <- trees.final[!is.na(as.double(trees.final$lat1)) | !is.na(trees.final$long1),] 

## Removing Spatial Duplicates ##
inv.data <- data.table(trees.final)
setkeyv(inv.data, "species.correct2")
inv.data[ , coord.string := as.character(paste(lat1, long1, sep="_")),]
inv.data[ , spat.dups := duplicated(coord.string), by = species.correct2]
inv.data <- inv.data[spat.dups == FALSE,,]
inv.data[,coord.string:=NULL] 
inv.data[,spat.dups:=NULL] 
trees.final <- data.frame(inv.data)

## Creating the colum of taxonomic confidence
trees.final$tax_check2 = FALSE

## Creating the conservation unit information
trees.final$UC <- trees.final$Unidade_de_conservacao
trees.final$UC[trees.final$UC %in% "private_in_APA"] <- "APA"
trees.final$UC[trees.final$UC %in% "UC Protecao Integral"] <- "PI"
trees.final$UC[trees.final$UC %in% c("UC Uso Sustentavel", "US","universities and research centers")] <- "US"
trees.final$UC[trees.final$UC %in% "Terra Indigena"] <- "TI"
trees.final$UC[trees.final$UC %in% c("", NA)] <- "unknown"
trees.final$UC[trees.final$UC %in% c("Military","other_public_lands")] <- "US"
trees.final$UC[trees.final$UC %in% c("private|UC Protecao Integral")] <- "unknown"

## Removing occurrences closest to 1km from herbarium data
source("R/dist.valid.points.R")
oc.data <- readRDS("data/threat_occ_data.rds")
oc.data <- data.frame(ddlat = as.double(oc.data$latitude.work1),
                      ddlon = as.double(oc.data$longitude.work1),
                      tax = oc.data$species.correct2,
                      stringsAsFactors = FALSE)
oc.data <- oc.data[!is.na(oc.data$ddlat),]
oc.data <- oc.data[!is.na(oc.data$ddlon),]
inv.data <- data.frame(ddlat = trees.final$lat1,
                       ddlon = trees.final$long1,
                       tax = trees.final$species.correct2,
                       stringsAsFactors = FALSE)
inv.data <- .dist.valid.points(inv.data, oc.data, NbeCores = 5)
trees.final$dist.herb.km <- round(inv.data$distance/1000,2)
trees.final <- trees.final[trees.final$dist.herb.km >= 2,]
rm(oc.data, inv.data)

## Removing unnecessary columns
cols = c("lat1","long1","species.correct2","family",
  "coletor","number","year_data","numTombo","DetBy","DetDate",
  "tax_check2","tax_check","notes","confiabilidade","domain","UC","forest_type",
  "country","state","county","locality","locality1",
  "species.correct","N","accession","record","year","Reference","ordem")
trees.final = trees.final[,cols]

## Saving ##
saveRDS(trees.final, file = "data/threat_inventory_data.rds", compress = "gzip")

#### OBTAINING THE NUMBER OF OCCURRENCES PER SPECIES FROM TREECO ####
tmp <- table(trees.final$species.correct2)
tmp1 <- merge(resultado, 
             data.frame(species.correct2 = names(tmp), treeco.occs = as.double(tmp)), 
             by = "species.correct2", all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species.correct2),]
#table(resultado$species.correct2 == tmp1$species.correct2)
resultado$treeco.occs <- tmp1$treeco.occs
rm(trees,trees1,trees.final)


#####################################################################################################################################################################H
#####################################################################################################################################################################H
#########################H
#### POPULATION SIZES ####
#########################H

pop.sizes <- readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//pop.size.est_nmx50_mxd5_idp2.rds")
pop.sizes.2018 <- pop.sizes$`2018`$mean 
tmp <- apply(pop.sizes.2018, 2, sum)
tmp <- data.frame(species.correct2 = gsub("_"," ",names(tmp)), 
                  pop.size.2018 = as.double(tmp), 
                  stringsAsFactors = FALSE)
tmp$species.correct2 <- gsub("\\.", "-", tmp$species.correct2)

tmp1 <- merge(resultado, tmp, by = "species.correct2", 
              all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species.correct2),]
table(resultado$species.correct2 == tmp1$species.correct2)
resultado$pop.size.2018 <- tmp1$pop.size.2018

## How many species in the checklist with pop.size estimates?
100*table(!is.na(resultado$pop.size.2018))/dim(resultado)[1]

## Pop size an number of occurrences are related?
plot(log(resultado$total.occs+1) ~ log(resultado$pop.size.2018+1))
abline(MASS::rlm(log(resultado$total.occs+1) ~ log(resultado$pop.size.2018+1)), lwd=2, col=3)

## Which specie swe have population sizes which are not in the checklist?
miss.sp <- tmp$species.correct2[!tmp$species.correct2 %in% resultado$species.correct2]
length(miss.sp) # 215 species; most are non-AF species
# tmp <- flora::get.taxa(miss.sp, life.form = TRUE, vegetation.type = TRUE, establishment = TRUE)
# tmp1 <- flora::get_domains(tmp)
# tmp2 <- tmp1[grepl("ata Atl", tmp1$domain),]
# spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
# tmp3 <- merge(tmp2, spp, by.x = "original.search", by.y = "SpeciesKTSA")
# write.csv(tmp3, "tmp.csv")

## Saving the data vailable per source of information
taxon_id <- paste0("sp", 1:dim(resultado)[1])
resultado <- cbind.data.frame(taxon_id, resultado, stringsAsFactors = FALSE)
saveRDS(resultado, "data/assess_iucn_spp.rds")


#### PREPARING THE DATA FOR THE ASSESSMENTS USING ConR ###
mean.pop.sizes <- sapply(pop.sizes, function(x) apply(x$mean, 2, sum, na.rm = TRUE))
low.pop.sizes <- sapply(pop.sizes, function(x) apply(x$low, 2, sum, na.rm = TRUE))
high.pop.sizes <- sapply(pop.sizes, function(x) apply(x$high, 2, sum, na.rm = TRUE))
table(rownames(mean.pop.sizes) == rownames(low.pop.sizes))
table(rownames(mean.pop.sizes) == rownames(high.pop.sizes))
rm(pop.sizes,pop.sizes.2018)

## OBTAINING POPULATION SIZES FOR MISSING YEARS ##
## Getting pop. sizes for all possible generation lengths
gen.lengths <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100)
pos.gen.length <- sort(unique(gen.lengths * rep(c(1:3), each=length(gen.lengths))), decreasing = TRUE)
miss.years <- 2018 - pos.gen.length
miss.years1 <- miss.years[miss.years < 1992]
miss.years2 <- miss.years[miss.years >= 1992]

## Looping the pop. decline trends for each species and time interval
require(snow)
require(doSNOW)
require(foreach)
require(ConR)

dados <- list(mean.pop.sizes, low.pop.sizes, high.pop.sizes)
results <- vector("list", length(dados))
for(i in 1:length(dados)) {
  #defining the object for the loop
  pop_data <- dados[[i]]
  
  #Setting the loop parameters
  cl <- snow::makeSOCKcluster(6)
  doSNOW::registerDoSNOW(cl)
  `%d%` <- foreach::`%dopar%`
  x <- NULL
  output <-
    foreach::foreach(
      x = 1:dim(pop_data)[1],
      #.combine = 'c',
      .options.snow = NULL
    ) %d% {
  
      res1 <- ConR:::pop.decline.fit(pop.size = pop_data[x,],
                              years = c(1850,1940,1945,1950,1955,1960,1965,1970,1975,1980,1985,1992),
                              #models = "all", 
                              models = c("linear", "quadratic", "exponential", 
                                         "logistic", "general_logistic"),
                              project.years = miss.years1,
                              plot.fit = FALSE)
      res2 <- ConR:::pop.decline.fit(pop.size = pop_data[x,],
                              years = c(1992,1995,2000,2005,2010,2015,2018),
                              #models = "all", 
                              models = c("linear", "quadratic", "exponential", 
                                         "logistic", "general_logistic"),
                              project.years = miss.years2,
                              #parsimony = TRUE,
                              max.count = 20,
                              plot.fit = FALSE)
      res1$data$Modelo <- attributes(res1$best.model)$best.model.name
      res2$data$Modelo <- head(attributes(res2$best.model)$best.model.name, 1)
      dados <- rbind.data.frame(head(res1$data, -1), res2$data, stringsAsFactors = FALSE) 
      #dados$Predicted[dados$Year < 1850] <- 
        #dados$Predicted[dados$Year == 1850]
      dados <- dados[,c("Year", "Observed", "Predicted", "Modelo")]
      dados
    }
  snow::stopCluster(cl)
  
  #Saving the predictions for each table (mean, low and high) 
  results[[i]] <- output 
}
  
#Setting species names
rownames(mean.pop.sizes) <- gsub("_"," ", rownames(mean.pop.sizes))
rownames(mean.pop.sizes) <- gsub("\\.","-", rownames(mean.pop.sizes))
for(i in 1:length(results)) names(results[[i]]) <- rownames(mean.pop.sizes)

## Putting data in the ConR format
years <- results[[1]][[1]]$Year
ncols <- length(years)
spp <- names(results[[1]])
nrows <- length(spp)
mean.pop.conR <- low.pop.conR <- high.pop.conR <- matrix(NA, ncol = ncols, nrow = nrows,
                                                            dimnames = list(spp, years))
for(x in 1:length(results[[1]])) {
  mean.pop.conR[x,] <- results[[1]][[x]]$Predicted
  low.pop.conR[x,] <- results[[2]][[x]]$Predicted
  high.pop.conR[x,] <- results[[3]][[x]]$Predicted
}

#Converting into data.frames
mean.pop.conR <- cbind.data.frame(species = rownames(mean.pop.conR), mean.pop.conR, 
                                   row.names = NULL, stringsAsFactors = FALSE)
low.pop.conR <- cbind.data.frame(species = rownames(low.pop.conR), low.pop.conR, 
                                  row.names = NULL, stringsAsFactors = FALSE)
high.pop.conR <- cbind.data.frame(species = rownames(high.pop.conR), high.pop.conR, 
                                   row.names = NULL, stringsAsFactors = FALSE)

#### SAVING ###
#Saving the estimated populations (ONLY "OBSERVED" POP SIZES)
saveRDS(mean.pop.sizes, "data/threat_mean_pop_sizes.rds")
saveRDS(low.pop.sizes, "data/threat_low_pop_sizes.rds")
saveRDS(high.pop.sizes, "data/threat_high_pop_sizes.rds")

#Saving the estimated and infered populations (BOTH "OBSERVED" AND ESTIMATED/INTERPOLATED POP SIZES)
saveRDS(results[[1]], "data/threat_mean_pop_sizes_infer.rds")
saveRDS(results[[2]], "data/threat_low_pop_sizes_infer.rds")
saveRDS(results[[3]], "data/threat_high_pop_sizes_infer.rds")

#Saving the estimated and infered populations in the ConR format
saveRDS(mean.pop.conR, "data/threat_mean_pop_sizes_for_ConR.rds")
saveRDS(low.pop.conR, "data/threat_low_pop_sizes_for_ConR.rds")
saveRDS(high.pop.conR, "data/threat_high_pop_sizes_for_ConR.rds")

#####################################################################################################################################################################H
#####################################################################################################################################################################H
#####################################################
#### SPECIES TAXONOMY, SYNONYMS AND COMMON NAMES ####
#####################################################
rm(list=ls())
require(taxize)
require(flora)
require(redlistr)

## Red List IUCN key: 814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25

## Seting the ncbi key
#API-key: c9366b59140a833ccd26fe9ac4f92e48d208 
options(ENTREZ_KEY = "c9366b59140a833ccd26fe9ac4f92e48d208")

## Getting NCBI full taxonomy ##
tax  <- readRDS("data/assess_iucn_spp.rds")[,1:3]

#higher taxonomy
fams <- sort(unique(tax$family.correct1))
tmp  <- sapply(fams, taxize::classification, "ncbi")
tmp1 <- lapply(tmp, function(x) x$name[x$rank %in% c("kingdom","phylum","class","order")]) 
tmp1 <- do.call(rbind.data.frame, tmp1)
tmp2 <- cbind.data.frame(tmp1, fams, stringsAsFactors = FALSE)
names(tmp2) <- c("kingdom","phylum","classname","ordername","family")
ncbi.taxonomy <- merge(tax, tmp2, by.x = "family.correct1", by.y = "family", all.x = TRUE, order = FALSE)
ncbi.taxonomy$genus <- sapply(strsplit(ncbi.taxonomy$species.correct2, " "), function(x) x[1])
ncbi.taxonomy <- ncbi.taxonomy[,c("kingdom","phylum","classname","ordername","family.correct1","genus","species.correct2")]
ncbi.taxonomy$kingdom <- "Plantae" #IUCN standard notation
ncbi.taxonomy$phylum <- "Tracheophyta" #IUCN standard notation
ncbi.taxonomy <- ncbi.taxonomy[order(ncbi.taxonomy$species.correct2),]

#lower taxonomy
# gens <- sort(unique(ncbi.taxonomy$genus))
# # res.gens <- vector("list", length(gens))
# # for(i in 1:length(gens)){
# #   cat(i,"\n")
# #   res.gens[[i]] <- try(children(gens[i], db = "ncbi"), TRUE)
# # }
# tmp  <- lapply(gens, function(x) try(taxize::children(x, "ncbi"), TRUE))
# tmp1 <- do.call(rbind.data.frame, tmp1)

## Getting Tropicos full taxonomy ##
#Tropicos API key: E8E538BC-0DAD-47EF-84CA-F6A663D9170A
options(ENTREZ_KEY = "E8E538BC-0DAD-47EF-84CA-F6A663D9170A")

#higher taxonomy
tmp  <- sapply(fams, taxize::classification, "tropicos", rows = 1)
tmp1 <- lapply(tmp, function(x)
    as.character(x$name[x$rank %in% c("class", "order")]))
tmp1 <- do.call(rbind.data.frame, tmp1)
tmp2 <- cbind.data.frame(tmp1, fams, stringsAsFactors = FALSE)
names(tmp2) <- c("classname","ordername","family")
tps.taxonomy <- merge(tax, tmp2, by.x = "family.correct1", by.y = "family", all.x = TRUE, order = FALSE)
tps.taxonomy$genus <- sapply(strsplit(tps.taxonomy$species.correct2, " "), function(x) x[1])
tps.taxonomy$kingdom <- "Plantae" #IUCN standard notation
tps.taxonomy$phylum <- "Tracheophyta" #IUCN standard notation
tps.taxonomy <- tps.taxonomy[,c("kingdom","phylum","classname","ordername","family.correct1","genus","species.correct2")]

#lower taxonomy
cl <- snow::makeSOCKcluster(6)
doSNOW::registerDoSNOW(cl)
`%d%` <- foreach::`%dopar%`
x <- NULL
output <- foreach::foreach(
  x = 1:length(tps.taxonomy$species.correct2),
  #.combine = 'c',
  .options.snow = NULL
) %d% {
  try(taxize::get_tpsid_(tps.taxonomy$species.correct2[x],
                         key = "E8E538BC-0DAD-47EF-84CA-F6A663D9170A"),
      TRUE)
}
snow::stopCluster(cl)
output1 <- sapply(output, function(x) x[[1]])
names(output1) <- sapply(output, function(x) names(x))
output2 <- lapply(1:length(output1), function(x)
  if (is.null(output1[[x]])) {
    c(nameid= NA, scientificname = names(output1[x]), rankabbreviation = "sp.", 
      nomenclaturestatusname = "Not found", author=NA, displayreference= NA, displaydate = NA,
      search.string = names(output1[x]))
  } else {
    out <- output1[[x]][,c("nameid",
              "scientificname",
              "rankabbreviation",
              "nomenclaturestatusname",
              "author",
              "displayreference",
              "displaydate")]
    if (dim(out)[1] > 1) {
      out <- out[out$rankabbreviation %in% "sp.", ]
      if (dim(out)[1] > 1 & any(out$scientificname %in% names(output1[x])))
        out <- out[out$scientificname %in% names(output1[x]),]
      if (dim(out)[1] > 1 & !all(out$nomenclaturestatusname %in% "Invalid"))
        out <- out[!out$nomenclaturestatusname %in% "Invalid",]
      if (dim(out)[1] > 1 & !all(out$nomenclaturestatusname %in% "Illegitimate"))
        out <- out[!out$nomenclaturestatusname %in% "Illegitimate",]
      if (dim(out)[1] > 1 & "Legitimate" %in% out$nomenclaturestatusname)
        out <- out[out$nomenclaturestatusname %in% "Legitimate",]
      if (dim(out)[1] > 1 & !all(out$displayreference %in% ""))
        out <- out[!out$displayreference %in% "",]
      if (dim(out)[1] > 1)
        out <- out[which.min(out$displaydate),]
    }
    cbind.data.frame(out, search.string = names(output1[x]), 
                     stringsAsFactors = FALSE)
  }
)
output3 <- do.call(rbind.data.frame, output2)
tps.taxonomy1 <- merge(tps.taxonomy, output3, by.x = "species.correct2", by.y = "search.string", all.x = TRUE, sort= FALSE)


## Getting FB 2020 taxonomy
toto <- flora::get.taxa(tax$species.correct2, replace.synonyms = FALSE, suggest.names = FALSE,
         life.form = TRUE, habitat = TRUE, vegetation.type = TRUE,
         vernacular = TRUE, states = TRUE, establishment = TRUE,
         drop = c(""),
         suggestion.distance = 0.9, parse = FALSE)
table(toto$accepted.name)
toto <- toto[c("family","genus","original.search","authorship","id","accepted.name","taxon.status","name.status","notes",
               "threat.status","life.form","habitat","vegetation.type","vernacular.name","occurrence","establishment")]

## Trying to get species ID from IUCN
## Red List IUCN key: 814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25
cl <- snow::makeSOCKcluster(6)
doSNOW::registerDoSNOW(cl)
`%d%` <- foreach::`%dopar%`
x <- NULL
iucn <- foreach::foreach(
  x = 1:length(tps.taxonomy$species.correct2),
  #.combine = 'c',
  .options.snow = NULL
) %d% {
  try(taxize::get_iucn(tax$species.correct2[x], 
                   key = "814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25"), TRUE)
}
snow::stopCluster(cl)

iucn1 <- lapply(1:length(iucn), function(x)
  if (is.null(iucn[[x]])) {
    c(iucnid= NA, match = "error", name = tax$species.correct2[x], uri = NA)
  } else {
    atributos <- unlist(attributes(iucn[[x]]))
    out <- c(iucnid= iucn[[x]][1], atributos[1:2]) 
  }
)
iucn2 <- do.call(rbind.data.frame, iucn1)
names(iucn2) <- c("Redlist_id","match","name")
for(i in 1:3) iucn2[,i] <- as.character(iucn2[,i])
toto1 <- merge(toto, iucn2[,!names(iucn2) %in% "match"], 
                  by.x = "original.search", by.y = "name", all.x = TRUE, sort= FALSE)

#### PUTTING ALL TAXONOMIC INFORMATION TOGETHER AND SAVING ####
full.tax <- merge(toto1, tps.taxonomy1[,!names(tps.taxonomy1) %in% "genus"], 
                            by.x = "original.search", by.y = "species.correct2", all.x = TRUE, sort= FALSE)
full.tax$source <- "the Brazilian Flora 2020 (http://floradobrasil.jbrj.gov.br)"
full.tax$source[is.na(full.tax$authorship)] <- "Tropicos (www.tropicos.org)" 
full.tax$id[is.na(full.tax$authorship)] <- 
  full.tax$nameid[is.na(full.tax$authorship)]
full.tax$authorship[is.na(full.tax$authorship)] <- 
  full.tax$author[is.na(full.tax$authorship)]
full.tax[full.tax$original.search %in% "Ocotea koscinskii",c("family","genus","authorship")] <-
  c("Lauraceae","Ocote","Baitello & Brotto")
full.tax$source[is.na(full.tax$authorship)] <- "GBIF (www.gbif.org)" 

#Final search using other databases (GBIF works better)
splist <- full.tax$original.search[is.na(full.tax$authorship)]
miss.tax <- get_gbifid_(splist)
miss.tax1 <- lapply(1:length(miss.tax), function(x) {
                  out <- miss.tax[[x]][,c("usagekey","kingdom","phylum","order","class","family","genus","canonicalname",
                                          "rank","status","matchtype","synonym","scientificname")]
                  if (dim(out)[1] > 1 & any(out$canonicalname %in% splist[x])) 
                    out <- out[out$canonicalname %in% splist[x], ]
                  if (dim(out)[1] > 1 & !all(out$status %in% "SYNONYM"))
                    out <- out[out$status %in% "ACCEPTED",]
                  out
})
miss.tax1 <- do.call(rbind.data.frame,miss.tax1)
miss.tax1$authorship <- stringr::str_trim(sapply(strsplit(miss.tax1$scientificname," "), function(x) paste(tail(x,-2), collapse = " ")))
cols <- c("family","genus","authorship","usagekey","status","synonym")
full.tax[full.tax$original.search %in% splist, 
         c("family","genus","authorship","id","taxon.status","notes")] <- miss.tax1[,cols]
full.tax$family[is.na(full.tax$family)] <- 
  full.tax$family.correct1[is.na(full.tax$family)]
full.tax$genus[is.na(full.tax$genus)] <- 
  sapply(strsplit(full.tax$original.search[is.na(full.tax$genus)], " "), function(x) x[1])

## Filtering and creating the final columns
full.tax1 <- full.tax[ ,c("kingdom","phylum","classname","ordername","family","genus","original.search","authorship","Redlist_id")]
                         #,"accepted.name","taxon.status","name.status","notes","nomenclaturestatusname","displayreference","displaydate"
                        #"life.form","habitat","vegetation.type","vernacular.name","occurrence"
#Taxonomic Reference
full.tax1$TaxonomicReference <-
  paste0("More taxonomic information on this species can be found at: http://servicos.jbrj.gov.br/flora/search/",
    gsub(" ", "_", full.tax$original.search),".")
full.tax1$TaxonomicReference[full.tax$source %in% "Tropicos (www.tropicos.org)"] <-
  paste0("More taxonomic information on this species can be found at: https://www.tropicos.org/name/", full.tax$id[full.tax$source %in% "Tropicos (www.tropicos.org)"],".")
full.tax1$TaxonomicReference[full.tax$source %in% "GBIF (www.gbif.org)"] <-
  paste0("More taxonomic information on this species can be found at: https://www.gbif.org/species/", full.tax$id[full.tax$source %in% "GBIF (www.gbif.org)"],".")
ids <- !(is.na(full.tax$displayreference) | full.tax$displayreference %in% c(""," ")) & !(is.na(full.tax$displaydate) | full.tax$displaydate %in% c(""," ")) 
full.tax1$TaxonomicReference[ids] <- 
  paste0(full.tax1$TaxonomicReference[ids], paste0(" According to Tropicos (www.tropicos.org), this species was published in: ",
                                                  paste0(full.tax$displayreference[ids]," (", full.tax$displaydate[ids], ").")))
ids <- !(is.na(full.tax$displayreference) | full.tax$displayreference %in% c(""," ")) & (is.na(full.tax$displaydate) | full.tax$displaydate %in% c(""," "))
full.tax1$TaxonomicReference[ids] <- 
  paste0(full.tax1$TaxonomicReference[ids], paste0(" According to Tropicos (www.tropicos.org), this species was published in: ",
                                                   full.tax$displayreference[ids]))

#Taxonomic Notes
full.tax1$TaxonomicNotes.value <- NA
ids <- full.tax$taxon.status %in% c("accepted","ACCEPTED") & (grepl("correct", full.tax$name.status) | is.na(full.tax$name.status))
full.tax1$TaxonomicNotes.value[ids] <- paste0("This taxon is accepted and is taken as being correct by ", full.tax$source[ids], ".")
ids <- is.na(full.tax1$TaxonomicNotes.value) & !is.na(full.tax$name.status)
full.tax1$TaxonomicNotes.value[ids] <- paste0("This taxon is taken as ",full.tax$name.status[ids], "by ", full.tax$source[ids], ".")
ids <- is.na(full.tax1$TaxonomicNotes.value) & full.tax$nomenclaturestatusname %in% c("No opinion","Legitimate") & full.tax$source %in% "Tropicos (www.tropicos.org)"
full.tax1$TaxonomicNotes.value[ids] <- "This taxon is accepted in Tropicos (www.tropicos.org)."
ids <- is.na(full.tax1$TaxonomicNotes.value)
full.tax1[ids,]
#### CHECK HERE?: ADD CONFLICTS BETWEEN DATABASES (ReFlora vs. GBIF) ####
#### CHECK HERE: OTHER ACCEPTED NAMES FROM REFLORA: e.g. Schinus terebinthifolius Raddi

#### Three new combintions for Calyptranthes and Marlierea:
#"Calyptranthes brasiliensis" -> ""
#"Marlierea eugenioides" -> "Myrcia eugenioides Cambess."
#"Marlierea laevigata" -> "Myrcia multipunctata Mazine"


## Final edits and saving
names(full.tax1)[names(full.tax1) %in% "original.search"] <- "species"
names(full.tax1)[names(full.tax1) %in% "authorship"] <- "taxonomicAuthority"
full.tax1 <- full.tax1[order(full.tax1$species),]
table(full.tax1$species == ncbi.taxonomy$species.correct2)
full.tax1$classname <- ncbi.taxonomy$classname

#Solving encoding problems
full.tax1$taxonomicAuthority <- gsub("Ã¼","ü",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã£","ã",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã©","é",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã¡","á",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã¨","è",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã³","ó",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã±","ñ",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã¶","ö",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã§","ç",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ãº","ú",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã¥","a",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã´","ô",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority <- gsub("Ã\u0081vila","Ávila",full.tax1$taxonomicAuthority)
full.tax1$taxonomicAuthority[grepl("Ã",full.tax1$taxonomicAuthority)] <- 
  gsub("A-","í",textclean::replace_non_ascii(full.tax1$taxonomicAuthority[grepl("Ã",full.tax1$taxonomicAuthority)]))

##Saving
taxon_id <- paste0("sp", 1:dim(full.tax1)[1])
full.tax2 <- cbind.data.frame(taxon_id, full.tax1, stringsAsFactors = FALSE)
write.csv(full.tax2, "data/threat_taxonomy.csv", row.names = FALSE, fileEncoding = "UTF-8")


#####################H
#### COMMON NAMES ####
#####################H

#Filtering and editing
taxon_id <- paste0("sp", 1:dim(full.tax1)[1])
common <- cbind.data.frame(taxon_id, full.tax[ ,c("original.search","vernacular.name")],
                           stringsAsFactors = FALSE)
table(common$taxon_id == full.tax2$taxon_id)
names(common) <- c("taxon_id","validName","name0")
common1 <- common[!is.na(common$name0),]

#Solving encoding problems
common1$name0 <- gsub("Ã¼", "ü", common1$name0)
common1$name0 <- gsub("Ã£", "ã", common1$name0)
common1$name0 <- gsub("Ã©", "é", common1$name0)
common1$name0 <- gsub("Ã¡", "á", common1$name0)
common1$name0 <- gsub("Ã¨", "è", common1$name0)
common1$name0 <- gsub("Ã³", "ó", common1$name0)
common1$name0 <- gsub("Ã±", "ñ", common1$name0)
common1$name0 <- gsub("Ã¶", "ö", common1$name0)
common1$name0 <- gsub("Ã§", "ç", common1$name0)
common1$name0 <- gsub("Ãº", "ú", common1$name0)
common1$name0 <- gsub("Ã¥", "a", common1$name0)
common1$name0 <- gsub("Ã´", "ô", common1$name0)
common1$name0 <- gsub("Ãª", "ê", common1$name0)
common1$name0 <- gsub("Ãµ", "õ", common1$name0)
common1$name0[grepl("Ã", common1$name0)] <- 
  gsub("A-", "í", textclean::replace_non_ascii(common1$name0[grepl("Ã", common1$name0)]))

#minor fixes
common1$name0 <- gsub("pao \\/ pau \\/ pao", "pau", common1$name0)
common1$name0 <- gsub("pao \\/ pau", "pau", common1$name0)
common1$name0 <- gsub("cumixá/cumichá", "cumixá", common1$name0)
common1$name0 <- gsub("Caixeta/Caixeto", "Caixeta", common1$name0)

#Breaking the info
tutu <- strsplit(common1$name0, "\\|")
tutu <- lapply(tutu, stringr::str_trim)
tutu <- lapply(tutu, function(x) t(sapply(strsplit(x, '\\/'), function(y) y[1:2])))
names(tutu) <- common1$validName
tutu1 <- do.call('rbind.data.frame', tutu)
tutu1 <- cbind.data.frame(validName = sapply(strsplit(rownames(tutu1), "\\."), function(x) x[1]),
                          language = as.character(tutu1[,2]),
                          name = stringr::str_to_title(as.character(tutu1[,1])),
                          stringsAsFactors = FALSE)
tutu1$language <- stringr::str_replace_all(tutu1$language,
                                           c("PORTUGUES" = "Portuguese","ESPANHOL" = "Spanish","INGLES" = "English", "KAXINAWA" = "South American Indian (Other)"))
## Saving
tutu2 <- merge(tutu1, common1[,1:2], by = "validName", all.x = TRUE)
tutu2 <- tutu2[order(tutu2$validName),]
tutu2 <- tutu2[,c("taxon_id","language","name")]
write.csv(tutu2, "data/threat_commonnames.csv", row.names = FALSE, fileEncoding = "UTF-8")

#################H
#### SYNONYMS ####
#################H

flora.syn <- lapply(tax$species.correct2, flora::get.synonyms)
flora.syn.authot <- lapply(1:length(flora.syn), function(x)
  if (is.null(flora.syn[[x]])) {
  } else {
    validName <- tax$species.correct2[x]
    out <- flora::get.taxa(flora.syn[[x]],replace.synonyms = FALSE, suggest.names = FALSE,
                           life.form = FALSE, habitat = FALSE, vegetation.type = FALSE,
                           vernacular = FALSE, states = FALSE, establishment = FALSE,
                           drop = c(""))
    name <- flora.syn[[x]]
    speciesName <- flora.syn[[x]]
    speciesAuthor <- out$authorship
    cbind.data.frame(validName, name, speciesName, speciesAuthor, 
                     stringsAsFactors = FALSE)
  }
)
synonyms <- do.call(rbind.data.frame, flora.syn.authot)

##Final Edits 
#Solving encoding problems
synonyms$speciesAuthor <- gsub("Ã¼","ü",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã£","ã",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã©","é",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã¡","á",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã¨","è",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã³","ó",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã±","ñ",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã¶","ö",synonyms$speciesAuthor)
synonyms$speciesAuthor <- gsub("Ã§","ç",synonyms$speciesAuthor)

#getting and removing infra-specific synonyms
synonyms$infraType <- NA
synonyms$infraType[grepl(" var\\. ",synonyms$speciesName)] <- "var."
synonyms$speciesName <- gsub(" var\\. "," ",synonyms$speciesName)
synonyms$infraType[grepl(" subsp\\. ",synonyms$speciesName)] <- "subsp."
synonyms$speciesName <- gsub(" subsp\\. "," ",synonyms$speciesName)
synonyms$infraType[grepl(" fo\\. ",synonyms$speciesName)] <- "fo."
synonyms$speciesName <- gsub(" fo\\. "," ",synonyms$speciesName)

#isolation the infra-specific author, name and cleaning speciesName
synonyms$infrarankAuthor <- NA
synonyms$infrarankAuthor[!is.na(synonyms$infraType)] <-
  synonyms$speciesAuthor[!is.na(synonyms$infraType)]
synonyms$speciesAuthor[!is.na(synonyms$infraType)] <- NA

synonyms$infrarankName <- NA
synonyms$infrarankName[!is.na(synonyms$infraType)] <-
  sapply(strsplit(synonyms$speciesName[!is.na(synonyms$infraType)], " "), function(x) x[3])
synonyms$speciesName <- 
  stringr::str_trim(sapply(strsplit(synonyms$speciesName, " "), function(x) paste(x[1],x[2], sep=" ")))

## Saving
synonyms1 <- merge(synonyms, tax[,c("taxon_id","species.correct2")],
                  by.x = "validName", by.y = "species.correct2", all.x = TRUE)
synonyms1 <- synonyms1[order(synonyms1$validName),]
synonyms1 <- synonyms1[,c("taxon_id","name","speciesName","infrarankName","infraType","speciesAuthor","infrarankAuthor")]
write.csv(synonyms, "data/threat_synonyms.csv", row.names = FALSE, fileEncoding = "UTF-8")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
##############################
#### HABITATS AND ECOLOGY ####
##############################
rm(list = ls())

## Getting Flora do BRasil information ##
tax  <- readRDS("data/assess_iucn_spp.rds")[,1:3]
toto <- flora::get.taxa(tax$species.correct2, replace.synonyms = FALSE, suggest.names = FALSE,
                        life.form = TRUE, habitat = TRUE, vegetation.type = TRUE,
                        vernacular = FALSE, states = FALSE, establishment = TRUE,
                        drop = c(""),
                        suggestion.distance = 0.9, parse = FALSE)
hab <- cbind.data.frame(tax, 
                         toto[,c("life.form","habitat","vegetation.type","establishment")],
                         stringsAsFactors = FALSE)

#### OBTAINING GOOD GUESSES OF GENERATION LENGTH BASED ON SPECIES ECOLOGICAL GROUPS AND MAXIMUM SIZE ####
path = "C:/Users/renato/Documents/raflima/Pos Doc/Databases/TreeCo Database Management"
source(paste(path,"trait_data_prep.R",sep="/"))

## List of species
hab$genus <- sapply(strsplit(hab$species.correct2," "), function(x) x[1])
hab$taxon.rank <- "species"
hab <- hab[,c("taxon_id","species.correct2","family.correct1","genus","species.correct2","taxon.rank","life.form","habitat","vegetation.type")]
names(hab) <- c("taxon_id", "Name_submitted","family","genus","species.correct","taxon.rank","life.form.reflora","habitat.reflora","vegetation.type.reflora")
hab$ordem.spp = 1:dim(hab)[1]

## South American tree species list
trees.sa = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv",as.is=TRUE,na.string=c(""," ",NA))
## Reading trait data per taxonomic level ####
trait.spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species traits//Atributos das espécies//traits.species.csv",as.is=T)
trait.gen = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species traits//Atributos das espécies//traits.genus.csv",as.is=T)
trait.fam = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species traits//Atributos das espécies//traits.family.csv",as.is=T)
traits1 = trait_data_prep(hab, trees.sa, trait.spp, trait.gen, trait.fam)$data_frame
table(hab$species.correct == tax$species.correct2)
traits1$wsg_gcm3 <- as.double(traits1$wsg_gcm3)
traits1$SeedMass_g <- as.double(traits1$SeedMass_g)


## Getting the final table
aux <- merge(hab[,c("family","species.correct")], trees.sa[,c(1:20,74:75,138:140,144:146,170:177)], 
             by.x = "species.correct", by.y = "SpeciesKTSA", all.x = TRUE)
aux <- aux[order(aux$species.correct),]

# Establishment in respect to the AF
hab$establishment <- traits1$establishment.AF
hab$establishment[is.na(hab$establishment)] <- traits1$establishment.BR[is.na(hab$establishment)]
#hab$establishment[is.na(hab$establishment)] <- "check"
hab$establishment[hab$species.correct %in% c("Ormosia monosperma")] <- "native"
hab$establishment[hab$establishment %in% "native" & grepl("not in B",aux$establishment.BR) & !is.na(aux$Distribution.VCS)] <- "AF native but not in Brazil"

## Flagging species according to their habit
hab$habito <- traits1$habito
hab$habito[hab$habito %in% c("1","1?")] = "tree"  
hab$habito[hab$habito %in% c("0.5","0.5?")] = "shrub"  

## Flagging species according to their growth form
#1 Tree - size unknown Tree (any size); 2 Tree - large Large tree (>15 m); 3 Tree - small Small tree (1-15 m)
#4 Shrub - size unknown; 5 Shrub (>1m); 6 Shrub - small (<1m)
#16 Succulent - form unknown; 18 Succulent - shrub (generally <1 m); 19 Succulent - tree (generally >1 m); #20 Fern
hab$life.form <- traits1$life.form
hab$life.form[hab$life.form %in% "palm" & !is.na(aux$Growth.habits)] <- 
  aux$Growth.habits[hab$life.form %in% "palm" & !is.na(aux$Growth.habits)]
hab$life.form[hab$life.form %in% "succulent_tree" & !is.na(aux$Growth.habits) & grepl("Succulent", aux$Growth.habits)] <- 
  aux$Growth.habits[hab$life.form %in% "succulent_tree" & !is.na(aux$Growth.habits) & grepl("Succulent", aux$Growth.habits)]
hab$life.form[hab$life.form %in% "succulent_tree"] <- "Succulent tree"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & grepl("^Tree$", aux$Growth.habits)] <- "Woody tree (Tree)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & grepl("^Shrub or treelet$", aux$Growth.habits)] <- "Woody tree (Shrub or treelet)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & grepl("^Tree or shrub$", aux$Growth.habits)] <- "Woody tree (Tree or shrub)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & grepl("^Hemiepiphytic", aux$Growth.habits)] <- "Woody tree (Hemiepiphytic)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & grepl("^Lianescent", aux$Growth.habits)] <- "Woody tree (Lianescent)"
hab$life.form[hab$life.form %in% "woody_vines_and_subshrubs"] <- "Woody tree (Lianescent)"
hab$life.form[hab$life.form %in% "palm"] <- "Palm or palmoid"
hab$life.form[hab$life.form %in% "palmoids"] <- "Palm or palmoid"
hab$life.form[hab$life.form %in% "tree_fern"] <- "Tree fern"
hab$life.form[hab$life.form %in% "woody_bamboo"] <- "Woody bamboo"
hab$life.form[hab$life.form %in% "woody_tree" & hab$habito %in% "shrub"] <- "Woody tree (Shrub or treelet)"
hab$life.form[hab$life.form %in% "woody_tree" & hab$habito %in% "tree"] <- "Woody tree (Tree)"
hab$life.form[hab$life.form %in% "woody_tree" & hab$life.form.reflora %in% "Árvore"] <- "Woody tree (Tree)"
hab[hab$species.correct %in% "Mimosa setosa",c("habito","life.form")] <- c("shrub", "Woody tree (Shrub or treelet)")
table(hab$life.form)

## Flagging species according to their mean maximum size
hab$MaxHeight <- traits1$MaxHeight_m
hab$MaxHeight[is.na(hab$MaxHeight)] <- as.double(aux$Potential.height[is.na(hab$MaxHeight)])
#Some missing values retrieved manually
hab$MaxHeight[hab$Name_submitted %in% "Ephedranthus dimerus"] = 20
hab$MaxHeight[hab$Name_submitted %in% "Myrcia costeira"] = 9
hab$MaxHeight[hab$Name_submitted %in% "Psychotria pubigera"] = 4
#some final generalizations
small.shrubs = c("Athenaea","Capsicum","Piper","Eumachia","Chiococca","Stachytarpheta","Schweiggeria")
hab$MaxHeight[grepl(paste(small.shrubs,collapse="|"),hab$Name_submitted) & is.na(hab$MaxHeight)] = 3
taller.shrubs = c("Erythroxylum","Chamaecrista","Faramea","Palicourea","Psychotria","Conchocarpus","Strychnos","Cybianthus","Chomelia")
hab$MaxHeight[grepl(paste(taller.shrubs,collapse="|"),hab$Name_submitted) & is.na(hab$MaxHeight)] = 4
table(is.na(hab$MaxHeight))
# hab[hab$habito %in% "shrub" & !is.na(hab$MaxHeight) & hab$MaxHeight >15,]
# hab[hab$habito %in% "tree" & !is.na(hab$MaxHeight) & hab$MaxHeight <5,]


## Classifying species into the growth from classifications
hab$GF <- findInterval(hab$MaxHeight, c(0,5,15))
hab$GF[hab$GF %in% 1 & !is.na(aux$Potential.height) & aux$Potential.height >8] <- 2
hab$GF[hab$GF %in% 2 & hab$life.form %in% "Woody tree (Shrub or treelet)" & !is.na(aux$Potential.height) & aux$Potential.height <5] <- 1
hab$GF[hab$GF %in% 2 & hab$habito %in% "shrub"  & !is.na(aux$Potential.height) & aux$Potential.height <5] <- 1
hab$GF[hab$GF %in% 3 & hab$life.form %in% "Woody tree (Shrub or treelet)" & !is.na(aux$Potential.height) & aux$Potential.height <5] <- 2
hab$GF[hab$GF %in% 3 & hab$habito %in% "shrub"] <- 2
hab$GF[is.na(hab$GF) & hab$habito %in% "shrub" & hab$life.form  %in% "Woody tree (Shrub or treelet)" & 
         hab$life.form.reflora %in% c("Arbusto","Arbusto|Subarbusto","Arbusto|Erva|Subarbusto","Erva|Subarbusto","Subarbusto")] <- 1
hab$GF[is.na(hab$GF) & hab$life.form.reflora %in% c("Arbusto","Subarbusto") & hab$life.form  %in% "Woody tree (Shrub or treelet)"] <- 1
hab$GF[is.na(hab$GF) & hab$life.form  %in% "Woody tree (Tree)" & 
         hab$life.form.reflora %in% c("Arbusto", "Arbusto|Árvore|Liana/volúvel/trepadeira")] <- 2
hab$GF[is.na(hab$GF) & hab$life.form  %in% "Woody tree (Shrub or treelet)" & 
         grepl("Árvore", hab$life.form.reflora)] <- 2
hab$GF[is.na(hab$GF) & hab$life.form  %in% "Woody tree (Shrub or treelet)" & 
         grepl("Arbusto", hab$life.form.reflora)] <- 1
hab$GF[is.na(hab$GF) & (hab$life.form  %in% "Woody tree (Tree)" | 
         grepl("Árvore", hab$life.form.reflora))] <- 4
hab$GF <- stringr::str_replace_all(hab$GF, c("1" = "large_shrub", "2" = "small_tree", "3" = "large_tree", "4" = "tree_unknown"))
table(hab$GF, useNA = "always")

## Wood density
hab$wsg <- traits1$wsg_gcm3
#Some generalizations for missing data
hab$wsg[hab$family %in% "Berberidaceae"] = 0.593 # value from DRYAD for Berberis jobii
hab$wsg[hab$family %in% "Gentianaceae"] = 0.5 # value from Brown, S. for Anthocleista grandiflora
hab$wsg[hab$family %in% "Onagraceae"] = 0.59 # value from IFMG for Ludwigia elegans
hab$wsg[hab$family %in% "Schoepfiaceae"] = 0.757 # value from Magnago for Schoepfia brasiliensis
hab$wsg[hab$family %in% "Trigoniaceae"] = 0.6967 # value from Magnago for Trigoniodendron_spiritusanctense
hab$wsg[hab$genus %in% "Acanthocladus"] = 0.70 # value from Magnago for Acanthocladus pulcherrimus
hab$wsg[hab$family %in% "Polygalaceae"] = 0.66 # family average for different references
hab[is.na(hab$wsg),]

## Flagging species according to their ecological group
hab$ecol.group <- traits1$ecological.group
hab$ecol.group[hab$ecol.group %in% 1.5] <- 1
hab$ecol.group[hab$ecol.group %in% 2.5] <- 3
hab$ecol.group[hab$ecol.group %in% 3.5] <- 4
hab$ecol.group[is.na(hab$ecol.group) & !hab$family.correct1 %in% c("Arecaceae","Cyatheaceae","Poaceae","Cactaceae","Anacardiaceae","Sapindaceae","Lauraceae","Moraceae","Myrtaceae","Symplocaceae") & 
                  !is.na(hab$wsg) & hab$wsg < 0.4 & !is.na(traits1$SeedMass_g) & traits1$SeedMass_g < 0.1] <- 1
hab$ecol.group[is.na(hab$ecol.group) & !hab$family.correct1 %in% c("Arecaceae","Cyatheaceae","Poaceae","Cactaceae") & 
                  !is.na(hab$wsg) & hab$wsg < 0.55 & !is.na(traits1$SeedMass_g) & traits1$SeedMass_g < 0.25] <- 2
hab$ecol.group[is.na(hab$ecol.group) & !hab$family.correct1 %in% c("Arecaceae","Cyatheaceae","Poaceae","Cactaceae") & 
                  !is.na(hab$wsg) & hab$wsg > 0.85 & !is.na(traits1$SeedMass_g) & traits1$SeedMass_g > 1] <- 4
hab$ecol.group[is.na(hab$ecol.group) & !hab$family.correct1 %in% c("Arecaceae","Cyatheaceae","Poaceae","Cactaceae") & 
                  !is.na(hab$wsg) & hab$wsg > 0.7 & !is.na(traits1$SeedMass_g) & traits1$SeedMass_g > 0.75] <- 3
hab$ecol.group <- stringr::str_replace_all(hab$ecol.group, c("1" = "pioneer", "2" = "early_secondary", "3" = "late_secondary", "4" = "climax"))
hab$ecol.group[is.na(hab$ecol.group)] <- "unknown"
table(hab$ecol.group)

## Including other traits that may be important
table(hab$Name_submitted == traits1$Name_submitted)
hab$SeedMass_g <- traits1$SeedMass_g
hab$dispersal.syndrome <- traits1$dispersal.syndrome

##################################################
#### GENERATION LENGHT AND MATURE INDIVIDUALS ####
##################################################

#### SETTING THE PROXIES OF GENERATION LENGTH FOR EACH SPECIES ####
## Defining the generation lengths per combination of GF and ecol.group
g1 <- sort(paste(c(1,3,2,4), sort(unique(hab$GF)),sep="."))
g2 <- sort(paste(c(4,2,3,1,5), sort(unique(hab$ecol.group)),sep = "."))
combo <- expand.grid(EG = g2, GF = g1)
combo[,1] <- gsub('[0-9:]\\.',"", combo[,1])
combo[,2] <- gsub('[0-9:]\\.',"", combo[,2])

#Generation lengths
combo$GL <- c(10, 20, 25, 35, 25, # for shrubs
              20, 40, 50, 65, 45, # for small trees
              30, 50, 65, 80, 60, # for large trees
              25, 45, 50, 60, 50) # for trees unknown
#Diameter at onset of maturity
combo$DBH <- c(5, 5, 5, 5, 5, # for shrubs
               6, 7, 8, 9, 7.5, # for small trees
               10, 12.5, 15, 20, 12.5, # for large trees
               7.5, 10, 12.5, 15, 10) # for trees unknown
#Prop. mature
combo$p.est <- c(1, 1, 1, 1, 1, # for shrubs
                 0.7213, 0.5998, 0.4927, 0.2868, 0.5122, # for small trees
                 0.5842, 0.4498, 0.3047, 0.2533, 0.3342, # for large trees
                 0.6371, 0.4529, 0.3476, 0.2786, 0.4529) # for trees unknown
combo$p.ci <- c("0.99-1.1", "0.975-1.075", "0.95-1.05", "0.90-1.025", "0.95-1.05", # for shrubs
                 "0.5838-0.8587", "0.5212-0.6783", "0.4304-0.555", "0.1829-0.3906", "0.4737-0.5507", # for small trees
                 "0.5141-0.6544", "0.417-0.4825", "0.2756-0.3338", "0.1827-0.3238", "0.3161-0.3524", # for large trees
                 "0.6206-0.6536", "0.4352-0.4707", "0.3306-0.3646", "0.2624-0.2949", "0.4352-0.4707") # for trees unknown
#Merging the info
combo$str.match <- paste(combo$EG, combo$GF, sep = "_") 
hab$str.match <- paste(hab$ecol.group, hab$GF, sep = "_")
hab1 <- dplyr::left_join(hab, combo[,c("str.match","GL","DBH","p.est","p.ci")], by= "str.match") 

#Adding Proportion of mature individuals in the population (from code '')
p.est <- read.csv("data/prop_mature.csv", as.is = TRUE)
hab1 <- dplyr::left_join(hab1, p.est[,c("species.correct","N","dbh.crit","p")], by= "species.correct") 


##################################################H
#### PUTTING ALL FIELDS IN THE IUCN SIS FORMAT ####
##################################################H

#HabitatDocumentation.narrative: Habitat Information.
hab1$HabitatDocumentation.narrative <- NA

#GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup:	(Coding Option, e.g. 3.4)
## 1. Forest & Woodland: 
# 1.5 Subtropical/Tropical Dry Forest: includes Deciduous forests and Chaco
# 1.6 Subtropical/Tropical Moist Lowland Forest: below 1.200 m
# 1.8 Subtropical/Tropical Swamp Forest: typically flooded for at least part of the year and dependent on this flooding for its existence
# 1.9 Subtropical/Tropical Moist Montane Forest: above 1.200 m (include cloud forests)
## 2. Savanna
# 2.1 Dry Savanna: includes cerrado/campos/caatinga (Brazil, Guiana), chaco seco (Argentina/Uruguay)
# 2.2 Moist Savanna: includes pantanal (Brazil/Bolivia/Paraguay), chaco húmedo (Paraguay/Bolivia/Argentina), and llanos (Venezuela/Colombia).
## 3. Shrubland
# 3.5 Subtropical/Tropical Dry Shrubland: includes restingas (coastal scrub in Brazil); xerophyllous scrub; chapparal (Mexico and SW US); matorral sammofilo, matorral seco (Argentina, Paraguay, Uruguay).
# 3.7 Subtropical/Tropical High Altitude Shrubland: included Campo Rupestre ares by in Jung et al. (2020)
## 4. Grassland
# 4.7 Subtropical/Tropical High Altitude Grassland: included here as Campos Rupestres, but Jung et al. 2020 also classify them as 3.7
hab1$vt <- hab1$vegetation.type.reflora
forests <- 'Floresta Ombrófila \\(\\= Floresta Pluvial\\)|Floresta Estacional Semidecidual|Floresta Ombrófila Mista|Floresta de Terra Firme|Floresta Ciliar ou Galeria|Floresta Estacional Perenifólia'
hab1$vt <- gsub(forests, "1.6", hab1$vt, perl = TRUE)
savanas <- 'Cerrado \\(lato sensu\\)|Caatinga \\(stricto sensu\\)|Campo Limpo|Savana Amazônica'
hab1$vt <- gsub(savanas, "2.1", hab1$vt, perl = TRUE)
hab1$vt <- gsub("Carrasco|Restinga|Campinarana", "3.5", hab1$vt)
hab1$vt <- gsub("Campo rupestre|Campo de Altitude|Campo de Altitude", "4.7", hab1$vt)
hab1$vt <- gsub("Floresta Estacional Decidual|Chaco", "1.5", hab1$vt)
hab1$vt <- gsub("Floresta de Várzea|Floresta de Igapó", "1.8", hab1$vt)
hab1$vt <- gsub("Vegetação Sobre Afloramentos Rochosos", "6", hab1$vt)
hab1$vt <- gsub("Área Antrópica", "14.5", hab1$vt)
hab1$vt <- gsub("Manguezal", "1.7", hab1$vt)
hab1$vt <- gsub("Campo de Várzea", "4.6", hab1$vt)
hab1$vt <- gsub("Palmeiral", "2.1", hab1$vt)
hab1$vt <- gsub("Vegetação Aquática", "4.6", hab1$vt)
hab1$vt <- sapply(strsplit(hab1$vt,"\\|"), 
                  function(x) paste(sort(as.double(unique(x)), na.last = T), collapse = "|", sep=""))
hab1$vt[hab1$vt %in% "NA"] <- "1.6"
names(hab1)[grepl("^vt", names(hab1))] <- "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup"
  
#GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName:	(e.g. Habitat description, e.g. Shrubland ->Shrubland - Temperate)
hab1$vt <- hab1$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup
hab1$vt <- gsub("^1.5$", "Subtropical/Tropical Dry Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1.6$", "Subtropical/Tropical Moist Lowland Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1.7$", "Subtropical/Tropical Mangrove Forest Vegetation Above High Tide Level", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1.8$", "Subtropical/Tropical Swamp Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^2.1$", "Dry Savanna", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^3.5$", "Subtropical/Tropical Dry Shrubland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^4.6$", "Subtropical/Tropical Seasonally Wet/Flooded Lowland Grassland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^4.7$", "Subtropical/Tropical High Altitude Grassland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^6$", "Inland Rocky Areas", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^14.5$", "Urban Areas", hab1$vt, perl = TRUE)
names(hab1)[grepl("^vt", names(hab1))] <- "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName"


#GeneralHabitats.GeneralHabitatsSubfield.suitability:	Suitability - defaults to suitable if not provided

#GenerationLength.range
names(hab1)[grepl("^GL", names(hab1))] <- "GenerationLength.range"

#FemaleMaturitySize.size:	Size at Maturity (in cms): Female
#MaleMaturitySize.size:	Size at Maturity (in cms): Male
hab1$FemaleMaturitySize.size <- hab1$MaleMaturitySize.size <- hab1$DBH

#MaxSize.size	Maximum Size (in cms)
names(hab1)[grepl("^MaxHeight", names(hab1))] <- "MaxSize.size"
hab1$MaxSize.size <- hab1$MaxSize.size * 100

#System.value:	System (can have more than one value, separated by pipes ("|"))
names(hab1)[grepl("^habitat.reflora", names(hab1))] <- "System.value"
hab1$System.value <- gsub("Terrícola|Hemiparasita|Rupícola|Epífita|Hemiepífita|Hemiparasita", "Terrestrial", hab1$System.value)
hab1$System.value <- gsub("Aquática", "Freshwater (=Inland waters)", hab1$System.value)
hab1$System.value <- gsub("Desconhecido|NA", "Unknown", hab1$System.value)
hab1$System.value <- sapply(strsplit(hab1$System.value,"\\|"), 
                            function(x) paste(sort(unique(x), na.last = T, decreasing = TRUE), collapse = "|", sep=""))
hab1$System.value[hab1$System.value %in% "Freshwater (=Inland waters)"] <- "Terrestrial|Freshwater (=Inland waters)"


#PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup: i.e. the Code for the Plants e.g. for Tree - Large , its TL.
hab1$pgf <- hab1$GF
hab1$pgf <- gsub("tree_unknown", "1", hab1$pgf)
hab1$pgf <- gsub("large_tree", "2", hab1$pgf)
hab1$pgf <- gsub("small_tree", "3", hab1$pgf)
hab1$pgf <- gsub("large_shrub", "5", hab1$pgf)
hab1$pgf[hab1$life.form %in% c("Succulent tree","Succulent treelet")] <- "19"
hab1$pgf[hab1$life.form %in% "Tree fern"] <- "20"

#PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName:	Provide the description here if you have it e.g. Tree - Large
hab1$pgf1 <- hab1$pgf
hab1$pgf1 <- gsub("^1$", "Tree - size unknown", hab1$pgf1)
hab1$pgf1 <- gsub("^2$", "Tree - large", hab1$pgf1)
hab1$pgf1 <- gsub("^3$", "Tree - small", hab1$pgf1)
hab1$pgf1 <- gsub("^5$", "Shrub - large", hab1$pgf1)
hab1$pgf1 <- gsub("^19$", "Succulent - tree", hab1$pgf1)
hab1$pgf1 <- gsub("^20$", "Fern", hab1$pgf1)
names(hab1)[grepl("^pgf$", names(hab1))] <- "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup"
names(hab1)[grepl("^pgf1$", names(hab1))] <- "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName"

#GenerationLength.justification
standard.texts <- data.frame(GL1 = "Field measurements of generation length are largely missing for Neotropical trees. Therefore, an approach based on tree functional types was used to set proxies of species generation lengths based on their growth form and ecological group. Values of the proxies ranged from 10 (for pioneer shrubs) to 80 years (large, climax trees), which are reasonable for tropical trees obtained using IUCN approximative method number 3 (IUCN 2019, p.29). This species was classified as ",
                             GL2 = "More details on this approximate method to define species generation lengths can be found in Lima et al. (XXXX).",
                             stringsAsFactors = FALSE)
hab1$GenerationLength.justification <- paste0("standard_text_GL1", hab1$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName,
                                                 ", regarding its growth form (maximum height of ", hab1$MaxSize.size/100," m), and as ",
                                                 hab1$ecol.group,", regarding its ecological group. standard_text_GL2")

## Any missing info?
#Any missing GL?
hab1$GenerationLength.range[is.na(hab1$GenerationLength.range)] <- 50
#Any missing p?
hab1$p.ci[is.na(hab1$p.est)] <- "0.44-0.45"
hab1$p.est[is.na(hab1$p.est)] <- 0.45


##Organizing and saving
hab2 <- hab1[order(hab1$species.correct),]
cols <- c("taxon_id","Name_submitted",
          "HabitatDocumentation.narrative",
          "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup",
          "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName",
          #"GeneralHabitats.GeneralHabitatsSubfield.suitability",
          "GenerationLength.range",
          "GenerationLength.justification",
          "MaleMaturitySize.size",                                           
          "FemaleMaturitySize.size", 
          "MaxSize.size",
          "System.value",
          "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup",
          "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName",
          "p.est","p.ci",
          #some extra info tha my be useful for species descriptions
          "habito","life.form","DBH","wsg","ecol.group","dispersal.syndrome","SeedMass_g"
)
hab2 <- hab2[, cols]
write.csv(standard.texts, "data/threat_standard_texts.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hab2, "data/threat_habitats.csv", row.names = FALSE, fileEncoding = "UTF-8")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
#############################H
#### PREVIOUS ASSESSMENTS ####
#############################H
rm(list=ls())


## Red List IUCN key: 814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25
??rredlist

## Getting the list of species in the Atlantic Forest
prev.assess <- readRDS("data/assess_iucn_spp.rds")
prev.assess <- prev.assess[,1:3]

## Getting the IUCN assessments for Brazil (CNCFlora) - National level
tmp = flora::get.taxa(prev.assess$species.correct2, replace.synonyms = FALSE, life.form = TRUE)
tmp$search.str = prev.assess$species.correct2
tmp1 <- merge(prev.assess, tmp, by.x = "species.correct2", by.y = "search.str", all.x = TRUE, sort=FALSE)
#table(prev.assess$species.correct2 == tmp1$species.correct2)
prev.assess$status.reflora <- tmp1$threat.status

## Getting the global IUCN assessments (IUCN) 
#Citation: IUCN 2020. IUCN Red List of Threatened Species. Version 2020-1 <www.iucnredlist.org>
iucn <- read.csv("IUCN_2020_assessments.csv", as.is = TRUE, na.strin = c(NA,""," "))
iucn$species.correct2 <- sapply(strsplit(iucn$scientificName," "), function(x) paste(x[1], x[2],sep=" "))
tmp <- merge(prev.assess, iucn,  by= "species.correct2", all.x = TRUE)
tmp <- tmp[,c("taxon_id","species.correct2","assessmentId","internalTaxonId","redlistCategory","redlistCriteria","yearPublished","assessmentDate","criteriaVersion",
          "language","rationale","habitat","threats","population","populationTrend","range","useTrade","systems","conservationActions",
          "realm","yearLastSeen","possiblyExtinct","possiblyExtinctInTheWild","scopes")]
table(prev.assess$species.correct2 == tmp$species.correct2)
prev.assess <- cbind.data.frame(prev.assess, tmp,
                                stringsAsFactors = FALSE)
saveRDS(prev.assess, "data/previous_assessments_spp.rds")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
##################################################
#### RANGE DESCRIPTION AND COUNTRY OCCURRENCE ####
##################################################

?red::countries


#####################################################################################################################################################################H
#####################################################################################################################################################################H
########################################################H
#### USE AND TREAD, THREATS AND CONSERVATION ACTIONS ####
########################################################H

##See TreeCo use database

