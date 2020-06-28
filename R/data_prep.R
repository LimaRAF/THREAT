#####################################################
#### PREPARING DATA FROM THREAT IUCN ASSESSMENTS ####
#####################################################
#rm(list=ls())

#### LOADING PACKAGES ###
library(data.table)
library(rgeos)
library(dplyr)
source("R/suggestions_for_ConR.r")

########################################
#### SIPLIFIED NEOTROPICAL COUNTOUR ####
########################################
neotrop <- readRDS("E:/ownCloud/W_GIS/Am_Lat_ADM_GADM_v3.6/gadm36_Neotrop_0_sp_simplified.rds")
neotrop <- gBuffer(neotrop, byid=TRUE, width=0)

#projecting and simplifying the neotropical contours
neotrop.simp <- gSimplify(neotrop,tol=0.05)
neotrop.simp <- gBuffer(neotrop.simp, byid=TRUE, width=0)
saveRDS(neotrop.simp, file = "data//Contour_Neotrop_simplified.rds")


#####################################################################################################################################################################
#####################################################################################################################################################################
########################
#### OCCURRENCE DATA ###
########################

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

#################################
#### MANAGING OCCURRENCE DATA ###
#################################

#### REMOVING SPECIES NOT IN THE ATLANTIC FOREST CHECKLIST ####
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
             "Aspidosperma quirandy","Cochlospermum vitifolium","Ficus pakkensis","Handroanthus selachidentatus",
             "Handroanthus spongiosus","Ilex cognata","Luetzelburgia purpurea","Mezilaurus synandra","Mimosa caesalpiniifolia",
             "Ocotea grandifructa","Persea punctata","Piptadenia killipii","Piranhea securinega","Quillaja lancifolia",
             "Tovomita longifolia","Trischidium decipiens")
oc.data <- oc.data[species.correct2 %in% c(spp.af.true, miss.spp),] # occurrences inside the AF
# oc.data[,uniqueN(species.correct2)] # 5183 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING DUPLICATES ####
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
toto1 - dim(oc.data)[1]; 100*(toto1 - dim(oc.data)[1])/toto1 ## 732,113 (23.64%) removed
# oc.data[,uniqueN(species.correct2)] # 5183 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING DATA NOT GEOGRAPHICALLY VALIDATED #### 
## But keeping records for species without coordinates but know to be in the Atlantic Forest
gc()
toto2 = dim(oc.data)[1]
oc.data <- oc.data[geo.check1 %like% "ok_county|ok_locality" | 
                     (geo.check1 %like% "ok_state" & af.check2 == TRUE)]
# oc.data[,uniqueN(species.correct2)] # 5165 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING MISSING COORDINATES ####
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
toto2 - dim(oc.data)[1]; 100*(toto2 - dim(oc.data)[1])/toto2 ## 970,170 (41.01% of the non-duplicated) removed
100*(toto2 - dim(oc.data)[1])/toto1 ## extra 31.32% in respect to all records
#oc.data[,uniqueN(species.correct2)] # 5183 species
#table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING SPATIAL DUPLICATES ####
gc()
toto3 = dim(oc.data)[1]
setkeyv(oc.data, "species.correct2") ## setting 'species.correct2' as key to the data.table (makes computations faster)
oc.data[ , coord.string := as.character(paste(latitude.work1, longitude.work1, sep="_")),]
oc.data[ , spat.dups := duplicated(coord.string, incomparables = "NA_NA"), by = species.correct2]
oc.data <- oc.data[spat.dups == FALSE,,]
oc.data[,coord.string:=NULL] 
oc.data[,spat.dups:=NULL] 
toto3 - dim(oc.data)[1]; 100*(toto3 - dim(oc.data)[1])/toto2 ## 573,557 (24.25% of remaining records) removed
100*(toto3 - dim(oc.data)[1])/toto1 ## extra 18.52% in respect to all records
# oc.data[,uniqueN(species.correct2)] # 5183 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING SPATIAL OUTLIERS  ####
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
# oc.data[,uniqueN(species.correct2)] # 5183 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

## REMOVING PROBABLE OUTLIERS ##
gc()
#2- remover probable outliers para cultivated species:  classic Mahalanobis distance >2.5 & robust (Mahalanobis) distance>11.5
tmp <- data.table(oc.data[species.correct2 %in% cult])
#table(tmp$probable.out, useNA = "always")
tmp1 <- data.table(tmp[probable.out %in% TRUE])
oc.data <- oc.data[!numTombo %in% tmp1$numTombo]
toto4 - dim(oc.data)[1]; 100*(toto4 - dim(oc.data)[1])/toto3 ## 3117+41 (0.23% of remaining records) removed
100*(toto4 - dim(oc.data)[1])/toto1 ## extra 0.1% in respect to all records
# oc.data[,uniqueN(species.correct2)] # 5183 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#3- remover non-core specimens for heavily cultivated native species (e.g. Araucaria angustifolia): rob.out.99 == FALSE
## REMOVE OCCURRENCES OUTSIDE THE NEOTROPICS, EXCEPT FOR PANTROPICAL SPECIES?
## REMOVE OCCURRENCE DATA FROM WELL-KNOW BOTANICAL GARDENS GROWING SOUTH AMERICA SPECIES?

#### REMOVING NON-NEOTROPICAL OCCURRENCES ####
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
toto5 - dim(oc.data)[1]; 100*(toto5 - dim(oc.data)[1])/toto4 ## 489 (0.11% of remaining records) removed
100*(toto5 - dim(oc.data)[1])/toto1 ## extra 0.02% in respect to all records
#oc.data[,uniqueN(species.correct2)] # 5183 species
#table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### GETTING SPECIES INFORMATION: TBC, geographical range and cultivation #### Getting species information and generating the vectors for the groups of species
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

#### SEPARATING DATA IDENTIFIED AND NOT IDENTIFIED BY TAXONOMISTS #### 
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
#oc.data[,uniqueN(species.correct2)] # 5165 species
#table(is.na(oc.data$latitude.work1))

#some last minute validation to avoid the removal of species from the checklist
oc.data[species.correct2 %in% "Agonandra brasiliensis" & determinador.name1 %like% c("Hiepko","Heipko","Marquete","Groppo"), tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Aiouea bracteata" & determinador.name1 %like% c("Baitello","Lorea"), tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Quillaja lancifolia" & determinador.name1 %in% c("Luebert, F."), tax.check2 := "TRUE" ]
oc.data[species.correct2 %like% "Handroanthus" & determinador.name1  %like% c("Santo, F.S.E.", "Espirito Santo, F.S.", "Santo, F.E.","Silva-Castro","ilva, M.M."), tax.check2 := "TRUE" ]

#### TRYING TO OBTAIN INFORMATION ON CONSERVATION UNITS FROM LOCALITY INFO ####
#oc.data <- readRDS("threat_occ_data.rds")
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

#### SAVING THE ESSENTIAL DATA FOR THE ASSESSMENTS ####

## Filtering only the necessary columns
oc.data1 <- oc.data[,c("latitude.work1","longitude.work1","species.correct2","family.correct1",
         "coletor.name","coletor.number","ano","numTombo","dup.ID1",
         "typeStatus","determinador.name1","ano.det1","tax.check2",
         "geo.check1","origin.coord1","af.check2","UC"),]
saveRDS(oc.data1, file = "threat_occ_data.rds", compress = "gzip")

#### OBTAINING THE NAME AND NUMBER OF OCCURRENCES (TOTAL AND SPATIALLY-UNIQUE) PER SPECIES ####
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
rm(oc.data,oc.data1,MyData)

#####################################################################################################################################################################
#####################################################################################################################################################################
########################
#### INVENTORY DATA ####
########################
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
tombos <- readRDS("threat_occ_data.rds")
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
oc.data <- readRDS("threat_occ_data.rds")
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
saveRDS(trees.final, file = "threat_inventory_data.rds", compress = "gzip")

#### OBTAINING THE NUMBER OF OCCURRENCES PER SPECIES FROM TREECO ####
tmp <- table(trees.final$species.correct2)
tmp1 <- merge(resultado, 
             data.frame(species.correct2 = names(tmp), treeco.occs = as.double(tmp)), 
             by = "species.correct2", all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species.correct2),]
#table(resultado$species.correct2 == tmp1$species.correct2)
resultado$treeco.occs <- tmp1$treeco.occs
rm(trees,trees1,trees.final)

#####################################################################################################################################################################
#####################################################################################################################################################################
##########################
#### POPULATION SIZES ####
##########################
pop.sizes <- readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//pop.size.est_nmx50_mxd5_idp2.rds")
pop.sizes.2018 <- pop.sizes$`2018`$mean 
rm(pop.sizes)
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
abline(lm(log(resultado$total.occs+1) ~ log(resultado$pop.size.2018+1)), lwd=2, col=2)

## Which specie swe have population sizes which are not in the checklist?
miss.sp <- tmp$species.correct2[!tmp$species.correct2 %in% resultado$species.correct2]
length(miss.sp) # 215 species; most are non-AF species
# tmp <- flora::get.taxa(miss.sp, life.form = TRUE, vegetation.type = TRUE, establishment = TRUE)
# tmp1 <- flora::get_domains(tmp)
# tmp2 <- tmp1[grepl("ata Atl", tmp1$domain),]
# spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
# tmp3 <- merge(tmp2, spp, by.x = "original.search", by.y = "SpeciesKTSA")
# write.csv(tmp3, "tmp.csv")

saveRDS(resultado, "assess_iucn_spp.rds")
#####################################################################################################################################################################
#####################################################################################################################################################################
##############################
#### PREVIOUS ASSESSMENTS ####
##############################
rm(list=ls())

## Getting the list of species in the Atlantic Forest
prev.assess <- readRDS("assess_iucn_spp.rds")
prev.assess <- prev.assess[,1:2]

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
tmp <- tmp[,c("assessmentId","redlistCategory","redlistCriteria","yearPublished","assessmentDate","criteriaVersion",
          "language","rationale","habitat","threats","population","populationTrend","range","useTrade","systems","conservationActions",
          "realm","yearLastSeen","possiblyExtinct","possiblyExtinctInTheWild","scopes")]
#table(prev.assess$species.correct2 == tmp$species.correct2)
prev.assess <- cbind.data.frame(prev.assess, tmp,
                                stringsAsFactors = FALSE)
saveRDS(prev.assess, "previous_assessments_spp.rds")
