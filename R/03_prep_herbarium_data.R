##################################################################
##################################################################
#### PREPARING HERBARIUM DATA FOR THE THREAT IUCN ASSESSMENTS ####
##################################################################
##################################################################
rm(list=ls())
gc()

#### LOADING PACKAGES ####
require(data.table)
require(plantR)
source("R/other_functions.R")

# #### LOADING THREAT HERBARIUM DATA ###
# ## Loading the paths to occurrence data 
# ## >200 csv files that make almost 4 GB in size, so the output of the codes below was saved in the data folder 
# paths <- list.files("./data/data-raw/occurrence_data/", full.names = TRUE)
# 
# ## Loading the different collections
# lista <- vector("list",length(paths))
# for (i in 1:length(paths)) {
#   path.csv <- paths[i]
#   #tmp = read.csv(path.csv,as.is=TRUE,na.strings=c(""," ","NA"))
#   tmp <- fread(path.csv, na.strings=c("","NA"))
#   lista[[i]] <- tmp
# }
# # Altering the class of the date vector to characters
# for (i in 1:length(lista)) {
#   lista.i <- lista[[i]]  
#   classes <- sapply(lista.i, class)
#   if (any(classes[[15]] != "character")) {
#     lista.i[ , dateIdentified := as.character(dateIdentified)]
#     lista[[i]] <- lista.i
#   } else { next }
# }
# ## Merging the different family records
# oc.data <- rbindlist(lista) #merging data from all csvs
# rm(lista, lista.i, tmp, paths, path.csv, i, classes)
# saveRDS(oc.data, "./data/data-raw/edited_occurrence_data.rds")
# gc()

#################################H
#### MANAGING OCCURRENCE DATA ####
#################################H

## Loading the combined/merged raw occurrence data (over 6 million initial records)
oc.data <- readRDS("./data/data-raw/edited_occurrence_data.rds")
gc()

#### REMOVING SPECIES NOT IN THE ATLANTIC FOREST CHECKLIST ###
## Creating a vector with the names of specimens determined at species and infra-specific levels
oc.data[,species.correct2 := sapply(strsplit(species.correct1," "), 
                                    function(x) paste(x[1], x[2],sep=" ")),]
setkeyv(oc.data, "species.correct2") ## setting 'species.correct2' as key to the data.table (makes computations faster)

#### SOME LAST MINUTE, NAMES TO BE REPLACED BEFORE FILTERING ###
## Changes not taken into account in th occurrence data editing by Lima et al. (2020)
oc.data[loc.correct %like% "brazil" & species.correct2 %in% "Ixora grandifolia", 
        species.correct2 := "Ixora muelleri"]
oc.data[species.correct2 %like% "Zanthoxylum" & numTombo %in% c("jpb_54464","ase_6935","ase_29975","ase_30740","ase_35130"), 
        species.correct2 := "Zanthoxylum unifoliolatum"]
oc.data[species.correct2 %in% "Ocotea lucida" & coletor.last.name %in% "gardner", 
        species.correct2 := "Ocotea brachybotrya"]
locais <- c("brazil_espirito santo|brazil_parana|brazil_rio janeiro|brazil_santa catarina|brazil_sao paulo")
oc.data[species.correct2 %in% "Parinari excelsa" &  loc.correct %like% locais, 
        species.correct2 := "Parinari brasiliensis"]
oc.data[species.correct2 %in% "Solanum pseudocapsicum" & 
          species.correct %in% "Solanum pseudoquina", species.correct2 := "Solanum pseudoquina"]

#### CHECAR AQUI QUANDO FOR RODAR PELA ULTIMA VEZ ####
# add year col.year (1877) e det.year (2010) for ny_375712 
oc.data[species.correct2 %in% "Pouteria stenophylla" & numTombo %in% "ny_375712",
        ano := "1877"]
oc.data[species.correct2 %in% "Pouteria stenophylla" & numTombo %in% "ny_375712",
        year := 1877]
oc.data[species.correct2 %in% "Pouteria stenophylla" & numTombo %in% "ny_375712",
        ano.det := "2010"]
oc.data[species.correct2 %in% "Pouteria stenophylla" & numTombo %in% "ny_375712",
        ano.det1 := "2010"]
#complete det. and col. year, and det Name for the dup.ID: f_72195f|g_439562|us_702351
oc.data[numTombo %in% c("f_72195f", "g_439562", "us_702351"), dup.ID := "f_72195f|g_439562|us_702351"]
oc.data[dup.ID %in% "f_72195f|g_439562|us_702351", ano := "1823"]
oc.data[dup.ID %in% "f_72195f|g_439562|us_702351", year := 1823]
oc.data[dup.ID %in% "f_72195f|g_439562|us_702351", determinador.name1 := "Pennington, T.D."]
oc.data[dup.ID %in% "f_72195f|g_439562|us_702351", ano.det1 := "1986"]
gc()

## Filtering for the AF preliminary list of species names
spp.af <- readRDS("data/threat_af_spp_list_preliminary.rds")
# spp.af$species.correct2[!spp.af$species.correct2 %in% unique(oc.data$species.correct2)]
# oc.data[, uniqueN(species.correct2)] # 21643 species; now 21793
oc.data <- oc.data[species.correct2 %in% spp.af$species.correct2,] # only occurrences for species inside the AF
oc.data[, uniqueN(species.correct2)] # before 5107 species; now 5081
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

##List of new synonyms (not taken into account by Lima et al. (2020))
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace") {
    oc.data[species.correct2 %in% sp.i, status := "ok"]
    oc.data[species.correct2 %in% sp.i, species.correct2 := rpl.i]
  }
}
#oc.data[,uniqueN(species.correct2)] # Before 5068 species; now 4953


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
toto1 - dim(oc.data)[1]; 100*(toto1 - dim(oc.data)[1])/toto1 ## 728,411 (23.65%) removed
# oc.data[,uniqueN(species.correct2)] # 5095 species; now 4953
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING DATA NOT GEOGRAPHICALLY VALIDATED ###
## But keeping records for species without coordinates but know to be in the Atlantic Forest
gc()
paises.adm0 <- paste(c("anguilla", "aruba", "curacao", "falkland islands",
                       "saint-martin", "sint maarten", "st\\. barthelemy"), collapse = "|")
paises.adm1 <- paste(c("antigua & barbuda", "bahamas","barbados","belize",
                       "bermuda", "caribbean netherlands", "british virgin islands", 
                       "cayman islands", "dominica", "grenada","jamaica", "martinique",
                       "montserrat", "puerto rico","st\\. kitts & nevis",
                       "st\\. lucia","st\\. vincent & grenadines",
                       "trinidad & tobago","turks & caicos islands",
                       "u\\.s\\. virgin islands"), collapse = "|")
toto2 = dim(oc.data)[1]
oc.data <- oc.data[geo.check1 %like% "ok_county|ok_locality" | 
                     (geo.check1 %like% "ok_state" & af.check2 == TRUE) |
                     (geo.check1 %like% "ok_state" & loc.correct %like% paises.adm1) |
                     (geo.check1 %like% "ok_country" & loc.correct %like% paises.adm0)]
#oc.data[,uniqueN(species.correct2)] # 5095 species; now 4953
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
toto2 - dim(oc.data)[1]; 100*(toto2 - dim(oc.data)[1])/toto2 ## 955,427 (40.64% of the non-duplicated) removed
100*(toto2 - dim(oc.data)[1])/toto1 ## extra 31.02% in respect to all records
#oc.data[,uniqueN(species.correct2)] # 5095 species; now 4953
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
toto3 - dim(oc.data)[1]; 100*(toto3 - dim(oc.data)[1])/toto2 ## 575,655 (24.48% of remaining records) removed
100*(toto3 - dim(oc.data)[1])/toto1 ## extra 18.69% in respect to all records
# oc.data[,uniqueN(species.correct2)] # 5095 species; now 4953
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### REMOVING SPATIAL OUTLIERS  ###
## GETTING SPECIES INFORMATION: TBC, geographical range and cultivation #### Getting species information and generating the vectors for the groups of species
spp = read.csv("data/data-raw/DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
tmp = unique(oc.data[,species.correct2,])
spp = spp[spp$SpeciesKTSA %in% tmp,]
uso = read.csv("data/data-raw/plant_uses.csv", na= c(""," ",NA))
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
oc.data <- oc.data[is.na(true.out) | true.out %in% FALSE] # removing 3142 true outliers
# oc.data[,uniqueN(species.correct2)] # 5095 species; now 4953
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

## REMOVING PROBABLE OUTLIERS ##
gc()
#2- remover probable outliers para cultivated species:  classic Mahalanobis distance >2.5 & robust (Mahalanobis) distance>11.5
tmp <- data.table(oc.data[species.correct2 %in% cult])
#table(tmp$probable.out, useNA = "always")
tmp1 <- data.table(tmp[probable.out %in% TRUE])
oc.data <- oc.data[!numTombo %in% tmp1$numTombo]
toto4 - dim(oc.data)[1]; 100*(toto4 - dim(oc.data)[1])/toto3 ## 3181 (0.23% of remaining records) removed
100*(toto4 - dim(oc.data)[1])/toto1 ## extra 0.1% in respect to all records
# oc.data[,uniqueN(species.correct2)] # 5095 species; now 4953
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

## SAVING THE DATA FOR EXTRACTING THE COUNTRY OF OCCURRENCE AND ALTITUDE
tmp.data <- oc.data[,c("latitude.work1", "longitude.work1", 
                       "loc.correct","species.correct2"),]
tmp.data[ , pais := gsub("_.*", "", loc.correct, perl = TRUE)]
f <- function (x) { sapply(strsplit(x, "_", fixed = TRUE), function(x) x[2]) }
tmp.data[ , estado := f(loc.correct)]
tmp.data[ , loc.correct := paste(pais, estado, sep = "_")]
tmp.data[ , c("pais", "estado") := NULL]
dim(tmp.data) # 816677 records
saveRDS(tmp.data, file = "data/threat_species_by_country.rds", compress = "gzip")
rm(tmp.data)

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
# oc.data[,uniqueN(species.correct2)] # 5095 species
# table(is.na(oc.data$latitude.work1)) # TRUE: with missing coordinates (but confirmed in the AF) for some species

#### GETTING SPECIES INFORMATION: TBC, geographical range and cultivation ### 
# Getting species information and generating the vectors for the groups of species
spp = read.csv("data/data-raw//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
spp = spp[spp$SpeciesKTSA %in% oc.data$species.correct2,]

## Preliminary editing of the list
spp$nomes = spp$species.reflora
#replacing names not found in ReFlora by the TNRS/Hans taxonomy
spp$nomes[is.na(spp$nomes)] = spp$Accepted_name[is.na(spp$nomes)]
#replacing names not found in ReFlora and TNRS by the Plant List
spp$nomes[is.na(spp$nomes)] = spp$NewTaxon.TPL[is.na(spp$nomes)]

## Getting the list of species per group
tbc = spp$nomes[spp$TBC %in% "TBC"] # taxa with low taxonomic complexity
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
#### CHECK HERE: what should be the best cut-off to set the most frequent taxonomists? Leaving the 20 that was fixed before... ####
taxonomists = names(tmp[tmp >= 20]) # quantile(tmp, prob = seq(0.2,0.5,0.05)) 
taxonomists = taxonomists[!taxonomists %in% c("Anonimo","Anonymous","Semdeterminador")]
oc.data[determinador.name %in% taxonomists & tax.check2 %in% "FALSE", tax.check2 := "TRUE_OTHER",]

#Validating all occurrences of TBC
oc.data[species.correct1 %in% tbc & tax.check2 %in% c("FALSE","cannot_check"), tax.check2 := "TRUE_TBC",]
# oc.data[,uniqueN(species.correct2)] # 5095 species
# table(is.na(oc.data$latitude.work1))

#some last minute validation to avoid the removal of species from the checklist
#Edits not previously taken into account by Lima et al. (2020)
oc.data[species.correct2 %in% "Agonandra brasiliensis" & determinador.name1 %like% c("Hiepko|Heipko|Marquete|Groppo"), tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Aiouea bracteata" & determinador.name1 %like% c("Baitello|Lorea"), tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Quillaja lancifolia" & determinador.name1 %in% c("Luebert, F."), tax.check2 := "TRUE" ]
oc.data[species.correct2 %like% "Handroanthus" & determinador.name1  %like% c("Santo, F.S.E.|Espirito Santo, F.S.|Santo, F.E.|Silva-Castro|ilva, M.M."), tax.check2 := "TRUE" ]
oc.data[species.correct2 %in% "Pouteria stenophylla" & determinador.name1 %like% "Palazzo", tax.check2 := "TRUE"]
oc.data[species.correct2 %in% "Pradosia glaziovii" & determinador.name1 %like% "Terra-Araujo", tax.check2 := "TRUE"]

## VALIDATING VOUCHER INFO FROM REFLORA
vouchers <- readRDS("data/vouchers_reflora.rds")

#### CHECK HERE ####
#By the 'Número de tombo'
for (i in 1:length(unique(vouchers$tax))) {
  spi <- unique(vouchers$tax)[i]
  vouc.list <- vouchers$numTombo[vouchers$tax %in% spi]
  vouc.list <- vouc.list[!is.na(vouc.list)]
  if (length(vouc.list) > 0)
    oc.data[species.correct2 %in% spi & 
              numTombo %in% vouc.list, tax.check2 := "TRUE" ]
}

#By coletor last name and number
oc.data[ , tmp.collect := paste(coletor.last.name, coletor.number, sep="_"),]
oc.data[coletor.number %in% "s.n." , tmp.collect := NA_character_,]
for(i in 1:length(unique(vouchers$tax))) {
  spi <- unique(vouchers$tax)[i]
  lst.name <- tolower(vouchers$recordBy[vouchers$tax %in% spi])
  col.numb <- vouchers$colNumber[vouchers$tax %in% spi]
  comb.collect <- paste(lst.name, col.numb, sep = "_")
  oc.data[species.correct2 %in% spi & 
            tmp.collect %in% comb.collect, tax.check2 := "TRUE" ]
}

#### TRYING TO OBTAIN INFORMATION ON CONSERVATION UNITS FROM LOCALITY INFO ###
#Editing the localities
locals <- oc.data$locality
# unwanted_array = list('Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 
#                       'Ç'='C', 'È'='E', 'É'='E', 'Ê'='E', 'Ë'='E', 
#                       'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 
#                       'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 
#                       'Ù'='U', 'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 
#                       'Þ'='B', 'ß'='S', 
#                       'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 
#                       'ç'='c', 'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 
#                       'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 
#                       'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
#                       'ö'='o', 'ø'='o', 'ü'='u', 'ù'='u', 'ú'='u', 
#                       'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )
# locals <- chartr(paste(names(unwanted_array), collapse=''), paste(unwanted_array, collapse=''), locals)
locals <- plantR:::rmLatin(locals)
locals <- gsub('\\.\\.',".", locals, perl = TRUE)
locals <- gsub(' - | -|- ',"-", locals, perl = TRUE)
locals <- gsub('^-- ',"", locals, perl = TRUE)
locals <- gsub(' --$',"", locals, perl = TRUE)
locals <- gsub('^\\.|^,',"", locals, perl = TRUE)
locals <- gsub('\\.$|,$',"", locals, perl = TRUE)
locals <- gsub('^-\\. ','', locals, perl = TRUE)
locals <- gsub('-\\.$','', locals, perl = TRUE)
locals <- gsub('^-','', locals, perl = TRUE)
locals <- gsub('-$','', locals, perl = TRUE)
locals <- stringr::str_squish(locals)
oc.data$locality1 <- locals

## Generating the file with the conservation unit information
## Some manual editing was performed here, so we are using the edited version loaded below
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
ucs <- readRDS("data/UCs.rds")
oc.data <- dplyr::left_join(oc.data, ucs, by= "locality1")

#### SAVING THE ESSENTIAL DATA FOR THE ASSESSMENTS ###
## Filtering only the necessary columns
oc.data1 <- oc.data[,c("latitude.work1","longitude.work1","species.correct2","family.correct1",
                       "coletor.name","coletor.number","ano","numTombo","dup.ID1",
                       "typeStatus","determinador.name1","ano.det1","tax.check2",
                       "geo.check1","origin.coord1","af.check2","UC"),]
saveRDS(oc.data1, file = "data/threat_occ_data.rds", compress = "gzip")

#### OBTAINING THE NAME AND NUMBER OF OCCURRENCES (TOTAL AND SPATIALLY-UNIQUE) PER SPECIES ###
resultado <- oc.data[!duplicated(oc.data$species.correct2), 
                     c("family.correct1", "species.correct2")]
#table(resultado$species.correct2 == names(table(oc.data$species.correct2)))
resultado$total.occs <- as.double(table(oc.data$species.correct2))
tmp <- merge(resultado, 
             data.frame(species.correct2 = names(extra.neotrop.spp), 
                        extra.neotrop.occurs = as.double(extra.neotrop.spp)), 
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
resultado$total.occs <- 
  as.double(apply(resultado[,c("total.occs","non.neotrop.occs")], 
                  1, sum, na.rm = TRUE))

## SAVING THE FINAL LIST OF SPECIES ##
taxon_id <- paste0("sp", 1:dim(resultado)[1])
resultado <- cbind.data.frame(internal_taxon_id = taxon_id, resultado, stringsAsFactors = FALSE)
saveRDS(resultado, "data/herbarium_spp_data.rds")
