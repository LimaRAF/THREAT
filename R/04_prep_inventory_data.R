##############################################################
##############################################################
#### PREPARING INVENTORY DATA FOR THREAT IUCN ASSESSMENTS ####
##############################################################
##############################################################
rm(list=ls())

#### LOADING PACKAGES ####
library(data.table)

#######################H
#### INVENTORY DATA ####
#######################H

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

## List of synonyms 
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$ï..status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$ï..status[i]
  
  if (st.i == "replace")
    trees$species.correct[trees$species.correct %in% sp.i] <- rpl.i
}

## READING AND EDITING THE ABUNDANCE AND PLOT DATA ##
#trees1 <- abundance_data_prep(trees, surveys1, linhas= sites)$abundances_per_ordem
trees1 <- aggregate(trees$N, list(trees$ordem, trees$species.correct), sum, na.rm=TRUE)
names(trees1)[1:3] = c("ordem","species.correct", "N")
trees1$RecordID <- aggregate(trees$RecordID, list(trees$ordem, trees$species.correct), function(x) x[1])$x
tmp <- aggregate(trees$family, list(trees$species.correct), unique)
names(tmp) = c("species.correct","family")
trees1 <- dplyr::left_join(trees1, tmp, by = "species.correct")
trees1 <- trees1[,c("RecordID","ordem","species.correct","family","N")]

##Binding other important info
abund = trees[,c("RecordID","taxon.rank","tax_check","notes","accession","coletor","number","record","speciesLink","status.spLink","DetBy","DetDate"),]
trees1 <- dplyr::left_join(trees1, abund, by = "RecordID")

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
factors.idx = c(1, which(sapply(1:ncol(sites2), function(x) is.character(sites2[,x])) == T)) 
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
dados2 = merge(dados, surveys.final, by="ordem",all.x=TRUE)

##Multiple-site plots
dados3 = merge(dados1, surveys2,by="ordem",all.x=TRUE)

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
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]
# spp.af = readRDS("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/spp.af.rds")
# spp.af.true <- spp.af[af.check2 %in% "TRUE", species.correct2]
# #Adding some missing names (don't know why they were excluded; probably due to missing valid determinations or coordinates)
# miss.spp = c("Agonandra brasiliensis","Aiouea bracteata","Aspidosperma brasiliense",
#              "Aspidosperma quirandy","Cochlospermum vitifolium","Ficus pakkensis","Handroanthus selachidentatus",
#              "Handroanthus spongiosus","Ilex cognata","Luetzelburgia purpurea","Mezilaurus synandra","Mimosa caesalpiniifolia",
#              "Ocotea grandifructa","Persea punctata","Piptadenia killipii","Piranhea securinega","Quillaja lancifolia",
#              "Piptocarpha regnellii","Pleroma echinatum","Pleroma vimineum", 
#              "Tovomita longifolia","Trischidium decipiens")
trees.final <- trees.final[trees.final$species.correct2 %in% tax$species.correct2,]

## Removing records already in the herbarium data
tombos <- readRDS("data/threat_occ_data.rds")
tombos <- tombos$dup.ID1
tombos <- as.character(unlist(strsplit(tombos, "\\|")))
trees.final <- trees.final[!trees.final$numTombo %in% tombos,]

## Removing sites with poor coordinates ## 
trees.final <- 
  trees.final[trees.final$confiabilidade %in% c("precisa","exata","boa","boa|precisa","precisa|exata","precisa|exata|boa"),]

## Removing missing coordinates ##
trees.final <- 
  trees.final[!is.na(as.double(trees.final$lat1)) | !is.na(trees.final$long1),] 

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

#Creating the typeStatus column
trees.final$typeStatus <- NA
trees.final$typeStatus[grepl("paratipo|paratype|type specimen", trees.final$notes)] <-
  trees.final$notes[grepl("paratipo|paratype|type specimen", trees.final$notes)]

## Saving ##
saveRDS(trees.final, file = "data/threat_inventory_data.rds", compress = "gzip")

#### OBTAINING THE NUMBER OF OCCURRENCES PER SPECIES FROM TREECO ####
tmp <- table(trees.final$species.correct2)
tmp1 <- data.frame(species.correct2 = names(tmp), treeco.occs = as.double(tmp), 
                   stringsAsFactors = FALSE)
resultado <- dplyr::left_join(tax, tmp1)
resultado <- resultado[order(resultado$species.correct2),]
saveRDS(resultado, "data/inventory_spp_data.rds")
rm(trees,trees1,trees.final)
