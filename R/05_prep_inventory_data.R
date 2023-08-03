##############################################################
##############################################################
#### PREPARING INVENTORY DATA FOR THREAT IUCN ASSESSMENTS ####
##############################################################
##############################################################
rm(list=ls())

#### LOADING PACKAGES AND FUNCTIONS ####
require(data.table)
source("R/99_functions.R")


####################################################
#### ADDING SPECIES RECORDS FROM INVENTORY DATA ####
####################################################

# Loading TreeCo species records
trees.final <- readRDS("data/treeco_inventory_records.rds")

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
         "species.correct","accession","record","year","ordem") 
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
rm(list = ls())
