##############################################################
#### OBTAINING THE LATEST FBO INFORMATION FOR ALL SPECIES ####
##############################################################
rm(list = ls())

## Loading the THREAT species list
tax <- readRDS("data/threat_af_spp_list_preliminary.rds")
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
tax <- rbind(tax, syn.br[, c("family", "search.str")], use.names = FALSE)
tax <- tax[!duplicated(species.correct2), ]

################################################################################################################################################H
################################################################################################################################################H
##################################H
#### INFO FROM THE FLORA DWC-A ####
##################################H

## Loading the ReFlora DWC-A file (version 393.285, downloaded the 06/05/2021)
# Citation: Brazil Flora Group (2021): Brazilian Flora 2020 project - Projeto Flora do 
#Brasil 2020. v393.274. Instituto de Pesquisas Jardim Botanico do Rio de Janeiro. 
#Dataset/Checklist. doi:10.15468/1mtkaw
path = here::here()
file = "dwca-lista_especies_flora_brasil-v393.285.zip"
path.dwc <- file.path(path, "data", "data-raw", file)
finch::dwca_cache$delete_all()
fbo <- finch::dwca_read(path.dwc, read=TRUE, encoding="UTF-8")$data

## Finding FBO taxonomic information
taxon <- fbo$taxon.txt
taxon$species.correct2 <- 
  stringr::str_squish(paste(taxon$genus, taxon$specificEpithet))
# taxon1 <- taxon[taxon$kingdom %in% "Plantae" &
#                   taxon$taxonRank %in% c("ESPECIE","VARIEDADE","SUB_ESPECIE","FORMA"),]
taxon1 <- taxon[taxon$species.correct2 %in% tax$species.correct2 &
                  taxon$taxonRank %in% "ESPECIE",]
# Getting the best names for homonyms
dups <- (duplicated(taxon1$species.correct2) |  
           duplicated(taxon1$species.correct2, fromLast = TRUE))
dup.spp <- unique(taxon1$species.correct2[dups])
taxon1$remove <- TRUE
for (i in 1:length(dup.spp)) {
  sp.i <- dup.spp[i]
  check_ids <- taxon1$species.correct2 %in% sp.i
  taxon1.i <- taxon1[check_ids,]
  if (length(unique(taxon1.i$nomenclaturalStatus)) > 1) {
    taxon1$remove[check_ids] <- 
      ifelse(taxon1.i$nomenclaturalStatus == "NOME_CORRETO", TRUE, FALSE)
    if (all(taxon1$remove[check_ids] %in% FALSE))
      taxon1$remove[check_ids] <-
        ifelse(taxon1.i$nomenclaturalStatus == "", TRUE, FALSE)
  } 
}
taxon1$remove[taxon1$scientificName == "Myrcia apiocarpa (O.Berg) N.Silveira"] <- 
  TRUE
taxon1 <- taxon1[taxon1$remove,]

## Life forms and vegetation types
prof <- fbo$speciesprofile.txt
prof$lifeForm <- gsub("\\\"", "", prof$lifeForm, perl = TRUE)
prof$lifeForm <- gsub("\\{|\\}", "", prof$lifeForm, perl = TRUE)
lista <- strsplit(prof$lifeForm, "(?<=\\]),", perl = TRUE)
df <- as.data.frame(stringi::stri_list2matrix(lista, byrow = TRUE), stringsAsFactors=FALSE)
names(df) <- c("life.form", "habitat", "vegetation.type")
#Organizing missing data
df[[3]][grepl("vegetationType", df[[2]]) & !grepl("habitat", df[[2]]) & !is.na(df[[2]])] <-
  df[[2]][grepl("vegetationType", df[[2]]) & !grepl("habitat", df[[2]]) & !is.na(df[[2]])]
df[[2]][grepl("vegetationType", df[[2]]) & !grepl("habitat", df[[2]]) & !is.na(df[[2]])] <- NA
df[[3]][grepl("vegetationType", df[[1]]) & !grepl("lifeForm", df[[1]]) & !is.na(df[[1]])] <-
  df[[1]][grepl("vegetationType", df[[1]]) & !grepl("lifeForm", df[[1]]) & !is.na(df[[1]])]
df[[1]][grepl("vegetationType", df[[1]]) & !grepl("lifeForm", df[[1]]) & !is.na(df[[1]])] <- NA

df[[2]][grepl("habitat", df[[1]]) & !grepl("lifeForm", df[[1]]) & !is.na(df[[1]])] <-
  df[[1]][grepl("habitat", df[[1]]) & !grepl("lifeForm", df[[1]]) & !is.na(df[[1]])]
df[[1]][grepl("habitat", df[[1]]) & !grepl("lifeForm", df[[1]]) & !is.na(df[[1]])] <- NA

#Cleaning
df$life.form <- gsub('^lifeForm\\:\\[|\\]$', "", df$life.form, perl = TRUE)
df$habitat <- gsub('^habitat\\:\\[|\\]$', "", df$habitat, perl = TRUE)
df$vegetation.type <- gsub('^vegetationType\\:\\[|\\]$', "", df$vegetation.type, perl = TRUE)

df$life.form <- gsub(',', '|', df$life.form)
df$habitat <- gsub(',', '|', df$habitat)
df$vegetation.type <- gsub(',', '|', df$vegetation.type)

#Binding
prof <- cbind.data.frame(id = prof[,1], df)
taxon1 <- dplyr::left_join(taxon1, prof, by = "id", suffix = c("", ".y"))


## Distribution within Brazil, Endemism and Establishment
dist <- fbo$distribution.txt
dist$endemism <- gsub(',.*', "", dist$occurrenceRemarks)
#Cleaning
dist$endemism[!grepl("endemism", dist$endemism)] <- "{}"
dist$endemism <- gsub('"endemism":', "", dist$endemism)
dist$endemism <- gsub('\\{\\"|\\"\\}|\\"$', "", dist$endemism)
dist$endemism <- gsub('\\{\\}', NA_character_, dist$endemism)

dist$phytogeographicDomain <- gsub('.*\\"phytogeographicDomain\\"\\:', "", dist$occurrenceRemarks)
#Cleaning
dist$phytogeographicDomain <- gsub('^\\[|\\]\\}$', '', dist$phytogeographicDomain)
dist$phytogeographicDomain[grepl("endemism", dist$phytogeographicDomain)] <- "{}"
dist$phytogeographicDomain <- gsub('^\\"|\\"$', '', dist$phytogeographicDomain)
dist$phytogeographicDomain <- gsub('\\",\\"', '|', dist$phytogeographicDomain)
dist$phytogeographicDomain <- gsub('\\{\\}', NA_character_, dist$phytogeographicDomain)

dist1 <- aggregate(cbind(locationID, establishmentMeans, endemism, phytogeographicDomain) ~ id, data = dist, 
                   function(x) paste(unique(x), collapse = "|")) 
names(dist1) <- c("id", "occurrence", "establishment", "endemism", "phytogeographicDomain")
taxon1 <- dplyr::left_join(taxon1, dist1, by = "id", suffix = c("", ".y"))

## References
refs <- fbo$reference.txt
# Removing dates from the citations, if present
refs$date[refs$date %in% c(""," ", NA)] <- "[no date]"
check_these <- is.na(refs$date) & 
  grepl(",[0-9][0-9][0-9][0-9],$| \\([0-9][0-9][0-9][0-9],$", refs$bibliographicCitation)
n.chars <- nchar(refs$bibliographicCitation[check_these])
refs$date[check_these] <- 
  substring(refs$bibliographicCitation[check_these], n.chars - 4)
refs$bibliographicCitation[refs$date == refs$bibliographicCitation] <- ""
citation <- mapply(function(x, y) { gsub(paste0(",",x), ".", y, fixed = TRUE) },
                   refs$date, refs$bibliographicCitation)
ids <- grepl("[0-9]$", refs$bibliographicCitation) & !grepl("\\D", refs$date)
citation[ids] <- mapply(function(x, y) { gsub(paste0(x, "$"), ".", y, perl = TRUE) },
                        refs$date[ids], citation[ids])
citation <- gsub("\\.\\.$", ".", citation)
#Validating dates
refs$date <- gsub("\\.$", "", refs$date)
refs$date <- gsub("\\)$|^\\(", "", refs$date)
good_dates <- (as.double(refs$date) >1000 & as.double(refs$date) <= 2021) |
  refs$date %in% "[no date]" | refs$date %in% "no prelo"
good_dates[is.na(good_dates)] <- FALSE
bad_dates <- sapply(strsplit(refs$date[!good_dates], "-|, |\\/| \\[", perl = TRUE),
                    function(x) any((as.double(x) >1000 & as.double(x) <= 2021)))
bad_dates[is.na(bad_dates)] <- FALSE
good_dates[!good_dates][bad_dates] <- TRUE
#Criando as citações  
refs$reference <- NA
check_these <- !refs$creator %in% c(""," ", NA) & good_dates
refs$reference[check_these] <- 
  paste0(refs$creator[check_these], " ", refs$date[check_these],
         ". ", citation[check_these])
check_these <- !refs$creator %in% c(""," ", NA) & !good_dates
refs$reference[check_these] <- 
  paste0(refs$creator[check_these], " ", refs$bibliographicCitation[check_these])
check_these <- refs$creator %in% c(""," ", NA) & good_dates
refs$reference[check_these] <- 
  paste0(refs$bibliographicCitation[check_these],".")
#Cleaning
refs$reference <- gsub(",$",".", refs$reference, perl = TRUE)
refs$reference <- gsub(",\\.$",".", refs$reference, perl = TRUE)
refs$reference <- stringr::str_squish(refs$reference)
refs$reference[!grepl("\\.$",refs$reference)] <-
  paste0(refs$reference[!grepl("\\.$",refs$reference)],".")
refs$reference <- gsub("\\.\\.$", ".", refs$reference, perl = TRUE)
hist(as.double(refs$date[refs$type %in% "original description"]), 
     nclass = 320, xlim = c(1700, 2022))

#Getting all references for each ID
refs1 <- aggregate(cbind(reference, type) ~ id, data = refs, 
                   function(x) paste(unique(x), collapse = "|")) 
refs1[[2]] <- gsub("^\\||\\|$", "", refs1[[2]], perl = TRUE)
refs1[[3]] <- gsub("^\\||\\|$", "", refs1[[3]], perl = TRUE)
#Binding
taxon1 <- dplyr::left_join(taxon1, refs1, by = "id", suffix = c("", ".y"))

## Vernacular names
vernac <- fbo$vernacularname.txt
vernac$vernacularName <- stringr::str_to_title(vernac$vernacularName)
vernac$vernacular.name <- apply(vernac[,2:4], 1, paste, collapse = "/")
vernac$vernacular.name <- gsub("^\\/|\\/$", "", vernac$vernacular.name, perl = TRUE)
vernac1 <- aggregate(vernacular.name ~ id, data = vernac, 
                     function(x) paste(unique(x), collapse = " | ")) 
#Binding
taxon1 <- dplyr::left_join(taxon1, vernac1, by = "id", suffix = c("", ".y"))


## Editing final info
rplc <- c("ORDEM" = "order", "SUB_FAMILIA" = "subfamily", "FAMILIA" = "family", 
          "GENERO" = "genus", "SUB_ESPECIE" = "subspecies", "ESPECIE" = "species", 
          "VARIEDADE" = "variety",  "CLASSE" = "class", "TRIBO" = "tribe",  
          "DIVISAO" = "division",  "FORMA" = "form")
taxon1$taxon.rank <- stringr::str_replace_all(taxon1$taxonRank, rplc)

taxon1$taxon.status <- gsub("NOME_ACEITO", "accepted", taxon1$taxonomicStatus)
taxon1$taxon.status <- gsub("SINONIMO", "synonym", taxon1$taxon.status)

rplc <- c("NOME_CORRETO_VIA_CONSERVACAO" = "correct via conservation", 
          "NOME_CORRETO" = "correct", 
          "NOME_APLICACAO_INCERTA" = "uncertain application",
          "NOME_LEGITIMO_MAS_INCORRETO" = "legitimate, but incorrect", 
          "NOME_ILEGITIMO" = "illegitimate", 
          "NOME_NAO_VALIDAMENTE_PUBLICADO" = "not validly published", 
          "VARIANTE_ORTOGRAFICA" = "orthographic variant",
          "NOME_NAO_EFETIVAMENTE_PUBLICADO" = "not effectively published", 
          "NOME_MAL_APLICADO" = "misapplied",
          "NOME_REJEITADO" = "rejected")
taxon1$name.status <- stringr::str_replace_all(taxon1$nomenclaturalStatus, rplc)

rplc <- c("NATIVA" = "native", 
          "NATURALIZADA" = "naturalized", 
          "CULTIVADA" = "cultivated")
taxon1$establishment <- stringr::str_replace_all(taxon1$establishment, rplc)

rplc <- c("Endemica" = "endemic", 
          "Não endemica" = "not endemic")
taxon1$endemism <- stringr::str_replace_all(taxon1$endemism, rplc)


## Saving
taxon2 <- taxon1[, c("id", "family", "genus", "specificEpithet", "infraspecificEpithet",
                     "scientificNameAuthorship",  "scientificName", "acceptedNameUsage",
                     "taxon.rank", "taxon.status", "name.status", 
                     #"notes", "threat.status",
                     "life.form","habitat","vegetation.type","phytogeographicDomain",
                     "occurrence","establishment", "endemism","reference","vernacular.name", 
                     "species.correct2")]
saveRDS(taxon2, "data/threat_fbo_tax_info.rds")

################################################################################################################################################H
################################################################################################################################################H
#######################################
#### VOUCHERS FROM THE FLORA DWC-A ####
#######################################
## Finding all synonyms and basionyms
relat <- fbo$resourcerelationship.txt
relat <- relat[relat$id %in% taxon1$id, ]
relat.syn <- relat[!(duplicated(relat$relatedResourceID) | 
                       duplicated(relat$relatedResourceID, fromLast = TRUE)),]
toto <- cbind.data.frame(id = taxon1$id,
                         relatedResourceID = taxon1$id,
                         relationshipOfResource = "nome VALIDO")
relat.syn <- rbind.data.frame(relat.syn, toto)

#### VOUCHERS FROM REFLORA ####
vouchers <- fbo$typesandspecimen.txt
vouchers <- vouchers[vouchers$id %in% relat.syn$relatedResourceID, ]
vouchers <- vouchers[order(vouchers$id), ]
vouchers.syn <- merge(vouchers, relat.syn, by.x = "id",
                      by.y = "relatedResourceID", all.x = TRUE, sort = FALSE,
                      suffixes = c(".x",".valid"))
vouchers.syn <- vouchers.syn[order(vouchers.syn$id), ]
vouchers.syn1 <- merge(vouchers.syn, taxon1[,c("id","species.correct2")],
                       by.x = "id.valid", by.y = "id", all.x = TRUE, sort = FALSE)

## Isolating and editing voucher information from Reflora
vouchers.syn1$recordedBy1 <-
  gsub(" [0-9].*| s\\.n\\..*","", vouchers.syn1$recordedBy)
vouchers.syn1$colNumber <-
  gsub(".*[^0-9]","", vouchers.syn1$recordedBy)
vouchers.syn1$colNumber[grepl(" s\\.n\\.$", vouchers.syn1$recordedBy)] <- 
  "s.n."
vouchers.syn1$lastName <- plantR::lastName(vouchers.syn1$recordedBy1)
vouchers.syn1$recordBy.init <- plantR::lastName(vouchers.syn1$recordedBy1, 
                                                invert = TRUE, initials = TRUE)

vouchers.syn1$numTombo <- paste(tolower(vouchers.syn1$collectionCode),
                                as.double(vouchers.syn1$catalogNumber), sep="_")

##Tipos
vouchers.syn1$typus <- FALSE
vouchers.syn1$typus[vouchers.syn1$typeStatus %in% "typus"] <- TRUE

##Tombo
vouchers.syn1$numTombo <- paste(tolower(vouchers.syn1$collectionCode),
                                as.double(vouchers.syn1$catalogNumber),
                                sep = "_")
vouchers.syn1$numTombo[grepl("_NA$", vouchers.syn1$numTombo)] <- NA 
vouchers.syn1$voucher <- NA

##Merging with previous data and saving
old.vouchers <- readRDS("data/vouchers_reflora_old.rds")
old.vouchers$tax <- gsub("[0-9]$", "", old.vouchers$tax)
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$ï..status %in% c("replace"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  old.vouchers$tax[old.vouchers$tax %in% sp.i] <- rpl.i
}
new.vouchers <- vouchers.syn1[, c("species.correct2","voucher","recordedBy",
                                  "typus", "lastName","recordBy.init",
                                  "colNumber","numTombo")]
names(new.vouchers) <- names(old.vouchers)
all.vouchers <- rbind.data.frame(old.vouchers, new.vouchers)
saveRDS(all.vouchers, "data/vouchers_reflora.rds")


################################################################################################################################################H
################################################################################################################################################H
##################################H
#### INFO FROM REFLORA WEBSITE ####
##################################H
# DONE ONCE SEMI-AUTOMATICALLY AND SAVED IN THE DATA FOLDER (CODES NEED TO BE UPDATED)
# source("C:/Users/renato/Documents/raflima/R_packages/plantR/R/getFlora.R")
# tax  <- readRDS("data/assess_iucn_spp.rds")[,1:3]
# flora.info <- getFlora(tax$species.correct2)
# tmp <- cbind.data.frame(tax, flora.info)

# ## Editing the free descriptions
# ids <- !is.na(tmp$descricaoLivrePT)
# #spliting the strings
# lista <- strsplit(tmp$descricaoLivrePT[ids], '\\"\\>|\\</span\\>|\\<span|\\</font\\>', perl = TRUE)
# #finding the text descriptions
# patt <- "rbusto|rvore|indumento|folha|pares|flor|fruto|semente|glabr|altura|ramos|maduro|nervura|axilar"
# no.patt <- "-TRAD"
# descr <- sapply(lista, function(x) grepl(patt, x, perl = TRUE) & 
#                   !grepl(no.patt, x, perl = TRUE))
# descr <- lapply(descr, function(x) {
#   y <- length(x)
#   x[y] <- TRUE
#   return(x)
# })
# lista.limpa <- mapply(function(x, y) { x[y] }, lista, descr)
# length.lista <- lengths(lista.limpa) > 1
# lista.limpa[length.lista] <- sapply(lista.limpa[length.lista], 
#                                     paste, collapse = " ")
# 
# #editing the text descriptions
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\<b\\>|\\</b\\>', "", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\<p\\>|\\</p\\>', "", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\<i\\>|\\</i\\>', "", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\<div\\>|\\</div\\>', "", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\<o\\:p\\>|\\</o\\:p\\>', "", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\<br\\>|\\</br\\>', "", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, function(x) gsub('\\\n', " ", x, perl = TRUE))
# lista.limpa <- sapply(lista.limpa, stringr::str_squish)
# check <- sapply(lista.limpa, function(x) x %in% "" | nchar(x) < 5)
# lista.limpa[check] <- NA
# #saving back to the file
# tmp$descricaoLivrePT[ids] <- as.character(unlist(lista.limpa))
# 
# # ## Saving
# # tmp$collectors[tmp$species.correct2 %in% "Miconia flammea"] <-
# #   gsub("Silva, M.F.O. \\| Prieto, P.V.", "Silva, M.F.O.; Prieto, P.V.", tmp$collectors[tmp$species.correct2 %in% "Miconia flammea"])
# 
# # saveRDS(tmp, "data/threat_flora_info.rds")
# tmp <- readRDS("data/threat_flora_info.rds")
# 
# ## Isolating and editing voucher information from Reflora
# tmp1 <- strsplit(tmp$tombos, "\\|")
# names(tmp1) <- tax$species.correct2
# tmp1 <- unlist(tmp1)
# 
# tmp2 <- strsplit(tmp$collectors, "\\|")
# tmp2 <- sapply(tmp2, stringr::str_trim)
# names(tmp2) <- tax$species.correct2
# tmp2 <- unlist(tmp2)
# 
# vouchers <- data.frame(tax = names(tmp1), 
#                        voucher = tmp1,
#                        collect = tmp2,
#                        stringsAsFactors = FALSE, row.names = NULL)
# vouchers$tax <- gsub("[0-9]$","", vouchers$tax)
# ##Tombo
# vouchers$voucher <- sapply(strsplit(vouchers$voucher, ","), function(x) x[1])
# ##Tipos
# vouchers$typus <- FALSE
# vouchers$typus[grepl('><b>Typus</b></font',vouchers$collect)] <- TRUE
# vouchers$collect <- gsub(", <font color='red'><b>Typus</b></font>", "", vouchers$collect)
# ##Coletor
# vouchers$recordBy <- stringr::str_trim(sapply(strsplit(vouchers$collect, ","), function(x) x[1]))
# vouchers$recordBy <- gsub("   ", " ", vouchers$recordBy)
# vouchers$recordBy <- gsub("  ", " ", vouchers$recordBy)
# tmp3 <- stringr::str_trim(sapply(strsplit(vouchers$collect, ","), function(x) x[2]))
# tmp4 <- stringr::str_trim(sapply(strsplit(vouchers$collect, ","), function(x) x[3]))
# ##Iniciais do Coletors
# vouchers$recordBy.init <- NA
# vouchers$recordBy.init[grepl("\\D", tmp3)] <-
#   tmp3[grepl('\\D', tmp3)]
# vouchers$recordBy.init[is.na(vouchers$recordBy.init) & grepl("\\D", tmp4)] <-
#   tmp4[is.na(vouchers$recordBy.init) & grepl('\\D', tmp4)]
# vouchers$recordBy.init[vouchers$recordBy.init %in% "s.n."] <- NA
# ##Número de coleta
# vouchers$colNumber <- NA
# vouchers$colNumber[grepl("\\d", tmp3) | tmp3 %in% "s.n."] <-
#   tmp3[grepl('\\d', tmp3) | tmp3 %in% "s.n."]
# vouchers$colNumber[is.na(vouchers$colNumber) & (grepl("\\d", tmp4) | tmp4 %in% "s.n.")] <-
#   tmp4[is.na(vouchers$colNumber) & (grepl('\\d', tmp4) | tmp4 %in% "s.n.")]
# 
# ##Final edits
# vouchers$recordBy[grepl(" & ", vouchers$recordBy)] <- 
#   sapply(strsplit(vouchers$recordBy[grepl(" & ", vouchers$recordBy)], " & "), 
#          function(x) x[1])
# vouchers$recordBy <- gsub(" \\(Brother\\) ", " ",vouchers$recordBy)
# vouchers$recordBy[grepl(" \\(", vouchers$recordBy)] <- 
#   sapply(strsplit(vouchers$recordBy[grepl(" \\(", vouchers$recordBy)], " \\("), 
#          function(x) x[1])
# vouchers$voucher[grepl("[A-Z] [0-9][0-9]", vouchers$recordBy.init) &
#                    vouchers$voucher %in% "NA"] <- 
#   vouchers$recordBy.init[grepl("[A-Z] [0-9][0-9]", vouchers$recordBy.init) &
#                            vouchers$voucher %in% "NA"]
# vouchers$recordBy.init[vouchers$recordBy.init %in% vouchers$colNumber] <- NA
# ids <- is.na(vouchers$recordBy.init) & grepl("\\. ", vouchers$recordBy)
# vouchers$recordBy.init[ids] <- 
#   sapply(strsplit(vouchers$recordBy[ids], "\\. "), function(x) paste0(paste(head(x, -1), collapse="."),".")) 
# vouchers$recordBy[ids] <- 
#   sapply(strsplit(vouchers$recordBy[ids], "\\. "), function(x) tail(x,1)) 
# vouchers$recordBy <- gsub(" et al\\.$", "", vouchers$recordBy)
# 
# ##PlantR columns
# #Tombo
# vouchers$numTombo <- sapply(strsplit(vouchers$voucher, " "), 
#                             function(x) paste0(tolower(x[1]),"_",x[2]))
# vouchers$numTombo[grepl("_NA$", vouchers$numTombo)] <- NA 
# #Last name
# 
# ##Saving
# saveRDS(vouchers, "data/vouchers_reflora_old.rds")
rm(list = ls())

