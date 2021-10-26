##################################################################H
##################################################################H
#### PREPARING SPECIES INFORMATION FOR THREAT IUCN ASSESSMENTS ####
##################################################################H
##################################################################H
rm(list=ls())
require(taxize)
require(flora)
require(redlistr)

######################H
#### FULL TAXONOMY ####
######################H
## Red List IUCN key: 814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25

## Reading the THREAT species list ##
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]
# tax  <- readRDS("data/assess_iucn_spp.rds")[,1:3]

#### NCBI ####
## Seting the ncbi key
#API-key: c9366b59140a833ccd26fe9ac4f92e48d208 
options(ENTREZ_KEY = "c9366b59140a833ccd26fe9ac4f92e48d208")

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
saveRDS(ncbi.taxonomy, "data/threat_ncbi_taxonomy.rds")

#lower taxonomy
# gens <- sort(unique(ncbi.taxonomy$genus))
# # res.gens <- vector("list", length(gens))
# # for(i in 1:length(gens)){
# #   cat(i,"\n")
# #   res.gens[[i]] <- try(children(gens[i], db = "ncbi"), TRUE)
# # }
# tmp  <- lapply(gens, function(x) try(taxize::children(x, "ncbi"), TRUE))
# tmp1 <- do.call(rbind.data.frame, tmp1)

#### TROPICOS ####
## Getting Tropicos full taxonomy ##
#Tropicos API key: E8E538BC-0DAD-47EF-84CA-F6A663D9170A
options(TROPICOS_KEY = "E8E538BC-0DAD-47EF-84CA-F6A663D9170A")

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
tps.taxonomy1 <- merge(tps.taxonomy, output3, by.x = "species.correct2", 
                       by.y = "search.string", all.x = TRUE, sort= FALSE)
saveRDS(tps.taxonomy1, "data/threat_tropicos_taxonomy.rds")

#### IUCN ####
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
    out <- c(iucnid= iucn[[x]][1], atributos[c(2:4)]) 
  }
)
iucn2 <- do.call(rbind.data.frame, iucn1)
names(iucn2) <- c("Redlist_id","match","name","uri")
for(i in 1:3) iucn2[,i] <- as.character(iucn2[,i])
saveRDS(iucn2, "data/threat_iucn_taxonomy.rds")

#### BRAZILIAN FLORA INFORMATION ####
## Getting FB 2020 taxonomy
fbo.info <- readRDS("data/threat_fbo_tax_info.rds")
toto <- dplyr::left_join(tax, fbo.info, by = "species.correct2")
# toto <- flora::get.taxa(tax$species.correct2, replace.synonyms = FALSE, suggest.names = FALSE,
#                         life.form = TRUE, habitat = TRUE, vegetation.type = TRUE,
#                         vernacular = TRUE, states = TRUE, establishment = TRUE,
#                         drop = c(""),
#                         suggestion.distance = 0.9, parse = FALSE)

toto <- toto[c("family","genus","species.correct2","scientificNameAuthorship",
               "id","acceptedNameUsage","taxon.rank","taxon.status",
               "name.status", #"notes", "threat.status",
               "life.form","habitat","vegetation.type","phytogeographicDomain",
               "occurrence","establishment", "reference", "vernacular.name",
               "internal_taxon_id")]


#### PUTTING ALL TAXONOMIC INFORMATION TOGETHER AND SAVING ####
ncbi.taxonomy <- readRDS("data/threat_ncbi_taxonomy.rds")
tps.taxonomy1 <- readRDS("data/threat_tropicos_taxonomy.rds")
iucn2 <- readRDS("data/threat_iucn_taxonomy.rds")

# Merging flora and IUCN species codes
toto1 <- merge(toto, iucn2[,!names(iucn2) %in% "match"], 
               by.x = "species.correct2", by.y = "name", all.x = TRUE, sort= FALSE)
full.tax <- merge(toto1, tps.taxonomy1[,!names(tps.taxonomy1) %in% "genus"], 
                  by = "species.correct2", all.x = TRUE, sort= FALSE)
full.tax$source <- "Flora do Brasil 2020 (http://floradobrasil.jbrj.gov.br)"
full.tax$source[is.na(full.tax$scientificNameAuthorship)] <- "Tropicos (www.tropicos.org)" 
full.tax$id[is.na(full.tax$scientificNameAuthorship)] <- 
  full.tax$nameid[is.na(full.tax$scientificNameAuthorship)]
full.tax$scientificNameAuthorship[is.na(full.tax$scientificNameAuthorship)] <- 
  full.tax$author[is.na(full.tax$scientificNameAuthorship)]
full.tax$source[is.na(full.tax$scientificNameAuthorship)] <- "GBIF (www.gbif.org)" 

#Final search using other databases (GBIF works better)
splist <- full.tax$species.correct2[is.na(full.tax$scientificNameAuthorship)]
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
miss.tax1$scientificNameAuthorship <- stringr::str_trim(sapply(strsplit(miss.tax1$scientificname," "), function(x) paste(tail(x,-2), collapse = " ")))
cols <- c("family","genus","scientificNameAuthorship","usagekey","status")
full.tax[full.tax$species.correct2 %in% splist, 
         c("family","genus","scientificNameAuthorship","id","taxon.status")] <- miss.tax1[,cols]
full.tax$family[is.na(full.tax$family)] <- 
  full.tax$family.correct1[is.na(full.tax$family)]
full.tax$genus[is.na(full.tax$genus)] <- 
  sapply(strsplit(full.tax$species.correct2[is.na(full.tax$genus)], " "), function(x) x[1])

##Taxonomic Reference
full.tax$TaxonomicReference1 <- ## Brazilian Flora
  paste0("Flora do Brasil 2020 em construção. Jardim Botânico do Rio de Janeiro. Available at: http://servicos.jbrj.gov.br/flora/search/",
         gsub(" ", "_", full.tax$species.correct2),".")
full.tax$TaxonomicReference1[full.tax$source %in% "Tropicos (www.tropicos.org)"] <- NA_character_

full.tax$TaxonomicReference2 <- # Tropicos
  paste0("Tropicos.org. Missouri Botanical Garden. Available at: https://www.tropicos.org/name/", full.tax$nameid,".")
full.tax$TaxonomicReference2[is.na(full.tax$nameid)] <- NA_character_

full.tax$TaxonomicReference3 <- NA_character_ # GBIF
full.tax$TaxonomicReference3[full.tax$source %in% "GBIF (www.gbif.org)"] <-
  paste0("GBIF Secretariat. GBIF Backbone Taxonomy. Available at: https://www.gbif.org/species/", full.tax$id[full.tax$source %in% "GBIF (www.gbif.org)"],".")

## Saving full tax info for THREAT
saveRDS(full.tax, "data/threat_full.taxonomy.rds")

## Filtering and creating the final columns
final.cols <- c("kingdom","phylum","classname","ordername","family","genus","species.correct2","scientificNameAuthorship")
full.tax1 <- full.tax[ ,final.cols]
#,"accepted.name","taxon.status","name.status","notes","nomenclaturestatusname","displayreference","displaydate"
#"life.form","habitat","vegetation.type","vernacular.name","occurrence"

## Taxonomic Notes
full.tax1$TaxonomicNotes.value <- NA

# Name publication
ids <- !is.na(full.tax$reference) 
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0("De acordo com a Flora do Brasil 2020 (http://floradobrasil.jbrj.gov.br), esta espécie foi descrita em: ",
            full.tax$reference[ids])
ids <- is.na(full.tax$reference) 
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0("De acordo com Tropicos (www.tropicos.org), esta espécie foi descrita em: ",
         paste0(full.tax$displayreference[ids],", ", full.tax$displaydate[ids], "."))
table(is.na(full.tax1$TaxonomicNotes.value)) # ok! all names have a citation (all false)

# Adding other taxonomic notes
ids <- full.tax$taxon.status %in% c("accepted","ACCEPTED") & (grepl("correct", full.tax$name.status) | full.tax$name.status %in% c("", NA))
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0(full.tax1$TaxonomicNotes.value[ids], " O nome desta espécie é considerado como correto e aceito por: ", full.tax$source[ids], ".")
ids <- is.na(full.tax$taxon.status) & (grepl("correct", full.tax$name.status) | full.tax$name.status %in% c("", NA))
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0(full.tax1$TaxonomicNotes.value[ids], " O nome desta espécie é considerado como correto por: ", full.tax$source[ids], ".")
ids <- full.tax$taxon.status %in% "" & grepl("published", full.tax$name.status)
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0(full.tax1$TaxonomicNotes.value[ids], " O nome desta espécie é considerado como não efetivamente ou validamente publicado por: ", full.tax$source[ids], ".")
ids <- full.tax$taxon.status %in% "" & grepl("uncertain application", full.tax$name.status)
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0(full.tax1$TaxonomicNotes.value[ids], " O nome desta espécie é considerado como de aplicação incerta ou incertae saedis por: ", full.tax$source[ids], ".")
ids <- full.tax$taxon.status %in% "" & grepl("orthographic variant", full.tax$name.status)
full.tax1$TaxonomicNotes.value[ids] <-  paste0(full.tax1$TaxonomicNotes.value[ids], 
         " O nome desta espécie é considerado como uma variante ortográfica de ", 
         full.tax$acceptedNameUsage[ids],
         " por: ", full.tax$source[ids], ".")
ids <- full.tax$taxon.status %in% "" & grepl("legitimate, but incorrect", full.tax$name.status)
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0(full.tax1$TaxonomicNotes.value[ids], 
         " O nome desta espécie é considerado como legítimo porém incorreto por: ", full.tax$source[ids], ".")
ids <- full.tax$species.correct2 %in% "Ocotea rariflora"
full.tax1$TaxonomicNotes.value[ids] <- 
  "This taxon was taken as a synonym of Ocotea daphnifolia by Rohwer (1986), but Baitello (2003) [Fl. Fanerog. Estado São Paulo 3: 203] disagrees. The name is maked as unresolved in the Brazilian Flora 2020 (http://floradobrasil.jbrj.gov.br)."
# ids <- is.na(full.tax1$TaxonomicNotes.value) & !is.na(full.tax$accepted.name)
# full.tax1$TaxonomicNotes.value[ids] <- paste0("This taxon is taken as being correct by the ", full.tax$source[ids], 
#                                               ", but under a different accepeted name: ", full.tax$accepted.name[ids], ".")

# Adding species descriptions
spp.info <- readRDS("data/threat_flora_info.rds")
spp.info1 <- dplyr::left_join(full.tax1, spp.info, by = "species.correct2")
ids <- !is.na(spp.info1$descricaoLivrePT)
full.tax1$TaxonomicNotes.value[ids] <- 
  paste0(full.tax1$TaxonomicNotes.value[ids], 
         "\nA Flora do Brasil 2020 (http://floradobrasil.jbrj.gov.br) fornece a seguinte descrição para essa espécie: ", spp.info1$descricaoLivrePT[ids], ".")
full.tax1$TaxonomicNotes.value[ids] <- 
  gsub("\\.\\.$", ".", full.tax1$TaxonomicNotes.value[ids])

## Final edits and saving
names(full.tax1)[names(full.tax1) %in% "species.correct2"] <- "species"
names(full.tax1)[names(full.tax1) %in% "scientificNameAuthorship"] <- "taxonomicAuthority"
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

## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/taxonomy.csv")
head(sample, 3)
apply(sample, 2, unique)
names(sample)

## Getting the missing columns
full.tax1$infraType <- NA
full.tax1$infra_authority <- NA
full.tax1$infra_name <- NA
full.tax2 <- merge(full.tax1, tax[,c("internal_taxon_id", "species.correct2")],
                   by.x = "species", by.y = "species.correct2", all.x = TRUE, sort = FALSE)
full.tax3 <- merge(full.tax2, iucn2[,c("Redlist_id", "name")],
                   by.x = "species", by.y = "name", all.x = TRUE, sort = FALSE)
full.tax3 <- full.tax3[, match(names(sample), names(full.tax3), nomatch = 0)]

##Editing the species column
full.tax3$species <- 
  sapply(strsplit(full.tax3$species, " "), function(x) x[2])

## Saving
full.tax3 <- full.tax3[order(as.double(gsub("sp", "", full.tax3$internal_taxon_id))), ] 
write.csv(full.tax3, "data/sis_connect/taxonomy_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


#####################H
#### COMMON NAMES ####
#####################H

#Filtering and editing
common <- merge(full.tax[ ,c("species.correct2","vernacular.name")], tax[,c("internal_taxon_id", "species.correct2")],
                   by = "species.correct2", all.x = TRUE, sort = FALSE)
table(common$internal_taxon_id == full.tax2$internal_taxon_id)
names(common) <- c("validName", "name0", "internal_taxon_id")
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
common1$name0 <- gsub("pao \\/ pau \\/ pao", "Pau", common1$name0, ignore.case = TRUE)
common1$name0 <- gsub("pao \\/ pau", "Pau", common1$name0, ignore.case = TRUE)
common1$name0 <- gsub("cumixá/cumichá", "Cumixá", common1$name0, ignore.case = TRUE)
common1$name0 <- gsub("Caixeta/Caixeto", "Caixeta", common1$name0, ignore.case = TRUE)

#Breaking the info
tutu <- strsplit(common1$name0, "\\|")
tutu <- lapply(tutu, stringr::str_squish)
tutu <- lapply(tutu, function(x) t(sapply(strsplit(x, '\\/'), function(y) y[1:2])))
names(tutu) <- common1$validName
tutu1 <- do.call('rbind.data.frame', tutu)
tutu1 <- cbind.data.frame(validName = sapply(strsplit(rownames(tutu1), "\\."), function(x) x[1]),
                          language = as.character(tutu1[,2]),
                          name = stringr::str_to_title(as.character(tutu1[,1])),
                          stringsAsFactors = FALSE)
tutu1$language <- stringr::str_replace_all(tutu1$language,
                                           c("PORTUGUES" = "Portuguese",
                                             "ESPANHOL" = "Spanish; Castilian",
                                             "INGLES" = "English", 
                                             "KAXINAWA" = "South American Indian (Other)"))


## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/commonnames.csv")
apply(sample, 2, unique)
names(sample)
head(sample, 3)

## Saving
commonnames <- merge(tutu1, common1[,1:3], by = "validName", all.x = TRUE)
commonnames <- commonnames[order(commonnames$validName), ]
commonnames <- commonnames[,match(names(sample), names(commonnames), nomatch = 0)]
commonnames$primary <- "false"
write.csv(commonnames, "data/sis_connect/commonnames_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


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

## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/synonyms.csv")
apply(sample, 2, unique)
names(sample)
head(sample, 3)

#Getting the final columns in the same order and format as the sample
synonyms1 <- merge(synonyms, tax[,c("internal_taxon_id", "species.correct2")],
                   by.x = "validName", by.y = "species.correct2", all.x = TRUE)
synonyms1 <- synonyms1[order(synonyms1$validName),]
synonyms2 <- synonyms1[,match(names(sample), names(synonyms1), nomatch = 0)]
synonyms3 <- synonyms2[!is.na(synonyms2$name), ]
ids <- synonyms3$infraType %in% c(NA, "")
synonyms3$name[ids] <- paste(synonyms3$name[ids], synonyms3$speciesAuthor[ids])
ids <- synonyms3$infraType %in% c("var.", "subsp.", "fo.")
synonyms3$name[ids] <- paste(synonyms3$name[ids], synonyms3$infrarankAuthor[ids])
head(sample, 3)
head(synonyms3, 3)

## Saving
write.csv(synonyms3, "data/sis_connect/synonyms_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
#############################H
#### HABITATS AND ECOLOGY ####
#############################H
rm(list = ls())

## Getting Flora do BRasil information ##
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]
# tax  <- readRDS("data/assess_iucn_spp.rds")[,1:3]

fbo.info <- readRDS("data/threat_fbo_tax_info.rds")
toto <- dplyr::left_join(tax, fbo.info, by = "species.correct2")
hab <- toto[, c(names(tax), "life.form","habitat","vegetation.type","establishment")]
# toto <- flora::get.taxa(tax$species.correct2, replace.synonyms = FALSE, suggest.names = FALSE,
#                         life.form = TRUE, habitat = TRUE, vegetation.type = TRUE,
#                         vernacular = TRUE, states = TRUE, establishment = TRUE,
#                         drop = c(""),
#                         suggestion.distance = 0.9, parse = FALSE)
# hab <- cbind.data.frame(tax, 
#                         toto[,c("life.form","habitat","vegetation.type","establishment")],
#                         stringsAsFactors = FALSE)

#### OBTAINING GOOD GUESSES OF GENERATION LENGTH BASED ON SPECIES ECOLOGICAL GROUPS AND MAXIMUM SIZE ###
path = "C:/Users/renato/Documents/raflima/Pos Doc/Databases/TreeCo Database Management"
source(paste(path,"trait_data_prep.R",sep="/"))

## List of species
hab$genus <- sapply(strsplit(hab$species.correct2," "), function(x) x[1])
hab$taxon.rank <- "species"
hab <- hab[,c("internal_taxon_id","species.correct2","family.correct1","genus","species.correct2","taxon.rank","life.form","habitat","vegetation.type")]
names(hab) <- c("internal_taxon_id", "Name_submitted","family","genus","species.correct","taxon.rank","life.form.reflora","habitat.reflora","vegetation.type.reflora")
hab$ordem.spp = 1:dim(hab)[1]

## South American tree species list
trees.sa = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv",as.is=TRUE,na.string=c(""," ",NA))
## Reading trait data per taxonomic level ##
trait.spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species traits//Atributos das espécies//traits.species.csv",as.is=T)
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$ï..status %in% c("replace"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$ï..status[i]
  
  if (st.i == "replace") {
    trait.spp$TAX[trait.spp$TAX %in% sp.i] <- rpl.i
  }
}
trait.spp <- trait.spp[!duplicated(trait.spp$TAX),]
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
table(hab$Name_submitted == aux$species.correct)

# Establishment in respect to the AF
hab$establishment <- traits1$establishment.AF
hab$establishment[is.na(hab$establishment)] <- traits1$establishment.BR[is.na(hab$establishment)]
#hab$establishment[is.na(hab$establishment)] <- "check"
hab$establishment[hab$species.correct %in% c("Ormosia monosperma")] <-
  "native"
hab$establishment[hab$establishment %in% "native" & 
                    grepl("not in B",aux$establishment.BR) & !is.na(aux$Distribution.VCS)] <- 
  "AF native but not in Brazil"

## Flagging species according to their habit
hab$habito <- traits1$habito
hab$habito[hab$habito %in% c("1","1?")] = "tree"  
hab$habito[hab$habito %in% c("0.5","0.5?")] = "shrub"  
hab$habito[hab$habito %in% c("0")] = "non_tree"  
table(hab$habito, useNA = "always")

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
hab$life.form[hab$life.form %in% "woody_tree" & hab$life.form.reflora %in% c("Arbusto|Subarbusto", "Subarbusto")] <- "Woody tree (Shrub or treelet)"
hab[hab$species.correct %in% "Mimosa setosa",c("habito","life.form")] <- c("shrub", "Woody tree (Shrub or treelet)")
table(hab$life.form, useNA = "always")

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
table(is.na(hab$MaxHeight), useNA = "always")
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
hab$GF[is.na(hab$GF) & hab$life.form  %in% c("Woody tree (Shrub or treelet)", "Palm or palmoid", "Woody tree (Lianescent)") & 
         hab$habito %in% "shrub"] <- 1
hab$GF[is.na(hab$GF) & hab$life.form  %in% "Tree fern" & 
         hab$habito %in% "tree"] <- 2
hab$GF[is.na(hab$GF) & (hab$life.form  %in% "Woody tree (Tree)" | 
                          grepl("Árvore", hab$life.form.reflora))] <- 4
hab$GF <- stringr::str_replace_all(hab$GF, c("1" = "large_shrub", "2" = "small_tree", "3" = "large_tree", "4" = "tree_unknown"))
# hab$GF[is.na(hab$GF)] <- "unknown"
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
table(hab$ecol.group, useNA = "always")

## Including other traits that may be important
table(hab$Name_submitted == traits1$Name_submitted)
hab$SeedMass_g <- traits1$SeedMass_g
hab$dispersal.syndrome <- traits1$dispersal.syndrome

## Saving preliminar habitat data
write.csv(hab, "data/threat_habitats_preliminar.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


#################################################H
#### GENERATION LENGHT AND MATURE INDIVIDUALS ####
#################################################H

#### SETTING THE PROXIES OF GENERATION LENGTH FOR EACH SPECIES ###
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
               #6, 7, 8, 9, 7.5, # for small trees
               7, 8, 9, 10, 8, # for small trees
               10, 12.5, 15, 20, 12.5, # for large trees
               8, 10, 12.5, 15, 10) # for trees unknown
#Prop. mature
combo$p.est <- c(1, 1, 1, 1, 1, # for shrubs
                 #0.7213, 0.5998, 0.4927, 0.2868, 0.5122, # for small trees (old dbh.crit: 6, 7, 8, 9, 7.5cm)
                 0.6125, 0.4661, 0.4143, 0.2385, 0.4669, # for small trees (new dbh.crit: 7, 8, 9, 10, 8cm)
                 0.5931, 0.4485, 0.3028, 0.2588, 0.4160, # for large trees
                 0.6414, 0.5089, 0.3346, 0.2470, 0.4540) # for trees unknown
combo$p.ci <- c("0.99-1.1", "0.975-1.075", "0.95-1.05", "0.90-1.025", "0.95-1.05", # for shrubs
                "0.4731-0.7518", "0.3827-0.5495", "0.3494-0.4792", "0.1356-0.3415", "0.4285-0.5054", # for small trees
                "0.5224-0.6637", "0.4152-0.4817", "0.2736-0.3321", "0.1874-0.3332", "0.3966-0.4353", # for large trees
                "0.5761-0.7066", "0.4778-0.5399", "0.3075-0.3617", "0.1875-0.3301", "0.4362-0.4718") # for trees unknown
#Merging the info
combo$str.match <- paste(combo$EG, combo$GF, sep = "_") 
hab$str.match <- paste(hab$ecol.group, hab$GF, sep = "_")
hab1 <- dplyr::left_join(hab, combo[,c("str.match","GL","DBH","p.est","p.ci")], by= "str.match") 

#Adding Proportion of mature individuals in the population (from code '')
p.est <- read.csv("data/prop_mature.csv", as.is = TRUE)
hab1 <- dplyr::left_join(hab1, p.est[,c("species.correct","N","dbh.crit","p")], by= "species.correct") 


#### PUTTING ALL FIELDS IN THE IUCN SIS FORMAT ###

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
hab1$vt <- gsub("^1\\.5$", "Subtropical/Tropical Dry Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1\\.5\\|", "Subtropical/Tropical Dry Forest|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.5$", "|Subtropical/Tropical Dry Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.5\\|", "|Subtropical/Tropical Dry Forest|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^1\\.6$", "Subtropical/Tropical Moist Lowland Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1\\.6\\|", "Subtropical/Tropical Moist Lowland Forest|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.6$", "|Subtropical/Tropical Moist Lowland Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.6\\|", "|Subtropical/Tropical Moist Lowland Forest|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^1\\.7$", "Subtropical/Tropical Mangrove Forest Vegetation Above High Tide Level", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1\\.7\\|", "Subtropical/Tropical Mangrove Forest Vegetation Above High Tide Level|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.7$", "|Subtropical/Tropical Mangrove Forest Vegetation Above High Tide Level", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.7\\|", "|Subtropical/Tropical Mangrove Forest Vegetation Above High Tide Level|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^1\\.8$", "Subtropical/Tropical Swamp Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^1\\.8\\|", "Subtropical/Tropical Swamp Forest|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.8$", "|Subtropical/Tropical Swamp Forest", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|1\\.8\\|", "|Subtropical/Tropical Swamp Forest|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^2\\.1$", "Dry Savanna", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^2\\.1\\|", "Dry Savanna|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|2\\.1$", "|Dry Savanna", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|2\\.1\\|", "|Dry Savanna|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^3\\.5$", "Subtropical/Tropical Dry Shrubland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^3\\.5\\|", "Subtropical/Tropical Dry Shrubland|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|3\\.5$", "|Subtropical/Tropical Dry Shrubland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|3\\.5\\|", "|Subtropical/Tropical Dry Shrubland|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^4\\.6$", "Subtropical/Tropical Seasonally Wet/Flooded Lowland Grassland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^4\\.6\\|", "Subtropical/Tropical Seasonally Wet/Flooded Lowland Grassland|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|4\\.6$", "|Subtropical/Tropical Seasonally Wet/Flooded Lowland Grassland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|4\\.6\\|", "|Subtropical/Tropical Seasonally Wet/Flooded Lowland Grassland|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^4\\.7$", "Subtropical/Tropical High Altitude Grassland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^4\\.7\\|", "Subtropical/Tropical High Altitude Grassland|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|4\\.7$", "|Subtropical/Tropical High Altitude Grassland", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|4\\.7\\|", "|Subtropical/Tropical High Altitude Grassland|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^6$", "Inland Rocky Areas", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^6\\|", "Inland Rocky Areas|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|6$", "|Inland Rocky Areas", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|6\\|", "|Inland Rocky Areas|", hab1$vt, perl = TRUE)

hab1$vt <- gsub("^14\\.5$", "Urban Areas", hab1$vt, perl = TRUE)
hab1$vt <- gsub("^14\\.5\\|", "Urban Areas|", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|14\\.5$", "|Urban Areas", hab1$vt, perl = TRUE)
hab1$vt <- gsub("\\|14\\.5\\|", "|Urban Areas|", hab1$vt, perl = TRUE)
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
hab1$System.value <- gsub("Terricola|Terrícola|Hemiparasita|Rupícola|Epífita|Hemiepífita|Hemiparasita", "Terrestrial", hab1$System.value)
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
hab1$p.ci[is.na(hab1$p.est)] <- "0.44-0.47"
hab1$p.est[is.na(hab1$p.est)] <- 0.454


## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/habitats.csv")
head(sample, 3)
apply(sample, 2, unique)
names(sample)


##Organizing and saving
hab2 <- hab1[order(hab1$species.correct),]
cols <- c("internal_taxon_id","Name_submitted",
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

## Saving the SIS connect file ##
lookup <- strsplit(hab2$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup, "\\|")
names(lookup) <- hab2$internal_taxon_id
lookup1 <- unlist(lookup)
taxa_id <- rep(names(lookup), lengths(lookup))
hab3 <- data.frame(
  GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup = lookup1,
  GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName = NA,
  GeneralHabitats.GeneralHabitatsSubfield.majorImportance = NA,
  GeneralHabitats.GeneralHabitatsSubfield.season = "resident",
  GeneralHabitats.GeneralHabitatsSubfield.suitability = "Suitable",
  internal_taxon_id = taxa_id, stringsAsFactors = FALSE
)
head(sample, 3)
head(hab3, 3)
write.csv(hab3, "data/sis_connect/habitats_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

#####################################################################################################################################################################H
#####################################################################################################################################################################H
#######################H
#### PLANT SPECIFIC ####
#######################H

## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/plantspecific.csv")
head(sample, 3)
apply(sample, 2, unique)
names(sample)

## Filtering and editing
cols <- c("internal_taxon_id",
          "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup",
          "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName")
hab4 <- hab2[, cols]
hab4 <- hab4[, match(names(sample), names(hab4), nomatch = 0)]
hab4[[1]][ hab4[[2]] %in% "Tree - large"] <- "TL"
hab4[[1]][ hab4[[2]] %in% "Tree - small"] <- "TS"
hab4[[1]][ hab4[[2]] %in% "Shrub - large"] <- "SL"
hab4[[1]][ hab4[[2]] %in% "Tree - size unknown"] <- "T"
hab4[[1]][ hab4[[2]] %in% "Succulent - tree"] <- "ST"
hab4[[1]][ hab4[[2]] %in% "Fern"] <- "PT"
head(sample)
head(hab4, 3)

## Saving
write.csv(hab4, "data/sis_connect/plantspecific_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

#####################################################################################################################################################################H
#####################################################################################################################################################################H
###########################H
#### COUNTRY OCCURRENCE ####
###########################H

## Generating the SIS CONNECT file "countries.csv"
## Sample from SIS CONNECT v6.1
head(read.csv("SIS_sample_6_1/countries.csv"))

## Getting info from reflora and formatting
taxon2 <- readRDS("data/threat_fbo_tax_info.rds")
taxon2 <- taxon2[,c("occurrence", "species.correct2")]
taxon2 <- taxon2[!is.na(taxon2$occurrence),]
br.states <- strsplit(gsub("BR-", "", taxon2$occurrence), "\\|")
br.states.trans <- data.table::transpose(br.states, fill = NA)
br.states.table <- br.states.trans
for(i in 1:length(br.states.trans)) {
  df <- cbind(species.correct2 = taxon2$species.correct2,
              country = "brazil",
              stateProvince = br.states.trans[[i]])
  br.states.table[[i]] <- df
}
taxon2.loc <- do.call(rbind, br.states.table) 
taxon2.loc <- as.data.frame(taxon2.loc[!is.na(taxon2.loc[,3]), ])
taxon2.loc <- plantR::fixLoc(taxon2.loc, loc.levels = c("country", "stateProvince"))
taxon2.loc$municipality.new <- taxon2.loc$locality.new <- NA_character_
taxon2.loc$loc.correct <- plantR::strLoc(taxon2.loc)$loc.string
states.fbo <- data.table::data.table(
  taxon2.loc[,c("loc.correct", "species.correct2")])
states.fbo[ , N := 1L]
data.table::setkeyv(states.fbo, "species.correct2")

## Loading species herbarium data
occs <- readRDS("data/threat_species_by_country.rds")

## Getting the country codes  (CountryOccurrenceLookup)
# ISO2 (countries)
countries <- occs[ , .N , by = c("loc.correct", "species.correct2")]
rm(occs)
countries$loc.correct <- gsub("_NA$", "", countries$loc.correct, perl = TRUE)
countries$loc.correct <- gsub("u\\.s\\. virgin islands", "united states virgin islands", 
                              countries$loc.correct, perl = TRUE)
countries$loc.correct <- gsub("cocos", "cocos islands", 
                              countries$loc.correct, perl = TRUE)
countries$loc.correct <- plantR::prepCountry(countries$loc.correct)
countries$loc.correct[is.na(countries$loc.correct)] <- "no_loc"

#merging occurrence and FBO information on species states of occurrence in BR
countries <- rbind(countries, states.fbo)
tmp <- paste(countries$loc.correct, countries$species.correct2, sep = "_")
countries <- countries[!duplicated(tmp),]
data.table::setkeyv(countries, "species.correct2")

#getting ADM names for the localities
tmp <- plantR::getAdmin(as.data.frame(countries))[,1:3]
paises.ADM1 <- c("Argentina", "Brazil", "Chile", "Mexico")
tmp$NAME_1[!tmp$NAME_0 %in% paises.ADM1] <- NA
tmp$iso2 <- countrycode::countrycode(tmp$NAME_0, "country.name", "iso2c")
tmp$code <- tmp$iso2

# SIS CONNECT codes (states/provinces)
lookup.codes <- read.csv("SIS_sample_6_1/Countryoccurrence.countryoccurrencelookup.csv", 
                         encoding = "UTF-8")
for (i in 1:length(paises.ADM1)) {
  pais.i <- countrycode::countrycode(paises.ADM1[i], "country.name", "iso2c")
  lookup.codes.i <- lookup.codes[lookup.codes$iso2c %in% pais.i, ]
  replace_these <- tmp$iso2 %in% pais.i
  tmp$code[replace_these] <- 
    lookup.codes.i$code[match(tmp$NAME_1[replace_these], lookup.codes.i$description)]
}
tmp$code[is.na(tmp$code)] <- tmp$iso2[is.na(tmp$code)]
# Building the file
countries[ , CountryOccurrenceLookup := tmp$code]
countries <- countries[!is.na(CountryOccurrenceLookup), ]
countries[ , CountryOccurrenceName := NA_character_]
countries[ , formerlyBred := NA_character_]
#Origin (originlookup: Native - 0, Reintroduced - 1, Introduced - 2, Prehistorically Introduced - 3, Vagrant - 4, Uncertain - 5)
countries[ , origin := "Native"]
#Presence (presencelookup: Extant - 0, PossiblyExtinct - 1, Extinct - 2, Presence Uncertain - 3)
countries[ , presence := "Extant"]
countries[ , seasonality := "Resident"]
countries[ , c("loc.correct", "N") := NULL]
head(countries)

# Obtaining the internal_taxon_id information
res <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
countries <- dplyr::left_join(countries, res)
countries[ , species.correct2 := NULL]

## Saving
names(countries)[1:6] <- paste0("CountryOccurrence.CountryOccurrenceSubfield.", 
                                names(countries)[1:6])
head(countries,3)
head(read.csv("SIS_sample_6_1/countries.csv"),3)
write.csv(countries, "data/sis_connect/countries_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
#############################H
#### PREVIOUS ASSESSMENTS ####
#############################H
# rm(list=ls())

## Red List IUCN key: 814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25
??rredlist

## Getting the list of species in the Atlantic Forest
prev.assess <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]
# prev.assess <- readRDS("data/assess_iucn_spp.rds")
# prev.assess <- prev.assess[,1:3]

## Getting the IUCN assessments for Brazil (CNCFlora) - National level
tmp = flora::get.taxa(prev.assess$species.correct2, replace.synonyms = FALSE, life.form = TRUE)
tmp$search.str = prev.assess$species.correct2
tmp1 <- merge(prev.assess, tmp, by.x = "species.correct2", by.y = "search.str", all.x = TRUE, sort=FALSE)
#table(prev.assess$species.correct2 == tmp1$species.correct2)
prev.assess$status.reflora <- tmp1$threat.status

## Getting the global IUCN assessments (IUCN) 
#Citation: IUCN 2020. IUCN Red List of Threatened Species. Version 2020-1 <www.iucnredlist.org>
# iucn <- read.csv("IUCN_2020_assessments.csv", as.is = TRUE, na.string = c(NA,""," "))
#Citation: IUCN 2021. IUCN Red List of Threatened Species. Version 2021-1 <www.iucnredlist.org>
iucn <- read.csv("IUCN_2021_v1_assessments.csv", as.is = TRUE, na.string = c(NA,""," "))

iucn$species.correct2 <- sapply(strsplit(iucn$scientificName," "), function(x) paste(x[1], x[2],sep=" "))
tmp <- merge(prev.assess, iucn,  by= "species.correct2", all.x = TRUE)
tmp <- tmp[,c("internal_taxon_id","species.correct2","assessmentId","internalTaxonId","redlistCategory","redlistCriteria","yearPublished","assessmentDate","criteriaVersion",
              "language","rationale","habitat","threats","population","populationTrend","range","useTrade","systems","conservationActions",
              "realm","yearLastSeen","possiblyExtinct","possiblyExtinctInTheWild","scopes")]
table(prev.assess$species.correct2 == tmp$species.correct2)
prev.assess <- cbind.data.frame(prev.assess, tmp,
                                stringsAsFactors = FALSE)
prev.assess <- prev.assess[order(prev.assess$species.correct2),]
saveRDS(prev.assess, "data/sis_connect/prev_assessments_threat.rds")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
##################################################################H
#### SPECIES ENDEMISM LEVELS IN RESPECT TO THE ATLANTIC FOREST ####
##################################################################H

## Species endemism levels from Lima et al. 2020
end <- read.csv("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/AppendixF_endemism_levels.csv", as.is=TRUE)
end <- end[,c(1,2,7,9,11)]

## Replacing the synonym 
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$ï..status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$ï..status[i]
  
  if (st.i == "replace") {
    end$species[end$species %in% sp.i] <- rpl.i
  }
}

##Getting mean endemism for new synonyms
tmp <- aggregate(cbind(endemism.level..validated.taxonomy.only., endemism.level..validated.and.probably.validated.taxonomy.) ~ 
                   family + species, FUN = mean, data=end)
tmp1 <- aggregate(endemism.accepted.by.the.BF.2020 ~ family + species, 
                  FUN = function(x) paste(unique(x), collapse = "|"), data=end)
end <- dplyr::left_join(tmp, tmp1)
fbo <- readRDS("data/threat_fbo_tax_info.rds")[,c("species.correct2","phytogeographicDomain","endemism")]
names(fbo)[1] <- "species"
fbo$endemism <- stringr::str_to_title(fbo$endemism)
fbo <- fbo[!is.na(duplicated(fbo$species)),]
end <- dplyr::left_join(end, fbo)
end$endemism.accepted.by.the.BF.2020[!is.na(end$phytogeographicDomain) & 
                                       !grepl("Mata Atlântica", end$phytogeographicDomain)] <- "Not in the AF"
end$endemism.accepted.by.the.BF.2020[!is.na(end$phytogeographicDomain) & 
                                       end$phytogeographicDomain %in% "Mata Atlântica" &
                                       end$endemism %in% "Endemic"] <- "Endemic"
end$endemism.accepted.by.the.BF.2020[!is.na(end$phytogeographicDomain) & 
                                       grepl("Mata Atlântica", end$phytogeographicDomain) &
                                       !end$phytogeographicDomain %in% "Mata Atlântica"] <- "Not endemic"
end$endemism.accepted.by.the.BF.2020[end$species %in% "Myrcia glomerata"] <- "Not endemic"


##Classification
end$endemic <- NA
end$endemic[end$endemism.level..validated.taxonomy.only. >= 85 & end$endemism.level..validated.and.probably.validated.taxonomy. >= 85 & 
              !end$endemism.accepted.by.the.BF.2020 %in% "Not in the AF"] <- "endemic"
end$endemic[is.na(end$endemic) & end$endemism.level..validated.taxonomy.only. >= 50 & end$endemism.level..validated.and.probably.validated.taxonomy. >= 50 &
              end$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "endemic"
end$endemic[end$endemism.level..validated.taxonomy.only. < 15 & end$endemism.level..validated.and.probably.validated.taxonomy. < 15 &
              !end$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "occasional"
end$endemic[is.na(end$endemic) & end$endemism.level..validated.taxonomy.only. >= 50 & 
              end$endemism.level..validated.and.probably.validated.taxonomy. >= 50] <- "widespread_common"
end$endemic[is.na(end$endemic) & end$endemism.level..validated.taxonomy.only. >= 15 & 
              end$endemism.level..validated.and.probably.validated.taxonomy. >= 15] <- "widespread_sparse"
end$endemic[is.na(end$endemic) & 
              end$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "widespread_common"
end$endemic[is.na(end$endemic) & 
              !end$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "widespread_sparse"
table(end$endemic, useNA = "always")
#Probable occurrences
end.prob <- read.csv("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/AppendixD_probable_occurrences.csv", as.is=TRUE)[,1:5]
#end.prob$species  <- sapply(end.prob$scientific.name, flora::remove.authors)
end.prob$species  <- sapply(end.prob$scientific.name, function(x) paste(strsplit(x, " ")[[1]][1:2], collapse = " "))
end.prob <- end.prob[end.prob$status %in% "no records found but cited in BF.2020",]
end.prob <- end.prob[!end.prob$species %in% end$species,]

tmp <- flora::get.taxa(end.prob$species)
tmp1 <- flora::get_domains(tmp)
tmp2 <- flora::get_endemism(tmp)
end.prob$domain <- tmp1$domain 
end.prob$end.br <- tmp2$endemism 
end.prob$endemic <- NA
end.prob$endemic[end.prob$domain %in% "Mata Atlântica" & 
                   end.prob$end.br %in%  "Endemica"] <- "endemic"
end.prob$endemic[is.na(end.prob$endemic) & grepl("Mata Atlântica", end.prob$domain)] <- 
  "not_endemic"
end.prob$endemic[is.na(end.prob$endemic) & !grepl("Mata Atlântica", end.prob$domain)] <- 
  "not in the AF"
table(end.prob$endemic, useNA = "always")

#Merging the two data.frames
end1 <- end[,c("family","species",
               "endemism.level..validated.and.probably.validated.taxonomy.",
               "endemism.level..validated.taxonomy.only.",
               "endemic")]
names(end1)[3:4] <- c("endemism.level.1","endemism.level.2") 

end.prob$endemism.level.1 <- NA
end.prob$endemism.level.2 <- NA
end.prob1 <- end.prob[,c("family","species",
                         "endemism.level.1",
                         "endemism.level.2",
                         "endemic")]
end.final <- rbind.data.frame(end1, end.prob1,
                              stringsAsFactors = FALSE)
end.final <- end.final[order(end.final$species),]

##Getting the endemism in Brazil (no CNCFlora 'jurisidction')
fbo.info <- readRDS("data/threat_fbo_tax_info.rds")
names(fbo.info)[which(names(fbo.info) == "species.correct2")] <- "species" 
tmp0 <- dplyr::left_join(end.final, fbo.info, by = "species")
tmp <- flora::get.taxa(end.final$species)
tmp1 <- flora::get_endemism(tmp)
tmp0$endemism[is.na(tmp0$endemism)] <- 
  tmp1$endemism[is.na(tmp0$endemism)]
table(tmp0$species == end.final$species)
end.final$endemic.BR <- tmp0$endemism 
end.final$endemic.BR[end.final$endemic.BR %in% "Endemica"] <- 
  "endemic"
end.final$endemic.BR[end.final$endemic.BR %in% "Não endemica"] <- 
  "not_endemic"
end.final$endemic.BR[end.final$endemic.BR %in% "not endemic"] <- 
  "not_endemic"

## Saving
saveRDS(end.final, "data/sis_connect/endemism_threat.rds")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
###################H
#### CITES LIST ####
###################H
# 
# ## Paralelizing to get all 166 CITES pages faster
# cl <- snow::makeSOCKcluster(6)
# doSNOW::registerDoSNOW(cl)
# `%d%` <- foreach::`%dopar%`
# x <- NULL
# cites0 <- foreach::foreach(
#   x = 1:166,
#   #.combine = 'c',
#   .options.snow = NULL
# ) %d% {
#   try(rcites::spp_taxonconcept(query_taxon = '', pages = x, 
#                                token = "r5BFEg69iKmfx8OCqm3tIgtt"), TRUE)
# }
# snow::stopCluster(cl)
# # for(i in 1:6) {
# #   ele.i <- c(161:166)[i]
# #   cites[[ele.i]] <- cites0[[i]]
# # }
# 
# ## Binding all pages together
# cites1 <- vector("list", 7)
# names(cites1) <- c("all_id", "general", "higher_taxa", "accepted_names", 
#                    "common_names", "synonyms", "cites_listings")
# for(i in 1:7) {
#   cites1.i <- lapply(cites, function(x) x[[i]])
#   non.empty <- sapply(cites1.i, function(x) dim(x)[1] >0 ) 
#   cites1[[i]] <- as.data.frame(dplyr::bind_rows(cites1.i[non.empty]))
#   cat(i, "\n")
# }
# saveRDS(cites1, "data/cites_data.rds") # DONE IN 07/05/2021
cites1 <- readRDS("data/threat_cites.rds")

# Getting the best names for homonyms
tax.cites <- cites1$all_id
dups <- (duplicated(tax.cites$full_name) |  
           duplicated(tax.cites$full_name, fromLast = TRUE))
dup.spp <- unique(tax.cites$full_name[dups])
tax.cites$remove <- TRUE
for (i in 1:length(dup.spp)) {
  sp.i <- dup.spp[i]
  check_ids <- tax.cites$full_name %in% sp.i
  taxon1.i <- tax.cites[check_ids,]
  tax.cites$remove[check_ids] <- 
    ifelse(taxon1.i$active %in% TRUE, TRUE, FALSE)
  check_ids1 <- tax.cites$full_name %in% sp.i & tax.cites$remove
  taxon2.i <- tax.cites[check_ids1,]
  if (length(unique(taxon2.i$name_status)) > 1) 
    tax.cites$remove[check_ids1] <-
        ifelse(taxon2.i$name_status %in% "A", TRUE, FALSE)
}
tax.cites <- tax.cites[tax.cites$remove,]


## Crossing with THREAT species
tax$genus <- gsub(" .*", "", tax$species.correct2)
tax <- as.data.frame(tax)
#Species level
ids_spp <- merge(tax, tax.cites, by.y = "full_name",
                 by.x = "species.correct2", all.x = TRUE, sort = FALSE)
ids_spp <- ids_spp[order(ids_spp$species.correct2),]

#Genus level
ids_gen <- merge(tax, tax.cites, by.y = "full_name",
                 by.x = "genus", all.x = TRUE, sort = FALSE)
ids_gen <- ids_gen[match(ids_spp$species.correct2, ids_gen$species.correct2),]
table(ids_spp$species.correct2 == ids_gen$species.correct2)
ids <- is.na(ids_spp$id) & !is.na(ids_gen$id)
cols <- c("id", "rank", "name_status", "updated_at", "active", "author_year")
ids_spp[ids,  cols] <- ids_gen[ids, cols]

#family level
ids_fam <- merge(tax, cites1$all_id, by.y = "full_name",
                 by.x = "family.correct1", all.x = TRUE, sort = FALSE)
ids_fam <- ids_fam[match(ids_spp$species.correct2, ids_fam$species.correct2),]
table(ids_spp$species.correct2 == ids_fam$species.correct2)
ids <- is.na(ids_spp$id) & !is.na(ids_fam$id)
ids_spp[ids, cols] <- ids_fam[ids, cols]
ids_spp <- ids_spp[ , -(which(names(ids_spp) == "remove"))]


## Getting valid CITES IDs for synonyms
syn.cites <- cites1$synonyms
syn.cites1 <- syn.cites[syn.cites$id...2 %in% 
                          ids_spp$id[ids_spp$name_status %in% "S"],]
ids_spp$id[match(syn.cites1$id...2, ids_spp$id)] <- 
  syn.cites1$id...1

## General table
gen.cites <- cites1$general
cites_spp <- dplyr::left_join(ids_spp, 
                              gen.cites[,c("id", "full_name", "cites_listing")],
                              by = "id")
## CITES annotations
list.cites <- cites1$cites_listings
all_cites_spp <- dplyr::left_join(cites_spp, list.cites, by = "id")
all_cites_spp$hash_annotation[all_cites_spp$hash_annotation %in% c("NA", NA)] <- 
  NA_character_
all_cites_spp$code_annotation <- 
  substr(all_cites_spp$hash_annotation, start = 1, stop = 3)
all_cites_spp$code_annotation <- gsub("A$|L$", "", all_cites_spp$code_annotation)
all_cites_spp$hash_annotation <- 
  gsub("^#10|^#15|^#4|^#5|^#6", "", all_cites_spp$hash_annotation)
all_cites_spp$hash_annotation <- 
  stringr::str_squish(gsub("\\\n|\\\r", " ", all_cites_spp$hash_annotation))

## EU legislation
# taxon_ids <- unique(all_cites_spp$id[!is.na(all_cites_spp$id)])
# # EU <- rcites::spp_eu_legislation(taxon_id = taxon_ids,  # website error....
# #                            token = "izCSojTRrM62oMgkkyAAmgtt")
# EU <- vector("list", length(taxon_ids))
# names(EU) <- taxon_ids
# for (i in 1:length(taxon_ids)) {
#   EU[[i]] <- try(rcites::spp_eu_legislation(taxon_id = taxon_ids[i], 
#                                  token = "izCSojTRrM62oMgkkyAAmgtt"), TRUE) 
# }
# rerun <- sapply(EU, class) %in% "try-error"
# for (i in 1:length(taxon_ids[rerun])) {
#   EU[rerun][[i]] <- try(rcites::spp_eu_legislation(taxon_id = taxon_ids[rerun][i], 
#                                             token = "izCSojTRrM62oMgkkyAAmgtt"), TRUE) 
# }
# 
# EU1 <- vector("list", 2)
# names(EU1) <- c("eu_listings", "eu_decisions")
# for(i in 1:2) {
#   EU.i <- lapply(EU, function(x) x[[i]])
#   non.empty <- sapply(EU.i, function(x) dim(x)[1] >0)
#   non.empty[lengths(EU.i) <= 1] <- FALSE
#   EU1[[i]] <- as.data.frame(dplyr::bind_rows(EU.i[unlist(non.empty)]))
# }
# # saveRDS(EU1, "data/EU_listings_data.rds") # DONE IN 10/05/2021
EU1 <- readRDS("data/EU_listings_data.rds")

# EU listings and decisions
list.EU <- EU1$eu_listings
names(list.EU)[4:5] <- paste0(names(list.EU)[4:5], ".EU")
all_cites.EU_spp <- merge(all_cites_spp, list.EU[, c("taxon_concept_id","annex.EU","change_type.EU")],
                       by.x = "id", by.y = "taxon_concept_id",
                       all.x = TRUE, sort = FALSE)
decis.EU <- EU1$eu_decisions
decis.EU <- decis.EU[, -(which(names(decis.EU) %in% 
                                 c("id","start_date","is_current","eu_decision_type.description","term.code","term.name")))]
decis.EU1 <- aggregate(. ~ taxon_concept_id, data = decis.EU, 
                   function(x) paste(x, collapse = "|")) 
all_cites.EU_spp1 <- merge(all_cites.EU_spp,  decis.EU1,
                          by.x = "id", by.y = "taxon_concept_id",
                          all.x = TRUE, sort = FALSE)
## Saving 
all_cites.EU_spp1 <- all_cites.EU_spp1[order(all_cites.EU_spp1$species.correct2),]
saveRDS(all_cites.EU_spp1, "data/threat_cites_EU.rds")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
######################H
#### USE AND TRADE ####
######################H

## Generating the SIS CONNECT file "usetrade.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/usetrade.csv")
head(sample, 3)

##Loading the species list
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]

##Loading and Preparing TreeCo use database
path = "C:/Users/renato/Documents/raflima/Pos Doc/Databases/TreeCo Database Management"
source(paste(path,"uses_data_prep.R",sep="/"))
usos <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species Uses//plant_uses.csv", 
                 na.strings = c(""," ",NA), encoding = "UTF-8")
usos1 <- uses_data_prep(usos, spp = tax$species.correct2)

##Crossing TreeCo and IUCN use codes and nomenclature
lookup.codes <- read.csv("SIS_sample_6_1/Utenduse.utenduselookup.csv", 
                        encoding = "UTF-8")
usos1$UTEndUseLookup <- NA
usos1$UTEndUseLookup <- 
  lookup.codes$code[match(usos1$uses, lookup.codes$treeco)]

##Creating the other required columns
usos1$UTEndUseName <- NA
usos1$international <- NA
usos1$national <- NA
usos1$other <- NA
usos1$subsistence <- NA

##Editing the other required columns
usos1$international[grepl("Agroforestree|Mark, J\\.", usos1$source)] <- "true"
usos1$international[!grepl("Agroforestree|Mark, J\\.", usos1$source)] <- "false"

usos1$national[grepl("Brazilian|Nacional|nativas", usos1$source)] <- "true"
usos1$national[!grepl("Brazilian|Nacional|nativas", usos1$source)] <- "false"

usos1$subsistence[usos1$UTEndUseLookup %in% c(1, 2, 3)] <- "true"
usos1$subsistence[!usos1$UTEndUseLookup %in% c(1, 2, 3)] <- "false"

##Adding the species internal id
usos1 <- merge(usos1, tax, by.x = "Name_submitted", by.y = "species.correct2", 
               all.x = TRUE, sort = FALSE)

##Adding missing sources
usos1$source[is.na(usos1$source)] <- 
  usos1$link[is.na(usos1$source)]

##Filtering the selected uses for the SIS CONNECT file
usos2 <- usos1[!is.na(usos1$UTEndUseLookup),]
usos2 <- usos2[ , c(names(usos2)[which(names(usos2) == "UTEndUseLookup"): dim(usos2)[2]], "source")]
names(usos2)[1:6] <- paste0("UTEndUse.UTEndUseSubfield.", 
                            names(usos2)[1:6])
head(usos2)

##Saving the SIS Connet file
write.csv(usos2, "data/sis_connect/usetrade_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


##Saving the file with the exploited species
usos3 <- usos1[usos1$uses %in% c("timber") |
                 usos1$Name_submitted %in% c("Euterpe edulis") , 
               names(usos1)[1:which(names(usos1) == "uses")]]

check_these <- !usos1$Name_submitted %in% usos3$Name_submitted
cols <- c("essential_oils", "fences", "fibres","fodder","food","resins_gums","rubber_latex", "fuel")
tab <- table(usos1$Name_submitted[check_these], usos1$uses[check_these])
linhas <- apply(tab[,cols], 1, sum) > 0 
tab[linhas, cols] # ok in 30/04/2021

#How many times were they cited as being used for timber
usos3 <- usos3[!grepl("Coradin", usos3$source) , ] 
usos3$obs <- gsub("list provided by GTA. Note: I haven’t updated it for taxonomy so there might be some issues", NA, usos3$obs)
usos4 <- as.data.frame(table(usos3$Name_submitted))
names(usos4)[1:2] <- c("species.correct2", "times.cites")
usos4$sources <- aggregate(usos3$source, list(usos3$Name_submitted), 
                           function (x) paste(x, collapse = " | "))$x
usos4$obs <- aggregate(usos3$obs, list(usos3$Name_submitted), 
                           function (x) paste(x, collapse = " | "))$x
write.csv(usos4, "data/threat_exploited_timber_spp.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


#####################################################################################################################################################################H
#####################################################################################################################################################################H
################H
#### THREATS ####
################H

## Generating the SIS CONNECT file "threats.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/threats.csv")
head(sample, 3)
apply(sample, 2, unique)

## CREATING THE TYPES OF THREATS FOR THE ATALNTIC FOREST ##
## Residential & commercial development
# 1.1 Housing & urban areas (generalizado para a Mata Atlântica)
res1.1 <- c(stress = "1.1|1.3|2.1|2.2", stressdesc = NA, ThreatsLookup = "1.1", 
            ThreatsName = NA, ancestry = NA, ias = NA,
            internationalTrade = NA, scope = "Majority (50-90%)", 
            severity = "Slow, Significant Declines", text = NA, 
            timing = "Ongoing", virus = NA)
# 1.2 Commercial & industrial areas  (e.g. proximo de grandes centros urbanos)
res1.2 <- c(stress = "1.1|1.3|2.1|2.2", stressdesc = NA, ThreatsLookup = "1.2", 
            ThreatsName = NA, ancestry = NA, ias = NA,
            internationalTrade = NA, scope = "Minority (<50%)", 
            severity = "Slow, Significant Declines", text = NA, 
            timing = "Ongoing", virus = NA)
# 1.3 Tourism & recreation areas (e.g. proximo ao litoral)
res1.3 <- c(stress = "1.2|2.2", stressdesc = NA, ThreatsLookup = "1.3", 
            ThreatsName = NA, ancestry = NA, ias = NA,
            internationalTrade = NA, scope = "Minority (<50%)", 
            severity = "Causing/Could cause fluctuations", text = NA, 
            timing = "Ongoing", virus = NA)


## Agriculture & aquaculture
# 2.1	Annual & perennial non-timber crops
#2.1.2	Small-holder farming (e.g. Araucaria)
agr2.1.2 <- c(stress = "1.1|1.2|2.1|2.2|2.3", stressdesc = NA, ThreatsLookup = "2.1.2", 
              ThreatsName = NA, ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Majority (50-90%)", 
              severity = "Slow, Significant Declines", text = NA, 
              timing = "Ongoing", virus = NA) # Can be "Past, Unlikely to Return" for some species
#2.1.3	Agro-industry farming (e.g. Araucaria, Ilex paraguariensis, Euterpe edulis)
agr2.1.3 <- c(stress = "1.1|1.2|2.1|2.2|2.3", stressdesc = NA, ThreatsLookup = "2.1.3", 
              ThreatsName = NA, ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Majority (50-90%)", 
              severity = "Rapid Declines", text = NA, 
              timing = "Ongoing", virus = NA)
#2.1.4  Scale Unknown/Unrecorded (mais usado que os anteriored na AmSul)
agr2.1.4 <- c(stress = "1.1|1.2|2.1|2.2|2.3", stressdesc = NA, ThreatsLookup = "2.1.4", 
              ThreatsName = NA, ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Majority (50-90%)", 
              severity = "Rapid Declines", text = NA, 
              timing = "Ongoing", virus = NA)


# 2.2	Wood & pulp plantations
#2.2.2	Agro-industry plantations (e.g. região de Aracruz)
agr2.2.2 <- c(stress = "1.1|2.1|2.3", stressdesc = NA, ThreatsLookup = "2.2.2", 
              ThreatsName = NA, ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Minority (<50%)", 
              severity = "Causing/Could cause fluctuations", text = NA, 
              timing = "Ongoing", virus = NA)

#2.2.3	Scale Unknown/Unrecorded
agr2.2.3 <- c(stress = "1.1|2.1|2.3", stressdesc = NA, ThreatsLookup = "2.2.3", 
              ThreatsName = NA, ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Minority (<50%)", 
              severity = "Causing/Could cause fluctuations", text = NA, 
              timing = "Ongoing", virus = NA)


# 2.3 Livestock farming & ranching
#2.3.3	Agro-industry grazing, ranching or farming
live2.3.3 <- c(stress = "1.1|1.2|2.1|2.2|2.3", stressdesc = NA, ThreatsLookup = "2.3.3", 
               ThreatsName = NA, ancestry = NA, ias = NA,
               internationalTrade = NA, scope = "Majority (50-90%)", 
               severity = "Rapid Declines", text = NA, 
               timing = "Ongoing", virus = NA)
#2.3.4	Scale Unknown/Unrecorded (mais seguro/usado que o anterior)
live2.3.4 <- c(stress = "1.1|1.2|2.1|2.2|2.3", stressdesc = NA, ThreatsLookup = "2.3.4", 
               ThreatsName = NA, ancestry = NA, ias = NA,
               internationalTrade = NA, scope = "Majority (50-90%)", 
               severity = "Rapid Declines", text = NA, 
               timing = "Ongoing", virus = NA)


## 3	Energy production & mining
# 3.2	Mining & quarrying
mine3.2 <- c(stress = "1.1|1.3|2.1|2.2", stressdesc = NA, ThreatsLookup = "3.2", 
             ThreatsName = NA, ancestry = NA, ias = NA,
             internationalTrade = NA, scope = "Minority (<50%)", 
             severity = "Causing/Could cause fluctuations", text = NA, 
             timing = "Ongoing", virus = NA)


## 4	Transportation & service corridors
# 4.1	Roads & railroads
road4.1 <- c(stress = "1.1|1.2", stressdesc = NA, ThreatsLookup = "3.2", 
             ThreatsName = NA, ancestry = NA, ias = NA,
             internationalTrade = NA, scope = "Minority (<50%)", 
             severity = "Causing/Could cause fluctuations", text = NA, 
             timing = "Ongoing", virus = NA) 


## 5	Biological resource use
# 5.2	Gathering terrestrial plants
#5.2.1	Intentional use (species is the target) e.g. Araucaria
bio5.2.1 <- c(stress = "2.2|2.3", stressdesc = NA, ThreatsLookup = "5.2.1", 
              ThreatsName = NA, ancestry = NA, ias = NA, 
              internationalTrade = "No", scope = "Unknown", 
              severity = "Causing/Could cause fluctuations", text = NA, 
              timing = "Ongoing", virus = NA)

# 5.3	Logging & wood harvesting
#5.3.1	Intentional use: (subsistence/small scale) [harvest]
wood5.3.1 <- c(stress = "2.1|2.3", stressdesc = NA, ThreatsLookup = "5.3.1", 
               ThreatsName = NA, ancestry = NA, ias = NA, 
               internationalTrade = "No", scope = "Minority (<50%)", 
               severity = "Causing/Could cause fluctuations", text = NA, 
               timing = "Ongoing", virus = NA)
#5.3.2	Intentional use: (large scale) [harvest]
wood5.3.2 <- c(stress = "2.1|2.3", stressdesc = NA, ThreatsLookup = "5.3.2", 
               ThreatsName = NA, ancestry = NA, ias = NA, 
               internationalTrade = "No", scope = "Majority (50-90%)", 
               severity = "Rapid Declines", text = NA, 
               timing = "Ongoing", virus = NA)
#5.3.5	Motivation Unknown/Unrecorded (mais seguro/usado que o anterior)
wood5.3.5 <- c(stress = "2.1|2.3", stressdesc = NA, ThreatsLookup = "5.3.5", 
               ThreatsName = NA, ancestry = NA, ias = NA, 
               internationalTrade = "No", scope = "Unknown", 
               severity = "Slow, Significant Declines", text = NA, 
               timing = "Ongoing", virus = NA)


# ## 7	Natural system modifications
# # 7.1	Fire & fire suppression
# #7.1.1	Increase in fire frequency/intensity
# fire7.1.1
# #7.1.3	Trend Unknown/Unrecorded
# fire7.1.3


## 8	Invasive and other problematic species, genes & diseases
# 8.1	Invasive non-native/alien species/diseases
#8.1.1	Unspecified species (mais usado para AmSul) - competição com gramíneas exóticas
invade8.1.1 <- c(stress = "2.3.2", stressdesc = NA, ThreatsLookup = "8.1.1", 
                 ThreatsName = NA, ancestry = NA, ias = NA, 
                 internationalTrade = NA, scope = "Unknown",  
                 severity = "Causing/Could cause fluctuations", text = NA, 
                 timing = "Ongoing", virus = NA)

# ## 11	Climate change & severe weather
# # 11.1. Habitat shifting & alteration (e.g. Mangue)
# change11.1

## CREATING COMBINED THREATS FOR ALL SPECIES ##
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
n.spp <- dim(tax)[1]
geral <- rbind.data.frame(res1.1, res1.2, # res1.3, # urban, industrial and tourism 
                          agr2.1.4, # agr2.1.2, agr2.1.3, # crops
                          # agr2.2.2, agr2.2.3, # wood & pulp
                          live2.3.4, # live2.3.3, # livestock
                          # mine3.2, # mining 
                          road4.1, # roads
                          # bio5.2.1, # Species intentional use (species is the target)
                          # wood5.3.1, wood5.3.2, wood5.3.5,
                          # invade8.1.1,
                          stringsAsFactors = FALSE)
names(geral) <- names(res1.1)
threats <- do.call(rbind.data.frame, 
                   replicate(n.spp, geral, simplify = FALSE))
threats$internal_taxon_id <- rep(tax$internal_taxon_id, each = dim(geral)[1]) 

## SPECIES INTENTIONAL USE ##
# uses: "food/alcohol", "fodder", "resins_gums/rubber_latex", "essential_oils", "fuel", "fibres"
# uses classes: 1, 2, 5, 6, 7, 8
used.spp <- unique(
  usos2$internal_taxon_id[usos2$UTEndUse.UTEndUseSubfield.UTEndUseLookup %in%
                            c(1,2,5,6,8)])
bio.uses <- do.call(rbind.data.frame,
                    replicate(length(used.spp), bio5.2.1, simplify = FALSE))
names(bio.uses) <- names(bio5.2.1)
bio.uses$internal_taxon_id <- used.spp 


## LOGGING ##
# uses: "timber", "fences"
# uses classes: 9
wood.spp <- unique( 
  usos1[usos1$uses %in% c("timber", "fences"), c("internal_taxon_id", "Name_submitted")])
# Separating species with actual from possible/potential uses
wood.spp.actual <- wood.spp[wood.spp$Name_submitted %in% usos4$species.correct2, ]
wood.spp.potential <- wood.spp[!wood.spp$Name_submitted %in% usos4$species.correct2, ]
wood.uses <- c(replicate(length(wood.spp.actual$internal_taxon_id), wood5.3.2, simplify = FALSE),
                   replicate(length(wood.spp.potential$internal_taxon_id), wood5.3.1, simplify = FALSE))
wood.uses <- do.call(rbind.data.frame, wood.uses)
names(wood.uses) <- names(wood5.3.2)
wood.uses$internal_taxon_id <- c(wood.spp.actual$internal_taxon_id, 
                                 wood.spp.potential$internal_taxon_id)

## Past large-scale harvesting woods (timing: "Past, Unlikely to Return")
#Fontes: 
#https://pt.wikipedia.org/wiki/Madeira_de_lei
#https://snif.florestal.gov.br/pt-br/especies-florestais: IN 06/2008, Decreto 5.975/2006, 6.472/2008 
#http://www.esalq.usp.br/trilhas/lei/maplei.php
madeiras.proibidas <- unique(c(
  #Wikipedia
  "Cedrela fissilis", "Carapa guianensis", "Dinizia excelsa",
  "Anadenanthera sp.", "Parapiptadenia sp.", "Piptadenia sp.",
  "Centrolobium tomentosum", "Ocotea porosa", "Tabebuia sp.", 
  "Handroanthus sp.", "Zeyheria tuberculosa", "Jacaranda sp", 
  "Dalbergia nigra", "Calophyllum brasiliense", "Hymenaea courbaril",
  "Swietenia macrophylla",  "Caesalpinia echinata", "Paubrasilia echinata", 
  "Caesalpinia ferrea","Libidibia ferrea", 
  "Platycyamus regnellii", "Aspidosperma polyneuron", 
  #SNIF
  "Bertholletia excelsa", "Swietenia macrophylla",
  "Araucaria angustifolia", "Paratecoma peroba", "Amburana acreana", 
  "Apuleia leiocarpa", "Hymenaea parvifolia",  "Hymenolobium excelsum", 
  "Melanoxylon brauna", "Peltogyne maranhensis",  "Ocotea catharinensis", 
  "Ocotea odorifera", "Mezilaurus itauba",  "Cariniana legalis", 
  "Cedrela odorata", "Virola bicuhyba", "Virola surinamensis", 
  "Euxylophora paraensis",
  #Tilha ESALQ
  "Aspidosperma cylindrocarpon", "Astronium graveolens", "Machaerium villosum", 
  "Anadenanthera colubrina", "Balfourodendron riedelianum", 
  "Copaifera langsdorffii", "Colubrina glandulosa", 
  "Tabebuia aurea", "Myracrodruon urundeuva", "Astronium urundeuva",
  "Aspidosperma pyrifolium", "Calycophyllum spruceanum",
  "Patagonula americana", "Cordia americana", "Aspidosperma ramiflorum",
  "Tabebuia heptaphylla", "Handroanthus heptaphyllus",
  "Securinega guaraiuva", "Savia dictyocarpa",
  "Myroxylon peruiferum", "Holocalyx balansae", "Machaerium scleroxylon",
  "Lafoensia glyptocarpa", "Esenbeckia leiocarpa"))

proibidas.id <- wood.spp$internal_taxon_id[wood.spp$Name_submitted %in% madeiras.proibidas]
wood.uses$timing[wood.uses$internal_taxon_id %in% proibidas.id] <- 
  "Past, Unlikely to Return"
  

## TOURISM ?? ##

## INDUSTRIAL PLANTATIONS ?? ##

## MINING ?? ##

## COMPETION WITH EXOTIC GRASSES ?? ##

## COMBINING ALL THREATS 
threats <- rbind.data.frame(threats, bio.uses, wood.uses, 
                            stringsAsFactors = FALSE)
threats <- threats[order(as.double(gsub("sp","",threats$internal_taxon_id)), 
                         threats$ThreatsLookup), ] 
#Renaming
names(threats)[1:2] <- paste0("Threats.ThreatsSubfield.StressesSubfield.", 
                            names(threats)[1:2])
names(threats)[3:12] <- paste0("Threats.ThreatsSubfield.", 
                              names(threats)[3:12])
table(names(threats) == names(sample))

## Saving the SIS Connet file
write.csv(threats, "data/sis_connect/threats_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


#### CONTINUAR DAQUI OS CSVs ####

#####################################################################################################################################################################H
#####################################################################################################################################################################H
#####################H
#### CONSERVATION ####
#####################H

## Generating the SIS CONNECT file "conservationneeded.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/conservationneeded.csv")
head(sample, 3)
apply(sample, 2, unique)

## CREATING THE TYPES OF CONSERVATION NEEDS FOR THE ATLANTIC FOREST ##
# 1.1 Site/area protection (generalizado para a Mata Atlântica)
con1.1 <- c(ConservationActionsLookup = "1.1", ConservationActionsName = NA, 
            note = "Criar ou expandir unidades de conservação (Categorias IUCN  I-VI) dentro da área de distribuição da espécie")
# 1.2	Resource & habitat protection
con1.2 <- c(ConservationActionsLookup = "1.2", ConservationActionsName = NA, 
            note = "Proteção dos habitats aos quais a espécie está associada, mesmo que fora de Unidades de Conservação")
# 2.3	Habitat & natural process restoration
con2.3 <- c(ConservationActionsLookup = "2.3", ConservationActionsName = NA, 
            note = "Recuperar a degradação e/ou restaurar habitats ou funções do ecossistema aos quais a espécie está associada")
# 3.1.1	Harvest management
con3.1.1 <- c(ConservationActionsLookup = "3.1.1", ConservationActionsName = NA, 
            note = "Regulamentar o volume e frequência da colheita de produtos florestais não-madeireiros associados à espécie")
# 3.1.2	Trade management
con3.1.2 <- c(ConservationActionsLookup = "3.1.2", ConservationActionsName = NA, 
            note = "Regulamentar o comércio de produtos florestais não-madeireiros associados à espécie")
# 3.3.1	Reintroduction
con3.3.1 <- c(ConservationActionsLookup = "3.3.1", ConservationActionsName = NA, 
            note = "Reintrodução (i.e. plantio) de indivíduos em áreas onde a espécie foi extirpada, seguindo as orientações da IUCN")
# 3.4	Ex-situ conservation
con3.4 <- c(ConservationActionsLookup = "3.4", ConservationActionsName = NA, 
            note = "Proteção ex-situ da espécie via plantios em jardins botânicos e incorporação em bancos de germoplasma e/ou bases de dados genéticos (e.g. GenBank)")
# 4.2	Training
con4.2 <- c(ConservationActionsLookup = "4.2", ConservationActionsName = NA, 
            note = "Capacitação de profissionais que possam atuar diretamente em ações de conservação da espécie (e.g. gestores de unidades de conservação, agentes de fiscalização, produtores de semente), aprimorando conhecimentos, habilidades e troca de informações entre esses profissionais (e.g. melhorando habilidades de identificação da espécie)")
# 4.3	Awareness & communications
con4.3 <- c(ConservationActionsLookup = "4.3", ConservationActionsName = NA, 
            note = "Sensibilização ambiental da sociedade em geral sobre o grau de ameaça da espécie via diferentes meios de comunicação")
# 5.1.1	International level
con5.1.1 <- c(ConservationActionsLookup = "5.1.1", ConservationActionsName = NA, 
            note = "Criar, implementar ou adaptar leis nacionais que oficializem a proteção da espécie")
# 5.1.2	National level
con5.1.2 <- c(ConservationActionsLookup = "5.1.2", ConservationActionsName = NA, 
            note = "Criar, implementar ou adaptar dispositivos inter-nacionais que facilitem a proteção da espécie")
# 5.2	Policies and regulations
con5.2 <- c(ConservationActionsLookup = "5.2", ConservationActionsName = NA, 
            note = "Criar, implementar ou adaptar políticas públicas que facilitem a proteção da espécie")
# 5.4.1	International level
con5.4.1 <- c(ConservationActionsLookup = "5.4.1", ConservationActionsName = NA, 
            note = "Reforçar e garantir a implementação de tratados e regulamentações internacionais que envolvem a espécie, incluindo a implementação das sanções e punições legais cabíveis")
# 6.2	Substitution
con6.2 <- c(ConservationActionsLookup = "6.2", ConservationActionsName = NA, 
            note = "Promover o uso de produtos oriundos de espécies manejadas e que causem menores pressões sobre as populações naturais da espécie")
# 6.3	Market forces
con6.3 <- c(ConservationActionsLookup = "6.3", ConservationActionsName = NA, 
            note = "Promover o uso de mecanismos de mercado que promovam mudanças de comportamento em relação ao uso da espécie (e.g. certificação florestal)")
# 7.1	Institutional & Civil Society Development
con7.1 <- c(ConservationActionsLookup = "7.1", ConservationActionsName = NA, 
            note = "Dar suporte e capacitar organizações, empresas e comunidades que possam atuar diretamente em ações para a conservação da espécie")
# 7.2	Alliance & Partnership Development
con7.2 <- c(ConservationActionsLookup = "7.2", ConservationActionsName = NA, 
            note = "Formar redes de colaboração multi-setoriais para divulgar e promover a conservação da espécie")


## CREATING COMBINED THREATS FOR ALL SPECIES ##
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
# end <- 
# cites <- 
# assess <-   
n.spp <- dim(tax)[1]

## LAW & POLICY ##
#5.1.1	International level	depend (if threatened and not endemic)

#5.1.2	National level	depend (if threatened and Brazilian endemic)
# national.spp <- unique()
national <- do.call(rbind.data.frame,
                    replicate(length(national.spp), con5.1.2, simplify = FALSE))
names(national) <- names(con5.1.2)
national$internal_taxon_id <- national.spp 

## INTERNATIONAL TRADE ##
#5.4.1 International level (if CITES)


## SPECIES ASSESSED AS CR or CR (PE) ##
# CR.PE.spp <- 
CR.PE <- rbind.data.frame(con1.1,
                          con1.2,
                          con2.3,
                          con3.4,
                          con4.2,
                          con4.3,
                          con7.1,
                          con7.2,
                          stringsAsFactors = FALSE)
names(CR.PE) <- names(con1.1)
CR <- do.call(rbind.data.frame, 
                   replicate(length(CR.PE.spp), geral, simplify = FALSE))
CR$internal_taxon_id <- rep(CR.PE.spp, each = length(CR.PE.spp))

## THREATENED TIMBER SPECIES ##

## THREATENED NOT WITHIN CONSERVATION UNITS ##

## THREATENED OVER-EXPLOITED SPECIES ##





## RENAMING
names(conservation)[1:3] <- paste0("ConservationActions.ConservationActionsSubfield.", 
                                   names(conservation)[1:3])
table(names(conservation) == names(sample))


## Saving the SIS Connet file
write.csv(conservation, "data/sis_connect/conservationneeded_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

#####################################################################################################################################################################H
#####################################################################################################################################################################H
########################H
#### RESEARCH NEEDED ####
########################H

## Generating the SIS CONNECT file "researchneeded.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/researchneeded.csv")
head(sample, 3)
apply(sample, 2, unique)


## CREATING THE TYPES OF RESEARCH NEEDS FOR THE ATLANTIC FOREST ##
# 1	Research					
#1.1	Taxonomy depend (if most records are old)				
res1.1 <- c(ResearchLookup = "1.1", ResearchName = NA, 
            note = "Espécie conhecida apenas para a localidade onde o espécimen tipo foi coletado e que não tenha sido re-coletado nos últimos 50 anos. Pode se tratar de uma sinônimo de outra espécie que não tenha recebido o tratamento taxonômico necessário")
#1.2 Population size, distribution & trends	depend (if no plot data is available)
res1.2 <- c(ResearchLookup = "1.2", ResearchName = NA, 
            note = "Espécie não registrada em nehum dos mais de 1100 inventários florestais compilados para a Mata Atlântica, portanto sem estimativas de seu tamanho populacional e possíveis tendências de declínio. Provavelmente representa uma espécie que ocorre em baixas densidades na floresta ou uma espécie associada a um tipo de habitat específico e sub-representado nos inventários florestais compilados")
#1.3	Life history & ecology					
res1.3 <- c(ResearchLookup = "1.3", ResearchName = NA, note = NA)
#1.4	Harvest, use & livelihoods					
res1.4 <- c(ResearchLookup = "1.4", ResearchName = NA, note = NA)
#1.5	Threats					
res1.5 <- c(ResearchLookup = "1.5", ResearchName = NA, note = NA)
#1.6	Conservation Actions				Research to determine how to mitigate particular threats (e.g., how to mitigate impacts of long-lining to reduce bycatch), or whether ex-situ breeding is possible	
res1.6 <- c(ResearchLookup = "1.6", ResearchName = NA, note = NA)

# 2	Conservation Planning					
#2.1	Species Action/Recovery Plan	depend (if CR or CR PE)			
res2.1 <- c(ResearchLookup = "2.1", ResearchName = NA, 
            note = "Estabelecer um Plano de Ação em escala nacional para recuperar as populações da espécie")
#2.2	Area-based Management Plan					
res2.2 <- c(ResearchLookup = "2.2", ResearchName = NA, note = NA)
#2.3	Harvest & Trade Management Plan					
res2.3 <- c(ResearchLookup = "2.3", ResearchName = NA, note = NA)

# 3	Monitoring					
#3.1	Population trends	depend (if threatened based on C or D)	
res3.1 <- c(ResearchLookup = "3.1", ResearchName = NA, 
            note = "Espécie com indicativos de tamanho populacional pequeno e em declínio ou com tamanho populacionais muito pequenos, o que é raro entre as espécies de árvores com ao menos um registro nos inventários florestais compilados. Estudos populacionais mais refinados são necessários para melhor estimar qual o tamanho populacional da espécie ao longo da sua área de distribuição")
#3.2	Harvest level trends	depend (if overharvested)
res3.2 <- c(ResearchLookup = "3.2", ResearchName = NA, 
            note = "Espécie com histórico de exploração na Mata Atlântica. Apesar das estimativas do tamanho populacional da espécie estarem disponíveis (baseadas em perda de área florestal), pouco se conhece sobre a intensidade da exploração histórica e de seu impacto no tamanho populacional efetivo da espécie ao longo do tempo. Estudos populacionais mais refinados são necessários para melhor estimar qual a real situação das populacões ao longo da área de distribuição da espécie")
#3.3	Trade trends	depend (if traded)?				Espécie de uso comercial comum na Mata Atlântica. Contudo, pouco se conhece sobre o impacto desse uso sobre a espécie.
res3.3 <- c(ResearchLookup = "3.3", ResearchName = NA, note = NA)
#3.4	Habitat trends					
res3.4 <- c(ResearchLookup = "3.4", ResearchName = NA, note = NA)


## ONLY OLD RECORDS or LESS THAN 3##

## NO PLOT DATA AVAILABLE ##

## Species without plot data 

#This species has no records in forest plot data for the Atlantic Forest, being probably rare in the region.
#Research is recommended to estimate the population size and declines within its are of distribution.


##Specific comments for EX and EW species
#Campomanesia lundiana: only one recent record (2018) is available for Monte Cabrão, Santos, São Paulo, Brasil. 
#This record was collected and identified by I. Villaça and Mara Magenta, which are not Myrtaceae specialists.
#Thus, this new record may be a misidentification of Campomaneisa phaea, since the species listed in the associated study
#only list this species as those used by the tradional local communities for fruit consumption.

# Pouteria stenophylla: Species mainly known from its type locality (probably Serra dos Orgãos, Rio de Janeiro). Only one recent (2010) and taxonomically validated record was found for this species (Palazzo, F.M.A. 40) 
#for the restingas of Rio das Ostras county, Rio de Janeiro. Another recent record for Prado, Bahia (Rezende, S.G. 1833) was also found but its identification was not yet validated by an specialist.
#Recent records from rocky outcrops in Minas Gerais state (Ouro Preto and MAriana) are also available but again identification by family specialists are still pending.
#ver tb: PALAZZO, F.M.A.; DIAS-NETO, A.O.; MONTEIRO, M.H.D.A.; ANDREATA, R.H.P. 2010. Sinopse comentada de Sapotaceae no município de Rio das Ostras (RJ, BRASIL). Pesquisas (Botânica) 61: 293-306.


## Saving the SIS Connet file
write.csv(research, "data/sis_connect/researchneeded_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



#####################################################################################################################################################################H
#####################################################################################################################################################################H
#################H
#### CREDITS  ####
#################H
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/credits.csv")
apply(sample, 2, unique)

## Creating the assessor(s) vector
assess1 <- data.frame(Order = 1L, 
                      affiliation = "Naturalis Biodiversity Center & IUCN SSC Global Tree Specialist Group", 
                      credit_type = "Assessor", email = "renato.lima@naturalis.nl",
                      firstName = "Renato A. Ferreira de",
                      initials = "R.A.F.", lastName = "Lima", nickname = "",
                      user_id = "") 

## Creating the compiler(s) vector
compila1 <- data.frame(Order = 1L, 
                      affiliation = "Naturalis Biodiversity Center & IUCN SSC Global Tree Specialist Group", 
                      credit_type = "Compiler", email = "renato.lima@naturalis.nl",
                      firstName = "Renato A. Ferreira de",
                      initials = "R.A.F.", lastName = "Lima", nickname = "",
                      user_id = "") 


## Creating the reviewer(s) vector
review1 <- data.frame(Order = 1L, affiliation = "",  credit_type = "Reviewer", 
                      email = "", firstName = "", initials = "", lastName = "", 
                      nickname = "", user_id = "") 

## Creating the reviewer(s) vector
contrib1 <- data.frame(Order = 1L, 
                       affiliation = "Centro Nacional de Conservação da Flora (CNCFlora) & IUCN SSC Brazil Plant Red List Authority",
                       credit_type = "Contributor", 
                       email = "", firstName = "", initials = "", lastName = "", 
                       nickname = "", user_id = "") 

## Combining all info
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
n.spp <- dim(tax)[1]
temp <- rbind.data.frame(assess1,  review1,  contrib1, 
                         stringsAsFactors = FALSE)
credits <- do.call(rbind.data.frame, 
                   replicate(n.spp, temp, simplify = FALSE))
credits$internal_taxon_id <- rep(tax$internal_taxon_id, each = dim(temp)[1]) 

#Checking the format
sample[sample$internal_taxon_id %in% "22823",]
credits[credits$internal_taxon_id %in% "sp1",]

# Saving
credits <- credits[ , match(names(sample), names(credits))]
write.csv(credits, "data/sis_connect/credits_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



#####################################################################################################################################################################H
#####################################################################################################################################################################H
###################H
#### REFERENCES ####
###################H

## Generating the SIS CONNECT file "references.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/references.csv")
head(sample, 3)
apply(sample[,1:10], 2, unique)

## Taxonomic references 
full.tax <- readRDS("data/threat_full.taxonomy.rds")

## Brazilian Flora
tax1 <- data.frame(Reference_type = "Taxonomy", 
                   access_date = "6 May 2021", alternate_title = NA, 
                   author = "Flora do Brasil 2020 em construção",
                   date = NA, edition = NA, externalBibCode = NA, isbnissn = NA,
                   keywords = NA, number = NA, number_of_volumes = NA, 
                   pages = NA, place_published = NA, publisher = NA, 
                   secondary_author = NA, secondary_title = NA, section = NA, 
                   short_title = NA, submission_type = NA, 
                   title = "Jardim Botânico do Rio de Janeiro", 
                   type = "electronic source", volume = NA, year = NA)
spp <- full.tax$internal_taxon_id[!is.na(full.tax$TaxonomicReference1)]
tax1 <- do.call(rbind.data.frame, 
                   replicate(length(spp), tax1, simplify = FALSE))
tax1$internal_taxon_id <- spp
url <- full.tax$TaxonomicReference1[!is.na(full.tax$TaxonomicReference1)]
url <- gsub(".*Available at\\: ", "", url, perl = TRUE)
url <- gsub("\\.$", "", url, perl = TRUE)
tax1$url <- url
tax1 <- tax1[ , match(names(sample), names(tax1))]

## Tropicos
tax2 <- data.frame(Reference_type = "Taxonomy", 
                   access_date = "6 May 2021", alternate_title = NA, 
                   author = "Tropicos.org",
                   date = NA, edition = NA, externalBibCode = NA, isbnissn = NA,
                   keywords = NA, number = NA, number_of_volumes = NA, 
                   pages = NA, place_published = NA, publisher = NA, 
                   secondary_author = NA, secondary_title = NA, section = NA, 
                   short_title = NA, submission_type = NA, 
                   title = "Missouri Botanical Garden", 
                   type = "electronic source", volume = NA, year = NA)
spp2 <- full.tax$internal_taxon_id[!is.na(full.tax$TaxonomicReference2)]
tax2 <- do.call(rbind.data.frame, 
                replicate(length(spp2), tax2, simplify = FALSE))
tax2$internal_taxon_id <- spp2
url <- full.tax$TaxonomicReference2[!is.na(full.tax$TaxonomicReference2)]
url <- gsub(".*Available at\\: ", "", url, perl = TRUE)
url <- gsub("\\.$", "", url, perl = TRUE)
tax2$url <- url
tax2 <- tax2[ , match(names(sample), names(tax2))]

## GBIF
tax3 <- data.frame(Reference_type = "Taxonomy", 
                   access_date = "6 May 2021", alternate_title = NA, 
                   author = "GBIF Secretariat",
                   date = NA, edition = NA, externalBibCode = NA, isbnissn = NA,
                   keywords = NA, number = NA, number_of_volumes = NA, 
                   pages = NA, place_published = NA, publisher = NA, 
                   secondary_author = NA, secondary_title = NA, section = NA, 
                   short_title = NA, submission_type = NA, 
                   title = "GBIF Backbone Taxonomy", 
                   type = "electronic source", volume = NA, year = NA)
spp3 <- full.tax$internal_taxon_id[!is.na(full.tax$TaxonomicReference3)]
tax3 <- do.call(rbind.data.frame, 
                replicate(length(spp3), tax3, simplify = FALSE))
tax3$internal_taxon_id <- spp3
url <- full.tax$TaxonomicReference3[!is.na(full.tax$TaxonomicReference3)]
url <- gsub(".*Available at\\: ", "", url, perl = TRUE)
url <- gsub("\\.$", "", url, perl = TRUE)
tax3$url <- url
tax3 <- tax3[ , match(names(sample), names(tax3))]


#####################################################################################################################################################################H
#####################################################################################################################################################################H
####################H
#### ASSESSMENTS ####
####################H

## Generating the SIS CONNECT file "assessments.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/assessments.csv")
head(sample[,8:17], 4)
# apply(sample, 2, unique)

#BiogeographicRealm.realm (Afrotropical, Antarctic, Australasian, Indomalayan, Nearctic, Neotropical, Oceanian, Palearctic)



#####################################################################################################################################################################H
#####################################################################################################################################################################H
###################H
#### ALL FIELDS ####
###################H

## Generating the SIS CONNECT file "allfields.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/allfields.csv")
head(sample, 3)
# apply(sample, 2, unique)



?red::countries
