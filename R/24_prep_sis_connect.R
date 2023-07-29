################################################H
#### BUILDING THE IUCN SIS CONNECT CSV FILES ####
################################################H
rm(list = ls())

# READING NECESSARY FILES ------------------------------------------------------
source("./R/98_internal_functions.R")

## Reading the THREAT species list
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]

## Reading THREAT full taxonomy and taxonomic references
full.tax <- readRDS("data/threat_full.taxonomy.rds")
ncbi.taxonomy <- readRDS("data/threat_ncbi_taxonomy.rds")
# tps.taxonomy1 <- readRDS("data/threat_tropicos_taxonomy.rds")
iucn2 <- readRDS("data/threat_iucn_taxonomy.rds")

## Reading THREAT species descriptions
spp.info <- readRDS("data/threat_flora_info.rds")

## Reading THREAT species habitat data
hab1 <- read.csv("data/threat_habitats_preliminar1.csv")

## Reading THREAT from the Brazilian Flora Online (FBO)
taxon2 <- readRDS("data/threat_fbo_tax_info.rds")
fbo.info <- taxon2
 
## Replacing new synonyms from the Brazilian Flora Online (FBO) 
syn.br <- read.csv("data/new_synonyms_floraBR.csv", 
                  na.strings = c(""," ",NA), as.is = TRUE)

## Reading THREAT species herbarium data
occs <- readRDS("data/threat_species_by_country.rds")

## Reading the global IUCN assessments (IUCN) 
#Citation: IUCN 2020. IUCN Red List of Threatened Species. Version 2020-1 <www.iucnredlist.org>
# iucn <- read.csv("IUCN_2020_assessments.csv", as.is = TRUE, na.string = c(NA,""," "))
#Citation: IUCN 2021. IUCN Red List of Threatened Species. Version 2021-1 <www.iucnredlist.org>
# iucn <- read.csv("IUCN_2021_v1_assessments.csv", as.is = TRUE, na.string = c(NA,""," "))
#Citation: IUCN 2023. IUCN Red List of Threatened Species. Version 2022-2 <www.iucnredlist.org>
iucn <- readRDS("IUCN_2022_v2_assessments_THREAT.rds")
iucn.sis <- readRDS("IUCN_2022_v2_sis_connect_THREAT.rds")

## Reading species endemism levels in respect to the AF (from a)
# af.list <- read.csv("C:/Users/renat/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/AppendixC_af_checklist.csv", as.is=TRUE)
# end <- read.csv("C:/Users/renat/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/AppendixF_endemism_levels.csv", as.is=TRUE)
af.list <- read.csv("./data/outside_data/AppendixC_af_checklist.csv", as.is=TRUE)
end <- read.csv("./data/outside_data/AppendixF_endemism_levels.csv", as.is=TRUE)

## Reading probable occurrences of tree species in the AF (from a)
# end.prob <- read.csv("C:/Users/renat/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/AppendixD_probable_occurrences.csv", as.is=TRUE)[,1:5]
end.prob <- read.csv("./data/outside_data/AppendixD_probable_occurrences.csv", as.is=TRUE)[,1:5]

## Reading the CITES/EU listings
#LOADING BELOW, DIRECTLY WITH THE CODES
# cites1 <- readRDS("data/threat_cites.rds")
# EU1 <- readRDS("data/EU_listings_data.rds")

## Reading TreeCo species use database and associated functions
# path = "C:/Users/renat/Documents/raflima/Pos Doc/Databases/TreeCo Database Management"
# source(paste(path,"uses_data_prep.R",sep="/"))
# usos <- read.csv("C://Users//renat//Documents//raflima//Pos Doc//Databases//Species Uses//plant_uses.csv", 
#                  na.strings = c(""," ",NA), encoding = "UTF-8")

source("./data/outside_data/uses_data_prep.R")
usos <- read.csv("./data/outside_data/plant_uses.csv", 
                 na.strings = c(""," ",NA), encoding = "UTF-8")

## Reading THREAT assessments
all.crit <- readRDS("data/all.criteria.rds")
all.crit <- all.crit[, -142]

## CITES
cites <- readRDS("data/threat_cites_EU.rds")

## CREATING COMBINED INFO FOR ALL SPECIES ##
all.spp <- merge(all.crit, tax,
                 by.x = "species", by.y = "species.correct2", all.x = TRUE, sort = FALSE)
all.spp <- merge(all.spp, cites[ ,c("internal_taxon_id", "cites_listing", "appendix", "annotation")],
                 by.x = "internal_taxon_id", 
                 by.y = "internal_taxon_id", all.x = TRUE, sort = FALSE)
all.spp <- all.spp[order(as.double(gsub("sp", "", all.spp$internal_taxon_id))), ]
saveRDS(all.spp, "data/all_spp_crits_tax_cites.rds")



###############################################################################H
###############################################################################H
# TAXONOMY ---------------------------------------------------------------------
## Filtering and creating the final columns
final.cols <- c("kingdom","phylum","classname","ordername","family","genus",
                "species.correct2","scientificNameAuthorship")
full.tax1 <- full.tax[ ,final.cols]
#,"accepted.name","taxon.status","name.status","notes","nomenclaturestatusname",
#"displayreference","displaydate",
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
full.tax1$TaxonomicNotes.value[ids] <- paste0(full.tax1$TaxonomicNotes.value[ids], 
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
saveRDS(full.tax1, "data/threat_final_taxonomy.rds")

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

## Replacing NA by empty characters
for (i in seq_len(dim(full.tax3)[2])) {
  check_these <- is.na(full.tax3[,i])
  if (any(check_these))
    full.tax3[check_these, i] <- ""
}

## Comparing with previous assessments
cols_names <- c("genus", "species", "infraType", "infra_name")
csv.sp.names  <- .build.name(full.tax3, cols_names)
full.tax4 <- .merge_sis_connect(full.tax3, type = "taxonomy", 
                                iucn, x.spp = csv.sp.names, replace = TRUE)

## Replacing NA by empty characters for the new data frame
for (i in seq_len(dim(full.tax4)[2])) {
  check_these <- full.tax4[,i] %in% c("NA", "", NA)  
  if (any(check_these))
    full.tax4[check_these, i] <- ""
}

## Saving
full.tax3 <- full.tax3[order(as.double(gsub("sp", "", full.tax3$internal_taxon_id))), ] 
write.csv(full.tax3, "data/sis_connect/taxonomy_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

full.tax4 <- full.tax4[order(as.double(gsub("sp", "", full.tax4$internal_taxon_id))), ] 
write.csv(full.tax4, "data/sis_connect/taxonomy_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



###############################################################################H
###############################################################################H
# COMMON NAMES -----------------------------------------------------------------

## Filtering and editing
common <- merge(full.tax[ ,c("species.correct2","vernacular.name")], tax[,c("internal_taxon_id", "species.correct2")],
                by = "species.correct2", all.x = TRUE, sort = FALSE)
table(common$internal_taxon_id == full.tax2$internal_taxon_id)
names(common) <- c("validName", "name0", "internal_taxon_id")
common1 <- common[!is.na(common$name0),]

## Solving encoding problems
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

### minor fixes
common1$name0 <- gsub("pao \\/ pau \\/ pao", "Pau", common1$name0, ignore.case = TRUE)
common1$name0 <- gsub("pao \\/ pau", "Pau", common1$name0, ignore.case = TRUE)
common1$name0 <- gsub("cumixá/cumichá", "Cumixá", common1$name0, ignore.case = TRUE)
common1$name0 <- gsub("Caixeta/Caixeto", "Caixeta", common1$name0, ignore.case = TRUE)

## Breaking the info
tutu <- strsplit(common1$name0, "\\|")
# tutu <- lapply(tutu, stringr::str_squish)
tutu <- lapply(tutu, function(x) gsub("\\s+", " ", x, perl = TRUE))
tutu <- lapply(tutu, function(x) gsub("^ | $", "", x, perl = TRUE))
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

## Organizing 
commonnames <- merge(tutu1, common1[,1:3], by = "validName", all.x = TRUE)
commonnames <- commonnames[order(commonnames$validName), ]
commonnames1 <- commonnames[,match(names(sample), names(commonnames), nomatch = 0)]
commonnames1$primary <- "false"


## Adding common names from previous assessments
iucn.sis.common <- iucn.sis[iucn.sis$sic.connect.file %in% "common_names", ]
commonnames2 <- .merge_sis_connect(commonnames, type = "commonnames", tax.table = tax,
                                   iucn.sis.common, x.spp = commonnames$validName, replace = TRUE)
# commonnames2 <- commonnames2[order(commonnames2$validName), ]
commonnames2$primary <- "false"
commonnames2 <- 
  commonnames2[,c("internal_taxon_id", "language", "name", "primary", "source.info")]

## Saving
write.csv(commonnames1, "data/sis_connect/commonnames_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(commonnames2, "data/sis_connect/commonnames_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


###############################################################################H
###############################################################################H
# SYNONYMS ---------------------------------------------------------------------

## Filtering and editing
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

## Replacing NA by empty characters
for (i in seq_len(dim(synonyms3)[2])) {
  check_these <- is.na(synonyms3[,i])
  if (any(check_these))
    synonyms3[check_these, i] <- ""
}

## Adding common names from previous IUCN assessments
iucn.sis.synonyms <- iucn.sis[iucn.sis$sic.connect.file %in% "synonyms", ]
synonyms3.1 <- dplyr::left_join(synonyms3, unique(synonyms1[,c("internal_taxon_id", "validName")]))
synonyms4 <- .merge_sis_connect(synonyms3, type = "synonyms", tax.table = tax,
                                iucn.sis.synonyms, x.spp = synonyms3.1$validName, replace = TRUE)

## Replacing NA by empty characters for the new data frame
for (i in seq_len(dim(synonyms4)[2])) {
  check_these <- is.na(synonyms4[,i])
  if (any(check_these))
    synonyms4[check_these, i] <- ""
}
synonyms4 <- synonyms4[, c("infraType", "infrarankAuthor", "infrarankName", 
                           "internal_taxon_id", "name", "speciesAuthor", 
                           "speciesName", "source.info")]

## Saving
write.csv(synonyms3, "data/sis_connect/synonyms_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(synonyms4, "data/sis_connect/synonyms_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



###############################################################################H
###############################################################################H
# HABITATS ---------------------------------------------------------------------

## Habitat Information (field 'HabitatDocumentation.narrative')
hab1$HabitatDocumentation.narrative <- NA

#GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup:	(Coding Option, e.g. 3.4)
## 1. Forest & Woodland: 
# 1.4 Temperate Forest
# 1.5 Subtropical/Tropical Dry Forest: includes Deciduous forests and Chaco
# 1.6 Subtropical/Tropical Moist Lowland Forest: below 1.200 m
# 1.8 Subtropical/Tropical Swamp Forest: typically flooded for at least part of the year and dependent on this flooding for its existence
# 1.9 Subtropical/Tropical Moist Montane Forest: above 1.200 m (include cloud forests)
## 2. Savanna
# 2.1 Dry Savanna: includes cerrado/campos/caatinga (Brazil, Guiana), chaco seco (Argentina/Uruguay)
# 2.2 Moist Savanna: includes pantanal (Brazil/Bolivia/Paraguay), chaco húmedo (Paraguay/Bolivia/Argentina), and llanos (Venezuela/Colombia).
## 3. Shrubland
# 3.5 Subtropical/Tropical Dry Shrubland: includes restingas (coastal scrub in Brazil); xerophyllous scrub; chapparal (Mexico and SW US); matorral sammofilo, matorral seco (Argentina, Paraguay, Uruguay).
# 3.7 Subtropical/Tropical High Altitude Shrubland: included Campo Rupestre areas by in Jung et al. (2020)
## 4. Grassland
# 4.7 Subtropical/Tropical High Altitude Grassland: included here as Campos Rupestres, but Jung et al. 2020 also classify them as 3.7

## Loading the table with the conversions
hab_equiv <- read.csv("./SIS_sample_6_1/IUCN_habitat_classification_and_equivalencies.csv", as.is = TRUE)

hab_equiv$vegetation.type.reflora <- 
  gsub("(", "\\(", hab_equiv$vegetation.type.reflora, fixed = TRUE)
hab_equiv$vegetation.type.reflora <- 
  gsub(")", "\\)", hab_equiv$vegetation.type.reflora, fixed = TRUE)
hab_equiv$vegetation.type.reflora <- 
  gsub("=", "\\=", hab_equiv$vegetation.type.reflora, fixed = TRUE)
hab_equiv1 <- hab_equiv[!is.na(hab_equiv$vegetation.type.reflora), ] 

hab1$vt <- hab1$vegetation.type.reflora

for (i in seq_along(hab_equiv1$vegetation.type.reflora)) {
  patt <- hab_equiv1$vegetation.type.reflora[i]
  rplc <- hab_equiv1$code[i]
  hab1$vt <- gsub(patt, rplc, hab1$vt, perl = TRUE)
}
# forests <- 'Floresta Ombrófila \\(\\= Floresta Pluvial\\)|Floresta Estacional Semidecidual|Floresta Ombrófila Mista|Floresta de Terra Firme|Floresta Ciliar ou Galeria|Floresta Estacional Perenifólia'
# hab1$vt <- gsub(forests, "1.6", hab1$vt, perl = TRUE)
# savanas <- 'Cerrado \\(lato sensu\\)|Caatinga \\(stricto sensu\\)|Campo Limpo|Savana Amazônica'
# hab1$vt <- gsub(savanas, "2.1", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("Carrasco|Restinga|Campinarana", "3.5", hab1$vt)
# hab1$vt <- gsub("Campo rupestre|Campo de Altitude|Campo de Altitude", "4.7", hab1$vt)
# hab1$vt <- gsub("Floresta Estacional Decidual|Chaco", "1.5", hab1$vt)
# hab1$vt <- gsub("Floresta de Várzea|Floresta de Igapó", "1.8", hab1$vt)
# hab1$vt <- gsub("Vegetação Sobre Afloramentos Rochosos", "6", hab1$vt)
# hab1$vt <- gsub("Área Antrópica", "14.5", hab1$vt)
# hab1$vt <- gsub("Manguezal", "1.7", hab1$vt)
# hab1$vt <- gsub("Campo de Várzea", "4.6", hab1$vt)
# hab1$vt <- gsub("Palmeiral", "2.1", hab1$vt)
# hab1$vt <- gsub("Vegetação Aquática", "4.6", hab1$vt)
hab1$vt <- sapply(strsplit(hab1$vt,"\\|"), 
                  function(x) paste(sort(as.double(unique(x)), na.last = T), collapse = "|", sep=""))
hab1$vt[hab1$vt %in% "NA"] <- "1.6"
names(hab1)[grepl("^vt", names(hab1))] <- "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup"

#GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName:	(e.g. Habitat description, e.g. Shrubland ->Shrubland - Temperate)
hab1$vt <- hab1$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup
for (i in seq_along(hab_equiv$code)) {
  patt <- hab_equiv$code[i]
  patt <- gsub(".", "\\.", patt, fixed = TRUE)
  rplc <- hab_equiv$name[i]
  hab1$vt <- gsub(paste0("^", patt, "$"), rplc, hab1$vt, perl = TRUE)
  hab1$vt <- gsub(paste0("^", patt, "\\|"), paste0(rplc, "|"), hab1$vt, perl = TRUE)
  hab1$vt <- gsub(paste0("\\|", patt, "$"), paste0("|", rplc), hab1$vt, perl = TRUE)
  hab1$vt <- gsub(paste0("\\|", patt, "\\|"), paste0("|", rplc, "|"), hab1$vt, perl = TRUE)
}

# hab1$vt <- gsub("^1\\.5$", "Forest - Subtropical/Tropical Dry", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^1\\.5\\|", "Forest - Subtropical/Tropical Dry|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.5$", "|Forest - Subtropical/Tropical Dry", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.5\\|", "|Forest - Subtropical/Tropical Dry|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^1\\.6$", "Forest - Subtropical/Tropical Moist Lowland", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^1\\.6\\|", "Forest - Subtropical/Tropical Moist Lowland|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.6$", "|Forest - Subtropical/Tropical Moist Lowland", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.6\\|", "|Forest - Subtropical/Tropical Moist Lowland|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^1\\.7$", "Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^1\\.7\\|", "Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.7$", "|Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.7\\|", "|Forest - Subtropical/Tropical Mangrove Vegetation Above High Tide Level|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^1\\.8$", "Forest - Subtropical/Tropical Swamp", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^1\\.8\\|", "Forest - Subtropical/Tropical Swamp|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.8$", "|Forest - Subtropical/Tropical Swamp", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|1\\.8\\|", "|Forest - Subtropical/Tropical Swamp|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^2\\.1$", "Savanna - Dry", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^2\\.1\\|", "Savanna - Dry|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|2\\.1$", "|Savanna - Dry", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|2\\.1\\|", "|Savanna - Dry|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^3\\.5$", "Shrubland - Subtropical/Tropical Dry", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^3\\.5\\|", "Shrubland - Subtropical/Tropical Dry|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|3\\.5$", "|Shrubland - Subtropical/Tropical Dry", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|3\\.5\\|", "|Shrubland - Subtropical/Tropical Dry|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^4\\.6$", "Grassland - Subtropical/Tropical Seasonally Wet/Flooded", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^4\\.6\\|", "Grassland - Subtropical/Tropical Seasonally Wet/Flooded|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|4\\.6$", "|Grassland - Subtropical/Tropical Seasonally Wet/Flooded", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|4\\.6\\|", "|Grassland - Subtropical/Tropical Seasonally Wet/Flooded|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^4\\.7$", "Grassland - Subtropical/Tropical High Altitude", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^4\\.7\\|", "Grassland - Subtropical/Tropical High Altitude|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|4\\.7$", "|Grassland - Subtropical/Tropical High Altitude", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|4\\.7\\|", "|Grassland - Subtropical/Tropical High Altitude|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^6$", "Rocky areas (eg. inland cliffs, mountain peaks)", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^6\\|", "Rocky areas (eg. inland cliffs, mountain peaks)|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|6$", "|Rocky areas (eg. inland cliffs, mountain peaks)", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|6\\|", "|Rocky areas (eg. inland cliffs, mountain peaks)|", hab1$vt, perl = TRUE)
# 
# hab1$vt <- gsub("^14\\.5$", "Artificial/Terrestrial - Urban Areas", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("^14\\.5\\|", "Artificial/Terrestrial - Urban Areas|", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|14\\.5$", "|Artificial/Terrestrial - Urban Areas", hab1$vt, perl = TRUE)
# hab1$vt <- gsub("\\|14\\.5\\|", "|Artificial/Terrestrial - Urban Areas|", hab1$vt, perl = TRUE)
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
hab1$System.value <- gsub("TerrÃ­cola", "Terricola", hab1$System.value)


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
standard.texts_en <- data.frame(GL1 = "Field measurements of generation length are largely missing for Neotropical trees. 
                             Therefore, an approach based on tree functional types was used to set proxies of species generation lengths based on their growth form and ecological group. 
                             Values of the proxies ranged from 10 (for pioneer shrubs) to 80 years (large, climax trees), which are reasonable values for tropical trees obtained using IUCN approximative method number 3 (IUCN 2019, p.29) and based on empirical evidence of tree species development in forested conditions. 
                             This species was classified as ",
                             GL2 = "More details on this approximate method to define species generation lengths can be found in Lima et al. (in prep.).",
                             stringsAsFactors = FALSE)
texto.pt1 <- c("Medidas em campo do comprimento ou tempo de geração (sensu IUCN 2019) não existem para a maioria das espécies arbóreas Neotropicais. 
Portanto, um abordagem baseada em tipos funcionais foi usado para definir valores aproximados do tempo de geração em função da forma de vida e do grupo ecológico das espécies. 
Esses valores variaram entre 10 (para arbustos pioneiros) e 80 anos (árvores altas e de final de sucessão), que são valores razoáveis considerando tempos de geração de árvores tropicais obtidos usando o método aproximativo número 3 da IUCN (IUCN 2019, p.29) e baseado nas evidências de campo sobre crescimento de espécies arbóreas tropicais em condições florestais. 
A espécie foi classificada como '")
texto.pt2 <- c("Mais detalhes sobre os métodos usados para definir o comprimento de geração aproximado podem ser encontrados em Lima et al. (in prep.).")
hab1$GenerationLength.justification <- paste0(texto.pt1, hab1$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName,
                                              "', em relação a sua forma de crescimento (altura máxima de ", hab1$MaxSize.size/100," m), e como '",
                                              hab1$ecol.group,"', em relação ao seu grupo ecológico. ", texto.pt2)

hab1$GenerationLength.justification <- 
  gsub(" \\(altura máxima de NA m\\)", "", hab1$GenerationLength.justification, perl = TRUE)
hab1$GenerationLength.justification <- 
  gsub("early_secondary", "Secundária Inicial", hab1$GenerationLength.justification, perl = TRUE)
hab1$GenerationLength.justification <- 
  gsub("late_secondary", "Secundária Tardia", hab1$GenerationLength.justification, perl = TRUE)
hab1$GenerationLength.justification <- 
  gsub("climax", "Clímax", hab1$GenerationLength.justification, perl = TRUE)
hab1$GenerationLength.justification <- 
  gsub("pioneer", "Pioneira", hab1$GenerationLength.justification, perl = TRUE)
hab1$GenerationLength.justification <- 
  gsub("unknown", "Desconhecido", hab1$GenerationLength.justification, perl = TRUE)

hab1$GenerationLength.justification <- 
  gsub("\\\n", "", hab1$GenerationLength.justification)
hab1$GenerationLength.justification <- 
  gsub("\\s+", " ", hab1$GenerationLength.justification)
hab1$GenerationLength.justification <- 
  gsub("^ | $", "", hab1$GenerationLength.justification)

## Any missing info?
#Any missing GL?
hab1$GenerationLength.range[is.na(hab1$GenerationLength.range)] <- 50
#Any missing p?
hab1$p.ci[is.na(hab1$p.est)] <- "0.44-0.47"
hab1$p.est[is.na(hab1$p.est)] <- 0.454


## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/habitats.csv")
head(sample, 3)
names(sample)
apply(sample, 2, unique)

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

## Manual updates
hab2$ecol.group[hab2$Name_submitted %in% "Ximenia americana"] <- "early_secondary"

## Saving the supporting files
# write.csv(standard.texts, "data/threat_standard_texts.csv", 
#           row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hab2, "data/threat_habitats.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

## Saving the SIS connect file ##
lookup <- strsplit(hab2$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup, "\\|")
lookup.name <- strsplit(hab2$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName, "\\|")

names(lookup) <- hab2$internal_taxon_id
names(lookup.name) <- hab2$internal_taxon_id

lookup1 <- unlist(lookup)
lookup.name1 <- unlist(lookup.name)
taxa_id <- rep(names(lookup), lengths(lookup))
hab3 <- data.frame(
  GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup = lookup1,
  GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName = lookup.name1,
  GeneralHabitats.GeneralHabitatsSubfield.majorImportance = NA,
  GeneralHabitats.GeneralHabitatsSubfield.season = "Resident",
  GeneralHabitats.GeneralHabitatsSubfield.suitability = "Suitable",
  internal_taxon_id = taxa_id, stringsAsFactors = FALSE
)
head(sample, 3)
head(hab3, 6)

## Adding Major importance for species with only one or multiple habitats
dup.IDs <- hab3$internal_taxon_id[duplicated(hab3$internal_taxon_id)]
rep_these <- hab3$internal_taxon_id %in% dup.IDs
hab3$GeneralHabitats.GeneralHabitatsSubfield.majorImportance[!rep_these] <- "Yes"

for (i in seq_along(unique(dup.IDs))) {
  id.i <- unique(dup.IDs)[i]
  check_these <- hab3$internal_taxon_id %in% id.i
  habs.i <- 
    hab3$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName[check_these]
  rep_these <- grepl("^Grassland|^Artificial|^Wetland", habs.i, perl = TRUE) & 
                !grepl("^Forest|^Savanna|^Rocky Areas|^Shrubland", habs.i, perl = TRUE)
  if (any(rep_these)) {
    hab3$GeneralHabitats.GeneralHabitatsSubfield.suitability[check_these][rep_these] <- "Marginal"
    hab3$GeneralHabitats.GeneralHabitatsSubfield.majorImportance[check_these][rep_these] <- "No"
  }
}

## Replacing NA by empty characters
hab4 <- hab3[!hab3$GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup %in% 
               c("NA", "", NA),]
for (i in seq_len(dim(hab4)[2])) {
  check_these <- hab4[,i] %in% c("NA", "", NA)  
  if (any(check_these))
    hab4[check_these, i] <- ""
}

## Adding common names from previous IUCN assessments
iucn.sis.habits <- iucn.sis[iucn.sis$sic.connect.file %in% "habitats", ]
hab4.1 <- dplyr::left_join(hab4, tax)
hab5 <- .merge_sis_connect(hab4, type = "habitats", tax.table = tax,
                                iucn.sis.habits, x.spp = hab4.1$species.correct2, replace = TRUE)

## Replacing NA by empty characters for the new data frame
for (i in seq_len(dim(hab5)[2])) {
  check_these <- hab5[,i] %in% c("NA", "", NA)  
  if (any(check_these))
    hab5[check_these, i] <- ""
}
hab5$validName <- NULL
hab5 <- hab5[ , c(names(hab4), "source.info")]

## Saving
write.csv(hab4, "data/sis_connect/habitats_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hab5, "data/sis_connect/habitats_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



###############################################################################H
###############################################################################H
# PLANT SPECIFIC ---------------------------------------------------------------

## Sample of the SIS CONNECT taxonomy v6.1
sample <- read.csv("SIS_sample_6_1/plantspecific.csv")
head(sample, 3)
apply(sample, 2, unique)
names(sample)

## Filtering and editing
cols <- c("internal_taxon_id",
          "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup",
          "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName")
hab6 <- hab2[, cols]
hab6 <- hab6[, match(names(sample), names(hab6), nomatch = 0)]
hab6[[1]][ hab6[[2]] %in% "Tree - large"] <- "TL"
hab6[[1]][ hab6[[2]] %in% "Tree - small"] <- "TS"
hab6[[1]][ hab6[[2]] %in% "Shrub - large"] <- "SL"
hab6[[1]][ hab6[[2]] %in% "Tree - size unknown"] <- "T"
hab6[[1]][ hab6[[2]] %in% "Succulent - tree"] <- "ST"
hab6[[1]][ hab6[[2]] %in% "Fern"] <- "PT"
head(sample)
head(hab6, 2)


## Adding common names from previous IUCN assessments
iucn.sis.plant.spec <- iucn.sis[iucn.sis$sic.connect.file %in% "plant_specific", ]
hab7 <- .merge_sis_connect(hab6, type = "plant_specific", tax.table = tax,
                           iucn.sis.plant.spec, x.spp = hab2$Name_submitted, replace = TRUE)


## Saving
write.csv(hab6, "data/sis_connect/plantspecific_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hab7, "data/sis_connect/plantspecific_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



###############################################################################H
###############################################################################H
# COUNTRY OF OCCURRENCE --------------------------------------------------------

## Generating the SIS CONNECT file "countries.csv"
### sample from SIS CONNECT v6.1
head(read.csv("SIS_sample_6_1/countries.csv"))

## Getting info from reflora and formatting
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

## Getting the country codes  (CountryOccurrenceLookup)
# ISO2 (countries)
countries <- occs[ , .N , by = c("loc.correct", "species.correct2")]
# rm(occs)
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
countries0 <- cbind.data.frame(countries, tmp[,-1])
saveRDS(countries0, "data/countries_per_species.rds")

# Building the file
countries[ , CountryOccurrenceLookup := tmp$code]
countries[ , CountryOccurrenceName := tmp$NAME_0]
countries <- countries[!is.na(CountryOccurrenceLookup), ]
countries[ , formerlyBred := NA_character_]
#Origin (originlookup: Native - 0, Reintroduced - 1, Introduced - 2, Prehistorically Introduced - 3, Vagrant - 4, Uncertain - 5)
countries[ , origin := "Native"]
#Presence (presencelookup: Extant - 0, PossiblyExtinct - 1, Extinct - 2, Presence Uncertain - 3)
countries[ , presence := "Extant"]
countries[ , seasonality := "Resident"]
countries[ , c("loc.correct", "N") := NULL]
head(countries)

# Obtaining the internal_taxon_id information
res <- tax[ , c(1,3)]
countries <- dplyr::left_join(countries, res)
names(countries)[2:7] <- paste0("CountryOccurrence.CountryOccurrenceSubfield.", 
                                names(countries)[2:7])
countries1 <- as.data.frame(unique(countries))

## Adding common names from previous IUCN assessments
iucn.sis.countries <- iucn.sis[iucn.sis$sic.connect.file %in% "countries", ]
countries2 <- .merge_sis_connect(countries1, type = "countries", tax.table = tax,
                           iucn.sis.countries, x.spp = countries$species.correct2, replace = TRUE)

## Saving
countries1$species.correct2 <- NULL
head(countries,3)
head(read.csv("SIS_sample_6_1/countries.csv"),3)
table(names(countries1) == names(read.csv("SIS_sample_6_1/countries.csv")))
write.csv(countries1, "data/sis_connect/countries_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(countries2, "data/sis_connect/countries_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



###############################################################################H
###############################################################################H
# PREVIOUS ASSESSMENTS ---------------------------------------------------------

## Red List IUCN key: 814b5b4ee7bc21205e9aec80f441b7e0ff8f957a326f2c76152ac9991b521e25
#??rredlist

## Getting the list of species in the Atlantic Forest
prev.assess <- tax
# prev.assess <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]

## Getting the IUCN assessments for Brazil (CNCFlora) - National level
tmp = flora::get.taxa(prev.assess$species.correct2, 
                      replace.synonyms = FALSE, life.form = TRUE)
tmp$search.str = prev.assess$species.correct2
tmp1 <- merge(prev.assess, tmp, 
              by.x = "species.correct2", by.y = "search.str", 
              all.x = TRUE, sort=FALSE)
#table(prev.assess$species.correct2 == tmp1$species.correct2)
prev.assess$status.reflora <- tmp1$threat.status

## Getting the global IUCN assessments (IUCN) 
iucn$species.correct2 <- sapply(strsplit(iucn$scientificName," "), 
                                function(x) paste(x[1], x[2],sep=" "))
tmp <- merge(prev.assess, iucn,  by= "species.correct2", all.x = TRUE)
tmp <- tmp[,c("internal_taxon_id","species.correct2","assessmentId","internalTaxonId",
              "redlistCategory","redlistCriteria","yearPublished","assessmentDate","criteriaVersion",
              "language","rationale","habitat","threats","population","populationTrend","range",
              "useTrade","systems","conservationActions",
              "realm","yearLastSeen","possiblyExtinct","possiblyExtinctInTheWild","scopes", 
              "AOO.range", "EOO.range", "NoThreats.noThreats", 
              "ElevationLower.limit", "ElevationUpper.limit", 
              "PopulationSize.range", "ThreatsUnknown.value", 
              "LocationsNumber.range", "GenerationLength.range", 
              "MovementPatterns.pattern", "SubpopulationNumber.range", 
              "AreaRestricted.isRestricted", "CropWildRelative.isRelative", 
              "YearOfPopulationEstimate.value", "InPlaceEducationControlled.value", 
              "SevereFragmentation.isFragmented", "InPlaceResearchRecoveryPlan.value", 
              "InPlaceLandWaterProtectionInPA.value", "InPlaceSpeciesManagementExSitu.value",
              "InPlaceResearchMonitoringScheme.value", "InPlaceEducationSubjectToPrograms.value",
              "InPlaceSpeciesManagementHarvestPlan.value", "InPlaceLandWaterProtectionAreaPlanned.value",
              "InPlaceEducationInternationalLegislation.value", "InPlaceLandWaterProtectionInvasiveControl.value",
              "InPlaceLandWaterProtectionSitesIdentified.value", "InPlaceLandWaterProtectionPercentProtected.value")]

table(prev.assess$species.correct2 == tmp$species.correct2)
prev.assess <- cbind.data.frame(prev.assess, tmp[, !names(tmp) %in% names(prev.assess)],
                                stringsAsFactors = FALSE)
prev.assess <- prev.assess[order(prev.assess$species.correct2),]


## Saving
saveRDS(prev.assess, "data/sis_connect/prev_assessments_threat.rds")



###############################################################################H
###############################################################################H
# ENDEMISM LEVELS --------------------------------------------------------------

## Species endemism levels in respect to the AF from 
end1 <- end[,c(1,2,7,9,11)]

## Getting species statuses from the AF checklist
af.list$scientific.name <- gsub("M\xfcll", "Müll", af.list$scientific.name, fixed = TRUE)
af.list1 <- plantR::fixSpecies(af.list, tax.name = "scientific.name")
fix.spp <- af.list1[!af.list1$scientificNameStatus %in% "name_w_authors", "scientificName.new"]
af.list1$scientificName.new[!af.list1$scientificNameStatus %in% "name_w_authors"] <-
  as.character(sapply(strsplit(fix.spp, " "), function(x) paste(x[1], x[2])))
names(af.list1)[9] <- "species"
end1 <- dplyr::left_join(end1, 
                         af.list1[ ,c("family", "species", "status")])

## Replacing the synonym 
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace") {
    end1$species[end1$species %in% sp.i] <- rpl.i
  }
}

##Getting mean endemism for new synonyms
tmp <- aggregate(cbind(endemism.level..validated.taxonomy.only., 
                       endemism.level..validated.and.probably.validated.taxonomy.) ~ 
                   family + species, FUN = mean, data=end1)
tmp1 <- aggregate(endemism.accepted.by.the.BF.2020 ~ family + species, 
                  FUN = function(x) paste(unique(x), collapse = "|"), data=end1)
tmp2 <- aggregate(status ~ family + species, 
                  FUN = function(x) paste(unique(x), collapse = "|"), data=end1)
tmp2$status[grepl("confirmed", tmp2$status)] <- "confirmed"
tmp2$status[grepl("new to AF", tmp2$status)] <- "new to AF"
tmp1 <- dplyr::left_join(tmp1, tmp2)
end1 <- dplyr::left_join(tmp, tmp1)

fbo <- readRDS("data/threat_fbo_tax_info.rds")[,c("species.correct2","phytogeographicDomain","endemism")]
names(fbo)[1] <- "species"
fbo$endemism <- stringr::str_to_title(fbo$endemism)
fbo <- fbo[!is.na(duplicated(fbo$species)),]
end1 <- dplyr::left_join(end1, fbo)
end1$endemism.accepted.by.the.BF.2020[!is.na(end1$phytogeographicDomain) & 
                                       !grepl("Mata Atlântica", end1$phytogeographicDomain)] <- "Not in the AF"
end1$endemism.accepted.by.the.BF.2020[!is.na(end1$phytogeographicDomain) & 
                                       end1$phytogeographicDomain %in% "Mata Atlântica" &
                                       end1$endemism %in% "Endemic"] <- "Endemic"
end1$endemism.accepted.by.the.BF.2020[!is.na(end1$phytogeographicDomain) & 
                                       grepl("Mata Atlântica", end1$phytogeographicDomain) &
                                       !end1$phytogeographicDomain %in% "Mata Atlântica"] <- "Not endemic"
end1$endemism.accepted.by.the.BF.2020[end1$species %in% "Myrcia glomerata"] <- "Not endemic"


##Classification
end1$endemic <- NA
end1$endemic[end1$endemism.level..validated.taxonomy.only. >= 85 & end1$endemism.level..validated.and.probably.validated.taxonomy. >= 85 & 
              !end1$endemism.accepted.by.the.BF.2020 %in% "Not in the AF"] <- "endemic"
end1$endemic[is.na(end1$endemic) & end1$endemism.level..validated.taxonomy.only. >= 50 & end1$endemism.level..validated.and.probably.validated.taxonomy. >= 50 &
              end1$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "endemic"
end1$endemic[end1$endemism.level..validated.taxonomy.only. < 15 & end1$endemism.level..validated.and.probably.validated.taxonomy. < 15 &
              !end1$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "occasional"
end1$endemic[is.na(end1$endemic) & end1$endemism.level..validated.taxonomy.only. >= 50 & 
              end1$endemism.level..validated.and.probably.validated.taxonomy. >= 50] <- "widespread_common"
end1$endemic[is.na(end1$endemic) & end1$endemism.level..validated.taxonomy.only. >= 15 & 
              end1$endemism.level..validated.and.probably.validated.taxonomy. >= 15] <- "widespread_sparse"
end1$endemic[is.na(end1$endemic) & 
              end1$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "widespread_common"
end1$endemic[is.na(end1$endemic) & 
              !end1$endemism.accepted.by.the.BF.2020 %in% "Endemic"] <- "widespread_sparse"
table(end1$endemic, useNA = "always")

#Probable occurrences
#end.prob$species  <- sapply(end.prob$scientific.name, flora::remove.authors)
end.prob$species  <- sapply(end.prob$scientific.name, function(x) paste(strsplit(x, " ")[[1]][1:2], collapse = " "))
end.prob1 <- end.prob[end.prob$status %in% "no records found but cited in BF.2020",]
end.prob1 <- end.prob1[!end.prob1$species %in% end1$species,]

tmp <- flora::get.taxa(end.prob1$species)
tmp1 <- flora::get_domains(tmp)
tmp2 <- flora::get_endemism(tmp)
end.prob1$domain <- tmp1$domain 
end.prob1$end.br <- tmp2$endemism 
end.prob1$endemic <- NA
end.prob1$endemic[end.prob1$domain %in% "Mata Atlântica" & 
                   end.prob1$end.br %in%  "Endemica"] <- "endemic"
end.prob1$endemic[is.na(end.prob1$endemic) & grepl("Mata Atlântica", end.prob1$domain)] <- 
  "not_endemic"
end.prob1$endemic[is.na(end.prob1$endemic) & !grepl("Mata Atlântica", end.prob1$domain)] <- 
  "not in the AF"
table(end.prob1$endemic, useNA = "always")

#Merging the two data.frames
end2 <- end1[,c("family","species",
               "endemism.level..validated.and.probably.validated.taxonomy.",
               "endemism.level..validated.taxonomy.only.",
               "endemic",
               "status")]
names(end2)[3:4] <- c("endemism.level.1","endemism.level.2") 

end.prob1$endemism.level.1 <- NA
end.prob1$endemism.level.2 <- NA
end.prob2 <- end.prob1[,c("family","species",
                         "endemism.level.1",
                         "endemism.level.2",
                         "endemic", "status")]
end.final <- rbind.data.frame(end2, end.prob2,
                              stringsAsFactors = FALSE)
end.final <- end.final[order(end.final$species),]

##Getting the endemism in Brazil (no CNCFlora 'jurisidction')
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

## Getting the endemism from other AF countries
end.final1 <- dplyr::left_join(end.final, 
                               all.crit[,c("species", "endemism.ARG")])

## Saving
saveRDS(end.final1, "data/sis_connect/endemism_threat.rds")



###############################################################################H
###############################################################################H
# CITES LIST -------------------------------------------------------------------

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
all_cites_spp$code_annotation <- gsub("A$|L$", "", all_cites_spp$code_annotation, perl = TRUE)
all_cites_spp$hash_annotation <- 
  gsub("^#10|^#15|^#4|^#5|^#6", "", all_cites_spp$hash_annotation, perl = TRUE)
# all_cites_spp$hash_annotation <- 
#   stringr::str_squish(gsub("\\\n|\\\r", " ", all_cites_spp$hash_annotation))
all_cites_spp$hash_annotation <- 
  gsub("\\s+", " ", gsub("\\\n|\\\r", " ", all_cites_spp$hash_annotation, perl = TRUE), perl = TRUE)
all_cites_spp$hash_annotation <- 
  gsub("^ | $", "", all_cites_spp$hash_annotation, perl = TRUE)

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
                                 c("id","start_date","is_current",
                                   "eu_decision_type.description","term.code","term.name")))]
decis.EU1 <- aggregate(. ~ taxon_concept_id, data = decis.EU, 
                       function(x) paste(x, collapse = "|")) 
all_cites.EU_spp1 <- merge(all_cites.EU_spp,  decis.EU1,
                           by.x = "id", by.y = "taxon_concept_id",
                           all.x = TRUE, sort = FALSE)
## Saving 
all_cites.EU_spp1 <- 
  all_cites.EU_spp1[order(all_cites.EU_spp1$species.correct2),]
saveRDS(all_cites.EU_spp1, "data/threat_cites_EU.rds")



###############################################################################H
###############################################################################H
# USE AND TRADE ----------------------------------------------------------------

## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/usetrade.csv")
head(sample, 3)

##Loading the species list
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]

##Loading and Preparing TreeCo use database
usos1 <- uses_data_prep(usos, spp = tax$species.correct2)

##Crossing TreeCo and IUCN use codes and nomenclature
lookup.codes <- read.csv("SIS_sample_6_1/Utenduse.utenduselookup.csv", 
                         encoding = "UTF-8")
usos1$UTEndUseLookup <- NA
usos1$UTEndUseLookup <- 
  lookup.codes$code[match(usos1$uses, lookup.codes$treeco)]

##Creating the other required columns
usos1$UTEndUseName <- NA
usos1$UTEndUseName <- 
  lookup.codes$description[match(usos1$uses, lookup.codes$treeco)]

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

##Filtering and saving the species with some type of use
usos1 <- usos1[!is.na(usos1$UTEndUseLookup),]
write.csv(usos1, "data/threat_species_uses.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

##Filtering the selected uses for the SIS CONNECT file
usos2 <- usos1[ , c(names(usos1)[which(names(usos1) == "UTEndUseLookup"): dim(usos1)[2]], "source")]
names(usos2)[1:6] <- paste0("UTEndUse.UTEndUseSubfield.", 
                            names(usos2)[1:6])
head(sample, 3)
head(usos2, 3)

## Replacing NA by empty characters
for (i in seq_len(dim(usos2)[2])) {
  check_these <- is.na(usos2[,i])
  if (any(check_these))
    usos2[check_these, i] <- ""
}

## Adding common names from previous IUCN assessments
iucn.sis.usetrade <- iucn.sis[iucn.sis$sic.connect.file %in% "usetrade", ]
usos3 <- .merge_sis_connect(usos2, type = "usetrade", tax.table = tax,
                                 iucn.sis.usetrade, x.spp = usos1$Name_submitted, replace = TRUE)

## Replacing NA by empty characters for the new data frame
for (i in seq_len(dim(usos3)[2])) {
  check_these <- is.na(usos3[,i])
  if (any(check_these))
    usos3[check_these, i] <- ""
}

##Saving the SIS Connet file
write.csv(usos2, "data/sis_connect/usetrade_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(usos3, "data/sis_connect/usetrade_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

##Saving the file with the exploited species
usos4 <- usos1[usos1$uses %in% c("timber") |
                 usos1$Name_submitted %in% c("Euterpe edulis") , 
               names(usos1)[1:which(names(usos1) == "uses")]]

check_these <- !usos1$Name_submitted %in% usos4$Name_submitted
cols <- c("essential_oils", "fences", "fibres","fodder","food","resins_gums","rubber_latex", "fuel")
tab <- table(usos1$Name_submitted[check_these], usos1$uses[check_these])
linhas <- apply(tab[,cols], 1, sum) > 0 
tab[linhas, cols] # ok in 30/04/2021

#How many times were they cited as being used for timber
usos4 <- usos4[!grepl("Coradin", usos3$source) , ] # Aqui Euterpe edulis está sendo removido. ok?
usos4$obs <- gsub("list provided by GTA. Note: I haven’t updated it for taxonomy so there might be some issues", NA, usos4$obs)
usos5 <- as.data.frame(table(usos4$Name_submitted))
names(usos5)[1:2] <- c("species.correct2", "times.cites")
usos5$sources <- aggregate(usos4$source, list(usos4$Name_submitted), 
                           function (x) paste(x, collapse = " | "))$x
usos5$obs <- aggregate(usos4$obs, list(usos4$Name_submitted), 
                       function (x) paste(x, collapse = " | "))$x
write.csv(usos5, "data/threat_exploited_timber_spp.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")



###############################################################################H
###############################################################################H
# THREATS ----------------------------------------------------------------------

## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/threats.csv")
head(sample, 3)
apply(sample, 2, unique)

#Reading uses data
usos1 <- read.csv("data/threat_species_uses.csv", as.is = TRUE)


## CREATING THE TYPES OF THREATS FOR THE ATALNTIC FOREST ##
## Residential & commercial development
# 1.1 Housing & urban areas (generalizado para a Mata Atlântica)
res1.1 <- c(stress = "1.1|1.3|2.1|2.2", 
            stressdesc = "Ecosystem Conversion|Indirect Ecosystem Effects|Species Mortality|Species Disturbance", 
            ThreatsLookup = "1.1", ThreatsName = "Housing & Urban Areas", ancestry = NA, ias = NA,
            internationalTrade = NA, scope = "Majority (50-90%)", 
            severity = "Slow, Significant Declines", text = NA, 
            timing = "Ongoing", virus = NA)
# 1.2 Commercial & industrial areas  (e.g. proximo de grandes centros urbanos)
res1.2 <- c(stress = "1.1|1.3|2.1|2.2", 
            stressdesc = "Ecosystem Conversion|Indirect Ecosystem Effects|Species Mortality|Species Disturbance", 
            ThreatsLookup = "1.2", ThreatsName = "Commercial & Industrial Areas", ancestry = NA, ias = NA,
            internationalTrade = NA, scope = "Minority (<50%)", 
            severity = "Slow, Significant Declines", text = NA, 
            timing = "Ongoing", virus = NA)
# 1.3 Tourism & recreation areas (e.g. proximo ao litoral)
res1.3 <- c(stress = "1.2|2.2", 
            stressdesc = "Ecosystem Degradation|Species Disturbance", 
            ThreatsLookup = "1.3", ThreatsName = "Tourism & Recreation Areas", ancestry = NA, ias = NA,
            internationalTrade = NA, scope = "Minority (<50%)", 
            severity = "Causing/Could cause fluctuations", text = NA, 
            timing = "Ongoing", virus = NA)


## Agriculture & aquaculture
# 2.1	Annual & perennial non-timber crops
#2.1.2	Small-holder farming (e.g. Araucaria)
agr2.1.2 <- c(stress = "1.1|1.2|2.1|2.2|2.3", 
              stressdesc = "Ecosystem Conversion|Ecosystem Degradation|Species Mortality|Species Disturbance|Indirect Species Effects", 
              ThreatsLookup = "2.1.2", ThreatsName = "Small-holder Farming", ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Majority (50-90%)", 
              severity = "Slow, Significant Declines", text = NA, 
              timing = "Ongoing", virus = NA) # Can be "Past, Unlikely to Return" for some species
#2.1.3	Agro-industry farming (e.g. Araucaria, Ilex paraguariensis, Euterpe edulis)
agr2.1.3 <- c(stress = "1.1|1.2|2.1|2.2|2.3", 
              stressdesc = "Ecosystem Conversion|Ecosystem Degradation|Species Mortality|Species Disturbance|Indirect Species Effects", 
              ThreatsLookup = "2.1.3", ThreatsName = "Agro-industry Farming", ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Majority (50-90%)", 
              severity = "Rapid Declines", text = NA, 
              timing = "Ongoing", virus = NA)
#2.1.4  Scale Unknown/Unrecorded (mais usado que os anteriored na AmSul)
agr2.1.4 <- c(stress = "1.1|1.2|2.1|2.2|2.3", 
              stressdesc = "Ecosystem Conversion|Ecosystem Degradation|Species Mortality|Species Disturbance|Indirect Species Effects", 
              ThreatsLookup = "2.1.4", ThreatsName = "Scale Unknown/Unrecorded", ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Majority (50-90%)", 
              severity = "Rapid Declines", text = NA, 
              timing = "Ongoing", virus = NA)

# 2.2	Wood & pulp plantations
#2.2.2	Agro-industry plantations (e.g. região de Aracruz)
agr2.2.2 <- c(stress = "1.1|2.1|2.3", 
              stressdesc = "Ecosystem Conversion|Species Mortality|Indirect Species Effects", 
              ThreatsLookup = "2.2.2", 
              ThreatsName = "Agro-industry Plantations", ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Minority (<50%)", 
              severity = "Causing/Could cause fluctuations", text = NA, 
              timing = "Ongoing", virus = NA)

#2.2.3	Scale Unknown/Unrecorded
agr2.2.3 <- c(stress = "1.1|2.1|2.3", 
              stressdesc = "Ecosystem Conversion|Species Mortality|Indirect Species Effects", 
              ThreatsLookup = "2.2.3", ThreatsName = "Scale Unknown/Unrecorded", ancestry = NA, ias = NA,
              internationalTrade = NA, scope = "Minority (<50%)", 
              severity = "Causing/Could cause fluctuations", text = NA, 
              timing = "Ongoing", virus = NA)


# 2.3 Livestock farming & ranching
#2.3.3	Agro-industry grazing, ranching or farming
live2.3.3 <- c(stress = "1.1|1.2|2.1|2.2|2.3", 
               stressdesc = "Ecosystem Conversion|Ecosystem Degradation|Species Mortality|Species Disturbance|Indirect Species Effects", 
               ThreatsLookup = "2.3.3", ThreatsName = "Agro-industry Grazing, Ranching or Farming", ancestry = NA, ias = NA,
               internationalTrade = NA, scope = "Majority (50-90%)", 
               severity = "Rapid Declines", text = NA, 
               timing = "Ongoing", virus = NA)
#2.3.4	Scale Unknown/Unrecorded (mais seguro/usado que o anterior)
live2.3.4 <- c(stress = "1.1|1.2|2.1|2.2|2.3", 
               stressdesc = "Ecosystem Conversion|Ecosystem Degradation|Species Mortality|Species Disturbance|Indirect Species Effects", 
               ThreatsLookup = "2.3.4", ThreatsName = "Scale Unknown/Unrecorded", ancestry = NA, ias = NA,
               internationalTrade = NA, scope = "Majority (50-90%)", 
               severity = "Rapid Declines", text = NA, 
               timing = "Ongoing", virus = NA)


## 3	Energy production & mining
# 3.2	Mining & quarrying
mine3.2 <- c(stress = "1.1|1.3|2.1|2.2", 
             stressdesc = "Ecosystem Conversion|Indirect Ecosystem Effects|Species Mortality|Species Disturbance", 
             ThreatsLookup = "3.2", ThreatsName = "Mining & Quarrying", ancestry = NA, ias = NA,
             internationalTrade = NA, scope = "Minority (<50%)", 
             severity = "Causing/Could cause fluctuations", text = NA, 
             timing = "Ongoing", virus = NA)


## 4	Transportation & service corridors
# 4.1	Roads & railroads
road4.1 <- c(stress = "1.1|1.2", 
             stressdesc = "Ecosystem Conversion|Ecosystem Degradation", 
             ThreatsLookup = "3.2", ThreatsName = "Roads & Railroads", ancestry = NA, ias = NA,
             internationalTrade = NA, scope = "Minority (<50%)", 
             severity = "Causing/Could cause fluctuations", text = NA, 
             timing = "Ongoing", virus = NA) 


## 5	Biological resource use
# 5.2	Gathering terrestrial plants
#5.2.1	Intentional use (species is the target) e.g. Araucaria
bio5.2.1 <- c(stress = "2.2|2.3", 
              stressdesc = "Species Disturbance|Indirect Species Effects", 
              ThreatsLookup = "5.2.1", 
              ThreatsName = "Intentional Use (species being assessed is the target)", 
              ancestry = NA, ias = NA, 
              internationalTrade = "No", scope = "Unknown", 
              severity = "Causing/Could cause fluctuations", text = NA, 
              timing = "Ongoing", virus = NA)

# 5.3	Logging & wood harvesting
#5.3.1	Intentional use: (subsistence/small scale) [harvest]
wood5.3.1 <- c(stress = "2.1|2.3", 
               stressdesc = "Species Mortality|Indirect Species Effects", 
               ThreatsLookup = "5.3.1", 
               ThreatsName = "Intentional Use: subsistence/small scale (species being assessed is the target) [harvest]", ancestry = NA, ias = NA, 
               internationalTrade = "No", scope = "Minority (<50%)", 
               severity = "Causing/Could cause fluctuations", text = NA, 
               timing = "Ongoing", virus = NA)
#5.3.2	Intentional use: (large scale) [harvest]
wood5.3.2 <- c(stress = "2.1|2.3", 
               stressdesc = "Species Mortality|Indirect Species Effects", 
               ThreatsLookup = "5.3.2", 
               ThreatsName = "Intentional Use: large scale (species being assessed is the target) [harvest]", 
               ancestry = NA, ias = NA, 
               internationalTrade = "No", scope = "Majority (50-90%)", 
               severity = "Rapid Declines", text = NA, 
               timing = "Ongoing", virus = NA)
#5.3.5	Motivation Unknown/Unrecorded (mais seguro/usado que o anterior)
wood5.3.5 <- c(stress = "2.1|2.3", 
               stressdesc = "Species Mortality|Indirect Species Effects", 
               ThreatsLookup = "5.3.5",  ThreatsName = "Motivation Unknown/Unrecorded", 
               ancestry = NA, ias = NA, 
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
invade8.1.1 <- c(stress = "2.3.2", stressdesc = "Competition", ThreatsLookup = "8.1.1", 
                 ThreatsName = "Unspecified Species", ancestry = NA, ias = NA, 
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
wood.spp.actual <- wood.spp[wood.spp$Name_submitted %in% usos5$species.correct2, ]
wood.spp.potential <- wood.spp[!wood.spp$Name_submitted %in% usos5$species.correct2, ]
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

## Replacing NA by empty characters
for (i in seq_len(dim(threats)[2])) {
  check_these <- is.na(threats[,i])
  if (any(check_these))
    threats[check_these, i] <- ""
}

## Adding common names from previous IUCN assessments
iucn.sis.threats <- iucn.sis[iucn.sis$sic.connect.file %in% "threats", ]
threats.1 <- dplyr::left_join(threats, tax)
threats1 <- .merge_sis_connect(threats, type = "threats", tax.table = tax,
                            iucn.sis.threats, x.spp = threats.1$species.correct2, replace = TRUE)

## Replacing NA by empty characters for the new data frame
for (i in seq_len(dim(threats1)[2])) {
  check_these <- is.na(threats1[,i])
  if (any(check_these))
    threats1[check_these, i] <- ""
}

## Saving the SIS Connet file
write.csv(threats, "data/sis_connect/threats_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(threats1, "data/sis_connect/threats_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


###############################################################################H
###############################################################################H
# CONSERVATION NEEDED ----------------------------------------------------------

## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/conservationneeded.csv")
# head(sample, 3)
# apply(sample, 2, unique)

## CREATING THE TYPES OF CONSERVATION NEEDS FOR THE ATLANTIC FOREST ##
# 1.1 Site/area protection (generalizado para a Mata Atlântica)
con1.1 <- c(ConservationActionsLookup = "1.1", ConservationActionsName = "Site/area Protection", 
            note = "Criar ou expandir unidades de conservação (Categorias IUCN I-VI) dentro da área de distribuição da espécie")
# 1.2	Resource & habitat protection
con1.2 <- c(ConservationActionsLookup = "1.2", ConservationActionsName = "Resource & Habitat Protection", 
            note = "Proteção dos recursos naturais aos quais a espécie está associada, mesmo que fora de Unidades de Conservação")
# 2.3	Habitat & natural process restoration
con2.3 <- c(ConservationActionsLookup = "2.3", ConservationActionsName = "Habitat & Natural Process Restoration", 
            note = "Recuperar a degradação ambiental e/ou restaurar os habitats, corredores ecológicos e funções do ecossistema aos quais a espécie está associada")
# 3.1.1	Harvest management
con3.1.1 <- c(ConservationActionsLookup = "3.1.1", ConservationActionsName = "Harvest Management", 
              note = "Regulamentar o volume e frequência da colheita de produtos florestais não-madeireiros associados à espécie")
# 3.1.2	Trade management
con3.1.2 <- c(ConservationActionsLookup = "3.1.2", ConservationActionsName = "Trade Management", 
              note = "Regulamentar o comércio de produtos florestais não-madeireiros associados à espécie")
# 3.3.1	Reintroduction
con3.3.1 <- c(ConservationActionsLookup = "3.3.1", ConservationActionsName = "Species Re-Introduction", 
              note = "Reintrodução (i.e. plantio) de indivíduos em áreas onde a espécie foi extirpada, seguindo as orientações da IUCN")
# 3.4	Ex-situ conservation
con3.4 <- c(ConservationActionsLookup = "3.4", ConservationActionsName = "Ex-situ Conservation", 
            note = "Proteção ex-situ da espécie via plantios em jardins botânicos e incorporação em bancos de germoplasma e/ou bases de dados genéticos (e.g. GenBank)")
# 4.2	Training
con4.2 <- c(ConservationActionsLookup = "4.2", ConservationActionsName = "Training", 
            note = "Capacitação de profissionais que possam atuar diretamente em ações de conservação da espécie (e.g. gestores de unidades de conservação, agentes de fiscalização, produtores de semente), aprimorando conhecimentos, habilidades e troca de informações entre esses profissionais (e.g. melhorando habilidades de identificação da espécie)")
# 4.3	Awareness & Communications
con4.3 <- c(ConservationActionsLookup = "4.3", ConservationActionsName = "Awareness & Communications", 
            note = "Sensibilização ambiental da sociedade em geral sobre o grau de ameaça da espécie via diferentes meios de comunicação")
# 5.1.1	International Level
con5.1.1 <- c(ConservationActionsLookup = "5.1.1", ConservationActionsName = "International Level", 
              note = "Criar, implementar ou adaptar dispositivos inter-nacionais que facilitem a proteção da espécie")
# 5.1.2	National Level
con5.1.2 <- c(ConservationActionsLookup = "5.1.2", ConservationActionsName = "National Level", 
              note = "Criar, implementar ou adaptar leis nacionais que oficializem a proteção da espécie")
# 5.2	Policies & Regulations
con5.2 <- c(ConservationActionsLookup = "5.2", ConservationActionsName = "Policies & Regulations", 
            note = "Criar, implementar ou das suporte a políticas públicas que facilitem a proteção da espécie")
# 5.4.1	International Level
con5.4.1 <- c(ConservationActionsLookup = "5.4.1", ConservationActionsName = "International level", 
              note = "Reforçar e garantir a implementação de tratados e regulamentações internacionais que envolvam a espécie (e.g. CITES), incluindo a implementação das sanções e punições legais cabíveis")
# 6.2	Substitution
con6.2 <- c(ConservationActionsLookup = "6.2", ConservationActionsName = "Substitution", 
            note = "Promover o uso de produtos oriundos de espécies manejadas e que causem menores pressões sobre as populações naturais da espécie")
# 6.3	Market Forces
con6.3 <- c(ConservationActionsLookup = "6.3", ConservationActionsName = "Market Forces", 
            note = "Promover o uso de mecanismos de mercado que promovam mudanças de comportamento em relação ao uso da espécie (e.g. certificação florestal)")
# 7.1	Institutional & Civil Society Development
con7.1 <- c(ConservationActionsLookup = "7.1", ConservationActionsName = "Institutional & Civil Society Development", 
            note = "Dar suporte e capacitar organizações, empresas e comunidades que possam atuar diretamente em ações para a conservação da espécie")
# 7.2	Alliance & Partnership Development
con7.2 <- c(ConservationActionsLookup = "7.2", ConservationActionsName = "Alliance & Partnership Development", 
            note = "Formar redes de colaboração multi-setoriais para divulgar e promover a conservação da espécie")

## CREATING COMBINED INFO FOR ALL SPECIES ##
# tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
# all.spp <- merge(all.crit, tax,
#                    by.x = "species", by.y = "species.correct2", all.x = TRUE, sort = FALSE)
# cites <- readRDS("data/threat_cites_EU.rds")
# all.spp <- merge(all.spp, cites[ ,c("internal_taxon_id", "cites_listing", "appendix", "annotation")],
#                  by.x = "internal_taxon_id", 
#                  by.y = "internal_taxon_id", all.x = TRUE, sort = FALSE)
# all.spp <- all.spp[order(as.double(gsub("sp", "", all.spp$internal_taxon_id))), ]
# saveRDS(all.spp, "data/all_spp_crits_tax_cites.rds")
n.spp <- dim(all.spp)[1]

#species uses
used.spp1 <- read.csv("data/threat_species_uses.csv")
used.spp2 <- read.csv("data/sis_connect/usetrade_threat.csv")
used.spp4 <- read.csv("data/threat_exploited_timber_spp.csv")

## 1 - LAND PROTECTION ##
threat.cats <- c("EN","VU","CR")
#1.1 Site/Area Protection
# species: threatened, endemic and with low % of occs (<10%) and EOO (<5%) in Protected Areas
# Threatened endemics with low occurrence in or not within conservation units
unprotected.spp <- unique( 
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" & 
              ((all.spp$prop.EOO.in.StrictUCs < 5 & 
                all.spp$protected < 10) | 
                  (all.spp$prop.EOO.in.StrictUCs < 1 | all.spp$protected < 1)),   
          c("internal_taxon_id", "species")])
unprotected <- do.call(rbind.data.frame,
                    replicate(dim(unprotected.spp)[1], con1.1, simplify = FALSE))
names(unprotected) <- names(con1.1)
unprotected$internal_taxon_id <- unprotected.spp$internal_taxon_id

## 2 - LAND MANAGEMENT ##
# 2.3	Habitat & natural process restoration 
# species: threatened, endemic and with low AOH (<20%) and less than EOO (<50%) in Protected Areas
restore.spp <- unique( 
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" & 
              all.spp$prop.EOO.forest < 20 &
                all.spp$protected < 50,   
          c("internal_taxon_id", "species")])
restore <- do.call(rbind.data.frame,
                       replicate(dim(restore.spp)[1], con2.3, simplify = FALSE))
names(restore) <- names(con2.3)
restore$internal_taxon_id <- restore.spp$internal_taxon_id

## 3 - SPECIES MANAGEMENT ##
# 3.1.1	Harvest management
# uses: "timber", "fuel", "fences"
# uses classes: 7, 9
# species: Threatened or not, endemic or common, with known timber-related uses
wood.sp <- used.spp1[used.spp1$UTEndUseLookup %in% c(7, 9),
            c("Name_submitted", "lifeform", "origin",
              "uses", "UTEndUseLookup", "subsistence","internal_taxon_id", "source")]
#Filtering only known timber species (sources including both actual and potential uses)
wood.sp.actual <- wood.sp[!grepl("Coradin|Agroforestree database", wood.sp$source),]
wood.spp <- unique( 
  all.spp[all.spp$internal_taxon_id %in% unique(wood.sp$internal_taxon_id) & 
           all.spp$cat.reg.clean %in% c(threat.cats, "NT") &
            all.spp$endemic %in% c("endemic", "widespread_common"),
            c("internal_taxon_id", "species")])
harvest <- do.call(rbind.data.frame,
                   replicate(dim(wood.spp)[1], con3.1.1, simplify = FALSE))
names(harvest) <- names(con3.1.1)
harvest$internal_taxon_id <- wood.spp$internal_taxon_id

# 3.1.2	Trade management
# Species intentional uses: "food/alcohol", "fodder", "resins_gums/rubber_latex", "essential_oils", "fibres"
# uses classes: 1, 2, 5, 6, 8
# species: Threatened or not, endemic or common, with known non-timber uses
used.sp <- unique(
  used.spp1[used.spp1$UTEndUseLookup %in% c(1,2,5,6,8),
            c("Name_submitted", "lifeform", "origin",
              "uses", "UTEndUseLookup", "subsistence","internal_taxon_id")])
#Filtering only known non-timber species
used.spp <- unique( 
  all.spp[all.spp$internal_taxon_id %in% used.sp$internal_taxon_id & 
            all.spp$cat.reg.clean %in% c(threat.cats, "NT") &
              all.spp$endemic %in% c("endemic", "widespread_common"), 
          c("internal_taxon_id", "species")])
used <- do.call(rbind.data.frame,
                   replicate(dim(used.spp)[1], con3.1.2, simplify = FALSE))
names(used) <- names(con3.1.2)
used$internal_taxon_id <- used.spp$internal_taxon_id


# 3.3.1	Reintroduction
# Threatened timber species or over-exploited species
wood.use.sp <- used.spp1[used.spp1$internal_taxon_id %in% wood.uses$internal_taxon_id,
                     c("Name_submitted", "lifeform", "origin",
                       "uses", "UTEndUseLookup", "subsistence","internal_taxon_id", "source")]
exploited.sp <- used.spp1[used.spp1$Name_submitted %in% used.spp4$species.correct2,
                         c("Name_submitted", "lifeform", "origin",
                           "uses", "UTEndUseLookup", "subsistence","internal_taxon_id", "source")]
proibidas.sp <- used.spp1[used.spp1$Name_submitted %in% madeiras.proibidas,
                          c("Name_submitted", "lifeform", "origin",
                            "uses", "UTEndUseLookup", "subsistence","internal_taxon_id", "source")]

all.explo.sp <- rbind.data.frame(wood.use.sp, exploited.sp, proibidas.sp)
all.explo.sp <- all.explo.sp[!duplicated(all.explo.sp$internal_taxon_id), ]
  
#Filtering over-exploited species
explo.spp <- unique( 
  all.spp[all.spp$internal_taxon_id %in% all.explo.sp$internal_taxon_id & 
            all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% c("endemic", "widespread_common"), 
          c("internal_taxon_id", "species")])
exploited <- do.call(rbind.data.frame,
                replicate(dim(all.explo.sp)[1], con3.3.1, simplify = FALSE))
names(exploited) <- names(con3.3.1)
exploited$internal_taxon_id <- all.explo.sp$internal_taxon_id


# 3.4	Ex-situ conservation
# Similar to 1 - Land protection
# species: threatened, endemic and with low % of occs (<20%) and EOO (<10%) in Protected Areas
# All threatened endemics with not many occurrence in conservation units
ex.situ.spp <- unique( 
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" & 
            ((all.spp$prop.EOO.in.StrictUCs < 10 & 
                all.spp$protected < 20) | 
               (all.spp$prop.EOO.in.StrictUCs < 5 | all.spp$protected < 5)),   
          c("internal_taxon_id", "species")])
ex.situ <- do.call(rbind.data.frame,
                       replicate(dim(ex.situ.spp)[1], con3.4, simplify = FALSE))
names(ex.situ) <- names(con3.4)
ex.situ$internal_taxon_id <- ex.situ.spp$internal_taxon_id


## 4 - EDUCATION & AWARNESS ##
# species: all CR endemics
CR.end.spp <- unique( 
  all.spp[all.spp$cat.reg.clean %in% "CR" &
            all.spp$endemic %in% "endemic",   
          c("internal_taxon_id", "species")])
# 4.2	Training
training <- do.call(rbind.data.frame,
                   replicate(dim(CR.end.spp)[1], con4.2, simplify = FALSE))
names(training) <- names(con4.2)
training$internal_taxon_id <- CR.end.spp$internal_taxon_id

# 4.3	Awareness & communications
aware <- do.call(rbind.data.frame,
                    replicate(dim(CR.end.spp)[1], con4.3, simplify = FALSE))
names(aware) <- names(con4.3)
aware$internal_taxon_id <- CR.end.spp$internal_taxon_id


## 5 - LAW & POLICY ##
#5.1.1	International level	depend (if threatened and not endemic)
international.spp <- unique(
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" &
            !all.spp$endemismo.reflora %in% "endemic",   
          c("internal_taxon_id", "species")])
international <- do.call(rbind.data.frame,
                    replicate(dim(international.spp)[1], con5.1.1, simplify = FALSE))
names(international) <- names(con5.1.1)
international$internal_taxon_id <- international.spp$internal_taxon_id 


#5.1.2	National level depend (if threatened and Brazilian endemic)
national.spp <- unique(
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" &
            all.spp$endemismo.reflora %in% "endemic",   
          c("internal_taxon_id", "species")])
national <- do.call(rbind.data.frame,
                    replicate(dim(national.spp)[1], con5.1.2, simplify = FALSE))
names(national) <- names(con5.1.2)
national$internal_taxon_id <- national.spp$internal_taxon_id 

## INTERNATIONAL TRADE ##
#5.4.1 International level (if CITES)
cites.spp <- unique(
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" &
            all.spp$cites_listing %in% c("I","II","III"),   
          c("internal_taxon_id", "species", "family",
            "cites_listing", "appendix", "annotation")])
cites.spp <- cites.spp[!grepl("Senna|Pavonia", cites.spp$species), ]
cites.spp <- cites.spp[!is.na(cites.spp$appendix), ]
inter.trade.spp <- unique(
  all.spp[all.spp$cat.reg.clean %in% threat.cats &
            all.spp$endemic %in% "endemic" &
            all.spp$species %in% cites.spp$species,
          c("internal_taxon_id", "species")])
international.trade <- do.call(rbind.data.frame,
                         replicate(dim(inter.trade.spp)[1], con5.4.1, simplify = FALSE))
names(international.trade) <- names(con5.4.1)
international.trade$internal_taxon_id <- inter.trade.spp$internal_taxon_id 


## 6 - LIVELIHOOD, ECONOMIC & OTHER INCENTIVES ##
# Euterpe edulis, Araucaria angustifolia?
# 6.1 Linked Enterprises & Livelihood Alternatives
# 6.2 Substitution 
# 6.3 Market forces 
# 6.4 Conservation Payments
CR.EN.explo.spp <- unique( 
  all.spp[all.spp$cat.reg.clean %in% c("EN","CR") &
            all.spp$endemic %in% c("endemic", "widespread_common") &
              (all.spp$internal_taxon_id %in% all.explo.sp$internal_taxon_id | 
                 all.spp$internal_taxon_id %in% "sp1780"),   
          c("internal_taxon_id", "species")])
substitution <- do.call(rbind.data.frame,
                    replicate(dim(CR.EN.explo.spp)[1], con6.2, simplify = FALSE))
names(substitution) <- names(con6.2)
substitution$internal_taxon_id <- CR.EN.explo.spp$internal_taxon_id

market <- do.call(rbind.data.frame,
                        replicate(dim(CR.EN.explo.spp)[1], con6.3, simplify = FALSE))
names(market) <- names(con6.3)
market$internal_taxon_id <- CR.EN.explo.spp$internal_taxon_id


## SPECIES ASSESSED AS CR or CR (PE) ##
CR.PE.spp <- all.spp[all.spp$cat.reg.clean %in% "CR" &
              all.spp$endemic %in% c("endemic"), c("internal_taxon_id","species")]
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
              replicate(dim(CR.PE.spp)[1], CR.PE, simplify = FALSE))
CR$internal_taxon_id <- rep(CR.PE.spp$internal_taxon_id, each = dim(CR.PE)[1])


## COMBINING ALL THREATS 
conservation <- rbind.data.frame(unprotected, restore, harvest, used, 
                                 exploited, ex.situ, training, aware,
                                 international, national, international.trade,
                                 substitution, market, CR, 
                            stringsAsFactors = FALSE)
#Re-ordering
conservation <- conservation[order(as.double(gsub("sp","", 
                                    conservation$internal_taxon_id)), 
                                   conservation$ConservationActionsLookup), ] 
#Removig possible duplicates
dups <- duplicated(paste(conservation$internal_taxon_id,  
                         conservation$ConservationActionsLookup, sep = "_"))
conservation <- conservation[!dups, ]

#Renaming
head(sample, 3)
head(conservation, 3)
names(conservation)[1:3] <- paste0("ConservationActions.ConservationActionsSubfield.", 
                                   names(conservation)[1:3])
table(names(conservation) == names(sample))

## Adding conservation info from previous IUCN assessments
iucn.sis.conservation <- iucn.sis[iucn.sis$sic.connect.file %in% "conservation_needed", ]
conservation.1 <- dplyr::left_join(conservation, tax)
conservation1 <- .merge_sis_connect(conservation, type = "conservation_needed", tax.table = tax,
                               iucn.sis.conservation, x.spp = conservation.1$species.correct2, replace = TRUE)

#Adding descriptions to the new categories from the previous assessments
new.cats <- c("2.1", "2.2", "3.2", "3.4.1", "3.4.2", "5.2", "5.4.2", "5.4.3")
names(new.cats) <- c("Implementação de atividades de manejo para fins de conservação em unidades de conservação e outras áreas protegidas onde a espécie ocorre naturalmente",
                     "Controle e/ou prevenção de espécies invasoras, patógenos e outras espécies problemas que afetem a sobrevivência dos indivíduos da espécie",
                     "Manejo, melhoramento e recuperação de populações da espécie (e.g. polinização manual, plantios de enriquecimento",
                     "Criação em cativeiro (e.g. plantios em jardins botânicos ou arboretos), propagação vegetativa e/ou produção de sementes e mudas",
                     "Incorporação em bancos de germoplasma e/ou bases de dados genéticos (e.g. GenBank)",
                     "Criar, implementar ou das suporte a políticas públicas que facilitem a proteção da espécie",
                     "Reforçar e garantir a implementação de tratados e regulamentações nacionais que envolvam a espécie (e.g. CITES), incluindo a implementação das sanções e punições legais cabíveis",
                     "Reforçar e garantir a implementação de tratados e regulamentações estaduais e municipais que envolvam a espécie (e.g. CITES), incluindo a implementação das sanções e punições legais cabíveis")

for (i in seq_along(new.cats)) {
  conservation1$ConservationActions.ConservationActionsSubfield.note[conservation1$ConservationActions.ConservationActionsSubfield.ConservationActionsLookup %in% 
                  new.cats[i]] <- names(new.cats)[i]
}


## Saving the SIS Connet file
write.csv(conservation, "data/sis_connect/conservationneeded_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(conservation1, "data/sis_connect/conservationneeded_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


###############################################################################H
###############################################################################H
# RESEARCH NEEDED --------------------------------------------------------------

## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/researchneeded.csv")
head(sample, 3)
apply(sample, 2, unique)


## CREATING THE TYPES OF RESEARCH NEEDS FOR THE ATLANTIC FOREST ##
# 1	Research					
#1.1	Taxonomy depend (if most records are old)		
res1.1 <- c(ResearchLookup = "1.1", ResearchName = "Taxonomy", 
            note = "Espécie conhecida apenas para a localidade onde o espécimen tipo foi coletado e/ou que não tenha sido re-coletado nos últimos 50 anos. Pode se tratar de uma sinônimo de outra espécie que não tenha recebido o tratamento taxonômico necessário. Caso contrário, pode se tratar de uma espécie possivelmente extinta na natureza")
#1.2 Population size, Distribution and Past Trends (if no plot data is available)
res1.2 <- c(ResearchLookup = "1.2", ResearchName = "Population size, Distribution and Past Trends", 
            note = "Espécie não registrada em nenhum dos mais de 1100 inventários florestais compilados para a Mata Atlântica, portanto sem estimativas de seu tamanho populacional e tendências de declínio. Adicionalmente, a espécie foi registrada em poucas localidades (menos que 10). Provavelmente trata-se de uma espécie que ocorre em baixas densidades na floresta (i.e. espécie rara) ou uma espécie de ocrrência restrita à um tipo de habitat específico e sub-representado nos registros de herbário e inventários florestais compilados")
#1.3	Life History and Ecology
res1.3 <- c(ResearchLookup = "1.3", ResearchName = "Life History and Ecology", 
            note = "Espécie classificada como criticamente ameaçada, porém com pouca ou nenhuma informação sobre sua distribuição, ecologia e história de vida")
#1.4	Harvest, Use and Livelihoods			
res1.4 <- c(ResearchLookup = "1.4", ResearchName = "Harvest, Use and Livelihoods", 
            note = "Espécie classificada como criticamente ameaçada, porém sem nenhuma informação sobre usos, colheita e seu papel para a subsistência de populações locais")
#1.5	Threats		
res1.5 <- c(ResearchLookup = "1.5", ResearchName = "Threats", 
            note = "Espécie classificada como criticamente ameaçada, porém sem nenhuma informação específica sobre as ameaças à manutenção a longo prazo de suas populações na natureza")
# #1.6	Conservation Actions	Research to determine how to mitigate particular threats (e.g., how to mitigate impacts of long-lining to reduce bycatch), or whether ex-situ breeding is possible	
# res1.6 <- c(ResearchLookup = "1.6", ResearchName = NA, note = NA)

# 2	Conservation Planning					
#2.1	Species Action/Recovery Plan	depend (if CR or CR PE)			
res2.1 <- c(ResearchLookup = "2.1", ResearchName = "Species Action/Recovery Plan", 
            note = "Estabelecer um Plano de Ação em escala nacional para recuperar as populações da espécie")
#2.2	Area-based Management Plan					
res2.2 <- c(ResearchLookup = "2.2", ResearchName = "Area-based Management Plan", note = "")
#2.3	Harvest and Trade Management Plan
res2.3 <- c(ResearchLookup = "2.3", ResearchName = "Harvest and Trade Management Plan", 
            note = "Espécie criticamente ameaçada com algum tipo de uso conhecido, real ou potencial, para a obtenção de produtos madereiros ou não-madereiro")

# 3	Monitoring					
#3.1	Population Trends (if threatened based on C or D)	
res3.1 <- c(ResearchLookup = "3.1", ResearchName = "Population Trends", 
            note = "Espécie com estimativas ou indícios de tamanho populacional pequeno e em declínio ou com tamanho populacionais muito pequenos. Estudos populacionais mais refinados são necessários para melhor estimar qual é o real tamanho populacional da espécie ao longo da sua área de distribuição")
#3.2	Harvest level trends depend (if overharvested)
res3.2 <- c(ResearchLookup = "3.2", ResearchName = "Harvest level Trends", 
            note = "Espécie com histórico de exploração na Mata Atlântica. Apesar das estimativas do tamanho populacional da espécie estarem disponíveis (baseadas em perda de área florestal), pouco se conhece sobre a intensidade da exploração histórica e de seu impacto no tamanho populacional efetivo da espécie ao longo do tempo. Estudos populacionais mais refinados são necessários para melhor estimar qual a real situação das populacões ao longo da área de distribuição da espécie")
#3.3	Trade Trends	depend (if traded)?	Espécie de uso comercial comum na Mata Atlântica. Contudo, pouco se conhece sobre o impacto desse uso sobre a espécie.
res3.3 <- c(ResearchLookup = "3.3", ResearchName = "Trade Trends", note = "")
#3.4	Habitat Trends
res3.4 <- c(ResearchLookup = "3.4", ResearchName = "Habitat Trends", note = "")


## CREATING COMBINED INFO FOR ALL SPECIES ##
# tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
# all.spp <- merge(all.crit, tax,
#                  by.x = "species", by.y = "species.correct2", all.x = TRUE, sort = FALSE)
# cites <- readRDS("data/threat_cites_EU.rds")
# all.spp <- merge(all.spp, cites[ ,c("internal_taxon_id", "cites_listing", "appendix", "annotation")],
#                  by.x = "internal_taxon_id", 
#                  by.y = "internal_taxon_id", all.x = TRUE, sort = FALSE)
# all.spp <- all.spp[order(as.double(gsub("sp", "", all.spp$internal_taxon_id))), ] 
n.spp <- dim(all.spp)[1]

#species uses
used.spp1 <- read.csv("data/threat_species_uses.csv")
used.spp2 <- read.csv("data/sis_connect/usetrade_threat.csv")
used.spp4 <- read.csv("data/threat_exploited_timber_spp.csv")


threat.cats <- c("EN","VU","CR")
## 1 - RESEARCH ##
# 1.1	Taxonomy (types only or only old records)
# Only old records
old.sp <- all.spp[((!is.na(all.spp$last.record) & 
                     all.spp$last.record < (2018-50)) |
                      (all.spp$only.type %in% TRUE &
                         !is.na(all.spp$last.record) &
                            all.spp$last.record < (2018-50)) |
                              grepl("^sp566$|^sp3874$", all.spp$internal_taxon_id)) &
                        # all.spp$cat.reg.clean %in% threat.cats &
                          all.spp$endemic %in% c("endemic"),
                  c("internal_taxon_id","species")]
                  # c("internal_taxon_id","species","cat.reg.clean","endemic","last.record","only.type")]
taxonomy <- do.call(rbind.data.frame,
                        replicate(dim(old.sp)[1], res1.1, simplify = FALSE))
names(taxonomy) <- names(res1.1)
taxonomy$internal_taxon_id <- old.sp$internal_taxon_id


#1.2 Population size, distribution & trends	(if no plot data is available)
herb.sp <- all.spp[(is.na(all.spp$reduction_A12) &
                     (is.na(all.spp$nbe_loc_total) |
                        all.spp$nbe_loc_total < 10) | 
                          (is.na(all.spp$nbe_loc_total) | 
                             all.spp$Nbe_occs < 3)) &
                              all.spp$endemic %in% c("endemic") & 
                     grepl("tree", hab1$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName, ignore.case = TRUE),
                  #c("internal_taxon_id","species")]
                  c("internal_taxon_id","species","endemic","cat.reg.clean","main.criteria",
                    "Nbe_occs","nbe_loc_total","pop.size","pop.size.low")]

population <- do.call(rbind.data.frame,
                    replicate(dim(herb.sp)[1], res1.2, simplify = FALSE))
names(population) <- names(res1.2)
population$internal_taxon_id <- herb.sp$internal_taxon_id

#1.3	Life history & ecology &					
# Rare populations classified as CR
herb.CR.sp <- herb.sp[herb.sp$cat.reg.clean %in% "CR",]
ecology <- do.call(rbind.data.frame, 
                   replicate(dim(herb.CR.sp)[1], res1.3, simplify = FALSE))
names(ecology) <- names(res1.3)
ecology$internal_taxon_id <- herb.CR.sp$internal_taxon_id

#1.4	Harvest, use & livelihoods & 					
# Rare populations classified as CR, but without known uses
livelihood.sp <- herb.sp[herb.sp$cat.reg.clean %in% "CR" & 
                           !herb.sp$internal_taxon_id %in% used.spp1$internal_taxon_id ,]
livelihood <- do.call(rbind.data.frame, 
                   replicate(dim(livelihood.sp)[1], res1.4, simplify = FALSE))
names(livelihood) <- names(res1.4)
livelihood$internal_taxon_id <- livelihood.sp$internal_taxon_id

#1.5	Threats
# Rare populations classified as CR
threat.sp <- herb.sp[herb.sp$cat.reg.clean %in% "CR" & 
                           !herb.sp$internal_taxon_id %in% used.spp1$internal_taxon_id ,]
threat <- do.call(rbind.data.frame, 
                      replicate(dim(threat.sp)[1], res1.5, simplify = FALSE))
names(threat) <- names(res1.5)
threat$internal_taxon_id <- threat.sp$internal_taxon_id


#1.6	Conservation Actions	Research to determine how to mitigate particular threats (e.g., how to mitigate impacts of long-lining to reduce bycatch), or whether ex-situ breeding is possible	
#Nothing for now

## 2 - CONSERVATION PLANNING ##
#	Species Action/Recovery Plan depend (if CR or CR PE)			
plan.sp <- all.spp[all.spp$cat.reg.clean %in% "CR" &
                       all.spp$endemic %in% c("endemic"), c("internal_taxon_id","species")]
plan <- do.call(rbind.data.frame, 
                  replicate(dim(plan.sp)[1], res2.1, simplify = FALSE))
names(plan) <- names(res2.1)
plan$internal_taxon_id <- plan.sp$internal_taxon_id


# #2.2	Area-based Management Plan					
#Nothing for now

#2.3	Harvest & Trade Management Plan					
harvest.sp <- herb.sp[herb.sp$cat.reg.clean %in% "CR" & 
                           herb.sp$internal_taxon_id %in% used.spp1$internal_taxon_id ,]
# harvest <- do.call(rbind.data.frame, 
#                       replicate(dim(harvest.sp)[1], res2.3, simplify = FALSE))
# names(harvest) <- names(res2.3)
# harvest$internal_taxon_id <- harvest.sp$internal_taxon_id
# NENHUMA ESPÉCIE COM ESSES CRITÉRIOS

## 3 - MONITORING ##
#3.1	Population trends	depend (if threatened or near-threatened based on C or D)	
monitor.sp <- all.spp[((all.spp$cat.reg.clean %in% c(threat.cats) & 
                        grepl("C|D", all.spp$main.criteria)) | 
                         (all.spp$cat.reg.clean %in% c(threat.cats, "NT") &
                            !is.na(all.spp$pop.size) & 
                              (all.spp$pop.size < 10000 | all.spp$pop.size.low < 2500))) & 
                                all.spp$endemic %in% c("endemic") & 
                                  grepl("tree", hab1$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName, ignore.case = TRUE),
                      c("internal_taxon_id","species","endemic",
                        "cat.reg.clean","main.criteria",
                        "reduction_A12","pop.size","pop.size.low"),]
monitor <- do.call(rbind.data.frame,
                      replicate(dim(monitor.sp)[1], res3.1, simplify = FALSE))
names(monitor) <- names(res3.1)
monitor$internal_taxon_id <- monitor.sp$internal_taxon_id


#3.2	Harvest level trends depend (if overharvested/exploited)
explo.sp <- unique( 
            all.spp[all.spp$internal_taxon_id %in% all.explo.sp$internal_taxon_id & 
              all.spp$cat.reg.clean %in% threat.cats &
                all.spp$endemic %in% c("endemic", "widespread_common"), 
          c("internal_taxon_id", "species")])
exploited <- do.call(rbind.data.frame,
                     replicate(dim(explo.sp)[1], res3.2, simplify = FALSE))
names(exploited) <- names(res3.2)
exploited$internal_taxon_id <- explo.sp$internal_taxon_id


#3.3	Trade trends	depend (if traded)?	Espécie de uso comercial comum na Mata Atlântica. Contudo, pouco se conhece sobre o impacto desse uso sobre a espécie.
res3.3
#Nothing for now

#3.4	Habitat trends					
res3.4
#Nothing for now

## COMBINING ALL THREATS 
research <- rbind.data.frame(taxonomy, population, ecology, 
                             livelihood, threat, plan,
                             monitor, exploited, 
                                 stringsAsFactors = FALSE)
#Re-ordering
research <- research[order(as.double(gsub("sp","", 
                                                  research$internal_taxon_id)), 
                                   research$ResearchLookup), ] 
#Removig possible duplicates
dups <- duplicated(paste(research$internal_taxon_id,  
                         research$ResearchLookup, sep = "_"))
research <- research[!dups, ]

#Renaming
head(sample, 3)
head(research, 3)
names(research)[1:3] <- paste0("Research.ResearchSubfield.", 
                                   names(research)[1:3])
table(names(research) == names(sample))

## ADDING SPECIFIC NOTES 
##Specific comments for EX and EW species
#Campomanesia lundiana: 
text <- c("Existe apenas um registro recente (2018) para a espécie, feito no Monte Cabrão, Santos, São Paulo, Brasil. 
Esse registro foi feito e identificado por I. Villaça and Mara Magenta, que não são especialistas da família. 
Assim, esse registro pode ser um erro de identificação, em especial visto que o estudo no qual a espécie foi listada inclui
 apenas espécies usadas para o consumo de frutos pela população local (provavél confusão com Campomaneisa phaea?)")

replace_these <- grepl("sp566", research$internal_taxon_id) & 
                    grepl("1.1", research$Research.ResearchSubfield.ResearchLookup)
research$Research.ResearchSubfield.note[replace_these] <-
  paste(research$Research.ResearchSubfield.note[replace_these],
        stringr::str_squish(gsub("\\\n", "", text)),
        sep = ". ")

# Pouteria stenophylla: 
text <- c("Espécie conhecida basicamente para a sua localidade tipo (provavelmente Serra dos Órgãos, Rio de Janeiro, Brasil). 
Apenas um registro recente (2010) e validado taxonomicamente foi encontrado para esta espécie (Palazzo, F.M.A. 40)
 para as restingas do município de Rio das Ostras, Rio de Janeiro. Há outro registro recente foi feito em Prado, Bahia (Rezende, S.G. 1833), 
 mas sua identificação ainda não foi validada por um especialista. Há também registros recentes feitos em afloramentos rochosos no estado de Minas Gerais 
 (Ouro Preto e Mariana), mas novamente a validação por um especialistas da família ainda está pendente. 
 ver tb: PALAZZO, F.M.A.; DIAS-NETO, A.O.; MONTEIRO, M.H.D.A.; ANDREATA, R.H.P. 2010. Sinopse comentada de Sapotaceae no município de Rio das Ostras (RJ, BRASIL). Pesquisas (Botânica) 61: 293-306.")

replace_these <- grepl("sp3874", research$internal_taxon_id) & 
                    grepl("1.1", research$Research.ResearchSubfield.ResearchLookup)
research$Research.ResearchSubfield.note[replace_these] <-
  paste(research$Research.ResearchSubfield.note[replace_these],
        stringr::str_squish(gsub("\\\n", "", text)),
        sep = ". ")

## Adding conservation info from previous IUCN assessments
iucn.sis.research <- iucn.sis[iucn.sis$sic.connect.file %in% "research_needed", ]
research.1 <- dplyr::left_join(research, tax)
research1 <- .merge_sis_connect(research, type = "research_needed", tax.table = tax,
                                    iucn.sis.research, x.spp = research.1$species.correct2, replace = TRUE)

#Adding descriptions to the new categories from the previous assessments
new.cats <- c("1.6", "2.2", "2.3", "3.3", "3.4")
names(new.cats) <- c("Pesquisa é necessária para determinar como mitigar ameaças específicas (e.g. extração ilegal, fogo, espécies invasoras), ou se a reprodução em cativeiro/viveiro é viável",
                     "Pesquisa necessária para avaliar a necessidade e viabilidade de manejo in situ das populações naturais da espécie visando a sua conservação",
                     "Espécie com algum tipo de uso conhecido, real ou potencial, para a obtenção de produtos madereiros ou não-madereiros. Pesquisa necessária para determinar a intensidade e extensão do uso da espécie",
                     "Espécie de uso comercial comum na Mata Atlântica. Contudo, pouco se conhece sobre o impacto desse uso sobre as populações da espécie",
                     "Pesquisa é necessária para detectar os habitats preferenciais da espécie e para avaliar a evolução da quantidade e qualidade desses habitats ao longo do tempo")

for (i in seq_along(new.cats)) {
    research1$Research.ResearchSubfield.note[research1$Research.ResearchSubfield.ResearchLookup %in% 
                                                                       new.cats[i]] <- names(new.cats)[i]
}
research1$Research.ResearchSubfield.note.iucn_rl[
  is.na(research1$Research.ResearchSubfield.note.iucn_rl)
] <- ""

## Saving the SIS Connet file
write.csv(research, "data/sis_connect/researchneeded_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(research1, "data/sis_connect/researchneeded_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

###############################################################################H
###############################################################################H
# ASSESSMENTS ------------------------------------------------------------------

## Generating the SIS CONNECT file "assessments.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/assessments.csv")
# head(sample[,8:17], 4)
# apply(sample, 2, unique)


## Biogeographic Realm
#realms: Afrotropical, Antarctic, Australasian, Indomalayan, Nearctic, Neotropical, Oceanian, Palearctic
countries0 <- readRDS("data/countries_per_species.rds")
# regions <- as.data.frame(readxl::read_xlsx("E:/ownCloud/W_GIS/WO_Biogeographical_divisions/country_list.xlsx"))
regions <- as.data.frame(readxl::read_xlsx("./data/outside_data//country_list.xlsx"))
countries0 <- merge(countries0, regions[, c("Country_Names", "WWF_Realm", "NAME_PT")],
                    by.x = "NAME_0", by.y = "Country_Names", sort = FALSE, all.x = TRUE)
# countries0$UN_Bioregion <- gsub(" Region", "", countries0$UN_Bioregion)
sort(unique(countries0$WWF_Realm))
tmp <- as.data.frame.matrix(table(countries0$species.correct2, countries0$WWF_Realm))
tmp <- tmp[, -which(names(tmp) == "NA")]
tmp1 <- apply(tmp, 1, function(x) 
                        paste0(names(tmp)[x > 0], collapse = "|"))
# tmp1 <- table(countries0$species.correct2)
# assess <- data.frame(species.correct2 = names(tmp1),
#                      BiogeographicRealm.realm = as.double(tmp1))
assess <- data.frame(species.correct2 = names(tmp1),
                     BiogeographicRealm.realm = tmp1)
assess <- assess[assess$species.correct2 %in% all.spp$species, , drop = FALSE]
assess <- assess[match(assess$species.correct2, all.spp$species), , drop = FALSE]
table(assess$species.correct2 == tax$species.correct2)
# toto <- read.csv("data/sis_connect/assessments_threat.csv", encoding = "UTF-8")
# table(toto$species.correct2 == assess$species.correct2)
# assess$BiogeographicRealm.realm <- toto$BiogeographicRealm.realm
assess$internal_taxon_id <- tax$internal_taxon_id

## Population Trend
table(assess$species.correct2 == all.spp$species)
# cont.decline: overall pop trend based on fitted statistical models (last 3 generations or last period if fit is piecewise - Criterio C1)
# any.decline: any rate of mean pop change (no stat fit), but only between 2000-2018 (Criterio C2)
# declineB: continuing decline based on habitat loss between 2000 and 2015, corrected for pioneer species (Criterio B)
# head(all.spp[,c("cont.decline", "any.decline", "declineB")])

tmp <- strsplit(all.spp$cont.decline, "\\|")
tmp <- sapply(tmp, function(x) gsub("\\(|\\)|[0-9]", "", x))
tmp <- sapply(tmp, function(x) gsub(" -", "", x))

# Priorizando o ajuste para o ultimo periodo (em geral, 1995-2018)
est.decline <- sapply(tmp, tail, 1) # usar o ajuste para o ultimo periodo apenas? 
table(est.decline, useNA = "always")
table(est.decline, all.spp[ ,c("any.decline")], useNA = "always")
table(est.decline, all.spp[ ,c("declineB")], useNA = "always")

# est.decline[!est.decline %in% "Decreasing" & all.spp$any.decline %in% "Decreasing"] <- "Stable"
est.decline[est.decline %in% "Increasing" & all.spp$declineB %in% "Decreasing"] <- "Stable" 
est.decline[est.decline %in% "Decreasing" & all.spp$any.decline %in% "Increasing"] <- "Stable" 

# Se declinio não está disponível, usar o declinio no forest cover 2000-2015
est.decline[is.na(est.decline)] <- gsub("Not Decreasing", "Stable",
  all.spp[ ,c("declineB")][is.na(est.decline)]) 

# Se declinio ainda não está disponível, usar o declinio para o ultimo assessment
est.decline[is.na(est.decline)] <- all.spp$populationTrend[is.na(est.decline)]

# Se declinio ainda não está disponível, fixar como "Unknown"
est.decline[is.na(est.decline)] <- "Unknown"

table(est.decline, useNA = "always")
assess$PopulationTrend.value <- est.decline

## System.value
assess$System.value <- "Terrestrial"

## version of the IUCN criteria
assess$RedListCriteria.critVersion <- 3.1

## Categories and Criteria 
assess$RedListCriteria.isManual <- "true"
assess$RedListCriteria.manualCategory <- all.spp$cat.reg.clean

crits <- all.spp$main.criteria
fill_these <- assess$RedListCriteria.manualCategory %in% c("EN","VU","CR")
crits[!fill_these] <- ""
crits <- gsub("A1+A2","A2b", crits, fixed = TRUE)
crits <- gsub("^A2$","A2b", crits, perl = TRUE)
crits <- gsub("^A2\\+","A2b+", crits, perl = TRUE)
crits <- gsub("B1+B2","B1ab(iii)_2ab(iii)", crits, fixed = TRUE)
crits <- gsub("B2","B2ab(iii)", crits, perl = TRUE)
crits <- gsub("B1$","B1ab(iii)", crits, perl = TRUE)
# crits <- gsub("B1","B1a,B1biii", crits, fixed = TRUE)
# crits <- gsub("B2","B2a,B2biii", crits, fixed = TRUE)
crits <- gsub("C2","C2a(ii)", crits, fixed = TRUE)
crits <- gsub("C1+C2","C1_2", crits, fixed = TRUE)
crits <- gsub("D","D1", crits, fixed = TRUE)
crits <- gsub("\\+","; ", crits, perl = TRUE)
crits <- gsub("_","+", crits, perl = TRUE)

# crits <- gsub("B1ab(iii)+B2ab(iii)","B1ab(iii)+2ab(iii)", crits, fixed = TRUE)
exploit <- grepl("Extra reduction", all.spp$basis_d)
crits[exploit] <- gsub("A2b","A2bd", crits[exploit], fixed = TRUE)
assess$RedListCriteria.manualCriteria <- crits
assess$RedListCriteria.manualCriteriaString <- ""

## Possibly extinct?
assess$RedListCriteria.possiblyExtinct <- "false"
assess$RedListCriteria.possiblyExtinctCandidate <- "false"
replace_these <- grepl("_PE", all.spp$category.regional) & 
                  assess$RedListCriteria.manualCategory %in% c("EN","VU","CR") & 
                    all.spp$endemic %in% "endemic" & all.spp$endemic.BR %in% "endemic"
assess$RedListCriteria.possiblyExtinctCandidate[replace_these] <- "true"
assess$RedListCriteria.yearLastSeen <- ""
# assess$RedListCriteria.yearLastSeen[replace_these] <- 
#   all.spp$last.record[replace_these]

## Data deficient 
assess$RedListCriteria.dataDeficientReason <- ""
replace_these <- assess$RedListCriteria.manualCategory %in% "DD" 
assess$RedListCriteria.dataDeficientReason[replace_these] <- "Taxonomic"

## Date of the Assessment (dd/mm/yyyy)
assess$RedListAssessmentDate.value <- "28/05/2021"

## CREATING COMBINED INFO FOR ALL SPECIES ##
# tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 3)]
# all.spp <- merge(all.crit, tax,
#                  by.x = "species", by.y = "species.correct2", all.x = TRUE, sort = FALSE)
# cites <- readRDS("data/threat_cites_EU.rds")
# all.spp <- merge(all.spp, cites[ ,c("internal_taxon_id", "cites_listing", "appendix", "annotation")],
#                  by.x = "internal_taxon_id", 
#                  by.y = "internal_taxon_id", all.x = TRUE, sort = FALSE)
# all.spp <- all.spp[order(as.double(gsub("sp", "", all.spp$internal_taxon_id))), ] 

## Re-assessments
prev.assess <- readRDS("data/sis_connect/prev_assessments_threat.rds")
all.spp0 <- all.spp[, c("internal_taxon_id", names(all.spp)[!names(all.spp) %in% names(prev.assess)])]
all.spp1 <- dplyr::left_join(all.spp0, prev.assess, by = "internal_taxon_id", suffix = c("", ".iucn_rl"))
all.spp1$assessmentYear <- plantR::getYear(all.spp1$assessmentDate, noYear = "")
reasses <- !is.na(all.spp1$redlistCategory) & 
              !all.spp1$cat.reg.clean %in% c("NA", NA)

# Reformatting categories
redlistCategory <- all.spp1$redlistCategory
redlistCategory[redlistCategory %in% "Critically Endangered"] <- "CR"
redlistCategory[redlistCategory %in% "Data Deficient"] <- "DD"
redlistCategory[redlistCategory %in% "Endangered"] <- "EN"
redlistCategory[redlistCategory %in% "Extinct"] <- "EX"
redlistCategory[redlistCategory %in% c("Least Concern", "Lower Risk/least concern")] <- "LC"
redlistCategory[redlistCategory %in% c("Lower Risk/near threatened", "Lower Risk/conservation dependent", "Near Threatened")] <- "NT"
redlistCategory[redlistCategory %in% "Vulnerable"] <- "VU"

# Type of change: "No change", "Genuine Change", "Nongenuine change"
change <- reasses & all.spp1$cat.reg.clean != redlistCategory
assess$RedListReasonsForChange.type <- ""
assess$RedListReasonsForChange.type[change] <- "Nongenuine change"
no.change <- reasses & all.spp1$cat.reg.clean == redlistCategory
assess$RedListReasonsForChange.type[no.change] <- "No change"

assess$RedListReasonsForChange.timeframe <- ""
assess$RedListReasonsForChange.timeframe[change] <- "Since first assessment"

assess$RedListReasonsForChange.changeReasons <- ""
assess$RedListReasonsForChange.changeReasons[change] <- "New information"
# assess$RedListReasonsForChange.changeReasons[change & !is.na(all.spp1$category_A)] <- "New information"
# assess$RedListReasonsForChange.changeReasons[change & is.na(all.spp1$category_A)] <- "Other"
assess$RedListReasonsForChange.otherReason <- ""

#Change in the criteria
crits <- all.spp1$main.criteria
fill_these <- assess$RedListCriteria.manualCategory %in% c("EN","VU","CR")
crits[!fill_these] <- ""
crits.iucn <- all.spp1$redlistCriteria
crits.iucn <- gsub("\\(|\\)", "", crits.iucn)
crits.iucn <- gsub("i|v|[a-z]", "", crits.iucn)
crits.iucn <- gsub(",,,", ",", crits.iucn)
crits.iucn <- gsub(",,", ",", crits.iucn)
crits.iucn <- gsub(",$", "", crits.iucn)
crits.iucn <- gsub(",\\+", "+", crits.iucn)
crits.iucn <- gsub(",; ", "+", crits.iucn)
crits.iucn <- gsub(", ", "+", crits.iucn)
crits.iucn <- gsub("; ", "+", crits.iucn)
crits.iucn <- gsub("B1\\+2", "B1+B2", crits.iucn)
crits.iucn <- gsub("A1\\+2", "A1+A2", crits.iucn)
crits.iucn <- gsub("A2\\+3", "A2+A3", crits.iucn)
crits.iucn <- gsub("A3\\+4", "A3+A4", crits.iucn)
crits.iucn <- gsub("A2\\+4", "A2+A4", crits.iucn)
crits.iucn <- gsub("D1", "D", crits.iucn)
change.crit <- reasses & crits != crits.iucn & !is.na(crits.iucn) &
                  assess$RedListReasonsForChange.type %in% "No change"
no.change.crit <- reasses & crits == crits.iucn & !is.na(crits.iucn) &
                    assess$RedListReasonsForChange.type %in% "No change"
assess$RedListReasonsForChange.catCritChanges <- ""
assess$RedListReasonsForChange.catCritChanges[change.crit] <- 
  "Same category but change in criteria"
assess$RedListReasonsForChange.catCritChanges[no.change.crit] <- 
  "Same category and criteria"

# Reasons for change
# assess$ReasonForChangeJustification.narrative <- ""

## Language
assess$language.value <- "Portuguese"

## Range description (all assessments but LC and NA)
# describe <- all.spp1$cat.reg.clean %in% c("EN", "VU", "CR", "NT", "DD")
describe <- all.spp1$cat.reg.clean %in% c("EN", "VU", "CR", "NT", "DD", "LC") # precisa pra todas se não o SIS Connect dá erro
endemismo <- all.spp1$endemic
replace_these <- endemismo %in% "endemic" & 
                    all.spp1$endemism.level.1 < 99.5 & all.spp1$endemism.level.2 < 99.5 &
                      !is.na(all.spp1$endemism.level.1) & !is.na(all.spp1$endemism.level.2)
endemismo[replace_these] <- "near_endemic"

assess$RangeDocumentation.narrative <- "Dados de distribuição não reportados para a espécie, pois ela é considerada ocasional na região estudada (i.e. categoria 'NA')"
assess$RangeDocumentation.narrative[describe] <-
  c("Espécie _ENDEMISMO_. Possui ocorrências confirmadas para _DISTRIBUICAO_")

# Endemicas do BR
br <- all.spp1$endemic.BR %in% "endemic"
# Endemicas puras do BR e da MA
texto <- "endêmica do Brasil (Flora do Brasil 2020) e da Mata Atlântica (Lima et al. 2020a)"
replace_these <- endemismo %in% "endemic" & br
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])
  
# Endemicas do BR e sub-endemicas da MA
texto <- "endêmica do Brasil (Flora do Brasil 2020) e sub-endêmica da Mata Atlântica (Lima et al. 2020a)"
replace_these <- endemismo %in% "near_endemic" & br 
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])

# Endemicas do BR e com maiorias das ocorrências na MA
texto <- "endêmica do Brasil (Flora do Brasil 2020) e com a maioria de seus registros concentrados na Mata Atlântica (Lima et al. 2020a)"
replace_these <- endemismo %in% "widespread_common" & br 
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])

# Demais endêmicas da MA
texto <- "endêmica da Mata Atlântica (Lima et al. 2020a), porém não é endêmica do Brasil (Flora do Brasil 2020)"
replace_these <- endemismo %in% "endemic" & !br 
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])

# Demais endêmicas da MA
texto <- "sub-endêmica da Mata Atlântica (Lima et al. 2020a), porém não é endêmica do Brasil (Flora do Brasil 2020)"
replace_these <- endemismo %in% "near_endemic" & !br 
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])

# Demais espécies
texto <- "com uma certa parte de seus registros ocorrendo dentro da Mata Atlântica (Lima et al. 2020a)"
replace_these <- endemismo %in% "widespread_common" & !br 
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])
texto <- "com uma pequena parte de seus registros ocorrendo dentro da Mata Atlântica (Lima et al. 2020a)"
replace_these <- endemismo %in% c("widespread_sparse", "occasional")
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("_ENDEMISMO_", texto, assess$RangeDocumentation.narrative[replace_these])

# DISTRIBUIÇÃO
texto <- "Espécie com distribuição Neotropical, "
replace_these <- assess$BiogeographicRealm.realm %in% "Neotropical"
assess$RangeDocumentation.narrative[replace_these] <- 
  gsub("^Espécie ", texto, assess$RangeDocumentation.narrative[replace_these])

#paises e estados
tmp <- as.data.frame.matrix(table(countries0$species.correct2, countries0$NAME_0))
tmp1 <- apply(tmp, 1, function(x) 
  paste0(names(tmp)[x > 0], collapse = "|"))
paises <- data.frame(species.correct2 = names(tmp1),
                     paises = tmp1)
tmp <- as.data.frame.matrix(table(countries0$species.correct2, countries0$NAME_1))
tmp1 <- apply(tmp, 1, function(x) 
  paste0(names(tmp)[x > 0], collapse = "|"))
paises$estados <- as.character(tmp1)
paises <- paises[paises$species.correct2 %in% all.spp1$species,]
paises <- paises[match(paises$species.correct2, all.spp1$species),]

# Replacing country names to Portuguese
paises.pt <- unique(unlist(strsplit(paises$paises, "\\|")))
paises.pt <- paises.pt[order(nchar(paises.pt), decreasing = TRUE)]
paises.pt.org <- paises.pt
paises.pt <- as.character(countrycode::countryname(paises.pt, destination = "cldr.name.pt"))
paises.pt[paises.pt.org %in% "Micronesia"] <- "Micronésia"
names(paises.pt) <- paises.pt.org 
paises$paises.pt <- paises$paises
replace_these <- assess$language.value %in% "Portuguese"
paises$paises.pt[replace_these] <-
  stringr::str_replace_all(paises$paises.pt[replace_these], paises.pt)
paises$paises.pt <- gsub("Franças", "Frances", paises$paises.pt)

# Replacing the endemism and distribution information
for (i in 1:dim(paises)[1]) {
  if (!br[i]) {
    repos <- paises$paises.pt[i]
    if (grepl("\\|", repos)) {
      repos <- paste0("os seguintes países: ", gsub("\\|", ", ", repos, perl = TRUE), ".")
    } else {
      repos <- paste0("o seguinte país: ", repos,".")
    }
    assess$RangeDocumentation.narrative[i] <-  
      gsub("_DISTRIBUICAO_",  repos, 
           assess$RangeDocumentation.narrative[i], fixed = TRUE)
    
    replace_these <- grepl("^Espécie com distribuição Neotropical, com",
                           assess$RangeDocumentation.narrative[i], perl = TRUE)
    if (replace_these)
      assess$RangeDocumentation.narrative[i] <-  
        gsub("^Espécie com distribuição Neotropical, com",
             "Espécie com distribuição Neotropical, não endêmica do Brasil (Flora do Brasil 2020), com", 
             assess$RangeDocumentation.narrative[i], perl = TRUE)
    
  } else {
    repos <- paises$estados[i]
    if (grepl("\\|", repos)) {
      repos <- paste0("os estados de: ", gsub("\\|", ", ", repos, perl = TRUE), ".")
    } else {
      repos <- paste0("o estado de ", repos,".")
    }
    assess$RangeDocumentation.narrative[i] <-  
      gsub("_DISTRIBUICAO_",  repos, 
           assess$RangeDocumentation.narrative[i], fixed = TRUE)
  }
}

## Map status (single select)??
assess$MapStatus.status <- "false"
check_these <- assess$species.correct2 %in% readLines("data/sis_connect/spp_with_shapes.txt") 
assess$MapStatus.status[check_these] <- "Done"

## Population description (all assessments but LC and NA)
info.threat <- readRDS("data/assess_iucn_spp.rds")
table(info.threat$species.correct2 == all.spp1$species)

assess$PopulationDocumentation.narrative <- "Dados populacionais inexistentes ou não reportados para a espécie, pois ela é considerada ocasional na região estudada (i.e. categoria 'NA')"
replace_these <- 
assess$PopulationDocumentation.narrative[describe] <-
  c("_HERBARIO_ para a espécie. ") 

# Espécie com pop estimates acima do limiar
info.threat$treeco.occs[is.na(info.threat$treeco.occs) & !is.na(all.spp1$reduction_A12)] <- 1 
add_these <- !is.na(info.threat$treeco.occs) & !grepl("D", all.spp1$main.criteria)
texto <- c("Adicionalmente, _TREECO_ entre os inventários florestais da Mata Atlântica. 
           As estimativas do tamanho populacional da espécie ficaram acima do limiar de 10,000 indíviduos adultos definidos pela IUCN (IUCN 2019), o que significa que as populações da espécie são relativamente grandes na Mata Atlântica. ")
assess$PopulationDocumentation.narrative[describe & add_these] <- 
  paste0(assess$PopulationDocumentation.narrative[describe & add_these], texto)

# Espécie com pop estimates abaixo do limiar
add_these <- !is.na(info.threat$treeco.occs) & grepl("D", all.spp1$main.criteria)
texto <- c("Adicionalmente, _TREECO_ entre os inventários florestais da Mata Atlântica. 
           As estimativas do tamanho populacional da espécie ficaram abaixo do limiar de 10,000 indíviduos adultos definidos pela IUCN (IUCN 2019), o que significa que as populações da espécie são pequenas na Mata Atlântica. ")
assess$PopulationDocumentation.narrative[describe & add_these] <- 
  paste0(assess$PopulationDocumentation.narrative[describe & add_these], texto)

add_these <- !is.na(info.threat$treeco.occs) &  !grepl("D", all.spp1$main.criteria) & 
                grepl("A", all.spp1$main.criteria)
texto <- c("Porém, foi estimado que essas populações sofreram uma redução de _REDUCAO_ nas últimas três gerações.")
assess$PopulationDocumentation.narrative[describe & add_these] <- 
  paste0(assess$PopulationDocumentation.narrative[describe & add_these], texto)

add_these <- is.na(info.threat$treeco.occs)
texto <- c("Não foi encontrado nenhum registro para a espécie entre as centenas de inventários florestais compilados para a Mata Atlântica, 
           sugerindo que se trata de uma espécie que ocorre em baixa densidade na Mata Atlântica (tamanho populacional pequeno) e/ou associada a algum habitat mal-representado entre os inventários compilados. 
           Foi estimado um tamanho populacional para a espécie baseado na relação entre as métricas de distribuição geográfica (i.e. EOO and AOO) e o tamanho populacional para as espécies da Mata Atlântica (_POPMEAN_LOW_). 
           Contudo, houve grande incerteza ao redor dessa estimativa (Lima et al. in prep.). Portanto, essa estimativa não foi usada para a avaliação de ameaça da espécie.")
assess$PopulationDocumentation.narrative[describe & add_these] <- 
  paste0(assess$PopulationDocumentation.narrative[describe & add_these], texto)

## Substituindo os valores
for (i in 1:dim(assess)[1]) {
  
  texto <- assess$PopulationDocumentation.narrative[i]

  if (!texto %in% "") {

    # Registros de herbário     
    # repos <- all.spp1$Nbe_occs[i]
    repos <- as.double(info.threat$non.dup.occs[i])
    if (repos == 1)
      repos <- "Foi encontrado apenas 1 registro de herbário"
    if (is.double(repos) & repos > 1 & repos < 6)
      repos <- paste0("Foram encontrados apenas ", repos, " registros de herbário espacialmente únicos")
    if (is.double(repos) & repos >= 6)
      repos <- paste0("Foram encontrados ", repos, " registros de herbário espacialmente únicos")
    texto <- gsub("_HERBARIO_", repos, texto, fixed = TRUE)
    
    # Registros nos inventários florestais     
    repos <- as.double(info.threat$treeco.occs[i])
    if (!is.na(repos)) {
      if (repos == 1)
        repos <- "foi encontrado apenas 1 registro da espécie"
      if (is.double(repos) & repos > 1)
        repos <- paste0("foram encontrados ", repos, " registros da espécie")
      texto <- gsub("_TREECO_", repos, texto, fixed = TRUE)
    }
    
    # Limiar geral de redução populacional     
    repos <- round(all.spp1$reduction_A12[i], 1)
    if (!is.na(repos)) {
      if (is.double(repos) & repos >= 90)
        repos <- "90% ou mais"
      if (is.double(repos) & repos >= 80)
        repos <- "80% ou mais"
      if (is.double(repos) & repos >= 70)
        repos <- "70% ou mais"
      if (is.double(repos) & repos >= 50)
        repos <- "50% ou mais"
      if (is.double(repos) & repos >= 30)
        repos <- "30% ou mais"
      texto <- gsub("_REDUCAO_", repos, texto, fixed = TRUE)
    }  
    
    # Pop size estimates based on AOO and EOO Limiar geral de redução populacional
    if (!is.na(all.spp1$pop.size[i])) {
      repos <- paste0("Predição média: ", round(all.spp1$pop.size[i], 0), 
                      "; Limite inferior do intervalo de confiança ao redor da estimativa: ",
                      round(all.spp1$pop.size.low[i], 0)," indivíduos adultos")
      texto <- gsub("_POPMEAN_LOW_", repos, texto, fixed = TRUE)
    } else {
      pat <- "Foi estimado um tamanho populacional para a espécie baseado na relação entre as métricas de distribuição geográfica (i.e. EOO and AOO) e o tamanho populacional para as espécies da Mata Atlântica (_POPMEAN_LOW_). "
      texto <- gsub(pat, "", texto, fixed = TRUE)
      pat <- "Contudo, houve grande incerteza ao redor dessas estimativas (Lima et al. in prep.). Portanto, elas não foram usadas para as avaliações de ameaça da espécie."
      texto <- gsub(pat, "", texto, fixed = TRUE)
      
    }
    
    # Fixing possible problems and replacing
    # texto <- stringr::str_squish(gsub("\\\n", "", texto))
    texto <- gsub("\\s+", " ", (gsub("\\\n", "", texto)), perl = TRUE)
    texto <- gsub("^ | $", "", texto, perl = TRUE)
    assess$PopulationDocumentation.narrative[i] <- texto
  }
}

## Fine-tuning the text
replace_these <- as.double(info.threat$non.dup.occs) > 15 &
                  grepl("apenas 1 registro da espécie entre", assess$PopulationDocumentation.narrative) &
                    grepl("relativamente grandes", assess$PopulationDocumentation.narrative) &
                      !grepl("sofreram uma redução", assess$PopulationDocumentation.narrative)
assess$PopulationDocumentation.narrative[replace_these] <- 
  gsub("\\. Adicionalmente,", ". Contudo,", 
       assess$PopulationDocumentation.narrative[replace_these], perl = TRUE)
assess$PopulationDocumentation.narrative[replace_these] <- 
  gsub("da espécie são relativamente grandes", "da espécie provavelmente não são pequenas", 
       assess$PopulationDocumentation.narrative[replace_these], perl = TRUE)
assess$PopulationDocumentation.narrative[replace_these] <- 
  gsub("pequenas na Mata Atlântica", "pequenas na Mata Atlântica (i.e. ocorrência baixa, mas com abundância local relativamente alta)", 
       assess$PopulationDocumentation.narrative[replace_these], perl = TRUE)


replace_these <- as.double(info.threat$non.dup.occs) <= 15 &
  grepl("apenas 1 registro da espécie entre", assess$PopulationDocumentation.narrative) &
  grepl("relativamente grandes", assess$PopulationDocumentation.narrative) &
  !grepl("sofreram uma redução", assess$PopulationDocumentation.narrative)
assess$PopulationDocumentation.narrative[replace_these] <- 
  gsub(" As estimativas do tamanho", " Contudo, as estimativas do tamanho", 
       assess$PopulationDocumentation.narrative[replace_these], perl = TRUE)
assess$PopulationDocumentation.narrative[replace_these] <- 
  gsub("da espécie são relativamente grandes", "da espécie provavelmente não são pequenas", 
       assess$PopulationDocumentation.narrative[replace_these], perl = TRUE)
assess$PopulationDocumentation.narrative[replace_these] <- 
  gsub("pequenas na Mata Atlântica", "pequenas na Mata Atlântica (i.e. ocorrência baixa, mas com abundância local relativamente alta)", 
       assess$PopulationDocumentation.narrative[replace_these], perl = TRUE)


## Habitat information (all assessments but LC and NA)
hab2 <- read.csv("data/threat_habitats.csv", fileEncoding = "UTF-8")
# head(hab2, 2)
names(hab2)[2] <- "species.correct2"
hab2.1 <- dplyr::left_join(hab2, spp.info)
# spp.info0 <- spp.info[spp.info$species.correct2 %in% all.spp1$species, ]
# spp.info0 <- spp.info0[match(spp.info0$species.correct2, all.spp1$species),]

table(hab2.1$species.correct2 == all.spp1$species)
# table(spp.info0$species.correct2 == all.spp1$species)

# assess$HabitatDocumentation.narrative <- ""
# assess$HabitatDocumentation.narrative[describe] <-
#   c("_HABITO_ atingindo até _ALTURA_. No Brasil, ocorre na _DOMINIO_, em particular na _VEGET_ (Flora do Brasil 2020).") 
assess$HabitatDocumentation.narrative <-
  c("_HABITO_ atingindo até _ALTURA_. No Brasil, ocorre na _DOMINIO_, em particular na _VEGET_ (Flora do Brasil 2020).") 

## Substituindo os valores
for (i in 1:dim(assess)[1]) {
  
  texto <- assess$HabitatDocumentation.narrative[i]
  
  if (!texto %in% "") {
    
    # Habito
    repos <- hab2.1$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName[i]
    if (grepl("Tree - small", repos))
      repos <- "Árvore pequena,"
    if (grepl("Tree - large", repos))
      repos <- "Árvore grande,"
    if (grepl("Shrub - large", repos))
      repos <- "Arbusto alto ou arvoreta,"
    if (grepl("Shrub - large", repos))
      repos <- "Arbusto alto ou arvoreta,"
    if (grepl("Tree - size unknown", repos))
      repos <- "Árvore,"
    if (grepl("Fern", repos))
      repos <- "Samambaia arborescente,"
    if (grepl("Succulent - tree", repos))
      repos <- "Árvore ou arbusto suculento,"
    texto <- gsub("_HABITO_", repos, texto, fixed = TRUE)
    
    # Registros nos inventários florestais     
    repos <- round(hab2.1$MaxSize.size[i]/100, 0)
    if (is.double(repos)) {
      repos <- paste0(repos, " metros de altura")
      texto <- gsub("_ALTURA_", repos, texto, fixed = TRUE)
    }
    
    # Domínios     
    repos <- gsub("\\|", ", ", hab2.1$dominioFitogeografico[i])
    texto <- gsub("_DOMINIO_", repos, texto, fixed = TRUE)
    
    # Vegetação
    repos <- gsub("\\|", ", ", hab2.1$tipoVegetacao[i])
    texto <- gsub("_VEGET_", repos, texto, fixed = TRUE)
    
    # Fixing possible problems and replacing
    # texto <- stringr::str_squish(gsub("\\\n", "", texto))
    texto <- gsub("\\s+", " ", gsub("\\\n", "", texto), perl = TRUE)
    texto <- gsub("^ | $", "", texto, perl = TRUE)
    
    assess$HabitatDocumentation.narrative[i] <- texto
  }
}
assess$HabitatDocumentation.narrative <-
  gsub("na Cerrado", "no Cerrado", assess$HabitatDocumentation.narrative)

## Final fixes after CNCFlora quality check
check_these <- grepl("_", assess$HabitatDocumentation.narrative)
if (any(check_these))
  hab1$HabitatDocumentation.narrative[check_these] <- 
    "Não foram encontradas informações sobre o habitat da espécie."

## Use and trade
#Reading uses data
usos1 <- read.csv("data/threat_species_uses.csv", as.is = TRUE)

# Preparing to replace info
replace_these <- grepl("wikipedia", usos1$source)
usos1$source[replace_these] <- "Wikipedia"

rep.usos <- c("alimentação (bebidas)", "produção de óleos essenciais", 
                     "contrução de cercas", "obtenção de fibras", "forragem para animais",
                     "alimentação (comidas)", "lenha ou carvão (combustíveis)",
                     "fonte de lipídeos", "uso medicinal", "ornamentação",
                     "produção de venenos", "produção de resinas ou gomas",
                     "produção de borracha ou látex", "produção de madeira")
names(rep.usos) <- sort(unique(usos1$uses))
usos1$uses1 <- stringr::str_replace_all(usos1$uses, rep.usos)

rep.refs <- c("Agroforestree database", "Brazilian Woods database", 
              "Carvalho 1998", "Coradin et al. 2011", "Coradin et al. 2018",
              "Lorenzi 1992", "Mark et al. 2014",
              "Personal comm.", "Sistema Nacional de Informações Florestais", 
              "Vieira et al. 2016", "Wikipedia")
names(rep.refs) <- sort(unique(usos1$source))

#Summary per species
tmp <- as.data.frame.matrix(table(usos1$Name_submitted, usos1$uses))
tmp <- apply(tmp, 1, function(x) 
       paste0(sort(names(tmp)[x > 0]), collapse = "|"))
tmp1 <- as.data.frame.matrix(table(usos1$Name_submitted, usos1$source))
tmp1 <- apply(tmp1, 1, function(x) 
  paste0(names(tmp1)[x > 0], collapse = "|"))
uses.trade <- data.frame(species = names(tmp),
                          usos = as.character(tmp),
                         refs = as.character(tmp1))

uses.trade$usos <- stringr::str_replace_all(uses.trade$usos, rep.usos)
uses.trade$refs <- stringr::str_replace_all(uses.trade$refs, stringr::coll(rep.refs))

## Substituindo os valores
assess$UseTradeDocumentation.value <- ""

for (i in 1:dim(assess)[1]) {
  spp <- assess$species.correct2[i]
  
  if (spp %in% uses.trade$species) {
    info <- uses.trade[uses.trade$species %in% spp, ]
    
    if (grepl("\\|", info$usos)) {
      texto <- paste0("A espécie possui os seguintes usos: ", 
                      gsub("\\|", ", ", info$usos), 
                      " (",
                      gsub("\\|", ", ", info$refs), ").")
    } else {
      texto <- paste0("A espécie é usada para ", info$usos, " (",
                      gsub("\\|", ", ", info$refs), ").")
    }
    
  } else {
    texto <- "Não foram encontrados registros sobre usos efetivos ou potenciais."
  }
  # assess$UseTradeDocumentation.value[i] <- stringr::str_squish(texto)
  texto <- gsub("\\s+", " ", texto, perl = TRUE)
  assess$UseTradeDocumentation.value[i] <- gsub("^ | $", "", texto, perl = TRUE)
  
}  

assess$UseTradeDocumentation.value <- 
  gsub(", Personal comm\\.\\)", ")", assess$UseTradeDocumentation.value, perl = TRUE)


## Threats information (all assessments but LC and NA)
assess$ThreatsDocumentation.value <- ""
texto <- c("O principal e mais comum tipo de ameaça para as espécies arbóreas da Mata Atlântica 
é a conversão de áreas de vegetação nativa para a agricultura e pastagem de larga-escala. 
Apenas as áreas mais montanhosas próximas ao litoral brasileiro e a região de Misiones na Argentina foram poupadas dessa conversão em larga-escala. 
Mais localmente, a Mata Atlântica também foi convertida para áreas urbanas, estradas e mineração. 
Esse processo foi mais intenso até os anos 2000, resultando na redução de mais de 80% da área original da Mata Atlântica. 
Mas ele continua em andamento, com taxas de desmatamento acima dos 20 mil hectares ao ano, e com taxas mais elevadas no Paraguay e 
em alguns estados brasileiros (Argentina & WWF 2017; SOS Mata Atlântica e INPE 2018; Rosa et al. 2021). 
Além do desmatamento (i.e. perda de habitat), a degradação florestal (e.g. corte seletivo, fogo, espécies invasoras) 
é um componente importante na Mata Atlântica, havendo indicações de que 80% ou mais dos remanescentes florestais 
da Mata Atlântica apresentam perdas variando entre 20 e 40% na sua biodiversidade (Lima et al. 2020b). 
Por fim, apenas 1% da área original da Mata Atlântica corresponde a Unidades de Conservação (Ribeiro et al. 2009). 
Como não existem informações sistematizadas e espacializadas sobre ameaças ao longo de toda a Mata Atlântica, 
não foi possível extrair informações sobre o tipo e magnitude das ameaças para as espécies individualmente.")
assess$ThreatsDocumentation.value[describe] <- texto

# Adding info from previous assessments
all.spp1$threats <- plantR:::fixEncoding(all.spp1$threats)
all.spp1$threats <- gsub("<em>|<\\/em>", "", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub("<p>|<\\/p>", "", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub('\\\"', "", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub('&#160;', "", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub('<span lang=en-US>', "", all.spp1$threats)
all.spp1$threats <- gsub('<br\\/><br\\/>', "", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub('Ã\u0081', "á", all.spp1$threats, fixed = TRUE)
all.spp1$threats <- gsub(' Ã ', " à ", all.spp1$threats, fixed = TRUE)
all.spp1$threats <- gsub("<span style=>|<\\/span>|<span>", "", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub("<strong><\\/strong>", "", all.spp1$threats, perl = TRUE)
# all.spp1$threats <- stringr::str_squish(all.spp1$threats)
all.spp1$threats <- gsub("\\s+", " ", all.spp1$threats, perl = TRUE)
all.spp1$threats <- gsub("^ | $", "", all.spp1$threats, perl = TRUE)

# add_these <- (!(grepl("^Dados publicados recentemente", all.spp1$threats, perl = TRUE) & 
#                        grepl("Atlântica e INPE 2018\\)\\.$", all.spp1$threats, perl = TRUE))) |
#               !grepl("taxon is not considered to be ", all.spp1$threats) |
#               !grepl("species is not known to be threatened", all.spp1$threats) |
#               !grepl("species is not threatened", all.spp1$threats) |
#               !grepl("taxon is not currently considered to be", all.spp1$threats)
# assess$ThreatsDocumentation.value.PrevAssessment <- ""
# assess$ThreatsDocumentation.value.PrevAssessment[describe & add_these] <- 
#   all.spp1$threats[describe & add_these]
# 
# patt <- c("Dados publicados recentemente (Fundação SOS Mata Atlântica e INPE 2018) apontam para uma redução maior que 85% da área originalmente 
# coberta com Mata Atlântica e ecossistemas associados no Brasil. De acordo com o relatório, cerca de 12,4% de vegetação original ainda resistem. 
# Embora a taxa de desmatamento tenha diminuído nos últimos anos, ainda está em andamento, e a qualidade e extensão de áreas florestais encontram-se 
# em declínio contínuo há pelo menos 30 anos (Fundação SOS Mata Atlântica e INPE 2018).")
# assess$ThreatsDocumentation.value.PrevAssessment <- 
#   gsub(patt, "", assess$ThreatsDocumentation.value.PrevAssessment, fixed = TRUE)
# assess$ThreatsDocumentation.value.PrevAssessment <-
#   gsub("\\s+", " ", assess$ThreatsDocumentation.value.PrevAssessment, perl = TRUE)
# assess$ThreatsDocumentation.value.PrevAssessment <-
#   gsub("^ | $", "", assess$ThreatsDocumentation.value.PrevAssessment, perl = TRUE)
# 
# replace_these <- !assess$ThreatsDocumentation.value.PrevAssessment %in% c("", NA)
# assess$ThreatsDocumentation.value[describe & replace_these] <- 
#   paste(assess$ThreatsDocumentation.value[describe & replace_these],
#         " Contudo, podem haver informações mais específicas sobre as ameaças às quais a espécie 
# está submetida na descrição feita para as avaliações anteriores. Abaixo, segue a transcrição literal 
# desta descrição (sem referências bibliográficas):_X_", 
#         assess$ThreatsDocumentation.value.PrevAssessment[describe & replace_these])
# 
# assess$ThreatsDocumentation.value <- gsub("\\\n", "", assess$ThreatsDocumentation.value, perl = TRUE)
# assess$ThreatsDocumentation.value <- gsub("_X_", "\n", assess$ThreatsDocumentation.value, perl = TRUE)
# # assess$ThreatsDocumentation.value <- stringr::str_squish(assess$ThreatsDocumentation.value)
# assess$ThreatsDocumentation.value <- gsub("\\s+", " ", assess$ThreatsDocumentation.value, perl = TRUE)
# assess$ThreatsDocumentation.value <- gsub("^ | $", "", assess$ThreatsDocumentation.value, perl = TRUE)


## Conservation Actions information
assess$ConservationActionsDocumentation.narrative <- ""

# With recors in UCs
protect <- !is.na(all.spp1$protected) & all.spp1$protected > 0 
# Endemicas puras do BR e da MA
texto <- "Há _XX_ para a espécie em Unidades de Conservação de Proteção Integral (categorias I, II ou III do 'IUCN Protected Area Management Categories' - www.iucn.org/theme/protected-areas), 
o que representa _YY_% das ocorrência confirmadas para a espécie. 
Adicionalmente, foi estimado que _ZZ_% da área do polígono que representa a Extensão de Ocorrência da espécie (EEO) corresponde a Unidades de Conservação de Proteção Integral. 
Não foram compiladas informações sobre a conservação ex situ da espécie."
assess$ConservationActionsDocumentation.narrative[protect] <- texto 

unprotect <- !is.na(all.spp1$protected) & all.spp1$protected == 0 
# Endemicas puras do BR e da MA
texto <- "Não há registros confirmados para a espécie em unidades de conservação. 
Foi estimado que _ZZ_% da área do polígono que representa a Extensão de Ocorrência da espécie (EEO) corresponde a Unidades de Conservação de Proteção Integral. 
Não foram compiladas informações sobre a conservação ex situ da espécie."
assess$ConservationActionsDocumentation.narrative[unprotect] <- texto 

texto <- "A espécie não possui registros suficientes para avaliar sua ocorrência em áreas protegidas. Não foram compiladas informações sobre a conservação ex situ da espécie."
assess$ConservationActionsDocumentation.narrative[is.na(all.spp1$protected)] <- texto


#Replacing text
for (i in 1:dim(assess)[1]) {
  
  spp <- assess$species.correct2[i]
  texto <- assess$ConservationActionsDocumentation.narrative[i]
  
  if (grepl("_XX_", texto)) {
    info <- all.spp1[all.spp1$species %in% spp, ]
    info1 <-info.threat[info.threat$species.correct2 %in% spp, ]
    occs1 <- round((info$protected/100) * info$Nbe_occs, 0)
    prop <- round(info$protected, 0)
    prop.eoo <- round(info$prop.EOO.in.StrictUCs, 0)
    
    if (occs1 == 1) {
      texto <- gsub("_XX_", "apenas um registro confirmado", texto, fixed = TRUE)
    } else {
      texto <- gsub("_XX_", paste0(occs1, " registros confirmados"), 
                    texto, fixed = TRUE)
                      
    }

    texto <- gsub("_YY_", prop, texto)
    if (!is.na(prop.eoo))
      texto <- gsub("_ZZ_", prop.eoo, texto, fixed = TRUE)

  } else {
    
    info <- all.spp1[all.spp1$species %in% spp, ]
    prop.eoo <- round(info$prop.EOO.in.StrictUCs, 0)
    
    if (!is.na(prop.eoo))
      texto <- gsub("_ZZ_", prop.eoo, texto, fixed = TRUE)
    
  }
  # assess$ConservationActionsDocumentation.narrative[i] <- 
  #   stringr::str_squish(gsub("\\\n", "", texto))
  texto <- gsub("\\s+", " ", gsub("\\\n", "", texto, perl = TRUE), perl = TRUE)
  assess$ConservationActionsDocumentation.narrative[i] <- 
    gsub("^ | $", "", texto, perl = TRUE)
}  


## Rationale for the Red List Assessment
assess$RedListRationale.value <- ""

for (i in 1:dim(assess)[1]) {
  
  spp <- assess$species.correct2[i]
  info <- all.spp1[all.spp1$species %in% spp, ]
  avalia <- assess[assess$species.correct2 %in% spp,]
  ameaca <- threats[threats$internal_taxon_id %in% info$internal_taxon_id,]
  conserva <- conservation[conservation$internal_taxon_id %in% info$internal_taxon_id,]
  pesquisa <- research[research$internal_taxon_id %in% info$internal_taxon_id,]
  
  texto <- ""
  
  if (info$cat.reg.clean %in% "NA") {
    texto <- "Esta espécie foi classificada como 'ocasional' para Mata Atlântica (Lima et al. 2020a), 
baseado na distribuição de seus registros de herbário. Por esse motivo, a espécie 
foi considerada como 'vagrant taxa' sensu IUCN (2012) durante as avaliações de ameaça. 
Seguindo as recomendações da IUCN para avaliações regionais (IUCN 2012), essa espécie 
não foi avaliada quanto ao seu estado de conservação e teve a sua categoria alterada 
para ‘Not Applicable’ (NA)."
  }

  if (info$cat.reg.clean %in% "DD") {
    texto <- "A espécie possui ocorrência confirmada para a Mata Atlântica, 
porém não foi encontrado nenhum registro de herbário e nenhuma dado populacional 
para a espécie. Ou seja, atualmente não há dados suficientes que permitam a avaliação 
dos critérios A, B, C ou D. Essa decisão foi baseada nas recomendações da IUCN 
que sugere “o uso de quaisquer informações disponíveis” e “de classificar táxons como 
Dados Deficientes apenas quando realmente não nehuma outra alternativa” 
(IUCN Standards and Petitions Committee 2019)."
  }

  if (info$cat.reg.clean %in% "NT") {
    texto <- "A espécie foi inicialmente classificada como Menor Preocupação (LC) após a valiação 
dos critérios A, B, C ou D. Mas, como um ou mais indicadores da distribuição espacial e de 
tamanho/declínio/fragmentação populacional da espécie estavam próximos de serem classificado como 
Vulnerável (VU), segundo alguns dos casos particulares definidos pela IUCN Standards and Petitions Committee 
(2019, p 76-77), a espécie foi finalmente classificada como Quase Ameaçada (NT) neste momento."
  }
  
  if (info$cat.reg.clean %in% "LC") {
    
    if (is.na(info$reduction_A12)) {
      
      if (info$AOO > 2000 & info$EOO > 20000) {
        texto <- "A espécie apresenta ampla distribuição geográfica (EOO > 20000 km2) e 
grande área de ocupação (AOO > 2000 km2). Por esse motivo, a espécie foi 
classificada como Menor Preocupação (LC) neste momento. Contudo, não foram 
foram encontrados dados para avaliar o tamanho e declínio das populações da 
espécie. Portanto, estudos populacionais específicos, assim como uma avaliação 
detalhada das pressões e ameaças que afetem a perpetuação da espécie na natureza, 
devem ser conduzidos no futuro para subsidiar uma avaliação de risco de extinção mais robusta."
      }

      if (info$AOO < 2000 & info$EOO > 20000) {
        if (!info$declineB %in% "Decreasing")
          texto <- "A espécie apresenta ampla distribuição geográfica (EOO > 20000 km2), 
mas uma área de ocupação relativamente pequena (AOO < 2000 km2). Porém, a espécie não 
apresentou evidências de declínio contínuo da quantidade de habitat dentro da sua extensão de ocorrência,
uma condição necessária para classificar espécies como ameaçada de extinção sob o critério B. 
Por esse motivo, a espécie foi classificada como Menor Preocupação (LC) neste momento. 
Contudo, não foram foram encontrados dados para avaliar o tamanho e declínio das populações da 
espécie. Portanto, estudos populacionais específicos, assim como uma avaliação 
detalhada das pressões e ameaças que afetem a perpetuação da espécie na natureza, 
devem ser conduzidos no futuro para subsidiar uma avaliação de risco de extinção mais robusta."

        if (info$declineB %in% "Decreasing" & (info$sever.frag | info$nbe_loc_total > 10))
          texto <- "A espécie apresenta ampla distribuição geográfica (EOO > 20000 km2), 
mas uma área de ocupação relativamente pequena (AOO < 2000 km2). Porém, a espécie não 
apresentou populações severamente fragmentadas ou ocorrências em um número de locais < 11, 
condições necessárias para classificar espécies como ameaçada de extinção sob o critério B. 
Por esse motivo, a espécie foi classificada como Menor Preocupação (LC) neste momento. 
Contudo, não foram foram encontrados dados para avaliar o tamanho e declínio das populações da 
espécie. Portanto, estudos populacionais específicos, assim como uma avaliação 
detalhada das pressões e ameaças que afetem a perpetuação da espécie na natureza, 
devem ser conduzidos no futuro para subsidiar uma avaliação de risco de extinção mais robusta."
        
      }

      if (info$AOO < 2000 & info$EOO < 20000) {
        texto <- "A espécie apresenta distribuição geográfica e área de ocupação pequenas 
(EOO < 20000 km2 e AOO < 2000 km2). Porém, a espécie não apresentou uma ou mais  
condições necessárias para classificar espécies como ameaçada de extinção sob o critério B 
(fragmentação severa, número de localidades pequena e/ou sem declínuo contíbuo de habitat). 
Por esse motivo, a espécie foi classificada como Menor Preocupação (LC) neste momento. 
Contudo, não foram foram encontrados dados para avaliar o tamanho e declínio das populações da 
espécie. Portanto, estudos populacionais específicos, assim como uma avaliação 
detalhada das pressões e ameaças que afetem a perpetuação da espécie na natureza, 
devem ser conduzidos no futuro para subsidiar uma avaliação de risco de extinção mais robusta."
      }
      
    } else {

      if (info$AOO > 2000 & info$EOO > 20000) {
        texto <- "A espécie apresenta ampla distribuição geográfica (EOO > 20000 km2) e 
grande área de ocupação (AOO > 2000 km2). A espécies também não apresentou 
evidências de fragmentação severa e ela ocorre em um grande número de localidades. 
Adicionalmente, não há evidencias de declínio populacional importante nas três últimas 
gerações (> 30%) ou de um tamanho populacional pequeno (<10.000 indivíduos adultos). 
Ou seja, a espécie foi avaliada contra os critérios A, B, C ou D e não foi encontrada
nenhuma indicação de ameaça. Por esse motivo, a espécie foi classificada, sem dúvidas, 
como Menor Preocupação (LC)."
      }
      
      if (info$AOO < 2000 & info$EOO > 20000) {
        texto <- "A espécie apresenta ampla distribuição geográfica (EOO > 20000 km2), 
mas uma área de ocupação relativamente pequena (AOO < 2000 km2). Porém, a espécie não 
apresentou evidências de fragmentação severa ou de ocorrência em poucas localidades (<11),
condições necessárias para classificar espécies como ameaçada de extinção sob o critério B. 
Adicionalmente, não há evidencias de declínio populacional importante nas três últimas 
gerações (> 30%) ou de um tamanho populacional pequeno (<10.000 indivíduos adultos). 
Ou seja, a espécie foi avaliada contra os critérios A, B, C ou D e não foi encontrada 
nenhuma indicação de ameaça. Por esse motivo, a espécie foi classificada como Menor Preocupação (LC)."
      }
      
      if (info$AOO < 2000 & info$EOO < 20000) {
        texto <- "A espécie apresenta distribuição geográfica e área de ocupação pequenas 
(EOO < 20000 km2 e AOO < 2000 km2). Porém, a espécie não apresentou uma ou mais  
condições necessárias para classificar espécies como ameaçada de extinção sob o critério B 
(fragmentação severa, número de localidades pequena e/ou sem declínuo contíbuo de habitat). 
Adicionalmente, não há evidencias de declínio populacional importante nas três últimas 
gerações (> 30%) ou de um tamanho populacional pequeno (<10.000 indivíduos adultos). 
Ou seja, a espécie foi avaliada contra os critérios A, B, C ou D e não foi encontrada 
nenhuma indicação de ameaça. Por esse motivo, a espécie foi classificada como Menor Preocupação (LC)."
      }
      
    }
    
  }
  
  if (info$cat.reg.clean %in% c("VU", "EN", "CR")) {
    
    critA <- grepl("A", info$main.criteria)
    critB <- grepl("B", info$main.criteria)
    critC <- grepl("C", info$main.criteria)
    critD <- grepl("D", info$main.criteria)

    ## Criando os textos base
    if (critA & !critB & !critC & !critD) { # only A
      texto <- c("Apesar das estimativas de tamanho populacional acima do limiar de 10,000 indíviduos adultos, 
foi estimada uma redução de _REDUCAO_ no tamanho das populações da espécie nas últimas três gerações. _EARLY_. _EXPLOIT_. 
_DECL_. _PROT_. _TYPE_. Assumindo, que as causas da perda de habitat na Mata Atlântica não cessaram, essa redução 
populacional permite classificar a espécie como _CAT_, sob _CRIT_. _AUX_. _REGION_.")
    }
    
    if (!critA & critB & !critC & !critD) { # only B
      
      texto <- c("A espécie apresentou _EOO_, _AOO_ e as condições necessárias 
para classificar espécies como ameaçada de extinção  sob o critério B (fragmentação 
severa ou ocorrência em poucos locais e declínuo contínuo). _EARLY_. _EXPLOIT_. _PROT_. _TYPE_. 
Assim, a espécie foi classificada como _CAT_, sob _CRIT_. _AUX_. _REGION_. _POPDEC_.") 

    }
    
    if (critA & critB & !critC & !critD) { # A and B
      texto <- c("A espécie apresentou _EOO_, _AOO_ e as condições necessárias 
para classificar espécies como ameaçada de extinção  sob o critério B (fragmentação 
severa ou ocorrência em poucos locais e declínuo contínuo). Adicionalmente, foi 
estimada uma redução de _REDUCAO_ no tamanho das populações da espécie nas últimas três gerações. 
_EARLY_. _EXPLOIT_. _PROT_. _TYPE_. Assim, a espécie foi classificada como _CAT_, sob _CRIT_. _AUX_. _REGION_.") 
      
    }
    
    if (critA & !critB & critC & !critD) { # A and C
      texto <- c("A espécie apresentou _EOO_, _AOO_ e as condições necessárias 
para classificar espécies como ameaçada de extinção  sob o critério B (fragmentação 
severa ou ocorrência em poucos locais e declínuo contínuo). Adicionalmente, foi 
estimado que a espécie possui uma população pequena (_POP_) e em declínio. _EARLY_. _EXPLOIT_. 
_DECL_. _PROT_. _TYPE_. Assim, a espécie foi classificada como _CAT_, sob _CRIT_. _AUX_. _REGION_.") 
      
    }
    
    if (critA & !critB & critC & critD) { # A, C and D
      texto <- c("Foi estimada uma redução de _REDUCAO_ no tamanho das populações da espécie nas 
últimas três gerações. Adicionalmente, foi estimado que a espécie possui uma população 
muito pequena (_POP_) e em declínio. _EARLY_. _EXPLOIT_. _DECL_. _PROT_. _TYPE_. Assim, a 
espécie foi classificada como _CAT_, sob _CRIT_. _AUX_. _REGION_.")
      
    }
    
    if (critA & critB & critC & !critD) { # A, B and C
      texto <- c("A espécie apresentou _EOO_, _AOO_ e as condições necessárias 
para classificar espécies como ameaçada de extinção  sob o critério B (fragmentação 
severa ou ocorrência em poucos locais e declínuo contínuo). Adicionalmente, foi 
estimada uma redução de _REDUCAO_ no tamanho das populações da espécie nas últimas três gerações. 
Também foi estimado que a espécie possui uma população pequena (_POP_) e em declínio. 
_EARLY_. _EXPLOIT_. _PROT_. _TYPE_. Assim, a espécie foi classificada como _CAT_, sob _CRIT_. _AUX_. _REGION_.") 
      
    }
    
    if (!critA & critB & critC & !critD) { # B and C (none in THREAT)
      texto <- ""
    }

    if (critA & critB & critC & critD) { # all criteria (none in THREAT)
      texto <- ""
    }

    ## Substituindo os campos específicos/variáveis
    
    # Espécies sem dados dos inventários
    if (grepl("_POPDEC_", texto)) {
      rep <- c("Não foi encontrado nenhum registro para a espécie entre as centenas de inventários 
florestais compilados para a Mata Atlântica, sugerindo a espécie seja relativamente rara ao longo se sua distribuição. 
As altas taxas de desmatamento e de degradação florestal às quais a Mata Atlântica foi submetida, 
também sugerem que as populações dessa espécie possam ter sofridos declínios importantes. 
Sendo assim, estudos populacionais específicos são necessários para melhor avaliar o risco 
e a categoria de ameaça de extinção para essa espécie")
      texto <- gsub("_POPDEC_", rep, texto)
    }
    
    # Population decline
    if (grepl("_REDUCAO_", texto)) {
      repos <- round(info$reduction_A12, 1)
      if (!is.na(repos)) {
        
        texto <- gsub("_REDUCAO_", paste0(repos, "%"), texto, fixed = TRUE)
      }  else {
        texto <- texto
      }
    }  

    # Correction in population decline for early successional species
    if (grepl("Correction of", info$reduction_A12_obs)) {

        repos <- c(" Espécie _EG_, cujas populações provavelmente se beneficiam das alterações na qualidade do habitat (e.g. efeitos de borda),
contrabalanceando os declínios de tamanho populacional causados pela perda de habitat. 
Este efeito da mudança na qualidade do habitat remanescente foi incorporado nas estimativas do declínio do tamanho populacional da espécie, mas de maneira 
simples e preliminar, havendo a necessidade de estudos mais detalhados 
para melhor prever a relação entre perda de habitat e declínio populacional ao longo da área de distribuição da espécie")
        
        EG <- ifelse(info$ecol.group %in% "early_secondary", 
                     "'Secundária Inicial'", "'Pioneira'")
        
        repos <- gsub("_EG_", EG, repos, fixed = TRUE)

        texto <- gsub("_EARLY_", repos, texto, fixed = TRUE)
        
    } else {
        texto <- gsub("_EARLY_\\.", "", texto, perl = TRUE)
    }

        
    # Uses and exploitation
    repos <- avalia$UseTradeDocumentation.value
    if (repos %in% c("", NA, "Não foram encontrados registros sobre usos efetivos ou potenciais.")) {
      # texto <- gsub(" _EXPLOIT_\\.", " Não foram econtrados registros sobre usos efetivos ou potenciais, como corte seletivo ou outros usos que afetem diretamente as populações da espécie.", 
      #               texto, perl = TRUE)
      texto <- gsub("_EXPLOIT_\\.", "", texto, perl = TRUE)
    } else {
      
      if (grepl("Extra", info$basis_d) | "3.2" %in% pesquisa[,1] | "3.3.1" %in% conserva[,1]) {
        repos <- c(" Espécie com histórico de exploração na Mata Atlântica. Esta exploração 
foi incorporada na avaliação do declínio do tamanho populacional da espécie, mas de maneira 
bastante preliminar e conservadora, havendo a necessidade de estudos populacionais detalhados 
para estimar qual é o real impacto dessa exploração sobre a diminuição e recuperação das 
populações da espécie ao longo da área de distribuição da espécie")
        texto <- gsub("_EXPLOIT_", repos, texto, fixed = TRUE)
      } else {
        repos <- c("A espécie possui registros de um ou mais usos pelo homem que podem afetar 
diretamente as populações da espécie")
        texto <- gsub("_EXPLOIT_", repos, texto, fixed = TRUE)
      }
    }
    
    # Threat category
    categ <- info$cat.reg.clean
    if (categ == "CR")
      texto <- gsub("_CAT_", "Criticamente em Perigo (CR) de extinção", texto, perl = TRUE)
    if (categ == "EN")
      texto <- gsub("_CAT_", "Em Perigo (EN) de extinção", texto, perl = TRUE)
    if (categ == "VU")
      texto <- gsub("_CAT_", "Vulnerável (VU) de extinção", texto, perl = TRUE)
    
    # IUCN criteria
    crits <- info$main.criteria
    crits <- gsub("A1+A2","A2", crits, fixed = TRUE)
    crits <- gsub("^A2$","A2", crits, perl = TRUE)
    crits <- gsub("D","D1", crits, fixed = TRUE)
    crits <- gsub("\\+",", ", crits, perl = TRUE)
    
    if (grepl(", ", crits)) {
      texto <- gsub("_CRIT_", paste0("os critérios ", crits), texto, fixed = TRUE)
    } else {
      texto <- gsub("_CRIT_", paste0("o critério ", crits), texto, fixed = TRUE)
    }
    
    # Auxiliary criteria
    aux.crit <- info$aux.criteria %in% c("", NA, " ")
    if (aux.crit) {
      texto <- gsub(" _AUX_\\.", "", texto, perl = TRUE)
    } else {
      rep <- paste0("A espécie também foi avaliada sob outros critérios. A categoria e critérios resultantes desta avaliação foram: ",
                    gsub("\\:", " sob o(s) critério(s)", 
                         gsub("\\+",", ", info$aux.criteria, perl = TRUE), perl = TRUE))
      texto <- gsub("_AUX_", rep, texto, perl = TRUE)
    }
      
    # Spatial indicators
    rep <- info$AOO
    if (!is.na(rep) & grepl("_AOO_", texto)) {
      texto <- gsub("_AOO_", paste0(rep, " km2 de AOO"), texto, fixed = TRUE)
    } else {
      texto <- gsub(", _AOO_", "", texto, fixed = TRUE)
    }
    
    rep <- info$EOO
    if (!is.na(rep) & grepl("_EOO_", texto)) {
      texto <- gsub("_EOO_", paste0(rep, " km2 de EOO"), texto, fixed = TRUE)
    } else {
      texto <- gsub("_EOO_,", "", texto, fixed = TRUE)
    }
    
    # Population size
    rep <- round(info$pop.size, 0)
    if (!is.na(rep) & grepl("_POP_", texto))
      texto <- gsub("_POP_", paste0(rep, " indivíduos adultos"), texto, fixed = TRUE)
    
    # Regional assessment
    categ <- info$cat.reg.clean
    categ.reg <- gsub("_PE","", info$category.regional)
    if (categ == categ.reg) {
      texto <- gsub(" _REGION_\\.", "", texto, perl = TRUE)
    } else {
      rep <- "Esta espécie possui distribuição mais ampla, ocorrendo com frequência 
em domínios vizinhos à Mata Atlântica (e.g. Caatinga, Cerrado ou Pantanal) ou 
preferencialmente em tipos de vegetacão não florestais dentro da Mata Atlântica. 
Sendo assim, a categoria de ameaça das populações da espécie na MAta Atlântica foi 
rebaixada em uma categoria, seguindo as recomendações da IUCN para avaliações regionais (IUCN 2012)."
      texto <- gsub("_REGION_", rep, texto, perl = TRUE)
    }

    # Declinio continuo
    rep <- info$declineB
    if (!is.na(rep) & grepl("_DECL_", texto)) {
      if (rep == "Decreasing")
        texto <- gsub("_DECL_", 
                      "Foi inferido que esteja havendo um declínio contínuo na área, extensão e qualidade do habitat ao qual a espécies está associada", 
                      texto, fixed = TRUE)
      if (rep == "Not Decreasing")
        texto <- gsub("_DECL_\\.", "", texto, perl = TRUE)
      
    } else {
      texto <- gsub("_DECL_\\.", "", texto, perl = TRUE)
    }
    
    # Ausente em Unidade de Proteção Integrals
    rep <- info$protected
    if (!is.na(rep) & rep == 0) {
      texto <- gsub("_PROT_", "Até a presente avaliação, não foram encontrados 
registros válidos da ocorrência da espécie em unidades de conservação de proteção integral, 
demandando portanto um planejamento de ações de conservação que visem garantir a 
perpetuação dessa espécie em seu habitat natural", texto, perl = TRUE)
    } else {
      texto <- gsub("_PROT_\\.", "", texto, perl = TRUE)
    }
    
    # Só tipos ou só registros velhos
    tipo <- info$only.type
    ultimo <- 2019 - info$last.record
    if ((!is.na(tipo) & tipo == TRUE) | (!is.na(ultimo) & ultimo >= 50)) {
      if ((!is.na(tipo) & tipo == TRUE) & (!is.na(ultimo) & ultimo >= 50))
        repos <- "Espécie conhecida apenas pelo seu espécimen tipo, que foi coletado há mais de 50 anos atrás."

      if (!(!is.na(tipo) & tipo == TRUE) & (!is.na(ultimo) & ultimo >= 50))
        repos <- "Espécie não registrada na natureza há mais de 50 anos atrás."

      if ((!is.na(tipo) & tipo == TRUE) & !(!is.na(ultimo) & ultimo >= 50))
        repos <- ""
      
      texto <- gsub("_TYPE_\\.", repos, texto, perl = TRUE)
      
    } else { 
      texto <- gsub("_TYPE_\\.", "", texto, perl = TRUE)
    }    
  }

  # assess$RedListRationale.value[i] <- stringr::str_squish(gsub("\\\n", "", texto))
  texto <- gsub("\\s+", " ", gsub("\\\n", "", texto, perl = TRUE), perl = TRUE)
  assess$RedListRationale.value[i] <- 
    gsub("^ | $", "", texto, perl = TRUE)
  
}

table(assess$RedListRationale.value %in% "", all.spp1$cat.reg.clean)
table(grepl("_", assess$RedListRationale.value))
assess$RedListRationale.value[grepl("_", assess$RedListRationale.value)]

## Reorganizing columns
names(sample)[!names(sample) %in% names(assess)]
# assess1 <- assess[, c(names(sample), "ThreatsDocumentation.value.PrevAssessment")]
assess1 <- assess[, c(names(sample))]
table(names(assess1) == names(sample))

head(sample, 2)
head(assess1, 2)

## Replacing NA by empty characters
for (i in seq_len(dim(assess1)[2])) {
  check_these <- assess1[,i] %in% c("", NA, NA_character_)
  if (any(check_these))
    assess1[check_these, i] <- ""
}

## Adding same columns from previous IUCN assessments
assess2 <- .merge_sis_connect(assess1, type = "assessments", tax.table = tax,
                              prev.assess, x.spp = prev.assess$species.correct2, 
                              replace = TRUE)

## Replacing NA by empty characters for the merged data frame
for (i in seq_len(dim(assess2)[2])) {
  check_these <- assess2[,i] %in% c("", NA, NA_character_)
  if (any(check_these))
    assess2[check_these, i] <- ""
}

## Saving
write.csv(assess1, "data/sis_connect/assessments_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(assess2, "data/sis_connect/assessments_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

###############################################################################H
###############################################################################H
# ALL FIELDS -------------------------------------------------------------------

## Generating the SIS CONNECT file "allfields.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/allfields.csv")
# head(sample, 3)
# apply(sample, 2, unique)

all.spp <- readRDS("data/all_spp_crits_tax_cites.rds")
allfields <- data.frame(internal_taxon_id = all.spp$internal_taxon_id)

## AOO
allfields$AOO.range <- all.spp$AOO
allfields$AOO.justification <- "A AOO foi estimada como sendo o número mínimo de células ocupadas em um total de 30 grids de 2×2 km grid com posições iniciais estabelecidas de maneira aleatória (Dauby et al. 2017)."
allfields$AOODetails.breedingRange <- NA
allfields$AOODetails.breedingRangeJustification <- NA
allfields$AOODetails.nonbreedingRange <- NA
allfields$AOODetails.nonbreedingRangeJustification <- NA
#Continuing decline in AOO: "1" (meaning true), "0" (meaning false), or "U" (meaning “unknown”)
allfields$AOOContinuingDecline.isInContinuingDecline <- ""
allfields$AOOContinuingDecline.qualifier <- "Unknown"
allfields$AOOContinuingDecline.justification <- ""
# allfields$AOOContinuingDecline.justification <- 
#   "Não foram feitas inferências sobre o declínuo contínuo da Área de Ocupação (AOO) da espécie. Porém, foi estimado o declínio contínuo da Área de Habitat (AOH) dentro a área de extensão terrestre da espécie."
allfields$AOOExtremeFluctuation.isFluctuating <- NA
allfields$AOOExtremeFluctuation.justification <- NA

## EOO
allfields$EOO.range <- all.spp$EOO
allfields$EOO.justification <- 
  "EOO foi definida como a área do polígono convexo mínimo que contém todas as ocorrências conhecidas para a espécies, incluindo áreas não-terrestres, como recomendado pelos critérios e categorias da IUCN"
# se EOO não for gerado ou for menor que AOO
check_these <- !is.na(allfields$EOO.range) &
                  (allfields$EOO.range %in% 0 | allfields$EOO.range < allfields$AOO.range)
allfields$EOO.range[check_these] <- allfields$AOO.range[check_these]
allfields$EOO.justification[check_these] <- 
  "EOO foi definida como sendo igual a AOO, pois o EOO estimado teve área igual a zero ou área menor que o AOO"

allfields$EOOContinuingDecline.isContinuingDecline <- ""
allfields$EOOContinuingDecline.qualifier <- "Unknown"
allfields$EOOContinuingDecline.justification <- ""
# allfields$EOOContinuingDecline.justification <- 
#   "Não foram feitas inferências sobre o declínuo contínuo da Área de Extensão (EOO) da espécie. Porém, foi estimado o declínio contínuo da Área de Habitat (AOH) dentro a área de extensão terrestre da espécie."
allfields$EOOExtremeFluctuation.isFluctuating <- NA
allfields$EOOExtremeFluctuation.justification <- NA
# allfields$EOODetails.breedingRange <- NA
# allfields$EOODetails.nonbreedingRange <- NA

## Locations
allfields$LocationsNumber.range <- all.spp$nbe_loc_total
allfields$LocationsNumber.justification <- 
  "O número de locais no qual a espécies ocorre foi estimado como a soma entre o número de Unidades de Conservação de Proteção Integral no qual a espécies ocorre e o número de células ocupadas for de Unidades de Conservação em um grid de 10 km de resolução, assumindo 10 km como uma escala espacial razoável da influêcia das ameaças à espécie (Dauby et al. 2017; Stévart et al. 2019)."
allfields$LocationContinuingDecline.inDecline <- ""
allfields$LocationContinuingDecline.qualifier <- "Unknown"
allfields$LocationContinuingDecline.justification <- ""
# allfields$LocationContinuingDecline.justification <- 
#   "Não foram feitas inferências sobre o declínuo contínuo do número de locais nos quais a espécie ocorre. Porém, foi estimado o declínio contínuo da Área de Habitat (AOH) dentro a área de extensão terrestre da espécie."
allfields$LocationExtremeFluctuation.isFluctuating <- NA
allfields$LocationExtremeFluctuation.justification <- NA

## Very Restricted
allfields$AreaRestricted.isRestricted <- 
  ifelse(all.spp$AOO <20, 1, ifelse(all.spp$nbe_loc_total <5, 1, 0)) 
allfields$AreaRestricted.isRestricted[is.na(allfields$AreaRestricted.isRestricted)] <- "U"
allfields$AreaRestricted.justification <- NA

## Elevation/Depth
occs <- readRDS("data/threat_species_by_country_environmental.rds")
data.table::setkeyv(occs, "species.correct2")
spp.ranges <- occs[ , range(ALT, na.rm = TRUE), by = species.correct2]
ranges <- as.data.frame.matrix(aggregate(spp.ranges$V1, list(spp.ranges$species.correct2), min))
names(ranges) <- c("species", "ElevationLower")
ranges$ElevationLower[!is.na(ranges$ElevationLower) & 
                        !is.infinite(ranges$ElevationLower) & 
                          ranges$ElevationLower < 0] <- 0
ranges$ElevationLower[is.infinite(ranges$ElevationLower)] <- NA
ranges$ElevationUpper <- aggregate(spp.ranges$V1, list(spp.ranges$species.correct2), max)$x
ranges$ElevationUpper[is.infinite(ranges$ElevationUpper)] <- NA
ranges <- dplyr::left_join(all.spp, ranges)
table(all.spp$species == ranges$species)

# Rounding
ranges$ElevationLowerClass <- cut(ranges$ElevationLower, breaks = seq(0,5000, by = 25))
levels(ranges$ElevationLowerClass) <- c(0, rep(seq(50,4950, by = 50), each = 2), 5000)
ranges$ElevationUpperClass <- cut(ranges$ElevationUpper, breaks = seq(0,5000, by = 25))
levels(ranges$ElevationUpperClass) <- c(0, rep(seq(50,4950, by = 50), each = 2), 5000)

replace_these <- !is.na(ranges$ElevationLowerClass) & !is.na(ranges$ElevationUpperClass) &
                    !is.na(ranges$nbe_loc_total) & ranges$nbe_loc_total > 10 & 
                      !ranges$ElevationLowerClass == ranges$ElevationUpperClass &
                        ranges$cat.reg.clean %in% c("LC","NT","VU","EN","CR")

allfields$ElevationLower.limit <- NA
allfields$ElevationLower.limit[replace_these] <- ranges$ElevationLowerClass[replace_these]
allfields$ElevationUpper.limit <- NA
allfields$ElevationUpper.limit[replace_these] <- ranges$ElevationUpperClass[replace_these]
allfields$DepthLower.limit <- NA
allfields$DepthUpper.limit <- NA
allfields$DepthZone.depthZone <- NA

## Population
allfields$PopulationSize.range <- all.spp$pop.size

allfields$CurrentTrendDataDerivation.value <- NA
replace_these <- !is.na(all.spp$reduction_A12)
allfields$CurrentTrendDataDerivation.value[replace_these] <- "Inferred"
replace_these <- is.na(all.spp$reduction_A12) & !is.na(all.spp$pop.size)
allfields$CurrentTrendDataDerivation.value[replace_these] <- "Suspected"

allfields$YearOfPopulationEstimate.value <- 2018
allfields$PopulationExtremeFluctuation.isFluctuating <- NA
allfields$PopulationExtremeFluctuation.justification <- NA

allfields$SevereFragmentation.isFragmented <- ifelse(all.spp$sever.frag, 1, 0) 
allfields$SevereFragmentation.justification <- 
  "Mais da metade da área de ocupação de uma determinada espécie sendo representada por subpopulações isoladas foi considerada como indicação de uma distribuição severamente fragmentada."

allfields$PopulationContinuingDecline.isDeclining <- ifelse(all.spp$any.decline == "Decreasing", 1, 0)
allfields$PopulationContinuingDecline.isDeclining[is.na(allfields$PopulationContinuingDecline.isDeclining)] <- "U"
allfields$PopulationContinuingDecline.qualifier <- ""
allfields$PopulationContinuingDecline.justification <- ""
replace_these <- !allfields$PopulationContinuingDecline.isDeclining %in% "U"
allfields$PopulationContinuingDecline.qualifier[replace_these] <- "Inferred" 
# allfields$PopulationContinuingDecline.justification <- 
#   "O tamanho da população foi considerada como em declínio continuo se houve uma diminuição de 1% ou mais no tamanho da população entre 2000 e 2018."

# critC <- readRDS("data\criterionC_optmim_params.rds")
critC <- readRDS("data/criterionC_all_prop_mature.rds")
critC <- dplyr::left_join(all.spp, 
                          critC[, c("species", 
                                    #"reduction_3gen", "reduction_2gen", "reduction_1gen",
                                    "prop.subpop.size", "max.subpop.size")])

reduc <- round(as.double(critC$reduction_1gen), 1)
reduc[!is.na(reduc) & reduc <= 0] <- 0
add_these <- !is.na(reduc)
allfields$PopulationDeclineGenerations1.range <- ""
allfields$PopulationDeclineGenerations1.range[add_these] <- reduc[add_these]
allfields$PopulationDeclineGenerations1.qualifier <- ""
allfields$PopulationDeclineGenerations1.justification <- ""
allfields$PopulationDeclineGenerations1.justification[add_these] <- "Inferred"

reduc <- round(as.double(critC$reduction_2gen), 1)
reduc[!is.na(reduc) & reduc <= 0] <- 0
add_these <- !is.na(reduc)
allfields$PopulationDeclineGenerations2.range <- ""
allfields$PopulationDeclineGenerations2.range[add_these] <- reduc[add_these]
allfields$PopulationDeclineGenerations2.qualifier <- ""
allfields$PopulationDeclineGenerations2.justification <- ""
allfields$PopulationDeclineGenerations2.justification[add_these] <- "Inferred"

reduc <- round(as.double(critC$reduction_3gen), 1)
reduc[!is.na(reduc) & reduc <= 0] <- 0
add_these <- !is.na(reduc)
allfields$PopulationDeclineGenerations3.range <- ""
allfields$PopulationDeclineGenerations3.range[add_these] <- reduc[add_these]
allfields$PopulationDeclineGenerations3.qualifier <- ""
allfields$PopulationDeclineGenerations3.justification <- ""
allfields$PopulationDeclineGenerations2.justification[add_these] <- "Inferred"

allfields$SubpopulationExtremeFluctuation.isFluctuating <- NA
allfields$SubpopulationExtremeFluctuation.justification <- NA
allfields$SubpopulationContinuingDecline.isDeclining <- NA
allfields$SubpopulationContinuingDecline.qualifier <- NA
allfields$SubpopulationContinuingDecline.justification <- NA

allfields$SubpopulationSingle.value <- ifelse(critC$prop.subpop.size == 100, 1, 0)

allfields$MaxSubpopulationSize.range <- round(as.double(critC$max.subpop.size))
allfields$SubpopulationNumber.range <- all.spp$Nbe_subPop
allfields$SubpopulationNumber.justification <- 
  "O número de subpopulações foi calculado usando o método de buffer circular e o quantil 90% da distribuição de 1/10 das distâncias máximas entre as ocorrências da espécie como um proxy da distância de dispersão média da espécie (Rivers et al. 2010; Dauby et al. 2017)"
allfields$MatureIndividualsSubpopulation.value <- NA

allfields$ExtinctionProbabilityGenerations3.justification <- NA
allfields$ExtinctionProbabilityGenerations3.range <- NA
allfields$ExtinctionProbabilityGenerations5.justification <- NA
allfields$ExtinctionProbabilityGenerations5.range <- NA
allfields$ExtinctionProbabilityYears100.justification <- NA
allfields$ExtinctionProbabilityYears100.range <- NA

## Past Reduction 
allfields$PopulationReductionPast.range <- abs(all.spp$reduction_A12)
allfields$PopulationReductionPast.direction <- 
  ifelse(all.spp$reduction_A12 > 0, "Reduction", "Increase")
allfields$PopulationReductionPast.qualifier <- ""
allfields$PopulationReductionPast.justification <- ""
check_these <- !is.na(allfields$PopulationReductionPast.direction)
allfields$PopulationReductionPast.justification[check_these] <- "Inferred"

allfields$PopulationReductionPastBasis.value <- NA
allfields$PopulationReductionPastBasis.value[!is.na(all.spp$reduction_A12)] <- 
  "b) an index of abundance appropriate for the taxon"
allfields$PopulationReductionPastBasis.value[grepl("Extra", all.spp$basis_d)] <- 
  "b) an index of abundance appropriate for the taxon|d) actual or potential levels of exploitation"

allfields$PopulationReductionPastBasis.detail <- NA
#as opções do SIS são AOO, EOO, Quality of Habitat?? Deixando vazio
# allfields$PopulationReductionPastBasis.detail[!is.na(all.spp$reduction_A12)] <- 
#   "Número de indivíduos adultos na população"

allfields$PopulationReductionPastReversible.value <- NA
allfields$PopulationReductionPastReversible.value[!is.na(all.spp$reduction_A12)] <- "U"
allfields$PopulationReductionPastUnderstood.value <- NA
allfields$PopulationReductionPastCeased.value <- NA
allfields$PopulationReductionPastCeased.value[!is.na(all.spp$reduction_A12)] <- "0"
check_these <- all.spp$declineB %in% "Not Decreasing" & allfields$PopulationReductionPastCeased.value %in% "0"
allfields$PopulationReductionPastCeased.value[check_these] <- "1"

## Future Reduction
allfields$PopulationReductionFuture.range <- NA
allfields$PopulationReductionFuture.direction <- NA
allfields$PopulationReductionFuture.qualifier <- NA
allfields$PopulationReductionFuture.justification <- NA
allfields$PopulationReductionFutureBasis.value <- NA
allfields$PopulationReductionFutureBasis.detail <- NA

## Ongoing Reduction
allfields$PopulationReductionPastandFuture.range <- NA
allfields$PopulationReductionPastandFuture.direction <- NA
allfields$PopulationReductionPastandFuture.numYears <- NA
allfields$PopulationReductionPastandFuture.qualifier <- NA
allfields$PopulationReductionPastandFuture.justification <- NA
allfields$PopulationReductionPastandFutureBasis.value <- NA
allfields$PopulationReductionPastandFutureBasis.detail <- NA
allfields$PopulationReductionPastandFutureReversible.value <- NA
allfields$PopulationReductionPastandFutureUnderstood.value <- NA
allfields$PopulationReductionPastandFutureCeased.value <- NA

## Habitat decline
allfields$HabitatContinuingDecline.isDeclining <- ifelse(all.spp$declineB == "Decreasing", "YES", "NO")
check_these <- is.na(all.spp$declineB)
allfields$HabitatContinuingDecline.isDeclining[check_these] <- "U"
allfields$HabitatContinuingDecline.qualifier <- ""
allfields$HabitatContinuingDecline.qualifier[!check_these] <- "Estimated"
allfields$HabitatContinuingDecline.justification <- ""
# allfields$HabitatContinuingDecline.justification[!check_these] <- "Estimated"
# allfields$HabitatContinuingDecline.justification <- 
#   "O declínuo contínuo foi estimado como a perda de área de habitat dentro da parte terrestre da EOO da espécie (i.e. extensão de habitat adequado) para o período entre 2000 e 2018."

## Life history
hab2 <- read.csv("data/threat_habitats.csv", as.is = TRUE, encoding = "UTF-8")
table(all.spp$species == hab2$Name_submitted)

allfields$GenerationLength.range <- hab2$GenerationLength.range
allfields$GenerationLength.justification <- hab2$GenerationLength.justification
allfields$GenerationLength.quality <- NA # good? what are the other options?
allfields$FemaleMaturityAge.age <- NA
allfields$FemaleMaturityAge.units <- NA
allfields$MaleMaturityAge.age <- NA
allfields$MaleMaturityAge.units <- NA
allfields$FemaleMaturitySize.size <- hab2$FemaleMaturitySize.size
allfields$MaleMaturitySize.size <- hab2$MaleMaturitySize.size
allfields$Longevity.longevity <- NA
allfields$Longevity.units <- NA
allfields$AvgReproductiveAge.age <- NA
allfields$AvgReproductiveAge.units <- NA
allfields$MaxSize.size <- hab2$MaxSize.size
allfields$BirthSize.size <- NA
allfields$GestationTime.time <- NA
allfields$GestationTime.units <- NA
allfields$ReproductivePeriodicity.value <- NA
allfields$AvgAnnualFecundity.fecundity <- NA
allfields$PopulationIncreaseRate.narrative <- NA
allfields$NaturalMortality.value <- NA
allfields$EggLaying.layEggs <- NA
allfields$LiveBirth.liveBirth <- NA
allfields$Parthenogenesis.exhibitsParthenogenesis <- NA
allfields$FreeLivingLarvae.hasStage <- NA
allfields$WaterBreeding.value <- NA

## Movement Patterns
allfields$MovementPatterns.pattern <- NA
allfields$Congregatory.value <- NA

## Use and Trade Information
allfields$CropWildRelative.isRelative <- "0"
crops <- c("Hevea brasiliensis", "Manihot esculenta", "Schinus terebinthifolia",
           "Psidium guajava", "Pouteria caimito", "Ilex paraguariensis",
           "Swietenia macrophylla", "Anacardium occidentale", 
           "Araucaria angustifolia")
allfields$CropWildRelative.isRelative[all.spp$species %in% crops] <- "1"
allfields$NotUtilized.isNotUtilized <- NA

usos2 <- read.csv("data/sis_connect/usetrade_threat.csv", encoding = "UTF-8")
allfields$UseTradeNoInformation.value <- "false"
check_these <- all.spp$internal_taxon_id %in% unique(usos2$internal_taxon_id)
allfields$UseTradeNoInformation.value[check_these] <- "true"

allfields$UTLocalLivelihood.subsistence <- 3
check_these <- all.spp$internal_taxon_id %in% unique(usos2$internal_taxon_id)
allfields$UTLocalLivelihood.subsistence[check_these] <- 2
check_these <- all.spp$internal_taxon_id %in% unique(usos2$internal_taxon_id[usos2$UTEndUse.UTEndUseSubfield.subsistence %in% "true"])
allfields$UTLocalLivelihood.subsistence[check_these] <- 1

allfields$UTLocalLivelihood.subsistenceRationale <- NA
allfields$UTLocalLivelihood.localcommercial <- NA
allfields$UTLocalLivelihood.localcommercialDetail <- NA

allfields$UTNatlCommercialValue.value <- NA
check_these <- all.spp$internal_taxon_id %in% unique(usos2$internal_taxon_id[usos2$UTEndUse.UTEndUseSubfield.national %in% "true"])
allfields$UTNatlCommercialValue.value[check_these] <- 1

allfields$UTIntlCommercialValue.value <- NA
check_these <- all.spp$internal_taxon_id %in% unique(usos2$internal_taxon_id[usos2$UTEndUse.UTEndUseSubfield.international %in% "true"])
allfields$UTIntlCommercialValue.value[check_these] <- 1

allfields$UTCaptiveHarvest.value <- NA
allfields$TrendInWildOfftake.value <- NA
allfields$TrendInDomesticOfftake.value <- NA
allfields$UTHarvestTrendComments.narrative <- NA
allfields$NonConsumptiveUse.isNonConsumptiveUse <- NA
allfields$NonConsumptiveUseDescription.narrative <- NA
allfields$Livelihoods.noInfo <- NA

## Coded Threats
threats <- read.csv("data/sis_connect/threats_threat.csv", encoding = "UTF-8")

allfields$NoThreats.noThreats <- "false"
check_these <- !all.spp$internal_taxon_id %in% unique(threats$internal_taxon_id)
allfields$NoThreats.noThreats[check_these] <- "true"
allfields$ThreatsUnknown.value <- "false"

## In-Place Conservation
allfields$InPlaceResearchRecoveryPlan.value <- NA # ver com CNCFlora
allfields$InPlaceResearchRecoveryPlan.note <- NA
allfields$InPlaceResearchMonitoringScheme.value <- NA # ver com CNCFlora
allfields$InPlaceResearchMonitoringScheme.note <- NA
allfields$InPlaceLandWaterProtectionSitesIdentified.value <- NA
allfields$InPlaceLandWaterProtectionSitesIdentified.note <- NA
#Occur in at least one PA
allfields$InPlaceLandWaterProtectionInPA.value <- 
  ifelse (all.spp$protected > 0, "Yes", "No")
allfields$InPlaceLandWaterProtectionInPA.note <- 
  "Presença em áreas protegidas se referem apenas aos registros válidos em Unidades de Conservação de Proteção Integral (categorias I, II ou III do 'IUCN Protected Area Management Categories' - www.iucn.org/theme/protected-areas)"
#Percentage of population protected by PAs (0-100)
allfields$InPlaceLandWaterProtectionPercentProtected.value <- round(all.spp$protected, 1)
allfields$InPlaceLandWaterProtectionPercentProtected.note <- 
  "A percentagem se refere ao número de registros 'de ocorrência dentro de Unidades de Conservação de Proteção Integral em relação total de registros da espécie, e não necessariamente à percentagem da população da espécie"
allfields$InPlaceLandWaterProtectionAreaPlanned.value <- NA
allfields$InPlaceLandWaterProtectionAreaPlanned.note <- NA
allfields$InPlaceLandWaterProtectionInvasiveControl.value <- NA
allfields$InPlaceLandWaterProtectionInvasiveControl.note <- NA
allfields$InPlaceSpeciesManagementHarvestPlan.value <- NA # ver com CNCFlora
allfields$InPlaceSpeciesManagementHarvestPlan.note <- NA
allfields$InPlaceSpeciesManagementReintroduced.value <- NA
allfields$InPlaceSpeciesManagementReintroduced.note <- NA
allfields$InPlaceSpeciesManagementExSitu.value <- NA
allfields$InPlaceSpeciesManagementExSitu.note <- NA
allfields$InPlaceEducationSubjectToPrograms.value <- NA
allfields$InPlaceEducationSubjectToPrograms.note <- NA
allfields$InPlaceEducationInternationalLegislation.value <- NA
allfields$InPlaceEducationInternationalLegislation.note <- NA

allfields$InPlaceEducationControlled.value <- "No"
spp <- conservation$internal_taxon_id[conservation$ConservationActions.ConservationActionsSubfield.ConservationActionsLookup %in% "5.4.1"]
check_these <- all.spp$internal_taxon_id %in% unique(spp)
allfields$InPlaceEducationControlled.value[check_these] <- "Yes"

allfields$InPlaceEducationControlled.note <- NA
allfields$InPlaceEducationControlled.note[check_these] <- 
  "Espécie incluída nos Anexos do CITES"

## Description of the species
# Take description from the taxonomy file?
table(spp.info$species.correct2 == all.spp$species)
allfields$IdentificationDescription.text <- NA

## Reorganizing columns
names(sample)[!names(sample) %in% names(allfields)]
names(allfields)[!names(allfields) %in% names(sample)]
allfields <- allfields[ , names(sample)]
table(names(allfields) == names(sample))

## Replacing NA by empty characters
for (i in seq_len(dim(allfields)[2])) {
  check_these <- is.na(allfields[,i])
  if (any(check_these))
    allfields[check_these, i] <- ""
}

## Adding same columns from previous IUCN assessments
allfields1 <- .merge_sis_connect(allfields, type = "allfields", tax.table = tax,
                              prev.assess, x.spp = prev.assess$species.correct2, 
                              replace = TRUE)

## Replacing NA by empty characters for the merged data frame
for (i in seq_len(dim(allfields1)[2])) {
  check_these <- allfields1[,i] %in% c("", NA, NA_character_)
  if (any(check_these))
    allfields1[check_these, i] <- ""
}


## Saving the SIS Connet files
write.csv(allfields, "data/sis_connect/allfields_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(allfields1, "data/sis_connect/allfields_threat_merged_iucn.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

###############################################################################H
###############################################################################H
# REFERENCES ----------------------------------------------------------------------

## Generating the SIS CONNECT file "references.csv"
## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/references.csv", encoding = "UTF-8")
head(sample, 3)
apply(sample[,1:10], 2, unique)

## Finding citations in each SIS Connect file
my.files <- list.files("data/sis_connect/", full.names = TRUE)
# files with references in the notes
keep <- "allfields|assessments|prev_assessment|researchneeded" 
# maybe depending on the content of the 'note' or 'desc' fields
maybe <- "conservationneeded|taxonomy|threats"
my.files1 <- my.files[grepl(keep, my.files, perl = TRUE)]

result.all <- vector("list", length(my.files1))
for (i in seq_len(length(my.files1))) {
  if (grepl("csv", my.files1[i])) { 
    arquivo <- read.csv(my.files1[i], encoding = "UTF-8")
  } else {
    arquivo <- readRDS(my.files1[i])
  } 
  
  planilha <- gsub('data/sis_connect/', '', my.files1[i], fixed = TRUE)
  planilha <- gsub('_threat.csv|_threat.rds', '', planilha, perl = TRUE)

  result1.spp <- vector("list", dim(arquivo)[1])
  names(result1.spp) <- arquivo$internal_taxon_id
    
  patt <- "et al\\.|in prep\\.|[1-2][0-9][0-9][0-9]\\)"
    
  for(j in 1:dim(arquivo)[1]) {
      linha <- arquivo[j, ]
      sp.i <- linha[, "internal_taxon_id"]
      if (any(grepl(patt, linha, perl = TRUE))) {
        indices <- which(grepl(patt, linha, perl = TRUE))
        res <- linha[, indices]
        res.split <- strsplit(as.character(res), "(?<=[()])", perl = TRUE)
        res.split1 <- lapply(res.split, function (x) x[grepl(patt, x, perl = TRUE)])
        
        res.split2 <- strsplit(as.character(res), " ", perl = TRUE)[[1]]
        indices2 <- which(grepl("\\([1-2][0-9][0-9][0-9]\\)", res.split2, perl = TRUE))
        if (length(indices2) > 0) {
          res.split2 <- res.split2[sort(c(indices2, indices2 - 1))]
          res.split1[[1]] <- c(res.split1[[1]], res.split2)
        }
        
        # et al.
        patt1 <- "et al\\. \\($"
        correcao1 <- 
          sapply(res.split1, function(x) any(grepl("et al\\. \\($", x, perl = TRUE)))
        if (any(correcao1)) {
          
          for (w in seq_len(sum(correcao1))) {
            input <- res.split1[which(correcao1)[w]][[1]]
            correcao.i <- unlist(lapply(input, function(x) grepl(patt1, x, perl = TRUE)))
            input[correcao.i] <- lapply(input[correcao.i], 
                                        function (x) substring(x, head(tail(gregexpr(" ", x)[[1]], 4),1)))
            # input <- list(stringr::str_squish(unlist(input)))
            input <- list(gsub("^ | $", "", gsub("\\s+", " ", unlist(input), perl = TRUE), perl = TRUE))
            
            patt2 <- "in prep\\.\\)$|^[1-2][0-9][0-9][0-9]\\)$|^\\([1-2][0-9][0-9][0-9]\\)$"
            if (any(sapply(input, function(x) grepl(patt2, x, perl = TRUE)))) {
              input2 <- lapply(input, 
                               function(x) paste0(x[which(grepl(patt2, x, perl = TRUE))-1], 
                                                  x[grepl(patt2, x, perl = TRUE)]))
              input2 <- lapply(input2, 
                               function(x) gsub("(\\))([1-2][0-9][0-9][0-9]\\)$)", "\\1", x, perl = TRUE))
              input <- append(input, input2)
            }
            
            input <- list(gsub("^ | $", "", gsub("\\s+", " ", unlist(input), perl = TRUE), perl = TRUE))
            res.split1[which(correcao1)[w]] <- input
            # res.split1[which(correcao1)[w]] <- list(stringr::str_squish(unlist(input)))
            
          }  
        }
        
        # autor (data)
        patt1 <- "^[1-2][0-9][0-9][0-9]\\)$"
        correcao2 <- 
          sapply(res.split1, function(x) any(!grepl("[0-9]", x, perl = TRUE)))
        if (any(correcao2)) {
          
          for (w in seq_len(sum(correcao2))) {
            input <- res.split1[which(correcao2)[w]][[1]]
            # correcao.i <- unlist(lapply(input, function(x) grepl(patt1, x, perl = TRUE)))
            # input[correcao.i] <- lapply(input[correcao.i], 
            #                             function (x) substring(x, head(tail(gregexpr(" ", x)[[1]], 4),1)))
            # input <- list(stringr::str_squish(unlist(input)))
            input <- list(gsub("^ | $", "", gsub("\\s+", " ", unlist(input), perl = TRUE), perl = TRUE))
            
            
            patt2 <- "in prep\\.\\)$|^[1-2][0-9][0-9][0-9]\\)$|^\\([1-2][0-9][0-9][0-9]\\)$"
            if (any(sapply(input, function(x) grepl(patt2, x, perl = TRUE)))) {
              input2 <- lapply(input, 
                               function(x) paste(x[which(grepl(patt2, x, perl = TRUE))-1], 
                                                  x[grepl(patt2, x, perl = TRUE)]))
              input2 <- lapply(input2, 
                               function(x) gsub("(\\))( [1-2][0-9][0-9][0-9]\\)$)", "\\1", x, perl = TRUE))
              input <- append(input, input2)
            }
            
            input <- list(gsub("^ | $", "", gsub("\\s+", " ", unlist(input), perl = TRUE), perl = TRUE))
            res.split1[which(correcao2)[w]] <- input
            # res.split1[which(correcao2)[w]] <- list(stringr::str_squish(unlist(input)))
            
          }  
        }
  
        #Cleaning
        res.split1 <- lapply(res.split1, 
                             function (x) x[!grepl("et al\\. \\($", x, perl = TRUE)])
        res.split1 <- lapply(res.split1, 
                             function (x) x[!x %in% "in prep.)"])
        res.split1 <- lapply(res.split1, 
                             function (x) x[!grepl("^[1-2][0-9][0-9][0-9]\\)$", x, perl = TRUE)])
        res.split1 <- lapply(res.split1, 
                             function (x) x[!grepl("^[1-2][0-9][0-9][0-9]$", x, perl = TRUE)])
        
        res.split1 <- lapply(res.split1, 
                             function (x) gsub("\\)$", "", x, perl = TRUE))
        res.split1 <- lapply(res.split1, 
                             function (x) gsub(" \\(", " ", x, perl = TRUE))
        res.split1 <- lapply(res.split1, 
                             function (x) gsub("sensu ", "", x, fixed = TRUE))
        res.split1 <- lapply(res.split1, 
                             function (x) gsub(", p\\.[0-9][0-9]", "", x, perl = TRUE))
        res.split1 <- lapply(res.split1, 
                             function (x) gsub("^in prep\\.\\) ", "", x, perl = TRUE))
        res.split1 <- lapply(res.split1, 
                             function (x) x[grepl("[0-9]|in prep", x, perl = TRUE)])
        res.split1 <- lapply(res.split1, 
                             function (x) x[!grepl("^\\([1-2][0-9][0-9][0-9]$", x, perl = TRUE)])
        
        #splitting and removing duplicates
        res.split1 <- lapply(res.split1, strsplit, "; |, ", perl = TRUE)
        
        # res2 <- paste0(unique(sort(stringr::str_squish(unlist(res.split1)))),
        #                collapse = "|")
        tmp <- gsub("^ | $", "", 
                    gsub("\\s+", " ", unlist(res.split1), perl = TRUE), perl = TRUE)
        res2 <- paste0(unique(tmp), collapse = "|")
        
        resultado <- cbind(planilha, internal_taxon_id = sp.i, refs = res2)
        result1.spp[[j]] <- resultado
      }
    }  
      
  result1 <- result1.spp[lengths(result1.spp) != 0]
  result1.1 <- do.call("rbind", result1)
  names(result.all)[i] <- planilha
  result.all[[i]] <- result1.1   
}  
lapply(result.all, head, 2)


# ASSESSMENTS
plan <- "allfields"
unique(result.all[[plan]][,"refs"])
#"Dauby et al. 2017|IUCN 2019|Lima et al. in prep.|Stévart et al. 2019|Rivers et al. 2010"
 
ref1 <- data.frame(Reference_type = "Assessment", 
                   author = "Dauby, G., Stévart, T., Droissart, V., Cosiaux, A., Deblauwe, V., Simo-Droissart, M., Sosef, M.S.M., Lowry, P.P., Schatz, G.E., Gereau, R.E., & Couvreur, T.L.P.",
                   year = "2017", volume = "7", number = "24", pages = "11292–11303",
                   title = "ConR: An R package to assist large-scale multispecies preliminary conservation assessments using distribution data", 
                   secondary_title = "Ecology and Evolution", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Dauby et al. 2017", result.all[[plan]][,"refs"])]
ref1 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref1, simplify = FALSE))
ref1$internal_taxon_id <- spp

ref2 <- data.frame(Reference_type = "Assessment", 
                   author = "IUCN Standards and Petitions Committee",
                   year = "2019", pages = "113p.",
                   title = "Guidelines for Using the IUCN Red List Categories and Criteria (Version 14)", 
                   secondary_title = "Prepared by the Standards and Petitions Committee", type = "book")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("IUCN 2019", result.all[[plan]][,"refs"])]
ref2 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref2, simplify = FALSE))
ref2$internal_taxon_id <- spp

ref3 <- data.frame(Reference_type = "Assessment", 
                   author = "Lima, R.A.F., Dauby, G., de Gasper, A.L., Vibrans, A.C., Oliveira, A.A., Prado, P.I., Souza, V.C., Siqueira, M.F., ter Steege, H.",
                   year = "in prep.", volume = NA, number = NA, pages = NA,
                   title = "The Atlantic Forest trees: a flora on the verge of extinction", 
                   secondary_title = NA, type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Lima et al. in prep.", result.all[[plan]][,"refs"])]
ref3 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref3, simplify = FALSE))
ref3$internal_taxon_id <- spp

ref4 <- data.frame(Reference_type = "Assessment", 
                   author = "Rivers, M.C., Bachman, S.P., Meagher, T.R., Nic Lughadha, E., Brummitt, N.A.",
                   year = "2010", volume = "19", number = "7", pages = "2071–2085",
                   title = "Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data", 
                   secondary_title = "Biodiversity and Conservation", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Rivers et al. 2010", result.all[[plan]][,"refs"])]
ref4 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref4, simplify = FALSE))
ref4$internal_taxon_id <- spp

ref5 <- data.frame(Reference_type = "Assessment", 
                  author = "Stévart, T., Dauby, G., Lowry, P.P., Blach-Overgaard, A., Droissart, V., Harris, D.J., Mackinder, B.A., Schatz, G.E., Sonké, B., Sosef, M.S.M., Svenning, J.-C., Wieringa, J.J., Couvreur, T.L.P.",
                  year = "2019", volume = "5", number = "11", pages = "eaax9444",
                  title = "A third of the tropical African flora is potentially threatened with extinction", 
                  secondary_title = "Science Advances", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
          grepl("Stévart et al. 2019", result.all[[plan]][,"refs"])]
ref5 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref5, simplify = FALSE))
ref5$internal_taxon_id <- spp


# ASSESSMENTS
plan <- "assessments"
head(unique(unlist(sapply(result.all[[plan]][,"refs"], strsplit, "\\|"))), 10)
#"IUCN 2012","Flora do Brasil 2020", "Lima et al. in prep.", 
#"Argentina & WWF 2017","SOS Mata Atlântica e INPE 2018",
#"Rosa et al. 2021","Ribeiro et al. 2009"

#"Lima et al. 2020a
ref6 <- data.frame(Reference_type = "Assessment", 
                   author = "Lima, R.A.F., Souza, V.C., de Siqueira, M.F. & ter Steege, H.",
                   year = "2020a", volume = "252", number = NA, pages = "108825",
                   title = "Defining endemism levels for biodiversity conservation: Tree species in the Atlantic Forest hotspot", 
                   secondary_title = "Biological Conservation", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Lima et al. 2020a", result.all[[plan]][,"refs"])]
ref6 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref6, simplify = FALSE))
ref6$internal_taxon_id <- spp

#Lima et al. 2020b
ref7 <- data.frame(Reference_type = "Assessment", 
                   author = "Lima, R.A.F., Oliveira, A.A., Pitta, G.R., de Gasper, A.L., Vibrans, A.C., Chave, J., ter Steege, H., Prado, P.I.",
                   year = "2020b", volume = "11", number = "1", pages = "6347",
                   title = "The erosion of biodiversity and biomass in the Atlantic Forest biodiversity hotspot", 
                   secondary_title = "Nature Communications", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Lima et al. 2020b", result.all[[plan]][,"refs"])]
ref7 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref7, simplify = FALSE))
ref7$internal_taxon_id <- spp

#"IUCN 2012"
ref8 <- data.frame(Reference_type = "Assessment", 
                   author = "IUCN",
                   year = "2012", volume = NA, number = NA, pages = NA,
                   title = "Guidelines for Application of Iucn Red List Criteria At Regional and National Levels", 
                   secondary_title = "IUCN, International Union for Conservation of Nature", type = "book")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("IUCN 2012", result.all[[plan]][,"refs"])]
ref8 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref8, simplify = FALSE))
ref8$internal_taxon_id <- spp

#"Flora do Brasil 2020"
ref9 <- data.frame(Reference_type = "Assessment", 
                   author = "Brazil Flora Group",
                   year = "2020", volume = NA, number = NA, pages = "doi:10.15468/1mtkaw",
                   title = "Brazilian Flora 2020 project (v393.285)", 
                   secondary_title = "Instituto de Pesquisas Jardim Botanico do Rio de Janeiro", type = "electronic source")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Flora do Brasil 2020", result.all[[plan]][,"refs"])]
ref9 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref9, simplify = FALSE))
ref9$internal_taxon_id <- spp

#"Lima et al. in prep.", 
ref10 <- data.frame(Reference_type = "Assessment", 
                    author = "Lima, R.A.F., Dauby, G., de Gasper, A.L., Vibrans, A.C., Oliveira, A.A., Prado, P.I., Souza, V.C., Siqueira, M.F., ter Steege, H.",
                    year = "in prep.", volume = NA, number = NA, pages = NA,
                    title = "The Atlantic Forest trees: a flora on the verge of extinction", 
                    secondary_title = NA, type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Lima et al. in prep.", result.all[[plan]][,"refs"])]
ref10 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref10, simplify = FALSE))
ref10$internal_taxon_id <- spp

#"Argentina & WWF 2017"
ref11 <- data.frame(Reference_type = "Assessment", 
                   author = "Fundación Vida Silvestre Argentina & WWF",
                   year = "2017", volume = NA, number = NA, pages = NA,
                   title = "State of the Atlantic Forest: Three countries, 148 million people, one of the richest forests on Earth", 
                   secondary_title = "Fundación Vida Silvestre Argentina, WWF-Brazil and WWF-Paraguay, Puerto Iguazú.", type = "report")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Argentina & WWF 2017", result.all[[plan]][,"refs"])]
ref11 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref11, simplify = FALSE))
ref11$internal_taxon_id <- spp

#"SOS Mata Atlântica e INPE 2018"
ref12 <- data.frame(Reference_type = "Assessment", 
                   author = "Fundação SOS Mata Atlântica & Instituto Nacional de Pesquisas Espaciais (INPE)",
                   year = "2018", volume = NA, number = NA, pages = NA,
                   title = "Atlas dos remanescentes florestais da Mata Atlântica: período 2016–2017", 
                   secondary_title = "Fundação SOS Mata Atlântica, São Paulo", 
                   type = "report")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("SOS Mata Atlântica e INPE 2018", result.all[[plan]][,"refs"])]
ref12 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref12, simplify = FALSE))
ref12$internal_taxon_id <- spp

#"Rosa et al. 2021"
ref13 <- data.frame(Reference_type = "Assessment", 
                   author = "Rosa, M.R., Brancalion, P.H.S., Crouzeilles, R., Tambosi, L.R., Piffer, P.R., Lenti, F.E.B., Hirota, M.M., Santiami, E., Metzger, J.P.",
                   year = "2021", volume = "7", pages = "eabc4547",
                   title = "Hidden destruction of older forests threatens Brazil’s Atlantic Forest and challenges restoration programs", 
                   secondary_title = "Science Advances", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Rosa et al. 2021", result.all[[plan]][,"refs"])]
ref13 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref13, simplify = FALSE))
ref13$internal_taxon_id <- spp

#"Ribeiro et al. 2009"
ref14 <- data.frame(Reference_type = "Assessment", 
                   author = "Ribeiro, M.C., Metzger, J.P., Martensen, A.C., Ponzoni, F.J., Hirota, M.M.",
                   year = "2009", volume = "142", number = "6", pages = "1141–1153",
                   title = "The Brazilian Atlantic Forest: How much is left, and how is the remaining forest distributed? Implications for conservation", 
                   secondary_title = "Biological Conservation", type = "journal article")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("Ribeiro et al. 2009", result.all[[plan]][,"refs"])]
ref14 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref14, simplify = FALSE))
ref14$internal_taxon_id <- spp

#IUCN 2019
ref15 <- data.frame(Reference_type = "Assessment", 
                   author = "IUCN Standards and Petitions Committee",
                   year = "2019", pages = "113p.",
                   title = "Guidelines for Using the IUCN Red List Categories and Criteria (Version 14)", 
                   secondary_title = "Prepared by the Standards and Petitions Committee", type = "book")
spp <- result.all[[plan]][,"internal_taxon_id"][
  grepl("IUCN 2019", result.all[[plan]][,"refs"])]
ref15 <- do.call(rbind.data.frame, 
                replicate(length(spp), ref15, simplify = FALSE))
ref15$internal_taxon_id <- spp


# PREVIOUS ASSESSMENTS
plan <- "prev_assessments"
head(unique(unlist(sapply(result.all[[plan]][,"refs"], strsplit, "\\|"))), 10)
#Sem citações para presente avaliação 

# RESEARCH NEEDED
plan <- "researchneeded"
unique(result.all[[plan]][,"refs"])
#Sem citações para presente avaliação 

# USES 
ref.usos <-  read.csv(my.files[grepl("usetrade", my.files)][1], encoding = "UTF-8")
sort(unique(ref.usos$source))
head(ref.usos)

#wikipedia
wiki <- data.frame(Reference_type = "Assessment", 
                   author = "Wikipedia",
                   year = "2019", #volume = 142, number = 6, pages = "1141–1153",
                   title = ": The free encyclopedia", 
                   secondary_title = "https://en.wikipedia.org/wiki/Main_Page", 
                   type = "electronic source")
spp <- ref.usos$internal_taxon_id[grepl("wikipedia", ref.usos$source)]
wiki <- do.call(rbind.data.frame, 
                 replicate(length(spp), wiki, simplify = FALSE))
wiki$internal_taxon_id <- spp

#Agroforestree Database
agro <- data.frame(Reference_type = "Assessment", 
                   author = "Orwa, C., Mutua, A., Kindt, R., Jamnadass, R., Simons, A.",
                   year = "2009", #volume = 142, number = 6, pages = "1141–1153",
                   title = "Agroforestree Database: a tree reference and selection guide", 
                   secondary_title = "World Agroforestry Centre, Kenya.", 
                   type = "electronic source")
spp <- ref.usos$internal_taxon_id[grepl("Agroforestree", ref.usos$source)]
agro <- do.call(rbind.data.frame, 
                replicate(length(spp), agro, simplify = FALSE))
agro$internal_taxon_id <- spp

#Brazilian Woods database
br.woods <- data.frame(Reference_type = "Assessment", 
                          author = "Coordenação do Laboratório de Produtos Florestais",
                          year = "2020", #volume = 142, number = 6, pages = "1141–1153",
                          title = "Brazilian Woods Database", 
                          secondary_title = "Available at: https://lpf.florestal.gov.br/en-us/brazilian-woods", 
                          type = "electronic source")
spp <- ref.usos$internal_taxon_id[grepl("Brazilian Woods database", ref.usos$source)]
br.woods <- do.call(rbind.data.frame, 
                replicate(length(spp), br.woods, simplify = FALSE))
br.woods$internal_taxon_id <- spp
      
#Carvalho, P.E.R. (1998)
car98 <- data.frame(Reference_type = "Assessment", 
                       author = "Carvalho, P.E.R.",
                       year = "1998", #volume = 142, number = 6, pages = "1141–1153",
                       title = "Especies nativas para fins produtivos", 
                       secondary_title = "Seminário sobre especies não tradicionais. Curitiba", 
                       type = "report")
spp <- ref.usos$internal_taxon_id[grepl("Especies nativas para fins produtivos", ref.usos$source)]
car98 <- do.call(rbind.data.frame, 
                    replicate(length(spp), car98, simplify = FALSE))
car98$internal_taxon_id <- spp

#Coradin et al. 2011
cor11 <- data.frame(Reference_type = "Assessment", 
                    author = "Coradin, L, Siminski, A., Reis, A. (eds,)",
                    year = "2011", #volume = 142, number = 6, 
                    pages = "934p.",
                    title = "Espécies nativas da flora brasileira de valor econômico atual ou potencial: plantas para o futuro – Região Sul.", 
                    secondary_title = "Ministério do Meio Ambiente, Brasília", 
                    type = "book")
spp <- ref.usos$internal_taxon_id[grepl("– Região Sul", ref.usos$source)]
cor11 <- do.call(rbind.data.frame, 
                 replicate(length(spp), cor11, simplify = FALSE))
cor11$internal_taxon_id <- spp

#Coradin et al. 2018
cor18 <- data.frame(Reference_type = "Assessment", 
                    author = "Coradin, L., Camillo, J., Pareyn, F.G.C. (eds.)",
                    year = "2018", #volume = 142, number = 6, 
                    pages = "1311p.",
                    title = "Espécies nativas da flora brasileira de valor econômico atual ou potencial: plantas para o futuro: região Nordeste", 
                    secondary_title = "Ministério do Meio Ambiente, Brasília", 
                    type = "book")
spp <- ref.usos$internal_taxon_id[grepl(" região Nordeste", ref.usos$source)]
cor18 <- do.call(rbind.data.frame, 
                 replicate(length(spp), cor18, simplify = FALSE))
cor18$internal_taxon_id <- spp

#Lorenzi 1992
lor92 <- data.frame(Reference_type = "Assessment", 
                    author = "Lorenzi, H.",
                    year = "1992", volume = "1", pages = "384p.",
                    title = "Árvores brasileiras: manual de identificação e cultivo de plantas arbóreas nativas do Brasil", 
                    secondary_title = "Editora Plantarum, Nova Odessa", 
                    type = "book")
spp <- ref.usos$internal_taxon_id[grepl("LORENZI, H", ref.usos$source)]
lor92 <- do.call(rbind.data.frame, 
                 replicate(length(spp), lor92, simplify = FALSE))
lor92$internal_taxon_id <- spp

# Mark et al. 2014
mark14 <- data.frame(Reference_type = "Assessment", 
                    author = "Mark, J., Newton, A.C, Oldfield, S. & Rivers, M.",
                    year = "2014", pages = "56p.",
                    title = "The International Timber Trade: A Working List of Commercial Timber Tree Species", 
                    secondary_title = "Botanic Gardens Conservation International, Richmond, UK", 
                    type = "report")
spp <- ref.usos$internal_taxon_id[grepl("The International Timber Trade", ref.usos$source)]
mark14 <- do.call(rbind.data.frame, 
                 replicate(length(spp), mark14, simplify = FALSE))
mark14$internal_taxon_id <- spp


#Vieira et al. 2016
vie16 <- data.frame(Reference_type = "Assessment", 
                    author = "Vieira, R.F., Camillo, J., Coradin, L. (eds.)",
                    year = "2016", 
                    pages = "1160p.",
                    title = "Espécies nativas da flora brasileira de valor econômico atual ou potencial: Plantas para o Futuro: Região Centro-Oeste", 
                    secondary_title = "Ministério do Meio Ambiente, Brasília", 
                    type = "book")
spp <- ref.usos$internal_taxon_id[grepl(" Região Centro-Oeste", ref.usos$source)]
vie16 <- do.call(rbind.data.frame, 
                 replicate(length(spp), vie16, simplify = FALSE))
vie16$internal_taxon_id <- spp

#Sistema Nacional de Informações Florestais
snif <- data.frame(Reference_type = "Assessment", 
                       author = "Sistema Nacional de Informações Florestais",
                       year = "2019",
                       title = "Espécies Florestais", 
                       secondary_title = "Available at: https://snif.florestal.gov.br/pt-br/especies-florestais", 
                       type = "electronic source")
spp <- ref.usos$internal_taxon_id[grepl("Sistema Nacional de Informações Florestais", ref.usos$source)]
snif <- do.call(rbind.data.frame, 
                    replicate(length(spp), snif, simplify = FALSE))
snif$internal_taxon_id <- spp

# TAXONOMY 
## Brazilian Flora
tax1 <- data.frame(Reference_type = "Taxonomy", 
                   access_date = "6 May 2021",
                   author = "Flora do Brasil 2020 em construção",
                   title = "Jardim Botânico do Rio de Janeiro", 
                   type = "electronic source")
spp <- full.tax$internal_taxon_id[!is.na(full.tax$TaxonomicReference1)]
tax1 <- do.call(rbind.data.frame, 
                replicate(length(spp), tax1, simplify = FALSE))
tax1$internal_taxon_id <- spp
url <- full.tax$TaxonomicReference1[!is.na(full.tax$TaxonomicReference1)]
url <- gsub(".*Available at\\: ", "", url, perl = TRUE)
url <- gsub("\\.$", "", url, perl = TRUE)
tax1$url <- url

## Tropicos
tax2 <- data.frame(Reference_type = "Taxonomy", 
                   access_date = "6 May 2021",
                   author = "Tropicos.org",
                   title = "Missouri Botanical Garden", 
                   type = "electronic source")
spp2 <- full.tax$internal_taxon_id[!is.na(full.tax$TaxonomicReference2)]
tax2 <- do.call(rbind.data.frame, 
                replicate(length(spp2), tax2, simplify = FALSE))
tax2$internal_taxon_id <- spp2
url <- full.tax$TaxonomicReference2[!is.na(full.tax$TaxonomicReference2)]
url <- gsub(".*Available at\\: ", "", url, perl = TRUE)
url <- gsub("\\.$", "", url, perl = TRUE)
tax2$url <- url

## GBIF
tax3 <- data.frame(Reference_type = "Taxonomy", 
                   access_date = "6 May 2021", 
                   author = "GBIF Secretariat",
                   title = "GBIF Backbone Taxonomy", 
                   type = "electronic source") 
spp3 <- full.tax$internal_taxon_id[!is.na(full.tax$TaxonomicReference3)]
tax3 <- do.call(rbind.data.frame, 
                replicate(length(spp3), tax3, simplify = FALSE))
tax3$internal_taxon_id <- spp3
url <- full.tax$TaxonomicReference3[!is.na(full.tax$TaxonomicReference3)]
url <- gsub(".*Available at\\: ", "", url, perl = TRUE)
url <- gsub("\\.$", "", url, perl = TRUE)
tax3$url <- url

## Binding everything together
ref.template <- as.data.frame(matrix(NA, ncol = dim(sample)[2], nrow = 1,
                                     dimnames = list(NULL, names(sample))))
refs.assess <- dplyr::bind_rows(ref1, ref2, ref3, ref4, ref5, ref6, ref7,
                                ref8, ref9, ref10, ref11, ref12, ref13,
                                ref14, ref15)
refs.usos <- dplyr::bind_rows(agro, br.woods, car98, cor11, cor18, lor92,
                              mark14, vie16, snif, wiki)
refs.tax <- rbind.data.frame(tax1, tax2, tax3)

references <- as.data.frame(dplyr::bind_rows(ref.template, 
                                             refs.assess,
                                             refs.usos,
                                             refs.tax))

# Final edits
references <- references[ , match(names(sample), names(references))]
references <- references[-1 ,]
combo <- apply(references[,c("internal_taxon_id", "title")], 1, paste0, collapse = "_")
references1 <- references[!duplicated(combo), ]
references1 <- references1[order(references1$internal_taxon_id), ]

# Removing NAs
for (i in seq_len(dim(references1)[2])) {
  check_these <- references1[,i] %in% c("NA", "", NA)  
  if (any(check_these))
    references1[check_these, i] <- ""
}

## Saving
write.csv(references1, "data/sis_connect/references_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


###############################################################################H
###############################################################################H
# CREDITS ----------------------------------------------------------------------

## Sample from SIS CONNECT v6.1
sample <- read.csv("SIS_sample_6_1/credits.csv")
apply(sample, 2, unique)

## Creating the assessor(s) vector
assess1 <- data.frame(Order = 1L, 
                      affiliation = "Naturalis Biodiversity Center, Leiden, The Netherlands & IUCN SSC Global Tree Specialist Group", 
                      credit_type = "Assessor", email = "raflima@usp.br",
                      firstName = "Renato A.",
                      initials = "R.", lastName = "Ferreira de Lima", nickname = "",
                      user_id = "") 

## Creating the compiler(s) vector
compila1 <- data.frame(Order = 1L, 
                       affiliation = "Naturalis Biodiversity Center, Leiden, The Netherlands & IUCN SSC Global Tree Specialist Group", 
                       credit_type = "Compiler", email = "raflima@usp.br",
                       firstName = "Renato A.",
                       initials = "R.", lastName = "Ferreira de Lima", nickname = "",
                       user_id = "") 

## Creating the reviewer(s) vector
review1 <- data.frame(Order = 1L, affiliation = "",  credit_type = "Reviewer", 
                      email = "", firstName = "", initials = "", lastName = "", 
                      nickname = "", user_id = "") 

## Creating the contributors(s) vector
contrib1 <- data.frame(Order = 1L, 
                       affiliation = "Naturalis Biodiversity Center, Leiden, The Netherlands",
                       credit_type = "Contributor", 
                       email = "hans.tersteege@naturalis.nl", 
                       firstName = "Hans", initials = "H.", lastName = "ter Steege", 
                       nickname = "", user_id = "") 
contrib2 <- data.frame(Order = 2L, 
                       affiliation = "Centro Nacional de Conservação da Flora (CNCFlora) & IUCN SSC Brazil Plant Red List Authority",
                       credit_type = "Contributor", 
                       email = "", firstName = "", initials = "", lastName = "", 
                       nickname = "", user_id = "") 

## Combining all info
n.spp <- dim(tax)[1]
temp <- rbind.data.frame(assess1,  compila1, review1,  contrib1, contrib2,
                         stringsAsFactors = FALSE)
credits <- do.call(rbind.data.frame, 
                   replicate(n.spp, temp, simplify = FALSE))
credits$internal_taxon_id <- rep(tax$internal_taxon_id, each = dim(temp)[1]) 

#Checking the format
sample[sample$internal_taxon_id %in% "22823",]
credits[credits$internal_taxon_id %in% "sp1",]

## Saving
credits <- credits[ , match(names(sample), names(credits))]
write.csv(credits, "data/sis_connect/credits_threat.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")


###############################################################################H
###############################################################################H
# OCCURRENCES ------------------------------------------------------------------

sample <- read.csv("SIS_sample_6_1/ex_maps.csv", sep = "\t")
oc.data <- readRDS("data/threat_occ_data_final.rds")
# oc.data.loc <- readRDS("data/threat_species_by_country.rds")
oc.data1 <- readRDS("data/threat_occ_data.rds")
names(oc.data1)[grepl("numTombo", names(oc.data1))] <- "vouchers"
oc.data <- dplyr::left_join(oc.data, 
                            oc.data1[, c("vouchers", "coletor.name", 
                                         "coletor.number", "dup.ID1", "geo.check1")])
oc.data <- oc.data[!is.na(oc.data$ddlat) | !is.na(oc.data$ddlon), ]
oc.data$SpatialRef <- "WGS84"
oc.data$Compiler <- "RAF Lima"
oc.data$Year <- "2019"
cols <- c("tax", "ddlat", "ddlon", "Compiler", "SpatialRef", "source", "coly", "vouchers",
          #"dup.ID1", 
          "coletor.name", "coletor.number", "detBy", "dety", 
          #"tax.check2", 
          "UC", "tax.check.final", "geo.check1", "typeStatus")
oc.data2 <- oc.data[, .SD, .SDcols = c(cols)]
head(sample, 3)
head(oc.data2, 3)
colnames(oc.data2) <- c("Binomial", "Dec_Lat", "Dec_Long", "Compiler", "SpatialRef",
                        "BasisOfRec", "Event_Year", "collectID", "recordedBy", 
                        "recordNo", "detBy", "detYear", "UC", 
                        "taxCheck", "coordCheck", "typeStatus")
data.table::fwrite(oc.data2, "data/sis_connect/occurrences_threat.csv.zip", compress = "gzip")


###############################################################################H
###############################################################################H
# EOO SHAPEFILES ---------------------------------------------------------------

shapes <- readRDS("data/spp.convex.hull.polys_sf_uncropped.rds")
shapes <- sf::st_transform(shapes, crs = 3857)
shapes <- sf::st_transform(shapes, crs = 4326)
sf::sf_use_s2(FALSE)
shp.area <- sf::st_area(shapes)
small.area <- !as.double(shp.area) <= (2010*2010) # EOO smaller than the minimum AOO
shapes1 <- shapes[small.area, ]
writeLines(shapes1$tax, "data/sis_connect/spp_with_shapes.txt")
toto <- full.tax[ , c("species.correct2", "internal_taxon_id")]
names(shapes1)[1] <- "species.correct2"
shapes2 <- dplyr::left_join(shapes1, toto)
saveRDS(shapes2, "data/sis_connect/polygons_threat.rds", compress = TRUE)

