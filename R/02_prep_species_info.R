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


################################################################################H
################################################################################H
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
hab1 <- dplyr::left_join(hab, combo[,c("str.match","GL","DBH","p.est","p.ci")], 
                         by= "str.match") 

#Adding Proportion of mature individuals in the population (from code '')
p.est <- read.csv("data/prop_mature.csv", as.is = TRUE)
hab1 <- dplyr::left_join(hab1, p.est[,c("species.correct","N","dbh.crit","p")], 
                         by= "species.correct") 
## Saving preliminar habitat data
write.csv(hab1, "data/threat_habitats_preliminar1.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
