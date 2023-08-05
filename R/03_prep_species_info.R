##################################################################H
##################################################################H
#### PREPARING SPECIES INFORMATION FOR THREAT IUCN ASSESSMENTS ####
##################################################################H
##################################################################H
rm(list=ls())
require(taxize)
require(flora)
require(redlistr)
source("R/99_functions.R")

##############################
#### THREAT FULL TAXONOMY ####
##############################
# Based on different taxonomic backbones
# Note: running the codes here need API keys 
# So, the resulting objects were already saved in the data folder and are loaded down here

## Reading the THREAT species list ##
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]

# #### NCBI ####
# ## Seting the ncbi key
# options(ENTREZ_KEY = "")
# 
# #higher taxonomy
# fams <- sort(unique(tax$family.correct1))
# # tmp  <- sapply(fams, taxize::classification, "ncbi") ## was having Bad Request Errors (HTTP 400)
# tmp <- vector("list", length = length(fams))
# for (i in seq_along(tmp)) {
#   res <- try(taxize::classification(fams[i], db = "ncbi"), TRUE)
#   
#   if(class(res) == "try-error") {
#     Sys.sleep(5)
#     res <- try(taxize::classification(fams[i], db = "ncbi"), TRUE)
#     
#     if(class(res) == "try-error") {
#       Sys.sleep(15)
#       res <- try(taxize::classification(fams[i], db = "ncbi"), TRUE)
#     } 
#   }
#   
#   if(class(res) == "try-error")
#     warning(fams[i])
#   
#   tmp[i] <- res
# }  
# tmp1 <- lapply(tmp, function(x) x$name[x$rank %in% c("kingdom","phylum","class","order")]) 
# tmp1 <- do.call(rbind.data.frame, tmp1)
# tmp2 <- cbind.data.frame(tmp1, fams, stringsAsFactors = FALSE)
# names(tmp2) <- c("kingdom","phylum","classname","ordername","family")
# ncbi.taxonomy <- merge(tax, tmp2, by.x = "family.correct1", by.y = "family", all.x = TRUE, order = FALSE)
# ncbi.taxonomy$genus <- sapply(strsplit(ncbi.taxonomy$species.correct2, " "), function(x) x[1])
# ncbi.taxonomy <- ncbi.taxonomy[,c("kingdom","phylum","classname","ordername","family.correct1","genus","species.correct2")]
# ncbi.taxonomy$kingdom <- "Plantae" #IUCN standard notation
# ncbi.taxonomy$phylum <- "Tracheophyta" #IUCN standard notation
# ncbi.taxonomy <- ncbi.taxonomy[order(ncbi.taxonomy$species.correct2),]
# saveRDS(ncbi.taxonomy, "data/threat_ncbi_taxonomy.rds")

# #### TROPICOS ####
# ## Getting Tropicos full taxonomy ##
# options(TROPICOS_KEY = "")
# 
# #higher taxonomy
# tmp  <- sapply(fams, taxize::classification, "tropicos", rows = 1)
# tmp1 <- lapply(tmp, function(x)
#   as.character(x$name[x$rank %in% c("class", "order")]))
# tmp1 <- do.call(rbind.data.frame, tmp1)
# tmp2 <- cbind.data.frame(tmp1, fams, stringsAsFactors = FALSE)
# names(tmp2) <- c("classname","ordername","family")
# tps.taxonomy <- merge(tax, tmp2, by.x = "family.correct1", by.y = "family", all.x = TRUE, order = FALSE)
# tps.taxonomy$genus <- sapply(strsplit(tps.taxonomy$species.correct2, " "), function(x) x[1])
# tps.taxonomy$kingdom <- "Plantae" #IUCN standard notation
# tps.taxonomy$phylum <- "Tracheophyta" #IUCN standard notation
# tps.taxonomy <- tps.taxonomy[,c("kingdom","phylum","classname","ordername","family.correct1","genus","species.correct2")]
# 
# #lower taxonomy
# cl <- snow::makeSOCKcluster(6)
# doSNOW::registerDoSNOW(cl)
# `%d%` <- foreach::`%dopar%`
# x <- NULL
# output <- foreach::foreach(
#   x = 1:length(tps.taxonomy$species.correct2),
#   #.combine = 'c',
#   .options.snow = NULL
# ) %d% {
#   try(taxize::get_tpsid_(tps.taxonomy$species.correct2[x],
#                          key = "E8E538BC-0DAD-47EF-84CA-F6A663D9170A"),
#       TRUE)
# }
# snow::stopCluster(cl)
# output1 <- sapply(output, function(x) x[[1]])
# names(output1) <- sapply(output, function(x) names(x))
# output2 <- lapply(1:length(output1), function(x)
#   if (is.null(output1[[x]])) {
#     c(nameid= NA, scientificname = names(output1[x]), rankabbreviation = "sp.", 
#       nomenclaturestatusname = "Not found", author=NA, displayreference= NA, displaydate = NA,
#       search.string = names(output1[x]))
#   } else {
#     out <- output1[[x]][,c("nameid",
#                            "scientificname",
#                            "rankabbreviation",
#                            "nomenclaturestatusname",
#                            "author",
#                            "displayreference",
#                            "displaydate")]
#     if (dim(out)[1] > 1) {
#       out <- out[out$rankabbreviation %in% "sp.", ]
#       if (dim(out)[1] > 1 & any(out$scientificname %in% names(output1[x])))
#         out <- out[out$scientificname %in% names(output1[x]),]
#       if (dim(out)[1] > 1 & !all(out$nomenclaturestatusname %in% "Invalid"))
#         out <- out[!out$nomenclaturestatusname %in% "Invalid",]
#       if (dim(out)[1] > 1 & !all(out$nomenclaturestatusname %in% "Illegitimate"))
#         out <- out[!out$nomenclaturestatusname %in% "Illegitimate",]
#       if (dim(out)[1] > 1 & "Legitimate" %in% out$nomenclaturestatusname)
#         out <- out[out$nomenclaturestatusname %in% "Legitimate",]
#       if (dim(out)[1] > 1 & !all(out$displayreference %in% ""))
#         out <- out[!out$displayreference %in% "",]
#       if (dim(out)[1] > 1)
#         out <- out[which.min(out$displaydate),]
#     }
#     cbind.data.frame(out, search.string = names(output1[x]), 
#                      stringsAsFactors = FALSE)
#   }
# )
# output3 <- do.call(rbind.data.frame, output2)
# tps.taxonomy1 <- merge(tps.taxonomy, output3, by.x = "species.correct2", 
#                        by.y = "search.string", all.x = TRUE, sort= FALSE)
# saveRDS(tps.taxonomy1, "data/threat_tropicos_taxonomy.rds")

# #### IUCN ####
# ## Trying to get species ID from IUCN
# cl <- snow::makeSOCKcluster(6)
# doSNOW::registerDoSNOW(cl)
# `%d%` <- foreach::`%dopar%`
# x <- NULL
# iucn <- foreach::foreach(
#   x = 1:length(tps.taxonomy$species.correct2),
#   #.combine = 'c',
#   .options.snow = NULL
# ) %d% {
#   try(taxize::get_iucn(tax$species.correct2[x], 
#                        key = ""), TRUE)
# }
# snow::stopCluster(cl)
# 
# iucn1 <- lapply(1:length(iucn), function(x)
#   if (is.null(iucn[[x]])) {
#     c(iucnid= NA, match = "error", name = tax$species.correct2[x], uri = NA)
#   } else {
#     atributos <- unlist(attributes(iucn[[x]]))
#     out <- c(iucnid= iucn[[x]][1], atributos[c(2:4)]) 
#   }
# )
# iucn2 <- do.call(rbind.data.frame, iucn1)
# names(iucn2) <- c("Redlist_id","match","name","uri")
# for(i in 1:3) iucn2[,i] <- as.character(iucn2[,i])
# saveRDS(iucn2, "data/threat_iucn_taxonomy.rds")

#### BRAZILIAN FLORA INFORMATION ####
## Getting FB 2020 taxonomy
fbo.info <- readRDS("data/threat_fbo_tax_info.rds")
toto <- dplyr::left_join(tax, fbo.info, by = "species.correct2")
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
miss.tax <- taxize::get_gbifid_(splist)
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


#############################H
#### HABITATS AND ECOLOGY ####
#############################H

## Getting Flora do BRasil information ##
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]

fbo.info <- readRDS("data/threat_fbo_tax_info.rds")
toto <- dplyr::left_join(tax, fbo.info, by = "species.correct2")
hab <- toto[, c(names(tax), "life.form","habitat","vegetation.type","establishment")]

#### OBTAINING AND GENERALISING SPECIES INFORMATION ####
## List of species
hab$genus <- sapply(strsplit(hab$species.correct2," "), function(x) x[1])
hab$taxon.rank <- "species"
hab <- hab[,c("internal_taxon_id","species.correct2","family.correct1","genus",
              "species.correct2","taxon.rank","life.form","habitat","vegetation.type")]
names(hab) <- c("internal_taxon_id", "Name_submitted","family","genus",
                "species.correct","taxon.rank","life.form.reflora","habitat.reflora",
                "vegetation.type.reflora")
hab$ordem.spp = 1:dim(hab)[1]

## South American tree species list
trees.sa = read.csv("data/data-raw/DomainsKnownTreesNeotropics.csv",as.is=TRUE,na.string=c(""," ",NA))
## Reading trait data per taxonomic level ##
traits.treeco <- readRDS("data/data-raw/traits_treeco.rds")
trait.spp <- traits.treeco$species
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$status %in% c("replace"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace") {
    trait.spp$TAX[trait.spp$TAX %in% sp.i] <- rpl.i
  }
}
trait.spp <- trait.spp[!duplicated(trait.spp$TAX),]
trait.gen <- traits.treeco$genus
trait.fam <- traits.treeco$family
rm(traits.treeco)
traits1 <- trait_data_prep(hab, trees.sa, trait.spp, trait.gen, trait.fam)$data_frame
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
hab$establishment[is.na(hab$establishment)] <- 
  traits1$establishment.BR[is.na(hab$establishment)]
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
hab$life.form[hab$life.form %in% "succulent_tree" & !is.na(aux$Growth.habits) 
              & grepl("Succulent", aux$Growth.habits)] <- 
  aux$Growth.habits[hab$life.form %in% "succulent_tree" & !is.na(aux$Growth.habits) 
                    & grepl("Succulent", aux$Growth.habits)]
hab$life.form[hab$life.form %in% "succulent_tree"] <- "Succulent tree"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & 
                grepl("^Tree$", aux$Growth.habits)] <- "Woody tree (Tree)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & 
                grepl("^Shrub or treelet$", aux$Growth.habits)] <- "Woody tree (Shrub or treelet)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & 
                grepl("^Tree or shrub$", aux$Growth.habits)] <- "Woody tree (Tree or shrub)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & 
                grepl("^Hemiepiphytic", aux$Growth.habits)] <- "Woody tree (Hemiepiphytic)"
hab$life.form[hab$life.form %in% "woody_tree" & !is.na(aux$Growth.habits) & 
                grepl("^Lianescent", aux$Growth.habits)] <- "Woody tree (Lianescent)"
hab$life.form[hab$life.form %in% "woody_vines_and_subshrubs"] <- "Woody tree (Lianescent)"
hab$life.form[hab$life.form %in% "palm"] <- "Palm or palmoid"
hab$life.form[hab$life.form %in% "palmoids"] <- "Palm or palmoid"
hab$life.form[hab$life.form %in% "tree_fern"] <- "Tree fern"
hab$life.form[hab$life.form %in% "woody_bamboo"] <- "Woody bamboo"
hab$life.form[hab$life.form %in% "woody_tree" & 
                hab$habito %in% "shrub"] <- "Woody tree (Shrub or treelet)"
hab$life.form[hab$life.form %in% "woody_tree" & 
                hab$habito %in% "tree"] <- "Woody tree (Tree)"
hab$life.form[hab$life.form %in% "woody_tree" & 
                hab$life.form.reflora %in% "Árvore"] <- "Woody tree (Tree)"
hab$life.form[hab$life.form %in% "woody_tree" & 
                hab$life.form.reflora %in% c("Arbusto|Subarbusto", "Subarbusto")] <- "Woody tree (Shrub or treelet)"
hab[hab$species.correct %in% "Mimosa setosa",c("habito","life.form")] <- 
  c("shrub", "Woody tree (Shrub or treelet)")
table(hab$life.form, useNA = "always")

## Flagging species according to their mean maximum size
# There are different sources of potential heigt and 2 ways of summaring these info (mean and max quantile)
# We were using at the beggining max quantile plus Ary's info for missing species
# Now, to avoid the effect of outliers (mostly due to misidentifications), we are using the mean between both (didn't changed much)
# Anyways, initial classifications are fine tunned using the available info of growth from from the BFO
# hab$MaxHeight <- round(traits1$MaxHeight_m, 1)
# hab$MaxHeight <- round(traits1$MeanMaxHeight_m, 1)
hab$MaxHeight <- round(apply(traits1[,c("MeanMaxHeight_m","MaxHeight_m")], 1, mean, na.rm = T), 1)
hab$MaxHeight[is.na(hab$MaxHeight)] <- 
  as.double(aux$Potential.height[is.na(hab$MaxHeight)])

#some final generalizations
small.shrubs = c("Athenaea","Capsicum","Piper","Eumachia","Chiococca","Stachytarpheta","Schweiggeria")
hab$MaxHeight[grepl(paste(small.shrubs,collapse="|"),hab$Name_submitted) & is.na(hab$MaxHeight)] = 3
taller.shrubs = c("Erythroxylum","Chamaecrista","Faramea","Palicourea","Psychotria","Conchocarpus","Strychnos","Cybianthus","Chomelia")
hab$MaxHeight[grepl(paste(taller.shrubs,collapse="|"),hab$Name_submitted) & is.na(hab$MaxHeight)] = 4
table(is.na(hab$MaxHeight), useNA = "always")
# hab[hab$habito %in% "shrub" & !is.na(hab$MaxHeight) & hab$MaxHeight >15,]
# hab[hab$habito %in% "tree" & !is.na(hab$MaxHeight) & hab$MaxHeight <5,]

## Classifying species into the growth form classifications
# based on available info on growth form and potential adult height
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
hab$GF <- stringr::str_replace_all(hab$GF, c("1" = "large_shrub", 
                                             "2" = "small_tree", 
                                             "3" = "large_tree", 
                                             "4" = "tree_unknown"))
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
saveRDS(hab, "data/threat_habitats_preliminar.rds")


#################################################H
#### GENERATION LENGHT AND MATURE INDIVIDUALS ####
#################################################H
## Reading the preliminary habitat data
hab <- readRDS("data/threat_habitats_preliminar.rds")

#### SETTING THE PROXIES OF GENERATION LENGTH FOR EACH SPECIES ###
## Defining the generation lengths per combination of GF and ecol.group
g1 <- sort(paste(c(1,3,2,4), sort(unique(hab$GF)),sep="."))
g2 <- sort(paste(c(4,2,3,1,5), sort(unique(hab$ecol.group)),sep = "."))
combo <- expand.grid(EG = g2, GF = g1)
combo[,1] <- gsub('[0-9:]\\.',"", combo[,1])
combo[,2] <- gsub('[0-9:]\\.',"", combo[,2])

#Generation lengths
combo$GL <- c(
  # 10, 20, 25, 35, 25, # for shrubs; first values
  7, 15, 25, 35, 20, # for shrubs; new values
  # 20, 40, 50, 65, 45, # for small trees; first values
  15, 30, 50, 65, 40, # for small trees; new values
  # 30, 50, 65, 80, 60, # for large trees; first values
  20, 40, 65, 80, 55, # for large trees; new values
  # 25, 45, 50, 60, 50) # for trees unknown; first values
  18, 35, 55, 60, 45) # for trees unknown; new values
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

## LOADING THE OF PROP. OF MATURE INDIVIDUALS EXPANDED FOR ALL SPECIES
# Values of the prop. of mature individuals in the population was estimated from tree-by-tree measurements in the TreeCo  database
p.est <- read.csv("data/prop_mature.csv", as.is = TRUE)

#Adding Proportion of mature individuals in the population
hab1 <- dplyr::left_join(hab1, p.est[,c("species.correct","N","dbh.crit","p")], 
                         by= "species.correct") 
## Saving preliminar habitat data
saveRDS(hab1, "data/threat_habitats.rds")


########################H
#### ENDEMISM LEVELS ####
########################H

## Reading species endemism levels in respect to the AF (from the SOM of Lima et al. 2020 Biological Conservation 252: 108825)
af.list <- read.csv("data/AppendixC_af_checklist.csv", as.is=TRUE)
end <- read.csv("data/AppendixF_endemism_levels.csv", as.is=TRUE)

## Reading probable occurrences of tree species in the AF (from the same source as above)
end.prob <- read.csv("data/AppendixD_probable_occurrences.csv", as.is=TRUE)[,1:5]

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
syn.br <- read.csv("data/new_synonyms_floraBR.csv", 
                   na.strings = c(""," ",NA), as.is = TRUE)
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

## Getting the endemism from other AF countries
# end.final1 <- dplyr::left_join(end.final, 
#                                all.crit[,c("species", "endemism.ARG")])

## Saving
saveRDS(end.final, "data/threat_endemism.rds")




###################H
#### CITES LIST ####
###################H


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
cites1 <- readRDS("data/cites_data.rds")

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
tax <- readRDS("data/herbarium_spp_data.rds")[, c(1, 2, 3)]

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
decis.EU1 <- aggregate(. ~ taxon_concept_id, data = na.omit(decis.EU), 
                       function(x) paste(x, collapse = "|")) 
all_cites.EU_spp1 <- merge(all_cites.EU_spp,  decis.EU1,
                           by.x = "id", by.y = "taxon_concept_id",
                           all.x = TRUE, sort = FALSE)
## Saving 
all_cites.EU_spp1 <- 
  all_cites.EU_spp1[order(all_cites.EU_spp1$species.correct2),]
saveRDS(all_cites.EU_spp1, "data/threat_cites_EU.rds")
rm(list=ls())
