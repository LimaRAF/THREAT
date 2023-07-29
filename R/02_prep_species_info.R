##################################################################H
##################################################################H
#### PREPARING SPECIES INFORMATION FOR THREAT IUCN ASSESSMENTS ####
##################################################################H
##################################################################H
rm(list=ls())
require(taxize)
require(flora)
require(redlistr)

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


#################################################H
#### GENERATION LENGHT AND MATURE INDIVIDUALS ####
#################################################H
rm(list = ls())

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
  7.5, 15, 25, 35, 25, # for shrubs; new values
  # 20, 40, 50, 65, 45, # for small trees; first values
  15, 30, 50, 65, 40, # for small trees; new values
  # 30, 50, 65, 80, 60, # for large trees; first values
  25, 40, 65, 80, 55, # for large trees; new values
  # 25, 45, 50, 60, 50) # for trees unknown; first values
  20, 35, 55, 60, 50) # for trees unknown; new values
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
saveRDS(hab1, "data/threat_habitats_preliminar1.rds")