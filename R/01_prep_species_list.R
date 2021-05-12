###################################################################
###################################################################
#### PREPARING THE LIST OF SPECIES FOR THREAT IUCN ASSESSMENTS ####
###################################################################
###################################################################
rm(list=ls())

#### LOADING THE AF SPECIES LIST ####
spp.af <- readRDS("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/spp.af.rds")

## REMOVING SPECIES THAT SHOULD NOT BE IN THE LIST (EXOTICS, CULTIVATED OR NOT IN THE AF) ###
rm.spp <- c("Annona calophylla","Bauhinia galpinii","Bunchosia glandulifera",
            "Cereus repandus","Eugenia acapulcensis",      
            "Garcinia longifolia","Grabowskia boerhaaviifolia","Humiriastrum cuspidatum",
            "Ixora grandifolia","Maytenus macrocarpa","Miconia pausana",           
            "Mimosa caesalpiniifolia","Mimosa rufescens","Myrceugenia obtusa",
            "Myrsine ligustrina","Ocotea megaphylla", 
            "Prunus oleifolia","Randia obovata","Rauvolfia tetraphylla",
            "Schinus fasciculata","Senegalia fiebrigii",
            "Sesbania macroptera","Sloanea multiflora","Solanum confusum",
            "Solanum lanceifolium","Tarenaya spinosa","Tricerma vitis-idaeum",
            "Vachellia macracantha","Zanthoxylum schreberi")
spp.af1 <- spp.af[!species.correct2 %in% rm.spp,] 

## REMOVING Unresolved, Uncertain application and etc
spp.af2 <- spp.af1[!TreeCo_status %in% c("remove", "remove?", "remove? Unresolved",
                                          "Incorrect", "Missapplied", "Illegitimate",
                                          "Illegitimate name", "Uncertain application",
                                          "Unresolved", "not in the AF"),]
#names with possible problems 
sort(table(spp.af2$TreeCo_status))

## REMOVING ALL WOODY BAMBOS - No abudance data
spp.af3 <- spp.af2[!family.correct1 %in% "Poaceae", ]

## REMOVING NAMES NOT CONFIRMED FOR THE ATLANTIC FOREST EXCEPT FOR SOME MISSING NAMES
## List of missing names (don't know why they were excluded; probably due to missing valid determinations or coordinates)
miss.spp = c("Agonandra brasiliensis","Aiouea bracteata","Alsophila dicomatolepis","Alsophila dichromatolepis",
             "Aspidosperma brasiliense", "Aspidosperma quirandy","Brasiliopuntia schulzii",
             "Cochlospermum vitifolium","Derris leucanthus","Ficus pakkensis","Handroanthus selachidentatus",
             "Handroanthus spongiosus","Ilex cognata","Luetzelburgia purpurea",
             "Mezilaurus synandra", "Ocotea grandifructa", "Ossaea loligomorpha",
             "Palicourea subcuspidata","Palicourea didymocarpos","Palicourea didymocarpa","Persea punctata",
             "Piptadenia killipii","Piranhea securinega", "Piptocarpha regnellii",
             "Pleroma echinatum", "Pleroma vimineum",
             "Quillaja brasiliensis", "Quillaja lancifolia", "Quillaja lanceolata", "Quillaja sellowiana", 
             "Tovomita longifolia", "Trischidium decipiens","Zanthoxylum unifoliolatum")
## Other names, synonyms or possible synonyms not even on the current version of the list
other.spp = c("Brasiliopuntia schulzii", "Derris leucanthus","Muellera campestris",
              "Palicourea subcuspidata","Palicourea didymocarpa","Palicourea didymocarpos",
              "Quillaja brasiliensis", "Quillaja lancifolia","Quillaja lanceolata","Quillaja sellowiana",
              "Ixora grandifolia", "Ixora muelleri", "Misanteca duartei",
              "Calyptranthes grandifolia","Calyptranthes brasiliensis",
              "Ossaea loligomorpha","Miconia loligomorpha", # synonym
              "Plinia brachybotrya","Plinia pseudodichasiantha",
              "Metternichia princeps", "Metternichia principis", 
              "Lonchocarpus guillemineanus", "Lonchocarpus cultratus",
              "Marlierea eugenioides", "Myrcia eugenioides", # synonym
              "Marlierea laevigata", "Myrcia multipunctata", # synonym
              "Piptocarpha regnelii", "Piptocarpha regnellii",
              "Pleroma echinata", "Pleroma echinatum", # synonym
              "Pleroma viminea", "Pleroma vimineum", # synonym
              "Alsophila dicomatolepis","Trichipteris dicomatolepis", "Cyathea dichromatolepis", # synonym
              "Ocotea lucida", "Ocotea brachybotrya",
              "Myrcia oreioeca", "Myrcia aethusa") # synonym
spp.af4 <- spp.af3[af.check2 %in% "TRUE" | 
                     species.correct2 %in% miss.spp | 
                     species.correct2 %in% other.spp, ]
spp.af5 <- spp.af4[ , c("family.correct1", "species.correct2")]

## ADDING ANY MISSING NAME (VALID OR RECENT SYNONYMIZATIONS)
still.miss.spp <- miss.spp[!miss.spp %in% spp.af4$species.correct2]
add.spp <- flora::get.taxa(still.miss.spp)[ , c("family", "original.search")]
names(add.spp) <- c("family.correct1", "species.correct2")
spp.af6 <- rbind.data.frame(spp.af5, add.spp) 
spp.af7 <- spp.af6[order(spp.af6$species.correct2), ]

## SAVING THE SPECIES LIST
dim(spp.af7)
saveRDS(spp.af7, "data/threat_af_spp_list_preliminary.rds")

# ## CHECKING ANY POSSIBLE NEW SYNONYMAZTIONS
# toto <- flora::get.taxa(spp.af6$species.correct2, suggestion.distance = 0.85)
# # toto <- flora::get.taxa(unique(oc.data$species.correct2), suggestion.distance = 0.85)
# toto1 <- toto[!toto$notes %in% c("", "not found"), ] #ok! included in 28/04/2021
# write.csv(toto1,"tmp2.csv", fileEncoding = "UTF-8")
# toto2 <- toto[toto$notes %in% c("not found"), ] # ok! checked in 28/04/2021
# 
# ##Final check for some possible changes in previous synonymizations
# oc.data$concatena <- paste(oc.data$species.correct, oc.data$species.correct1, sep="_")
# temp <- oc.data[!duplicated(oc.data$concatena), ]
# toto3 <- flora::get.taxa(temp$species.correct, suggestion.distance = 0.85)
# toto3$species.correct1 <- temp$species.correct1 
# toto3$species.correct2 <- temp$species.correct2 
# toto3$dup.ID <- temp$dup.ID 
# toto3 <- toto3[order(toto3$dup.ID, na.last = FALSE),]
# toto4 <- toto3[!toto3$notes %in% c("", "not found"), ]
# toto4 <- toto4[!duplicated(toto4$original.search), ]
# write.csv(toto4, "tmp4.csv")
