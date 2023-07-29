###################################################################
###################################################################
#### PREPARING THE LIST OF SPECIES FOR THREAT IUCN ASSESSMENTS ####
###################################################################
###################################################################
### DATA ARE BASICALLY THE SAME USED BY LIMA ET AL. 2020 TO GENERATE THE AF CHECKLIST
### HERE WE USE A MOST UP-TO-DATE TAXONOMY
rm(list=ls())

#### LOADING PACKAGES  ####
require(data.table)

#### LOADING THE THE AF CHECKLIST ####
## LIST OF ALL VALID NAMES CONFIRMED FOR THE ATLANTIC FOREST
spp.all <- readRDS("data/all_spp_atlantic_forest.rds")

## Removing exotic, cultivated and naturalized species
spp.af <- spp.all[!spp.all$establishment.BR %in% 
                    c("check","cultivated","cultivated?","exotic","naturalized",
                      "not in Brazil","not in Brazil?"),]

## Species not in Brazil but maybe in the AF?
tmp <- data.table(spp.all[spp.all$establishment.BR %in% "not in Brazil",]) 
tmp <- tmp[is.na(tmp$Distribution.VCS) | tmp$Distribution.VCS %in% c("PAY","PAR","BOL|PAY","PAY|BOL","ARG|BOL","ARG|BOL|PAR")]
tmp <- tmp[is.na(tmp$establishment.AF) & is.na(tmp$checagem_VCS)]
tmp <- tmp[!is.na(tmp$PotentialAtlanticForest) | !is.na(tmp$VCS.jan2019..establishment.)]
spp.af <- data.table(rbind(spp.af, tmp))
spp.af <- spp.af[!establishment.AF %in% c("exotic","naturalized","cultivated",
                                          "cultivated in the AF",
                                          "cultivated?","not in the AF","not in the AF?",
                                          "probable error in locality description",
                                          "probably cultivated in JBRJ"),]

## Removing Unresolved, Uncertain application and etc
spp.af <- spp.af[!TreeCo_status %in% 
                   c("check","remove","remove?","remove? Unresolved","Incorrect",
                     "Missapplied","Illegitimate","Illegitimate name",
                     "Uncertain application","Uncertain Application",
                     "Unresolved","not in the AF"),]

## Checking the results
table(spp.af$establishment.BR, useNA = "always")
table(spp.af$establishment.AF, useNA = "always")
table(spp.af$habito, useNA = "always")
sort(table(spp.af$domain.reflora[!grepl("Atlantic",spp.af$domain.reflora)], useNA = "always"))
table(spp.af$notes.reflora, useNA = "always")
sort(table(spp.af$checagem_VCS))

#names with possible problems 
sort(table(spp.af$TreeCo_status, useNA = "always"))
spp.af <- spp.af[TreeCo_status %like% "^ok",]

#names that are synonyms that are still in the mix
#tmp = get.taxa(spp.af$species.correct1)
#tmp1 = tmp[!tmp$notes %in% "",]

## How many valid names in the Atlantic Forest checklist? Disregarding the 'valid_in_Reflora': they are included only as potential occurrences.
dim(spp.af[af.check2 %in% "TRUE" & 
             spp.group %in% c("valid","valid_not_in_reflora","valid_TBC"),])[1] # before:???; then: 5166 spp; now: 5079 spp; new 5047 spp

## REMOVING SPECIES THAT SHOULD NOT BE IN THE LIST (EXOTICS, CULTIVATED OR NOT IN THE AF) ###
# spp.af <- readRDS("C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo AF checklist/data analysis/spp.af.rds")
setkeyv(spp.af, "species.correct2")
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
# No longer necessary (feito acima)
# spp.af2 <- spp.af1[!TreeCo_status %in% c("remove", "remove?", "remove? Unresolved",
#                                           "Incorrect", "Missapplied", "Illegitimate",
#                                           "Illegitimate name", "Uncertain application",
#                                           "Unresolved", "not in the AF"),]
# #names with possible problems 
# sort(table(spp.af2$TreeCo_status, useNA = "always"))
spp.af2 <- spp.af1

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
