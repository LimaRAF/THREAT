###################################################################
###################################################################
#### PREPARING THE LIST OF SPECIES FOR THREAT IUCN ASSESSMENTS ####
###################################################################
###################################################################
### INPUT DATA IS BASICALLY THE ONE BY LIMA ET AL. 2020 TO GENERATE THE AF CHECKLIST
### HERE WE USE A MOST UP-TO-DATE TAXONOMY
rm(list=ls())

#### LOADING PACKAGES  ####
require(data.table)

#### LOADING THE HERBARIUM DATA #### -----------------------------------------

## Loading the combined/merged raw occurrence data (over 6 million initial records)
oc.data <- readRDS("./data/data-raw/edited_occurrence_data.rds")

#### REMOVING SPATIAL OUTLIERS  ####
## Getting species information: TBC, geographical range and cultivation 
## Generating the vectors for the groups of species
spp <- read.csv("data/data-raw/DomainsKnownTreesNeotropics.csv", 
               as.is=T, na.string=c(NA,""," "))
tmp <- unique(oc.data[,species.correct2,])
tmp0 <- unique(oc.data[,species.correct1,])
spp <- spp[spp$SpeciesKTSA %in% unique(c(tmp, tmp0)),]
uso <- read.csv("data/data-raw/plant_uses.csv", na= c(""," ",NA))
uso <- uso[as.character(uso$Name_submitted) %in% unique(c(tmp, tmp0)),]

## Preliminary editing of the list
spp$nomes = spp$species.reflora
#replacing names not found in ReFlora by the TNRS/Hans taxonomy
spp$nomes[is.na(spp$nomes)] = spp$Accepted_name[is.na(spp$nomes)]
#replacing names not found in ReFlora and TNRS by the Plant List
spp$nomes[is.na(spp$nomes)] = spp$NewTaxon.TPL[is.na(spp$nomes)]

## Getting the list of species per group
#tbc = spp$nomes[spp$TBC %in% "TBC"]
cult <- 
  as.character(uso$Name_submitted[as.character(uso$group_renato) %in% "cultivated"])

## REMOVING TRUE OUTLIERS ##
#1- remover true outliers para (todas?) as especies: classic Mahalanobis distance >3 & robust (Mahalanobis) distance>16
#table(oc.data$true.out, useNA = "always")
toto = dim(oc.data)[1]
oc.data <- data.table(oc.data[is.na(true.out) | true.out %in% FALSE])
toto1 = dim(oc.data)[1]
toto-toto1; 100*(toto-toto1)/toto # 49325 occs excluded or 0.807%

#### REMOVING PROBABLE OUTLIERS ####
#2- remover probable outliers para cultivated species:  classic Mahalanobis distance >2.5 & robust (Mahalanobis) distance>11.5
tmp <- data.table(oc.data[species.correct2 %in% cult])
#table(tmp$probable.out, useNA = "always")
tmp1 <- data.table(tmp[probable.out %in% TRUE])
oc.data <- oc.data[!numTombo %in% tmp1$numTombo]
toto2 = dim(oc.data)[1]
toto1-toto2; 100*(toto1-toto2)/dim(tmp)[1] # 2974 occs excluded or 1.43%
rm(tmp, tmp1, toto, toto1, toto2)

#### REMOVING DUPLICATES  ####
## Ordering the dataset to more-easily remove duplicates
#first by dup.ID, then decreasing dup.numb, then decreasing dup.prop, then decreasing source name (preferenc to splink, than jabot, then splink_new, then gbif)
oc.data[, source1 := as.double(as.character(factor(source, levels = c("gbif","sndb","jabot","splink","splink_new"), labels = c(5,4,2,1,3)))), ]
setorder(oc.data,  dup.ID, -dup.numb, -dup.prop, source1)               ## re-ordering the data.table, using setorder()
oc.data[,source1:=NULL] #removing the extra column created for ranking

#creating the unique accession number and multiple-acession number for each specimen
oc.data[, dup.ID1 := dup.ID]
oc.data[is.na(dup.ID1), dup.ID1 := numTombo, by = numTombo]
#removing the duplicates
oc.data1 = unique(oc.data, by = "dup.ID1")
#removing the extra column created for ranking
oc.data[,dup.ID1:=NULL] 
oc.data1[,dup.ID1:=NULL]

## How many non-duplicates occurrences?
dim(oc.data1) # 4,628,942 previously; now: 4,859,649

## How many non-duplicates occurrences inside the AF?
oc.data.af = data.table(oc.data1[oc.data1$af.check2 %in% c("TRUE","TRUE_in_transition"),]) # occurrences inside the AF
dim(oc.data.af) # 810,803, previously; now: 808,176 
dim(oc.data.af[af.check2 %in% c("TRUE"),])[1] #then 763,796; now 758,209 truly inside the AF

## How many non-duplicates, taxonomic-validated occurrences inside the AF?
oc.data.af.tax = oc.data.af[oc.data.af$tax.check1 == "TRUE",] # occurrences inside the AF
dim(oc.data.af.tax)
# then: 300,895 (279,424 TRUE and 21,471 TRUE_in_transition)
# now: 299,070 (277,649 TRUE and 21421  TRUE_in_transition)
rm(oc.data1)


#### LIST OF ALL NAMES CITED FOR THE ATLANTIC FOREST ----------------------
## How many names for the records inside the AF (valid, invalid, synonyms, etc.)?
# re-ordering the data.table, using setorder(): af.check == TRUE on top
setorder(oc.data.af.tax,  af.check2)
# valid names
tmp = oc.data.af.tax[!duplicated(oc.data.af.tax$species.correct1),]
tmp$spp.group = "valid"
# extra input names (that will be part of the excluded names from the checklist)
tmp1 = oc.data.af.tax[!oc.data.af.tax$scientificName %in% tmp$species.correct1,]  
tmp1 = tmp1[!duplicated(tmp1$scientificName),]
tmp1 = tmp1[!stringr::str_count(tmp1$scientificName," ") %in% 2 & !stringr::str_count(tmp1$species.correct1," ") %in% 3,]
tmp1$spp.group = "replaced"
# any other name removed by duplicated homogeneization?
tmp2 = oc.data.af.tax[!oc.data.af.tax$species.correct %in% tmp$species.correct1,]
tmp2 = tmp2[!tmp2$species.correct %in% tmp1$scientificName,]
tmp2 = tmp2[!stringr::str_count(tmp2$scientificName," ") %in% 2 & !stringr::str_count(tmp2$species.correct," ") %in% 3,]
tmp2$spp.group = "replaced"
# binding all names
tmp = rbindlist(list(tmp,tmp1,tmp2)) # merging the names from all tables
tmp = tmp[,c("status","scientificName","species.correct","family.correct1","species.correct1","af.check","af.check2","spp.group")]
names(tmp)[2] = "SpeciesKTSA"
rm(tmp1,tmp2)

#number of valid names inside the AF
length(unique(tmp$species.correct1)) # 6027 spp previously now: 5991 spp

## Any valid names not in the entry name list? No!
spp = read.csv("data/data-raw/DomainsKnownTreesNeotropics.csv", 
               as.is=T, na.string=c(NA,""," "))
tmp$species.correct1[!tmp$species.correct1 %in% spp$SpeciesKTSA] #ok!

## Crossing the list of names with the list of valid names found for the AF
spp1= merge(tmp, spp, by.x= "species.correct1", by.y="SpeciesKTSA",all.x=TRUE)

#preliminary edits before replacement
spp1[is.na(spp1$species.reflora)]$spp.group <-
  paste(spp1[is.na(spp1$species.reflora)]$spp.group,"_not_in_reflora",sep="")
spp1$scientific.name.tnrs = paste(spp1$Accepted_name, spp1$Accepted_name_author,sep=" ")
spp1$scientific.name.tpl = paste(spp1$NewTaxon.TPL,spp1$New.Authority.TPL,sep=" ")
#replacing names not found in ReFlora by the TNRS/Hans taxonomy
spp1[is.na(spp1$species.reflora),c("scientific.name.reflora","family.reflora","genus.reflora","species.reflora","name.status.reflora")] =
  spp1[is.na(spp1$species.reflora),c("scientific.name.tnrs","Accepted_name_family","Accepted_name_genus","Accepted_name","Taxonomic_status")]
#replacing names not found in ReFlora and TNRS by the Plant List
spp1[is.na(spp1$species.reflora),c("scientific.name.reflora","family.reflora","genus.reflora","species.reflora")] =
  spp1[is.na(spp1$species.reflora),c("scientific.name.tpl","Family.TPL","New.Genus.TPL","NewTaxon.TPL")]
#removing the extra column created for replacing full species names
spp1[,c("scientific.name.tnrs","scientific.name.tpl"):=NULL]

## Any TBC missing from the list? Yes!
tmp1 = spp[spp$TBC %in% "TBC",]$species.reflora[!spp[spp$TBC %in% "TBC",]$species.reflora %in% spp1$species.reflora]
setorder(oc.data.af,  af.check2) # re-ordering the data.table, using setorder(): af.check == TRUE on top
tmp2 = oc.data.af[oc.data.af$species.correct1 %in% tmp1,]
tmp2 = tmp2[!duplicated(tmp2$species.correct1),]
tmp2$spp.group = "valid_TBC"
tmp2 = tmp2[,c("status","scientificName","species.correct","family.correct1","species.correct1","af.check","af.check2","spp.group")]
spp2= merge(tmp2, spp, by.x= "species.correct1", by.y="SpeciesKTSA",all.x=TRUE)

## Any tree names cited as AF in Reflora missing from our specimens? Yes!
tmp1 = spp[!spp$species.reflora %in% spp1$species.reflora,]
tmp1 = tmp1[grepl("AtlanticForest",tmp1$domain.reflora),]
tmp1 = tmp1[!tmp1$habito %in% c("0","0.5","0.5?") & tmp1$establishment.BR %in% "native",]
tmp2 = oc.data.af[oc.data.af$species.correct1 %in% tmp1$species.reflora,]
tmp2 = tmp2[!duplicated(tmp2$species.correct1),]
tmp2$spp.group = "valid_in_reflora"
tmp2 = tmp2[,c("status","scientificName","species.correct","family.correct1","species.correct1","af.check","af.check2","spp.group")]
spp3= merge(tmp2, spp, by.x= "species.correct1", by.y="SpeciesKTSA",all.x=TRUE)

## Merging names from all "sources": specimens, TBC and missing
spp.all = rbindlist(list(spp1,spp2,spp3),use.names=FALSE) #merging the unicate and duplicated tables
rm(tmp,tmp1,tmp2,spp,spp1,spp2,spp3)
# filtering the final table
spp.all= spp.all[,c("order","checagem_VCS","TreeCo_status","SpeciesKTSA","family.correct1","species.correct1",
                    "af.check", "af.check2", "spp.group",
                    "establishment.AF","establishment.BR","habito","habito.reflora",
                    "scientific.name.reflora","family.reflora","genus.reflora","species.reflora",
                    "taxon.rank.reflora","taxon.status.reflora","name.status.reflora","notes.reflora", 
                    "life.form.reflora","domain.reflora","establishment.reflora",
                    "vegtype.reflora","domain.reflora","endemism.Reflora",
                    "AtlanticForest","Domain","Endemism",
                    "PotentialAtlanticForest","AtlanticForestChecklist_preliminary","treeco",
                    "VCS.jan2019..nome.","VCS.jan2019..habito.","VCS.jan2019..establishment.","Distribution.VCS","TBC"
)]

## Any names left that are known to occur only under cultivation in the AF? Verified by RAFL in 03/01/2019
toto = c("Antrocaryon amazonicum","Attalea cohune","Annona salicifolia","Acoelorraphe wrightii",
         "Brownea rosa-de-monte","Diospyros nigra","Mimozyganthus carinatus","Calliandra calothyrsus",
         "Himatanthus attenuatus","Tabernaemontana heterophylla","Euterpe oleracea",
         "Socratea exorrhiza","Aiphanes aculeata","Dipteryx rosea","Elizabetha paraensis",
         "Inga cinnamomea","Inga pilosula","Vatairea guianensis","Pseudobombax munguba",
         "Plinia pinnata","Erisma calcaratum","Bactris brongniartii","Heterocoma lanuginosa",
         "Clusia insignis","Amburana acreana","Dipteryx punctata","Elizabetha duckei","Elizabetha durissima",
         "Parkia nitida","Allantoma lineata","Guadua superba","Triplaris weigeltiana","Centrolobium paraense",
         "Cecropia obtusa","Cecropia peltata","Mimosa glutinosa","Inga sapindoides","Bertholletia excelsa",
         "Libidibia paraguariensis","Bowdichia nitida","Diplotropis brasiliensis","Swartzia polyphylla",
         "Macrosamanea pubiramea","Macrolobium bifolium","Parkia multijuga",
         "Cereus repandus","Rauvolfia tetraphylla","Capraria biflora")
spp.all[species.reflora %in% toto, establishment.AF := "cultivated"]
toto = c("Bowdichia nitida","Diplotropis brasiliensis","Swartzia polyphylla","Macrosamanea pubiramea",
         "Elizabetha princeps","Isertia coccinea","Psychotria poeppigiana",
         "Garcinia longifolia","Humiriastrum cuspidatum")
spp.all[!establishment.AF %in% "cultivated" & species.reflora %in% toto, establishment.AF := "probably cultivated in JBRJ"]

## Any names left due to spurius locality/synonym assigning or other problems?
toto1 = c("Terminalia amazonia","Dalbergia amazonica","Quiina amazonica","Vachellia albicorticata",
          "Ocotea porphyria","Trichilia trifolia","Platymiscium trinitatis",
          "Grabowskia boerhaaviifolia","Maytenus macrocarpa","Mimosa rufescens",
          "Myrceugenia obtusa","Ocotea megaphylla","Sloanea multiflora","Sloanea parviflora",
          "Solanum confusum","Solanum lanceifolium","Zanthoxylum schreberi")
spp.all[species.reflora %in% toto1, establishment.AF := "probable error in locality description"]

##Removing unecessary objects
rm(oc.data.af, oc.data.af.tax)

#### THE AF CHECKLIST - LIST OF ALL VALID NAMES CONFIRMED --------------------
## Separating synonyms and valid duplicated valid names
spp.synonym <- data.table(spp.all[duplicated(spp.all$species.correct1),])
spp.all <- data.table(spp.all[!duplicated(spp.all$species.correct1),])
spp.synonym <- spp.synonym[!SpeciesKTSA %in% spp.all$species.correct1]
spp.synonym <- spp.synonym[grepl("replaced", spp.synonym$spp.group)]

## Removing names at infraspecific level already cited in the 
tmp <- spp.all$species.correct1[stringr::str_count(spp.all$species.correct1," ") %in% 3]
tmp1 <- sapply(strsplit(tmp," "), function(x) paste(x[1],x[2],sep=" "))
tmp2 <- tmp[tmp1 %in% spp.all$species.correct1]
tmp[!tmp1 %in% spp.all$species.correct1] # species only found at infra-specific level (check synonyms?)
#spp2 = spp2[str_count(spp2$species.correct1," ") %in% 1,]
# spp.all <- spp.all[!spp.all$species.correct1 %in% tmp2,]

## Creating a vector with the names of specimens determined at species and infra-specific levels
spp.all[,species.correct2 := sapply(strsplit(species.correct1," "), function(x) paste(x[1],x[2],sep=" ")),]
spp.synonym[,species.correct2 := sapply(strsplit(species.correct1," "), function(x) paste(x[1],x[2],sep=" ")),]
spp.synonym[,SpeciesKTSA1 := sapply(strsplit(SpeciesKTSA," "), function(x) paste(x[1],x[2],sep=" ")),]
spp.synonym <- spp.synonym[!SpeciesKTSA1 %in% spp.all$species.correct2]

## Removing herbs and lianas (only shrubs, treelets and trees)
spp.all <- spp.all[!habito %in% c("0","check"),]
spp.synonym <- spp.synonym[!habito %in% c("0","check"),]

## Saving the entire Atlantic Forest checklist
setkeyv(spp.all, "species.correct2")
saveRDS(spp.all, "data/all_spp_atlantic_forest.rds")


#### LOADING THE THE AF CHECKLIST #### ---------------------------------------
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
             spp.group %in% c("valid","valid_not_in_reflora","valid_TBC"),])[1] # before:: 5079 spp; new 5044 spp

## REMOVING SPECIES THAT SHOULD NOT BE IN THE LIST (EXOTICS, CULTIVATED OR NOT IN THE AF) ###
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
dim(spp.af7) # now 5087
saveRDS(spp.af7, "data/threat_af_spp_list_preliminary.rds")
rm(list = ls())
