##############################################################
#### EXTRAPOLATING THE RESULTS TO OTHER TROPICAL HOTSPOTS ####
##############################################################
rm(list = ls())

### LOADING PACKAGES AND FUNCTIONS ###
require(dplyr)
require(writexl)
source("R/99_functions.R")

## In 21/09/2021: 15,160 of the 58,173 tree species were globally 
# threatened (26%) - https://www.bgci.org/resources/bgci-databases/globaltree-portal/global-overview/

## In 05/08/2023: 16,551 of the 58,013 tree species were globally 
# threatened (28.5%) - https://www.bgci.org/resources/bgci-databases/globaltree-portal/global-overview/

# 20-30% das plantas são árvores: florestas  neotropicais  (Gentry  &  Dodson
# 1987b,  Hammel  1990,  Foster  1990,  Foster  &  Hubbell  1990)

### Reading the hotspot information ###
info <- read.csv("data/hotspots_info.csv")
lu_hotspots <- read.csv("data/esa_2018_lc_per_hotspot.csv")
info <- dplyr::left_join(info, lu_hotspots[,c("hotspot.region","area_km2",
                                              "ForestCover",
                                              "Mosaic_ForestCover_GrassLand",
                                              "Mosaic_GrassLand_ForestCover",
                                              "OpenForestCover",
                                              "ShrubLand","ShrubLand_Caatinga.",
                                              "ShrubLand_Cerrado","Shrubland_Restinga")])

## getting the necessary vector to extrapolate threat
#habitat loss
info$habitat_loss <- I(100 - info$ForestCover)
add_open_forest <- info$hotspot.region %in% c("Coastal Forests of Eastern Africa", 
                                               "Eastern Afromontane",
                                               "Guinean Forests of West Africa", 
                                               "Madagascar and the Indian Ocean Islands", 
                                               "Central Africa")
info$habitat_loss[add_open_forest] <- I(100 - (info$ForestCover[add_open_forest] + 
                                          info$OpenForestCover[add_open_forest]/2))
add_mosaic <- info$hotspot.region %in% c("Western Ghats and Sri Lanka")
info$habitat_loss[add_mosaic] <- I(100 - (info$ForestCover[add_mosaic] + 
                                          info$Mosaic_ForestCover_GrassLand[add_mosaic]/2))
# correcting areas which estimates were much lower than previous ones (Schmitt et al. 2009 apud Sloan et al. 2014)
add_correction <- info$hotspot.region %in% c("Caribbean Islands", "New Caledonia")
info$habitat_loss[add_correction] <- I(100 - (info$ForestCover[add_correction] -  
                                            0.25*info$ForestCover[add_correction]))
#endemism
end_ratio <- info$number_of_endemic_plants_perc_Habel_etal_2019
info$number_of_endemic_plants_Habel_etal_2019[grepl("Central Africa", info$hotspot.region)] <- 2600
end_richness <- info$number_of_endemic_plants_Habel_etal_2019
info$plant_species_Myers_etal_2000[grepl("Central Africa", info$hotspot.region)] <- 9500
info$tree_richness <- as.double(sapply(strsplit(info$tree_species," "), function(x) x[1]))
info$tree_richness[grepl("Central Africa", info$hotspot.region)] <- 2000

trop.hotspots <- c("Atlantic Forest", "Caribbean Islands", 
                   "Coastal Forests of Eastern Africa", #"East Melanesian Islands",
                   "Eastern Afromontane", "Guinean Forests of West Africa",
                   "Indo-Burma", "Madagascar and the Indian Ocean Islands",
                   "Mesoamerica", "New Caledonia", "Philippines", 
                   "Sundaland", "Tropical Andes",
                   "Tumbes-Choco-Magdalena", "Wallacea", 
                   "Western Ghats and Sri Lanka",
                   "Amazon", "Central Africa", "New Guinea")
tropical_hotspots <- info$hotspot.region %in% trop.hotspots
info$obs_prop_trees <- info$tree_richness/
                      as.double(sapply(strsplit(info$plant_species_Myers_etal_2000," "), 
                                       function(x) x[1]))
summary(info$obs_prop_trees[tropical_hotspots], 
        na.rm = TRUE) # median 0.290, mean 0.313 [0.25-0.43] for tropical forests; 2nd revision med 0,261 and mean 0.2866
info$tree_end_richness <- end_richness * 0.3
richness <- as.double(sapply(strsplit(info$plant_species_Myers_etal_2000," "), function(x) x[1]))
info$tree_richness <- richness * 0.3

#Generating other important columns 
info$obs_threat_or_ext_endemics_trees <-
  as.double(info$observed_threatened_or_extinct_endemics_Brooks_etal_2002) * 0.3
info$est_threat_or_ext_endemics_trees <-
  info$estimated_threatened_or_extinct_endemic_plants_Brooks_etal_2002 * 0.3

## Filtering to the target regions and info
cols <- c("hotspot.region", "continent", "area_km2", #"domain",
          "plant_species_Myers_etal_2000", "Endemic_plants_perc_Myers_etal_2000",
          "number_of_endemic_plants_Habel_etal_2019",
          "number_of_endemic_plants_perc_Habel_etal_2019",
          "ForestCover", "OpenForestCover",
          "habitat_loss", "obs_prop_trees", "tree_richness", "tree_end_richness",
          "obs_threat_or_ext_endemics_trees", "est_threat_or_ext_endemics_trees") 
info_TF <- info[tropical_hotspots, cols]

## extrapolating threat
#reading the spline fits to the empirical AF threat-habitat loss relationship
splines <- readRDS("data/spline_threat.per.loss.clumped.rds")
splines.VU <- readRDS("data/spline_threat.per.loss.clumped_VU.rds")
splines.EN <- readRDS("data/spline_threat.per.loss.clumped_EN.rds")
splines.CR <- readRDS("data/spline_threat.per.loss.clumped_CR.rds")
threatened_per_hotspot <- 
  cbind.data.frame(info_TF, 
                   get_threat(info_TF$habitat_loss, info_TF$tree_end_richness,
                              all = splines, VU = splines.VU, 
                              EN = splines.EN, CR = splines.CR))
#head(threatened_per_hotspot, 3)

#### Sumary stats ####
#BGCI 2021: total = 58497
#BGCI 2021: without Neartic, Paleartic and Oceania = 49469
58497 - (1432 + 5994 + 1602)
#BGCI 2021: without Neartic, Paleartic and Oceania + species and Australia and New Zealand endemcis = 48140  
58497 - (1432 + 5994 + 2730 + 201) # or 48,140
23631 + 9237 + 13739 + (7442 - (2730 + 201)) # 51,118


sum(threatened_per_hotspot$area_km2) / 1000000 # area cover in million km2
round(sum(threatened_per_hotspot$area_km2) / 1000, 2) # area cover in thousand km2
round(apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", 
                                       "high_all_spp")], 2, sum, na.rm = TRUE), 0) # total number of threatened species

100 * round(apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", 
                                      "high_all_spp")], 2, sum, na.rm = TRUE), 0)/
  sum(threatened_per_hotspot$tree_end_richness)# threat proportion of endemic tropical trees threatened: 57-64-69%

apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", 
                                "high_all_spp")], 2, sum, na.rm = TRUE)/584.97 # 34.9-38.9-42.3% total number of trees in the world; 2nd revision 35.1-39.0-42.6%
apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", 
                                "high_all_spp")], 2, sum, na.rm = TRUE)/481.40 # 42.1-47.0-51.1% total number of tropical trees in the world; 2nd revision 42.6-47.4.0-51.8%
apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", 
                                "high_all_spp")], 2, sum, na.rm = TRUE)/
sum(threatened_per_hotspot$tree_end_richness) # 56.5-68.5% of the total number of trees in these forests; 2nd revision 57.1-63.5-69.4%
sum(threatened_per_hotspot$tree_end_richness)/584.97 #61.7% of total number of trees in the world; 2nd revision 61.4%
sum(threatened_per_hotspot$tree_end_richness)/481.40 #76.6% of total number of tropical trees in the world; 2nd revision 74.6%

sum(threatened_per_hotspot$med_all_spp)/584.97 #38.9% of endemic threat in respect to the total number of trees in the world; 2nd revision 39.0%
sum(threatened_per_hotspot$med_all_spp)/481.40 #47.0% of endemic threat in respect to the number of tropical trees in the world; 2nd revision 47.4%

rem_habitat_area <- ((100 - threatened_per_hotspot$habitat_loss)/100) * threatened_per_hotspot$area_km2
sum(rem_habitat_area)/sum(threatened_per_hotspot$area_km2) # overall 57.9% of lost in these regions; 2nd revision the same  

100*round(sum(threatened_per_hotspot$med_all_spp), 0)/
  sum(threatened_per_hotspot$tree_end_richness) # 62.9% of the endemic tree species are threatened; 2nd revision 63.5%

# Simples regra de três: quantas árvores no total devem estar ameaçadas no mundo)
other.spp <- 58497 - sum(threatened_per_hotspot$tree_end_richness)
sum(threatened_per_hotspot$med_all_spp)/584.97 + 
  (0.25*other.spp)/584.97 # considering a 25% theat level for non-endemic tropical trees; 2nd revision 48.6%

saveRDS(threatened_per_hotspot, "data/hotspots_results.rds")

#### MAKING THE TABLES FOR THE TEXT ####
# Table 2 (was removed from the main text and merged with Table S6)
# cols1 <- c("hotspot.region","habitat_loss",
#            "number_of_endemic_plants_Habel_etal_2019","tree_end_richness",
#            "est_threat_or_ext_endemics_trees",
#            "low_all_spp", "med_all_spp", "high_all_spp")
# result3 <- threatened_per_hotspot[, cols1]
# result3$est_range <- apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", "high_all_spp")], 
#                            1, function(x) paste(ceiling(x), collapse =  "-"))
# result4 <- result3[,!names(result3) %in% c("low_all_spp", "med_all_spp", "high_all_spp")]
# # result4$area_km2 <- ceiling(result4$area_km2)
# # result4$ForestCover <- round(result4$ForestCover, 1)
# result4$tree_end_richness <- ceiling(result4$tree_end_richness)
# result4$est_threat_or_ext_endemics_trees <- 
#   round(result4$est_threat_or_ext_endemics_trees, 0)
# # result4$HabitatLoss <- 100 - result4$ForestCover
# 
# names(result4) <- c("Biodiversity Hotspot", #"Total extent (km2)", 
#                     "Habitat Loss (%)a", "Endemic plantsb",
#                     "Endemic treesc", "Estimated threate", "Estimated threatd")
# result4 <- result4[, c("Biodiversity Hotspot", #"Total extent (km2)", 
#                        "Habitat Loss (%)a", "Endemic plantsb",
#                        "Endemic treesc", "Estimated threatd", "Estimated threate")] 
# write.csv(result4, "data//Table2.csv")

## Table S6 and S7
cols1 <- c("hotspot.region","area_km2",
           # "plant_species_Myers_etal_2000",
           "number_of_endemic_plants_Habel_etal_2019",
           "tree_end_richness",
           "habitat_loss",
           # "ForestCover", "OpenForestCover",
           "low_all_spp", "med_all_spp", "high_all_spp",
           # "number_of_endemic_plants_perc_Habel_etal_2019",
           "low_all", "med_all", "high_all", 
           "low_VU", "med_VU", "high_VU", 
           "low_EN", "med_EN", "high_EN", 
           "low_CR", "med_CR", "high_CR",
           "low_rli", "median_rli", "high_rli",
           "est_threat_or_ext_endemics_trees"
)
result3.1 <- threatened_per_hotspot[, cols1]
result3.1$all <- apply(threatened_per_hotspot[,c("low_all", "med_all", "high_all")], 
                           1, function(x) paste0(x[2]," [",paste(x[c(1,3)], collapse =  "-"),"]"))
result3.1$all_spp <- apply(threatened_per_hotspot[,c("low_all_spp", "med_all_spp", "high_all_spp")], 
                       1, function(x) paste0(round(x[2],0)," [",paste(round(x[c(1,3)],0), collapse =  "-"),"]"))
result3.1$VU <- apply(threatened_per_hotspot[,c("low_VU", "med_VU", "high_VU")],
                       1, function(x) paste0(x[2]," [",paste(x[c(1,3)], collapse =  "-"),"]"))
result3.1$EN <- apply(threatened_per_hotspot[,c("low_EN", "med_EN", "high_EN")],
                       1, function(x) paste0(x[2]," [",paste(x[c(1,3)], collapse =  "-"),"]"))
result3.1$CR <- apply(threatened_per_hotspot[,c("low_CR", "med_CR", "high_CR")],
                       1, function(x) paste0(x[2]," [",paste(x[c(1,3)], collapse =  "-"),"]"))
result3.1$rli <- apply(threatened_per_hotspot[,c("low_rli", "median_rli", "high_rli")],
                       1, function(x) paste0(x[2]," [",paste(round(x[c(1,3)],3), collapse =  "-"),"]"))
result4.1 <- result3.1[,!names(result3.1) %in% c("low_all", "med_all", "high_all", 
                                           "low_VU", "med_VU", "high_VU", 
                                           "low_EN", "med_EN", "high_EN", 
                                           "low_CR", "med_CR", "high_CR",
                                           "low_rli", "median_rli", "high_rli")]
result4.1$area_km2 <- round(ceiling(result4.1$area_km2)/1000, 2)
# result4.1$ForestCover <- round(result4.1$ForestCover, 1)
# result4.1$OpenForestCover <- round(result4.1$OpenForestCover, 1)
head(result4.1)
# Table S6
result5.1 <- result4.1[, c("hotspot.region","area_km2",
                           "number_of_endemic_plants_Habel_etal_2019",
                           "tree_end_richness",
                           "habitat_loss",
                           "all", "all_spp", 
                           "est_threat_or_ext_endemics_trees")]
names(result5.1) <- c("topical.forest", "area_K_km2", "endemic_plants", 
                      "endemic_trees", "habitat_loss",
                      "overall_threat", "spp_threat", "previous_threat")
writexl::write_xlsx(result5.1, "tables/TableS6.xlsx")

# Table S7
result6.1 <- result3.1[, c("hotspot.region",
                           "habitat_loss", 
                           "VU", "EN", "CR", "rli")]
names(result6.1) <- c("topical.forest", "habitat_loss", 
                      "VU", "EN", "CR", "Red List Index")
writexl::write_xlsx(result6.1, "tables/TableS7.xlsx")
