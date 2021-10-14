##############################################################
#### EXTRAPOLATING THE RESULTS TO OTHER TROPICAL HOTSPOTS ####
##############################################################
rm(list = ls())

## In 21/09/2021: 15,160 of the 58,173 tree species were globally 
# threatened (26%) - https://www.bgci.org/resources/bgci-databases/globaltree-portal/global-overview/

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

## Theorectical curves 
y <- c(0,0.001, 0.999, 1)
x <- c(0,0.001, 0.999, 1)
lm(y ~ x)
logis <- try(stats::nls(y ~ exp(a + b*x)/(1 + exp(a + b*x)),
                        start = list(a=1,b=0.1)), TRUE)
# mm <- try(stats::nls(y ~ x/(b + x), 
#                      start = list(b=0.08)), TRUE)
# pow <- try(stats::nls(y[-1] ~ a*x[-1]^b, 
#                      start = list(b=0.01)), TRUE)

par(mfrow = c(1,1))
plot(c(0,1) ~  I(c(0,1)), xlim=c(0,1), ylim = c(0,1), col = "white",
     xlab = "Habitat loss", ylab = "Threatened species")
abline(v = 0.3, lty = 3)
curve(0 + 1*x, from = 0, to = 1, add = TRUE, lwd = 2)
# curve(exp(coef(logis)[1] + coef(logis)[2]*x)/
#         (1 + exp(coef(logis)[1] + coef(logis)[2]*x)), from = 0, to = 1,
#       lwd=2, col = 2, add= TRUE, lty = 2)
curve(exp(-5.5 + 11*x)/ (1 + exp(-5.5 + 11*x)), from = 0, to = 1,
      lwd=2, col = 2, add= TRUE)
# a = -1
# curve(a*x^2 +(1 − a/2)*x, from = 0, to = 1,
#       lwd=2, col = 6, add= TRUE)
curve(exp(-(log(1/x))^0.6), from = 0, to = 1,
      lwd=2, col = 6, add= TRUE)
curve(x^0.25, add= TRUE, lwd = 1, col=4)
curve(x^0.5, add= TRUE, lwd = 2, col=4)
# curve(x/(coef(mm)[1] + x), from = 0, to = 1,
#       lwd=2, col = 3, add= TRUE)
curve(x^2, add= TRUE, lwd = 2, col=4, lty = 2)
selected <- info$hotspot.region == "Atlantic Forest"
threat.all <- I(info$observed_percentage_threat_all_criteria[selected]/100)
threat.end <- I(info$observed_percentage_threat_endemic_all_criteria[selected]/100)
loss <- I((100 - info$ForestCover[selected])/100)
points(threat.all ~ loss, pch = 21)
points(threat.end ~ loss, pch = 19)
points(0.114, 0.25, pch = 17) # Amazon from ter Steege et al. 2015
points(0.114, 0.09, pch = 17) # Amazon from ter Steege et al. 2015, criteria A

## function to obtain threat from habitat loss
get_threat <- function(habitat_loss = NULL, spp_richness = NULL, 
                       endemism_ratio = NULL, pow1 = 0.5, pow2 = 2, 
                       logis.a = -5.5, logis.b = 11) {
  # getting the proportion of loss
  if (any(habitat_loss[!is.na(habitat_loss)] > 1))
    habitat_loss <- habitat_loss/100
    
  # getting the threat proportion
  lin_threat <- 0 + habitat_loss*1
  pow1_threat <- habitat_loss^pow1
  pow2_threat <- habitat_loss^pow2
  logis_threat <- exp(logis.a + logis.b*habitat_loss)/
                    (1 + exp(logis.a + logis.b*habitat_loss))
  
  threat <- cbind(lin_threat, pow1_threat, pow2_threat, logis_threat)
  
  # getting the number of threatened species
  if (!is.null(spp_richness)) {
    threat_spp <- round(spp_richness * threat, 1)
    colnames(threat_spp) <- paste0(colnames(threat_spp),"_spp") 
      
      # getting the number of threatened endemic species
      if (!is.null(endemism_ratio)) {
        if (any(endemism_ratio[!is.na(endemism_ratio)] > 1))
          endemism_ratio <- endemism_ratio/100
        
        threat_end_spp <- threat_spp * endemism_ratio
        colnames(threat_end_spp) <- gsub("_spp", "end_spp", colnames(threat_end_spp)) 
        threat <- cbind(round(100*threat, 1), threat_spp, threat_end_spp)
        
      } else {
        threat <- cbind(round(100*threat, 1), threat_spp)    
      }
  } else {
    threat <- round(100*threat, 1)
  }
  return(threat)
}


## getting the necessary vector to extrapolate threat
info$habitat_loss <- I(100 - info$ForestCover)

end_ratio <- info$number_of_endemic_plants_perc_Habel_etal_2019
end_richness <- info$number_of_endemic_plants_Habel_etal_2019
info$tree_richness <- as.double(sapply(strsplit(info$tree_species," "), function(x) x[1]))

trop.hotspots <- c("Atlantic Forest", "Caribbean Islands", 
                   "Coastal Forests of Eastern Africa", #"East Melanesian Islands",
                   "Eastern Afromontane", "Guinean Forests of West Africa",
                   "Indo-Burma", "Madagascar and the Indian Ocean Islands",
                   "Mesoamerica", "New Caledonia", "Philippines", 
                   "Sundaland", "Tropical Andes",
                   "Tumbes-Choco-Magdalena", "Wallacea", "Western Ghats and Sri Lanka")
tropical_hotspots <- info$hotspot.region %in% trop.hotspots
summary(info$tree_richness[tropical_hotspots]/
                      as.double(sapply(strsplit(info$plant_species_Myers_etal_2000," "), 
                                       function(x) x[1]))[tropical_hotspots],
                  na.rm = TRUE) # median 0.290, mean 0.313 [0.25-0.43]
info$tree_end_richness <- end_richness*0.3
richness <- as.double(sapply(strsplit(info$plant_species_Myers_etal_2000," "), function(x) x[1]))
info$tree_richness <- richness*0.3


## extrapolating threat
head(get_threat(info$habitat_loss, info$tree_richness, end_ratio), 3)
head(get_threat(info$habitat_loss, info$tree_end_richness), 3)
info$obs_threat_or_ext_endemics_trees <-
  as.double(info$observed_threatened_or_extinct_endemics_Brooks_etal_2002)*0.3
info$est_threat_or_ext_endemics_trees <-
  info$estimated_threatened_or_extinct_endemic_plants_Brooks_etal_2002*0.3
threatened_per_hotspot <- cbind.data.frame(info, get_threat(info$habitat_loss, info$tree_end_richness))
cols <- c("hotspot.region", "continent", "area_km2", #"domain",
          "plant_species_Myers_etal_2000", #"Endemic_plants_perc_Myers_etal_2000",
          "number_of_endemic_plants_Habel_etal_2019",
          "habitat_loss","tree_end_richness",
          "obs_threat_or_ext_endemics_trees", "est_threat_or_ext_endemics_trees", 
          "lin_threat_spp", "pow1_threat_spp", "pow2_threat_spp", "logis_threat_spp" )
threatened_per_hotspot_tropic <- threatened_per_hotspot[tropical_hotspots, cols]
head(threatened_per_hotspot_tropic, 3)

apply(threatened_per_hotspot_tropic[,c("lin_threat_spp", "pow1_threat_spp", 
                                       "pow2_threat_spp", "logis_threat_spp")], 2, sum, na.rm = TRUE)

cols1 <- c("hotspot.region","area_km2","ForestCover",
           "number_of_endemic_plants_Habel_etal_2019","tree_end_richness",
           "est_threat_or_ext_endemics_trees", 
           "lin_threat_spp","pow1_threat_spp","pow2_threat_spp","logis_threat_spp")
result3 <- threatened_per_hotspot[tropical_hotspots, cols1]
result3$est_range <- apply(apply(threatened_per_hotspot_tropic[,c("lin_threat_spp", "pow1_threat_spp", 
                                       "pow2_threat_spp", "logis_threat_spp")], 1, range, na.rm = TRUE),
      2, function(x) paste(ceiling(x), collapse =  "-"))
result4 <- result3[,!names(result3) %in% c("lin_threat_spp", "pow1_threat_spp", 
                       "pow2_threat_spp", "logis_threat_spp")]
result4$area_km2 <- ceiling(result4$area_km2)
result4$ForestCover <- round(result4$ForestCover, 1)
result4$tree_end_richness <- ceiling(result4$tree_end_richness)
result4$est_threat_or_ext_endemics_trees <- 
  round(result4$est_threat_or_ext_endemics_trees, 0)
result4$HabitatLoss <- 100 - result4$ForestCover

names(result4) <- c("Biodiversity Hotspot", "Total extent (km2)", 
                    "Forest Cover (%)a", "Endemic plantsb",
                    "Endemic treesc", "Estimated threate", "Estimated threatd")
result4 <- result4[, c("Biodiversity Hotspot", "Total extent (km2)", 
                       "Forest Cover (%)a", "Endemic plantsb",
                       "Endemic treesc", "Estimated threatd", "Estimated threate")] 
write.csv(result4, "data//Table2.csv")
