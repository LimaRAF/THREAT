#############################
#### APPLYING CRITERIA B ####
#############################
rm(list = ls())
gc()

#### LOADING PACKAGES ###
require(ConR)
require(dplyr)
require(circlize)

#### LOADING ACCESSORY FUNCTIONS ###
source("./R/99_functions.R")


#Getting the estimates for criterion B
critB_opt <- readRDS("data/criteriaB_metrics_optim_confidence.rds")
critB_opt$tax <- as.character(critB_opt$tax)
critB_low <- readRDS("data/criteriaB_metrics_low_confidence.rds")
critB_low$tax <- as.character(critB_low$tax)
critB_high <- readRDS("data/criteriaB_metrics_high_confidence.rds")
critB_high$tax <- as.character(critB_high$tax)

#Getting estimates of species continuing decline based on the index of abundance (Criterion C)
critC <- readRDS("data/criterionC_all_prop_mature.rds")
est.decline <- sapply(strsplit(critC$cont.decline,"\\|"), tail, 1)  
est.decline <- gsub("\\(|\\)|[0-9]", "", est.decline)
est.decline <- gsub(" -", "", est.decline)
table(est.decline, critC$any.decline, useNA = "always")
critC$declineC <- critC$any.decline
critC$declineC[!critC$declineC %in% "Decreasing" & est.decline %in% "Decreasing"] <- "Decreasing" 
critC$declineC[critC$declineC %in% "Increasing" | critC$declineC %in% "Stable"] <- "Not Decreasing" 

#Getting estimates of species continuing decline based on habitat loss (Criterion B)
eoo.decline <- readRDS("data/EOO_hab_loss_2000_2015_uncropped.rds")
eoo.decline$declineB <- eoo.decline$rel.loss
eoo.decline$declineB[!is.na(eoo.decline$declineB) & eoo.decline$declineB >= 1] <- 1
eoo.decline$declineB[!is.na(eoo.decline$declineB) & as.double(eoo.decline$declineB) < 1 & as.double(eoo.decline$declineB) >= 0] <- 0
eoo.decline$declineB[is.nan(eoo.decline$declineB)] <- NA
eoo.decline$declineB[eoo.decline$declineB %in% 1] <- "Decreasing"
eoo.decline$declineB[eoo.decline$declineB %in% 0] <- "Not Decreasing"
table(eoo.decline$declineB, useNA = "always")

#Correcting the decrease for pioneer species
hab <- readRDS("data/threat_habitats.rds")
tmp <- merge(eoo.decline, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
tmp <- tmp[order(tmp$tax),]
table(tmp$tax == eoo.decline$tax)
eoo.decline$net.loss <- eoo.decline$recover - eoo.decline$loss 
hist(eoo.decline$net.loss, nclass=80)
eoo.decline$declineB[eoo.decline$net.loss >=0.1 & tmp$ecol.group %in% "pioneer"] <- "Not Decreasing"
table(eoo.decline$declineB, useNA = "always")

#Merging any decline info with the criterion B metrics
critB_opt <- merge(critB_opt, critC[,c("species","declineC")],
                   by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_opt <- merge(critB_opt, eoo.decline[,c("tax","declineB")],
                   by = "tax", all.x = TRUE, sort = FALSE)
critB_opt <- critB_opt[order(critB_opt$tax),]

critB_low <- merge(critB_low, critC[,c("species","declineC")],
                   by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_low <- merge(critB_low, eoo.decline[,c("tax","declineB")],
                   by = "tax", all.x = TRUE, sort = FALSE)
critB_low <- critB_low[order(critB_low$tax),]

critB_high <- merge(critB_high, critC[,c("species","declineC")],
                    by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB_high <- merge(critB_high, eoo.decline[,c("tax","declineB")],
                    by = "tax", all.x = TRUE, sort = FALSE)
critB_high <- critB_high[order(critB_high$tax),]

#Comparing the decline from criterion C and B and correcting if necessary (priority for abudance decline over EOO decline) For assuming, that all other species (rarer) are decreasing as well
hab <- readRDS("data/threat_habitats.rds")
tmp <- merge(critB_opt, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
table(critB_opt$tax == tmp$tax)
critB_opt$declineB[critB_opt$declineB %in% "Not Decreasing" & critB_opt$declineC %in% "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
#table(critB_opt$declineB, critB_opt$declineC, useNA = "always")
tmp <- merge(critB_low, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
table(critB_low$tax == tmp$tax)
critB_low$declineB[critB_low$declineB %in% "Not Decreasing" & critB_low$declineC %in% "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
#table(critB_low$declineB, critB_low$declineC, useNA = "always")
tmp <- merge(critB_high, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
table(critB_high$tax == tmp$tax)
critB_high$declineB[critB_high$declineB %in% "Not Decreasing" & critB_high$declineC %in% "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
#table(critB_high$declineB, critB_high$declineC, useNA = "always")

#Combining the info on number of localities and % in PAs
critB_opt <- 
  critB_opt %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)
critB_low <- 
  critB_low %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)
critB_high <- 
  critB_high %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)

#Creating the severely fragmented colum
critB_opt$sever.frag <- 100 * critB_opt$Nbe_subPop/ (critB_opt$AOO / 4) > 50
critB_low$sever.frag <- 100 * critB_low$Nbe_subPop/ (critB_low$AOO / 4) > 50
critB_high$sever.frag <- 100 * critB_high$Nbe_subPop/ (critB_high$AOO / 4) > 50


#### PERFORMING THE ASSESSMENTES OF CRITERION B ####
results_Cb_opt <- cat_criterion_b(EOO = critB_opt$EOO,
                                  AOO = critB_opt$AOO,
                                  locations = critB_opt$nbe_loc_total,
                                  sever.frag = critB_opt$sever.frag,
                                  #protected = critB_opt$protected, 
                                  decline = critB_opt$declineB,
                                  protected.threshold = 100 
)
sum((100 * table(results_Cb_opt$ranks_B, useNA = "always")/dim(critB_opt)[1])[c(1,2,4)]) #17.25%; now 15.65%

results_Cb_low <- cat_criterion_b(EOO = critB_low$EOO,
                                  AOO = critB_low$AOO,
                                  locations = critB_low$nbe_loc_total,
                                  sever.frag = critB_low$sever.frag,
                                  #protected = critB_low$protected, 
                                  decline = critB_low$declineB,
                                  protected.threshold = 100
)
sum((100 * table(results_Cb_low$ranks_B, useNA = "always")/dim(critB_low)[1])[c(1,2,4)]) #14.68%; now 13.57%

results_Cb_high <- cat_criterion_b(EOO = critB_high$EOO,
                                   AOO = critB_high$AOO,
                                   locations = critB_high$nbe_loc_total,
                                   sever.frag = critB_high$sever.frag,
                                   #protected = critB_high$protected, 
                                   decline = critB_high$declineB,
                                   protected.threshold = 100 
)
sum((100 * table(results_Cb_high$ranks_B, useNA = "always")/dim(critB_high)[1])[c(1,2,4)]) #26.95%; now 24.99%


#Saving the results
results_Cb_opt <- do.call(cbind.data.frame, c(results_Cb_opt, stringsAsFactors = FALSE))
critB_opt.all <- critB_opt[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_opt.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_opt.all[, c("category_B", "category_B_code","category_B1","category_B2")] <-
  results_Cb_opt
critB_opt.all[is.na(critB_opt.all$AOO),]

results_Cb_low <- do.call(cbind.data.frame, c(results_Cb_low, stringsAsFactors = FALSE))
critB_low.all <- critB_low[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_low.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_low.all[, c("category_B", "category_B_code","category_B1","category_B2")] <-
  results_Cb_low
critB_low.all[is.na(critB_low.all$AOO),]

results_Cb_high <- do.call(cbind.data.frame, c(results_Cb_high, stringsAsFactors = FALSE))
critB_high.all <- critB_high[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_high.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_high.all[, c("category_B", "category_B_code","category_B1","category_B2")] <-
  results_Cb_high
critB_high.all[is.na(critB_high.all$AOO),]

#Saving
saveRDS(as.data.frame(critB_opt.all), "data/criterionB_all_optim_tax_confidence.rds")
saveRDS(as.data.frame(critB_low.all), "data/criterionB_all_low_tax_confidence.rds")
saveRDS(as.data.frame(critB_high.all), "data/criterionB_all_high_tax_confidence.rds")
rm(list = ls())



