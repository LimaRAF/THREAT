##############################################################
##############################################################
#### COMBINING ASSESSMENTS FROM ALL CRITERIA - A, B, C, D ####
##############################################################
##############################################################
rm(list=ls())

## Loading packages ##


#####################
#### PRE-EDITING ####
#####################
## Loading and merging the assessments ##
#Criterion A
critA <- readRDS("data/criterionA_all_GLs.rds")
critA <- critA[,c("species","assessment.period","reduction_A12","A1","A2","category_A","category_A_code","reduction_A12.25ys","A1.25ys","A2.25ys")]
critA$reduction_A12 <- round(as.double(critA$reduction_A12),2)
critA$reduction_A12.25ys <- round(as.double(critA$reduction_A12.25ys),2)

#Criterion B
critB_high <- readRDS("data/criterionB_all_high_tax_confidence.rds")
#critB_high <- critB_high[,c("tax","EOO","AOO","Nbe_subPop","nbe_loc_total","Nbe_occs","B1","B2","category_B_code","category_B")]
critB_high <- critB_high[,c("tax","EOO","AOO","Nbe_subPop","nbe_loc_total","Nbe_occs","category_B_code","category_B")]
names(critB_high)[1] <- "species"

critB_low <- readRDS("data/criterionB_all_low_tax_confidence.rds")
#critB_low <- critB_low[,c("tax","EOO","AOO","Nbe_subPop","nbe_loc_total","Nbe_occs","B1","B2","category_B_code","category_B")]
critB_low <- critB_low[,c("tax","EOO","AOO","Nbe_subPop","nbe_loc_total","Nbe_occs","category_B_code","category_B")]
names(critB_low)[1] <- "species"
#critB <- merge(critB_high, critB_low[,c("species","EOO","AOO","Nbe_subPop","nbe_loc_total","B1","B2","category_B_code","category_B")], by = "species", all = TRUE, sort= FALSE)
critB <- merge(critB_high, critB_low[,c("species","EOO","AOO","Nbe_subPop","nbe_loc_total","category_B_code","category_B")], by = "species", all = TRUE, sort= FALSE)
names(critB)[grepl("\\.x", names(critB))] <- gsub("\\.x","",names(critB)[grepl("\\.x", names(critB))])
names(critB)[grepl("\\.y", names(critB))] <- gsub("\\.y",".low.tax",names(critB)[grepl("\\.y", names(critB))])

#Criterion C
critC <- readRDS("data/criterionC_all_prop_mature.rds")
critC <- critC[,c("species","any.decline","cont.decline","C1","C2","category_C","category_C_code","C1.p0.45")]

#Criterion D
critD <- readRDS("data/criterionD_all_prop_mature.rds")
critD <- critD[,c("species","pop.size","D","D.AOO","D2.Loc","category_D","category_D_code","D.p0.45")]
critD$pop.size <- round(as.double(critD$pop.size),1)

## Merging all assessments
all.crit <- merge(critA, critB, by = "species", all = TRUE)
all.crit <- merge(all.crit, critC, by = "species", all = TRUE)
all.crit <- merge(all.crit, critD, by = "species", all = TRUE)
all.crit <- all.crit[order(all.crit$species),]

#########################################
#### GETTING SPECIES ENDEMISM LEVELS ####
#########################################

## Getting endemism levels from Lima et al. (2020)
endemism <- readRDS("data/threat_endemism_levels.rds") # endemism levels from Lima et al. (2020)
all.crit <- merge(all.crit, endemism, by = "species", all.x = TRUE)
#Getting missing endemism status
SA.spp <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
tmp <- all.crit[is.na(all.crit$endemic),]
tmp1 <- merge(tmp[,c("species","endemic")] , 
             SA.spp[,c("SpeciesKTSA","geographical_distribution","domain.reflora","endemism.Reflora",
                       "Domain.Ary","Endemism.Ary","AtlanticForest","Domain","Endemism")],
                        by.x="species", by.y= "SpeciesKTSA", all.x = TRUE)
tmp1$Domain[is.na(tmp1$Domain)] <- 
  tmp1$domain.reflora[is.na(tmp1$Domain)]
tmp1$Domain[tmp1$Domain %in% "no_records"] <- 
  tmp1$domain.reflora[tmp1$Domain %in% "no_records"]
#Homogenizing the informations for missing endemism levels
tmp1$endemic <- tmp1$Endemism
tmp1$endemic[(is.na(tmp1$endemic) | tmp1$endemic %in% "endemic?") &
               tmp1$Domain %in% "AtlanticForest"] <- "endemic"
tmp1$endemic[(is.na(tmp1$endemic) | tmp1$endemic %in% "not_endemic") &
               !grepl("AtlanticForest", tmp1$Domain)] <- "occasional"
tmp1$endemic[(is.na(tmp1$endemic) | tmp1$endemic %in% "not_endemic") & 
               grepl("AtlanticForest", tmp1$Domain) & 
               stringr::str_count(tmp1$Domain, "\\|") >2] <- "widespread_sparse"
tmp1$endemic[(is.na(tmp1$endemic) | tmp1$endemic %in% "not_endemic") & 
               grepl("AtlanticForest", tmp1$Domain) & 
               stringr::str_count(tmp1$Domain, "\\|") >=1] <- "widespread_common"
tmp1$endemic[tmp1$Domain %in% c("Amazonia","Cerrado","Caatinga")] <- "occasional"
tmp1$endemic[!tmp1$endemic %in% "endemic" & 
               tmp1$geographical_distribution %in% c("regional_endemic","local_endemic")] <- "widespread_common"
table(tmp1$endemic, useNA = "always")
all.crit$endemic[is.na(all.crit$endemic)] <- tmp1$endemic
  
#Homogenizing all categories
table(all.crit$endemic, useNA = "always")
all.crit$endemic[all.crit$endemic %in% "not in the AF"] <- "occasional"
all.crit$endemic[all.crit$endemic %in% "not_endemic"] <- "widespread_sparse"


###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
#########################
#### NEAR THREATENED ####
#########################



###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
##############################
#### CONSENSUS ASSESSMENT ####
##############################

#tmp <- all.crit[,c("species", "A1", "A2", "B1","B2","C1", "C2", "D")]
tmp <- all.crit[,c("species","A1", "A2", "category_B","C1", "C2", "D")]
tmp[] <- lapply(tmp, gsub, pattern = "^LC or NT$", replacement = 0)
tmp[] <- lapply(tmp, gsub, pattern = "^LC$", replacement = 0)
tmp[] <- lapply(tmp, gsub, pattern = "^NT$", replacement = 1)
tmp[] <- lapply(tmp, gsub, pattern = "^VU$", replacement = 2)
tmp[] <- lapply(tmp, gsub, pattern = "^EN$", replacement = 3)
tmp[] <- lapply(tmp, gsub, pattern = "^CR$", replacement = 4)
tmp[,2:7] <- apply(tmp[,2:7], 2, as.double)
tmp$cat <- as.character(apply(tmp[,2:7], 1, max, na.rm=TRUE))
tmp$crit.main <- 
  apply(tmp[,2:7], 1, 
        function(x) paste(names(x)[which(x == max(x, na.rm = TRUE))], collapse="+"))
rpl.cds <- c("CR", "EN", "VU", "NT", "LC", "LC or NT")
names(rpl.cds) <- c("4", "3", "2", "1", "0", "0")
tmp$crit.aux <- 
  apply(tmp[,2:7], 1, function(x) 
                      if(length(unique(x)) > 1) {
                        crits <- sort(unique(as.double(x)), decreasing = TRUE)[-1]
                        cats <- sapply(crits, function(y) paste(names(tmp[,2:7])[x==y], collapse = "+"))
                        paste(paste(stringr::str_replace_all(crits, rpl.cds),": ", cats, sep=""), collapse = "; ")
                      } else {
                        ""
                      }  
      )
tmp$crit.aux <- gsub("^: $", "", tmp$crit.aux)
tmp$cat <- stringr::str_replace_all(tmp$cat,  rpl.cds)
head(tmp)



###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
#########################
#### DATA DEFICIENT  ####
#########################

## All species with < 3 occurrences and no abundance data at all
## We also classified as DD species known only from its type locality

###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
##############################
#### PREVIOUS ASSESSMENTS ####
##############################

2*2