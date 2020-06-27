##############################################################
#### EFFECT OF TAXONOMIC CONFIDENCE OVER THE EOO ESTIMATE ####
##############################################################
rm(list=ls())

#devtools::install_github("gdauby/ConR", ref = "master", force = TRUE)
library("ConR")
require(rgeos)
source("./R/suggestions_for_ConR.r")
source("C://Users//renato//Documents//raflima//R_packages//ConR//R//EOO.sensitivity.R")
#source("C://Users//renato//Documents//raflima//R_packages//ConR//R//over.valid.poly.R")

## Reading the Neotropics shapefile ##
# neotrop.simp <- readRDS("data/Contour_Neotrop_simplified.rds")

## Reading herbarium data
# oc.data <- readRDS("data/threat_occ_data.rds")

#Putting data in the format demanded by the package
# MyData <- cbind.data.frame(ddlat = as.double(oc.data$latitude.work1),
#                            ddlon = as.double(oc.data$longitude.work1),
#                            tax = as.character(oc.data$species.correct2),
#                            tax.check2 = oc.data$tax.check2,
#                            stringsAsFactors = FALSE)
# rm(oc.data)

#### CALCULATING EOO USING DIFFERENT CONFIDENCE LEVELS OF THE OCCURRENCES ####
# teste <- MyData[grepl("Persea",MyData$tax), ] 
# sens <- EOO.sensitivity (MyData, 
#                          levels.order = c("FALSE", "cannot_check", "TRUE_TBC","TRUE_OTHER", "TRUE"),
#                             occ.based = TRUE,
#                             exclude.area = TRUE,
#                             country_map = neotrop.simp,
#                             method.range = "convex.hull",
#                             parallel = TRUE,
#                             NbeCores = 5,
#                             show_progress = TRUE,
#                             proj_user = 5641,
#                             value = "dist")
# sapply(sens, head)
# oc.data <- readRDS("data/threat_occ_data.rds")
# table(oc.data$species.correct2 == sens[[2]]$tax)
# oc.data$dist.eoo <- sens[[2]]$prop.dist.eoo
# saveRDS(oc.data,"data/threat_occ_data_new.rds")
# saveRDS(sens$EOO.change,"data/eoo.change_preliminar.rds")
# rm(oc.data, sens, MyData, neotrop.simp)

## Reading previously saved files
#resultado <- readRDS("data/eoo.change_preliminar.rds")

#### CLASSYFYING SPECIES ACCORDING TO 3 CRITERIA ####

### NEW CLASSIFICATION ###
tmp <- findInterval(resultado$Occs.level.5, c(5,15,30,75)) # number of occs indentied by family specialists
tmp0 <- findInterval(resultado$Occs.level.4, c(5,15,30,75)) # number of occs indentied by any specialists
tmp1 <- findInterval(resultado$Occs.level.5 / resultado$Occs.level.1, c(0.25,0.5,0.75,0.9)) # proportion of all occs indentied by faily specialists
tmp2 <- as.double(as.character(cut(resultado$EOO.increase.1, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens

## an inex of data quality
index <- apply(cbind(tmp,tmp1,tmp2), 1, sum, na.rm = TRUE)
table(index, useNA = "always")
resultado$index.conf <- index

## Assigning classes of taxonomic confidence
resultado$tax.conf <- NA
#class 1: best taxonomy (>90), many occurrences (> 75), little or no changes in EOO (do nothing)
resultado$tax.conf[tmp == 4] <- "class1"
#class 2: good taxonomy (75), good amount of occurrences (>30)
resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 %in% 4 & tmp1>=3] <- "class2"
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 3 & tmp1>=3] <- "class2"
#class 3: low taxonomy (>50), feasible amount of occurrences (>15)
resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 >= 3 & tmp1>=2] <- "class3"
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & tmp1>=2] <- "class3"
#class 4: low taxonomy (>50), bad amount of occurrences (<15)
resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 >= 2 & tmp1>=1] <- "class4"
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & tmp1>=1] <- "class4"
#class 5: wort taxonomy (<25)
resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 >= 2 & tmp1 %in% 0] <- "class5"
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & tmp1 %in% 0] <- "class5"
#class 6: very few (<5) or no occurrences validated by an specialist
resultado$tax.conf[is.na(resultado$tax.conf) & is.na(tmp)] <- "class6"
resultado$tax.conf[!is.na(tmp) & tmp %in% 0] <- "class6"

##Inspecting
table(resultado$tax.conf, useNA = "always")
par(mfrow=c(2,2), mar=c(3,3,2,1))
boxplot(index ~ resultado$tax.conf,
        varwidth = TRUE, notch = TRUE, main = "index x classes")
boxplot(resultado$Occs.level.5 / resultado$Occs.level.1 ~ resultado$tax.conf,
        varwidth = TRUE, notch = TRUE, main = "prop. true x classes")
boxplot(jitter(tmp0) ~ resultado$tax.conf,
        varwidth = TRUE, notch = TRUE, main = "prop. true_other x classes")
boxplot(log(resultado$EOO.increase.1+1) ~ resultado$tax.conf,
        varwidth = TRUE, notch = TRUE, main = "EOO increase x classes")

## Which species we can include all TRUE_TBC automatically:
oc.data <- readRDS("data/threat_occ_data_new.rds")
tbc <- unique(oc.data$species.correct2[oc.data$tax.check2 %in% "TRUE_TBC"])
resultado$tax.conf[resultado$Species %in% tbc & !is.na(tmp) & tmp<4 & tmp0 >= 3] = "class2"
resultado$tax.conf[resultado$Species %in% tbc & !is.na(tmp) & tmp<3 & tmp0 %in% 2] = "class3"
resultado$tax.conf[resultado$Species %in% tbc & !is.na(tmp) & tmp<2 & tmp0 %in% 1] = "class4"
saveRDS(resultado, "data/eoo.change_preliminar_new.rds")

#### PERFORMING THE FINAL CLASSIFICATIONS OF TAXONOMIC CONFIDENCE ####
require(data.table)
rm(list=ls()); gc()

#Loading the data and the species classifications and merging
oc.data <- data.table(readRDS("data/threat_occ_data_new.rds"))
resultado <- data.table(readRDS("data/eoo.change_preliminar_new.rds"))
oc.data <- merge(oc.data, resultado[,c("Species","tax.conf")], 
                 by.x = "species.correct2", by.y = "Species", all.x = TRUE, sort = FALSE)

#Creating the final column for the assessments by class of taxonomical confidence (see colorful table in the appendix)
oc.data[, tax.check.final := tax.check2]
oc.data[tax.check.final %in% "cannot_check", tax.check.final := "FALSE"]

table(oc.data[,tax.check.final])
oc.data[!tax.conf %in% "class1" & 
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC") & 
          dist.eoo %in% 0, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])
oc.data[!tax.conf %in% "class1" & 
          tax.conf %in% c("class3", "class4", "class5", "class6") & 
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC", "FALSE") & 
          dist.eoo %in% 0, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])
oc.data[!tax.conf %in% "class1" & 
          tax.conf %in% c("class4", "class5", "class6") & 
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC", "FALSE") & 
          dist.eoo < 0.05, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])
oc.data[!tax.conf %in% "class1" & 
          tax.conf %in% c("class5", "class6") & 
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC", "FALSE") & 
          dist.eoo < 0.1, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])
oc.data[tax.conf %in% "class6" & 
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC", "FALSE") & 
          dist.eoo < 0.25, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])
oc.data[tax.conf %in% "class6" & 
          tax.check2 %in% c("TRUE_OTHER") & 
          dist.eoo < 0.5, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])
oc.data[tax.conf %in% "class6" & 
          tax.check2 %in% c("TRUE_OTHER") & 
          is.na(dist.eoo), tax.check.final := "TRUE"]
tmp <- table(oc.data[tax.conf %in% "class6" & tax.check.final %in% "TRUE", species.correct2])
spp <- names(tmp[tmp<5])
oc.data[tax.conf %in% "class6" & 
          species.correct2 %in% spp &
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC", "FALSE") & 
          dist.eoo < 1, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])

tbc <- unique(oc.data$species.correct2[oc.data$tax.check2 %in% "TRUE_TBC"])
tmp <- table(oc.data[species.correct2 %in% tbc & tax.check.final %in% "TRUE", species.correct2])
spp <- names(tmp[tmp<30])
oc.data[species.correct2 %in% spp &
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC"), tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])

## Final edits and saving the data with the new column of tax.check
oc.data[tax.check.final %in% c("TRUE_OTHER","TRUE_TBC"), tax.check.final := "medium"]
oc.data[tax.check.final %in% c("FALSE"), tax.check.final := "low"]
oc.data[tax.check.final %in% c("TRUE"), tax.check.final := "high"]
saveRDS(oc.data,"data/threat_occ_data_final.rds")

### OLD CLASSIFICATION ###
# tmp <- findInterval(resultado$Occs.level.5, c(15,30,75)) # number of occs indentied by family specialists
# tmp0 <- findInterval(resultado$Occs.level.4, c(15,30,75)) # number of occs indentied by any specialists
# tmp1 <- findInterval(resultado$Occs.level.5 / resultado$Occs.level.1, c(0.5,0.75,0.9)) # proportion of all occs indentied by faily specialists
# tmp2 <- as.double(as.character(cut(resultado$EOO.increase.4, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens
# tmp3 <- as.double(as.character(cut(resultado$EOO.increase.2, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens
# tmp4 <- as.double(as.character(cut(resultado$EOO.increase.1, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens
# table(tmp, tmp1)
# table(tmp, tmp3)
# 
# boxplot(jitter(tmp2) ~ tmp1, varwidth = TRUE, notch = TRUE, 
#         xlab = "Tax/all", ylab = "EOO change (other tax or better)") # increase from high to medium in respect to the prop. of occs identfied by specialists
# boxplot(jitter(tmp3) ~ tmp1, varwidth = TRUE, notch = TRUE,
#         xlab = "Tax/all", ylab = "EOO change(cannot check or better)") # increase from high to low in respect to the prop. of occs identfied by specialists
# boxplot(jitter(tmp4) ~ tmp1, varwidth = TRUE, notch = TRUE,
#         xlab = "Tax/all", ylab = "EOO change (all occs)") # increase from high to bad in respect to the prop. of occs identfied by specialists
# 
# 
# ## an inex of data quality
# index <- apply(cbind(tmp,tmp1,tmp4), 1, sum, na.rm = TRUE)
# index1 <- round(apply(cbind(tmp,tmp1,apply(cbind(tmp2,tmp3,tmp4),1,mean,na.rm=TRUE)), 1, sum, na.rm = TRUE),1)
# plot(jitter(index1) ~ index)
# table(index)
# 
# ## Which species we can include all TRUE_OTHER automatically:
# resultado$add.other <- FALSE
# resultado$add.other[!is.na(tmp2) & (tmp2 >= 2) ] <- TRUE
# resultado$add.other[!is.na(tmp2) & (tmp1 >= 2) ] <- TRUE
# resultado$add.other[is.na(tmp2)] <- TRUE
# oc.data <- readRDS("data/threat_occ_data.rds")
# tbc <- unique(oc.data$species.correct2[oc.data$tax.check2 %in% "TRUE_TBC"])
# resultado$add.other[resultado$tax %in% tbc & (tmp3 >= 1) ] <- TRUE
# rm(oc.data)
# 
# ## Which species we can include all TRUE_TBC automatically:
# resultado$add.tbc <- FALSE
# resultado$add.tbc[resultado$tax %in% tbc & (tmp<3 | is.na(tmp3)) & (tmp1 >=1 | tmp3 >=1)] <- TRUE
# 
# ## Assigning classes of taxonomic confidence
# resultado$tax.conf <- NA
# resultado$action <- NA
# #class 1: best taxonomy, many occurrences (> 75), little or no changes in EOO
# resultado$tax.conf[tmp == 3] <- "class1"
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 %in% 3 & (tmp1 %in% 3 | tmp2 %in% 3)] <- "class1.1" # TRUE + TRUE_OTHER >=75 + best tax or small EOO change
# resultado$action[resultado$tax.conf == "class1"] <- "nothing"
# #class 2: good taxonomy
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & (tmp1>=3 | tmp4>=3)] <- "class2"   # good occs + best tax or small EOO change
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & (tmp1>=2 & tmp4>=2)] <- "class2.1" # good occs + good tax + low EOO change
# #class 3
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & (tmp1>=2 | tmp4>=2)] <- "class3"   # good occs + good tax or low EOO change
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & (tmp1>=2 & tmp4>=2)] <- "class3.1" # low occs + good tax + low EOO change
# #class 4 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & (tmp1>=2 | tmp4>=2)] <- "class4.1"   # low occs + good tax or low EOO change
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & (tmp1>=1 & tmp4>=1)] <- "class4.2" # low occs + low tax + considerable EOO change 
# #class 5 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2] <- "class5" # all other good
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=3 | tmp4>=3)] <- "class5.1"   # bad occs + best tax or small EOO change 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=2 & tmp4>=2)] <- "class5.2" # bad occs + good tax + low EOO change 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=2 | tmp4>=2)] <- "class5.3"   # bad occs + good tax or small EOO change 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=1 & tmp4>=1)] <- "class5.4"  # bad occs + low tax + considerable EOO change 
# #class 6
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1] <- "class6"   # all other low occs 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=1 | tmp4>=1)] <- "class6.1"   # bad occs + low tax or considerable EOO change 
# #class 7
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp == 0 & resultado$true >5] <- "class7"   # bad occs + low tax or considerable EOO change 
# resultado$tax.conf[is.na(resultado$tax.conf) & tmp == 0] <- "class7.1" # all other bad occs
# #class 8: no occurrences validated by a family taxonomist
# resultado$tax.conf[is.na(resultado$tax.conf) & is.na(tmp)] <- "class8"



