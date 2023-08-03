##############################################################
#### EFFECT OF TAXONOMIC CONFIDENCE OVER THE EOO ESTIMATE ####
##############################################################
rm(list=ls())
gc()
require(rgeos)
require(sf)
require(fields)
require(dplyr)
require(data.table)
source("R/99_functions.R")

## Reading herbarium data
oc.data <- readRDS("data/threat_occ_data.rds")
## Reading/editing inventory data
inv.data <- readRDS("data/threat_inventory_data.rds")
inv.data$year_data[is.na(inv.data$year_data)] <-
  inv.data$year[is.na(inv.data$year_data)]
inv.data$numTombo[is.na(inv.data$numTombo)] <-
  paste0("ordem_",inv.data$ordem[is.na(inv.data$numTombo)])
inv.data$DetDate <- gsub("\\\n","",inv.data$DetDate)
inv.data$DetDate <- gsub("\\|s\\.d\\.$|^s\\.d\\.\\|","",inv.data$DetDate)
inv.data$DetDate <- plantR::getYear(inv.data$DetDate)

familias <- inv.data$family
multiple.fam <- lengths(familias) > 1
familias[multiple.fam] <- lapply(familias[multiple.fam], function(x) x[!is.na(x)])
familias <- unlist(familias)

#Putting data in the format demanded by the ConR package
MyData <- cbind.data.frame(ddlat = as.double(c(oc.data$latitude.work1, inv.data$lat1)),
                           ddlon = as.double(c(oc.data$longitude.work1, inv.data$long1)),
                           tax = as.character(c(oc.data$species.correct2, inv.data$species.correct2)),
                           higher.tax.rank = c(oc.data$family.correct1, familias),
                           coly = as.double(c(oc.data$ano, inv.data$year_data)),
                           vouchers = c(oc.data$dup.ID1, inv.data$numTombo),
                           typeStatus = c(oc.data$typeStatus, inv.data$typeStatus),
                           detBy = c(oc.data$determinador.name, inv.data$DetBy),
                           dety = c(oc.data$ano.det, inv.data$DetDate),
                           tax.check2 = c(oc.data$tax.check2, inv.data$tax_check2),
                           UC = c(oc.data$UC, inv.data$UC),
                           source = c(rep("herbaria",dim(oc.data)[1]), rep("treeco",dim(inv.data)[1])),
                           stringsAsFactors = FALSE)
rm(oc.data, inv.data)

#### CALCULATING EOO USING DIFFERENT CONFIDENCE LEVELS OF THE OCCURRENCES ####
#teste <- MyData[grepl("Persea",MyData$tax), ]
sens <- my.EOO.sensitivity (
                            #teste[,c("ddlat", "ddlon", "tax", "tax.check2")],
                            MyData[,c("ddlat", "ddlon", "tax", "tax.check2")],
                            levels.order = c("FALSE", "cannot_check", "TRUE_TBC","TRUE_OTHER", "TRUE"),
                            occ.based = TRUE,
                            exclude.area = FALSE,
                            country_map = NULL,
                            method.range = "convex.hull",
                            parallel = TRUE,
                            NbeCores = 6,
                            show_progress = TRUE,
                            proj_user = 5641,
                            value = "dist")
sapply(sens, head)
table(MyData$tax == sens[[2]]$tax)
MyData$dist.eoo <- sens[[2]]$prop.dist.eoo
saveRDS(MyData, "data/threat_occ_data_new.rds", compress = "gzip")
saveRDS(sens$EOO.change, "data/eoo.change_preliminar_uncropped.rds", compress = "gzip")
rm(sens, MyData)

## Reading previously saved files
resultado <- readRDS("data/eoo.change_preliminar_uncropped.rds")

#### CLASSIFYING SPECIES ACCORDING TO 3 CRITERIA ####

## NEW CLASSIFICATION ##
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
resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 >= 4 & tmp1>=3] <- "class2"
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
# table(resultado$tax.conf, useNA = "always")
# par(mfrow=c(2,2), mar=c(3,3,2,1))
# boxplot(index ~ resultado$tax.conf,
#         varwidth = TRUE, notch = TRUE, main = "index x classes")
# boxplot(resultado$Occs.level.5 / resultado$Occs.level.1 ~ resultado$tax.conf,
#         varwidth = TRUE, notch = TRUE, main = "prop. true x classes")
# boxplot(jitter(tmp0) ~ resultado$tax.conf,
#         varwidth = TRUE, notch = TRUE, main = "prop. true_other x classes")
# boxplot(log(resultado$EOO.increase.1+1) ~ resultado$tax.conf,
#         varwidth = TRUE, notch = TRUE, main = "EOO increase x classes")

## Which species we can include all TRUE_TBC automatically:
oc.data <- readRDS("data/threat_occ_data_new.rds")
tbc <- unique(oc.data$tax[oc.data$tax.check2 %in% "TRUE_TBC"])
resultado$tax.conf[resultado$Species %in% tbc & !is.na(tmp) & tmp<4 & tmp0 >= 3] = "class2"
resultado$tax.conf[resultado$Species %in% tbc & !is.na(tmp) & tmp<3 & tmp0 %in% 2] = "class3"
resultado$tax.conf[resultado$Species %in% tbc & !is.na(tmp) & tmp<2 & tmp0 %in% 1] = "class4"
saveRDS(resultado, "data/eoo.change_preliminar_new.rds")

#### PERFORMING THE FINAL CLASSIFICATIONS OF TAXONOMIC CONFIDENCE ####
rm(list=ls()); gc()

#Loading the data and the species classifications and merging
oc.data <- data.table(readRDS("data/threat_occ_data_new.rds"))
resultado <- data.table(readRDS("data/eoo.change_preliminar_new.rds"))
oc.data <- merge(oc.data, resultado[,c("Species","tax.conf")], 
                 by.x = "tax", by.y = "Species", all.x = TRUE, sort = FALSE)

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
tmp <- table(oc.data[tax.conf %in% "class6" & tax.check.final %in% "TRUE", tax])
spp <- names(tmp[tmp<5])
oc.data[tax.conf %in% "class6" & 
          tax %in% spp &
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC", "FALSE") & 
          dist.eoo < 1, tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])

tbc <- unique(oc.data$tax[oc.data$tax.check2 %in% "TRUE_TBC"])
tmp <- table(oc.data[tax %in% tbc & tax.check.final %in% "TRUE", tax])
spp <- names(tmp[tmp<30])
oc.data[tax %in% spp &
          tax.check2 %in% c("TRUE_OTHER", "TRUE_TBC"), tax.check.final := "TRUE"]
table(oc.data[,tax.check.final])

## Final edits and saving the data with the new column of tax.check
oc.data[tax.check.final %in% c("TRUE_OTHER","TRUE_TBC"), tax.check.final := "medium"]
oc.data[tax.check.final %in% c("FALSE"), tax.check.final := "low"]
oc.data[tax.check.final %in% c("TRUE"), tax.check.final := "high"]
saveRDS(oc.data,"data/threat_occ_data_final.rds", compress = "gzip")
rm(list = ls())
