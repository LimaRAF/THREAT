##############################################################
##############################################################
#### COMBINING ASSESSMENTS FROM ALL CRITERIA - A, B, C, D ####
##############################################################
##############################################################
rm(list=ls())

## Loading packages ##
# detach("package:ConR", unload=TRUE)
# install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
#                  repos = NULL, 
#                  type = "source")
library("ConR")
require(data.table)

#####################
#### PRE-EDITING ####
#####################
## Loading and merging the assessments ##
#Criterion A
critA <- readRDS("data/criterionA_all_GLs.rds")
critA <- critA[,c("species","assessment.period","reduction_A12","category_A","category_A_code",
                  "reduction_A12_obs","basis_d",
                  "reduction_A12.10ys","reduction_A12.15ys","reduction_A12.20ys",
                  "reduction_A12.25ys","reduction_A12.30ys","reduction_A12.35ys",
                  "reduction_A12.40ys","reduction_A12.45ys","reduction_A12.50ys", 
                  "A2.20ys", "A2.25ys")]
critA$reduction_A12 <- round(as.double(critA$reduction_A12), 2)
critA$reduction_A12.10ys <- round(as.double(critA$reduction_A12.10ys), 2)
critA$reduction_A12.15ys <- round(as.double(critA$reduction_A12.15ys), 2)
critA$reduction_A12.20ys <- round(as.double(critA$reduction_A12.20ys), 2)
critA$reduction_A12.25ys <- round(as.double(critA$reduction_A12.25ys), 2)
critA$reduction_A12.30ys <- round(as.double(critA$reduction_A12.30ys), 2)
critA$reduction_A12.35ys <- round(as.double(critA$reduction_A12.35ys), 2)
critA$reduction_A12.40ys <- round(as.double(critA$reduction_A12.40ys), 2)
critA$reduction_A12.45ys <- round(as.double(critA$reduction_A12.45ys), 2)
critA$reduction_A12.50ys <- round(as.double(critA$reduction_A12.50ys), 2)
critA$A2 <- critA$category_A

#Criterion B
critB_opt <- readRDS("data/criterionB_all_optim_tax_confidence.rds")
critB_opt <- critB_opt[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag","category_B1","category_B2","category_B","category_B_code")]

critB_low <- readRDS("data/criterionB_all_low_tax_confidence.rds")
critB_low <- critB_low[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag","category_B1","category_B2","category_B","category_B_code")]

critB_high <- readRDS("data/criterionB_all_high_tax_confidence.rds")
critB_high <- critB_high[,c("tax","Nbe_occs","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag","category_B1","category_B2","category_B","category_B_code")]

names(critB_high)[1] <- names(critB_opt)[1] <- names(critB_low)[1] <- "species"
names(critB_high)[10:11] <- names(critB_opt)[10:11] <- names(critB_low)[10:11] <- c("B1","B2")

critB <- merge(critB_opt, critB_low[,c("species","EOO","AOO","B1","B2","category_B_code","category_B")], by = "species", all = TRUE, sort= FALSE)
names(critB)[grepl("\\.x", names(critB))] <- gsub("\\.x","",names(critB)[grepl("\\.x", names(critB))])
names(critB)[grepl("\\.y", names(critB))] <- gsub("\\.y",".low.tax",names(critB)[grepl("\\.y", names(critB))])
critB <- merge(critB, critB_high[,c("species","EOO","AOO","B1","B2","category_B_code","category_B")], by = "species", all = TRUE, sort= FALSE)
names(critB)[grepl("\\.x", names(critB))] <- gsub("\\.x","",names(critB)[grepl("\\.x", names(critB))])
names(critB)[grepl("\\.y", names(critB))] <- gsub("\\.y",".high.tax",names(critB)[grepl("\\.y", names(critB))])


#Criterion C
critC <- readRDS("data/criterionC_all_prop_mature.rds")
# critC <- critC[,c("species","any.decline","cont.decline","C1","C2","category_C","category_C_code","C1.p0.45")]
critC <- critC[,c("species","any.decline","cont.decline","assess.pop.size",
                  "reduction_3gen","reduction_2gen","reduction_1gen","reduction_obs",
                  "C1","C2","category_C","category_C_code",
                  "C2.p0.18","C2.p0.25","C2.p0.33","C2.p0.45", "C2.p0.51")]
replace_these <- critC$category_C_code %in% "C1" & critC$C2 %in% "LC or NT" 
critC$category_C[replace_these] <- "LC or NT"
critC$category_C_code[replace_these] <- ""

replace_these <- critC$category_C_code %in% "C1" & !critC$C2 %in% "LC or NT" & as.double(critC$assess.pop.size) < 1000
critC$category_C[replace_these] <- critC$C2[replace_these]
critC$category_C_code[replace_these] <- "C2ai"

replace_these <- critC$category_C_code %in% "C1" & !critC$C2 %in% "LC or NT" & as.double(critC$assess.pop.size) > 1000
critC$category_C[replace_these] <- "LC or NT"
critC$category_C_code[replace_these] <- ""
critC$category_C_code <- gsub("C1\\+", "", critC$category_C_code)
critC$assess.pop.size <- NULL

#Criterion D
critD <- readRDS("data/criterionD_all_prop_mature.rds")
critD <- critD[,c("species","pop.size","pop.size.low","D","D.AOO","D2.Loc",
                  "category_D","category_D_code",
                  "D.p0.18","D.p0.25","D.p0.33","D.p0.45", "D.p0.51")]
critD$pop.size <- round(as.double(critD$pop.size),1)
critD$pop.size.low <- round(as.double(critD$pop.size.low),1)

#Estimated pop. sizes based on AOO, taxonomy, life-form and endemism level (not used in the assessments)
est.pop <- readRDS("data/estimated_pop_size_from_AOO.rds")
critD <- merge(critD, est.pop[,c("species","pred","pred.low")], by = "species", all = TRUE, sort = FALSE)
critD$pop.size[is.na(critD$pop.size)] <- 
  round(as.double(critD$pred[is.na(critD$pop.size)]),1)
critD$pop.size.low[is.na(critD$pop.size.low)] <- 
  round(as.double(critD$pred.low[is.na(critD$pop.size.low)]),1)
critD <- critD[order(critD$species),]
critD <- critD[,c("species","pop.size","pop.size.low","D","D.AOO","D2.Loc",
                  "category_D","category_D_code",
                  "D.p0.18","D.p0.25","D.p0.33","D.p0.45", "D.p0.51")]

#Desregarding assessments based on criteria D2 only (no species-specific info on threats)
critD$category_D_code[!critD$category_D %in% "LC or NT" & 
                    critD$D %in% "LC or NT"] <- "" 
critD$category_D[!critD$category_D %in% "LC or NT" & 
                    critD$D %in% "LC or NT"] <- "LC or NT" 

## Merging all assessments
all.crit <- merge(critA, critB, by = "species", all = TRUE)
all.crit <- merge(all.crit, critC, by = "species", all = TRUE)
all.crit <- merge(all.crit, critD, by = "species", all = TRUE)
all.crit <- all.crit[order(all.crit$species),]
rm(critA, critB, critB_high, critB_low, critB_opt, critC, critD, est.pop)

###############################################################################################################################################################################################H
#########################
#### NEAR THREATENED ####
#########################
subcriteria <- c("A2", "B1", "B2", "C2", "D")
for(i in 1:length(subcriteria)) {
  all.crit[,subcriteria[i]] <- 
    ConR::near.threatened(all.crit[,subcriteria[i]],
                          all.crit$EOO, all.crit$AOO, all.crit$declineB,
                          all.crit$reduction_A12, all.crit$pop.size, all.crit$pop.size.low,
                          all.crit$nbe_loc_total, all.crit$sever.frag, ext.fluct = NULL,
                          subpop = all.crit$Nbe_subPop, subcriteria = subcriteria[i])
}
rm(subcriteria)

###############################################################################################################################################################################################H
##############################
#### CONSENSUS ASSESSMENT ####
##############################

assess.df <- all.crit[,c("species", "A2", "B1", "B2", "C2", "D")]
tmp <- ConR:::cat_mult_criteria(assess.df)
table(tmp$species == all.crit$species)
all.crit <- cbind.data.frame(all.crit, 
                             tmp[,c("category","main.criteria","aux.criteria")], stringsAsFactors = FALSE)
rm(assess.df)

###############################################################################################################################################################################################H
#########################
#### DATA DEFICIENT  ####
#########################

#### LOADING THREAT OCCURRENCE DATA (HERBARIUM + TREECO) ###
oc.data <- readRDS("data/threat_occ_data_final.rds")

#Organizing the data
MyData <- oc.data[, c("higher.tax.rank","tax",
                      "ddlat","ddlon",
                      "coly","vouchers",
                      "detBy","dety","typeStatus",
                      "tax.check2","tax.check.final",
                      #"UC","dist.eoo",
                      "tax.conf","source")]
rm(oc.data)

##Taxonomic confidence of the optimum scheme
tmp0 <- MyData[tax.check2 %in% TRUE & tax.check.final %in% c("high"), uniqueN(vouchers), by = tax]
tmp0.1 <- MyData[tax.check.final %in% c("high"), uniqueN(vouchers), by = tax]
tmp0 <- merge(tmp0.1, tmp0, by = "tax", all = TRUE)
tmp0$V1.y[is.na(tmp0$V1.y)] <- 0
tmp0$prop.high.tax <- tmp0$V1.y / tmp0$V1.x 
quantile(tmp0$prop.high.tax, prob = seq(0,1,0.05))
boxplot(tmp0$prop.high.tax, notch = TRUE)
plot(tmp0$prop.high.tax ~ log10(tmp0$V1.x))
names(tmp0)[1:3] <- c("species", "n.occs","n.occs.high")

## Getting the year of last colection date
MyData[coly %in% "2997", coly := 1997] 
MyData[coly %in% "2992", coly := 1992]
MyData[coly %in% "2188", coly := 1988]
MyData[coly %in% "2097", coly := 1997]
MyData[coly %in% "1012", coly := 2012] 
MyData[coly %in% "1394", coly := 1934]
MyData[coly %in% "1089", coly := 1989]
MyData[coly %in% "1076", coly := 1976]
MyData[coly %in% "1014", coly := 2014]
MyData[coly %in% "1190", coly := 1990]
MyData[coly %in% "1235", coly := 2012]
MyData[coly %in% "1197", coly := 1997]
MyData[coly %in% "1192", coly := 1997]
MyData[coly %in% "1071", coly := 1971]

tmp <- aggregate(MyData$coly, list(MyData$tax), max, na.rm = TRUE)
tmp$x[is.infinite(tmp$x)] <- NA
tmp.1 <- aggregate(MyData$vouchers[MyData$tax.check.final %in% c("high")], list(MyData$tax[MyData$tax.check.final %in% c("high")]), length)
tmp.2 <- aggregate(MyData$coly[MyData$tax.check.final %in% c("high")], list(MyData$tax[MyData$tax.check.final %in% c("high")]), function(x) (sum(is.na(x))))
tmp <- merge(tmp, tmp.1, by = "Group.1", all = TRUE)
tmp <- merge(tmp, tmp.2, by = "Group.1", all = TRUE)
names(tmp) <- c("species", "last.record","n.occs","n.occs.coly.na")
tmp$prop.coly.na <- tmp$n.occs.coly.na/tmp$n.occs

#tmp[!is.na(tmp$last.record) & !is.na(tmp$prop.coly.na) & tmp$last.record < 1970 & tmp$prop.coly.na <20, ]
table(tmp$last.record < 1970) # species without new in situ records in the past 50 years

##Known only from type locality
MyData$typeStatus[MyData$typeStatus %in% c("not a type","NOTATYPE","Não idetificado","Mendonça, J.O.","Barbosa, A.V.G.","Campos, C.J.","1")] <- NA
MyData$typeStatus[grepl("possible|probable", MyData$typeStatus, ignore.case = TRUE)] <- NA

tmp1 <- aggregate(MyData$typeStatus[MyData$tax.check.final %in% c("high")], list(MyData$tax[MyData$tax.check.final %in% c("high")]), function(x) (sum(!is.na(x))))
tmp1$n.occs <- aggregate(MyData$typeStatus[MyData$tax.check.final %in% c("high")], list(MyData$tax[MyData$tax.check.final %in% c("high")]), length)$x
tmp1$only.type <- FALSE
tmp1$only.type[tmp1$x == tmp1$n.occs] <- TRUE
names(tmp1)[1:2] <- c("species", "n.types")

#Merging all info
tmp <- merge(tmp[,c("species",  "n.occs", "last.record", "prop.coly.na")], 
             tmp0[,c("species",  "prop.high.tax")], by = "species", all = TRUE, sort = FALSE)
tmp <- merge(tmp, tmp1[,c("species",  "only.type")], by = "species", all = TRUE, sort = FALSE)
tmp <- tmp[order(tmp$species),]

#Binding with the rest of the information
table(tmp$species == all.crit$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp[,c("last.record","prop.coly.na","prop.high.tax","only.type")], stringsAsFactors = FALSE)

### DATA DEFICIENCY ###
## No taxonomic valid occurrence (unknown provenance)
all.crit$category[all.crit$category %in% c("", " ", NA)] <- "DD"

rm(tmp, tmp.1, tmp.2, tmp0, tmp0.1, tmp1)
###############################################################################################################################################################################################H
###################################################
#### Critically Endangered (Possibly Extinct)  ####
###################################################

## We also classified as DD species known only from its (old) type locality
all.crit$category[!is.na(all.crit$only.type) & !is.na(all.crit$last.record) &
           all.crit$only.type & all.crit$last.record < (2018-50) & all.crit$category %in% "CR"] <- "CR_PE"
rm(MyData); gc()


###############################################################################################################################################################################################H
#########################################
#### GETTING SPECIES ENDEMISM LEVELS ####
#########################################

## Getting endemism levels from Lima et al. (2020)
endemism <- readRDS("data/sis_connect/endemism_threat.rds") # endemism levels from Lima et al. (2020)
all.crit <- merge(all.crit, endemism, by = "species", all.x = TRUE)
#Getting missing endemism status
# SA.spp <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
SA.spp <- read.csv("C://Users//renato//Desktop//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
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
all.crit$endemic[all.crit$endemic %in% "not in the AF"] <- "occasional"
all.crit$endemic[all.crit$endemic %in% "not_endemic"] <- "widespread_sparse"
table(all.crit$endemic, useNA = "always")

#Saving the domain information
tmp1 <- merge(all.crit[,c("species"), drop = FALSE] , 
              SA.spp[,c("SpeciesKTSA","domain.reflora",
                        "Domain.Ary","Domain","geographical_distribution","Grand.Total","vegtype.reflora")],
              by.x="species", by.y= "SpeciesKTSA", all.x = TRUE)
tmp1$Domain[is.na(tmp1$Domain)] <- 
  tmp1$domain.reflora[is.na(tmp1$Domain)]
tmp1$Domain[tmp1$Domain %in% "no_records"] <- 
  tmp1$domain.reflora[tmp1$Domain %in% "no_records"]
tmp1$Domain[is.na(tmp1$Domain)] <- 
  tmp1$Domain.Ary[is.na(tmp1$Domain)]
table(tmp1$species == all.crit$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp1[,c("Domain","geographical_distribution","vegtype.reflora","Grand.Total"), drop = FALSE], stringsAsFactors = FALSE)
names(all.crit)[dim(all.crit)[2]] <- "Amazonian.Occurrences"
rm(SA.spp, endemism, tmp, tmp1)

##Fixing some final problems related to endemism (Pantropical or Neotropical species)
all.crit$endemic[all.crit$species %in% "Alsophila capensis"] <- "widespread_sparse"
all.crit$endemic[all.crit$species %in% "Lycianthes rantonnetii"] <- "widespread_sparse"
all.crit$endemic[all.crit$endemic %in% "endemic" & all.crit$geographical_distribution %in% c("neotropical") &
           all.crit$endemism.level.1 <85 & all.crit$endemism.level.2 <85] <- "widespread_common"
all.crit$endemic[all.crit$endemic %in% "endemic" & all.crit$geographical_distribution %in% c("south_american") &
                   all.crit$endemism.level.1 <85 & all.crit$endemism.level.2 <85] <- "widespread_common"
all.crit$endemic[all.crit$species %in% 
                   c("Schinus terebinthifolia","Acnistus arborescens","Cupania oblongifolia",
                     "Nectandra oppositifolia","Pereskia aculeata","Piper mollicomum",
                     "Quararibea turbinata","Styrax glabratus","Acca sellowiana")] <- "widespread_common"

### THE SPECIES IS ENDEMIC TO BRAZIL?
reflora <- readRDS("data/threat_fbo_tax_info.rds") # information from ReFlora
tmp <- merge(all.crit, reflora[,c("species.correct2","endemism","phytogeographicDomain")], 
             by.x = "species", by.y = "species.correct2", all.x = TRUE)
table(tmp$species == all.crit$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp[,c("endemism","phytogeographicDomain"), drop = FALSE], stringsAsFactors = FALSE)
names(all.crit)[dim(all.crit)[2]-1] <- "endemismo.reflora"
names(all.crit)[dim(all.crit)[2]] <- "dominio.reflora"


###############################################################################################################################################################################################H
##########################################
#### GETTING SPECIES EOO DESCRIPTIONS ####
##########################################

## % of the EOO inside the AF
tmp <- readRDS("data/EOO_AtlanticForest_cropped.rds")
names(tmp)[3] <- "prop.EOO.in.AF"
tmp1 <- merge(all.crit, tmp, by.x = "species", by.y = "tax", all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species),]
table(all.crit$species == tmp1$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp1[,c("prop.EOO.in.AF"), drop = FALSE], stringsAsFactors = FALSE)
plot(all.crit$endemism.level.1 ~ all.crit$prop.EOO.in.AF)

## % of the EOO inside protected areas
tmp <- readRDS("data/EOO_StrictUCs_cropped.rds")
names(tmp)[3] <- "prop.EOO.in.StrictUCs"
tmp1 <- merge(all.crit, tmp, by.x = "species", by.y = "tax", all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species),]
table(all.crit$species == tmp1$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp1[,c("prop.EOO.in.StrictUCs"), drop = FALSE], stringsAsFactors = FALSE)
plot(all.crit$protected ~ all.crit$prop.EOO.in.StrictUCs)

## mean value of the Human Influence Index
tmp <- readRDS("data/EOO_HII_cropped.rds")
names(tmp)[4] <- "Median.HII"
tmp1 <- merge(all.crit, tmp, by.x = "species", by.y = "tax", all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species),]
table(all.crit$species == tmp1$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp1[,c("Median.HII"), drop = FALSE], stringsAsFactors = FALSE)

## Forest cover
tmp <- readRDS("data/EOO_hab_loss_2000_2015_cropped.rds")
names(tmp)[2] <- "prop.EOO.forest"
#names(tmp)[10] <- "habitat.quality" # -1 a 2
tmp1 <- merge(all.crit, tmp, by.x = "species", by.y = "tax", all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species),]
table(all.crit$species == tmp1$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp1[,c("prop.EOO.forest"), drop = FALSE], stringsAsFactors = FALSE)

#### CITES ####
cites <- readRDS("data/threat_cites_EU.rds")
cites <- cites[!is.na(cites$cites_listing),]
cites <- cites[!grepl("NC",cites$cites_listing),]

# cites <- read.csv("Index_of_CITES_Species_[CUSTOM]_2020-09-28 15_34.csv", as.is = TRUE, na.strings = c(NA, "", " "))
# cites.syns <- cites[!is.na(cites$SynonymsWithAuthors),]
# tmp <- strsplit(cites.syns$SynonymsWithAuthors,",")
# tmp <- sapply(tmp, stringr::str_trim)
# names(tmp) <- cites.syns$FullName
# tmp1 <- unlist(tmp)
# names(tmp1) <- gsub('[0-9]$','',names(tmp1))
# tmp2 <- sapply(tmp1, flora::remove.authors)
# tmp3 <- data.frame(FullName = names(tmp1),
#                    Synonym = tmp2,
#                    stringsAsFactors = FALSE)
# tmp4 <- merge(tmp3, cites.syns, by = "FullName", all.x = TRUE)
# tmp4 <- tmp4[,c("TaxonId","Kingdom","Family","Synonym","CurrentListing","AnnotationEnglish","X.AnnotationSymbol")]
# names(tmp4)[4] <- "FullName" 
# #tmp4 <- tmp4[!duplicated(tmp4$FullName),]
# tmp5 <- rbind.data.frame(cites[,c("TaxonId","Kingdom","Family","FullName","CurrentListing","AnnotationEnglish","X.AnnotationSymbol")],
#                          tmp4, stringsAsFactors = FALSE)
tmp6 <- merge(all.crit, 
              cites, by.x = "species", by.y = "species.correct2", all.x = TRUE, sort= FALSE)
tmp6 <- tmp6[!duplicated(tmp6$species),]
tmp6 <- tmp6[order(tmp6$species),]
table(tmp6$species == all.crit$species)
all.crit <- cbind.data.frame(all.crit,
                             tmp6[,c("cites_listing", "annotation", "code_annotation"), drop = FALSE], stringsAsFactors = FALSE)
names(all.crit)[(dim(all.crit)[2]-2):dim(all.crit)[2]] <- 
  c("CITES.CurrentListing","CITES.AnnotationEnglish","CITES.AnnotationCode")
rm(cites, cites.syns, tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)


###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
#####################################################
#### ADDING GL AND UNKNOWN ECOLOGICAL GROUP INFO ####
####################################################

hab <- read.csv("data/threat_habitats.csv")
hab$ecol.group[hab$Name_submitted %in% "Ximenia americana"] <- "early_secondary"
hab$ecol.group[hab$Name_submitted %in% "Tococa guianensis"] <- "early_secondary"
tmp <- hab[,c("GenerationLength.range", "ecol.group", 
              "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName",
              "Name_submitted"), drop = FALSE]
names(tmp) <- c("gen.length", "ecol.group", "growth.form", "Name_submitted")
table(all.crit$species == tmp$Name_submitted)
all.crit <- cbind.data.frame(all.crit,
                             tmp[,c(1, 2, 3)], stringsAsFactors = FALSE)



###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
##############################
#### PREVIOUS ASSESSMENTS ####
##############################

## Loading previous assessments from IUCN, Argentina and Paraguay
prev.assess <- readRDS("IUCN_2022_v2_assessments_THREAT.rds")
# prev.assess <- readRDS("data/sis_connect/prev_assessments_threat.rds")
prev.assess.PY <- read.csv("especies_ameacadas_Paraguay.csv", as.is = TRUE, na.strings = c("", " ", NA))
prev.assess.AR <- read.csv("especies_ameacadas_Argentina.csv", as.is = TRUE, na.strings = c("", " ", NA))
end.AR <- read.csv("especies_endemicas_Argentina.csv", as.is = TRUE, na.strings = c("", " ", NA))

## Merging IUCN with THREAT assessments
prev.spp.name <- "scientificName"
prev.cols <- c("redlistCategory", "redlistCriteria", "yearPublished",
               "assessmentDate", "criteriaVersion","language", "rationale", "habitat",
               "threats","population","populationTrend","range","useTrade",
               "conservationActions", "systems", "realm", 
               "yearLastSeen", "possiblyExtinct", "possiblyExtinctInTheWild",
               "scopes",
               "AOO.range", "EOO.range",
               "Congregatory.value", "NoThreats.noThreats", 
               "ElevationLower.limit", "ElevationUpper.limit", 
               "PopulationSize.range", 
               "ThreatsUnknown.value",
               "LocationsNumber.range", "GenerationLength.range", 
               "SubpopulationNumber.range", 
               "AreaRestricted.isRestricted",
               "CropWildRelative.isRelative",
               "YearOfPopulationEstimate.value", 
               "InPlaceEducationControlled.value", 
               "SevereFragmentation.isFragmented", 
               "InPlaceResearchRecoveryPlan.value", 
               "InPlaceLandWaterProtectionInPA.value", 
               "InPlaceSpeciesManagementExSitu.value", 
               "InPlaceResearchMonitoringScheme.value", 
               "InPlaceSpeciesManagementHarvestPlan.value", 
               "InPlaceLandWaterProtectionAreaPlanned.value",     
               "InPlaceEducationInternationalLegislation.value",
               "InPlaceLandWaterProtectionInvasiveControl.value",
               "InPlaceLandWaterProtectionSitesIdentified.value",
               "InPlaceLandWaterProtectionPercentProtected.value")

tmp <- merge(all.crit[,c("species","category")],
             # prev.assess[,c("species.correct2","status.reflora","redlistCategory","redlistCriteria","yearPublished","population","populationTrend","threats","rationale")],
             # by.x = "species", by.y="species.correct2", all.x = TRUE, sort = FALSE)
              prev.assess[, c(prev.spp.name, prev.cols)],
              by.x = "species", by.y = prev.spp.name, all.x = TRUE, sort = FALSE)
tmp$redlistCategory[tmp$redlistCategory %in% "Extinct"] <- "EX"
tmp$redlistCategory[tmp$redlistCategory %in% "Extinct in the Wild"] <- "EW"
tmp$redlistCategory[tmp$redlistCategory %in% "Critically Endangered"] <- "CR"
tmp$redlistCategory[tmp$redlistCategory %in% "Endangered"] <- "EN"
tmp$redlistCategory[tmp$redlistCategory %in% "Vulnerable"] <- "VU"
tmp$redlistCategory[tmp$redlistCategory %in% c("Near Threatened","Lower Risk/near threatened","Lower Risk/conservation dependent")] <- "NT"
tmp$redlistCategory[tmp$redlistCategory %in% "Data Deficient"] <- "DD"
tmp$redlistCategory[tmp$redlistCategory %in% c("Least Concern","Lower Risk/least concern")] <- "LC"
tmp <- tmp[order(tmp$species),]
table(all.crit$species == tmp$species)
all.crit <- cbind.data.frame(all.crit, 
                             # tmp[,c("status.reflora","redlistCategory","redlistCriteria","yearPublished","population","populationTrend","threats","rationale")], 
                             tmp[, prev.cols], 
                             stringsAsFactors = FALSE)

## Getting CNCFlora previous assessments with THREAT assessments
tmp <- flora::get.taxa(all.crit$species, replace.synonyms = FALSE,
                       drop = c())
table(all.crit$species == tmp$original.search)
all.crit$status.reflora <- tmp$threat.status

## Merging Argentina and Paraguay with THREAT assessments
prev.assess.AR$species <- stringr::str_trim(prev.assess.AR$species)
prev.assess.AR$categoria <- stringr::str_trim(prev.assess.AR$categoria)
prev.assess.PY$species <- stringr::str_trim(prev.assess.PY$species)
prev.assess.PY$categoria <- stringr::str_trim(prev.assess.PY$categoria)
tmp1 <- merge(prev.assess.AR[,c("species","categoria")], 
              prev.assess.PY[,c("species","categoria")],
              by = "species", all = TRUE)
names(tmp1)[2:3] <- c("category.ARG","category.PAY")
tmp1 <- tmp1[!duplicated(tmp1$species),]
tmp1 <- merge(tmp1, end.AR[,c("species","categoria","obs")],
              by = "species", all = TRUE)
names(tmp1)[4:5] <- c("endemism.ARG","obs.endemism.ARG")

#Replacing some categories
tmp1$category.ARG[tmp1$category.ARG %in% "Endangered"] <- "EN"
tmp1$category.ARG[tmp1$category.ARG %in% "Vulnerable"] <- "VU"
tmp1$category.PAY[tmp1$category.PAY %in% "en_peligro"] <- "EN"
tmp1$category.PAY[tmp1$category.PAY %in% "amenazada"] <- "VU"

#Solving some typos and synonymizations
tmp1.1 <- flora::get.taxa(stringr::str_trim(tmp1$species), drop=c("scientific.name","authorship","threat.status"))
tmp1.1$species <- paste(tmp1.1$genus, tmp1.1$specific.epiteth, sep = " ")
tmp1.1$accepted.name[!is.na(tmp1.1$accepted.name)] <- 
  as.character(sapply(tmp1.1$accepted.name[!is.na(tmp1.1$accepted.name)], flora::remove.authors))
tmp1.1$notes[tmp1.1$species %in% "Solanum robustum"] <- ""
tmp1$species[tmp1.1$notes %in% c("was misspelled","was misspelled|replaced synonym")] <- 
  tmp1.1$species[tmp1.1$notes %in% c("was misspelled","was misspelled|replaced synonym")]
tmp1$species[tmp1.1$notes %in% c("check no accepted name")] <- 
  tmp1.1$accepted.name[tmp1.1$notes %in% c("check no accepted name")]
tmp1$species[tmp1.1$notes %in% c("replaced synonym")] <- 
  tmp1.1$species[tmp1.1$notes %in% c("replaced synonym")]
tmp1$category.ARG[!is.na(tmp1$category.ARG) & nchar(tmp1$category.ARG) == 2 & tmp1.1$notes %in% c("replaced synonym")] <-
  paste(tmp1$category.ARG[!is.na(tmp1$category.ARG) & nchar(tmp1$category.ARG) == 2 & tmp1.1$notes %in% c("replaced synonym")],
        " (under ", 
        tmp1.1$original.search[!is.na(tmp1$category.ARG) & nchar(tmp1$category.ARG) == 2 & tmp1.1$notes %in% c("replaced synonym")],
        ")",sep="")
tmp1$category.PAY[!is.na(tmp1$category.PAY) & nchar(tmp1$category.PAY) == 2 & tmp1.1$notes %in% c("replaced synonym")] <-
  paste(tmp1$category.PAY[!is.na(tmp1$category.PAY) & nchar(tmp1$category.PAY) == 2 & tmp1.1$notes %in% c("replaced synonym")],
        " (under ", 
        tmp1.1$original.search[!is.na(tmp1$category.PAY) & nchar(tmp1$category.PAY) == 2 & tmp1.1$notes %in% c("replaced synonym")],
        ")",sep="")
spp <- tmp1$species[duplicated(tmp1$species)]
tmp1$category.PAY[tmp1$species %in% spp[1]] <- "VU (under Baccharis cylindrica and Baccharis trimera)"
tmp1$category.ARG[tmp1$species %in% spp[2]] <- "VU"
tmp1$category.ARG[tmp1$species %in% spp[3]] <- "EN (under Cyathea vestita)"
tmp1$category.ARG[tmp1$species %in% spp[4]] <- "Indeterminada (maybe LC)"
tmp1 <- tmp1[!duplicated(tmp1$species),]

#Some final edits
tmp1$category.ARG[tmp1$category.ARG %in% "EN (under Rheedia brasiliensis)"] <- "EN"
tmp1$category.ARG[tmp1$category.ARG %in% "EN (under Tabebuia caraiba)"] <- "EN"
tmp1$category.ARG[tmp1$category.ARG %in% "EN (under Trichipteris atrovirens)"] <- "EN"
tmp1$category.ARG[tmp1$category.ARG %in% "LC (under Patagonula americana)"] <- "LC"

tmp1$category.PAY[tmp1$category.PAY %in% "EN (under Astronium balansae)"] <- "EN"
tmp1$category.PAY[tmp1$category.PAY %in% "EN (under Tabebuia heptaphylla)"] <- "EN"
tmp1$category.PAY[tmp1$category.PAY %in% "EN (under Tabebuia pulcherrima)"] <- "EN"
tmp1$category.PAY[tmp1$category.PAY %in% "EN (under Trichipteris atrovirens)"] <- "EN"
tmp1$category.PAY[tmp1$category.PAY %in% "VU (under Tabebuia alba)"] <- "VU"

#Crossing with THREAT
tmp2 <- merge(all.crit[,c("species","category")], tmp1, 
             by = "species", all.x = TRUE, sort = FALSE)
tmp2 <- tmp2[order(tmp2$species),]
table(all.crit$species == tmp2$species)
all.crit <- cbind.data.frame(all.crit, 
                             tmp2[,c("category.ARG","category.PAY","endemism.ARG")], 
                             stringsAsFactors = FALSE)
rm(prev.assess, prev.assess.AR, prev.assess.PY, end.AR, spp, tmp, tmp1, tmp1.1, tmp2)

###############################################################################################################################################################################################H
###############################################################################################################################################################################################H
##############################
#### REGIONAL ASSESSMENTS ####
##############################
### FROM IUCN (2012) ###
all.crit$category.regional <- all.crit$category

## 1- Vagrant taxa (i.e found only occasionally within the boundaries of a region) should NOT be assessed ##
#I include an execption to criteria B which uses all occurrence data available (inside and outside) 
all.crit$category.regional[all.crit$endemic %in% "occasional"] <- "NA"
all.crit$category.regional[all.crit$endemic %in% "occasional" &
                             !all.crit$category_B %in% c("LC or NT")] <-
  all.crit$category_B[all.crit$endemic %in% "occasional" &
                               !all.crit$category_B %in% c("LC or NT")]  
all.crit$main.criteria[all.crit$endemic %in% "occasional" &
                             !all.crit$category_B %in% c("LC or NT")] <-
  all.crit$category_B_code[all.crit$endemic %in% "occasional" &
                        !all.crit$category_B %in% c("LC or NT")]
all.crit$main.criteria <- gsub("ab", "", all.crit$main.criteria)

## 2 - Downlistings ##
##If there are no conspecific populations in neighbouring regions or if propagules are unable
#to disperse to the region, the regional population behaves as an endemic and the category should be left unchanged.
## If there is not enough suitable habitat and if current conservation measures are not leading
#to an improvement in the quality and/or quantity of habitat within the foreseeable future, (...)
#immigration from outside the region will not decrease extinction risk and the category should be left unchanged.
##Are conditions outside or within the region deteriorating? If yes the taxon will decrease and thus
#category should remain unchanged

###EXAMPLES OF DOWNLISTING CASES
##Based on the most probable number of mature individuals, the species meets the
#criterion for CR D. Since there is an obvious probability of recolonisation from neighbouring
#countries, the category is downlisted to EN° D.

##Bird in Sweden: Depending on the values of population size span are used, the criteria 
#meet is VU D1 or EN D, of which EN D is the most plausible. Since immigration from neighbouring 
#countries is possible, the risk of extinction is probably less than if the subpopulation were isolated.
#For instance, the Norwegian subpopulation is stable at 1000-3000 pairs. Accordingly, the
#category is downlisted to VU° D1.

## Taxon was listed preliminarily under CR A2acd; C2a(ii). Conditions are deteriorating within 
#Viet Nam but there is uncertainty about conditions outside the region (e.g. in Cambodia); 
#the global population is also in decline. The preliminary assessment category is therefore unchanged.

##The North African populations meet the criteria for Endangered under B2.
#However, the species is classified as Least Concern in Europe and the Mediterranean
#and it is readily transported (e.g. by ducks), so a rescue effect from European populations
#is expected. Therefore, the assessment is downlisted to Vulnerable (VU° B2ab(iii)).

## Decision: Dowlisting was performed only in some cases: 
#Case 1: for Cerrado, Caatinga, Pantanal and Pampa 'specialists':
tmp <- all.crit$category.regional[grepl("widespread", all.crit$endemic) &
                             !all.crit$category.regional %in% c("NT","LC","LC or NT") &
                             !all.crit$geographical_distribution %in% c("local_endemic","regional_endemic") &
                             all.crit$Domain %in% c(
                               "Caatinga","Caatinga|Cerrado","Caatinga|Cerrado|Chaco","Caatinga|Pantanal",
                               "Cerrado","Cerrado|Pantanal","Pantanal","Pampa"
                             )]
tmp1 <- ConR:::cat_downlist(tmp)
#tmp1 <- cat_downlist(tmp)

all.crit$category.regional[grepl("widespread", all.crit$endemic) &
                             !all.crit$category.regional %in% c("NT","LC","LC or NT") &
                             !all.crit$geographical_distribution %in% c("local_endemic","regional_endemic") &
                             all.crit$Domain %in% c(
                               "Caatinga","Caatinga|Cerrado","Caatinga|Cerrado|Chaco","Caatinga|Pantanal",
                               "Cerrado","Cerrado|Pantanal","Pantanal","Pampa"
                             )] <- tmp1

#Case 1.1: for specialists of exclusive of Cerrado, Caatinga, Pantanal and Pampa vegetation types:
veg.types <- 
  unique(all.crit$vegtype.reflora[!grepl("Floresta|Restinga|Altitude|rupestre|Rochosos|Manguezal|Amaz|Campinarana",all.crit$vegtype.reflora)])[-1]
  
tmp <- all.crit$category.regional[grepl("widespread", all.crit$endemic) &
                                    !all.crit$category.regional %in% c("NT","LC","LC or NT","VUo","ENo","NTo") &
                                    !all.crit$geographical_distribution %in% c("local_endemic","regional_endemic") &
                                    all.crit$vegtype.reflora %in% veg.types]
tmp1 <- ConR:::cat_downlist(tmp)
#tmp1 <- cat_downlist(tmp)

all.crit$category.regional[grepl("widespread", all.crit$endemic) &
                             !all.crit$category.regional %in% c("NT","LC","LC or NT","VUo","ENo","NTo") &
                             !all.crit$geographical_distribution %in% c("local_endemic","regional_endemic") &
                             all.crit$vegtype.reflora %in% veg.types] <- tmp1

#Case 2: widespread, multi-domain species (high dispersal ability?)
## DOING NOTHING FOT THE MOMENT

#Case 3: Taxa listed at a higher category at the global level than for the Atlantic Forest 
## DOING NOTHING FOT THE MOMENT

#Clean column for data analyses
all.crit$cat.reg.clean <- all.crit$category.regional
all.crit$cat.reg.clean <- gsub("_PE|o$", "", all.crit$cat.reg.clean)
all.crit$status.reflora <- gsub("DD\\|DD", "DD", all.crit$status.reflora)

#### SAVING THE FINAL RESULTS TABLE ####
all.crit <- all.crit[order(all.crit$species),]
saveRDS(all.crit, "data/all.criteria.rds")
