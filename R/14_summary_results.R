#################################
#### SUMMARIZING THE RESULTS ####
#################################
rm(list = ls())

##Loading the data
all.crit <- readRDS("data/all.criteria.rds")

##Loading packages
require(webr)
require(ggplot2)
library(cowplot)
require(multcompView)
# library(ggarrange)
# require(ggforce)
# require(grid)
# require(moonBook)

##Species for CNCFlora
# tmp <- c("Kuhlmanniodendron macrocarpum", "Xylopia langsdorffiana","Protium atlanticum",
#          "Licania belemii","Paubrasilia echinata","Calyptranthes ubatubana",
#          "Myrcia brunea")
# tmp1 <- c("family","species","assessment.period","reduction_A12","A1","A2","category_A","category_A_code",
#           "Nbe_occs","EOO","EOO.high.tax","AOO","AOO.high.tax","Nbe_subPop","nbe_loc_total","declineB","sever.frag",            
#           "protected","B1","B1.high.tax","B2","B2.high.tax","category_B","category_B_code",
#           "category_B_code.high.tax","category_B.high.tax",
#           "any.decline","cont.decline","C1","C2","category_C","category_C_code",
#           "pop.size","pop.size.low","D","D.AOO","D2.Loc","category_D","category_D_code",
#           "category","main.criteria","aux.criteria","endemic",
#           "last.record","only.type","prop.EOO.in.AF","prop.EOO.in.StrictUCs", "Median.HII","prop.EOO.forest",
#           "status.reflora", "redlistCategory","redlistCriteria","yearPublished")
# (tmp2 <- all.crit[all.crit$species %in% tmp, tmp1])
# write.csv(tmp2, "threat.assesments.to.cncflora.csv")

#### REMOVE ANY SPECIES BEFORE FINAL ANALYSES: UNRESOLVED NAMES, SMALL SHRUBS? ####


######################################################
#### THE CONSERVATION STATUS OF THE AF TREE FLORA ####
######################################################

## HOW MANY SPECIES ASSESSED?
dim(all.crit)[1] # 5094 before, now 4953

## HOW MANY HERBARIUM RECORDS USED?
oc.data <- readRDS("data/threat_occ_data_final.rds")
table(oc.data$source) # 816,192 herbarium records + 91463 additional records from TreeCo

## HOW MANY TREE COUNTS?
oc.data <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//species_vs_sites_table.csv")
dim(oc.data)[1]
sum(oc.data[,-1], na.rm = TRUE) # 1,348,987 trees on 1133 inventories
oc.data <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//plot_metadata.csv")
sum(oc.data$effort_ha, na.rm = TRUE)  # 1154.36 ha
rm(oc.data)

## HOW MANY THREATENED SPECIES?
## Overall
sum(table(all.crit$cat.reg.clean)[c("CR", "EN", "VU")]) # 3408 threatened; now 3303
round(100 * sum(table(all.crit$cat.reg.clean)[c("CR", "EN", "VU")]) / dim(all.crit)[1], 1) # 66.9 threatened; now 66.7% 
## Endemic
round(100 * sum(table(all.crit$endemic)["endemic"]) / dim(all.crit)[1], 1) # 50.8 are pure and near endemics; 49.7%
round(100 * sum(table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])[c("CR", "EN", "VU")]) / 
        dim(all.crit[all.crit$endemic %in% "endemic",])[1], 1) # 84.4% of the endemics are threatened; now 84.5%
sum(table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])[c("CR", "EN", "VU")]) #2186 species; now 2083

## THE RLI FOR THE ATLANTIC FOREST TREE FLORA 
red::rli(all.crit$cat.reg.clean, boot = TRUE, runs = 50000)
red::rli(all.crit$cat.reg.clean[!all.crit$cat.reg.clean %in% "NA"], boot = TRUE, runs = 50000)

## MORE FREQUENT CATEGORIES AND CRITERIA
round(100 * table(all.crit$cat.reg.clean)[c("CR", "EN", "VU")] / dim(all.crit)[1], 1) # 43.1% of EN, now 42%
round(100 * table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])[c("CR", "EN", "VU")] / 
        dim(all.crit[all.crit$endemic %in% "endemic",])[1], 1) # 55.9% of EN; now 55.6%

## How many near threatened
round(100 * sum(table(all.crit$cat.reg.clean)[c("NT")]) / dim(all.crit)[1], 1) #1.4% of NT; 1.3%

## SOME EXAMPLES OF THREATENED SPECIES WITH IMPORTANT USES
cols <- c("species","reduction_A12","reduction_A12.25ys","EOO","AOO","nbe_loc_total","sever.frag","declineB","any.decline","pop.size",
          "category","cat.reg.clean","main.criteria","endemism.level.1","endemism.level.2","endemic","last.record","only.type",
          "status.reflora","redlistCategory","redlistCriteria")
iconic <- c("Araucaria angustifolia","Euterpe edulis","Ilex paraguariensis","Paubrasilia echinata")
all.crit[all.crit$species %in% iconic, cols]
timber <- c("Cariniana legalis","Dalbergia nigra","Melanoxylon brauna",
            "Myrocarpus frondosus","Ocotea odorifera", 
            "Ocotea porosa","Paratecoma peroba")
all.crit[all.crit$species %in% timber, cols]
# food.orn <- c("Eugenia brasiliensis","Handroanthus chrysotrichus","Ilex paraguariensis")
# all.crit[all.crit$species %in% food.orn, cols]


#########################
#### IUCN CATEGORIES ####
#########################
tmp0 <- table(all.crit$main.criteria[!all.crit$cat.reg.clean %in% c("LC","NA","DD","NT")], useNA = "always")
tmp <- sort(round(100 *  tmp0 / 
                    dim(all.crit[!all.crit$cat.reg.clean %in% c("LC","NA","DD","NT"),])[1], 4))
sum(tmp[grepl("A2", names(tmp))]) #74.9% of the threatened species had A2 within its main criteria; now 75.9%
sum(tmp[grepl("B2", names(tmp))]) #28.4% of the threatened species had B2 within its main criteria; now 27.3%
sum(tmp[grepl("C", names(tmp))]) #0.24% of the threatened species had C within its main criteria; now 0.24%
sum(tmp[grepl("D", names(tmp))]) #0.03% of the threatened species had D within its main criteria; now 0.03%
sum(tmp0[grepl("C", names(tmp0))]) # now 8
sum(tmp0[grepl("D", names(tmp0))]) # now 1

## MEAN POPULATION REDUCTIONS
summary(all.crit$reduction_A12); hist(all.crit$reduction_A12, nclass=40)
plot(all.crit$reduction_A12 ~ all.crit$reduction_A12.25ys); abline(0,1); abline(v=0,h=0,lty=2)
confint(lm(all.crit$reduction_A12 ~ 1))
table(all.crit$reduction_A12>=90)
table(all.crit$reduction_A12>=80)
pop.decline.thres <- all.crit$reduction_A12>=30
toto <- table(pop.decline.thres, all.crit$endemic, useNA = "always")
100 * toto[rownames(toto) %in% "TRUE", !colnames(toto) %in% "occasional"]/
  sum(toto[, !colnames(toto) %in% "occasional"]) # 62.6% com declínio acima do cutoff (não ocasionais)
100 * toto[rownames(toto) %in% "FALSE", !colnames(toto) %in% "occasional"]/
  sum(toto[, !colnames(toto) %in% "occasional"]) # 3.56% com declínio abaixo do cutoff (não ocasionais)
100 * toto[rownames(toto) %in% "TRUE", colnames(toto) %in% "endemic"]/
  sum(toto[, colnames(toto) %in% "endemic"]) # 59.38% com declínio acima do cutoff (só endêmicas)
100 * toto[rownames(toto) %in% "FALSE", colnames(toto) %in% "endemic"]/
  sum(toto[, colnames(toto) %in% "endemic"]) # 4.99% com declínio acima do cutoff (só endêmicas)

# AF species with large declines (>80%)
100 * dim(all.crit[!is.na(all.crit$reduction_A12) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$reduction_A12 >= 80,])[1] / dim(all.crit[!is.na(all.crit$reduction_A12) & 
                                                             all.crit$endemic %in% "endemic",])[1] # 9.33%
#Common AF species with very large declines
all.crit[!is.na(all.crit$reduction_A12) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$reduction_A12 >= 90 &
           all.crit$pop.size >= 2e6, cols]
#Species with increasing populations (0.2% of the endemics)
100 * dim(all.crit[!is.na(all.crit$reduction_A12) & 
                     all.crit$endemic %in% "endemic" & 
                     all.crit$reduction_A12<0 &
                     all.crit$reduction_A12.25ys<0,])[1] / dim(all.crit[!is.na(all.crit$reduction_A12) & 
                                                                       all.crit$endemic %in% "endemic",])[1] 
all.crit[!is.na(all.crit$reduction_A12) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$reduction_A12<0 &
           all.crit$reduction_A12.25ys<0, cols]


## SMALL EOO AND AOO
summary(all.crit$AOO); hist(log10(all.crit$AOO), nclass=40); abline(v=log10(c(10,500,2000)), col=c("red","darkorange","gold"))
confint(lm(all.crit$AOO ~ 1))
summary(all.crit$nbe_loc_total); hist(log10(all.crit$nbe_loc_total), nclass=40); abline(v=log10(c(1,5,10)), col=c("red","darkorange","gold"))
confint(lm(all.crit$nbe_loc_total ~ 1))
100*table(all.crit$nbe_loc_total>10)/dim(all.crit)[1]
100*table(all.crit$sever.frag)/dim(all.crit)[1]

plot(log10(all.crit$AOO) ~ log10(all.crit$nbe_loc_total)); abline(h=log10(c(10,500,2000)), lty=2, col=c("red","darkorange","gold")); abline(v=log10(c(1,5,10)),lty=2, col=c("red","darkorange","gold"))
all.crit[!is.na(all.crit$reduction_A12) & !is.na(all.crit$EOO) & !is.na(all.crit$AOO) & 
           all.crit$endemic %in% "endemic" & all.crit$cat.reg.clean %in% "CR" &
           all.crit$reduction_A12>=90 & 
           #all.crit$EOO<30000 &
           all.crit$AOO<100 & 
           all.crit$pop.size<3000, ] 

## EXAMPLES OF HIGHLY THREATENED SPECIES
#Endemics classified as threatened under all criteria
tmp <- all.crit[all.crit$category_A %in% c("CR", "EN", "VU") &
                  all.crit$category_B %in% c("CR", "EN", "VU") & 
                  all.crit$category_C %in% c("CR", "EN", "VU") &
                  all.crit$category_D %in% c("CR", "EN", "VU"), ]
tmp$species[tmp$endemic %in% "endemic"] # there were examples before but not anymore

#Endemics with the largest declines and smaller distributions
all.crit[!is.na(all.crit$reduction_A12) & !is.na(all.crit$EOO) & !is.na(all.crit$AOO) & 
           all.crit$endemic %in% "endemic" & all.crit$cat.reg.clean %in% "CR" &
           all.crit$reduction_A12>85 &
           #all.crit$EOO<30000 &
           all.crit$AOO<100 &
           all.crit$pop.size<3000, cols]


##############################################
#### COMPARISON WITH PREVIOUS ASSESSMENTS ####
##############################################

## OVERALL ASSESSMENTS 
#How many species with IUCN assessments
table(all.crit$redlistCategory)
100*sum(table(all.crit$redlistCategory))/dim(all.crit)[1] # 28.9% with previous IUCN asses.; now 49.9%
table(all.crit$cat.reg.clean, all.crit$redlistCategory)

#How many species with CNCFlora assessments
table(all.crit$status.reflora)
100*sum(table(all.crit$status.reflora))/dim(all.crit)[1]  # 19.7% with previous CNCFlora asses.; remained the same 

#How many species with national assessments
tmp <- all.crit[!is.na(all.crit$status.reflora) | 
                  !is.na(all.crit$category.ARG) |
                  !is.na(all.crit$category.PAY),]
100 * dim(tmp)[1]/dim(all.crit)[1] # 20.3%; now 20.2%

## ENDEMIC SPECIES - proportions not show in the new version of the main text
# #How may species with IUCN assessments
# table(all.crit$redlistCategory[all.crit$endemic %in% "endemic"])
# 100*sum(table(all.crit$redlistCategory[all.crit$endemic %in% "endemic"])) / 
#   dim(all.crit[all.crit$endemic %in% "endemic",])[1] # 20.01% with previous IUCN asses.; now 44.6%
# table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"], 
#       all.crit$redlistCategory[all.crit$endemic %in% "endemic"])
# 
# #How may species with CNCFlora assessments
# table(all.crit$status.reflora[all.crit$endemic %in% "endemic"])
# 100*sum(table(all.crit$status.reflora[all.crit$endemic %in% "endemic"])) / 
#   dim(all.crit[all.crit$endemic %in% "endemic",])[1]  # 24.5% with previous CNCFlora asses.; now 24.6% 


## THEAT NEW ASSESSMENTS
tmp <- all.crit[is.na(all.crit$redlistCategory) &
                  is.na(all.crit$status.reflora) &
                  is.na(all.crit$category.ARG) &
                  is.na(all.crit$category.PAY),]
dim(tmp)[1] # 2959 species (58.09%) without assessments; now 1923 assessments (38.8%)
100 * dim(tmp)[1] / dim(all.crit)[1] # before 58.09%; now 38.8%
dim(tmp[tmp$endemic %in% "endemic",])[1] # before 1632 endemic species; now 1011
100 * dim(tmp[tmp$endemic %in% "endemic",])[1] / dim(tmp)[1] # before 55.12% without assessments; now: 52.6%


##Species EX or EW
ex.spp <- all.crit$species[all.crit$redlistCategory %in% c("EX","EW")]
all.crit[all.crit$species %in% ex.spp, cols]
all.crit[all.crit$species %in% ex.spp, c("species","population")]
oc.data <- readRDS("data/threat_occ_data_final.rds")
oc.data <- oc.data[oc.data$tax %in% ex.spp,]
#Date of last collection
table(oc.data$coly, oc.data$tax)
table(oc.data$coly>=1998, oc.data$tax)
table(oc.data$coly[is.na(oc.data$typeStatus)]>=1998, oc.data$tax[is.na(oc.data$typeStatus)])
oc.data[oc.data$tax %in% "Pouteria stenophylla",]
#Identifications confirmed by specialists
oc.data[oc.data$tax.check2 %in% "TRUE", c("coly","tax")]
oc.data[oc.data$tax %in% c("Pouteria stenophylla"), 
        c("coly","tax","tax.check2","source","detBy","dety","vouchers")]

paths = dir("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//occurrence_data",full.names=TRUE)
paths = paths[grepl('merged_outliers.csv',paths) & grepl('Sapotaceae',paths)]
tmp = data.table::fread(paths,na.strings=c(""," ","NA"))
tmp[tmp$species.correct2 %in% "Pouteria stenophylla",]
rm(tmp,paths,oc.data)
#Data with abundance records
spp.data <- readRDS("data/assess_iucn_spp.rds")
spp.data[spp.data$species %in% ex.spp,]


## PREVIOUS VS. NEW RLI
#overall
red::rli(all.crit$redlistCategory[!is.na(all.crit$redlistCategory)], boot = TRUE, runs = 50000)
red::rli(all.crit$cat.reg.clean[!is.na(all.crit$redlistCategory)], boot = TRUE, runs = 50000)
#endemics: presenting only this one in the main text
red::rli(all.crit$redlistCategory[!is.na(all.crit$redlistCategory) & all.crit$endemic %in% "endemic"], boot = TRUE, runs = 50000) #0.736
red::rli(all.crit$cat.reg.clean[!is.na(all.crit$redlistCategory) & all.crit$endemic %in% "endemic"], boot = TRUE, runs = 50000) #0.478

## COMPARISON OF ALL CRITERIA VS. CRITERIA B
#not in the main text anymore
# tmp <- table(all.crit$main.criteria[all.crit$endemic %in% "endemic"], 
#              all.crit$redlistCriteria[all.crit$endemic %in% "endemic"])
# tmp <- table(all.crit$redlistCriteria[all.crit$endemic %in% "endemic"])
# 100 * sum(tmp[grepl("B", names(tmp))]) / sum(tmp) # prop. of species evaluated under criterion B in previous and new assessments
# 
# tmp1 <- all.crit[grepl("B",all.crit$redlistCriteria) &
#                    all.crit$endemic %in% "endemic",]
# tmp1 <- all.crit[(grepl("B",all.crit$redlistCriteria) | 
#                     all.crit$redlistCategory %in% "LC") & 
#                    all.crit$endemic %in% "endemic",]
# tmp1$redlistCriteria <- gsub("i|v|a|b|c|d|\\,|\\(|\\)", "", tmp1$redlistCriteria)
# tmp1$redlistCriteria <- stringr::str_squish(tmp1$redlistCriteria)
# tmp1.igual <- table(tmp1$main.criteria[tmp1$cat.reg.clean == tmp1$redlistCategory], 
#                     tmp1$redlistCriteria[tmp1$cat.reg.clean == tmp1$redlistCategory])
# tmp1.diff <- table(tmp1$main.criteria[tmp1$cat.reg.clean != tmp1$redlistCategory], 
#                    tmp1$redlistCriteria[tmp1$cat.reg.clean != tmp1$redlistCategory])
# 
# sum(tmp1.diff)/dim(tmp1)[1] # before: ?; now 52%
# sum(tmp1.igual)/dim(tmp1)[1] # before: ?;  now 48%
# 
# 100 * sum(tmp1.diff[grepl("A", row.names(tmp1.diff)),]) / sum(tmp1.diff) #68% das diferenças foram para espécies avaliadas usando o criterio A, entre outros; agora 50% 
# 100 * sum(tmp1.diff[c("A1+A2","A2"),]) / sum(tmp1.diff) #58% das diferenças foram para espécies avaliadas usando apenas o criterio A; agora 38% 
# 100 * sum(tmp1.diff[!grepl("A", row.names(tmp1.diff)),]) / sum(tmp1.diff) # before: 32%; now: 50%
# 
# 100 * sum(tmp1.igual[grepl("B", row.names(tmp1.igual)),]) / sum(tmp1.igual) #32% das igualdades usaram o criterio B entre outros; now 63%
# 100 * sum(tmp1.igual[c("B1","B1+B2","B2"),]) / sum(tmp1.igual) #25% apenas o criterio B; now 55%
# 100 * sum(tmp1.igual[!grepl("B", row.names(tmp1.igual)),]) / sum(tmp1.igual) #68% das igualidades vieram do critério A; now 37%
# 
# sum(
#   table(tmp1$main.criteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean != tmp1$redlistCategory], 
#         tmp1$redlistCriteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean != tmp1$redlistCategory])[c("A1+A2","A2"),]
# )
# 
# tmp1.all.cats <- table(tmp1$cat.reg.clean, tmp1$redlistCategory)


## COMPARISON OF PREVIOUS AND NEW ASSESS. RLI ONLY FOR CRITERIA B
## Another approach: comparing directly the classification using only criteria B
tmp2 <- all.crit[(grepl("B",all.crit$redlistCriteria) & !grepl("A|C|D",all.crit$redlistCriteria) |
                    all.crit$redlistCategory %in% "LC") & all.crit$endemic %in% "endemic",]
tmp2$remove <- FALSE
tmp2$remove[grepl("This species has a large population size", tmp2$population)] <- TRUE
tmp2$remove[grepl("population size is estimated to be large", tmp2$population)] <- TRUE
tmp2 <- tmp2[!tmp2$remove, ]
tmp2$redlistCriteria <- gsub("i|v|a|b|c|d|\\,|\\(|\\)", "", tmp2$redlistCriteria)
table(tmp2$redlistCriteria, useNA = "always")

##Consensus based on B only
assess.df <- tmp2[,c("species", "B1","B2")]
tmp <- ConR:::cat_mult_criteria(assess.df)
table(tmp$species == tmp2$species)
result <- tmp[,c("category","main.criteria","aux.criteria")]
names(result) <- c("category_B_consensus","main.criteria_B","aux.criteria_B")
tmp2 <- cbind.data.frame(tmp2, result, stringsAsFactors = FALSE)
rm(assess.df)

100 * dim(tmp2)[1] / 
  dim(all.crit[all.crit$endemic %in% "endemic" &
                 !is.na(all.crit$redlistCategory),])[1] #83% of previous global assessments were made only using criterion B

tmp2.igual <- table(tmp2$category_B_consensus[tmp2$category_B_consensus == tmp2$redlistCategory], 
                    tmp2$redlistCategory[tmp2$category_B_consensus == tmp2$redlistCategory])
tmp2.diff <- table(tmp2$category_B_consensus[tmp2$category_B_consensus != tmp2$redlistCategory], 
                   tmp2$redlistCategory[tmp2$category_B_consensus != tmp2$redlistCategory])
sum(tmp2.diff)/dim(tmp2)[1] # 31.9%
sum(tmp2.igual)/dim(tmp2)[1] # 68.1%

table(tmp2$category_B_consensus, tmp2$redlistCategory, 
      dnn = c("new", "previous"))[c(1,2,5,4,3),c(1,2,5,4,3)]
round(100*table(tmp2$category_B_consensus, tmp2$redlistCategory, 
      dnn = c("new", "previous"))[c(1,2,5,4,3),c(1,2,5,4,3)]/dim(tmp2)[1],2)

#comparing the RLI fro criteria B only
red::rli(tmp2$redlistCategory, boot = TRUE, runs = 50000) #0.753
red::rli(tmp2$category_B_consensus, boot = TRUE, runs = 50000) #0.829


## LARGE CHANGES BETWEEN CATEGORIES (TWO STESP OR MORE UP OR DOWN)
#same matrix used to build Figure2a:
tmp <- table(paste0(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],"_new"), 
             all.crit$redlistCategory[all.crit$endemic %in% "endemic"])
#LC to CR
100 * sum(tmp[c("CR_new"), c("LC","NT")]) / sum(tmp) # now: 4.9%
#Threatened to not Threatened
100 * sum(tmp[c("LC_new","NT_new"), c("CR","EN","EW","EX","VU")]) / sum(tmp) # now: 2.8%

## EXPLAINING CHANGES BETWEEN PREVIOUS AND NEW CRITERIA
#Which species were up-listed from not threatened to CR?
up.spp <- all.crit$species[all.crit$cat.reg.clean %in% c("CR") & 
                             all.crit$redlistCategory %in% c("LC","NT") &
                             all.crit$endemic %in% "endemic"]
tmp0 <- all.crit[all.crit$species %in% up.spp, ]
tmp0$redlistCriteria <- gsub("i|v|a|b|c|d|\\,|\\(|\\)", "", tmp0$redlistCriteria)
tmp0$redlistCriteria <- stringr::str_squish(tmp0$redlistCriteria)
tmp0$main.criteria <- gsub("B1\\+B2$", "B1+2", tmp0$main.criteria, perl = TRUE)
tmp0[, cols]
table(tmp0$EOO > 20000, useNA = "always")
table(tmp0$AOO > 2000, useNA = "always")
tail(sort(table(tmp0$rationale)))


#Which species were down-listed for being not threatened
down.spp <- all.crit$species[all.crit$cat.reg.clean %in% c("LC","NT") & 
                               all.crit$redlistCategory %in% c("CR","EN","EW","EX","VU") &
                               all.crit$endemic %in% "endemic"]
tmp0 <- all.crit[all.crit$species %in% down.spp, ]
tmp0$redlistCriteria <- gsub("i|v|a|b|c|d|\\,|\\(|\\)", "", tmp0$redlistCriteria)
tmp0$redlistCriteria <- stringr::str_squish(tmp0$redlistCriteria)
tmp0$main.criteria <- gsub("B1\\+B2$", "B1+2", tmp0$main.criteria, perl = TRUE)
tmp0[, cols]
tmp0[, c("species","population")]
tmp0[tmp0, c("species","redlistCategory","redlistCriteria","yearPublished")]
#Why species were downlisted: nos pop decline over 30% and/or no severe fragmentation/ <10 locations
toto <- 100 * table(tmp0[,c("redlistCriteria")])/dim(tmp0)[1]
table(tmp0[,c("redlistCriteria", "reduction_A12")])
table(tmp0$redlistCriteria, tmp0$EOO < 20000)
table(tmp0$redlistCriteria, tmp0$AOO <= 2000, useNA = "always")
table(tmp0$redlistCriteria, tmp0$reduction_A12 < 30, useNA = "always")


###################################################
#### THE INFLUENCE OF USING DIFFERENT CRITERIA ####
###################################################
##How many threatened species classified in IUCN 2021 using criteria B only?
iucn <- read.csv("IUCN_2021_v1_assessments.csv", as.is = TRUE)
iucn$redlistCriteria <- gsub("i|v|a|b|c|d|e|\\,|\\(|\\)", "", iucn$redlistCriteria)
iucn$redlistCriteria <- stringr::str_squish(iucn$redlistCriteria)
onlyA <- grepl("A", iucn$redlistCriteria) & 
  !grepl("B|C|D", iucn$redlistCriteria) & !iucn$redlistCriteria %in% ""
onlyB <- grepl("B", iucn$redlistCriteria) & 
  !grepl("A|C|D", iucn$redlistCriteria) & !iucn$redlistCriteria %in% ""
onlyC <- grepl("C", iucn$redlistCriteria) & 
  !grepl("A|B|D", iucn$redlistCriteria) & !iucn$redlistCriteria %in% ""
onlyD <- grepl("D", iucn$redlistCriteria) & 
  !grepl("A|C|B", iucn$redlistCriteria) & !iucn$redlistCriteria %in% ""
# threatened species: 23860; not threatened: 34069   
table(iucn$redlistCriteria[onlyA]) # 2389
table(iucn$redlistCriteria[onlyB]) # 15496 -> 64.95% of all threatened
table(iucn$redlistCriteria[onlyC]) # 553
table(iucn$redlistCriteria[onlyD]) # 3310
rm(iucn)

##Species classified under all criteria vs. criterion A and criterion B
tmp <- all.crit[!is.na(all.crit$reduction_A12) &  # Assessed under criterion A
                  !is.na(all.crit$AOO) & # Assessed under criterion B
                  !all.crit$cat.reg.clean %in% "NA",]
dim(tmp)[1] ## 2757 species; now 2698 species
100*dim(tmp)[1]/dim(all.crit)[1]
# tmp$cat.reg.clean[tmp$cat.reg.clean %in% c("LC","NT")] <- "LC+NT"

subcriteria <- c("category_A", "category_B","category_C","category_D")
for(i in 1:length(subcriteria)) {
  tmp[,subcriteria[i]] <- 
    ConR::near.threatened(tmp[,subcriteria[i]],
                          tmp$EOO, tmp$AOO, tmp$declineB,
                          tmp$reduction_A12, tmp$pop.size, tmp$pop.size.low,
                          tmp$nbe_loc_total, tmp$sever.frag, ext.fluct = NULL,
                          subpop = tmp$Nbe_subPop, subcriteria = subcriteria1[i])
}

## ALL SPECIES ##
mat1 <- table(tmp[, c("cat.reg.clean","category_A")]) ## All vs. A
mat2 <- table(tmp[, c("cat.reg.clean","category_B")]) ## All vs. B
mat3 <- table(tmp[, c("cat.reg.clean","category_C")]) ## All vs. C
mat4 <- table(tmp[, c("cat.reg.clean","category_D")]) ## All vs. D
mat5 <- table(tmp[, c("category_A","category_B")]) ## A vs. B

#RLI
rl0 <- round(red::rli(tmp$cat.reg.clean, boot = TRUE, runs = 50000),4)
rl1 <- round(red::rli(tmp$category_A, boot = TRUE, runs = 50000),4)
rl2 <- round(red::rli(tmp$category_B, boot = TRUE, runs = 50000),4)
rl3 <- round(red::rli(tmp$category_C, boot = TRUE, runs = 50000),4)
rl4 <- round(red::rli(tmp$category_D, boot = TRUE, runs = 50000),4)

#Number of threatened species
ts0 <- round(100*sum(mat1[-3,]) / sum(mat1),2) #94.8% of threatened species; now 95.89%
ts1 <- round(100*sum(mat1[,-3]) / sum(mat1),2) #94.5% of threatened species; now 95.59%
ts2 <- round(100*sum(mat2[,-3]) / sum(mat2),2) #10.2% of threatened species; now 9.97%
ts3 <- round(100*sum(mat3[,-3]) / sum(mat3),2) #3.2% of threatened species; now 4.52%
ts4 <- round(100*sum(mat4[,-3]) / sum(mat4),2) #1.7% of threatened species; now 2.22%

#Categories of threat
cg1 <- round(100*sum(diag(mat1))/sum(mat1),2) #95.7% of congruence in the assessments; now 95.7%
cg2 <- round(100*sum(diag(mat2))/sum(mat2),2) #10.9% of congruence in the assessments; 9.93%
cg3 <- round(100*sum(diag(mat3))/sum(mat3),2) #5.6% of congruence in the assessments; now 5.3%
cg4 <- round(100*sum(diag(mat4))/sum(mat4),2) #5.3% of congruence in the assessments; now 4.93%
cg5 <- round(100*sum(diag(mat5))/sum(mat5),2) #9.3% of congruence in the assessments; now 8.82%


## ENDEMIC SPECIES ##
tmp1 <- tmp[tmp$endemic %in% "endemic",]
dim(tmp1)[1] ## 1646 species; now 1586
100*dim(tmp1)[1]/dim(all.crit)[1] # 32.3% of the total; now 32.02%

mat1.1 <- table(tmp1[, c("cat.reg.clean","category_A")]) ## All vs. A
mat2.1 <- table(tmp1[, c("cat.reg.clean","category_B")]) ## All vs. B
mat3.1 <- table(tmp1[, c("cat.reg.clean","category_C")]) ## All vs. C
mat4.1 <- table(tmp1[, c("cat.reg.clean","category_D")]) ## All vs. D
mat5.1 <- table(tmp1[, c("category_A","category_B")]) ## A vs. B

#RLI
rl0.1 <- round(red::rli(tmp1$cat.reg.clean, boot = TRUE, runs = 50000),4)
rl1.1 <- round(red::rli(tmp1$category_A, boot = TRUE, runs = 50000),4)
rl2.1 <- round(red::rli(tmp1$category_B, boot = TRUE, runs = 50000),4)
rl3.1 <- round(red::rli(tmp1$category_C, boot = TRUE, runs = 50000),4)
rl4.1 <- round(red::rli(tmp1$category_D, boot = TRUE, runs = 50000),4)

#Number of threatened species
ts0.1 <- round(100*sum(mat1.1[-3,]) / sum(mat1.1),2) #93.0% of threatened species; now 94.26%
ts1.1 <- round(100*sum(mat1.1[,-3]) / sum(mat1.1),2) #92.2% of threatened species; now 93.88%
ts2.1 <- round(100*sum(mat2.1[,-3]) / sum(mat2.1),2) #15.0% of threatened species; now 15.13%
ts3.1 <- round(100*sum(mat3.1[,-3]) / sum(mat3.1),2) #2.6% of threatened species; now 4.67%
ts4.1 <- round(100*sum(mat4.1[,-2]) / sum(mat4.1),2) #2.3% of threatened species; now 2.96%

#Categories of threat
cg1.1 <- round(100*sum(diag(mat1.1))/sum(mat1.1),2) #96.9% of congruence in the assessments; now 96.97%
cg2.1 <- round(100*sum(diag(mat2.1))/sum(mat2.1),2) #15.3% of congruence in the assessments; now 14.75%
cg3.1 <- round(100*sum(diag(mat3.1))/sum(mat3.1),2) #7.4% of congruence in the assessments; now 7.38%
cg4.1 <- round(100*sum(diag(mat4.1[c(1,3,4),c(1,2,3)]))/sum(mat4.1),2) #7.1% of congruence in the assessments; now 7.12%
cg5.1 <- round(100*sum(diag(mat5.1))/sum(mat5.1),2) #9.3% of congruence in the assessments; now 13.11%

## Table (S)1? ##
tabS1 <- matrix(c(ts1, paste0(rl1[2]," [",paste0(round(rl1[c(1,3)],3), collapse = "-"),"]"), ts1.1, paste0(rl1.1[2]," [",paste0(round(rl1.1[c(1,3)],3), collapse = "-"),"]"), #cg1.1,
                  ts2, paste0(rl2[2]," [",paste0(round(rl2[c(1,3)],3), collapse = "-"),"]"), ts2.1, paste0(rl2.1[2]," [",paste0(round(rl2.1[c(1,3)],3), collapse = "-"),"]"), #cg2.1,
                  ts3, paste0(rl3[2]," [",paste0(round(rl3[c(1,3)],3), collapse = "-"),"]"), ts3.1, paste0(rl3.1[2]," [",paste0(round(rl3.1[c(1,3)],3), collapse = "-"),"]"), #cg3.1,
                  ts4, paste0(rl4[2]," [",paste0(round(rl4[c(1,3)],3), collapse = "-"),"]"), ts4.1, paste0(rl4.1[2]," [",paste0(round(rl4.1[c(1,3)],3), collapse = "-"),"]"),#, cg4.1), 
                  ts0, paste0(rl0[2]," [",paste0(round(rl0[c(1,3)],3), collapse = "-"),"]"), ts0.1, paste0(rl0.1[2]," [",paste0(round(rl0.1[c(1,3)],3), collapse = "-"),"]")),
                nrow = 5, ncol = 4, byrow = TRUE,
                dimnames = list(c("catA","catB","catC","catD","all"), c("Threatened","RLI","Threatened","RLI")))
write.csv(tabS1,"TableS1.csv")



## How many threatened species for criteria A and B?
# All populations
sum(table(all.crit$category_A)[c("CR", "EN", "VU")]) # 2980 threatened
round(100 * sum(table(all.crit$category_A)[c("CR", "EN", "VU")]) / dim(all.crit)[1], 1) # 60.2% threatened 
sum(table(all.crit$category_B)[c("CR", "EN", "VU")]) # 775 threatened
round(100 * sum(table(all.crit$category_B)[c("CR", "EN", "VU")]) / dim(all.crit)[1], 1) # 15.6% threatened 

# Endemics
tmp.end <- all.crit[all.crit$endemic %in% "endemic",]
sum(table(tmp.end$category_A)[c("CR", "EN", "VU")]) # 1463 threatened
round(100 * sum(table(tmp.end$category_A)[c("CR", "EN", "VU")]) / dim(tmp.end)[1], 1) # 59.4% threatened 
sum(table(tmp.end$category_B)[c("CR", "EN", "VU")]) # 614 threatened
round(100 * sum(table(tmp.end$category_B)[c("CR", "EN", "VU")]) / dim(tmp.end)[1], 1) # 24.9% threatened 



#### INFLUENCES OF DATA UNCERTAINTY ####

## GL AND P ##
#Already in the criterion C and D codes (only figures)

## SPECIES IDENTIFICATION
summary(all.crit$prop.high.tax) # mean 68%; median 72%
100 * table(all.crit$prop.high.tax >= 1, useNA = "always") /
  sum(table(all.crit$prop.high.tax >= 1, useNA = "always")) # 34% with tax.conf >= 1

#How many classified as threatened
tmp <- all.crit[!is.na(all.crit$category_B.high.tax) &
                  !is.na(all.crit$category_B),]
100 * table(tmp$category_B)/sum(table(tmp$category_B)) # 15.6% threatened
100 * table(tmp$category_B.high.tax)/sum(table(tmp$category_B)) # 25.2% threatened

subcriteria <- c("category_B","category_B.high.tax")
for(i in 1:length(subcriteria)) {
  tmp[,subcriteria[i]] <- 
    ConR::near.threatened(tmp[,subcriteria[i]],
                          tmp$EOO, tmp$AOO, tmp$declineB,
                          tmp$reduction_A12, tmp$pop.size, tmp$pop.size.low,
                          tmp$nbe_loc_total, tmp$sever.frag, ext.fluct = NULL,
                          subpop = tmp$Nbe_subPop, subcriteria = NULL)
}

## RLI ##
(opt.rli <- round(red::rli(tmp$category_B, boot = TRUE, runs = 50000),4))
(high.rli <- round(red::rli(tmp$category_B.high.tax, boot = TRUE, runs = 50000),4))

## Finding a confidence level cutoff 
cut <- seq(0.05,0.95, by=0.05)
res <- matrix(NA, nrow = length(cut), ncol = 6)
rownames(res) <- cut
for (i in 1:length(cut)) {
  res[i,1:3] <- round(red::rli(tmp$category_B[tmp$prop.high.tax >= cut[i]], boot = TRUE, runs = 50000),4)
  res[i,4:6] <- round(red::rli(tmp$category_B.high.tax[tmp$prop.high.tax >= cut[i]], boot = TRUE, runs = 50000),4)
}
res
names(res) <- c("Low.Opt", "Median.Opt", "High.Opt", 
                "Low.High", "Median.High", "High.High")
saveRDS(res, "data/tax_conf_effect_on_RLI.rds")
## See Figure SU

###################################
#### RECORDS IN TIME AND SPACE ####
###################################

#How many records were used in the assessments
#Overall number of records retrieved for the AF tree flora: 3,079,350
oc.data <- readRDS("data/threat_occ_data_final.rds")
table(oc.data$source) # 816,192 herbarium records + 91463 additional records from TreeCo
100*table(oc.data$source) / 3079350 # total 26.51%
100*table(oc.data$source[oc.data$tax.check2 %in% "TRUE"]) / 
  sum(table(oc.data$tax.check2[oc.data$source %in% "herbaria"])) # herbarium records with valid identification 37.69%

#Correcting problematic years
oc.data[coly %in% "2997", coly := 1997] 
oc.data[coly %in% "2992", coly := 1992]
oc.data[coly %in% "2188", coly := 1988]
oc.data[coly %in% "2097", coly := 1997]
oc.data[coly %in% "1012", coly := 2012] 
oc.data[coly %in% "1394", coly := 1934]
oc.data[coly %in% "1089", coly := 1989]
oc.data[coly %in% "1076", coly := 1976]
oc.data[coly %in% "1014", coly := 2014]
oc.data[coly %in% "1190", coly := 1990]
oc.data[coly %in% "1235", coly := 2012]
oc.data[coly %in% "1197", coly := 1997]
oc.data[coly %in% "1192", coly := 1992]
oc.data[coly %in% "1071", coly := 1971]

#How many valid occurrences
dim(oc.data)[1] # before 904,776; now 907,655
#How many valid occurrences used in the assessments
table(oc.data$tax.check.final %in% "high") # 446,077; now: 444,834
100 * table(oc.data$tax.check.final %in% "high")/dim(oc.data)[1]

## Most of species occurrences come from which period?
tmp1 <- oc.data[oc.data$tax.check.final %in% "high",] 
tmp2 <- table(tmp1$coly)

100 * table(is.na(tmp1$coly)) / 
  sum(table(tmp1$coly, useNA = "always"))
1 - round(cumsum(table(oc.data$coly)) / 
            sum(table(oc.data$coly)),2)["2000"]
1 - round(cumsum(tmp2) / 
            sum(tmp2),2)["1990"]
1 - round(cumsum(tmp2) / 
            sum(tmp2),2)["1980"]
1 - round(cumsum(tmp2) / 
            sum(tmp2),2)["1970"]

plot(cumsum(tmp2) / 
       sum(tmp2), xaxt = "n")
axis(1, at = which(names(tmp2) %in% seq(1840,2020,10)), label= seq(1840,2020,10))

quantile(tmp1$coly, probs = seq(0,1,by=0.05), na.rm = TRUE)
diff(quantile(tmp1$coly, probs = seq(0,1,by=0.05), na.rm = TRUE))


## How many species not recorded in the last 50 years
dim(all.crit[!is.na(all.crit$last.record) & all.crit$last.record < (2018-50) &
               all.crit$endemic %in% "endemic", ])[1] #41 endemics without records in the past 50 years
all.crit[!is.na(all.crit$last.record) & 
           all.crit$last.record < (2018-50) &
           all.crit$endemic %in% "endemic", cols]
dim(all.crit[!is.na(all.crit$last.record) & 
               all.crit$last.record < (2018-50) &
               !is.na(all.crit$endemic) &
               all.crit$endemic %in% "endemic" &
               all.crit$only.type %in% TRUE, cols])[1] # 14 of the only know from the type specimen
all.crit[all.crit$endemic %in% "endemic" &
           all.crit$category %in% "CR_PE", cols]


## Spatial distribution across the AF: see codes 'data/14.5_grid_summaries' and 'data/15_maps.R'


## Occurrence insided protected areas
table(all.crit$protected[!is.na(all.crit$protected)] > 0) / # 89.3% have at least one occurrence inside UCs de PI
  dim(all.crit)[1] 
table(all.crit$protected[!is.na(all.crit$protected) & all.crit$endemic %in% "endemic"]>0) / 
  dim(all.crit[all.crit$endemic %in% "endemic",])[1] # 84.3% of the endemics have at least one occurrence inside UCs de PI


table(all.crit$protected[!is.na(all.crit$protected) & 
                           all.crit$endemic %in% "endemic" & 
                           all.crit$cat.reg.clean %in% c("CR","EN","VU")] > 0) / 
  dim(all.crit[all.crit$endemic %in% "endemic"& 
                 all.crit$cat.reg.clean %in% c("CR","EN","VU"),])[1] # before 82%; now 82.7% of the endangered endemics have at least one occurrence inside UCs de PI

summary(all.crit$protected)
summary(all.crit$protected[all.crit$endemic %in% "endemic"])
table(all.crit$protected[all.crit$endemic %in% "endemic"] < 25)/
  dim(all.crit[all.crit$endemic %in% "endemic",])[1] # 74.2% of the endemics have less than 25% of records inside UCS of PI

summary(all.crit$prop.EOO.in.StrictUCs)
summary(all.crit$prop.EOO.in.StrictUCs[all.crit$endemic %in% "endemic"])
table(all.crit$prop.EOO.in.StrictUCs[all.crit$endemic %in% "endemic"] < 10)/
  dim(all.crit[all.crit$endemic %in% "endemic",])[1] # 84.3% of the endemics have less than 10% of their EOO inside UC of PI


## COMPARING % PROTECT BETWEEN THREAT CATEGORIES
## SEE FIGURE 4
tmp <- na.omit(all.crit[all.crit$endemic %in% "endemic",
                        c("protected","prop.EOO.in.StrictUCs","cat.reg.clean")])
tmp$cats.f <- factor(tmp$cat.reg.clean,
                     levels = c("LC","NT","VU","EN","CR"))
levels(tmp$cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")
## proportion of records inside PAs
model <- lm(log(protected+1) ~ cats.f,
            data = all.crit[all.crit$endemic %in% "endemic",])
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")
LABELS[c("CR", "EN", "VU", "LC+NT"),]
summary(model)
car::Anova(model)
plot(TUKEY)

## proportion of terrestrial EOO within PAs
model <- lm(log(prop.EOO.in.StrictUCs+1) ~ cats.f,
            data = all.crit[all.crit$endemic %in% "endemic",])
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")
LABELS[c("CR", "EN", "VU", "LC+NT"),]
summary(model)
car::Anova(model)
plot(TUKEY)

##MULTINOMIAL REGRESSIONS (NOT USING SO FAR)
# library(brant)
# library(nnet)
# 
# # ordered logit model
# mod0 <- MASS::polr(cats.f ~ 1, 
#                   data = tmp, Hess=TRUE)
# mod <- MASS::polr(cats.f ~ protected + prop.EOO.in.StrictUCs, 
#             data = tmp, Hess=TRUE)
# summary(mod)
# anova(mod0, mod)
# car::Anova(mod)
# pr <- profile(mod)
# confint(pr)
# plot(pr)
# pairs(pr)
# 
# brant(mod)
# 
# # Run the model
# model0 <- nnet::multinom(cats.f ~ 1, 
#                         data =  tmp, Hess = TRUE)
# model <- nnet::multinom(cats.f ~ protected + prop.EOO.in.StrictUCs, 
#                         data =  tmp, Hess = TRUE)
# summary(model)
# bbmle::AICtab(model0, model)
# coef(model)
# # Calculate z-values
# zvalues <- summary(model)$coefficients / summary(model)$standard.errors
# zvalues
# #p-values
# pnorm(abs(zvalues), lower.tail=FALSE)*2


### THE REMAINING PROPORTION OF FOREST WITHIN EOO
summary(all.crit$prop.EOO.forest[all.crit$endemic %in% "endemic" ]) # median 21.1%; mean 29.5%
summary(all.crit$prop.EOO.forest[all.crit$endemic %in% "endemic"  & # median 21.0%; mean 29.7% 
                                   all.crit$cat.reg.clean %in% c("CR","EN","VU")])
hist(all.crit$prop.EOO.forest[all.crit$endemic %in% "endemic"  &  # very long tail!
                                   all.crit$cat.reg.clean %in% c("CR","EN","VU")],
     nclass = 40)

table(all.crit$prop.EOO.forest[!is.na(all.crit$prop.EOO.forest) & 
                                 all.crit$endemic %in% "endemic"] < 40) / 
  dim(all.crit[all.crit$endemic %in% "endemic",])[1] # 73.9% of the endemics have less at least one occurrence inside UCs de PI
table(all.crit$prop.EOO.forest[!is.na(all.crit$prop.EOO.forest) & 
                           all.crit$endemic %in% "endemic" & 
                           all.crit$cat.reg.clean %in% c("CR","EN","VU")] < 40) / 
  dim(all.crit[all.crit$endemic %in% "endemic"& 
                 all.crit$cat.reg.clean %in% c("CR","EN","VU"),])[1] # 73.0% of the endemics have less at least one occurrence inside UCs de PI

cats.f <- factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],
                 levels = c("LC","NT","VU","EN","CR"))
levels(cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")
# boxplot(prop.EOO.forest ~ cats.f,
#         data = all.crit[all.crit$endemic %in% "endemic",],
#         notch = TRUE, varwidth = TRUE, col = c("yellowgreen","gold","darkorange","red"))
# abline(h=c(25,40),lty=2)

tmp <- all.crit[all.crit$endemic %in% "endemic",]
model <- lm(log(prop.EOO.forest+1) ~ cats.f, 
          data = tmp)
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")
LABELS[c("CR", "EN", "VU", "LC+NT"),]
summary(model)
car::Anova(model)
plot(TUKEY)


aggregate(all.crit$prop.EOO.forest[all.crit$endemic %in% "endemic"],
          list(cats.f), summary)


### THE AVERAGE HII WITHIN SPECIES EOO
summary(all.crit$Median.HII[all.crit$endemic %in% "endemic"])
quantile(all.crit$Median.HII[!is.na(all.crit$Median.HII) & all.crit$endemic %in% "endemic"],
         prob = c(0.1,0.5,0.9))

boxplot(Median.HII ~ cat.reg.clean, 
        data = all.crit[all.crit$endemic %in% "endemic",], 
        notch = TRUE)
cats.f <- factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],
                 levels = c("LC","NT","VU","EN","CR"))
levels(cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")
boxplot(Median.HII ~ cats.f, ylim=c(0,65),
        data = all.crit[all.crit$endemic %in% "endemic",], 
        notch = TRUE, varwidth = TRUE, col = c("yellowgreen","gold","darkorange","red"))

tmp <- all.crit[all.crit$endemic %in% "endemic",]
mod <- lm(log(Median.HII+1) ~ cats.f, 
          data = tmp)
summary(mod); car::Anova(mod)
confint(lm(Median.HII ~ cats.f - 1, data = tmp))
out = TukeyHSD(aov(mod))
plot(out)

##Examples high HII
all.crit[!is.na(all.crit$Median.HII) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$Median.HII>35,c("species","prop.EOO.in.StrictUCs","Median.HII","prop.EOO.forest","redlistCategory","main.criteria","cat.reg.clean")]

##Examples low HII
all.crit[!is.na(all.crit$Median.HII) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$Median.HII<5,c("species","prop.EOO.in.StrictUCs","Median.HII","prop.EOO.forest","redlistCategory","main.criteria","cat.reg.clean")]


## BOTH
plot(Median.HII ~ prop.EOO.forest, data = tmp)
plot(Median.HII ~ prop.EOO.forest, data = tmp, log = "xy")







