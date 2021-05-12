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

######################################################
#### THE CONSERVATION STATUS OF THE AF TREE FLORA ####
######################################################

## HOW MANY SPECIES ASSESSED?
dim(all.crit)[1] # 5094

## HOW MANY THREATENED SPECIES?
## Overall
sum(table(all.crit$cat.reg.clean)[c("CR", "EN", "VU")]) # 3408 threatened
round(100 * sum(table(all.crit$cat.reg.clean)[c("CR", "EN", "VU")]) / dim(all.crit)[1], 1) # 66.9 threatened 
## Endemic
round(100 * sum(table(all.crit$endemic)["endemic"]) / dim(all.crit)[1], 1) # 50.8 are pure and near endemics
round(100 * sum(table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])[c("CR", "EN", "VU")]) / 
        dim(all.crit[all.crit$endemic %in% "endemic",])[1], 1) # 84.4% of the endemics are threatened
sum(table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])[c("CR", "EN", "VU")]) #2186 species

## THE RLI FOR THE ATLANTIC FOREST TREE FLORA 
red::rli(all.crit$cat.reg.clean, boot = TRUE, runs = 50000)
red::rli(all.crit$cat.reg.clean[!all.crit$cat.reg.clean %in% "NA"], boot = TRUE, runs = 50000)

## MORE FREQUENT CATEGORIES AND CRITERIA
round(100 * table(all.crit$cat.reg.clean)[c("CR", "EN", "VU")] / dim(all.crit)[1], 1) # 43.1% of EN
round(100 * table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])[c("CR", "EN", "VU")] / 
        dim(all.crit[all.crit$endemic %in% "endemic",])[1], 1) # 55.9% of EN

## How many near threatened
round(100 * sum(table(all.crit$cat.reg.clean)[c("NT")]) / dim(all.crit)[1], 1) #1.4% of NT

## SOME EXAMPLES OF THREATENED SPECIES WITH IMPORTANT USES
cols <- c("species","reduction_A12","reduction_A12.25ys","EOO","AOO","nbe_loc_total","sever.frag","declineB","any.decline","pop.size",
          "category","cat.reg.clean","main.criteria","endemism.level.1","endemism.level.2","endemic","last.record","only.type",
          "status.reflora","redlistCategory")
iconic <- c("Araucaria angustifolia","Euterpe edulis","Ilex paraguariensis","Paubrasilia echinata")
all.crit[all.crit$species %in% iconic, cols]
timber <- c("Cariniana legalis","Dalbergia nigra","Melanoxylon brauna","Myrocarpus frondosus","Ocotea porosa","Parapiptadenia rigida","Paratecoma peroba")
all.crit[all.crit$species %in% timber, cols]
# food.orn <- c("Eugenia brasiliensis","Handroanthus chrysotrichus","Ilex paraguariensis")
# all.crit[all.crit$species %in% food.orn, cols]



#########################
#### IUCN CATEGORIES ####
#########################
tmp0 <- table(all.crit$main.criteria[!all.crit$cat.reg.clean %in% c("LC","NA","DD","NT")], useNA = "always")
tmp <- sort(round(100 *  tmp0 / 
                    dim(all.crit[!all.crit$cat.reg.clean %in% c("LC","NA","DD","NT"),])[1], 4))
sum(tmp[grepl("A2", names(tmp))]) #74.9% of the threatened species had A2 within its main criteria
sum(tmp[grepl("B2", names(tmp))]) #28.4% of the threatened species had B2 within its main criteria
sum(tmp[grepl("C", names(tmp))]) #0.24% of the threatened species had C within its main criteria
sum(tmp[grepl("D", names(tmp))]) #0.03% of the threatened species had D within its main criteria
sum(tmp0[grepl("C", names(tmp0))]) #80.24% of the threatened species had C within its main criteria
sum(tmp0[grepl("D", names(tmp0))]) #0.03% of the threatened species had D within its main criteria

## MEAN POPULATION REDUCTIONS
summary(all.crit$reduction_A12); hist(all.crit$reduction_A12, nclass=40)
plot(all.crit$reduction_A12 ~ all.crit$reduction_A12.25ys); abline(0,1); abline(v=0,h=0,lty=2)
confint(lm(all.crit$reduction_A12 ~ 1))
table(all.crit$reduction_A12>=90)
table(all.crit$reduction_A12>=80)
table(all.crit$reduction_A12>=50)
table(all.crit$reduction_A12>=30)
#Common AF species with large declines
all.crit[!is.na(all.crit$reduction_A12) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$reduction_A12>=90 &
           all.crit$pop.size>=2e6, cols]
#Species with increasing populations
all.crit[!is.na(all.crit$reduction_A12) & 
           all.crit$endemic %in% "endemic" & 
           all.crit$reduction_A12<0 &
           all.crit$reduction_A12.25ys<0, cols]


## SMALL EOO AND AOO
summary(all.crit$AOO); hist(log10(all.crit$AOO), nclass=40); abline(v=log10(c(10,500,2000)), col=c("red","darkorange","gold"))
confint(lm(all.crit$AOO ~ 1))
summary(all.crit$nbe_loc_total); hist(log10(all.crit$nbe_loc_total), nclass=40); abline(v=log10(c(1,5,10)), col=c("red","darkorange","gold"))
confint(lm(all.crit$nbe_loc_total ~ 1))
table(all.crit$nbe_loc_total>10)/50.94
table(all.crit$sever.frag>10)/50.94

plot(log10(all.crit$AOO) ~ log10(all.crit$nbe_loc_total)); abline(h=log10(c(10,500,2000)), lty=2, col=c("red","darkorange","gold")); abline(v=log10(c(1,5,10)),lty=2, col=c("red","darkorange","gold"))
all.crit[!is.na(all.crit$reduction_A12) & !is.na(all.crit$EOO) & !is.na(all.crit$AOO) & 
           all.crit$endemic %in% "endemic" & all.crit$cat.reg.clean %in% "CR" &
           all.crit$reduction_A12>85 &
           #all.crit$EOO<30000 &
           all.crit$AOO<100 &
           all.crit$pop.size<2500, cols]

## EXAMPLES OF HIGHLY THREATENED SPECIES
#Endemics classified as threatened under all criteria
tmp <- all.crit[all.crit$category_A %in% c("CR", "EN", "VU") &
                  all.crit$category_B %in% c("CR", "EN", "VU") & 
                  all.crit$category_C %in% c("CR", "EN", "VU") &
                  all.crit$category_D %in% c("CR", "EN", "VU"), ]
tmp$species[tmp$endemic %in% "endemic"]

#Endemics with the largest declines and smaller distributions
all.crit[!is.na(all.crit$reduction_A12) & !is.na(all.crit$EOO) & !is.na(all.crit$AOO) & 
           all.crit$endemic %in% "endemic" & all.crit$cat.reg.clean %in% "CR" &
           all.crit$reduction_A12>85 &
           #all.crit$EOO<30000 &
           all.crit$AOO<100 &
           all.crit$pop.size<2500, cols]


##############################################
#### COMPARISON WITH PREVIOUS ASSESSMENTS ####
##############################################

## OVERALL ASSESSMENTS 
#How many species with IUCN assessments
table(all.crit$redlistCategory)
100*sum(table(all.crit$redlistCategory))/dim(all.crit)[1] # 28.9% with previous IUCN asses.
table(all.crit$cat.reg.clean, all.crit$redlistCategory)

#How many species with CNCFlora assessments
table(all.crit$status.reflora)
100*sum(table(all.crit$status.reflora))/dim(all.crit)[1]  # 19.7% with previous CNCFlora asses. 

#How many species with national assessments
tmp <- all.crit[!is.na(all.crit$status.reflora) | 
                  !is.na(all.crit$category.ARG) |
                  !is.na(all.crit$category.PAY),]
100 * dim(tmp)[1]/dim(all.crit)[1] # 20.3%

## ENDEMIC SPECIES
#How may species with IUCN assessments
table(all.crit$redlistCategory[all.crit$endemic %in% "endemic"])
100*sum(table(all.crit$redlistCategory[all.crit$endemic %in% "endemic"])) / 
  dim(all.crit[all.crit$endemic %in% "endemic",])[1] # 20.01% with previous IUCN asses.
table(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"], 
      all.crit$redlistCategory[all.crit$endemic %in% "endemic"])

#How may species with CNCFlora assessments
table(all.crit$status.reflora[all.crit$endemic %in% "endemic"])
100*sum(table(all.crit$status.reflora[all.crit$endemic %in% "endemic"])) / 
  dim(all.crit[all.crit$endemic %in% "endemic",])[1]  # 24.5% with previous CNCFlora asses. 


## THEAT NEW ASSESSMENTS
tmp <- all.crit[is.na(all.crit$redlistCategory) &
			is.na(all.crit$status.reflora) &
			is.na(all.crit$category.ARG) &
			is.na(all.crit$category.PAY),]
dim(tmp)[1] # 2959 species (58.09%) without assessments
100 * dim(tmp)[1] / dim(all.crit)[1] # 2959 species (58.09%) without assessments
100 * dim(tmp[tmp$endemic %in% "endemic",])[1] / dim(tmp)[1] # 1631 endemic species (55.12%) without assessments


##Species EX or EW
ex.spp <- all.crit$species[all.crit$redlistCategory %in% c("EX","EW")]
all.crit[all.crit$species %in% ex.spp, cols]
all.crit[all.crit$species %in% ex.spp, c("species","population")]
oc.data <- readRDS("data/threat_occ_data_final.rds")
oc.data <- oc.data[oc.data$tax %in% ex.spp,]
table(oc.data$coly, oc.data$tax)
table(oc.data$coly>=1998, oc.data$tax)
table(oc.data$coly[is.na(oc.data$typeStatus)]>=1998, oc.data$tax[is.na(oc.data$typeStatus)])
oc.data[oc.data$tax %in% "Pouteria stenophylla",]
paths = dir("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//occurrence_data",full.names=TRUE)
paths = paths[grepl('merged_outliers.csv',paths) & grepl('Sapotaceae',paths)]
tmp = data.table::fread(paths,na.strings=c(""," ","NA"))
tmp[tmp$species.correct2 %in% "Pouteria stenophylla",]
rm(tmp,paths,oc.data)

## MAIN CHANGES BETWEEN CATEGORIES
tmp <- table(paste0(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],"_new"), 
      all.crit$redlistCategory[all.crit$endemic %in% "endemic"])
#LC to CR
100 * sum(tmp[c("CR_new"), c("LC","NT")]) / sum(tmp)
#Threatened to not Threatened
100 * sum(tmp[c("LC_new","NT_new"), c("CR","EN","EW","EX","VU")]) / sum(tmp)

#Which species were down-listed for being not threatened
down.spp <- all.crit$species[all.crit$cat.reg.clean %in% c("LC","NT") & 
                               all.crit$redlistCategory %in% c("CR","EN","EW","EX","VU") &
                               all.crit$endemic %in% "endemic"]
all.crit[all.crit$species %in% down.spp, cols]
all.crit[all.crit$species %in% down.spp, c("species","population")]
all.crit[all.crit$species %in% down.spp, c("species","redlistCategory","redlistCriteria","yearPublished")]

## PREVIOUS VS. NEW RLI
#overall
red::rli(all.crit$redlistCategory[!is.na(all.crit$redlistCategory)], boot = TRUE, runs = 50000)
red::rli(all.crit$cat.reg.clean[!is.na(all.crit$redlistCategory)], boot = TRUE, runs = 50000)
#endemics
red::rli(all.crit$redlistCategory[!is.na(all.crit$redlistCategory) & all.crit$endemic %in% "endemic"], boot = TRUE, runs = 50000)
red::rli(all.crit$cat.reg.clean[!is.na(all.crit$redlistCategory) & all.crit$endemic %in% "endemic"], boot = TRUE, runs = 50000)

## COMPARISON OF ALL CRITERIA VS. CRITERIA B
tmp <- table(all.crit$main.criteria[all.crit$endemic %in% "endemic"], 
             all.crit$redlistCriteria[all.crit$endemic %in% "endemic"])
tmp <- table(all.crit$redlistCriteria[all.crit$endemic %in% "endemic"])
100 * sum(tmp[grepl("B",names(tmp))]) / sum(tmp)

tmp1 <- all.crit[grepl("B",all.crit$redlistCriteria) &
                   all.crit$endemic %in% "endemic",]
tmp1$redlistCriteria <- gsub("i|\\,|\\(|\\)", "", tmp1$redlistCriteria)
tmp1.igual <- table(tmp1$main.criteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean == tmp1$redlistCategory], 
      tmp1$redlistCriteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean == tmp1$redlistCategory])
tmp1.diff <- table(tmp1$main.criteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean != tmp1$redlistCategory], 
                    tmp1$redlistCriteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean != tmp1$redlistCategory])

sum(tmp1.diff)/dim(tmp1)[1]
sum(tmp1.igual)/dim(tmp1)[1]

100 * sum(tmp1.diff[grepl("A", row.names(tmp1.diff)),]) / sum(tmp1.diff) #68% das diferenças foram para espécies avaliadas usando o criterio A, entre outros 
100 * sum(tmp1.diff[c("A1+A2","A2"),]) / sum(tmp1.diff) #58% das diferenças foram para espécies avaliadas usando apenas o criterio A 
100 * sum(tmp1.diff[!grepl("A", row.names(tmp1.diff)),]) / sum(tmp1.diff)

100 * sum(tmp1.igual[grepl("B", row.names(tmp1.igual)),]) / sum(tmp1.igual) #32% das igualdades usaram o criterio B entre outros
100 * sum(tmp1.igual[c("B1","B1+B2","B2"),]) / sum(tmp1.igual) #25% apenas o criterio B
100 * sum(tmp1.igual[!grepl("B", row.names(tmp1.igual)),]) / sum(tmp1.igual) #68% das igualidades vieram do critério A


sum(
  table(tmp1$main.criteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean != tmp1$redlistCategory], 
        tmp1$redlistCriteria[tmp1$endemic %in% "endemic" & tmp1$cat.reg.clean != tmp1$redlistCategory])[c("A1+A2","A2"),]
)


###################################
#### RECORDS IN TIME AND SPACE ####
###################################

#How many records were used in the assessments
oc.data <- readRDS("data/threat_occ_data_final.rds")
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
dim(oc.data)[1] # 904,776
#How many valid occurrences used in the assessments
table(oc.data$tax.check.final %in% "high") # 446,077
100 * table(oc.data$tax.check.final %in% "high")/dim(oc.data)[1]

#Most of species occurrences come from which period?
tmp1 <- oc.data[oc.data$tax.check.final %in% "high",] 
100 * table(is.na(tmp1$coly)) / 
  sum(table(tmp1$coly, useNA = "always"))
1 - round(cumsum(table(oc.data$coly)) / 
        sum(table(oc.data$coly)),2)["1989"]
1 - round(cumsum(table(tmp1$coly)) / 
            sum(table(tmp1$coly)),2)["1989"]
1 - round(cumsum(table(tmp1$coly)) / 
            sum(table(tmp1$coly)),2)["1980"]

plot(cumsum(table(tmp1$coly)) / 
        sum(table(tmp1$coly)))
quantile(tmp1$coly, probs = seq(0,1,by=0.05), na.rm = TRUE)
diff(quantile(tmp1$coly, probs = seq(0,1,by=0.05), na.rm = TRUE))


#How many species not recorded in the last 50 years
dim(all.crit[!is.na(all.crit$last.record) & all.crit$last.record < (2018-50) &
               all.crit$endemic %in% "endemic", ])[1] #54 endemics without records in the past 50 years
all.crit[!is.na(all.crit$last.record) & 
           all.crit$last.record < (2018-50) &
           all.crit$endemic %in% "endemic", cols]
dim(all.crit[!is.na(all.crit$last.record) & 
           all.crit$last.record < (2018-50) &
           !is.na(all.crit$endemic) &
           all.crit$endemic %in% "endemic" &
           all.crit$only.type %in% TRUE, cols])[1] # 21 of the only know from the type specimen
all.crit$species[!is.na(all.crit$last.record) & 
           all.crit$last.record < (2018-50) &
           !is.na(all.crit$endemic) &
           all.crit$endemic %in% "endemic" &
           all.crit$only.type %in% TRUE]


# Occurrence insided protected areas
table(all.crit$protected>0) / dim(all.crit)[1]
table(all.crit$protected[all.crit$endemic %in% "endemic"]>0) / 
  dim(all.crit[all.crit$endemic %in% "endemic",])[1]
table(all.crit$protected[all.crit$endemic %in% "endemic" & 
                           all.crit$cat.reg.clean %in% c("CR","EN","VU")]>0) / 
  dim(all.crit[all.crit$endemic %in% "endemic"& 
                 all.crit$cat.reg.clean %in% c("CR","EN","VU"),])[1]

summary(all.crit$protected)
summary(all.crit$protected[all.crit$endemic %in% "endemic"])
summary(all.crit$prop.EOO.in.StrictUCs)
summary(all.crit$prop.EOO.in.StrictUCs[all.crit$endemic %in% "endemic"])


## COMPARING % PROTECT BETWEEN THREAT CATEGORIES
## SEE FIGURE 3
tmp <- na.omit(all.crit[all.crit$endemic %in% "endemic",c("protected","prop.EOO.in.StrictUCs","cat.reg.clean")])
tmp$cats.f <- factor(tmp$cat.reg.clean,
                     levels = c("LC","NT","VU","EN","CR"))
levels(tmp$cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")

mod <- lm(protected ~ cats.f, data = tmp)
summary(mod); car::Anova(mod)
confint(lm(protected ~ cats.f - 1, data = tmp))
out = TukeyHSD(aov(mod))
plot(out)

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
summary(all.crit$prop.EOO.forest[all.crit$endemic %in% "endemic"])

boxplot(prop.EOO.forest ~ cat.reg.clean, 
        data = all.crit[all.crit$endemic %in% "endemic",], 
        notch = TRUE, varwidth = TRUE)
cats.f <- factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],
                 levels = c("LC","NT","VU","EN","CR"))
levels(cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")
boxplot(prop.EOO.forest ~ cats.f, 
        data = all.crit[all.crit$endemic %in% "endemic",], 
        notch = TRUE, varwidth = TRUE, col = c("yellowgreen","gold","darkorange","red"))
abline(h=c(25,40),lty=2)

tmp <- all.crit[all.crit$endemic %in% "endemic",]
mod <- lm(log(prop.EOO.forest+1) ~ cats.f, 
          data = tmp)
summary(mod); car::Anova(mod)
confint(lm(Median.HII ~ cats.f - 1, data = tmp))
out = TukeyHSD(aov(mod))
plot(out)

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

###################################################
#### THE INFLUENCE OF USING DIFFERENT CRITERIA ####
###################################################
##Species classified under all criteria vs. criterion A and criterion B
tmp <- all.crit[!is.na(all.crit$reduction_A12) &  # Assessed under criterion A
                  !is.na(all.crit$AOO) & # Assessed under criterion B
                  !all.crit$cat.reg.clean %in% "NA",]
dim(tmp)[1] ## 2757 species
100*dim(tmp)[1]/dim(all.crit)[1]
tmp$cat.reg.clean[tmp$cat.reg.clean %in% c("LC","NT")] <- "LC+NT"

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
ts0 <- round(100*sum(mat1[-3,]) / sum(mat1),2) #94.8% of threatened species
ts1 <- round(100*sum(mat1[,-3]) / sum(mat1),2) #94.5% of threatened species
ts2 <- round(100*sum(mat2[,-3]) / sum(mat2),2) #10.2% of threatened species
ts3 <- round(100*sum(mat3[,-3]) / sum(mat3),2) #3.2% of threatened species
ts4 <- round(100*sum(mat4[,-3]) / sum(mat4),2) #1.7% of threatened species

#Categories of threat
cg1 <- round(100*sum(diag(mat1))/sum(mat1),2) #95.7% of congruence in the assessments
cg2 <- round(100*sum(diag(mat2))/sum(mat2),2) #10.9% of congruence in the assessments
cg3 <- round(100*sum(diag(mat3))/sum(mat3),2) #5.6% of congruence in the assessments
cg4 <- round(100*sum(diag(mat4))/sum(mat4),2) #5.3% of congruence in the assessments
cg5 <- round(100*sum(diag(mat5))/sum(mat5),2) #9.3% of congruence in the assessments


## ENDEMIC SPECIES ##
tmp1 <- tmp[tmp$endemic %in% "endemic",]
dim(tmp1)[1] ## 1646 species
100*dim(tmp1)[1]/dim(all.crit)[1] # 32.3% of the total

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
ts0.1 <- round(100*sum(mat1.1[-3,]) / sum(mat1.1),2) #93.0% of threatened species
ts1.1 <- round(100*sum(mat1.1[,-3]) / sum(mat1.1),2) #92.2% of threatened species
ts2.1 <- round(100*sum(mat2.1[,-3]) / sum(mat2.1),2) #15.0% of threatened species
ts3.1 <- round(100*sum(mat3.1[,-3]) / sum(mat3.1),2) #2.6% of threatened species
ts4.1 <- round(100*sum(mat4.1[,-2]) / sum(mat4.1),2) #2.3% of threatened species

#Categories of threat
cg1.1 <- round(100*sum(diag(mat1.1))/sum(mat1.1),2) #96.9% of congruence in the assessments
cg2.1 <- round(100*sum(diag(mat2.1))/sum(mat2.1),2) #15.3% of congruence in the assessments
cg3.1 <- round(100*sum(diag(mat3.1))/sum(mat3.1),2) #7.4% of congruence in the assessments
cg4.1 <- round(100*sum(diag(mat4.1[c(1,3,4),c(1,2,3)]))/sum(mat4.1),2) #7.1% of congruence in the assessments
cg5.1 <- round(100*sum(diag(mat5.1))/sum(mat5.1),2) #9.3% of congruence in the assessments

## Table (S)1? ##
tabS1 <- matrix(c(ts1, paste0(rl1[2]," [",paste0(round(rl1[c(1,3)],3), collapse = "-"),"]"), ts1.1, paste0(rl1.1[2]," [",paste0(round(rl1.1[c(1,3)],3), collapse = "-"),"]"), #cg1.1,
         ts2, paste0(rl2[2]," [",paste0(round(rl2[c(1,3)],3), collapse = "-"),"]"), ts2.1, paste0(rl2.1[2]," [",paste0(round(rl2.1[c(1,3)],3), collapse = "-"),"]"), #cg2.1,
         ts3, paste0(rl3[2]," [",paste0(round(rl3[c(1,3)],3), collapse = "-"),"]"), ts3.1, paste0(rl3.1[2]," [",paste0(round(rl3.1[c(1,3)],3), collapse = "-"),"]"), #cg3.1,
         ts4, paste0(rl4[2]," [",paste0(round(rl4[c(1,3)],3), collapse = "-"),"]"), ts4.1, paste0(rl4.1[2]," [",paste0(round(rl4.1[c(1,3)],3), collapse = "-"),"]"),#, cg4.1), 
         ts0, paste0(rl0[2]," [",paste0(round(rl0[c(1,3)],3), collapse = "-"),"]"), ts0.1, paste0(rl0.1[2]," [",paste0(round(rl0.1[c(1,3)],3), collapse = "-"),"]")),
         nrow = 5, ncol = 4, byrow = TRUE,
       dimnames = list(c("catA","catB","catC","catD","all"), c("Threatened","RLI","Threatened","RLI")))
write.csv(tabS1,"TableS1.csv")

#### INFLUENCES OF DATA UNCERTAINTY ####


##################################################################################################################################################################################################################################H
##################################################################################################################################################################################################################################H
#################
#### FIGURES ####
#################
source("R/my.PieDonut.R")

## Colors for each category
cores <- c(CR_PE = "darkred", CR = "red", EN = "darkorange", VU = "gold", NT = "yellowgreen", `NA` = "grey", DD = "grey", LC = "forestgreen")
#LC = "green4"?; NT = "palegreen4"?

#### FIGURE 1 ####
#Creating and organizing the data frame to be plotted
pie.df <- all.crit[,c("cat.reg.clean","main.criteria","endemic")]
names(pie.df)[1] <- "category"
pie.df$category[pie.df$category %in% "CR_PE"] <- "CR" 
pie.df$main.criteria[pie.df$main.criteria %in% "A1+A2+B1+B2+C1+C2+D"] <- "all" 
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+B1+B2+C1","A1+A2+B2+C1")] <- "other" #"A+B+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+C1","A1+A2+C1+C2","A1+A2+C1")] <- "other" #"A+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+B2","A2+B1+B2","A2+B2","A1+A2+B1+B2")] <- "A2, B2"
pie.df$main.criteria[pie.df$main.criteria %in% c("B1")] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+C2+D")] <- "other" #"A+C+D"
pie.df$main.criteria[is.na(pie.df$main.criteria)] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2")] <- "A2"
pie.df$main.criteria[pie.df$category %in% "NA"] <- ""
pie.df <- pie.df[order(match(pie.df$category, names(cores)[-1])),]
pie.df <- pie.df[!pie.df$category %in% "DD",]
pie.df$main.criteria[pie.df$main.criteria %in% c("B1+B2")] <- "B1+2"

## Geral
dados <- pie.df
toto <- rep(TRUE, length(table(dados$category, dados$main.criteria)))
toto[c(10,20,32)] <- FALSE
pie.all <- my.PieDonut(dados, aes(category, main.criteria),
                   ratioByGroup=FALSE, showPieName = FALSE,
                   start=75,r0=0, r1=1,r2=1.3,
                   showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
                   main.colors = cores[c(2:5,8,6)],
                   #donut.colors = NULL,
                   pieAlpha = 0.75,
                   correct.legend.x = c(0.00,0.00,-0.075,-0.100,-0.15,-0.10),
                   correct.legend.y = c(0.15,0.00,-0.200,-0.075, 0.00, 0.10),
                   tidy.legend.donut = toto,
                   tidy.legend.pie = c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE),
                   plot.piedonut = FALSE,
                   return = "pie", pieLabelSize = 4.5, donutLabelSize = 3)
donut.all <- my.PieDonut(dados, aes(category, main.criteria),
                       ratioByGroup=FALSE, showPieName = FALSE,
                       start=75,r0=0, r1=1,r2=1.3,
                       showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
                       main.colors = cores[c(2:5,8,6)],
                       #donut.colors = NULL,
                       pieAlpha = 0.75,
                       correct.legend.x = c(0.00,0.00,-0.075,-0.100,-0.15,-0.10),
                       correct.legend.y = c(0.15,0.00,-0.200,-0.075, 0.00, 0.10),
                       tidy.legend.donut = toto,
                       tidy.legend.pie = c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE),
                       plot.piedonut = FALSE,
                       return = "donut", pieLabelSize = 4.5, donutLabelSize = 3)

## Endemicas
dados <- pie.df[pie.df$endemic %in% "endemic",]
toto <- rep(TRUE, length(table(dados$category, dados$main.criteria)))
toto[c(16)] <- FALSE
pie.end <- my.PieDonut(dados, aes(category, main.criteria),
         ratioByGroup=FALSE, showPieName = FALSE,
         start=75,r0=0, r1=1,r2=1.3,
         showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
         #labelposition = 2, selected=c(1:15), title="labelposition=1",
         #explode = c(1:3), explodeDonut = TRUE,
         main.colors = cores[c(2:5,8)],
         #donut.colors = NULL,
         correct.legend.x = c(0,0,-0.15,-0.1,-0.1),
         correct.legend.y = c(0.15,0,-0.05,0.05,0.1),
         tidy.legend.donut = toto,
         tidy.legend.pie = c(TRUE,TRUE,TRUE,FALSE,TRUE),
         plot.piedonut = TRUE,
         return = "pie", pieLabelSize = 4.5, donutLabelSize = 3)
donut.end <- my.PieDonut(dados, aes(category, main.criteria),
          ratioByGroup=FALSE, showPieName = FALSE,
          start=75,r0=0, r1=1,r2=1.3,
          showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
          #labelposition = 2, selected=c(1:15), title="labelposition=1",
          #explode = c(1:3), explodeDonut = TRUE,
          main.colors = cores[c(2:5,8)],
          #donut.colors = NULL,
          correct.legend.x = c(0,0,-0.15,-0.1,-0.1),
          correct.legend.y = c(0.15,0,-0.05,0.05,0.1),
          tidy.legend.donut = toto,
          tidy.legend.pie = c(TRUE,TRUE,TRUE,FALSE,TRUE),
          plot.piedonut = TRUE,
          return = "donut", pieLabelSize = 4.5, donutLabelSize = 3)

## Construting the plots
pie.donut.all <- ggdraw(pie.all) + draw_plot(donut.all)
pie.donut.end <- ggdraw(pie.end) + draw_plot(donut.end)

fig1 <- plot_grid(pie.donut.all, pie.donut.end,
  labels = c('A - All populations', 'B - Endemics'),
  align="v", vjust = 6
)

path <- "C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo Extincao na MA/data analysis/figures"
ggsave2("Figure1.jpg", fig1, "jpeg", path,
        width = 50, height = 20, units = "cm", dpi = 320)


#### FIGURE 2 ####
# Previous vs. new assessments
require(circlize)

jpeg(filename = "figures/Figure2.jpg", width = 3750, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(1,2))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)

## NEW VS. PREVIOUS (GLOBAL)
mat <- as.matrix(table(paste0(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],"_new"), 
					paste0(all.crit$redlistCategory[all.crit$endemic %in% "endemic"],"_prev")))
mat <- mat[c(1,3,6,5,4,2), rev(c(5,4,1,3,9,8,6,2))]

#Defining the colors of tracks and links
grid.col = c(EX_prev = "black", EW_prev = "purple", CR_prev = "red", EN_prev = "darkorange", VU_prev = "gold", NT_prev = "yellowgreen", LC_prev = "forestgreen", DD_prev = "grey",
		   CR_new = "red", EN_new = "darkorange", VU_new = "gold", NT_new = "yellowgreen", LC_new = "forestgreen", DD_new = "grey")
col_mat = rep(rev(c("black", "purple", "red", "darkorange", "gold", "yellowgreen", "forestgreen", "grey")), each=6)
#col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
#col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 10 & mat >= 5] = mat[mat < 10 & mat >= 5] + 1
mat[mat > 0 & mat < 5] = mat[mat > 0 & mat < 5] + 2
transp <- rep(0.3, length(col_mat))
transp[col_mat %in% "forestgreen"] <- 0.6

#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible[,-c(1:2)]) = FALSE
chordDiagram(mat, big.gap = 10, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
             grid.col = grid.col, col = col_mat,
		         transparency = transp ,
             #self.link = 1, #link.visible = visible,
             #h=1.1, #w=0.5,
             #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
             link.lwd = 4,
             h.ratio = 0.9,
             #reduce_to_mid_line = TRUE,
             w2=0.5, rou=0.2
             #point1 = rep(0,16)
)
#Putting legends on
sec.ind <- c("CR","EN","VU","NT","LC","DD","DD","LC","NT","VU","EN","CR")#,"EW","EX")
for(si in get.all.sector.index()) {
  lab <- sec.ind[which(si == get.all.sector.index())]
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#  if(si == "VU_hi") {
#    circos.text(mean(xlim), mean(ylim), labels = "VU", 
#                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
#                facing = "bending", niceFacing = FALSE, col = "black")
#    
#  } else {
    circos.text(mean(xlim), mean(ylim), labels = lab, 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = TRUE, col = "black")
#  }  
}
legend("topleft","Previous", bty="n", cex=1.3, adj=c(-.5,2.5))
legend("topright","New assess.", bty="n", cex=1.3, adj=c(0.25,2.5))
legend("topleft",legend=expression(bold("A - Global")),
	bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.3)
text(-0.205,1.02,"EW", cex=0.9)
text(-0.105,1.04,"EX", cex=0.9)

## NEW VS. PREVIOUS (NATIONAL)
all.crit$redlistCategory.national <- all.crit$status.reflora
all.crit$redlistCategory.national[is.na(all.crit$redlistCategory.national) & !is.na(all.crit$category.PAY)] <-
	all.crit$category.PAY[is.na(all.crit$redlistCategory.national) & !is.na(all.crit$category.PAY)]
all.crit$redlistCategory.national[is.na(all.crit$redlistCategory.national) & !is.na(all.crit$category.ARG)] <-
	all.crit$category.ARG[is.na(all.crit$redlistCategory.national) & !is.na(all.crit$category.ARG)]
all.crit$redlistCategory.national[all.crit$redlistCategory.national %in% "Indeterminada"] <- "DD"
all.crit$redlistCategory.national[all.crit$redlistCategory.national %in% "Decreasing in the Cordoba province"] <- "NT"
all.crit$redlistCategory.national <- sapply(strsplit(all.crit$redlistCategory.national," \\("), function(x) x[1])

mat <- as.matrix(table(paste0(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],"_new"), 
					paste0(all.crit$redlistCategory.national[all.crit$endemic %in% "endemic"],"_prev")))
mat <- mat[c(1,3,6,5,4,2), c(2,4,6,7,3,1)]

#Defining the colors of tracks and links
grid.col = c(CR_prev = "red", EN_prev = "darkorange", VU_prev = "gold", NT_prev = "yellowgreen", LC_prev = "forestgreen", DD_prev = "grey",
		   CR_new = "red", EN_new = "darkorange", VU_new = "gold", NT_new = "yellowgreen", LC_new = "forestgreen", DD_new = "grey")
col_mat = rep(rev(c("red", "darkorange", "gold", "yellowgreen", "forestgreen", "grey")), each=5)
# mat[mat < 15] = mat[mat < 15]*1.25
# mat[mat > 0 & mat < 5] = 5
mat[mat < 10 & mat >= 5] = mat[mat < 10 & mat >= 5] + 1
mat[mat > 0 & mat < 5] = mat[mat > 0 & mat < 5] + 2
transp <- rep(0.3, length(col_mat))
transp[col_mat %in% "forestgreen"] <- 0.6

#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
chordDiagram(mat, big.gap = 10, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
             grid.col = grid.col, col = col_mat,
		    transparency = transp ,
             link.lwd = 4,
             h.ratio = 0.9,
             w2=0.5, rou=0.2
)
#Putting legends on
sec.ind <- c("CR","EN","VU","NT","LC","DD","LC","NT","VU","EN","CR")
for(si in get.all.sector.index()) {
  lab <- sec.ind[which(si == get.all.sector.index())]
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), labels = lab, 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = TRUE, col = "black")
}
legend("topleft","Previous", bty="n", cex=1.3, adj=c(-.5,2.5))
legend("topright","New assess.", bty="n", cex=1.3, adj=c(0.25,2.5))
legend("topleft",legend=expression(bold("B - Regional")),
	bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.2)
dev.off()


#### FIGURE 3 ####
#Proportion of protected areas

jpeg(filename = "figures/Figure3.jpg", width = 2500, height = 2250, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

cexs <- pchs <- cores1 <- as.factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])
levels(cores1) <- c("red", "grey", "darkorange", "forestgreen", "yellowgreen", "gold")
levels(pchs) <- c(19,1,15,1,1,17)
levels(cexs) <- c(1,1,1,1,1,1.2)

layout.matrix <- matrix(c(2, 1, 0, 3), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(0.5, 2), # Heights of the two rows
       widths = c(2, 0.5)) # Widths of the two columns
# Plot 1: Scatterplot
par(mar=c(3,3,0.5,0.5), mgp=c(1.5,0.25,0),tcl=-0.2,las=1)
plot(log(all.crit$protected[all.crit$endemic %in% "endemic"] + 1) ~ 
       log(all.crit$prop.EOO.in.StrictUCs[all.crit$endemic %in% "endemic"] + 1),
     xlab = "Proportion of the EOO", ylab = "Proportion of the occurrences",
     cex.lab  =1.2,
     xaxt = "n", yaxt= "n",
     cex = as.double(as.character(cexs)),
     pch = as.double(as.character(pchs)),
     col = alpha(as.character(cores1), alpha = 0.6),
     xlim = log(c(0,100)+1), ylim = log(c(0,100)+1)
)
axis(1, at=log(c(1,5,10,25,50)+1), labels = c(1,5,10,25,50))
axis(2, at=log(c(1,5,10,25,50)+1), labels = c(1,5,10,25,50))
#abline(0,1,lty=2)
linhas <- log(c(10,25,50)+1)
# abline(v = linhas, lty=3)
# abline(h = linhas, lty=3)
abline(v = linhas[1], lty=3)
abline(h = linhas[2], lty=3)

legend(log(45),log(4), c("CR","EN","VU","LC+NT"),
       pch=c(19,15,17,1), col=c("red", "darkorange",  "gold", "forestgreen"), 
       bty= "n", cex=1.1)

# Plot 2: Top (height) boxplot
cats.f <- factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],
                 levels = c("LC","NT","VU","EN","CR"))
levels(cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")
par(mar = c(0, 3, 0, 0.5), mgp=c(1,0,0))
pos <- c(0.2,0.75,1.5,2.2)
a <- boxplot(log(prop.EOO.in.StrictUCs+1) ~ cats.f,
        data = all.crit[all.crit$endemic %in% "endemic",], 
        notch = TRUE, varwidth = TRUE, frame = FALSE, horizontal = TRUE,
        at = pos, xaxt = "n", yaxt = "n", 
        ylim = log(c(0,100)+1), xlim=c(0,2.5),
        col = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4))
axis(2, pos, col = NA, col.ticks = NA, 
     labels = c("LC+NT","VU","EN","CR"), cex.axis = 0.7)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
model <- lm(log(prop.EOO.in.StrictUCs+1) ~ cats.f,
            data = all.crit[all.crit$endemic %in% "endemic",])
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")

#Add the legend
par(xpd=T)
text( log(120), pos, rev(LABELS[,1])  , col=1 )
par(xpd=F)

# Plot 3: Right (weight) boxplot
par(mar = c(3, 0, 0.5, 0))
a <- boxplot(log(protected+1) ~ cats.f,
        data = all.crit[all.crit$endemic %in% "endemic",], 
        notch = TRUE, varwidth = TRUE, frame = FALSE, horizontal = FALSE,
        at = c(0.2,0.75,1.5,2.2), yaxt = "n", xaxt = "n", 
        ylim = log(c(0,100)+1), xlim=c(0,2.5),
        col = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4))
axis(1, pos, col = NA, col.ticks = NA, 
     labels = c("LC+NT","VU","EN","CR"),  cex.axis = 0.7)

model <- lm(log(protected+1) ~ cats.f,
            data = all.crit[all.crit$endemic %in% "endemic",])
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")

#Add the legend
par(xpd=T)
text( pos, log(120),  rev(LABELS[,1])  , col=1 )
par(xpd=F)

dev.off()

#### FIGURE X ####

par(mar = c(3, 3, 0.5, 0.5))
a <- boxplot(log10(pop.size) ~ cats.f,
             data = all.crit[all.crit$endemic %in% "endemic",], 
             notch = TRUE, varwidth = TRUE, # frame = FALSE, horizontal = FALSE,
             #at = c(0.2,0.75,1.5,2.2), yaxt = "n", xaxt = "n", 
             #ylim = log(c(0,100)+1), xlim=c(0,2.5),
             col = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4))
# axis(1, 1:4, col = NA, col.ticks = NA, 
#      labels = c("LC+NT","VU","EN","CR"),  cex.axis = 0.7)

model <- lm(log(pop.size) ~ cats.f,
            data = all.crit[all.crit$endemic %in% "endemic",])
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")

#Add the legend
par(xpd=T)
text(1:4, max(log10(all.crit$pop.size)+0.2, na.rm = TRUE), rev(LABELS[,1])  , col=1 )
par(xpd=F)
