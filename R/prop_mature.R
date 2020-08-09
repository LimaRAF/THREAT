####################################################################################
#### GETTING MEAN VALUES OF PROP. OF MATURE INDIVIDUALS FROM TREECO RAW DATASET ####
####################################################################################


##Loading TreeCo individual-tree data
rawdata <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Databases//Species abundances//dados_brutos_treeco.csv", na.strings=c(""," ","NA"),as.is=T)

#Removing non Atlantic Forest sites
rawdata <- rawdata[!rawdata$SiteCode %in% c("SPeea1", "SPeea2"),] 

#Removing individuals not identified to the species level
rawdata <- rawdata[grepl(" ", rawdata$Latin),]
rawdata <- rawdata[!grepl(" sp\\.", rawdata$Latin),]
rawdata <- rawdata[!grepl(" aff\\.", rawdata$Latin),]
rawdata$Latin <- gsub(" cf\\. "," ", rawdata$Latin)
rawdata$Latin <- sapply(strsplit(rawdata$Latin," "), function(x) paste(x[1], x[2],sep=" "))

#Fixing some main identification problems
rawdata$Latin <- stringr::str_replace_all(rawdata$Latin,
                                         c("Xylopia langsdorfiana" = "Xylopia langsdorffiana",  "Marlierea tomentosa" = "Myrcia strigipes",
                                           "Maytenus gonoclada" = "Monteverdia gonoclada", "Eugenia cuprea" = "Eugenia expansa",
                                           "Ocotea elegans" = "Ocotea indecora", "Quiina glazovii" = "Quiina glaziovii",
                                           "Marlierea racemosa" = "Myrcia vellozoi", "Eugenia rostrata" = "Eugenia zuccarinii",
                                           "Marlierea eugeniopsoides" = "Myrcia eugeniopsoides", "Myrcia pulchra" = "Myrcia subcordata",
                                           "Almeidea coerulea" = "Conchocarpus coeruleus", "Erythroxylum amplifolium" = "Erythroxylum umbu",
                                           "Marlierea reitzii" = "Myrcia reitzii", "Quiina magallano-gomesii" = "Quiina glaziovii",
                                           "Tetragastris catuaba" = "Protium catuaba", "Laplacea fructicosa" = "Laplacea fruticosa",
                                           "Maytenus distichophylla" = "Monteverdia distichophylla", "Maytenus obtusifolia" = "Monteverdia obtusifolia",
                                           "Marlierea obversa" = "Myrcia obversa", "Maytenus aquifolia" = "Monteverdia aquifolia",
                                           "Marlierea silvatica" = "Myrcia ferruginosa", "Tibouchina stenocarpa" = "Pleroma stenocarpum",
                                           "Eugenia bocainensis" = "Eugenia puberula", "Henriettea saldanhaei" = "Henriettea saldanhae",
                                           "Sloanea stipitata" = "Sloanea guianensis", "Tibouchina trichopoda" = "Pleroma trichopodum",
                                           "Bathysa cuspidata" = "Schizocalyx cuspidatus", "Tibouchina elegans" = "Pleroma elegans",
                                           "Simaba cedron" = "Homalolepis cedron", "Maytenus brasiliensis" = "Monteverdia brasiliensis",
                                           "Marlierea glazioviana" = "Myrcia insigniflora", "Faramea glaziovii" = "Faramea martiana",
                                           "Tovomita excelsa" = "Tovomita choisyana", "Caesalpinia echinata" = "Paubrasilia echinata",
                                           "Marlierea glabra" = "Myrcia neoglabra", "Andradaea floribunda" = "Andradea floribunda",
                                           "Calyptranthes grandifolia" = "Calyptranthes brasiliensis", "Crepidospermum atlanticum" = "Protium atlanticum",
                                           "Mollinedia micrantha" = "Mollinedia elegans", "Pouteria hispida" = "Pouteria guianensis",
                                           "Psychotria myriantha" = "Palicourea mamillaris", "Sloanea usurpatrix" = "Sloanea sinemariensis",
                                           "Tibouchina fissinervia" = "Pleroma fissinervium", "Davilla kunthii" = "Davilla nitida",
                                           "Eugenia convexinervia" = "Eugenia supraaxillaris", "Ixora grandifolia" = "Ixora muelleri",
                                           "Lonchocarpus campestris" = "Muellera campestris", "Margaritopsis cephalantha" = "Eumachia cephalantha",
                                           "Pradosia verrucosa" = "Pradosia glaziovii", "Psychotria racemosa" = "Palicourea racemosa",
                                           "Psychotria stellaris" = "Palicourea octocuspis", "Sebastiania serrata" = "Gymnanthes serrata",
                                           "Simaba floribunda" = "Homalolepis floribunda", "Tovomita bahiensis" = "Tovomita choisyana",
                                           "Psidium cattleianum" = "Psidium cattleyanum", "Maytenus ilicifolia" = "Monteverdia ilicifolia",
                                           "Mollinedia marqueteana" = "Mollinedia lamprophylla", "Maytenus ardisiaefolia" = "Monteverdia ardisiifolia",
                                           "Xylosma glaberrimum" = "Xylosma glaberrima", "Xylosma tweedianum" = "Xylosma tweediana",
                                           "Cinnamomum hirsutum" = "Aiouea hirsuta", "Miconia theizans" = "Miconia theaezans",
                                           "Ocotea vaccinioides" = "Ocotea daphnifolia", "Swartzia alternatifolia" = "Swartzia alternifoliolata"))

#Getting the dbh
rawdata$dbh <- rawdata$PAP/pi

#Removing DBH<5 cm
rawdata <- rawdata[rawdata$PAP>=15,]

#Keeping only the largest dbh for multiple-stemmed trees
rawdata <- rawdata[order(rawdata$ordem),]
rawdata <- rawdata[order(rawdata$stemID),]
rawdata <- rawdata[order(rawdata$dbh, decreasing = TRUE),]
rawdata <- rawdata[!duplicated(rawdata$ordem),]

#Getting the number of individuals per species per dbh interval
length(unique(rawdata$SiteCode)) # how many surveys?
dim(rawdata)[1] # how many tree measurements?

result <- aggregate(rawdata$ordem, list(rawdata$Latin), length)
names(result) <- c("species.correct", "N")
result$dbh.6 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=6))$x/result$N
result$dbh.7 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=7))$x/result$N
result$dbh.7.5 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=7.5))$x/result$N
result$dbh.8 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=8))$x/result$N
result$dbh.9 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=9))$x/result$N
result$dbh.10 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=10))$x/result$N
result$dbh.12.5 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=12.5))$x/result$N
result$dbh.15 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=15))$x/result$N
result$dbh.20 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=20))$x/result$N
#result$dbh.25 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=25))$x/result$N
#result$dbh.30 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=30))$x/result$N
#result$dbh.40 <- aggregate(rawdata$dbh, list(rawdata$Latin), function(x) sum(x>=40))$x/result$N

##Merging with MaxHeight, GRowth form and ecological group data
result1 <- merge(result, hab[,c("species.correct","habito","life.form","MaxHeight","GF","wsg","ecol.group")]
                 , by = "species.correct", all.x = TRUE)

## Summary statistics
result2 <- result1[result1$N>15,]

#Overall p for all species (i.e. unknown)
apply(result1[!result1$GF %in% "large_shrub", 3:11], 2 , mean, na.rm = TRUE)
apply(result1[!result1$GF %in% "large_shrub", 3:11], 2 , function(x) confint(lm(x~1)))

apply(result2[!result2$GF %in% "large_shrub", 3:11], 2 , mean, na.rm = TRUE)
apply(result2[!result2$GF %in% "large_shrub", 3:11], 2 , function(x) confint(lm(x~1)))

apply(result1[result1$GF %in% "small_tree", 3:11], 2 , mean, na.rm = TRUE)
apply(result1[result1$GF %in% "small_tree", 3:11], 2 , function(x) confint(lm(x~1)))

apply(result1[result2$GF %in% "large_tree", 3:11], 2 , mean, na.rm = TRUE)
apply(result1[result2$GF %in% "large_tree", 3:11], 2 , function(x) confint(lm(x~1)))


#p for each GF and each EG
aggregate(result1[,3:11], list(result1$GF), mean, na.rm = TRUE)
aggregate(result1[,3:11], list(result1$ecol.group), mean, na.rm = TRUE)


#Overall p for each GF and EG combination
af.p <- af.p1 <- matrix(NA, ncol = 18, nrow = 8)
rownames(af.p) <- rownames(af.p1) <- rep(c("pioneer","early","late","climax"), each = 2)
colnames(af.p) <- colnames(af.p1) <- rep(c("dbh.6","dbh.7","dbh.7.5","dbh.8","dbh.9","dbh.10","dbh.12.5","dbh.15","dbh.20"), 2)
egs <- rep(c("pioneer","early_secondary","late_secondary","climax"), each = 2)
gfs <- c("small_tree", "large_tree")

for(i in c(1,3,5,7)) {
  eg <- egs[i]
  data1 <- result1[result1$ecol.group %in% eg,]
  for(j in 1:2) {
    gf <- gfs[j]
    data2 <- data1[data1$GF %in% gf,]
    if(j==1) {
      af.p[i,1:9] <- round(apply(data2[,3:11], 2 , mean, na.rm = TRUE), 4) 
      af.p[i+1,1:9] <- apply(round(apply(data2[,3:11], 2 , function(x) confint(lm(x~1))), 4), 2, paste, collapse = "-")  
    }
    if(j==2) {
      af.p[i,10:18] <- round(apply(data2[,3:11], 2 , mean, na.rm = TRUE), 4) 
      af.p[i+1,10:18] <- apply(round(apply(data2[,3:11], 2 , function(x) confint(lm(x~1))), 4), 2, paste, collapse = "-")  
    }
  }
}
af.p <- af.p[,c(1,2,4,5,3,15:18)]

for(i in c(1,3,5,7)) {
  eg <- egs[i]
  data1 <- result2[result2$ecol.group %in% eg,]
  for(j in 1:2) {
    gf <- gfs[j]
    data2 <- data1[data1$GF %in% gf,]
    if(j==1) {
      af.p1[i,1:9] <- round(apply(data2[,3:11], 2 , mean, na.rm = TRUE), 4) 
      af.p1[i+1,1:9] <- apply(round(apply(data2[,3:11], 2 , function(x) confint(lm(x~1))), 4), 2, paste, collapse = "-")  
    }
    if(j==2) {
      af.p1[i,10:18] <- round(apply(data2[,3:11], 2 , mean, na.rm = TRUE), 4) 
      af.p1[i+1,10:18] <- apply(round(apply(data2[,3:11], 2 , function(x) confint(lm(x~1))), 4), 2, paste, collapse = "-")  
    }
  }
}
af.p1 <- af.p1[,c(1,2,4,5,3,15:18)]

#Selecting the info for each Dcrit
res <- res1 <- matrix(NA, ncol = 4, nrow = 2)
colnames(res) <- colnames(res1) <- c("pioneer","early","late","climax")
rownames(res) <- rownames(res1) <- c("small_tree", "large_tree")

res[1,1] <- paste0(af.p[1,1]," [", af.p[2,1],"]"); res1[1,1] <- paste0(af.p1[1,1]," [", af.p1[2,1],"]")
res[1,2] <- paste0(af.p[3,2]," [", af.p[4,2],"]"); res1[1,2] <- paste0(af.p1[3,2]," [", af.p1[4,2],"]")
res[1,3] <- paste0(af.p[5,3]," [", af.p[6,3],"]"); res1[1,3] <- paste0(af.p1[5,3]," [", af.p1[6,3],"]")
res[1,4] <- paste0(af.p[7,4]," [", af.p[8,4],"]"); res1[1,4] <- paste0(af.p1[7,4]," [", af.p1[8,4],"]")

res[2,1] <- paste0(af.p[1,6]," [", af.p[2,6],"]"); res1[2,1] <- paste0(af.p1[1,6]," [", af.p1[2,6],"]")
res[2,2] <- paste0(af.p[3,7]," [", af.p[4,7],"]"); res1[2,2] <- paste0(af.p1[3,7]," [", af.p1[4,7],"]")
res[2,3] <- paste0(af.p[5,8]," [", af.p[6,8],"]"); res1[2,3] <- paste0(af.p1[5,8]," [", af.p1[6,8],"]")
res[2,4] <- paste0(af.p[7,9]," [", af.p[8,9],"]"); res1[2,4] <- paste0(af.p1[7,9]," [", af.p1[8,9],"]")

res # usando todas as espécies para o artigo...
res1


## Getting prop. of mature for each group
result1$p <- NA
result1$dbh.crit <- NA
result1$p[result1$GF %in% "large_shrub"] <- 1
result1$p[result1$GF %in% "small_tree" & result1$ecol.group %in% "pioneer"] <- 
  result1$dbh.6[result1$GF %in% "small_tree" & result1$ecol.group %in% "pioneer"]
result1$p[result1$GF %in% "small_tree" & result1$ecol.group %in% "early_secondary"] <- 
  result1$dbh.7[result1$GF %in% "small_tree" & result1$ecol.group %in% "early_secondary"]
result1$p[result1$GF %in% "small_tree" & result1$ecol.group %in% "late_secondary"] <- 
  result1$dbh.8[result1$GF %in% "small_tree" & result1$ecol.group %in% "late_secondary"]
result1$p[result1$GF %in% "small_tree" & result1$ecol.group %in% "climax"] <- 
  result1$dbh.9[result1$GF %in% "small_tree" & result1$ecol.group %in% "climax"]

result1$p[result1$GF %in% "large_tree" & result1$ecol.group %in% "pioneer"] <- 
  result1$dbh.10[result1$GF %in% "large_tree" & result1$ecol.group %in% "pioneer"]
result1$p[result1$GF %in% "large_tree" & result1$ecol.group %in% "early_secondary"] <- 
  result1$dbh.12.5[result1$GF %in% "large_tree" & result1$ecol.group %in% "early_secondary"]
result1$p[result1$GF %in% "large_tree" & result1$ecol.group %in% "late_secondary"] <- 
  result1$dbh.15[result1$GF %in% "large_tree" & result1$ecol.group %in% "late_secondary"]
result1$p[result1$GF %in% "large_tree" & result1$ecol.group %in% "climax"] <- 
  result1$dbh.20[result1$GF %in% "large_tree" & result1$ecol.group %in% "climax"]

result1$dbh.crit[result1$GF %in% "large_shrub"] <- 5
result1$dbh.crit[result1$GF %in% "small_tree" & result1$ecol.group %in% "pioneer"] <- 6
result1$dbh.crit[result1$GF %in% "small_tree" & result1$ecol.group %in% "early_secondary"] <- 7
result1$dbh.crit[result1$GF %in% "small_tree" & result1$ecol.group %in% "late_secondary"] <- 8
result1$dbh.crit[result1$GF %in% "small_tree" & result1$ecol.group %in% "climax"] <- 9

result1$dbh.crit[result1$GF %in% "large_tree" & result1$ecol.group %in% "pioneer"] <- 10
result1$dbh.crit[result1$GF %in% "large_tree" & result1$ecol.group %in% "early_secondary"] <- 12.5
result1$dbh.crit[result1$GF %in% "large_tree" & result1$ecol.group %in% "late_secondary"] <- 15
result1$dbh.crit[result1$GF %in% "large_tree" & result1$ecol.group %in% "climax"] <- 20

result2 <- result1[result1$N>15,]
summary(result1$p); summary(result2$p)
boxplot(result1$p ~ result1$GF + result1$ecol.group, notch = TRUE, varwidth = TRUE)
boxplot(result2$p ~ result2$GF + result2$ecol.group, notch = TRUE, varwidth = TRUE)

##Getting average proportions per groups
tmp <- lme4::lmer(p ~ MaxHeight + (1|ecol.group), data = result2)
tmp1 <- lme4::lmer(p ~ MaxHeight + I(MaxHeight^2) + (1|ecol.group), data = result2)
tmp2 <- lm(p ~ MaxHeight + wsg, data = result2)
tmp3 <- lme4::lmer(p ~ MaxHeight + wsg + (1|ecol.group), data = result2)
tmp4 <- lme4::lmer(p ~ MaxHeight + I(MaxHeight^2) + wsg + (1|ecol.group), data = result2)
bbmle::AICtab(tmp, tmp1, tmp2, tmp3, tmp4)
lme4::fixef(tmp4)
lme4::ranef(tmp4)
car::Anova(tmp4)

##Resultados meio estranho... (não usando pra nada no final das contas)
(lme4::fixef(tmp4)[1] - lme4::ranef(tmp4)[[1]][1,1]) + 
  lme4::fixef(tmp4)[2]*c(7.5,10,12.5,15,20) +
  lme4::fixef(tmp4)[3]*(c(7.5,10,12.5,15,20)^2) +
  lme4::fixef(tmp4)[4]*0.75 # climax
(lme4::fixef(tmp4)[1] - lme4::ranef(tmp4)[[1]][3,1]) + 
  lme4::fixef(tmp4)[2]*c(7.5,10,12.5,15,20) +
  lme4::fixef(tmp4)[3]*(c(7.5,10,12.5,15,20)^2) +
  lme4::fixef(tmp4)[4]*0.65 # late
(lme4::fixef(tmp4)[1] - lme4::ranef(tmp4)[[1]][2,1]) + 
  lme4::fixef(tmp4)[2]*c(7.5,10,12.5,15,20) +
  lme4::fixef(tmp4)[3]*(c(7.5,10,12.5,15,20)^2) +
  lme4::fixef(tmp4)[4]*0.6 # early
(lme4::fixef(tmp4)[1] - lme4::ranef(tmp4)[[1]][4,1]) + 
  lme4::fixef(tmp4)[2]*c(7.5,10,12.5,15,20) +
  lme4::fixef(tmp4)[3]*(c(7.5,10,12.5,15,20)^2) +
  lme4::fixef(tmp4)[4]*0.5 # pioneer


## RESULTS AVERAGE FOR THE MORE ABUNDANT SPECIES
ids <- !result1$GF %in% "large_shrub"
aggregate(result1$p[ids], list(result1$GF[ids], result1$ecol.group[ids]), mean, na.rm=TRUE)
aggregate(result1$p[ids], list(result1$GF[ids]), mean, na.rm=TRUE)
aggregate(result1$p[ids], list(result1$ecol.group[ids]), mean, na.rm=TRUE)


############################
#### Merging and Saving ####
############################
g1 <- sort(paste(c(1,3,2,4), sort(unique(result1$GF)),sep="."))
g2 <- sort(paste(c(4,2,3,1,5), sort(unique(result1$ecol.group)),sep = "."))
combo <- expand.grid(EG = g2, GF = g1)
combo[,1] <- gsub('[0-9:]\\.',"", combo[,1])
combo[,2] <- gsub('[0-9:]\\.',"", combo[,2])

#Prop. mature
combo$p.est <- c(1, 1, 1, 1, 1, # for shrubs
             af.p[1,1], af.p[3,2], af.p[5,3], af.p[7,4], 0.5122, # for small trees
             af.p[1,6], af.p[3,7], af.p[5,8], af.p[7,9], 0.3342, # for large trees
             0.6371, 0.4529, 0.3476, 0.2786, 0.4529) # for trees unknown
combo$str.match <- paste(combo$EG,combo$GF,sep = "_") 
result1$str.match <- paste(result1$ecol.group,result1$GF,sep = "_")
result.final <- dplyr::left_join(result1, combo[,c("str.match","p.est")], by= "str.match") 
result.final <- result.final[,c("species.correct","N","dbh.crit","p","p.est","str.match")]

## Some species with known p
result.final[result.final$species.correct %in% "Euterpe edulis",] # field p: 0.42
result.final$p[result.final$species.correct %in% "Euterpe edulis"] <- 0.42  # field p: 0.42
result.final[result.final$species.correct %in% "Hyeronima alchorneoides",] # field p in Panama: 0.61
#result.final[result.final$species.correct %in% "Araucaria angustifolia",] # dbh.crit == 20
result.final[result.final$species.correct %in% "Tapirira guianensis",] # dbh.crit == 10
result.final$p[result.final$species.correct %in% "Tapirira guianensis"] <- 0.5842
result.final[result.final$species.correct %in% "Cecropia glaziovii",] # dbh.crit == 10
plot(result.final$p.est[result.final$N>30] ~ result.final$p[result.final$N>30])
abline(lm(result.final$p.est[result.final$N>30] ~ result.final$p[result.final$N>30]), lwd=2, col=2)

## Saving...
write.csv(result.final,"data/prop_mature.csv", row.names = FALSE)
