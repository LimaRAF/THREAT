#######################################################
#### ASSESSING VERY SMALL POPULATIONS - CRITERIA D ####
#######################################################
rm(list=ls())

#### LOADING PACKAGES ###
#devtools::install_github("gdauby/ConR", ref = "master", force = TRUE) # old version
#devtools::install_github("gdauby/ConR@devel") # new version on GitHub
# detach("package:ConR", unload=TRUE)
# install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
#  repos = NULL, 
#  type = "source")
library("ConR")
library("red")
library("circlize")

#### LOADING THREAT POPULATION SIZE DATA (TREECO) ###
#Already with pop. sizes estimated for all necessary years
res.means <- readRDS("data/threat_mean_pop_sizes_infer.rds")
mean.pop.sizes <- readRDS("data/threat_mean_pop_sizes_for_ConR.rds")
low.pop.sizes <- readRDS("data/threat_low_pop_sizes_for_ConR.rds")
high.pop.sizes <- readRDS("data/threat_high_pop_sizes_for_ConR.rds") 

#Putting data in the ConR format
spp <- names(res.means)
spp <- gsub("_"," ", spp)
spp <- gsub("\\.","-", spp)
nrows <- length(spp)
decline.models <- matrix(NA, ncol = 2, nrow = nrows,
                         dimnames = list(spp, c("Before 1992", "After 1992")))
for(x in 1:length(res.means)) {
  decline.models[x, 1] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year<1990])
  decline.models[x, 2] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year>1995])
}

#### LOADING THREAT HABITAT AND ECOLOGY DATA ####
## Includes species info on Generation Length and Proportion of matrue individuals
hab <- read.csv("data/threat_habitats.csv", as.is = TRUE)
PopData <- merge(decline.models, hab, by.x= "row.names", by.y = "Name_submitted", all.x = TRUE)
PopData <- PopData[!is.na(PopData$internal_taxon_id),]
names(PopData)[1] <- "species"
PopData <- PopData[order(PopData$species),]

## Filtering the populational datasets
ids <- mean.pop.sizes$species %in% PopData$species
mean.pop.sizes <- mean.pop.sizes[ids,]
low.pop.sizes <- low.pop.sizes[ids,]
high.pop.sizes <- high.pop.sizes[ids,]
table(mean.pop.sizes$species == PopData$species)
rm(res.means)

#### LOADING THREAT EXPLOITED SPECIES INFORMATION  ####
explo <- read.csv("data/threat_exploited_timber_spp.csv", as.is = TRUE, encoding = "UTF-8")
explo$commercial <- grepl("interesse comercial|construction", explo$obs) |
  grepl("International Timber Trade", explo$sources) |
  grepl("Especies nativas para fins produtivos", explo$sources)
PopData <- merge(PopData, explo[,c("species.correct2", "times.cites","commercial")], 
                 by.x= "species", by.y = "species.correct2", all.x = TRUE)
PopData <- PopData[order(PopData$species),]
table(mean.pop.sizes$species == PopData$species)
PopData$timber <- !is.na(PopData$times.cites)
PopData$timber <- ifelse(PopData$timber == FALSE, 0, 10)
PopData$timber[!is.na(PopData$commercial) & PopData$commercial == FALSE] <- 5
PopData$timber[PopData$species %in% "Euterpe edulis"] <- 10
PopData$timber[PopData$species %in% "Dicksonia sellowiana"] <- 10
PopData$timber <- as.numeric(PopData$timber)

#### GETTING THREAT ECOLOGICAL GROUP INFORMATION  ####
hab1 <- hab[match(PopData$species, hab$Name_submitted), ]
table(hab1$Name_submitted == PopData$species)
early.sucession <- hab1$ecol.group
early.sucession[hab1$Name_submitted %in% "Ximenia americana"] <- "early_secondary"
early.sucession[early.sucession %in% c("late_secondary", "climax", "unknown")] <- 1L
early.sucession[early.sucession %in% c("early_secondary")] <- 0.9
early.sucession[early.sucession %in% c("pioneer")] <- 0.8
early.sucession <- as.numeric(early.sucession)

#########################################################################################################################################################H
#########################################################################################################################################################H
#############################
#### APPLYING CRITERIA D ####
#############################
#Reading the files needed from species distributions
spp <- readRDS("data/assess_iucn_spp.rds")
EOO <- readRDS("data/EOO.convex.hull_uncropped.rds")
EOO$Species <- as.character(EOO$Species)
AOO <- readRDS("data/AOO.rds")
AOO$tax <- as.character(AOO$tax)
localities <- readRDS("data/number_localities.rds")
localities$tax <- gsub("_", " ", as.character(localities$tax)) 
localities$tax <- gsub("\\.", "-", localities$tax) 
localities$Loc.level1 <- apply(localities[,2:3], 1, sum) 
localities$Loc.level2 <- apply(localities[,4:5], 1, sum) 
spp1<- merge(spp, EOO, by.x = "species.correct2", by.y = "Species", all.x = TRUE, sort= FALSE)
spp1<- merge(spp1, AOO, by.x = "species.correct2", by.y = "tax", all.x = TRUE, sort= FALSE)
spp1<- merge(spp1, localities[,c("tax","Loc.level1","Loc.level2")], by.x = "species.correct2", by.y = "tax", all.x = TRUE, sort= FALSE)

spp1$EOO.level.1[!is.na(spp1$EOO.level.1) & spp1$EOO.level.1 < spp1$AOO.level.1] <-
  spp1$AOO.level.1[!is.na(spp1$EOO.level.1) & spp1$EOO.level.1 < spp1$AOO.level.1]
spp1$EOO.level.2[!is.na(spp1$EOO.level.2) & spp1$EOO.level.2 < spp1$AOO.level.2] <-
  spp1$AOO.level.2[!is.na(spp1$EOO.level.2) & spp1$EOO.level.2 < spp1$AOO.level.2]
spp1$EOO.level.3[!is.na(spp1$AOO.level.3) &!is.na(spp1$EOO.level.3) & spp1$EOO.level.3 < spp1$AOO.level.3] <-
  spp1$AOO.level.3[!is.na(spp1$AOO.level.3) & !is.na(spp1$EOO.level.3) & spp1$EOO.level.3 < spp1$AOO.level.3]

#Creating and merging the df for analysis
df <- data.frame(species = mean.pop.sizes$species,
                 pop.size = mean.pop.sizes$`2018`,
                 pop.size.low = low.pop.sizes$`2018`,
                 pop.size.high = high.pop.sizes$`2018`,
                 p = PopData$p.est, 
                 exploitation = PopData$timber,
                 stringsAsFactors = FALSE)
#Taking into account exploitation over adut population
df$p1 <- df$p - (df$p * df$exploitation/100) 
df1 <- merge(df, spp1, by.x = "species", by.y = "species.correct2",
             all.x = TRUE, sort = FALSE)


## Should we apply the corrections to early-successional species here as well?

### ASSESSMENTS ###
#Optimal params
critD <- ConR::criterion_D(pop.size = df1$pop.size, 
                     Name_Sp = df1$species, 
                     AOO = df1$AOO.level.2,
                     n.Locs = df1$Loc.level2,
                     prop.mature = df1$p1, # p1 accounts for exploitation of commercial species
                     subcriteria = c("D", "D2"),
                     D.threshold = c(1000, 250, 50), 
                     AOO.threshold = 16, Loc.threshold = 2)
critD$D.low <- ConR::criterion_D(pop.size = df1$pop.size.low,
                           Name_Sp = df1$species,
                           prop.mature = df1$p1, subcriteria = c("D"),
                           D.threshold = c(1000, 250, 50))$D
critD$D.high <- ConR::criterion_D(pop.size = df1$pop.size.high,
                            Name_Sp = df1$species,
                            prop.mature = df1$p1, subcriteria = c("D"),
                            D.threshold = c(1000, 250, 50))$D

table(critD$D)
table(critD$D.low)
table(critD$D.high)

#Varying values of p
ps <- sort(c(1,.85,.72,.60,.49,.51,.58,.31,.25,.33,.64,.45,.35,.28,0.18,0.4))
critD.all <- cbind.data.frame(critD,
                              D = as.character(criterion_D(df1$pop.size, Name_Sp = df1$species, subcriteria = c("D"),
                                                           prop.mature = ps[1])[,c("D")]), stringsAsFactors = FALSE)
for(i in 2:length(ps)){
  critD.all <- cbind.data.frame(critD.all,
                                D= as.character(criterion_D(df1$pop.size, Name_Sp = df1$species, subcriteria = c("D"),
                                                            prop.mature = ps[i])[,c("D")]), stringsAsFactors = FALSE)
}

## Calculating the Red List Index for subcriterion A1 and A2
#Optmimal params
all.GL2 <- critD.all[,c(1:4,6:9,5,10:27)]
for(i in 9:27) all.GL2[,i] <- as.character(all.GL2[,i])
for(i in 9:27) all.GL2[,i] <- gsub("LC or NT", "LC", all.GL2[,i])


## Calculating the Red List Index for subcriterion A2
rli.all2 <- apply(all.GL2[,9:27], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL2[,9:27], 2, table)
rli.all.opt2 <- red::rli(gsub("LC or NT", "LC", as.character(all.GL2[,"category_D"])), boot = TRUE, runs = 4999)
rli.all.opt.low2 <- red::rli(gsub("LC or NT", "LC", as.character(all.GL2[,"D.low"])), boot = TRUE, runs = 4999)
rli.all.opt.high2 <- red::rli(gsub("LC or NT", "LC", as.character(all.GL2[,"D.high"])), boot = TRUE, runs = 4999)


## Renaming columns
names(all.GL2)[grepl("D\\.[0-9]", names(all.GL2))] <- 
  paste0("D.p", ps, sep = "")

## Renaming the LC category
all.GL2[] <- lapply(all.GL2, gsub, pattern = "^LC$", replacement = "LC or NT")

## Adding the low population estimates
table(all.GL2$species == df1$species)
all.GL2$pop.size.low <- df1$pop.size.low * df1$p1

#### Saving ####
saveRDS(all.GL2, "data/criterionD_all_prop_mature.rds")


###################
#### FIGURE SZ ####
###################

jpeg(filename = "figures/Figure_SZ.jpg", width = 2500, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
#par(mfrow=c(1,2))
par(mar=c(4,4,0.75,0.5), mgp=c(2.5,0.25,0),tcl=-0.2,las=1)
#optimum GL
plot(rev(rli.all2[2, grepl("D\\.[0-9]", colnames(rli.all2))]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch=19, ylim = c(0.99,1))
#axis(1, at=rev(ps), cex.axis = 1)
axis(1, at=seq(0.2,1,0.1), cex.axis = 1)
#axis(2, at=c(60,80,100,120,140,160,180,200,220), cex.axis = 1)
arrows(x0=rev(ps), y0 = rev(rli.all2[1,grepl("D\\.[0-9]", colnames(rli.all2))]), #using mean CIs
       y1 = rev(rli.all2[3,grepl("D\\.[0-9]", colnames(rli.all2))]),
       code = 3, angle = 90, length = 0.05)
#arrows(x0=rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][17:32]), #using low and high CIs
#       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][33:48]),
#       code = 3, angle = 90, length = 0.05, col=2)
#legend("topleft", expression(bold(A)), bty="n", cex=1.3)
abline(h=rli.all2[2,1], lty = 2)
legend("bottomright", c("Group-specific", "Fixed values"),
       lty = c(3,0), pch=c(NA,19),
       bty = "n", lwd=2)
dev.off()


#########################################################################################################################################################H
#########################################################################################################################################################H
#### COMPARING POP SIZE AND SPECIES DISTRIBUTION ####
#####################################################

df1$mature.pop.size <- df1$pop.size * df1$p
df1$mature.pop.size.exploited <- df1$pop.size * df1$p1
#Species information
df1 <- merge(df1, PopData[,c("species","life.form","dispersal.syndrome","SeedMass_g","habito")],
             by = "species", all.x = TRUE, sort= FALSE)
#Taxonomic information
taxonomy <- read.csv("data/sis_connect/taxonomy_threat.csv", as.is = TRUE)
taxonomy$species <- paste(taxonomy$genus, taxonomy$species)
df1 <- merge(df1, taxonomy[,c("species","classname","ordername")],
             by = "species", all.x = TRUE, sort= FALSE)
df1$taxonomy <- df1$classname
df1$taxonomy[df1$taxonomy %in% "Magnoliopsida"] <-
  df1$ordername[df1$taxonomy %in% "Magnoliopsida"]
df1$taxonomy[grepl("palm",df1$life.form)] <- "Palms"
df1$taxonomy[grepl("ucculent",df1$life.form)] <- "Cactus"
#Species endemism
end <- readRDS("data/sis_connect/endemism_threat.rds")
end$endemic[end$endemic %in% "not in the AF"] <- "occasional"
end$endemic[end$endemic %in% "not_endemic"] <- "widespread_sparse"
end <- end[!duplicated(end$species),]
df1 <- merge(df1, end[,c("species","endemism.level.1","endemism.level.2","endemic")],
             by = "species", all.x = TRUE, sort= FALSE)

## Exploring patterns
par(mfrow = c(2,2))
par(mar=c(4,4,0.75,0.5), mgp=c(2.5,0.25,0),tcl=-0.2,las=1)
plot(log(mature.pop.size) ~ log(treeco.occs), data = df1, col = (df1$endemic %in% "endemic")+1) # as expected, since one was generated from the other
abline(lm(log(mature.pop.size) ~ log(treeco.occs), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(non.dup.occs), data = df1, col = (df1$endemic %in% "endemic")+1)
abline(lm(log(mature.pop.size) ~ log(non.dup.occs), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ EOO.level.1, data = df1, col = (df1$endemic %in% "endemic")+1)
abline(lm(log(mature.pop.size) ~ EOO.level.1, data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(AOO.level.1), data = df1, col = (df1$endemic %in% "endemic")+1)
abline(lm(log(mature.pop.size) ~ log(AOO.level.1), data = df1), lwd=2,col=2)

#Number of occurences, non-duplicated occurrences, n.localities, EOO or AOO? AOO!
pairs(log(mature.pop.size) ~ EOO.level.1 + log(AOO.level.1), data = df1, col = (df1$endemic %in% "endemic")+1)
bbmle::AICtab(lm(log(mature.pop.size) ~ log(AOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(non.dup.occs), data = df1),
              lm(log(mature.pop.size) ~ log(total.occs), data = df1),
              lm(log(mature.pop.size) ~ log(Loc.level1), data = df1),
              lm(log(mature.pop.size) ~ log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(AOO.level.1) + log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(AOO.level.1) * log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(Loc.level1) * log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(non.dup.occs) * log(EOO.level.1), data = df1))

df2 <- na.omit(df1[,c("mature.pop.size","EOO.level.1","AOO.level.1","taxonomy","habito","endemic")])
mod <- lm(log(mature.pop.size) ~ log(EOO.level.1) * 
            log(AOO.level.1), data = df2)
options(na.action = "na.fail")
dd <- MuMIn::dredge(mod, rank = "AIC")
subset(dd, delta < 4)
options(na.action = "na.omit")
car::Anova(mod)
car::avPlots(mod)
car::vif(mod) ## vif very high for the interaction
car::vif(lm(log(mature.pop.size) ~ log(EOO.level.1) + log(AOO.level.1), data = df1)) #ok

## Taking into account unequal variance and random effects

#Getting the data.frame in the good format for predictions
df2 <- na.omit(df1[,c("mature.pop.size","EOO.level.1","AOO.level.1","taxonomy","habito","endemic")])
df2$PopSize <- log(df2$mature.pop.size)
df2$EOO <- log(df2$EOO.level.1)
df2$AOO <- log(df2$AOO.level.1)
df2$taxonomy <- as.factor(df2$taxonomy)
df2$habito <- as.factor(df2$habito)
df2$endemic <- as.factor(df2$endemic)
df2 <- df2[,c("PopSize","AOO","EOO","taxonomy","habito","endemic")]
# df2$cor <- factor(df2$endemic, levels = sort(unique(df2$endemic)),
#                                              labels = c("darkred","yellow","red","orange")) 

##Fitting the models
mod.ols <- stats::lm(PopSize ~ EOO + #Ordinary Least Squares
                       AOO, data = df2)
mod.rob <- MASS::rlm(PopSize ~ EOO + #Robust lm
                       AOO, data = df2)
mod.gls <- nlme::gls(PopSize ~ EOO + #Gen. Least Squares 
                       AOO, data = df2)
mod.gls.pow <- nlme::gls(PopSize ~ EOO + #Gen. Least Squares 
                           AOO, data = df2, weights = nlme::varPower())
mod.gls.exp <- nlme::gls(PopSize ~ EOO + #Gen. Least Squares 
                           AOO, data = df2, weights = nlme::varExp())
mod.gls.fix <- nlme::gls(PopSize ~ EOO + #Gen. Least Squares 
                           AOO, data = df2, weights = nlme::varFixed(~AOO))
mod.gls.fix1 <- nlme::gls(PopSize ~ EOO + #Gen. Least Squares 
                            AOO, data = df2, weights = nlme::varFixed(~EOO))
mod.mix <- nlme::lme(PopSize ~ EOO + #Gen. Least Squares 
                            AOO, random = (~1|taxonomy), data = df2, weights = nlme::varFixed(~EOO))
mod.mix1 <- nlme::lme(PopSize ~ EOO + #Gen. Least Squares 
                       AOO + habito, random = (~1|taxonomy), data = df2, weights = nlme::varFixed(~AOO))
mod.mix2 <- nlme::lme(PopSize ~ AOO + habito, 
                      random = (~1|taxonomy), data = df2, weights = nlme::varPower())
mod.mix3 <- nlme::lme(PopSize ~ EOO + #Gen. Least Squares 
                        AOO + habito, random = (~1|endemic), data = df2, weights = nlme::varFixed(~AOO))
mod.mix4 <- nlme::lme(PopSize ~ EOO + #Gen. Least Squares 
                        AOO + habito, random = list(~1|taxonomy, ~1|endemic), data = df2, weights = nlme::varFixed(~AOO))
mod.mix5 <- nlme::lme(PopSize ~ EOO + #Gen. Least Squares 
                        AOO + habito, random = (~AOO|endemic), data = df2, weights = nlme::varFixed(~AOO))
mod.gls.fix2 <- nlme::gls(PopSize ~ AOO + endemic +
                            habito, data = df2, weights = nlme::varPower())
mod.gls.fix3 <- nlme::gls(PopSize ~ AOO * endemic +
                            habito, data = df2, weights = nlme::varPower())
mod.gls.fix4 <- nlme::gls(PopSize ~ EOO + AOO * endemic +
                            habito, data = df2, weights = nlme::varPower())
mod.mix6 <- lme4::lmer(PopSize ~ AOO + (AOO|taxonomy) + 
                        endemic + habito, data = df2, REML = TRUE)
mod.mix7 <- lme4::lmer(PopSize ~ AOO + (AOO|taxonomy) + 
                         EOO + endemic + habito, data = df2, REML = TRUE)
mod.mix8 <- lme4::lmer(PopSize ~ AOO + (AOO|endemic) +
                         (AOO|taxonomy) + EOO + habito, data = df2, REML = TRUE)
mod.mix9 <- lme4::lmer(PopSize ~ AOO + (AOO|taxonomy) + 
                         EOO + endemic + habito +
                         AOO:endemic +
                         AOO:habito, data = df2, REML = TRUE)
mod.mix10 <- lme4::lmer(PopSize ~ AOO + (AOO|taxonomy) + 
                         endemic + habito +
                         #AOO:endemic +
                         AOO:habito, data = df2, REML = TRUE)

bbmle::AICtab(mod.ols, mod.rob, mod.gls, 
              mod.gls.pow, mod.gls.exp, mod.gls.fix, mod.gls.fix1, mod.gls.fix2, mod.gls.fix3, mod.gls.fix4,
              mod.mix, mod.mix1, mod.mix2, mod.mix3, mod.mix4, mod.mix5, mod.mix6, mod.mix7, mod.mix8, mod.mix9, mod.mix10)
car::Anova(mod.mix7)
car::vif(mod.mix7)
plot(mod.mix7)
#car::avPlot(mod.mix7)
summary(mod.mix7)
r2glmm::r2beta(mod.mix7, partial = FALSE)
piecewiseSEM::rsquared(mod.mix7)$Marginal
piecewiseSEM::rsquared(mod.mix7)$Conditional
lmerTest::ranova(mod.mix7)

## Inspectingthe results of the regression using jtools
require(ggplot2)
require(jtools)
require(ggstance)
require(broom)
require(broom.mixed)
#Tidy summary
jtools::summ(mod.mix7, confint = TRUE, vifs = TRUE,
             scale= FALSE, digits = 2, data = df2)
#Size effects
jtools::plot_summs(mod.mix7, scale = TRUE)
#Comparing size effects between two or more models
jtools::plot_summs(mod.mix7, mod.mix8, scale = TRUE, model.names = c("mix7", "mix8"))

#Partial residual plots
jtools::effect_plot(mod.mix7, pred = AOO, 
                         interval = TRUE, plot.points = TRUE, partial.residuals = TRUE, 
                         int.type = "confidence",
                         #int.type = "prediction",
                         data = df2)
jtools::effect_plot(mod.mix7, pred = EOO, 
            interval = TRUE, plot.points = TRUE, partial.residuals = TRUE, data = df2)
jtools::effect_plot(mod.mix7, pred = habito, 
            interval = TRUE, plot.points = TRUE, partial.residuals = TRUE, 
            data = df2, jitter = c(0.2,0))
jtools::effect_plot(mod.mix7, pred = endemic, 
            interval = TRUE, plot.points = TRUE, partial.residuals = TRUE, 
            data = df2, jitter = c(0.2,0))


## A model only for endemic species...
toto <- df2[df2$endemic %in% "endemic",]
mod.mix7.end <- lme4::lmer(PopSize ~ AOO + (AOO|taxonomy) + 
                         EOO + habito, data = toto, REML = TRUE)
summary(mod.mix7.end)
stats::coef(mod.mix7.end)
lme4::fixef(mod.mix7.end)
lme4::ranef(mod.mix7.end)

# par(mfrow=c(1,1))
# plot(PopSize ~ AOO, data = df2, xlim=c(0,10),ylim = c(0,20))
# abline(lm(PopSize ~ AOO - 1, data = df2), lwd=2, col=6)
# abline(v=log(4), lty=3); abline(h=log(1000), lty=3) 
# curve(5.56 + 1.19*x, add=TRUE, lwd=2, col="red") #gls fixed 2
# curve(5.44 + 1.21*x, add=TRUE, lwd=2, col="green") #gls fixed 3
# curve(6.57 + 1.45*x, add=TRUE, lwd=2, col="blue") #gls fix 4
# curve(6.38 + 1.128*x, add=TRUE, lwd=2, col="cyan") # mixed 6 (average per order)
# curve(7.448 + 1.28*x, add=TRUE, lwd=2, col="orange") # mixed 7
# curve(7.25 + 1.42*x, add=TRUE, lwd=2, col="purple") # mixed 8
# ## EVEN AT THE MINIMUM AOO (i.e log(4)) THE MODELS PREDICT MORE THEN 1000 INDIVIDUALS... 


# ## Inspecting the predicted values of each model
# par(mfrow=c(2,3))
# plot(exp(predict(mod.mix4)) ~ df2$mature.pop.size,
#      log = "xy", xlab = "Obs", ylab = "Pred", ylim = c(1000,1e+08)); abline(0,1,lwd=2,lty=2); abline(h=1000,lwd=2,lty=3,col=2); abline(v=1000,lwd=2,lty=3,col=2)
# plot(exp(predict(mod.mix6)) ~ df2$mature.pop.size,
#      log = "xy", xlab = "Obs", ylab = "Pred", ylim = c(1000,1e+08)); abline(0,1,lwd=2,lty=2); abline(h=1000,lwd=2,lty=3,col=2); abline(v=1000,lwd=2,lty=3,col=2)
# plot(exp(predict(mod.mix7)) ~ df2$mature.pop.size,
#      log = "xy", xlab = "Obs", ylab = "Pred", ylim = c(1000,1e+08)); abline(0,1,lwd=2,lty=2); abline(h=1000,lwd=2,lty=3,col=2); abline(v=1000,lwd=2,lty=3,col=2)
# plot(exp(predict(mod.gls.fix2)) ~ df2$mature.pop.size,
#      log = "xy", xlab = "Obs", ylab = "Pred", ylim = c(1000,1e+08)); abline(0,1,lwd=2,lty=2); abline(h=1000,lwd=2,lty=3,col=2); abline(v=1000,lwd=2,lty=3,col=2)
# plot(exp(predict(mod.gls.fix3)) ~ df2$mature.pop.size,
#      log = "xy", xlab = "Obs", ylab = "Pred", ylim = c(1000,1e+08)); abline(0,1,lwd=2,lty=2); abline(h=1000,lwd=2,lty=3,col=2); abline(v=1000,lwd=2,lty=3,col=2)
# plot(exp(predict(mod.gls.fix4)) ~ df2$mature.pop.size,
#      log = "xy", xlab = "Obs", ylab = "Pred", ylim = c(1000,1e+08)); abline(0,1,lwd=2,lty=2); abline(h=1000,lwd=2,lty=3,col=2); abline(v=1000,lwd=2,lty=3,col=2)


#### MAKING THE PREDICTIONS OF POP. SIZE FOR SPECIES WITHOUT ABUNDANCE DATA
df3 <- taxonomy[,c("classname","ordername","family","species")]
df3 <- merge(df3, spp1, by.x = "species", by.y = "species.correct2",
             all.x = TRUE, sort = FALSE)
df3 <- merge(df3, hab[,c("Name_submitted","life.form","dispersal.syndrome","SeedMass_g","habito")],
             by.x = "species", by.y = "Name_submitted", all.x = TRUE, sort= FALSE)
df3 <- merge(df3, end[,c("species","endemic")],
             by = "species", all.x = TRUE, sort= FALSE)
df3$taxonomy <- df3$classname
df3$taxonomy[df3$taxonomy %in% "Magnoliopsida"] <-
  df3$ordername[df3$taxonomy %in% "Magnoliopsida"]
df3$taxonomy[grepl("palm",df3$life.form)] <- "Palms"
df3$taxonomy[grepl("ucculent",df3$life.form)] <- "Cactus"

##Getting the predictions
new.dat <- na.omit(df3[,c("species","EOO.level.1","AOO.level.1","habito","endemic","taxonomy")])
new.dat$EOO <- log(new.dat$EOO.level.1)
new.dat$AOO <- log(new.dat$AOO.level.1)
new.dat$taxonomy <- as.factor(new.dat$taxonomy)
new.dat$habito <- as.factor(new.dat$habito)
new.dat$endemic <- as.factor(new.dat$endemic)

ids <- new.dat$endemic %in% "endemic"
summary(exp(predict(mod.mix6, new.dat[ids,], re.form=NULL))) # dAIC: 61.7
summary(exp(predict(mod.mix7, new.dat[ids,], re.form=NULL))) # dAIC: 34.5 (media: 122,775)
summary(exp(predict(mod.mix8, new.dat[ids,], re.form=NULL))) # best model (media: 126,255)
## mod.mix8 ganhou, mas como END tem poucos niveis para ser randomica, usaremos o mod.mix7, que tb teve os mesnores VIFs

best.model <- mod.mix7
foreach::registerDoSEQ()
preds <- exp(merTools::predictInterval(best.model, new.dat, include.resid.var = TRUE,
                                       fix.intercept.variance = TRUE,
                                       which = "full", level = 0.9, n.sims = 5000))
new.dat1 <- cbind.data.frame(new.dat, preds, 
                             stringsAsFactors = FALSE)
names(new.dat1)[(dim(new.dat1)[2]-2):dim(new.dat1)[2]] <- c("pred","pred.high","pred.low")
# preds1 <- exp(merTools::predictInterval(best.model, new.dat, include.resid.var = FALSE,
#                                         fix.intercept.variance = TRUE,
#                                        which = "full", level = 0.9, n.sims = 5000))
# new.dat1 <- cbind.data.frame(new.dat1, preds1, 
#                              stringsAsFactors = FALSE)
# names(new.dat1)[(dim(new.dat1)[2]-2):dim(new.dat1)[2]] <- c("pred1","pred1.high","pred1.low")

# #create design matrix
# Designmat <- model.matrix(eval(eval(best.model$call$fixed)[-2]), new.dat[,2:5])
# 
# #compute standard error for predictions
# predvar <- diag(Designmat %*% best.model$varFix %*% t(Designmat))
# new.dat$SE <- sqrt(predvar) 
# new.dat$SE2 <- sqrt(predvar + best.model$sigma^2)
# new.dat$low <- new.dat$pred - 2*new.dat$SE
# new.dat$high <- new.dat$pred + 2*new.dat$SE
# new.dat$low2 <- new.dat$pred - 2*new.dat$SE2
# new.dat$high2 <- new.dat$pred + 2*new.dat$SE2

# #re-tranforming
# new.dat$pred <- exp(new.dat$pred)
# new.dat$low <- exp(new.dat$low)
# new.dat$high <- exp(new.dat$high)
# new.dat$low2 <- exp(new.dat$low2)
# new.dat$high2 <- exp(new.dat$high2)

## SUMMARY ##
sumario <- jtools::summ(best.model, confint = TRUE, vifs = TRUE,
             scale= FALSE, digits = 2, data = df2)
anova.best <- anova(best.model, lm(PopSize ~1, data = df2))
preds.obs <- merTools::predictInterval(best.model, include.resid.var = TRUE,
                                       fix.intercept.variance = TRUE,
                                       which = "full", level = 0.9, n.sims = 5000) 

jpeg(filename = "figures/Figure_ZZ.jpg", width = 2500, height = 2500, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(2,2))
par(mar=c(4,4,0.75,0.5), mgp=c(2,0.4,0),tcl=-0.2,las=1)
## Partial-regression plots: AOO
p <- jtools::effect_plot(best.model, pred = AOO, 
                         interval = TRUE, plot.points = TRUE, partial.residuals = TRUE, 
                         int.type = "confidence",
                         #int.type = "prediction",
                         data = df2)
pg <- ggplot2::ggplot_build(p)
plot(pg[[1]][[1]][df2$habito %in% "tree",1:2], 
     pch =19, col = adjustcolor("red", alpha.f = 0.4),
     xlim=c(1,10),ylim = c(2.5,20), xlab = "log(AOO)", ylab = "log(Pop. size)")
points(pg[[1]][[1]][df2$habito %in% "shrub",1:2],  
       col = adjustcolor("blue", alpha.f = 0.4), pch = 15)
lines(pg[[1]][[2]][,1:2], col=1, lwd = 2)
lines(pg[[1]][[3]][,c(1,3)], col=1, lwd = 1, lty = 3)
lines(pg[[1]][[3]][,c(1,4)], col=1, lwd = 1, lty = 3)
#abline(v=log(4), lty=3) 
abline(h=log(1000), lty=3) 
legend("topleft", c("Trees", "Shrubs"), 
       pch = c(19,15), col=c("red", "blue"), cex=1.1, bty="n")

## Partial-regression plots: EOO
p <- jtools::effect_plot(best.model, pred = EOO, 
                         interval = TRUE, plot.points = TRUE, partial.residuals = TRUE, 
                         int.type = "confidence",
                         #int.type = "prediction",
                         data = df2)
pg <- ggplot2::ggplot_build(p)
par(mar=c(4,4,0.75,0.5), mgp=c(2,0.4,0),tcl=-0.2,las=1)
plot(pg[[1]][[1]][df2$habito %in% "tree",1:2], 
     pch =19, col = adjustcolor("red", alpha.f = 0.4),
     xlim=c(1,18),ylim = c(2.5,20), xlab = "log(EOO)", ylab = "log(Pop. size)")
points(pg[[1]][[1]][df2$habito %in% "shrub",1:2],  
       col = adjustcolor("blue", alpha.f = 0.4), pch = 15)
lines(pg[[1]][[2]][,1:2], col=1, lwd = 2)
lines(pg[[1]][[3]][,c(1,3)], col=1, lwd = 1, lty = 3)
lines(pg[[1]][[3]][,c(1,4)], col=1, lwd = 1, lty = 3)
abline(h=log(1000), lty=3) 
#abline(v=log(4), lty=3)
legend("bottomleft", c("Trees", "Shrubs"), 
       pch = c(19,15), col=c("red", "blue"), cex=1.1, bty="n")

##Predicted vs. observed
plot(df2$PopSize[df2$habito %in% "tree"],
     preds.obs$fit[df2$habito %in% "tree"],
     pch =19, col = adjustcolor("red", alpha.f = 0.4),
     xlim=c(3,20),ylim = c(8,19), xlab = "Obs. Pop. Size", ylab = "Pred. Pop. Size")
points(df2$PopSize[df2$habito %in% "shrub"],
       preds.obs$fit[df2$habito %in% "shrub"],
       col = adjustcolor("blue", alpha.f = 0.4), pch = 15)
abline(0,1, lty=2)
legend("topleft",legend=do.call(expression, list(
  bquote({Marginal~italic(R)^2}~"="~.(round(attributes(sumario)$rsqs[[1]]*100,1))~"%"),
  bquote({Conditional~italic(R)^2}~"="~.(round(attributes(sumario)$rsqs[[2]]*100,1))~"%"),
  bquote("Number of obs. = "~.(attributes(sumario)$n)),
  bquote(italic(chi)^2~"="~.(round(anova.best$Chisq[2],1)))
  #bquote(italic(F)~"- test="~.(round(summary(toto)$fstatistic[1],2)))
)),
bty="n",horiz=F,cex=1.1,adj=c(0.1,0.1)
)
legend("bottomright", c("Trees", "Shrubs"), 
       pch = c(19,15), col=c("red", "blue"), cex=1.1, bty="n")

##Mean and low Predictions
hist(df2$PopSize, nclass=40, probability = FALSE, ylim=c(0,310), 
     col= adjustcolor("black", alpha.f = 0.3), main="",xlab="log(Population size)")
# ids <- new.dat1$species %in% df1$species
# hist(log(new.dat1$pred[ids]), nclass=20, add=TRUE, border="green", probability = TRUE,
#      col= adjustcolor("green", alpha.f = 0.3))
abline(v = log(1000), lty=3)
ids1 <- !new.dat1$species %in% df1$species
hist(log(new.dat1$pred[ids1]), nclass=20, add=TRUE, border="green", probability = FALSE,
     col= adjustcolor("green", alpha.f = 0.35))
hist(log(new.dat1$pred.low[ids1]), nclass=20, add=TRUE, border="darkorange", probability = FALSE,
      col= adjustcolor("darkorange", alpha.f = 0.3))
abline(v = log(1000), lty=3)
#abline(v = log(c(50,250,1000)), lty=3)
100*table(new.dat1$pred.low[ids1] < 1000)/length(new.dat1$pred.low[ids1])
table(new.dat1$pred.low[ids1] < 250)
table(new.dat1$pred.low[ids1] < 50)
legend("topright", c("Observed","Mean Predict.","Low Predict."),
       pch=15, col=c("darkgrey","green","orange"), bty= "n", cex=1.1)
dev.off()

#### Saving ####
saveRDS(new.dat1, "data/estimated_pop_size_from_AOO.rds")
