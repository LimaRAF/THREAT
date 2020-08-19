#######################################################
#### ASSESSING VERY SMALL POPULATIONS - CRITERIA D ####
#######################################################
rm(list=ls())

#### LOADING PACKAGES ###
#devtools::install_github("gdauby/ConR", ref = "master", force = TRUE) # old version
#devtools::install_github("gdauby/ConR@devel") # new version on GitHub
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
 repos = NULL, 
 type = "source")
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
PopData <- PopData[!is.na(PopData$taxon_id),]
names(PopData)[1] <- "species"
PopData <- PopData[order(PopData$species),]

## Filtering the populational datasets
ids <- mean.pop.sizes$species %in% PopData$species
mean.pop.sizes <- mean.pop.sizes[ids,]
low.pop.sizes <- low.pop.sizes[ids,]
high.pop.sizes <- high.pop.sizes[ids,]
table(mean.pop.sizes$species == PopData$species)
rm(res.means)


#########################################################################################################################################################H
#########################################################################################################################################################H
#############################
#### APPLYING CRITERIA D ####
#############################
#Reading the files needed from species distributions
spp <- readRDS("data/assess_iucn_spp.rds")
EOO <- readRDS("data/EOO.convex.hull.rds")
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
spp1$EOO.level.3[!is.na(spp1$EOO.level.3) & spp1$EOO.level.3 < spp1$AOO.level.3] <-
  spp1$AOO.level.3[!is.na(spp1$EOO.level.3) & spp1$EOO.level.3 < spp1$AOO.level.3]


#Creating and merging the df for analysis
df <- data.frame(species = mean.pop.sizes$species,
                 pop.size = mean.pop.sizes$`2018`,
                 pop.size.low = low.pop.sizes$`2018`,
                 pop.size.high = high.pop.sizes$`2018`,
                 p = PopData$p.est, stringsAsFactors = FALSE)
df1 <- merge(df, spp1, by.x = "species", by.y = "species.correct2",
             all.x = TRUE, sort = FALSE)

### ASSESSMENTS ###
#Optimal params
critD <- criterion_D(pop.size = df1$pop.size, 
                     Name_Sp = df1$species, 
                     AOO = df1$AOO.level.2,
                     n.Locs = df1$Loc.level2,
                     prop.mature = df1$p,
                     subcriteria = c("D", "D2"),
                     D.threshold = c(1000, 250, 50), 
                     AOO.threshold = 16, Loc.threshold = 2)
critD$D.low <- criterion_D(pop.size = df1$pop.size.low,
                           Name_Sp = df1$species,
                           prop.mature = df1$p, subcriteria = c("D"),
                           D.threshold = c(1000, 250, 50))$D
critD$D.high <- criterion_D(pop.size = df1$pop.size.high,
                            Name_Sp = df1$species,
                            prop.mature = df1$p, subcriteria = c("D"),
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
all.GL2 <- critD.all[,c(1:4,6:9,5,10:25)]
for(i in 9:25) all.GL2[,i] <- as.character(all.GL2[,i])
for(i in 9:25) all.GL2[,i] <- gsub("LC or NT", "LC", all.GL2[,i])

rli.all2 <- apply(all.GL2[,9:25], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL2[,9:25], 2, table)

###################
#### FIGURE SZ ####
###################

jpeg(filename = "figures/Figure_SZ.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
#par(mfrow=c(1,2))
par(mar=c(4,4,0.75,0.5), mgp=c(2.5,0.25,0),tcl=-0.2,las=1)
#optimum GL
plot(rev(rli.all2[2, grepl("D\\.", colnames(rli.all2))]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch=19, ylim = c(0.99,1))
#axis(1, at=rev(ps), cex.axis = 1)
axis(1, at=seq(0.2,1,0.1), cex.axis = 1)
#axis(2, at=c(60,80,100,120,140,160,180,200,220), cex.axis = 1)
arrows(x0=rev(ps), y0 = rev(rli.all2[1,grepl("D\\.", colnames(rli.all2))]), #using mean CIs
       y1 = rev(rli.all2[3,grepl("D\\.", colnames(rli.all2))]),
       code = 3, angle = 90, length = 0.05)
#arrows(x0=rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][17:32]), #using low and high CIs
#       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][33:48]),
#       code = 3, angle = 90, length = 0.05, col=2)
#legend("topleft", expression(bold(A)), bty="n", cex=1.3)
abline(h=rli.all2[2,1], lty = 2)
legend("bottomright", c("Group-specific", "Fixed"),
       lty = c(3,0), pch=c(NA,19),
       bty = "n", lwd=2)
dev.off()


#########################################################################################################################################################H
#########################################################################################################################################################H
#### COMPARING POP SIZE AND SPECIES DISTRIBUTION ####
#####################################################

df1$mature.pop.size <- df1$pop.size * df1$p
#Species information
df1 <- merge(df1, PopData[,c("species","life.form","dispersal.syndrome","SeedMass_g","habito")],
             by = "species", all.x = TRUE, sort= FALSE)
#Taxonomic information
taxonomy <- read.csv("data/threat_taxonomy.csv", as.is = TRUE)
df1 <- merge(df1, taxonomy[,c("species","classname","ordername")],
             by = "species", all.x = TRUE, sort= FALSE)
df1$taxonomy <- df1$classname
df1$taxonomy[df1$taxonomy %in% "Magnoliopsida"] <-
  df1$ordername[df1$taxonomy %in% "Magnoliopsida"]
df1$taxonomy[grepl("palm",df1$life.form)] <- "Palms"
df1$taxonomy[grepl("ucculent",df1$life.form)] <- "Cactus"
#Species endemism
end <- readRDS("data/threat_endeism_levels.rds")
end$endemic[end$endemic %in% "not in the AF"] <- "occasional"
end$endemic[end$endemic %in% "not_endemic"] <- "widespread_sparse"
df1 <- merge(df1, end[,c("species","endemism.level.1","endemism.level.2","endemic")],
             by = "species", all.x = TRUE, sort= FALSE)

## Exploring patterns
par(mfrow = c(2,2))
par(mar=c(4,4,0.75,0.5), mgp=c(2.5,0.25,0),tcl=-0.2,las=1)
plot(log(mature.pop.size) ~ log(treeco.occs), data = df1, col = (df1$endemic %in% "endemic")+1) # as expected, since one was generated from the other
abline(lm(log(mature.pop.size) ~ log(treeco.occs), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(non.dup.occs), data = df1, col = (df1$endemic %in% "endemic")+1)
abline(lm(log(mature.pop.size) ~ log(non.dup.occs), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(EOO.level.1), data = df1, col = (df1$endemic %in% "endemic")+1)
abline(lm(log(mature.pop.size) ~ log(EOO.level.1), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(AOO.level.1), data = df1, col = (df1$endemic %in% "endemic")+1)
abline(lm(log(mature.pop.size) ~ log(AOO.level.1), data = df1), lwd=2,col=2)

#Number of occurences, non-duplicated occurrences, n.localities, EOO or AOO? AOO!
pairs(log(mature.pop.size) ~ log(EOO.level.1) + log(AOO.level.1), data = df1, col = (df1$endemic %in% "endemic")+1)
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
mod.ols <- stats::lm(log(mature.pop.size) ~ log(EOO.level.1) + #Ordinary Least Squares
                       log(AOO.level.1), data = df2)
mod.rob <- MASS::rlm(log(mature.pop.size) ~ log(EOO.level.1) + #Robust lm
                       log(AOO.level.1), data = df2)
mod.gls <- nlme::gls(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                       log(AOO.level.1), data = df2)
mod.gls.pow <- nlme::gls(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                           log(AOO.level.1), data = df2, weights = nlme::varPower())
mod.gls.exp <- nlme::gls(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                           log(AOO.level.1), data = df2, weights = nlme::varExp())
mod.gls.fix <- nlme::gls(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                           log(AOO.level.1), data = df2, weights = nlme::varFixed(~log(AOO.level.1)))
mod.gls.fix1 <- nlme::gls(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                            log(AOO.level.1), data = df2, weights = nlme::varFixed(~log(EOO.level.1)))
mod.mix <- nlme::lme(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                            log(AOO.level.1), random = (~1|taxonomy), data = df2, weights = nlme::varFixed(~log(EOO.level.1)))
mod.mix1 <- nlme::lme(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                       log(AOO.level.1) + as.factor(habito), random = (~1|taxonomy), data = df2, weights = nlme::varFixed(~log(AOO.level.1)))
mod.mix2 <- nlme::lme(log(mature.pop.size) ~ log(AOO.level.1) + as.factor(habito), 
                      random = (~1|taxonomy), data = df2, weights = nlme::varPower())
mod.mix3 <- nlme::lme(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                        log(AOO.level.1) + as.factor(habito), random = (~1|endemic), data = df2, weights = nlme::varFixed(~log(AOO.level.1)))
mod.mix4 <- nlme::lme(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                        log(AOO.level.1) + as.factor(habito), random = list(~1|taxonomy, ~1|endemic), data = df2, weights = nlme::varFixed(~log(AOO.level.1)))
mod.mix5 <- nlme::lme(log(mature.pop.size) ~ log(EOO.level.1) + #Gen. Least Squares 
                        log(AOO.level.1) + as.factor(habito), random = (~log(AOO.level.1)|endemic), data = df2, weights = nlme::varFixed(~log(AOO.level.1)))
mod.gls.fix2 <- nlme::gls(log(mature.pop.size) ~ log(AOO.level.1) + as.factor(endemic) +
                            as.factor(habito), data = df2, weights = nlme::varPower())
mod.gls.fix3 <- nlme::gls(log(mature.pop.size) ~ log(AOO.level.1) * as.factor(endemic) +
                            as.factor(habito), data = df2, weights = nlme::varPower())
mod.gls.fix4 <- nlme::gls(log(mature.pop.size) ~ log(EOO.level.1) + log(AOO.level.1) * as.factor(endemic) +
                            as.factor(habito), data = df2, weights = nlme::varPower())
mod.mix6 <- lme4::lmer(log(mature.pop.size) ~ log(AOO.level.1) + (log(AOO.level.1)|taxonomy) + 
                        as.factor(endemic) + as.factor(habito), data = df2, REML = TRUE)
mod.mix7 <- lme4::lmer(log(mature.pop.size) ~ log(AOO.level.1) + (log(AOO.level.1)|taxonomy) + 
                         log(EOO.level.1) + as.factor(endemic) + as.factor(habito), data = df2, REML = TRUE)
mod.mix8 <- lme4::lmer(log(mature.pop.size) ~ log(AOO.level.1) + (log(AOO.level.1)|endemic) +
                         (log(AOO.level.1)|taxonomy) + log(EOO.level.1) + as.factor(habito), data = df2, REML = TRUE)

bbmle::AICtab(mod.ols, mod.rob, mod.gls, 
              mod.gls.pow, mod.gls.exp, mod.gls.fix, mod.gls.fix1, mod.gls.fix2, mod.gls.fix3, mod.gls.fix4,
              mod.mix, mod.mix1, mod.mix2, mod.mix3, mod.mix4, mod.mix5, mod.mix6, mod.mix7, mod.mix8)
car::Anova(mod.mix7)
car::vif(mod.mix7)
plot(mod.mix7)
car::avPlot(mod.mix7)
summary(mod.mix7)

par(mfrow=c(1,1))
plot(log(mature.pop.size) ~ log(AOO.level.1), data = df1, xlim=c(0,10),ylim = c(0,20))
abline(lm(log(mature.pop.size) ~ log(AOO.level.1) - 1, data = df1), lwd=2, col=6)
abline(v=log(4), lty=3); abline(h=log(1000), lty=3) 
curve(5.56 + 1.19*x, add=TRUE, lwd=2, col="red") #gls fixed 2
curve(5.44 + 1.21*x, add=TRUE, lwd=2, col="green") #gls fixed 3
curve(6.57 + 1.45*x, add=TRUE, lwd=2, col="blue") #gls fix 4
curve(6.38 + 1.128*x, add=TRUE, lwd=2, col="cyan") # mixed 6 (average per order)
curve(7.448 + 1.28*x, add=TRUE, lwd=2, col="orange") # mixed 7
curve(7.25 + 1.42*x, add=TRUE, lwd=2, col="purple") # mixed 8
## EVEN AT THE MINIMUM AOO (i.e log(4)) THE MODELS PREDICT MORE THEN 1000 INDIVIDUALS... 


## Inspecting the predicted values of each model
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
df.miss <- taxonomy[!taxonomy$species %in% df1$species,c("classname","ordername","family","species")]
df.miss <- merge(df.miss, spp1, by.x = "species", by.y = "species.correct2",
             all.x = TRUE, sort = FALSE)
df.miss <- merge(df.miss, hab[,c("Name_submitted","life.form","dispersal.syndrome","SeedMass_g","habito")],
             by.x = "species", by.y = "Name_submitted", all.x = TRUE, sort= FALSE)
df.miss <- merge(df.miss, end[,c("species","endemic")],
             by = "species", all.x = TRUE, sort= FALSE)
df.miss$taxonomy <- df.miss$classname
df.miss$taxonomy[df.miss$taxonomy %in% "Magnoliopsida"] <-
  df.miss$ordername[df.miss$taxonomy %in% "Magnoliopsida"]
df.miss$taxonomy[grepl("palm",df.miss$life.form)] <- "Palms"
df.miss$taxonomy[grepl("ucculent",df.miss$life.form)] <- "Cactus"

##Getting the predictions
new.dat <- na.omit(df.miss[,c("species","EOO.level.1","AOO.level.1","habito","endemic","taxonomy")])
summary(exp(predict(mod.gls.fix2, new.dat, re.form=NULL))) # dAIC: 117 (media: 94,000)
summary(exp(predict(mod.gls.fix3, new.dat, re.form=NULL))) # dAIC: 96 (media: 106,690) # high vifs...
summary(exp(predict(mod.gls.fix4, new.dat, re.form=NULL))) # dAIC: 62 (media: 114,859) # high vifs...
summary(exp(predict(mod.mix6, new.dat, re.form=NULL))) # dAIC: 55.6 (media: 110,185)
summary(exp(predict(mod.mix7, new.dat, re.form=NULL))) # dAIC: 32.1 (media: 122,775)
summary(exp(predict(mod.mix8, new.dat, re.form=NULL))) # best model (media: 126,255)
best.model <- mod.mix8
preds <- exp(merTools::predictInterval(tmp1, new.dat1, include.resid.var = TRUE,
                                    which = "full", level = 0.9, n.sims = 5000))
new.dat1 <- cbind.data.frame(new.dat, preds, 
                             stringsAsFactors = FALSE)
names(new.dat1)[(dim(new.dat1)[2]-2):dim(new.dat1)[2]] <- c("pred","pred.high","pred.low")

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
par(mfrow=c(1,2))
hist(log(df2$mature.pop.size), nclass=40, probability = TRUE, ylim=c(0,0.5), 
     col= adjustcolor("grey", alpha.f = 0.3), main="",xlab="log(Population size)")
hist(log(new.dat1$pred), nclass=20, add=TRUE, border=4, probability = TRUE,
     col= adjustcolor("blue", alpha.f = 0.3))
hist(log(new.dat1$pred.low), nclass=20, add=TRUE, border=2, probability = TRUE,
     col= adjustcolor("red", alpha.f = 0.3))
abline(v = log(c(50,250,1000)), lty=2)
table(new.dat1$pred.low < 1000)
table(new.dat1$pred.low < 250)
table(new.dat1$pred.low < 50)

boxplot(log(new.dat1$pred) ~ new.dat1$endemic, varwidth = TRUE, notch= TRUE)


#### Saving ####
saveRDS(all.GL2, "data/criterionD_all_prop_mature.rds")
saveRDS(new.dat1, "data/estimated_pop_size_no_abundance.rds")
