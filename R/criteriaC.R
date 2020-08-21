###################################################
###################################################
#### ASSESSING POPULATION DECLINE - CRITERIA C ####
###################################################
###################################################
rm(list=ls())

#### LOADING PACKAGES ###
#devtools::install_github("gdauby/ConR", ref = "master", force = TRUE) # old version
#devtools::install_github("gdauby/ConR@devel") # new version on GitHub
detach("package:ConR", unload=TRUE)
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


#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA C ####
#############################

#### NEED TO RUN FOR low and high POP SIZES?

## Getting info on the number of subpopulations for each species
subpops <- readRDS("data/number.subpops.rds")
subpops <- merge(mean.pop.sizes[,c("species","2018")], subpops,
                 by.x = "species", by.y = "tax", all.x = TRUE, sort= FALSE)
subpops <- subpops[order(subpops$species),]
table(mean.pop.sizes$species == subpops$species)
subpop.sizes <- vector("list", dim(subpops)[1])
names(subpop.sizes) <- subpops$species
for(i in 1:length(subpop.sizes)) {
  subpop.sizes[[i]] <- rep(subpops$`2018`[i]/subpops$Number.subpops.level.2[i],
                           subpops$Number.subpops.level.2[i])
}
#Testing
table(round(mean.pop.sizes[, which(names(mean.pop.sizes) == 2018)],0) == round(sapply(subpop.sizes, sum),0))

##Running the criterion C for the optimal parameters (GL and p)
system.time(
critC <- criterion_C(x = mean.pop.sizes,
                     assess.year = 2018,
                     project.years = NULL,
                     project = FALSE,
                     ignore.years = c(1718,1748,1778,1793,1808,1818,1823,1838),
                     subcriteria = c("C1", "C2"),
                     generation.time = PopData$GenerationLength.range,
                     prop.mature = PopData$p.est,
                     subpop.size = subpop.sizes,
                     parallel = TRUE,
                     NbeCores = 7)
) ## took 24.3 min with all statistical models, 5130 spp and 7 cores 
saveRDS(critC, "data/criterionC_optmim_params.rds")
critC <- readRDS("data/criterionC_optmim_params.rds")

##Running the criterion C for GL = 25 ys and the optimal p
critC.gl25 <- criterion_C(mean.pop.sizes,
                       assess.year = 2018,
                       project.years = NULL,
                       project = FALSE,
                       ignore.years = c(1718,1748,1778,1793,1808,1818,1823,1838),
                       subcriteria = c("C1","C2"),
                       generation.time = 25,
                       prop.mature = PopData$p.est,
                       subpop.size = subpop.sizes,
                       parallel = TRUE,
                       NbeCores = 7)
saveRDS(critC.gl25, "data/criterionC_GL_25ys.rds")
critC.gl25 <- readRDS("data/criterionC_GL_25ys.rds")

### RUNNING ASSESSMENTS FOR DIFFERENT VALUES OF p ###
ps <- sort(c(1,.85,.72,.60,.49,.51,.58,.31,.25,.33,.64,.45,.35,.28,0.18,0.4), decreasing = TRUE)

# Optimal params
df <- cbind.data.frame(critC[,c("assess.pop.size", "cont.decline")],
                       critC[,grepl("reduction", names(critC))],
                       stringsAsFactors = FALSE)
df$assess.pop.size <- mean.pop.sizes$`2018`
df$assess.pop.size.low <- low.pop.sizes$`2018`
df$assess.pop.size.high <- high.pop.sizes$`2018`

res <- res.low <- res.high <- vector("list",length(ps))
names(res) <- names(res.low) <- names(res.high) <- ps
for(i in 1:length(ps)){
  C1 <- df
  C1$assess.pop.size <- df$assess.pop.size * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res[[i]] <- as.character(all_ranks$ranks_C)

  C1$assess.pop.size <- df$assess.pop.size.low * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.low[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.high * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.high[[i]] <- as.character(all_ranks$ranks_C)
}
res <- do.call(cbind.data.frame, res)
colnames(res) <- paste0("C1.p",colnames(res))
res.low <- do.call(cbind.data.frame, res.low)
colnames(res.low) <- paste0("C1.p",colnames(res.low), ".low")
res.high <- do.call(cbind.data.frame, res.high)
colnames(res.high) <- paste0("C1.p",colnames(res.high), ".high")
critC.all <- cbind.data.frame(critC, res, res.low, res.high)

# GenLeng = 25 years
df <- cbind.data.frame(critC.gl25[,c("assess.pop.size", "cont.decline")],
                       critC.gl25[,grepl("reduction", names(critC.gl25))],
                       stringsAsFactors = FALSE)
df$assess.pop.size <- mean.pop.sizes$`2018`
df$assess.pop.size.low <- low.pop.sizes$`2018`
df$assess.pop.size.high <- high.pop.sizes$`2018`

res <- res.low <- res.high <- vector("list",length(ps))
names(res) <- names(res.low) <- names(res.high) <- ps
for(i in 1:length(ps)){
  C1 <- df
  C1$assess.pop.size <- df$assess.pop.size * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.low * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.low[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.high * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.high[[i]] <- as.character(all_ranks$ranks_C)
}
res <- do.call(cbind.data.frame, res)
colnames(res) <- paste0("C1.p",colnames(res))
res.low <- do.call(cbind.data.frame, res.low)
colnames(res.low) <- paste0("C1.p",colnames(res.low), ".low")
res.high <- do.call(cbind.data.frame, res.high)
colnames(res.high) <- paste0("C1.p",colnames(res.high), ".high")
critC.all.gl25 <- cbind.data.frame(critC.gl25, res, res.low, res.high)


## Calculating the Red List Index for subcriterion A1 and A2
#Optmimal params
all.GL1 <- critC.all[,c(1:15,18:19,16,17,20:67)]
for(i in 18:67) all.GL1[,i] <- as.character(all.GL1[,i])
for(i in 18:67) all.GL1[,i] <- gsub("LC or NT", "LC", all.GL1[,i])

rli.all1 <- apply(all.GL1[,18:67], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL1[,18:67], 2, table)[1:2]

#Gl 25ys
all.GL1.gl25 <- critC.all.gl25[,c(1:15,18:19,16,17,20:67)]
for(i in 18:67) all.GL1.gl25[,i] <- as.character(all.GL1.gl25[,i])
for(i in 18:67) all.GL1.gl25[,i] <- gsub("LC or NT", "LC", all.GL1.gl25[,i])

rli.all1.gl25 <- apply(all.GL1.gl25[,18:67], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL1.gl25[,18:67], 2, table)[1:2]

###################
#### FIGURE SY ####
###################

jpeg(filename = "figures/Figure_SY.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(1,2))
par(mar=c(3,3.5,0.75,0.5), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
#optimum GL
plot(rev(rli.all1[2, grepl("\\.p", colnames(rli.all1))][1:16]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch=19, ylim = c(0.95,1))
#axis(1, at=rev(ps), cex.axis = 1)
axis(1, at=seq(0.2,1,0.1), cex.axis = 1)
#axis(2, at=c(60,80,100,120,140,160,180,200,220), cex.axis = 1)
arrows(x0=rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][1:16]), #using mean CIs
       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][1:16]),
       code = 3, angle = 90, length = 0.05)
#arrows(x0=rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][17:32]), #using low and high CIs
#       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][33:48]),
#       code = 3, angle = 90, length = 0.05, col=2)
legend("topleft", expression(bold(A)), bty="n", cex=1.3)
abline(h=rli.all1[2,1], lty = 2)

#GL: 25ys
plot(rev(rli.all1.gl25[2, grepl("\\.p", colnames(rli.all1.gl25))][1:16]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch=19, ylim = c(0.95,1))
#axis(1, at=rev(ps), cex.axis = 1)
axis(1, at=seq(0.2,1,0.1), cex.axis = 1)
#axis(2, at=c(60,80,100,120,140,160,180,200,220), cex.axis = 1)
arrows(x0=rev(ps), y0 = rev(rli.all1.gl25[1,grepl("\\.p", colnames(rli.all1.gl25))][1:16]), #using mean CIs
       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1.gl25))][1:16]),
       code = 3, angle = 90, length = 0.05)
# arrows(x0=rev(ps), y0 = rev(rli.all1.gl25[1,grepl("\\.p", colnames(rli.all1.gl25))][17:32]), #using low and high CIs
#        y1 = rev(rli.all1.gl25[3,grepl("\\.p", colnames(rli.all1.gl25))][33:48]),
#        code = 3, angle = 90, length = 0.05, col=2)
legend("topleft", expression(bold(B)), bty="n", cex=1.3)
abline(h=rli.all1.gl25[2,1], lty = 2)
legend("topright", c("Group-specific", "Fixed"),
       lty = c(3,0), pch=c(NA,19),
       bty = "n", lwd=2)
dev.off()


#### Saving ####
saveRDS(all.GL1, "data/criterionC_all_prop_mature.rds")
