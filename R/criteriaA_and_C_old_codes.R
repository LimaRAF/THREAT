########################################################
#### ASSESSING POPULATION DECLINE - CRITERIA A AND C ###
########################################################
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


#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA A ####
#############################

## ASSESSMENTS USING SPECIES-SPECIFIC PROXIES OF GENERATION LENGTH ##
critA <- criterion_A(mean.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A1","A2"),
                     generation.time = PopData$GenerationLength.range)
critA.low <- criterion_A(low.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A1","A2"),
                     generation.time = PopData$GenerationLength.range)
critA.high <- criterion_A(high.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A1","A2"),
                     generation.time = PopData$GenerationLength.range)
table(critA$A1)
table(critA.low$A1)
table(critA.high$A1)
critA[critA$A1 %in% "VU" & critA.low$A1 %in% "LC or NT",]
critA[critA$A1 %in% "VU" & critA.high$A1 %in% "EN",]

## ASSESSMENTS of the influence of different GENERATION LENGTHS (FIXED FOR ALL SPECIES) ##
all.GL <- cbind.data.frame(critA,
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 10)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 20)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 25)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 30)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 35)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 40)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 45)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 50)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"),
                                       generation.time = 55)[,c("reduction_A12","A1","A2")],
                           stringsAsFactors = FALSE, deparse.level = 0
                           )
all.GL <- all.GL[,c(1:4,6,7,
                    5,10,13,16,19,22,25,28,31,34,
                    8,11,14,17,20,23,26,29,32,35,
                    9,12,15,18,21,24,27,30,33,36)]
for(i in 17:36) all.GL[,i] <- as.character(all.GL[,i])
for(i in 17:36) all.GL[,i] <- gsub("LC or NT", "LC", all.GL[,i])

## Calculating the Red List Index for subcriterion A1 and A2
rli.all <- apply(all.GL[,17:36], 2, red::rli, boot = TRUE, runs = 4999)

###################
#### FIGURE SX ####
###################

jpeg(filename = "figures/Figure_SX.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
gls = c(10,20,25,30,35,40,45,50,55)
par(mfrow=c(1,2))
par(mar=c(3,3.5,0.75,0.5), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
#subcriterion A1
plot(rli.all[2,grepl("A1\\.", colnames(rli.all))] ~ gls, #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Generation lenght (years)", ylab = "Red List Index", 
     pch=19, ylim = c(0.39,1))
axis(1, at=gls, cex.axis = 1)
#axis(2, at=c(60,80,100,120,140,160,180,200,220), cex.axis = 1)
arrows(x0=gls, y0 = rli.all[1,grepl("A1\\.", colnames(rli.all))], 
       y1 = rli.all[3,grepl("A1\\.", colnames(rli.all))],
       code = 3, angle = 90, length = 0.05)
legend("topright", expression(bold(A)), bty="n", cex=1.3)
abline(h=rli.all[2,1], lty = 2)

#subcriterion A2
plot(rli.all[2,grepl("A2\\.", colnames(rli.all))] ~ gls, #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n", #yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Generation lenght (years)", ylab = "Red List Index", 
     pch=19, ylim = c(0.39,1))
axis(1, at=gls, cex.axis = 1)
arrows(x0=gls, y0 = rli.all[1,grepl("A2\\.", colnames(rli.all))], 
       y1 = rli.all[3,grepl("A2\\.", colnames(rli.all))],
       code = 3, angle = 90, length = 0.05)
legend("topright", expression(bold(B)), bty="n", cex=1.3)
abline(h=rli.all[2,11], lty = 2)
legend("top", c("Group-specific", "Fixed"),
       lty = c(2,0), pch=c(NA,19),
       bty = "n", lwd=2)
dev.off()

#Where are the changes related to low generation lengths (25 years)
critA.gl25 <- criterion_A(mean.pop.sizes,
                     assess.year = 2018,
                     project.years = NULL,
                     subcriteria = c("A1","A2"),
                     generation.time = 25)

critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "EN",] # 56 cases, 21 for 35 ys
critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "VU",] # No cases
critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "LC or NT",] # No cases
critA.gl25[critA$A1 %in% "EN" & critA.gl25$A1 %in% "VU",] #390 cases, 109 for 35 ys
critA.gl25[critA$A1 %in% "EN" & critA.gl25$A1 %in% "LC or NT",] #1 cases, 0 for 35 ys
critA.gl25[critA$A1 %in% "VU" & critA.gl25$A1 %in% "LC or NT",] #393 cases, 96 for 35 ys

#Inpecting some classic examples
critA[critA$species %in% "Araucaria angustifolia",]
critA.gl25[critA$species %in% "Araucaria angustifolia",]
critA[critA$species %in% "Euterpe edulis",]
critA.gl25[critA$species %in% "Euterpe edulis",]
critA[critA$species %in% "Cedrela fissilis",]
critA.gl25[critA$species %in% "Cedrela fissilis",]
critA[critA$species %in% "Cecropia pachystachya",]
critA.gl25[critA$species %in% "Cecropia pachystachya",]

#### CHORD DIAGRAM ####
mat <- as.matrix(table(paste0(critA.gl25$A1,"_25"), paste0(critA$A1,"_opt")))
mat <- mat[c(1,2,4,3), c(3,4,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
grid.col = c(CR_25 = "red", EN_25 = "darkorange", VU_25 = "gold", LCorNT_25 = "khaki",
             CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
col_mat[mat >= 10] = adjustcolor(col_mat[mat >= 10], alpha.f = 0.5)
col_mat[mat < 10] = adjustcolor(col_mat[mat < 10], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 10] = mat[mat < 10]*2
mat[mat > 0 & mat < 5] = 10

#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible) = FALSE
chordDiagram(mat, big.gap = 20, annotationTrack = "grid", annotationTrackHeight = mm_h(4),
             grid.col = grid.col, col = col_mat,
             self.link = 1, link.visible = visible,
             #h=0.9,
             #w=1,
             #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
             link.lwd = 4
             #h.ratio = 0.7
             #reduce_to_mid_line = FALSE,
             #w2=0.5,
             #rou=0.2
             #point1 = rep(0,16)
             )
#Putting legends on
sec.ind <- c("CR","EN","VU","LC or NT","LC or NT","VU","EN","CR")
for(si in get.all.sector.index()) {
  lab <- sec.ind[which(si == get.all.sector.index())]
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  if(si == "VU_25") {
    circos.text(mean(xlim), mean(ylim), labels = "VU", 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = FALSE, col = "black")
    
  } else {
    circos.text(mean(xlim), mean(ylim), labels = lab, 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = TRUE, col = "black")
  }  
}
legend("topleft","Optimal gen. lenght", bty="n")
legend("topright","Gen. lenght = 25 years", bty="n")

## Renaming columns
names(all.GL)[grepl("\\.1$", names(all.GL))] <- gsub("\\.1$", ".10ys", names(all.GL)[grepl("\\.1$", names(all.GL))])
names(all.GL)[grepl("\\.2$", names(all.GL))] <- gsub("\\.2$", ".20ys", names(all.GL)[grepl("\\.2$", names(all.GL))])
names(all.GL)[grepl("\\.3$", names(all.GL))] <- gsub("\\.3$", ".25ys", names(all.GL)[grepl("\\.3$", names(all.GL))])
names(all.GL)[grepl("\\.4$", names(all.GL))] <- gsub("\\.4$", ".30ys", names(all.GL)[grepl("\\.4$", names(all.GL))])
names(all.GL)[grepl("\\.5$", names(all.GL))] <- gsub("\\.5$", ".35ys", names(all.GL)[grepl("\\.5$", names(all.GL))])
names(all.GL)[grepl("\\.6$", names(all.GL))] <- gsub("\\.6$", ".40ys", names(all.GL)[grepl("\\.6$", names(all.GL))])
names(all.GL)[grepl("\\.7$", names(all.GL))] <- gsub("\\.7$", ".45ys", names(all.GL)[grepl("\\.7$", names(all.GL))])
names(all.GL)[grepl("\\.8$", names(all.GL))] <- gsub("\\.8$", ".50ys", names(all.GL)[grepl("\\.8$", names(all.GL))])
names(all.GL)[grepl("\\.9$", names(all.GL))] <- gsub("\\.9$", ".55ys", names(all.GL)[grepl("\\.9$", names(all.GL))])

## Saving
saveRDS(all.GL, "data/criterionA_all_GLs.rds")


#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA C ####
#############################

#### NEED TO RUN FOR low and high POP SIZES?

##Running the criterion C for the optimal parameters (GL and p)
system.time(
critC <- criterion_C(mean.pop.sizes,
                     assess.year = 2018,
                     project.years = NULL,
                     project = FALSE,
                     ignore.years = c(1718,1748,1778,1793,1808,1818,1823,1838),
                     subcriteria = c("C1"),
                     generation.time = PopData$GenerationLength.range,
                     prop.mature = PopData$p.est,
                     parallel = TRUE,
                     NbeCores = 7)
) ## took 18.3 min with all statistical models, 5130 spp and 7 cores 
#saveRDS(critC, "data/criterionC_optmim_params.rds")

##Running the criterion C for GL = 25 ys and the optimal p
critC.gl25 <- criterion_C(mean.pop.sizes,
                       assess.year = 2018,
                       project.years = NULL,
                       project = FALSE,
                       ignore.years = c(1718,1748,1778,1793,1808,1818,1823,1838),
                       subcriteria = c("C1"),
                       generation.time = 25,
                       prop.mature = PopData$p.est,
                       parallel = TRUE,
                       NbeCores = 7)
#saveRDS(critC.gl25, "data/criterionC_GL_25ys.rds")

### RUNNINGG ASSESSMENTS FOR DIFFERENT VALUES OF p ###
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

# Optimal params
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
all.GL1 <- critC.all[,c(1:11,13:14,12,15:62)]
for(i in 14:62) all.GL1[,i] <- as.character(all.GL1[,i])
for(i in 14:62) all.GL1[,i] <- gsub("LC or NT", "LC", all.GL1[,i])

rli.all1 <- apply(all.GL1[,14:62], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL1[,14:62], 2, table)[1:2]

#Gl 25ys
all.GL1.gl25 <- critC.all.gl25[,c(1:11,13:14,12,15:62)]
for(i in 14:62) all.GL1.gl25[,i] <- as.character(all.GL1.gl25[,i])
for(i in 14:62) all.GL1.gl25[,i] <- gsub("LC or NT", "LC", all.GL1.gl25[,i])

rli.all1.gl25 <- apply(all.GL1.gl25[,14:62], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL1.gl25[,14:62], 2, table)[1:2]

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


## Saving
saveRDS(all.GL1, "data/criterionC_all_prop_mature.rds")


#########################################################################################################################################################
#########################################################################################################################################################
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
apply(all.GL2[,9:25], 2, table)[1:2]

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

#### COMPARING MAT. POP SIZE AND SPECIES DISTRIBUTION ####
df1$mature.pop.size <- df1$pop.size * df1$p

par(mfrow = c(2,2))
par(mar=c(4,4,0.75,0.5), mgp=c(2.5,0.25,0),tcl=-0.2,las=1)
plot(log(mature.pop.size) ~ log(treeco.occs), data = df1) # as expected, since one was generated from the other
abline(lm(log(mature.pop.size) ~ log(treeco.occs), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(non.dup.occs), data = df1)
abline(lm(log(mature.pop.size) ~ log(non.dup.occs), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(EOO.level.1), data = df1)
abline(lm(log(mature.pop.size) ~ log(EOO.level.1), data = df1), lwd=2,col=2)
plot(log(mature.pop.size) ~ log(AOO.level.1), data = df1)
abline(lm(log(mature.pop.size) ~ log(AOO.level.1), data = df1), lwd=2,col=2)

# number of occurences, non-duplicated occurrences, n.localities, EOO or AOO? AOO!
pairs(log(mature.pop.size) ~ log(non.dup.occs) + 
        log(EOO.level.1) + log(AOO.level.1) + log(Loc.level1), data = df1)
bbmle::AICtab(lm(log(mature.pop.size) ~ log(AOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(non.dup.occs), data = df1),
              lm(log(mature.pop.size) ~ log(total.occs), data = df1),
              lm(log(mature.pop.size) ~ log(Loc.level1), data = df1),
              lm(log(mature.pop.size) ~ log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(AOO.level.1) + log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(AOO.level.1) * log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(Loc.level1) * log(EOO.level.1), data = df1),
              lm(log(mature.pop.size) ~ log(non.dup.occs) * log(EOO.level.1), data = df1))

df2 <- na.omit(df1[,c("mature.pop.size","EOO.level.1","AOO.level.1")])
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

## Takink into account unequal variance
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
bbmle::AICtab(mod.ols, mod.rob, mod.gls, 
              mod.gls.pow, mod.gls.exp, mod.gls.fix, mod.gls.fix1)
summary(mod.rob)$sigma
summary(mod.gls)$sigma
summary(mod.ols)$sigma
summary(mod.gls.fix)$sigma
summary(mod.gls.fix1)$sigma
summary(mod.gls.pow)$sigma
summary(mod.gls.exp)$sigma

car::Anova(mod.gls.fix1)
plot(mod)
car::avPlots(mod)
summary(mod.gls.fix)
summary(mod.gls.fix1)
