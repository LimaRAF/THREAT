###################################################
###################################################
#### ASSESSING POPULATION DECLINE - CRITERIA A ####
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
timber <- !is.na(PopData$times.cites)
timber <- ifelse(timber == FALSE, 0, 10)
timber[!is.na(PopData$commercial) & PopData$commercial == FALSE] <- 5
timber[PopData$species %in% "Euterpe edulis"] <- 10

#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA A ####
#############################

## ASSESSMENTS USING SPECIES-SPECIFIC PROXIES OF GENERATION LENGTH ##
critA <- criterion_A(mean.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A1","A2"),
                     generation.time = PopData$GenerationLength.range,
                     exploitation = timber)
critA.low <- criterion_A(low.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A1","A2"),
                     generation.time = PopData$GenerationLength.range,
                     exploitation = timber)
critA.high <- criterion_A(high.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A1","A2"),
                     generation.time = PopData$GenerationLength.range,
                     exploitation = timber)
table(critA$A1)
table(critA.low$A1)
table(critA.high$A1)
critA[critA$A1 %in% "VU" & critA.low$A1 %in% "LC or NT",]
critA[critA$A1 %in% "VU" & critA.high$A1 %in% "EN",]

## ASSESSMENTS of the influence of different GENERATION LENGTHS (FIXED FOR ALL SPECIES) ##
all.GL <- cbind.data.frame(critA,
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 10)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 20)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 25)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 30)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 35)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 40)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 45)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 50)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A1","A2"), exploitation = timber,
                                       generation.time = 55)[,c("reduction_A12","A1","A2")],
                           stringsAsFactors = FALSE, deparse.level = 0
                           )
all.GL <- all.GL[,c(1:4,6,7,10,
                    5,11,14,17,20,23,26,29,32,35, #Reductions A12
                    8,12,15,18,21,24,27,30,33,36, #A1
                    9,13,16,19,22,25,28,31,34,37)] #A2
for(i in 18:37) all.GL[,i] <- as.character(all.GL[,i])
for(i in 18:37) all.GL[,i] <- gsub("LC or NT", "LC", all.GL[,i])

## Calculating the Red List Index for subcriterion A1 and A2
rli.all <- apply(all.GL[,18:37], 2, red::rli, boot = TRUE, runs = 4999)

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
                     generation.time = 25,
                     exploitation = timber)

critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "EN",] # 66 cases, 21 for 35 ys
critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "VU",] # No cases
critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "LC or NT",] # No cases
critA.gl25[critA$A1 %in% "EN" & critA.gl25$A1 %in% "VU",] #383 cases, 109 for 35 ys
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

## Renaming the LC category
all.GL[] <- lapply(all.GL, gsub, pattern = "^LC$", replacement = "LC or NT")

#### Saving ####
saveRDS(all.GL, "data/criterionA_all_GLs.rds")


#### CHORD DIAGRAMS ####
require(circlize)

# SUBCRITERIA A1
jpeg(filename = "figures/Figure_SW.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

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
par(mfrow=c(1,1))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible) = FALSE
chordDiagram(mat, big.gap = 15, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
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
legend("topleft","Optimal gen. lenght", bty="n", cex=1.2)
legend("topright","Gen. lenght = 25 years", bty="n", cex=1.2)
dev.off()


# SUBCRITERIA A2
jpeg(filename = "figures/Figure_SWb.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

mat <- as.matrix(table(paste0(critA.gl25$A2,"_25"), paste0(critA$A2,"_opt")))
mat <- mat[c(1,2,4,3), c(3,4,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
grid.col = c(CR_25 = "red", EN_25 = "darkorange", VU_25 = "gold", LCorNT_25 = "khaki",
             CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 15] = mat[mat < 15]*2
mat[mat > 0 & mat < 5] = 10

#plotting the diagram
par(mfrow=c(1,1))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible) = FALSE
chordDiagram(mat, big.gap = 15, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
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
                facing = "bending", niceFacing = TRUE, col = "black")
    
  } else {
    circos.text(mean(xlim), mean(ylim), labels = lab, 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = TRUE, col = "black")
  }  
}
legend("topleft","Optimal gen. lenght", bty="n", cex=1.2)
legend("topright","Gen. lenght = 25 years", bty="n", cex=1.2)
dev.off()
