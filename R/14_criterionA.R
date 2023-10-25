###################################################
###################################################
#### ASSESSING POPULATION DECLINE - CRITERIA A ####
###################################################
###################################################
rm(list=ls())

#### LOADING PACKAGES ###
require("ConR")
require("red")
require("circlize")


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
hab <- readRDS("data/threat_habitats.rds")
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
timber[PopData$species %in% "Euterpe edulis"] <- 10 # non-timber heavily exploited species
timber[PopData$species %in% "Dicksonia sellowiana"] <- 10 # non-timber heavily exploited species
timber <- as.numeric(timber)

#### GETTING THREAT ECOLOGICAL GROUP INFORMATION  ####
hab1 <- hab[match(PopData$species, hab$Name_submitted), ]
table(hab1$Name_submitted == PopData$species)

## Taking into account habitat change in the population decline of early-successional species
early.sucession <- hab1$ecol.group
early.sucession[early.sucession %in% c("late_secondary", "climax", "unknown")] <- 1L
early.sucession[early.sucession %in% c("early_secondary")] <- 0.9
early.sucession[early.sucession %in% c("pioneer")] <- 0.8
early.sucession <- as.numeric(early.sucession)

## Taking into account non-forested habitats for "ruderal" small-statured pioneers
ruderals <- grepl("Antr", hab1$vegetation.type.reflora) & 
  hab1$GF %in% c("large_shrub", "small_tree") &
  hab1$ecol.group %in% "pioneer"
early.sucession[ruderals] <- 0.5


#########################################################################################################################################################
#########################################################################################################################################################
#############################
#### APPLYING CRITERIA A ####
#############################

## ASSESSMENTS USING SPECIES-SPECIFIC PROXIES OF GENERATION LENGTH ##
critA <- criterion_A(mean.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A2"),
                     generation.time = PopData$GL,
                     exploitation = timber,
                     correction = early.sucession)
critA.low <- criterion_A(low.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A2"),
                     generation.time = PopData$GL,
                     exploitation = timber,
                     correction = early.sucession)
names(critA.low) <- paste(names(critA.low), ".low", sep= "")
critA.high <- criterion_A(high.pop.sizes, assess.year = 2018,
                     project.years = NULL, subcriteria = c("A2"),
                     generation.time = PopData$GL,
                     exploitation = timber,
                     correction = early.sucession)
names(critA.high) <- paste(names(critA.high), ".high", sep= "")

critA[critA$category_A %in% "VU" & critA.low$category_A.low %in% "LC or NT",]
critA[critA$category_A %in% "VU" & critA.high$category_A.high %in% "EN",]
critA[critA$category_A %in% "EN" & critA.high$category_A.high %in% "CR",]

## ASSESSMENTS of the influence of different GENERATION LENGTHS (FIXED FOR ALL SPECIES) ##
all.GL <- cbind.data.frame(critA, 
                           critA.low[, c("category_A.low"), drop = FALSE],
                           critA.high[, c("category_A.high"), drop = FALSE],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 10)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 15)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 20)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 25)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 30)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 35)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 40)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 45)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 50)[,c("reduction_A12","category_A")],
                           criterion_A(mean.pop.sizes, assess.year = 2018, subcriteria = c("A2"), exploitation = timber, correction = early.sucession,
                                       generation.time = 55)[,c("reduction_A12","category_A")],
                           stringsAsFactors = FALSE, deparse.level = 0
                           )

#### REVER AQUI APÓS A INCLUSÃO DOS 15 ANOS DE GL ####
all.GL <- all.GL[,c(1:12,
                    13,15,17,19,21,23,25,27,29,31, #Reductions A12
                    14,16,18,20,22,24,26,28,30,32 #A2
                    )] 
for(i in 23:32) all.GL[,i] <- as.character(all.GL[,i])
for(i in 23:32) all.GL[,i] <- gsub("LC or NT", "LC", all.GL[,i])

## Calculating the Red List Index for subcriterion A2
rli.all <- apply(all.GL[,23:32], 2, red::rli, boot = TRUE, runs = 4999)
rli.all.opt <- red::rli(gsub("LC or NT", "LC", as.character(all.GL[,"category_A"])), boot = TRUE, runs = 4999)
rli.all.opt.low <- red::rli(gsub("LC or NT", "LC", as.character(all.GL[,"category_A.low"])), boot = TRUE, runs = 4999)
rli.all.opt.high <- red::rli(gsub("LC or NT", "LC", as.character(all.GL[,"category_A.high"])), boot = TRUE, runs = 4999)


#Where are the changes related to low generation lengths (25 years)
critA.gl25 <- criterion_A(mean.pop.sizes,
                     assess.year = 2018,
                     project.years = NULL,
                     subcriteria = c("A2"),
                     generation.time = 25,
                     exploitation = timber, correction = early.sucession)

critA.gl25[critA$category_A %in% "CR" & critA.gl25$category_A %in% "EN",] # 66 cases, 21 for 35 ys; revision 127 cases
critA.gl25[critA$category_A %in% "CR" & critA.gl25$category_A %in% "VU",] # No cases; revision No cases
critA.gl25[critA$category_A %in% "CR" & critA.gl25$category_A %in% "LC or NT",] # No cases; revision No cases
critA.gl25[critA$category_A %in% "EN" & critA.gl25$category_A %in% "VU",] #383 cases, 109 for 35 ys; revision 335 cases 
critA.gl25[critA$category_A %in% "EN" & critA.gl25$category_A %in% "LC or NT",] #1 cases, 0 for 35 ys; revision No cases
critA.gl25[critA$category_A %in% "VU" & critA.gl25$category_A %in% "LC or NT",] #393 cases, 96 for 35 ys; revision 101 cases

#Inspecting some classic examples
critA[critA$species %in% "Araucaria angustifolia",]
critA.gl25[critA$species %in% "Araucaria angustifolia",] # EN para VU
critA[critA$species %in% "Euterpe edulis",]
critA.gl25[critA$species %in% "Euterpe edulis",] # no change
critA[critA$species %in% "Cedrela fissilis",]
critA.gl25[critA$species %in% "Cedrela fissilis",] # no change
critA[critA$species %in% "Cecropia pachystachya",]
critA.gl25[critA$species %in% "Cecropia pachystachya",] # no change

#Where are the changes related to low generation lengths (20 years)
critA.gl20 <- criterion_A(mean.pop.sizes,
                          assess.year = 2018,
                          project.years = NULL,
                          subcriteria = c("A2"),
                          generation.time = 20,
                          exploitation = timber, correction = early.sucession)

critA.gl20[critA$category_A %in% "CR" & critA.gl20$category_A %in% "EN",] # revision 202 cases
critA.gl20[critA$category_A %in% "CR" & critA.gl20$category_A %in% "VU",] # revision No cases
critA.gl20[critA$category_A %in% "CR" & critA.gl20$category_A %in% "LC or NT",] # revision No cases
critA.gl20[critA$category_A %in% "EN" & critA.gl20$category_A %in% "VU",] #revision 661 cases 
critA.gl20[critA$category_A %in% "EN" & critA.gl20$category_A %in% "LC or NT",] #revision 18 cases
critA.gl20[critA$category_A %in% "VU" & critA.gl20$category_A %in% "LC or NT",] #revision 231 cases


## Renaming columns
names(all.GL)[grepl("\\.1$", names(all.GL))] <- gsub("\\.1$", ".10ys", names(all.GL)[grepl("\\.1$", names(all.GL))])
names(all.GL)[grepl("\\.2$", names(all.GL))] <- gsub("\\.2$", ".15ys", names(all.GL)[grepl("\\.2$", names(all.GL))])
names(all.GL)[grepl("\\.3$", names(all.GL))] <- gsub("\\.3$", ".20ys", names(all.GL)[grepl("\\.3$", names(all.GL))])
names(all.GL)[grepl("\\.4$", names(all.GL))] <- gsub("\\.4$", ".25ys", names(all.GL)[grepl("\\.4$", names(all.GL))])
names(all.GL)[grepl("\\.5$", names(all.GL))] <- gsub("\\.5$", ".30ys", names(all.GL)[grepl("\\.5$", names(all.GL))])
names(all.GL)[grepl("\\.6$", names(all.GL))] <- gsub("\\.6$", ".35ys", names(all.GL)[grepl("\\.6$", names(all.GL))])
names(all.GL)[grepl("\\.7$", names(all.GL))] <- gsub("\\.7$", ".40ys", names(all.GL)[grepl("\\.7$", names(all.GL))])
names(all.GL)[grepl("\\.8$", names(all.GL))] <- gsub("\\.8$", ".45ys", names(all.GL)[grepl("\\.8$", names(all.GL))])
names(all.GL)[grepl("\\.9$", names(all.GL))] <- gsub("\\.9$", ".50ys", names(all.GL)[grepl("\\.9$", names(all.GL))])
names(all.GL)[grepl("\\.10$", names(all.GL))] <- gsub("\\.10$", ".55ys", names(all.GL)[grepl("\\.10$", names(all.GL))])
names(all.GL)[grepl("category_A\\.", names(all.GL))] <- gsub("category_A", "A2", names(all.GL)[grepl("category_A\\.", names(all.GL))])

## Renaming the LC category
all.GL[] <- lapply(all.GL, gsub, pattern = "^LC$", replacement = "LC or NT")

#### Saving ####
saveRDS(all.GL, "data/criterionA_all_GLs.rds")


#### CHORD DIAGRAMS ####

# SUBCRITERIA A1
# jpeg(filename = "figures/Figure_SW.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
#      res = 300, family = "sans", type="cairo", bg="white")
# 
# mat <- as.matrix(table(paste0(critA.gl25$A1,"_25"), paste0(critA$A1,"_opt")))
# mat <- mat[c(1,2,4,3), c(3,4,2,1)]
# colnames(mat) <- gsub(" ", "", colnames(mat))
# rownames(mat) <- gsub(" ", "", rownames(mat))
# 
# #Defining the colors of tracks and links
# grid.col = c(CR_25 = "red", EN_25 = "darkorange", VU_25 = "gold", LCorNT_25 = "khaki",
#              CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
# col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
# col_mat[mat >= 10] = adjustcolor(col_mat[mat >= 10], alpha.f = 0.5)
# col_mat[mat < 10] = adjustcolor(col_mat[mat < 10], alpha.f = 0.9)
# #col_mat[mat < 5] = "#00000000"
# mat[mat < 10] = mat[mat < 10]*2
# mat[mat > 0 & mat < 5] = 10
# 
# #plotting the diagram
# par(mfrow=c(1,1))
# par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
# circos.clear()
# circos.par(start.degree = 90)
# visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
# #diag(visible) = FALSE
# lava::revdiag(visible) = FALSE
# chordDiagram(mat, big.gap = 15, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
#              grid.col = grid.col, col = col_mat,
#              self.link = 1, link.visible = visible,
#              #h=0.9,
#              #w=1,
#              #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
#              link.lwd = 4
#              #h.ratio = 0.7
#              #reduce_to_mid_line = FALSE,
#              #w2=0.5,
#              #rou=0.2
#              #point1 = rep(0,16)
# )
# #Putting legends on
# sec.ind <- c("CR","EN","VU","LC or NT","LC or NT","VU","EN","CR")
# for(si in get.all.sector.index()) {
#   lab <- sec.ind[which(si == get.all.sector.index())]
#   xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#   ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#   if(si == "VU_25") {
#     circos.text(mean(xlim), mean(ylim), labels = "VU", 
#                 sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
#                 facing = "bending", niceFacing = FALSE, col = "black")
#     
#   } else {
#     circos.text(mean(xlim), mean(ylim), labels = lab, 
#                 sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
#                 facing = "bending", niceFacing = TRUE, col = "black")
#   }  
# }
# legend("topleft","Optimal gen. lenght", bty="n", cex=1.2)
# legend("topright","Gen. lenght = 25 years", bty="n", cex=1.2)
# dev.off()


# SUBCRITERIA A2 - comparison with 25 years
jpeg(filename = "figures/Figure_SWb_25y.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

mat <- as.matrix(table(paste0(critA.gl25$category_A,"_25"), paste0(critA$category_A,"_opt")))
mat <- mat[c(1,2,4,3), c(3,4,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
# grid.col = c(CR_25 = "red", EN_25 = "darkorange", VU_25 = "gold", LCorNT_25 = "khaki",
#              CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
# col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
grid.col = c(CR_25 = "#CD0000", EN_25 = "#E69F00", VU_25 = "#FFFF00", LCorNT_25 = "khaki",
             CR_opt = "#CD0000", EN_opt = "#E69F00", VU_opt = "#FFFF00", LCorNT_opt = "khaki")
col_mat = rep(rev(c("#CD0000","#E69F00","#FFFF00","khaki")), each=4)
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
legend("topleft","Group-specific gen. lenght", bty="n", cex=1.2)
legend("topright","Gen. lenght = 25 years", bty="n", cex=1.2)
dev.off()


# SUBCRITERIA A2 - comparison with 20 years
jpeg(filename = "figures/Figure_SWb_20y.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

mat <- as.matrix(table(paste0(critA.gl20$category_A,"_20"), paste0(critA$category_A,"_opt")))
mat <- mat[c(1,2,4,3), c(3,4,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
# grid.col = c(CR_20= "red", EN_20 = "darkorange", VU_20 = "gold", LCorNT_20 = "khaki",
#              CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
# col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
grid.col = c(CR_20= "#CD0000", EN_20 = "#E69F00", VU_20 = "#FFFF00", LCorNT_20 = "khaki",
             CR_opt = "#CD0000", EN_opt = "#E69F00", VU_opt = "#FFFF00", LCorNT_opt = "khaki")
col_mat = rep(rev(c("#CD0000","#E69F00","#FFFF00","khaki")), each=4)
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
legend("topleft","Group-specific gen. lenght", bty="n", cex=1.2)
legend("topright","Gen. lenght = 20 years", bty="n", cex=1.2)
dev.off()


#####################################################################
#### NEW VERSION OF FIGURE SX USING RLI AND % OF THREAT PER CATS ####
#####################################################################

jpeg(filename = "figures/Figure_SXab.jpg", width = 3750, height = 2000, 
     units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
gls = c(10,15,20,25,30,35,40,45,50,55)
par(mfrow=c(1,2))
par(mar=c(3,3.5,0.75,0.5), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)

# Panel A: RLI results
#subcriterion A2
plot(rli.all[2,] ~ gls, #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n", #yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Generation lenght (years)", ylab = "Red List Index", 
     pch=19, ylim = c(0.39,1))
axis(1, at=gls, cex.axis = 1)
arrows(x0=gls, y0 = rli.all[1,], y1 = rli.all[3,],
       code = 3, angle = 90, length = 0.05)
legend("topleft", expression(bold(A)), bty="n", cex=1.3) #, x.intersp=-0.1,y.intersp=-0.1)
# abline(h=rli.all[2,11], lty = 2)
abline(h=rli.all.opt[2], lty = 2)
# abline(h=rli.all.opt[c(1,3)], lty = 3, col = "grey")
legend("topright", c("Group-specific", "Fixed values"),
       lty = c(2,0), pch=c(NA,19),
       bty = "n", lwd=2)

# Panel B: results per threat category
pct <- t(rbind(table(all.GL$A2.10ys),
               table(all.GL$A2.15ys),
               table(all.GL$A2.20ys),
               table(all.GL$A2.25ys),
               table(all.GL$A2.30ys),
               table(all.GL$A2.35ys),
               table(all.GL$A2.40ys),
               table(all.GL$A2.45ys),
               table(all.GL$A2.50ys),
               table(all.GL$A2.55ys)))
pct <- pct[match(row.names(pct), c("CR", "EN", "VU", "LC or NT")),]
gl.vals <- c("10","15","20","25","30","35","40","45","50","55")
colnames(pct) <- gl.vals
pct <- round(100 * pct/dim(all.GL)[1], 1)

par(mar=c(3,3.5,0.75, 0.25), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
barplot(pct, horiz = FALSE, 
        # col = c("red","darkorange","gold","khaki"),
        col = c("#CD0000","#E69F00","#FFFF00","khaki"),
        legend.text = c("CR", "EN", "VU", "LC or NT"),
        xlab = "Generation length (years)",  
        args.legend = list(ncol=4, x.intersp = 1, y.intersp = 1, yjust = -1,
                           adj = c(0, 0.5)), 
        #xaxt = "n", 
        cex.lab = 1.2,
        ylab = "IUCN Red List Categories (%)", 
)
axis(1, at=seq(0.725, 12.7, 1.2), cex.axis = 1, labels = FALSE)
legend("topleft",legend=expression(bold("B")),
       bty="n",horiz=F,cex=1.5,y.intersp=0.4,x.intersp=0.4)

par(xpd = TRUE)
legend(10.65,12, expression(bold(CR)), bty="n", cex=0.9)
legend(10.65,55, expression(bold(EN)), bty="n", cex=0.9)
legend(10.65,92.5, expression(bold(VU)), bty="n", cex=0.9)
legend(10.27,101.7, bquote(bold('LC+NT')), bty="n", cex=0.9)
par(xpd = FALSE)
dev.off()


# pct <- t(rbind(table(all.GL$A2.10ys),
#                table(all.GL$A2.15ys),
#                table(all.GL$A2.20ys),
#                table(all.GL$A2.25ys),
#                table(all.GL$A2.30ys),
#                table(all.GL$A2.35ys),
#                table(all.GL$A2.40ys),
#                table(all.GL$A2.45ys),
#                table(all.GL$A2.50ys),
#                table(all.GL$A2.55ys)))
# pct <- pct[match(row.names(pct), c("CR", "EN", "VU", "LC or NT")),]
# gl.vals <- c("10","15","20","25","30","35","40","45","50","55")
# colnames(pct) <- gl.vals
# pct <- round(100 * pct/dim(all.GL)[1], 1)
# 
# jpeg(filename = "figures/Figure_SXb.jpg", width = 2500, height = 2000, units = "px", pointsize = 12,
#      res = 300, family = "sans", type="cairo", bg="white")
# par(mfrow=c(1,1))
# par(mar=c(3.5,3.5,0.75, 2), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
# barplot(pct, horiz = FALSE, col = c("red","darkorange","gold","khaki"),
#         legend.text = c("CR", "EN", "VU", "LC or NT"),
#         xlab = "Generation length (years)",  
#         args.legend = list(ncol=4, x.intersp = 1, y.intersp = 1, yjust = -1,
#                            adj = c(0, 0.5)), 
#         #xaxt = "n", 
#         cex.lab = 1.2,
#         ylab = "IUCN Red List Categories (%)", 
# )
# axis(1, at=seq(0.725, 12.7, 1.2), cex.axis = 1, labels = FALSE)
# par(xpd = TRUE)
# legend(10.75,12, expression(bold(CR)), bty="n", cex=1)
# legend(10.75,55, expression(bold(EN)), bty="n", cex=1)
# legend(10.75,92.5, expression(bold(VU)), bty="n", cex=1)
# legend(10.485,101.75, bquote(bold('LC+NT')), bty="n", cex=1)
# par(xpd = FALSE)
# dev.off()
