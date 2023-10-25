##################################################################################################################################################################################################################################H
##################################################################################################################################################################################################################################H
#################
#### FIGURES ####
#################
### NOT ALL FIGURES ARE BEING GENERATED HERE
### MAPS ARE BEING GENERATED IN SCRIPT '22_maps.R'
### SOME OTHER FIGURES ARE BEING GENERATED WITHIN THE SCRIPTS...
rm(list = ls())
gc()

## Loading packages
require(circlize)
require(webr)
require(ggplot2)
library(cowplot)
require(multcompView)
require(lava)
require(car)
source("R/99_functions.R")

## Loading the data
all.crit <- readRDS("data/all.criteria.rds")

## Defining the colors for each category
# cores <- c(CR_PE = "darkred", CR = "red", EN = "darkorange", VU = "gold", 
#            NT = "yellowgreen", `NA` = "grey", DD = "grey", LC = "forestgreen")
cores <- c(CR_PE = "#8B0000", CR = "#CD0000", EN = "#E69F00", VU = "#F0E442", 
                NT = "#C1FF73", `NA` = "#A0A0A0", DD = "#A0A0A0", LC = "#10E210")
crit_cores <- c("#A0A0A0", "#648FFF", "#2455D8", "#88CCEE", "cyan", "#E280DC", "#EE72AF", "#785EF0", "black")

crit_cores <- grey.colors(9, start = 0.3, end = 0.9)


color.f <- colorRampPalette(c("cyan","blue"))
color.f <- colorRampPalette(c("cyan","#8000FF"))
crit_cores <- color.f(9)
crit_cores <- c("#A0A0A0", crit_cores[-1])

#### FIGURE 1 ####
#Creating and organizing the data frame to be plotted
pie.df <- all.crit[,c("cat.reg.clean","main.criteria","endemic")]
names(pie.df)[1] <- "category"
pie.df$category[pie.df$category %in% "CR_PE"] <- "CR" 
pie.df$main.criteria[pie.df$main.criteria %in% "A2+B1+B2+C2+D"] <- "All" 
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+B2")] <- "A2, B1+2"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+B1+B2")] <- "A2, B1+2"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+B1+B2+C2","A2+B2+C2")] <- "other" #"A+B+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+C1","A2+C1+C2","A2+C2")] <- "other" #"A+C"
# pie.df$main.criteria[pie.df$main.criteria %in% c("B1")] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+C2+D")] <- "other" #"A+C+D"
pie.df$main.criteria[is.na(pie.df$main.criteria)] <- "other"
# pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2")] <- "A2"
pie.df$main.criteria[pie.df$category %in% "NA"] <- ""
pie.df <- pie.df[order(match(pie.df$category, names(cores)[-1])),]
pie.df <- pie.df[!pie.df$category %in% "DD",]
pie.df$main.criteria[pie.df$main.criteria %in% c("B1+B2")] <- "B1+2"

## Geral
dados <- pie.df
toto <- rep(TRUE, length(table(dados$category, dados$main.criteria)))
horiz <- c(8, 25, 26)  
toto[horiz] <- FALSE
pie.all <- my.PieDonut(dados, aes(category, main.criteria),
                       ratioByGroup=FALSE, showPieName = FALSE,
                       start=75,r0=0, r1=1,r2=1.3,
                       showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
                       main.colors = cores[c(2:5,8,6)],
                       #donut.colors = NULL,
                       pieAlpha = 0.9,
                       correct.legend.x = c(0.00,0.00,-0.02,-0.075,-0.15,-0.10),
                       correct.legend.y = c(0.2,0.00,-0.200,-0.075, 0.00, 0.10),
                       tidy.legend.donut = toto,
                       tidy.legend.pie = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                       plot.piedonut = FALSE,
                       return = "pie", pieLabelSize = 6, donutLabelSize = 4.5)
                       # return = "pie", pieLabelSize = 6, donutLabelSize = 4.5,
                       # explode = 1:6,
                       # explode.cor.pie = 0.5, explode.cor.donut = 0.7)
crit_cores1 <-  
donut.all <- my.PieDonut(dados, aes(category, main.criteria, colour = main.criteria),
                         ratioByGroup=FALSE, 
                         showPieName = FALSE,
                         start=75, r0=0, r1=1, r2=1.3,
                         showRatioThreshold = 0.01, 
                         labelpositionThreshold =  0.028,
                         main.colors = cores[c(2:5,8,6)],
                         sub.colors = crit_cores,
                         #donut.colors = NULL,
                         pieAlpha = 0.75,
                         correct.legend.x = c(0.00,0.00,-0.045,-0.100,-0.15,-0.10),
                         correct.legend.y = c(0.2,0.00,-0.200,-0.075, 0.00, 0.10),
                         tidy.legend.donut = toto,
                         tidy.legend.pie = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                         plot.piedonut = FALSE,
                         return = "donut", pieLabelSize = 6, donutLabelSize = 4.5)
                         # return = "donut", pieLabelSize = 6, donutLabelSize = 4.5,
                         # explode = 1:6,
                         # explode.cor.pie = 0.5, explode.cor.donut = 0.7)
ggdraw(pie.all) + draw_plot(donut.all)

## Endemicas
dados <- pie.df[pie.df$endemic %in% "endemic",]
dim(dados)
toto <- rep(TRUE, length(table(dados$category, dados$main.criteria)))
toto[c(19)] <- FALSE
pie.end <- my.PieDonut(dados, aes(category, main.criteria),
                       ratioByGroup=FALSE, showPieName = FALSE,
                       start=75,r0=0, r1=1,r2=1.3,
                       showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
                       #labelposition = 2, selected=c(1:15), title="labelposition=1",
                       #explode = c(1:3), explodeDonut = TRUE,
                       main.colors = cores[c(2:5,8)],
                       pieAlpha = 0.9,
                       #donut.colors = NULL,
                       correct.legend.x = c(0,0,-0.15,-0.1,-0.1),
                       correct.legend.y = c(0.2,0,-0.05,0.05,0.1),
                       tidy.legend.donut = toto,
                       tidy.legend.pie = c(TRUE,TRUE,TRUE,FALSE,TRUE),
                       plot.piedonut = TRUE,
                       return = "pie", pieLabelSize = 6, donutLabelSize = 4.5)
                       # return = "pie", pieLabelSize = 6, donutLabelSize = 4.5,
                       # explode = 1:5,
                       # explode.cor.pie = 0.5, explode.cor.donut = 0.7)
donut.end <- my.PieDonut(dados, aes(category, main.criteria),
                         ratioByGroup=FALSE, showPieName = FALSE,
                         start=75,r0=0, r1=1,r2=1.3,
                         showRatioThreshold = 0.01, labelpositionThreshold =  0.03,
                         #labelposition = 2, selected=c(1:15), title="labelposition=1",
                         #explode = c(1:3), explodeDonut = TRUE,
                         main.colors = cores[c(2:5,8)],
                         sub.colors = crit_cores[-c(1,4)],
                         #donut.colors = NULL,
                         correct.legend.x = c(0,0,-0.15,-0.1,-0.1),
                         correct.legend.y = c(0.2,0,-0.05,0.05,0.1),
                         tidy.legend.donut = toto,
                         tidy.legend.pie = c(TRUE,TRUE,TRUE,FALSE,TRUE),
                         plot.piedonut = TRUE,
                         return = "donut", pieLabelSize = 6, donutLabelSize = 4.5)
                         # return = "donut", pieLabelSize = 6, donutLabelSize = 4.5,
                         # explode = 1:5,
                         # explode.cor.pie = 0.5, explode.cor.donut = 0.7)

## Constructing the plots
pie.donut.all <- ggdraw(pie.all) + draw_plot(donut.all)
pie.donut.end <- ggdraw(pie.end) + draw_plot(donut.end)


## Combining the plots into one image and saving
# fig1 <- plot_grid(pie.donut.all, pie.donut.end,
#                   labels = c('A - All populations', 'B - Endemics'),
#                   align="v", vjust = 7, label_size = 16
# )
fig1 <- ggdraw() + 
  draw_plot(pie.donut.all + guides(color=F), x=-0.05, y=0, width = 0.65, height = 1.07) + 
  draw_plot(pie.donut.end + guides(color=F), x=0.45, y=0, width = 0.65, height = 1.07) + 
  draw_label(label = 'A - All populations',x =0.095, y =0.95, size = 22, fontface = "bold") +
  draw_label(label = 'B - Endemics',x =0.57, y =0.95, size = 22, fontface = "bold")
fig1

path <- "./figures"
ggsave2("Figure1_AB_new.jpg", fig1, "jpeg", path,
        width = 45, height = 22, units = "cm", dpi = 320)
ggsave2("Figure1_AB_new.pdf", fig1, "pdf", path,
        width = 45, height = 22, units = "cm", dpi = 320)

## Creating fake legends
df <- data.frame(x = rep(c(1,2,3,4,5,6), each =7), y = c(letters[1:7]))
p <- ggplot(data = df, aes(y = y, fill = as.factor(y))) + 
  geom_bar() +
  theme_minimal() +
  scale_fill_manual(labels = c("A2", "A2, B1+2", "A2+D", "All", "B1", "B1+2", "B2"),
                    values = crit_cores[-1]) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 28)) + 
  # labs(fill = NULL) + 
  labs(fill = "Criteria") + 
  guides(fill = guide_legend(nrow=1))
#guides(colour = "colorbar", size = "legend", shape = "legend")
lp = get_legend(p)

cores1 <- unname(cores[-1][c(1,2,3,4,7,5)])
df1 <- data.frame(x = rep(c(1,2,3,4,5,6), each =6), y = c(letters[1:6]))
p1 <- ggplot(data = df1, aes(y = y, fill = as.factor(y))) + 
  geom_bar() +
  theme_minimal() +
  scale_fill_manual(labels = c("CR", "EN", "VU", "NT", "LC", "NA"),
                    values = cores1) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 28)) + 
  labs(fill = NULL) + 
  labs(fill = "Categories") + 
  guides(fill = guide_legend(nrow=1))
lp1 = get_legend(p1)


## Combining the plots and legends into one image and saving
fig2 <- ggdraw() + 
  draw_plot(pie.donut.all + guides(color=F), x=-0.05, y=0, width = 0.65, height = 1.07) + 
  draw_plot(pie.donut.end + guides(color=F), x=0.45, y=0, width = 0.65, height = 1.07) + 
  draw_label(label = 'A - All populations',x =0.095, y =0.95, size = 22, fontface = "bold") +
  draw_label(label = 'B - Endemics',x =0.57, y =0.95, size = 22, fontface = "bold") + 
  draw_grob(lp, x=0.025, y=-0.455) +
  draw_grob(lp1, x=0.02, y=-0.4)

ggsave2("Figure1_AB_new1.jpg", fig2, "jpg", path,
        width = 45, height = 22, units = "cm", dpi = 320)


# fig1A <- plot_grid(pie.donut.all,
#                    labels = c('A - All populations'),
#                    align="v", vjust = 7, label_size = 18
# )
# path <- "./figures"
# ggsave2("Figure1_A.pdf", fig1A, "pdf", path,
#         width = 25, height = 22, units = "cm", dpi = 320)
# 
# fig1B <- plot_grid(pie.donut.end,
#                      labels = c('B - Endemics'),
#                    align="v", vjust = 7, label_size = 18
# )
# path <- "./figures"
# ggsave2("Figure1_B.pdf", fig1B, "pdf", path,
#         width = 25, height = 22, units = "cm", dpi = 320)


#### FIGURE 2 ####
# Previous vs. new assessments

# jpeg(filename = "figures/Figure2.jpg", width = 3750, height = 2000, units = "px", pointsize = 12,
#      res = 300, family = "sans", type="cairo", bg="white")
pdf(file = "figures/Figure2.pdf", width = 13, height = 7, pointsize = 12,
    family = "sans", bg="white")
par(mfrow=c(1,2))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)

## NEW VS. PREVIOUS (GLOBAL)
mat <- as.matrix(table(paste0(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],"_new"), 
                       paste0(all.crit$redlistCategory[all.crit$endemic %in% "endemic"],"_prev")))
mat <- mat[c(1,3,6,5,4,2), rev(c(4,1,3,8,7,5,2))]

#how many LC actually remained as LC?
100*mat[5,2]/sum(mat[,2]) # 14.94%; 2n revision 17.6%

#Defining the colors of tracks and links
# grid.col = c(EX_prev = "black", CR_prev = "red", EN_prev = "darkorange", VU_prev = "gold", NT_prev = "yellowgreen", LC_prev = "forestgreen", DD_prev = "grey",
#              CR_new = "red", EN_new = "darkorange", VU_new = "gold", NT_new = "yellowgreen", LC_new = "forestgreen", DD_new = "grey")
# col_mat = rep(rev(c("black", "red", "darkorange", "gold", "yellowgreen", "forestgreen", "grey")), each=6)
grid.col = c(EX_prev = "#000000", CR_prev = "#CD0000", EN_prev = "#E69F00", VU_prev = "#F0E442", NT_prev = "#C1FF73", LC_prev = "#10E210", DD_prev = "#A0A0A0",
             CR_new = "#CD0000", EN_new = "#E69F00", VU_new = "#F0E442", NT_new = "#C1FF73", LC_new = "#10E210", DD_new = "#A0A0A0")
col_mat = rep(rev(c("#000000", "#CD0000", "#E69F00", "#F0E442", "#C1FF73", "#10E210", "#A0A0A0")), each=6)

#col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
#col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 10 & mat >= 5] = mat[mat < 10 & mat >= 5] + 1
mat[mat > 0 & mat < 5] = mat[mat > 0 & mat < 5] + 2
transp <- rep(0.3, length(col_mat))
# transp[col_mat %in% "forestgreen"] <- 0.6
transp[col_mat %in% "#10E210"] <- 0.6

#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
lava::revdiag(visible[,-c(7)]) = FALSE
chordDiagram(mat, big.gap = 10, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
             grid.col = grid.col, col = col_mat,
             transparency = 0.3 ,
             # transparency = transp ,
             #self.link = 1, 
             # link.visible = visible,
             #h=1.1, #w=0.5,
             #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
             link.lwd = 4,
             h.ratio = 0.9,
             #reduce_to_mid_line = TRUE,
             w2=0.5, rou=0.2
             #point1 = rep(0,16)
)
#Putting legends on
sec.ind <- c("CR","EN","VU","NT","LC","","DD","LC","NT","VU","EN","CR")#,"EW","EX")
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
legend("topleft","Previous assess.", bty="n", cex=1.3, adj=c(.25,2.5))
legend("topright","New assess.", bty="n", cex=1.3, adj=c(0.25,2.5))
legend("topleft",legend=expression(bold("A - Global")),
       bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.3)
# text(-0.16,1.02,"EW", cex=0.9)
text(-0.07,1.04,"EX", cex=1.1)
text(0.1,-1.025,"DD", cex=1.1)


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
mat <- mat[c(1,3,6,5,4,2), rev(c(1,3,7,6,4,2))]

#how many LC actually remained as LC?
100*mat[5,2]/sum(mat[,2]) # 14.94%; 2n revision 17.5%

#Defining the colors of tracks and links
grid.col = c(CR_prev = "#CD0000", EN_prev = "#E69F00", VU_prev = "#F0E442", NT_prev = "#C1FF73", LC_prev = "#10E210", DD_prev = "#A0A0A0",
             CR_new = "#CD0000", EN_new = "#E69F00", VU_new = "#F0E442", NT_new = "#C1FF73", LC_new = "#10E210", DD_new = "#A0A0A0")
col_mat = rep(rev(c("#CD0000", "#E69F00", "#F0E442", "#C1FF73", "#10E210", "#A0A0A0")), each=6)

# mat[mat < 15] = mat[mat < 15]*1.25
# mat[mat > 0 & mat < 5] = 5
mat[mat < 10 & mat >= 5] = mat[mat < 10 & mat >= 5] + 1
mat[mat > 0 & mat < 5] = mat[mat > 0 & mat < 5] + 2
transp <- rep(0.3, length(col_mat))
transp[col_mat %in% "#10E210"] <- 0.6

#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
chordDiagram(mat, big.gap = 10, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
             grid.col = grid.col, 
             col = col_mat,
             transparency = 0.3,
             #transparency = transp ,
             link.lwd = 4,
             h.ratio = 0.9,
             w2=0.5, rou=0.2
)
#Putting legends on
sec.ind <- c("CR","EN","VU","NT","LC","","DD","LC","NT","VU","EN","CR")
for(si in get.all.sector.index()) {
  lab <- sec.ind[which(si == get.all.sector.index())]
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), labels = lab, 
              sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
              facing = "bending", niceFacing = TRUE, col = "black")
}
legend("topleft","Previous assess.", bty="n", cex=1.3, adj=c(.25,2.5))
legend("topright","New assess.", bty="n", cex=1.3, adj=c(0.25,2.5))
legend("topleft",legend=expression(bold("B - Regional")),
       bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.2)
text(0.1,-1.025,"DD", cex=1.1)
dev.off()


#### FIGURE 3 ####
#Spatial distribution of threatened species
## see code 22_maps.R


#### NOW FIGURE XX ####
#Proportion of protected areas

jpeg(filename = "figures/Figure_XX.jpg", width = 2500, height = 2250, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

cexs <- pchs <- cores1 <- as.factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"])
# levels(cores1) <- c("red", "grey", "darkorange", "forestgreen", "yellowgreen", "gold")
levels(cores1) <- c("#CD0000", "#A0A0A0", "#E69F00", "#10E210", "#C1FF73", "#F0E442")
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
     xlab = "EOO in protected areas (%)", ylab = "Occurrences in protected areas (%)",
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
       pch=c(19,15,17,1), 
       # col=c("red", "darkorange",  "gold", "forestgreen"), 
       col=c("#CD0000", "#E69F00", "#F0E442", "#10E210"), 
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
             # col = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4),
             # outcol = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4),
             col = alpha(rev(c("#CD0000", "#E69F00", "#F0E442", "#10E210")),0.4),
             outcol = alpha(rev(c("#CD0000", "#E69F00", "#F0E442", "#10E210")),0.4),
             outpch = 19,
             ylab = "", xlab = "")
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
car::Anova(model)
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")
LABELS <- LABELS[c("CR", "EN", "VU", "LC+NT"),]

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
             # col = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4),
             # outcol = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4),
             col = alpha(rev(c("#CD0000", "#E69F00", "#F0E442", "#10E210")),0.4),
             outcol = alpha(rev(c("#CD0000", "#E69F00", "#F0E442", "#10E210")),0.4),
             outpch = 19,
             ylab = "", xlab = "")
axis(1, pos, col = NA, col.ticks = NA, 
     labels = c("LC+NT","VU","EN","CR"),  cex.axis = 0.7)

model <- lm(log(protected+1) ~ cats.f,
            data = all.crit[all.crit$endemic %in% "endemic",])
car::Anova(model)
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")

#Add the legend
par(xpd=T)
text( pos, log(120),  rev(LABELS[,1])  , col=1 )
par(xpd=F)
dev.off()


#### FIGURE 4 ####
#The Red List Index for other tropical forests
## see code 22_maps.R


#### FIGURE SU: OPTIMUM VS. HIGH CONFIDENCE LEVEL ####
critB_opt.all <- readRDS("data/criterionB_all_optim_tax_confidence.rds")
critB_high.all <- readRDS("data/criterionB_all_high_tax_confidence.rds")
res <- readRDS("data/tax_conf_effect_on_RLI.rds")

jpeg(filename = "figures/Figure_SU.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(1,2))
par(mar=c(1,1,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)

## OPTIMUM VS. HIGH
mat <- as.matrix(table(paste0(critB_opt.all$category_B,"_opt"), 
                       paste0(critB_high.all$category_B,"_hi")))
mat <- mat[c(1,2,5,3), c(3,5,2,1)]
colnames(mat) <- gsub(" ", "", colnames(mat))
rownames(mat) <- gsub(" ", "", rownames(mat))

#Defining the colors of tracks and links
# grid.col = c(CR_hi = "red", EN_hi = "darkorange", VU_hi = "gold", LCorNT_hi = "khaki",
#              CR_opt = "red", EN_opt = "darkorange", VU_opt = "gold", LCorNT_opt = "khaki")
# col_mat = rep(rev(c("red","darkorange","gold","khaki")), each=4)
grid.col = c(CR_hi = "#CD0000", EN_hi = "#E69F00", VU_hi = "#F0E442", LCorNT_hi = "#F0E68C",
             CR_opt = "#CD0000", EN_opt = "#E69F00", VU_opt = "#F0E442", LCorNT_opt = "#F0E68C")
col_mat = rep(rev(c("#CD0000","#E69F00","#F0E442","#F0E68C")), each=4)


col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"
mat[mat < 15] = mat[mat < 15]*2
mat[mat > 0 & mat < 5] = 10

#plotting the diagram
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
  if(si == "VU_hi") {
    circos.text(mean(xlim), mean(ylim), labels = "VU", 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = FALSE, col = "black")
    
  } else {
    circos.text(mean(xlim), mean(ylim), labels = lab, 
                sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
                facing = "bending", niceFacing = TRUE, col = "black")
  }  
}
legend("topleft","High confidence", bty="n", cex=1.2)
legend("topright","Optim. confidence", bty="n", cex=1.2)
legend("topleft",legend=expression(bold("A")),
       bty="n",horiz=F,cex=1.5,x.intersp=-1,y.intersp=-0.3)


## Assessing the confidence level cutoff 
par(mar=c(3,3.5,1.5,0.5))
cut <- seq(0.05,0.95, by=0.05)
ids <- cut < 0.85
par(mar=c(3,3.5,0.75,0.5), mgp=c(2,0.25,0), tcl=-0.2, las=1)
plot(res[,2][ids] ~ cut[ids], ylim = c(0.85, .92),
     xlab = "Tax. confidence level", ylab = "Red List Index", 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2, pch=19)
axis(1, at=cut, cex.axis = 1, labels = cut * 100)
lines(res[,1][ids] ~ cut[ids] , lty = 2)
lines(res[,3][ids] ~ cut[ids] , lty = 2)
# abline(h = opt.rli[2])
# abline(h = high.rli[2], col = 2)
points(res[,5][ids] ~ cut[ids], col = "red", pch=19)
lines(res[,4][ids] ~ cut[ids], col = "red", lty = 2)
lines(res[,6][ids] ~ cut[ids], col = "red", lty = 2)
# abline(v=c(0.6,0.75))
legend("topleft", expression(bold(B)), bty="n", cex=1.3,
       x.intersp=-0.7, y.intersp=0.1)
dev.off()



#### Figure WW ####
jpeg(filename = "figures/Figure_WW.jpg", width = 2500, height = 2250, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")

cats.f <- factor(all.crit$cat.reg.clean[all.crit$endemic %in% "endemic"],
                 levels = c("LC","NT","VU","EN","CR"))
levels(cats.f) <- c("LC+NT","LC+NT","VU","EN","CR")
par(mar = c(3, 3, 0.5, 0.5), mgp=c(1.8,0.5,0), las = 1, tcl = -0.25)
pos <- c(0.2,0.75,1.5,2.2)
a <- boxplot(log(prop.EOO.forest+1) ~ cats.f,
             data = all.crit[all.crit$endemic %in% "endemic",], 
             notch = TRUE, varwidth = TRUE, frame = TRUE, horizontal = FALSE,
             at = pos, #xaxt = "n", 
             yaxt = "n", 
             ylim = log(c(1,110)), xlim=c(0,2.5),
             # col = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4),
             # outcol = alpha(rev(c("red", "darkorange",  "gold", "forestgreen")),0.4),
             col = alpha(rev(c("#CD0000", "#E69F00", "#F0E442", "#10E210")),0.4),
             outcol = alpha(rev(c("#CD0000", "#E69F00", "#F0E442", "#10E210")),0.4),
             outpch = 19, cex.lab = 1.2,
             ylab = "Area of Habitat (%)", xlab = "")
# axis(1, pos, col = NA, col.ticks = NA, 
#      labels = c("LC+NT","VU","EN","CR"), cex.axis = 0.7)
sequencia <- seq(0,100,20)+1
axis(2, log(sequencia), tick = TRUE,
     labels = seq(0,100,20), cex.axis = 1)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
tmp <- all.crit[all.crit$endemic %in% "endemic",]
model <- lm(log(prop.EOO.forest+1) ~ cats.f, 
            data = tmp)
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'cats.f', conf.level=0.95)
LABELS <- generate_label_df(TUKEY , "cats.f")
LABELS[c("CR", "EN", "VU", "LC+NT"),]
summary(model)
car::Anova(model)

#Add the legend
par(xpd=T)
text( pos, log(120),  rev(LABELS[,1])  , col=1 )
par(xpd=F)
dev.off()
