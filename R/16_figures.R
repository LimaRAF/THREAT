##################################################################################################################################################################################################################################H
##################################################################################################################################################################################################################################H
#################
#### FIGURES ####
#################
rm(list = ls())

##Loading the data
all.crit <- readRDS("data/all.criteria.rds")

##Loading packages
require(circlize)
require(webr)
require(ggplot2)
library(cowplot)
require(multcompView)
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
# pie.df$main.criteria[pie.df$main.criteria %in% c("B1")] <- "other"
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
toto[c(11,23,36)] <- FALSE
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
                         showRatioThreshold = 0.01, labelpositionThreshold =  0.028,
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
toto[c(19)] <- FALSE
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
                  align="v", vjust = 7, label_size = 16
)

path <- "C:/Users/renato/Documents/raflima/Pos Doc/Manuscritos/Artigo Extincao na MA/data analysis/figures"
ggsave2("Figure1.jpg", fig1, "jpeg", path,
        width = 50, height = 22, units = "cm", dpi = 320)


#### FIGURE 2 ####
# Previous vs. new assessments

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
text(-0.16,1.02,"EW", cex=0.9)
text(-0.06,1.04,"EX", cex=0.9)

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
