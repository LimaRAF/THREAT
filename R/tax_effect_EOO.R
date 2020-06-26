##########################################################################
#### EFFECT OF TAXONOMIC CONFIDENCE OVER THE EOO ESTIMATE ####
##########################################################################
require(ConR)
require(rgeos)
source("suggestions_for_ConR.r")


## Reading the Neotropics shapefile ##
neotrop <- readRDS("E:/ownCloud/W_GIS/Am_Lat_ADM_GADM_v3.6/gadm36_Neotrop_0_sp_simplified.rds")
neotrop <- gBuffer(neotrop, byid=TRUE, width=0)

#projecting and simplifying the neotropical contours
neotrop.simp <- gSimplify(neotrop,tol=0.05)
neotrop.simp <- gBuffer(neotrop.simp, byid=TRUE, width=0)
neotrop.simp.proj <- spTransform(neotrop.simp, CRS("+init=epsg:5641"))  # define projection system of our data

## Reading herbarium data
oc.data <- readRDS("threat_occ_data.rds")

#Putting data in the format demanded by the package
MyData <- cbind.data.frame(ddlat = as.double(oc.data$latitude.work1),
                           ddlon = as.double(oc.data$longitude.work1),
                           tax = as.character(oc.data$species.correct2),
                           higher.tax.rank = oc.data$family.correct1,
                           tax.check2 = oc.data$tax.check2,
                           stringsAsFactors = FALSE)
rm(oc.data)

#### CALCULATING EOO USING DIFFERENT CONFIDENCE LEVELS OF THE OCCURRENCES ####
#teste <- MyData[grepl("Rutaceae",MyData$higher.tax.rank), c(1:3,5)] 
t1 <- Sys.time()
sens <- EOO.sensitivity (MyData[,c(1:3,5)], 
                         levels.order = c("FALSE", "cannot_check", "TRUE_TBC","TRUE_OTHER", "TRUE"),
                            occ.based = TRUE,
                            exclude.area = TRUE,
                            country_map = neotrop.simp,
                            parallel = TRUE,
                            NbeCores = 5,
                            show_progress = TRUE,
                            proj_user = 5641,
                            value = "dist")
t2 <- Sys.time()
t2 - t1
sapply(sens, head)
resultado <- sens$EOO.change

#### CLASSYFYING SPECIES ACCORDING TO 3 CRITERIA ####
tmp <- findInterval(resultado$Occs.level.5, c(15,30,75)) # number of occs indentied by family specialists
tmp0 <- findInterval(resultado$Occs.level.4, c(15,30,75)) # number of occs indentied by any specialists
tmp1 <- findInterval(resultado$Occs.level.5 / resultado$Occs.level.1, c(0.5,0.75,0.9)) # proportion of all occs indentied by faily specialists
tmp2 <- as.double(as.character(cut(resultado$EOO.increase.4, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens
tmp3 <- as.double(as.character(cut(resultado$EOO.increase.2, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens
tmp4 <- as.double(as.character(cut(resultado$EOO.increase.1, c(-Inf,0.05001, 0.09999, 0.29999, Inf), labels = c(3:0)))) # proportion of change in EOO due to the inclusion of non-taxonomic vetted specimens
table(tmp, tmp1)
table(tmp, tmp3)

boxplot(jitter(tmp2) ~ tmp1, varwidth = TRUE, notch = TRUE, 
        xlab = "Tax/all", ylab = "EOO change high to medium") # increase from high to medium in respect to the prop. of occs identfied by specialists
boxplot(jitter(tmp3) ~ tmp1, varwidth = TRUE, notch = TRUE,
        xlab = "Tax/all", ylab = "EOO change high to low") # increase from high to low in respect to the prop. of occs identfied by specialists
boxplot(jitter(tmp4) ~ tmp1, varwidth = TRUE, notch = TRUE,
        xlab = "Tax/all", ylab = "EOO change high to bad") # increase from high to bad in respect to the prop. of occs identfied by specialists

#an inex of data quality
index <- apply(cbind(tmp,tmp1,tmp4), 1, sum, na.rm = TRUE)
index1 <- round(apply(cbind(tmp,tmp1,apply(cbind(tmp2,tmp3,tmp4),1,mean,na.rm=TRUE)), 1, sum, na.rm = TRUE),1)
plot(jitter(index1) ~ index)
table(index)

## Which species we can include all TRUE_OTHER automatically:
resultado$add.other <- FALSE
resultado$add.other[!is.na(tmp2) & (tmp2 >= 2) ] <- TRUE
resultado$add.other[!is.na(tmp2) & (tmp1 >= 2) ] <- TRUE
resultado$add.other[is.na(tmp2)] <- TRUE
tbc <- unique(MyData$tax[MyData$tax.check2 %in% "TRUE_TBC"])
resultado$add.other[resultado$tax %in% tbc & (tmp3 >= 1) ] <- TRUE

## Which species we can include all TRUE_TBC automatically:
resultado$add.tbc <- FALSE
resultado$add.tbc[resultado$tax %in% tbc & (tmp<3 | is.na(tmp3)) & (tmp1 >=1 | tmp3 >=1)] <- TRUE

## Assigning classes of taxonomic confidence
resultado$tax.conf <- NA
resultado$action <- NA
#class 1: best taxonomy, many occurrences (> 75), little or no changes in EOO
resultado$tax.conf[tmp == 3] <- "class1"
resultado$tax.conf[is.na(resultado$tax.conf) & tmp0 %in% 3 & (tmp1 %in% 3 | tmp2 %in% 3)] <- "class1.1" # TRUE + TRUE_OTHER >=75 + best tax or small EOO change
resultado$action[resultado$tax.conf == "class1"] <- "nothing"
#class 2: good taxonomy
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & (tmp1>=3 | tmp4>=3)] <- "class2"   # good occs + best tax or small EOO change
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & (tmp1>=2 & tmp4>=2)] <- "class2.1" # good occs + good tax + low EOO change
#class 3
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2 & (tmp1>=2 | tmp4>=2)] <- "class3"   # good occs + good tax or low EOO change
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & (tmp1>=2 & tmp4>=2)] <- "class3.1" # low occs + good tax + low EOO change
#class 4 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & (tmp1>=2 | tmp4>=2)] <- "class4.1"   # low occs + good tax or low EOO change
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1 & (tmp1>=1 & tmp4>=1)] <- "class4.2" # low occs + low tax + considerable EOO change 
#class 5 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 2] <- "class5" # all other good
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=3 | tmp4>=3)] <- "class5.1"   # bad occs + best tax or small EOO change 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=2 & tmp4>=2)] <- "class5.2" # bad occs + good tax + low EOO change 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=2 | tmp4>=2)] <- "class5.3"   # bad occs + good tax or small EOO change 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=1 & tmp4>=1)] <- "class5.4"  # bad occs + low tax + considerable EOO change 
#class 6
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 1] <- "class6"   # all other low occs 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp >= 0 & (tmp1>=1 | tmp4>=1)] <- "class6.1"   # bad occs + low tax or considerable EOO change 
#class 7
resultado$tax.conf[is.na(resultado$tax.conf) & tmp == 0 & resultado$true >5] <- "class7"   # bad occs + low tax or considerable EOO change 
resultado$tax.conf[is.na(resultado$tax.conf) & tmp == 0] <- "class7.1" # all other bad occs

##Inspecting
table(resultado$tax.conf, useNA = "always")
table(gsub("\\.1|\\.2|\\.3|\\.4","",resultado$tax.conf), useNA = "always")
boxplot(index1 ~ resultado$tax.conf, varwidth = TRUE)

par(mfrow=c(2,2), mar=c(4,4,1,1))
boxplot(index1 ~ gsub("\\.1|\\.2|\\.3|\\.4", "", resultado$tax.conf),
        varwidth = TRUE, notch = TRUE)
boxplot(resultado$Occs.level.5 / resultado$Occs.level.1 ~ gsub("\\.1|\\.2|\\.3|\\.4", "", resultado$tax.conf),
  varwidth = TRUE, notch = TRUE)
boxplot(jitter(tmp0) ~ gsub("\\.1|\\.2|\\.3|\\.4", "", resultado$tax.conf),
          varwidth = TRUE, notch = TRUE)
boxplot(log(resultado$EOO.increase.1+1) ~ gsub("\\.1|\\.2|\\.3|\\.4", "", resultado$tax.conf),
        varwidth = TRUE, notch = TRUE)

#class x: few occurrences, unknow taxonomy
resultado$tax.conf[resultado$true < 3] <- "classX"















#### CALCULATING EOO USING DIFFERENT CONFIDENCE LEVELS ####
resultado <- MyData[!duplicated(MyData$tax), c("higher.tax.rank", "tax")]
resultado$non.dup.occurs <- unique.coords(MyData)$non.dup.occurs
tmp <- as.data.frame.matrix(table(MyData$tax, MyData$tax.check2))
tmp$false.true = apply(tmp[,c("FALSE","TRUE","TRUE_OTHER","TRUE_TBC")], 1, sum)
tmp$true.true = apply(tmp[,c("TRUE","TRUE_OTHER")], 1, sum)
tmp$true = tmp[,c("TRUE")]
resultado <- cbind.data.frame(resultado, tmp[,6:8])

## Convex Hull method
EOO.best.tax <- EOO.computing(MyData[MyData$tax.check2 %in% "TRUE",], export_shp = FALSE,
                              exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                              write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                              write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                              parallel = FALSE, NbeCores = 6) # run on parallel? How many cores?
EOO.good.tax <- EOO.computing(MyData[MyData$tax.check2 %in% c("TRUE","TRUE_OTHER"),], export_shp = FALSE,
                              exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                              write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                              write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                              parallel = FALSE, NbeCores = 6) # run on parallel? How many cores?
EOO.low.tax <- EOO.computing(MyData[MyData$tax.check2 %in% c("FALSE","TRUE","TRUE_OTHER","TRUE_TBC"),], export_shp = FALSE,
                             exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                             write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                             write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                             parallel = FALSE, NbeCores = 6) # run on parallel? How many cores?
EOO.bad.tax <- EOO.computing(MyData, export_shp = FALSE,
                             exclude.area = TRUE, country_map = neotrop.simp, # If 'exclude.area' is TRUE, the EEO is cropped using the shapefile defined in the argument country_map
                             write_shp = FALSE, # If TRUE, a directory named 'shapesIUCN' is created to receive all the shapefiles created for computing EOO
                             write_results=FALSE, file.name = "EOO.hull", # If TRUE, a csv fiel is created in the working directory
                             parallel = FALSE, NbeCores = 6) # run on parallel? How many cores?

## Merging the results
tmp = merge(resultado, EOO.bad.tax, by.x = "tax", by.y = "row.names", all.x = TRUE)
tmp1 = merge(resultado, EOO.low.tax, by.x = "tax", by.y = "row.names", all.x = TRUE)
tmp2 = merge(resultado, EOO.good.tax, by.x = "tax", by.y = "row.names", all.x = TRUE)
tmp3 = merge(resultado, EOO.best.tax, by.x = "tax", by.y = "row.names", all.x = TRUE)
resultado = cbind.data.frame(resultado, tmp$EOO, tmp1$EOO, tmp2$EOO, tmp3$EOO)
names(resultado)[7:10] <- c("bad.tax", "low.tax", "good.tax","best.tax")

## Calculating the proportional increases
resultado$best.good <- resultado$best.low <- resultado$best.bad <-   NA
#resultado$occ.best.good <- resultado$occ.best.low <-  resultado$occ.best.bad <- NA
for (i in 1:dim(resultado)[1]) {
  x <- as.numeric(resultado[i,7:10])
  #x1 <- as.numeric(resultado[i,3:6])
  if(is.na(x[4])) { 
    resultado[i,11:13] <- c(NA, NA, NA) 
  } else {
    if(x[4] == 0) { 
      resultado[i,11:13] <- c(NA, NA, NA)
    } else {  
      diffs <- 100*((x[1:3] - x[4])/x[4])
      resultado[i,11:13] <- diffs
      #diffs1 <- 100*((x1[1:3] - x1[4])/x1[4])
      #resultado[i,14:16] <- diffs1
    }
  }
}
apply(resultado[,11:13], 2, summary)
## EOO can decrease when adding more occurrences, when the new occurrences adjust the convex hull sides towards the center of the species distribution
#apply(resultado[,14:16], 2, summary)
## Differences in the number of occurrences can never be negative...


