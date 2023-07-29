#######################################
#######################################
#### PREPARING THE NEOTROPICAL MAP ####
#######################################
#######################################
rm(list=ls())

#### LOADING PACKAGES ####
library(rgeos)

#### SIMPLIFIED NEOTROPICAL COUNTOUR ####
#Loading the files
neotrop <- readRDS("data/data-raw/gadm36_Neotrop_0_sp.rds")
neotrop <- gBuffer(neotrop, byid=TRUE, width=0)

#projecting and simplifying the neotropical contours
neotrop.simp <- gSimplify(neotrop, tol=0.0001, topologyPreserve = TRUE)
neotrop.simp <- gBuffer(neotrop.simp, byid=TRUE, width=0)
neotrop.simp1 <- gSimplify(neotrop, tol=0.001, topologyPreserve = TRUE)
neotrop.simp1 <- gBuffer(neotrop.simp1, byid=TRUE, width=0)
neotrop.simp2 <- gSimplify(neotrop, tol=0.005, topologyPreserve = TRUE)
neotrop.simp2 <- gBuffer(neotrop.simp2, byid=TRUE, width=0)
neotrop.simp3 <- gSimplify(neotrop, tol=0.01, topologyPreserve = TRUE)
neotrop.simp3 <- gBuffer(neotrop.simp3, byid=TRUE, width=0)

#how many polygons left?
length(neotrop)
length(neotrop.simp)
length(neotrop.simp1)
length(neotrop.simp2)
length(neotrop.simp3)

# any bad polygons remaining?
sum(gIsValid(neotrop.simp, byid=TRUE)==FALSE) #no!
sum(gIsValid(neotrop.simp1, byid=TRUE)==FALSE) #no!
sum(gIsValid(neotrop.simp2, byid=TRUE)==FALSE) #no!
sum(gIsValid(neotrop.simp3, byid=TRUE)==FALSE) #no!

#removing very small polygons
length(neotrop.simp2)
area <- gArea(neotrop.simp2, byid=TRUE) 
neotrop.simp4 <- neotrop.simp2[area > 0.004*0.004,] #removing polygons smaller than ~16 ha/0.16 km2
length(neotrop.simp4)


#Inspecting more closely...
plot(neotrop.simp1, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=1)
plot(neotrop.simp2, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=2)
plot(neotrop.simp3, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=3)
plot(neotrop.simp4, xlim=c(-45.5,-43), ylim=c(-23.5,-23), border=4)

#Saving
# saveRDS(neotrop.simp, file = "data//Contour_Neotrop_simplified_very_large.rds")
# saveRDS(neotrop.simp1, file = "data//Contour_Neotrop_simplified_tol_001.rds")
# saveRDS(neotrop.simp2, file = "data//Contour_Neotrop_simplified_tol_005.rds")
# saveRDS(neotrop.simp3, file = "data//Contour_Neotrop_simplified_tol_01.rds")
saveRDS(neotrop.simp4, file = "data//Contour_Neotrop_simplified_tol_005_no_small.rds") # We are currently using this one

