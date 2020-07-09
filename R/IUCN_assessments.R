#################################################
#### PERFORMIN THE IUCN CRITERIA ASSESSMENTS ####
#################################################
rm(list=ls())


#### LOADING PACKAGES AND FUNCTIONS ####
#devtools::install_github("gdauby/ConR", ref = "master", force = TRUE)
require(ConR)
source("./R/suggestions_for_ConR.r")
source("C://Users//renato//Documents//raflima//R_packages//ConR//R//cat_criterion_b.R")

#### OBTAINING THE NAME AND NUMBER OF OCCURRENCES (TOTAL AND SPATIALLY-UNIQUE) PER SPECIES ####
resultado <- readRDS("data/assess_iucn_spp.rds")
names(resultado)[1:2] <- c("family","tax")

# Olhar com mais cuidado a lista de species assigned para mim pelo GTA
renato.gta = read.csv("renato_assignments_GTA.csv", as.is= TRUE)

#### READING THE PARAMETERS FOR EACH CRITERION ####
critB_high <- readRDS("data/criteriaB_high_confidence.rds")
critB_high$tax <- as.character(critB_high$tax) 
critB_low <- readRDS("data/criteriaB_low_confidence.rds")
critB_low$tax <- as.character(critB_low$tax) 

##Convert Nbe_loc to Nbe_loc + Nbe_loc_PA?


#################################################################################################################################################################################H
#################################################################################################################################################################################H
##########################
#### IUCN ASSESSMENTS ####
##########################

#### CRITERION B ####
tmp <- cat_criterion_b(EOO = critB_high$EOO,
                AOO = critB_high$AOO,
                locations = as.double(apply(critB_high[,c("Nbe_loc","Nbe_loc_PA")], 1, sum, na.rm = TRUE)))
table(tmp$ranks_B12a)
table(tmp$cat_codes)
