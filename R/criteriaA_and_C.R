########################################################
#### ASSESSING POPULATION DECLINE - CRITERIA A AND C ###
########################################################
rm(list=ls())

#### LOADING PACKAGES ###
#devtools::install_github("gdauby/ConR", ref = "master", force = TRUE)
devtools::install_github("gdauby/ConR@devel")
require(ConR)
library("ConR")

#### LOADING THREAT POPULATION SIZE DATA (TREECO) ###
mean.pop.sizes <- readRDS("data/threat_mean_pop_sizes.rds")

#Putting data in the ConR format
sapply(1:length(res.high), function(x) {
  
  
})

PopData <- oc.data[, c("ddlat","ddlon",
                      "tax","higher.tax.rank",
                      "coly","vouchers",
                      "detBy","dety",
                      "tax.check.final","UC",
                      "dist.eoo","tax.conf","source")]
rm(oc.data)


