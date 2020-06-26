#### FITTING MODELS ####

## Loading packages ##
#install_version("BayesLogit", version = "0.5.1", repos = "http://cran.us.r-project.org")
library(Hmsc)
require(abind)
setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/InPreparation/Renato trees")

localDir = "."
dataDir = file.path(localDir, "data")
ModelDir = file.path(localDir, "models")

## Reading the data ##
vegtab = readRDS(file.path(dataDir, "data.rds"))$vegtab.ha
sites  = readRDS(file.path(dataDir, "data.rds"))$sites
traits = readRDS(file.path(dataDir, "data.rds"))$traits

## Filtering the data (most prevalent species) ##
Y = as.matrix(vegtab)
prev = colSums(Y>0)
Y = Y[,prev>=200]
XData = data.frame(temp = sites$temp, ppt = sites$ppt, ecoreg = as.factor(sites$ecoreg))
XData$ecoreg[XData$ecoreg=="SerraDoMar|BA_Interior"] = "SerraDoMar"
XData$ecoreg = as.factor(as.character(XData$ecoreg))
xy = as.matrix(cbind(sites$long1, sites$lat1))

studyDesign = data.frame(site = as.factor(sites$SiteCode))
rownames(xy) = studyDesign[,1]
rL = HmscRandomLevel(sData = xy)
XFormula = ~ temp + ppt + ecoreg
modeltype = 1

m = Hmsc(Y = if(modeltype==1){1*(Y>0)} else {Yabu}, XData = XData, XFormula = XFormula,
             distr = if(modeltype==1){"probit"} else {"normal"}, studyDesign = studyDesign,
             ranLevels = list(site = rL))

# Setting the default parameters
nChains = 1
thin = 1
samples = 10
transient = 5*thin
verbose = 5*thin

# Fit models to data

m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains,
                           verbose = verbose)

save(m,file=file.path(ModelDir, paste("model_",c("pa","abu")[modeltype],"_chains=", as.character(nChains),"_thin=", as.character(thin), "_","samples=",as.character(samples),".Rdata",sep = "")))

