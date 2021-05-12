################################################################
################################################################
#### PREPARING POPULATION SIZES FOR THREAT IUCN ASSESSMENTS ####
################################################################
################################################################
rm(list=ls())

## List of species with occurrence and inventory data (the final Atlantic Forest checklist)
herb_data <- readRDS("data/herbarium_spp_data.rds")
inv_data <- readRDS("data/inventory_spp_data.rds")
resultado <- dplyr::left_join(herb_data, inv_data)

#########################H
#### POPULATION SIZES ####
#########################H

# Reading the population estimates
pop.sizes <- 
  readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//pop.size.est_nmx50_mxd5_idp2.rds")
pop.sizes.2018 <- pop.sizes$`2018`$mean 
tmp <- apply(pop.sizes.2018, 2, sum)
tmp <- data.frame(species.correct2 = gsub("_"," ",names(tmp)), 
                  pop.size.2018 = as.double(tmp), 
                  stringsAsFactors = FALSE)
tmp$species.correct2 <- gsub("\\.", "-", tmp$species.correct2)

# Grouping data from the new synonyms
syn.br <- read.csv("data/new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$ï..status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$ï..status[i]
  
  if (st.i == "replace")
    tmp$species.correct2[tmp$species.correct2 %in% sp.i] <- rpl.i
}
tmp <- as.data.frame(xtabs(pop.size.2018 ~ species.correct2, tmp))
names(tmp)[2] <- "pop.size.2018"

#Joining with the original list
tmp1 <- merge(resultado, tmp, by = "species.correct2", 
              all.x = TRUE, sort = FALSE)
tmp1 <- tmp1[order(tmp1$species.correct2),]
table(resultado$species.correct2 == tmp1$species.correct2)
resultado$pop.size.2018 <- tmp1$pop.size.2018

## How many species in the checklist with pop.size estimates?
100*table(!is.na(resultado$pop.size.2018))/dim(resultado)[1] # 63%

## Pop size an number of occurrences are related?
# plot(log(resultado$total.occs+1) ~ log(resultado$pop.size.2018+1))
# abline(MASS::rlm(log(resultado$total.occs+1) ~ log(resultado$pop.size.2018+1)), lwd=2, col=3)

## Which specie swe have population sizes which are not in the checklist?
miss.sp <- tmp$species.correct2[!tmp$species.correct2 %in% resultado$species.correct2]
length(miss.sp) # 213 species; most are non-AF species
# tmp <- flora::get.taxa(miss.sp, life.form = TRUE, vegetation.type = TRUE, establishment = TRUE)
# tmp1 <- flora::get_domains(tmp)
# tmp2 <- tmp1[grepl("ata Atl", tmp1$domain),]
# spp = read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo AF checklist//data analysis//DomainsKnownTreesNeotropics.csv", as.is=T, na.string=c(NA,""," "))
# tmp3 <- merge(tmp2, spp, by.x = "original.search", by.y = "SpeciesKTSA")
# write.csv(tmp3, "tmp.csv")

## Saving the data vailable per source of information
saveRDS(resultado, "data/assess_iucn_spp.rds")


##########################################################################################################################################################################################################H
##########################################################################################################################################################################################################H
#### PREPARING THE DATA FOR THE ASSESSMENTS USING ConR ####
detach("package:ConR", unload=TRUE)
install.packages("C:/Users/renato/Documents/raflima/R_packages/ConR", # working version on Renato's local 
                 repos = NULL, 
                 type = "source")
library("ConR")

## LOADING AND EDITING THE POP DATA
# Obtaining the mean, low and high pop size estimates per year
mean.pop.sizes <- sapply(pop.sizes, function(x) apply(x$mean, 2, sum, na.rm = TRUE))
low.pop.sizes <- sapply(pop.sizes, function(x) apply(x$low, 2, sum, na.rm = TRUE))
high.pop.sizes <- sapply(pop.sizes, function(x) apply(x$high, 2, sum, na.rm = TRUE))
table(rownames(mean.pop.sizes) == rownames(low.pop.sizes))
table(rownames(mean.pop.sizes) == rownames(high.pop.sizes))

# Grouping data from the new synonyms
spp_names <- rownames(mean.pop.sizes)
spp_names <- gsub("_", " ", spp_names)
spp_names <- gsub("\\.", "-", spp_names)

for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$ï..status[i]
  
  if (st.i == "replace")
    spp_names[spp_names %in% sp.i] <- rpl.i
}
mean.pop.sizes <- rowsum(mean.pop.sizes, spp_names)
low.pop.sizes <- rowsum(low.pop.sizes, spp_names)
high.pop.sizes <- rowsum(high.pop.sizes, spp_names)

## OBTAINING POPULATION SIZES FOR MISSING YEARS ##
## Getting pop. sizes for all possible generation lengths
gen.lengths <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100)
pos.gen.length <- sort(unique(gen.lengths * rep(c(1:3), each=length(gen.lengths))), decreasing = TRUE)
miss.years <- 2018 - pos.gen.length
miss.years1 <- miss.years[miss.years < 1992]
miss.years2 <- miss.years[miss.years >= 1992]

## Looping the pop. decline trends for each species and time interval
require(snow)
require(doSNOW)
require(foreach)

dados <- list(mean.pop.sizes, low.pop.sizes, high.pop.sizes)
results <- vector("list", length(dados))
for(i in 1:length(dados)) {
  #defining the object for the loop
  pop_data <- dados[[i]]
  
  #Setting the loop parameters
  cl <- snow::makeSOCKcluster(6)
  doSNOW::registerDoSNOW(cl)
  `%d%` <- foreach::`%dopar%`
  x <- NULL
  output <-
    foreach::foreach(
      x = 1:dim(pop_data)[1],
      #.combine = 'c',
      .options.snow = NULL
    ) %d% {
      
      res1 <- ConR:::pop.decline.fit(pop.size = pop_data[x,],
                                     years = c(1850,1940,1945,1950,1955,1960,1965,1970,1975,1980,1985,1992),
                                     #models = "all", 
                                     models = c("linear", "quadratic", "exponential", 
                                                "logistic", "general_logistic"),
                                     project.years = miss.years1,
                                     plot.fit = FALSE)
      res2 <- ConR:::pop.decline.fit(pop.size = pop_data[x,],
                                     years = c(1992,1995,2000,2005,2010,2015,2018),
                                     #models = "all", 
                                     models = c("linear", "quadratic", "exponential", 
                                                "logistic", "general_logistic"),
                                     project.years = miss.years2,
                                     #parsimony = TRUE,
                                     max.count = 20,
                                     plot.fit = FALSE)
      res1$data$Modelo <- attributes(res1$best.model)$best.model.name
      res2$data$Modelo <- head(attributes(res2$best.model)$best.model.name, 1)
      dados <- rbind.data.frame(head(res1$data, -1), res2$data, stringsAsFactors = FALSE) 
      #dados$Predicted[dados$Year < 1850] <- 
      #dados$Predicted[dados$Year == 1850]
      dados <- dados[,c("Year", "Observed", "Predicted", "Modelo")]
      dados
    }
  snow::stopCluster(cl)
  
  #Saving the predictions for each table (mean, low and high) 
  results[[i]] <- output 
}

#Setting species names
# rownames(mean.pop.sizes) <- gsub("_"," ", rownames(mean.pop.sizes))
# rownames(mean.pop.sizes) <- gsub("\\.","-", rownames(mean.pop.sizes))
for(i in 1:length(results)) names(results[[i]]) <- rownames(mean.pop.sizes)

## Putting data in the ConR format
years <- results[[1]][[1]]$Year
ncols <- length(years)
spp <- names(results[[1]])
nrows <- length(spp)
mean.pop.conR <- low.pop.conR <- high.pop.conR <- matrix(NA, ncol = ncols, nrow = nrows,
                                                         dimnames = list(spp, years))
for(x in 1:length(results[[1]])) {
  mean.pop.conR[x,] <- results[[1]][[x]]$Predicted
  low.pop.conR[x,] <- results[[2]][[x]]$Predicted
  high.pop.conR[x,] <- results[[3]][[x]]$Predicted
}

#Converting into data.frames
mean.pop.conR <- cbind.data.frame(species = rownames(mean.pop.conR), mean.pop.conR, 
                                  row.names = NULL, stringsAsFactors = FALSE)
low.pop.conR <- cbind.data.frame(species = rownames(low.pop.conR), low.pop.conR, 
                                 row.names = NULL, stringsAsFactors = FALSE)
high.pop.conR <- cbind.data.frame(species = rownames(high.pop.conR), high.pop.conR, 
                                  row.names = NULL, stringsAsFactors = FALSE)

#### SAVING ###
#Saving the estimated populations (ONLY "OBSERVED" POP SIZES)
saveRDS(mean.pop.sizes, "data/threat_mean_pop_sizes.rds")
saveRDS(low.pop.sizes, "data/threat_low_pop_sizes.rds")
saveRDS(high.pop.sizes, "data/threat_high_pop_sizes.rds")

#Saving the estimated and infered populations (BOTH "OBSERVED" AND ESTIMATED/INTERPOLATED POP SIZES)
saveRDS(results[[1]], "data/threat_mean_pop_sizes_infer.rds")
saveRDS(results[[2]], "data/threat_low_pop_sizes_infer.rds")
saveRDS(results[[3]], "data/threat_high_pop_sizes_infer.rds")

#Saving the estimated and infered populations in the ConR format
saveRDS(mean.pop.conR, "data/threat_mean_pop_sizes_for_ConR.rds")
saveRDS(low.pop.conR, "data/threat_low_pop_sizes_for_ConR.rds")
saveRDS(high.pop.conR, "data/threat_high_pop_sizes_for_ConR.rds")
