#####################################################################################
#### SIMULATING THE EFFECT OF DIFFERENT PATTERNS OF HABITAT LOSS ON POP. DECLINE ####
#####################################################################################
rm(list = ls())

##Packages
require(sads)
require(spatstat)
require(mobsim)
require(Rcpp)
require(sp)
require(sf)
require(SpatialTools)

####################################################################
#### SIMULATING DIFFERENT AGGREGATION X DEFORESTATION SCENARIOS ####
####################################################################
# Setting the parameters for generating the random Log-normal SAD
iterations <- 499
mulog <- 5
sdlog <- 2
S <- 300
n_sim <- 15000
xlim <- 200
ylim <- 500
raio <- 20
area <- (xlim * ylim)
IUCN.threshold1 <- 0.3
IUCN.threshold2 <- 0.5
IUCN.threshold3 <- 0.8


# Generating the grid
sfc <- st_sfc(st_polygon(list(rbind(c(0,0), c(xlim,0), c(xlim,ylim), c(0,0)))))
rec_grid <- st_make_grid(sfc, cellsize = 10, square = TRUE)
rec_grid <- as_Spatial(rec_grid)
# gridded(rec_grid) <- TRUE
rec_grid.df <- as(rec_grid, "SpatialPolygonsDataFrame")
rec_grid.df@data$ID <- 
        as.character(sapply(slot(rec_grid, "polygons"), function(x) slot(x, "ID")))
rec_grid.df@data <- rec_grid.df@data[,-1, drop = FALSE]

# Starting the simulations
simula.result <- vector("list", iterations)
for (j in seq_len(iterations)){
        # Generating the random SAD
        rand.sad <- sim_sad(s_pool = S, n_sim = n_sim, sad_type = "lnorm",
                            sad_coef = list("meanlog" = mulog, "sdlog" = sdlog))
        
        # Species randomnly distributed
        rand.comm <- sim_poisson_coords(abund_vec = rand.sad, 
                                        xrange = c(0, xlim), yrange = c(0, ylim))$census
        
        # Species with conspecific aggregation
        clump.comm <- NULL
        for (i in seq_along(rand.sad)) {
                ni <- runif(1, 1, 6)
                x <- rThomas(kappa = rand.sad[i]/(ni*area), scale = 1.5, mu = ni, 
                             win = owin(c(0, xlim), c(0, ylim)))
                if (x$n > 0) {
                        sp.i <- cbind.data.frame(x = x$x, y = x$y, species = names(rand.sad[i]))
                        if(is.null(clump.comm)) clump.comm <- sp.i else clump.comm <- rbind.data.frame(clump.comm, sp.i)
                }
        }
        
        # Species with conspecific aggregation, and range restriction
        clump.comm1 <- NULL
        for (i in seq_along(rand.sad)) {
                ni <- runif(1, 1, 6)
                
                seqs <- 0:xlim
                probs <- dbeta(seq(0,1,length.out = length(seqs)), 0.5, 0.5)
                probs[is.infinite(probs)] <- 6
                x.range <- sort(c(sample(seqs, 1, prob = sqrt(probs)), 
                                  sample(seqs, 1, prob = sqrt(probs))))
                while(diff(x.range)<25)
                        x.range <- sort(c(sample(seqs, 1, prob = sqrt(probs)), 
                                          sample(seqs, 1, prob = sqrt(probs))))
                
                seqs <- 0:ylim
                probs <- dbeta(seq(0,1,length.out = length(seqs)), 0.5, 0.5)
                probs[is.infinite(probs)] <- 6
                y.range <- sort(c(sample(seqs, 1, prob = sqrt(probs)), 
                                  sample(seqs, 1, prob = sqrt(probs))))
                while(diff(y.range)<25)
                        y.range <- sort(c(sample(seqs, 1, prob = sqrt(probs)), 
                                          sample(seqs, 1, prob = sqrt(probs))))
                
                area.i <- (diff(x.range) * diff(y.range))
                
                x <- rThomas(kappa = rand.sad[i]/(ni*area.i), scale = 1.5, mu = ni, 
                             win = owin(x.range, y.range))
                if (x$n > 0) {
                        sp.i <- cbind.data.frame(x = x$x, y = x$y, species = names(rand.sad[i]))
                        if(is.null(clump.comm1)) clump.comm1 <- sp.i else clump.comm1 <- rbind.data.frame(clump.comm1, sp.i)
                }
        }
        # plot(clump.comm1[,1:2], col = as.factor(clump.comm1[,3]), pch = 19)
        
        # Overlaying the grid with the species distribution
        rand.comm.sp <- rand.comm 
        coordinates(rand.comm.sp) <- ~x + y
        rand.comm$cellID <- as.character(sp::over(rand.comm.sp, rec_grid.df)[,1])
        clump.comm.sp <- clump.comm 
        coordinates(clump.comm.sp) <- ~x + y
        clump.comm$cellID <- as.character(sp::over(clump.comm.sp, rec_grid.df)[,1])
        clump.comm.sp1 <- clump.comm1 
        coordinates(clump.comm.sp1) <- ~x + y
        clump.comm1$cellID <- as.character(sp::over(clump.comm.sp1, rec_grid.df)[,1])
        
        ## REMOVING CELLS AT RANDOM, BUT IN SEQUENCE ##
        seq.loss <- seq(25, 975, 25)
        random.loss <- vector("list", length(seq.loss))
        pool.ids <- rec_grid.df$ID
        for (i in seq_len(length(seq.loss))) {
                # random.pts <- spsample(rec_grid.df, seq.loss[i], 
                #                        type = "random", replace = FALSE)
                # random.IDS <- as.character(sp::over(random.pts, rec_grid.df)[,1])        
                # random.loss[[i]] <- random.IDS
                sample.IDs <- sample(pool.ids, 25)
                pool.ids <- pool.ids[!pool.ids %in% sample.IDs]
                if (i==1) {
                        random.loss[[i]] <- sample.IDs
                } else {
                        random.loss[[i]] <- c(random.loss[[i-1]], sample.IDs)
                }
        }
        
        # Calculating population declines
        result.loss <- result.loss1 <- result.loss2 <- 
                matrix(NA, ncol = 6, nrow = length(seq.loss),
                       dimnames = list(NULL, c("comm_type","loss", "VU", "EN", "CR", 
                                               "threat")))
        for(i in seq_len(length(seq.loss))) {
                remove_these <- random.loss[[i]] 
                
                #Random species distribution
                classes <- unique(rand.comm$species)
                sad.orig <- table(factor(rand.comm$species, levels = classes))
                comm.loss <- rand.comm[!rand.comm$cellID %in% remove_these,]
                sad.loss <- table(factor(comm.loss$species, levels = classes))
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                # threat <- (vu + en + cr)/length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/
                        ((xlim/10) * (ylim/10))
                result.loss[i,] <- cbind("random", loss, vu, en, cr, threat)
                
                #Clumped species distribution
                classes <- unique(clump.comm$species)
                sad.orig <- table(factor(clump.comm$species, levels = classes))
                comm.loss <- clump.comm[!clump.comm$cellID %in% remove_these,]
                sad.loss <- table(factor(comm.loss$species, levels = classes))
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/
                        ((xlim/10) * (ylim/10))
                result.loss1[i,] <- cbind("clumped", loss, vu, en, cr, threat)
                
                #Clumped species distribution with stratification
                classes <- unique(clump.comm1$species)
                sad.orig <- table(factor(clump.comm1$species, levels = classes))
                comm.loss <- clump.comm1[!clump.comm1$cellID %in% remove_these,]
                sad.loss <- table(factor(comm.loss$species, levels = classes))
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/
                        ((xlim/10) * (ylim/10))
                result.loss2[i,] <- cbind("clumped1", loss, vu, en, cr, threat)
        }
        # plot(result.loss[,2:3], ylim = c(0,1))
        # points(result.loss1[,2:3], col = 2)
        # points(result.loss2[,2:3], col = 3)
        random.losses <- cbind.data.frame(loss_type = "random", 
                                          rbind.data.frame(result.loss, result.loss1, result.loss2))
        
        ## REMOVING CELLS WITH AGGREGATION ##
        seq.loss <- seq(25, 3500, 25)
        coords <- sp::coordinates(rec_grid.df)
        closest <- RANN::nn2(data=coords, k=length(seq.loss))[[1]]
        row.names(closest) <- rec_grid.df$ID
        
        clump.loss <- vector("list", length(seq.loss))
        centers <- sample(row.names(closest), 25)
        for (i in seq_len(length(seq.loss))) {
                # clump.pts <- spsample(rec_grid.df, seq.loss1[i], 
                #                       type = "clustered", nclusters = 2+i)
                # clump.IDS <- as.character(sp::over(clump.pts, rec_grid.df)[,1])        
                sample.IDs <- paste0("ID", as.vector(closest[centers, seq_len(i)]))
                clump.loss[[i]] <- sample.IDs
        }
        
        # Calculating population declines
        result.loss3 <- result.loss4 <- result.loss5 <- 
                matrix(NA, ncol = 6, nrow = length(seq.loss),
                       dimnames = list(NULL, c("comm_type","loss", "VU", "EN", 
                                               "CR", "threat")))
        for(i in seq_len(length(seq.loss))) {
                remove_these <- clump.loss[[i]] 
                
                #Random species distribution
                classes <- unique(rand.comm$species)
                sad.orig <- table(factor(rand.comm$species, levels = classes))
                comm.loss <- rand.comm[!rand.comm$cellID %in% remove_these,]
                sad.loss <- table(factor(comm.loss$species, levels = classes))
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/
                        ((xlim/10) * (ylim/10))
                result.loss3[i,] <- cbind("random", loss, vu, en, cr, threat)
                
                #Clumped species distribution
                classes <- unique(clump.comm$species)
                sad.orig <- table(factor(clump.comm$species, levels = classes))
                comm.loss <- clump.comm[!clump.comm$cellID %in% remove_these,]
                sad.loss <- table(factor(comm.loss$species, levels = classes))
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/
                        ((xlim/10) * (ylim/10))
                result.loss4[i,] <- cbind("clumped", loss, threat)
                
                #Clumped species distribution with stratification
                classes <- unique(clump.comm1$species)
                sad.orig <- table(factor(clump.comm1$species, levels = classes))
                comm.loss <- clump.comm1[!clump.comm1$cellID %in% remove_these,]
                sad.loss <- table(factor(comm.loss$species, levels = classes))
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/
                        ((xlim/10) * (ylim/10))
                result.loss5[i,] <- cbind("clumped1", loss, threat)
        }
        # plot(result.loss3[,2:3], ylim = c(0,1), xlim=c(0,1))
        # points(result.loss4[,2:3], col = 2)
        # points(result.loss5[,2:3], col = 3)
        clump.losses <- cbind.data.frame(loss_type = "clumped", 
                                         rbind.data.frame(result.loss3, result.loss4, result.loss5))
        
        all.losses <- cbind(rbind.data.frame(random.losses, clump.losses), iter = j)
        simula.result[[j]] <- all.losses
}
saveRDS(simula.result, "data//simula_pop_decline.rds")

##################################
#### THE ATLANTIC FOREST CASE ####
##################################
af.grid <- readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//af_hex_grid_50km.rds")
pop.sizes <- readRDS("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//pop.size.est_nmx50_mxd5_idp2.rds")
pop.sizes <- pop.sizes$`1850`$mean
row.names(pop.sizes) <- paste0("ID", 1:dim(pop.sizes)[1])

#### ALL POPULATIONS ####
# Starting the simulations
simula.result.af <- vector("list", iterations)
for (j in seq_len(iterations)){

        cat(j,"\n")
        ## REMOVING CELLS AT RANDOM, BUT IN SEQUENCE ##
        seq.loss <- seq(30, 1140, 30)
        random.loss <- vector("list", length(seq.loss))
        pool.ids <- row.names(pop.sizes)
        for (i in seq_len(length(seq.loss))) {
                sample.IDs <- sample(pool.ids, 30)
                pool.ids <- pool.ids[!pool.ids %in% sample.IDs]
                if (i==1) {
                        random.loss[[i]] <- sample.IDs
                } else {
                        random.loss[[i]] <- c(random.loss[[i-1]], sample.IDs)
                }
        }
        
        # Calculating population declines
        result.loss <- matrix(NA, ncol = 6, nrow = length(seq.loss),
                              dimnames = list(NULL, c("loss_type","loss", 
                                                      "VU", "EN", "CR", "threat")))
        sad.orig <- apply(pop.sizes, 2, sum)
        
        for(i in seq_len(length(seq.loss))) {
                remove_these <- random.loss[[i]]
                
                comm.loss <- pop.sizes[!row.names(pop.sizes) %in% remove_these,]
                if (!is.null(dim(comm.loss))) {
                        sad.loss <- apply(comm.loss, 2, sum)        
                } else {
                        sad.loss <- comm.loss
                }
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/dim(pop.sizes)[1]
                result.loss[i, ] <- cbind("random", loss, vu, en, cr, threat)
        }
        
        ## REMOVING CELLS WITH AGGREGATION ##
        seq.loss <- seq(30, 1140*4, 30)
        coords <- as.data.frame(coordinates(af.grid))
        closest <- RANN::nn2(data=coords, k=length(seq.loss))[[1]]
        row.names(closest) <- row.names(pop.sizes)
        
        clump.loss <- vector("list", length(seq.loss))
        centers <- sample(row.names(closest), 30)
        for (i in seq_len(length(seq.loss))) {
                sample.IDs <- paste0("ID", as.vector(closest[centers, seq_len(i)]))
                clump.loss[[i]] <- sample.IDs
        }
        
        # Calculating population declines
        result.loss3 <- matrix(NA, ncol = 6, nrow = length(seq.loss),
                               dimnames = list(NULL, c("loss_type","loss", "VU", "EN", "CR", 
                                                       "threat")))
        for(i in seq_len(length(seq.loss))) {
                remove_these <- clump.loss[[i]] 
                
                comm.loss <- pop.sizes[!row.names(pop.sizes) %in% remove_these,]
                if (!is.null(dim(comm.loss))) {
                        sad.loss <- apply(comm.loss, 2, sum)        
                } else {
                        sad.loss <- comm.loss
                }
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/dim(pop.sizes)[1]
                result.loss3[i, ] <- cbind("clumped", loss, vu, en, cr, threat)
        }
        all.losses <- cbind(rbind.data.frame(result.loss, result.loss3), iter = j)
        simula.result.af[[j]] <- all.losses
}                
saveRDS(simula.result.af, "data//simula_pop_decline_AtlanticForest.rds")


#### ONLY ENDEMICS ####
all.crit <- readRDS("data/all.criteria.rds")
pop.sizes.end <- pop.sizes[, gsub("_", " ", colnames(pop.sizes)) %in% 
                                   all.crit$species[all.crit$endemic %in% "endemic"]]
af.grid.end <- af.grid[apply(pop.sizes.end, 1, sum) > 0, ]
pop.sizes.end <- pop.sizes.end[apply(pop.sizes.end, 1, sum) > 0,]

# Starting the simulations
simula.result.af.end <- vector("list", iterations)
for (j in seq_len(iterations)){
        cat(j,"\n")
        ## REMOVING CELLS AT RANDOM, BUT IN SEQUENCE ##
        seq.loss <- seq(30, 1140, 30)
        random.loss <- vector("list", length(seq.loss))
        pool.ids <- row.names(pop.sizes.end)
        for (i in seq_len(length(seq.loss))) {
                sample.IDs <- sample(pool.ids, 30)
                pool.ids <- pool.ids[!pool.ids %in% sample.IDs]
                if (i==1) {
                        random.loss[[i]] <- sample.IDs
                } else {
                        random.loss[[i]] <- c(random.loss[[i-1]], sample.IDs)
                }
        }
        
        # Calculating population declines
        result.loss <- matrix(NA, ncol = 6, nrow = length(seq.loss),
                              dimnames = list(NULL, c("loss_type","loss", 
                                                      "VU", "EN", "CR", "threat")))
        sad.orig <- apply(pop.sizes.end, 2, sum)
        
        for(i in seq_len(length(seq.loss))) {
                remove_these <- random.loss[[i]]
                
                comm.loss <- pop.sizes.end[!row.names(pop.sizes.end) %in% remove_these,]
                if (!is.null(dim(comm.loss))) {
                        sad.loss <- apply(comm.loss, 2, sum)        
                } else {
                        sad.loss <- comm.loss
                }
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/dim(pop.sizes.end)[1]
                result.loss[i, ] <- cbind("random", loss, vu, en, cr, threat)
        }
        
        ## REMOVING CELLS WITH AGGREGATION ##
        seq.loss <- seq(30, 1140*4, 30)
        coords <- as.data.frame(sp::coordinates(af.grid.end))
        closest <- RANN::nn2(data=coords, k=length(seq.loss))[[1]]
        row.names(closest) <- row.names(pop.sizes.end)
        
        clump.loss <- vector("list", length(seq.loss))
        centers <- sample(row.names(closest), 30)
        for (i in seq_len(length(seq.loss))) {
                sample.IDs <- paste0("ID", as.vector(closest[centers, seq_len(i)]))
                clump.loss[[i]] <- sample.IDs
        }
        
        # Calculating population declines
        result.loss3 <- matrix(NA, ncol = 6, nrow = length(seq.loss),
                               dimnames = list(NULL, c("loss_type","loss", "VU", "EN", "CR", 
                                                       "threat")))
        for(i in seq_len(length(seq.loss))) {
                remove_these <- clump.loss[[i]] 
                
                comm.loss <- pop.sizes.end[!row.names(pop.sizes.end) %in% remove_these,]
                if (!is.null(dim(comm.loss))) {
                        sad.loss <- apply(comm.loss, 2, sum)        
                } else {
                        sad.loss <- comm.loss
                }
                diff.loss <- (sad.orig - sad.loss)/sad.orig
                vu <- sum(diff.loss > IUCN.threshold1 & diff.loss <= IUCN.threshold2)/
                        length(sad.orig)
                en <- sum(diff.loss > IUCN.threshold2 & diff.loss <= IUCN.threshold3)/
                        length(sad.orig)
                cr <- sum((sad.orig - sad.loss)/sad.orig > IUCN.threshold3)/
                        length(sad.orig)
                threat <- sum(diff.loss > IUCN.threshold1)/length(sad.orig)
                loss <- length(unique(remove_these))/dim(pop.sizes.end)[1]
                result.loss3[i, ] <- cbind("clumped", loss, vu, en, cr, threat)
        }
        all.losses <- cbind(rbind.data.frame(result.loss, result.loss3), iter = j)
        simula.result.af.end[[j]] <- all.losses
}                
saveRDS(simula.result.af.end, "data//simula_pop_decline_AtlanticForest_endemics.rds")

###############################
#### ANALYZING THE RESULTS ####
###############################
simula.result <- readRDS("data//simula_pop_decline.rds")
simula.result.af <- readRDS("data//simula_pop_decline_AtlanticForest.rds")
simula.result.af.end <- readRDS("data//simula_pop_decline_AtlanticForest_endemics.rds")


## Simulations
all.data <- do.call(rbind, simula.result)
all.data$combo <- paste(all.data$loss_type, all.data$comm_type, sep= "_")
all.data$loss1 <- cut(as.double(all.data$loss), breaks = seq(0,1,0.025), )
levels(all.data$loss1) <- seq(0.025,1,0.025)
all.data$loss1 <- as.double(as.character(all.data$loss1))

## Atlantic Forest
all.data.af <- do.call(rbind, simula.result.af)
all.data.af$combo <- all.data.af$loss_type
all.data.af$loss1 <- cut(as.double(all.data.af$loss), breaks = seq(0,1,0.025), )
levels(all.data.af$loss1) <- seq(0.025,1,0.025)
all.data.af$loss1 <- as.double(as.character(all.data.af$loss1))

## Atlantic Forest endemics
all.data.af.end <- do.call(rbind, simula.result.af.end)
all.data.af.end$combo <- all.data.af.end$loss_type
all.data.af.end$loss1 <- cut(as.double(all.data.af.end$loss), breaks = seq(0,1,0.025), )
levels(all.data.af.end$loss1) <- seq(0.025,1,0.025)
all.data.af.end$loss1 <- as.double(as.character(all.data.af.end$loss1))
# getting the summary tables
tmp <- all.data.af.end[all.data.af.end$combo %in% "random", ]
threat.per.loss.random <- as.data.frame.list(
        aggregate(as.double(tmp$threat), list(tmp$loss1), summary))
tmp <- all.data.af.end[all.data.af.end$combo %in% "clumped", ]
threat.per.loss.clumped <- as.data.frame.list(
        aggregate(as.double(tmp$threat), list(tmp$loss1), summary))
threat.per.loss.clumped.VU <- as.data.frame.list(
        aggregate(as.double(tmp$VU), list(tmp$loss1), summary))
threat.per.loss.clumped.EN <- as.data.frame.list(
        aggregate(as.double(tmp$EN), list(tmp$loss1), summary))
threat.per.loss.clumped.CR <- as.data.frame.list(
        aggregate(as.double(tmp$CR), list(tmp$loss1), summary))


## Atlantic Forest remnants and threat level
info <- read.csv("data/hotspots_info.csv")
lu_hotspots <- read.csv("data/esa_2018_lc_per_hotspot.csv")
info <- dplyr::left_join(info, lu_hotspots[,c("hotspot.region","area_km2",
                                              "ForestCover",
                                              "Mosaic_ForestCover_GrassLand",
                                              "Mosaic_GrassLand_ForestCover",
                                              "OpenForestCover",
                                              "ShrubLand","ShrubLand_Caatinga.",
                                              "ShrubLand_Cerrado","Shrubland_Restinga")])
selected <- info$hotspot.region == "Atlantic Forest"
threat.end <- I(info$observed_percentage_threat_endemic_all_criteria[selected]/100)
loss <- I((100 - info$ForestCover[selected])/100)
losses <- I(100 - info$ForestCover)/100
trop.hotspots <- c("Atlantic Forest", "Caribbean Islands", 
                   "Coastal Forests of Eastern Africa", #"East Melanesian Islands",
                   "Eastern Afromontane", "Guinean Forests of West Africa",
                   "Indo-Burma", "Madagascar and the Indian Ocean Islands",
                   "Mesoamerica", "New Caledonia", "Philippines", 
                   "Sundaland", "Tropical Andes", "Tumbes-Choco-Magdalena", 
                   "Wallacea", "Western Ghats and Sri Lanka",
                   "Amazon", "Central Africa", "New Guinea")
trop_hotspots <- info$hotspot.region %in% trop.hotspots

## Atlantic Forest threat level only for endemics and criteria A
all.crit <- readRDS("data/all.criteria.rds")
all.crit.end <- all.crit[all.crit$endemic %in% "endemic", ]
all.crit.end.A <- all.crit.end[!is.na(all.crit.end$reduction_A12), ]
100*table(all.crit.end.A$category_A, useNA = "always")/
        dim(all.crit.end.A)[1]
threat.end.A.cats <- (table(all.crit.end.A$category_A, useNA = "always")/
                dim(all.crit.end.A)[1])[c(1,2,4)]
threat.end.A <- sum(threat.end.A.cats)

# # 2018 Forest Cover used for analyses
# esa.grid <- read.csv("C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Hyperdominance//ESA_LandUse_1992_2018_50km.csv", as.is=TRUE)
# esa.grid <- esa.grid[, grepl("ForestCover", names(esa.grid))]
# esa.grid <- esa.grid[, grepl("2018", names(esa.grid))]
# summary(100 - esa.grid)

# Getting the spline fit and predictions
spl.low <- smooth.spline(x= I(c(0, threat.per.loss.clumped$Group.1, 1)), 
              y = I(c(0, threat.per.loss.clumped$x.1st.Qu., 1)), spar = 0.001)
spl.med <- smooth.spline(x= I(c(0, threat.per.loss.clumped$Group.1, 1)), 
                         y = I(c(0, threat.per.loss.clumped$x.Median, 1)), spar = 0.001)
spl.high <- smooth.spline(x= I(c(0, threat.per.loss.clumped$Group.1, 1)), 
                         y = I(c(0, threat.per.loss.clumped$x.3rd.Qu., 1)), spar = 0.001)
splines <- list(spl.low, spl.med, spl.high)
names(splines) <- c("spline_threat.per.loss.clumped_1stQt",
                    "spline_threat.per.loss.clumped_Median",
                    "spline_threat.per.loss.clumped_3rdQt")
saveRDS(splines, "data/spline_threat.per.loss.clumped.rds")
# spl.pr.low <- predict(spl.low, x=losses[!is.na(losses)], na.action=na.exclude)
# spl.pr.med <- predict(spl.med, x=losses[!is.na(losses)], na.action=na.exclude)
# spl.pr.high <- predict(spl.high, x=losses[!is.na(losses)], na.action=na.exclude)

# x <- I(c(0, threat.per.loss.clumped$Group.1, 1))
# y <- I(c(0, threat.per.loss.clumped$x.Median, 1))
# plot(y ~ x)
# mod <- nls(y ~ SSlogis(x, Asym, xmid, scal))
# lines(predict(mod) ~ x, lwd = 2, col=2)
# f <- y ~ exp(a + b*x)/(1 + exp(a + b*x))
# mod1 <- try(stats::nls(f, 
#                         start = list(a=1,b=-0.1)), TRUE)
# lines(predict(mod1) ~ x, lwd = 2, col=2)
# lines(predict(spl.med, x = x)$y ~ x, lwd = 2, col=8)


#VU
spl.low.VU <- smooth.spline(x= I(c(0, threat.per.loss.clumped.VU$Group.1, 1)), 
                         y = I(c(0, threat.per.loss.clumped.VU$x.1st.Qu., 1)), spar = 0.001)
spl.med.VU <- smooth.spline(x= I(c(0, threat.per.loss.clumped.VU$Group.1, 1)), 
                         y = I(c(0, threat.per.loss.clumped.VU$x.Median, 1)), spar = 0.001)
spl.high.VU <- smooth.spline(x= I(c(0, threat.per.loss.clumped.VU$Group.1, 1)), 
                          y = I(c(0, threat.per.loss.clumped.VU$x.3rd.Qu., 1)), spar = 0.001)
splines.VU <- list(spl.low.VU, spl.med.VU, spl.high.VU)
names(splines.VU) <- c("spline_threat.per.loss.clumped_1stQt_VU",
                    "spline_threat.per.loss.clumped_Median_VU",
                    "spline_threat.per.loss.clumped_3rdQt_VU")
saveRDS(splines.VU, "data/spline_threat.per.loss.clumped_VU.rds")

#EN
spl.low.EN <- smooth.spline(x= I(c(0, threat.per.loss.clumped.EN$Group.1, 1)), 
                            y = I(c(0, threat.per.loss.clumped.EN$x.1st.Qu., 1)), spar = 0.001)
spl.med.EN <- smooth.spline(x= I(c(0, threat.per.loss.clumped.EN$Group.1, 1)), 
                            y = I(c(0, threat.per.loss.clumped.EN$x.Median, 1)), spar = 0.001)
spl.high.EN <- smooth.spline(x= I(c(0, threat.per.loss.clumped.EN$Group.1, 1)), 
                             y = I(c(0, threat.per.loss.clumped.EN$x.3rd.Qu., 1)), spar = 0.001)
splines.EN <- list(spl.low.EN, spl.med.EN, spl.high.EN)
names(splines.EN) <- c("spline_threat.per.loss.clumped_1stQt_EN",
                       "spline_threat.per.loss.clumped_Median_EN",
                       "spline_threat.per.loss.clumped_3rdQt_EN")
saveRDS(splines.EN, "data/spline_threat.per.loss.clumped_EN.rds")

#CR
spl.low.CR <- smooth.spline(x= I(c(0, threat.per.loss.clumped.CR$Group.1, 1)), 
                            y = I(c(0, threat.per.loss.clumped.CR$x.1st.Qu., 1)), spar = 0.001)
spl.med.CR <- smooth.spline(x= I(c(0, threat.per.loss.clumped.CR$Group.1, 1)), 
                            y = I(c(0, threat.per.loss.clumped.CR$x.Median, 1)), spar = 0.001)
spl.high.CR <- smooth.spline(x= I(c(0, threat.per.loss.clumped.CR$Group.1, 1)), 
                             y = I(c(0, threat.per.loss.clumped.CR$x.3rd.Qu., 1)), spar = 0.001)
splines.CR <- list(spl.low.CR, spl.med.CR, spl.high.CR)
names(splines.CR) <- c("spline_threat.per.loss.clumped_1stQt_CR",
                       "spline_threat.per.loss.clumped_Median_CR",
                       "spline_threat.per.loss.clumped_3rdQt_CR")
saveRDS(splines.CR, "data/spline_threat.per.loss.clumped_CR.rds")


#### FIGURE ST ####
jpeg(filename = "figures/Figure_ST.jpg", width = 4000, height = 4000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white")
par(mfrow=c(2,2))
par(mar=c(3,3,1,1), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)

# PANEL A - AF - endemics
plot(c(0,1) ~ c(1,0), col = "white", cex.lab = 1.4,
     xlab = "", ylab = "Threatened species")
tmp1 <- threat.per.loss.random
lines(c(0,tmp1$x.Median,1) ~  c(0,tmp1$Group.1,1), 
      col = "red", lwd = 2, lty = 1)
# segments(x0 = c(0,as.double(tmp1$Group.1),1)-0.001, x1 = c(0,as.double(tmp1$Group.1),1)-0.001,
#          y0 = c(0,as.double(tmp1$x.1st.Qu.),1), y1 = c(0,as.double(tmp1$x.3rd.Qu.),1), 
#          col = "red", lty = 1)
arrows(x0 = c(0,as.double(tmp1$Group.1),1)-0.001, 
       y0 = c(0,as.double(tmp1$x.1st.Qu.),1),
       y1 = c(0,as.double(tmp1$x.3rd.Qu.),1),
       col = "red", code = 3, angle = 90, length = 0.025)

tmp1 <- threat.per.loss.clumped
lines(c(0,tmp1$x.Median,1) ~  c(0,tmp1$Group.1,1), 
      col = "black", lwd = 2, lty = 1)
# segments(x0 = as.double(tmp1$Group.1)+0.001, x1 = as.double(tmp1$Group.1)+0.001,
#          y0 = as.double(tmp1$x.1st.Qu.), y1 = as.double(tmp1$x.3rd.Qu.), 
#          col = "black", lty = 1)
arrows(x0 = c(0,as.double(tmp1$Group.1),1)+0.001, 
       y0 = c(0,as.double(tmp1$x.1st.Qu.),1),
       y1 = c(0,as.double(tmp1$x.3rd.Qu.),1),
       col = "black", code = 3, angle = 90, length = 0.025)
# abline(0,1, lty = 3, col = "grey")
# points(threat.end ~ loss, pch = 19)
points(threat.end.A ~ loss, pch = 19, cex = 1.1)
points(0.114, 0.09, pch = 17, cex = 1.2) # Hans criterion A
legend("bottomright", c("Random", "Aggregate"),
       # col = seq_along(combos) + 1,
       lty = 1,
       lwd = 2, bty = "n",
       cex = 1.2, col = c("red", "black"))
legend(0.7, 0.125, expression(bold("Type of loss:")),
       #pch = c(19,17),
       bty = "n", #x.intersp= 1.85,
       cex = 1.2, col = "black")
legend(0.3, 0.125, expression(bold("Observed threat:")),
       #pch = c(19,17),
       bty = "n", #x.intersp= 1.85,
       cex = 1.2, col = "black")
legend("bottom", c("Atlantic Forest","Amazon"),
       pch = c(19,17),
       bty = "n", x.intersp= 1.5,
       cex = 1.2, col = "black")
legend("topleft",legend=expression(bold("A")),
       bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.1)

# PANEL B - Clumped loss + VU/EN/CR curves
plot(c(0,1) ~ c(1,0), col = "white", cex.lab = 1.4,
     xlab = "", ylab = "")
tmp1 <- threat.per.loss.clumped.VU
lines(c(0,tmp1$x.Median) ~  c(0,tmp1$Group.1), 
      col = "gold", lwd = 2, lty = 1)
tmp1 <- threat.per.loss.clumped.EN
lines(c(0,tmp1$x.Median) ~  c(0,tmp1$Group.1), 
      col = "darkorange", lwd = 2, lty = 1)
tmp1 <- threat.per.loss.clumped.CR
lines(c(0,tmp1$x.Median) ~  c(0,tmp1$Group.1), 
      col = "red2", lwd = 2, lty = 1)
tmp1 <- threat.per.loss.clumped
lines(c(0,tmp1$x.Median,1) ~  c(0,tmp1$Group.1,1), 
      col = "black", lwd = 2, lty = 1)
points(c(threat.end.A, threat.end.A.cats) ~ 
               rep(loss, 4), pch = 21,
       bg = c("black", "gold", "red", "darkorange"), cex = 1.1)
legend("right", c("Threatened", "VU", "EN", "CR"),
       # col = seq_along(combos) + 1,
       lty = 1,
       lwd = 2, bty = "n",
       cex = 1.2, col = c("black", "gold", "darkorange", "red"))
legend(0.8, 0.41, c("Obs. threat."),
       pch = c(21),
       bty = "n", x.intersp= 1.5,
       cex = 1.2, col = "black")
legend("topleft",legend=expression(bold("B")),
       bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.1)

## PANEL C - Comparing other simulated scenarios with one derived form the real AF data
# using only the clumped removal of grid cells
combos <- unique(all.data$combo)
combos <- combos[grepl("clumped_", combos)]
plot(c(0,1) ~ c(1,0), col = "white", las = 1, cex.lab = 1.4,
     xlab = "Habitat loss", ylab = "Threatened species")
legend(0.6, 0.175, expression(bold("Type of community:")),
       #pch = c(19,17),
       bty = "n", #x.intersp= 1.85,
       cex = 1.2, col = "black")
cores.tmp <- c("magenta", "darkgreen", "blue", "black")
legend("bottomright", 
       c("Random", "Aggregate", "Aggregate in blocks", "Atlantic Forest"), 
       col = cores.tmp, 
       lwd = c(2,2,2,1), lty = c(1,1,1,1),
       bty = "n", cex = 1.1)
tmp1 <- threat.per.loss.clumped
lines(c(0,tmp1$x.Median,1) ~  c(0,tmp1$Group.1,1), 
      col = "black", lwd = 1, lty = 1)
arrows(x0 = c(0,as.double(tmp1$Group.1),1)+0.001, 
       y0 = c(0,as.double(tmp1$x.1st.Qu.),1),
       y1 = c(0,as.double(tmp1$x.3rd.Qu.),1),
       col = "black", code = 3, angle = 90, length = 0.025)
# lines(c(0,tmp1$x.1st.Qu.,1) ~  c(0,tmp1$Group.1,1), 
#       col = "black", lwd = 1, lty = 2)
# lines(c(0,tmp1$x.3rd.Qu.,1) ~  c(0,tmp1$Group.1,1), 
#       col = "black", lwd = 1, lty = 2)

for(i in seq_along(combos)) {
        tmp <- all.data[all.data$combo %in% combos[i], ]
        tmp1 <- aggregate(as.double(tmp$threat), list(tmp$loss1), summary)
        lines(c(0,tmp1$x[,3],1) ~ c(0,tmp1$Group.1,1), col = cores.tmp[i], lwd = 2)
        # segments(x0 = as.double(tmp1$Group.1), x1 = as.double(tmp1$Group.1),
        #          y0 = as.double(tmp1$x[,2]), y1 = as.double(tmp1$x[,5]))
}
legend("topleft", expression(bold(C)), bty="n", cex=1.3,
       x.intersp=-0.7, y.intersp=0.1)

# PANEL D - Projections for other regions
hots <- info$hotspot.region[trop_hotspots]
hots <- gsub("West |Western ", "W. ", hots)
hots <- gsub("Eastern ", "E. ", hots)
hots[7] <- "Madagascar, Indian Oc. Isl."

continents <- info$continent[trop_hotspots]
continents <- gsub("North |South ", "", continents)
cores.cont <- factor(continents)
levels(cores.cont) <- c("darkgoldenrod", "forestgreen", "brown", "cadetblue1")
cores.cont <- as.character(cores.cont)
# plot(1:4, 1:4, bg = unique(cores.cont), pch = 21)
pch.hot <- rep(21, length(hots))
pch.hot[continents == "America"] <- c(21, 22, 23, 24, 25, 4)
pch.hot[continents == "Africa"] <- c(21, 22, 23, 24, 25)
pch.hot[continents == "Asia"] <- c(21, 22)
pch.hot[continents == "Oceania"] <- c(21, 22, 23, 24, 25)

spl.pr.low <- predict(spl.low, x=losses[trop_hotspots], na.action=na.exclude)
spl.pr.med <- predict(spl.med, x=losses[trop_hotspots], na.action=na.exclude)
spl.pr.high <- predict(spl.high, x=losses[trop_hotspots], na.action=na.exclude)

plot(c(0,1) ~ c(1,0), col = "white", las = 1, cex.lab = 1.4,
     xlab = "Habitat loss", ylab = "")
tmp1 <- threat.per.loss.clumped
lines(c(0,tmp1$x.Median,1) ~  c(0,tmp1$Group.1,1), 
      col = "black", lwd = 1, lty = 1)
segments(x0 = losses[trop_hotspots], x1 = losses[trop_hotspots],
         y0 = spl.pr.low$y, y1 = spl.pr.high$y, 
         lty = 1, col = cores.cont)
points(spl.pr.med, bg = cores.cont, pch = pch.hot, cex = 1.3)
legend("topleft", expression(bold(D)), bty="n", cex=1.3,
       x.intersp=-0.7, y.intersp=0.1)
legend(0.55, 0.63, expression(bold("Region:")),
       #pch = c(19,17),
       bty = "n", #x.intersp= 1.85,
       cex = 1.2, col = "black")
legend("bottomright", 
       hots, 
       pt.bg = cores.cont, 
       pch = pch.hot,
       bty = "n", cex = 1.1)
dev.off()


## Plotting the loss x community combinations
# par(mfrow=c(1,1), mar = c(5,5,0.5,0.5))
# #simulations
# combos <- unique(all.data$combo)
# plot(c(0,1) ~ c(1,0), col = "white", las = 1, cex.lab = 1.2,
#      xlab = "Habitat loss", ylab = "Threatned species")
# legend("topleft", gsub("_", " x ", combos), col = seq_along(combos) + 1, lwd = 2, bty = "n",
#        cex = 1.1)
# for(i in seq_along(combos)) {
#         tmp <- all.data[all.data$combo %in% combos[i], ]
#         tmp1 <- aggregate(as.double(tmp$threat), list(tmp$loss1), summary)
#         lines(c(0,tmp1$x[,3],1) ~ c(0,tmp1$Group.1,1), col = i + 1, lwd = 2)
#         # segments(x0 = as.double(tmp1$Group.1), x1 = as.double(tmp1$Group.1),
#         #          y0 = as.double(tmp1$x[,2]), y1 = as.double(tmp1$x[,5]))
# }
# # AF - all pops
# combos <- unique(all.data.af$combo)
# # plot(c(0,1) ~ c(1,0), col = "white", las = 1, cex.lab = 1.2,
# #      xlab = "Habitat loss", ylab = "Threatned species")
# legend("bottomright", gsub("_", " x ", combos), 
#        # col = seq_along(combos) + 1, 
#        lty = 1:2,
#        lwd = 2, bty = "n",
#        cex = 1.1)
# cores <- c("black","black")
# for(i in seq_along(combos)) {
#         tmp <- all.data.af[all.data.af$combo %in% combos[i], ]
#         tmp1 <- aggregate(as.double(tmp$threat), list(tmp$loss1), summary)
#         lines(c(0,tmp1$x[,3],1) ~ c(0,tmp1$Group.1,1), 
#               col = cores[i], 
#               lwd = 2,
#               lty = i)
#         segments(x0 = as.double(tmp1$Group.1), x1 = as.double(tmp1$Group.1),
#                   y0 = as.double(tmp1$x[,2]), y1 = as.double(tmp1$x[,5]))
# }
# # AF - endemics
# combos <- unique(all.data.af.end$combo)
# # plot(c(0,1) ~ c(1,0), col = "white", las = 1, cex.lab = 1.2,
# #      xlab = "Habitat loss", ylab = "Threatned species")
# legend("bottomright", gsub("_", " x ", combos),
#        # col = seq_along(combos) + 1,
#        lty = 1:2,
#        lwd = 2, bty = "n",
#        cex = 1.1)
# combos <- unique(all.data.af.end$combo)
# cores <- c("red","red")
# for(i in seq_along(combos)) {
#         tmp <- all.data.af.end[all.data.af.end$combo %in% combos[i], ]
#         tmp1 <- aggregate(as.double(tmp$threat), list(tmp$loss1), summary)
#         lines(c(0,tmp1$x[,3],1) ~ c(0,tmp1$Group.1,1), 
#               col = cores[i], 
#               lwd = 2,
#               lty = i)
#         segments(x0 = as.double(tmp1$Group.1), x1 = as.double(tmp1$Group.1),
#                  y0 = as.double(tmp1$x[,2]), y1 = as.double(tmp1$x[,5]))
# }
