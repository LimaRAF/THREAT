#### QUESTIONS/SUGGESTIONS/FUNCTIONS TO GILES - ConR ####

#### QUESTIONS ####
# 1- Why do you use Resol_sub_pop = 30 in function subpop.comp and 5 in the IUCN.eval? 
# 2- Edit hidden function .subpop.comp to cope with species.specific 'Resol_sub_pop'? See suggested function below.
# 3- Any suggestiong of a good rule-of-thumb for the argument 'nbe.rep.rast.AOO' in AOO.computing?
# 4- Add progress bar to the 'subpop.comp' function?
# 5- Progress bar in IUCN.eval "starts but don't finish" when using parallel=TRUE
# 6- Why in IUCN.eval we need to provide a 'spatialpolygondataframe' in country_map and not in 'EOO.computing'? Does it make a difference if it is a 'spatialpolygondataframe' or a 'spatialpolygon'?
# 7- Option/argument to specify in which folder the shapefiles/figures could be saved?
# 8- Maybe offer individual functions just for the generation of figures?
# 9- I am getting the following error in the foreach loop of the function '.EOO.comp' for a species with 2 occurrences and 1 unique.occur.
#The error emerges for certain subsets of the data and not others. Could not find the root of the problem:
#Error in { : 
#    task 581 failed - "Unable to parse: LINESTRING(-42.50158901 -22.31972091)
#  GEOS reported: "rgeos_readWKT: unable to read wkt""
#I looked into the loop and the problems seems to be that the 'nrow(unique(XY)) < 2' is not recognizing the
# two duplicated coordinates as being duplicated. Then I think that when it try 

#### SUGGESTED FUNCTIONS ####
# Return the number of unique spatial coordinates use to produce the EOO
unique.coords = function(XY, min.dist = 0) {
  if (any(is.na(XY[, c(1:2)]))) {
    print(paste("Skipping", length(which(rowMeans(is.na(XY[,1:2])) > 0)), "occurrences because of missing coordinates for", 
                paste(as.character(unique(XY[which(rowMeans(is.na(XY[,1:2])) > 0), 3])), collapse = " AND ")))
    XY <- XY[which(!is.na(XY[, 1])), ]
    XY <- XY[which(!is.na(XY[, 2])), ]
  }
  XY <- as.data.frame(XY)
  if (any(XY[, 2] > 180) || any(XY[, 2] < -180) || any(XY[,1] < -180) || any(XY[, 1] > 180)) 
    stop("coordinates are outside of expected range")  
  if (ncol(XY) > 2) {
    colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
    XY$tax <- as.character(XY$tax)
    list_data <- split(XY, f = XY$tax)
  } else {
    colnames(XY)[1:2] <- c("ddlat", "ddlon")
    list_data <- list(XY)
  }
  if (is.null(names(list_data))) {
    names_ <- rep(Name_Sp, length(list_data))
  } else { names_ <- names(list_data) }

  ## Total number of occurrences per species
  XY.dt <- data.table(XY)
  setkeyv(XY.dt, "tax") ## setting 'species.correct2' as key to the data.table (makes computations faster)
  all <- as.data.frame(XY.dt[ , .N, by = tax]) ## How many occurrences (i.e. number of rows) per species (by = "taxa")
  
  ## Total number of spatially-unique occurrences per species
  # Exact same coordinates
  non.dup <- as.data.frame(XY.dt[ , uniqueN(paste(ddlat, ddlon, sep="_")), by = tax])
  
  # Same coordinates until a certain distance defined by the user ## TO BE FINISHED
  #f = function(lat, lon, min) {
  #  d = as.double(dist(cbind(lon,lat)))
  #  n = sum(d <= min)
  #}
  #non.dup1 = XY.dt[ , f(ddlat,ddlon,min.dist) , by = tax]
  res = merge(all, non.dup, by = "tax")
  names(res) = c("tax","all.occurs","non.dup.occurs")
  return(res)
}

# Calculate the 1/10th maximum inter-point distance for each species
# quant.max is the upper-quantile of the inter-point distance to be considered as the threshold
subpop.radius = function(XY, factor.div = 10, quant.max = 1) {
  if (any(is.na(XY[, c(1:2)]))) {
    print(paste("Skipping", length(which(rowMeans(is.na(XY[,1:2])) > 0)), "occurrences because of missing coordinates for", 
                paste(as.character(unique(XY[which(rowMeans(is.na(XY[,1:2])) > 0), 3])), collapse = " AND ")))
    XY <- XY[which(!is.na(XY[, 1])), ]
    XY <- XY[which(!is.na(XY[, 2])), ]
  }
  XY <- as.data.frame(XY)
  if (any(XY[, 2] > 180) || any(XY[, 2] < -180) || any(XY[,1] < -180) || any(XY[, 1] > 180)) 
    stop("coordinates are outside of expected range")  
  if (ncol(XY) > 2) {
    colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
    XY$tax <- as.character(XY$tax)
    list_data <- split(XY, f = XY$tax)
  } else {
    colnames(XY)[1:2] <- c("ddlat", "ddlon")
    list_data <- list(XY)
  }
  if (is.null(names(list_data))) {
    names_ <- rep(Name_Sp, length(list_data))
  } else { names_ <- names(list_data) }
  
  ## Getting the maximum inter-point distance
  f = function(lat, lon, fact) {
    x = cbind.data.frame(lon, lat)
    if(dim(x)[1]>3) {
      projWSG = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      projEAC = crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
      #x <- as.data.frame(XY.dt[tax == "Acalypha diversifolia"])
      coordinates(x) <- ~ lon + lat 
      proj4string(x) <- projWSG
      coordEAC <- as.data.frame(coordinates(spTransform(x, projEAC)))
      d.inter = as.character(round(quantile(as.double(dist(coordEAC)/1000), prob=quant.max, na.rm = TRUE)/fact,2))
    } else {
      d.inter = NA_character_
    } 
    return(d.inter)
  }

  ## Getting the maximum inter-point distance using data.table
  XY.dt <- data.table(XY)
  setkeyv(XY.dt, "tax") ## setting 'species.correct2' as key to the data.table (makes computations faster)
  radius <- as.data.frame(XY.dt[ , f(ddlat, ddlon, factor.div) , by = tax]) 
  names(radius) = c("tax","radius")
  return(radius)
}

# Renato's edited version of '.subpop.comp' to include species-specific radius of circles in kilometres
my.subpop.comp = function (XY, Resol_sub_pop = NULL) {
  if (is.null(Resol_sub_pop)) 
    stop("Resol_sub_pop is missing, please provide a value")
  if (any(is.na(XY[, c(1, 2)]))) {
    length(which(rowMeans(is.na(XY[, 1:2])) > 0))
    unique(XY[which(rowMeans(is.na(XY[, 1:2])) > 0), 3])
    print(paste("Skipping", length(which(rowMeans(is.na(XY[, 1:2])) > 0)), "occurrences because of missing coordinates for", 
                paste(as.character(unique(XY[which(rowMeans(is.na(XY[, 1:2])) > 0), 3])), collapse = " AND ")))
    XY <- XY[which(!is.na(XY[, 1])), ]
    XY <- XY[which(!is.na(XY[, 2])), ]
  }
  if (any(XY[, 1] > 180) || any(XY[, 1] < -180) || any(XY[, 2] < -180) || any(XY[, 2] > 180)) 
    stop("coordinates are outside of expected range")
  colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
  XY$tax <- as.character(XY$tax)
  
  ### PART INCLUDED/EDITED BY RENATO ###
  if(class(Resol_sub_pop) == "data.frame") { 
    tmp = dplyr::left_join(XY, Resol_sub_pop, by= "tax")
    XY$spp.radius = as.double(tmp[,dim(tmp)[2]])
  } else {
    XY = XY  
  }
  list_data <- split(XY, f = XY$tax)
  if(class(Resol_sub_pop) == "data.frame") {
    OUTPUT <- lapply(list_data, function(x) ConR:::.subpop.comp(x, Resol_sub_pop = unique(x$spp.radius)))
  } else {
    OUTPUT <- lapply(list_data, function(x) ConR:::.subpop.comp(x, Resol_sub_pop = Resol_sub_pop))
  }
  ### END OF PART INCLUDED/EDITED BY RENATO ###
  
  if (length(OUTPUT) == 1) 
    OUTPUT <- OUTPUT[[1]]
  return(OUTPUT)
}

# Renato's version trying to fix an error I was having (not sure why)
my.locations.comp <- function (XY, method = "fixed_grid", nbe_rep = 0, protec.areas = NULL, 
                               Cell_size_locations = NULL, method_protected_area = "no_more_than_one", 
                               ID_shape_PA = "WDPA_PID", Rel_cell_size = 0.05, parallel = FALSE, 
                               NbeCores = 2, show_progress = TRUE) 
{
  if (!any(class(XY) == "data.frame")) 
    XY <- as.data.frame(XY)
  if (any(XY[, 2] > 180) || any(XY[, 2] < -180) || any(XY[, 
                                                          1] < -180) || any(XY[, 1] > 180)) 
    stop("coordinates are outside of expected range")
  if (method != "fixed_grid" & method != "sliding scale") 
    stop("method should be fixed_grid sliding scale")
  projEAC <- ConR:::.proj_crs()
  coordEAC <- data.frame(matrix(unlist(rgdal::project(as.matrix(XY[, 
                                                                   c(2, 1)]), proj = as.character(projEAC), inv = FALSE)), 
                                ncol = 2), tax = XY[, 3])
  if (any(is.na(coordEAC[, c(1:2)]))) {
    print(paste("Skipping", length(which(rowMeans(is.na(coordEAC[, 
                                                                 1:2])) > 0)), "occurrences because of missing coordinates for", 
                paste(as.character(unique(coordEAC[which(rowMeans(is.na(coordEAC[, 
                                                                                 1:2])) > 0), 3])), collapse = " AND ")))
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 1])), ]
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 2])), ]
  }
  coordEAC$tax <- as.character(coordEAC$tax)
  list_data <- split(coordEAC, f = coordEAC$tax)
  crs_proj <- projEAC
  
  if (is.null(protec.areas)) {
    # if (nrow(coordEAC) > 1) 
    #   pairwise_dist <- stats::dist(coordEAC[, 1:2], upper = F)
    if (any(method == "fixed_grid") & !is.null(Cell_size_locations)) { 
      Resolution <- Cell_size_locations
    } else { 
      Resolution <- 10
    }  
    if (any(method == "sliding scale")) {
      if (nrow(coordEAC) > 1) {
        pairwise_dist <- stats::dist(coordEAC[, 1:2], upper = F)
        Resolution <- max(pairwise_dist) * Rel_cell_size/1000
      }
      else {
        Resolution <- 10
      }
    }
    if (parallel) {
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      message("Parallel running with ", NbeCores, " cores")
      `%d%` <- foreach::`%dopar%`
    }
    else {
      `%d%` <- foreach::`%do%`
    }
    x <- NULL
    if (show_progress) {
      pb <- utils::txtProgressBar(min = 0, max = length(list_data), 
                                  style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, 
                                                       n)
      opts <- list(progress = progress)
    }
    else {
      opts <- NULL
    }
    output <- foreach::foreach(x = 1:length(list_data), .combine = "c", 
                               .options.snow = opts) %d% {
                                 if (!parallel & show_progress) 
                                   utils::setTxtProgressBar(pb, x)
                                 res <- .cell.occupied(size = Resolution, coord = list_data[[x]], 
                                                       nbe_rep = nbe_rep)
                                 names(res) <- c("spatial", "nbe_occ")
                                 res
                               }
    if (parallel) 
      snow::stopCluster(cl)
    if (show_progress) 
      close(pb)
    Locations <- unlist(output[names(output) == "nbe_occ"])
    r2 <- unlist(output[names(output) == "spatial"])
    names(Locations) <- names(r2) <- gsub(pattern = " ", 
                                          replacement = "_", names(list_data))
  }
  
  if (!is.null(protec.areas)) {
    DATA_SF <- as.data.frame(XY[, 1:2])
    colnames(DATA_SF) <- c("ddlat", "ddlon")
    sp::coordinates(DATA_SF) <- ~ddlon + ddlat
    raster::crs(DATA_SF) <- raster::crs(protec.areas)
    Links_NatParks <- sp::over(DATA_SF, protec.areas)
    coordEAC_pa <- coordEAC[!is.na(Links_NatParks[, 1]), ]
    coordEAC_pa <- cbind(coordEAC_pa, id_pa = Links_NatParks[which(!is.na(Links_NatParks[, 
                                                                                         1])), ID_shape_PA])
    LocNatParks <- vector(mode = "numeric", length = length(list_data))
    names(LocNatParks) <- gsub(pattern = " ", replacement = "_", 
                               names(list_data))
    if (nrow(coordEAC_pa) > 0) {
      if (method_protected_area == "no_more_than_one") {
        loc_pa <- by(coordEAC_pa[, c("tax", "id_pa")], 
                     coordEAC_pa[, "tax"], FUN = function(x) length(unique(x$id_pa)))
        names(loc_pa) <- gsub(pattern = " ", replacement = "_", 
                              names(loc_pa))
        LocNatParks[names(LocNatParks) %in% names(loc_pa)] <- loc_pa
        r2_PA <- NA
      }
      else {
        coordEAC_pa$tax <- as.character(coordEAC_pa$tax)
        list_data_pa <- split(coordEAC_pa, f = coordEAC_pa$tax)
        # if (nrow(coordEAC_pa) > 1) 
        #   pairwise_dist_pa <- stats::dist(coordEAC_pa[, 
        #                                               1:2], upper = F)
        if (any(method == "fixed_grid") & !is.null(Cell_size_locations)) { 
          Resolution <- Cell_size_locations
        } else { 
          Resolution <- 10
        } 
        
        if (any(method == "sliding scale")) {
          if (nrow(coordEAC_pa) > 1) {
            pairwise_dist_pa <- stats::dist(coordEAC_pa[, 
                                                        1:2], upper = F)
            Resolution <- max(pairwise_dist_pa)/1000 * 
              Rel_cell_size
          }
          else {
            Resolution <- 10
          }
        }
        if (parallel) {
          cl <- snow::makeSOCKcluster(NbeCores)
          doSNOW::registerDoSNOW(cl)
          message("Parallel running with ", NbeCores, 
                  " cores")
          `%d%` <- foreach::`%dopar%`
        }
        else {
          `%d%` <- foreach::`%do%`
        }
        x <- NULL
        if (show_progress) {
          pb <- utils::txtProgressBar(min = 0, max = length(list_data), 
                                      style = 3)
          progress <- function(n) utils::setTxtProgressBar(pb, 
                                                           n)
          opts <- list(progress = progress)
        }
        else {
          opts <- NULL
        }
        output <- foreach::foreach(x = 1:length(list_data_pa), 
                                   .combine = "c", .options.snow = opts) %d% {
                                     if (!parallel & show_progress) 
                                       utils::setTxtProgressBar(pb, x)
                                     res <- ConR:::.cell.occupied(size = Resolution, coord = list_data_pa[[x]], 
                                                           nbe_rep = nbe_rep)
                                     names(res) <- c("spatial", "nbe_occ")
                                     res
                                   }
        if (parallel) 
          snow::stopCluster(cl)
        if (show_progress) 
          close(pb)
        loc_pa <- unlist(output[names(output) == "nbe_occ"])
        r2_PA <- unlist(output[names(output) == "spatial"])
        names(loc_pa) <- names(r2_PA) <- gsub(pattern = " ", 
                                              replacement = "_", names(list_data_pa))
        LocNatParks[names(LocNatParks) %in% names(loc_pa)] <- loc_pa
      }
    }
    else {
      r2_PA <- NA
    }
    
    coordEAC_not_pa <- coordEAC[is.na(Links_NatParks[, 1]),]
    LocOutNatParks <- vector(mode = "numeric", length = length(list_data))
    names(LocOutNatParks) <- gsub(pattern = " ", replacement = "_", 
                                  names(list_data))
    if (nrow(coordEAC_not_pa) > 0) {
      coordEAC_not_pa$tax <- as.character(coordEAC_not_pa$tax)
      list_data_not_pa <- split(coordEAC_not_pa, f = coordEAC_not_pa$tax)
      # if (nrow(coordEAC_pa) > 1) 
      #   pairwise_dist_not_pa <- stats::dist(coordEAC_not_pa[, 
      #                                                       1:2], upper = F)
      if (any(method == "fixed_grid") & !is.null(Cell_size_locations)) { 
        Resolution <- Cell_size_locations
      } else {
        Resolution <- 10
      }  
      if (any(method == "sliding scale")) {
        if (nrow(coordEAC_pa) > 1) {
          pairwise_dist_not_pa <- stats::dist(coordEAC_not_pa[, 
                                                              1:2], upper = F)
          Resolution <- max(pairwise_dist_not_pa) * Rel_cell_size
        }
        else {
          Resolution <- 10
        }
      }
      if (parallel) {
        cl <- snow::makeSOCKcluster(NbeCores)
        doSNOW::registerDoSNOW(cl)
        message("Parallel running with ", NbeCores, " cores")
        `%d%` <- foreach::`%dopar%`
      }
      else {
        `%d%` <- foreach::`%do%`
      }
      x <- NULL
      if (show_progress) {
        pb <- utils::txtProgressBar(min = 0, max = length(list_data), 
                                    style = 3)
        progress <- function(n) utils::setTxtProgressBar(pb, 
                                                         n)
        opts <- list(progress = progress)
      }
      else {
        opts <- NULL
      }
      output <- foreach::foreach(x = 1:length(list_data_not_pa), 
                                 .combine = "c", .options.snow = opts) %d% {
                                   if (!parallel & show_progress) 
                                     utils::setTxtProgressBar(pb, x)
                                   res <- ConR:::.cell.occupied(size = Resolution, coord = list_data_not_pa[[x]], 
                                                         nbe_rep = nbe_rep)
                                   names(res) <- c("spatial", "nbe_occ")
                                   res
                                 }
      if (parallel) 
        snow::stopCluster(cl)
      if (show_progress & show_progress) 
        close(pb)
      loc_not_pa <- unlist(output[names(output) == "nbe_occ"])
      r2 <- unlist(output[names(output) == "spatial"])
      names(loc_not_pa) <- names(r2) <- gsub(pattern = " ", 
                                             replacement = "_", names(list_data_not_pa))
      LocOutNatParks[names(LocOutNatParks) %in% names(loc_not_pa)] <- loc_not_pa
    }
    else {
      r2 <- NA
    }
  }
  if (!is.null(protec.areas)) 
    return(list(r2, r2_PA, LocNatParks, LocOutNatParks))
  if (is.null(protec.areas)) 
    return(list(r2, Locations))
}


# Renato's version of .IUCN.comp to include pre-calculated inputs and with thresholds as arguments (Global Tree Assessments use slightly different threholds to define LC species)
# NOT CHECKED FOR: !is.null(protec.areas)

DATA <- critB_high
criteria.B = function (DATA, poly_borders = NULL, 
                        #Cell_size_AOO = 2, Cell_size_locations = 10, Resol_sub_pop = 5, method_locations = c("fixed_grid"), Rel_cell_size = 0.05, 
                        protec.areas = NULL, 
                        #exclude.area = FALSE, method_protected_area = "no_more_than_one", ID_shape_PA = "WDPA_PID", buff_width = 0.1, 
                        NamesSp = "species1", 
                        #write_shp = FALSE, 
                        file_name = NULL, add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                        #SubPop = TRUE, MinMax = c(min(list_data[[x]][, 2]), max(list_data[[x]][, 2]), min(list_data[[x]][, 1]), max(list_data[[x]][, 1])), 
                        #alpha = 1, buff.alpha = 0.1, method.range = "convex.hull", nbe.rep.rast.AOO = 0, 
                        #verbose = TRUE, showWarnings = TRUE,
                        EOO.threshold = c(20000, 5000, 100), AOO.threshold = c(2000, 500, 10), Loc.threshold = c(10, 5, 1)) {
  # if (is.null(poly_borders)) {
  #   data("land", package = "ConR", envir = environment())
  #   land <- get("land", envir = environment())
  #   poly_borders = land
  # }
  # if (DrawMap) {
  #   full_poly_borders <- poly_borders
  #   if (!is.null(poly_borders)) 
  #     poly_borders <- raster::crop(poly_borders, extent(MinMax) + 30)
  # }
  #projEAC <- crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
  
  if (is.null(protec.areas)) {
    tmp <- as.data.frame(matrix(NA, 4, dim(DATA)[2]))
    rownames(tmp) <- c("Category_CriteriaB", "Category_code","Category_AOO", "Category_EOO")
    colnames(tmp) <- colnames(DATA)
    Results = rbind.data.frame(DATA, tmp, make.row.names = TRUE, deparse.level = 2)
  } else {
    tmp <- as.data.frame(matrix(NA, 5, dim(DATA)[2]))
    rownames(tmp) <- c("Category_CriteriaB", "Category_code","Ratio_occ_within_PA", "Category_AOO", "Category_EOO")
    colnames(tmp) <- colnames(DATA)
    Results = rbind.data.frame(DATA, tmp, make.row.names = TRUE, deparse.level = 2)
  }
  

  if (any(Results[,"EOO"] < Results[,"AOO"]) )
    Results[,"EOO"][!is.na(Results[,"EOO"]) & Results[,"EOO"] < Results[,"AOO"]] <- 
      Results[,"AOO"][!is.na(Results[,"EOO"]) & Results[,"EOO"] < Results[,"AOO"]]

    if (!is.null(protec.areas)) {
      Results$Ratio_occ_within_PA <- NA
      Results[, "Ratio_occ_within_PA"] <- 
        round((Results[,"Nbe_loc_PA"]/apply(Results[,c("Nbe_loc","Nbe_loc_PA")], 1, sum, na.rm=TRUE)) * 100, 1)
    }  
   
    Rank_EOO = rep(4, dim(DATA)[2])
    Rank_EOO[Results[,"EOO"] < EOO.threshold[1]] <-3
    Rank_EOO[Results[,"EOO"] < EOO.threshold[2]] <-2
    Rank_EOO[Results[,"EOO"] < EOO.threshold[3]] <-1

    Rank_AOO = rep(4, dim(DATA)[2])
    Rank_AOO[Results[,"AOO"] < AOO.threshold[1]] <-3
    Rank_AOO[Results[,"AOO"] < AOO.threshold[2]] <-2
    Rank_AOO[Results[,"AOO"] < AOO.threshold[3]] <-1
    
    Rank_Loc = rep(4, dim(DATA)[2])
    Rank_Loc[apply(Results[,c("Nbe_loc","Nbe_loc_PA")], 1, sum, na.rm=TRUE) < Loc.threshold[1]] <-3
    Rank_Loc[apply(Results[,c("Nbe_loc","Nbe_loc_PA")], 1, sum, na.rm=TRUE) < Loc.threshold[2]] <-2
    Rank_Loc[apply(Results[,c("Nbe_loc","Nbe_loc_PA")], 1, sum, na.rm=TRUE) < Loc.threshold[3]] <-1
    
    Rank_B1a <- apply(rbind(Rank_EOO, Rank_Loc), 2, max, na.rm=TRUE)
    Rank_B2a <- apply(rbind(Rank_AOO, Rank_Loc), 2, max, na.rm=TRUE)
    Rank_CriteriaB <- apply(rbind(Rank_B1a, Rank_B2a), 2, min, na.rm=TRUE)

    Nbe_Loc <- as.numeric(Results[,"Nbe_loc"])
    Cat = rep("LC or NT", dim(DATA)[2])
    Cat[Rank_CriteriaB == 1] <- "CR"
    Cat[Rank_CriteriaB == 2] <- "EN"
    Cat[Rank_CriteriaB == 3 & Nbe_Loc > 0 & Nbe_Loc < 11] <- "VU"
    
    Cat_Code <- rep(NA, dim(DATA)[1])
    if(any(Rank_B1a > Rank_B2a))    
      Cat_Code[Rank_B1a > Rank_B2a] <- paste(Cat[Rank_B1a > Rank_B2a], "B2a")
    if(any(Rank_B1a < Rank_B2a)) 
      Cat_Code[Rank_B1a < Rank_B2a] <- paste(Cat[Rank_B1a < Rank_B2a], "B1a")
    if (any(Rank_B1a == Rank_B2a)) 
      Cat_Code[Rank_B1a == Rank_B2a & Rank_B1a != 4] <- paste(Cat[Rank_B1a == Rank_B2a & Rank_B1a != 4], "B1a+B2a") 

    
    if (!is.null(protec.areas)) {
      if (any(as.numeric(Results[,"Ratio_occ_within_PA"]) == 100)) {
        Cat[as.numeric(Results[, "Ratio_occ_within_PA"]) %in% 100] <- "LC or NT"
        Cat_Code[as.numeric(Results[, "Ratio_occ_within_PA"]) %in% 100] <- NA
      } 
    }   
    
    Results$Category_EOO <- NA
    Results$Category_AOO <- NA
    
    if (any(Rank_B2a == 1)) 
      Results[Rank_B2a == 1, "Category_AOO"] <- "CR"
    if (any(Rank_B2a == 2)) 
      Results[Rank_B2a == 2, "Category_AOO"] <- "EN"
    if (any(Rank_B2a == 3)) 
      Results[Rank_B2a == 3, "Category_AOO"] <- "VU"
    if (any(Rank_B2a > 3)) 
      Results[Rank_B2a > 3, "Category_AOO"] <- "LC or NT"
    if (any(Rank_B1a == 1)) 
      Results[Rank_B1a == 1, "Category_EOO"] <- "CR"
    if (any(Rank_B1a == 2)) 
      Results[Rank_B1a == 2, "Category_EOO"] <- "EN"
    if (any(Rank_B1a == 3)) 
      Results[Rank_B1a == 3, "Category_EOO"] <- "VU"
    if (any(Rank_B1a > 3)) 
      Results[Rank_B1a > 3, "Category_EOO"] <- "LC or NT"

    if (!is.null(protec.areas)) {
      if (any(as.numeric(Results[,"Ratio_occ_within_PA"]) %in% 100)) {
        Results[as.numeric(Results[,"Ratio_occ_within_PA"]) %in% 100, "Category_CriteriaB"] <- "LC or NT"
      } else {
        if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3])) {
          Results["Category_AOO", as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc",
                      ]) <= Loc.threshold[3]] <- Results["Category_CriteriaB", as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3]] <- "CR"
        } else {
          if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2])) {
            Results["Category_AOO", as.numeric(Results[,"AOO" 1]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", 
                      ]) <= Loc.threshold[2]] <- Results["Category_CriteriaB", as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2]] <- "EN"
          } else {
            if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1])) {
              Results["Category_AOO", as.numeric(Results[,"AOO" 1]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", 
                      ]) <= Loc.threshold[1]] <- Results["Category_CriteriaB", as.numeric(Results[,"AOO" ]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
            } else {
              Results["Category_AOO", ] <- Results["Category_CriteriaB", ] <- "LC or NT"
            }
          }
        }
      }
    } else {
      if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3])) {
        Results["Category_AOO", as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3]] <- "CR"
        Results["Category_CriteriaB", as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3]] <- "CR"
      } else {
        if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2])) {
          Results["Category_AOO", as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2]] <- "EN"
          Results["Category_CriteriaB", as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2]] <- "EN"
        } else {
          if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[1] &  as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1])) {
            Results["Category_AOO", as.numeric(Results[,"AOO" ]) < AOO.threshold[1] &  as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
            Results["Category_CriteriaB", as.numeric(Results[,"AOO" ]) < AOO.threshold[1] &  as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
          } else {
            #Results["Category_AOO", ] <- "LC or NT"
            #Results["Category_CriteriaB", ] <- "LC or NT"
          }
        }
      }
    }
    if (any(is.na(Results["Category_AOO", ]))) {
      if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[3])) {
        Results["Category_AOO", as.numeric(Results[,"AOO" ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[3]] <- "CR"
      } else {
        if (any(as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[2])) {
          Results["Category_AOO", as.numeric(Results[,"AOO" ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[2]] <- "EN"
        } else {
          if (any(as.numeric(Results[,"AOO" 1]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1])) {
            Results["Category_AOO", as.numeric(Results[,"AOO" 1]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
          } else {
            #Results["Category_AOO", ] <- "LC or NT"
          }
        }
      }
    }
    #Results["Category_code", ] <- paste(Results["Category_CriteriaB", ], "B2a")

    # if (showWarnings) 
    #   warning(paste("EOO statistic is not computed for", 
    #                 NamesSp, "because there is less than 3 records"))
  # }
  # if (DrawMap) {
  #   if (!map_pdf) {
  #     if (!is.null(file_name)) {
  #       NAME_FILE <- paste(file_name, gsub(" ", replacement = "_", as.character(NamesSp)), sep = "")
  #     } else {
  #       NAME_FILE <- paste("IUCN_", gsub(" ", replacement = "_",  as.character(NamesSp)), sep = "")
  #     }
  #     FILE_NAME <- ifelse(!is.null(file_name), file_name,  "IUCN_")
  #     dir.create(file.path(paste(getwd(), paste("/", FILE_NAME, "_results_map", sep = ""), sep = "")), showWarnings = FALSE)
  #   }
  #   if (!map_pdf) 
  #     png(paste(file.path(paste(getwd(), paste("/", FILE_NAME,  "_results_map", sep = ""), sep = "")), "/", NAME_FILE, 
  #               ".png", sep = ""), width = 2000, height = 2000)
  #   par(mar = c(10, 12, 10, 2), xpd = FALSE, las = 1)
  #   if (add.legend & !any(colnames(DATA) == "coly")) 
  #     nf <- layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), c(4, 1.5), c(4, 1.5))
  #   if (any(colnames(DATA) == "coly") & add.legend) 
  #     nf <- layout(matrix(c(1, 1, 1, 1, 1, 1, 2, 3, 4), 3, 3, byrow = TRUE), c(2, 1.5, 1.5), c(4, 1.5, 1.5))
  #   if (!is.null(protec.areas)) {
  #     if (LocOutNatParks == 0) {
  #       plot(poly_borders, xlim = c(range(XY[, 1])[1] - 1, range(XY[, 1])[2] + 1), ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1), axes = FALSE, 
  #            xlab = "", ylab = "")
  #     } else {
  #       if (LocOutNatParks == 1) {
  #         plot(r2, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.2), xlim = c(range(XY[, 1])[1] - 
  #                                                     1, range(XY[, 1])[2] + 1), ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1))
  #       } else {
  #         plot(r2, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.2), xlim = c(range(XY[, 1])[1] - 
  #                                                     1, range(XY[, 1])[2] + 1), ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1))
  #       }
  #     }
  #   } else {
  #     plot(r2, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.2), xlim = c(range(XY[, 1])[1] - 1, 
  #                                               range(XY[, 1])[2] + 1), ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1))
  #   }
  #   if (SubPop) 
  #     plot(SubPopPoly, add = T, border = "black", lwd = 2,  lty = 1)
  #   if (!is.null(protec.areas)) {
  #     if (LocNatParks > 0) {
  #       if (method_protected_area != "no_more_than_one") {
  #         plot(r2_PA, add = T, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
  #       }
  #     }
  #   }
  #   if (!is.null(p1) & draw.poly.EOO) 
  #     plot(p1, add = T, col = rgb(red = 0.2, green = 0.2, blue = 0.2, alpha = 0.1))
  #   plot(poly_borders, axes = FALSE, lty = 1, add = T, lwd = 1)
  #   if (!is.null(protec.areas)) 
  #     plot(protec.areas, add = T, col = rgb(red = 0.2, green = 0.2, blue = 0.2, alpha = 0.05), lty = 2)
  #   if (!is.null(protec.areas)) {
  #     colnames(XY) <- c("ddlon", "ddlat")
  #     XY_sp <- XY[which(is.na(Links_NatParks[, 1])), ]
  #     if (nrow(XY_sp) > 0) {
  #       coordinates(XY_sp) <- ~ddlon + ddlat
  #       plot(XY_sp, pch = 19, cex = 2, col = "black",  add = T)
  #     }
  #     XY_sp <- XY[which(!is.na(Links_NatParks[, 1])), ]
  #     if (nrow(XY_sp) > 0) {
  #       coordinates(XY_sp) <- ~ddlon + ddlat
  #       plot(XY_sp, pch = 19, cex = 2, col = "blue",  add = T)
  #     }
  #   }  else {
  #     colnames(XY) <- c("ddlon", "ddlat")
  #     XY_sp <- XY
  #     coordinates(XY_sp) <- ~ddlon + ddlat
  #     plot(XY_sp, pch = 19, cex = 2, col = "black", add = T)
  #   }
  #   axis(1, outer = FALSE, cex.axis = 3, tick = FALSE, line = 1.5)
  #   axis(1, outer = FALSE, labels = FALSE, cex.axis = 3, tick = TRUE, line = 0)
  #   axis(2, outer = FALSE, cex.axis = 3, tick = FALSE, line = 1.5)
  #   axis(2, outer = FALSE, labels = FALSE, cex.axis = 3, tick = TRUE, line = 0)
  #   box()
  #   if (Results["Nbe_loc", 1] > 1) {
  #     xlim <- par("xaxp")[1:2]
  #     xlim <- abs(xlim[1] - xlim[2])
  #     border_to_center <- as.data.frame(matrix(NA, 2, 2))
  #     border_to_center[, 1] <- c(xlim/10, 0)
  #     border_to_center[, 2] <- c(0, 0)
  #     scaleBAR <- round(matrix(unlist(rgdal::project(as.matrix(border_to_center), proj = as.character(projEAC), inv = F)), ncol = 2)/1000, 0)[1, 1]
  #   }  else {
  #     scaleBAR <- Resolution/1000
  #   }
  #   scalebar(scaleBAR, type = "bar", below = "kilometres", cex = 2.5)
  #   mtext(NamesSp, side = 3, cex = 3, line = 3)
  #   if (any(colnames(DATA) == "higher.tax.rank")) 
  #     mtext(DATA[which(DATA[, 3] == NamesSp), "higher.tax.rank"][1], side = 3, cex = 3, line = 0.4)
  #   if (add.legend) {
  #     par(mar = c(1, 1, 1, 1), xpd = T)
  #     plot(1:10, 1:10, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #     if (is.null(protec.areas)) {
  #       legend(1, 10, c(paste("EOO=", ifelse(!is.na(Results["EOO", 1]), round(as.numeric(Results["EOO", 1]), 1),
  #                                             NA), "km2"), paste("AOO (grid res.", Cell_size_AOO, "km)=", format(Results[,"AOO" 1], scientific = 5), 
  #                                             "km2"), paste("Number of unique occurrences=", Results["Nbe_unique_occ.", 1]), paste("Number of sub-populations (radius", 
  #                                             Resol_sub_pop, "km)=", Results["Nbe_subPop", 1]), paste("Number of locations (grid res.:", 
  #                                             round(Resolution/1000, 1), " km)", "=", Results["Nbe_loc",1]), paste("IUCN category according to criterion B:", 
  #                                             Results["Category_CriteriaB", 1])), cex = 3.5, bg = grey(0.9))
  #     }
  #     if (!is.null(protec.areas)) {
  #       legend(1, 10, c(paste("EOO=", ifelse(!is.na(Results["EOO", 1]), round(as.numeric(Results["EOO", 1]), 1), 
  #                                             NA), "km2"), paste("AOO (grid res.", Cell_size_AOO, "km)=", format(Results[,"AOO" 1], scientific = 5), 
  #                                             "km2"), paste("Number of unique occurrences=", Results["Nbe_unique_occ.", 1]), paste("Number of sub-populations (radius", 
  #                                             Resol_sub_pop, "km)=", Results["Nbe_subPop", 1]), paste("Number of locations (grid res.:", 
  #                                             round(Resolution/1000, 1), " km)", "=", Results["Nbe_loc",1]), paste("Number of occupied protected areas=", 
  #                                             Results["Nbe_loc_PA", 1]), paste("IUCN category according to criterion B:", 
  #                                             Results["Category_CriteriaB", 1]), paste("Proportion of occurences within protected areas"), 
  #                                             Results["Ratio_occ_within_PA", 1]), cex = 3.5, bg = grey(0.9))
  #     }
  #     par(mar = c(4, 1, 1, 1))
  #     plot(full_poly_borders, lty = 1, lwd = 1, axes = FALSE)
  #     points(XY[, 1], XY[, 2], pch = 8, cex = 2, col = "red")
  #   }
  #   if (any(colnames(DATA) == "coly") & add.legend) {
  #     par(mar = c(12, 6, 1, 2), las = 2, yaxs = "r", xpd = FALSE)
  #     subdata <- DATA[which(DATA[, "tax"] == NamesSp), "coly"]
  #     if ((sum(subdata, na.rm = T)) > 0) {
  #       plot(table(subdata), col = "grey", ylab = " ", xlab = " ", cex.lab = 4, cex.axis = 4, axes = F)
  #       axis(1, outer = FALSE, cex.axis = 3, tick = FALSE, line = 1.5)
  #       axis(1, outer = FALSE, labels = FALSE, cex.axis = 3, tick = TRUE, line = 0)
  #       axis(2, outer = FALSE, cex.axis = 3, tick = FALSE, line = 1.5)
  #       axis(2, outer = FALSE, labels = FALSE, cex.axis = 3, tick = TRUE, line = 0)
  #     }
  #   }
  #   if (!map_pdf) 
  #     dev.off()
  # }
  # if (write_shp) {
  #   dir.create(file.path(paste(getwd(), "/shapesIUCN", sep = "")), showWarnings = FALSE)
  #   if (!is.null(p1)) {
  #     if (length(list.files(paste(getwd(), "/shapesIUCN",  sep = ""))) > 0) {
  #       if (length(grep(paste(NamesSp, "_EOO_poly", sep = ""), 
  #                       unique(sub("....$", "", list.files(paste(getwd(), "/shapesIUCN", sep = "")))))) > 0) {
  #         FILES <- list.files(paste(getwd(), "/shapesIUCN", sep = ""), full.names = TRUE)
  #         file.remove(FILES[grep(paste(NamesSp, "_EOO_poly", sep = ""), FILES)])
  #       }
  #     }
  #     NAME <- names(p1)
  #     p1@polygons[[1]]@ID <- "1"
  #     ConvexHulls_poly_dataframe <- SpatialPolygonsDataFrame(p1, data = as.data.frame(names(p1)))
  #     colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = "")
  #     writeOGR(ConvexHulls_poly_dataframe, "shapesIUCN", paste(NamesSp, "_EOO_poly", sep = ""), driver = "ESRI Shapefile")
  #   }
  #   if (SubPop) {
  #     if (length(list.files(paste(getwd(), "/shapesIUCN", sep = ""))) > 0) {
  #       if (length(grep(paste(NamesSp, "_subpop_poly", sep = ""), unique(sub("....$", "", list.files(paste(getwd(), 
  #                       "/shapesIUCN", sep = "")))))) > 0) {
  #         FILES <- list.files(paste(getwd(), "/shapesIUCN", sep = ""), full.names = TRUE)
  #         file.remove(FILES[grep(paste(NamesSp, "_subpop_poly", sep = ""), FILES)])
  #       }
  #     }
  #     NAME <- names(SubPopPoly)
  #     SubPopPoly@polygons[[1]]@ID <- "1"
  #     ConvexHulls_poly_dataframe <- SpatialPolygonsDataFrame(SubPopPoly, data = as.data.frame(names(SubPopPoly)))
  #     colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = "")
  #     writeOGR(ConvexHulls_poly_dataframe, "shapesIUCN", paste(NamesSp, "_subpop_poly", sep = ""), driver = "ESRI Shapefile")
  #   }
  # }
  # if (SubPop) {
  #   OUTPUT <- list(Results, p1, SubPopPoly)
  #   names(OUTPUT) <- c("Results", "spatialPoly_EOO", "spatialPoly_subpop")
  # } else {
  #   OUTPUT <- list(Results, p1)
  #   names(OUTPUT) <- c("Results", "spatialPoly_EOO")
  # }
  # return(OUTPUT)
  Results1 = as.data.frame(t(Results), stringsAsFactors = FALSE)  
  return(Results1)
}

