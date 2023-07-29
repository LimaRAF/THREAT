#### INTERNAL FUNCTIONS ####

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
  
  res = merge(all, non.dup, by = "tax")
  names(res) = c("tax","all.occurs","non.dup.occurs")
  return(res)
}

# Calculate the 1/10th maximum inter-point distance for each species
# quant.max is the upper-quantile of the inter-point distance to be considered as the threshold
my.subpop.radius = function(XY, factor.div = 10, quant.max = 1) {
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
      #projEAC = crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
      projSIRGAS <- CRS("+init=epsg:5641")
      #x <- as.data.frame(XY.dt[tax == "Acalypha diversifolia"])
      coordinates(x) <- ~ lon + lat 
      proj4string(x) <- projWSG
      coordEAC <- as.data.frame(coordinates(spTransform(x, projSIRGAS)))
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

# Renato's version; trying to fix an error I was having (not sure why)
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
  projEAC <- my.proj_crs()
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
                                 
                                 source("./R/other_functions.R")
                                 
                                 res <- my.cell.occupied(size = Resolution, coord = list_data[[x]], 
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
                                     
                                     source("./R/other_functions.R")
                                     
                                     res <- my.cell.occupied(size = Resolution, coord = list_data_pa[[x]], 
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
                                   
                                   source("./R/other_functions.R")
                                   
                                   res <- my.cell.occupied(size = Resolution, coord = list_data_not_pa[[x]], 
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


## Renato's first version of the function now in ConR
my.over.valid.poly <- function(poly, points, proj_user = NULL, min.dist = 0.1, value = "dist") {
  
  if (all(!poly$tax %in% unique(points$tax)))
    return(rep(NA, nrow(points)))
  
  poly <- poly[poly$tax %in% unique(points$tax),]
  
  if (is.null(proj_user)) {
    proj_user <- 3857
    warning("no projected coordinate reference system provided by the user: assuming WGS 84 Pseudo-Mercator (see https://epsg.io)"
    )
  }
  
  poly_sf <- sf::st_as_sf(poly)
  if (is.na(sf::st_crs(poly_sf)[[1]]))
    sf::st_crs(poly_sf) <- proj_user
  
  points_sf <- sf::st_as_sf(points, coords = c("ddlon", "ddlat"))
  if (is.na(sf::st_crs(points_sf)[[1]]))
    sf::st_crs(points_sf) <- sf::st_crs(poly_sf)
  
  points_sf <- sf::st_transform(points_sf, crs = proj_user)
  poly_sf <- sf::st_transform(poly_sf, crs = proj_user)
  
  # Solution for multiple species at once# solution for single species at a time
  if (nrow(poly_sf) == 1) {
    
    dist_poly <- sf::st_distance(points_sf, poly_sf)
    dists <- round(as.double(dist_poly[,1]), 2)
    
    
    # Solution for multiple species EOO at once
  } else {
    dist_poly <- sf::st_distance(points_sf, poly_sf)
    dist_poly <- matrix(dist_poly, nrow = dim(dist_poly)[1], ncol = dim(dist_poly)[2])
    colnames(dist_poly) <- as.character(poly$tax)
    row.names(dist_poly) <- as.character(points$tax)
    j <- match(rownames(dist_poly), colnames(dist_poly))
    dists <- round(diag(dist_poly[,j]),2)
  }
  
  flag <- dists < min.dist
  
  if(value == "flag") {
    return(flag)
  }
  
  if(value == "dist") {
    
    coords <- as.data.frame(sf::st_coordinates(points_sf))
    true.ids <- points_sf$classes >= max(points_sf$classes, na.rm = TRUE)
    coords.true <- coords[true.ids,]
    coords.true$tax <- points_sf$tax[true.ids]
    coords.true <- coords.true[!is.na(coords.true[, 1]), ]
    
    max_dists <- tapply(1:nrow(coords.true), coords.true$tax, 
                        function(x) quantile(fields::rdist(coords.true[x, 1:2], compact = TRUE), prob=0.95 , na.rm = TRUE))
    tmp <- dplyr::left_join(data.frame(dist = dists, tax = points_sf$tax, stringsAsFactors = FALSE),
                            data.frame(
                              #med_dists = as.double(med_dists),
                              inter_dists = as.double(max_dists),
                              tax = names(max_dists), stringsAsFactors = FALSE), 
                            by = 'tax')
    dist <- round(tmp$d / tmp$inter_dists, 5)
    return(dist)
  }
}  

## Version from the package's previous version
my.EOO.comp <-  function(XY,
                         exclude.area = FALSE,
                         buff_width = 0.1,
                         country_map = NULL,
                         Name_Sp = "tax",
                         method.range = "convex.hull",
                         alpha = 1,
                         buff.alpha = 0.1,
                         method.less.than3 = "not comp",
                         mode = "spheroid",
                         proj_type = "cea"
) {
  # , verbose=TRUE
  
  XY <- 
    my.coord.check(XY = XY, listing = FALSE)
  
  ### Getting by default land map if poly_borders is not provided
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  }else{
    
    if(any(grepl('sf', class(country_map))))
      country_map <- 
        as(country_map, "Spatial")
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
    country_map <- 
      as(country_map, "sf")
  }
  
  ### Checking if the method of calculating EOO has been chosen
  # if (!convex.hull & !alpha.hull)
  #   stop("alpha.hull and convex.hull are both FALSE, one should TRUE")
  
  if (nrow(unique(XY)) > 1)
    if (max(dist(XY[, 2]), na.rm = T) >= 180)
      warning(
        paste(
          "Occurrences spans more than 180 degrees longitude for species",
          as.character(Name_Sp),
          ". EOO unlikely reliable, check the projection for a proper estimation and used a 'planar' (projected) mode"
        )
      )
  
  # if(alpha.hull) {
  #   convex.hull <- FALSE
  # }
  
  ## Check if there are less than 3 unique occurrences
  if (nrow(XY) < 3) {
    ## if there is only one occurrence, EOO is NA
    if (nrow(XY) < 2) {
      
      EOO <- NA
      message(
        paste(
          "\nEOO parameter cannot be estimated for",
          as.character(Name_Sp),
          "because there is only 1 unique occurrence"
        )
      )
      
    } else {
      if (method.less.than3 == "arbitrary") {
        
        projEAC <- my.proj_crs(proj_type = proj_type)
        
        coordEAC <-
          sf_project(
            from = sf::st_crs(4326),
            to =
              sf::st_crs(projEAC),
            pts = XY[, c(1, 2)]
          )
        
        EOO <-
          as.numeric(dist(coordEAC / 1000) * 0.1 * dist(coordEAC / 1000))
      }
      
      if (method.less.than3 == "not comp") {
        ## if there are two unique occurences, EOO is not computed neither
        message(
          paste(
            "\nEOO parameter cannot be estimated for",
            as.character(Name_Sp),
            "because there is less than 3 unique occurrences"
          )
        )
        EOO <- NA
      }
    }
    
    OUTPUT <- list(EOO = EOO, spatial.polygon = NA)
    
  } else {
    
    ### Checking if all occurrences are on a straight line
    if (length(unique(XY[, 1])) == 1 ||
        length(unique(XY[, 2])) == 1 ||
        round(abs(cor(XY[, 1], XY[, 2])), 6) == 1) {
      ## If so, a straight line is built and a buffer of buff_width is added
      
      message(
        paste(
          "\nOccurrences of",
          as.character(Name_Sp),
          "follow a straight line, thus EOO is based on an artificial polygon using buff_width"
        )
      )
      
      hpts <- unique(XY[, c(2, 1)])
      POLY <- "LINESTRING("
      for (Z in 1:dim(hpts)[1]) {
        POLY <- paste(POLY, hpts[Z, 1], " ", hpts[Z, 2], sep = "")
        if (Z != dim(hpts)[1])
          POLY <- paste(POLY, ", ", sep = "")
        if (Z == dim(hpts)[1])
          POLY <- paste(POLY, ")", sep = "")
      }
      p1 <- rgeos::readWKT(POLY)
      
      p1 <- rgeos::readWKT(POLY)
      sp::proj4string(p1) <- CRS(SRS_string='EPSG:4326')
      
      # crs <- CRS("+proj=longlat +datum=WGS84")
      # crs(p1) <- crs
      
      p1 <-
        suppressWarnings(geosphere::makeLine(p1)) ### Add vertices to line
      
      p1 <-
        suppressWarnings(rgeos::gBuffer(p1, width = buff_width)) ### Add buffer to line
      
      ## If exclude.area is TRUE
      if (exclude.area) {
        
        p1_sf <- as(p1, "sf")
        
        p1 <-
          suppressWarnings(suppressMessages(st_union(
            st_intersection(p1_sf, country_map)
          )))
        
        sf::st_crs(p1) <-
          "+proj=longlat +datum=WGS84"
        
        p1 <- as(p1, "Spatial")
        
        EOO <- 
          suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
        
      } else {
        
        EOO <- suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
        
      }
      
    } else {
      
      if (method.range == "alpha.hull") {
        
        ### work around to avoid bug appearing randomly
        # cont_try <- TRUE
        # while(cont_try) {
        
        p1 <-
          alpha.hull.poly(
            XY = XY[, c(2, 1)],
            alpha = alpha,
            buff = buff.alpha,
            exclude.area = exclude.area,
            poly_exclude = as(country_map, "sf"),
            mode = mode,
            proj_type = proj_type)
        
        #   if (!grepl("trye-error", class(p1)))
        #     cont_try <- FALSE
        # }
        
      }
      
      
      if (method.range == "convex.hull") {
        p1 <- try(
          Convex.Hull.Poly(XY = XY[, c(2, 1)],
                           mode = mode,
                           proj_type = proj_type, 
                           exclude.area = exclude.area,
                           poly_exclude = country_map), TRUE)
        if (class(p1) == "try-error") {
          sf::sf_use_s2(FALSE)
          country_map1 <- sf::st_make_valid(country_map)
          p1 <- try(
            Convex.Hull.Poly(XY = XY[, c(2, 1)],
                             mode = mode,
                             proj_type = proj_type, 
                             exclude.area = exclude.area,
                             poly_exclude = country_map1), TRUE)
          sf::sf_use_s2(TRUE)
          
          if (class(p1) == "try-error")
            return(list(EOO = NA, spatial.polygon = NA))
        }
      }
      
      if(any(class(p1) == "SpatialPolygons") | any(class(p1) == "sfc")) {
        if (mode == "spheroid") {
          EOO <-
            suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
        }
        
        if (mode == "planar") {
          
          EOO <-
            as.numeric(st_area(p1)) / 1000000
          
          p1 <- 
            as(st_transform(p1, 4326), "Spatial")
        }
        
      } else  {
        
        EOO <- NA
        
      }
    }
    
    OUTPUT <- list(EOO = EOO, spatial.polygon = p1)
  }
  
  digits <-
    c(6, 5, 4, 3, 2, 1, 0)[findInterval(OUTPUT$EOO, 
                                        c(0, 0.0001, 0.01, 0.1, 1, 10, 30000, Inf))]
  
  OUTPUT$EOO <- round(OUTPUT$EOO, digits)
  
  # if(verbose) cat(" ",paste(Name_Sp,"EOO comp."))
  
  return(OUTPUT)
}


## Renato's first version of the function now in ConR
my.EOO.sensitivity <- function(XY,
                               levels.order = NULL,
                               occ.based = TRUE,
                               value = "dist",
                               proj_user = NULL,
                               exclude.area = FALSE,
                               country_map = NULL,
                               alpha = 1,
                               buff.alpha = 0.1,
                               method.range = "convex.hull",
                               #Name_Sp = "species1",
                               buff_width = 0.1,
                               method.less.than3 = "not comp",
                               write_results = FALSE,
                               file.name = "EOO.sensitivity.results",
                               parallel = FALSE,
                               NbeCores = 2,
                               show_progress = TRUE
){ 
  
  XY$recordID = 1:dim(XY)[1]
  
  if (any(is.na(XY[, c(1:2)]))) {
    print(paste(
      "Skipping",
      length(which(rowMeans(is.na(
        XY[, 1:2]
      )) > 0)) ,
      "occurrences because of missing coordinates for",
      # if(verbose)
      paste(as.character(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                           0), 3])), collapse = " AND ")
    ))
    XY1 <- XY[which(!is.na(XY[, 1])), ]
    XY1 <- XY1[which(!is.na(XY1[, 2])), ]
    XY1 <- as.data.frame(XY1)
  } else {
    XY1 <- as.data.frame(XY)
  }
  
  if (exclude.area & is.null(country_map))
    stop("exclude.area is TRUE but no country_map is provided")
  
  if (buff_width > 80)
    stop("buff_width has an unrealistic value")
  
  if (any(XY1[, 2] > 180) ||
      any(XY1[, 2] < -180) ||
      any(XY1[, 1] < -180) ||
      any(XY1[, 1] > 180))
    stop("coordinates are outside of the expected range")
  
  if(length(unique(XY1[,4])) <2)
    stop("there is only one class of confidence level")
  
  if(!is.null(country_map))
    country_map <- suppressWarnings(rgeos::gBuffer(country_map, byid=TRUE, width=0))
  
  if (method.range == "convex.hull") {
    convex.hull = TRUE
    alpha.hull = FALSE
  }
  
  if (method.range == "alpha.hull") {
    convex.hull = FALSE
    alpha.hull = TRUE
  }
  
  if (ncol(XY1) > 4) {
    colnames(XY1)[1:4] <- c("ddlat", "ddlon", "tax", "valid")
    XY1$valid <- as.character(XY1$valid)
    XY1$tax <- as.character(XY1$tax)
    XY1 <- XY1[, c("ddlat", "ddlon", "tax", "valid", "recordID")]
  } else{
    colnames(XY1)[1:3] <- c("ddlat", "ddlon", "valid")
    XY1$valid <- as.character(XY1$valid)
    XY1$tax <- "Species 1"
    XY1 <- XY1[, c("ddlat", "ddlon", "tax", "valid", "recordID")]
  }
  
  levels.order <- levels.order[levels.order %in% unique(XY1$valid)]
  XY1 <- XY1[XY1$valid %in% levels.order, ]
  n.levels <- length(levels.order)
  
  if(occ.based) {
    export_shp = c(rep(FALSE, n.levels - 1), TRUE)
  } else {
    export_shp = c(rep(FALSE, n.levels - 1), FALSE) 
  }
  
  ## Obtaining the data for each class of confidence level
  XY1$classes <- as.double(factor(XY1$valid, levels = levels.order, labels = 1:n.levels))
  
  XY.list <- sp_names <- vector("list", n.levels)
  for(i in 1:n.levels) {
    tmp <- XY1[XY1$classes >= i, ]
    XY.list[[i]] <- tmp
    sp_names[[i]] <- sort(unique(tmp$tax))
    rm(tmp)
  }
  names(XY.list) <- names(sp_names) <- paste0("level.", 1:n.levels)
  
  ## Obtaining EOO for each taxon and class of confidence level
  result <- vector("list", n.levels) 
  names(result) <- paste0("level.", 1:n.levels)  
  cat("Starting the EOO analysis for each species and confidence levels...", sep= "\n")
  for(i in 1:length(result)) {
    result[[i]] <- my.EOO.computing(XY = XY.list[[i]],
                                    exclude.area = exclude.area,
                                    buff_width = buff_width,
                                    country_map = country_map,
                                    write_shp = FALSE,
                                    #Name_Sp = names_list[[i]],
                                    method.range = method.range,
                                    alpha = 1,
                                    buff.alpha = 0.1,
                                    method.less.than3 = "not comp",
                                    #alpha.hull = alpha.hull,
                                    #convex.hull = convex.hull,
                                    write_results = FALSE,
                                    export_shp = export_shp[i],
                                    parallel = parallel,
                                    NbeCores = NbeCores,
                                    show_progress = show_progress)
  }
  
  if(occ.based) {
    eoos <- do.call(rbind.data.frame, result[[n.levels]][grepl("EOO", names(result[[n.levels]]))])
    dimnames(eoos) <- list(sp_names[[n.levels]], "EOO")
    shps <- result[[n.levels]][grepl("spatial.polygon", names(result[[n.levels]]))]
    for (i in 1:length(shps)) {
      shp.i <- shps[[i]]
      if (class(shp.i) == "SpatialPolygons")
        methods::slot(methods::slot(shps[[i]], "polygons")[[1]], "ID") <- 
          sp_names[[n.levels]][i]
    }
    shps <- shps[!is.na(shps)]
    shps <- do.call(rbind, shps)
    shps_df <- sp::SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
    shps_df$tax <- as.character(shps_df$tax) 
    result[[n.levels]] <- eoos
    rm(eoos, shps)
  }
  
  for (i in 1:length(result))
    result[[i]]$n.occs <- as.double(table(XY.list[[i]]$tax)) 
  
  Results_short <- merge(result[[1]], result[[2]], by="row.names", 
                         all = TRUE, suffixes = c(".conf1",".conf2"))
  if(n.levels > 2) {
    for(i in 3:n.levels)
      Results_short <- merge(Results_short, result[[i]], 
                             all = TRUE, by.x="Row.names", by.y="row.names", suffixes = c("", paste0(".conf",i)))
  }
  names(Results_short) <- c("Species", paste0(rep(c("EOO.level.","Occs.level."),n.levels), rep(1:n.levels, each = 2)))
  
  result1 <- Results_short[, grepl("EOO", names(Results_short))]
  result1 <- do.call(rbind.data.frame,
                     lapply(1:nrow(result1),  
                            function(x) round((result1[x, 1:(n.levels-1)] - result1[x, n.levels]) / 
                                                result1[x, n.levels], 5)))
  names(result1) <- paste0("EOO.increase.", 1:(n.levels-1))
  Results_short <- cbind.data.frame(Results_short,
                                    result1, stringsAsFactors = FALSE)
  Results_short <- Results_short[,c(1, grep("Occs", names(Results_short)),
                                    grep("EOO.level", names(Results_short)),
                                    grep("EOO.increase", names(Results_short)))]
  
  if(occ.based) {
    cat("Starting the occurrence-based analysis...", sep= "\n")
    
    list_data <- split(XY.list[[1]], f = XY.list[[1]]$tax)
    
    if (parallel) {
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      
      message('Parallel running with ',
              NbeCores, ' cores')
      
      `%d%` <- foreach::`%dopar%`
    } else {
      `%d%` <- foreach::`%do%`
    }
    
    if(show_progress) {
      pb <- txtProgressBar(min = 0,
                           max = length(list_data),
                           style = 3)
      
      progress <- function(n)
        setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }
    
    x <- NULL
    output <-
      foreach::foreach(
        x = 1:length(list_data),
        .combine = 'c',
        .options.snow = opts
      ) %d% {
        
        if (!parallel & show_progress)
          setTxtProgressBar(pb, x)
        
        source("./R/other_functions.R")
        library(sf)
        library(sp)
        
        res <- my.over.valid.poly(shps_df, list_data[[x]], 
                                  proj_user = proj_user, value = value)
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    XY.list[[1]]$prop.dist.eoo <- output
    Results_long <- dplyr::left_join(XY,
                                     XY.list[[1]][, c("recordID", "classes","prop.dist.eoo")],
                                     by = "recordID")
    Results_long$prop.dist.eoo[is.na(Results_long$prop.dist.eoo) &
                                 Results_long$classes >= max(Results_long$classes, na.rm = TRUE)] <- 0
    Results_long <- Results_long[order(Results_long$recordID), ]
    Results_long <- 
      Results_long[, -which(names(Results_long) %in% c("recordID", "classes"))]
  }
  
  if (write_results)
    write.csv(Results_short, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (occ.based) {
    output <- list(Results_short,
                   Results_long)
    names(output) = c("EOO.change", "Occ.influence")
  } else {
    output <- Results_short
  }
  
  cat("Returning the results.", sep= "\n") 
  output
}

## Version of the first package version
my.EOO.computing <- function(XY,
                             exclude.area = FALSE,
                             country_map = NULL,
                             export_shp = FALSE,
                             write_shp = FALSE,
                             alpha = 1,
                             buff.alpha = 0.1,
                             method.range = "convex.hull",
                             Name_Sp = "species1",
                             buff_width = 0.1,
                             method.less.than3 = "not comp",
                             write_results = TRUE,
                             file.name = "EOO.results",
                             parallel = FALSE,
                             NbeCores = 2,
                             show_progress = TRUE,
                             proj_type = "cea",
                             mode = "spheroid"
) {
  
  
  list_data <- my.coord.check(XY = XY)
  
  
  ### Getting by default land map if poly_borders is not provided
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  } else {
    
    if(any(grepl('sf', class(country_map))))
      country_map <- 
        as(country_map, "Spatial")
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
    country_map <- 
      as(country_map, "sf")
  }
  
  if (buff_width > 80)
    stop("buff_width has unrealistic value")
  
  
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
  
  
  if (is.null(names(list_data))) {
    names_ <-
      rep(Name_Sp, length(list_data))
  } else {
    names_ <- names(list_data)
  }
  
  x <- NULL
  if (show_progress) {
    pb <-
      txtProgressBar(min = 0,
                     max = length(list_data),
                     style = 3)
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }
  
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      source("./R/other_functions.R")
      library(sf)
      library(sp)
      
      res <- my.EOO.comp(
        XY = list_data[[x]],
        exclude.area = exclude.area,
        buff_width = buff_width,
        country_map = country_map,
        Name_Sp = names_[x],
        method.range = method.range,
        alpha = alpha,
        buff.alpha = buff.alpha,
        method.less.than3 = method.less.than3, 
        mode = mode, 
        proj_type = proj_type
      )
      
      names(res)[1] <-
        paste0(names(res)[1], "_" , x)
      if (length(res) > 1)
        names(res)[2] <-
        paste0(names(res)[2], "_" , x)
      
      res
      
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  Results_short <-
    data.frame(EOO = unlist(output[grep("EOO", names(output))]))
  row.names(Results_short) <- names_
  
  if (length(output) == 1)
    names(output) <- Name_Sp
  
  if(write_shp) {
    
    message("Writing EOO shapefiles in shapesIUCN directory")
    
    dir.create(file.path(paste(getwd(), "/shapesIUCN", sep = "")), showWarnings = FALSE)
    output_spatial <- unlist(output[grep("spatial", names(output))])
    output_spatial <- 
      output_spatial[unlist(lapply(output_spatial, function(x) !is.vector(x)))]
    id_spatial <-
      as.numeric(unlist(lapply(strsplit(
        names(output_spatial), "_"
      ), function(x)
        x[[2]])))
    
    
    for (i in 1:length(output_spatial)) {
      
      
      NAME <- names_[id_spatial[i]]
      NAME <- gsub(" ", "_", NAME)
      
      sf::write_sf(as(output_spatial[[i]], "sf"),
                   dsn = "shapesIUCN",
                   layer = paste(NAME, "_EOO_poly", sep = ""),
                   driver = "ESRI Shapefile",
                   overwrite = TRUE)
      
    }
  }
  
  if (write_results)
    write.csv(Results_short, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (!export_shp)
    output <- Results_short
  
  output
}

## Version of the first package version
my.coord.check <- function(XY, listing = TRUE, proj_type = NULL) {
  XY <- as.data.frame(XY)
  
  if (any(is.na(XY[, c(1:2)]))) {
    print(
      paste(
        "Skipping",
        length(which(rowMeans(is.na(
          XY[, 1:2]
        )) > 0)) ,
        "occurrences because of missing coordinates for",
        # if(verbose)
        length(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                 0), 3])), "taxa"
      )
    )
    XY <- XY[which(!is.na(XY[, 1])), ]
    XY <- XY[which(!is.na(XY[, 2])), ]
  }
  
  if (any(XY[, 2] > 180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 1] > 180))
    stop("coordinates are outside of expected range")
  
  if (!is.null(proj_type)) {
    
    XY_proj <-
      sf::sf_project(
        from = sf::st_crs(4326),
        to =
          sf::st_crs(proj_type),
        pts = XY[, c(2, 1)]
      )[, c(2, 1)]
    
    XY[, c(1, 2)] <- 
      XY_proj[, c(1, 2)]
    
  }
  
  
  if (listing) {
    if (ncol(XY) > 2) {
      colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
      XY$tax <- as.character(XY$tax)
      
      if(length(grep("[?]", XY[,3]))>0) XY[,3] <- gsub("[?]", "_", XY[,3])
      if(length(grep("[/]", XY[,3]))>0) XY[,3] <- gsub("[/]", "_", XY[,3])
      
      list_data <- split(XY, f = XY$tax)
    } else{
      colnames(XY)[1:2] <- c("ddlat", "ddlon")
      list_data <- list(XY)
    }
  } else{
    list_data <-
      XY
  }
  
  return(list_data)
}



## Version of the first package version
Convex.Hull.Poly <- function(XY,
                             mode = "spheroid",
                             exclude.area = FALSE,
                             poly_exclude = NULL,
                             proj_type = "cea") {
  
  if (exclude.area & is.null(poly_exclude))
    stop("exclude.area is TRUE but no shape provided")
  
  if(any(grepl('Spatial', class(poly_exclude))))
    poly_exclude <- as(poly_exclude, "sf")
  
  if (mode == "spheroid") {
    hpts <- grDevices::chull(x =  XY[, 1], y = XY[, 2])
    hpts <- c(hpts, hpts[1])
    coord <- matrix(NA, length(hpts), 2)
    POLY <- "POLYGON(("
    for (i in 1:length(hpts)) {
      POLY <- paste(POLY, XY[hpts[i], 1], " ", XY[hpts[i], 2], sep = "")
      if (i != length(hpts))
        POLY <- paste(POLY, ", ", sep = "")
      if (i == length(hpts))
        POLY <- paste(POLY, "))", sep = "")
      
      coord[i, 1] <- XY[hpts[i], 2]
      coord[i, 2] <- XY[hpts[i], 1]
      
    }
    
    p1 <- rgeos::readWKT(POLY) 
    
    ## Gilles, there is a warning here that may be a potential the problem... use sf to make poly?
    sp::proj4string(p1) <-
      sp::CRS("+proj=longlat +datum=WGS84", doCheckCRSArgs = TRUE)
    
    p1 <- suppressWarnings(geosphere::makePoly(p1)) 
    
    ##Gilles, why not define directly a polygon using 'sf'? See code examples commented below:
    #p1_sf <- sf::st_sfc(sf::st_polygon(list(coord[,2:1])))
    #sf::st_crs(p1_sf) <- sf::st_crs(poly_exclude)
    
    if (exclude.area) {
      p1_sf <- as(p1, "sf")
      
      p1 <-
        suppressWarnings(suppressMessages(sf::st_union(
          sf::st_intersection(p1_sf, poly_exclude)
        )))
      
      sf::st_crs(p1) <-
        "+proj=longlat +datum=WGS84"
      
      if(length(p1) == 0) {
        warning("After excluding areas, the convex hull is empty. EOO is NA.")
        
        p1 <- NA 
      } else {
        
        p1 <- 
          as(p1, "Spatial")
        
      }
      
    }
    
  }
  
  if (mode == "planar") {
    
    if(class(proj_type) != "CRS") {
      projEAC <- my.proj_crs(proj_type = proj_type)
    } else {
      projEAC <- proj_type
    }
    
    XY_sf_proj <-
      sf::sf_project(
        from = sf::st_crs(4326),
        to =
          sf::st_crs(projEAC),
        pts = XY[, c(1, 2)]
      )
    
    p1 <- sf::st_convex_hull(x = sf::st_multipoint(XY_sf_proj))
    # eoo <- st_area(p1)
    
    p1 <-
      sf::st_sfc(p1)
    
    sf::st_crs(p1) <- projEAC
    
    if (exclude.area) {
      
      poly_exclude_proj <-
        sf::st_transform(poly_exclude, crs = projEAC)
      
      p1 <-
        sf::st_union(sf::st_intersection(p1, poly_exclude_proj))
      
      if(length(p1) == 0) {
        warning("After excluding areas, the convex hull is empty. EOO is NA.")
        
        p1 <- NA 
      } else {
        
        p1 <- 
          as(p1, "Spatial")
        
      }
      
    }
    
    # p1 <- suppressWarnings(geosphere::makePoly(as(p1, "Spatial")))
    
  }
  
  return(p1)
}


## Version of the first package version
my.proj_crs <- function(proj_type = "cea") {
  
  ## https://epsg.io/54032
  # World Azimuthal Equidistant
  # proj <-
  #   "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  ## https://epsg.io/54002
  # World Azimuthal Equidistant
  # proj <-
  #   "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  
  if(!is.numeric(proj_type)) {
    
    # https://epsg.io/6933
    if(proj_type == "cea")
      proj <-
        "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    
    if(proj_type == "Antarctic")
      proj <-
        "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    if(proj_type == "Africa_eac")
      proj <-
        "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    
  } else {
    
    all_epsg <-
      rgdal::make_EPSG()
    
    proj <- 
      all_epsg[which(all_epsg$code == proj_type), "prj4"]
    
    if(length(proj) == 0)
      stop("No projection found given proj_type")
  }
  
  
  if (rgdal::rgdal_extSoftVersion()[1] >= "3.0.0")  {
    
    # if(proj_type == "cea")
    #   wkt_crs <- 
    #     'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["standard_parallel_1",0.0],UNIT["Meter",1.0]]'
    
    if(!is.numeric(proj_type)) {
      if(proj_type == "cea")
        wkt_crs <- 
          'EPSG:6933'
      
      if(proj_type == "Antarctic")
        wkt_crs <- 
          'EPSG:3031'
      
      if(proj_type == "Africa_eac")
        wkt_crs <- 
          'PROJCS["Africa_Albers_Equal_Area_Conic",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["longitude_of_center",25],PARAMETER["Standard_Parallel_1",20],PARAMETER["Standard_Parallel_2",-23],PARAMETER["latitude_of_center",0],UNIT["Meter",1],AUTHORITY["EPSG","102022"]]'
      
    } else {
      
      wkt_crs <- paste0("EPSG:", proj_type)
      
    }
    
    # wkt_crs <-
    #   rgdal::showWKT(
    #     proj
    #   )
    
    # rgdal::make_EPSG() %>% 
    #   dplyr::as_tibble() %>% 
    #   dplyr::filter(code == 3031)
    # 
    # 
    # rgdal::make_EPSG() %>% 
    #   dplyr::as_tibble() %>% 
    #   dplyr::filter(grepl("+proj=cea", prj4)) %>% 
    #   dplyr::pull(prj4)
    
    crs_proj <- 
      sp::CRS(projargs = proj, 
              SRS_string = wkt_crs, doCheckCRSArgs = TRUE)
    
    # crs_proj <- 
    #   sp::CRS(projargs = proj, 
    #                     SRS_string = wkt_crs, doCheckCRSArgs = TRUE)
  }
  
  if (rgdal::rgdal_extSoftVersion()[1] < "3.0.0")
    crs_proj <-
      sp::CRS(projargs = proj, doCheckCRSArgs = FALSE)
  
  return(crs_proj)
}


## Version of the first package version
alpha.hull.poly <-
  function(XY,
           alpha = 1,
           buff = 0.1,
           exclude.area = FALSE,
           poly_exclude = NULL,
           mode = "spheroid",
           proj_type = "cea")
  {
    
    
    if (mode == "planar") {
      
      if(class(proj_type) != "CRS") {
        projEAC <- my.proj_crs(proj_type = proj_type)
      } else {
        projEAC <- proj_type
      }
      
      XY <-
        sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(projEAC),
          pts = XY[, c(1, 2)]
        )
      
      if (alpha < median(dist(XY,  upper = F))/1000/10) {
        warning("alpha value is defined in kilometer given planar estimation and might be too small to get any polygon")
      }
      
      if (buff < quantile(dist(XY,  upper = F), probs = 0.01)/1000/5) {
        warning("buff value is defined in kilometer given planar estimation and might be too small to get any polygon")
      }
      
      
      buff <-
        buff*1000
      
      alpha <-
        alpha*1000
      
    }
    
    
    Used_data <- unique(XY)
    if (any(rownames(installed.packages()) == "alphahull")) {
      
      loadNamespace("alphahull")
      
      run_alpha <- TRUE
      while(run_alpha) {
        
        ahull.obj <-
          try(alphahull::ahull(Used_data[, c(1, 2)], alpha = alpha), silent = T)
        
        if(class(ahull.obj) != "try-error") {
          
          run_alpha <- FALSE
          
        } else {
          
          Used_data <- apply(Used_data, 2, jitter)
          
        }
      }
      
      y.as.spldf <- ahull_to_SPLDF(ahull.obj)
      y.as.spldf_buff <- rgeos::gBuffer(y.as.spldf, width = buff)
      
      NZp <- slot(y.as.spldf_buff, "polygons")
      holes <-
        lapply(NZp, function(x)
          sapply(slot(x, "Polygons"), slot,
                 "hole"))
      res <- lapply(1:length(NZp), function(i)
        slot(NZp[[i]],
             "Polygons")[!holes[[i]]])
      IDs <- row.names(y.as.spldf_buff)
      NZfill <- SpatialPolygons(
        lapply(1:length(res), function(i)
          Polygons(res[[i]], ID = IDs[i])),
        proj4string = sp::CRS(sp::proj4string(y.as.spldf_buff), doCheckCRSArgs = TRUE)
      )
      
      # crs <- sp::CRS("+proj=longlat +datum=WGS84", doCheckCRSArgs = TRUE)
      # raster::crs(NZfill) <- crs
      
      NZfill <-
        suppressWarnings(rgeos::gBuffer(NZfill, byid = TRUE, width = 0))
      
      if (mode == "planar") {
        
        NZfill <- as(NZfill, "sf")
        
        sf::st_crs(NZfill) <-
          projEAC
        
        # sp::proj4string(NZfill) <- projEAC
        
        if (exclude.area) {
          
          poly_exclude_proj <-
            st_transform(st_make_valid(poly_exclude), crs = projEAC)
          
          
          ### work around to avoid bug appearing randomly
          # cont_try <- TRUE
          # while(cont_try) {
          #   NZfill <-
          #     try(suppressWarnings(suppressMessages(st_union(
          #       st_union(st_intersection(st_make_valid(NZfill), poly_exclude_proj))
          #     ))))
          #
          #   if (class(NZfill) != "try-error")
          #     cont_try <- FALSE
          # }
          
          
          NZfill <-
            st_union(st_intersection(st_make_valid(NZfill), poly_exclude_proj))
          
          sf::st_crs(NZfill) <-
            projEAC
          
          if(length(NZfill) == 0) {
            warning("After excluding areas, the alpha hull is empty. EOO is NA.")
            
            NZfill <- NA
          }
          
          NZfill <- as(NZfill, "Spatial")
          
        }
        
      }
      
      if (mode == "spheroid") {
        
        sp::proj4string(NZfill) <- sp::CRS(SRS_string = 'EPSG:4326')
        
        NZfill_vertices <- 
          try(suppressWarnings(geosphere::makePoly(NZfill)), silent = TRUE)
        
        if(class(NZfill_vertices) == "try-error") {
          NZfill_vertices <- 
            NZfill
          
          warning("Failed to add vertices to polygon, be careful with EOO estimation if large EOO")
        }
        
        if (exclude.area) {
          
          NZfill_sf <- as(NZfill, "sf")
          
          # poly_exclude <-
          #   st_make_valid(poly_exclude)
          
          NZfill <-
            suppressWarnings(suppressMessages(st_union(
              st_intersection(st_make_valid(NZfill_sf), poly_exclude)
            )))
          
          # NZfill <- NZfill_sf
          
          if(length(NZfill) == 0) {
            warning("After excluding areas, the alpha hull is empty. EOO is NA.")
            
            NZfill <- NA
          } else {
            
            sf::st_crs(NZfill) <-
              "+proj=longlat +datum=WGS84"
            
            NZfill <- as(NZfill, "Spatial")
            
          }
          
        }
        
      }
      
      # raster::crs(NZfill) <- "+proj=longlat +datum=WGS84"
      
    } else{
      stop("The package alphahull is required for this procedure, please install it")
    }
    
    return(NZfill)
  }

## Version of the first package version
ahull_to_SPLDF <- function(x) {
  if (class(x) != 'ahull')
    stop('This function only works with `ahull` class objects')
  
  # convert ashape edges to DF
  x.ah.df <- as.data.frame(x$arcs)
  
  # convert each arc to a line segment
  l.list <- list()
  for (i in 1:nrow(x.ah.df))
  {
    # extract row i
    row_i <- x.ah.df[i, ]
    
    # extract elements for arc()
    v <- c(row_i$v.x, row_i$v.y)
    theta <- row_i$theta
    r <- row_i$r
    cc <- c(row_i$c1, row_i$c2)
    # from arc()
    angles <- alphahull::anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- cc[1] + r * cos(seqang)
    y <- cc[2] + r * sin(seqang)
    
    # convert to line segment
    l.list[[i]] <- Line(cbind(x, y))
  }
  
  # promote to Lines class, then to SpatialLines class
  l <- Lines(l.list, ID = 1)
  
  # copy over CRS data from original point data
  l.spl <-
    SpatialLines(list(l), proj4string = CRS(as.character(NA), doCheckCRSArgs=TRUE))
  
  # promote to SpatialLinesDataFrame, required for export to GRASS / OGR
  l.spldf <-
    SpatialLinesDataFrame(l.spl, data = data.frame(id = 1), match.ID =
                            FALSE)
  
  return(l.spldf)
}


## Version of the first package version
my.AOO.computing <- function(XY,
                             Cell_size_AOO = 2,
                             nbe.rep.rast.AOO = 0,
                             parallel = FALSE,
                             NbeCores = 2,
                             show_progress = TRUE,
                             export_shp = FALSE,
                             proj_type = "cea"
) {
  
  proj_type <- my.proj_crs(proj_type = proj_type)
  
  list_data <- 
    my.coord.check(XY = XY, proj_type = proj_type)
  
  if(parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    `%d%` <- foreach::`%dopar%`
  }else{
    `%d%` <- foreach::`%do%`
  }
  
  
  x <- NULL
  if(show_progress) {
    pb <-
      txtProgressBar(min = 0,
                     max = length(list_data),
                     style = 3)
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{opts <- NULL}
  
  
  # print(proj_type)
  
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c', .options.snow = opts
    ) %d% {
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      source("./R/other_functions.R")
      
      res <- my.AOO.estimation(
        coordEAC = list_data[[x]],
        cell_size = Cell_size_AOO,
        nbe_rep = nbe.rep.rast.AOO,
        export_shp = export_shp,
        proj_type = proj_type
      )
      
      if (export_shp)
        names(res) <- c("aoo", "spatial")
      
      res
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  if(!export_shp) {
    
    res <- unlist(output)
    names(res) <- names(list_data)
    
  }
  
  if(export_shp) {
    
    res <- unlist(output[names(output) == "aoo"])
    names(res) <- names(list_data)
    shapes <-  unlist(output[names(output) == "spatial"])
    names(shapes) <- names(list_data)
    
  }
  
  if(!export_shp) 
    return(res)
  if(export_shp) 
    return(list(AOO = res, AOO_poly = shapes))
  
}

## Version of the first package version
my.AOO.estimation <- function(coordEAC,
                              cell_size = 2,
                              nbe_rep = 0,
                              # poly_borders = NULL,
                              export_shp = FALSE,
                              proj_type = proj_type
) {
  
  res <-
    my.cell.occupied(
      nbe_rep = nbe_rep,
      size = cell_size,
      coord = coordEAC[,c(2, 1)],
      export_shp = export_shp,
      proj_type = proj_type
    )
  
  AOO <- res[[2]] * cell_size * cell_size  ### AOO
  
  if (export_shp)
    return(list(AOO = AOO, poly_AOO = res[[1]]))
  
  if (!export_shp)
    return(AOO)
  
}


## Version of the first package version
my.cell.occupied <-
  function(nbe_rep = 0,
           size = 4,
           coord,
           export_shp = TRUE,
           proj_type = NULL)  {
    
    crs_proj <-
      proj_type
    
    Corners <- rbind(c(min(coord[, 1]),
                       max(coord[, 1])),
                     c(min(coord[, 2]),
                       max(coord[, 2])))
    
    
    if (nbe_rep == 0) {
      
      Occupied_cells <- vector(mode = "numeric", length = 4)
      decal <- c(0, 1, 2, 3)
      
      
      for (h in decal) {
        
        ext <-
          raster::extent(
            floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
            floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
              size * 1000,
            floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
              size * 1000,
            floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
              size * 1000
          )
        
        r <-
          raster::raster(ext,
                         resolution = size * 1000,
                         crs = crs_proj)
        
        r2_ <-
          raster::rasterize(coord[, 1:2], r)
        
        OCC <-
          length(which(!is.na(raster::values(r2_))))
        
        Occupied_cells[h + 1] <- OCC
        
        ### If only one occupied cell, stop the production of raster
        if (OCC == 1)
          break
      }
    }
    
    if (nbe_rep > 0) {
      Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
      
      for (h in 1:nbe_rep) {
        rd.1 <- runif(1) * size * 1000
        rd.2 <- runif(1) * size * 1000
        
        ext = raster::extent(
          floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
          floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
          floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
          floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
        )
        r = raster::raster(ext, resolution = size * 1000, crs = crs_proj)
        # r
        r2_ <- raster::rasterize(coord[, 1:2], r)
        OCC <- length(which(!is.na(raster::values(r2_))))
        Occupied_cells[h] <- OCC
        # rd.1.vec <- c(rd.1.vec, rd.1)
        # rd.2.vec <- c(rd.2.vec, rd.2)
        if (OCC == 1)
          break
      }
      
    }
    
    Occupied_cells <- Occupied_cells[Occupied_cells > 0]
    Occupied_cells <- min(Occupied_cells)
    
    if (export_shp)
      r2_pol <-
      raster::rasterToPolygons(
        r2_,
        fun = NULL,
        n = 4,
        na.rm = TRUE,
        digits = 6,
        dissolve = FALSE
      )
    
    if (export_shp)
      return(list(r2_pol, Occupied_cells))
    if (!export_shp)
      return(list(NA, Occupied_cells))
    
  }


## Previous version from the version 1 of the package
my.get.patches <- function(XY, 
                           cell_size = NULL,
                           nbe_rep = 0,
                           AOO = NULL,
                           Resol_sub_pop = NULL,
                           subpop_poly = NULL,
                           dist_isolated = NULL,
                           proj_type = "cea",
                           export_shp = FALSE) {
  
  if (!is.null(AOO) & is.null(cell_size))
    stop("Please provide the size (in km) of the grid cells used for AOO")
  
  proj_type <- 
    proj_crs(proj_type = proj_type)
  
  ## Creating and projecting the points  
  XY_proj <-
    sf::st_as_sf(XY, coords = c("ddlon", "ddlat"))
  sf::st_crs(XY_proj) <- 4326   
  XY_proj <- sf::st_transform(XY_proj, crs = proj_type)
  XY_proj_coord <- sf::st_coordinates(XY_proj)
  XY_proj_coord <- as.data.frame(XY_proj_coord)
  
  
  ## if no AOO is provided, getting the AOO raster for the points
  if (is.null(AOO)) {
    res_aoo <-
      my.AOO.estimation(
        coordEAC = XY_proj_coord[, c(2, 1)],
        cell_size = cell_size,
        nbe_rep = nbe_rep,
        export_shp = TRUE,
        proj_type = proj_type
      )
    
    res_aoo_poly <-
      res_aoo$poly_AOO
    cells <- nrow(res_aoo_poly)
    
  } else {
    cells <- round(AOO / (cell_size * cell_size), 0)
  }  
  ## if no subpopulation sf provided, getting it
  if (is.null(subpop_poly)) {
    
    res_subpop <- 
      my.subpop.estimation(
        XY = XY_proj_coord,
        Resol_sub_pop = Resol_sub_pop,
        proj_type = proj_type,
        export_shp = TRUE
      )
    
    res_subpop_poly <-
      res_subpop$poly_subpop
    
  } else {
    
    res_subpop_poly <- 
      subpop_poly
    
  }
  
  ##Obtaining the distances between subpopulations above the dist_isolated threshold
  dist_btw_poly <- sf::st_distance(res_subpop_poly)
  dist_btw_poly <- matrix(dist_btw_poly,
                          nrow = nrow(dist_btw_poly), ncol = nrow(dist_btw_poly))
  dist_btw_poly[upper.tri(dist_btw_poly, diag = TRUE)] <- NA
  
  mat_above_thres <- dist_btw_poly > dist_isolated * 1000
  
  isolated_subpop1 <- 
    apply(mat_above_thres, 1, FUN = function(x) all(x, na.rm = T))
  isolated_subpop2 <- 
    apply(mat_above_thres, 2, FUN = function(x) all(x, na.rm = T))
  
  isolated_subpop_poly <-
    res_subpop_poly[isolated_subpop1 & isolated_subpop2, ]
  
  if (nrow(isolated_subpop_poly) > 0) {
    df <- data.frame(tax = as.character(unique(XY$tax)),
                     frag = "isolated", stringsAsFactors = FALSE)
    isolated_subpop_poly <-
      sf::st_as_sf(data.frame(sf::st_geometry(isolated_subpop_poly), df))
  }
  
  connected_subpop_poly <- res_subpop_poly[!(isolated_subpop1 & isolated_subpop2),]
  
  if (nrow(connected_subpop_poly) > 0) {
    df <- data.frame(tax = as.character(unique(XY$tax)),
                     frag = "connected")
    connected_subpop_poly <- 
      sf::st_as_sf(data.frame(sf::st_geometry(connected_subpop_poly), df))
  }
  
  
  #### GILLES: WHY THIS STEP IS NECESSARY? CAN'T YOU JUST DO: ####
  fraction_aoo <- 100 * nrow(isolated_subpop_poly)/cells 
  
  if(export_shp) {
    
    polys <- rbind(isolated_subpop_poly, connected_subpop_poly)
    
    return(list(fraction_aoo = fraction_aoo, subpop_poly = polys))
    
  } else {
    names(fraction_aoo) <- as.character(unique(XY$tax))
    return(fraction_aoo)
  }
}

## Previous version from the version 1 of the package
my.subpop.estimation <- function(XY,
                                 Resol_sub_pop, 
                                 export_shp = FALSE,
                                 proj_type) {
  
  if(utils::packageVersion("sf") < '0.9.0')
    stop('A more recent version of the "sf" package is needed, please update this package')
  
  points_sf <- sf::st_as_sf(XY, coords = c(1, 2))
  buff_sf <- sf::st_buffer(points_sf, dist = Resol_sub_pop * 1000)
  buff_sf <- sf::st_union(buff_sf) ### GILLES, WE NEED TO JOIN THE OVERLAPPING CIRCLES, RIGHT? PLEASE CHECK!
  # buff_sf_multipoly <- sf::st_cast(buff_sf, "MULTIPOLYGON")
  buff_sf <- sf::st_cast(buff_sf, "POLYGON")
  SubPopPoly <- sf::st_as_sf(data.frame(buff_sf))
  
  sf::st_crs(SubPopPoly) <- proj_type
  
  NbeSubPop <- nrow(SubPopPoly)
  
  if (export_shp) { ## IF/ELSE ADDED BY RENATO
    #SubPopPoly <- as(SubPopPoly, "Spatial") # not transforming to sp anymore
    OUTPUT <- list(NbeSubPop, SubPopPoly)
    names(OUTPUT) <- c("number_subpop", "poly_subpop")
    
  } else {
    OUTPUT <- NbeSubPop
    names(OUTPUT) <- c("number_subpop")
  }
  
  return(OUTPUT)
}

#' 
#' @title Remove Unwanted Spaces
#' 
#' @param x a character or vector
#'
#' @return the character `x` without trailing or double spaces
#'  
#' @keywords internal
#'
#' @noRd
#' 
.squish <- function (x) {
  x <- gsub("\\s\\s+", " ", as.character(x), perl = TRUE)
  # x <- gsub("  ", " ", x, perl = TRUE)
  x <- gsub("^ | $", "", x, perl = TRUE)
  return(x)
}
#' 
#' @title Build Organism Name
#' 
#' @description Combine diffent table columns with species name information (i.e.
#'   genus, epiteth, infra-epiteth) into a single organism name
#' 
#' @param x the data frame with the taxonomic information
#' @param col.names the name of the columns containing the information to be
#'   combined in the desired order
#'
#' @return a vector with the combined information
#' 
#'  
#' @keywords internal
#' 
#' @noRd
#' 
.build.name <- function(x, col.names = c("Genus_original", "Species_original"))
  {
  
  if (any(!col.names %in% colnames(x)))
    stop("One or more names in 'col.names' were not found in 'x'")
  
  # cols <- names(x)[names(x) %in% col.names]
  cols <- names(x)[match(col.names, names(x), nomatch = 0)]

  if (length(cols) > 1) {
    
    #organismName <- apply(x[, cols], 1, paste, collapse = " ")
    organismName <- do.call(paste, x[, cols])
    organismName[organismName %in% c("NA NA", "NA NA NA",
                                     "NULL NULL", "NULL NULL NULL")] <- NA_character_
    organismName <- gsub(" NA NA$| NULL NULL$", "", organismName, perl = TRUE)
    organismName <- gsub(" NA$| NULL$", "", organismName, perl = TRUE)
    organismName <- .squish(organismName)
    return(organismName)
  } else {
    warning("Less than two columns found; skipping...")
  }    
}
#' 
#' @title Complete SIS Connect Files
#' 
#' @description Combine the output SIS Connect files from THREAT assessments 
#'   with those from the available IUCN assessments
#' 
#' @param x the target data frame with the SIS connect file to be completed
#' @param type the type of SIS connect file 
#' @param sis the data frame with the SIS connect file from the IUCN Red List
#' @param x.spp vector with the full species names (binomials) in `x` for
#'   merging
#' @param iucn.spp vector with the full species names (binomials) in `sis` for
#'   merging
#' @param replace should the different info found in `sis` be merged, replaced
#'   or added to `x`? Defaults to TRUE.
#' 
#' @return the input data frame x with the extra columns resulting from the
#'   merge bewteen input and IUCN SIS connect files
#' 
#'  
#' @keywords internal
#' 
#' @noRd
#' 
.merge_sis_connect <- function(x = NULL, type = NULL,
                               sis = NULL, 
                               x.spp = NULL,
                               tax.table = NULL, 
                               replace = TRUE) {
  
  # x <- allfields
  # type <- "allfields"
  # sis <- prev.assess
  # x.spp <- prev.assess$species.correct2
  # tax.table <- tax
  # replace <- TRUE
  
  # Format standardizations
  sis[] <- lapply(sis, gsub, pattern = "&amp;", replacement = "&", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "\t$", replacement = "", perl = TRUE)
  sis[] <- lapply(sis, gsub, pattern = " \n$", replacement = "", perl = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "<em>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "</em>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "<i>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "</i>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, gsub, pattern = "<i>", replacement = "", fixed = TRUE)
  sis[] <- lapply(sis, .squish)
  
  if (type == "taxonomy") {
    x$scientificName <- x.spp
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("kingdomName", "phylumName", "className", "orderName", 
                      "familyName", "genusName", "speciesName", "infraType", 
                      "infraName", "infraAuthority", "authority","taxonomicNotes")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("kingdom", "phylum", "classname", "ordername", "family", 
                   "genus", "species", "infraType", "infra_name", 
                   "infra_authority", "taxonomicAuthority", 
                   "TaxonomicNotes.value")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
                               by = iucn.spp.col)
    rep.col <- paste0("merged.", type)
    x.temp[[rep.col]] <- FALSE
      
    if (replace) {
      
      vector.true <- rep(TRUE, length(nomes))
      names(vector.true) <- nomes
      no.rep <- c("infraType", "infra_authority", "infra_name")
      vector.true[names(vector.true) %in% no.rep] <- FALSE
      rep.tax <- nomes[vector.true]

      vector.action <- rep("replace", length(nomes))
      names(vector.action) <- nomes
      concatenar <- c("TaxonomicNotes.value")
      vector.action[names(vector.action) %in% concatenar] <- "concatenate"
      ignore <- c("taxonomicAuthority")
      vector.action[names(vector.action) %in% ignore] <- ""
      
      for (i in seq_along(rep.tax)) {
        col.i <- rep.tax[i]
        col.rep.i <- 
          colnames(x.temp)[which(grepl(col.i, colnames(x.temp), perl = TRUE))]
        flag.i <- !is.na(x.temp[[col.rep.i[2]]]) &
                    tolower(gsub(" ", "", x.temp[[col.rep.i[1]]], perl = TRUE)) != tolower(gsub(" ", "", x.temp[[col.rep.i[2]]], perl = TRUE))
        if (any(flag.i)) {
          vector.action.i <- vector.action[names(vector.action) %in% col.i]
          new.col.i <- paste0("merged.", col.i)
          
          if (vector.action.i == "replace") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            x.temp[flag.i, new.col.i] <- x.temp[flag.i, col.rep.i[2]]
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
          
          if (vector.action.i == "concatenate") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            
            f1 <- function(x) paste(unique(x), collapse = "|")
            x.temp[flag.i, new.col.i] <- 
              apply(x.temp[flag.i, col.rep.i], 1, f1)
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
        }
      }  
    }
    
    new.columns <- colnames(x.temp)[!colnames(x.temp) %in% colnames(x)]
    x1 <- cbind.data.frame(x, x.temp[ , new.columns])
    
    x2 <- x1[order(as.double(gsub("sp", "", x1$internal_taxon_id))), ]
    
    return(x2)
  }

  if (type == "commonnames") { 
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "language")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("name", "language")
    names(sis.prev) <- c("validName", nomes)
    sis.prev$combo <- 
      paste(sis.prev[["validName"]], sis.prev[[colunas[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[["validName"]], x[[colunas[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- "validName"
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "validName")
    x1 <- x[ , c("validName", nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    # sis.prev.miss2 <- sis.prev.miss1[,match(names(x), names(sis.prev.miss1), nomatch = 0)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
  }
  
  if (type == "countries") {
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "code", "presence", "origin", "seasonality", "formerlyBred")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("CountryOccurrence.CountryOccurrenceSubfield.", colunas)
    nomes[1:2] <- c("CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName", 
                    "CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceLookup")
    names(sis.prev) <- c("species.correct2", nomes)
    
    check_these <- !sis.prev[[nomes[1]]] %in% x$CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName
    diff <- sis.prev[[nomes[1]]][check_these]
    sort(table(diff))
    diff1 <- plantR::prepCountry(countrycode::countrycode(diff, "country.name", "iso3c"), 
                                 to.lower = FALSE, special.char = FALSE, rm.abbrev = FALSE)
    diff1 <- gsub("u\\.s\\. virgin islands", "United States Virgin Islands", 
                  diff1, perl = TRUE, ignore.case = TRUE)
    diff1 <- gsub("St. Barthelemy", "St. Barthélemy", 
                  diff1, perl = TRUE, ignore.case = TRUE)
    sis.prev[[nomes[1]]][check_these] <- diff1

    sis.prev$combo <- 
      paste(sis.prev[["species.correct2"]], sis.prev[[nomes[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[["species.correct2"]], x[[nomes[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(as.data.frame(x), sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    rep_these <- !is.na(tmp1[[paste0(nomes[4], ".y")]]) &
                  !tmp1[[paste0(nomes[4], ".y")]] %in% "Native" &
                  tmp1[[nomes[4]]] %in% "Native" 
    x[rep_these, nomes[c(4,6)]] <- tmp1[rep_these, paste0(nomes[c(4,6)], ".y")]
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "species.correct2")
    
    x1 <- x[ , c("species.correct2", nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    ## Taging countries with uncertain occurrences (probable undeclared cultivated individuals)
    # not doing it (need to know for which new countries may be uncertain e.g. Araucaria and which are not e.g. Cabralea canjerana)
    # for (i in seq_along(unique(x1$internal_taxon_id))) {
    #   id.i <- unique(x1$species.correct2)[i]
    #   x.i <- x1[x1$species.correct2 %in% id.i, ]
    #   if (id.i %in% sis.prev$species.correct2) {
    #     iucn.i <- sis.prev[sis.prev$species.correct2 %in% id.i, ]
    #     new.countries.i <- unique(x.i[[2]][!x.i[[2]] %in% iucn.i[[2]]])
    #     
    #     if (length(new.countries.i) > 0)
    #       x1[x1$species.correct2 %in% id.i &
    #            x1$CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName %in% new.countries.i, 
    #          "CountryOccurrence.CountryOccurrenceSubfield.origin"] <- ""
    #   }
    # }
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
    
  }
  
  if (type == "habitats") {
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "code", "majorImportance", "season", "suitability")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName", 
               "GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup",
               "GeneralHabitats.GeneralHabitatsSubfield.majorImportance",
               "GeneralHabitats.GeneralHabitatsSubfield.season",
               "GeneralHabitats.GeneralHabitatsSubfield.suitability")
    names(sis.prev) <- c("validName", nomes)
    
    sis.prev <- sis.prev[!sis.prev[[nomes[2]]] %in% 
                           c(paste0("9.", 1:10), 
                             paste0("10.", 1:4),
                             paste0("11.", 1:6),
                             paste0("12.", 1:7),
                             paste0("13.", 1:5),
                             paste0("15.", 1:13)),]
    
    # diferenca <- sis.prev[[nomes[2]]][!sis.prev[[nomes[2]]] %in%
    #                                             unique(x[[nomes[2]]])]
    # sort(table(diferenca))

    sis.prev$combo <- 
      paste(sis.prev[["validName"]], sis.prev[[nomes[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    x[["validName"]] <- x.spp
    combo.x <- paste(x[["validName"]], x[[nomes[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- gsub(" – ", " ", combo.x, perl = TRUE)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
      
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    rep_these <- !is.na(tmp1[[paste0(nomes[3], ".y")]]) &
                  tmp1[[nomes[3]]] %in% "" 
    x[rep_these, nomes[c(3,5)]] <- tmp1[rep_these, paste0(nomes[c(3,5)], ".y")]

    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- "validName"
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "validName")
    x1 <- x[ , c("validName", nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    x$combo <- NULL
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss1$combo <- NULL
    sis.prev.miss2 <- sis.prev.miss1[, names(x1)]
  
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    x3$validName <- NULL
    
    return(x3)
  }
  
  if (type == "plant_specific") {
    
    x$scientificName <- x.spp
    
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "code")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName", 
               "PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    # x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
    #                            by = iucn.spp.col)
    # rep.col <- paste0("merged.", type)
    # x.temp[[rep.col]] <- FALSE
    
    # Removing vines, epiphytes and herbs
    sis.prev <- sis.prev[grepl(" - ",iucn.sis.plant.spec$name),]
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    x[[iucn.spp.col]] <- x.spp
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- gsub(" – ", " ", combo.x, perl = TRUE)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    sis.prev.miss2 <- sis.prev.miss1[,match(names(x), names(sis.prev.miss1), nomatch = 0)]
    x1 <- x[ , c(iucn.spp.col, nomes, "internal_taxon_id")]
    x1$source.info <- "Brazilian_Flora_Online_v.393.285"
    sis.prev.miss2$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2$combo <- NULL
    sis.prev.miss3 <- sis.prev.miss2[,match(names(x1), names(sis.prev.miss2), nomatch = 0)]
    
    x2 <- rbind.data.frame(x1, sis.prev.miss3)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    x3[[iucn.spp.col]] <- NULL
    
    ## Removing duplicated info for the same species with lower resolution (i.e. "T" and "S")
    dup.IDs <- x3$internal_taxon_id[duplicated(x3$internal_taxon_id)]
    remove_these <- rep(TRUE, dim(x3)[1])

    for (i in seq_along(unique(dup.IDs))) {
      id.i <- unique(dup.IDs)[i]
      check_these <- x3$internal_taxon_id %in% id.i
      code.i <- 
        x3$PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup[check_these]
      rep_these <- nchar(code.i) == 1 
      if (any(rep_these) & !all(rep_these)) {
        remove_these[check_these][rep_these] <- FALSE
      }
    }
    
    x4 <- x3[remove_these, ]
    return(x4)
  }
  
  if (type == "synonyms") {
   
    iucn.spp.col <- "scientificName"
    
    colunas <- c("name", "speciesName", "speciesAuthor", "infraType", 
                 "infrarankAuthor")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- c("name", "speciesName", "speciesAuthor", "infraType", 
               "infrarankAuthor")
    names(sis.prev) <- c("species.correct2" , nomes)
    sis.prev$combo <- 
      paste(sis.prev[["species.correct2"]], sis.prev[[colunas[1]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE))
    sis.prev$combo <- .squish(tolower(gsub(" [orth. error]", "", sis.prev$combo, perl = TRUE)))
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    x1 <- dplyr::left_join(x, tax.table1, by = "internal_taxon_id")
    
    combo.x <- paste(x1[["species.correct2"]], x[[colunas[1]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    sis.prev.miss$combo <- NULL
    
    x2 <- x1[ , c("species.correct2", nomes, "infrarankName", "internal_taxon_id")]
    
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = "species.correct2")
    
    sis.prev.miss1$infrarankName <- NA
    sis.prev.miss2 <- sis.prev.miss1[,match(names(x2), names(sis.prev.miss1), nomatch = 0)]

    x2$source.info <- "Brazilian_Flora_Online_v.393.285"
    sis.prev.miss2$source.info <- "IUCN_Red_List_v.2022-2"
    x3 <- rbind.data.frame(x2, sis.prev.miss2)
    
    return(x3)
  }

  if (type == "conservation_needed") {
    
    iucn.spp.col <- "scientificName"
    
    x[[iucn.spp.col]] <- x.spp
    x$ConservationActions.ConservationActionsSubfield.note.iucn_rl <- ""
    
    colunas <- c("name", "code", "note")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("ConservationActions.ConservationActionsSubfield.", colunas)
    nomes[1:3] <- c("ConservationActions.ConservationActionsSubfield.ConservationActionsName", 
                    "ConservationActions.ConservationActionsSubfield.ConservationActionsLookup",
                    "ConservationActions.ConservationActionsSubfield.note.iucn_rl")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    sis.prev$ConservationActions.ConservationActionsSubfield.note <- ""
    
    sis.prev <-
      sis.prev[!sis.prev[[nomes[2]]] %in% c("4.1"), ]
    # sis.prev[[nomes[1]]] <- gsub("", "", sis.prev[[nomes[1]]], perl = TRUE)
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x[[nomes[2]]]
    # diff <- sis.prev[[nomes[2]]][check_these]
    # table(diff)
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))

    # Adding previous notes (not doing it - too difficult, using separate columns)
    # rep_these <- !is.na(tmp1[paste0(nomes[3], ".y")])
    # if (any(rep_these))
    #   x[rep_these, nomes[3]] <- tmp1[rep_these, paste0(nomes[3], ".y")]
    
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% x$combo, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, 
                 "ConservationActions.ConservationActionsSubfield.note", "internal_taxon_id")]
    x1$source.info <- "THREAT_project"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    x2$remove_these <- TRUE
    
    ## Replacing the categories Scale Unknown by categories with known threat scale
    # not doing it: ask CNCFlora what they prefer
    unknown <- c("3.4")
    for (i in seq_along(unique(x1$internal_taxon_id))) {
      id.i <- unique(x1$scientificName)[i]
      x.i <- x1[x1$scientificName %in% id.i, ]
      if (id.i %in% sis.prev$scientificName) {
        iucn.cats.i <- sis.prev[sis.prev$scientificName %in% id.i, ]
        new.cats.i <- unique(x.i[[3]][!x.i[[3]] %in% iucn.cats.i[[3]]])

        unknown.cats.i <- new.cats.i[new.cats.i %in% unknown]
        if (length(unknown.cats.i) > 0) {
          for (j in seq_along(unknown.cats.i)) {
            if (any(grepl(unknown.cats.i[j], iucn.cats.i[[3]])))
              x2[x2$scientificName %in% id.i &
                   x2[[3]] %in% unknown.cats.i[j],
                 "remove_these"] <- FALSE
          }
        }
      }
    }
    x3 <- x2[x2$remove_these, ]
    x3$remove_these <- NULL

    x3 <- x3[order(as.double(gsub("sp", "", x3$internal_taxon_id))), ]
    x3 <- x3[ , c("ConservationActions.ConservationActionsSubfield.ConservationActionsLookup",
                  "ConservationActions.ConservationActionsSubfield.ConservationActionsName",
                  "ConservationActions.ConservationActionsSubfield.note", 
                  "internal_taxon_id",
                  "ConservationActions.ConservationActionsSubfield.note.iucn_rl")]
    return(x3)
  }
  
  if (type == "research_needed") {
    
    iucn.spp.col <- "scientificName"
    
    x[[iucn.spp.col]] <- x.spp
    x$Research.ResearchSubfield.note.iucn_rl <- ""
    
    colunas <- c("name", "code", "note")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("Research.ResearchSubfield.", colunas)
    nomes[1:3] <- c("Research.ResearchSubfield.ResearchName", 
                    "Research.ResearchSubfield.ResearchLookup",
                    "Research.ResearchSubfield.note.iucn_rl")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    sis.prev$Research.ResearchSubfield.note <- ""
    
    sis.prev <-
      sis.prev[!sis.prev[[nomes[2]]] %in% c("4"), ]
    # sis.prev[[nomes[1]]] <- gsub("", "", sis.prev[[nomes[1]]], perl = TRUE)
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x[[nomes[2]]]
    # diff <- sis.prev[[nomes[2]]][check_these]
    # table(diff)
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    # sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    # tmp1 <- dplyr::left_join(x, sis.prev.common, 
    #                          by = "combo", suffix = c("", ".y"))
    # 
    # Adding previous notes (not doing it - too difficult, using separate columns)
    # rep_these <- !is.na(tmp1[paste0(nomes[3], ".y")]) &
    #                 !tmp1[paste0(nomes[3], ".y")] %in% "Population size, distribution & trends"
    # if (any(rep_these))
    #   x[rep_these, nomes[3]] <- tmp1[rep_these, paste0(nomes[3], ".y")]
    
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% x$combo, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, 
                 "Research.ResearchSubfield.note", "internal_taxon_id")]
    x1$source.info <- "THREAT_project"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    x3 <- x3[ , c("Research.ResearchSubfield.ResearchLookup",
                  "Research.ResearchSubfield.ResearchName",
                  "Research.ResearchSubfield.note", 
                  "internal_taxon_id",
                  "Research.ResearchSubfield.note.iucn_rl")]
    return(x3)
    
    
  }

  if (type == "threats") {
    
    iucn.spp.col <- "scientificName"
    
    x[[iucn.spp.col]] <- x.spp
    
    colunas <- c("name", "code", "stressCode", "stressName", "text",
                 "ancestry", "ias", "internationalTrade",
                 "scope", "severity", "timing", "virus")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("Threats.ThreatsSubfield.", colunas)
    nomes[1:2] <- c("Threats.ThreatsSubfield.ThreatsName", 
                    "Threats.ThreatsSubfield.ThreatsLookup")
    nomes[3:4] <- c("Threats.ThreatsSubfield.StressesSubfield.stress", 
                    "Threats.ThreatsSubfield.StressesSubfield.stressdesc")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    
    sis.prev <-
      sis.prev[!sis.prev$Threats.ThreatsSubfield.ThreatsLookup %in% c("1", "12.1"), ]
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x$Threats.ThreatsSubfield.ThreatsLookup
    # diff <- sis.prev[[nomes[2]]][check_these]
    # table(diff)
    
    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    tmp1[[paste0(nomes[4], ".y")]] <- 
      stringr::str_to_title(tmp1[[paste0(nomes[4], ".y")]])
    tmp1[[paste0(nomes[4], ".y")]] <- 
      gsub(" Of ", " of ",tmp1[[paste0(nomes[4], ".y")]], perl = TRUE)
    
    
    # Adding missing stress codes
    rep_these <- !is.na(tmp1[[paste0(nomes[3], ".y")]])
    if (any(rep_these)) {
      prev.i <- strsplit(tmp1[[paste0(nomes[3], ".y")]][rep_these], "|", fixed = TRUE)
      prev.names.i <- strsplit(tmp1[[paste0(nomes[4], ".y")]][rep_these], "|", fixed = TRUE)
      prev.ii <- mapply(function(x, y) { 
                        names(x) <- y 
                        return(x) }, 
                        prev.i, prev.names.i)
      new.i <- strsplit(tmp1[[nomes[3]]][rep_these], "|", fixed = TRUE)  
      new.names.i <- strsplit(tmp1[[nomes[4]]][rep_these], "|", fixed = TRUE)
      new.ii <- mapply(function(x, y) { 
        names(x) <- y 
        return(x) }, 
        new.i, new.names.i)

      combo <- mapply(function(x, y) { 
                        vector <- c(unlist(x), unlist(y))
                        vector1 <- vector[!duplicated(tolower(names(vector)))]
                        ordem <- order(as.numeric(gsub("\\.3\\.",".3",vector1, perl = TRUE)))
                        vector2 <- vector1[ordem]
                        return(vector2)
                        }, new.ii, prev.ii)
      new.iii <- sapply(combo, paste0, collapse = "|")
      tmp1[[nomes[3]]][rep_these] <- new.iii
      
      new.names.iii <- sapply(combo, function(x) paste0(names(x), collapse = "|"))
      tmp1[[nomes[4]]][rep_these] <- new.names.iii
    }

    for (i in c(5:12)) {
      rep_these <- !is.na(tmp1[[paste0(nomes[i], ".y")]]) &
                    !tmp1[[paste0(nomes[i], ".y")]] %in% c("") &
                    tmp1[[nomes[i]]] %in% c("", "Unknown")
      if (any(rep_these))
        x[rep_these, nomes[i]] <- tmp1[rep_these, paste0(nomes[i], ".y")]
    }
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% x$combo, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, "internal_taxon_id")]
    x1$source.info <- "THREAT_project"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    # x2$remove_these <- rep(TRUE, dim(x2)[1])
    # 
    ## Replacing the categories Scale Unknown by categories with known threat scale
    # not doing it: ask CNCFlora what they prefer
    # unknown <- c("2.1.4", "2.3.4", "5.2.4", "5.3.5", "5.4.6")
    # for (i in seq_along(unique(x1$internal_taxon_id))) {
    #   id.i <- unique(x1$scientificName)[i]
    #   x.i <- x1[x1$scientificName %in% id.i, ]
    #   if (id.i %in% sis.prev$scientificName) {
    #     iucn.cats.i <- sis.prev[sis.prev$scientificName %in% id.i, ]
    #     new.cats.i <- unique(x.i[[3]][!x.i[[3]] %in% iucn.cats.i[[3]]])
    #     
    #     unknown.cats.i <- new.cats.i[new.cats.i %in% unknown]
    #     if (length(unknown.cats.i) > 0) {
    #       for (j in seq_along(unknown.cats.i)) {
    #         if (any(grepl(substr(unknown.cats.i[j],1,4), iucn.cats.i[[3]])))
    #           x2[x2$scientificName %in% id.i &
    #                x2[[3]] %in% unknown.cats.i[j],
    #              "remove_these"] <- FALSE
    #       }
    #     }
    #   }
    # }
    # x3 <- x2[x2$remove_these, ]
    # x3$remove_these <- NULL
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
  }
  
  if (type == "usetrade") {
    
    iucn.spp.col <- "scientificName"
   
    x[[iucn.spp.col]] <- x.spp

    colunas <- c("name", "code", "international", "national", "other", "subsistence")
    sis.prev <- sis[, c(iucn.spp.col, colunas)]
    nomes <- paste0("UTEndUse.UTEndUseSubfield.", colunas)
    nomes[1:2] <- c("UTEndUse.UTEndUseSubfield.UTEndUseName", 
                    "UTEndUse.UTEndUseSubfield.UTEndUseLookup")
    names(sis.prev) <- c(iucn.spp.col, nomes)
    
    sis.prev <-
      sis.prev[!sis.prev$UTEndUse.UTEndUseSubfield.UTEndUseLookup %in% c("15", "18"), ]
    
    # check_these <- !sis.prev[[nomes[2]]] %in% x$UTEndUse.UTEndUseSubfield.UTEndUseLookup
    # diff <- sis.prev[[nomes[2]]][check_these]
    # sort(table(diff))

    sis.prev$combo <- 
      paste(sis.prev[[iucn.spp.col]], sis.prev[[nomes[2]]], sep = "_")
    sis.prev$combo <- plantR::rmLatin(sis.prev$combo)
    sis.prev$combo <- .squish(tolower(gsub("[[:punct:]]", " ", sis.prev$combo, perl = TRUE)))
    
    combo.x <- paste(x[[iucn.spp.col]], x[[nomes[2]]], sep = "_")
    combo.x <- plantR::rmLatin(combo.x)
    combo.x <- .squish(tolower(gsub("[[:punct:]]", " ", combo.x, perl = TRUE)))
    x[["combo"]] <- combo.x
    
    sis.prev.common <- sis.prev[sis.prev$combo %in% x$combo, ]
    tmp1 <- dplyr::left_join(x, sis.prev.common, 
                             by = "combo", suffix = c("", ".y"))
    
    for (i in 3:6) {
      rep_these <- !is.na(tmp1[[paste0(nomes[i], ".y")]]) &
                    tmp1[[paste0(nomes[i], ".y")]] %in% "true" &
                    tmp1[[nomes[i]]] %in% "false" 
      x[rep_these, nomes[i]] <- tmp1[rep_these, paste0(nomes[i], ".y")]
    }
    
    sis.prev.miss <- sis.prev[!sis.prev$combo %in% combo.x, ]
    sis.prev.miss$combo <- NULL
    x$combo <- NULL
    
    tax.table1 <- tax.table[, c("species.correct2", "internal_taxon_id")]
    names(tax.table1)[1] <- iucn.spp.col
    sis.prev.miss1 <- dplyr::left_join(sis.prev.miss, 
                                       tax.table1, 
                                       by = iucn.spp.col)
    
    x1 <- x[ , c(iucn.spp.col, nomes, "internal_taxon_id")]
    x1$source.info <- "TreeCo_database_v.4.1"
    
    sis.prev.miss1$source.info <- "IUCN_Red_List_v.2022-2"
    sis.prev.miss2 <- sis.prev.miss1[,names(x1)]
    x2 <- rbind.data.frame(x1, sis.prev.miss2)
    
    x3 <- x2[order(as.double(gsub("sp", "", x2$internal_taxon_id))), ]
    
    return(x3)  
  }

  if (type == "assessments") {
    
    iucn.spp.col <- "species.correct2"

    x[[iucn.spp.col]] <- x.spp
    
    colunas <- c("realm", 
                 "conservationActions", 
                 "habitat", 
                 #"MapStatus", 
                 "population",
                 "populationTrend",
                 "range",
                 "yearPublished",
                 "criteriaVersion",
                 #"RedListCriteria.dataDeficientReason",
                 #"RedListCriteria.isManual",
                 "redlistCategory",
                 "redlistCriteria",
                 #"RedListCriteria.manualCriteriaString",
                 "possiblyExtinct",
                 #"RedListCriteria.possiblyExtinctCandidate",
                 "yearLastSeen",
                 "rationale",
                 #"RedListReasonsForChange.catCritChanges",
                 #"RedListReasonsForChange.changeReasons",
                 #"RedListReasonsForChange.otherReason",
                 #"RedListReasonsForChange.timeframe",
                 #"RedListReasonsForChange.type",
                 "systems",
                 "threats",
                 "useTrade",
                 "language")
    sis.prev <- sis[, c("internal_taxon_id", iucn.spp.col, colunas)]
    
    nomes <- c("BiogeographicRealm.realm", 
                 "ConservationActionsDocumentation.narrative", 
                 "HabitatDocumentation.narrative", 
                 #"MapStatus.status", 
                 "PopulationDocumentation.narrative",
                 "PopulationTrend.value",
                 "RangeDocumentation.narrative",
                 "RedListAssessmentDate.value",
                 "RedListCriteria.critVersion",
                 #"RedListCriteria.dataDeficientReason",
                 #"RedListCriteria.isManual",
                 "RedListCriteria.manualCategory",
                 "RedListCriteria.manualCriteria",
                 #"RedListCriteria.manualCriteriaString",
                 "RedListCriteria.possiblyExtinct",
                 #"RedListCriteria.possiblyExtinctCandidate",
                 "RedListCriteria.yearLastSeen",
                 "RedListRationale.value",
                 #"RedListReasonsForChange.catCritChanges",
                 #"RedListReasonsForChange.changeReasons",
                 #"RedListReasonsForChange.otherReason",
                 #"RedListReasonsForChange.timeframe",
                 #"RedListReasonsForChange.type",
                 "System.value",
                 "ThreatsDocumentation.value",
                 "UseTradeDocumentation.value",
                 "language.value")
    
    names(sis.prev) <- c("internal_taxon_id", iucn.spp.col, nomes)

    x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
                               by = c("internal_taxon_id", iucn.spp.col))
    rep.col <- paste0("merged.", type)
    x.temp[[rep.col]] <- FALSE
    
    if (replace) {
      
      vector.true <- rep(TRUE, length(nomes))
      names(vector.true) <- nomes
      no.rep <- c("RedListAssessmentDate.value", "System.value", "language.value")
      vector.true[names(vector.true) %in% no.rep] <- FALSE
      rep.tax <- nomes[vector.true]
      
      vector.action <- rep("replace", length(nomes))
      names(vector.action) <- nomes
      concatenar <- c("BiogeographicRealm.realm", 
                      "ConservationActionsDocumentation.narrative",
                      "HabitatDocumentation.narrative",
                      "PopulationDocumentation.narrative",
                      "RangeDocumentation.narrative",
                      "RedListRationale.value",
                      "ThreatsDocumentation.value",
                      "UseTradeDocumentation.value")
      vector.action[names(vector.action) %in% concatenar] <- "concatenate"
      ignore <- c("PopulationTrend.value", "RedListCriteria.critVersion", 
                  "RedListCriteria.manualCategory", 
                  "RedListCriteria.manualCriteria", 
                  "RedListCriteria.possiblyExtinct",
                  "RedListCriteria.yearLastSeen")
      vector.action[names(vector.action) %in% ignore] <- ""
      
      for (i in seq_along(rep.tax)) {
        col.i <- rep.tax[i]
        col.rep.i <- 
          colnames(x.temp)[which(grepl(col.i, colnames(x.temp), perl = TRUE))]
        flag.i <- !is.na(x.temp[[col.rep.i[2]]]) &
          tolower(gsub(" ", "", x.temp[[col.rep.i[1]]], perl = TRUE)) != tolower(gsub(" ", "", x.temp[[col.rep.i[2]]], perl = TRUE))
        if (any(flag.i)) {
          vector.action.i <- vector.action[names(vector.action) %in% col.i]
          new.col.i <- paste0("merged.", col.i)
          
          if (vector.action.i == "replace") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            x.temp[flag.i, new.col.i] <- x.temp[flag.i, col.rep.i[2]]
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
          
          if (vector.action.i == "concatenate") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            
            f1 <- function(x) paste(unique(x), collapse = "|")
            x.temp[flag.i, new.col.i] <- 
              apply(x.temp[flag.i, col.rep.i], 1, f1)
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
        }
      }  
    }
    
    new.columns <- colnames(x.temp)[!colnames(x.temp) %in% colnames(x)]
    x1 <- cbind.data.frame(x, x.temp[ , new.columns])
    
    x2 <- x1[order(as.double(gsub("sp", "", x1$internal_taxon_id))), ]
    
    return(x2)
    
  }

  if (type == "allfields") {
    
    iucn.spp.col <- "species.correct2"
    
    x[[iucn.spp.col]] <- x.spp
    
    colunas <- c("AOO.range",
                 "EOO.range",
                 "NoThreats.noThreats",
                 "ElevationLower.limit",
                 "ElevationUpper.limit",
                 "PopulationSize.range",
                 "ThreatsUnknown.value",
                 "LocationsNumber.range",
                 "GenerationLength.range",
                 "MovementPatterns.pattern",
                 "SubpopulationNumber.range",
                 "AreaRestricted.isRestricted",
                 "CropWildRelative.isRelative",
                 "YearOfPopulationEstimate.value",
                 "InPlaceEducationControlled.value",
                 "SevereFragmentation.isFragmented",
                 "InPlaceResearchRecoveryPlan.value",
                 "InPlaceLandWaterProtectionInPA.value",
                 "InPlaceSpeciesManagementExSitu.value",
                 "InPlaceResearchMonitoringScheme.value",
                 "InPlaceEducationSubjectToPrograms.value",
                 "InPlaceSpeciesManagementHarvestPlan.value",
                 "InPlaceLandWaterProtectionAreaPlanned.value",
                 "InPlaceEducationInternationalLegislation.value",
                 "InPlaceLandWaterProtectionInvasiveControl.value",
                 "InPlaceLandWaterProtectionSitesIdentified.value",
                 "InPlaceLandWaterProtectionPercentProtected.value")
    sis.prev <- sis[, c("internal_taxon_id", iucn.spp.col, colunas)]
    nomes <- colunas
    
    x.temp <- dplyr::left_join(x, sis.prev, suffix = c("", ".iucn_rl"), 
                               by = c("internal_taxon_id", iucn.spp.col))
    rep.col <- paste0("merged.", type)
    x.temp[[rep.col]] <- FALSE
    
    if (replace) {
      
      vector.true <- rep(TRUE, length(nomes))
      names(vector.true) <- nomes
      no.rep <- c("MovementPatterns.pattern")
      vector.true[names(vector.true) %in% no.rep] <- FALSE
      rep.tax <- nomes[vector.true]
      
      vector.action <- rep("replace", length(nomes))
      names(vector.action) <- nomes
      concatenar <- c("NoThreats.noThreats", 
                      "ElevationLower.limit",
                      "ElevationUpper.limit",
                      "PopulationSize.range",
                      "ThreatsUnknown.value",
                      "RedListRationale.value",
                      "ThreatsDocumentation.value",
                      "UseTradeDocumentation.value",
                      "AreaRestricted.isRestricted", 
                      "CropWildRelative.isRelative", 
                      "SevereFragmentation.isFragmented")
      vector.action[names(vector.action) %in% concatenar] <- "concatenate"
      ignore <- c("AOO.range", "EOO.range", "LocationsNumber.range",
                  "SubpopulationNumber.range", "GenerationLength.range",
                  "SubpopulationNumber.range", "YearOfPopulationEstimate.value",
                  "InPlaceEducationControlled.value", 
                  "InPlaceResearchRecoveryPlan.value",
                  "InPlaceLandWaterProtectionInPA.value",
                  "InPlaceSpeciesManagementExSitu.value",
                  "InPlaceResearchMonitoringScheme.value",
                  "InPlaceEducationSubjectToPrograms.value",
                  "InPlaceSpeciesManagementHarvestPlan.value",
                  "InPlaceLandWaterProtectionAreaPlanned.value",
                  "InPlaceEducationInternationalLegislation.value",
                  "InPlaceLandWaterProtectionInvasiveControl.value",
                  "InPlaceLandWaterProtectionSitesIdentified.value",
                  "InPlaceLandWaterProtectionPercentProtected.value")
      vector.action[names(vector.action) %in% ignore] <- ""
      
      for (i in seq_along(rep.tax)) {
        col.i <- rep.tax[i]
        col.rep.i <- 
          colnames(x.temp)[which(grepl(col.i, colnames(x.temp), perl = TRUE))]
        flag.i <- !is.na(x.temp[[col.rep.i[2]]]) &
          tolower(gsub(" ", "", x.temp[[col.rep.i[1]]], perl = TRUE)) != tolower(gsub(" ", "", x.temp[[col.rep.i[2]]], perl = TRUE))
        if (any(flag.i)) {
          vector.action.i <- vector.action[names(vector.action) %in% col.i]
          new.col.i <- paste0("merged.", col.i)
          
          if (vector.action.i == "replace") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            x.temp[flag.i, new.col.i] <- x.temp[flag.i, col.rep.i[2]]
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
          
          if (vector.action.i == "concatenate") {
            x.temp[[new.col.i]] <- x.temp[[col.rep.i[1]]]
            
            f1 <- function(x) paste(unique(x), collapse = "|")
            x.temp[flag.i, new.col.i] <- 
              apply(x.temp[flag.i, col.rep.i], 1, f1)
            
            x.temp[[rep.col]][flag.i] <- TRUE
          }
        }
      }  
    }
    
    new.columns <- colnames(x.temp)[!colnames(x.temp) %in% colnames(x)]
    x1 <- cbind.data.frame(x, x.temp[ , new.columns])
    
    x2 <- x1[order(as.double(gsub("sp", "", x1$internal_taxon_id))), ]
    
    return(x2)
    
  }
  
  if (type == "references") {
    
  }
  
  if (type == "credits") {
    
  }
}
