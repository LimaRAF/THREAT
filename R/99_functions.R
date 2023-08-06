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
                                 
                                 source("./R/99_functions.R")
                                 
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
                                     
                                     source("./R/99_functions.R")
                                     
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
                                   
                                   source("./R/99_functions.R")
                                   
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
        
        source("./R/99_functions.R")
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
      source("./R/99_functions.R")
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
      source("./R/99_functions.R")
      
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

## Preparing trai data
trait_data_prep = function(trees,trees.sa,trait.spp,trait.gen,trait.fam,generaliza=list(c("wsg_gcm3","LeafType","SeedMass_g"),c("wsg_gcm3")),add.spp=NULL,rm.colunas=NULL,inc.colunas=NULL,...) {
  

  #### TRANFORMING THE VARIABLES ###
  trees$ordem = as.character(trees$ordem)
  trees$Name_submitted = as.character(trees$Name_submitted)
  
  #### CREATING THE DOCUMENT WITH THE LIST OF DECISIONS TAKES AT EACH STEP ####
  notas = NULL
  
  #### CREATING THE LIST OF TAXA INSIDE THE ABUNDANCE DATASET ####
  taxa = trees[,c("Name_submitted","family","genus","taxon.rank")] 
  taxa = taxa[!duplicated(taxa$Name_submitted),]
  taxa = taxa[order(taxa$Name_submitted),]
  taxa = taxa[!taxa$taxon.rank %in% "unidentified",]

  if(!is.null(add.spp)) {
   taxa1 = add.spp[,c("Name_submitted","family","genus","taxon.rank")]
   taxa = rbind.data.frame(taxa,taxa1)
  }
  
  ########################################################################
  #### TREE TRAITS AVAILABLE FROM THE SOUTH AMERICA TREE SPECIES LIST ####
  ########################################################################
  #Info: threat status, endemism, establishment, dispersal syndrome and habit
  names(trees.sa)[names(trees.sa)%in%"SpeciesKTSA"] = "Name_submitted"
  taxa = merge(taxa,trees.sa[,c("Name_submitted","establishment.AF","establishment.BR","habito","succesional.group",
                                "dispersal.syndrome","geographical_distribution","status_final","status_brazil_2005",
                                "status_iucn_2005","life.form.reflora","habitat.reflora","vegtype.reflora",
                                "domain.reflora")],by="Name_submitted",all.x = TRUE)
  
  #### TRAIT DATA EDITING ####
  ## Species establishment (endemic or exotic/cultivated) in respect to the Atlantic Forest and to Brazil ##
  taxa$establishment.AF[is.na(taxa$establishment.AF)&!is.na(taxa$establishment.BR)] = 
    taxa$establishment.BR[is.na(taxa$establishment.AF)&!is.na(taxa$establishment.BR)]
  taxa$establishment.AF[taxa$establishment.AF%in%"check"&taxa$establishment.BR%in%c("native","cultivated")] = 
    taxa$establishment.BR[taxa$establishment.AF%in%"check"&taxa$establishment.BR%in%c("native","cultivated")]
  taxa$establishment.AF[taxa$establishment.AF %in% c("cultivated","naturalized","not in Brazil")] = "exotic"
  notas[1] = "The native/exotic classification is in respect to the Atlantic Forest and more broadly to Brazil"
  
  ## Extinction threat level - Brazil + IUCN (global) ##
  taxa$status_final[is.na(taxa$status_final)] = "not_evaluated_or_not_threatened"
  taxa$extinction = as.factor(gsub("_iucn|_reflora|_author|_brasil","",taxa$status_final))
  #CR=4, DD=0.5, EN=3, LC=0, NT=1, not_evaluated=0, VU=2
  niveis = c(4,0.5,3,0,1,0,2)
  names(niveis) = c("critically_endangered_CR","data_deficient_DD","endangered_EN","least_concern_LC","near_threatened_NT","not_evaluated_or_not_threatened","vulnerable_VU")
  levels(taxa$extinction) = as.character(niveis[names(niveis) %in% levels(taxa$extinction)])
  taxa$extinction = as.double(as.character(taxa$extinction))
  taxa$extinction[taxa$taxon.rank %in% c("unidentified","genus","family")] = NA
  taxa$status_final[taxa$taxon.rank %in% c("unidentified","genus","family")] = NA
  notas = c(notas,paste("The values given for CR, DD, EN, LC, NT, NE, and VU were: ",paste(niveis, collapse=","),", respectively",sep=""))
  
  ## Endemism level - in respect to the Atlantic Forest ##
  #edtiting species geographical occurrence by genus
  taxa$geographical_distribution[is.na(taxa$geographical_distribution) & taxa$genus %in% 
                                   c("Artocarpus","Citrus","Coffea","Camellia","Corymbia","Cupressus",
                                     "Dracaena","Eucalyptus","Eriobotrya","Hovenia","Leucaena","Ligustrum",
                                     "Mangifera","Mimusops","Morus","Pinus","Pittosporum","Platanus",
                                     "Ricinus","Spathodea","Syzygium","Tamarindus","Tecoma")] = "exotic"
  taxa$geographical_distribution[is.na(taxa$geographical_distribution) & taxa$genus %in% 
                                   c("Neomitranthes","Salzmannia"
                                   )] = "eastern_south_america"
  taxa$geographical_distribution[is.na(taxa$geographical_distribution) & taxa$genus %in% 
                                   c("Bocagea","Hornschuchia","Tovomitopsis","Urbanodendron","Behuria",
                                     "Macrotorus","Andradea","Ramisia","Macrothumia","Tripterodendron"
                                   )] = "regional_endemic"
  taxa$geographical_distribution[is.na(taxa$geographical_distribution) & taxa$genus %in% 
                                   c("Kuhlmanniodendron","Paratecoma","Ophthalmoblapton","Arapatiella","Brodriguesia",
                                     "Grazielodendron","Harleyodendron","Hydrogaster","Curitiba","Riodocea","Melanopsidium",
                                     "Andreadoxa","Alatococcus","Trigoniodendron"
                                     )] = "local_endemic"
  #More concise classification: C/E/N/W SA=1, exotic=-1, local=3, pan./neotropical=0, regional=2, SA=0 
  taxa$endemism = as.factor(taxa$geographical_distribution)
  niveis = c(1,1,-1,3,0,1,0,2,0,1,1,1)
  names(niveis) = c("central_south_america","eastern_south_america","exotic","local_endemic","neotropical","northern_south_america",         
                    "pantropical","regional_endemic","south_american","south_and_central_south_america","southern_south_america","western_south_america")
  levels(taxa$endemism) = as.character(niveis[names(niveis) %in% levels(taxa$endemism)])
  taxa$endemism = as.double(as.character(taxa$endemism))
  notas = c(notas,paste("The more concise values given for CSA, ESA, Exotic, Local endemic, Neotrop., NSA, Pantrop., Regional endemic, SA, SCSA, SSA and WSA were: ",paste(niveis, collapse=","),", respectively",sep=""))
  #More detailed classification: C/E/N/W SA=2, exotic=-1, local=4, pan./neotropical=0, regional=3, SA=1
  taxa$endemism1 = as.factor(taxa$geographical_distribution)
  niveis = c(2,2,-1,4,0,2,0,3,0,1,2,2)
  levels(taxa$endemism) = as.character(niveis[names(niveis) %in% levels(taxa$endemism)])
  taxa$endemism = as.double(as.character(taxa$endemism))
  notas = c(notas,paste("The more detailed values given for CSA, ESA, Exotic, Local endemic, Neotrop., NSA, Pantrop., Regional endemic, SA, SCSA, SSA and WSA were: ",paste(niveis, collapse=","),", respectively",sep=""))
  #notas = c(notas,"")
  
  ## Ecological/successional group ##
  taxa$ecological.group = taxa$succesional.group
  taxa$ecological.group[taxa$ecological.group %in% "climax"] = 4
  taxa$ecological.group[taxa$ecological.group %in% "late_secondary"] = 3
  taxa$ecological.group[taxa$ecological.group %in% c("early_secondary","early_seconday","early_secondary?")] = 2
  taxa$ecological.group[taxa$ecological.group %in% c("pioneer")] = 1
  taxa$ecological.group[taxa$ecological.group %in% c("non_pioneer")] = 2.5
  taxa$ecological.group[taxa$ecological.group %in% c("light_demanding","pioneer?")] = 1.5
  taxa$ecological.group[taxa$ecological.group %in% c("shade_tolerant","climax shade_tolerant","climax?")] = 3.5
  taxa$ecological.group = suppressWarnings(as.double(taxa$ecological.group))
  notas = c(notas,"The ecological group classification for Pioneer, Early_secondary, Late_secondary, and Climax species were 1, 2, 3, 4, respectively. In addition: Light-demanding= 1.5, Non-pioneer= 2.5 and Shade-tolerant= 3.5")
  
  ## Automatic assignment of very common pioneer, early secondary, late secondary and climax genera
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$genus %in% 
         c("Trema","Celtis","Guazuma","Bixa","Dodonea","Heliocarpus","Tecoma","Ricinus","Leucaena",
           "Manihot","Vernonanthura","Urera","Vasconcellea","Carica","Cecropia","Vismia","Chiococca",
           "Cochlospermum","Commiphora","Curatella","Prosopis","Physocalymma","Wedelia","Chromolaena")] = 1
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$genus %in% 
         c("Lithraea","Vismia","Schinus","Moquiniastrum","Lithraea","Clethra","Aparisthmium","Schizolobium",
           "Acrocomia","Schinopsis","Acca","Ligustrum","Alchornea","Erythrina","Annona","Sapium")] = 1.5
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$genus %in% 
         c("Zanthoxylum","Cedrela","Maclura","Enterolobium","Lonchocarpus","Lacistema","Tapirira",
           "Cybistax","Prockia","Coutarea","Sparattosperma","Sapindus","Balfourodendron",
           "Crateva","Hura","Simarouba","Plathymenia","Dasyphyllum","Peltophorum","Rhamnidium",
           "Bowdichia","Basiloxylon","Pterogyne","Helietta","Tetrorchidium")] = 2
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$genus %in% 
         c("Blepharocalyx","Pimenta","Galipea","Senefeldera","Schoepfia","Tetrastylidium","Amphirrhox",
           "Leretia","Platycyamus","Melanoxylon","Diploon","Sweetia")] = 3
  notas = c(notas,"Missing classification of Ecological Groups for some species was carried at genus level for some genera typically pioneer (e.g. Trema, Celtis, Guazuma, Vernonanthura, Urera, Cecropia) or late_secondary (e.g. Blepharocalyx, Senefeldera, Schoepfia, Tetrastylidium, Amphirrhox, Platycyamus, Melanoxylon, Diploon, Sweetia)")

  ## Final edits for ecological groups
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$species.correct %in% 
                          c("Myrcia citrifolia","Myrceugenia ovalifolia","Aniba intermedia")] = 3
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$species.correct %in% 
                          c("Citrus limon","Citrus sinensis","Coffea arabica","Bougainvillea praecox")] = 2
  taxa$ecological.group[is.na(taxa$ecological.group) & taxa$species.correct %in% 
                          c("Crepidospermum atlanticum","Eugenia puberula")] = 2.5
  ## Inspecting species still missing ecological groups info
  #tmp = taxa[is.na(taxa$ecological.group)&!is.na(taxa$GS),]
  #tmp1 = aggregate(tmp$N,list(tmp$species.correct),sum)
  #tmp1$GS = aggregate(tmp$GS,list(tmp$species.correct),mean)$x
  #tmp1$SuccesionalGroup = as.character(aggregate(tmp$SuccesionalGroup,list(tmp$species.correct),unique)$x)
  #tmp1[order(tmp1$x),]
  #tmp1[,c("species.correct","SuccesionalGroup","GS")]
  #tmp1[tmp1$species.correct %in% c("Myrcia citrifolia","Myrceugenia ovalifolia","Crepidospermum atlanticum","Aniba intermedia"),c("species.correct","SuccesionalGroup","GS")]

  
  ## Creating column for different life-forms (palms, ferns, etc)
  taxa$life.form = NA
  taxa$life.form[taxa$family %in% c("Arecaceae")] = "palm"
  taxa$life.form[taxa$family %in% c("Araceae","Asparagaceae","Laxmanniaceae","Velloziaceae")] = "palmoids"
  taxa$life.form[taxa$family %in% c("Cyatheaceae","Dicksoniaceae","Blechnaceae")] = "tree_fern"
  taxa$life.form[taxa$family %in% c("Cactaceae")] = "succulent_tree"
  taxa$life.form[taxa$family %in% c("Poaceae")] = "woody_bamboo"
  taxa$life.form[taxa$family %in% c("Convolvulaceae","Cucurbitaceae","Menispermaceae","Smilacaceae",
                                    "Turneraceae","Vitaceae")] = "woody_vines_and_subshrubs"
  #"Convolvulaceae" => "Ipomoea carnea"  "Ipomoea pinnata"    
  #"Cucurbitaceae" => "Cayaponia tayuya" "Wilbrandia"   
  #"Menispermaceae" => "Abuta grandifolia" "Cissampelos ovalifolia"    
  #"Smilacaceae" => "Smilax fluminensis"
  #"Turneraceae" =>  "Turnera blanchetiana" "Turnera calyptrocarpa" "Turnera diffusa" "Turnera macrophylla"
  #"Vitaceae"=> "Cissus pulcherrima"
  taxa$life.form[is.na(taxa$life.form)] = "woody_tree"
  
  ## Dispersal syndromes ## 
  #by families
  taxa$dispersal.syndrome[!is.na(taxa$family) & is.na(taxa$dispersal.syndrome) & taxa$family %in% 
                            c("Myrtaceae","Lauraceae","Annonaceae","Chrysobalanaceae","Siparunaceae","Piperaceae",
                              "Aquifoliaceae","Erythroxylaceae","Arecaceae","Burseraceae","Clusiaceae")] = "zoochoric"
  taxa$dispersal.syndrome[!is.na(taxa$Name_submitted) & grepl("Eucalyptus|Corymbia",taxa$Name_submitted)] = "anemochoric"
  taxa$dispersal.syndrome[!is.na(taxa$family) & taxa$family %in% 
                            c("Cyatheaceae")] = "anemochoric" #Asteraceae,Bignoniaceae,Vochysiaceae??
  
  #by genera
  taxa$dispersal.syndrome[!is.na(taxa$genus) & is.na(taxa$dispersal.syndrome) & taxa$genus%in% 
                            c("Miconia","Leandra","Pisonia","Guapira","Myrsine","Pouteria","Solanum","Capsicum","Ficus","Mollinedia","Diospyros","Symplocos",
                              "Ardisia","Cabralea","Calophyllum","Cereus","Cheiloclinium","Diploon","Dipteryx","Lantana","Ossaea","Prockia","Sideroxylon","Anacardium","Emmotum","Hyeronima","Macropeplus","Meliosma",
                              "Pera","Pradosia","Pseudolmedia","Sacoglottis","Sessea","Simaba","Strychnos","Tapirira","Agonandra","Alchornea","Citronella","Heisteria","Henriettea","Jacaratia","Schinus","Tetrorchidium","Urera",
                              "Vantanea","Aureliana","Bunchosia","Cecropia","Celtis","Copaifera","Coussapoa","Humiriastrum","Pourouma","Prunus","Salacia","Sapium","Stylogyne","Talisia","Brunfelsia","Citharexylum","Dendropanax",
                              "Pilocarpus","Rauvolfia","Tabernaemontana","Zollernia","Cestrum","Chionanthus","Manilkara","Meriania","Micropholis","Ormosia","Virola","Lacistema","Sorocea","Vitex","Xylosma","Brosimum","Connarus",
                              "Euplassa","Guarea","Matayba","Schefflera","Aegiphila","Picramnia","Chrysophyllum","Daphnopsis","Styrax","Andira","Coccoloba","Cybianthus","Neea","Mouriri","Allophylus","Ouratea","Zanthoxylum",
                              "Cupania","Byrsonima","Casearia","Maytenus","Trichilia","Swartzia","Inga","Rubus","Citrus",
                              "Alibertia","Amaioua","Chiococca","Chomelia","Cordiera","Coussarea","Faramea","Genipa","Guettarda","Ixora","Kutchubaea","Palicourea","Posoqueria","Psychotria","Randia","Rudgea","Tocoyena",
                              "Vismia","Quiina","Quararibea","Phytolacca","Homalolepis","Simaba","Capparidastrum","Neocalyptrocalyx","Athenaea","Monteverdia"
                            )] = "zoochoric"
  taxa$dispersal.syndrome[!is.na(taxa$genus) & is.na(taxa$dispersal.syndrome) & taxa$genus%in% 
                            c("Gymnanthes","Galipea","Jatropha","Pachystroma","Sebastiania","Cnidoscolus","Calliandra","Cassia","Esenbeckia","Helicteres","Mabea","Metrodorea","Mimosa","Actinostemon","Manihot","Senegalia",
                              "Senna","Croton","Senefeldera","Algernonia","Brasiliocroton","Stillingia"
                            )] = "autochoric"
  taxa$dispersal.syndrome[!is.na(taxa$genus) & is.na(taxa$dispersal.syndrome) & taxa$genus%in% 
                            c("Adenocalymma","Pachira","Peltophorum","Plathymenia","Raulinoreitzia","Bougainvillea","Cedrela","Hymenolobium","Myrocarpus","Pinus","Pterodon","Ruprechtia","Alsophila","Lafoensia","Moquiniastrum","Seguieria",
                              "Symphyopappus","Triplaris","Callisthene","Centrolobium","Pterocarpus","Astronium","Cariniana","Lonchocarpus","Pseudobombax","Ceiba","Luehea","Tabebuia","Tachigali","Vernonanthura","Eriotheca","Terminalia",
                              "Handroanthus","Piptocarpha","Jacaranda","Cyathea","Dalbergia","Qualea","Tibouchina","Handroanthus","Vochysia","Baccharis","Kielmeyera","Aspidosperma","Machaerium","Vernonia","Eupatorium",
                              "Alseis","Bathysa","Coutarea","Rustia","Simira","Combretum","Luetzelburgia","Pleroma"
                            )] = "anemochoric"
  
  #tmp = taxa[is.na(taxa$dispersal.syndrome),]
  #tmp = tmp[!grepl(" ",tmp$Name_submitted),]
  #sort(table(tmp$Name_submitted))
  #sort(table(tmp$genus))
  notas = c(notas,"Missing classification of Dispersal Syndrome for some species or morphospecies was carried at genus or family level (e.g. all Miconia are zoochorich and all ferns are anemochoric")
  
  ####################################################################
  #### MEAN TREE TRAITS AVAILABLE FROM THE TREECO TRAIT DATABASE  ####
  ####################################################################
  
  #### MERGING TRAIT AT SPECIES LEVEL ####
  tmp = trait.spp[trait.spp$TAX %in% taxa$Name_submitted,]
  taxa = merge(taxa,tmp,by.x="Name_submitted",by.y="TAX",all.x=T)
  
  ##CHECK!!
  #Checkings possible problems
  #tax[!is.na(tax$wsg_gcm3)&tax$wsg_gcm3>1&tax$GS<2,c("species.correct","wsg","wsg_gcm3","SuccesionalGroup","GS")]
  #tax[!is.na(tax$wsg_gcm3)&tax$wsg_gcm3<0.6&tax$GS>=3.5,c("species.correct","wsg","wsg_gcm3","SuccesionalGroup","GS")]
  
  #### PRELIMINAR EDITING OF TRAIT DATA ####
  #### UPDATE THIS PART OF THE CODE FROM THE CODES OF THE PAPER SUBMITTED TO SCIENCE ####
  ## Seed mass ~ Seed volume
  taxa$SeedVolume = NA
  #for circular fruits
  taxa$SeedVolume[is.na(taxa$SeedLength_cm)&is.na(taxa$SeedWidth_cm)&!is.na(taxa$SeedDiameter_cm)] = 
    (4/3)*pi*(as.double(taxa$SeedDiameter_cm[is.na(taxa$SeedLength_cm)&is.na(taxa$SeedWidth_cm)&!is.na(taxa$SeedDiameter_cm)])/2)^3
  #for ellipsoidal fruits
  taxa$SeedVolume[!is.na(taxa$SeedLength_cm)&!is.na(taxa$SeedWidth_cm)&!is.na(taxa$SeedDiameter_cm)] = 
    (4/3)*pi*(as.double(taxa$SeedDiameter_cm[!is.na(taxa$SeedLength_cm)&!is.na(taxa$SeedWidth_cm)&!is.na(taxa$SeedDiameter_cm)])/2)*
    (as.double(taxa$SeedWidth_cm[!is.na(taxa$SeedLength_cm)&!is.na(taxa$SeedWidth_cm)&!is.na(taxa$SeedDiameter_cm)])/2)*
    (as.double(taxa$SeedLength_cm[!is.na(taxa$SeedLength_cm)&!is.na(taxa$SeedWidth_cm)&!is.na(taxa$SeedDiameter_cm)])/2)
  #plot(log(taxa$SeedMass_g) ~ log(taxa$SeedVolume)); abline(0,1,lwd=2,col=2)
  #text(x=log(taxa$SeedVolume), y=log(taxa$SeedMass_g), labels=taxa$Name_submitted,cex=0.5,pos=1:3)
  
  tmp = na.omit(taxa[,c("Name_submitted","SeedMass_g","SeedVolume")])
  mod = lm(log(as.double(tmp$SeedMass_g[!tmp$Name_submitted %in% c("Miconia prasina","Myrcia pubipetala")])) ~ log(as.double(tmp$SeedVolume[!tmp$Name_submitted %in% c("Miconia prasina","Myrcia pubipetala")])))
  #abline(mod)
  ## CHECK!!
  ##RE-CHECK THE EQUATION USING ONLY SPECIES WITH ALL SEED INFO (RAW TRAIT VALUES): length/width and/or diameter and seed mass
  taxa$est_SeedMass_g = exp(coef(mod)[1] + coef(mod)[2]*log(taxa$SeedVolume))
  #plot(taxa$est_SeedMass_g ~ taxa$SeedMass_g, log="xy"); abline(0,1,lwd=2,col=2)
  #taxa$SeedMass_g[is.na(taxa$SeedMass_g)&!is.na(taxa$est_SeedMass_g)] = 
  #  taxa$est_SeedMass_g[is.na(taxa$SeedMass_g)&!is.na(taxa$est_SeedMass_g)]
  
  #seed mass <- 0.8*fruits/kg
  #plot(SeedMass_g ~ FruitMass_g, data=trait.spp, log="xy",
  #     ylim=c(6.150288e-06, 1.157476e+02),
  #     xlim=c(0.00001, 850)); abline(0,1,lwd=2,col=2)
  #tmp = na.omit(trait.spp[,c("SeedMass_g","FruitMass_g")])
  #mod = lm(log(tmp$SeedMass_g) ~ log(tmp$FruitMass_g))
  #abline(mod)
  #abline(-0.5,0.85)
  #text(x=trait.taxa$FruitMass_g, y=trait.taxa$SeedMass_g, labels=trait.taxa$TAX,cex=0.4,pos=1:3)
  ##RE-CHECK THE EQUATION USING ONLY SPECIES WITH ALL SEED INFO (RAW TRAIT VALUES): fruit mass and seed mass
  taxa$est_SeedMass_g[is.na(taxa$est_SeedMass_g)] = exp(coef(mod)[1] + coef(mod)[2]*log(taxa$FruitMass_g[is.na(taxa$est_SeedMass_g)]))
  taxa$SeedMass_g[is.na(taxa$SeedMass_g)&!is.na(taxa$FruitMass_g)] = 
    taxa$est_SeedMass_g[is.na(taxa$SeedMass_g)&!is.na(taxa$FruitMass_g)]
  notas = c(notas,"Missing information of Seed Mass was also obtained from studies providing values of number of seeds per kilogram. In additin, for some monospermic species was based on a conversion formula using Fruit_mass as an input")
  
  #### MERGING TRAIT AT GENUS LEVEL ####
  gen = as.character(trait.gen$TAX[trait.gen$TAX %in% taxa$genus])
  trait = generaliza[[1]]  
  for(j in 1:length(gen)) {
    tmp = taxa[taxa$genus %in% gen[j],]
    for(i in 1:length(trait)) {
      tmp[is.na(tmp[,trait[i]]),trait[i]] = as.character(trait.gen[as.character(trait.gen$TAX) %in% gen[j],trait[i]])
    }
    taxa[taxa$genus %in% gen[j],trait] = tmp[,trait]
  }
  notas = c(notas,paste("Missing info at species level was replaced by genus-level averages for the following traits: ",paste(generaliza[[1]], collapse=","),".",sep=""))
  
  #### MERGING TRAIT AT FAMILY LEVEL ####
  fam = as.character(trait.fam$TAX[trait.fam$TAX %in% taxa$family])
  trait = generaliza[[2]] 
  for(j in 1:length(fam)) {
    tmp = taxa[taxa$family %in% fam[j],]
    for(i in 1:length(trait)) {
      tmp[is.na(tmp[,trait[i]]),trait[i]] = as.character(trait.fam[as.character(trait.fam$TAX) %in% fam[j],trait[i]])
    }
    taxa[taxa$family %in% fam[j],trait] = tmp[,trait]
  }
  notas = c(notas,paste("Missing info at species level was replaced by family-level averages for the following traits: ",paste(generaliza[[2]], collapse=","),".",sep=""))
  
  #### EDITING MISSING TRAIT DATA FOR ABUNDANT OR FREQUENT SPECIES IN TREECO ###
  ## Maximum Height (m) ##
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Vellozia variabilis"] = 2
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Vellozia tubiflora"] = 2
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Ipomoea carnea"] = 2.5
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Erythroxylum nummularium"] = 4
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Cratylia mollis"] = 4
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Hyptis pachyphylla"] = 3
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Croton adenodontus"] = 3
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Eriope exaltata"] = 6
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Simaba ferruginea"] = 13
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Homalolepis ferruginea"] = 13
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Erythroxylum pauferrense"] = 6
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Ephedranthus dimerus"] = 20
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Myrcia costeira"] = 9
  taxa$MaxHeight_m[!is.na(taxa$Name_submitted) & taxa$Name_submitted == "Psychotria pubigera"] = 4
 
  ## Getting leaf Length and Width from the speciesLink measurements ##
  tmp = read.csv("data/data-raw/leaf_size_speciesLink.csv",as.is=T)


  tmp = tmp[!is.na(tmp$Length),]
  tmp1 = cbind.data.frame(aggregate(tmp$Length,list(tmp$species),mean),
                          W=aggregate(tmp$Width,list(tmp$species),mean)$x)
  names(tmp1)[1:2] = c("species","L")
  tmp1$LA = pi*tmp1$L*tmp1$W/4
  notas = c(notas,"Missing information of Leaf Area was obtained from Leaf_length and Leaf_width using the formula of an ellipse: pi*L*W/4")
  
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Astrocaryum aculeatissimum"] = tmp1$LA[tmp1$species == "Astrocaryum aculeatissimum"]
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Euterpe edulis"] = tmp1$LA[tmp1$species == "Euterpe edulis"]
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Marlierea racemosa"] = tmp1$LA[tmp1$species == "Marlierea racemosa"]
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Myrcia aethusa"] = tmp1$LA[tmp1$species == "Myrcia aethusa"]
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Myrcia pubipetala"] = tmp1$LA[tmp1$species == "Myrcia pubipetala"]
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Pachystroma longifolium"] = tmp1$LA[tmp1$species == "Pachystroma longifolium"]
  taxa$LeafArea[!is.na(taxa$Name_submitted)&taxa$Name_submitted == "Astrocaryum aculeatissimum"] = tmp1$LA[tmp1$species == "Astrocaryum aculeatissimum"]
  
  ## Editing LeafType data ##
  taxa$LeafType[grepl("compostas",taxa$LeafType)&!grepl("simples",taxa$LeafType)] = "compostas"
  taxa$LeafType[!grepl("compostas",taxa$LeafType)&grepl("simples",taxa$LeafType)] = "simples"
  taxa$LeafType[!grepl("compostas",taxa$LeafType)&!grepl("simples",taxa$LeafType)] = NA
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Annona",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Bauhinia rufa",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Cheiloclinium",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Conchocarpus fontanesianus",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Conchocarpus gaudichaudianus",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Conchocarpus pentandrus",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Esenbeckia grandiflora",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Matayba",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Pseudobombax simplicifolium",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Roupala montana",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Swartzia simplex",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Tabebuia insignis",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Vitex mexiae",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Weinmannia humilis",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Aralia warmingiana",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Eriotheca pubescens",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Malouetia cestroides",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Joannesia princeps",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Malouetia cestroides",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Piptocarpha axillaris|Piptocarpha macropoda|Piptocarpha regnellii",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Protium heptaphyllum|Protium spruceanum|Protium widgrenii",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Schefflera calva|Schefflera macrocarpa",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Tabebuia roseoalba|Tabebuia aurea",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Vitex cymosa|Vitex megapotamica|Vitex polygama",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("Zeyheria montana|Zeyheria tuberculosa",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Actinostemon",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Allophylus",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Callisthene",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Celtis",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Dendropanax",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Faramea",taxa$Name_submitted)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Jacaranda",taxa$Name_submitted)] = "compostas"
  taxa$LeafType[!is.na(taxa$Name_submitted)&grepl("\\|",taxa$LeafType)&grepl("Sebastiania",taxa$Name_submitted)] = "simples"
  taxa$LeafType[grepl("compostas",taxa$LeafType)&grepl("simples",taxa$LeafType)] = NA
  taxa$LeafArea[is.infinite(taxa$LeafArea)] = NA

  #some final edits at genus or family level
  taxa$LeafType[!is.na(taxa$family)&is.na(taxa$LeafType)&grepl("Annonaceae|Melastomataceae|Myrtaceae",taxa$family)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted)&is.na(taxa$LeafType)&grepl("Beilschmiedia|Mezilaurus",taxa$Name_submitted)] = "simples"

  #some final edits at species level (looked up in herbarium vouchers)
  taxa$LeafType[!is.na(taxa$Name_submitted) & is.na(taxa$LeafType) & taxa$Name_submitted %in% 
			c("Bauhinia_longifolia","Cavanillesia_umbellata","Guapira_areolata","Pogonophora_schomburgkiana"
	)] = "simples"
  taxa$LeafType[!is.na(taxa$Name_submitted) & is.na(taxa$LeafType) & taxa$Name_submitted %in% 
			c("Caryocar_brasiliense","Cenostigma_bracteosum","Connarus_suberosus","Esenbeckia_febrifuga","Pachira_endecaphylla",
			"Galipea_jasminiflora","Protium_atlanticum","Protium_icicariba","Sparattosperma_leucanthum"
	)] = "compostas"
  notas = c(notas,"Missing information of Leaf Type was obtained at genus-level for some taxa. Species with both compound and simple leafs along their onthogeny (e.g. Matayba, Roupala) were treated as having compound leaves")
  
  ### Final edits: Seed size
  taxa$SeedMass_g[is.na(taxa$SeedMass_g) & taxa$genus %in% c("Actinostemon")] = mean(as.double(unique(taxa$SeedMass_g[taxa$genus %in% c("Sebastiania")])))
  taxa$SeedMass_g[is.na(taxa$SeedMass_g) & taxa$genus %in% c("Athenaea")] = mean(as.double(unique(taxa$SeedMass_g[taxa$genus %in% c("Capsicum")])))
  taxa$SeedMass_g[is.na(taxa$SeedMass_g) & taxa$genus %in% c("Anaxagorea")] = mean(as.double(unique(taxa$SeedMass_g[taxa$genus %in% c("Xylopia")])))
  
  ### Edits for very species rich genera without habit
  ##VERIFICAR A PROPORÇÃO DE ÁRVORES E ARBUSTOS POR FAMILIA E GENERO PARA FECHAR AS LISTAGEMS ABAIXO, PRINCIPALMENTE GENEROS
  #By family
  taxa$habito[!is.na(taxa$family) & is.na(taxa$habito) & taxa$family%in% 
                c("Annonaceae","Burseraceae","Chrysobalanaceae","Lauraceae","Lecythidaceae","Moraceae","Myristicaceae","Myrtaceae","Proteaceae","Sapotaceae"
                )] = 1
  taxa$habito[!is.na(taxa$family) & is.na(taxa$habito) & taxa$family%in% 
                c("Lacistemataceae","Loganiaceae","Pipeaceae","Siparunaceae","Violaceae"
                )] = 0.5
  #By genus
  taxa$habito[!is.na(taxa$genus) & is.na(taxa$habito) & taxa$genus%in% 
                c("Alchornea","Anadenanthera"
                )] = 1
  taxa$habito[!is.na(taxa$genus) & is.na(taxa$habito) & taxa$genus%in% 
                c("Athenaea","Brunfelsia","Capsicum","Piper","Eumachia","Chiococca","Stachytarpheta","Schweiggeria"
                )] = 0.5
    
  #### CHECKING COMMON SPECIES STILL WITHOUT TRAIT INFO ####
  tmp = c("establishment.AF","habito","extinction","endemism","ecological.group","dispersal.syndrome",
          "wsg_gcm3","MaxHeight_m","MaxDbh_cm",
          "SeedLength_cm","SeedWidth_cm","SeedDiameter_cm","SeedMass_g","FruitMass_g",
          "LeafArea","LeafType","SuccesionalGroup",
          "Deciduousness","LT","PericarpConsistency","SexualSystem","PollinationSyndrome","SLA_cm2_g.1")
  tmp1= NULL
  for(i in 1:length(tmp)) {
    vrvl = tmp[i]
    if(is.numeric(taxa[,vrvl])) { tmp2 = taxa$Name_submitted[is.na(taxa[,vrvl])] 
    } else { tmp2 = taxa$Name_submitted[is.na(taxa[,vrvl])|(grepl("NA:",taxa[,vrvl])&!grepl("\\|",taxa[,vrvl]))]
    }
    tmp2 = tmp2[!grepl(" sp\\.",tmp2)&!is.na(tmp2)]
    tmp2 = tmp2[grepl(" ",tmp2)]
    if(is.null(tmp1)) { tmp1 = c(length(table(tmp2)),100*sum(table(tmp2))/dim(taxa)[1]) 
    } else { tmp1 = rbind.data.frame(tmp1,c(length(table(tmp2)),100*sum(table(tmp2))/dim(taxa)[1])) }
  }
  names(tmp1) =c("Nsp_without","prop.sp.without")
  rownames(tmp1) = tmp
  print("Check the number and proportion of species (not genera or families) without trait information.")
  print("Proportion below 40% generally results in good coverage for all individuals at the community level.")
  print(tmp1[order(tmp1$prop.sp.without),])
  
  #### FILTERING THE SELECTED METADATA (COLUMNS) ####
  colunas = c("X","species.original","Habit","Bark","LeafMargin","LeafTip","AcumenLength_cm","LeafSubtype","Leaflets",
              "LeafHairiness_upper","LeafHairiness_lower","FruitLipids","FruitWater","FruitProtein","FruitCarbo",
              "PollinationSubsyndrome","Anthesis","FlowerType","CorollaType","CorollaColour","SepalColour","CorollaLength_cm","CorollaWidth_cm",
              "ChemicalCompounds","PhysicalCompounds","Allelopathy","sem.kg","Flowers","Fruits","Seeds","LeafN","obs1","obs2","obs3",
              "GeogDist","Domain","Vegetation","TaxNotes","MinLeafLength_cm","MaxLeafLength_cm","MinLeafWidth_cm","MaxLeafWidth_cm",
              "LeafThickness_micro","FruitLength_cm","FruitWidth_cm","FruitDiameter_cm","LeafChlorophyllAB","LeafCarbon","LeafWater",
              "dry_mass_fraction","WoodCarbon_gkg.1","LeafDryMatter_mg.g","ChlorophyllContent_FCI","LeafTough_N.mm","PithProp","XylProp","BarkProp",
              "BarkProp","VesselDensity_n.mm2","VesselDiameter_microm","StomataDensity_n.mm2","StomataLength_microm","Kpot_kg.m.s.Mpa",
              "SuccesionalGroup","wsg")
  if(is.null(rm.colunas)) {colunas1 = colunas}  else { colunas1 = c(colunas,rm.colunas) }
  if(is.null(inc.colunas)) {colunas1 = colunas1}  else { colunas1 = colunas1[!colunas1 %in% inc.colunas] }
  if(!is.null(colunas1)) taxa = taxa[,!names(taxa) %in% colunas1]
  
  #### PREPARING AND RETURNING THE FINAL RESULTS ####
  notas = gsub("\\\n","",notas)
  notas1 = paste(">>>NOTE ",1:length(notas),": ",notas,sep="")
  resulta = list(taxa, notas1)
  names(resulta) = c("data_frame","editing_notes")
  return(resulta)
}


#' Internal function
#'
#' Get Smallest Distance Between Occurrences
#'
#' @param points1 XY data frame (test points)
#' @param points2 XY data frame (reference points)
#' @param value output value: distance ("dist") or smaller/greater than `min.dist` 
#'   ("flag")?
#' @param min.dist minimum tolerated distance between points, in meters.
#'   Default to 1000 m.
#' 
#' @details The first two arguments are a XY data frame has the same structure
#'   as other XY objects within `ConR` with the three first columns being
#'   `ddlat`, `ddlon` (both in decimal degrees) and `tax`. In the context of
#'   species assessments, the third column `tax` refers to the name of a taxon
#'   but it practice it can be any grouping character. At least one of the
#'   values in `tax` must be common to the two data frames. If one of the taxa
#'   in `points1` in not in `points2`, distances cannot be computed and NAs are
#'   returned.
#'
#' @return The data frame in provided in `points1` with an aditional column with
#'   the smallest distance from the points in `points2` for each taxon.
#'   
#' @examples
#' 
#' mydf1 <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0),
#'                    ddlon = c(-42.2,-41.6,-45.3,-42.5,-42.3,-39.0),
#'                    tax = rep(c("a","b"), each = 3),
#'                    stringsAsFactors = FALSE)
#' mydf2 <- data.frame(ddlat = runif(20, -47, -41),
#'                    ddlon = runif(20, -46, -38),
#'                    tax = rep(c("a","b"), each = 10),
#'                    stringsAsFactors = FALSE)
#' plot(mydf2[,2:1], pch =19)
#' points(mydf1[,2:1], pch =19, col =2)
#' .dist.valid.points(mydf1, mydf2)
#' .dist.valid.points(mydf1, mydf2, min.dist = 100000, value = "flag")
#' 
#' @importFrom fields rdist.earth
#'
.dist.valid.points <- function(points1, points2, 
                               min.dist = 1000, 
                               value = "dist", 
                               parallel = TRUE,
                               NbeCores = 2) {
  
  if (all(!unique(points1$tax) %in% unique(points2$tax))) {
    warning("no match between characters in points 1 and points 2: returning NAs")
    return(cbind.data.frame(points1, 
                            distance = rep(NA, nrow(points1)),
                            stringsAsFactors = FALSE))
  }
  
  points2 <- points2[points2$tax %in% unique(points1$tax),]
  
  if(any(!unique(points1$tax) %in% unique(points2$tax))) {
    miss.spp <- unique(points1$tax)[!unique(points1$tax) %in% unique(points2$tax)]
    points2 <- rbind.data.frame(points2,
                                data.frame(ddlat = rep(NA, length(miss.spp)),
                                           ddlon = rep(NA, length(miss.spp)),
                                           tax = miss.spp,
                                           stringsAsFactors = FALSE)
                                )  
  }
    
  list_data1 <- split(points1, f = points1$tax)
  list_data2 <- split(points2, f = points2$tax)

  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
  
  x <- NULL
  output <-
    foreach::foreach(
      x = 1:length(list_data1),
      .combine = 'c',
      .options.snow = NULL
    ) %d% {
      #source("C://Users//renato//Documents//raflima//R_packages//ConR//R//over.valid.poly.R")

      # res <- RANN::nn2(list_data2[[x]][,c("ddlon","ddlat")],
      #                  list_data1[[x]][,c("ddlon","ddlat")],  k=1)$nn.dists
      # res <- apply(
      #           fields::rdist(list_data1[[x]][,c("ddlon","ddlat")], 
      #                       list_data2[[x]][,c("ddlon","ddlat")]), 
      #           1, min, na.rm = TRUE)
      res <- apply(
                fields::rdist.earth(list_data1[[x]][,c("ddlon","ddlat")], 
                          list_data2[[x]][,c("ddlon","ddlat")], miles = FALSE), 
                1, min, na.rm = TRUE)*1000
      # res <- apply(
      #           geosphere::distm(list_data1[[x]][,c("ddlon","ddlat")],
      #                           list_data2[[x]][,c("ddlon","ddlat")]), 
      #           1, min, na.rm = TRUE)
      res
    }
  
  if(parallel) snow::stopCluster(cl)

  output[is.infinite(output)] <- NA
  
  if(value == "flag") {
    flag <- output < min.dist
    points1$distance <- flag
    return(points1)
  }
  
  if(value == "dist") {
    points1$distance <- output
    return(points1)
  }
}  
#' @Title Get Occupied Cells and Patches

#' @param EOO.poly spatial object containig species EOO
#' @param hab.map raster, raster layer/stack or spatial polygons containing 
#'  the habitat spatial informaton 
#' @param hab.class classess of values in ```hab.map``` to be considered as 'habitat'
#' @param years numeric. Time interval between the first and last ```hab.map```
#'  if more than one raster is provided (e.g. if ```hab.map``` is a RasterStack object)
#' @param export_shp logical. Should the species habitat map be exported? Default to FALSE.
#' @param plot logical. Should the habitat change map be plotted? Default to FALSE. 
#' @param output_raster the output format in the case raster are quantities and not classes.
#'  Default to "summary".
#' @param proj_type character string, numeric or object of CRS class. Default to "cea".
#' @param parallel logical. Should computing run in parallel? Default to FALSE.
#' @param NbeCores integer. The number of cores for parallel execution. Default to 2.
#' @param show_progress logical. whether a bar showing progress in computation should be shown. 
#'  Default to TRUE.
#' 
#' @details The function compute the amount of habitat within the EOO of each species.
#' If two or more time intervals are provided, the function also returns the values
#' of habitat change. 
#'    
#' If a RasterStack is provided, the last stacked raster is taken as the most recent 
#' raster to calculate current habitat proportion and area. Similarly, the first and 
#' last rasters are used to calculate habitat loss.      
#'    
#'    
#' @examples #To be included
#' 
#' @importFrom raster crs crop extract mask
#' @importFrom sf st_crs st_transform st_intersection st_intersects st_as_sf st_geometry
#' 
#' 
#' 
EOO.habitat <- function(EOO.poly,
                        hab.map = NULL,
                        hab.class = NULL,
                        years = NULL,
                        export_shp = FALSE,
                        plot = FALSE,
                        output_raster = "summary",
                        proj_type = "cea",
                        parallel = FALSE,
                        NbeCores = 2,
                        show_progress = TRUE
) {

  proj_crs <- my.proj_crs(proj_type = proj_type)
  
  proj_hab <- raster::crs(hab.map)
  
  if(is.na(sf::st_crs(EOO.poly))) 
    sf::st_crs(EOO.poly) <- 4326
  
  EOO.poly <- sf::st_transform(EOO.poly, crs = proj_hab)
  
  # Croping EOO shapefiles to the raster extent
  EOO.poly.crop <- suppressWarnings(
    sf::st_intersection(EOO.poly, sf::st_as_sfc(sf::st_bbox(hab.map))))
  
  if (length(EOO.poly$geometry) > length(EOO.poly.crop$geometry))
    warning(paste0("The EOO of ",
                   length(EOO.poly$geometry) - length(EOO.poly.crop$geometry),
                   " species became empty after croping them to the habitat map extent"))
  
  EOO.poly.proj <- sf::st_transform(EOO.poly.crop, crs = proj_crs)
  EOO.poly.area <- as.double(sf::st_area(EOO.poly.proj))/1000000 # area in km2
  EOO.poly.name <- EOO.poly.crop
  sf::st_geometry(EOO.poly.name) <- NULL
  EOO.poly.name <- EOO.poly.name[,1] 
  
  if(any(grepl("SpatialPolygon", class(hab.map))))
    hab.map <- sf::st_as_sf(hab.map)
  
  if (is.null(hab.class) & any(grepl("Raster", class(hab.map))))
    classes <- seq(hab.map[[1]]@data@min, hab.map[[1]]@data@max, 1)

  # Extracting the map information for each EOO polygon 
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    `%d%` <- foreach::`%dopar%`
  } else{
    `%d%` <- foreach::`%do%`
  }
  
  x <- NULL
  if (show_progress) {
    pb <-
      txtProgressBar(
        min = 0,
        max = length(EOO.poly.crop$geometry),
        style = 3
      )
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else{
    opts <- NULL
  }
  
  output <-
    foreach::foreach(
      x = 1:length(EOO.poly.crop$geometry),
      #x = 4000:5000,
      #.combine = 'c',
      #.packages=c("raster","sf"),
      .options.snow = opts
    ) %d% {
      
    # if (!parallel & show_progress)
    #   setTxtProgressBar(pb, x)
    
    # for(x in 1:100) {
    # cat(x, "\n")  
    if(any(grepl("sf", class(hab.map)))) {

      crop1 <- sf::st_intersection(hab.map, EOO.poly.crop[x, ]) # cropping to the extent of the polygon
      crop.proj <- sf::st_transform(crop1, crs = proj_crs)
      areas.x <- try(sf::st_area(crop1), TRUE)
      if (class(areas.x) == "try-error") {
        sf::sf_use_s2(FALSE)
        crop1 <- sf::st_make_valid(crop1)
        areas.x <- try(sf::st_area(crop1), TRUE)
        sf::sf_use_s2(TRUE)
      }
      
      area.hab <- as.numeric(sum(areas.x, na.rm = TRUE))/1000000
      tmp <- c(as.numeric(sf::st_area(EOO.poly.crop[x, ]))/1000000 - area.hab, area.hab)
      tmp1 <- crop1
      sf::st_geometry(tmp1) <- NULL
      n.polys <- length(unique(as.character(tmp1[,1]))) ### CHECK HERE: FAIL TO GET THE NUMBER OF SHAPE IDs
      hab.mat <- matrix(c(NA,  n.polys, 100*tmp/EOO.poly.area[x], tmp), ncol = 3, nrow = 2, byrow = FALSE)
      row.names(hab.mat) <- c("non_habitat", "habitat")
      colnames(hab.mat) <- c("numb.polys","prop.EOO", "area.EOO")
      res <- cbind.data.frame(tax = EOO.poly.name[x], hab.mat,
                              stringsAsFactors = FALSE)["habitat",]

      if(export_shp) {

        ## Tried to combine the two sf objects but could not make it with the df info
        # toto <- c(st_geometry(crop1), st_geometry(EOO.poly.crop[x, ]))
        # merge(toto, crop1)
        # toto <- sf::st_multipolygon(st_geometry(eoo1), st_geometry(crop1))
        # crop1 <- st_cast(crop1, "POLYGON")
        # eoo1 <- st_cast(EOO.poly.crop[x,], "POLYGON")

        res <- list(res, crop1)

        if (plot) {
          par(mar = c(3,3,2,2), las=1, tcl=-0.25, mgp=c(2.5,0.5,0))
          plot(sf::st_geometry(EOO.poly.crop[x, ]),
               main = EOO.poly.name[x], bg=0)
          plot(sf::st_geometry(crop1), add = TRUE, col = 1)
        }

      }

    }

  # result <- NULL
   # for(x in 1001:2000) {
   #   cat(x, "\n")
    if(any(grepl("Raster", class(hab.map)))) {
      
      ext <- raster::extent(sf::st_bbox(EOO.poly.crop[x, ])[c(1,3,2,4)])
      crop1 <- raster::crop(hab.map, ext) # cropping raster to the extent of the polygon
      mask1 <- raster::mask(crop1, mask = as(EOO.poly.crop[x, ], "Spatial"))
      tmp  <- raster::getValues(mask1) # much faster than raster::extract

      if("numeric" %in% class(tmp) | "integer" %in% class(tmp))
        tmp <- tmp[!is.na(tmp)] # excluding pixels outside the EOO
      if("matrix" %in% class(tmp) | "data.frame" %in% class(tmp))
        tmp  <- tmp[!is.na(tmp[,1]), , drop = FALSE] # excluding pixels outside the EOO
      
      if (dim(tmp)[1] > 0) {
        #Getting the approximated area of the masked raster (in km2)
        mask.area <- raster::area(mask1[[1]], na.rm = TRUE)
        r.area <- raster::cellStats(mask.area, 'sum')
        
        if (!is.null(hab.class)) {
          
          # Getting habitat quantity (land-use maps)
          mask.dim <- dim(tmp)
          pixs <- mask.dim[1]
          last <- mask.dim[2]
          hab.mat <- as.matrix(table(factor(tmp[, last] %in% hab.class, levels = c(FALSE, TRUE))))
          row.names(hab.mat) <- c("non_habitat", "habitat")
          hab.mat <- cbind(hab.mat, 100 * hab.mat/pixs)
          hab.mat <- cbind(hab.mat, (hab.mat[,2]* r.area)/100)
          colnames(hab.mat) <- c("n.pixs", "prop.EOO", "area.EOO")
          
          if (class(crop1) == "RasterBrick") {
            
            # Getting the overal percentage of each transition
            transit.mat <- 100 * as.matrix(table(factor(tmp[, 1] %in% hab.class, levels = c(FALSE, TRUE)),
                                                 factor(tmp[, last] %in% hab.class, levels = c(FALSE, TRUE)))) / pixs
            colnames(transit.mat) <- c("non_habitat_t1", "habitat_t1")
            row.names(transit.mat) <- c("non_habitat_t0", "habitat_t0")
            hab.mat <- cbind(hab.mat, loss = c(transit.mat[1, 2], transit.mat[2, 1]))
            hab.mat <- cbind(hab.mat, recover = c(transit.mat[2, 1], transit.mat[1, 2]))
            
            # % of habitat loss and recovery
            hab_loss <- 100 * transit.mat[2, 1] / sum(transit.mat[2, ])
            non_hab_loss <- 100 * transit.mat[1, 2] / sum(transit.mat[1, ])
            hab_rec <- 100 * transit.mat[1, 2] / sum(transit.mat[, 2])
            non_hab_rec <- 100 * transit.mat[2, 1] / sum(transit.mat[, 1])
            hab.mat <- cbind(hab.mat, rel.loss = c(non_hab_loss, hab_loss))
            hab.mat <- cbind(hab.mat, rel.recover = c(non_hab_rec, hab_rec))
            
            # rate of loss and recovery
            if(is.null(years))
              years <- dim(hab.map)[3] - 1
            
            non_hab_loss_rate <- non_hab_loss / years
            hab_loss_rate <- hab_loss / years
            non_hab_rec_rate <- non_hab_rec / years
            hab_rec_rate <- hab_rec / years
            hab.mat <- cbind(hab.mat, rate.loss = c(non_hab_loss_rate, hab_loss_rate))
            hab.mat <- cbind(hab.mat, rate.recover = c(non_hab_rec_rate, hab_rec_rate))
            
            # defining habitat quality at time t+1
            transit <- matrix(NA, ncol = 1, nrow = pixs)
            transit[tmp[, 1] %in% hab.class &
                      tmp[, last] %in% hab.class, 1] <- 2
            transit[tmp[, 1] %in% hab.class &
                      !tmp[, last] %in% hab.class, 1] <- -2
            transit[!tmp[, 1] %in% hab.class &
                      tmp[, last] %in% hab.class, 1] <- 1
            transit[!tmp[, 1] %in% hab.class &
                      !tmp[, last] %in% hab.class, 1] <- 0
            mod <- stats::lm(base::jitter(transit, factor=0.1) ~ 1)
            ci <- round(stats::confint(mod),4)
            hab.mat <- cbind(hab.mat,
                             hab.quality = round(as.numeric(stats::coef(mod)),4),
                             hab.quality.lo = ci[1],
                             hab.quality.hi = ci[2])
            
            res <- cbind.data.frame(tax = EOO.poly.name[x], hab.mat[,-1],
                                    stringsAsFactors = FALSE)["habitat",]
            
            # location of decline and recover
            if (export_shp) {
              
              loc.loss <- data.frame(
                cell = 1:length(mask1[[1]]@data@values),
                raster::coordinates(mask1[[1]]),
                classes = NA_integer_
              )
              #loc.loss$classes[loc.loss$cell %in% tmp[, 1]] <- transit
              loc.loss$classes[!is.na(mask1[[1]]@data@values)] <- transit
              r <- raster::rasterFromXYZ(loc.loss[!is.na(loc.loss$classes), 2:4],
                                         crs = raster::crs(mask1[[1]]))
              
              if (plot) {
                par(mar = c(3,3,2,2), las=1, tcl=-0.25, mgp=c(2.5,0.5,0))
                raster::plot(r, col = c("red", "grey", "green", "darkgreen"),
                             main = EOO.poly.name[x])
                plot(sf::st_geometry(EOO.poly.crop[x,]), add = TRUE, bg = 0)
              }
              
              res <- list(res, r)
              
            }
          }
          
        } else {
          
          #Getting summary stats for the variable in the EOO (quantitative maps)
          if("numeric" %in% class(tmp) | "integer" %in% class(tmp)) {
            vals <- tmp # excluding pixels outside the EOO
            pixs <- length(vals)
          }
          if("matrix" %in% class(tmp) | "data.frame" %in% class(tmp)) {
            last <- dim(tmp)[2]
            vals <- tmp[, last] # excluding pixels outside the EOO
            pixs <- length(vals)
          }
          sumario <- matrix(unclass(summary(vals)), nrow = 1, ncol = 6, 
                            dimnames = list("", c("Min","1st_Qu","Median","Mean","3rd_Qu","Max")))
          
          if(output_raster %in% "summary") {
            res <- cbind.data.frame(tax = EOO.poly.name[x],
                                    sumario, # ci,
                                    stringsAsFactors = FALSE)
            
          } else {
            
            #Tabulating pixels
            tmp1 <- table(factor(vals, levels = classes))
            prop <-  matrix(round(100*unclass(tmp1)/pixs, 8), nrow = 1, ncol = length(classes),
                            dimnames = list("", classes))
            
            if(output_raster %in% "prop.table")
              res <- cbind.data.frame(tax = EOO.poly.name[x],
                                      prop,
                                      stringsAsFactors = FALSE)
            
            if(output_raster %in% "area.table")
              area <- round((prop * r.area)/100, 3)
            res <- cbind.data.frame(tax = EOO.poly.name[x],
                                    area,
                                    stringsAsFactors = FALSE)
            
          } 
        }
      } else {
        
        col.names <- c("prop.EOO", "area.EOO", "loss", "recover", "rel.loss", 
                       "rel.recover", "rate.loss", "rate.recover", "hab.quality",   
                       "hab.quality.lo", "hab.quality.hi")
        hab.mat <- matrix(NA, nrow = 1, ncol = length(col.names), 
                          dimnames = list("habitat", col.names))
        res <- cbind.data.frame(tax = EOO.poly.name[x], hab.mat,
                                stringsAsFactors = FALSE)
      } 
    }
      
    res
  }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  result <- do.call(rbind.data.frame, output)
  rownames(result) <- NULL
  return(result)
}

#### GETTING TRUE POLYGON CENTROIDS ####
polygonCenter <- function(x, type = "cm") {
  # how many polygons?
  n.polys <- length(x)
  
  # object that will store the results
  results <- NULL
  for (i in seq_len(n.polys)){
    # cat(i, "\n")
    
    # Getting the polygon with larger area for the cell
    poly.i <- x[i,]
    areas.i <- rgeos::gArea(poly.i, byid = TRUE)
    poly.i.max <- poly.i[which.max(areas.i),]
    
    # Calculating the centroid using rGeos
    centroid <- sp::coordinates(rgeos::gCentroid(poly.i.max))
    
    polys <- attr(poly.i.max, 'polygons')
    polys2 <- attr(polys[[which.max(areas.i)]], 'Polygons')
    coords <- sp::coordinates(polys2[[1]])
    cm <- pracma::poly_center(coords[,1], coords[,2]) # center of mass
    
    # plot(x, border = "grey")
    # plot(poly.i.max, add = TRUE, col = )
    # points(centroid, pch = 19)
    # points(t(cm), pch = 15, col=2)
    
    if (type == "cm")    
      if(is.null(results)) { 
        results <- t(cm) 
      } else { 
        results <- rbind(results, t(cm))
      }
    
    if (type == "centroid")    
      if(is.null(results)) { 
        results <- centroid 
      } else { 
        results <- rbind(results, centroid)
      }
  }
  
  return(results)
}

## Another option taken from: https://gis.stackexchange.com/a/265475/171796

#' find the center of mass / furthest away from any boundary
#' 
#' Takes as input a spatial polygon
#' @param pol One or more polygons as input
#' @param ultimate optional Boolean, TRUE = find polygon furthest away from
#'   centroid. False = ordinary centroid

centroid <- function(pol, ultimate=TRUE, iterations=5, initial_width_step=10){
  if (ultimate){
    new_pol <- pol
    # For every polygon do this:
    for (i in 1:length(pol)){
      width <- -initial_width_step
      area <- rgeos::gArea(pol[i,])
      centr <- pol[i,]
      wasNull <- FALSE
      for (j in 1:iterations){
        if (!wasNull){ # stop when buffer polygon was alread too small
          centr_new <- rgeos::gBuffer(centr,width=width)
          # if the buffer has a negative size:
          substract_width <- width/20
          while (is.null(centr_new)){ #gradually decrease the buffer size until it has positive area
            width <- width-substract_width
            centr_new <- rgeos::gBuffer(centr,width=width)
            wasNull <- TRUE
          }
          # if (!(is.null(centr_new))){
          #   plot(centr_new,add=T)
          # }
          new_area <- rgeos::gArea(centr_new)
          #linear regression:
          slope <- (new_area-area)/width
          #aiming at quarter of the area for the new polygon
          width <- (area/4-area)/slope
          #preparing for next step:
          area <- new_area
          centr<- centr_new
        }
      }
      #take the biggest polygon in case of multiple polygons:
      d <- disaggregate(centr)
      if (length(d)>1){
        biggest_area <- rgeos::gArea(d[1,])
        which_pol <- 1                             
        for (k in 2:length(d)){
          if (rgeos::gArea(d[k,]) > biggest_area){
            biggest_area <- rgeos::gArea(d[k,])
            which_pol <- k
          }
        }
        centr <- d[which_pol,]
      }
      #add to class polygons:
      new_pol@polygons[[i]] <- remove.holes(new_pol@polygons[[i]])
      new_pol@polygons[[i]]@Polygons[[1]]@coords <- centr@polygons[[1]]@Polygons[[1]]@coords
    }
    centroids <- rgeos::gCentroid(new_pol,byid=TRUE)
  }else{
    centroids <- rgeos::gCentroid(pol,byid=TRUE)  
  }  
  return(centroids)
}

#Given an object of class Polygons, returns
#a similar object with no holes

remove.holes <- function(Poly){
  # remove holes
  is.hole <- lapply(Poly@Polygons,function(P)P@hole)
  is.hole <- unlist(is.hole)
  polys <- Poly@Polygons[!is.hole]
  Poly <- Polygons(polys,ID=Poly@ID)
  # remove 'islands'
  max_area <- largest_area(Poly)
  is.sub <- lapply(Poly@Polygons,function(P)P@area<max_area)  
  is.sub <- unlist(is.sub)
  polys <- Poly@Polygons[!is.sub]
  Poly <- Polygons(polys,ID=Poly@ID)
  Poly
}

largest_area <- function(Poly){
  total_polygons <- length(Poly@Polygons)
  max_area <- 0
  for (i in 1:total_polygons){
    max_area <- max(max_area,Poly@Polygons[[i]]@area)
  }
  max_area
}

#' Make transparent theme
transparent=function(size=0){
  
  
  temp=theme(rect= element_rect(fill = 'transparent',size=size),
             panel.background=element_rect(fill = 'transparent'),
             panel.border=element_rect(size=size),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank())
  temp
}

#' Make default palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Make subcolors with main colors
makeSubColor=function(main,no=3){
  result=c()
  for(i in 1:length(main)){
    temp=ztable::gradientColor(main[i],n=no+2)[2:(no+1)]
    result=c(result,temp)
  }
  result
}


#'Draw a PieDonut plot
my.PieDonut=function(data,mapping,
                     start=getOption("PieDonut.start",0),
                     addPieLabel=TRUE,addDonutLabel=TRUE,
                     showRatioDonut=TRUE,showRatioPie=TRUE,
                     ratioByGroup=TRUE,
                     showRatioThreshold=getOption("PieDonut.showRatioThreshold",0.02),
                     labelposition=getOption("PieDonut.labelposition",2),
                     labelpositionThreshold=0.1,
                     r0=getOption("PieDonut.r0",0.3),
                     r1=getOption("PieDonut.r1",1.0),
                     r2=getOption("PieDonut.r2",1.2),
                     explode=NULL,
                     selected=NULL,
                     explodePos=0.1,
                     color="white",
                     pieAlpha=0.8,
                     donutAlpha=1.0,
                     maxx=NULL,
                     showPieName=FALSE,
                     showDonutName=FALSE,
                     title=NULL,
                     pieLabelSize=4,
                     donutLabelSize=3,
                     titlesize=5,
                     explodePie=TRUE,explodeDonut=FALSE,
                     use.label=TRUE,use.labels=TRUE,
                     family=getOption("PieDonut.family",""),
                     labels = NULL,
                     main.colors = NULL,
                     correct.legend.x = NULL,
                     correct.legend.y = NULL,
                     tidy.legend.donut = NULL,
                     tidy.legend.pie = NULL,
                     plot.piedonut = TRUE,
                     ratioAccuracy = c("pie" = 0.1, "donut" = 0.1),
                     return = "both"){
  
  # data = pie.df
  # mapping = aes(category, main.criteria)
  # ratioByGroup=FALSE
  # start=75
  # r0=0
  # r1=1
  # r2=1.3
  # labelposition = 2
  # family = ""
  # main.colors = cores[-1]
  
  (cols=colnames(data))
  if(use.labels) data=moonBook::addLabelDf(data,mapping)
  
  count<-NULL
  
  if("count" %in% names(mapping)) count <- moonBook::getMapping(mapping,"count")
  count
  
  pies<-donuts<-NULL
  (pies = moonBook::getMapping(mapping,"pies"))
  if(is.null(pies)) (pies = moonBook::getMapping(mapping,"pie"))
  if(is.null(pies)) (pies = moonBook::getMapping(mapping,"x"))
  
  (donuts = moonBook::getMapping(mapping,"donuts"))
  if(is.null(donuts)) (donuts = moonBook::getMapping(mapping,"donut"))
  if(is.null(donuts)) (donuts = moonBook::getMapping(mapping,"y"))
  
  if(!is.null(count)){
    
    df<-data %>% group_by(.data[[pies]]) %>%dplyr::summarize(Freq=sum(.data[[count]]))
    df
  } else{
    df=data.frame(table(data[[pies]]))
  }
  colnames(df)[1]=pies
  
  ##Putting the data frame into the order of the colors
  if(!is.null(main.colors)) {
    df <- df[order(match(df$category, names(main.colors))),]
    # df3 <- df3[order(match(df3$category, names(main.colors))),]
  }
  
  df$end=cumsum(df$Freq)
  df$start=dplyr::lag(df$end)
  df$start[1]=0
  total=sum(df$Freq)
  df$start1=df$start*2*pi/total
  df$end1=df$end*2*pi/total
  df$start1=df$start1+start
  df$end1=df$end1+start
  df$focus=0
  if(explodePie) df$focus[explode]=explodePos
  df$mid=(df$start1+df$end1)/2
  df$x=ifelse(df$focus==0,0,df$focus*sin(df$mid))
  df$y=ifelse(df$focus==0,0,df$focus*cos(df$mid))
  df$label=df[[pies]]
  df$ratio=df$Freq/sum(df$Freq)
  
  if(is.null(tidy.legend.pie))
    tidy.legend.pie = rep(TRUE, nrow(df))
  if(showRatioPie) {
    
    # df$label = ifelse(df$ratio >= showRatioThreshold,
    #                 paste0(df$label, "\n(", scales::percent(df$ratio), ")"),
    #                 as.character(df$label))
    df$label = paste0(df$label, " (", scales::percent(df$ratio, accuracy = ratioAccuracy["pie"]), ")")
    df$label[tidy.legend.pie] = gsub(" \\(", "\n(",df$label[tidy.legend.pie])
    
    
  }
  
  if(is.null(correct.legend.x)) {
    correct.legend.x = rep(0, nrow(df))
  }
  if(is.null(correct.legend.y)) {
    correct.legend.y = rep(0, nrow(df))
  }
  
  df$labelx=(r0+r1)/2*sin(df$mid)+df$x + correct.legend.x
  df$labely=(r0+r1)/2*cos(df$mid)+df$y + correct.legend.y
  if(!is.factor(df[[pies]])) df[[pies]]<-factor(df[[pies]])
  df
  
  if(!is.null(main.colors)) {
    mainCol = main.colors
  } else {
    mainCol = gg_color_hue(nrow(df))
  }
  
  df$radius=r1
  df$radius[df$focus!=0]=df$radius[df$focus!=0]+df$focus[df$focus!=0]
  df$hjust=ifelse((df$mid %% (2*pi))>pi,1,0)
  df$vjust=ifelse(((df$mid %% (2*pi)) <(pi/2))|(df$mid %% (2*pi) >(pi*3/2)),0,1)
  df$segx=df$radius*sin(df$mid)
  df$segy=df$radius*cos(df$mid)
  df$segxend=(df$radius+0.05)*sin(df$mid)
  df$segyend=(df$radius+0.05)*cos(df$mid)
  df
  
  if(!is.null(donuts)){
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    
    
    data
    if(!is.null(count)){
      
      df3 <- as.data.frame(data[c(donuts,pies,count)])
      colnames(df3)=c("donut","pie","Freq")
      df3
      df3<-eval(parse(text="complete(df3,donut,pie)"))
      
      # df3<-df3 %>% complete(donut,pie)
      df3$Freq[is.na(df3$Freq)]=0
      if(!is.factor(df3[[1]])) df3[[1]]=factor(df3[[1]])
      if(!is.factor(df3[[2]])) df3[[2]]=factor(df3[[2]])
      
      df3<-df3 %>% arrange(.data$pie,.data$donut)
      a<-df3 %>% spread(.data$pie,value=.data$Freq)
      # a<-df3 %>% spread(pie,value=Freq)
      a=as.data.frame(a)
      a
      rownames(a)=a[[1]]
      a=a[-1]
      a
      colnames(df3)[1:2]=c(donuts,pies)
      
      
      
    } else {
      df3 <- data.frame(table(data[[donuts]],data[[pies]]),stringsAsFactors = FALSE)
      colnames(df3)[1:2] <- c(donuts,pies)
      a = table(data[[donuts]],data[[pies]])
      
      #Putting the data frame into the order of the colors
      if(!is.null(main.colors)) {
        # df <- df[order(match(df$category, names(main.colors))),]
        df3 <- df3[order(match(df3$category, names(main.colors))),]
        a <- a[,order(match(colnames(a), names(main.colors)))]
      }
      
    }
    
    a
    df3
    df3$group = rep(colSums(a),each=nrow(a))
    df3$pie = rep(1:ncol(a),each=nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if(ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group, accuracy = ratioAccuracy["donut"])
    } else {
      df3$ratio <- scales::percent(df3$ratio1, accuracy = ratioAccuracy["donut"])
    }
    df3$end <- cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    
    df3$start1 = df3$start * 2 * pi / total
    df3$end1 = df3$end * 2 * pi / total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1) / 2
    df3$focus = 0
    
    if(!is.null(selected)){
      df3$focus[selected]=explodePos
    } else if(!is.null(explode)) {
      selected=c()
      for(i in 1:length(explode)){
        start=1+nrow(a)*(explode[i]-1)
        selected=c(selected,start:(start+nrow(a)-1))
      }
      selected
      df3$focus[selected]=explodePos
    }
    df3
    df3$x=0
    df3$y=0
    df
    
    if(!is.null(explode)){
      explode
      for(i in 1:length(explode)){
        
        xpos=df$focus[explode[i]]*sin(df$mid[explode[i]])
        ypos=df$focus[explode[i]]*cos(df$mid[explode[i]])
        
        df3$x[df3$pie==explode[i]]=xpos
        df3$y[df3$pie==explode[i]]=ypos
      }
    }
    df3$no=1:nrow(df3)
    df3$label=df3[[donuts]]
    
    if(is.null(tidy.legend.donut))
      tidy.legend.donut <- rep(TRUE, nrow(df3))
    
    if(showRatioDonut) {
      #### CHECK HERE ####
      if(max(nchar(levels(df3$label)))<=2) { 
        df3$label = paste0(df3$label,"(",df3$ratio,")")
      } else { 
        #df3$label = paste0(df3$label,"\n(",df3$ratio,")")
        df3$label = paste0(df3$label," (",df3$ratio,")")
        df3$label[tidy.legend.donut] <- gsub(" \\(","\n(",df3$label[tidy.legend.donut]) 
      }
    }
    df3$label[df3$ratio1==0]=""
    
    # if(labelposition==0)
    df3$label[df3$ratio1<showRatioThreshold]=""
    
    
    df3$hjust=ifelse((df3$mid %% (2*pi))>pi,1,0)
    df3$vjust=ifelse(((df3$mid %% (2*pi)) <(pi/2))|(df3$mid %% (2*pi) >(pi*3/2)),0,1)
    df3$no=factor(df3$no)
    df3
    # str(df3)
    labelposition
    if(labelposition>0){
      df3$radius=r2
      if(explodeDonut) df3$radius[df3$focus!=0]=df3$radius[df3$focus!=0]+df3$focus[df3$focus!=0]
      
      df3$segx=df3$radius*sin(df3$mid)+df3$x
      df3$segy=df3$radius*cos(df3$mid)+df3$y
      df3$segxend=(df3$radius+0.05)*sin(df3$mid)+df3$x
      df3$segyend=(df3$radius+0.05)*cos(df3$mid)+df3$y
      
      if(labelposition==2) df3$radius=(r1+r2)/2
      df3$labelx= (df3$radius)*sin(df3$mid)+df3$x
      df3$labely= (df3$radius)*cos(df3$mid)+df3$y
    } else{
      df3$radius=(r1+r2)/2
      if(explodeDonut) df3$radius[df3$focus!=0]=df3$radius[df3$focus!=0]+df3$focus[df3$focus!=0]
      df3$labelx=df3$radius*sin(df3$mid)+df3$x
      df3$labely=df3$radius*cos(df3$mid)+df3$y
    }
    df3$segx[df3$ratio1==0]=0
    df3$segxend[df3$ratio1==0]=0
    df3$segy[df3$ratio1==0]=0
    df3$segyend[df3$ratio1==0]=0
    if(labelposition==0){
      df3$segx[df3$ratio1<showRatioThreshold]=0
      df3$segxend[df3$ratio1<showRatioThreshold]=0
      df3$segy[df3$ratio1<showRatioThreshold]=0
      df3$segyend[df3$ratio1<showRatioThreshold]=0
    }
    df3
    
    del=which(df3$Freq==0)
    del
    if(length(del)>0) subColor<-subColor[-del]
    subColor
  }
  
  p <- ggplot() + ggforce::theme_no_axes() + coord_fixed()
  
  if(is.null(maxx)) {
    r3=r2+0.3
  } else{
    r3=maxx
  }
  
  p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y",
                                             r0 = as.character(r0), r = as.character(r1),
                                             start="start1",end="end1",
                                             fill = pies),alpha=pieAlpha,color=color, data = df) + 
    transparent()+
    scale_fill_manual(values=mainCol)+
    xlim(r3*c(-1.2,1.2)) + ylim(r3*c(-1,1)) + guides(fill=FALSE) +
    theme_void()
  
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy",
                                       xend="segxend",yend="segyend"),data=df)+
      geom_text(aes_string(x="segxend",y="segyend",label="label",hjust="hjust",vjust="vjust"),size=pieLabelSize,data=df,family=family)
    
  } else 
    if ((labelposition == 2) & (is.null(donuts))) {
      p1<-p1+ geom_segment(aes_string(x="segx",y="segy",
                                      xend="segxend",yend="segyend"),data=df[df$ratio<labelpositionThreshold,])+
        geom_text(aes_string(x="segxend",y="segyend",label="label",hjust="hjust",vjust="vjust"),size=pieLabelSize,data=df[df$ratio<labelpositionThreshold,],family=family)+
        geom_text(aes_string(x="labelx",y="labely",label="label"),size=pieLabelSize,data=df[df$ratio>=labelpositionThreshold,],family=family)
      
      
    } else{
      p1 <- p1 + geom_text(
        aes_string(x = "labelx", y = "labely", label = "label"),
        size = pieLabelSize,
        data = df,
        family = family
      )
    }
  
  if(showPieName) p1<-p1+annotate("text",x=0,y=0,label=pies,size=titlesize,family=family)
  
  p1 <- p1+theme(text=element_text(family=family))
  
  if(!is.null(donuts)){
    
    # donutAlpha=1.0;color="white"
    # explodeDonut=FALSE
    if(explodeDonut) {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r1),
                                                 r = as.character(r2), start="start1",end="end1",
                                                 fill = "no",explode="focus"),alpha=donutAlpha,color=color,
                                      data = df3)
    } else{
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r1),
                                                 r = as.character(r2), start="start1",end="end1",
                                                 fill = "no"),alpha=donutAlpha,color=color,
                                      data = df3)
    }
    
    p3 <- p3 + transparent()+
      scale_fill_manual(values=subColor)+
      xlim(r3*c(-1,1))+ylim(r3*c(-1,1))+guides(fill=FALSE)+
      theme_void()
    
    p3
    
    if(labelposition==1){
      p3<-p3+ geom_segment(aes_string(x="segx",y="segy",
                                      xend="segxend",yend="segyend"),data=df3)+
        geom_text(aes_string(x="segxend",y="segyend",
                             label="label",hjust="hjust",vjust="vjust"),size=donutLabelSize,data=df3,family=family)
    } else 
      if(labelposition==0){
        p3<-p3+geom_text(aes_string(x="labelx",y="labely",
                                    label="label"),size=donutLabelSize,data=df3,family=family)
      } else{
        p3 <- p3 + geom_segment(aes_string(x="segx",y="segy",
                                           xend="segxend",yend="segyend"),data=df3[df3$ratio1<labelpositionThreshold,])+
          geom_text(aes_string(x="segxend",y="segyend",
                               label="label",hjust="hjust",vjust="vjust"),size=donutLabelSize,data=df3[df3$ratio1<labelpositionThreshold,],family=family)+
          geom_text(aes_string(x="labelx",y="labely",
                               label="label"),size=donutLabelSize,data=df3[df3$ratio1>=labelpositionThreshold,],family=family)
        
      }
    
    if(!is.null(title)) 
      p3<-p3+annotate("text",x=0,y=r3,label=title,size=titlesize,family=family)
    if(showDonutName) 
      p3<-p3+annotate("text",x=(-1)*r3,y=r3,label=donuts,hjust=0,size=titlesize,family=family)
    
    p3 <- p3 + theme(text=element_text(family=family))# +
    #theme(axis.ticks.length = unit(.1, "cm"))
    #theme(axis.ticks = element_blank())
    #p3
    # grid::grid.newpage()
    # print(p1, vp = grid::viewport(height = 1, width = 1))
    # print(p3, vp = grid::viewport(height = 1, width = 1))
  } # else {
  #   p1
  # }
  
  if (plot.piedonut) {
    grid::grid.newpage()
    print(p1, vp = grid::viewport(height = 1, width = 1))
    print(p3, vp = grid::viewport(height = 1, width = 1))
  }
  
  if(return == "both"){
    return(list(p1, p3))
  } else{
    if(return == "pie") return(p1)
    if(return == "donut") return(p3)
  }
  
}

## function to obtain threat from habitat loss from empirical spline fits to AF derived relationships
get_threat <- function(habitat_loss = NULL, spp_richness = NULL,
                       endemism_ratio = NULL, all = NULL, 
                       VU = NULL, EN = NULL, CR = NULL,
                       correct.cats = TRUE) {
  # getting the proportion of loss
  if (any(habitat_loss[!is.na(habitat_loss)] > 1))
    habitat_loss <- habitat_loss/100
  
  perdas <- habitat_loss[!is.na(habitat_loss)]
  
  # getting the overall threat proportion
  low <- predict(all[[1]], x=perdas)$y
  med <- predict(all[[2]], x=perdas)$y
  high <- predict(all[[3]], x=perdas)$y
  threat <- cbind(low, med, high)
  
  low <- predict(VU[[1]], x=perdas)$y
  med <- predict(VU[[2]], x=perdas)$y
  high <- predict(VU[[3]], x=perdas)$y
  threat.VU <- cbind(low, med, high)
  
  low <- predict(EN[[1]], x=perdas)$y
  med <- predict(EN[[2]], x=perdas)$y
  high <- predict(EN[[3]], x=perdas)$y
  threat.EN <- cbind(low, med, high)
  
  low <- predict(CR[[1]], x=perdas)$y
  med <- predict(CR[[2]], x=perdas)$y
  high <- predict(CR[[3]], x=perdas)$y
  threat.CR <- cbind(low, med, high)
  
  
  # Averaging estimates of threat categories based on overall threat
  if (correct.cats) {
    low <- apply(cbind(threat.VU[,"low"], threat.EN[,"low"], threat.CR[,"low"]), 1, sum)
    med <- apply(cbind(threat.VU[,"med"], threat.EN[,"med"], threat.CR[,"med"]), 1, sum)
    high <- apply(cbind(threat.VU[,"high"], threat.EN[,"high"], threat.CR[,"high"]), 1, sum)
    threat.cats <- cbind(low, med, high)
    
    low.ratio <- threat[,"low"]/threat.cats[,"low"]
    threat.VU[,"low"] <- threat.VU[,"low"] * low.ratio
    threat.EN[,"low"] <- threat.EN[,"low"] * low.ratio
    threat.CR[,"low"] <- threat.CR[,"low"] * low.ratio
    
    med.ratio <- threat[,"med"]/threat.cats[,"med"]
    threat.VU[,"med"] <- threat.VU[,"med"] * med.ratio
    threat.EN[,"med"] <- threat.EN[,"med"] * med.ratio
    threat.CR[,"med"] <- threat.CR[,"med"] * med.ratio
    
    high.ratio <- threat[,"med"]/threat.cats[,"med"]
    threat.VU[,"high"] <- threat.VU[,"high"] * high.ratio
    threat.EN[,"high"] <- threat.EN[,"high"] * high.ratio
    threat.CR[,"high"] <- threat.CR[,"high"] * high.ratio
    
  }
  
  # Making sure that overall proportion is not above 100%
  dados <- cbind(VU = threat.VU[,1], EN = threat.EN[,1], CR = threat.CR[,1])
  dados <- round(100 * cbind(dados, LC = (1 - apply(dados, 1, sum))), 0)
  dados[apply(dados,1,sum) != 100, 4] <- 
    dados[apply(dados,1,sum) != 100, 4] + 
    (100 - apply(dados,1,sum)[apply(dados,1,sum) != 100])
  dados[,4][dados[,4] < 0 ] <- 0
  # Calculating the RLI
  high.RLI <- 
    sapply(seq_along(perdas), 
           function(i) red::rli(rep(colnames(dados), times = dados[i,])))
  dados <- cbind(VU = threat.VU[,2], EN = threat.EN[,2], CR = threat.CR[,2])
  dados <- round(100 * cbind(dados, LC = 1 - apply(dados, 1, sum)), 0)
  dados[apply(dados,1,sum) != 100, 4] <- 
    dados[apply(dados,1,sum) != 100, 4] + 
    (100 - apply(dados,1,sum)[apply(dados,1,sum) != 100])
  dados[,4][dados[,4] < 0 ] <- 0
  med.RLI <- 
    sapply(seq_along(perdas), 
           function(i) red::rli(rep(colnames(dados), times = dados[i,])))
  dados <- cbind(VU = threat.VU[,3], EN = threat.EN[,3], CR = threat.CR[,3])
  dados <- round(100 * cbind(dados, LC = 1 - apply(dados, 1, sum)), 0)
  dados[apply(dados, 1, sum) != 100, 4] <-
    dados[apply(dados, 1, sum) != 100, 4] +
    (100 - apply(dados,1,sum)[apply(dados,1,sum) != 100])
  dados[,4][dados[,4] < 0 ] <- 0
  low.RLI <- 
    sapply(seq_along(perdas), 
           function(i) red::rli(rep(colnames(dados), times = dados[i,])))
  
  # combining the results
  threat <- cbind(threat, threat.VU, threat.EN, threat.CR)
  colnames(threat) <- paste0(colnames(threat), 
                             rep(c("_all","_VU","_EN","_CR"), each = 3))
  RLI <- cbind(low_rli = low.RLI, median_rli = med.RLI, high_rli = high.RLI)
  
  # getting the number of threatened species
  if (!is.null(spp_richness)) {
    threat_spp <- round(spp_richness * threat, 1)
    colnames(threat_spp) <- paste0(colnames(threat_spp),"_spp")
    
    # getting the number of threatened endemic species
    if (!is.null(endemism_ratio)) {
      if (any(endemism_ratio[!is.na(endemism_ratio)] > 1))
        endemism_ratio <- endemism_ratio/100
      
      threat_end_spp <- threat_spp * endemism_ratio
      colnames(threat_end_spp) <- gsub("_spp", "end_spp", colnames(threat_end_spp))
      threat <- cbind(round(100 * threat, 1), threat_spp, threat_end_spp)
      
    } else {
      threat <- cbind(round(100 * threat, 1), threat_spp)
    }
  } else {
    threat <- round(100 * threat, 1)
  }
  
  threat <- cbind(threat, RLI)
  return(threat)
}

