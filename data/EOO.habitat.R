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

  proj_crs <- proj_crs(proj_type = proj_type)
  
  proj_hab <- raster::crs(hab.map)
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
  
  if (is.null(hab.class))
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
      #.combine = 'c',
      #.packages=c("raster","sf"),
      .options.snow = opts
    ) %d% {
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)

    if(any(grepl("sf", class(hab.map)))) {

      crop <- sf::st_intersection(hab.map, EOO.poly.crop[x, ]) # cropping raster to the extent of the polygon
      crop.proj <- sf::st_transform(crop, crs = proj_crs)
      area.hab <- as.numeric(sum(sf::st_area(crop), na.rm = TRUE))/1000000
      tmp <- c(as.numeric(sf::st_area(EOO.poly.crop[x, ]))/1000000 - area.hab, area.hab)
      hab.mat <- matrix(c(tmp/EOO.poly.area[x], tmp), ncol = 2, nrow = 2, byrow = FALSE)
      row.names(hab.mat) <- c("non_habitat", "habitat")
      colnames(hab.mat) <- c("prop.EOO", "area.EOO")
      res <- cbind.data.frame(name = EOO.poly.name[x], hab.mat,
                              stringsAsFactors = FALSE)["habitat",]


      if(export_shp) {

        ## Tried to combine the two sf objects but could not make it with the df info
        # toto <- c(st_geometry(crop), st_geometry(EOO.poly.crop[x, ]))
        # merge(toto, crop)
        # toto <- sf::st_multipolygon(st_geometry(eoo1), st_geometry(crop1))
        # crop1 <- st_cast(crop, "POLYGON")
        # eoo1 <- st_cast(EOO.poly.crop[x,], "POLYGON")

        res <- list(res, crop)

        if (plot) {
          par(mar = c(3,3,2,2), las=1, tcl=-0.25, mgp=c(2.5,0.5,0))
          plot(sf::st_geometry(EOO.poly.crop[x, ]),
               main = EOO.poly.name[x], bg=0)
          plot(sf::st_geometry(crop), add = TRUE, col = 1)
        }

      }

    }

    if(grepl("Raster", class(hab.map))) {
      
      ext <- raster::extent(sf::st_bbox(EOO.poly.crop[x, ])[c(1,3,2,4)])
      crop <- raster::crop(hab.map, ext) # cropping raster to the extent of the polygon
      mask <- raster::mask(crop, EOO.poly.crop[x, ])
      tmp  <- raster::getValues(mask) # much faster than raster::extract
      tmp  <- tmp[!is.na(tmp[,1]),] # excluding pixels outside the EOO
      #gc()

      if (!is.null(hab.class)) {

      # Getting habitat quantity (land-use maps)
        mask.dim <- dim(tmp)
        pixs <- mask.dim[1]
        last <- mask.dim[2]
        hab.mat <- as.matrix(table(factor(tmp[, last] %in% hab.class, levels = c(FALSE, TRUE))))
        row.names(hab.mat) <- c("non_habitat", "habitat")
        hab.mat <- cbind(hab.mat, 100 * hab.mat/pixs)
        hab.mat <- cbind(hab.mat, (hab.mat[,2]* EOO.poly.area[x])/100)
        colnames(hab.mat) <- c("n.pixs", "prop.EOO", "area.EOO")

        if (class(crop) == "RasterBrick") {

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

          res <- cbind.data.frame(name = EOO.poly.name[x], hab.mat[,-1],
                                  stringsAsFactors = FALSE)["habitat",]

          # location of decline and recover
          if (export_shp) {
            #mask <- raster::mask(crop, EOO.poly.crop[x,], cellnumbers = TRUE)
            loc.loss <- data.frame(
              cell = 1:length(mask[[1]]@data@values),
              raster::coordinates(mask[[1]]),
              classes = NA_integer_
            )
            #loc.loss$classes[loc.loss$cell %in% tmp[, 1]] <- transit
            loc.loss$classes[!is.na(mask[[1]]@data@values)] <- transit
            r <- raster::rasterFromXYZ(loc.loss[!is.na(loc.loss$classes), 2:4],
                                       crs = raster::crs(mask[[1]]))

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
        vals <- tmp[, last]
        sumario <- summary(vals)
        mod <- stats::glm(vals ~ 1, family = "poisson" )
        ci <- suppressMessages(round(exp(confint(mod)),4))
        
        #Tabulating pixels
        tmp1 <- table(factor(tmp[, last], levels = classes))
        prop <- round(100*tmp1/pixs, 8)
        area <- round((prop * EOO.poly.area[x])/100, 3)

        if(output_raster %in% "summary")
          res <- c(name = EOO.poly.name[x],
                   sumario,
                   ci)

        if(output_raster %in% "prop.table")
          res <- c(name = EOO.poly.name[x],
                   prop)

        if(output_raster %in% "area.table")
          res <- c(name = EOO.poly.name[x],
                   area)

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
  