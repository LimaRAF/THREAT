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
