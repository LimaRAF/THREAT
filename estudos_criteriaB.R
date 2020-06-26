
##Convex vs alpha-hull
x <- runif(20,-24,-19)
y <- runif(20,-44,-39)


# Delaunay triangulation
dv <- delvor(x,y)
plot(dv, wlines="del", col=1:3, xlim=c(-25,-18))

## convex hull
tmp <- grDevices::chull(x,y) # convex order
hull <- cbind(x[tmp], y[tmp]) # extract the row numbers of the boundary points, in convex order
hull.polygon <- SpatialPolygons(list(Polygons(list(Polygon(hull[,1:2])), ID = 1)))
# Add projection??: SpatialPolygons(proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),list(Polygons(list(Polygon(hull[,1:2])), ID = 1)))
plot(hull.polygon, lwd=2, add=TRUE)
points(x,y,pch=19)
#text(hull[,1],hull[,2],labels = tmp, pos=1:4)
# Save polygon??: writeOGR(hull.polygon,"chull", layer="chull", driver="ESRI Shapefile")

## alpha-hull
alfa <- 2
tmp <- ahull(x, y, alpha = alfa)
tmp1 = tmp$arcs[,"end1"] # convex order   
alpha <- cbind(x[tmp1], y[tmp1]) # extract the row numbers of the boundary points, in convex order
alpha.polygon <- SpatialPolygons(list(Polygons(list(Polygon(alpha[,1:2])), ID = 1)))
plot(alpha.polygon, add=TRUE, border="blue", lwd=2)
plot(tmp, add=TRUE, col=c("green",rep("black",5)), lwd=2)

