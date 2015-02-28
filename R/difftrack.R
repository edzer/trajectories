setClass("difftrack",
         slots=c(track1 ="Track", track2 = "Track", 
                 conns1 = "SpatialLinesDataFrame", conns2 = "SpatialLinesDataFrame"),
)

## plots a difftrack
plot.difftrack <- function(x, y, ..., axes = TRUE) {  
  plot(x@track1@sp, col="red", ..., axes = axes)
  points(x@track2@sp, col="blue")
  lines(x@conns1)
  lines(x@conns2)
}

setMethod("plot", "difftrack", plot.difftrack)



## stcube for difftrack
stcube.difftrack <- function(x, showMap = FALSE, mapType = "osm", ...) {
  stcube(x@track1, col = "red", showMap = showMap, mapType = mypType, ...)
  stcube(x@track2, col = "blue", add = TRUE, ...)
  
  lines1 <- x@conns1@lines
  lines2 <- x@conns2@lines
  time1 <- x@conns1@data$time
  time1 <- time1 - min(index(x@track1@time))
  time2 <- x@conns2@data$time
  time2 <- time2 - min(index(x@track2@time))
  
  sapply(lines1, function(l) {
    coords <- coordinates(l)
    id <- as.numeric(l@ID)
    z1 <- time1[id]
    z2 <- time2[id]
    rgl::lines3d(coords[[1]][,1], coords[[1]][,2], c(z2, z1), col = "red")    
  })
  sapply(lines2, function(l) {
    coords <- coordinates(l)
    id <- as.numeric(l@ID)
    z1 <- time2[id]
    z2 <- time1[id]
    rgl::lines3d(coords[[1]][,1], coords[[1]][,2], c(z1, z2), col = "blue")    
  })
}

setMethod("stcube", "difftrack", stcube.difftrack)