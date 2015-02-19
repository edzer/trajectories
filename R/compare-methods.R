## get distances between 2 tracks for each point in time where they overlap
setGeneric(
  name = "compare",
  def = function(track1, track2, ...) standardGeneric("compare")
)

compare.track <- function(track1, track2) {
  require(xts)
  if (!(first(track1@endTime) < last(track2@endTime) && first(track2@endTime) < last(track1@endTime)))
      stop("Time itervals don't overlap!")
  if (!identicalCRS(track1, track2))
      stop("CRS are not identical!")
  track1.df <- cbind(id = rep(1, dim(track1)), 
                     as.data.frame(track1)[c(coordnames(track1), "time")])
  track2.df <- cbind(id = rep(2, dim(track2)), 
                     as.data.frame(track2)[c(coordnames(track2), "time")])
  track1.df$iv <- findInterval(track1.df$time, track2.df$time) # intervals timestamps fall in
  track2.df$iv <- findInterval(track2.df$time, track1.df$time)
  # find points and create data frame
  dtrack <- rbind(findPoints(track1.df, track2.df), findPoints(track2.df, track1.df))
  dtrack <- dtrack[!is.na(dtrack$p.x),] # remove points falling outside the time intervals
  dtrack <- dtrack[order(dtrack$time),] # sort by timestamp
  dtrack$dist <- NA
  for (i in 1:nrow(dtrack)) { # distance at timestamp
    p1 <- SpatialPoints(cbind(dtrack[i,2], dtrack[i,3]), CRS(proj4string(track1)))
    p2 <- SpatialPoints(cbind(dtrack$p.x[i], dtrack$p.y[i]), CRS(proj4string(track1)))
    dtrack$dist[i] <- spDistsN1(p1,p2)
  }
  dtrack
}

setMethod("compare", signature("Track"), compare.track)


## finds corresponding points for track1 on track2
findPoints <- function(tr1, tr2) {
  tr1$p.x <- NA # projected point
  tr1$p.y <- NA
  for (i in 1:nrow(tr1)) {
    if (!tr1[i,]$iv == 0 && !tr1[i,]$iv == nrow(tr2)) {
      iv <- tr1[i,]$iv
      tdiff1 <- tr1$time[i] - tr2$time[iv] # diff between timestamp and start of the interval it falls in
      tdiff2 <- tr2$time[iv+1] - tr2$time[iv] # diff between timestamps (calculated here because it often varies)
      ratio <- as.numeric(tdiff1)/as.numeric(tdiff2) 
      x1 <- tr2[iv,2] # segment coordinates
      y1 <- tr2[iv,3]
      x2 <- tr2[iv+1,2]
      y2 <- tr2[iv+1,3]
      tr1$p.x[i] <- x1 + ratio * (x2 - x1) #find point 
      tr1$p.y[i] <- y1 + ratio * (y2 - y1) 
    }
  }
  tr1
}



## calculates frechet distance
setGeneric(
  name = "frechetDist",
  def = function(track1, track2, ...) standardGeneric("frechetDist")
)

frechetDist.track <- function(track1, track2) {
  if (!identicalCRS(track1, track2))
    stop("CRS are not identical!")
  dists <- spDists(track1@sp, track2@sp) #dists between all points
  dists[,1] <- Reduce(max, dists[,1], accumulate=TRUE) # cases where one of the trajectories is a point 
  dists[1,] <- Reduce(max, dists[1,], accumulate=TRUE)
  for (i in 2:nrow(dists)) { # build rest of frechet distance matrix
    for (j in 2:ncol(dists)) {
      dists[i,j] <- max(dists[i,j], min(dists[i-1,j], dists[i-1,j-1], dists[i,j-1]))
    }
  }
  last(last(dists))
}

setMethod("frechetDist", signature("Track"), frechetDist.track)



## downsamples a track to the length of another one
setGeneric(
  name = "downsample",
  def = function(track1, track2, ...) standardGeneric("downsample")
)

# track1: track that will be downsampled
# track2: to the dimension of track2
downsample.track <- function(track1, track2) {
  if (!identicalCRS(track1, track2))
    stop("CRS are not identical!")
  if (dim(track1) == dim(track2))
    stop("Dimensions are euqal!")
  tr <- track1
  xy <- coordinates(track1)
  time <- index(track1@time)
  crs <- CRS (proj4string(track1))
  while (dim(track1) > dim(track2)) {
    d1 <- track1$distance # distances
    n <- length (d1) - 1 # number of segments between every second point
    xy1 <- cbind (head (xy, n), tail (xy, n))    
    d2.long <- head (d1, n) + tail (d1, n)
    xy.new <- list ()
    for (i in 1:n) xy.new [[i]] <- rbind (head (xy, n) [i,], tail (xy, n) [i,])
    d2.short <- sapply (xy.new, function (x) spDists (x, longlat=TRUE)[1,2])
    remove <- which.min (d2.long - d2.short) + 1
    xy <- xy[- remove,]
    time <- time[- remove]
    stidf <- STIDF(SpatialPoints (xy, crs), time, data.frame(extraDat=rnorm(n)))
    tr  <- Track (stidf)
  }
  tr
}

setMethod("downsample", signature("Track"), downsample.track)