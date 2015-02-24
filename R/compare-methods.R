setClass("difftrack",
         slots=c(track1 ="Track", track2 = "Track", sl = "SpatialLines", data = "data.frame"),
)


## get distances between 2 tracks for each point in time where they overlap
## extend each track with these points
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
  crs <- CRS(proj4string(track1))
  track1.df <- cbind(as.data.frame(track1)[c(coordnames(track1), "time")])
  track2.df <- cbind(as.data.frame(track2)[c(coordnames(track2), "time")])  
  # intervals timestamps fall in
  ivs1 <- findInterval(track1.df$time, track2.df$time) 
  ivs2 <- findInterval(track2.df$time, track1.df$time)
  # find points and create new extended data frames
  newTrack1.df <- findPoints(track2.df, track1.df, ivs2)
  newTrack2.df <- findPoints(track1.df, track2.df, ivs1)
  # equal timestamps
  idx <- timeMatch(newTrack1.df$time, newTrack2.df$time)
  dists <- rep(NA, length(idx))
  lines <- list()
  for (i in 1:length(idx)) { # distance at timestamp
    if (!is.na(idx[i])) {
      coords1 <- cbind(newTrack1.df$x[i], newTrack1.df$y[i])
      coords2 <- cbind(newTrack2.df$x[idx[i]], newTrack2.df$y[idx[i]])
      p1 <- SpatialPoints(coords1, crs)
      p2 <- SpatialPoints(coords2, crs)
      dists[i] <- spDistsN1(p1,p2) # dists between tracks at matching timestamps
      lines <- c(lines, Line(rbind(coords1, coords2))) # corresponding line
    }
  }
  # extended tracks
  newTrack1 <- STIDF(SpatialPoints(cbind(newTrack1.df$x, newTrack1.df$y), crs), newTrack1.df$time, data.frame(1:nrow(newTrack1.df)))
  newTrack2 <- STIDF(SpatialPoints(cbind(newTrack2.df$x, newTrack2.df$y), crs), newTrack2.df$time, data.frame(1:nrow(newTrack2.df)))
  newTrack1 <- Track(newTrack1)
  newTrack2 <- Track(newTrack2)
  sl <- SpatialLines(list(Lines(lines, ID = 1)), crs) #lines between tracks
  new("difftrack", track1 = newTrack1, track2 = newTrack2, sl = sl, data = data.frame(idx, dists))
}

setMethod("compare", signature("Track"), compare.track)


## finds corresponding points for track1 on track2
findPoints <- function(tr1, tr2, ivs) {
  x <- tr2[,1]
  y <- tr2[,2]
  time <- tr2[,3]
  for (i in 1:nrow(tr1)) {
    if (!ivs[i] == 0 && !ivs[i] == nrow(tr2)) {
      iv <- ivs[i]
      tdiff1 <- tr1$time[i] - tr2$time[iv] # diff between timestamp and start of the interval it falls in
      tdiff2 <- tr2$time[iv+1] - tr2$time[iv] # diff between timestamps (calculated here because it often varies)
      ratio <- as.numeric(tdiff1)/as.numeric(tdiff2) 
      x1 <- tr2[iv,1] # segment coordinates
      y1 <- tr2[iv,2]
      x2 <- tr2[iv+1,1]
      y2 <- tr2[iv+1,2]
      x <- c(x, x1 + ratio * (x2 - x1)) #find point 
      y <- c(y, y1 + ratio * (y2 - y1))
      time <- c(time, tr1$time[i])
    }
  }
  newTrack <- data.frame(x, y, time)
  newTrack <- newTrack[!duplicated(newTrack),] # remove duplicates
  newTrack <- newTrack[order(newTrack$time),] # sort by timestamp
  newTrack
}


## plots a difftrack
plot.difftrack <- function(x, ...) {  
  plot(x@track1@sp, col="red")
  points(x@track2@sp, col="blue")
  lines(x@sl)
}

setMethod("plot", "difftrack", plot.difftrack)


## calculates frechet distance
setGeneric(
  name = "frechetDist",
  def = function(track1, track2, ...) standardGeneric("frechetDist")
)

frechetDist.track <- function(track1, track2) {
  if (!identicalCRS(track1, track2))
    stop("CRS are not identical!")
  dists <- spDists(track1@sp, track2@sp) #dists between all points
  dists[,1] <- cummax(dists[,1]) # cases where one of the trajectories is a point 
  dists[1,] <- cummax(dists[1,])
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
  if(!identicalCRS(track1, track2))
    stop("CRS are not identical!")
  if(dim(track1) == dim(track2))
    stop("Dimensions are euqal!")
  tr <- track1
  xy <- coordinates(track1)
  time <- index(track1@time)
  crs <- CRS (proj4string(track1))
  while(dim(track1) > dim(track2)) {
    d1 <- track1$distance # distances
    n <- length(d1) - 1 # number of segments between every second point
    xy1 <- cbind(head (xy, n), tail (xy, n))    
    d2.long <- head(d1, n) + tail(d1, n)
    xy.new <- list()
    for(i in 1:n) xy.new[[i]] <- rbind(head(xy, n)[i,], tail(xy, n)[i,])
    d2.short <- sapply (xy.new, function(x) spDists(x, longlat=TRUE)[1,2])
    remove <- which.min(d2.long - d2.short) + 1
    xy <- xy[- remove,]
    time <- time[- remove]
    stidf <- STIDF(SpatialPoints (xy, crs), time, data.frame(extraDat=rnorm(n)))
    tr  <- Track (stidf)
  }
  tr
}

setMethod("downsample", signature("Track"), downsample.track)