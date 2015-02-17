setGeneric(
  name = "compare",
  def = function(track1, ...) standardGeneric("compare")
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
  dtrack <- rbind(findPoints(track1.df, track2.df), findPoints(track2.df, track1.df))
  dtrack <- dtrack[!is.na(dtrack$p.x),]
  dtrack <- dtrack[order(dtrack$time),]
  dtrack$dist <- NA
  for (i in 1:nrow(dtrack)) { # distance at timestamp
    p1 <- SpatialPoints(cbind(dtrack[i,2], dtrack[i,3]), CRS(proj4string(track1)))
    p2 <- SpatialPoints(cbind(dtrack$p.x[i], dtrack$p.y[i]), CRS(proj4string(track1)))
    dtrack$dist[i] <- spDistsN1(p1,p2)
  }
  dtrack
}

setMethod("compare", signature("Track"), compare.track)



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