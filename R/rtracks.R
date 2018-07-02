rTrack <- function (n = 100, origin = c(0, 0), start = as.POSIXct("1970-01-01"), 
                    ar = 0.8, step = 60, sd0 = 1,bbox=bbox, transform=FALSE,nrandom=FALSE, ...){
  
  if(nrandom)  repeat{n <- rpois(1,n);if(!n==0) break()}
  if (missing(bbox) & transform) {
    xo <- runif(1)
    yo <- runif(1)
    origin <- c(xo,yo)
  }  
  if (!missing(bbox) & transform) {
    xo <- runif(1,bbox[1,1],bbox[1,2])
    yo <- runif(1,bbox[2,1],bbox[2,2])
    origin <- c(xo,yo)
  }
  if (length(ar) == 1 && ar == 0) 
    xy = cbind(x=cumsum(rnorm(n, sd = sd0)) + origin[1], y=cumsum(rnorm(n, 
                                                                    sd = sd0)) + origin[2])
  else {xy = cbind(x=origin[1] + cumsum(as.vector(arima.sim(list(ar = ar), 
                                                          n, sd = sd0, ...))), 
                   y=origin[2] + cumsum(as.vector(arima.sim(list(ar = ar), 
                                                          
                                                          n, sd = sd0, ...))))}
  if(transform) {
    if(missing(bbox))  bbox <- matrix(c(0,1,0,1),nrow = 2,byrow = T); colnames(bbox) <- c("min","max");rownames(bbox) <- c("x","y")
    
    xr <- max(xy[,1])-min(xy[,1])
    yr <- max(xy[,2])-min(xy[,2])
    
    xt <- (xy[,1]-min(xy[,1]))/xr
    yt <- (xy[,2]-min(xy[,2]))/yr
    
    xy <- cbind(xt,yt)
    xy <- cbind(x=xy[,1]*bbox[1,2],y=xy[,2]*bbox[2,2])
  }
  
  T = start + 0:(n - 1) * step
  sti = STI(SpatialPoints(xy), T)
  out <- Track(sti)
  if (transform) out@sp@bbox <- bbox
  return(out)
}



rTracks <- function (m = 20, start = as.POSIXct("1970-01-01"), delta = 7200, 
                     sd1 = 0, origin = c(0, 0), ...) 
  Tracks(lapply(0:(m - 1) * delta, function(x) rTrack(start = start + 
                                                        x, origin = origin + rnorm(2, sd = sd1), ...)))

rTracksCollection <- function (p = 10, sd2 = 0, ...) 
  TracksCollection(lapply(1:p, function(x) rTracks(origin = rnorm(2, 
                                                                  sd = sd2), ...)))


print.Track <- function(x,...){
  X = x
  if (class(X@sp)=="SpatialPoints") {
    cat("An object of class Track \n");
    cat(paste0(nrow(as.data.frame(X@sp)), "points"),"\n");
    }
  if (class(X@sp)=="SpatialLines") {
    cat("A generalized object of class Track \n");
    cat(paste0(length(X@sp@lines), "lines"),"\n"); 
  }
  cat(paste0("bbox:"),"\n");
  print(X@sp@bbox);
  cat(paste0("Time period: [",range(X@endTime)[1],", ", range(X@endTime)[2],"]"))
}

print.Tracks <- function(X){
  cat("An object of class Tracks" ,"\n");
  cat(paste0(length(X@tracks)), "tracks followed by a single object")
}

print.TracksCollection <- function(X){
  cat("An object of class TracksCollection" ,"\n");
  cat(paste0(length(X@tracksCollection))
      , "collection of tracks followed by", paste0(length(X@tracksCollection)), " object")
}
