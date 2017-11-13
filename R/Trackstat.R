# function as.Track  accepts a dataframe (let say X) with three columns "xcoor", "ycoor" and time (or any class that can be converted to a data.frame)
# and converts it to an object of class Track. It can also accepts covariates for the corresponding 
# locations, covariates must be a dataframe with some columns and length of each column is equal
# to length of each column in X.
as.Track <- function(X,covariate){
  stopifnot(nrow(X)>0)
  colnames(X) <- c("xcoor","ycoor","time")
  if(!is.data.frame(X)) X <- as.data.frame(X)
  sp <- cbind(x=X$xcoor,y=X$ycoor)
  sp <- SpatialPoints(sp)
  t <- as.POSIXct(paste(X$date,X$time))
  if(missing(covariate)) covariate <- data.frame(d=rep(NA,length(X$xcoor)))
  
  return(Track(STIDF(sp,time = t,data =covariate)))
}


# function reTrack accepts X as an object of class Track. Output is a reconstructed Track (an object of class Track), based on "timestamp".
# It only returns the interpolated points.
reTrack <- function(X,at=c("track","dfrm"),timestamp=timestamp){
  tsq <- tsqTracks(X,timestamp = timestamp)
  X <-  cbind(as.data.frame(X)[c(coordnames(X), "time")])
  Xrange <- range(X$time)
  xnew <- c()
  ynew <- c()
  
  x <- X$x
  y <- X$y
  
  time <- tsq[tsq<Xrange[2] & tsq>Xrange[1]]
  ivs <- findInterval(time,X$time)
  
  for (i in 1:length(ivs)) {
    if (!ivs[i] == 0 && !ivs[i] == nrow(X)) {
      
      iv <- ivs[i]
      tdiff1 <- difftime(time[i],X$time[iv],units = "sec") # diff between timestamp and start of the interval it falls in
      tdiff2 <- difftime(X$time[iv+1],X$time[iv],units = units(tdiff1)) # diff between timestamps (calculated here because it often varies)
      ratio <- as.numeric(tdiff1)/as.numeric(tdiff2)
      x1 <- X[iv,1] # segment coordinates
      y1 <- X[iv,2]
      x2 <- X[iv+1,1]
      y2 <- X[iv+1,2]
      xnew <- c(xnew, x1 + ratio * (x2 - x1)) #find point
      ynew <- c(ynew, y1 + ratio * (y2 - y1))
    }
  }
  
  newTrack <- data.frame(xnew, ynew, time)
  newTrack <- newTrack[!duplicated(newTrack),] # remove duplicates
  newTrack <- newTrack[order(newTrack$time),] # sort by timestamp
  colnames(newTrack) <- c("xcoor","ycoor","time")
  if (at=="dfrm") return(newTrack)
  return(as.Track(newTrack))
}

# rngTrack returns the timerange of an object of class Track
rngTrack <- function(X){
  Y <- cbind(as.data.frame(X)[c(coordnames(X), "time")])
  return(range(Y$time)) 
}

# tsqtracks returns a sequance of time based on a list of tracks (or a single object of class Track) and an argument timestamp
tsqTracks <- function(X,timestamp){
  
  if (!is.list(X)) timerange <- rngTrack(X)
  else timerange <- lapply(X, rngTrack)
  
  Trackrg <- range(timerange)
  class(Trackrg) <- c('POSIXt','POSIXct')
  # a seq from the range has been created every timestamp
  timeseq <- seq(from=as.POSIXct(strftime(Trackrg[1])),to=as.POSIXct(strftime(Trackrg[2])),by = timestamp)
  
  return(timeseq)
  
}