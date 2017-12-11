# function as.Track  accepts a dataframe (let say X) with three columns "xcoor", "ycoor" and time (or any class that can be converted to a data.frame)
# and converts it to an object of class Track. It can also accepts covariates for the corresponding 
# locations, covariates must be a dataframe with some columns and length of each column is equal
# to length of each column in X.
as.Track <- function(X,covariate){
  stopifnot(nrow(X)>0)
  # colnames(X) <- c("xcoor","ycoor","time")
  if(!is.data.frame(X)) X <- as.data.frame(X)
  sp <- cbind(x=X$xcoor,y=X$ycoor)
  sp <- SpatialPoints(sp)
  t <- as.POSIXct(paste(X$date,X$time))
  if(missing(covariate)) covariate <- data.frame(d=rep(1,length(X$xcoor)))
  
  return(Track(STIDF(sp,time = t,data =covariate)))
}


# function reTrack accepts X as an object of class Track. Output is a reconstructed Track (an object of class Track), based on "timestamp".
# It only returns the interpolated points.
reTrack <- function(X,at=c("track","dfrm"),timestamp=timestamp,tsq=NULL){
  
  if (missing(tsq)) tsq <- tsqTracks(X,timestamp = timestamp)
  if(missing(at)) at <- "track"
  Xrange <- rngTrack(X)
  X <-  cbind(as.data.frame(X)[c(coordnames(X), "time")])
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
  if (at=="dfrm") {attr(newTrack,"tsq") <-tsq;return(newTrack) }
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



# function avedistTrack accepts X as a list of tracks and reports the average distance between
# tracks over time, output is an object of class "distrack"
avedistTrack <- function(X,timestamp){
  
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp")  
  # calculate a sequance of time to interpolate tracks within this sequance
  timeseq <- tsqTracks(X,timestamp = timestamp)
  
  Y <- as.Track.ppp(X,timestamp)
  
  avedist <- lapply(X=1:length(Y), function(i){
    pd <-  pairdist(Y[[i]])
    mean(pd[pd>0])
  })
  
  avedist <- data.frame(timeseq[-1],unlist(avedist))
  colnames(avedist) <- c("timeseq","avedist")
  class(avedist) <- c("distrack")
  attr(avedist,"ppp") <- Y
  return(avedist)
}
print.distrack <- function(x){
  print(as.vector(x$avedist))
}

plot.distrack <- function(x,...){
  plot(x$timeseq,x$avedist,,xlab="time",ylab="average distance",...)
}


rmdupTrack <- function(X){
  X <-  cbind(as.data.frame(X)[c(coordnames(X), "time")])
  X <- X[!duplicated(X),]
  return(as.Track(X)) 
}


as.Track.ppp <- function(X,timestamp){
  
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  # calculate a sequance of time to interpolate tracks within this sequance
  timeseq <- tsqTracks(X,timestamp = timestamp)
  
  # reconstruct tracks in sequance timeseq
  Z <- lapply(X,reTrack,tsq = timeseq,at="dfrm")
  id <- rep(1:length(Z),sapply(Z, nrow))
  Z <- do.call("rbind",Z)
  Z <- cbind(Z,id)
  allZ <- split(Z,Z[,3])
  w <- owin(c(min(Z$xcoor)-0.001,max(Z$xcoor)+0.001),c(min(Z$ycoor)-0.001,max(Z$ycoor)+0.001))
  
  Tppp <- lapply(X=1:length(allZ), function(i){
    p <- as.ppp(allZ[[i]][,-c(3,4)],W=w)
    marks(p) <- allZ[[i]][,4]
    return(p)
  })
  return(Tppp)
}

density.Track <- function(X,timestamp,...){
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  
  p <- as.Track.ppp(X,timestamp)
  p <- p[!sapply(p, is.null)] 
  imlist <- lapply(p, density.ppp,...)
  out <- Reduce("+",imlist)/length(imlist)
  attr(out,"Tracksim") <- imlist
  attr(out,"ppps") <- p
  return(out)
}

as.Track.arrow <- function(X,timestamp,epsilon=epsilon){
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  
  Z <- as.Track.ppp(X,timestamp)
  Z <- Z[!sapply(Z, is.null)]
  wind <- Z[[1]]$window
  arrows <- list()
  Y <- list()
  for (i in 1:length(Z)) {
    if(i==length(Z)) break()
    j <- i+1
    m1 <- match(marks(Z[[i]]),marks(Z[[j]]))
    m2 <- match(marks(Z[[j]]),marks(Z[[i]]))
    m1 <- m1[!is.na(m1)]
    m2 <- m2[!is.na(m2)]
    x <- Z[[j]][m1]
    y <- Z[[i]][m2]
    l <- psp(y$x,y$y,x$x,x$y,window = wind)
    arrows[[i]] <- l
    center <- midpoints.psp(l)
    mark <- lengths.psp(l)
    marks(center) <- mark
    if (missing(epsilon)) epsilon <- 0
    Y[[i]] <- center[mark>epsilon]
  }
  class(Y) <- c("list","Trrow")
  attr(Y,"psp") <- arrows
  return(Y)  
}

print.Trrow <- function(x, ...) { 
  attributes(x) <- NULL 
  print(x) 
} 

Track.idw <- function(X,timestamp,epsilon=epsilon,...){
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp")
  
  Y <- as.Track.arrow(X,timestamp,epsilon=epsilon)
  Z <- lapply(Y, idw,...)
  meanIDW <- Reduce("+",Z)/length(Z)
  return(meanIDW)
}

avemove <- function(X,timestamp,epsilon=epsilon){
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  timeseq <- tsqTracks(X,timestamp = timestamp)
  Y <- as.Track.arrow(X,timestamp,epsilon=epsilon)
  Z <- attr(Y,"psp")
  preout <- lapply(X=1:length(Z), function(i){
    mean(lengths.psp(Z[[i]]))
  })
  out <- unlist(preout)
  class(out) <- c("numeric", "arwlen")
  attr(out,"tsq") <- timeseq
  return(out)
}

print.arwlen <- function(x){
  print(as.vector(x))
}

plot.arwlen <- function(x,...){
  tsq <- attr(x,"tsq")
  tsq <- tsq[c(-1,-length(tsq))]
  plot(tsq,x,xlab="time",ylab="average movement",...)
}

chimaps <- function(X,timestamp,rank,...){
  if (!is.numeric(rank)) stop("rank must be numeric")
  if (rank < 1 | rank >length(X)) stop("rank must be number between one and the number of Tracks")
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  timeseq <- tsqTracks(X,timestamp = timestamp)
  d <- density.Track(X,timestamp,...)
  imlist <- attr(d,"Tracksim")
  sumim <- Reduce("+",imlist)
  chi <- lapply(X=1:length(imlist),FUN = function(i){
    E1 <- sumim*sum(imlist[[i]]$v)/(sum(sumim$v))
    return((imlist[[i]]-E1)/sqrt(E1))
  })
  out <- chi[[rank]]
  attr(out,"ims") <- chi
  attr(out,"time") <- timeseq[rank]
  attr(out,"timevec") <- timeseq
  return(out)
}

Kinhom.Track <- function(X,timestamp,
                correction=c("border", "bord.modif", "isotropic", "translate"),q,
                sigma=c("bw.diggle","bw.ppl"," bw.scott"),...){
  
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  
  cor <- match.arg(correction,correction)
  bw <- match.arg(sigma,sigma)
  bw <- match.fun(bw)
  ZZ <- density.Track(X,timestamp,bw)
  
  Z <- attr(ZZ,"Tracksim")
  Y <- attr(ZZ,"ppps")
  W <- Y[[1]]$window
  ripley <- min(diff(W$xrange), diff(W$yrange))/4
  rr <- seq(0,ripley,length.out = 513)
  
  K <- lapply(X=1:length(Y), function(i){
      kk <- Kinhom(Y[[i]],lambda = Z[[i]],correction=cor,r=rr,...)
      return(as.data.frame(kk))
  })
  Kmat <- matrix(nrow = length(K[[1]]$theo),ncol = length(K))
  for (i in 1:length(K)) {
    Kmat[,i] <- K[[i]][,3]
  }
  # Kmat <- as.data.frame(K)
  lowk <- numeric()
  upk <- numeric()
  avek <- numeric()
  for (i in 1:nrow(Kmat)) {
    avek[i] <- mean(Kmat[i,])
    lowk[i] <- quantile(Kmat[i,],q)
    upk[i] <- quantile(Kmat[i,],1-q)
  }
  
  out <- data.frame(lowk=lowk,upk=upk,avek=avek,r=K[[1]]$r,theo=K[[1]]$theo)
  class(out) <- c("list","KTrack")
  attr(out,"out") <- out
  return(out)
}
print.KTrack <- function(x){
  print("variability area of K-function")
}

plot.KTrack <- function(x,type="l",col= "grey70",...){
  ylim <- c(min(x$lowk),max(x$upk))
  plot(x$r,x$lowk,ylim=ylim,xlab="r",ylab="K(r)",type=type,...)
  points(x$r,x$upk,type=type)
  polygon(c(x$r, rev(x$r)), c(x$upk, rev(x$lowk)),
          col = col, border = NA)
  points(x$r,x$theo,type=type,col=2)
  points(x$r,x$avek,type=type)
  legend(0,max(x$upk),col = c(2,1),legend=c("poisson","average"),lty=c(1,1))
}

pcfinhom.Track <- function(X,timestamp,
                           correction = c("translate", "Ripley"),q,...){
  
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  
  cor <- match.arg(correction,correction)
  
  ZZ <- density.Track(X,timestamp)
  
  Z <- attr(ZZ,"Tracksim")
  Y <- attr(ZZ,"ppps")
  
  g <- lapply(X=1:length(Y), function(i){
    gg <- pcfinhom(Y[[i]],lambda = Z[[i]],correction=cor,...)
    return(as.data.frame(gg))
  })
  gmat <- matrix(nrow = length(g[[1]]$theo),ncol = length(g))
  for (i in 1:length(g)) {
    gmat[,i] <- g[[i]][,3]
  }
  # Kmat <- as.data.frame(K)
  lowg <- numeric()
  upg <- numeric()
  aveg <- numeric()
  for (i in 1:nrow(gmat)) {
    aveg[i] <- mean(gmat[i,])
    lowg[i] <- quantile(gmat[i,],q)
    upg[i] <- quantile(gmat[i,],1-q)
  }
  
  out <- data.frame(lowg=lowg,upg=upg,aveg=aveg,r=g[[1]]$r,theo=g[[1]]$theo)
  class(out) <- c("list","gTrack")
  attr(out,"out") <- out
  return(out)
}


print.gTrack <- function(x){
  print("variability area of pair correlatio function")
}

plot.gTrack <- function(x,type="l",col= "grey70",...){
  ylim <- c(min(x$lowg),max(x$upg))
  plot(x$r,x$lowg,ylim=ylim,xlab="r",ylab="g",type=type,...)
  points(x$r,x$upg,type=type)
  polygon(c(x$r, rev(x$r)), c(x$upg, rev(x$lowg)),
          col = col, border = NA)
  points(x$r,x$theo,type=type,col=2)
  points(x$r,x$aveg,type=type)
}


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
    xy = cbind(cumsum(rnorm(n, sd = sd0)) + origin[1], cumsum(rnorm(n, 
                                                                    sd = sd0)) + origin[2])
  else {xy = cbind(origin[1] + cumsum(as.vector(arima.sim(list(ar = ar), 
                                                          n, sd = sd0, ...))), 
                   origin[2] + cumsum(as.vector(arima.sim(list(ar = ar), 
                                                          
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


