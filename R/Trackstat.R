as.Track <- function(x,y,t,covariate){
  stopifnot(length(x)>0 |  length(y)>0 | length(t)>0)
  sp <- cbind(x,y)
  sp <- SpatialPoints(sp)
  td <- as.POSIXct(paste(t))
  if(missing(covariate)) covariate <- data.frame(d=rep(1,length(x)))
  return(Track(STIDF(sp,time = td,data =covariate)))
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
  return(as.Track(newTrack[,1],newTrack[,2],newTrack[,3]))
}

# range.Track returns the timerange of an object of class Track
rngTrack <- function(X) {
  Y <- cbind(as.data.frame(X)[c(coordnames(X), "time")])
  return(range(Y$time)) 
}

# tsqtracks returns a sequance of time based on a list of tracks (or a single object of class Track) and an argument timestamp
tsqTracks <- function(X, timestamp){
  
  timerange = if (is.list(X)) 
    lapply(X, rngTrack)
  else 
  	rngTrack(X)
  
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
  
  avedist <- unlist(avedist)
  class(avedist) <- c("distrack","numeric")
  attr(avedist,"ppp") <- Y
  attr(avedist,"tsq") <- timeseq[-1]
  return(avedist)
}
print.distrack <- function(x, ...){
  print(as.vector(x), ...)
}

plot.distrack <- function(x,...){
  x = unclass(x)
  plot(attr(x,"tsq")[1:length(x)], x, xlab="time",ylab="average distance",...)
}


uniqueTrack <- function(X){
  X <-  cbind(as.data.frame(X)[c(coordnames(X), "time")])
  X <- X[!duplicated(X),]
  return(as.Track(X[,1],X[,2],X[,3])) 
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
  dx <- (max(Z$xcoor)-min(Z$xcoor))/1000
  dy <- (max(Z$ycoor)-min(Z$ycoor))/1000
  w <- owin(c(min(Z$xcoor)-dx,max(Z$xcoor)+dx),c(min(Z$ycoor)-dy,max(Z$ycoor)+dy))
  
  Tppp <- lapply(X=1:length(allZ), function(i){
    p <- as.ppp(allZ[[i]][,-c(3,4)],W=w)
    marks(p) <- allZ[[i]][,4]
    return(p)
  })
  return(Tppp)
}

density.Track <- function(x, ..., timestamp) {
  stopifnot(length(x)>1 & is.list(x))
  
  if (missing(timestamp)) stop("set timestamp") 
  
  p <- as.Track.ppp(x, timestamp)
  p <- p[!sapply(p, is.null)] 
  imlist <- lapply(p, density.ppp, ...)
  out <- Reduce("+", imlist) / length(imlist)
  attr(out, "Tracksim") <- imlist
  attr(out, "ppps") <- p
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
  print(x, ...) 
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

print.arwlen <- function(x, ...){
  print(as.vector(x), ...)
}

plot.arwlen <- function(x,...){
  tsq <- attr(x,"tsq")
  tsq <- tsq[-c(1,length(tsq))]
  plot(tsq,x,xlab="time",ylab="average movement",...)
}

chimaps <- function(X,timestamp,rank,...){
  if (!is.numeric(rank)) stop("rank must be numeric")
  if (rank < 1 | rank >length(X)) stop("rank must be number between one and the number of Tracks")
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  timeseq <- tsqTracks(X,timestamp = timestamp)
  d <- density.Track(X, timestamp = timestamp,...)
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
                sigma=c("default","bw.diggle","bw.ppl"," bw.scott"),...){
  
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  if (missing(q)) q <- 0
  
  cor <- match.arg(correction,correction)
  bw <- match.arg(sigma,sigma)
  if (bw == "default") {
    Y <- as.Track.ppp(X,timestamp = timestamp)
    W <- Y[[1]]$window
    ripley <- min(diff(W$xrange), diff(W$yrange))/4
    rr <- seq(0,ripley,length.out = 513)
    
    K <- lapply(X=1:length(Y), function(i){
      kk <- Kinhom(Y[[i]],correction=cor,r=rr,...)
      return(as.data.frame(kk))
    })
    Kmat <- matrix(nrow = length(K[[1]]$theo),ncol = length(K))
    for (i in 1:length(K)) {
      Kmat[,i] <- K[[i]][,3]
    }
  }
  else{
    bw <- match.fun(bw)
    ZZ <- density.Track(X, timestamp = timestamp, bw)
    
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
print.KTrack <- function(x, ...){
  print("variability area of K-function", ...)
}

plot.KTrack <- function(x,type="l",col= "grey70",...){
  ylim <- c(min(c(x$lowk,x$theo)),max(c(x$upk,x$theo)))
  plot(x$r,x$lowk,ylim=ylim,type=type,xlab="",ylab="",...)
  title(ylab=expression(K[inhom](r)),xlab="r", line=2.2, cex.lab=1.2)
  points(x$r,x$upk,type=type)
  polygon(c(x$r, rev(x$r)), c(x$upk, rev(x$lowk)),
          col = col, border = NA)
  points(x$r,x$theo,type=type,col=2)
  points(x$r,x$avek,type=type)
  legend(0,max(c(x$upk,x$theo)),col = c(2,0,1),legend=c(expression(K[inhom]^{pois}),"",expression(bar(K)[inhom])),lty=c(1,1))
}

pcfinhom.Track <- function(X,timestamp,
                           correction = c("translate", "Ripley"), q,
                           sigma=c("default", "bw.diggle", "bw.ppl", "bw.scott"), ...) {
  
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  if (missing(q)) q <- 0
  
  cor <- match.arg(correction,correction)
  bw <- match.arg(sigma,sigma)
  
  if (bw == "default"){
    Y <- as.Track.ppp(X,timestamp = timestamp)
    
    W <- Y[[1]]$window
    ripley <- min(diff(W$xrange), diff(W$yrange))/4
    rr <- seq(0,ripley,length.out = 513)
    
    g <- lapply(X=1:length(Y), function(i){
      gg <- pcfinhom(Y[[i]],correction=cor,r=rr,...)
      return(as.data.frame(gg))
    })
    gmat <- matrix(nrow = length(g[[1]]$theo),ncol = length(g))
    for (i in 1:length(g)) {
      gmat[,i] <- g[[i]][,3]
    }
  }
  else {
    
    bw <- match.fun(bw)
    ZZ <- density.Track(X, timestamp = timestamp, bw)
    
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
  }
  gmat <- gmat[-1,]
  
  lowg <- numeric()
  upg <- numeric()
  aveg <- numeric()
  for (i in 1:nrow(gmat)) {
    aveg[i] <- mean(gmat[i,])
    lowg[i] <- quantile(gmat[i,],q)
    upg[i] <- quantile(gmat[i,],1-q)
  }
  
  out <- data.frame(lowg=lowg,upg=upg,aveg=aveg,r=g[[1]]$r[-1],theo=g[[1]]$theo[-1])
  class(out) <- c("list","gTrack")
  attr(out,"out") <- out
  return(out)
}


print.gTrack <- function(x, ...){
  print("variability area of pair correlatio function", ...)
}

plot.gTrack <- function(x,type="l",col= "grey70",...){
  ylim <- c(min(x$lowg),max(x$upg))
  plot(x$r,x$lowg,ylim=ylim,xlab="r",ylab=expression(g[inhom](r)),type=type,...)
  points(x$r,x$upg,type=type)
  polygon(c(x$r, rev(x$r)), c(x$upg, rev(x$lowg)),
          col = col, border = NA)
  points(x$r,x$theo,type=type,col=2)
  points(x$r,x$aveg,type=type)
  legend(0.01*max(x$r),max(x$upg),col = c(2,0,1),
         legend=c(expression(g[inhom]^{pois}),"",expression(bar(g)[inhom])),lty=c(1,1))
}

auto.arima.Track <- function(X,...){
  if (! requireNamespace("forecast", quietly = TRUE))
    stop("package forecast required, please install it first")

  stopifnot(class(X)=="Track")
  xseries <- coordinates(X)[,1]
  yseries <- coordinates(X)[,2]
  
  xfit <- forecast::auto.arima(xseries,...)
  yfit <- forecast::auto.arima(yseries,...)
  
  out <- list(xfit,yfit)
  attr(out,"models") <- out
  class(out) <- c("ArimaTrack")
  return(out)
}

print.ArimaTrack <- function(x, ...){
  attributes(x) <- NULL
  cat("Arima model fitted to x-coordinate: ");
  cat(paste0(x[[1]]),"\n")
  cat("Arima model fitted to y-coordinate: ");
  cat(paste0(x[[2]]))
}
