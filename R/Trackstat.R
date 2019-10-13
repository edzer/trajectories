as.list.Tracks <- function(x,...){
  stopifnot(class(x)=="Tracks")
  return(as.list(x@tracks,...))
}

as.list.TracksCollection <- function(x,...){
  stopifnot(class(x)=="TracksCollection")
  out <-  lapply(X=1:length(x@tracksCollection), function(i){
    as.list.Tracks(x@tracksCollection[[i]],...)
  })
  return(unlist(out, recursive=FALSE))
}


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
  Xrange <- range.Track(X)
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
range.Track <- function(X,...) {
  out <- as.data.frame(X)
  return(range(out$time,...)) 
}

range.Tracks <- function(X,...) {
  out <- lapply(X@tracks,as.data.frame)
  out <- do.call(rbind,out)
  return(range(out$time,...))
}

range.TracksCollection <- function(X,...) {
  outf <- list()
  for (i in 1:length(X@tracksCollection)) {
    out <- lapply(X@tracksCollection[[i]]@tracks,as.data.frame)
    outf[[i]] <- do.call(rbind,out)
  }
  outf <- do.call(rbind,outf)
  
  return(range(outf$time,...))
}
# tsqtracks returns a sequance of time based on a list of tracks (or a single object of class Track) and an argument timestamp
tsqTracks <- function(X, timestamp,from=NULL,to=NULL){
  
  if(!(is.null(from) & is.null(to))) {
    return(seq(from=as.POSIXct(strftime(from)),to=as.POSIXct(strftime(to)),by = timestamp))
  }
  if(class(X)=="Track" | class(X)=="Tracks" | class(X)=="TracksCollection") timerange <- range(X)
  
  if(class(X)=="list") timerange <- lapply(X, range) ; timerange <- range(timerange)
  
  class(timerange) <- c('POSIXct','POSIXt')
  # a seq from the range has been created every timestamp
  if (missing(timestamp)) stop("set timestamp")
  timeseq <- seq(from=as.POSIXct(strftime(timerange[1])),to=as.POSIXct(strftime(timerange[2])),by = timestamp)
  
  return(timeseq)
  
}



# function avedistTrack accepts X as a list of tracks and reports the average distance between
# tracks over time, output is an object of class "distrack"
avedistTrack <- function(X,timestamp){
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  
  stopifnot(length(X)>1 & is.list(X))
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  
  if (missing(timestamp)) stop("set timestamp")  
  # calculate a sequance of time to interpolate tracks within this sequance
  timeseq <- tsqTracks(X,timestamp = timestamp)
  
  Y <- as.Track.ppp(X,timestamp)
  
  avedist <- lapply(X=1:length(Y), function(i){
    pd <-  spatstat::pairdist(Y[[i]])
    mean(pd[pd>0])
  })
  
  avedist <- unlist(avedist)
  class(avedist) <- c("distrack","numeric")
  attr(avedist,"ppp") <- Y
  attr(avedist,"tsq") <- attr(Y,"tsq")
  return(avedist)
}
print.distrack <- function(x){
  print(as.vector(x))
}

plot.distrack <- function(x,...){
  x = unclass(x)
  plot(attr(x,"tsq"), x, xlab="time",ylab=expression(italic(bar(D))),...)
}


unique.Track <- function(x,...){
  x <-  cbind(as.data.frame(x)[c(coordnames(x), "time")])
  x <- unique(x,...)
  return(as.Track(x[,1],x[,2],x[,3])) 
}


as.Track.ppp <- function(X,timestamp){
  
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  stopifnot(length(X)>1 & is.list(X))
  
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  
  
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
  w <- spatstat::owin(c(min(Z$xcoor)-dx,max(Z$xcoor)+dx),c(min(Z$ycoor)-dy,max(Z$ycoor)+dy))
  
  Tppp <- lapply(X=1:length(allZ), function(i){
    p <- spatstat::as.ppp(allZ[[i]][,-c(3,4)],W=w)
    p <- spatstat::`marks<-`(p, value = allZ[[i]][,4])
    return(p)
  })
  class(Tppp) <- c("list","ppplist")
  attr(Tppp,"tsq") <- as.POSIXlt.character(attributes(allZ)$names)
  return(Tppp)
}

print.ppplist <- function(x){
  attributes(x) <- NULL 
  print(x) 
}

density.list <- function(x, timestamp, method=c("kernel","Voronoi"), ...) {
  stopifnot(class(x)=="list" | class(x)=="Tracks" | class(x)=="TracksCollection")
  
  if(class(x)=="Tracks") x <- as.list.Tracks(x)
  if (class(x)=="TracksCollection") x <- as.list.TracksCollection(x)
  
  stopifnot(length(x)>1 & is.list(x))
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  
  if (missing(timestamp)) stop("set timestamp")
  if(missing(method)) method <- "kernel"
  
  p <- as.Track.ppp(x, timestamp)
  p <- p[!sapply(p, is.null)] 
  if(any(method == "kernel")){
    imlist <- lapply(p, spatstat::density.ppp, ...)  
  }
  else{
    imlist <- lapply(p, spatstat::densityVoronoi, ...)  
  }
  
  out <- Reduce("+", imlist) / length(imlist)
  attr(out, "Tracksim") <- imlist
  attr(out, "ppps") <- p
  return(out)
}

as.Track.arrow <- function(X,timestamp,epsilon=0){
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  stopifnot(length(X)>1 & is.list(X))
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  
  if (missing(timestamp)) stop("set timestamp") 
  if(missing(epsilon))  epsilon <- 0
  
  Z <- as.Track.ppp(X,timestamp)
  tsq <- attr(Z,"tsq")
  Z <- Z[!sapply(Z, is.null)]
  wind <- Z[[1]]$window
  arrows <- list()
  Y <- list()
  for (i in 1:length(Z)) {
    if(i==length(Z)) break()
    j <- i+1
    m1 <- match(spatstat::marks(Z[[i]]),spatstat::marks(Z[[j]]))
    m2 <- match(spatstat::marks(Z[[j]]),spatstat::marks(Z[[i]]))
    m1 <- m1[!is.na(m1)]
    m2 <- m2[!is.na(m2)]
    x <- Z[[j]][m1]
    y <- Z[[i]][m2]
    l <- spatstat::psp(y$x,y$y,x$x,x$y,window = wind)
    arrows[[i]] <- l
    center <- spatstat::midpoints.psp(l)
    mark <- spatstat::lengths.psp(l)
    center <- spatstat::`marks<-`(center, value = mark)
    if (missing(epsilon)) epsilon <- 0
    Y[[i]] <- center[mark>epsilon]
  }
  class(Y) <- c("list","Trrow")
  attr(Y, "psp") <- arrows
  attr(Y,"time") <- tsq[-length(tsq)]
  return(Y)  
}

print.Trrow <- function(x) { 
  attributes(x) <- NULL 
  print(x) 
} 

idw.Track <- function(X,timestamp,epsilon=0,...){
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp")
  if(missing(epsilon))  epsilon <- 0
  
  Y <- as.Track.arrow(X,timestamp,epsilon=epsilon)
  Z <- lapply(Y, spatstat::idw, ...)
  meanIDW <- Reduce("+",Z)/length(Z)
  return(meanIDW)
}

avemove <- function(X,timestamp,epsilon=0){
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  stopifnot(length(X)>1 & is.list(X))
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  
  if (missing(timestamp)) stop("set timestamp") 
  timeseq <- tsqTracks(X,timestamp = timestamp)
  if (missing(epsilon)) epsilon <- 0
  Y <- as.Track.arrow(X,timestamp,epsilon=epsilon)
  Z <- attr(Y,"psp")
  preout <- lapply(X=1:length(Z), function(i){
    mean(spatstat::lengths.psp(Z[[i]]))
  })
  out <- unlist(preout)
  class(out) <- c("numeric", "arwlen")
  attr(out,"time") <- attr(Y,"time")
  return(out)
}

print.arwlen <- function(x){
  print(as.vector(x))
}

plot.arwlen <- function(x,...){
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  x = unclass(x)
  tsq <- attr(x,"time")
  plot(tsq,x,xlab="time",ylab="average movement",...)
}

chimaps <- function(X,timestamp,rank,...){
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  stopifnot(length(X)>1 & is.list(X))
  
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  
  if(missing(rank)) rank <- 1
  if (!is.numeric(rank)) stop("rank must be numeric")
  if (rank < 1 | rank >length(X)) stop("rank must be number between one and the number of Tracks")
  stopifnot(length(X)>1 & is.list(X))
  
  if (missing(timestamp)) stop("set timestamp") 
  timeseq <- tsqTracks(X,timestamp = timestamp)
  d <- density.list(X, timestamp = timestamp,...)
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
                         correction=c("border", "bord.modif", "isotropic", "translate"),q=0,
                         sigma=c("default","bw.diggle","bw.ppl"," bw.scott"),...){
  
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
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
      kk <- spatstat::Kinhom(Y[[i]],correction=cor,r=rr,...)
      return(as.data.frame(kk))
    })
    Kmat <- matrix(nrow = length(K[[1]]$theo),ncol = length(K))
    for (i in 1:length(K)) {
      Kmat[,i] <- K[[i]][,3]
    }
  }
  else{
    bw <- match.fun(bw)
    ZZ <- density.list(X, timestamp = timestamp, bw)
    
    Z <- attr(ZZ,"Tracksim")
    Y <- attr(ZZ,"ppps")
    W <- Y[[1]]$window
    ripley <- min(diff(W$xrange), diff(W$yrange))/4
    rr <- seq(0,ripley,length.out = 513)
    
    K <- lapply(X=1:length(Y), function(i){
      kk <- spatstat::Kinhom(Y[[i]],lambda = Z[[i]],correction=cor,r=rr,...)
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
    avek[i] <- mean(Kmat[i,],na.rm = TRUE)
    lowk[i] <- quantile(Kmat[i,],q,na.rm = TRUE)
    upk[i] <- quantile(Kmat[i,],1-q,na.rm = TRUE)
  }
  
  out <- data.frame(lowk=lowk,upk=upk,avek=avek,r=K[[1]]$r,theo=K[[1]]$theo)
  class(out) <- c("list","KTrack")
  attr(out,"out") <- out
  return(out)
}
print.KTrack <- function(x){
  print("variability area of K-function")
}

plot.KTrack <- function(x,type="l",col= "grey70",cex=1,line=2.2,...){
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  ylim <- c(min(c(x$lowk,x$theo)),max(c(x$upk,x$theo)))
  plot(x$r,x$lowk,ylim=ylim,type=type,ylab="",xlab="r",...)
  title(ylab=expression(italic(K[inhom](r))),line = line,...)
  points(x$r,x$upk,type=type)
  polygon(c(x$r, rev(x$r)), c(x$upk, rev(x$lowk)),
          col = col, border = NA)
  points(x$r,x$theo,type=type,col=2)
  points(x$r,x$avek,type=type)
  legend(0,max(c(x$upk,x$theo)),col = c(2,1,"grey70","grey70"),
         legend=c(expression(italic(K[inhom]^{pois})),
                  expression(italic(bar(K)[inhom])),
                  expression(italic({hat(K)[inhom]^{high}}(r))),
                  expression(italic({hat(K)[inhom]^{low}}(r)))
         ),
         lty=c(1,1,1,1),cex = cex)
}

pcfinhom.Track <- function(X,timestamp,
                           correction = c("translate", "Ripley"), q,
                           sigma=c("default", "bw.diggle", "bw.ppl", "bw.scott"), ...) {
  
  stopifnot(class(X)=="list" | class(X)=="Tracks" | class(X)=="TracksCollection")
  
  if(class(X)=="Tracks") X <- as.list.Tracks(X)
  if (class(X)=="TracksCollection") X <- as.list.TracksCollection(X)
  
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
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
      gg <- spatstat::pcfinhom(Y[[i]],correction=cor,r=rr,...)
      return(as.data.frame(gg))
    })
    gmat <- matrix(nrow = length(g[[1]]$theo),ncol = length(g))
    for (i in 1:length(g)) {
      gmat[,i] <- g[[i]][,3]
    }
  }
  else {
    
    bw <- match.fun(bw)
    ZZ <- density.list(X, timestamp = timestamp, bw)
    
    Z <- attr(ZZ,"Tracksim")
    Y <- attr(ZZ,"ppps")
    
    g <- lapply(X=1:length(Y), function(i){
      gg <- spatstat::pcfinhom(Y[[i]],lambda = Z[[i]],correction=cor,...)
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
    aveg[i] <- mean(gmat[i,],na.rm = TRUE)
    lowg[i] <- quantile(gmat[i,],q,na.rm = TRUE)
    upg[i] <- quantile(gmat[i,],1-q,na.rm = TRUE)
  }
  
  out <- data.frame(lowg=lowg,upg=upg,aveg=aveg,r=g[[1]]$r[-1],theo=g[[1]]$theo[-1])
  class(out) <- c("list","gTrack")
  attr(out,"out") <- out
  return(out)
}


print.gTrack <- function(x){
  print("variability area of pair correlatio function")
}

plot.gTrack <- function(x,type="l",col= "grey70",cex=1,line=2.2,...){
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("spatstat required: install first?")
  ylim <- c(min(x$lowg),max(x$upg))
  plot(x$r,x$lowg,ylim=ylim,xlab="r",ylab="",type=type,...)
  title(ylab=expression(italic(g[inhom](r))),line = line,...)
  points(x$r,x$upg,type=type)
  polygon(c(x$r, rev(x$r)), c(x$upg, rev(x$lowg)),
          col = col, border = NA)
  points(x$r,x$theo,type=type,col=2)
  points(x$r,x$aveg,type=type)
  legend(0.01*max(x$r),max(x$upg),col = c(2,1,"grey70","grey70"),
         legend=c(expression(italic(g[inhom]^{pois})),
                  expression(italic(bar(g)[inhom])),
                  expression(italic({hat(g)[inhom]^{high}}(r))),
                  expression(italic({hat(g)[inhom]^{low}}(r)))
         ),
         lty=c(1,1,1,1),cex=cex)
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

print.ArimaTrack <- function(x){
  attributes(x) <- NULL
  cat("Arima model fitted to x-coordinate: ");
  cat(paste0(x[[1]]),"\n")
  cat("Arima model fitted to y-coordinate: ");
  cat(paste0(x[[2]]))
}
