### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("MASS")


###################################################
### code chunk number 2: article.Rnw:109-114
###################################################
do_all <- FALSE
if(do_all){
 install.packages("taxidata", 
 repos = "http://pebesma.staff.ifgi.de",type = "source") 
}


###################################################
### code chunk number 3: article.Rnw:119-133
###################################################
library(trajectories)
if(do_all){
  library(taxidata)
Beijing <- taxidata
Z <- lapply(X=1:length(Beijing), function(i){
  q <-  cut(Beijing[[i]], "day", touch = F)
  return(q@tracks[[3]])
})
plot(Z[[21]],xlim=c(420000,470000),ylim=c(4390000,4455000),lwd=2)
plot(Z[[26]],add=T,col="orange",lwd=2)
plot(Z[[20]],add=T,col=2,lwd=2)
plot(Z[[12]],add=T,col=3,lwd=2)
plot(Z[[15]],add=T,col=4,lwd=2)
}


###################################################
### code chunk number 4: article.Rnw:166-167
###################################################
library(trajectories)


###################################################
### code chunk number 5: article.Rnw:172-183
###################################################
set.seed(10)
t0 = as.POSIXct(as.Date("2013-09-30",tz="CET"))
x = c(7,6,5,5,4,3,3)
y = c(7,7,6,5,5,6,7)
n = length(x)
t = t0 + cumsum(runif(n) * 60)
crs = CRS("+proj=longlat +ellps=WGS84") # longlat
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, 
                data.frame(co2 = rnorm(n,mean = 10)))
A1 = Track(stidf)
A1


###################################################
### code chunk number 6: article.Rnw:187-188
###################################################
plot(A1)


###################################################
### code chunk number 7: article.Rnw:198-204
###################################################
x <- runif(10,0,1)
y <- runif(10,0,1)
date <- seq(as.POSIXct("2015-1-1 0:00"), as.POSIXct("2015-1-1 9:00"),
               by = "hour")
records <- as.data.frame(rpois(10,5))
as.Track(x,y,date,covariate = records)


###################################################
### code chunk number 8: article.Rnw:210-220
###################################################
x = c(7,6,6,7,7)
y = c(6,5,4,4,3)
n = length(x)
t = max(t) + cumsum(runif(n) * 60)
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, 
                data.frame(co2 = rnorm(n,mean = 10)))
A2 = Track(stidf)
# Tracks for person A:
A = Tracks(list(A1=A1,A2=A2))
A


###################################################
### code chunk number 9: article.Rnw:226-246
###################################################
# person B, track 1:
x = c(2,2,1,1,2,3)
y = c(5,4,3,2,2,3)
n = length(x)
t = max(t) + cumsum(runif(n) * 60)
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, 
                data.frame(co2 = rnorm(n,mean = 10)))
B1 = Track(stidf)
# person B, track 2:
x = c(3,3,4,3,3,4)
y = c(5,4,3,2,1,1)
n = length(x)
t = max(t) + cumsum(runif(n) * 60)
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, 
                data.frame(co2 = rnorm(n,mean = 10)))
B2 = Track(stidf)
# Tracks for person B:
B = Tracks(list(B1=B1,B2=B2))
Tr = TracksCollection(list(A=A,B=B))
Tr


###################################################
### code chunk number 10: article.Rnw:295-299
###################################################
dim(A1)
dim(B1)
stbox(A1)
downsample(A1,B1)


###################################################
### code chunk number 11: article.Rnw:303-304
###################################################
stplot(Tr, attr = "co2", arrows = TRUE, lwd = 3, by = "IDs")


###################################################
### code chunk number 12: article.Rnw:314-320
###################################################
set.seed(10)
x <- rTrack();x
y <- rTrack(transform = T);y
m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
w <- rTrack(bbox = m,transform = T);w
z <- rTrack(bbox = m,transform = T,nrandom = T);z


###################################################
### code chunk number 13: article.Rnw:324-327
###################################################
par(mfrow=c(2,2),mar=rep(2.2,4))
plot(x,lwd=2,main="x");plot(y,lwd=2,main="y")
plot(w,lwd=2,main="w");plot(z,lwd=2,main="z")


###################################################
### code chunk number 14: article.Rnw:341-343
###################################################
library(forecast)
auto.arima.Track(A1)


###################################################
### code chunk number 15: article.Rnw:376-389
###################################################
library(xts)
data(A3)
track2 <- A3
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

## create Tracks objects
tracks1 <- Tracks(list(A3, track2))
tracks2 <- Tracks(list(track2, A3))

## calculate distances
## Not run: 
dists(tracks1, tracks2,mean)


###################################################
### code chunk number 16: article.Rnw:404-413
###################################################
 if (do_all){
 meandist <- avedistTrack(Beijing,timestamp = "20 mins")
 plot(meandist,type="l",lwd=2)
 distinframe <- data.frame(tsq=attr(meandist,"tsq"),dist=meandist)
 dist3rd <- distinframe[substr(distinframe$tsq,start = 1,stop=10)==
                         "2008-02-03",]
 plot(dist3rd$tsq,dist3rd$dist,type="l",xlab="time",
       ylab="average distance",lwd=2)
            }


###################################################
### code chunk number 17: article.Rnw:449-453
###################################################
if(do_all){
b <- Track.idw(Beijing,timestamp = "20 mins",epsilon=1000)
plot(b,main="",ribwid=0.04,ribsep=0.02)
}


###################################################
### code chunk number 18: article.Rnw:465-475
###################################################
 if(do_all){
 q <- avemove(Beijing,timestamp = "20 mins",epsilon=1000)
  par(mfrow=c(1,2))
  plot(q,type="l",lwd=2)
  qdata <- data.frame(q,attr(q,"tsq")[-c(1,length(attr(q,"tsq")))])
  colnames(qdata) <- c("dist","startingtime")
  q3rd <- qdata[substr(qdata$startingtime,start = 1,stop=10)=="2008-02-03",]
  plot(q3rd$startingtime,q3rd$dist,type="l",xlab="time (hour)"
       ,ylab="average movement",lwd=2)  
 }


###################################################
### code chunk number 19: article.Rnw:517-535
###################################################
 if(do_all){
 d <- density.Track(Beijing,timestamp = "20 mins",bw.ppl)
 par(mfrow=c(1,2))
 plot(d,main="",ribwid=0.04,ribsep=0.02)
 #focus on the center
 w <- owin(c(440000,455000),c(4410000,4430000))
 pps <- attr(d,"ppps")
 npps <- lapply(X=1:length(pps),FUN = function(i){
   pps[[i]][w]
 })
 
 centerimg <- lapply(X=1:length(npps),FUN = function(i){
   density(npps[[i]],bw.ppl(npps[[i]]))
 })
 fcenterimg <- Reduce("+",centerimg)/length(centerimg)
 
 plot(fcenterimg,main="",ribwid=0.04,ribsep=0.02)
 }


###################################################
### code chunk number 20: article.Rnw:560-581
###################################################
 if(do_all){
ch <- chimaps(Beijing,timestamp = "20 mins",rank = 200)
 chall <- attr(ch,"ims")
 minmax <- lapply(X=1:length(chall),function(i){
     return(list(min(chall[[i]]$v),max(chall[[i]]$v)))
   })
 minmax <- do.call("rbind",minmax)
 col5 <- colorRampPalette(c('blue','white','red'))
 color_levels=200 
 par(mar=c(0,0,1,1))
 par(mfrow=c(1,3))
 plot(chall[[51]],zlim=c(-max(abs(unlist(minmax))),max(abs(unlist(minmax))))
           ,main=attr(ch,"timevec")[51],ribwid=0.04,ribsep=0.02,
           col=col5(n=color_levels))
 plot(chall[[75]],zlim=c(-max(abs(unlist(minmax))),max(abs(unlist(minmax))))
           ,main=attr(ch,"timevec")[75],ribwid=0.04,ribsep=0.02,
           col=col5(n=color_levels))
 plot(chall[[104]],zlim=c(-max(abs(unlist(minmax))),max(abs(unlist(minmax))))
           ,main=attr(ch,"timevec")[104],ribwid=0.04,ribsep=0.02,
           col=col5(n=color_levels))
 }


###################################################
### code chunk number 21: article.Rnw:633-641
###################################################
 if(do_all){
 K <- Kinhom.Track(Beijing,correction = "translate",
                       timestamp = "20 mins",q=0)
  par(mfrow=c(1,2))
  plot(K)
  g <- pcfinhom.Track(Beijing,timestamp = "20 mins",q=0)
  plot(g)
 }


