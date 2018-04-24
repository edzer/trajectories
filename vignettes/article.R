### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("MASS")


###################################################
### code chunk number 2: article.Rnw:140-141
###################################################
library(trajectories)


###################################################
### code chunk number 3: article.Rnw:146-157
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
### code chunk number 4: article.Rnw:161-162
###################################################
plot(A1)


###################################################
### code chunk number 5: article.Rnw:172-178
###################################################
x <- runif(10,0,1)
y <- runif(10,0,1)
date <- seq(as.POSIXct("2015-1-1 0:00"), as.POSIXct("2015-1-1 9:00"),
               by = "hour")
records <- as.data.frame(rpois(10,5))
as.Track(x,y,date,covariate = records)


###################################################
### code chunk number 6: article.Rnw:184-194
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
### code chunk number 7: article.Rnw:200-220
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
### code chunk number 8: article.Rnw:268-272
###################################################
dim(A1)
dim(B1)
stbox(A1)
downsample(A1,B1)


###################################################
### code chunk number 9: article.Rnw:276-277
###################################################
stplot(Tr, attr = "co2", arrows = TRUE, lwd = 3, by = "IDs")


###################################################
### code chunk number 10: article.Rnw:287-293
###################################################
set.seed(10)
x <- rTrack();x
y <- rTrack(transform = T);y
m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
w <- rTrack(bbox = m,transform = T);w
z <- rTrack(bbox = m,transform = T,nrandom = T);z


###################################################
### code chunk number 11: article.Rnw:297-300
###################################################
par(mfrow=c(2,2),mar=rep(2.2,4))
plot(x,lwd=2,main="x");plot(y,lwd=2,main="y")
plot(w,lwd=2,main="w");plot(z,lwd=2,main="z")


###################################################
### code chunk number 12: article.Rnw:314-317
###################################################
## EJP:
#data("Beijing")
load("/home/edzer/data/mehdi/taxi/Y.RData")
Beijing = Y

library(forecast)
auto.arima.Track(Beijing[[5]])


###################################################
### code chunk number 13: article.Rnw:350-353
###################################################
 tracks1 <- Tracks(list(Beijing[[1]],Beijing[[2]]))
 tracks2 <- Tracks(list(Beijing[[3]],Beijing[[4]]))
 dists(tracks1,tracks2,mean)


###################################################
### code chunk number 14: article.Rnw:368-375
###################################################
meandist <- avedistTrack(Beijing,timestamp = "20 mins")
plot(meandist,type="l",lwd=2)
distinframe <- data.frame(tsq=meandist$timeseq,dist=meandist$avedist)
dist3rd <- distinframe[substr(distinframe$tsq,start = 1,stop=10)==
                         "2008-02-03",]
plot(dist3rd$tsq,dist3rd$dist,type="l",xlab="time",
      ylab="average distance",lwd=2)


###################################################
### code chunk number 15: article.Rnw:412-414
###################################################
b <- Track.idw(Beijing,timestamp = "20 mins",epsilon=1000)
plot(b,main="",ribwid=0.04,ribsep=0.02)


###################################################
### code chunk number 16: article.Rnw:426-434
###################################################
 q <- avemove(Beijing,timestamp = "20 mins",epsilon=1000)
 par(mfrow=c(1,2))
 plot(q,type="l",lwd=2)
 qdata <- data.frame(q,attr(q,"tsq")[-c(1,length(attr(q,"tsq")))])
 colnames(qdata) <- c("dist","startingtime")
 q3rd <- qdata[substr(qdata$startingtime,start = 1,stop=10)=="2008-02-03",]
 plot(q3rd$startingtime,q3rd$dist,type="l",xlab="time (hour)"
      ,ylab="average movement",lwd=2)


###################################################
### code chunk number 17: article.Rnw:476-492
###################################################
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


###################################################
### code chunk number 18: article.Rnw:517-535
###################################################
ch <- chimaps(Beijing,timestamp = "20 mins",rank = 200)
chall <- attr(ch,"ims")
minmax <- mclapply(X=1:length(chall),function(i){
    return(list(min(chall[[i]]$v),max(chall[[i]]$v)))
  })
minmax <- do.call("rbind",minmax)
col5 <- colorRampPalette(c('blue','white','red'))
color_levels=200 
par(mar=c(0,0,1,1))
plot(chall[[51]],zlim=c(-max(abs(unlist(minmax))),max(abs(unlist(minmax))))
          ,main=attr(ch,"timevec")[51],ribwid=0.04,ribsep=0.02,
          col=col5(n=color_levels))
plot(chall[[75]],zlim=c(-max(abs(unlist(minmax))),max(abs(unlist(minmax))))
          ,main=attr(ch,"timevec")[75],ribwid=0.04,ribsep=0.02,
          col=col5(n=color_levels))
plot(chall[[104]],zlim=c(-max(abs(unlist(minmax))),max(abs(unlist(minmax))))
          ,main=attr(ch,"timevec")[104],ribwid=0.04,ribsep=0.02,
          col=col5(n=color_levels))


###################################################
### code chunk number 19: article.Rnw:587-591
###################################################
 K <- Kinhom.Track(Beijing,timestamp = "20 mins",q=0)
 plot(K)
 g <- pcfinhom.Track(Beijing,timestamp = "20 mins",q=0)
 plot(g)


