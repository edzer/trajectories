\name{approxTrack}
\alias{approxTrack}
\alias{approxTracks}
\alias{approxTracksCollection}

\title{
Interpolation for trajectories
}
\description{
Interpolating objects of classes \code{Track}, \code{Tracks} and \code{TracksCollection}.
}
\usage{
approxTrack(track, when, ..., n = 50, by, FUN = stats::approx, warn.if.outside = TRUE)
approxTracks(tr,...)
approxTracksCollection(tc,...)

}
\arguments{
\item{track,tr,tc}{
an object of class
\link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.
}
\item{when}{
a time sequance specifying where interpolation is to take place. 
}
\item{n}{
if when is not specified, interpolation takes place at n equally spaced times spanning the interval [min(index(track)), max(index(track))].
}
\item{by}{
if when is not specified, by can be used to create a sequance of equally spaced times spanning the interval [min(index(track)), max(index(track))].
}
\item{FUN}{
a function to be used for interpolation. 
}
\item{\dots}{
passed to FUN.
}
\item{warn.if.outside}{
warning if when has times outside the range of data.
}
}
\details{
Interpolating Track, Tracks, and TracksCollection based on a desirable function. Default is to use \link[stats]{approx} for a linear interpolation. See exapmle below for the use of other functions.
}
\value{
an object of class Track, Tracks or TracksCollection.
}
\author{
Edzer Pebesma, \email{edzer.pebesma@uni-muenster.de} 
}
\seealso{
\link[stats]{approx}, \link[stats]{spline}
}
\examples{
library(sp)
library(spacetime)

t0 = as.POSIXct(as.Date("2013-09-30",tz="CET"))
# person A, track 1:
x = c(7,6,5,5,4,3,3)
y = c(7,7,6,5,5,6,7)
n = length(x)
set.seed(131)
t = t0 + cumsum(runif(n) * 60)
require(rgdal)
crs = CRS("+proj=longlat +ellps=WGS84") # longlat
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, data.frame(co2 = rnorm(n)))
A1 = Track(stidf)

# person A, track 2:
x = c(7,6,6,7,7)
y = c(6,5,4,4,3)
n = length(x)
t = max(t) + cumsum(runif(n) * 60)
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, data.frame(co2 = rnorm(n)))
A2 = Track(stidf)

# Tracks for person A:
A = Tracks(list(A1=A1,A2=A2))

# person B, track 1:
x = c(2,2,1,1,2,3)
y = c(5,4,3,2,2,3)
n = length(x)
t = max(t) + cumsum(runif(n) * 60)
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, data.frame(co2 = rnorm(n)))
B1 = Track(stidf)

# person B, track 2:
x = c(3,3,4,3,3,4)
y = c(5,4,3,2,1,1)
n = length(x)
t = max(t) + cumsum(runif(n) * 60)
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, data.frame(co2 = rnorm(n)))
B2 = Track(stidf)

# Tracks for person B:
B = Tracks(list(B1=B1,B2=B2))

Tr = TracksCollection(list(A=A,B=B))

# approxTrack:
opar = par()
par(mfrow = c(1, 2))
plot(B2, ylim = c(.5, 6))
plot(B2, pch = 16, add = TRUE)
title("irregular time steps")
i = index(B2)
B3 = approxTrack(B2, seq(min(i), max(i), length.out = 50))
plot(B3, col = 'red', type = 'p', add = TRUE)

B4 = approxTrack(B2, seq(min(i), max(i), length.out = 50), FUN = spline)
plot(B4, col = 'blue', type = 'b', add = TRUE)

# regular time steps:
t = max(t) + (1:n) * 60 # regular
B2 = Track(STIDF(SpatialPoints(cbind(x,y),crs), t, data.frame(co2 = rnorm(n))))
plot(B2, ylim = c(.5, 6))
plot(B2, pch = 16, add = TRUE)
title("constant time steps")
i = index(B2)
B3 = approxTrack(B2)
plot(B3, type = 'p', col = 'red', add = TRUE)
B4 = approxTrack(B2, FUN = spline)
plot(B4, type = 'p', col = 'blue', add = TRUE)
# par(opar) # good to do, but would generate warnings
smth = function(x,y,xout,...) predict(smooth.spline(as.numeric(x), y), as.numeric(xout))
data(storms)
plot(storms, type = 'p')
storms.smooth = approxTracksCollection(storms, FUN = smth, n = 200)
plot(storms.smooth, add = TRUE, col = 'red')

}
