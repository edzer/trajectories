\name{tsqTracks}
\alias{tsqTracks}
\title{
sequance generation of time of a trajectory pattern
}
\description{
sequance generation of time of a trajectory pattern based on a given timestamp.
}
\usage{
tsqTracks(X,timestamp)
}

\arguments{
  \item{X}{an object of classes \link[trajectories]{Track}, \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{a timestamp to generate the time sequence based on it.}
}
\details{
This function generates a sequence of time based on a given timestamp for an object of classes \link[trajectories]{Track}, \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}

\value{
An object of class \link[base]{POSIXct} or \link[base]{POSIXt}.
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>}

\seealso{
\link[trajectories]{rTrack}, \link[base]{seq}
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
tsqTracks(A1,timestamp = "1 sec")
<<<<<<< HEAD
}
=======
}
>>>>>>> refs/remotes/edzer/master
