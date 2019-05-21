\name{reTrack}
\alias{reTrack}
\title{
Reconstruct objects of class Track.
}
\description{
Reconstruct objects of class \link[trajectories]{Track} according to a given timestamp.
}
\usage{
reTrack(X,at=c("track","dfrm"),timestamp=timestamp,tsq=NULL)
}
\arguments{
\item{X}{an object of class \link[trajectories]{Track}.}
\item{at}{to set the type of output as either an object of class \link[trajectories]{Track} or \link[base]{data.frame}.}
\item{timestamp}{a timestamp which \link[trajectories]{Track} be reconstructed based on.}
\item{tsq}{(optional) a time sequence to reconstruct X based on it. If this is not given, the function creates the time sequance based on timestamp.}
}

\details{
Sometimes tracks data are not collected according to a regular timestamp. In order to compare different tracks which share some time intervals, we might need to be aware of the locations in a regular timestamp. The function \link[trajectories]{reTrack} unables us to reconstruct an object of class \link[trajectories]{Track} based on a regular timestamp.

Time sequance can be given by user through argument tsq, if not \link[trajectories]{reTrack} first creates a regulare time sequance based on the given timestamp. 
}
\value{
an object of class \link[trajectories]{Track} or \link[base]{data.frame}.
}
\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{as.Track},
 \link[base]{as.POSIXct}, \link[trajectories]{compare}
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
reTrack(A1,timestamp = "1 sec")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
