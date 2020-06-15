\name{reTrack}
\alias{reTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct objects of class "Track"
}
\description{
Function reTrack accepts X as an object of class "Track". Output is a reconstructed Track (again an object of class Track), based on a regular "timestamp". It only returns the interpolated points.
}
\usage{
reTrack(X,at=c("track","dfrm"),timestamp=timestamp,tsq=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
an object of class Track
}
\item{at}{
to set the type of output as either an object of class "Track" or data.frame
}
\item{timestamp}{
timestamp which Track be reconstructed based on
}
\item{tsq}{
a time sequence to reconstruct Track X based on it. This is optional. If this is not given, the function creates the time sequance based on timestamp.
}
}
\details{
Sometimes tracks data are not collected according to a regular timestamp. In order to compare different tracks which share some time intervals, we might need to be aware of the locations in a regular timestamp. Function reTrack unables us to reconstruct an object of class "Track" based on a regular timestamp. Time sequance can be given by user, if not reTrack creates a regulare time sequance based on the given timestamp. 
}
\value{
Either an object of class "Track" or a data.frame
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{rTrack}, \link{as.Track},
 \link{as.POSIXct}, \link{compare}
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
crs = CRS("+proj=longlat +datum=WGS84") # longlat
stidf = STIDF(SpatialPoints(cbind(x,y),crs), t, data.frame(co2 = rnorm(n)))
A1 = Track(stidf)
reTrack(A1,timestamp = "1 sec")

}
