\name{tsqTracks}
\alias{tsqTracks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tsqTracks
}
\description{
tsqtracks returns a sequance of time based on a list of tracks (or a single object of class "Track"") and an argument timestamp.

}
\usage{
tsqTracks(X,timestamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{either an object of class "Track"" or a list of some objects of class "Track"
%%     ~~Describe \code{x} here~~
}
\item{timestamp}{a timestamp to create the time sequence based on it}
}
\details{
This creates a sequence of time based on a track or a list of tracks.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class "POSIXct" or "POSIXt".
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
rTrack
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
tsqTracks(A1,timestamp = "1 sec")
}
