\name{rngTrack}
\alias{rngTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rngTrack
}
\description{
rngTrack returns the timerange of an object of class "Track".
}
\usage{
rngTrack(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ an object of class Track
%%     ~~Describe \code{x} here~~
}
}
\details{
This calculates the range of an object of class Track
}
\value{
an object of "POSIXct" "POSIXt"
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mohammad Mehdi moradi <moradi@uji.es>}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
rTrack}
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
rngTrack(A1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
