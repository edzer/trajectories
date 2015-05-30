\name{storms}
\alias{storms}
\title{Storm trajectories}
\description{storm trajectories, 2009-2012, from http://weather.unisys.com/hurricane/atlantic/}
\usage{
data(storms)
}
\keyword{datasets}
\examples{
data(storms)
dim(storms)
plot(storms)
x = approxTracksCollection(storms, by = "30 min", FUN = spline)
plot(x, col = 'red', add = TRUE)
\dontrun{
demo(storms) # regenerates these data from their source
}
}
