\name{compare}

\docType{methods}

\alias{compare}
\alias{compare,Track-method}

\title{Compares objects of class \code{Track}}

\description{
Calculates distances between two trajectories for the overlapping time interval.   
}

\usage{%
\S4method{compare}{Track}(track1, track2)
}

\arguments{
\item{track1}{An object of class \code{Track}.}
\item{track2}{An object of class \code{Track}.}
}

\value{A difftrack object. Includes both tracks extended with additional points 
for the timestamps of the other track. Also includes SpatialLines representing
the distances between the tracks.}

\examples{
require(enviroCaR)
require(xts)
tracks <- importSingleTrack("https://envirocar.org/api/stable", "53f8c09be4b04c314e7ce221")
track1 <- tracks[1]
track2 <- track1
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

difftrack <- compare(track1, track2)
plot(difftrack)
}

\author{
Nikolai Gorte <n.gorte@gmail.com>
}

\keyword{compare}

