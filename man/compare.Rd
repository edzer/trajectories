\name{compare}

\docType{methods}

\alias{compare}
\alias{compare,Track-method}
\alias{compare,Tracks-method}

\title{Compares objects of class \code{Track}, \code{Tracks}}

\description{
Calculates distances between two trajectories for the overlapping time interval.   
}

\usage{%
\S4method{compare}{Track}(tr1, tr2)
\S4method{compare}{Tracks}(tr1, tr2)
}

\arguments{
\item{tr1}{An object of class \code{Track} or \code{Tracks}.}
\item{tr2}{An object of class \code{Track} or \code{Tracks}.}
}

\value{A difftrack object. Includes both tracks extended with additional points 
for the timestamps of the other track. Also includes SpatialLines representing
the distances between the tracks.}

\examples{
## example tracks
library(enviroCaR)
library(sp)
library(xts)
tracks <- importSingleTrack("https://envirocar.org/api/stable", "53f8c09be4b04c314e7ce221")
track1 <- tracks[1]
track2 <- track1
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

## compare and plot
difftrack <- compare(track1, track2)
plot(difftrack)
}

\author{
Nikolai Gorte <n.gorte@gmail.com>
}

\keyword{compare}

