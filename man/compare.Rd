\name{compare}

\docType{methods}

\alias{compare}
\alias{compare,Track-method}

\title{Compares objects of class \code{Track}}

\description{
Calculates distances between two tracks for the overlapping time interval.   
}

\usage{%
\S4method{compare}{Track}(tr1, tr2)
}

\arguments{
\item{tr1}{An object of class \code{Track}.}
\item{tr2}{An object of class \code{Track}.}
}

\value{A difftrack object. Includes both tracks extended with additional points 
for the timestamps of the other track. Also includes SpatialLines representing
the distances between the tracks.}

\examples{
## example tracks
library(sp)
library(xts)
data(A3)
track2 <- A3
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

## compare and plot
difftrack <- compare(A3, track2)
plot(difftrack)
}

\author{
Nikolai Gorte <n.gorte@gmail.com>
}

\keyword{compare}

