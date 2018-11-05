\name{dists
}
\docType{methods}

\alias{dists}
\alias{dists,Tracks,Tracks-method}

\title{Calculate distances between two \code{Tracks} objects}

\description{
Calculates a distance matrix with distances for each pair of tracks.
}

\usage{%
\S4method{dists}{Tracks,Tracks}(tr1, tr2, f, ...)
}

\arguments{
\item{tr1}{An object of class \code{Tracks}.}
\item{tr2}{An object of class \code{Tracks}.}
\item{f}{A function to calculate distances. Default is \code{mean}.}
\item{...}{Additional parameters passed to \code{f}.}
}

\value{A matrix with distances between each pair of tracks or \code{NA}
if they don't overlap in time.}

\details{
\code{f} can be any function applicable to a numerical vector or \link{frechetDist}.
}

\examples{
## example tracks
library(sp)
library(xts)
data(A3)
track2 <- A3
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

## create Tracks objects
tracks1 <- Tracks(list(A3, track2))
tracks2 <- Tracks(list(track2, A3))

## calculate distances
\dontrun{
dists(tracks1, tracks2)
dists(tracks1, tracks2, sum)
dists(tracks1, tracks2, frechetDist)
}
}

\keyword{dists}

