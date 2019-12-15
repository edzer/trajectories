\name{dists
}
\docType{methods}

\alias{dists}
\alias{dists,Tracks,Tracks-method}

\title{Calculate matrix of distances}

\description{
Calculates distances between two \link[trajectories]{Tracks} objects.
}

\usage{%
\S4method{dists}{Tracks,Tracks}(tr1, tr2, f, ...)
}

\arguments{
\item{tr1}{an object of class \link[trajectories]{Tracks}.}
\item{tr2}{an object of class \link[trajectories]{Tracks}.}
\item{f}{a function to calculate distances. Default is \code{mean}.}
\item{...}{additional parameters passed to \code{f}.}
}

\value{A matrix with distances between each pair of tracks or \code{NA}
if they don't overlap in time.}

\details{
This function makes use of fucntion \link[trajectories]{compare}, and 
\code{f} can be any function applicable to a numerical vector or \link{frechetDist}.
}

\seealso{
\link[trajectories]{compare}, \link[trajectories]{frechetDist}, \link[trajectories]{Tracks}
}


\examples{
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
dists(tracks1, tracks2)
dists(tracks1, tracks2, sum)
dists(tracks1, tracks2, frechetDist)
}

\keyword{dists}

