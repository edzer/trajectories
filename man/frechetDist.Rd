\name{frechetDist}
\docType{methods}
\alias{frechetDist}
\alias{frechetDist,Track-method}
\title{Frechet distance}
\description{
Compute the discrete Frechet distance between two \code{Track} objects.
}

\usage{%
\S4method{frechetDist}{Track}(track1, track2)
}

\arguments{
\item{track1}{An object of class \code{Track}.}
\item{track2}{An object of class \code{Track}.}
}


\value{A numeric value representing the Frechet distance.}

\references{ http://en.wikipedia.org/wiki/Fr\'echet_distance }

\author{
Nikolai Gorte <n.gorte@gmail.com>
}

\examples{
x <- rTrack()
y <- rTrack()
frechetDist(x,y)
}
\keyword{methods}

