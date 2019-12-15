\name{frechetDist}
\docType{methods}
\alias{frechetDist}
\alias{frechetDist,Track-method}
\title{Frechet distance}
\description{
Compute the discrete Frechet distance between two \link[trajectories]{Track} objects.
}

\usage{%
\S4method{frechetDist}{Track}(track1, track2)
}

\arguments{
\item{track1}{an object of class Track.}
\item{track2}{an object of class Track.}
}


\value{A numeric value representing the Frechet distance.}

\references{http://en.wikipedia.org/wiki/Fr\'echet_distance}

\author{
Nikolai Gorte <n.gorte@gmail.com>

Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}

}

\examples{
x <- rTrack(10)
y <- rTrack(10)
frechetDist(x,y)
}
\keyword{methods}

