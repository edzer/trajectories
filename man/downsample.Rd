\name{downsample}
\docType{methods}
\alias{downsample}
\alias{downsample,Track-method}
\title{Downsample a Track}

\description{
Downsamples a \link[trajectories]{Track} to the size (amount of points) of another Track.
}

\usage{%
\S4method{downsample}{Track}(track1, track2)
}

\arguments{
\item{track1}{Track that will be downsampled.}
\item{track2}{Reference Track.}
}

\value{A Track object. The downsampled track1.}

\author{
Nikolai Gorte <n.gorte@gmail.com>
}

\examples{
x <- rTrack(15)
dim(x)
y <- rTrack(20)
dim(y)
z <- downsample(y,x)
dim(z)
}

\keyword{downsample}

