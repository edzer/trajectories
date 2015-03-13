\name{downsample}

\docType{methods}

\alias{downsample}
\alias{downsample,Track-method}
\title{Downsample a \code{Track}}

\description{
Downsamples a \code{Track} to the size (amount of points) of another \code{Track}.
}

\usage{%
\S4method{downsample}{Track}(track1, track2)
}

\arguments{
\item{track1}{\code{Track} that will be downsampled.}
\item{track2}{Reference \code{Track}.}
}

\value{A \code{Track} object. The downsampled track1.}

\author{
Nikolai Gorte <n.gorte@gmail.com>
}

\keyword{downsample}

