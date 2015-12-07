\name{rTrack}

\alias{rTrack}
\alias{rTracks}
\alias{rTracksCollection}

\title{Generate random \code{Track}, \code{Tracks} or \code{TracksCollection} objects}
\description{Generate random \code{Track}, \code{Tracks} or \code{TracksCollection} objects}
\usage{
rTrack(n = 100, origin = c(0,0), start = as.POSIXct("1970-01-01"), ar = .8, step = 60, sd = 1, ...)
rTracks(m = 20, start = as.POSIXct("1970-01-01"), delta = 7200, sdo = 0, ...)
rTracksCollection(p = 10, ...)
}
\arguments{
\item{n}{number of points per Track}
\item{origin}{numeric, length two, indicating the origin of the Track}
\item{start}{POSIXct, indicating the start time of the Track}
\item{ar}{numeric vector, indicating the amound of correlation in the Track}
\item{step}{ numeric; time step(s) in seconds between Track fixes }
\item{sd}{standard deviation of the random steps}
\item{...}{rTrack: arguments passed on to \link[stats]{arima.sim}, rTracks: arguments
passed on to rTrack; rTracksCollection: arguments passed on to rTracks}
\item{m}{ number of Track objects to simulate}
\item{delta}{ time difference between consecutive Track start times}
\item{sdo}{ standard deviation of the consecutive Track origin values (using rnorm)}
\item{p}{ number of IDs with Tracks to generate}
}
\note{\code{ar} is passed on to \link[stats]{arima.sim} as \code{ar}
element, so may contain multiple AR coefficients. The generated
track is a \link[base]{cumsum} over the simulated AR values.
In case it has length 1 and value 0, random walk is created using \link[stats]{rnorm}.}
\value{An object of class \code{Track}, \code{Tracks} or \code{TracksCollection}.}
\keyword{random}
