\name{rTrack}

\alias{rTrack}
\alias{rTracks}
\alias{rTracksCollection}

\title{Generate random \code{Track}, \code{Tracks} or \code{TracksCollection} objects}
\description{Generate random \code{Track}, \code{Tracks} or \code{TracksCollection} objects}
\usage{
rTrack(n = 100, origin = c(0,0), start = as.POSIXct("1970-01-01"), ar = .8, 
	step = 60, sd0 = 1, bbox = bbox, transform = FALSE, nrandom = FALSE, ...)
rTracks(m = 20, start = as.POSIXct("1970-01-01"), delta = 7200, sd1 = 0, 
	origin = c(0,0), ...)
rTracksCollection(p = 10, sd2 = 0, ...)
}
\arguments{
\item{n}{number of points per Track}
\item{origin}{numeric, length two, indicating the origin of the Track}
\item{start}{POSIXct, indicating the start time of the Track}
\item{ar}{numeric vector, indicating the amound of correlation in the Track}
\item{step}{ numeric; time step(s) in seconds between Track fixes }
\item{sd0}{standard deviation of the random steps in a Track}
\item{sd1}{standard deviation of the consecutive Track origin values (using rnorm)}
\item{sd2}{standard deviation of the consecutive Tracks origin values (using rnorm)}
\item{bbox}{bbox object FIXME:fill in}
\item{transform}{logical; FIXME:fill in }
\item{nrandom}{logical; if \code{TRUE}, draw \code{n} from \code{rpois(n)}}
\item{...}{rTrack: arguments passed on to \link{arima.sim}, rTracks: arguments
passed on to rTrack; rTracksCollection: arguments passed on to rTracks}
\item{m}{ number of Track objects to simulate}
\item{delta}{ time difference between consecutive Track start times}
\item{p}{ number of IDs with Tracks to generate}
}
\details{\code{ar} is passed on to \link{arima.sim} as \code{ar}
element, and may contain multiple AR coefficients. The generated
track is a \link{cumsum} over the simulated AR values, for each dimension.
In case it has length 1 and value 0, random walk is created using \link{rnorm}. If bbox is given, the generated track will be transformed to bbox. If transform is TRUE and no bbox is given, it transforms the track to a unit box. If nrandom is TRUE, it generates a random number using \link{rpois} with parameter n as the number of locations per track.}
\value{An object of class \code{Track}, \code{Tracks} or \code{TracksCollection}.}
\author{
Edzer Pebesma <edzer.pebesma@uni-muenster.de>,
Mohammad Mehdi Moradi <moradi@uji.es>
}
\examples{
x = rTrack()
dim(x)
plot(x)
# x = rTracks(sd1 = 120)
# dim(x)
# plot(as(x, "SpatialLines"), col = 1:dim(x)[1], axes=TRUE)
# x = rTracksCollection() # star
# dim(x)
# plot(x)
x = rTracksCollection(sd2 = 200,p=4,m=10)
plot(x, col=1:dim(x)[1])
}
\keyword{random}
