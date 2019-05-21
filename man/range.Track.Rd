\name{range.Track}
\alias{range.Track}
\title{
range of a Track
}
\description{
Retrieves the range of an object of class \link[trajectories]{Track}.
}
\usage{
\method{range}{Track}(X,...)
}
\arguments{
\item{X}{an object of class Track.}
\item{...}{arguments passed to the \link[base]{range} function.}
}
\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{range}
}
\examples{
x <-  rTrack()
range(x)
}
