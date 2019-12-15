\name{range.Track}
\alias{range.Track}
\alias{range.Tracks}
\alias{range.TracksCollection}
\title{
range
}
\description{
Retrieves the range of objects of classes \link[trajectories]{Track},\link[trajectories]{Tracks} and \link[trajectories]{TracksCollection}.
}
\usage{
\method{range}{Track}(X,...)
\method{range}{Tracks}(X,...)
\method{range}{TracksCollection}(X,...)
}
\arguments{
\item{X}{an object of class \link[trajectories]{Track},\link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.}
\item{...}{arguments passed to the \link[base]{range} function.}
}
\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{range}
}
\examples{
x <-  rTrack()
range(x)
}
