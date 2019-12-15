\name{tsqTracks}
\alias{tsqTracks}
\title{
sequance generation of time of a trajectory pattern
}
\description{
sequance generation of time of a trajectory pattern based on a given timestamp.
}
\usage{
tsqTracks(X,timestamp)
}

\arguments{
\item{X}{
an object of classes \link[trajectories]{Track}, \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}
\item{timestamp}{
a given timestamp based on secs,mins, etc to generate the time sequence based on it.
}
}
\details{
this function generates a regular sequence of time based on a given timestamp for an object of classes \link[trajectories]{Track}, \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}

\value{
an object of class \link[base]{POSIXct} or \link[base]{POSIXt}.
}

\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}
\seealso{
\link[trajectories]{rTrack}, \link[base]{seq}
}

\examples{
X <- rTrack()
tsqTracks(X,timestamp = "120 sec")
}