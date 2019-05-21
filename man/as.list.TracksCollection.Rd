\name{as.list.TracksCollection}
\alias{as.list.TracksCollection}
\title{
TracksCollection methods
}
\description{
Convert an object of class \link[trajectories]{TracksCollection} to a list of tracks.
}

\usage{
\method{as.list}{TracksCollection}(x,...)
}

\arguments{
\item{x}{an object of class \link[trajectories]{TracksCollection}}
\item{...}{arguments passed on the \link[base]{as.list} fucntion.}
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{as.list}
}

\examples{
x <-  rTracksCollection()
as.list(x)
}
