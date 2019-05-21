\name{as.list.Tracks}
\alias{as.list.Tracks}
\title{
Tracks methods
}
\description{
Convert an object of class \link[trajectories]{Tracks} to a list of tracks.
}

\usage{
\method{as.list}{Tracks}(x,...)
}

\arguments{
\item{x}{an object of class "Tracks".}
\item{...}{arguments passed on the \link[base]{as.list} fucntion.}
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{as.list}
}

\examples{
x <-  rTracks()
as.list(x)
}
