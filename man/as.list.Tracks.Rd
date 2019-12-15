\name{as.list.Tracks}
\alias{as.list.Tracks}
\alias{as.list.TracksCollection}
\title{
Tracks methods
}
\description{
Convert an objects of classes \link[trajectories]{Tracks} and \link[trajectories]{TracksCollection} to a list of of objects of class \link[trajectories]{Track}.
}

\usage{
\method{as.list}{Tracks}(tr,...)
\method{as.list}{TracksCollection}(tc,...)
}

\arguments{
\item{tr,tc}{objects of classes \link[trajectories]{Tracks} and \link[trajectories]{TracksCollection}.
}
\item{...}{optional. additional arguments passed on the \link[base]{as.list} fucntion.}
}

\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{as.list}
}

\examples{
tr <-  rTracks(m=2,n=10)
as.list(tr)
tc <-  rTracksCollection(p=2,m=2,n=10)
as.list(tc)
}
