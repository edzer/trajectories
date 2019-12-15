\name{unique.Track}
\alias{unique.Track}
\title{
extract unique points from a trajectory pattern
}
\description{
Removing duplicated points in a \link[trajectories]{Track}.}
\usage{
\method{unique}{Track}(x,...)
}

\arguments{
\item{x}{
an object of class \link[trajectories]{Track}.
}
\item{...}{
arguments passed to \link[base]{unique}.
}
}

\details{
this function removes duplicated points in an object of class \link[trajectories]{Track}.
}
\value{
an object of class \link[trajectories]{Track} with  no duplicated point.
}
\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{unique}
}
\examples{
x <-  rTrack()
unique(x)
}
