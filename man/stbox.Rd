\name{stbox}
\alias{stbox}
\alias{stbox,Tracks-method} 
\alias{stbox,TracksCollection-method} 
\title{
obtain ranges of space and time coordinates
}
\description{
obtain ranges of space and time coordinates of objects of classes \link[trajectories]{Track}, \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection}.
}


\usage{ 
stbox(obj)
}
\arguments{
\item{obj}{
object of a class deriving from \link[trajectories]{Track}, \link[trajectories]{Tracks} and \link[trajectories]{TracksCollection}.
} 
}
\value{
\code{stbox} returns a \code{data.frame}, with three columns
representing x-, y- and time-coordinates, and two rows containing min
and max values. \code{bbox} gives a matrix with coordinate min/max values,
compatible to \link[sp]{bbox}.
}

\examples{
x <- rTrack()
stbox(x)
}
\keyword{dplot}
