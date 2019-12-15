\name{generalize}

\alias{generalize}
\alias{generalize,Track-method} 
\alias{generalize,Tracks-method} 
\alias{generalize,TracksCollection-method}  

\title{Generalize trajectories}

\description{Generalize objects of class \link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.}

\usage{%
\S4method{generalize}{Track}(t, FUN = mean, ...,  timeInterval, distance, n, tol, toPoints)
\S4method{generalize}{Tracks}(t, FUN = mean, ...)
\S4method{generalize}{TracksCollection}(t, FUN = mean, ...)}

\arguments{
\item{t}{an object of class \link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.
}
\item{FUN}{
the generalization method to be applied. Defaults to mean if none is passed.
}
\item{timeInterval}{
(lower limit) time interval to split Track into segments.
}
\item{distance}{
(lower limit) distance to split Track into segments.
}
\item{n}{number of points to form segments.
}
\item{tol}{
tolerance passed on to \link[rgeos]{gSimplify}, to generalize segments using the Douglas-Peucker algorithm. 
}
\item{toPoints}{
keep mid point rather than forming \link[sp]{SpatialLines} segments.
}
\item{...}{additional arguments passed to FUN.
}
}

\value{A generalized object of class \link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.}


\details{
Generalize object of classes \link[trajectories]{Track}, \link[trajectories]{Traks} or \link[trajectories]{TracksCollection} to \link[sp]{SpatialLines} based on number of spatial points, time or distance passed.
}

\examples{
x <- rTrack()
generalize(x,n=15)
generalize(x,distance=40)
generalize(x,timeInterval=50)
}

\keyword{generalize}
