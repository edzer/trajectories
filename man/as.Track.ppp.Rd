\name{as.Track.ppp}
\alias{as.Track.ppp}
\title{
Convert a trajectory pattern to a list of spatial point patterns
}

\description{
This function convert objects of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects to a list of spatial point patterns of class \link[spatstat]{ppp}.
}

\usage{
as.Track.ppp(X,timestamp)
}

\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs, mins,...}
}

\details{
The function \link[trajectories]{as.Track.ppp} converts objects of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects to a list of spatial point patterns of class \link[spatstat]{ppp} according to a given timestamp. Each obtained spatial point pattern represents the locatoin of all "moving objects" (e.g., taxis) in a fixed time.

Note that, tracks in X should overlap in time.
}

\value{
A list of spatial point patterns as objects of class \link[spatstat]{ppp}.
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{avedistTrack}, \link[spatstat]{as.ppp}
}

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Y <- as.Track.ppp(X,timestamp="120 secs")
}