\name{as.Track.arrow}
\alias{as.Track.arrow}
\title{
Convert a trajectory pattern to a list of marked point patterns
}

\description{
Convert objects of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects to a list of marked point patterns.
}

\usage{
as.Track.arrow(X,timestamp,epsilon=0)
}

\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs, mins,...}
\item{epsilon}{optional. movements with a length less than epsilon are not considered in the calculation.}
}
\details{
This function converts an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects to a list of marked point patterns. 
The function first converts X to a list of spatial point patterns based on a given timestamp and using the function \link[trajectories]{as.Track.ppp}. Then each two consecutive (according to their corresponding time) point patterns make a segment pattern of class \link[spatstat]{psp}. Finally, each segment pattern will be converted to a marked point pattern according to which points represent the midpoints of segments and marks are the length of segments. The segment patterns are kept as attributes.

Dealing with real data, sometimes moving objects such as taxies and pedestrians are not moving while being recorded their location. In order to deal with such a situation, argument epsilon is implemented in \link[trajectories]{as.Track.arrow} to ignore the movements with a length less than epsilon. epsilon=0 does not ignore any movement.

Note that, tracks in X should overlap in time.
}
\value{
a list of marked point patterns.
}
\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{as.Track.ppp}, \link[spatstat]{psp}, \link[spatstat]{midpoints.psp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Y <- as.Track.arrow(X,timestamp="120 secs")
}
