\name{avemove}
\alias{avemove}
\title{
Average length of movement of a trajectory pattern over time
}

\description{
This returns the average length of movements of an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects over time.
}

\usage{
avemove(X,timestamp,epsilon=0)
}

\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{a timestamp to calculate the length of movements accordingly.}
\item{epsilon}{(optional) movements with a length less than epsilon are not considered in the calculation.}
}

\details{

The function first converts X to a list of spatial point patterns based on a given timestamp and using the function \link[trajectories]{as.Track.ppp}. Then each two consecutive (according to their corresponding time) point patterns make a segment pattern of class \link[spatstat]{psp}.
Each segment shows the movement, connecting an origin point (the location at the previous time) to a destination point (the location at the current time).
Finally, for the obtained segment patterns, the length of each segment is calculated using \link[spatstat]{lengths.psp} and their average is considered as "average length of movement" at that time.

Dealing with real data, sometimes moving objects such as taxies and pedestrians are not moving while being recorded their location. In order to deal with such a situation, argument epsilon is implemented in \link[trajectories]{avemove} to ignore the movements with a length less than epsilon. epsilon=0 does not ignore any movement.

Note that, tracks in X should overlap in time.

}

\value{
a numeric vector represents the average length of movements that can be plotted with respect to time.
}
\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{as.Track.arrow}, \link[spatstat]{lengths.psp}
}

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
<<<<<<< HEAD
out <- avemove(X,timestamp = "30 secs")
plot(out,type="l")
}
=======
avemove(X,timestamp = "30 secs")
}
>>>>>>> refs/remotes/edzer/master
