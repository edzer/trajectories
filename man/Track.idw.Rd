\name{Track.idw}
\alias{Track.idw}
\title{
Movement smoothing of a trajectory pattern
}
\description{
Movement smoothing of an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}
\usage{
Track.idw(X,timestamp,epsilon=0,...)
}

\arguments{
  \item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs,mins, ... .}
\item{epsilon}{(optional) movements with a length less than epsilon are not considered in the calculation.}
\item{...}{arguments passed to \link[spatstat]{idw} function.}
}
\details{
Performs spatial smoothing to an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.

X will be converted to a list of marked point patterns using \link[trajectories]{as.Track.arrow} and then inverse-distance weighted smoothing will be applied to each obtained marked point pattern using \link[spatstat]{idw}. The average of all smoothed maps is considered as an inverse-distance weighted smoothing of X.

Note that, tracks in X should overlap in time.
}
\value{
an image of class \link[spatstat]{im}.
}
\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{as.Track.arrow}, \link[spatstat]{idw}
}

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Track.idw(X,timestamp="180 secs")
}
