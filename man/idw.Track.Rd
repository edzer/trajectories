\name{idw.Track}
\alias{idw.Track}
\title{
Movement smoothing of a trajectory pattern
}
\description{
Movement smoothing of an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}
\usage{
\method{idw}{Track}(X,timestamp,epsilon=0,Fun=mean,...)
}

\arguments{
  \item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs,mins, ... .}
\item{epsilon}{(optional) movements with a length less than epsilon are not considered in the calculation.}
\item{Fun}{a function to be applied to the values of each pixel, see details}
\item{...}{arguments passed to \link[spatstat]{idw} function.}
}
\details{
Performs spatial smoothing to an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.

X will be converted to a list of marked point patterns using \link[trajectories]{as.Track.arrow} where marks are the length of movements at corresponding points, and then inverse-distance weighted smoothing will be applied to each obtained marked point pattern using \link[spatstat]{idw}. This results in a list of images in which function FUN (default mean) will then be applied to the corresponding values of each pixel.

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
idw.Track(X,timestamp="180 secs")
}
