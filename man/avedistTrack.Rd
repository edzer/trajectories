\name{avedistTrack}
\alias{avedistTrack}

\title{
Average pairwise distance of a trajectory pattern over time
}

\description{
This measures the average of pairwise distances between tracks over time.
}

\usage{
avedistTrack(X,timestamp)
}

\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{a timestamp to calculate the pairwise distances between tracks accordingly.}
}

\details{
This function calculates the average pairwise distance for an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects according to a given timestamp.

The function first converts X to a list of spatial point patterns based on a given timestamp and using the function \link[trajectories]{as.Track.ppp}. For each obtained point pattern pairwise distances are calculated using \link[spatstat]{pairdist} and then their average is considered as "average pairwise distance" for their corresponding time. Distance from each point to itself is ignored.

Note that, tracks in X should overlap in time.
}

\value{
a numeric vector represents the average pairwise distances that can be plotted with respect to time.
}

\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}

\seealso{
\link[trajectories]{as.Track.ppp}, \link[spatstat]{pairdist}
}

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
out <- avedistTrack(X,timestamp = "120 secs")
plot(out,type="l")
}