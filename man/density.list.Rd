\name{density.list}
\alias{density.list}
\title{
Kernel estimate of intensity of a trajectory pattern
}

\description{
Estimating the intensity of an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}

\usage{
\method{density}{list}(x, timestamp,FUN=mean,...)
}
\arguments{
\item{x}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs, mins, ... .}
\item{FUN}{a function to be applied to the values of each pixel, see details}
\item{...}{arguments passed to  the \link[spatstat]{density.ppp}.}
}

\details{
This function converts a trajectory pattern to a list of spatial point patterns and then it calclulates the intensity of each obtained point pattern using \link[spatstat]{density.ppp}. This results in a list of images in which function FUN (default mean) will then be applied to the corresponding values of each pixel. Bandwidth selection methods such as \link[spatstat]{bw.diggle}, \link[spatstat]{bw.scott} and \link[spatstat]{bw.ppl} can be considered in \link[trajectories]{density.list}.

Note that, tracks in X should overlap in time.
}

\value{
an image of class \link[spatstat]{im}.
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{rTrack}, \link[spatstat]{density.ppp}, \link[spatstat]{bw.diggle}, \link[spatstat]{bw.scott}, \link[spatstat]{bw.ppl}
}

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
density(X, timestamp = "180 secs")
}
