\name{Kinhom.Track}
\alias{Kinhom.Track}
\title{
Inhomogeneous K-function for a trajectory pattern
}

\description{
Estimate the variability area of the K-function for an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects over time.}

\usage{
Kinhom.Track(X,timestamp,
                correction=c("border", "bord.modif", "isotropic", "translate"),q=0,
                sigma=c("default","bw.diggle","bw.ppl"," bw.scott"),...)
}

\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs,mins,... .}
\item{correction}{the type of correction to be used when computing K-function.}
\item{q}{(optional) a numeric value between 0 and 1. quantile to be applied to calculate the variability area.}
\item{sigma}{method to be considered for bandwidth selection in computing intensity function.}
\item{...}{arguments passed to the \link[spatstat]{Kinhom}.}
}

\details{
The type of interaction (clustering, repulsion, randomness) between tracks in X might vary over time. The function \link[trajectories]{Kinhom.Track} calculates the variability area of \link[spatstat]{Kinhom} over time. The variability area is calculated based on the K-functions of all obtained point patterns when converting a trajectory pattern to a list of point patterns. 

If sigma=default, it calculates the variability area using the defaults of \link[spatstat]{Kinhom}, otherwise it first estimates the intensity function using the given bandwidth selection method and then using the estimated intensity function it estimates the variability area.

The lower and upper border of the variability area is calculated based on the \link[stats]{quantile} of values of all obtained K-functions. q=0 means considering maximum and minimum.

Note that, tracks in X should overlap in time.
}

\value{
an object of class KTrack.
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{as.Track.ppp}, \link[spatstat]{Kinhom}
}

\examples{
library(spatstat)
X <- list()
for(i in 1:100){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
<<<<<<< HEAD
k <- Kinhom.Track(X,timestamp = "180 secs")
plot(k)
}
=======
Kinhom.Track(X,timestamp = "180 secs")

}
>>>>>>> refs/remotes/edzer/master
