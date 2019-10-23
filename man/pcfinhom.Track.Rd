\name{pcfinhom.Track}
\alias{pcfinhom.Track}
\title{
Pair correlation funcrion of a trajectory pattern
}

\description{
Estimate the variability area of the pair correlation funcrion for an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects over time.
}
\usage{
pcfinhom.Track(X,timestamp,correction = c("translate", "Ripley"),q=0,
                sigma=c("default","bw.diggle","bw.ppl","bw.scott"),...)
}
\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs,mins,... .}
\item{correction}{the type of correction to be used when computing pair correlation function.}
\item{q}{(optional) a numeric value between 0 and 1. quantile to be applied to calculate the variability area.}
\item{sigma}{method to be considered for bandwidth selection in computing intensity function.}
\item{...}{arguments passed to the \link[spatstat]{pcfinhom}.}
}
\details{
The type of interaction (clustering, repulsion, randomness) between tracks in X might vary over time. The function \link[trajectories]{pcfinhom.Track} calculates the variability area of \link[spatstat]{pcfinhom} over time. The variability area is calculated based on the pair correlation functions of all obtained point patterns when converting a trajectory pattern to a list of point patterns. 

If sigma=default, it calculates the variability area using the defaults of \link[spatstat]{pcfinhom}, otherwise it first estimates the intensity function using the given bandwidth selection method and then using the estimated intensity function it estimates the variability area.

The lower and upper border of the variability area is calculated based on the \link[stats]{quantile} of values of all obtained pair correlation functions. q=0 means considering maximum and minimum.

Note that, tracks in X should overlap in time.
}
\value{
an object of class gTrack.
}
\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{as.Track.ppp}, \link[spatstat]{pcfinhom}}
\examples{
X <- list()
for(i in 1:100){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
g <- pcfinhom.Track(X,timestamp = "180 sec")
plot(g)
<<<<<<< HEAD
}
=======
}
>>>>>>> refs/remotes/edzer/master
