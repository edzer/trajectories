\name{chimaps}
\alias{chimaps}
\title{
Chi-squared maps of a tarjectory pattern
}
\description{
Chi-squared maps of a tarjectory pattern.}
\usage{
chimaps(X,timestamp,rank=1,...)
}

\arguments{
\item{X}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.}
\item{timestamp}{based on secs, mins, ... .}
\item{rank}{a number between one and the length of corresponding time sequance which is created based on a given timestamp.}
\item{...}{arguments passed to \link[trajectories]{density.list}.}
}

\details{
This function converts a trajectory pattern to a list of spatial point patterns according to a given timestamp and calculates their intensities using \link[trajectories]{density.list}. Based on the estimated intensities, it calculates an expected intensity for a sequence of time which is already obtained based on a given timestamp. Having estimated and expected intensities, \link[trajectories]{chimaps} computes Chi-squared maps for a sequence of time as follows

\eqn{\chi=[\lambda_o - \lambda_e]/\sqrt\lambda_e}

where \eqn{\lambda_o} and \eqn{\lambda_e} are the estimated intensity and expected intensity respectively.

According to a given timestamp, a time sequence will be created based on the range of times of tracks so that Chi-squared maps will be computed for this time sequence. Argument rank means which one to be represented as output, however, all Chi-squared maps are still stored as attributes.

Note that, tracks in X should overlap in time.
}
 
\value{
an image of class \link[spatstat]{im}.
}

\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}

\seealso{
\link[trajectories]{density.list}, \link[spatstat]{density.ppp}
}

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
chimaps(X, timestamp = "180 secs",rank = 2)
}
