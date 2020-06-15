\name{Kinhom.Track}
\alias{Kinhom.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inhomogeneous K-function for trajectory pattern
}
\description{
Estimate the variability area of K-function of a list of tracks.}
\usage{
Kinhom.Track(X,timestamp,
                correction=c("border", "bord.modif", "isotropic", "translate"),q,
                sigma=c("default","bw.diggle","bw.ppl"," bw.scott"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{
A list of Track objects
}
\item{timestamp}{based on secs,mins,...}
\item{correction}{the type of correction to be used in computing K-function}
\item{q}{(optional) a numeric value between 0 and 1. quantile to be applied to calculate the variability area}
\item{sigma}{method to be used in computing intensity function}
\item{...}{passed to the arguments of Kinhom}
}
\details{
This calculates the variability area of K-function over time. If sigma=default, it calculates the variability area using the defaults of Kinhom, otherwise it first estimate the intensity function using the given sigma as bandwidth selection method and then using the estimated intensity function, it estimates the variability area.
}
\value{
an object of class "KTrack".
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{rTrack}, \link{as.Track.ppp}, \link[spatstat]{Kinhom}}
\examples{
library(spatstat)
X <- list()
for(i in 1:100){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Kinhom.Track(X,timestamp = "180 secs")

}
