\name{Track.idw}
\alias{Track.idw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Movement smoothing of trajectory pattern
}
\description{
Movement smoothing of trajectory pattern
}
\usage{
Track.idw(X,timestamp,epsilon=epsilon,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a list of objects of class "Track"
}
\item{timestamp}{
based on secs,mins, ...
}
\item{epsilon}{(optional) movements with length less than epsilon are not considered in the calculation}
\item{...}{passed to arguments of fucntion idw in spatstat}
}
\details{
Performs spatial smoothing to the movements of a list of tracks.
}
\value{
an image of class "im".
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{as.Track.arrow}, \link[spatstat]{idw}}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Track.idw(X,timestamp="180 secs")
}
