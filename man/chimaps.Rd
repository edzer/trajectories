\name{chimaps}
\alias{chimaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chimaps of tarjectory pattern.
}
\description{
Computes the chimaps, corresponding to a list of objects of class "Track". chimaps are based on the discrepancy between computed and expected intensity in a given location.}
\usage{
chimaps(X,timestamp,rank,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{
A list of Track objects
}
\item{timestamp}{based on secs,mins,...}
\item{rank}{a number between one and the length of corresponding time sequance which is created based on given timestamp.}
\item{...}{passed to arguments of density.Track}
}
\details{
[estimated intensity - expected intensity] / sqrt(expected intensity).
}
\value{
an image of class "im".
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{density.list}, \link[spatstat]{density.ppp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
chimaps(X, timestamp = "180 secs",rank = 2)
}
