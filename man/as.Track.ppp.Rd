\name{as.Track.ppp}
\alias{as.Track.ppp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conver trajectory pattern to a list of objects of class ppp
}
\description{
This function converts a list of Tracks to a list of point patterns (class "ppp")}
\usage{
as.Track.ppp(X,timestamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a list of Track objects
}
\item{timestamp}{based on secs, mins,...}
}
\details{
as.Track.ppp converts a list of Track objetcs to a list of ppp objetcs.}
\value{
A list of point patterns, objects of class "ppp".
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{avedistTrack}, \link[spatstat]{as.ppp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Y <- as.Track.ppp(X,timestamp="120 secs")
}
