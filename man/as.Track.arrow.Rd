\name{as.Track.arrow}
\alias{as.Track.arrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert trajectory pattern to a list of marked point patterns
}
\description{
Converting a list of Track objects to a list of marked point patterns. Each mark shows the length of movement.
}
\usage{
as.Track.arrow(X,timestamp,epsilon=epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A list of Track objects
}
\item{timestamp}{based on secs, mins,...}
\item{epsilon}{(optional) movements with length less than epsilon are not considered in the calculation}
}
\details{
Converting a list of Track objetcs to a list of marked point patterns. Marks show the length of movement with respect to the previous location.
}
\value{
a list of marked point patterns.
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{rTrack}, \link{as.Track.ppp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
Y <- as.Track.arrow(X,timestamp="120 secs")
}
