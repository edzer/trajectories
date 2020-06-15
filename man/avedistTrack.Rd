\name{avedistTrack}
\alias{avedistTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average pairwise distance of trajectory pattern over time
}
\description{
This measures the average of pairwise distances between tracks over time.
}
\usage{
avedistTrack(X,timestamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of some objects of class "Track"
%%     ~~Describe \code{x} here~~
}
\item{timestamp}{timestamp to calculate the pairwise distances between tarcks}
}
\details{
This function calculates the average pairwise distance between a list of tracks according to a given timestamp.
}
\value{
An object of class "distrack". It can be plotted over time.
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{as.Track.ppp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}

ave <- avedistTrack(X,timestamp = "120 secs")
plot(ave,type="l")


}
