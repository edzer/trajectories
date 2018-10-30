\name{density.list}
\alias{density.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel estimate of intensity of trajectory pattern
}
\description{
Estimating the intensity of a list of tracks.
}
\usage{
\method{density}{list}(X,..., timestamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a list of "Track" objects, an object of class "Tracks" or "TracksCollection"
}
\item{timestamp}{based on secs, mins, ...}
\item{...}{passed to arguments of density.ppp}
}
\details{
This estimate the average intensity function of moving objects over time. Bandwidth selection methods such as bw.diggle, bw.scott and bw.ppl can be passed to this density.list.
}
\value{
an image of class "im".
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[trajectories]{rTrack}, \link[spatstat]{density.ppp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
density(X, timestamp = "180 secs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
